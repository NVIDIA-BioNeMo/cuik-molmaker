// SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

//! @file CGR (Condensed Graph of Reaction) featurization.
//! Implements parse_rxn_side_mol, parse_reaction, reaction_mode_names_to_array,
//! and batch_reaction_featurizer (all declared in features.h).

#include <GraphMol/Atom.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/types.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "features.h"
#include "float_features.h"
#include "one_hot.h"

namespace py = pybind11;

static constexpr uint32_t NO_IDX = std::numeric_limits<uint32_t>::max();

// ---------------------------------------------------------------------------
// Internal: populate CompactAtom/Bond arrays from an already-parsed mol.
// Mirrors read_graph() in features.cpp but accepts a pre-built mol.
// ---------------------------------------------------------------------------
static void populate_graph_arrays(GraphData& gd) {
  const RDKit::ROMol& mol       = *gd.mol;
  const size_t        num_atoms = gd.num_atoms;
  const size_t        num_bonds = gd.num_bonds;

  gd.atoms = std::unique_ptr<CompactAtom[]>(new CompactAtom[num_atoms]);
  for (size_t i = 0; i < num_atoms; ++i) {
    const RDKit::Atom* a = mol.getAtomWithIdx(i);
    gd.atoms[i]          = CompactAtom{uint8_t(a->getAtomicNum()),
                              uint8_t(a->getTotalDegree()),
                              int8_t(a->getFormalCharge()),
                              uint8_t(a->getChiralTag()),
                              uint8_t(a->getTotalNumHs()),
                              uint8_t(a->getHybridization()),
                              a->getIsAromatic(),
                              float(a->getMass())};
  }

  gd.bonds                              = std::unique_ptr<CompactBond[]>(new CompactBond[num_bonds]);
  const RDKit::RingInfo* const ringInfo = mol.getRingInfo();
  for (size_t i = 0; i < num_bonds; ++i) {
    const RDKit::Bond* b = mol.getBondWithIdx(i);
    gd.bonds[i]          = CompactBond{uint8_t(b->getBondType()),
                              b->getIsConjugated(),
                              ringInfo->numBondRings(i) != 0,
                              uint8_t(b->getStereo()),
                              b->getBeginAtomIdx(),
                              b->getEndAtomIdx()};
  }
}

// ---------------------------------------------------------------------------
// Internal: bond lookup map — key = (min_atom_idx << 32) | max_atom_idx → bond_idx
// ---------------------------------------------------------------------------
static std::unordered_map<uint64_t, uint32_t> build_bond_lookup(const GraphData& gd) {
  std::unordered_map<uint64_t, uint32_t> lookup;
  lookup.reserve(gd.num_bonds);
  for (size_t i = 0; i < gd.num_bonds; ++i) {
    uint32_t a = gd.bonds[i].beginAtomIdx;
    uint32_t b = gd.bonds[i].endAtomIdx;
    if (a > b)
      std::swap(a, b);
    lookup[(uint64_t(a) << 32) | uint64_t(b)] = uint32_t(i);
  }
  return lookup;
}

// ---------------------------------------------------------------------------
// parse_rxn_side_mol
// ---------------------------------------------------------------------------
std::unique_ptr<RDKit::RWMol> parse_rxn_side_mol(const std::string& smiles, bool keep_h, bool add_h) {
  RDKit::SmilesParserParams params;
  params.removeHs = !keep_h;  // keep_h=true keeps explicit [H:n] atoms
  std::unique_ptr<RDKit::RWMol> mol{RDKit::SmilesToMol(smiles, params)};
  if (!mol)
    return mol;
  // Do NOT clearProp(molAtomMapNumber) — needed to build r2p_idx_map.
  // Do NOT reorder atoms — reactions preserve SMILES parse order.
  if (add_h)
    RDKit::MolOps::addHs(*mol);
  return mol;
}

// ---------------------------------------------------------------------------
// parse_reaction
// ---------------------------------------------------------------------------
CompactReaction parse_reaction(const std::string& reac_smi, const std::string& prod_smi, bool keep_h, bool add_h) {
  std::unique_ptr<RDKit::RWMol> reac_mol = parse_rxn_side_mol(reac_smi, keep_h, add_h);
  std::unique_ptr<RDKit::RWMol> prod_mol = parse_rxn_side_mol(prod_smi, keep_h, add_h);
  if (!reac_mol)
    throw std::runtime_error("Failed to parse reactant SMILES: " + reac_smi);
  if (!prod_mol)
    throw std::runtime_error("Failed to parse product SMILES: " + prod_smi);

  const size_t n_reac_atoms = reac_mol->getNumAtoms();
  const size_t n_prod_atoms = prod_mol->getNumAtoms();

  // Build mapno → product atom index
  std::unordered_map<int, uint32_t> mapno_to_prod_idx;
  mapno_to_prod_idx.reserve(n_prod_atoms);
  for (const auto* atom : prod_mol->atoms()) {
    int map_num = atom->getAtomMapNum();
    if (map_num > 0)
      mapno_to_prod_idx[map_num] = atom->getIdx();
  }

  // Build r2p / p2r maps; classify reactant atoms as matched or reac-only
  std::unordered_map<uint32_t, uint32_t> r2p_idx_map, p2r_idx_map;
  std::vector<uint32_t>                  reac_only_idxs;
  r2p_idx_map.reserve(n_reac_atoms);
  p2r_idx_map.reserve(n_reac_atoms);

  for (const auto* atom : reac_mol->atoms()) {
    uint32_t r_idx   = atom->getIdx();
    int      map_num = atom->getAtomMapNum();
    auto     it      = (map_num > 0) ? mapno_to_prod_idx.find(map_num) : mapno_to_prod_idx.end();
    if (it != mapno_to_prod_idx.end()) {
      r2p_idx_map[r_idx]      = it->second;
      p2r_idx_map[it->second] = r_idx;
    } else {
      reac_only_idxs.push_back(r_idx);
    }
  }

  // Classify product atoms: those without a reactant counterpart are product-only
  std::vector<uint32_t> prod_only_idxs;
  for (const auto* atom : prod_mol->atoms()) {
    if (p2r_idx_map.find(atom->getIdx()) == p2r_idx_map.end())
      prod_only_idxs.push_back(atom->getIdx());
  }

  // Build GraphData for each side (mol is moved in; arrays populated from it)
  GraphData reac_gd{n_reac_atoms,
                    std::unique_ptr<CompactAtom[]>(),
                    reac_mol->getNumBonds(),
                    std::unique_ptr<CompactBond[]>(),
                    std::move(reac_mol)};
  GraphData prod_gd{n_prod_atoms,
                    std::unique_ptr<CompactAtom[]>(),
                    prod_mol->getNumBonds(),
                    std::unique_ptr<CompactBond[]>(),
                    std::move(prod_mol)};

  populate_graph_arrays(reac_gd);
  populate_graph_arrays(prod_gd);

  auto reac_bond_lookup = build_bond_lookup(reac_gd);
  auto prod_bond_lookup = build_bond_lookup(prod_gd);

  return CompactReaction{std::move(reac_gd),
                         std::move(prod_gd),
                         std::move(r2p_idx_map),
                         std::move(p2r_idx_map),
                         std::move(reac_only_idxs),
                         std::move(prod_only_idxs),
                         std::move(reac_bond_lookup),
                         std::move(prod_bond_lookup)};
}

// ---------------------------------------------------------------------------
// reaction_mode_names_to_array
// ---------------------------------------------------------------------------
static const std::unordered_map<std::string, int64_t> rxn_mode_name_to_enum{
  {        std::string("REAC_DIFF"),         int64_t(ReactionMode::REAC_DIFF)},
  {        std::string("REAC_PROD"),         int64_t(ReactionMode::REAC_PROD)},
  {        std::string("PROD_DIFF"),         int64_t(ReactionMode::PROD_DIFF)},
  {std::string("REAC_DIFF_BALANCE"), int64_t(ReactionMode::REAC_DIFF_BALANCE)},
  {std::string("REAC_PROD_BALANCE"), int64_t(ReactionMode::REAC_PROD_BALANCE)},
  {std::string("PROD_DIFF_BALANCE"), int64_t(ReactionMode::PROD_DIFF_BALANCE)},
};

py::array_t<int64_t> reaction_mode_names_to_array(const std::vector<std::string>& modes) {
  const size_t               n = modes.size();
  std::unique_ptr<int64_t[]> out(new int64_t[n ? n : 1]);
  for (size_t i = 0; i < n; ++i) {
    auto it = rxn_mode_name_to_enum.find(modes[i]);
    out[i]  = (it != rxn_mode_name_to_enum.end()) ? it->second : int64_t(ReactionMode::UNKNOWN);
  }
  const int64_t dims[1] = {int64_t(n)};
  return py_array_from_array(std::move(out), dims, 1);
}

// ---------------------------------------------------------------------------
// Internal: fill all-atom feature temp array using existing all-atoms functions.
// Result: buf[atomIdx * single_fdim .. (atomIdx+1)*single_fdim - 1] = full feature vector for atomIdx.
// ---------------------------------------------------------------------------
static void fill_all_atom_features(const GraphData&            gd,
                                   std::vector<float>&         buf,
                                   size_t                      single_fdim,
                                   const py::array_t<int64_t>& atom_property_list_onehot,
                                   const py::array_t<int64_t>& atom_property_list_float,
                                   bool                        offset_carbon) {
  buf.assign(gd.num_atoms * single_fdim, 0.0f);
  if (gd.num_atoms == 0)
    return;

  const size_t   n_onehot = (atom_property_list_onehot.ndim() == 1) ? size_t(atom_property_list_onehot.shape(0)) : 0;
  const size_t   n_float  = (atom_property_list_float.ndim() == 1) ? size_t(atom_property_list_float.shape(0)) : 0;
  const int64_t* oh_ptr   = (n_onehot > 0) ? static_cast<const int64_t*>(atom_property_list_onehot.data()) : nullptr;
  const int64_t* fl_ptr   = (n_float > 0) ? static_cast<const int64_t*>(atom_property_list_float.data()) : nullptr;

  float* p = buf.data();
  for (size_t i = 0; i < n_onehot; ++i) {
    auto   feat = AtomOneHotFeature(oh_ptr[i]);
    size_t fsz  = get_one_hot_atom_feature<float>(gd, p, feat, single_fdim);
    p += fsz;
  }
  for (size_t i = 0; i < n_float; ++i) {
    get_atom_float_feature<float>(gd, p, AtomFloatFeature(fl_ptr[i]), single_fdim, offset_carbon);
    ++p;
  }
}

// ---------------------------------------------------------------------------
// Internal: build a num_only feature vector into out[0..single_fdim-1].
// All zeros except the atomic-number one-hot index for atomicNum, at offset 0.
// Mirrors chemprop's MultiHotAtomFeaturizer.num_only(atom).
// ---------------------------------------------------------------------------
static void build_num_only(uint8_t atomicNum, float* out, size_t single_fdim, AtomOneHotFeature first_onehot_feat) {
  std::fill(out, out + single_fdim, 0.0f);
  size_t idx = get_atomic_num_onehot_index(atomicNum, first_onehot_feat);
  out[idx]   = 1.0f;
}

// ---------------------------------------------------------------------------
// CGR atom feature computation (fills output for one reaction).
//
// CGR atom layout:
//   nodes 0..n_reac-1           = reactant atoms (in RDKit parse order)
//   nodes n_reac..n_cgr-1       = product-only atoms (in prod_only_idxs order)
//
// Feature layout per CGR node (cgr_atom_fdim = single_fdim + second_len):
//   [first_half: single_fdim values | second_half: second_len values]
//
// second_len = single_fdim - atomic_num_block_w
//   (strips the entire atomic-num one-hot block from the second side)
//
// Diff direction (for *_DIFF modes): ALWAYS prod_feats - reac_feats.
// ---------------------------------------------------------------------------
static void fill_cgr_atom_features(const CompactReaction&      rxn,
                                   ReactionMode                mode,
                                   float*                      cgr_out,
                                   size_t                      cgr_atom_fdim,
                                   size_t                      single_fdim,
                                   size_t                      atomic_num_block_w,
                                   const py::array_t<int64_t>& atom_property_list_onehot,
                                   const py::array_t<int64_t>& atom_property_list_float,
                                   bool                        offset_carbon) {
  const size_t n_reac     = rxn.reac.num_atoms;
  const size_t n_cgr      = n_reac + rxn.prod_only_idxs.size();
  const size_t second_len = single_fdim - atomic_num_block_w;

  // Mode flags
  const bool is_balance = (mode == ReactionMode::REAC_DIFF_BALANCE || mode == ReactionMode::REAC_PROD_BALANCE ||
                           mode == ReactionMode::PROD_DIFF_BALANCE);
  const bool prod_first = (mode == ReactionMode::PROD_DIFF || mode == ReactionMode::PROD_DIFF_BALANCE);
  const bool use_diff   = (mode == ReactionMode::REAC_DIFF || mode == ReactionMode::PROD_DIFF ||
                         mode == ReactionMode::REAC_DIFF_BALANCE || mode == ReactionMode::PROD_DIFF_BALANCE);

  // 1. Pre-compute full feature vectors for all reactant and product atoms.
  // reac_buf[r * single_fdim .. (r+1)*single_fdim-1] = full feature vector for reactant atom r.
  std::vector<float> reac_buf, prod_buf;
  fill_all_atom_features(rxn.reac,
                         reac_buf,
                         single_fdim,
                         atom_property_list_onehot,
                         atom_property_list_float,
                         offset_carbon);
  fill_all_atom_features(rxn.prod,
                         prod_buf,
                         single_fdim,
                         atom_property_list_onehot,
                         atom_property_list_float,
                         offset_carbon);

  // Determine the first onehot feature (needed for num_only).
  AtomOneHotFeature first_feat = AtomOneHotFeature::ATOMIC_NUM;
  if (atom_property_list_onehot.ndim() == 1 && atom_property_list_onehot.shape(0) > 0)
    first_feat = AtomOneHotFeature(static_cast<const int64_t*>(atom_property_list_onehot.data())[0]);

  // Scratch buffers for num_only vectors (reused each iteration).
  std::vector<float> num_only_reac(single_fdim, 0.0f);
  std::vector<float> num_only_prod(single_fdim, 0.0f);

  // 2. For each CGR node, mix reac and prod feature vectors.
  for (size_t u = 0; u < n_cgr; ++u) {
    float* atom_out = cgr_out + u * cgr_atom_fdim;

    // Determine reac_feats and prod_feats for this CGR node.
    const float* reac_feats = nullptr;
    const float* prod_feats = nullptr;

    if (u < n_reac) {
      // Reactant atom
      reac_feats    = reac_buf.data() + u * single_fdim;
      auto it       = rxn.r2p_idx_map.find(uint32_t(u));
      bool has_prod = (it != rxn.r2p_idx_map.end());
      if (has_prod) {
        prod_feats = prod_buf.data() + it->second * single_fdim;
      } else {
        // Reac-only: product side
        if (is_balance) {
          prod_feats = reac_feats;  // copy own feats
        } else {
          build_num_only(rxn.reac.atoms[u].atomicNum, num_only_prod.data(), single_fdim, first_feat);
          prod_feats = num_only_prod.data();
        }
      }
    } else {
      // Product-only atom
      uint32_t p_idx = rxn.prod_only_idxs[u - n_reac];
      prod_feats     = prod_buf.data() + p_idx * single_fdim;
      if (is_balance) {
        reac_feats = prod_feats;  // copy own feats
      } else {
        build_num_only(rxn.prod.atoms[p_idx].atomicNum, num_only_reac.data(), single_fdim, first_feat);
        reac_feats = num_only_reac.data();
      }
    }

    // Write first half: reac side for REAC_* modes, prod side for PROD_* modes.
    const float* first_src = prod_first ? prod_feats : reac_feats;
    std::copy(first_src, first_src + single_fdim, atom_out);

    // Write second half (length = second_len), starting at atom_out[single_fdim].
    // Strips atomic_num_block_w from the start of the second side.
    // Diff = prod - reac (always, regardless of which is "first").
    float* second_out = atom_out + single_fdim;
    if (use_diff) {
      for (size_t k = 0; k < second_len; ++k)
        second_out[k] = prod_feats[atomic_num_block_w + k] - reac_feats[atomic_num_block_w + k];
    } else {
      // REAC_PROD / REAC_PROD_BALANCE: second half = prod feats (atomic-num block stripped)
      std::copy(prod_feats + atomic_num_block_w, prod_feats + atomic_num_block_w + second_len, second_out);
    }
  }
}

// ---------------------------------------------------------------------------
// Internal: write one-hot + float features for a single bond into buf[0..single_bond_fdim-1].
// bond_idx == NO_IDX means the bond doesn't exist on that side (IS_NULL featurization).
// ---------------------------------------------------------------------------
static void fill_single_bond_feats(const GraphData&            gd,
                                   uint32_t                    bond_idx,
                                   float*                      buf,
                                   size_t                      single_bond_fdim,
                                   const py::array_t<int64_t>& bond_property_list) {
  const size_t   n_props = (bond_property_list.ndim() == 1) ? size_t(bond_property_list.shape(0)) : 0;
  const int64_t* props   = (n_props > 0) ? static_cast<const int64_t*>(bond_property_list.data()) : nullptr;
  const bool     is_null = (bond_idx == NO_IDX);

  float* p = buf;
  for (size_t i = 0; i < n_props; ++i) {
    auto feat = BondFeature(props[i]);
    switch (feat) {
      case BondFeature::IS_NULL:
        *p++ = is_null ? 1.0f : 0.0f;
        break;
      case BondFeature::TYPE_ONE_HOT: {
        // SINGLE→0, DOUBLE→1, TRIPLE→2, other(AROMATIC etc.)→3
        size_t fsz = get_one_hot_bond_feature_size(feat);
        std::fill(p, p + fsz, 0.0f);
        if (!is_null) {
          uint8_t bt = gd.bonds[bond_idx].bondType;
          size_t  slot;
          if (bt == 1)
            slot = 0;  // RDKit::Bond::SINGLE
          else if (bt == 2)
            slot = 1;  // RDKit::Bond::DOUBLE
          else if (bt == 3)
            slot = 2;  // RDKit::Bond::TRIPLE
          else
            slot = 3;  // AROMATIC or other
          p[slot] = 1.0f;
        }
        p += fsz;
        break;
      }
      case BondFeature::STEREO_ONE_HOT: {
        // STEREONONE=0, STEREOANY=1, STEREOZ=2, STEREOE=3, STEREOCIS=4, STEREOTRANS=5, other=6
        size_t fsz = get_one_hot_bond_feature_size(feat);
        std::fill(p, p + fsz, 0.0f);
        if (!is_null) {
          uint8_t st   = gd.bonds[bond_idx].stereo;
          size_t  slot = (size_t(st) < fsz) ? size_t(st) : fsz - 1;
          p[slot]      = 1.0f;
        }
        p += fsz;
        break;
      }
      case BondFeature::IN_RING:
        *p++ = (!is_null && gd.bonds[bond_idx].isInRing) ? 1.0f : 0.0f;
        break;
      case BondFeature::CONJUGATED:
        *p++ = (!is_null && gd.bonds[bond_idx].isConjugated) ? 1.0f : 0.0f;
        break;
      case BondFeature::TYPE_FLOAT:
        *p++ = is_null ? 0.0f : float(gd.bonds[bond_idx].bondType);
        break;
      default:
        *p++ = 0.0f;
        break;
    }
  }
  (void)single_bond_fdim;
}

// ---------------------------------------------------------------------------
// Internal: write CGR bond features for one undirected bond pair into out[0..cgr_bond_fdim-1].
// cgr_bond_fdim = 2 * single_bond_fdim.
// Diff direction: ALWAYS prod - reac.
// ---------------------------------------------------------------------------
static void write_cgr_bond_feats(const CompactReaction&      rxn,
                                 uint32_t                    b_reac_idx,
                                 uint32_t                    b_prod_idx,
                                 ReactionMode                mode,
                                 float*                      out,
                                 size_t                      single_bond_fdim,
                                 const py::array_t<int64_t>& bond_property_list) {
  const bool prod_first = (mode == ReactionMode::PROD_DIFF || mode == ReactionMode::PROD_DIFF_BALANCE);
  const bool use_diff   = (mode == ReactionMode::REAC_DIFF || mode == ReactionMode::PROD_DIFF ||
                         mode == ReactionMode::REAC_DIFF_BALANCE || mode == ReactionMode::PROD_DIFF_BALANCE);

  // Fill per-side feature vectors
  std::vector<float> reac_bf(single_bond_fdim, 0.0f);
  std::vector<float> prod_bf(single_bond_fdim, 0.0f);
  fill_single_bond_feats(rxn.reac, b_reac_idx, reac_bf.data(), single_bond_fdim, bond_property_list);
  fill_single_bond_feats(rxn.prod, b_prod_idx, prod_bf.data(), single_bond_fdim, bond_property_list);
  // Note: BALANCE adjustment is applied by the caller (enumerate_cgr_bonds) before calling this function.
  // The b_reac_idx / b_prod_idx passed in are already the post-balance-adjusted indices.
  // No BALANCE copying is done here.

  // First half
  const float* first_src = prod_first ? prod_bf.data() : reac_bf.data();
  std::copy(first_src, first_src + single_bond_fdim, out);

  // Second half: diff = prod - reac (always), or prod feats for REAC_PROD modes
  float* second_out = out + single_bond_fdim;
  if (use_diff) {
    for (size_t k = 0; k < single_bond_fdim; ++k)
      second_out[k] = prod_bf[k] - reac_bf[k];  // ALWAYS prod - reac
  } else {
    // REAC_PROD / REAC_PROD_BALANCE: second half = prod feats
    std::copy(prod_bf.begin(), prod_bf.end(), second_out);
  }
}

// ---------------------------------------------------------------------------
// Internal: bond enumeration result for one reaction.
// ---------------------------------------------------------------------------
struct BondEnumResult {
  size_t                     num_directed;    // total directed edges (undirected * 2)
  std::unique_ptr<float[]>   bond_feats;      // [num_directed, cgr_bond_fdim]
  std::unique_ptr<int64_t[]> edge_index;      // [2 * num_directed]: sources then dests
  std::unique_ptr<int64_t[]> rev_edge_index;  // [num_directed]: reverse edge index (local)
};

// ---------------------------------------------------------------------------
// Internal: enumerate CGR bonds for one reaction via the O(n²) scan that
// mirrors Python's CondensedGraphOfReactionFeaturizer ordering exactly.
// atom_offset: global atom index offset for this reaction (added to CGR node indices).
// ---------------------------------------------------------------------------
static BondEnumResult enumerate_cgr_bonds(const CompactReaction&      rxn,
                                          ReactionMode                mode,
                                          size_t                      cgr_bond_fdim,
                                          size_t                      single_bond_fdim,
                                          const py::array_t<int64_t>& bond_property_list,
                                          size_t                      atom_offset) {
  const size_t n_reac = rxn.reac.num_atoms;
  const size_t n_cgr  = n_reac + rxn.prod_only_idxs.size();

  // Helper: get product-side atom index for CGR node u (NO_IDX if u has no product atom)
  auto get_p = [&](size_t u) -> uint32_t {
    if (u >= n_reac)
      return rxn.prod_only_idxs[u - n_reac];
    auto it = rxn.r2p_idx_map.find(uint32_t(u));
    return (it != rxn.r2p_idx_map.end()) ? it->second : NO_IDX;
  };

  // Helper: look up a bond in a lookup map; returns NO_IDX if not found
  auto lookup_bond = [](const std::unordered_map<uint64_t, uint32_t>& lut, uint32_t a, uint32_t b) -> uint32_t {
    if (a == NO_IDX || b == NO_IDX)
      return NO_IDX;
    if (a > b)
      std::swap(a, b);
    auto it = lut.find((uint64_t(a) << 32) | uint64_t(b));
    return (it != lut.end()) ? it->second : NO_IDX;
  };

  // First pass: count undirected bonds
  size_t n_undirected = 0;
  for (size_t u = 0; u < n_cgr; ++u) {
    uint32_t r_u = (u < n_reac) ? uint32_t(u) : NO_IDX;
    uint32_t p_u = get_p(u);
    for (size_t v = u + 1; v < n_cgr; ++v) {
      uint32_t r_v = (v < n_reac) ? uint32_t(v) : NO_IDX;
      uint32_t p_v = get_p(v);
      if (lookup_bond(rxn.reac_bond_lookup, r_u, r_v) == NO_IDX &&
          lookup_bond(rxn.prod_bond_lookup, p_u, p_v) == NO_IDX)
        continue;
      ++n_undirected;
    }
  }

  const size_t               n_directed = 2 * n_undirected;
  std::unique_ptr<float[]>   bond_feats(new float[n_directed * cgr_bond_fdim]());
  std::unique_ptr<int64_t[]> edge_index(new int64_t[2 * n_directed]);
  std::unique_ptr<int64_t[]> rev_edge_index(new int64_t[n_directed]);

  const bool is_balance = (mode == ReactionMode::REAC_DIFF_BALANCE || mode == ReactionMode::REAC_PROD_BALANCE ||
                           mode == ReactionMode::PROD_DIFF_BALANCE);

  // Second pass: fill arrays
  size_t ep = 0;  // undirected bond index
  for (size_t u = 0; u < n_cgr; ++u) {
    uint32_t r_u = (u < n_reac) ? uint32_t(u) : NO_IDX;
    uint32_t p_u = get_p(u);
    for (size_t v = u + 1; v < n_cgr; ++v) {
      uint32_t r_v    = (v < n_reac) ? uint32_t(v) : NO_IDX;
      uint32_t p_v    = get_p(v);
      uint32_t b_reac = lookup_bond(rxn.reac_bond_lookup, r_u, r_v);
      uint32_t b_prod = lookup_bond(rxn.prod_bond_lookup, p_u, p_v);
      if (b_reac == NO_IDX && b_prod == NO_IDX)
        continue;

      // Apply BALANCE bond copying — mirrors Python's _get_bonds logic exactly:
      // * Both product-only (u>=n_reac && v>=n_reac): b_reac = b_prod
      // * Both reactant-only (u<n_reac, v<n_reac, neither mapped): b_prod = b_reac
      // * All other cases (including bonds between matched atoms where one side
      //   happens to be absent): no copying — missing side stays as null.
      if (is_balance) {
        bool u_prod_only = (u >= n_reac);
        bool v_prod_only = (v >= n_reac);
        if (u_prod_only && v_prod_only) {
          b_reac = b_prod;  // both product-only
        } else if (!u_prod_only && !v_prod_only) {
          bool u_matched = (p_u != NO_IDX);
          bool v_matched = (p_v != NO_IDX);
          if (!u_matched && !v_matched)
            b_prod = b_reac;  // both reactant-only
          // else: at least one is matched → no copy (null stays null)
        }
        // mixed case (one reactant atom, one product-only): no copy
      }

      // Directed edges: fwd (u→v) at position 2*ep, rev (v→u) at 2*ep+1
      size_t fwd = 2 * ep, rev = 2 * ep + 1;

      // Bond features (same for both directions)
      write_cgr_bond_feats(rxn,
                           b_reac,
                           b_prod,
                           mode,
                           bond_feats.get() + fwd * cgr_bond_fdim,
                           single_bond_fdim,
                           bond_property_list);
      std::copy(bond_feats.get() + fwd * cgr_bond_fdim,
                bond_feats.get() + fwd * cgr_bond_fdim + cgr_bond_fdim,
                bond_feats.get() + rev * cgr_bond_fdim);

      // edge_index: [sources | dests], interleaved fwd/rev
      edge_index[fwd]              = int64_t(u + atom_offset);  // fwd src
      edge_index[rev]              = int64_t(v + atom_offset);  // rev src
      edge_index[n_directed + fwd] = int64_t(v + atom_offset);  // fwd dst
      edge_index[n_directed + rev] = int64_t(u + atom_offset);  // rev dst

      // rev_edge_index: fwd and rev are each other's reverse (local indices)
      rev_edge_index[fwd] = int64_t(rev);
      rev_edge_index[rev] = int64_t(fwd);

      ++ep;
    }
  }

  return BondEnumResult{n_directed, std::move(bond_feats), std::move(edge_index), std::move(rev_edge_index)};
}

// ---------------------------------------------------------------------------
// batch_reaction_featurizer
// ---------------------------------------------------------------------------
std::vector<py::array> batch_reaction_featurizer(const std::vector<std::string>& reac_smiles_list,
                                                 const std::vector<std::string>& prod_smiles_list,
                                                 const py::array_t<int64_t>&     atom_property_list_onehot,
                                                 const py::array_t<int64_t>&     atom_property_list_float,
                                                 const py::array_t<int64_t>&     bond_property_list,
                                                 bool                            keep_h,
                                                 bool                            add_h,
                                                 bool                            offset_carbon,
                                                 ReactionMode                    mode) {
  if (reac_smiles_list.size() != prod_smiles_list.size())
    throw std::runtime_error("reac_smiles_list and prod_smiles_list must have the same length");
  const size_t n_rxns = reac_smiles_list.size();

  // Parse all reactions
  std::vector<CompactReaction> reactions;
  reactions.reserve(n_rxns);
  for (size_t i = 0; i < n_rxns; ++i)
    reactions.push_back(parse_reaction(reac_smiles_list[i], prod_smiles_list[i], keep_h, add_h));

  // Feature dimensions
  const size_t single_atom_fdim   = compute_atom_dim(atom_property_list_onehot, atom_property_list_float);
  // atomic_num_block_w = size of the atomic-num one-hot block (including "other" slot)
  // = get_one_hot_atom_feature_size(first_onehot_feature)
  size_t       atomic_num_block_w = 0;
  if (atom_property_list_onehot.ndim() == 1 && atom_property_list_onehot.shape(0) > 0) {
    auto f             = AtomOneHotFeature(static_cast<const int64_t*>(atom_property_list_onehot.data())[0]);
    atomic_num_block_w = get_one_hot_atom_feature_size(f);
  }
  const size_t second_atom_len  = single_atom_fdim - atomic_num_block_w;
  const size_t cgr_atom_fdim    = single_atom_fdim + second_atom_len;
  const size_t single_bond_fdim = compute_bond_dim(bond_property_list);
  const size_t cgr_bond_fdim    = 2 * single_bond_fdim;

  // Total CGR atom count across all reactions
  size_t total_cgr_atoms = 0;
  for (const auto& rxn : reactions)
    total_cgr_atoms += rxn.reac.num_atoms + rxn.prod_only_idxs.size();

  // Allocate atom feature array and batch array
  std::unique_ptr<float[]>   atom_data(new float[total_cgr_atoms * cgr_atom_fdim]());
  std::unique_ptr<int64_t[]> batch_data(new int64_t[total_cgr_atoms ? total_cgr_atoms : 1]);

  size_t atom_offset = 0;
  for (size_t i = 0; i < n_rxns; ++i) {
    const CompactReaction& rxn   = reactions[i];
    const size_t           n_cgr = rxn.reac.num_atoms + rxn.prod_only_idxs.size();
    fill_cgr_atom_features(rxn,
                           mode,
                           atom_data.get() + atom_offset * cgr_atom_fdim,
                           cgr_atom_fdim,
                           single_atom_fdim,
                           atomic_num_block_w,
                           atom_property_list_onehot,
                           atom_property_list_float,
                           offset_carbon);
    for (size_t k = 0; k < n_cgr; ++k)
      batch_data[atom_offset + k] = int64_t(i);
    atom_offset += n_cgr;
  }

  // Enumerate bonds per reaction
  std::vector<BondEnumResult> bond_results;
  bond_results.reserve(n_rxns);
  size_t total_directed = 0;
  atom_offset           = 0;
  for (size_t i = 0; i < n_rxns; ++i) {
    const CompactReaction& rxn = reactions[i];
    bond_results.push_back(
      enumerate_cgr_bonds(rxn, mode, cgr_bond_fdim, single_bond_fdim, bond_property_list, atom_offset));
    total_directed += bond_results.back().num_directed;
    atom_offset += rxn.reac.num_atoms + rxn.prod_only_idxs.size();
  }

  // Assemble global bond_feats, edge_index, rev_edge_index
  std::unique_ptr<float[]>   bond_data(new float[total_directed * cgr_bond_fdim + 1]());
  std::unique_ptr<int64_t[]> edge_index(new int64_t[2 * total_directed + 1]);
  std::unique_ptr<int64_t[]> rev_edge_index(new int64_t[total_directed + 1]);

  size_t bond_offset = 0;
  for (size_t i = 0; i < n_rxns; ++i) {
    const BondEnumResult& br = bond_results[i];
    const size_t          n  = br.num_directed;
    // Bond features
    std::copy(br.bond_feats.get(),
              br.bond_feats.get() + n * cgr_bond_fdim,
              bond_data.get() + bond_offset * cgr_bond_fdim);
    // Source indices
    std::copy(br.edge_index.get(), br.edge_index.get() + n, edge_index.get() + bond_offset);
    // Dest indices (second half of local array → second half of global array)
    std::copy(br.edge_index.get() + n, br.edge_index.get() + 2 * n, edge_index.get() + total_directed + bond_offset);
    // rev_edge_index: shift local indices by bond_offset
    for (size_t k = 0; k < n; ++k)
      rev_edge_index[bond_offset + k] = br.rev_edge_index[k] + int64_t(bond_offset);
    bond_offset += n;
  }

  // Wrap in py::array_t
  const int64_t atom_dims[2]  = {int64_t(total_cgr_atoms), int64_t(cgr_atom_fdim)};
  const int64_t bond_dims[2]  = {int64_t(total_directed), int64_t(cgr_bond_fdim)};
  const int64_t ei_dims[2]    = {int64_t(2), int64_t(total_directed)};
  const int64_t rev_dims[1]   = {int64_t(total_directed)};
  const int64_t batch_dims[1] = {int64_t(total_cgr_atoms)};

  return {
    py_array_from_array<float>(std::move(atom_data), atom_dims, 2),
    py_array_from_array<float>(std::move(bond_data), bond_dims, 2),
    py_array_from_array<int64_t>(std::move(edge_index), ei_dims, 2),
    py_array_from_array<int64_t>(std::move(rev_edge_index), rev_dims, 1),
    py_array_from_array<int64_t>(std::move(batch_data), batch_dims, 1),
  };
}
