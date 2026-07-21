// SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

//! @file This file defines generic feature-related functions,
//!       some of which are declared in features.h for exporting to Python.

#include <cstdint>
#define DEBUG_LOGGING 0

#include <GraphMol/MolOps.h>  // For RDKit's addHs

#include "features.h"
#include "float_features.h"
#include "one_hot.h"
// RDkit headers
#include <GraphMol/Atom.h>
#include <GraphMol/Canon.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/types.h>

#include <unordered_map>

// This is called by `mol_featurizer` and `batch_mol_featurizer` to parse the SMILES string into an RWMol and
// cache some data about the atoms and bonds.
static GraphData read_graph(const std::string& smiles_string, bool explicit_H, bool ordered) {
  std::unique_ptr<RDKit::RWMol> mol{parse_mol(smiles_string, explicit_H, ordered)};

  if (!mol) {
    return GraphData{0, std::unique_ptr<CompactAtom[]>(), 0, std::unique_ptr<CompactBond[]>(), std::move(mol)};
  }

  const size_t num_atoms = mol->getNumAtoms();
  const size_t num_bonds = mol->getNumBonds();
#if DEBUG_LOGGING
  printf("# atoms = %zu\n# bonds = %zu\n", num_atoms, num_bonds);
#endif
#if REPORT_STATS
  ++statsMolAtomCounts[(num_atoms >= STATS_NUM_MOL_ATOM_COUNTS) ? (STATS_NUM_MOL_ATOM_COUNTS - 1) : num_atoms];
  ++statsMolBondCounts[(num_bonds >= STATS_NUM_MOL_BOND_COUNTS) ? (STATS_NUM_MOL_BOND_COUNTS - 1) : num_bonds];
  statsTotalNumAtoms += num_atoms;
  statsTotalNumBonds += num_bonds;
#endif

#if ORDER_ATOMS
  // Determine a canonical ordering of the atoms, if desired.
  std::vector<unsigned int> atomOrder;
  atomOrder.reserve(num_atoms);
  RDKit::Canon::rankMolAtoms(*mol, atomOrder);
  assert(atomOrder.shape() == num_atoms);
#endif

  // Allocate an array of atom data, and fill it from the RDKit atom data.
  std::unique_ptr<CompactAtom[]> atoms(new CompactAtom[num_atoms]);
  for (size_t atomIdx = 0; atomIdx < num_atoms; ++atomIdx) {
    const RDKit::Atom* const             atom          = mol->getAtomWithIdx(atomIdx);
    auto                                 atomicNum     = atom->getAtomicNum();
    auto                                 totalDegree   = atom->getTotalDegree();
    auto                                 formalCharge  = atom->getFormalCharge();
    const RDKit::Atom::ChiralType        chiralType    = atom->getChiralTag();
    auto                                 totalNumHs    = atom->getTotalNumHs();
    const RDKit::Atom::HybridizationType hybridization = atom->getHybridization();

    const bool isAromatic = atom->getIsAromatic();
#if REPORT_STATS
    ++statsElementCounts[(atomicNum < 0 || atomicNum >= STATS_NUM_ELEMENTS) ? (STATS_NUM_ELEMENTS - 1) : atomicNum];
    ++statsDegreeCounts[(totalDegree < 0 || totalDegree >= STATS_NUM_DEGREES) ? (STATS_NUM_DEGREES - 1) : totalDegree];
    size_t formalChargeIndex = formalCharge + int(STATS_CHARGE_OFFSET);
    if (formalCharge < -int(STATS_CHARGE_OFFSET)) {
      formalChargeIndex = 0;
    } else if (formalCharge > int(STATS_CHARGE_OFFSET)) {
      formalChargeIndex = STATS_NUM_CHARGES - 1;
    }

    ++statsChargeCounts[formalChargeIndex];
    ++statsChiralityCounts[(size_t(chiralType) >= STATS_NUM_CHIRALITIES) ? (STATS_NUM_CHIRALITIES - 1) :
                                                                           size_t(chiralType)];
    ++statsHCounts[(totalNumHs < 0 || totalNumHs >= STATS_NUM_HS) ? (STATS_NUM_HS - 1) : totalNumHs];
    ++statsHybridizationCounts[(size_t(hybridization) >= STATS_NUM_HYBRIDIZATIONS) ? (STATS_NUM_HYBRIDIZATIONS - 1) :
                                                                                     size_t(hybridization)];
    statsAromaticAtomCount += (isAromatic ? 1 : 0);
#endif
    const double mass = atom->getMass();

#if ORDER_ATOMS
    const size_t destAtomIdx = atomOrder[atomIdx];
#else
    const size_t destAtomIdx = atomIdx;
#endif
    atoms[destAtomIdx] = CompactAtom{uint8_t(atomicNum),
                                     uint8_t(totalDegree),
                                     int8_t(formalCharge),
                                     uint8_t(chiralType),
                                     uint8_t(totalNumHs),
                                     uint8_t(hybridization),
                                     isAromatic,
                                     float(mass)};
#if DEBUG_LOGGING
    printf("atom[%zu] = {%zu, %u, %d, %u, %u, %u, %s, %f}\n",
           destAtomIdx,
           int(atomicNum),
           int(totalDegree),
           int(formalCharge),
           int(chiralType),
           int(totalNumHs),
           int(hybridization),
           isAromatic ? "true" : "false",
           mass);
#endif
  }

  // Allocate an array of bond data, and fill it from the RDKit bond data.
  std::unique_ptr<CompactBond[]> bonds(new CompactBond[num_bonds]);
  const RDKit::RingInfo* const   ringInfo = mol->getRingInfo();
  for (size_t bondIdx = 0; bondIdx < num_bonds; ++bondIdx) {
    const RDKit::Bond* const      bond         = mol->getBondWithIdx(bondIdx);
    const RDKit::Bond::BondType   bondType     = bond->getBondType();
    const bool                    isConjugated = bond->getIsConjugated();
    // TODO: Verify that it's the same index as bond->getIdx()
    const bool                    isInRing     = (ringInfo->numBondRings(bondIdx) != 0);
    const RDKit::Bond::BondStereo stereo       = bond->getStereo();

#if REPORT_STATS
    ++statsBondTypeCounts[(size_t(bondType) >= STATS_NUM_BOND_TYPES) ? (STATS_NUM_BOND_TYPES - 1) : size_t(bondType)];
    ++statsBondStereoCounts[(size_t(stereo) >= STATS_NUM_BOND_STEREOS) ? (STATS_NUM_BOND_STEREOS - 1) : size_t(stereo)];
    statsConjugatedBondCount += (isConjugated ? 1 : 0);
    statsBondInRingCount += (isInRing ? 1 : 0);
#endif

    auto beginAtomIdx = bond->getBeginAtomIdx();
    auto endAtomIdx   = bond->getEndAtomIdx();
#if ORDER_ATOMS
    beginAtomIdx = atomOrder[beginAtomIdx];
    endAtomIdx   = atomOrder[endAtomIdx];
#endif
    bonds[bondIdx] = CompactBond{uint8_t(bondType), isConjugated, isInRing, uint8_t(stereo), beginAtomIdx, endAtomIdx};
#if DEBUG_LOGGING
    printf("bond[%zu] = {%u, %s, %s, %u, {%u, %u}}\n",
           bondIdx,
           int(bondType),
           isConjugated ? "true" : "false",
           isInRing ? "true" : "false",
           int(stereo),
           beginAtomIdx,
           endAtomIdx);
#endif
  }

  // Return a GraphData structure, taking ownership of the atom and bond data arrays.
  return GraphData{num_atoms, std::move(atoms), num_bonds, std::move(bonds), std::move(mol)};
}

// This is a structure for managing the adjacency data (CSR format)
struct NeighborData {
  // This owns the data of all 3 arrays, which are actually a single, contiguous allocation.
  std::unique_ptr<uint32_t[]> deleter;

  // This is an array of indices into the other two arrays, indicating where
  // each atom's neighbors start, including the first entry being 0 for the start of
  // atom 0, and the num_atoms entry being 2*num_bonds (2x because each bond is on 2 atoms),
  // so there are num_atoms+1 entries.  The number of neighbors of an atom i is
  // neighbor_starts[i+1]-neighbor_starts[i]
  const uint32_t* neighbor_starts;

  // The neighbor atom for each bond, with each atom having an entry for each of
  // its neighbors, so each bond occurs twice.
  const uint32_t* neighbors;

  // This is in the same order as neighbors, but indicates the index of the bond.
  // Each bond occurs twice, so each number occurs twice.
  const uint32_t* bond_indices;
};

// Construct a NeighborData structure representing the molecule's graph in CSR format.
static NeighborData construct_neighbors(const GraphData& graph) {
  const uint32_t              num_atoms = graph.num_atoms;
  const uint32_t              num_bonds = graph.num_bonds;
  // Do a single allocation for all 3 arrays.
  std::unique_ptr<uint32_t[]> deleter(new uint32_t[num_atoms + 1 + 4 * num_bonds]);

  uint32_t* neighbor_starts = deleter.get();
  for (uint32_t i = 0; i <= num_atoms; ++i) {
    neighbor_starts[i] = 0;
  }

  // First, get atom neighbor counts
  const CompactBond* const bonds = graph.bonds.get();
  for (uint32_t i = 0; i < num_bonds; ++i) {
    uint32_t a = bonds[i].beginAtomIdx;
    uint32_t b = bonds[i].endAtomIdx;
    // NOTE: +1 is because first entry will stay zero.
    ++neighbor_starts[a + 1];
    ++neighbor_starts[b + 1];
  }

  // Find the starts by partial-summing the neighbor counts.
  // NOTE: +1 is because first entry will stay zero.
  std::partial_sum(neighbor_starts + 1, neighbor_starts + 1 + num_atoms, neighbor_starts + 1);

  // Fill in the neighbors and bond_indices arrays.
  uint32_t* neighbors    = neighbor_starts + num_atoms + 1;
  uint32_t* bond_indices = neighbors + 2 * num_bonds;
  for (uint32_t i = 0; i < num_bonds; ++i) {
    uint32_t a = bonds[i].beginAtomIdx;
    uint32_t b = bonds[i].endAtomIdx;

    uint32_t ai      = neighbor_starts[a];
    neighbors[ai]    = b;
    bond_indices[ai] = i;
    ++neighbor_starts[a];

    uint32_t bi      = neighbor_starts[b];
    neighbors[bi]    = a;
    bond_indices[bi] = i;
    ++neighbor_starts[b];
  }

  // Shift neighbor_starts forward one after incrementing it.
  uint32_t previous = 0;
  for (uint32_t i = 0; i < num_atoms; ++i) {
    uint32_t next      = neighbor_starts[i];
    neighbor_starts[i] = previous;
    previous           = next;
  }

  // NeighborData takes ownership of the memory.
  return NeighborData{std::move(deleter), neighbor_starts, neighbors, bond_indices};
}

// Maps float atom feature name strings to `AtomFloatFeature` enum values
static const std::unordered_map<std::string, int64_t> atom_float_name_to_enum{
  {         std::string("atomic-number"),          int64_t(AtomFloatFeature::ATOMIC_NUMBER)},
  {                  std::string("mass"),                   int64_t(AtomFloatFeature::MASS)},
  {                std::string("weight"),                   int64_t(AtomFloatFeature::MASS)},
  {               std::string("valence"),                int64_t(AtomFloatFeature::VALENCE)},
  {         std::string("total-valence"),                int64_t(AtomFloatFeature::VALENCE)},
  {      std::string("implicit-valence"),       int64_t(AtomFloatFeature::IMPLICIT_VALENCE)},
  {         std::string("hybridization"),          int64_t(AtomFloatFeature::HYBRIDIZATION)},
  {             std::string("chirality"),              int64_t(AtomFloatFeature::CHIRALITY)},
  {              std::string("aromatic"),               int64_t(AtomFloatFeature::AROMATIC)},
  {                  std::string("ring"),                int64_t(AtomFloatFeature::IN_RING)},
  {               std::string("in-ring"),                int64_t(AtomFloatFeature::IN_RING)},
  {              std::string("min-ring"),               int64_t(AtomFloatFeature::MIN_RING)},
  {              std::string("max-ring"),               int64_t(AtomFloatFeature::MAX_RING)},
  {              std::string("num-ring"),               int64_t(AtomFloatFeature::NUM_RING)},
  {                std::string("degree"),                 int64_t(AtomFloatFeature::DEGREE)},
  {      std::string("radical-electron"),       int64_t(AtomFloatFeature::RADICAL_ELECTRON)},
  {         std::string("formal-charge"),          int64_t(AtomFloatFeature::FORMAL_CHARGE)},
  {                 std::string("group"),                  int64_t(AtomFloatFeature::GROUP)},
  {                std::string("period"),                 int64_t(AtomFloatFeature::PERIOD)},
  {           std::string("single-bond"),            int64_t(AtomFloatFeature::SINGLE_BOND)},
  {         std::string("aromatic-bond"),          int64_t(AtomFloatFeature::AROMATIC_BOND)},
  {           std::string("double-bond"),            int64_t(AtomFloatFeature::DOUBLE_BOND)},
  {           std::string("triple-bond"),            int64_t(AtomFloatFeature::TRIPLE_BOND)},
  {             std::string("is-carbon"),              int64_t(AtomFloatFeature::IS_CARBON)},
  {   std::string("hydrogen-bond-donor"),    int64_t(AtomFloatFeature::HYDROGEN_BOND_DONOR)},
  {std::string("hydrogen-bond-acceptor"), int64_t(AtomFloatFeature::HYDROGEN_BOND_ACCEPTOR)},
  {                std::string("acidic"),                 int64_t(AtomFloatFeature::ACIDIC)},
  {                 std::string("basic"),                  int64_t(AtomFloatFeature::BASIC)},
};

// This is called from Python to list atom float features in a format that will be faster
// to interpret inside `mol_featurizer`, passed in the `atom_property_list_float` parameter.
// See the declaration in features.h for more details.
py::array_t<int64_t> atom_float_feature_names_to_array(const std::vector<std::string>& features) {
  const size_t               num_features = features.size();
  std::unique_ptr<int64_t[]> feature_enum_values(new int64_t[num_features]);
  for (size_t i = 0; i < num_features; ++i) {
    auto it = atom_float_name_to_enum.find(features[i]);
    if (it != atom_float_name_to_enum.end()) {
      feature_enum_values[i] = it->second;
    } else {
      feature_enum_values[i] = int64_t(AtomFloatFeature::UNKNOWN);
    }
  }
  const int64_t dims[1] = {int64_t(num_features)};
  return py_array_from_array(std::move(feature_enum_values), dims, 1);
}

std::vector<std::string> list_all_atom_float_features() {
  std::vector<std::string> names;
  names.reserve(atom_float_name_to_enum.size());
  for (const auto& pair : atom_float_name_to_enum) {
    names.push_back(pair.first);
  }
  return names;
}

// Maps one-hot atom feature name strings to `AtomOneHotFeature` enum values
static const std::unordered_map<std::string, int64_t> atom_onehot_name_to_enum{
  {         std::string("atomic-number"),             int64_t(AtomOneHotFeature::ATOMIC_NUM)},
  {  std::string("atomic-number-common"),      int64_t(AtomOneHotFeature::ATOMIC_NUM_COMMON)},
  { std::string("atomic-number-organic"),     int64_t(AtomOneHotFeature::ATOMIC_NUM_ORGANIC)},
  {                std::string("degree"),                 int64_t(AtomOneHotFeature::DEGREE)},
  {          std::string("total-degree"),           int64_t(AtomOneHotFeature::TOTAL_DEGREE)},
  {               std::string("valence"),                int64_t(AtomOneHotFeature::VALENCE)},
  {         std::string("total-valence"),                int64_t(AtomOneHotFeature::VALENCE)},
  {      std::string("implicit-valence"),       int64_t(AtomOneHotFeature::IMPLICIT_VALENCE)},
  {         std::string("hybridization"),          int64_t(AtomOneHotFeature::HYBRIDIZATION)},
  {std::string("hybridization-expanded"), int64_t(AtomOneHotFeature::HYBRIDIZATION_EXPANDED)},
  { std::string("hybridization-organic"),  int64_t(AtomOneHotFeature::HYBRIDIZATION_ORGANIC)},
  {             std::string("chirality"),              int64_t(AtomOneHotFeature::CHIRALITY)},
  {                 std::string("group"),                  int64_t(AtomOneHotFeature::GROUP)},
  {                std::string("period"),                 int64_t(AtomOneHotFeature::PERIOD)},
  {         std::string("formal-charge"),          int64_t(AtomOneHotFeature::FORMAL_CHARGE)},
  {         std::string("num-hydrogens"),          int64_t(AtomOneHotFeature::NUM_HYDROGENS)},
  {             std::string("ring-size"),              int64_t(AtomOneHotFeature::RING_SIZE)},
};

// This is called from Python to list atom one-hot features in a format that will be faster
// to interpret inside `mol_featurizer` and `batch_mol_featurizer`, passed in the
// `atom_property_list_onehot` parameter.
py::array_t<int64_t> atom_onehot_feature_names_to_array(const std::vector<std::string>& features) {
  const size_t               num_features = features.size();
  std::unique_ptr<int64_t[]> feature_enum_values(new int64_t[num_features]);
  for (size_t i = 0; i < num_features; ++i) {
    auto it = atom_onehot_name_to_enum.find(features[i]);
    if (it != atom_onehot_name_to_enum.end()) {
      feature_enum_values[i] = it->second;
    } else {
      feature_enum_values[i] = int64_t(AtomOneHotFeature::UNKNOWN);
    }
  }
  const int64_t dims[1] = {int64_t(num_features)};
  return py_array_from_array(std::move(feature_enum_values), dims, 1);
}

std::vector<std::string> list_all_atom_onehot_features() {
  std::vector<std::string> names;
  names.reserve(atom_onehot_name_to_enum.size());
  for (const auto& pair : atom_onehot_name_to_enum) {
    names.push_back(pair.first);
  }
  return names;
}

// Maps bond feature name strings to `BondFeature` enum values
static const std::unordered_map<std::string, int64_t> bond_name_to_enum{
  {              std::string("is-null"),               int64_t(BondFeature::IS_NULL)},
  {     std::string("bond-type-onehot"),          int64_t(BondFeature::TYPE_ONE_HOT)},
  {      std::string("bond-type-float"),            int64_t(BondFeature::TYPE_FLOAT)},
  {               std::string("stereo"),        int64_t(BondFeature::STEREO_ONE_HOT)},
  {              std::string("in-ring"),               int64_t(BondFeature::IN_RING)},
  {           std::string("conjugated"),            int64_t(BondFeature::CONJUGATED)},
  {std::string("conformer-bond-length"), int64_t(BondFeature::CONFORMER_BOND_LENGTH)},
  {std::string("estimated-bond-length"), int64_t(BondFeature::ESTIMATED_BOND_LENGTH)},
};

std::vector<std::string> list_all_bond_features() {
  std::vector<std::string> names;
  names.reserve(bond_name_to_enum.size());
  for (const auto& pair : bond_name_to_enum) {
    names.push_back(pair.first);
  }
  return names;
}

// This is called from Python to list bond features in a format that will be faster
// to interpret inside `mol_featurizer`, passed in the `bond_property_list` parameter.
// See the declaration in features.h for more details.
py::array_t<int64_t> bond_feature_names_to_array(const std::vector<std::string>& features) {
  const size_t               num_features = features.size();
  std::unique_ptr<int64_t[]> feature_enum_values(new int64_t[num_features]);
  for (size_t i = 0; i < num_features; ++i) {
    auto it = bond_name_to_enum.find(features[i]);
    if (it != bond_name_to_enum.end()) {
      feature_enum_values[i] = it->second;
    } else {
      feature_enum_values[i] = int64_t(BondFeature::UNKNOWN);
    }
  }
  const int64_t dims[1] = {int64_t(num_features)};
  return py_array_from_array(std::move(feature_enum_values), dims, 1);
}

// Maps feature level strings to `FeatureLevel` enum values
static const std::unordered_map<std::string, int64_t> feature_level_to_enum{
  {    std::string("node"),     int64_t(FeatureLevel::NODE)},
  {    std::string("edge"),     int64_t(FeatureLevel::EDGE)},
  {std::string("nodepair"), int64_t(FeatureLevel::NODEPAIR)},
  {   std::string("graph"),    int64_t(FeatureLevel::GRAPH)},
};

// This is called by `mol_featurizer` and `batch_mol_featurizer` to create the atom (node) features pybind array.
template <typename T>
void create_atom_features(const GraphData&            graph,
                          const py::array_t<int64_t>& atom_property_list_onehot,
                          const py::array_t<int64_t>& atom_property_list_float,
                          size_t                      single_atom_float_count,
                          bool                        offset_carbon,
                          T*                          current_atom_data) {
  const size_t num_onehot_properties =
    (atom_property_list_onehot.dtype() == py::dtype("int64") && atom_property_list_onehot.ndim() == 1) ?
      atom_property_list_onehot.shape(0) :
      0;
  const int64_t* const property_list_onehot =
    (num_onehot_properties != 0) ? static_cast<const int64_t*>(atom_property_list_onehot.data()) : nullptr;
  const size_t num_float_properties =
    (atom_property_list_float.dtype() == py::dtype("int64") && atom_property_list_float.ndim() == 1) ?
      atom_property_list_float.shape(0) :
      0;
  const int64_t* const property_list_float =
    (num_float_properties != 0) ? static_cast<const int64_t*>(atom_property_list_float.data()) : nullptr;

  for (size_t i = 0; i < num_onehot_properties; ++i) {
    const int64_t property = property_list_onehot[i];
    current_atom_data +=
      get_one_hot_atom_feature(graph, current_atom_data, AtomOneHotFeature(property), single_atom_float_count);
  }
  for (size_t i = 0; i < num_float_properties; ++i) {
    const int64_t property = property_list_float[i];
    get_atom_float_feature(graph,
                           current_atom_data,
                           AtomFloatFeature(property),
                           single_atom_float_count,
                           offset_carbon);
    ++current_atom_data;
  }
}

// This is called by `mol_featurizer` and `batch_mol_featurizer` to create the bond (edge) features pybind array.
template <typename T>
void create_bond_features(const GraphData&            graph,
                          const py::array_t<int64_t>& bond_property_list,
                          size_t                      single_bond_float_count,
                          const bool                  duplicate_edges,
                          bool                        add_self_loop,
                          T*                          current_bond_data) {
  const size_t num_properties = (bond_property_list.dtype() == py::dtype("int64") && bond_property_list.ndim() == 1) ?
                                  bond_property_list.shape(0) :
                                  0;
  const int64_t* const property_list =
    (num_properties != 0) ? static_cast<const int64_t*>(bond_property_list.data()) : nullptr;

  // add_self_loop is only supported if duplicating edges
  add_self_loop = add_self_loop && duplicate_edges;

  // This is the stride length (in floats) for each unique bond
  const size_t duplicated_bond_float_count = duplicate_edges ? (2 * single_bond_float_count) : single_bond_float_count;

  // Make a copy of current_bond_data pointer for use in duplicate_edges
  T* current_bond_data_copy = current_bond_data;
  for (size_t i = 0; i < num_properties; ++i) {
    const int64_t property = property_list[i];
    if (BondFeature(property) == BondFeature::TYPE_ONE_HOT || BondFeature(property) == BondFeature::STEREO_ONE_HOT) {
      current_bond_data +=
        get_one_hot_bond_feature(graph, current_bond_data, BondFeature(property), duplicated_bond_float_count);
    } else {
      get_bond_float_feature(graph, current_bond_data, BondFeature(property), duplicated_bond_float_count);
      ++current_bond_data;
    }
  }
  if (duplicate_edges) {
    // Duplicate the data for each bond
    for (size_t i = 0; i < graph.num_bonds; ++i) {
      for (size_t j = 0; j < single_bond_float_count; ++j) {
        current_bond_data_copy[i * duplicated_bond_float_count + j + single_bond_float_count] =
          current_bond_data_copy[i * duplicated_bond_float_count + j];
      }
      // current_bond_data += duplicated_bond_float_count;
    }
    if (add_self_loop) {
      // Self loops don't have valid bond data, but don't treat them as NaNs.
      // Fill with zeros, instead.
      memset(current_bond_data, 0, graph.num_atoms * graph.num_atoms);
    }
  }
}

// Computes the total dimension of atom features based on the property lists
size_t compute_atom_dim(const py::array_t<int64_t>& atom_property_list_onehot,
                        const py::array_t<int64_t>& atom_property_list_float) {
  const size_t num_onehot_properties =
    (atom_property_list_onehot.dtype() == py::dtype("int64") && atom_property_list_onehot.ndim() == 1) ?
      atom_property_list_onehot.shape(0) :
      0;
  const int64_t* const property_list_onehot =
    (num_onehot_properties != 0) ? static_cast<const int64_t*>(atom_property_list_onehot.data()) : nullptr;
  const size_t num_float_properties =
    (atom_property_list_float.dtype() == py::dtype("int64") && atom_property_list_float.ndim() == 1) ?
      atom_property_list_float.shape(0) :
      0;
  const int64_t* const property_list_float =
    (num_float_properties != 0) ? static_cast<const int64_t*>(atom_property_list_float.data()) : nullptr;

  size_t single_atom_float_count = num_float_properties;
  for (size_t i = 0; i < num_onehot_properties; ++i) {
    const int64_t property = property_list_onehot[i];
    single_atom_float_count += get_one_hot_atom_feature_size(AtomOneHotFeature(property));
  }
  return single_atom_float_count;
}

// Computes the total dimension of bond features based on the property list
size_t compute_bond_dim(const py::array_t<int64_t>& bond_property_list) {
  const size_t num_properties = (bond_property_list.dtype() == py::dtype("int64") && bond_property_list.ndim() == 1) ?
                                  bond_property_list.shape(0) :
                                  0;
  const int64_t* const property_list =
    (num_properties != 0) ? static_cast<const int64_t*>(bond_property_list.data()) : nullptr;
  size_t single_bond_float_count = 0;
  for (size_t i = 0; i < num_properties; ++i) {
    const int64_t property = property_list[i];
    if (BondFeature(property) == BondFeature::TYPE_ONE_HOT || BondFeature(property) == BondFeature::STEREO_ONE_HOT) {
      single_bond_float_count += get_one_hot_bond_feature_size(BondFeature(property));
    } else {
      ++single_bond_float_count;
    }
  }
  return single_bond_float_count;
}

std::vector<py::array> mol_featurizer(const std::string&          smiles_string,
                                      const py::array_t<int64_t>& atom_property_list_onehot,
                                      const py::array_t<int64_t>& atom_property_list_float,
                                      const py::array_t<int64_t>& bond_property_list,
                                      bool                        explicit_H,
                                      bool                        offset_carbon,
                                      bool                        duplicate_edges,
                                      bool                        add_self_loop,
                                      bool                        ordered) {
  return batch_mol_featurizer(std::vector{smiles_string},
                              atom_property_list_onehot,
                              atom_property_list_float,
                              bond_property_list,
                              explicit_H,
                              offset_carbon,
                              duplicate_edges,
                              add_self_loop,
                              ordered);
}

std::vector<py::array> batch_mol_featurizer(const std::vector<std::string>& smiles_list,
                                            const py::array_t<int64_t>&     atom_property_list_onehot,
                                            const py::array_t<int64_t>&     atom_property_list_float,
                                            const py::array_t<int64_t>&     bond_property_list,
                                            bool                            explicit_H,
                                            bool                            offset_carbon,
                                            bool                            duplicate_edges,
                                            bool                            add_self_loop,
                                            bool                            ordered) {
  const size_t n_smiles = smiles_list.size();

  // Create graphs
  std::vector<GraphData> graph_list;
  graph_list.reserve(smiles_list.size());

  size_t total_num_atoms = 0, total_num_bonds = 0;

  for (const auto& smiles : smiles_list) {
    GraphData igraph = read_graph(smiles, explicit_H, ordered);
    total_num_atoms += igraph.num_atoms;
    total_num_bonds += igraph.num_bonds;
    graph_list.push_back(std::move(igraph));
  }

  // Compute atom dimension
  size_t       single_atom_float_count = compute_atom_dim(atom_property_list_onehot, atom_property_list_float);
  const size_t atom_float_count        = single_atom_float_count * total_num_atoms;
  std::unique_ptr<float[]>   atom_data(new float[atom_float_count]);
  std::unique_ptr<int64_t[]> batch(new int64_t[total_num_atoms]);
  // Create atom features
  float*                     current_atom_data = atom_data.get();

  for (size_t igraph = 0; igraph < n_smiles; ++igraph) {
    const GraphData& graph = graph_list[igraph];

    create_atom_features(graph,
                         atom_property_list_onehot,
                         atom_property_list_float,
                         single_atom_float_count,
                         offset_carbon,
                         current_atom_data);
    current_atom_data += graph.num_atoms * single_atom_float_count;
  }
  const int64_t      dims[2]                 = {int64_t(total_num_atoms), int64_t(single_atom_float_count)};
  py::array_t<float> atom_features_array     = py_array_from_array<float>(std::move(atom_data), dims, 2);
  // Compute bond dimension
  size_t             single_bond_float_count = compute_bond_dim(bond_property_list);

  // add_self_loop is only supported if duplicating edges
  add_self_loop = add_self_loop && duplicate_edges;

  if (duplicate_edges) {
    total_num_bonds = 2 * total_num_bonds + size_t(add_self_loop);
  }
  const size_t             bond_float_count = single_bond_float_count * total_num_bonds;
  // Create bond features
  std::unique_ptr<float[]> bond_data(new float[bond_float_count]);
  float*                   current_bond_data = bond_data.get();
  // This is the stride length (in floats) for each unique bond
  const size_t duplicated_bond_float_count = duplicate_edges ? (2 * single_bond_float_count) : single_bond_float_count;
  for (size_t igraph = 0; igraph < n_smiles; ++igraph) {
    const GraphData& graph = graph_list[igraph];

    create_bond_features(graph,
                         bond_property_list,
                         single_bond_float_count,
                         duplicate_edges,
                         add_self_loop,
                         current_bond_data);
    current_bond_data += graph.num_bonds * duplicated_bond_float_count;  // move to position
  }
  const int64_t      bond_dims[2]        = {int64_t(total_num_bonds), int64_t(single_bond_float_count)};
  py::array_t<float> bond_features_array = py_array_from_array<float>(std::move(bond_data), bond_dims, 2);

  // Create batch
  size_t mol_offset     = 0;
  size_t n_atoms_offset = 0;

  for (size_t igraph = 0; igraph < n_smiles; ++igraph) {
    const GraphData& graph = graph_list[igraph];
    for (size_t iatom = 0; iatom < graph.num_atoms; ++iatom) {
      batch[n_atoms_offset + iatom] = igraph;
    }
    mol_offset += total_num_bonds;
    n_atoms_offset += graph.num_atoms;
  }
  const int64_t              batch_dims[1] = {int64_t(total_num_atoms)};
  py::array_t<int64_t>       batch_array   = py_array_from_array<int64_t>(std::move(batch), batch_dims, 1);
  // Create edge_index
  std::unique_ptr<int64_t[]> edge_index(new int64_t[2 * total_num_bonds]);

  mol_offset     = 0;
  n_atoms_offset = 0;
  for (size_t igraph = 0; igraph < n_smiles; ++igraph) {
    const GraphData& graph = graph_list[igraph];
    for (size_t i = 0; i < graph.num_bonds; ++i) {
      if (duplicate_edges) {
        // PyG has all directed edge begin indices followed by all end indices.
        edge_index[mol_offset + 2 * i]                       = graph.bonds[i].beginAtomIdx + n_atoms_offset;
        edge_index[mol_offset + 2 * i + 1]                   = graph.bonds[i].endAtomIdx + n_atoms_offset;
        edge_index[mol_offset + total_num_bonds + 2 * i]     = graph.bonds[i].endAtomIdx + n_atoms_offset;
        edge_index[mol_offset + total_num_bonds + 2 * i + 1] = graph.bonds[i].beginAtomIdx + n_atoms_offset;
      } else {
        // PyG has all directed edge begin indices followed by all end indices.
        edge_index[mol_offset + i]                   = graph.bonds[i].beginAtomIdx + n_atoms_offset;
        edge_index[mol_offset + total_num_bonds + i] = graph.bonds[i].endAtomIdx + n_atoms_offset;
      }
    }
    if (duplicate_edges) {
      mol_offset += graph.num_bonds * 2;
    } else {
      mol_offset += graph.num_bonds;
    }
    n_atoms_offset += graph.num_atoms;
  }

  int64_t              edge_index_dims[2] = {int64_t(2), int64_t(total_num_bonds)};
  py::array_t<int64_t> edge_index_array   = py_array_from_array<int64_t>(std::move(edge_index), edge_index_dims, 2);

  // Create rev_edge_index
  mol_offset = 0;

  size_t                     duplication_factor = duplicate_edges ? 1 : 2;  // account for lack of doubling
  std::unique_ptr<int64_t[]> rev_edge_index(new int64_t[duplication_factor * total_num_bonds]);

  for (size_t igraph = 0; igraph < n_smiles; ++igraph) {
    const GraphData& graph = graph_list[igraph];
    for (size_t i = 0; i < graph.num_bonds * 2; i += 2) {
      rev_edge_index[mol_offset + i]     = i + 1 + mol_offset;
      rev_edge_index[mol_offset + i + 1] = i + mol_offset;
    }
    mol_offset += graph.num_bonds * 2;
  }

  int64_t              rev_edge_index_dims[1] = {int64_t(duplication_factor * total_num_bonds)};
  py::array_t<int64_t> rev_edge_index_array =
    py_array_from_array<int64_t>(std::move(rev_edge_index), rev_edge_index_dims, 1);
  // Prepare features for return
  std::vector<py::array> feature_arrays;

  feature_arrays.reserve(5);
  feature_arrays.push_back(atom_features_array);
  feature_arrays.push_back(bond_features_array);
  feature_arrays.push_back(edge_index_array);
  feature_arrays.push_back(rev_edge_index_array);
  feature_arrays.push_back(batch_array);
  return feature_arrays;
}

// Creates an RWMol from a SMILES string.
// See the declaration in features.h for more details.
std::unique_ptr<RDKit::RWMol> parse_mol(const std::string& smiles_string, bool explicit_H, bool ordered) {
  // Parse SMILES string with default options
  RDKit::SmilesParserParams     params;
  std::unique_ptr<RDKit::RWMol> mol{RDKit::SmilesToMol(smiles_string, params)};
  if (!mol) {
    return mol;
  }

  if (ordered) {
    // Do not order atoms to the canonical order.
    // Order them based only on the atom map, and only
    // if they indicate a valid order.
    const unsigned int        num_atoms = mol->getNumAtoms();
    std::vector<unsigned int> atom_order(num_atoms);
    unsigned int              map_num_min = std::numeric_limits<unsigned int>::max();
    unsigned int              map_num_max = 0;
    for (unsigned int i = 0; i < num_atoms; ++i) {
      RDKit::Atom* atom = mol->getAtomWithIdx(i);
      if (!atom->hasProp(RDKit::common_properties::molAtomMapNumber)) {
        ordered = false;
        // Don't break, because the property needs to be cleared
        // from any following atoms that might have it.
      } else {
        atom_order[i] = (unsigned int)atom->getAtomMapNum();
        if (atom_order[i] < map_num_min) {
          map_num_min = atom_order[i];
        }
        if (atom_order[i] > map_num_max) {
          map_num_max = atom_order[i];
        }

        // Clear the property, so that any equivalent molecules will
        // get the same canoncial order.
        atom->clearProp(RDKit::common_properties::molAtomMapNumber);
      }
    }

    // Must be 0-based or 1-based, and must be a complete permutation
    if (ordered && num_atoms != 0 && (map_num_min > 1 || map_num_max - map_num_min + 1 != num_atoms)) {
      ordered = false;
    }

    if (ordered) {
      // Invert the order
      // Use max value as a "not found yet" value
      constexpr unsigned int    not_found_value = std::numeric_limits<unsigned int>::max();
      std::vector<unsigned int> inverse_order(num_atoms, not_found_value);
      for (unsigned int i = 0; i < num_atoms; ++i) {
        unsigned int index = atom_order[i] - map_num_min;
        // Can't have the same index twice
        if (inverse_order[index] != not_found_value) {
          ordered = false;
          break;
        }
        inverse_order[index] = i;
      }

      if (ordered) {
        // Reorder the atoms to the explicit order
        mol.reset(static_cast<RDKit::RWMol*>(RDKit::MolOps::renumberAtoms(*mol, inverse_order)));
      }
    }
  }
  if (explicit_H) {
    RDKit::MolOps::addHs(*mol);
  } else {
    // Default params for SmilesToMol already calls removeHs,
    // and calling it again shouldn't have any net effect.
    // RDKit::MolOps::removeHs(*mol);
  }
  return mol;
}