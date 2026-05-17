// SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

//! @file This header file declares feature-related enums, functions, and structs,
//!       some of which are defined in features.cpp and exported to Python.

#pragma once

#include <stdint.h>

#include <cmath>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

// RDKit headers
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

// PyBind headers
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "export.h"

namespace py = pybind11;

//! Levels at which features or labels can be associated
//! String names are in `feature_level_to_enum` in features.cpp
enum class FeatureLevel {
  NODE,      //!< Values for each node (atom)
  EDGE,      //!< Values for each edge (bond)
  NODEPAIR,  //!< Values for each pair of nodes (pair of atoms), even if no edge (bond)
  GRAPH      //!< Values for whole molecule
};

//! Features for use by `get_atom_float_feature` in float_features.cpp
//! String names are in `atom_float_name_to_enum` in features.cpp
enum class AtomFloatFeature {
  ATOMIC_NUMBER,
  MASS,
  VALENCE,
  IMPLICIT_VALENCE,
  HYBRIDIZATION,
  CHIRALITY,
  AROMATIC,
  IN_RING,
  MIN_RING,
  MAX_RING,
  NUM_RING,
  DEGREE,
  RADICAL_ELECTRON,
  FORMAL_CHARGE,
  GROUP,
  PERIOD,
  SINGLE_BOND,
  AROMATIC_BOND,
  DOUBLE_BOND,
  TRIPLE_BOND,
  IS_CARBON,
  HYDROGEN_BOND_DONOR,
  HYDROGEN_BOND_ACCEPTOR,
  ACIDIC,
  BASIC,
  UNKNOWN
};

//! Features for use by `get_one_hot_atom_feature` in one_hot.cpp
//! String names are in `atom_onehot_name_to_enum` in features.cpp
enum class AtomOneHotFeature {
  ATOMIC_NUM,              //!< All atomic numbers from 1 - 100
  ATOMIC_NUM_COMMON,       //!< First 4 rows of periodic table and Iodine
  ATOMIC_NUM_ORGANIC,      //!< Organic chemistry elements only
  DEGREE,                  //!< Number of explicit neighboring atoms
  TOTAL_DEGREE,            //!< Number of neighboring atoms including hydrogens
  VALENCE,                 //!< Total valence of the atom
  IMPLICIT_VALENCE,        //!< Implicit valence of the atom
  HYBRIDIZATION,           //!< Hybridizations specified in `hybridizationList` in one_hot.cpp
  HYBRIDIZATION_EXPANDED,  //!< Hybridizations expanded to include all possible values
  HYBRIDIZATION_ORGANIC,   //!< Hybridizations of organic elements
  CHIRALITY,               //!< "R", anything other value ("S") or no value, and an extra
                           //!< chirality-related value (independent of the other two, so can
                           //!< have a 2nd one value)
  GROUP,                   //!< Specified by `atomicNumToGroupTable` in float_features.h
  PERIOD,                  //!< Specified by `atomicNumToPeriodTable` in float_features.h
  FORMAL_CHARGE,           //!< Formal charge on atom
  NUM_HYDROGENS,           //!< Total number of hydrogens (explicit and implicit) on an atom
  RING_SIZE,               //!< Number of rings the atom is in
  UNKNOWN                  //!< Sentinel value.  Do not use.
};

//! Features for use by `get_one_hot_bond_feature` in one_hot.cpp (if ends in `ONE_HOT`), and
//! `get_bond_float_feature` in float_features.cpp
//! String names are in `bond_name_to_enum` in features.cpp
enum class BondFeature {
  IS_NULL,                //!< 1 if Bond is nullptr else 0
  TYPE_FLOAT,             //!< Bond type as a float, e.g. 2.0 for double, 1.5 for aromatic
  TYPE_ONE_HOT,           //!< Selected bond types specified in `bondTypeList` in one_hot.cpp
  IN_RING,                //!< 1.0 if the bond is in at least one ring, else 0.0
  CONJUGATED,             //!< 1.0 if the bond is conjugated, else 0.0
  STEREO_ONE_HOT,         //!< Selected bond stereo values specified in `bondStereoList` in
                          //!< one_hot.cpp
  CONFORMER_BOND_LENGTH,  //!< Length of the bond from a conformer (either first or computed)
  ESTIMATED_BOND_LENGTH,  //!< Length of the bond estimated with a fast heuristic
  UNKNOWN                 //!< Sentinel value.  Do not use.
};

//! Options for handling NaN or infinite values, passed from Python to `mol_featurizer` in
//! features.cpp.  Masking is done in `mask_nans` in features.h
enum class MaskNaNStyle {
  NONE,    //!< Ignore (keep) NaN values
  REPORT,  //!< (default behaviour) Count NaN values and report that with the index of the
           //!< first tensor that contained NaNs
  REPLACE  //!< Replace NaN values with a specific value (defaults to zero)
};

//! Class to help supporting `int16_t` as if it's a 16-bit floating-point (FP16) type,
//! while still supporting `float` (FP32) and `double` (FP64).
template <typename T> struct FeatureValues {};

//! Explicit instantiation of `FeatureValues` for `int16_t` as if it's a 16-bit
//! floating-point (FP16) type.
template <> struct FeatureValues<int16_t> {
  static constexpr int16_t zero      = 0x0000;
  static constexpr int16_t one       = 0x3C00;
  static constexpr int16_t nan_value = 0x7C01;

  template <typename T> static int16_t convertToFeatureType(T inputType) {
    static_assert(std::is_floating_point_v<T>);
    return float32_to_fp16_ieee(float(inputType));
  }

  static constexpr bool is_finite(int16_t v) {
    // If the exponent bits are the maximum value, v is infinite or NaN
    return (v & 0x7C00) != 0x7C00;
  }

  using MathType = float;

 private:
  // Manual IEEE 754 half-precision conversion (replacing c10 dependency)
  static int16_t float32_to_fp16_ieee(float f) {
    union {
      float    f;
      uint32_t i;
    } u = {f};

    uint32_t sign     = (u.i >> 16) & 0x8000;               // Sign bit
    int32_t  exp      = int32_t((u.i >> 23) & 0xff) - 127;  // Exponent (remove bias)
    uint32_t mantissa = u.i & 0x7fffff;                     // Mantissa (23 bits)

    // Handle special cases
    if (exp == 128) {  // Infinity or NaN
      return static_cast<int16_t>(sign | 0x7c00 | (mantissa ? 0x0200 : 0));
    }

    if (exp > 15) {  // Overflow -> Infinity
      return static_cast<int16_t>(sign | 0x7c00);
    }

    if (exp < -14) {  // Underflow -> Zero or denormal
      if (exp < -25)
        return static_cast<int16_t>(sign);  // Too small -> zero

      // Denormalized number
      mantissa |= 0x800000;  // Add implicit leading 1
      int shift = -14 - exp;

      // Apply proper rounding considering all discarded bits from both shifts
      // Total shift needed: shift + 13 (to go from 23-bit to 10-bit mantissa)
      int      total_shift  = shift + 13;
      uint32_t rounding_bit = (mantissa >> (total_shift - 1)) & 1;         // Bit at position total_shift-1
      uint32_t lsb          = (mantissa >> total_shift) & 1;               // LSB of result
      uint32_t sticky_bits  = mantissa & ((1u << (total_shift - 1)) - 1);  // All bits below rounding bit

      // Round-to-even logic
      uint32_t round_up = rounding_bit && (lsb || (sticky_bits != 0));

      return static_cast<int16_t>(sign | ((mantissa >> total_shift) + round_up));
    }

    // Normal number
    return static_cast<int16_t>(
      sign | ((exp + 15) << 10) |
      ((mantissa >> 13)  // Shift by 13 bits
                         // Bit 12 is rounding bit. If bit 12 is 0, we round down (no change to other bits)
                         // If bit 12 is 1, we are exactly halfway and rounding may be needed.

       /* Truth table for rounding:
       | bit 13 | sticky bits | Round up? | Explanation |
       -------------------------------------------------
       |  0     |   non-zero  |    Yes    | More than half way |
       |  0     |   zero      |    No     | Exactly halfway - so no round to even |
       |  1     |   non-zero  |    Yes    | FP16 result is odd, so round up   |
       |  1     |   zero      |    Yes     | Exactly halfway, round up to even |
       This truth table is same as OR operation between bit 13 and non-zero sticky bits
       Code for non-zero sticky bits is ((mantissa & 0xFFF) != 0)
       */
       + (((mantissa >> 12) & 1) &&  // Round if bit 12 is 1
          (((mantissa >> 13) & 1) || ((mantissa & 0xFFF) != 0)))));
  }
};
//! Explicit instantiation of `FeatureValues` for `float` (FP32)
template <> struct FeatureValues<float> {
  static constexpr float zero      = 0.0f;
  static constexpr float one       = 1.0f;
  static constexpr float nan_value = std::numeric_limits<float>::quiet_NaN();

  template <typename T> static float convertToFeatureType(T inputType) {
    static_assert(std::is_floating_point_v<T>);
    return float(inputType);
  }

  static bool is_finite(float v) { return std::isfinite(v); }

  using MathType = float;
};
//! Explicit instantiation of `FeatureValues` for `double` (FP64)
template <> struct FeatureValues<double> {
  static constexpr double zero      = 0.0;
  static constexpr double one       = 1.0;
  static constexpr double nan_value = std::numeric_limits<double>::quiet_NaN();

  template <typename T> static double convertToFeatureType(T inputType) {
    static_assert(std::is_floating_point_v<T>);
    return double(inputType);
  }

  // Note: std::isfinite is not constexpr in MSVC, so we can't make this constexpr
  static inline bool is_finite(double v) { return std::isfinite(v); }

  using MathType = double;
};

//! Handling for NaN or infinite values in an array, `data`,  of `n` values.
//! @see MaskNaNStyle
template <typename T> constexpr int64_t mask_nans(T* data, size_t n, MaskNaNStyle style, T value) {
  if (style == MaskNaNStyle::NONE) {
    return 0;
  }
  if (style == MaskNaNStyle::REPLACE) {
    for (size_t i = 0; i < n; ++i) {
      if (!FeatureValues<T>::is_finite(data[i])) {
        data[i] = value;
      }
    }
    return 0;
  }

  // assert(mask_nan_style == MaskNaNStyle::REPORT);
  int64_t num_nans = 0;
  for (size_t i = 0; i < n; ++i) {
    num_nans += (!FeatureValues<T>::is_finite(data[i]));
  }
  return num_nans;
}

// This is just a function to provide to torch, so that we don't have to copy
// the tensor data to put it in a torch tensor, and torch can delete the data
// when it's no longer needed.
template <typename T> void deleter(void* p) {
  delete[] (T*)p;
}

//! Helper function to construct a pybind11 `py::array_t` from a C++ array.
//! The `py::array_t` takes ownership of the memory owned by `source`.
template <typename T>
py::array_t<T> py_array_from_array(std::unique_ptr<T[]>&& source, const int64_t* dims, size_t num_dims) {
  // Convert dims to vector of sizes for pybind11
  std::vector<size_t> shape(num_dims);
  for (size_t i = 0; i < num_dims; ++i) {
    shape[i] = static_cast<size_t>(dims[i]);
  }

  // Create capsule that will handle deletion when array is destroyed
  T*          raw_ptr = source.release();
  py::capsule cleanup(raw_ptr, [](void* ptr) { delete[] static_cast<T*>(ptr); });

  // Create py::array_t with shape, data pointer, and cleanup capsule
  return py::array_t<T>(shape, raw_ptr, cleanup);
}

//! Most of the data needed about an atom
struct CompactAtom {
  uint8_t atomicNum;
  uint8_t totalDegree;
  int8_t  formalCharge;
  uint8_t chiralTag;
  uint8_t totalNumHs;
  uint8_t hybridization;
  bool    isAromatic;
  float   mass;
};

//! Most of the data needed about a bond
struct CompactBond {
  uint8_t  bondType;
  bool     isConjugated;
  bool     isInRing;
  uint8_t  stereo;
  uint32_t beginAtomIdx;
  uint32_t endAtomIdx;
};

//! Data representing a molecule before featurization
struct GraphData {
  const size_t                   num_atoms;
  std::unique_ptr<CompactAtom[]> atoms;
  const size_t                   num_bonds;
  std::unique_ptr<CompactBond[]> bonds;

  std::unique_ptr<RDKit::RWMol> mol;
};

//! Condensed Graph of Reaction featurization modes, matching chemprop's RxnMode enum
enum class ReactionMode {
  REAC_DIFF,          //!< First half = reactant feats; second half = prod - reac diff
  REAC_PROD,          //!< First half = reactant feats; second half = product feats
  PROD_DIFF,          //!< First half = product feats;  second half = prod - reac diff
  REAC_DIFF_BALANCE,  //!< Like REAC_DIFF but unmatched atoms copy own feats (diff = 0)
  REAC_PROD_BALANCE,  //!< Like REAC_PROD but unmatched atoms copy own feats
  PROD_DIFF_BALANCE,  //!< Like PROD_DIFF but unmatched atoms copy own feats
  UNKNOWN
};

//! Data representing a reaction (two molecules + atom correspondence) before CGR featurization.
//! Both GraphData members retain their RDKit mol pointers — required by one_hot.cpp features.
struct CompactReaction {
  GraphData reac;  //!< Reactant side (owns RDKit mol + CompactAtom/Bond caches)
  GraphData prod;  //!< Product side

  //! Atom mapping: reactant atom index → product atom index (built from atom-map numbers)
  std::unordered_map<uint32_t, uint32_t> r2p_idx_map;
  //! Inverse: product atom index → reactant atom index
  std::unordered_map<uint32_t, uint32_t> p2r_idx_map;

  //! Reactant atoms with no matching product atom (map num absent on product side)
  std::vector<uint32_t> reac_only_idxs;
  //! Product atoms with no matching reactant atom; these become CGR nodes n_reac..n_cgr-1
  std::vector<uint32_t> prod_only_idxs;

  //! Bond lookup for O(1) cross-referencing. Key = (min_atom_idx << 32) | max_atom_idx
  //! using the *side-local* (reactant or product) atom indices.
  std::unordered_map<uint64_t, uint32_t> reac_bond_lookup;
  std::unordered_map<uint64_t, uint32_t> prod_bond_lookup;
};

//! Computes the total dimension of atom features based on the property lists
CUIK_EXPORT size_t compute_atom_dim(const py::array_t<int64_t>& atom_property_list_onehot,
                                    const py::array_t<int64_t>& atom_property_list_float);

//! Computes the total dimension of bond features based on the property list
CUIK_EXPORT size_t compute_bond_dim(const py::array_t<int64_t>& bond_property_list);

//! This is called from Python to list atom one-hot features in a format that will be faster
//! to interpret inside `mol_featurizer`, passed in the `atom_property_list_onehot` parameter.
//! Implemented in features.cpp, but declared here so that cuik_molmaker_cpp.cpp can expose them to
//! Python via pybind.
CUIK_EXPORT py::array_t<int64_t> atom_onehot_feature_names_to_array(const std::vector<std::string>& features);

//! This is called from Python to list all atom one-hot features.
//! Implemented in features.cpp, but declared here so that cuik_molmaker_cpp.cpp can expose them to
//! Python via pybind.
CUIK_EXPORT std::vector<std::string> list_all_atom_onehot_features();

//! This is called from Python to list atom float features in a format that will be faster
//! to interpret inside `mol_featurizer`, passed in the `atom_property_list_float` parameter.
//! Implemented in features.cpp, but declared here so that cuik_molmaker_cpp.cpp can expose them to
//! Python via pybind.
CUIK_EXPORT py::array_t<int64_t> atom_float_feature_names_to_array(const std::vector<std::string>& features);

//! This is called from Python to list all atom float features.
//! Implemented in features.cpp, but declared here so that cuik_molmaker_cpp.cpp can expose them to
//! Python via pybind.
CUIK_EXPORT std::vector<std::string> list_all_atom_float_features();

//! This is called from Python to list bond features in a format that will be faster
//! to interpret inside `mol_featurizer`, passed in the `bond_property_list` parameter.
//! Implemented in features.cpp, but declared here so that cuik_molmaker_cpp.cpp can expose them to
//! Python via pybind.
CUIK_EXPORT py::array_t<int64_t> bond_feature_names_to_array(const std::vector<std::string>& features);

//! This is called from Python to list all bond features.
//! Implemented in features.cpp, but declared here so that cuik_molmaker_cpp.cpp can expose them to
//! Python via pybind.
CUIK_EXPORT std::vector<std::string> list_all_bond_features();

//! `mol_featurizer` is called from Python to get feature arrays for `smiles_string`.
//!
//! @param smiles_string SMILES string of the molecule to featurize
//! @param atom_property_list_onehot NumPy/pybind array returned by
//!                                  `atom_onehot_feature_names_to_array` representing the
//!                                  list of one-hot atom features to create.
//! @param atom_property_list_float NumPy/pybind array returned by
//!                                 `atom_float_feature_names_to_array` representing the
//!                                 list of float atom features to create.
//! @param bond_property_list NumPy/pybind array returned by `bond_feature_names_to_array`
//!                           representing the list of bond features to create.
//! @param explicit_H If true, implicit hydrogen atoms will be added explicitly
//!                   before featurizing.
//! @param duplicate_edges If true, bond features will have values stored for
//!                        both edge directions.
//! @param add_self_loop If true, bond features will have values stored for
//!                      self-edges.
//! @param offset_carbon If true, some atom float features will subtract a
//!                      value representing carbon, so that carbon atoms would have value zero.
//! @return A vector of torch NumPy/pybind arrays for the features.  The first array is the atom features
//!         array, `num_atoms` by the number of values required for all one-hot and float atom
//!         features.  The second array is the bond features array, `num_edges` (or
//!         `2*num_edges` if `duplicate_edges` is true) by the number of values required for all
//!         bond features. The third array is a 2 by `num_edges` (or `2*num_edges` if
//!         `duplicate_edges` is true) representing the indices of the nodes each edge is
//!         connected to. The fourth array is a 1 by `num_edges` array representing the reverse
//!         of the third array. The fifth array is a 1 by `num_atoms` array containing 0s.
CUIK_EXPORT std::vector<py::array> mol_featurizer(const std::string&          smiles_string,
                                                  const py::array_t<int64_t>& atom_property_list_onehot,
                                                  const py::array_t<int64_t>& atom_property_list_float,
                                                  const py::array_t<int64_t>& bond_property_list,
                                                  bool                        explicit_H,
                                                  bool                        offset_carbon,
                                                  bool                        duplicate_edges,
                                                  bool                        add_self_loop);

//! Creates an RWMol from a SMILES string.
//!
//! If `ordered` is true, and the string contains atom classes, called "bookmarks" in RDKit,
//! that form a complete (0-based) ordering of the atoms, the atoms will be reordered according
//! to this explicit order, and the bookmarks will be removed, so that canonical orders
//! can be correctly compared later.
//!
//! This is implemented in cuik_molmaker_cpp.cpp, but is declared in this header so
//! that both labels.cpp and features.cpp can call it.
std::unique_ptr<RDKit::RWMol> parse_mol(const std::string& smiles_string, bool explicit_H, bool ordered = true);

//! `batch_mol_featurizer` is called from Python to get feature arrays for `smiles_list`.
//!
//! @param smiles_list List of SMILES strings of molecules to featurize
//! @param atom_property_list_onehot NumPy/pybind array returned by
//!                                  `atom_onehot_feature_names_to_array` representing the
//!                                  list of one-hot atom features to create.
//! @param atom_property_list_float NumPy/pybind array returned by
//!                                 `atom_float_feature_names_to_array` representing the
//!                                 list of float atom features to create.
//! @param bond_property_list NumPy/pybind array returned by `bond_feature_names_to_array`
//!                           representing the list of bond features to create.
//! @param explicit_H If true, implicit hydrogen atoms will be added explicitly
//!                   before featurizing.
//! @param duplicate_edges If true, bond features will have values stored for
//!                        both edge directions.
//! @param offset_carbon If true, some atom float features will subtract a
//!                      value representing carbon, so that carbon atoms would have value zero.
//! @param add_self_loop If true, bond features will have values stored for
//!                      self-edges.
//! @return A vector of NumPy/pybind arrays for the features.  The first array is the atom features
//!         array, total number of atoms  by the number of values required for all one-hot and
//!         float atom features.  The second array is the bond features array, total number of
//!         edges (or `2*total_num_edges` if `duplicate_edges` is true) by the number of values
//!         required for all bond features. The third array is a 2 by `total_num_edges` (or
//!         `2*total_num_edges` if `duplicate_edges` is true) representing the indices of the
//!         nodes each edge is connected to. The fourth array is a 1 by `total_num_edges`
//!         array representing the reverse of the third array. The fifth array is a 1 by
//!         `total_num_atoms` array containing the index of the molecule each atom belongs to.
CUIK_EXPORT std::vector<py::array> batch_mol_featurizer(const std::vector<std::string>& smiles_list,
                                                        const py::array_t<int64_t>&     atom_property_list_onehot,
                                                        const py::array_t<int64_t>&     atom_property_list_float,
                                                        const py::array_t<int64_t>&     bond_property_list,
                                                        bool                            explicit_H,
                                                        bool                            offset_carbon,
                                                        bool                            duplicate_edges,
                                                        bool                            add_self_loop);

//! Parses one side of a reaction SMILES into an RWMol, preserving atom-map numbers.
//! Unlike parse_mol, this function does NOT clear atom-map numbers and does NOT reorder atoms.
//! @param keep_h  If true, SmilesParserParams.removeHs = false (retains explicit [H:n] atoms)
//! @param add_h   If true, RDKit::MolOps::addHs is called after parsing (adds unmapped Hs)
std::unique_ptr<RDKit::RWMol> parse_rxn_side_mol(const std::string& smiles, bool keep_h, bool add_h);

//! Parses a reaction SMILES pair into a CompactReaction (atom correspondence + both GraphData).
//! Both reac_smi and prod_smi must contain atom-map numbers.
//! keep_h / add_h semantics match chemprop's _ReactionDatapointMixin.from_smi exactly.
CUIK_EXPORT CompactReaction parse_reaction(const std::string& reac_smi,
                                           const std::string& prod_smi,
                                           bool               keep_h,
                                           bool               add_h);

//! Converts reaction mode name strings to a NumPy int64 array (mirrors atom_onehot_feature_names_to_array).
CUIK_EXPORT py::array_t<int64_t> reaction_mode_names_to_array(const std::vector<std::string>& modes);

//! Featurizes a batch of reactions as Condensed Graphs of Reaction (CGR).
//! Mirrors batch_mol_featurizer in interface and return convention (5 arrays).
//! @param reac_smiles_list  Reactant SMILES (atom-mapped); parallel to prod_smiles_list
//! @param prod_smiles_list  Product SMILES (atom-mapped)
//! @param keep_h  If true, retain explicit mapped [H:n] atoms (required for E2/SN2)
//! @param add_h   If true, add unmapped Hs via RDKit::MolOps::addHs (after parsing)
//! @param mode    CGR featurization mode (which combination of reac/prod/diff)
//! @return 5 arrays: [atom_feats, bond_feats, edge_index, rev_edge_index, batch]
CUIK_EXPORT std::vector<py::array> batch_reaction_featurizer(const std::vector<std::string>& reac_smiles_list,
                                                             const std::vector<std::string>& prod_smiles_list,
                                                             const py::array_t<int64_t>&     atom_property_list_onehot,
                                                             const py::array_t<int64_t>&     atom_property_list_float,
                                                             const py::array_t<int64_t>&     bond_property_list,
                                                             bool                            keep_h,
                                                             bool                            add_h,
                                                             bool                            offset_carbon,
                                                             ReactionMode                    mode);
