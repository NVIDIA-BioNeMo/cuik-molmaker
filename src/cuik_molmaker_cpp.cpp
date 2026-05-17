// SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

//! @file This file specifies which functions are exported to Python,
//!       as well as defining `parse_mol` and `get_canonical_atom_order`,
//!       declared in features.h and called from features.cpp and labels.cpp

#include "features.h"

// C++ standard library headers
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include <filesystem>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

// RDKit headers
#include <GraphMol/Atom.h>
#include <GraphMol/Canon.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/types.h>

// PyBind headers for use by library to be imported by Python
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

// This is necessary to export Python functions in a Python module named cuik_molmaker.
PYBIND11_MODULE(cuik_molmaker_cpp, m) {
  m.doc() = "Cuik MolMaker C++ plugin";  // Python module docstring

  // Functions in features.cpp
  m.def("atom_onehot_feature_names_to_array",
        &atom_onehot_feature_names_to_array,
        "Accepts feature names and returns a NumPy array representing them as integers");
  m.def("atom_float_feature_names_to_array",
        &atom_float_feature_names_to_array,
        "Accepts feature names and returns a NumPy array representing them as integers");
  m.def("bond_feature_names_to_array",
        &bond_feature_names_to_array,
        "Accepts feature names and returns a NumPy array representing them as integers");
  m.def(
    "mol_featurizer",
    &mol_featurizer,
    "Accepts a SMILES string and returns a list of NumPy arrays representing atom and bond features of the molecule.");
  m.def(
    "batch_mol_featurizer",
    &batch_mol_featurizer,
    "Accepts a list of SMILES strings and returns a list of NumPy arrays representing atom and bond features of the molecules.");
  m.def("list_all_atom_onehot_features",
        &list_all_atom_onehot_features,
        "Returns a list of all atom one-hot features.");

  m.def("list_all_atom_float_features", &list_all_atom_float_features, "Returns a list of all atom float features.");

  m.def("list_all_bond_features", &list_all_bond_features, "Returns a list of all bond features.");

  // Reaction featurization (CGR)
  m.def("reaction_mode_names_to_array",
        &reaction_mode_names_to_array,
        "Convert reaction mode name strings to int64 array (mirrors atom_onehot_feature_names_to_array).");
  m.def(
    "batch_reaction_featurizer",
    [](const std::vector<std::string>& reac_smiles_list,
       const std::vector<std::string>& prod_smiles_list,
       const py::array_t<int64_t>&     atom_property_list_onehot,
       const py::array_t<int64_t>&     atom_property_list_float,
       const py::array_t<int64_t>&     bond_property_list,
       bool                            keep_h,
       bool                            add_h,
       bool                            offset_carbon,
       int64_t                         mode_int) {
      return batch_reaction_featurizer(reac_smiles_list,
                                       prod_smiles_list,
                                       atom_property_list_onehot,
                                       atom_property_list_float,
                                       bond_property_list,
                                       keep_h,
                                       add_h,
                                       offset_carbon,
                                       ReactionMode(mode_int));
    },
    "Featurize a batch of reactions (CGR) and return 5 NumPy arrays.");
}
