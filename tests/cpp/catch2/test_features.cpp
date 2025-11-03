// SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

#include <pybind11/embed.h>  // Add this for py::scoped_interpreter
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../../../src/features.h"
#include "../../../src/one_hot.h"

namespace py = pybind11;

// Helper function to create a sequential array like numpy.arange
py::array_t<int64_t> create_sequential_array(int64_t size) {
  std::unique_ptr<int64_t[]> data = std::make_unique<int64_t[]>(size);
  for (int64_t i = 0; i < size; ++i) {
    data[i] = i;
  }
  const int64_t dims[1] = {size};
  return py_array_from_array(std::move(data), dims, 1);
}

template <typename T> T sum_py_array(const py::array_t<T>& arr) {
  T    sum = T{0};                         // Initialize to zero of type T
  auto buf = arr.template unchecked<1>();  // Add 'template' keyword
  for (ssize_t i = 0; i < buf.shape(0); ++i) {
    sum += buf(i);
  }
  return sum;
}

template <typename T> py::array_t<T> abs_py_array(const py::array_t<T>& arr) {
  auto                 buf  = arr.template unchecked<1>();
  ssize_t              n    = buf.shape(0);
  std::unique_ptr<T[]> data = std::make_unique<T[]>(n);
  for (ssize_t i = 0; i < n; ++i) {
    data[i] = std::abs(buf(i));
  }
  const int64_t dims[1] = {static_cast<int64_t>(n)};
  return py_array_from_array(std::move(data), dims, 1);
}

TEST_CASE("Molecule parsing works correctly", "[utils]") {
  std::string smiles = R"(CC1=C(CC(=O)O)c2cc(F)ccc2/C1=C\c1ccc([S+](C)[O-])cc1)";

  SECTION("With implicit hydrogens") {
    bool                          explicit_H = false;
    std::unique_ptr<RDKit::RWMol> mol        = parse_mol(smiles, explicit_H);

    REQUIRE(mol != nullptr);
    CHECK(mol->getNumAtoms() == 25);
    CHECK(mol->getNumBonds() == 27);
  }

  SECTION("With explicit hydrogens") {
    bool                          explicit_H = true;
    std::unique_ptr<RDKit::RWMol> mol_with_h = parse_mol(smiles, explicit_H);

    REQUIRE(mol_with_h != nullptr);
    CHECK(mol_with_h->getNumAtoms() == 42);
  }
}

TEST_CASE("Feature name to array conversion", "[features]") {
  py::scoped_interpreter guard{};
  SECTION("Atom float feature names") {
    std::vector<std::string> atom_float_feature_names = {
      "atomic-number",
      "mass",
      "valence",
      "implicit-valence",
      "hybridization",
      "chirality",
      "aromatic",
      "in-ring",
      "min-ring",
      "max-ring",
      "num-ring",
      "degree",
      "radical-electron",
      "formal-charge",
      "group",
      "period",
      "single-bond",
      "aromatic-bond",
      "double-bond",
      "triple-bond",
      "is-carbon",
      "hydrogen-bond-donor",
      "hydrogen-bond-acceptor",
      "acidic",
      "basic",
      "unknown-placeholder-feature"  // represents unknown feature
    };

    py::array_t<int64_t> atom_float_feature_array = atom_float_feature_names_to_array(atom_float_feature_names);

    py::array_t<int64_t> atom_float_feature_array_ref = create_sequential_array(26);
    int64_t              atom_float_diff =
      sum_py_array(abs_py_array(py::array_t<int64_t>(atom_float_feature_array - atom_float_feature_array_ref)));

    CHECK(atom_float_diff == 0);
  }

  SECTION("Atom onehot feature names") {
    std::vector<std::string> atom_onehot_feature_names = {
      "atomic-number",
      "atomic-number-common",
      "atomic-number-organic",
      "degree",
      "total-degree",
      "valence",
      "implicit-valence",
      "hybridization",
      "hybridization-expanded",
      "hybridization-organic",
      "chirality",
      "group",
      "period",
      "formal-charge",
      "num-hydrogens",
      "ring-size",
      "unknown-placeholder-feature"  // represents unknown feature
    };

    py::array_t<int64_t> atom_onehot_feature_array     = atom_onehot_feature_names_to_array(atom_onehot_feature_names);
    py::array_t<int64_t> atom_onehot_feature_array_ref = create_sequential_array(17);
    int                  atom_onehot_diff =
      sum_py_array(abs_py_array(py::array_t<int64_t>(atom_onehot_feature_array - atom_onehot_feature_array_ref)));

    CHECK(atom_onehot_diff == 0);
  }

  SECTION("Bond feature names") {
    std::vector<std::string> bond_feature_names = {
      "is-null",
      "bond-type-float",
      "bond-type-onehot",
      "in-ring",
      "conjugated",
      "stereo",
      "conformer-bond-length",
      "estimated-bond-length",
      "unknown-placeholder-feature"  // represents unknown feature
    };

    py::array_t<int64_t> bond_feature_array     = bond_feature_names_to_array(bond_feature_names);
    py::array_t<int64_t> bond_feature_array_ref = create_sequential_array(9);
    int                  bond_feature_diff =
      sum_py_array(abs_py_array(py::array_t<int64_t>(bond_feature_array - bond_feature_array_ref)));

    CHECK(bond_feature_diff == 0);
  }
}

TEST_CASE("Feature size calculation", "[features]") {
  SECTION("Atom onehot feature sizes") {
    // Test atomic number feature size
    size_t atomic_num_size = get_one_hot_atom_feature_size(AtomOneHotFeature::ATOMIC_NUM);
    CHECK(atomic_num_size == 101);

    // Test degree feature size
    size_t degree_size = get_one_hot_atom_feature_size(AtomOneHotFeature::DEGREE);
    CHECK(degree_size == 6);

    // Test formal charge feature size
    size_t formal_charge_size = get_one_hot_atom_feature_size(AtomOneHotFeature::FORMAL_CHARGE);
    CHECK(formal_charge_size == 6);
  }

  SECTION("Bond onehot feature sizes") {
    // Test bond type feature size
    size_t bond_type_size = get_one_hot_bond_feature_size(BondFeature::TYPE_ONE_HOT);
    CHECK(bond_type_size == 4);

    // Test bond stereo feature size
    size_t bond_stereo_size = get_one_hot_bond_feature_size(BondFeature::STEREO_ONE_HOT);
    CHECK(bond_stereo_size == 7);
  }
}

TEST_CASE("Multiple feature sizes", "[features]") {
  py::scoped_interpreter guard{};
  SECTION("Atom dimension calculation") {
    // Setup for atom dimension calculation
    std::vector<std::string> atom_onehot_feature_names =
      {"atomic-number", "total-degree", "formal-charge", "chirality", "num-hydrogens", "hybridization"};
    py::array_t<int64_t>     atom_onehot_feature_array = atom_onehot_feature_names_to_array(atom_onehot_feature_names);
    std::vector<std::string> atom_float_feature_names  = {"aromatic", "mass"};
    py::array_t<int64_t>     atom_float_feature_array  = atom_float_feature_names_to_array(atom_float_feature_names);

    // Calculate and check the dimension
    size_t atom_dim = compute_atom_dim(atom_onehot_feature_array, atom_float_feature_array);
    CHECK(atom_dim == 133);
  }

  SECTION("Bond dimension calculation") {
    // Setup for bond dimension calculation
    std::vector<std::string> bond_feature_names = {"is-null", "bond-type-onehot", "conjugated", "in-ring", "stereo"};
    py::array_t<int64_t>     bond_feature_array = bond_feature_names_to_array(bond_feature_names);

    // Calculate and check the dimension
    size_t bond_dim = compute_bond_dim(bond_feature_array);
    CHECK(bond_dim == 14);
  }
}
