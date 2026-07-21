# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved. # noqa: E501
# SPDX-License-Identifier: Apache-2.0

"""Tests for the atom-map based reordering performed by `parse_mol`.

Atom maps ("bookmarks" in RDKit) are honored only when they form a complete
permutation of the atoms, numbered either 0-based or 1-based. Anything else is
ignored, leaving the atoms in the order the SMILES parser produced them.
"""

import numpy as np
import pytest

import cuik_molmaker


def feature_arrays():
    atom_onehot_property_list = cuik_molmaker.atom_onehot_feature_names_to_array(
        ["atomic-number-common"]
    )
    atom_float_property_list = cuik_molmaker.atom_float_feature_names_to_array(
        ["atomic-number"]
    )
    bond_property_list = cuik_molmaker.bond_feature_names_to_array(["bond-type-float"])
    return atom_onehot_property_list, atom_float_property_list, bond_property_list


def atom_features(smiles, **kwargs):
    """Atom feature rows for `smiles`, one row per atom, in featurizer order.

    Extra keyword arguments are forwarded to `mol_featurizer`, so passing
    nothing exercises the default value of `ordered`.
    """
    explicit_H, offset_carbon, duplicate_edges, add_self_loop = (
        False,
        False,
        False,
        False,
    )
    atom_feats = cuik_molmaker.mol_featurizer(
        smiles,
        *feature_arrays(),
        explicit_H,
        offset_carbon,
        duplicate_edges,
        add_self_loop,
        **kwargs,
    )[0]
    return atom_feats


def batch_atom_features(smiles_list, **kwargs):
    """Same as `atom_features`, but through `batch_mol_featurizer`."""
    explicit_H, offset_carbon, duplicate_edges, add_self_loop = (
        False,
        False,
        False,
        False,
    )
    atom_feats = cuik_molmaker.batch_mol_featurizer(
        smiles_list,
        *feature_arrays(),
        explicit_H,
        offset_carbon,
        duplicate_edges,
        add_self_loop,
        **kwargs,
    )[0]
    return atom_feats


def test_atom_order_is_observable():
    """Guard for the other tests: writing the atoms in a different order in the
    SMILES really does produce differently ordered features."""
    assert not np.array_equal(atom_features("CO"), atom_features("OC"))


@pytest.mark.parametrize(
    "smiles",
    [
        "[CH3:0][OH:1]",  # 0-based, already in order
        "[OH:1][CH3:0]",  # 0-based, reversed
        "[CH3:1][OH:2]",  # 1-based, already in order
        "[OH:2][CH3:1]",  # 1-based, reversed
    ],
)
def test_complete_map_reorders_atoms(smiles):
    """A complete 0- or 1-based atom map puts the atoms in the mapped order,
    regardless of the order they appear in the SMILES."""
    np.testing.assert_allclose(atom_features(smiles), atom_features("CO"))


def test_zero_and_one_based_maps_agree():
    np.testing.assert_allclose(
        atom_features("[OH:1][CH3:0]"), atom_features("[OH:2][CH3:1]")
    )


@pytest.mark.parametrize(
    "smiles",
    [
        "[CH3:3][NH:2][OH:1]",  # 1-based, fully reversed
        "[CH3:2][NH:1][OH:0]",  # 0-based, fully reversed
    ],
)
def test_complete_map_applies_arbitrary_permutation(smiles):
    """Three atoms, mapped so the featurizer order is O, N, C."""
    np.testing.assert_allclose(atom_features(smiles), atom_features("ONC"))


@pytest.mark.parametrize(
    "smiles",
    [
        "OC",  # no atom maps at all
        "[OH:1]C",  # only some atoms mapped
        "[OH:4][CH3:3]",  # complete permutation, but neither 0- nor 1-based
        "[OH:1][CH3:1]",  # duplicated map number
        "[OH:1][CH3:5]",  # 1-based but not contiguous
        "[OH:0][CH3:5]",  # 0-based but not contiguous
    ],
)
def test_invalid_map_leaves_atoms_in_parsed_order(smiles):
    """Maps that are not a complete 0- or 1-based permutation are ignored, so
    the atoms stay in the order the SMILES string lists them (O, then C)."""
    np.testing.assert_allclose(atom_features(smiles), atom_features("OC"))


def test_map_numbers_do_not_leak_into_features():
    """The map numbers are stripped, so a mapped molecule and the equivalent
    unmapped one featurize identically."""
    np.testing.assert_allclose(atom_features("[CH3:1][OH:2]"), atom_features("CO"))


@pytest.mark.parametrize("smiles", ["[OH:1][CH3:0]", "[OH:2][CH3:1]"])
def test_ordered_false_disables_reordering(smiles):
    """With ordered=False the atom maps are ignored entirely, leaving the atoms
    in the order the SMILES string lists them (O, then C)."""
    np.testing.assert_allclose(
        atom_features(smiles, ordered=False), atom_features("OC")
    )


@pytest.mark.parametrize("smiles", ["[OH:1][CH3:0]", "[OH:2][CH3:1]"])
def test_ordered_defaults_to_true(smiles):
    np.testing.assert_allclose(
        atom_features(smiles), atom_features(smiles, ordered=True)
    )


def test_ordered_accepted_positionally():
    """`ordered` is the ninth positional argument, after add_self_loop."""
    np.testing.assert_allclose(
        cuik_molmaker.mol_featurizer(
            "[OH:2][CH3:1]", *feature_arrays(), False, False, False, False, False
        )[0],
        atom_features("OC"),
    )


@pytest.mark.parametrize("ordered", [True, False])
def test_batch_featurizer_honors_ordered(ordered):
    expected = "CO" if ordered else "OC"
    np.testing.assert_allclose(
        batch_atom_features(["[OH:2][CH3:1]"], ordered=ordered),
        atom_features(expected),
    )
