# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved. # noqa: E501
# SPDX-License-Identifier: Apache-2.0

"""Tests for the `keep_h` and `ignore_stereo` parsing options on
`mol_featurizer` / `batch_mol_featurizer`, and `ignore_stereo` on
`batch_reaction_featurizer`.
"""

import numpy as np
import pytest

import cuik_molmaker


def featurize(
    smiles,
    explicit_H=False,
    keep_h=False,
    ignore_stereo=False,
    ordered=True,
    atom_onehot=("atomic-number-common",),
    atom_float=(),
    bond=(),
):
    """mol_featurizer output for `smiles`, with the new options exposed as
    keyword arguments and everything else defaulted."""
    oh = cuik_molmaker.atom_onehot_feature_names_to_array(list(atom_onehot))
    fl = cuik_molmaker.atom_float_feature_names_to_array(list(atom_float))
    bp = cuik_molmaker.bond_feature_names_to_array(list(bond))
    return cuik_molmaker.mol_featurizer(
        smiles,
        oh,
        fl,
        bp,
        explicit_H=explicit_H,
        offset_carbon=False,
        duplicate_edges=False,
        add_self_loop=False,
        ordered=ordered,
        keep_h=keep_h,
        ignore_stereo=ignore_stereo,
    )


# ---------------------------------------------------------------------------
# keep_h
# ---------------------------------------------------------------------------


def test_keep_h_defaults_to_false():
    """Without keep_h, explicit [H] atoms written in the SMILES are folded
    back into implicit valence on the heavy atoms (C, O), leaving 2 atoms."""
    feats = featurize("[H]C([H])([H])O")
    assert feats[0].shape[0] == 2


def test_keep_h_true_retains_explicit_atoms():
    """With keep_h, the 3 explicit [H] on carbon are kept as real atoms; the
    implicit O-H is not (nothing adds it), giving 5 atoms total."""
    feats = featurize("[H]C([H])([H])O", keep_h=True)
    assert feats[0].shape[0] == 5


def test_keep_h_and_explicit_H_converge():
    """Whether the SMILES starts with some explicit Hs (keep_h) or none, once
    explicit_H adds every remaining implicit H, methanol always has 6 atoms."""
    default_then_add = featurize("[H]C([H])([H])O", explicit_H=True)
    keep_then_add = featurize("[H]C([H])([H])O", explicit_H=True, keep_h=True)
    assert default_then_add[0].shape[0] == 6
    assert keep_then_add[0].shape[0] == 6


def test_batch_mol_featurizer_honors_keep_h():
    oh = cuik_molmaker.atom_onehot_feature_names_to_array(["atomic-number-common"])
    fl = cuik_molmaker.atom_float_feature_names_to_array([])
    bp = cuik_molmaker.bond_feature_names_to_array([])
    feats = cuik_molmaker.batch_mol_featurizer(
        ["[H]C([H])([H])O"],
        oh,
        fl,
        bp,
        explicit_H=False,
        offset_carbon=False,
        duplicate_edges=False,
        add_self_loop=False,
        keep_h=True,
    )
    assert feats[0].shape[0] == 5


# ---------------------------------------------------------------------------
# ignore_stereo
# ---------------------------------------------------------------------------


def chirality_features(smiles, ignore_stereo):
    return featurize(
        smiles,
        atom_onehot=("chirality",),
        atom_float=("chirality",),
        ignore_stereo=ignore_stereo,
    )[0]


def test_ignore_stereo_clears_onehot_and_float_chirality():
    """The stereocenter is atom index 1 (N, C@@H, CH3, C(=O), O, O). Clearing
    stereo must change BOTH the one-hot chirality block (reads GetChiralTag)
    and the float chirality feature (reads the _CIPCode property) -- this is
    exactly the gap a manual port of Chemprop's SetChiralTag-only loop would
    miss, which is why ignore_stereo is implemented with
    RDKit::MolOps::removeStereochemistry instead."""
    with_stereo = chirality_features("N[C@@H](C)C(=O)O", ignore_stereo=False)
    without_stereo = chirality_features("N[C@@H](C)C(=O)O", ignore_stereo=True)
    plain = chirality_features("NC(C)C(=O)O", ignore_stereo=False)

    assert not np.allclose(with_stereo[1], without_stereo[1])
    np.testing.assert_allclose(without_stereo, plain)


def bond_stereo_features(smiles, ignore_stereo):
    return featurize(smiles, bond=("stereo",), ignore_stereo=ignore_stereo)[1]


def test_ignore_stereo_clears_bond_stereo():
    """cis and trans difluoroethylene must differ in bond features by
    default, and become identical (and match the unspecified form) once
    ignore_stereo strips E/Z information."""
    cis = bond_stereo_features(r"F/C=C\F", ignore_stereo=False)
    trans = bond_stereo_features("F/C=C/F", ignore_stereo=False)
    cis_cleared = bond_stereo_features(r"F/C=C\F", ignore_stereo=True)
    trans_cleared = bond_stereo_features("F/C=C/F", ignore_stereo=True)
    plain = bond_stereo_features("FC=CF", ignore_stereo=False)

    assert not np.allclose(cis, trans)
    np.testing.assert_allclose(cis_cleared, trans_cleared)
    np.testing.assert_allclose(cis_cleared, plain)


# ---------------------------------------------------------------------------
# ignore_stereo on the reaction (CGR) path
# ---------------------------------------------------------------------------


def test_batch_reaction_featurizer_ignore_stereo_clears_chirality():
    """Identity reaction (product == reactant) on an atom-mapped molecule
    with one stereocenter (atom map 2). ignore_stereo must change the CGR
    atom features, since it strips the stereocenter on both reactant and
    product sides."""
    reac = ["[NH2:1][C@@H:2]([CH3:3])[C:4](=[O:5])[OH:6]"]
    prod = ["[NH2:1][C@@H:2]([CH3:3])[C:4](=[O:5])[OH:6]"]
    oh = cuik_molmaker.atom_onehot_feature_names_to_array(
        ["atomic-number-common", "chirality"]
    )
    fl = cuik_molmaker.atom_float_feature_names_to_array([])
    bp = cuik_molmaker.bond_feature_names_to_array(["is-null"])
    mode = cuik_molmaker.reaction_mode_to_int("REAC_DIFF")

    with_stereo = cuik_molmaker.batch_reaction_featurizer(
        reac, prod, oh, fl, bp, False, False, False, mode, ignore_stereo=False
    )
    without_stereo = cuik_molmaker.batch_reaction_featurizer(
        reac, prod, oh, fl, bp, False, False, False, mode, ignore_stereo=True
    )

    assert not np.allclose(with_stereo[0], without_stereo[0])
