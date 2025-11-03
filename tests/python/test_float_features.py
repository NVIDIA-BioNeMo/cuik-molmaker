# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved. # noqa: E501
# SPDX-License-Identifier: Apache-2.0

import numpy as np

import cuik_molmaker


def test_hydrogen_bond_donor_feature():
    smiles = "C(C(=O)O)C(=O)O"
    atom_onehot_features_names = []
    atom_float_features_names = [
        "hydrogen-bond-donor",
        "hydrogen-bond-acceptor",
        "acidic",
        "basic",
    ]
    atom_onehot_features_array = cuik_molmaker.atom_onehot_feature_names_to_array(
        atom_onehot_features_names
    )
    atom_float_features_array = cuik_molmaker.atom_float_feature_names_to_array(
        atom_float_features_names
    )
    bond_property_array = cuik_molmaker.bond_feature_names_to_array(
        ["is-null", "bond-type-onehot", "conjugated", "in-ring", "stereo"]
    )

    explicit_H, offset_carbon, duplicate_edges, add_self_loop = (
        False,
        False,
        True,
        False,
    )
    atom_feats, bond_feats, edge_index, rev_edge_index, _ = (
        cuik_molmaker.mol_featurizer(
            smiles,
            atom_onehot_features_array,
            atom_float_features_array,
            bond_property_array,
            explicit_H,
            offset_carbon,
            duplicate_edges,
            add_self_loop,
        )
    )

    atom_feats_ref = np.array(
        [
            [
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            [
                0.0,
                0.0,
                1.0,
                0.0,
            ],
            [
                0.0,
                1.0,
                0.0,
                0.0,
            ],
            [
                1.0,
                0.0,
                0.0,
                0.0,
            ],
            [
                0.0,
                0.0,
                1.0,
                0.0,
            ],
            [
                0.0,
                1.0,
                0.0,
                0.0,
            ],
            [
                1.0,
                0.0,
                0.0,
                0.0,
            ],
        ]
    )
    np.testing.assert_allclose(
        atom_feats,
        atom_feats_ref,
        err_msg=f"atom feats diff: {np.abs(atom_feats - atom_feats_ref).sum()}",
    )


def test_ring_size_feature():
    smiles = "Nc1nc(N2CCC(n3cc(C(=O)O)cn3)CC2)nc2ccccc12"

    atom_onehot_features_names = ["ring-size"]
    atom_float_features_names = []
    atom_onehot_features_array = cuik_molmaker.atom_onehot_feature_names_to_array(
        atom_onehot_features_names
    )
    atom_float_features_array = cuik_molmaker.atom_float_feature_names_to_array(
        atom_float_features_names
    )
    bond_property_array = cuik_molmaker.bond_feature_names_to_array(
        ["is-null", "bond-type-onehot", "conjugated", "in-ring", "stereo"]
    )

    explicit_H, offset_carbon, duplicate_edges, add_self_loop = (
        False,
        False,
        True,
        False,
    )
    atom_feats, bond_feats, edge_index, rev_edge_index, _ = (
        cuik_molmaker.mol_featurizer(
            smiles,
            atom_onehot_features_array,
            atom_float_features_array,
            bond_property_array,
            explicit_H,
            offset_carbon,
            duplicate_edges,
            add_self_loop,
        )
    )

    atom_feats_ref = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        ]
    )
    np.testing.assert_allclose(
        atom_feats,
        atom_feats_ref,
        err_msg=f"atom feats diff: {np.abs(atom_feats - atom_feats_ref).sum()}",
    )
