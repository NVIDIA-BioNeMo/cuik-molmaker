# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved. # noqa: E501
# SPDX-License-Identifier: Apache-2.0

from rdkit import Chem

import cuik_molmaker


def test_import_cuik_molmaker():
    print(f"cuik_molmaker.__file__: {cuik_molmaker.__file__}")
    assert cuik_molmaker is not None


def test_mol_featurizer_dims(sample_smiles1):
    atom_onehot_property_list = cuik_molmaker.atom_onehot_feature_names_to_array(
        [
            "atomic-number",
            "total-degree",
            "formal-charge",
            "chirality",
            "num-hydrogens",
            "hybridization",
        ]
    )
    atom_float_property_list = cuik_molmaker.atom_float_feature_names_to_array(
        ["aromatic", "mass"]
    )
    bond_property_list = cuik_molmaker.bond_feature_names_to_array(
        ["is-null", "bond-type-onehot", "conjugated", "in-ring", "stereo"]
    )

    explicit_H, offset_carbon, duplicate_edges, add_self_loop = (
        False,
        False,
        True,
        False,
    )

    all_feats = cuik_molmaker.mol_featurizer(
        sample_smiles1,
        atom_onehot_property_list,
        atom_float_property_list,
        bond_property_list,
        explicit_H,
        offset_carbon,
        duplicate_edges,
        add_self_loop,
    )
    atom_feats, bond_feats, edge_index, rev_edge_index, batch = all_feats

    mol = Chem.MolFromSmiles(sample_smiles1)
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()

    assert (
        atom_feats.shape[0] == num_atoms
    ), f"atom_feats.shape[0]: {atom_feats.shape[0]}, num_atoms: {num_atoms}"
    assert (
        bond_feats.shape[0] == num_bonds * 2
    ), f"bond_feats.shape[0]: {bond_feats.shape[0]}, num_bonds: {num_bonds}"
    assert (
        edge_index.shape[1] == num_bonds * 2
    ), f"edge_index.shape[1]: {edge_index.shape[1]}, num_bonds: {num_bonds}"
    assert (
        rev_edge_index.shape[0] == num_bonds * 2
    ), f"rev_edge_index.shape[0]: {rev_edge_index.shape[0]}, num_bonds: {num_bonds}"
    assert (
        batch.shape[0] == num_atoms
    ), f"batch.shape[0]: {batch.shape[0]}, num_atoms: {num_atoms}"


def test_mol_featurizer_variations_dims(sample_smiles1):
    atom_onehot_property_list = cuik_molmaker.atom_onehot_feature_names_to_array(
        [
            "atomic-number",
            "total-degree",
            "formal-charge",
            "chirality",
            "num-hydrogens",
            "hybridization",
        ]
    )
    atom_float_property_list = cuik_molmaker.atom_float_feature_names_to_array(
        ["aromatic", "mass"]
    )
    bond_property_list = cuik_molmaker.bond_feature_names_to_array(
        ["is-null", "bond-type-onehot", "conjugated", "in-ring", "stereo"]
    )
    explicit_H, offset_carbon, duplicate_edges, add_self_loop = (
        False,
        False,
        False,
        False,
    )

    all_feats = cuik_molmaker.mol_featurizer(
        sample_smiles1,
        atom_onehot_property_list,
        atom_float_property_list,
        bond_property_list,
        explicit_H,
        offset_carbon,
        duplicate_edges,
        add_self_loop,
    )
    atom_feats, bond_feats, edge_index, rev_edge_index, batch = all_feats

    mol = Chem.MolFromSmiles(sample_smiles1)
    num_atoms_ref = mol.GetNumAtoms()
    num_bonds_ref = mol.GetNumBonds()

    assert (
        atom_feats.shape[0] == num_atoms_ref
    ), f"atom_feats.shape[0]: {atom_feats.shape[0]}, num_atoms: {num_atoms_ref}"
    assert (
        bond_feats.shape[0] == num_bonds_ref
    ), f"bond_feats.shape[0]: {bond_feats.shape[0]}, num_bonds: {num_bonds_ref}"
    assert (
        edge_index.shape[1] == num_bonds_ref
    ), f"edge_index.shape[1]: {edge_index.shape[1]}, num_bonds: {num_bonds_ref}"
    assert (
        rev_edge_index.shape[0] == num_bonds_ref * 2
    ), f"rev_edge_index.shape[0]: {rev_edge_index.shape[0]}, num_bonds: {num_bonds_ref}"
    assert (
        batch.shape[0] == num_atoms_ref
    ), f"batch.shape[0]: {batch.shape[0]}, num_atoms: {num_atoms_ref}"


def test_batch_mol_featurizer_dims(smiles_list_100):

    atom_onehot_property_list = cuik_molmaker.atom_onehot_feature_names_to_array(
        [
            "atomic-number",
            "total-degree",
            "formal-charge",
            "chirality",
            "num-hydrogens",
            "hybridization",
        ]
    )
    atom_float_property_list = cuik_molmaker.atom_float_feature_names_to_array(
        ["aromatic", "mass"]
    )
    bond_property_list = cuik_molmaker.bond_feature_names_to_array(
        ["is-null", "bond-type-onehot", "conjugated", "in-ring", "stereo"]
    )
    explicit_H, offset_carbon, duplicate_edges, add_self_loop = (
        False,
        False,
        True,
        False,
    )

    all_feats = cuik_molmaker.batch_mol_featurizer(
        smiles_list_100,
        atom_onehot_property_list,
        atom_float_property_list,
        bond_property_list,
        explicit_H,
        offset_carbon,
        duplicate_edges,
        add_self_loop,
    )
    atom_feats, bond_feats, edge_index, rev_edge_index, batch = all_feats

    # Prepare reference data
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list_100]
    total_num_atoms = sum(mol.GetNumAtoms() for mol in mol_list)
    total_num_bonds = sum(mol.GetNumBonds() for mol in mol_list)

    assert (
        atom_feats.shape[0] == total_num_atoms
    ), f"atom_feats.shape[0]: {atom_feats.shape[0]}, total_num_atoms: {total_num_atoms}"
    assert (
        bond_feats.shape[0] == total_num_bonds * 2
    ), f"bond_feats.shape[0]: {bond_feats.shape[0]}, total_num_bonds: {total_num_bonds}"
    assert (
        edge_index.shape[1] == total_num_bonds * 2
    ), f"edge_index.shape[1]: {edge_index.shape[1]}, total_num_bonds: {total_num_bonds}"
    assert (
        rev_edge_index.shape[0] == total_num_bonds * 2
    ), f"rev_edge_index.shape[0]: {rev_edge_index.shape[0]}, "
    f"total_num_bonds: {total_num_bonds}"
    assert (
        batch.shape[0] == total_num_atoms
    ), f"batch.shape[0]: {batch.shape[0]}, total_num_atoms: {total_num_atoms}"
