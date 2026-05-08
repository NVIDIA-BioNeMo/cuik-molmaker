# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved. # noqa: E501
# SPDX-License-Identifier: Apache-2.0

import json
from importlib.resources import files

import numpy as np
import pytest
from rdkit import Chem

from cuik_molmaker import MoleculeFeaturizer
from cuik_molmaker.utils.descriptor_normalization import (
    DESCRIPTASTORUS_DESC_LIST,
    get_normalization_functions,
)


def test_rdkit2D_some_props(smiles_list_100, test_data_path):
    num_props = 50
    featurizer = MoleculeFeaturizer(
        molecular_descriptor_type="rdkit2D",
        rdkit2D_descriptor_list=DESCRIPTASTORUS_DESC_LIST[:num_props],
    )
    desc = featurizer.featurize(smiles_list_100)
    assert (
        desc.shape[1] == num_props
    ), "RDKit 2D descriptor cannot generate subset of properties"


@pytest.mark.parametrize("normalization_type", ["none", "fast", "best"])
def test_rdkit2D(smiles_list_100, normalization_type):

    if normalization_type == "none":
        featurizer = MoleculeFeaturizer(molecular_descriptor_type="rdkit2D")
    else:
        featurizer = MoleculeFeaturizer(
            molecular_descriptor_type="rdkit2D",
            rdkit2D_normalization_type=normalization_type,
        )

    desc_list = featurizer.rdkit2D_descriptor_list
    rows = []
    for smi in smiles_list_100:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rows.append([np.nan] * len(desc_list))
        else:
            rows.append([func(mol) for _, func in desc_list])
    desc_ref = np.array(rows)

    if normalization_type != "none":
        # Load the same normalization params and apply the same per-column
        # transform as MoleculeFeaturizer._normalize_rdkit2D_descriptors.
        norm_params_path = files("cuik_molmaker").joinpath(
            "data", f"{normalization_type}_normalization_params.json"
        )
        with norm_params_path.open("r", encoding="utf-8") as f:
            norm_params = json.load(f)
        desc_names = [d[0] for d in desc_list]
        norm_fn_dict = get_normalization_functions(desc_names, norm_params)
        norm_desc = np.zeros_like(desc_ref, dtype=np.float64)
        for i, name in enumerate(desc_names):
            norm_desc[:, i] = norm_fn_dict[name](desc_ref[:, i])
        desc_ref = norm_desc

    desc = featurizer.featurize(smiles_list_100)

    np.testing.assert_allclose(
        desc_ref, desc, atol=1e-4
    ), f"RDKit 2D descriptor generation and normalization type {normalization_type} "
    "do not match reference"


def test_rdkit2D_descriptastorus(smiles_list_100):

    normalization_type = "descriptastorus"
    featurizer = MoleculeFeaturizer(
        molecular_descriptor_type="rdkit2D",
        rdkit2D_descriptor_list=DESCRIPTASTORUS_DESC_LIST,
        rdkit2D_normalization_type=normalization_type,
    )

    # Build the reference live from RDKit by mirroring the steps in
    # MoleculeFeaturizer.compute_rdkit2D_descriptors (src/mol_features.py).
    desc_list = featurizer.rdkit2D_descriptor_list
    rows = []
    for smi in smiles_list_100:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rows.append([np.nan] * len(desc_list))
        else:
            rows.append([func(mol) for _, func in desc_list])
    desc_ref = np.array(rows)

    norm_params_path = files("cuik_molmaker").joinpath(
        "data", f"{normalization_type}_normalization_params.json"
    )
    with norm_params_path.open("r", encoding="utf-8") as f:
        norm_params = json.load(f)
    desc_names = [d[0] for d in desc_list]
    norm_fn_dict = get_normalization_functions(desc_names, norm_params)
    norm_desc = np.zeros_like(desc_ref, dtype=np.float64)
    for i, name in enumerate(desc_names):
        norm_desc[:, i] = norm_fn_dict[name](desc_ref[:, i])
    desc_ref = norm_desc

    desc = featurizer.featurize(smiles_list_100)

    np.testing.assert_allclose(
        desc_ref, desc, atol=1e-4
    ), f"RDKit 2D descriptor generation and normalization type {normalization_type} do "
    "not match reference"
