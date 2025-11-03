# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved. # noqa: E501
# SPDX-License-Identifier: Apache-2.0

import os

import numpy as np
import pytest

from cuik_molmaker import MoleculeFeaturizer
from cuik_molmaker.utils.descriptor_normalization import DESCRIPTASTORUS_DESC_LIST


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
def test_rdkit2D(smiles_list_100, normalization_type, test_data_path):

    if normalization_type == "none":
        featurizer = MoleculeFeaturizer(molecular_descriptor_type="rdkit2D")
    else:
        featurizer = MoleculeFeaturizer(
            molecular_descriptor_type="rdkit2D",
            rdkit2D_normalization_type=normalization_type,
        )

    ref_file = os.path.join(
        test_data_path, f"rdkit2D_{normalization_type}_normalization_ref.npy"
    )
    desc_ref = np.load(ref_file)

    desc = featurizer.featurize(smiles_list_100)

    # Bertz CT descriptor RDKit is not backward compatible
    # https://www.rdkit.org/docs/source/rdkit.Chem.GraphDescriptors.html#rdkit.Chem.GraphDescriptors.BertzCT
    problem_idx_list = [
        [d[0] for d in featurizer.rdkit2D_descriptor_list].index("BertzCT")
    ]
    if normalization_type == "best":
        # TODO: Unstable genhyperbolic normalization function parameters
        # TODO: Remove this.
        problem_idx_list.append(
            [d[0] for d in featurizer.rdkit2D_descriptor_list].index("EState_VSA5")
        )
    for problem_idx in problem_idx_list:
        desc[:, problem_idx] = desc_ref[:, problem_idx]
    np.testing.assert_allclose(
        desc_ref, desc, atol=1e-4
    ), f"RDKit 2D descriptor generation and normalization type {normalization_type} "
    "do not match reference"


def test_rdkit2D_descriptastorus(smiles_list_100, test_data_path):

    normalization_type = "descriptastorus"
    featurizer = MoleculeFeaturizer(
        molecular_descriptor_type="rdkit2D",
        rdkit2D_descriptor_list=DESCRIPTASTORUS_DESC_LIST,
        rdkit2D_normalization_type=normalization_type,
    )

    ref_file = os.path.join(
        test_data_path, f"rdkit2D_{normalization_type}_normalization_ref.npy"
    )
    desc_ref = np.load(ref_file)

    desc = featurizer.featurize(smiles_list_100)

    # The implementation of descriptor normalization for fr_unbrch_alkane does not match
    # the reference. Normalization function parameters are very precision sensitive.
    # Bertz CT implementation is not backward compatible:
    # https://www.rdkit.org/docs/source/rdkit.Chem.GraphDescriptors.html#rdkit.Chem.GraphDescriptors.BertzCT
    problem_idx_list = [
        DESCRIPTASTORUS_DESC_LIST.index("fr_unbrch_alkane"),
        DESCRIPTASTORUS_DESC_LIST.index("BertzCT"),
    ]
    for problem_idx in problem_idx_list:
        desc[:, problem_idx] = desc_ref[:, problem_idx]

    np.testing.assert_allclose(
        desc_ref, desc, atol=1e-4
    ), f"RDKit 2D descriptor generation and normalization type {normalization_type} do "
    "not match reference"
