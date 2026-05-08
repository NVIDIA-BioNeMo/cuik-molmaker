# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved. # noqa: E501
# SPDX-License-Identifier: Apache-2.0

import argparse
import json
import time
from importlib.resources import files

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

from cuik_molmaker.utils.descriptor_normalization import (
    DESCRIPTASTORUS_DESC_LIST,
    get_normalization_functions,
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Path to the file containing SMILES strings",
    )
    parser.add_argument(
        "--output_file", type=str, required=True, help="Path to the output file"
    )
    parser.add_argument(
        "--normalization_params",
        type=str,
        choices=("best", "fast", "descriptastorus"),
        help="Path to the JSON file containing normalization parameters",
    )
    parser.add_argument(
        "--smiles_column",
        type=str,
        required=True,
        help="Name of the column containing SMILES strings",
    )
    parser.add_argument(
        "--rdkit2D_descriptors_nan_sub_value",
        type=float,
        help="Value to substitute for NaN values in rdkit2D descriptors",
    )
    return parser.parse_args()


class MoleculeFeaturizer:
    """A class for computing molecular descriptors from SMILES strings.

    Currently, this class provides functionality for computing RDKit 2D molecular
    descriptors.
    """

    def __init__(self, molecular_descriptor_type: str, **kwargs):

        if molecular_descriptor_type == "rdkit2D":
            self.molecular_descriptor_type = "rdkit2D"

            # Set up rdkit2D descriptors params
            all_desc_fn_dict = {desc[0]: desc[1] for desc in Descriptors._descList}
            if "rdkit2D_descriptor_list" in kwargs:
                self.rdkit2D_descriptor_list = [
                    (desc, all_desc_fn_dict[desc])
                    for desc in kwargs["rdkit2D_descriptor_list"]
                ]
            else:
                # Use all rdkit2D descriptors
                self.rdkit2D_descriptor_list = Descriptors._descList

            if "rdkit2D_normalization_type" in kwargs:
                self.rdkit2D_normalization_type = kwargs["rdkit2D_normalization_type"]
                assert self.rdkit2D_normalization_type in (
                    "best",
                    "fast",
                    "descriptastorus",
                    None,
                )
            else:
                self.rdkit2D_normalization_type = None

            if self.rdkit2D_normalization_type is not None:
                resource = files("cuik_molmaker").joinpath(
                    "data",
                    f"{self.rdkit2D_normalization_type}_normalization_params.json",
                )
                with resource.open("r", encoding="utf-8") as f:
                    norm_params = json.load(f)
                rdkit2D_desc_name_list = [d[0] for d in self.rdkit2D_descriptor_list]
                self.rdkit2D_normalization_fn_dict = get_normalization_functions(
                    rdkit2D_desc_name_list, norm_params
                )

        else:
            raise ValueError(
                f"Invalid molecular descriptor type: {molecular_descriptor_type}."
                " Only 'rdkit2D' is supported."
            )

        if "rdkit2D_descriptors_nan_sub_value" in kwargs:
            self.rdkit2D_descriptors_nan_sub_value = kwargs[
                "rdkit2D_descriptors_nan_sub_value"
            ]
        else:
            self.rdkit2D_descriptors_nan_sub_value = None

    def _compute_rdkit2D_descriptors(self, smi: str) -> np.ndarray:
        """Computes RDKit 2D molecular descriptors for a single SMILES string.

        This function converts a SMILES string to an RDKit molecule object and computes
        all configured 2D molecular descriptors. The molecular descriptors are returned
        as a numpy array. If the SMILES string is invalid or cannot be parsed, it
        returns an array of NaN values.

        Args:
            smi: SMILES string representing the molecule to compute descriptors for.

        Returns:
            numpy.ndarray: Array of computed descriptor values. If the SMILES is
            invalid, returns an array of NaN values with length equal to the number of
                configured descriptors.
        """

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return [np.nan] * len(self.rdkit2D_descriptor_list)
        return [func[1](mol) for func in self.rdkit2D_descriptor_list]

    def compute_rdkit2D_descriptors(self, smi_list: list[str]) -> np.ndarray:
        """Computes RDKit 2D molecular descriptors for a list of SMILES strings.

        This function processes a list of SMILES strings and computes 2D molecular
        descriptors for each molecule. The function performs the following steps:
        1. Computes descriptors for each SMILES  using _compute_rdkit2D_descriptors
        2. Applies normalization if normalization_type is configured
        3. Substitutes NaN values with a specified value if
        rdkit2D_descriptors_nan_sub_value is set

        Args:
            smi_list: List of SMILES strings representing molecules to compute
            descriptors for.

        Returns:
            numpy.ndarray: 2D array where each row contains the computed descriptors for
                one molecule. Shape is (n_molecules, n_descriptors). If normalization is
                enabled, the descriptors are normalized. NaN values are replaced if
                rdkit2D_descriptors_nan_sub_value is specified.

        """

        desc = [self._compute_rdkit2D_descriptors(smi) for smi in smi_list]
        desc = np.array(desc)
        if self.rdkit2D_normalization_type is not None:
            desc = self._normalize_rdkit2D_descriptors(desc)

        # Set nans to a value
        if self.rdkit2D_descriptors_nan_sub_value is not None:
            desc[np.isnan(desc)] = self.rdkit2D_descriptors_nan_sub_value

        return desc

    def featurize(self, smi_list: list[str]) -> np.ndarray:
        """Public interface for computing molecular descriptors from SMILES strings.

        This is the primary method for converting a list of SMILES strings into
        molecular descriptors. It delegates to the appropriate descriptor computation
        method based on the configured molecular_descriptor_type.

        Args:
            smi_list: List of SMILES strings representing molecules to featurize.

        Returns:
            numpy.ndarray: 2D array of molecular descriptors where each row corresponds
                to one molecule and each column corresponds to one descriptor.
                Shape is (n_molecules, n_descriptors).
        """
        return self.compute_rdkit2D_descriptors(smi_list)

    def _normalize_rdkit2D_descriptors(self, desc: np.ndarray) -> np.ndarray:

        norm_desc = np.zeros_like(desc, dtype=np.float64)
        for i, desc_name in enumerate(self.rdkit2D_descriptor_list):
            norm_fn = self.rdkit2D_normalization_fn_dict[desc_name[0]]
            norm_desc[:, i] = norm_fn(desc[:, i])

        return norm_desc


def main():
    """Main entry point for the molecular features computation script."""
    args = parse_args()
    smi_list = pd.read_csv(args.input_file)[args.smiles_column].to_list()

    desc_name_fn_list = Descriptors._descList
    if args.normalization_params == "descriptastorus":
        print(
            f"Using {len(DESCRIPTASTORUS_DESC_LIST)} descriptors from Descriptastorus"
            " package as descriptastorus normalization parameters is specified."
        )
        desc_name_list = DESCRIPTASTORUS_DESC_LIST
    else:
        desc_name_list = [d[0] for d in desc_name_fn_list]
    desc_start_time = time.time()
    featurizer = MoleculeFeaturizer(
        molecular_descriptor_type="rdkit2D",
        rdkit2D_descriptor_list=desc_name_list,
        rdkit2D_normalization_type=args.normalization_params,
    )
    desc = featurizer.featurize(smi_list)
    desc_end_time = time.time()
    print(f"Descriptor computation time: {(desc_end_time - desc_start_time):.4f}s")

    output_df = pd.DataFrame(desc, columns=desc_name_list)
    output_df.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main()
