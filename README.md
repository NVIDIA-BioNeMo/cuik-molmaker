# cuik-molmaker
`cuik-molmaker` is a specialized package designed for molecular featurization, converting chemical structures into formats that can be effectively used as inputs for deep learning models, particularly graph neural networks (GNNs).

`cuik-molmaker` is built as a hybrid package, leveraging both C++ and Python to deliver high performance and ease of use. The core featurization logic is implemented in C++ for maximum speed and efficiency, while the Python interface provides a user-friendly API that integrates seamlessly with modern GNN training and inference workflows. This design combines the computational power of C++ with the flexibility and accessibility of Python, making `cuik-molmaker` both fast and intuitive for researchers and developers. As `cuik-molmaker` interfaces with the C++ API of `rdkit`, the produced features are expected to be identical to those produced by `rdkit`.


## Quick start
#### Setup conda environment
```bash
# Set environment variables
export PYTHON_VERSION=3.11
export RDKIT_VERSION=2025.03.2

conda create -n cuik_molmaker_env python=${PYTHON_VERSION} conda-forge::rdkit==${RDKIT_VERSION} conda-forge::pybind11==2.13.6 conda-forge::libboost-devel==1.86.0 conda-forge::libboost-python-devel==1.86.0 

conda activate cuik_molmaker_env
```
This step is optional if you already have a conda environment with the required dependencies.

#### Install wheel from [NVIDIA PyPI](https://pypi.nvidia.com)
We provide a handy script to install the wheel from NVIDIA PyPI based on your OS and other dependencies.
```bash
python scripts/check_and_install_cuik_molmaker.py
```

#### Usage: Computing atom and bond features
```python
import cuik_molmaker
import numpy as np

# List all available atom onehot features
print(cuik_molmaker.list_all_atom_onehot_features())

# Compute atom (atomic number, number of hydrogen, chirality) and bond (bond type) features for acetic acid
acetic_acid_smiles = "CC(=O)O"

# Get atom onehot feature names as NumPy array
atom_onehot_feature_array = cuik_molmaker.atom_onehot_feature_names_to_array(['atomic-number', 'num-hydrogens', 'chirality'])

# Get bond feature names as NumPy array
bond_feature_array = cuik_molmaker.bond_feature_names_to_array(['bond-type-onehot'])

# Set parameters for featurization
explicit_h, offset_carbon, duplicate_edges, add_self_loop = False, False, True, False

# Featurize
all_features =cuik_molmaker.mol_featurizer(acetic_acid_smiles, atom_onehot_feature_array, np.array([]), bond_feature_array, explicit_h, offset_carbon, duplicate_edges, add_self_loop)

# This returns a list of NumPy arrays.
# First index contains atom features
print(all_features[0].shape)

# Second index contains bond features
print(all_features[1].shape)

# Third index contains edge indices in COO format
print(all_features[2].shape)
```

#### Usage: Computing molecular descriptors
```python
from cuik_molmaker.mol_features import MoleculeFeaturizer

featurizer = MoleculeFeaturizer(molecular_descriptor_type="rdkit2D", rdkit2D_normalization_type="fast")

smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O", # aspirin
               "CN(C)CCOC(C1=CC=CC=C1)C1=CC=CC=C1", # diphenhydramine
]
rdkit2D_descriptors = featurizer.featurize(smiles_list)

# Print the shape of the descriptors
print(rdkit2D_descriptors.shape)

```

## Source of acceleration
The hybrid C++/Python design of `cuik-molmaker` allows for the core featurization logic to be implemented in C++ and reduces the python overhead. Another source of acceleration is the creation of features for the entire minibatch of SMILES at once, which saves the overhead of creating memory allocation and concatenation. 

## Additional Documentation

| File | Description |
|------|-------------|
| [USAGE.md](docs/USAGE.md) | Examples and instructions for using `cuik-molmaker` to featurize molecules including batching. |
| [FEATURES.md](docs/FEATURES.md) | Detailed list and explanation of all atom and bond features available for featurization. |
| [BUILD.md](docs/BUILD.md) | Step-by-step instructions for building `cuik-molmaker` from source, including prerequisites and troubleshooting. |
| [TESTING.md](docs/TESTING.md) | Guidelines and commands for running the test suite to verify installation and functionality. |

## Hardware Requirements
`cuik-molmaker` is designed to run on any CPU-based system.

## Adoption
`cuik-molmaker` has currently been integrated into the following projects:
- [Chemprop](https://github.com/chemprop/chemprop): `cuik-molmaker` is available for use with conda and Docker installations of Chemprop. It can be enabled by setting `--use-cuikmolmaker-featurization` flag in the command line with all use cases: training, prediction, fingerprinting, and hyperparameter optimization.

