# cuik-molmaker
`cuik-molmaker` is a specialized package designed for molecular featurization, converting chemical structures into formats that can be effectively used as inputs for deep learning models, particularly graph neural networks (GNNs).
## Setup conda environment
```
# Set environment variables
export PYTHON_VERSION=3.11
export RDKIT_VERSION=2025.03.2
export TORCH_VERSION=2.6.0

# Create conda env
conda create -n cuik_molmaker_build python=${PYTHON_VERSION} conda-forge::rdkit==${RDKIT_VERSION} conda-forge::pybind11==2.13.6 conda-forge::pytorch-cpu==${TORCH_VERSION} conda-forge::libboost-devel==1.86.0 conda-forge::libboost-python-devel==1.86.0 

# Activate conda env
conda activate cuik_molmaker_build
```

## Build C++ extension and pip wheels
```
# Build C++ extension
TORCH_VERSION=${TORCH_VERSION} RDKIT_VERSION=${RDKIT_VERSION}   PYTHON_VERSION=${PYTHON_VERSION}  python setup.py build_ext --inplace


# Build pip wheels
pip install build
TORCH_VERSION=${TORCH_VERSION} RDKIT_VERSION=${RDKIT_VERSION}   PYTHON_VERSION=${PYTHON_VERSION}  python -m build  --outdir conda_rdkit  --wheel

```

#### HACK: This is to accommodate NumPy<=2.0. In NumPy 2.0, `numpy/core/include` was moved to `numpy/_core/include`
ln -s $CONDA_PREFIX/lib/python3.11/site-packages/numpy/_core/include $CONDA_PREFIX/lib/python3.11/site-packages/numpy/core/include


## Install from wheel and test
```
# Install from wheel
pip install /path/to/wheel_dir/cuik_molmaker-0.1-*py311*-none-manylinux2014_x86_64.whl

# Test that installation is working correctly
pip install pytest
pytest -s tests/python/test_featurize_dims.py
```

## Usage
```
python
>>> import torch
>>> import cuik_molmaker
>>> atom_props_onehot = ["atomic-number", "total-degree", "formal-charge", "chirality", "num-hydrogens", "hybridization"]
>>> atom_property_list_onehot = cuik_molmaker.atom_onehot_feature_names_to_tensor(atom_props_onehot)
>>> print(f"{atom_property_list_onehot}")
tensor([ 0,  2,  9,  6, 10,  5])
```

## Minimal conda env for import and running
```
# Set environment variables
export PYTHON_VERSION=3.11
export RDKIT_VERSION=2025.03.2
export TORCH_VERSION=2.6.0

# Create minimal conda env
conda create -n cuik_molmaker_import python=${PYTHON_VERSION} conda-forge::rdkit==${RDKIT_VERSION} conda-forge::pytorch==${TORCH_VERSION}

conda activate cuik_molmaker_import

# Install from wheel
pip install /path/to/wheel_dir/cuik_molmaker-0.1-*py311*-none-manylinux2014_x86_64.whl

# Test that installation is working correctly
python -c "import cuik_molmaker; print(dir(cuik_molmaker))"
```

## Testing
### Running C++ tests using Catch2
```
# Step 1: Build with test flag set to on
cmake -DCUIKMOLMAKER_BUILD_TESTS=ON -DCMAKE_PREFIX_PATH="$CONDA_PREFIX/lib/python3.11/site-packages/torch/share/cmake;$CONDA_PREFIX" ..

# Optional: List tests/tags
./catch2_tests --list-tests
./catch2_tests --list-tags

# Step 2: Run all C++ tests
cd /path/to/build
./catch2_tests
```

### Running python tests using pytest
```
# Step 1: Install cuik-molmaker using pip
pip install /path/to/wheel_dir/cuik_molmaker-0.1-*py311*-none-manylinux2014_x86_64.whl

# Step 2: Run pytest
cd path/to/repo
pytest -s tests/python
```