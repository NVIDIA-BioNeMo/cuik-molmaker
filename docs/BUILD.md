# Build cuik_molmaker
Follow the instructions below to build cuik_molmaker from source.

## Setup linux environment
#### Install eigen
```bash
sudo apt-get install libeigen3-dev
```
#### Verify eigen is installed correctly
```bash
dpkg -s libeigen3-dev
```
#### Install cmake
```bash
sudo apt install cmake
# Check cmake version
cmake --version
```

## Setup conda environment
```bash
# Set desired python and RDKit versions as environment variables. 
# Here is an example for python 3.11, RDKit 2025.03.2.
export PYTHON_VERSION=3.11
export RDKIT_VERSION=2025.03.2

# Create conda env
conda create -n cuik_molmaker_build python=$PYTHON_VERSION conda-forge::rdkit==$RDKIT_VERSION conda-forge::pybind11==2.13.6 conda-forge::librdkit-dev==$RDKIT_VERSION conda-forge::pybind11==2.13.6  conda-forge::libboost-devel==1.86.0 conda-forge::libboost-python-devel==1.86.0

# Activate conda env
conda activate cuik_molmaker_build
```
> **Note:** If you are building against an older version of RDKit, you may need to use an older version of Boost that is compatible with that RDKit release. Check the RDKit documentation for the recommended Boost version for your RDKit version.

## Build cuik_molmaker
```bash
# Clone cuik_molmaker repo
git clone https://github.com/NVIDIA-Digital-Bio/cuik-molmaker
cd cuik_molmaker

# Build C++ extension for cuik_molmaker
RDKIT_VERSION=$RDKIT_VERSION PYTHON_VERSION=$PYTHON_VERSION python setup.py build_ext --inplace

# Install build
pip install build

# Build wheel
RDKIT_VERSION=$RDKIT_VERSION PYTHON_VERSION=$PYTHON_VERSION python -m build --outdir path/to/output/directory --wheel
```

## Install cuik_molmaker from wheel
```bash
pip install path/to/output/directory/cuik_molmaker*.whl
```

## Test installation
```bash
pip install pytest
pytest tests/python/
```
If the installation was successful, you should see that all tests pass. If not, please refer to the Troubleshooting section below.

## Troubleshooting
- If you see a CMake error about `RDKit::GraphMol` including the non-existent path `$CONDA_PREFIX/include/python3.9` even though you are using a different python version, symlink it to the correct path as shown:
```bash
ln -s $CONDA_PREFIX/lib/python$PYTHON_VERSION $CONDA_PREFIX/lib/python3.9
ln -s $CONDA_PREFIX/include/python$PYTHON_VERSION $CONDA_PREFIX/include/python3.9
```


- If any of the build steps fails, delete the `build`, `cuik-molmaker.egg-info`, `cuik-molmaker` and wheel output directories before trying again.