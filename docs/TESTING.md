# Testing
The cuik-molmaker package has both C++ and Python tests. C++ code is tested using Catch2, and Python code is tested using pytest.

## Clone repository
```bash
git clone https://github.com/NVIDIA-Digital-Bio/cuik-molmaker
cd cuik-molmaker
```

## Running C++ tests using Catch2

### Build cuik-molmaker with testing enabled
#### Configure cmake with test flag set to on
```bash
mkdir build && cd build
cmake -DCUIKMOLMAKER_BUILD_TESTS=ON -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" ..
```
#### Build and install `cuik-molmaker`
```bash
make -j4 
make install 
mkdir lib && cp libcuik_molmaker_core.so lib/
```

### Run tests with Catch2
```bash
# List tests/tags
LD_LIBRARY_PATH=./lib ./catch2_tests --list-tests
LD_LIBRARY_PATH=./lib ./catch2_tests --list-tags

# Run all C++ tests
LD_LIBRARY_PATH=./lib ./catch2_tests
```


## Running python tests using pytest
#### Install `cuik-molmaker` from pre-built wheel.
Ensure prerequisites are installed using conda. See [README.md](../README.md) for more details.
```bash
python scripts/check_and_install_cuik_molmaker.py
```
Alternatively, you can build and install `cuik-molmaker` from source. See [BUILD.md](BUILD.md) for more details.

#### Run pytest
```bash
pytest tests/python
```