# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0

from setuptools import setup
import os
import shutil
import sys
import configparser
from setuptools.command.build_ext import build_ext
import subprocess
import sysconfig

class CMakeBuild(build_ext):
    def run(self):
        # Ensure build directory exists
        build_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "build")
        os.makedirs(build_dir, exist_ok=True)

        # Prepare CMake command
        cmake_args = [
            "cmake",
            f"-DCMAKE_PREFIX_PATH={os.environ['CONDA_PREFIX']}/lib/python{PYTHON_VERSION}/site-packages/torch/share/cmake;{os.environ['CONDA_PREFIX']}",
            os.path.abspath(os.path.dirname(__file__))
        ]

        # Run CMake
        print("Running CMake:", " ".join(cmake_args))
        subprocess.check_call(cmake_args, cwd=build_dir)

        # Run Make
        print("Running make -j4")
        subprocess.check_call(["make", "-j4"], cwd=build_dir)

        # Call the original build_ext to copy .so files, etc.
        super().run()

# Get torch and rdkit versions with priority: config settings > env vars > default
TORCH_VERSION = os.environ.get('TORCH_VERSION')
RDKIT_VERSION = os.environ.get('RDKIT_VERSION')
PYTHON_VERSION = os.environ.get('PYTHON_VERSION')


if TORCH_VERSION is None:
    print("Error: PyTorch version is not set.")
    print("Please specify it as follows:")
    print("Using environment variables: TORCH_VERSION=2.6.0 RDKIT_VERSION=2024.03.4 PYTHON_VERSION=3.11 pip install .")
    sys.exit(1)

if RDKIT_VERSION is None:
    print("Error: RDKit version is not set.")
    print("Please specify it as follows:")
    print("Using environment variables: TORCH_VERSION=2.6.0 RDKIT_VERSION=2024.03.4 PYTHON_VERSION=3.11 pip install .")
    sys.exit(1)

if PYTHON_VERSION is None:
    print("Error: Python version is not set.")
    print("Please specify it as follows:")
    print("Using environment variables: TORCH_VERSION=2.6.0 RDKIT_VERSION=2024.03.4 PYTHON_VERSION=3.11 pip install .")
    sys.exit(1)


# Update setup.cfg with the Python tag
PYTHON_DIGIT_ONLY_VERSION = PYTHON_VERSION.replace('.', '')

config = configparser.ConfigParser()
config.read('setup.cfg')
if 'bdist_wheel' not in config:
    config['bdist_wheel'] = {}
config['bdist_wheel']['python-tag'] = f"py{PYTHON_DIGIT_ONLY_VERSION}"
config['bdist_wheel']['plat_name'] = sysconfig.get_platform()
with open('setup.cfg', 'w') as f:
    config.write(f)

print(f"Environment TORCH_VERSION: {os.environ.get('TORCH_VERSION')}")
print(f"Environment RDKIT_VERSION: {os.environ.get('RDKIT_VERSION')}")

print(f"Building with TORCH_VERSION={TORCH_VERSION}, RDKIT_VERSION={RDKIT_VERSION}, PYTHON_VERSION={PYTHON_VERSION}")


# Create package directory structure first
dest_dir = os.path.join('cuik_molmaker')
lib_dir = os.path.join(dest_dir, 'lib')
os.makedirs(dest_dir, exist_ok=True)
os.makedirs(lib_dir, exist_ok=True)


# Check if .so file exists, display helpful message if not
so_file = os.path.join('build', f'cuik_molmaker.cpython-{PYTHON_DIGIT_ONLY_VERSION}-x86_64-linux-gnu.so')
print(f"Looking for compiled extension at: {so_file}")

if os.path.exists(so_file):
    print(f"Found compiled extension, copying to {dest_dir}")
    shutil.copy2(so_file, dest_dir)
else:
    print("WARNING: Compiled extension not found. You need to build the C++ extension first.")
    print("Try running: python setup.py build_ext --inplace")
    # Uncomment to abort if .so is missing:
    # sys.exit(1)

# Check for the shared library
lib_file = os.path.join('build', 'libcuik_molmaker_core.so')
if os.path.exists(lib_file):
    print(f"Found shared library, copying to {lib_dir}")
    shutil.copy2(lib_file, lib_dir)
else:
    print("WARNING: Shared library not found. You need to build the C++ extension first.")
    print("Try running: python setup.py build_ext --inplace")

# Ensure __init__.py exists
init_file = os.path.join(dest_dir, '__init__.py')
if not os.path.exists(init_file):
    print(f"Creating {init_file}")
    with open(init_file, 'w') as f:
        f.write("# Import compiled extension\n")
        f.write("from pathlib import Path\n")
        f.write("import os\n")
        f.write("import sys\n")
        f.write("\n")
        f.write("# Find the .so file in this directory\n")
        f.write("_module_dir = Path(__file__).parent\n")
        f.write("for file in os.listdir(_module_dir):\n")
        f.write("    if file.endswith('.so') and 'cpython' in file:\n")
        f.write("        # Add the extension module directly\n")
        f.write("        from importlib.machinery import ExtensionFileLoader\n")
        f.write("        from importlib.util import spec_from_loader, module_from_spec\n")
        f.write("        \n")
        f.write("        _loader = ExtensionFileLoader('cuik_molmaker', str(_module_dir / file))\n") 
        f.write("        _spec = spec_from_loader('cuik_molmaker', _loader)\n")
        f.write("        _module = module_from_spec(_spec)\n")
        f.write("        _loader.exec_module(_module)\n")
        f.write("        \n")
        f.write("        # Import all attributes from the module\n")
        f.write("        for attr in dir(_module):\n")
        f.write("            if not attr.startswith('_'):\n") 
        f.write("                globals()[attr] = getattr(_module, attr)\n")
        f.write("        break\n")


setup(
    name="cuik_molmaker",
    version="0.1",
    author="S. Veccham",
    author_email="sveccham@nvidia.com",
    description="C++ module for featurizing molecules",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    license="Apache 2.0",
    # Explicitly list the package instead of using find_packages()
    packages=["cuik_molmaker", "cuik_molmaker.lib"],
    package_data={
        'cuik_molmaker': ['*.so'],  # Include Python extension
        'cuik_molmaker.lib': ['*.so'],  # Include shared libraries
    },
    cmdclass={
        'build_ext': CMakeBuild,
    },
    install_requires=[
        f'rdkit=={RDKIT_VERSION}',
        f'torch=={TORCH_VERSION}',
    ],
    build_requires=[
        f'rdkit=={RDKIT_VERSION}',
        f'torch=={TORCH_VERSION}',
    ],
    tests_require=['pytest'],
    python_requires=f'=={PYTHON_VERSION}.*',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        f'Programming Language :: Python :: {PYTHON_VERSION}',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
) 