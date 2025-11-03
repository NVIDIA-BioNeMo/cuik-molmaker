# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved. # noqa: E501
# SPDX-License-Identifier: Apache-2.0

import configparser
import glob
import os
import platform
import shutil
import subprocess
import sys
import sysconfig

from setuptools import find_packages, setup
from setuptools.command.build_ext import build_ext

# Set global vars
RDKIT_VERSION = os.environ.get("RDKIT_VERSION")
PYTHON_VERSION = os.environ.get("PYTHON_VERSION")
CXX11_ABI = os.environ.get("CUIKMOLMAKER_CXX11_ABI")
SYSTEM = platform.system()


class CMakeBuild(build_ext):
    def run(self):

        # Detect if we're doing an install against pip rdkit
        cmake_extra_args = []
        cuikmolmaker_build_against_pip = os.getenv(
            "CUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT"
        )
        cmake_extra_args.extend(
            [
                f"-DCUIKMOLMAKER_CXX11_ABI={CXX11_ABI}",
            ]
        )
        if cuikmolmaker_build_against_pip:
            CUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR = os.getenv(
                "CUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR"
            )
            if not CUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR:
                raise ValueError(
                    "CUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR must be set when "
                    "building against pip rdkit"
                )
            CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR = os.getenv(
                "CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR"
            )
            if not CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR:
                raise ValueError(
                    "CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR must be set when "
                    "building against pip rdkit"
                )
            CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR = os.getenv(
                "CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR"
            )
            if not CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR:
                raise ValueError(
                    "CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR must be set "
                    "when building against pip rdkit"
                )
            cmake_extra_args.extend(
                [
                    "-DCUIKMOLMAKER_BUILD_AGAINST_PIP_RDKIT=ON",
                    f"-DCUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR="
                    f"{CUIKMOLMAKER_BUILD_AGAINST_PIP_LIBDIR}",
                    f"-DCUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR="
                    f"{CUIKMOLMAKER_BUILD_AGAINST_PIP_INCDIR}",
                    f"-DCUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR="
                    f"{CUIKMOLMAKER_BUILD_AGAINST_PIP_BOOSTINCLUDEDIR}",
                ]
            )
        # Prepare platform-specific CMake command
        platform_name = platform.system()

        if platform_name == "Windows":
            cmake_prefix_path = os.environ["CONDA_PREFIX"]
            cmake_extra_args.extend(
                ["-DCMAKE_BUILD_TYPE=Release", "-G", "Ninja", "-S", ".", "-B", "build"]
            )

            cmake_cmd = [
                "cmake",
                f"-DCMAKE_PREFIX_PATH={cmake_prefix_path}",
            ] + cmake_extra_args

            # Run CMake configure
            print("Running CMake configure command:", " ".join(cmake_cmd))
            subprocess.check_call(
                cmake_cmd, cwd=os.path.abspath(os.path.dirname(__file__))
            )

            # Run CMake build
            cmake_build_cmd = ["cmake", "--build", "build", "-j", "4"]
            print("Running CMake build command:", " ".join(cmake_build_cmd))
            subprocess.check_call(
                cmake_build_cmd, cwd=os.path.abspath(os.path.dirname(__file__))
            )

        elif platform_name in ("Linux", "Darwin"):

            # Ensure build directory exists
            build_dir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)), "build"
            )
            os.makedirs(build_dir, exist_ok=True)

            # Prepare cmake command
            cmake_args = [
                "cmake",
                f"-DCMAKE_PREFIX_PATH={os.environ['CONDA_PREFIX']}",
            ]
            cmake_args.extend(cmake_extra_args)
            cmake_args.append(os.path.abspath(os.path.dirname(__file__)))

            # Run CMake configure
            print("Running CMake:", " ".join(cmake_args))
            subprocess.check_call(cmake_args, cwd=build_dir)

            # Run make
            print("Running make -j4")
            subprocess.check_call(["make", "-j4"], cwd=build_dir)
        else:
            raise ValueError(f"Unsupported platform: {platform_name}")

        # Call the original build_ext to copy .so files, etc.
        super().run()


if RDKIT_VERSION is None:
    print("Error: RDKit version is not set.")
    print("Please specify it as follows using environment variables:")
    print(
        "RDKIT_VERSION=2024.03.4 PYTHON_VERSION=3.11 CUIKMOLMAKER_CXX11_ABI=ON "
        "python setup.py build_ext --inplace"
    )
    sys.exit(1)

if PYTHON_VERSION is None:
    print("Error: Python version is not set.")
    print("Please specify it as follows using environment variables:")
    print(
        "RDKIT_VERSION=2024.03.4 PYTHON_VERSION=3.11 CUIKMOLMAKER_CXX11_ABI=ON "
        "python setup.py build_ext --inplace"
    )
    sys.exit(1)

if CXX11_ABI is None:
    print("Error: CXX11_ABI is not set.")
    print("CXX11_ABI can be either ON or OFF.")
    print("Please specify it as follows using environment variables:")
    print(
        "RDKIT_VERSION=2024.03.4 PYTHON_VERSION=3.11 CUIKMOLMAKER_CXX11_ABI=ON "
        "python setup.py build_ext --inplace"
    )
    sys.exit(1)

# Update setup.cfg with the Python tag
PYTHON_DIGIT_ONLY_VERSION = PYTHON_VERSION.replace(".", "")

config = configparser.ConfigParser()
config.read("setup.cfg")
if "bdist_wheel" not in config:
    config["bdist_wheel"] = {}
config["bdist_wheel"]["python-tag"] = f"py{PYTHON_DIGIT_ONLY_VERSION}"
config["bdist_wheel"]["plat_name"] = sysconfig.get_platform()
with open("setup.cfg", "w") as f:
    config.write(f)

print(
    f"Building with RDKIT_VERSION={RDKIT_VERSION}, "
    f"PYTHON_VERSION={PYTHON_VERSION}, "
    f"CXX11_ABI={CXX11_ABI}"
)

# Create package directory structure first
dest_dir = os.path.join("cuik_molmaker")
utils_dir = os.path.join(dest_dir, "utils")
data_dir = os.path.join(dest_dir, "data")

# Set appropriate file extensions based on system
if SYSTEM == "Darwin":  # macOS
    so_suffix = f"cpython-{PYTHON_DIGIT_ONLY_VERSION}-darwin.so"
    lib_extension = "dylib"
    lib_dir = os.path.join(dest_dir, "lib")
    lib_file = os.path.join("build", f"libcuik_molmaker_core.{lib_extension}")
elif SYSTEM == "Linux":
    so_suffix = f"cpython-{PYTHON_DIGIT_ONLY_VERSION}-x86_64-linux-gnu.so"
    lib_extension = "so"
    lib_dir = os.path.join(dest_dir, "lib")
    lib_file = os.path.join("build", f"libcuik_molmaker_core.{lib_extension}")
elif SYSTEM == "Windows":
    machine = platform.machine().lower()  # like AMD64
    so_suffix = f"cp{PYTHON_DIGIT_ONLY_VERSION}-win_{machine}.pyd"
    lib_extension = "dll"
    lib_file = os.path.join("build", f"cuik_molmaker_core.{lib_extension}")
    # On Windows, DLLs should be in the same directory or on PATH
    lib_dir = dest_dir
else:
    raise ValueError(f"Unsupported platform: {SYSTEM}")

os.makedirs(dest_dir, exist_ok=True)
os.makedirs(lib_dir, exist_ok=True)
os.makedirs(utils_dir, exist_ok=True)
os.makedirs(data_dir, exist_ok=True)

# Copy all Python files from src/ to package directory
src_py_files = glob.glob(os.path.join("src", "**", "*.py"), recursive=True)
for src_file in src_py_files:
    # Get relative path from src/
    rel_path = os.path.relpath(src_file, "src")
    # Create destination path
    dest_path = os.path.join(dest_dir, rel_path)
    # Create parent directories if they don't exist
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    # Copy the file
    print(f"Copying {src_file} to {dest_path}")
    shutil.copy2(src_file, dest_path)

# Copy data files
data_files = [
    "best_normalization_params.json",
    "fast_normalization_params.json",
    "descriptastorus_normalization_params.json",
    "README.md",
]
for data_file in data_files:
    src_path = os.path.join("data", data_file)
    dest_path = os.path.join(data_dir, data_file)
    if os.path.exists(src_path):
        print(f"Copying {src_path} to {dest_path}")
        shutil.copy2(src_path, dest_path)
    else:
        print(f"WARNING: {src_path} not found")


# Create an empty __init__.py in the lib directory to make it a package
lib_init_file = os.path.join(lib_dir, "__init__.py")
if not os.path.exists(lib_init_file):
    print(f"Creating {lib_init_file}")
    with open(lib_init_file, "w") as f:
        f.write("# This file makes the lib directory a Python package\n")


setup(
    name="cuik_molmaker",
    version="0.2",
    author="S. Veccham",
    author_email="sveccham@nvidia.com",
    description="C++ module for featurizing molecules",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    license="Apache 2.0",
    # Include both src directory and cuik_molmaker_py package
    packages=find_packages(
        include=[
            "cuik_molmaker",
            "cuik_molmaker.lib",
            "cuik_molmaker.utils",
            "cuik_molmaker.data",
        ]
    ),
    package_data={
        "cuik_molmaker": [
            "*.so",
            "*.pyd",
            "*.dll",
            "*.py",
            "data/*.json",
            "data/*.md",
        ],  # Include Python extension and Python files
        "cuik_molmaker.lib": [
            "*.so",
            "__init__.py",
            "*.dll",
            "*.dylib",
        ],  # Include shared libraries and __init__.py
        "cuik_molmaker.utils": ["*.py"],  # Include Python files
    },
    include_package_data=True,
    cmdclass={
        "build_ext": CMakeBuild,
    },
    install_requires=[
        f"rdkit=={RDKIT_VERSION}",
        "scipy",
        "pandas",
    ],
    build_requires=[
        f"rdkit=={RDKIT_VERSION}",
    ],
    tests_require=["pytest"],
    extras_require={
        "dev": [
            "black>=24.2.0",
            "flake8>=7.3.0",
            "isort>=5.13.2",
            "pre-commit>=3.6.0",
        ],
    },
    python_requires=f"=={PYTHON_VERSION}.*",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        f"Programming Language :: Python :: {PYTHON_VERSION}",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    entry_points={
        "console_scripts": [
            "cuik-molmaker-fit-distribution=cuik_molmaker.utils.fit_distribution:main",
            "cuik-molmaker-mol-features=cuik_molmaker.mol_features:main",
        ],
    },
)

# Check if compiled extension exists, display helpful message if not
so_file = os.path.join("build", f"cuik_molmaker_cpp.{so_suffix}")
print(f"Looking for compiled extension at: {so_file}")

if os.path.exists(so_file):
    print(f"Found compiled extension, copying to {dest_dir}")
    shutil.copy2(so_file, dest_dir)
else:
    print(
        "WARNING: Compiled extension not found."
        "You need to build the C++ extension first."
    )
    sys.exit(1)

# Check for the shared library
if os.path.exists(lib_file):
    print(f"Found shared library, copying to {lib_dir}")
    shutil.copy2(lib_file, lib_dir)
else:
    print(
        "WARNING: Shared library not found. You need to build the C++ extension first."
    )
    print("Run cmake and make in the build directory first.")
    sys.exit(1)
