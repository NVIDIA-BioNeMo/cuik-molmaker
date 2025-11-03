# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved. # noqa: E501
# SPDX-License-Identifier: Apache-2.0

import os
import platform

# Import compiled extension
from pathlib import Path

from .mol_features import MoleculeFeaturizer
from .utils.descriptor_normalization import get_normalization_functions
from .utils.fit_distribution import best_fit_distribution, get_fast_distribution

platform_name = platform.system()
if platform_name == "Windows":
    extension = ".pyd"
    string_to_match = "cp"

    # On Windows, add DLL directories to search path
    # This is needed so that cuik_molmaker_core.dll can find its RDKit dependencies

    # At runtime, we only need DLLs from rdkit.libs (PyPI RDKit)
    # The .lib files from header environment are only needed at build/link time
    try:
        import rdkit

        rdkit_path = Path(rdkit.__file__).parent

        # Windows pip-installed RDKit stores DLLs in site-packages/rdkit.libs
        # These are the same DLLs we linked against at build time
        rdkit_libs = rdkit_path.parent / "rdkit.libs"
        if rdkit_libs.exists():
            if hasattr(os, "add_dll_directory"):
                os.add_dll_directory(str(rdkit_libs))
            else:
                os.environ["PATH"] = (
                    str(rdkit_libs) + os.pathsep + os.environ.get("PATH", "")
                )
    except ImportError:
        # rdkit not found, continue anyway (will fail later if actually needed)
        pass
else:
    extension = ".so"
    string_to_match = "cpython"

# Find the .so file in this directory
_module_dir = Path(__file__).parent
for file in os.listdir(_module_dir):
    if file.endswith(extension) and string_to_match in file:
        # Add the extension module directly
        from importlib.machinery import ExtensionFileLoader
        from importlib.util import module_from_spec, spec_from_loader

        _loader = ExtensionFileLoader("cuik_molmaker_cpp", str(_module_dir / file))
        _spec = spec_from_loader("cuik_molmaker_cpp", _loader)
        _module = module_from_spec(_spec)
        _loader.exec_module(_module)

        # Import all attributes from the module
        for attr in dir(_module):
            if not attr.startswith("_"):
                globals()[attr] = getattr(_module, attr)
        break


__all__ = [
    # mol_features functions
    "MoleculeFeaturizer",
    # fit_distribution functions
    "best_fit_distribution",
    "get_fast_distribution",
    # descriptor_normalization functions
    "get_normalization_functions",
]
