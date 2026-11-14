# Available Builds

`cuik-molmaker` is published through three channels. The channels have different package names, version schemes, and compatibility constraints, so each has its own build matrix.


| Channel                                                                                  | Package             | Version scheme                 | Choose this channel when                                                             |
| ---------------------------------------------------------------------------------------- | ------------------- | ------------------------------ | ------------------------------------------------------------------------------------ |
| [conda-forge](https://anaconda.org/channels/conda-forge/packages/cuik_molmaker/overview) | `cuik_molmaker`     | cuik-molmaker semantic version | You install RDKit and cuik-molmaker with conda.                                      |
| [PyPI](https://pypi.org/project/cuik-molmaker-pin/)                                      | `cuik-molmaker-pin` | RDKit date version             | You install with pip and want RDKit to be pinned automatically.                      |
| [NVIDIA PyPI](https://pypi.nvidia.com/rdkit-2026.09.2/cuik-molmaker/)                    | `cuik_molmaker`     | cuik-molmaker semantic version | You use RDKit from PyPI and you want a specific semantic version of `cuik-molmaker`. |




## Conda-forge

Conda-forge package versions use cuik-molmaker semantic versioning. Each entry
below lists every published RDKit build variant for the cuik-molmaker release.


| cuik-molmaker | RDKit builds                                                                                                            | Linux x86_64                 | Linux arm64            | macOS x86_64                      | macOS arm64            | Windows x86_64               |
| ------------- | ----------------------------------------------------------------------------------------------------------------------- | ---------------------------- | ---------------------- | --------------------------------- | ---------------------- | ---------------------------- |
| 0.3.1         | 2025.09.1, 2025.09.3, 2025.09.6, 2026.03.1, 2026.03.2, 2026.03.3, 2026.03.4, 2026.03.5, 2026.03.6, 2026.09.1, 2026.09.2 | 3.11, 3.12, 3.13, 3.14       | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 (see note) | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14       |
| 0.3.0         | 2025.09.6, 2026.03.1, 2026.03.2, 2026.03.3, 2026.03.4                                                                   | 3.11, 3.12, 3.13, 3.14       | —                      | 3.11, 3.12, 3.13, 3.14            | —                      | 3.11, 3.12, 3.13, 3.14       |
| 0.2.1         | 2025.09.2, 2025.09.3, 2025.09.4, 2025.09.6, 2026.03.1, 2026.03.2, 2026.03.3                                             | 3.10, 3.11, 3.12, 3.13, 3.14 | —                      | 3.10, 3.11, 3.12, 3.13, 3.14      | —                      | 3.10, 3.11, 3.12, 3.13, 3.14 |
| 0.2           | 2025.03.5, 2025.03.6, 2025.09.1                                                                                         | 3.10, 3.11, 3.12, 3.13       | —                      | 3.10, 3.11, 3.12, 3.13            | —                      | 3.10, 3.11, 3.12, 3.13       |

Note: the macOS x86_64 / Python 3.14 build for RDKit 2026.09.2 is not available yet, because conda-forge has no `rdkit` 2026.09.2 package for that combination.


Install the latest compatible build:

```bash
conda install -c conda-forge cuik_molmaker
```



## PyPI: pinned RDKit builds

`cuik-molmaker-pin` uses the RDKit release version as its own version and
requires that exact `rdkit` version. For example,
`cuik-molmaker-pin==2026.9.2` requires `rdkit==2026.09.2`.


| Pin / RDKit version | Linux x86_64           | Linux arm64            | macOS x86_64     | macOS arm64                           | Windows x86_64         |
| ------------------- | ---------------------- | ---------------------- | ---------------- | ------------------------------------- | ---------------------- |
| 2026.9.2            | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 | —                | 3.11, 3.12, 3.13, 3.14                | 3.11, 3.12, 3.13, 3.14 |
| 2026.9.1            | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 | —                | 3.11, 3.12, 3.13, 3.14                | 3.11, 3.12, 3.13, 3.14 |
| 2026.3.6            | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 | —                | 3.11, 3.12, 3.13, 3.14                | 3.11, 3.12, 3.13, 3.14 |
| 2026.3.5            | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 | —                | 3.11, 3.12, 3.13, 3.14                | 3.11, 3.12, 3.13, 3.14 |
| 2026.3.4            | 3.11, 3.12, 3.13, 3.14 | —                      | —                | 3.11, 3.12, 3.13, 3.14                | 3.11, 3.12, 3.13, 3.14 |
| 2026.3.3            | 3.11, 3.12, 3.13       | —                      | —                | 3.11, 3.12, 3.13 (`.post1`, `.post2`) | 3.11, 3.12, 3.13       |
| 2026.3.2            | 3.11, 3.12, 3.13       | —                      | —                | —                                     | 3.11, 3.12, 3.13       |
| 2026.3.1            | 3.11, 3.12, 3.13       | —                      | —                | —                                     | 3.11, 3.12, 3.13       |
| 2025.9.6            | 3.11, 3.12, 3.13       | —                      | —                | —                                     | 3.11, 3.12, 3.13       |
| 2025.9.4            | 3.11, 3.12, 3.13       | —                      | —                | —                                     | 3.11, 3.12, 3.13       |
| 2025.9.3            | 3.11, 3.12, 3.13       | —                      | 3.11, 3.12, 3.13 | —                                     | 3.11, 3.12, 3.13       |
| 2025.9.2            | 3.11, 3.12, 3.13       | —                      | 3.11, 3.12, 3.13 | —                                     | 3.11, 3.12, 3.13       |


Install a pinned build:

```bash
python -m pip install "cuik-molmaker-pin==2026.9.2"
```



## NVIDIA PyPI

NVIDIA PyPI groups wheels by the compatible RDKit version in different subdirectories in the [NVIDIA PyPI index](https://pypi.nvidia.com/).


| NVIDIA PyPI subdirectory                                                  | cuik-molmaker | Linux x86_64           | Linux arm64            | macOS x86_64     | macOS arm64            | Windows x86_64         |
| ------------------------------------------------------------------------- | ------------- | ---------------------- | ---------------------- | ---------------- | ---------------------- | ---------------------- |
| [rdkit-2026.09.2](https://pypi.nvidia.com/rdkit-2026.09.2/cuik-molmaker/) | 0.3.1         | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 | —                | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 |
| [rdkit-2026.09.1](https://pypi.nvidia.com/rdkit-2026.09.1/cuik-molmaker/) | 0.3.1         | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 | —                | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 |
| [rdkit-2026.03.6](https://pypi.nvidia.com/rdkit-2026.03.6/cuik-molmaker/) | 0.3.1         | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 | —                | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 |
| [rdkit-2026.03.5](https://pypi.nvidia.com/rdkit-2026.03.5/cuik-molmaker/) | 0.3.1         | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 | —                | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 |
| [rdkit-2026.03.4](https://pypi.nvidia.com/rdkit-2026.03.4/cuik-molmaker/) | 0.3.1         | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 | —                | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 |
| [rdkit-2026.03.4](https://pypi.nvidia.com/rdkit-2026.03.4/cuik-molmaker/) | 0.3.0         | 3.11, 3.12, 3.13, 3.14 | —                      | —                | 3.11, 3.12, 3.13, 3.14 | 3.11, 3.12, 3.13, 3.14 |
| [rdkit-2026.03.2](https://pypi.nvidia.com/rdkit-2026.03.2/cuik-molmaker/) | 0.2.1         | 3.11, 3.12, 3.13       | —                      | —                | —                      | 3.11, 3.12, 3.13       |
| [rdkit-2026.03.1](https://pypi.nvidia.com/rdkit-2026.03.1/cuik-molmaker/) | 0.2.1         | 3.11, 3.12, 3.13       | —                      | —                | —                      | 3.11, 3.12, 3.13       |
| [rdkit-2025.09.6](https://pypi.nvidia.com/rdkit-2025.09.6/cuik-molmaker/) | 0.2.1         | 3.11, 3.12, 3.13       | —                      | —                | —                      | 3.11, 3.12, 3.13       |
| [rdkit-2025.09.4](https://pypi.nvidia.com/rdkit-2025.09.4/cuik-molmaker/) | 0.2.1         | 3.11, 3.12, 3.13       | —                      | —                | —                      | 3.11, 3.12, 3.13       |
| [rdkit-2025.09.3](https://pypi.nvidia.com/rdkit-2025.09.3/cuik-molmaker/) | 0.2.1         | 3.11, 3.12, 3.13       | —                      | 3.11, 3.12, 3.13 | —                      | 3.11, 3.12, 3.13       |
| [rdkit-2025.09.3](https://pypi.nvidia.com/rdkit-2025.09.3/cuik-molmaker/) | 0.2           | 3.11, 3.12, 3.13       | —                      | —                | —                      | —                      |
| [rdkit-2025.09.2](https://pypi.nvidia.com/rdkit-2025.09.2/cuik-molmaker/) | 0.2.1         | 3.11, 3.12, 3.13       | —                      | 3.11, 3.12, 3.13 | —                      | 3.11, 3.12, 3.13       |
| [rdkit-2025.09.2](https://pypi.nvidia.com/rdkit-2025.09.2/cuik-molmaker/) | 0.2           | 3.11, 3.12, 3.13       | —                      | —                | —                      | —                      |
| [rdkit-2025.09.1](https://pypi.nvidia.com/rdkit-2025.09.1/cuik-molmaker/) | 0.2           | 3.11, 3.12, 3.13       | —                      | 3.11, 3.12, 3.13 | 3.11, 3.12, 3.13       | 3.11, 3.12, 3.13       |
| [rdkit-2025.03.6](https://pypi.nvidia.com/rdkit-2025.03.6/cuik-molmaker/) | 0.2           | 3.11, 3.12, 3.13       | —                      | 3.11, 3.12, 3.13 | 3.11, 3.12, 3.13       | 3.11, 3.12, 3.13       |
| [rdkit-2025.03.5](https://pypi.nvidia.com/rdkit-2025.03.5/cuik-molmaker/) | 0.2           | 3.11, 3.12, 3.13       | —                      | 3.11, 3.12, 3.13 | 3.11, 3.12, 3.13       | 3.11, 3.12, 3.13       |


For example, to install `cuik-molmaker` `0.3.1` compatible with `rdkit-2026.09.2`, use the following command:

```bash
python -m pip install "cuik-molmaker==0.3.1" --extra-index-url https://pypi.nvidia.com/rdkit-2026.09.2/
```

