# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved. # noqa: E501
# SPDX-License-Identifier: Apache-2.0

import numpy as np
import pytest

from cuik_molmaker import best_fit_distribution
from cuik_molmaker.utils import get_discrete_cdf


def test_get_discrete_cdf():
    x = np.arange(1, 6)
    y = np.array([0.1, 0.2, 0.3, 0.2, 0.2])
    cdf = get_discrete_cdf(x, y)

    discrete_cdf_ref = np.array([0.1, 0.3, 0.6, 0.8, 1.0])
    np.testing.assert_allclose(np.array(cdf[1]), discrete_cdf_ref)

    # Case where there are missing values in the middle in x
    x = np.array([1, 2, 3, 4, 5, 8, 9, 10])
    y = np.array([0.1, 0.2, 0.3, 0.05, 0.05, 0.1, 0.1, 0.1])
    cdf = get_discrete_cdf(x, y)
    discrete_cdf_ref = np.array([0.1, 0.3, 0.6, 0.65, 0.7, 0.7, 0.7, 0.8, 0.9, 1.0])
    np.testing.assert_allclose(np.array(cdf[1]), discrete_cdf_ref)


@pytest.mark.slow
def test_best_fit_distribution_continuous(rdkit2D_desc_df):

    data = rdkit2D_desc_df["SPS"].values
    best_dist = best_fit_distribution(data, "continuous")

    sse_array = np.array([dist[2] for dist in best_dist])
    assert np.all(
        np.diff(sse_array) >= 0
    ), f"SSE array is not monotonic increasing: {sse_array}"

    assert (
        best_dist[0][0] == "laplace_asymmetric"
    ), f"Best distribution is not laplace_asymmetric: {best_dist[0][0]}"
    assert (
        best_dist[1][0] == "exponnorm"
    ), f"Second best distribution is not exponnorm: {best_dist[1][0]}"
    assert (
        best_dist[2][0] == "skewnorm"
    ), f"Second best distribution is not skewnorm: {best_dist[2][0]}"


@pytest.mark.slow
def test_best_fit_distribution_discrete(rdkit2D_desc_df):

    data = rdkit2D_desc_df["fr_Ar_N"].values
    best_dist = best_fit_distribution(data, "discrete")
    best_dist_param_ref = np.array(
        [
            0.2813,
            0.5066999999999999,
            0.7699999999999999,
            0.9027,
            0.9794999999999999,
            0.9955999999999999,
            0.9994999999999999,
            0.9999999999999999,
        ]
    )
    np.testing.assert_allclose(
        best_dist[1],
        best_dist_param_ref,
        err_msg="Parameters do not match reference: "
        f"{best_dist[1]} != {best_dist_param_ref}",
    )
