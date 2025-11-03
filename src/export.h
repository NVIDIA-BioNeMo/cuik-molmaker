// SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#ifdef _WIN32
#ifdef cuik_molmaker_core_EXPORTS
// We are building the DLL
#define CUIK_EXPORT __declspec(dllexport)
#else
// We are using the DLL
#define CUIK_EXPORT __declspec(dllimport)
#endif
#else
// On non-Windows platforms, no special export needed
#define CUIK_EXPORT __attribute__((visibility("default")))
#endif