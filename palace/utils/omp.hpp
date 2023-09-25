// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_OMP_HPP
#define PALACE_UTILS_OMP_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_OPENMP)
#define PalacePragmaOmpHelper(x) _Pragma(#x)
#define PalacePragmaOmp(x) PalacePragmaOmpHelper(omp x)
#else
#define PalacePragmaOmp(x)
#endif

#endif  // PALACE_UTILS_OMP_HPP
