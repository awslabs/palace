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

namespace palace::utils
{

// Set the number of OpenMP threads to be used for parallel regions.
void SetNumThreads(int nt);

// Return maximum number of OpenMP threads.
int GetMaxThreads();

// Return number of active OpenMP threads.
int GetNumActiveThreads();

// Return the current thread ID.
int GetThreadNum();

// Return whether or not the current scope is inside a parallel OpenMP region.
int InParallel();

// Set and return the number of OpenMP threads depending on OMP_NUM_THREADS.
int ConfigureOmp();

}  // namespace palace::utils

#endif  // PALACE_UTILS_OMP_HPP
