// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "omp.hpp"

#if defined(MFEM_USE_OPENMP)
#include <omp.h>
#endif

namespace palace::utils
{

void SetNumThreads(int nt)
{
#if defined(MFEM_USE_OPENMP)
  omp_set_num_threads(nt);
#endif
}

int GetMaxThreads()
{
#if defined(MFEM_USE_OPENMP)
  return omp_get_max_threads();
#else
  return 1;
#endif
}

int GetNumActiveThreads()
{
#if defined(MFEM_USE_OPENMP)
  return omp_get_num_threads();
#else
  return 1;
#endif
}

int GetThreadNum()
{
#if defined(MFEM_USE_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}

int InParallel()
{
#if defined(MFEM_USE_OPENMP)
  return omp_in_parallel();
#else
  return 0;
#endif
}

}  // namespace palace::utils
