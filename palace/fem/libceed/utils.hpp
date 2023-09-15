// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_UTILS_HPP
#define PALACE_LIBCEED_UTILS_HPP

#include <functional>
#include <string>
#include <vector>
#include <ceed.h>
#include <mfem.hpp>

namespace palace::ceed
{

namespace internal
{

extern std::vector<Ceed> ceed;

}  // namespace internal

#define PalaceCeedCall(ceed, ...)      \
  do                                   \
  {                                    \
    int ierr_ = __VA_ARGS__;           \
    if (ierr_ != CEED_ERROR_SUCCESS)   \
    {                                  \
      const char *msg;                 \
      CeedGetErrorMessage(ceed, &msg); \
      MFEM_ABORT(msg);                 \
    }                                  \
  } while (0)

#define PalaceCeedCallBackend(...)                      \
  do                                                    \
  {                                                     \
    int ierr_ = __VA_ARGS__;                            \
    if (ierr_ != CEED_ERROR_SUCCESS)                    \
    {                                                   \
      MFEM_ABORT("libCEED encountered a fatal error!"); \
    }                                                   \
  } while (0)

#define PalaceQFunctionRelativePath(path) strstr(path, "qfunctions")

#if defined(MFEM_USE_OPENMP)
#define PalacePragmaOmpHelper(x) _Pragma(#x)
#define PalacePragmaOmp(x) PalacePragmaOmpHelper(omp x)
#else
#define PalacePragmaOmp(x)
#endif

// Call libCEED's CeedInit for the given resource. The specific device to use is set prior
// to this using mfem::Device.
void Initialize(const char *resource, const char *jit_source_dir);

// Finalize libCEED with CeedDestroy.
void Finalize();

// Get the configured libCEED backend.
std::string Print();

// Wrapper for std::hash.
template <typename T>
inline std::size_t CeedHash(const T key)
{
  return std::hash<T>{}(key);
}

// Effective way to combine hashes (from libCEED).
inline std::size_t CeedHashCombine(std::size_t seed, std::size_t hash)
{
  // See https://doi.org/10.1002/asi.10170 or https://dl.acm.org/citation.cfm?id=759509.
  return seed ^ (hash + (seed << 6) + (seed >> 2));
}

// Initialize a CeedVector from an mfem::Vector
void InitCeedVector(const mfem::Vector &v, Ceed ceed, CeedVector *cv);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_UTILS_HPP
