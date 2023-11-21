// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_UTILS_HPP
#define PALACE_LIBCEED_UTILS_HPP

#include <string>
#include <vector>
#include <ceed.h>
#include <mfem.hpp>

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

namespace palace::ceed
{

// Call libCEED's CeedInit for the given resource. The specific device to use is set prior
// to this using mfem::Device.
void Initialize(const char *resource, const char *jit_source_dir);

// Finalize libCEED with CeedDestroy.
void Finalize();

// Get the configured libCEED backend.
std::string Print();

// Initialize a CeedVector from an mfem::Vector.
void InitCeedVector(const mfem::Vector &v, Ceed ceed, CeedVector *cv);

namespace internal
{

// Access the Ceed objects initialized by CeedInit.
const std::vector<Ceed> &GetCeedObjects();

}  // namespace internal

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_UTILS_HPP
