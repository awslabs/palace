// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_CEED_HPP
#define PALACE_LIBCEED_CEED_HPP

#include <string>
#include <unordered_map>
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

// Useful alias templates for libCEED objects specific to a specific Ceed context and
// element geometry type.
template <typename T>
using GeometryObjectMap = std::unordered_map<mfem::Geometry::Type, T>;
template <typename T>
using CeedObjectMap = std::unordered_map<Ceed, GeometryObjectMap<T>>;

// Call libCEED's CeedInit for the given resource. The specific device to use is set prior
// to this using mfem::Device.
void Initialize(const char *resource, const char *jit_source_dir);

// Finalize libCEED with CeedDestroy.
void Finalize();

// Get the configured libCEED backend.
std::string Print();

// Initialize a CeedVector from an mfem::Vector. When init is false, expects the CeedVector
// has already been initialized and just sets the data pointer.
void InitCeedVector(const mfem::Vector &v, Ceed ceed, CeedVector *cv, bool init = true);

// Convert an MFEM geometry type to a libCEED one.
CeedElemTopology GetCeedTopology(mfem::Geometry::Type geom);

// Convert a libCEED geometry type to an MFEM one.
mfem::Geometry::Type GetMfemTopology(CeedElemTopology geom);

namespace internal
{

// Access the Ceed objects initialized by CeedInit.
const std::vector<Ceed> &GetCeedObjects();

// Convenience method for number of ceeds.
std::size_t NumCeeds();

}  // namespace internal

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_CEED_HPP
