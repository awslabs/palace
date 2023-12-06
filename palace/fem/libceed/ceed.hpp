// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_CEED_HPP
#define PALACE_LIBCEED_CEED_HPP

#include <string>
#include <unordered_map>
#include <utility>
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

// Base case for combining hashes.
inline void CeedHashCombine(std::size_t &seed) {}

// See for example https://onlinelibrary.wiley.com/doi/abs/10.1002/asi.10170, the source
// of https://www.boost.org/doc/libs/1_35_0/doc/html/boost/hash_combine_id241013.html.
template <typename T, typename... U>
inline void CeedHashCombine(std::size_t &seed, const T &v, const U &...args)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  (CeedHashCombine(seed, args), ...);
}

// Hash function for the CeedObjectMap type.
struct CeedObjectHash
{
  std::size_t operator()(const std::pair<Ceed, mfem::Geometry::Type> &k) const
  {
    std::size_t hash = 0;
    CeedHashCombine(hash, k.first, k.second);
    return hash;
  }
};

// Useful alias template for libCEED objects specific to a specific Ceed context and element
// geometry type.
template <typename T>
using CeedObjectMap =
    std::unordered_map<std::pair<Ceed, mfem::Geometry::Type>, T, CeedObjectHash>;

// Call libCEED's CeedInit for the given resource. The specific device to use is set prior
// to this using mfem::Device.
void Initialize(const char *resource, const char *jit_source_dir);

// Finalize libCEED with CeedDestroy.
void Finalize();

// Get the configured libCEED backend.
std::string Print();

// Initialize a CeedVector from an mfem::Vector.
void InitCeedVector(const mfem::Vector &v, Ceed ceed, CeedVector *cv);

// Convert an MFEM geometry type to a libCEED one.
CeedElemTopology GetCeedTopology(mfem::Geometry::Type geom);

namespace internal
{

// Access the Ceed objects initialized by CeedInit.
const std::vector<Ceed> &GetCeedObjects();

}  // namespace internal

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_OPERATOR_HPP
