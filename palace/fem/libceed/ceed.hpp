// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_CEED_HPP
#define PALACE_LIBCEED_CEED_HPP

#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>
#include <mfem.hpp>

// XX TODO COMBINE WITH libceed/utils.hpp AND SKIP FORWARD DECLARATIONS OF LIBCEED TYPES?

// Forward declarations of libCEED objects.
typedef struct Ceed_private *Ceed;
typedef struct CeedVector_private *CeedVector;
typedef struct CeedElemRestriction_private *CeedElemRestriction;
typedef struct CeedBasis_private *CeedBasis;
typedef struct CeedOperator_private *CeedOperator;

namespace palace::ceed
{

// XX TODO: Do we need to store qw * |J| separately in each quadrature data? Is it
//          significantly worse if we just use multiple inputs to the QFunction for the
//          different quantities?

// XX TODO: Can we skip adjugate storage and just compute from J on the fly? Or, probably
//          better, can we skip J storage and compute from adj(adj(J)/|J|) = J?

// XX TODO RENAME? CeedMeshData and SetUpCeedMeshData?

// Data structure for geometry information stored at quadrature points. Jacobian matrix is
// dim x space_dim, the adjugate is space_dim x dim, column-major storage by component.
struct CeedGeomFactorData_private
{
  // Dimension and space dimension for this element topology.
  int dim, space_dim;

  // Element indices from the mfem::Mesh used to construct Ceed objects with these geometry
  // factors.
  std::vector<int> indices;

  mfem::Vector wdetJ;  // qw * |J|, for H1 conformity with quadrature weights
  mfem::Vector adjJt;  // adj(J)^T / |J|, for H(curl) conformity
  mfem::Vector J;      // J / |J|, for H(div) conformity
  mfem::Vector attr;   // Mesh element attributes

  // Objects for libCEED interface to the quadrature data.
  CeedVector wdetJ_vec, adjJt_vec, J_vec, attr_vec;
  CeedElemRestriction wdetJ_restr, adjJt_restr, J_restr, attr_restr;
  Ceed ceed;

  CeedGeomFactorData_private(Ceed ceed)
    : dim(0), space_dim(0), wdetJ_vec(nullptr), adjJt_vec(nullptr), J_vec(nullptr),
      attr_vec(nullptr), wdetJ_restr(nullptr), adjJt_restr(nullptr), J_restr(nullptr),
      attr_restr(nullptr), ceed(ceed)
  {
  }
  ~CeedGeomFactorData_private();
};

using CeedGeomFactorData = std::unique_ptr<CeedGeomFactorData_private>;

inline CeedGeomFactorData CeedGeomFactorDataCreate(Ceed ceed)
{
  return std::make_unique<CeedGeomFactorData_private>(ceed);
}

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

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_OPERATOR_HPP
