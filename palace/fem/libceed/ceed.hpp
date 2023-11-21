// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_CEED_HPP
#define PALACE_LIBCEED_CEED_HPP

#include <unordered_map>
#include <utility>
#include <mfem.hpp>

// Forward declarations of libCEED objects.
typedef struct Ceed_private *Ceed;
typedef struct CeedVector_private *CeedVector;
typedef struct CeedElemRestriction_private *CeedElemRestriction;
typedef struct CeedBasis_private *CeedBasis;
typedef struct CeedOperator_private *CeedOperator;

namespace palace::ceed
{

// Useful alias template for libCEED objects specific to a specific Ceed context and element
// geometry type.
template <typename T>
using CeedObjectMap = std::unordered_map<std::pair<Ceed, mfem::Geometry::Type>, T>;

// XX TODO: Do we need to store qw * |J| separately in each? Is it significantly worse if we
//          just use multiple inputs to the QFunction for the different quantities?
// XX TODO: Can we skip adjugate storage and just compute from J on the fly?

// Data structure for geometry information stored at quadrature points. Jacobian matrix is
// dim x space_dim, the adjugate is space_dim x dim, column-major storage by component.
struct CeedGeomFactorData
{
  mfem::Vector wdetJ;  // qw * |J|, H1 conformity
  mfem::Vector J;      // J / |J|, H(curl) conformity
  mfem::Vector adjJt;  // adj(J)^T / |J|, H(div) conformity
  mfem::Vector attr;   // Mesh element attributes

  // Objects for libCEED interface to the quadrature data.
  CeedVector wdetJ_vec, J_vec, adjJt_vec, attr_vec;
  CeedElemRestriction wdetJ_restr, J_restr, adjJt_restr, attr_restr;

  CeedGeomFactorData()
    : wdetJ_vec(nullptr), J_vec(nullptr), adjJt_vec(nullptr), attr_vec(nullptr),
      wdetJ_restr(nullptr), J_restr(nullptr), adjJt_restr(nullptr), attr_restr(nullptr)
  {
  }
};

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_OPERATOR_HPP
