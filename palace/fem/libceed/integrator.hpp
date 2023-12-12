// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_INTEGRATOR_HPP
#define PALACE_LIBCEED_INTEGRATOR_HPP

#include <string>
#include <vector>
#include "fem/libceed/ceed.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"

namespace palace::ceed
{

// Geometry factor information as quadrature data.
enum GeomFactorInfo : unsigned int
{
  Determinant = 1 << 0,
  Adjugate = 1 << 1,
  Weight = 1 << 2
};

// Evaluation modes for CeedOperator fields for various integrators.
enum EvalMode : unsigned int
{
  None = 1 << 0,
  Interp = 1 << 1,
  Grad = 1 << 2,
  Div = 1 << 3,
  Curl = 1 << 4
};

// Data structure for CeedOperator construction for various integrators.
struct IntegratorInfo
{
  // QFunctions for operator construction and application.
  CeedQFunctionUser apply_qf;

  // Path and name of the QFunctions for operator construction and application.
  std::string apply_qf_path;

  // Geometry factor required as QFunction inputs.
  unsigned int geom_info;

  // Evaluation modes for the test and trial basis.
  unsigned int trial_ops, test_ops;

  // Control whether or not to pre-assemble the quadrature data or compute it during
  // operator application in true matrix-free fashion.
  bool assemble_qdata;

  IntegratorInfo()
    : apply_qf(nullptr), apply_qf_path(""), geom_info(0), trial_ops(0), test_ops(0),
      assemble_qdata(false)
  {
  }
};

// Create libCEED quadrature data and element restriction for use in a partially assembled
// libCEED operator.
void AssembleCeedGeometryData(Ceed ceed, CeedElemRestriction mesh_restr,
                              CeedBasis mesh_basis, const Vector &mesh_nodes,
                              CeedGeomFactorData &geom_data);

// Create libCEED operator using the given quadrature data, element restriction, and basis
// objects.
void AssembleCeedOperator(const IntegratorInfo &info, void *ctx, std::size_t ctx_size,
                          const CeedGeomFactorData &geom_data, Ceed ceed,
                          CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                          CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op);

// Construct libCEED operators for interpolation operations and their transpose between
// the two spaces. Note that contributions for shared degrees of freedom are added, so the
// output of the operator application must be scaled by the inverse multiplicity.
void AssembleCeedInterpolator(Ceed ceed, CeedElemRestriction trial_restr,
                              CeedElemRestriction test_restr, CeedBasis interp_basis,
                              CeedOperator *op, CeedOperator *op_t);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_INTEGRATOR_HPP
