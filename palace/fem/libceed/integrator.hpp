// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_INTEGRATOR_HPP
#define PALACE_LIBCEED_INTEGRATOR_HPP

#include <string>
#include <vector>
#include "fem/libceed/ceed.hpp"

namespace palace::ceed
{

// Evaluation modes for CeedOperator fields for various integrators.
enum EvalMode : unsigned int
{
  Weight = 1 << 0,
  None = 1 << 1,
  Interp = 1 << 2,
  Grad = 1 << 3,
  Div = 1 << 4,
  Curl = 1 << 5
};

// Data structure for CeedOperator construction for various integrators.
struct CeedQFunctionInfo
{
  // QFunctions for operator construction and application.
  CeedQFunctionUser apply_qf;

  // Path and name of the QFunctions for operator construction and application.
  std::string apply_qf_path;

  // Evaluation modes for the test and trial basis.
  unsigned int trial_ops, test_ops;

  // Control whether or not to pre-assemble the quadrature data or compute it during
  // operator application in true matrix-free fashion.
  bool assemble_q_data;

  CeedQFunctionInfo()
    : apply_qf(nullptr), apply_qf_path(""), trial_ops(0), test_ops(0),
      assemble_q_data(false)
  {
  }
};

// Helper function to get the geometry space dimension.
int CeedGeometryDataGetSpaceDimension(CeedElemRestriction geom_data_restr, CeedInt dim,
                                      CeedInt *space_dim);

// Assemble libCEED mesh geometry factor quadrature data for use in a partially assembled
// libCEED operator.
void AssembleCeedGeometryData(Ceed ceed, CeedElemRestriction mesh_restr,
                              CeedBasis mesh_basis, CeedVector mesh_nodes,
                              CeedElemRestriction attr_restr, CeedBasis attr_basis,
                              CeedVector elem_attr, CeedVector geom_data,
                              CeedElemRestriction geom_data_restr);

// Construct libCEED operator using the given quadrature data, element restriction, and
// basis objects.
void AssembleCeedOperator(const CeedQFunctionInfo &info, void *ctx, std::size_t ctx_size,
                          Ceed ceed, CeedElemRestriction trial_restr,
                          CeedElemRestriction test_restr, CeedBasis trial_basis,
                          CeedBasis test_basis, CeedVector geom_data,
                          CeedElemRestriction geom_data_restr, CeedOperator *op);

// Construct libCEED operators for interpolation operations and their transpose between
// the two spaces. Note that contributions for shared degrees of freedom are added, so the
// output of the operator application must be scaled by the inverse multiplicity.
void AssembleCeedInterpolator(Ceed ceed, CeedElemRestriction trial_restr,
                              CeedElemRestriction test_restr, CeedBasis interp_basis,
                              CeedOperator *op, CeedOperator *op_t);

// Construct a libCEED operator which integrates the squared difference between two
// functions over every element.
void AssembleCeedElementErrorIntegrator(
    const CeedQFunctionInfo &info, void *ctx, std::size_t ctx_size, Ceed ceed,
    CeedVector input1, CeedVector input2, CeedElemRestriction input1_restr,
    CeedElemRestriction input2_restr, CeedBasis input1_basis, CeedBasis input2_basis,
    CeedElemRestriction mesh_elem_restr, CeedVector geom_data,
    CeedElemRestriction geom_data_restr, CeedOperator *op);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_INTEGRATOR_HPP
