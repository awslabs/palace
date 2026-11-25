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

  // Whether the assembled quadrature data couples a pair of 3-component active fields
  // across, (u₁, u₂) ↦ (qd₁ u₂, qd₂ u₁), rather than each field to itself (see
  // f_apply_33_cross).
  bool q_data_cross;

  CeedQFunctionInfo()
    : apply_qf(nullptr), apply_qf_path(""), trial_ops(0), test_ops(0),
      assemble_q_data(false), q_data_cross(false)
  {
  }
};

// Number of components of the geometry factor quadrature data {attr, w |J|, adj(J)ᵀ / |J|},
// followed with coords by the physical coordinates of the quadrature points (see
// Mesh::GetCeedGeomFactorData).
constexpr CeedInt CeedGeometryDataSize(CeedInt space_dim, CeedInt dim, bool coords = false)
{
  return 2 + space_dim * dim + (coords ? space_dim : 0);
}

// Helper function to get the geometry space dimension, for geometry factor quadrature data
// without coordinates.
int CeedGeometryDataGetSpaceDimension(CeedElemRestriction geom_data_restr, CeedInt dim,
                                      CeedInt *space_dim);

// Assemble libCEED mesh geometry factor quadrature data for use in a partially assembled
// libCEED operator, with the quadrature point coordinates if coords (see
// CeedGeometryDataSize).
void AssembleCeedGeometryData(Ceed ceed, CeedElemRestriction mesh_restr,
                              CeedBasis mesh_basis, CeedVector mesh_nodes,
                              CeedElemRestriction attr_restr, CeedBasis attr_basis,
                              CeedVector elem_attr, CeedVector geom_data,
                              CeedElemRestriction geom_data_restr, bool coords = false);

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

// Construct a libCEED operator for the derivative of a form a(u, w) with respect to the
// positions of the mesh nodes. The fields u and w are passive inputs with the given
// restriction and basis (u is omitted if info.trial_ops has no evaluation mode), and the
// output is on the mesh nodal restriction and basis.
void AssembleCeedShapeIntegrator(const CeedQFunctionInfo &info, void *ctx,
                                 std::size_t ctx_size, Ceed ceed, CeedVector u,
                                 CeedVector w, CeedElemRestriction restr, CeedBasis basis,
                                 CeedElemRestriction mesh_restr, CeedBasis mesh_basis,
                                 CeedVector geom_data, CeedElemRestriction geom_data_restr,
                                 CeedOperator *op);

// Construct a libCEED operator which integrates the squared difference between two
// functions over every element.
void AssembleCeedElementErrorIntegrator(
    const CeedQFunctionInfo &info, void *ctx, std::size_t ctx_size, Ceed ceed,
    CeedVector input1, CeedVector input2, CeedElemRestriction input1_restr,
    CeedElemRestriction input2_restr, CeedBasis input1_basis, CeedBasis input2_basis,
    CeedElemRestriction mesh_elem_restr, CeedVector geom_data,
    CeedElemRestriction geom_data_restr, CeedOperator *op);

// Construct a libCEED operator which integrates a function of several fields, each a
// passive input named u_1, u_2, ... with its own restriction and basis, over every element
// (for example a complex flux error from the real and imaginary parts of two fields). The
// evaluation modes of every input are info.trial_ops. If ctx_out is non-null, the
// QFunction context handle is returned (ownership transferred to the caller, which must
// destroy it) for in-place runtime updates of its values.
struct CeedErrorIntegratorInput
{
  CeedVector vec;
  CeedElemRestriction restr;
  CeedBasis basis;
};
void AssembleCeedElementErrorIntegrator(const CeedQFunctionInfo &info, void *ctx,
                                        std::size_t ctx_size, Ceed ceed,
                                        const std::vector<CeedErrorIntegratorInput> &inputs,
                                        CeedElemRestriction mesh_elem_restr,
                                        CeedVector geom_data,
                                        CeedElemRestriction geom_data_restr,
                                        CeedOperator *op, CeedQFunctionContext *ctx_out);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_INTEGRATOR_HPP
