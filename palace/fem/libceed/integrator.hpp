// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_INTEGRATOR_HPP
#define PALACE_LIBCEED_INTEGRATOR_HPP

#include <string>
#include <ceed/backend.h>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"
#include "linalg/vector.hpp"

#include "fem/qfunctions/geom_qf.h"

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
};

// Create libCEED quadrature data and element restriction for use in a partially assembled
// libCEED operator.
inline void AssembleCeedGeometryData(Ceed ceed, CeedElemRestriction mesh_restr,
                                     CeedBasis mesh_basis, const Vector &mesh_nodes,
                                     CeedGeomFactorData &geom_data)
{
  CeedInt ne, dim, space_dim, nqpts;
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(mesh_restr, &ne));
  PalaceCeedCall(ceed, CeedBasisGetDimension(mesh_basis, &dim));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(mesh_basis, &space_dim));
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &nqpts));

  // Create the QFunction that builds the operator (i.e. computes its quadrature data).
  CeedQFunction build_qf;
  switch (10 * space_dim + dim)
  {
    case 22:
      PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                               ceed, 1, f_build_geom_factor_22,
                               PalaceQFunctionRelativePath(f_build_geom_factor_22_loc),
                               &build_qf));
      break;
    case 33:
      PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                               ceed, 1, f_build_geom_factor_33,
                               PalaceQFunctionRelativePath(f_build_geom_factor_33_loc),
                               &build_qf));
      break;
    case 21:
      PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                               ceed, 1, f_build_geom_factor_21,
                               PalaceQFunctionRelativePath(f_build_geom_factor_21_loc),
                               &build_qf));
      break;
    case 32:
      PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                               ceed, 1, f_build_geom_factor_32,
                               PalaceQFunctionRelativePath(f_build_geom_factor_32_loc),
                               &build_qf));
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << dim << ", " << space_dim << ") for geometry factor quadrature data!");
      build_qf = nullptr;  // Silence compiler warning
  }

  // Inputs
  PalaceCeedCall(
      ceed, CeedQFunctionAddInput(build_qf, "grad_x", space_dim * dim, CEED_EVAL_GRAD));
  PalaceCeedCall(ceed, CeedQFunctionAddInput(build_qf, "w", 1, CEED_EVAL_WEIGHT));

  // Outputs
  {
    CeedInt qdata_size = 1;
    geom_data->wdetJ.SetSize(ne * nqpts * qdata_size);
    InitCeedVector(geom_data->wdetJ, ceed, &geom_data->wdetJ_vec);
    PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(
                             ceed, ne, nqpts, qdata_size, ne * nqpts * qdata_size,
                             CEED_STRIDES_BACKEND, &geom_data->wdetJ_restr));
    PalaceCeedCall(ceed,
                   CeedQFunctionAddOutput(build_qf, "w_det_J", qdata_size, CEED_EVAL_NONE));
  }
  {
    CeedInt qdata_size = space_dim * dim;
    geom_data->adjJt.SetSize(ne * nqpts * qdata_size);
    InitCeedVector(geom_data->adjJt, ceed, &geom_data->adjJt_vec);
    PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(
                             ceed, ne, nqpts, qdata_size, ne * nqpts * qdata_size,
                             CEED_STRIDES_BACKEND, &geom_data->adjJt_restr));
    PalaceCeedCall(ceed,
                   CeedQFunctionAddOutput(build_qf, "adj_Jt", qdata_size, CEED_EVAL_NONE));
  }

  // Create the operator that builds the quadrature data.
  CeedOperator build_op;
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, build_qf, nullptr, nullptr, &build_op));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&build_qf));

  PalaceCeedCall(ceed, CeedOperatorSetField(build_op, "grad_x", mesh_restr, mesh_basis,
                                            CEED_VECTOR_ACTIVE));
  PalaceCeedCall(ceed, CeedOperatorSetField(build_op, "w", CEED_ELEMRESTRICTION_NONE,
                                            mesh_basis, CEED_VECTOR_NONE));

  {
    PalaceCeedCall(ceed, CeedOperatorSetField(build_op, "w_det_J", geom_data->wdetJ_restr,
                                              CEED_BASIS_NONE, geom_data->wdetJ_vec));
  }
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(build_op, "adj_Jt", geom_data->adjJt_restr,
                                              CEED_BASIS_NONE, geom_data->adjJt_vec));
  }

  PalaceCeedCall(ceed, CeedOperatorCheckReady(build_op));

  // Compute the quadrature data for the operator (all outputs are passive).
  CeedVector nodes_vec;
  InitCeedVector(mesh_nodes, ceed, &nodes_vec);

  PalaceCeedCall(ceed, CeedOperatorApply(build_op, nodes_vec, CEED_VECTOR_NONE,
                                         CEED_REQUEST_IMMEDIATE));

  PalaceCeedCall(ceed, CeedVectorDestroy(&nodes_vec));
  PalaceCeedCall(ceed, CeedOperatorDestroy(&build_op));
}

// Create libCEED operator using the given quadrature data and element restriction.
template <typename IntegratorContext>
inline void AssembleCeedOperator(const IntegratorInfo &info, const IntegratorContext &ctx,
                                 const CeedGeomFactorData &geom_data, Ceed ceed,
                                 CeedElemRestriction trial_restr,
                                 CeedElemRestriction test_restr, CeedBasis trial_basis,
                                 CeedBasis test_basis, CeedOperator *op)
{
  // XX TODO: Add quadrature data assembly option which computes once all the quadrature
  //          data and creates the "simple" operator to just multiply

  // Create the QFunction that defines the action of the operator.
  CeedQFunction apply_qf;
  PalaceCeedCall(ceed, CeedQFunctionCreateInterior(ceed, 1, info.apply_qf,
                                                   info.apply_qf_path.c_str(), &apply_qf));

  CeedQFunctionContext apply_ctx;
  PalaceCeedCall(ceed, CeedQFunctionContextCreate(ceed, &apply_ctx));
  PalaceCeedCall(ceed,
                 CeedQFunctionContextSetData(apply_ctx, CEED_MEM_HOST, CEED_COPY_VALUES,
                                             sizeof(ctx), (void *)&ctx));
  PalaceCeedCall(ceed, CeedQFunctionSetContext(apply_qf, apply_ctx));
  PalaceCeedCall(ceed, CeedQFunctionContextDestroy(&apply_ctx));

  // Inputs
  CeedInt trial_ncomp, test_ncomp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_ncomp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_ncomp));

  if (info.geom_info & GeomFactorInfo::Determinant)
  {
    CeedInt qdata_size;
    PalaceCeedCall(
        ceed, CeedElemRestrictionGetNumComponents(geom_data->wdetJ_restr, &qdata_size));
    PalaceCeedCall(ceed,
                   CeedQFunctionAddInput(apply_qf, "w_det_J", qdata_size, CEED_EVAL_NONE));
  }
  if (info.geom_info & GeomFactorInfo::Adjugate)
  {
    CeedInt qdata_size;
    PalaceCeedCall(
        ceed, CeedElemRestrictionGetNumComponents(geom_data->adjJt_restr, &qdata_size));
    PalaceCeedCall(ceed,
                   CeedQFunctionAddInput(apply_qf, "adj_Jt", qdata_size, CEED_EVAL_NONE));
  }
  if (info.geom_info & GeomFactorInfo::Weight)
  {
    PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "w", 1, CEED_EVAL_WEIGHT));
  }

  // Last non-active input is always the element attribute data.
  PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "attr", 1, CEED_EVAL_INTERP));

  if (info.trial_ops & EvalMode::None)
  {
    PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "u", trial_ncomp, CEED_EVAL_NONE));
  }
  if (info.trial_ops & EvalMode::Interp)
  {
    CeedInt qcomp;
    PalaceCeedCall(
        ceed, CeedBasisGetNumQuadratureComponents(trial_basis, CEED_EVAL_INTERP, &qcomp));
    PalaceCeedCall(
        ceed, CeedQFunctionAddInput(apply_qf, "u", trial_ncomp * qcomp, CEED_EVAL_INTERP));
  }
  if (info.trial_ops & EvalMode::Grad)
  {
    CeedInt qcomp;
    PalaceCeedCall(
        ceed, CeedBasisGetNumQuadratureComponents(trial_basis, CEED_EVAL_GRAD, &qcomp));
    PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "grad_u", trial_ncomp * qcomp,
                                               CEED_EVAL_GRAD));
  }
  if (info.trial_ops & EvalMode::Div)
  {
    CeedInt qcomp;
    PalaceCeedCall(ceed,
                   CeedBasisGetNumQuadratureComponents(trial_basis, CEED_EVAL_DIV, &qcomp));
    PalaceCeedCall(
        ceed, CeedQFunctionAddInput(apply_qf, "div_u", trial_ncomp * qcomp, CEED_EVAL_DIV));
  }
  if (info.trial_ops & EvalMode::Curl)
  {
    CeedInt qcomp;
    PalaceCeedCall(
        ceed, CeedBasisGetNumQuadratureComponents(trial_basis, CEED_EVAL_CURL, &qcomp));
    PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "curl_u", trial_ncomp * qcomp,
                                               CEED_EVAL_CURL));
  }

  // Outputs
  if (info.test_ops & EvalMode::None)
  {
    PalaceCeedCall(ceed, CeedQFunctionAddOutput(apply_qf, "v", test_ncomp, CEED_EVAL_NONE));
  }
  if (info.test_ops & EvalMode::Interp)
  {
    CeedInt qcomp;
    PalaceCeedCall(
        ceed, CeedBasisGetNumQuadratureComponents(test_basis, CEED_EVAL_INTERP, &qcomp));
    PalaceCeedCall(
        ceed, CeedQFunctionAddOutput(apply_qf, "v", test_ncomp * qcomp, CEED_EVAL_INTERP));
  }
  if (info.test_ops & EvalMode::Grad)
  {
    CeedInt qcomp;
    PalaceCeedCall(ceed,
                   CeedBasisGetNumQuadratureComponents(test_basis, CEED_EVAL_GRAD, &qcomp));
    PalaceCeedCall(ceed, CeedQFunctionAddOutput(apply_qf, "grad_v", test_ncomp * qcomp,
                                                CEED_EVAL_GRAD));
  }
  if (info.test_ops & EvalMode::Div)
  {
    CeedInt qcomp;
    PalaceCeedCall(ceed,
                   CeedBasisGetNumQuadratureComponents(test_basis, CEED_EVAL_DIV, &qcomp));
    PalaceCeedCall(
        ceed, CeedQFunctionAddOutput(apply_qf, "div_v", test_ncomp * qcomp, CEED_EVAL_DIV));
  }
  if (info.test_ops & EvalMode::Curl)
  {
    CeedInt qcomp;
    PalaceCeedCall(ceed,
                   CeedBasisGetNumQuadratureComponents(test_basis, CEED_EVAL_CURL, &qcomp));
    PalaceCeedCall(ceed, CeedQFunctionAddOutput(apply_qf, "curl_v", test_ncomp * qcomp,
                                                CEED_EVAL_CURL));
  }

  // Create the operator.
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, apply_qf, nullptr, nullptr, op));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&apply_qf));

  if (info.geom_info & GeomFactorInfo::Determinant)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "w_det_J", geom_data->wdetJ_restr,
                                              CEED_BASIS_NONE, geom_data->wdetJ_vec));
  }
  if (info.geom_info & GeomFactorInfo::Adjugate)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "adj_Jt", geom_data->adjJt_restr,
                                              CEED_BASIS_NONE, geom_data->adjJt_vec));
  }
  if (info.geom_info & GeomFactorInfo::Weight)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "w", CEED_ELEMRESTRICTION_NONE,
                                              trial_basis, CEED_VECTOR_NONE));
  }

  PalaceCeedCall(ceed, CeedOperatorSetField(*op, "attr", geom_data->attr_restr,
                                            geom_data->attr_basis, geom_data->attr_vec));

  if (info.trial_ops & EvalMode::None)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "u", trial_restr, CEED_BASIS_NONE,
                                              CEED_VECTOR_ACTIVE));
  }
  if (info.trial_ops & EvalMode::Interp)
  {
    PalaceCeedCall(
        ceed, CeedOperatorSetField(*op, "u", trial_restr, trial_basis, CEED_VECTOR_ACTIVE));
  }
  if (info.trial_ops & EvalMode::Grad)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "grad_u", trial_restr, trial_basis,
                                              CEED_VECTOR_ACTIVE));
  }
  if (info.trial_ops & EvalMode::Div)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "div_u", trial_restr, trial_basis,
                                              CEED_VECTOR_ACTIVE));
  }
  if (info.trial_ops & EvalMode::Curl)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "curl_u", trial_restr, trial_basis,
                                              CEED_VECTOR_ACTIVE));
  }

  if (info.test_ops & EvalMode::None)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "v", test_restr, CEED_BASIS_NONE,
                                              CEED_VECTOR_ACTIVE));
  }
  if (info.test_ops & EvalMode::Interp)
  {
    PalaceCeedCall(
        ceed, CeedOperatorSetField(*op, "v", test_restr, test_basis, CEED_VECTOR_ACTIVE));
  }
  if (info.test_ops & EvalMode::Grad)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "grad_v", test_restr, test_basis,
                                              CEED_VECTOR_ACTIVE));
  }
  if (info.test_ops & EvalMode::Div)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "div_v", test_restr, test_basis,
                                              CEED_VECTOR_ACTIVE));
  }
  if (info.test_ops & EvalMode::Curl)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "curl_v", test_restr, test_basis,
                                              CEED_VECTOR_ACTIVE));
  }

  PalaceCeedCall(ceed, CeedOperatorCheckReady(*op));
}

// Construct libCEED operators for interpolation operations and their transpose between
// the two spaces. The operation for interpolation is decided by the conformity of the trial
// and test spaces.
inline void AssembleCeedInterpolator(Ceed ceed, CeedElemRestriction trial_restr,
                                     CeedElemRestriction test_restr, CeedBasis interp_basis,
                                     CeedOperator *op, CeedOperator *op_t)
{
  // Create the QFunction that defines the action of the operator (only an identity as
  // element dof multiplicity is handled outside of libCEED).
  CeedQFunction apply_qf, apply_qf_t;
  PalaceCeedCall(ceed, CeedQFunctionCreateIdentity(ceed, 1, CEED_EVAL_INTERP,
                                                   CEED_EVAL_NONE, &apply_qf));
  PalaceCeedCall(ceed, CeedQFunctionCreateIdentity(ceed, 1, CEED_EVAL_NONE,
                                                   CEED_EVAL_INTERP, &apply_qf_t));

  // Create the operator.
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, apply_qf, nullptr, nullptr, op));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&apply_qf));

  PalaceCeedCall(ceed, CeedOperatorSetField(*op, "input", trial_restr, interp_basis,
                                            CEED_VECTOR_ACTIVE));
  PalaceCeedCall(ceed, CeedOperatorSetField(*op, "output", test_restr, CEED_BASIS_NONE,
                                            CEED_VECTOR_ACTIVE));

  PalaceCeedCall(ceed, CeedOperatorCheckReady(*op));

  // Create the transpose operator.
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, apply_qf_t, nullptr, nullptr, op_t));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&apply_qf_t));

  PalaceCeedCall(ceed, CeedOperatorSetField(*op_t, "input", test_restr, CEED_BASIS_NONE,
                                            CEED_VECTOR_ACTIVE));
  PalaceCeedCall(ceed, CeedOperatorSetField(*op_t, "output", trial_restr, interp_basis,
                                            CEED_VECTOR_ACTIVE));

  PalaceCeedCall(ceed, CeedOperatorCheckReady(*op_t));
}

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_INTEGRATOR_HPP
