// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_INTEGRATOR_HPP
#define PALACE_LIBCEED_INTEGRATOR_HPP

#include <string>
#include <vector>
#include <ceed.h>
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/utils.hpp"
#include "linalg/vector.hpp"

#include "fem/qfunctions/geom_qf.h"

// XX TODO WIP: FOR NOW, NO COEFFICIENTS IN ASSEMBLY

namespace palace::ceed
{

// Geometry factor information as quadrature data.
enum class GeomFactorInfo : unsigned int
{
  Determinant = 1 << 0,
  Jacobian = 1 << 1,
  Adjugate = 1 << 2,
  Weight = 1 << 3
};

// Evaluation modes for CeedOperator fields for various integrators.
enum class EvalMode : unsigned int
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
  CeedQFunctionUser build_qf, apply_qf;

  // Path and name of the QFunctions for operator construction and application.
  std::string build_qf_path, apply_qf_path;

  // Geometry factors required as QFunction inputs.
  unsigned int geom_data;

  // Evaluation modes for the test and trial basis.
  unsigned int trial_ops, test_ops;
};

// Create libCEED quadrature data and element restriction for use in a partially assembled
// libCEED operator.
inline void AssembleCeedGeometryData(GeomFactorInfo info, Ceed ceed,
                                     CeedElemRestriction mesh_restr, CeedBasis mesh_basis,
                                     const Vector &mesh_nodes, Vector &qdata,
                                     CeedVector *qdata_vec,
                                     CeedElemRestriction *qdata_restr)
{
  CeedInt ne, dim, space_dim, nqpts;
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(mesh_restr, &ne));
  PalaceCeedCall(ceed, CeedBasisGetDimension(mesh_basis, &dim));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(mesh_basis, &space_dim));
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &nqpts));

  // Create storage for quadrature point data.
  CeedInt qdata_size;
  switch (info)
  {
    case GeomFactorInfo::Determinant:
      qdata_size = 1;
      break;
    case GeomFactorInfo::Jacobian:
    case GeomFactorInfo::Adjugate:
      qdata_size = space_dim * dim;
      break;
    case GeomFactorInfo::Weight:
      MFEM_ABORT(
          "GeomFactorInfo::Weight is not a valid input for AssembleCeedGeometryData!");
      break;
  }
  qdata.SetSize(ne * nqpts * qdata_size);
  InitCeedVector(qdata, ceed, qdata_vec);
  PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(ceed, ne, nqpts, qdata_size,
                                                        ne * nqpts * qdata_size,
                                                        CEED_STRIDES_BACKEND, qdata_restr));

  // Create the QFunction that builds the operator (i.e. computes its quadrature data).
  CeedQFunction build_qf;
  CeedQFunctionUser build_qf_func;
  std::string build_qf_path;
  switch (info)
  {
    case GeomFactorInfo::Determinant:
      switch (10 * space_dim + dim)
      {
        case 22:
          build_qf_func = f_build_geom_factor_detJ22;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_detJ22_loc);
          break;
        case 21:
          build_qf_func = f_build_geom_factor_detJ21;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_detJ21_loc);
          break;
        case 33:
          build_qf_func = f_build_geom_factor_detJ33;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_detJ33_loc);
          break;
        case 32:
          build_qf_func = f_build_geom_factor_detJ32;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_detJ32_loc);
          break;
        default:
          MFEM_ABORT("Invalid value of (dim, space_dim) = (" << dim << ", " << space_dim
                                                             << ")!");
      }
      break;
    case GeomFactorInfo::Jacobian:
      switch (10 * space_dim + dim)
      {
        case 22:
          build_qf_func = f_build_geom_factor_J22;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_J22_loc);
          break;
        case 21:
          build_qf_func = f_build_geom_factor_J21;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_J21_loc);
          break;
        case 33:
          build_qf_func = f_build_geom_factor_J33;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_J33_loc);
          break;
        case 32:
          build_qf_func = f_build_geom_factor_J32;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_J32_loc);
          break;
        default:
          MFEM_ABORT("Invalid value of (dim, space_dim) = (" << dim << ", " << space_dim
                                                             << ")!");
      }
      break;
    case GeomFactorInfo::Adjugate:
      {
        case 22:
          build_qf_func = f_build_geom_factor_adjJt22;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_adjJt22_loc);
          break;
        case 21:
          build_qf_func = f_build_geom_factor_adjJt21;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_adjJt21_loc);
          break;
        case 33:
          build_qf_func = f_build_geom_factor_adjJt33;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_adjJt33_loc);
          break;
        case 32:
          build_qf_func = f_build_geom_factor_adjJt32;
          build_qf_path = PalaceQFunctionRelativePath(f_build_geom_factor_adjJt32_loc);
          break;
        default:
          MFEM_ABORT("Invalid value of (dim, space_dim) = (" << dim << ", " << space_dim
                                                             << ")!");
      }
      break;
  }
  PalaceCeedCall(ceed, CeedQFunctionCreateInterior(ceed, 1, build_qf_func,
                                                   build_qf_path.c_str(), &build_qf));

  // Inputs
  PalaceCeedCall(
      ceed, CeedQFunctionAddInput(build_qf, "grad_x", dim * space_dim, CEED_EVAL_GRAD));
  if (info == GeomFactorInfo::Determinant)
  {
    PalaceCeedCall(ceed, CeedQFunctionAddInput(build_qf, "w", 1, CEED_EVAL_WEIGHT));
  }

  // Output
  PalaceCeedCall(ceed,
                 CeedQFunctionAddOutput(build_qf, "qdata", qdata_size, CEED_EVAL_NONE));

  // Create the operator that builds the quadrature data for the actual operator.
  CeedOperator build_op;
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, build_qf, nullptr, nullptr, &build_op));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&build_qf));

  PalaceCeedCall(ceed, CeedOperatorSetField(build_op, "grad_x", mesh_restr, mesh_basis,
                                            CEED_VECTOR_ACTIVE));
  if (info == GeomFactorInfo::Determinant)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(build_op, "w", CEED_ELEMRESTRICTION_NONE,
                                              mesh_basis, CEED_VECTOR_NONE));
  }
  PalaceCeedCall(ceed, CeedOperatorSetField(build_op, "qdata", *qdata_restr,
                                            CEED_BASIS_NONE, CEED_VECTOR_ACTIVE));

  PalaceCeedCall(ceed, CeedOperatorCheckReady(build_op));

  // Compute the quadrature data for the operator.
  CeedVector nodes_vec;
  InitCeedVector(mesh_nodes, ceed, &nodes_vec);

  PalaceCeedCall(
      ceed, CeedOperatorApply(build_op, nodes_vec, *qdata_vec, CEED_REQUEST_IMMEDIATE));

  PalaceCeedCall(ceed, CeedVectorDestroy(&nodes_vec));
  PalaceCeedCall(ceed, CeedOperatorDestroy(&build_op));
}

// Create libCEED operator using the given quadrature data and element restriction.
template <typename CeedIntegratorInfo>
inline void AssembleCeedOperator(const CeedIntegratorInfo &info,
                                 const CeedGeomFactorData &geom_data, Ceed ceed,
                                 CeedElemRestriction trial_restr,
                                 CeedElemRestriction test_restr, CeedBasis trial_basis,
                                 CeedBasis test_basis, CeedOperator *op)
{
  // Create the QFunction that defines the action of the operator.
  CeedQFunction apply_qf;
  PalaceCeedCall(ceed, CeedQFunctionCreateInterior(ceed, 1, info.apply_qf,
                                                   info.apply_qf_path.c_str(), &apply_qf));

  CeedQFunctionContext apply_ctx;
  PalaceCeedCall(ceed, CeedQFunctionContextCreate(ceed, &apply_ctx));
  PalaceCeedCall(ceed,
                 CeedQFunctionContextSetData(apply_ctx, CEED_MEM_HOST, CEED_COPY_VALUES,
                                             sizeof(info.ctx), (void *)&info.ctx));
  PalaceCeedCall(ceed, CeedQFunctionSetContext(apply_qf, apply_ctx));
  PalaceCeedCall(ceed, CeedQFunctionContextDestroy(&apply_ctx));

  // Inputs
  CeedInt trial_ncomp, test_ncomp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_ncomp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_ncomp));

  if (info.geom_data & GeomFactorInfo::Determinant)
  {
    CeedInt qdata_size;
    PalaceCeedCall(ceed,
                   CeedElemRestrictionGetNumComponents(geom_data.wdetJ_restr, &qdata_size));
    PalaceCeedCall(ceed,
                   CeedQFunctionAddInput(apply_qf, "w_det_J", qdata_size, CEED_EVAL_NONE));
  }
  if (info.geom_data & GeomFactorInfo::Jacobian)
  {
    CeedInt qdata_size;
    PalaceCeedCall(ceed,
                   CeedElemRestrictionGetNumComponents(geom_data.J_restr, &qdata_size));
    PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "J", qdata_size, CEED_EVAL_NONE));
  }
  if (info.geom_data & GeomFactorInfo::Adjugate)
  {
    CeedInt qdata_size;
    PalaceCeedCall(ceed,
                   CeedElemRestrictionGetNumComponents(geom_data.adjJt_restr, &qdata_size));
    PalaceCeedCall(ceed,
                   CeedQFunctionAddInput(apply_qf, "adj_Jt", qdata_size, CEED_EVAL_NONE));
  }
  if (info.geom_data & GeomFactorInfo::Weight)
  {
    PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "w", 1, CEED_EVAL_WEIGHT));
  }

  if (info.trial_op & EvalMode::None):
    {
      PalaceCeedCall(ceed,
                     CeedQFunctionAddInput(apply_qf, "u", trial_ncomp, CEED_EVAL_NONE));
    }
  if (info.trial_op & EvalMode::Interp):
    {
      CeedInt qcomp;
      PalaceCeedCall(
          ceed, CeedBasisGetNumQuadratureComponents(trial_basis, CEED_EVAL_INTERP, &qcomp));
      PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "u", trial_ncomp * qcomp,
                                                 CEED_EVAL_INTERP));
    }
  if (info.trial_op & EvalMode::Grad):
    {
      CeedInt qcomp;
      PalaceCeedCall(
          ceed, CeedBasisGetNumQuadratureComponents(trial_basis, CEED_EVAL_GRAD, &qcomp));
      PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "grad_u", trial_ncomp * qcomp,
                                                 CEED_EVAL_GRAD));
    }
  if (info.trial_op & EvalMode::Div):
    {
      CeedInt qcomp;
      PalaceCeedCall(
          ceed, CeedBasisGetNumQuadratureComponents(trial_basis, CEED_EVAL_DIV, &qcomp));
      PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "div_u", trial_ncomp * qcomp,
                                                 CEED_EVAL_DIV));
    }
  if (info.trial_op & EvalMode::Curl):
    {
      CeedInt qcomp;
      PalaceCeedCall(
          ceed, CeedBasisGetNumQuadratureComponents(trial_basis, CEED_EVAL_CURL, &qcomp));
      PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "curl_u", trial_ncomp * qcomp,
                                                 CEED_EVAL_CURL));
    }

  // Output
  if (info.test_op & EvalMode::None):
    {
      PalaceCeedCall(ceed,
                     CeedQFunctionAddOutput(apply_qf, "v", test_ncomp, CEED_EVAL_NONE));
    }
  if (info.test_op & EvalMode::Interp):
    {
      CeedInt qcomp;
      PalaceCeedCall(
          ceed, CeedBasisGetNumQuadratureComponents(test_basis, CEED_EVAL_INTERP, &qcomp));
      PalaceCeedCall(ceed, CeedQFunctionAddOutput(apply_qf, "v", test_ncomp * qcomp,
                                                  CEED_EVAL_INTERP));
    }
  if (info.test_op & EvalMode::Grad):
    {
      CeedInt qcomp;
      PalaceCeedCall(
          ceed, CeedBasisGetNumQuadratureComponents(test_basis, CEED_EVAL_GRAD, &qcomp));
      PalaceCeedCall(ceed, CeedQFunctionAddOutput(apply_qf, "grad_v", test_ncomp * qcomp,
                                                  CEED_EVAL_GRAD));
    }
  if (info.test_op & EvalMode::Div):
    {
      CeedInt qcomp;
      PalaceCeedCall(
          ceed, CeedBasisGetNumQuadratureComponents(test_basis, CEED_EVAL_DIV, &qcomp));
      PalaceCeedCall(ceed, CeedQFunctionAddOutput(apply_qf, "div_v", test_ncomp * qcomp,
                                                  CEED_EVAL_DIV));
    }
  if (info.test_op & EvalMode::Curl):
    {
      CeedInt qcomp;
      PalaceCeedCall(
          ceed, CeedBasisGetNumQuadratureComponents(test_basis, CEED_EVAL_CURL, &qcomp));
      PalaceCeedCall(ceed, CeedQFunctionAddOutput(apply_qf, "curl_v", test_ncomp * qcomp,
                                                  CEED_EVAL_CURL));
    }

  // Create the operator.
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, apply_qf, nullptr, nullptr, op));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&apply_qf));

  if (info.geom_data & GeomFactorInfo::Determinant)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "w_det_J", geom_data.wdetJ_restr,
                                              CEED_BASIS_NONE, geom_data.wdetJ_vec));
  }
  if (info.geom_data & GeomFactorInfo::Jacobian)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "J", geom_data.J_restr, CEED_BASIS_NONE,
                                              geom_data.J_vec));
  }
  if (info.geom_data & GeomFactorInfo::Adjugate)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "adj_Jt", geom_data.adjJt_restr,
                                              CEED_BASIS_NONE, geom_data.adjJt_vec));
  }
  if (info.geom_data & GeomFactorInfo::Weight)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "w", CEED_ELEMRESTRICTION_NONE,
                                              trial_basis, CEED_VECTOR_NONE));
  }

  if (info.trial_op & EvalMode::None):
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(apply_qf, "u", trial_restr, CEED_BASIS_NONE,
                                                CEED_VECTOR_ACTIVE));
    }
  if (info.trial_op & EvalMode::Interp):
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(apply_qf, "u", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
    }
  if (info.trial_op & EvalMode::Grad):
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(apply_qf, "grad_u", trial_restr,
                                                trial_basis, CEED_VECTOR_ACTIVE));
    }
  if (info.trial_op & EvalMode::Div):
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(apply_qf, "div_u", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
    }
  if (info.trial_op & EvalMode::Curl):
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(apply_qf, "curl_u", trial_restr,
                                                trial_basis, CEED_VECTOR_ACTIVE));
    }

  if (info.test_op & EvalMode::None):
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(apply_qf, "v", test_restr, CEED_BASIS_NONE,
                                                CEED_VECTOR_ACTIVE));
    }
  if (info.test_op & EvalMode::Interp):
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(apply_qf, "v", test_restr, test_basis,
                                                CEED_VECTOR_ACTIVE));
    }
  if (info.test_op & EvalMode::Grad):
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(apply_qf, "grad_v", test_restr, test_basis,
                                                CEED_VECTOR_ACTIVE));
    }
  if (info.test_op & EvalMode::Div):
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(apply_qf, "div_v", test_restr, test_basis,
                                                CEED_VECTOR_ACTIVE));
    }
  if (info.test_op & EvalMode::Curl):
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(apply_qf, "curl_v", test_restr, test_basis,
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
