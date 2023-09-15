// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_INTEGRATOR_HPP
#define PALACE_LIBCEED_INTEGRATOR_HPP

#include <string>
#include <vector>
#include <ceed.h>
#include <mfem.hpp>
#include "basis.hpp"
#include "coefficient.hpp"
#include "restriction.hpp"
#include "utils.hpp"

namespace palace::ceed
{

// Evaluation modes for CeedOperator fields for various integrators.
enum class EvalMode
{
  None,
  Interp,
  Grad,
  Div,
  Curl,
  InterpAndGrad,
  InterpAndDiv,
  InterpAndCurl
};

// Data structure for CeedOperator construction for various integrators.
struct IntegratorInfo
{
  // QFunctions for operator construction and application.
  CeedQFunctionUser build_qf, apply_qf;

  // Path and name of the QFunctions for operator construction and application.
  std::string build_qf_path, apply_qf_path;

  // Evaluation modes for the test and trial basis.
  EvalMode trial_op, test_op;

  // Size of the data at each quadrature point.
  int qdata_size;
};

// Helper function which combines quadrature data assembly and operator assembly in a single
// method.
template <typename CeedIntegratorInfo>
inline void AssembleCeedOperator(const CeedIntegratorInfo &info,
                                 const mfem::FiniteElementSpace &trial_fespace,
                                 const mfem::FiniteElementSpace &test_fespace,
                                 const mfem::IntegrationRule &ir,
                                 const std::vector<int> &indices, const bool use_bdr,
                                 const std::vector<QuadratureCoefficient> &Q, Ceed ceed,
                                 CeedOperator *op, CeedOperator *op_t)
{
  // Assemble quadrature data.
  CeedVector qdata;
  CeedElemRestriction qdata_restr;
  AssembleCeedQuadratureData(info, trial_fespace, test_fespace, ir, indices, use_bdr, Q,
                             ceed, &qdata, &qdata_restr);

  // Assemble the operator (no transpose).
  AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, qdata,
                       qdata_restr, ceed, op);
  *op_t = nullptr;

  // Cleanup (these are now owned by the operator).
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&qdata_restr));
  PalaceCeedCall(ceed, CeedVectorDestroy(&qdata));
}

// Create libCEED quadrature data and element restriction for use in a partially assembled
// libCEED operator.
template <typename CeedIntegratorInfo>
inline void
AssembleCeedQuadratureData(const CeedIntegratorInfo &info,
                           const mfem::FiniteElementSpace &trial_fespace,
                           const mfem::FiniteElementSpace &test_fespace,
                           const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                           const bool use_bdr, const std::vector<QuadratureCoefficient> &Q,
                           Ceed ceed, CeedVector *qdata, CeedElemRestriction *qdata_restr)
{
  MFEM_VERIFY(trial_fespace.GetMesh() == test_fespace.GetMesh(),
              "Trial and test finite element spaces must correspond to the same mesh!");
  mfem::Mesh &mesh = *trial_fespace.GetMesh();
  MFEM_VERIFY(mesh.GetNodes(), "The mesh has no nodal FE space!");
  const mfem::FiniteElementSpace &mesh_fespace = *mesh.GetNodalFESpace();

  CeedInt ne = static_cast<CeedInt>(indices.size());
  CeedInt dim = mesh.Dimension() - use_bdr;
  CeedInt space_dim = mesh.SpaceDimension();

  CeedElemRestriction mesh_restr;
  CeedBasis mesh_basis;
  CeedInt nqpts, qdata_size = info.qdata_size;
  InitRestriction(mesh_fespace, indices, use_bdr, ceed, &mesh_restr);
  InitBasis(mesh_fespace, ir, indices, use_bdr, ceed, &mesh_basis);
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &nqpts));

  // Strided restrictions are cheap to construct and not stored in the global cache.
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, ne * nqpts * qdata_size, qdata));
  PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(ceed, ne, nqpts, qdata_size,
                                                        ne * nqpts * qdata_size,
                                                        CEED_STRIDES_BACKEND, qdata_restr));

  // Create the QFunction that builds the operator (i.e. computes its quadrature data).
  CeedQFunction build_qf;
  PalaceCeedCall(ceed, CeedQFunctionCreateInterior(ceed, 1, info.build_qf,
                                                   info.build_qf_path.c_str(), &build_qf));

  CeedQFunctionContext build_ctx;
  PalaceCeedCall(ceed, CeedQFunctionContextCreate(ceed, &build_ctx));
  PalaceCeedCall(ceed,
                 CeedQFunctionContextSetData(build_ctx, CEED_MEM_HOST, CEED_COPY_VALUES,
                                             sizeof(info.ctx), (void *)&info.ctx));
  PalaceCeedCall(ceed, CeedQFunctionSetContext(build_qf, build_ctx));
  PalaceCeedCall(ceed, CeedQFunctionContextDestroy(&build_ctx));

  // Inputs
  for (std::size_t i = 0; i < Q.size(); i++)
  {
    std::string name = "coeff" + std::to_string(i + 1);
    const CeedInt ncomp = Q[i].ncomp;
    PalaceCeedCall(ceed,
                   CeedQFunctionAddInput(build_qf, name.c_str(), ncomp, CEED_EVAL_NONE));
  }
  PalaceCeedCall(ceed,
                 CeedQFunctionAddInput(build_qf, "dx", dim * space_dim, CEED_EVAL_GRAD));
  PalaceCeedCall(ceed, CeedQFunctionAddInput(build_qf, "weights", 1, CEED_EVAL_WEIGHT));

  // Output
  PalaceCeedCall(ceed,
                 CeedQFunctionAddOutput(build_qf, "qdata", qdata_size, CEED_EVAL_NONE));

  // Create the operator that builds the quadrature data for the actual operator.
  CeedOperator build_op;
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, build_qf, nullptr, nullptr, &build_op));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&build_qf));

  for (std::size_t i = 0; i < Q.size(); i++)
  {
    std::string name = "coeff" + std::to_string(i + 1);
    const CeedInt ncomp = Q[i].ncomp;
    CeedInt strides[3] = {ncomp, 1, ncomp * nqpts};
    CeedElemRestriction coeff_restr;
    CeedVector coeff_vector;

    PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(ceed, ne, nqpts, ncomp,
                                                          ne * nqpts * ncomp, strides,
                                                          &coeff_restr));
    InitCeedVector(Q[i].data, ceed, &coeff_vector);

    PalaceCeedCall(ceed, CeedOperatorSetField(build_op, name.c_str(), coeff_restr,
                                              CEED_BASIS_NONE, coeff_vector));

    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&coeff_restr));
    PalaceCeedCall(ceed, CeedVectorDestroy(&coeff_vector));
  }
  PalaceCeedCall(ceed, CeedOperatorSetField(build_op, "dx", mesh_restr, mesh_basis,
                                            CEED_VECTOR_ACTIVE));
  PalaceCeedCall(ceed, CeedOperatorSetField(build_op, "weights", CEED_ELEMRESTRICTION_NONE,
                                            mesh_basis, CEED_VECTOR_NONE));
  PalaceCeedCall(ceed, CeedOperatorSetField(build_op, "qdata", *qdata_restr,
                                            CEED_BASIS_NONE, CEED_VECTOR_ACTIVE));

  PalaceCeedCall(ceed, CeedOperatorCheckReady(build_op));

  // Compute the quadrature data for the operator.
  CeedVector nodes;
  InitCeedVector(*mesh.GetNodes(), ceed, &nodes);

  PalaceCeedCall(ceed, CeedOperatorApply(build_op, nodes, *qdata, CEED_REQUEST_IMMEDIATE));

  PalaceCeedCall(ceed, CeedVectorDestroy(&nodes));
  PalaceCeedCall(ceed, CeedOperatorDestroy(&build_op));
}

// Create libCEED operator using the given quadrature data and element restriction.
template <typename CeedIntegratorInfo>
inline void AssembleCeedOperator(const CeedIntegratorInfo &info,
                                 const mfem::FiniteElementSpace &trial_fespace,
                                 const mfem::FiniteElementSpace &test_fespace,
                                 const mfem::IntegrationRule &ir,
                                 const std::vector<int> &indices, const bool use_bdr,
                                 CeedVector qdata, CeedElemRestriction qdata_restr,
                                 Ceed ceed, CeedOperator *op)
{
  MFEM_VERIFY(trial_fespace.GetMesh() == test_fespace.GetMesh(),
              "Trial and test finite element spaces must correspond to the same mesh!");
  mfem::Mesh &mesh = *trial_fespace.GetMesh();

  CeedInt dim = mesh.Dimension() - use_bdr;
  CeedInt curl_dim = (dim < 3) ? 1 : dim;
  CeedInt trial_vdim = trial_fespace.GetVDim();
  CeedInt test_vdim = test_fespace.GetVDim();
  bool trial_vectorfe =
      (trial_fespace.FEColl()->GetRangeType(dim) == mfem::FiniteElement::VECTOR);
  bool test_vectorfe =
      (test_fespace.FEColl()->GetRangeType(dim) == mfem::FiniteElement::VECTOR);

  CeedElemRestriction trial_restr, test_restr;
  CeedBasis trial_basis, test_basis;
  InitRestriction(trial_fespace, indices, use_bdr, ceed, &trial_restr);
  InitRestriction(test_fespace, indices, use_bdr, ceed, &test_restr);
  InitBasis(trial_fespace, ir, indices, use_bdr, ceed, &trial_basis);
  InitBasis(test_fespace, ir, indices, use_bdr, ceed, &test_basis);

  CeedInt trial_nqpts, test_nqpts, mesh_nqpts, qdata_size;
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(trial_basis, &trial_nqpts));
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(test_basis, &test_nqpts));
  PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(qdata_restr, &mesh_nqpts));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(qdata_restr, &qdata_size));
  MFEM_VERIFY(trial_nqpts == test_nqpts && trial_nqpts == mesh_nqpts,
              "Trial and test basis must have the same number of quadrature points!");

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
  switch (info.trial_op)
  {
    case EvalMode::None:
      PalaceCeedCall(ceed,
                     CeedQFunctionAddInput(apply_qf, "u", trial_vdim, CEED_EVAL_NONE));
      break;
    case EvalMode::Interp:
      PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "u",
                                                 trial_vdim * (trial_vectorfe ? dim : 1),
                                                 CEED_EVAL_INTERP));
      break;
    case EvalMode::Grad:
      MFEM_VERIFY(!trial_vectorfe, "EvalMode::Grad is not intended for vector FE!");
      PalaceCeedCall(
          ceed, CeedQFunctionAddInput(apply_qf, "gu", trial_vdim * dim, CEED_EVAL_GRAD));
      break;
    case EvalMode::Div:
      PalaceCeedCall(ceed,
                     CeedQFunctionAddInput(apply_qf, "du", trial_vdim, CEED_EVAL_DIV));
      break;
    case EvalMode::Curl:
      PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "cu", trial_vdim * curl_dim,
                                                 CEED_EVAL_CURL));
      break;
    case EvalMode::InterpAndGrad:
      MFEM_VERIFY(!trial_vectorfe,
                  "EvalMode::InterpAndGrad is not intended for vector FE!");
      PalaceCeedCall(ceed,
                     CeedQFunctionAddInput(apply_qf, "u", trial_vdim, CEED_EVAL_INTERP));
      PalaceCeedCall(
          ceed, CeedQFunctionAddInput(apply_qf, "gu", trial_vdim * dim, CEED_EVAL_GRAD));
      break;
    case EvalMode::InterpAndDiv:
      MFEM_VERIFY(trial_vectorfe, "EvalMode::InterpAndDiv is only intended for vector FE!");
      PalaceCeedCall(
          ceed, CeedQFunctionAddInput(apply_qf, "u", trial_vdim * dim, CEED_EVAL_INTERP));
      PalaceCeedCall(ceed,
                     CeedQFunctionAddInput(apply_qf, "du", trial_vdim, CEED_EVAL_DIV));
      break;
    case EvalMode::InterpAndCurl:
      MFEM_VERIFY(trial_vectorfe,
                  "EvalMode::InterpAndCurl is only intended for vector FE!");
      PalaceCeedCall(
          ceed, CeedQFunctionAddInput(apply_qf, "u", trial_vdim * dim, CEED_EVAL_INTERP));
      PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "cu", trial_vdim * curl_dim,
                                                 CEED_EVAL_CURL));
      break;
  }
  PalaceCeedCall(ceed,
                 CeedQFunctionAddInput(apply_qf, "qdata", qdata_size, CEED_EVAL_NONE));

  // Output
  switch (info.test_op)
  {
    case EvalMode::None:
      PalaceCeedCall(ceed,
                     CeedQFunctionAddOutput(apply_qf, "v", test_vdim, CEED_EVAL_NONE));
      break;
    case EvalMode::Interp:
      PalaceCeedCall(ceed, CeedQFunctionAddOutput(apply_qf, "v",
                                                  test_vdim * (test_vectorfe ? dim : 1),
                                                  CEED_EVAL_INTERP));
      break;
    case EvalMode::Grad:
      MFEM_VERIFY(!test_vectorfe, "EvalMode::Grad is not intended for vector FE!");
      PalaceCeedCall(
          ceed, CeedQFunctionAddOutput(apply_qf, "gv", test_vdim * dim, CEED_EVAL_GRAD));
      break;
    case EvalMode::Div:
      PalaceCeedCall(ceed,
                     CeedQFunctionAddOutput(apply_qf, "dv", test_vdim, CEED_EVAL_DIV));
      break;
    case EvalMode::Curl:
      PalaceCeedCall(ceed, CeedQFunctionAddOutput(apply_qf, "cv", test_vdim * curl_dim,
                                                  CEED_EVAL_CURL));
      break;
    case EvalMode::InterpAndGrad:
      MFEM_VERIFY(!test_vectorfe, "EvalMode::InterpAndGrad is not intended for vector FE!");
      PalaceCeedCall(ceed,
                     CeedQFunctionAddOutput(apply_qf, "v", test_vdim, CEED_EVAL_INTERP));
      PalaceCeedCall(
          ceed, CeedQFunctionAddOutput(apply_qf, "gv", test_vdim * dim, CEED_EVAL_GRAD));
      break;
    case EvalMode::InterpAndDiv:
      MFEM_VERIFY(test_vectorfe, "EvalMode::InterpAndDiv is only intended for vector FE!");
      PalaceCeedCall(
          ceed, CeedQFunctionAddOutput(apply_qf, "v", test_vdim * dim, CEED_EVAL_INTERP));
      PalaceCeedCall(ceed,
                     CeedQFunctionAddOutput(apply_qf, "dv", test_vdim, CEED_EVAL_DIV));
      break;
    case EvalMode::InterpAndCurl:
      MFEM_VERIFY(test_vectorfe, "EvalMode::InterpAndCurl is only intended for vector FE!");
      PalaceCeedCall(
          ceed, CeedQFunctionAddOutput(apply_qf, "v", test_vdim * dim, CEED_EVAL_INTERP));
      PalaceCeedCall(ceed, CeedQFunctionAddOutput(apply_qf, "cv", test_vdim * curl_dim,
                                                  CEED_EVAL_CURL));
      break;
  }

  // Create the operator.
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, apply_qf, nullptr, nullptr, op));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&apply_qf));

  switch (info.trial_op)
  {
    case EvalMode::None:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "u", trial_restr, CEED_BASIS_NONE,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::Interp:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "u", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::Grad:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "gu", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::Div:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "du", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::Curl:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "cu", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::InterpAndGrad:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "u", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "gu", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::InterpAndDiv:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "u", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "du", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::InterpAndCurl:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "u", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "cu", trial_restr, trial_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
  }
  PalaceCeedCall(ceed,
                 CeedOperatorSetField(*op, "qdata", qdata_restr, CEED_BASIS_NONE, qdata));
  switch (info.test_op)
  {
    case EvalMode::None:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "v", test_restr, CEED_BASIS_NONE,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::Interp:
      PalaceCeedCall(
          ceed, CeedOperatorSetField(*op, "v", test_restr, test_basis, CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::Grad:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "gv", test_restr, test_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::Div:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "dv", test_restr, test_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::Curl:
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "cv", test_restr, test_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::InterpAndGrad:
      PalaceCeedCall(
          ceed, CeedOperatorSetField(*op, "v", test_restr, test_basis, CEED_VECTOR_ACTIVE));
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "gv", test_restr, test_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::InterpAndDiv:
      PalaceCeedCall(
          ceed, CeedOperatorSetField(*op, "v", test_restr, test_basis, CEED_VECTOR_ACTIVE));
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "dv", test_restr, test_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
    case EvalMode::InterpAndCurl:
      PalaceCeedCall(
          ceed, CeedOperatorSetField(*op, "v", test_restr, test_basis, CEED_VECTOR_ACTIVE));
      PalaceCeedCall(ceed, CeedOperatorSetField(*op, "cv", test_restr, test_basis,
                                                CEED_VECTOR_ACTIVE));
      break;
  }

  PalaceCeedCall(ceed, CeedOperatorCheckReady(*op));
}

// Construct libCEED operators for interpolation operations and their transpose between
// the two spaces. The operation for interpolation is decided by the conformity of the trial
// and test spaces.
inline void AssembleCeedInterpolator(const mfem::FiniteElementSpace &trial_fespace,
                                     const mfem::FiniteElementSpace &test_fespace,
                                     const std::vector<int> &indices, Ceed ceed,
                                     CeedOperator *op, CeedOperator *op_t)
{
  CeedInt trial_vdim = trial_fespace.GetVDim();
  CeedInt test_vdim = test_fespace.GetVDim();
  MFEM_VERIFY(trial_vdim == 1 && test_vdim == 1,
              "AssembleCeedInterpolator does not support spaces with vdim > 1!");

  CeedElemRestriction trial_restr, test_restr;
  CeedBasis basis_ctof;
  InitRestriction(trial_fespace, indices, false, true, false, ceed, &trial_restr);
  InitRestriction(test_fespace, indices, false, true, true, ceed, &test_restr);
  InitInterpolatorBasis(trial_fespace, test_fespace, indices, ceed, &basis_ctof);

  // Create the QFunction that defines the action of the operator (only an identity as
  // element dof multiplicity is handled outside of libCEED).
  CeedQFunction apply_qf, apply_qf_t;
  PalaceCeedCall(ceed, CeedQFunctionCreateIdentity(ceed, trial_vdim, CEED_EVAL_INTERP,
                                                   CEED_EVAL_NONE, &apply_qf));
  PalaceCeedCall(ceed, CeedQFunctionCreateIdentity(ceed, trial_vdim, CEED_EVAL_NONE,
                                                   CEED_EVAL_INTERP, &apply_qf_t));

  // Create the operator.
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, apply_qf, nullptr, nullptr, op));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&apply_qf));

  PalaceCeedCall(ceed, CeedOperatorSetField(*op, "input", trial_restr, basis_ctof,
                                            CEED_VECTOR_ACTIVE));
  PalaceCeedCall(ceed, CeedOperatorSetField(*op, "output", test_restr, CEED_BASIS_NONE,
                                            CEED_VECTOR_ACTIVE));

  PalaceCeedCall(ceed, CeedOperatorCheckReady(*op));

  // Create the transpose operator.
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, apply_qf_t, nullptr, nullptr, op_t));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&apply_qf_t));

  PalaceCeedCall(ceed, CeedOperatorSetField(*op_t, "input", test_restr, CEED_BASIS_NONE,
                                            CEED_VECTOR_ACTIVE));
  PalaceCeedCall(ceed, CeedOperatorSetField(*op_t, "output", trial_restr, basis_ctof,
                                            CEED_VECTOR_ACTIVE));

  PalaceCeedCall(ceed, CeedOperatorCheckReady(*op_t));
}

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_INTEGRATOR_HPP
