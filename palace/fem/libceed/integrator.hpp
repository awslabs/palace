// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_INTEGRATOR_HPP
#define PALACE_LIBCEED_INTEGRATOR_HPP

#include <string>
#include <vector>
#include <ceed/backend.h>
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

// Construct libCEED operators for interpolation operations and their transpose between
// the two spaces. The operation for interpolation is decided by the conformity of the trial
// and test spaces.
void AssembleCeedInterpolator(Ceed ceed, CeedElemRestriction trial_restr,
                              CeedElemRestriction test_restr, CeedBasis interp_basis,
                              CeedOperator *op, CeedOperator *op_t);

namespace internal
{

void AddQFunctionActiveInputsOutputs(const IntegratorInfo &info, Ceed ceed,
                                     CeedBasis trial_basis, CeedBasis test_basis,
                                     CeedQFunction qf);

void AddOperatorActiveFields(const IntegratorInfo &info, Ceed ceed,
                             CeedElemRestriction trial_restr,
                             CeedElemRestriction test_restr, CeedBasis trial_basis,
                             CeedBasis test_basis, CeedOperator op);

std::vector<CeedInt> QuadratureDataSetup(const IntegratorInfo &info, Ceed ceed,
                                         CeedElemRestriction trial_restr,
                                         CeedBasis trial_basis, CeedVector *qdata_vec,
                                         CeedElemRestriction *qdata_restr);

void QuadratureDataAssembly(const std::vector<CeedInt> &qf_active_sizes,
                            const IntegratorInfo &info, Ceed ceed,
                            CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                            CeedBasis trial_basis, CeedBasis test_basis,
                            CeedVector qdata_vec, CeedElemRestriction qdata_restr,
                            CeedOperator *op);

}  // namespace internal

// Create libCEED operator using the given quadrature data and element restriction.
template <typename IntegratorContext>
inline void AssembleCeedOperator(const IntegratorInfo &info, const IntegratorContext &ctx,
                                 const CeedGeomFactorData &geom_data, Ceed ceed,
                                 CeedElemRestriction trial_restr,
                                 CeedElemRestriction test_restr, CeedBasis trial_basis,
                                 CeedBasis test_basis, CeedOperator *op)
{
  // If we are going to be assembling the quadrature data, construct the storage vector for
  // it (to be owned by the operator).
  CeedVector qdata_vec = nullptr;
  CeedElemRestriction qdata_restr = nullptr;
  std::vector<CeedInt> qf_active_sizes;
  if (info.assemble_qdata)
  {
    qf_active_sizes = internal::QuadratureDataSetup(info, ceed, trial_restr, trial_basis,
                                                    &qdata_vec, &qdata_restr);
  }

  // Create the QFunction that defines the action of the operator (or its setup).
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

  // Active inputs/outputs
  if (!info.assemble_qdata)
  {
    internal::AddQFunctionActiveInputsOutputs(info, ceed, trial_basis, test_basis,
                                              apply_qf);
  }
  else
  {
    CeedInt qdata_size;
    PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(qdata_restr, &qdata_size));
    PalaceCeedCall(ceed,
                   CeedQFunctionAddOutput(apply_qf, "qdata", qdata_size, CEED_EVAL_NONE));
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

  if (!info.assemble_qdata)
  {
    internal::AddOperatorActiveFields(info, ceed, trial_restr, test_restr, trial_basis,
                                      test_basis, *op);
  }
  else
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "qdata", qdata_restr, CEED_BASIS_NONE,
                                              CEED_VECTOR_ACTIVE));
  }

  PalaceCeedCall(ceed, CeedOperatorCheckReady(*op));

  // Assemble the quadrature data and create the actual operator.
  if (info.assemble_qdata)
  {
    internal::QuadratureDataAssembly(qf_active_sizes, info, ceed, trial_restr, test_restr,
                                     trial_basis, test_basis, qdata_vec, qdata_restr, op);
  }

  // Cleanup (these are now owned by the operator).
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&qdata_restr));
  PalaceCeedCall(ceed, CeedVectorDestroy(&qdata_vec));
}

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_INTEGRATOR_HPP
