// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "functional.hpp"

#include <ceed/backend.h>
#include <mfem.hpp>

namespace palace::ceed
{

namespace
{

void AddQFunctionFieldInput(const CeedFunctionalFieldInput &input, Ceed ceed,
                            CeedQFunction qf)
{
  if (input.ops & EvalMode::None)
  {
    // EvalMode::None inputs may have no basis (e.g. quadrature point data passed
    // through directly), so get the component count from the restriction.
    CeedInt num_comp;
    PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(input.restr, &num_comp));
    PalaceCeedCall(ceed,
                   CeedQFunctionAddInput(qf, input.name.c_str(), num_comp, CEED_EVAL_NONE));
  }
  if (input.ops & EvalMode::Interp)
  {
    CeedInt num_comp, q_comp;
    PalaceCeedCall(ceed, CeedBasisGetNumComponents(input.basis, &num_comp));
    PalaceCeedCall(
        ceed, CeedBasisGetNumQuadratureComponents(input.basis, CEED_EVAL_INTERP, &q_comp));
    PalaceCeedCall(ceed, CeedQFunctionAddInput(qf, input.name.c_str(), num_comp * q_comp,
                                               CEED_EVAL_INTERP));
  }
  MFEM_VERIFY(
      !(input.ops & (EvalMode::Grad | EvalMode::Div | EvalMode::Curl | EvalMode::Weight)),
      "Unsupported evaluation mode for surface functional field input '" << input.name
                                                                         << "'!");
}

void AddOperatorFieldInput(const CeedFunctionalFieldInput &input, Ceed ceed,
                           CeedOperator op)
{
  if (input.ops & EvalMode::None)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(op, input.name.c_str(), input.restr,
                                              CEED_BASIS_NONE, input.vec));
  }
  if (input.ops & EvalMode::Interp)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(op, input.name.c_str(), input.restr,
                                              input.basis, input.vec));
  }
}

}  // namespace

void AssembleCeedSurfaceFunctional(
    const CeedQFunctionInfo &info, void *ctx, std::size_t ctx_size, Ceed ceed,
    const std::vector<CeedFunctionalFieldInput> &inputs, CeedVector face_geom_data,
    CeedElemRestriction face_geom_data_restr, CeedVector vol_geom_data,
    CeedElemRestriction vol_geom_data_restr, CeedInt num_out_comp,
    CeedElemRestriction out_restr, CeedOperator *op)
{
  MFEM_VERIFY(!info.assemble_q_data,
              "Surface functional integrator does not support quadrature data assembly!");
  MFEM_VERIFY((vol_geom_data == nullptr) == (vol_geom_data_restr == nullptr),
              "Mismatched volume geometry data and restriction for surface functional!");

  // Create basis for summing contributions from all quadrature points on each element
  // for each output component (all-ones interpolation matrix).
  CeedInt num_qpts;
  CeedBasis sum_basis;
  {
    CeedBasis qpts_basis = nullptr;
    for (const auto &input : inputs)
    {
      if (input.basis)
      {
        qpts_basis = input.basis;
        break;
      }
    }
    if (qpts_basis)
    {
      PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(qpts_basis, &num_qpts));
    }
    else
    {
      // No field inputs with a basis: deduce the number of quadrature points from the
      // face geometry data restriction (one set of geometry data per quadrature point).
      CeedInt num_elem;
      CeedSize l_size;
      CeedInt num_geom_comp;
      PalaceCeedCall(ceed,
                     CeedElemRestrictionGetNumElements(face_geom_data_restr, &num_elem));
      PalaceCeedCall(
          ceed, CeedElemRestrictionGetNumComponents(face_geom_data_restr, &num_geom_comp));
      PalaceCeedCall(ceed,
                     CeedElemRestrictionGetLVectorSize(face_geom_data_restr, &l_size));
      num_qpts = static_cast<CeedInt>(l_size / (num_elem * num_geom_comp));
    }
    mfem::Vector Bt(num_qpts), Gt(num_qpts), qX(num_qpts), qW(num_qpts);
    Bt = 1.0;
    Gt = 0.0;
    qX = 0.0;
    qW = 0.0;
    // Note: ceed::GetCeedTopology(CEED_TOPOLOGY_LINE) == 1.
    PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, num_out_comp, 1,
                                           num_qpts, Bt.GetData(), Gt.GetData(),
                                           qX.GetData(), qW.GetData(), &sum_basis));
  }

  // Create the QFunction that defines the action of the operator.
  CeedQFunction apply_qf;
  PalaceCeedCall(ceed, CeedQFunctionCreateInterior(ceed, 1, info.apply_qf,
                                                   info.apply_qf_path.c_str(), &apply_qf));

  if (ctx && ctx_size > 0)
  {
    CeedQFunctionContext apply_ctx;
    PalaceCeedCall(ceed, CeedQFunctionContextCreate(ceed, &apply_ctx));
    PalaceCeedCall(ceed, CeedQFunctionContextSetData(apply_ctx, CEED_MEM_HOST,
                                                     CEED_COPY_VALUES, ctx_size, ctx));
    PalaceCeedCall(ceed, CeedQFunctionSetContext(apply_qf, apply_ctx));
    PalaceCeedCall(ceed, CeedQFunctionContextDestroy(&apply_ctx));
  }

  // Inputs/outputs.
  {
    CeedInt geom_data_size;
    PalaceCeedCall(
        ceed, CeedElemRestrictionGetNumComponents(face_geom_data_restr, &geom_data_size));
    PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "face_geom_data", geom_data_size,
                                               CEED_EVAL_NONE));
  }
  if (vol_geom_data)
  {
    CeedInt geom_data_size;
    PalaceCeedCall(
        ceed, CeedElemRestrictionGetNumComponents(vol_geom_data_restr, &geom_data_size));
    PalaceCeedCall(ceed, CeedQFunctionAddInput(apply_qf, "vol_geom_data", geom_data_size,
                                               CEED_EVAL_NONE));
  }
  for (const auto &input : inputs)
  {
    AddQFunctionFieldInput(input, ceed, apply_qf);
  }
  PalaceCeedCall(ceed,
                 CeedQFunctionAddOutput(apply_qf, "v", num_out_comp, CEED_EVAL_INTERP));

  // Create the operator.
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, apply_qf, nullptr, nullptr, op));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&apply_qf));

  PalaceCeedCall(ceed, CeedOperatorSetField(*op, "face_geom_data", face_geom_data_restr,
                                            CEED_BASIS_NONE, face_geom_data));
  if (vol_geom_data)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(*op, "vol_geom_data", vol_geom_data_restr,
                                              CEED_BASIS_NONE, vol_geom_data));
  }
  for (const auto &input : inputs)
  {
    AddOperatorFieldInput(input, ceed, *op);
  }
  PalaceCeedCall(ceed,
                 CeedOperatorSetField(*op, "v", out_restr, sum_basis, CEED_VECTOR_ACTIVE));

  PalaceCeedCall(ceed, CeedOperatorCheckReady(*op));

  // Cleanup (this is now owned by the operator).
  PalaceCeedCall(ceed, CeedBasisDestroy(&sum_basis));
}

void AssembleCeedPointEvaluator(const CeedQFunctionInfo &info, void *ctx,
                                std::size_t ctx_size, Ceed ceed,
                                const std::vector<CeedFunctionalFieldInput> &inputs,
                                CeedVector geom_data, CeedElemRestriction geom_data_restr,
                                CeedInt num_out_comp, CeedElemRestriction out_restr,
                                CeedOperator *op)
{
  MFEM_VERIFY(!info.assemble_q_data,
              "Point evaluator does not support quadrature data assembly!");

  // Create the QFunction that defines the action of the operator.
  CeedQFunction apply_qf;
  PalaceCeedCall(ceed, CeedQFunctionCreateInterior(ceed, 1, info.apply_qf,
                                                   info.apply_qf_path.c_str(), &apply_qf));

  if (ctx && ctx_size > 0)
  {
    CeedQFunctionContext apply_ctx;
    PalaceCeedCall(ceed, CeedQFunctionContextCreate(ceed, &apply_ctx));
    PalaceCeedCall(ceed, CeedQFunctionContextSetData(apply_ctx, CEED_MEM_HOST,
                                                     CEED_COPY_VALUES, ctx_size, ctx));
    PalaceCeedCall(ceed, CeedQFunctionSetContext(apply_qf, apply_ctx));
    PalaceCeedCall(ceed, CeedQFunctionContextDestroy(&apply_ctx));
  }

  // Inputs/outputs.
  {
    CeedInt geom_data_size;
    PalaceCeedCall(ceed,
                   CeedElemRestrictionGetNumComponents(geom_data_restr, &geom_data_size));
    PalaceCeedCall(
        ceed, CeedQFunctionAddInput(apply_qf, "geom_data", geom_data_size, CEED_EVAL_NONE));
  }
  for (const auto &input : inputs)
  {
    AddQFunctionFieldInput(input, ceed, apply_qf);
  }
  PalaceCeedCall(ceed, CeedQFunctionAddOutput(apply_qf, "v", num_out_comp, CEED_EVAL_NONE));

  // Create the operator.
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, apply_qf, nullptr, nullptr, op));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&apply_qf));

  PalaceCeedCall(ceed, CeedOperatorSetField(*op, "geom_data", geom_data_restr,
                                            CEED_BASIS_NONE, geom_data));
  for (const auto &input : inputs)
  {
    AddOperatorFieldInput(input, ceed, *op);
  }
  PalaceCeedCall(
      ceed, CeedOperatorSetField(*op, "v", out_restr, CEED_BASIS_NONE, CEED_VECTOR_ACTIVE));

  PalaceCeedCall(ceed, CeedOperatorCheckReady(*op));
}

}  // namespace palace::ceed
