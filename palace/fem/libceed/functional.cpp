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
  if (input.ops & EvalMode::Grad)
  {
    CeedInt num_comp, q_comp;
    PalaceCeedCall(ceed, CeedBasisGetNumComponents(input.basis, &num_comp));
    PalaceCeedCall(
        ceed, CeedBasisGetNumQuadratureComponents(input.basis, CEED_EVAL_GRAD, &q_comp));
    PalaceCeedCall(ceed, CeedQFunctionAddInput(qf, ("grad_" + input.name).c_str(),
                                               num_comp * q_comp, CEED_EVAL_GRAD));
  }
  if (input.ops & EvalMode::Weight)
  {
    PalaceCeedCall(ceed,
                   CeedQFunctionAddInput(qf, input.name.c_str(), 1, CEED_EVAL_WEIGHT));
  }
  MFEM_VERIFY(!(input.ops & (EvalMode::Div | EvalMode::Curl)),
              "Unsupported evaluation mode for pointwise operator field input '"
                  << input.name << "'!");
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
  if (input.ops & EvalMode::Grad)
  {
    PalaceCeedCall(ceed, CeedOperatorSetField(op, ("grad_" + input.name).c_str(),
                                              input.restr, input.basis, input.vec));
  }
  if (input.ops & EvalMode::Weight)
  {
    PalaceCeedCall(ceed,
                   CeedOperatorSetField(op, input.name.c_str(), CEED_ELEMRESTRICTION_NONE,
                                        input.basis, CEED_VECTOR_NONE));
  }
}

}  // namespace

namespace
{

void CreatePointEvaluatorQFunction(const CeedQFunctionInfo &info, void *ctx,
                                   std::size_t ctx_size, Ceed ceed,
                                   const std::vector<CeedFunctionalFieldInput> &inputs,
                                   CeedInt num_out_comp, CeedQFunction *apply_qf)
{
  MFEM_VERIFY(!info.assemble_q_data,
              "Point evaluator does not support quadrature data assembly!");
  PalaceCeedCall(ceed, CeedQFunctionCreateInterior(ceed, 1, info.apply_qf,
                                                   info.apply_qf_path.c_str(), apply_qf));

  if (ctx && ctx_size > 0)
  {
    CeedQFunctionContext apply_ctx;
    PalaceCeedCall(ceed, CeedQFunctionContextCreate(ceed, &apply_ctx));
    PalaceCeedCall(ceed, CeedQFunctionContextSetData(apply_ctx, CEED_MEM_HOST,
                                                     CEED_COPY_VALUES, ctx_size, ctx));
    PalaceCeedCall(ceed, CeedQFunctionSetContext(*apply_qf, apply_ctx));
    PalaceCeedCall(ceed, CeedQFunctionContextDestroy(&apply_ctx));
  }

  for (const auto &input : inputs)
  {
    AddQFunctionFieldInput(input, ceed, *apply_qf);
  }
  PalaceCeedCall(ceed,
                 CeedQFunctionAddOutput(*apply_qf, "v", num_out_comp, CEED_EVAL_NONE));
}

}  // namespace

void AssembleCeedPointEvaluator(const CeedQFunctionInfo &info, void *ctx,
                                std::size_t ctx_size, Ceed ceed,
                                const std::vector<CeedFunctionalFieldInput> &inputs,
                                CeedInt num_out_comp, CeedElemRestriction out_restr,
                                CeedOperator *op)
{
  CeedQFunction apply_qf;
  CreatePointEvaluatorQFunction(info, ctx, ctx_size, ceed, inputs, num_out_comp, &apply_qf);

  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, apply_qf, nullptr, nullptr, op));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&apply_qf));

  for (const auto &input : inputs)
  {
    AddOperatorFieldInput(input, ceed, *op);
  }
  PalaceCeedCall(
      ceed, CeedOperatorSetField(*op, "v", out_restr, CEED_BASIS_NONE, CEED_VECTOR_ACTIVE));

  PalaceCeedCall(ceed, CeedOperatorCheckReady(*op));
}

}  // namespace palace::ceed
