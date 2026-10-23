// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/h1_qf.h"
#include "fem/qfunctions/l2_qf.h"
#include "fem/qfunctions/l2h1_qf.h"

namespace palace
{

using namespace ceed;

void MassIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                              CeedElemRestriction test_restr, CeedBasis trial_basis,
                              CeedBasis test_basis, CeedVector geom_data,
                              CeedElemRestriction geom_data_restr, CeedOperator *op) const
{
  CeedQFunctionInfo info;
  info.assemble_q_data = assemble_q_data;

  // Set up QFunctions.
  CeedInt trial_num_comp, test_num_comp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_num_comp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_num_comp));
  MFEM_VERIFY(
      trial_num_comp == test_num_comp,
      "MassIntegrator requires test and trial spaces with same number of components!");
  // L2 spaces with INTEGRAL map type (u = 1 / det(J) ̂u) use the L2 QFunctions, which need
  // the quadrature weights.
  const bool trial_l2 = (trial_map_type == mfem::FiniteElement::INTEGRAL),
             test_l2 = (test_map_type == mfem::FiniteElement::INTEGRAL);
  if (trial_l2 && test_l2)
  {
    switch (trial_num_comp)
    {
      case 1:
        info.apply_qf = assemble_q_data ? f_build_l2_1 : f_apply_l2_1;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_l2_1_loc : f_apply_l2_1_loc);
        break;
      case 2:
        info.apply_qf = assemble_q_data ? f_build_l2_2 : f_apply_l2_2;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_l2_2_loc : f_apply_l2_2_loc);
        break;
      case 3:
        info.apply_qf = assemble_q_data ? f_build_l2_3 : f_apply_l2_3;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_l2_3_loc : f_apply_l2_3_loc);
        break;
      default:
        MFEM_ABORT("Invalid value of num_comp = " << trial_num_comp
                                                  << " for MassIntegrator!");
    }
  }
  else if (trial_l2 || test_l2)
  {
    MFEM_VERIFY(trial_num_comp == 1,
                "MassIntegrator only supports scalar mixed L2-H1 spaces!");
    info.apply_qf = assemble_q_data ? f_build_l2h1_1 : f_apply_l2h1_1;
    info.apply_qf_path = PalaceQFunctionRelativePath(assemble_q_data ? f_build_l2h1_1_loc
                                                                     : f_apply_l2h1_1_loc);
  }
  else
  {
    switch (trial_num_comp)
    {
      case 1:
        info.apply_qf = assemble_q_data ? f_build_h1_1 : f_apply_h1_1;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_h1_1_loc : f_apply_h1_1_loc);
        break;
      case 2:
        info.apply_qf = assemble_q_data ? f_build_h1_2 : f_apply_h1_2;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_h1_2_loc : f_apply_h1_2_loc);
        break;
      case 3:
        info.apply_qf = assemble_q_data ? f_build_h1_3 : f_apply_h1_3;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_h1_3_loc : f_apply_h1_3_loc);
        break;
      default:
        MFEM_ABORT("Invalid value of num_comp = " << trial_num_comp
                                                  << " for MassIntegrator!");
    }
  }
  info.trial_ops = EvalMode::Interp;
  if (trial_l2 || test_l2)
  {
    info.trial_ops |= EvalMode::Weight;
  }
  info.test_ops = EvalMode::Interp;

  // Set up the coefficient and assemble.
  auto ctx = PopulateCoefficientContext(trial_num_comp, Q, transpose);
  AssembleCeedOperator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed,
                       trial_restr, test_restr, trial_basis, test_basis, geom_data,
                       geom_data_restr, op);
}

}  // namespace palace
