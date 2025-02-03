// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/h1_qf.h"

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
  switch (trial_num_comp)
  {
    case 1:
      info.apply_qf = assemble_q_data ? f_build_h1_1 : f_apply_h1_1;
      info.apply_qf_path = PalaceQFunctionRelativePath(assemble_q_data ? f_build_h1_1_loc
                                                                       : f_apply_h1_1_loc);
      break;
    case 2:
      info.apply_qf = assemble_q_data ? f_build_h1_2 : f_apply_h1_2;
      info.apply_qf_path = PalaceQFunctionRelativePath(assemble_q_data ? f_build_h1_2_loc
                                                                       : f_apply_h1_2_loc);
      break;
    case 3:
      info.apply_qf = assemble_q_data ? f_build_h1_3 : f_apply_h1_3;
      info.apply_qf_path = PalaceQFunctionRelativePath(assemble_q_data ? f_build_h1_3_loc
                                                                       : f_apply_h1_3_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of num_comp = " << trial_num_comp
                                                << " for MassIntegrator!");
  }
  info.trial_ops = EvalMode::Interp;
  info.test_ops = EvalMode::Interp;

  // Set up the coefficient and assemble.
  auto ctx = PopulateCoefficientContext(trial_num_comp, Q, transpose);
  AssembleCeedOperator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed,
                       trial_restr, test_restr, trial_basis, test_basis, geom_data,
                       geom_data_restr, op);
}

}  // namespace palace
