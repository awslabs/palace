// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/hdivmass_qf.h"

namespace palace
{

using namespace ceed;

void CurlCurlMassIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                      CeedElemRestriction test_restr, CeedBasis trial_basis,
                                      CeedBasis test_basis, CeedVector geom_data,
                                      CeedElemRestriction geom_data_restr,
                                      CeedOperator *op) const
{
  CeedQFunctionInfo info;
  info.assemble_q_data = assemble_q_data;

  // Set up QFunctions.
  CeedInt dim, space_dim, trial_num_comp, test_num_comp;
  PalaceCeedCall(ceed, CeedBasisGetDimension(trial_basis, &dim));
  PalaceCeedCall(ceed, CeedGeometryDataGetSpaceDimension(geom_data_restr, dim, &space_dim));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_num_comp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_num_comp));
  MFEM_VERIFY(
      trial_num_comp == test_num_comp && trial_num_comp == 1,
      "CurlCurlMassIntegrator requires test and trial spaces with a single component!");
  switch (10 * space_dim + dim)
  {
    case 22:
      info.apply_qf = assemble_q_data ? f_build_hdivmass_22 : f_apply_hdivmass_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_q_data ? f_build_hdivmass_22_loc : f_apply_hdivmass_22_loc);
      break;
    case 33:
      info.apply_qf = assemble_q_data ? f_build_hdivmass_33 : f_apply_hdivmass_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_q_data ? f_build_hdivmass_33_loc : f_apply_hdivmass_33_loc);
      break;
    case 32:
      info.apply_qf = assemble_q_data ? f_build_hdivmass_32 : f_apply_hdivmass_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_q_data ? f_build_hdivmass_32_loc : f_apply_hdivmass_32_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << dim << ", " << space_dim << ") for CurlCurlMassIntegrator!");
  }
  info.trial_ops = EvalMode::Curl | EvalMode::Interp;
  info.test_ops = EvalMode::Curl | EvalMode::Interp;
  if (dim < 3)
  {
    info.trial_ops |= EvalMode::Weight;
  }

  // Set up the coefficient and assemble. Mass goes first.
  auto ctx = PopulateCoefficientContext(space_dim, Q_mass, (dim < 3) ? 1 : dim, Q,
                                        transpose_mass, transpose);
  AssembleCeedOperator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed,
                       trial_restr, test_restr, trial_basis, test_basis, geom_data,
                       geom_data_restr, op);
}

}  // namespace palace
