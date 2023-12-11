// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/hcurlh1d_build_qf.h"
#include "fem/qfunctions/hcurlh1d_qf.h"

namespace palace
{

void GradientIntegrator::Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                                  CeedElemRestriction trial_restr,
                                  CeedElemRestriction test_restr, CeedBasis trial_basis,
                                  CeedBasis test_basis, CeedOperator *op) const
{
  ceed::IntegratorInfo info;
  info.assemble_qdata = assemble_qdata;

  // Set up geometry factor quadrature data.
  MFEM_VERIFY(geom_data->wdetJ_vec && geom_data->wdetJ_restr && geom_data->adjJt_vec &&
                  geom_data->adjJt_restr,
              "Missing geometry factor quadrature data for GradientIntegrator!");
  info.geom_info = ceed::GeomFactorInfo::Determinant | ceed::GeomFactorInfo::Adjugate;

  // Set up QFunctions.
  CeedInt trial_ncomp, test_ncomp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_ncomp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_ncomp));
  MFEM_VERIFY(trial_ncomp == 1 && test_ncomp == geom_data->space_dim,
              "GradientIntegrator requires trial space with a single component and test "
              "space with space_dim components!");
  switch (10 * geom_data->space_dim + geom_data->dim)
  {
    case 22:
      info.apply_qf = assemble_qdata ? f_build_hcurlh1d_22 : f_apply_hcurlh1d_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurlh1d_22_loc : f_apply_hcurlh1d_22_loc);
      break;
    case 33:
      info.apply_qf = assemble_qdata ? f_build_hcurlh1d_33 : f_apply_hcurlh1d_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurlh1d_33_loc : f_apply_hcurlh1d_33_loc);
      break;
    case 21:
      info.apply_qf = assemble_qdata ? f_build_hcurlh1d_21 : f_apply_hcurlh1d_21;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurlh1d_21_loc : f_apply_hcurlh1d_21_loc);
      break;
    case 32:
      info.apply_qf = assemble_qdata ? f_build_hcurlh1d_32 : f_apply_hcurlh1d_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurlh1d_32_loc : f_apply_hcurlh1d_32_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = (" << geom_data->dim << ", "
                                                         << geom_data->space_dim
                                                         << ") for GradientIntegrator!");
  }
  info.trial_ops = ceed::EvalMode::Grad;
  info.test_ops = ceed::EvalMode::Interp;

  // Set up the coefficient and assemble.
  switch (geom_data->space_dim)
  {
    case 2:
      {
        auto ctx = ceed::PopulateCoefficientContext2(Q);
        ceed::AssembleCeedOperator(info, ctx, geom_data, ceed, trial_restr, test_restr,
                                   trial_basis, test_basis, op);
      }
      break;
    case 3:
      {
        auto ctx = ceed::PopulateCoefficientContext3(Q);
        ceed::AssembleCeedOperator(info, ctx, geom_data, ceed, trial_restr, test_restr,
                                   trial_basis, test_basis, op);
      }
      break;
  }
}

}  // namespace palace
