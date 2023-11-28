// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/integrator.hpp"
#include "fem/libceed/utils.hpp"

#include "fem/qfunctions/hcurlh1d_qf.h"

namespace palace
{

namespace
{

struct GradientIntegratorInfo : public ceed::IntegratorInfo
{
  bool ctx;  // XX TODO WIP COEFFICIENTS
};

}  // namespace

void GradientIntegrator::Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                                  CeedElemRestriction trial_restr,
                                  CeedElemRestriction test_restr, CeedBasis trial_basis,
                                  CeedBasis test_basis, CeedOperator *op)
{
  GradientIntegratorInfo info;

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
      info.apply_qf = f_apply_hcurlh1d_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlh1d_22_loc);
      break;
    case 33:
      info.apply_qf = f_apply_hcurlh1d_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlh1d_33_loc);
      break;
    case 21:
      info.apply_qf = f_apply_hcurlh1d_21;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlh1d_21_loc);
      break;
    case 32:
      info.apply_qf = f_apply_hcurlh1d_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlh1d_32_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = (" << geom_data->dim << ", "
                                                         << geom_data->space_dim
                                                         << ") for GradientIntegrator!");
  }
  info.trial_ops = ceed::EvalMode::Grad;
  info.test_ops = ceed::EvalMode::Interp;

  ceed::AssembleCeedOperator(info, geom_data, ceed, trial_restr, test_restr, trial_basis,
                             test_basis, op);
}

}  // namespace palace
