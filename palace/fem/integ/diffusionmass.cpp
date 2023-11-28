// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/hcurlmass_qf.h"

namespace palace
{

void DiffusionMassIntegrator::Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                                       CeedElemRestriction trial_restr,
                                       CeedElemRestriction test_restr,
                                       CeedBasis trial_basis, CeedBasis test_basis,
                                       CeedOperator *op)
{
  ceed::IntegratorInfo info;

  // Set up geometry factor quadrature data.
  MFEM_VERIFY(geom_data->wdetJ_vec && geom_data->wdetJ_restr && geom_data->adjJt_vec &&
                  geom_data->adjJt_restr,
              "Missing geometry factor quadrature data for DiffusionMassIntegrator!");
  info.geom_info = ceed::GeomFactorInfo::Determinant | ceed::GeomFactorInfo::Adjugate;

  // Set up QFunctions.
  CeedInt trial_ncomp, test_ncomp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_ncomp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_ncomp));
  MFEM_VERIFY(
      trial_ncomp == test_ncomp && trial_ncomp == 1,
      "DiffusionMassIntegrator requires test and trial spaces with a single component!");
  switch (10 * geom_data->space_dim + geom_data->dim)
  {
    case 22:
      info.apply_qf = f_apply_hcurlmass_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlmass_22_loc);
      break;
    case 33:
      info.apply_qf = f_apply_hcurlmass_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlmass_33_loc);
      break;
    case 21:
      info.apply_qf = f_apply_hcurlmass_21;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlmass_21_loc);
      break;
    case 32:
      info.apply_qf = f_apply_hcurlmass_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlmass_32_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << geom_data->dim << ", " << geom_data->space_dim
                 << ") for DiffusionMassIntegrator!");
  }
  info.trial_ops = ceed::EvalMode::Grad | ceed::EvalMode::Interp;
  info.test_ops = ceed::EvalMode::Grad | ceed::EvalMode::Interp;

  // Set up the coefficient and assemble.
  switch (geom_data->space_dim)
  {
    case 2:
      {
        MatCoeffPairContext21 ctx{ceed::PopulateCoefficientContext2(),
                                  ceed::PopulateCoefficientContext1()};
        ceed::AssembleCeedOperator(info, ctx, geom_data, ceed, trial_restr, test_restr,
                                   trial_basis, test_basis, op);
      }
      break;
    case 3:
      {
        MatCoeffPairContext31 ctx{ceed::PopulateCoefficientContext3(),
                                  ceed::PopulateCoefficientContext1()};
        ceed::AssembleCeedOperator(info, ctx, geom_data, ceed, trial_restr, test_restr,
                                   trial_basis, test_basis, op);
      }
      break;
  }
}

}  // namespace palace
