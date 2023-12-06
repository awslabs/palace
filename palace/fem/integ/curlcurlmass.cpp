// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/hdivmass_qf.h"

namespace palace
{

void CurlCurlMassIntegrator::Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                                      CeedElemRestriction trial_restr,
                                      CeedElemRestriction test_restr, CeedBasis trial_basis,
                                      CeedBasis test_basis, CeedOperator *op) const
{
  ceed::IntegratorInfo info;

  // Set up geometry factor quadrature data.
  MFEM_VERIFY(geom_data->wdetJ_vec && geom_data->wdetJ_restr && geom_data->adjJt_vec &&
                  geom_data->adjJt_restr,
              "Missing geometry factor quadrature data for CurlCurlIntegrator!");
  info.geom_info = ceed::GeomFactorInfo::Determinant | ceed::GeomFactorInfo::Adjugate;
  if (geom_data->dim == 3)
  {
    MFEM_VERIFY(geom_data->J_vec && geom_data->J_restr,
                "Missing geometry factor quadrature data for CurlCurlIntegrator!");
    info.geom_info |= ceed::GeomFactorInfo::Jacobian;
  }
  else
  {
    // Curl in 2D has a single component.
    info.geom_info |= ceed::GeomFactorInfo::Weight;
  }

  // Set up QFunctions.
  CeedInt trial_ncomp, test_ncomp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_ncomp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_ncomp));
  MFEM_VERIFY(
      trial_ncomp == test_ncomp && trial_ncomp == 1,
      "CurlCurlMassIntegrator requires test and trial spaces with a single component!");
  switch (10 * geom_data->space_dim + geom_data->dim)
  {
    case 22:
      info.apply_qf = f_apply_hdivmass_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdivmass_22_loc);
      break;
    case 33:
      info.apply_qf = f_apply_hdivmass_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdivmass_33_loc);
      break;
    case 32:
      info.apply_qf = f_apply_hdivmass_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdivmass_32_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << geom_data->dim << ", " << geom_data->space_dim
                 << ") for CurlCurlMassIntegrator!");
  }
  info.trial_ops = ceed::EvalMode::Curl | ceed::EvalMode::Interp;
  info.test_ops = ceed::EvalMode::Curl | ceed::EvalMode::Interp;

  // Set up the coefficient and assemble.
  switch (geom_data->space_dim)
  {
    case 2:
      {
        MatCoeffPairContext21 ctx{ceed::PopulateCoefficientContext2(Q_mass),
                                  ceed::PopulateCoefficientContext1(Q)};
        ceed::AssembleCeedOperator(info, ctx, geom_data, ceed, trial_restr, test_restr,
                                   trial_basis, test_basis, op);
      }
      break;
    case 3:
      {
        MatCoeffPairContext33 ctx{ceed::PopulateCoefficientContext3(Q_mass),
                                  ceed::PopulateCoefficientContext3(Q)};
        ceed::AssembleCeedOperator(info, ctx, geom_data, ceed, trial_restr, test_restr,
                                   trial_basis, test_basis, op);
      }
      break;
  }
}

}  // namespace palace
