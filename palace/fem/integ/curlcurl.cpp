// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"
#include "fem/libceed/utils.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/hdiv_qf.h"
#include "fem/qfunctions/l2_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

void CurlCurlIntegrator::Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                                  CeedElemRestriction trial_restr,
                                  CeedElemRestriction test_restr, CeedBasis trial_basis,
                                  CeedBasis test_basis, CeedOperator *op)
{
  ceed::IntegratorInfo info;

  // Set up geometry factor quadrature data.
  MFEM_VERIFY(geom_data->wdetJ_vec && geom_data->wdetJ_restr,
              "Missing geometry factor quadrature data for CurlCurlIntegrator!");
  info.geom_info = ceed::GeomFactorInfo::Determinant;
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
  MFEM_VERIFY(trial_ncomp == test_ncomp && trial_ncomp == 1,
              "CurlCurlIntegrator requires test and trial spaces with a single component!");
  switch (10 * geom_data->space_dim + geom_data->dim)
  {
    case 22:
      info.apply_qf = f_apply_l2_1;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_l2_1_loc);
      break;
    case 33:
      info.apply_qf = f_apply_hdiv_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdiv_33_loc);
      break;
    case 32:
      info.apply_qf = f_apply_l2_1;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_l2_1_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = (" << geom_data->dim << ", "
                                                         << geom_data->space_dim
                                                         << ") for CurlCurlIntegrator!");
  }
  info.trial_ops = ceed::EvalMode::Curl;
  info.test_ops = ceed::EvalMode::Curl;

  // Set up the coefficient and assemble.
  switch (geom_data->dim)
  {
    case 2:
      {
        auto ctx = ceed::PopulateCoefficientContext1();
        ceed::AssembleCeedOperator(info, ctx, geom_data, ceed, trial_restr, test_restr,
                                   trial_basis, test_basis, op);
      }
      break;
    case 3:
      {
        auto ctx = ceed::PopulateCoefficientContext3();
        ceed::AssembleCeedOperator(info, ctx, geom_data, ceed, trial_restr, test_restr,
                                   trial_basis, test_basis, op);
      }
      break;
  }
}

}  // namespace palace
