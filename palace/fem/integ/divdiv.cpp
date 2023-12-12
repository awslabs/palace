// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/l2_build_qf.h"
#include "fem/qfunctions/l2_qf.h"

namespace palace
{

void DivDivIntegrator::Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                                CeedElemRestriction trial_restr,
                                CeedElemRestriction test_restr, CeedBasis trial_basis,
                                CeedBasis test_basis, CeedOperator *op) const
{
  ceed::IntegratorInfo info;
  info.assemble_qdata = assemble_qdata;

  // Set up geometry factor quadrature data.
  MFEM_VERIFY(geom_data->wdetJ_vec && geom_data->wdetJ_restr,
              "Missing geometry factor quadrature data for DivDivIntegrator!");
  info.geom_info = ceed::GeomFactorInfo::Determinant | ceed::GeomFactorInfo::Weight;

  // Set up QFunctions.
  CeedInt trial_ncomp, test_ncomp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_ncomp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_ncomp));
  MFEM_VERIFY(
      trial_ncomp == test_ncomp,
      "DivDivIntegrator requires test and trial spaces with same number of components!");
  switch (trial_ncomp)
  {
    case 1:
      info.apply_qf = assemble_qdata ? f_build_l2_1 : f_apply_l2_1;
      info.apply_qf_path =
          PalaceQFunctionRelativePath(assemble_qdata ? f_build_l2_1_loc : f_apply_l2_1_loc);
      break;
    case 2:
      info.apply_qf = assemble_qdata ? f_build_l2_2 : f_apply_l2_2;
      info.apply_qf_path =
          PalaceQFunctionRelativePath(assemble_qdata ? f_build_l2_2_loc : f_apply_l2_2_loc);
      break;
    case 3:
      info.apply_qf = assemble_qdata ? f_build_l2_3 : f_apply_l2_3;
      info.apply_qf_path =
          PalaceQFunctionRelativePath(assemble_qdata ? f_build_l2_3_loc : f_apply_l2_3_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of ncomp = " << trial_ncomp << " for DivDivIntegrator!");
  }
  info.trial_ops = ceed::EvalMode::Div;
  info.test_ops = ceed::EvalMode::Div;

  // Set up the coefficient and assemble.
  auto ctx = [&]()
  {
    switch (trial_ncomp)
    {
      case 1:
        return ceed::PopulateCoefficientContext<1>(Q);
      case 2:
        return ceed::PopulateCoefficientContext<2>(Q);
      case 3:
        return ceed::PopulateCoefficientContext<3>(Q);
    }
    return std::vector<CeedIntScalar>();
  }();
  ceed::AssembleCeedOperator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar),
                             geom_data, ceed, trial_restr, test_restr, trial_basis,
                             test_basis, op);
}

}  // namespace palace
