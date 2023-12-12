// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/hcurl_build_qf.h"
#include "fem/qfunctions/hcurl_qf.h"

namespace palace
{

void MixedVectorGradientIntegrator::Assemble(const ceed::CeedGeomFactorData &geom_data,
                                             Ceed ceed, CeedElemRestriction trial_restr,
                                             CeedElemRestriction test_restr,
                                             CeedBasis trial_basis, CeedBasis test_basis,
                                             CeedOperator *op) const
{
  ceed::IntegratorInfo info;
  info.assemble_qdata = assemble_qdata;

  // Set up geometry factor quadrature data.
  MFEM_VERIFY(geom_data->wdetJ_vec && geom_data->wdetJ_restr && geom_data->adjJt_vec &&
                  geom_data->adjJt_restr,
              "Missing geometry factor quadrature data for MixedVectorGradientIntegrator!");
  info.geom_info = ceed::GeomFactorInfo::Determinant | ceed::GeomFactorInfo::Adjugate;

  // Set up QFunctions.
  CeedInt trial_ncomp, test_ncomp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_ncomp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_ncomp));
  MFEM_VERIFY(trial_ncomp == test_ncomp && trial_ncomp == 1,
              "MixedVectorGradientIntegrator requires test and trial spaces with a single "
              "component!");
  switch (10 * geom_data->space_dim + geom_data->dim)
  {
    case 22:
      info.apply_qf = assemble_qdata ? f_build_hcurl_22 : f_apply_hcurl_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurl_22_loc : f_apply_hcurl_22_loc);
      break;
    case 33:
      info.apply_qf = assemble_qdata ? f_build_hcurl_33 : f_apply_hcurl_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurl_33_loc : f_apply_hcurl_33_loc);
      break;
    case 21:
      info.apply_qf = assemble_qdata ? f_build_hcurl_21 : f_apply_hcurl_21;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurl_21_loc : f_apply_hcurl_21_loc);
      break;
    case 32:
      info.apply_qf = assemble_qdata ? f_build_hcurl_32 : f_apply_hcurl_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurl_32_loc : f_apply_hcurl_32_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << geom_data->dim << ", " << geom_data->space_dim
                 << ") for MixedVectorGradientIntegrator!");
  }
  info.trial_ops = ceed::EvalMode::Grad;
  info.test_ops = ceed::EvalMode::Interp;

  // Set up the coefficient and assemble.
  auto ctx = [&]()
  {
    switch (geom_data->space_dim)
    {
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

void MixedVectorWeakDivergenceIntegrator::Assemble(
    const ceed::CeedGeomFactorData &geom_data, Ceed ceed, CeedElemRestriction trial_restr,
    CeedElemRestriction test_restr, CeedBasis trial_basis, CeedBasis test_basis,
    CeedOperator *op) const
{
  ceed::IntegratorInfo info;
  info.assemble_qdata = assemble_qdata;

  // Set up geometry factor quadrature data.
  MFEM_VERIFY(
      geom_data->wdetJ_vec && geom_data->wdetJ_restr && geom_data->adjJt_vec &&
          geom_data->adjJt_restr,
      "Missing geometry factor quadrature data for MixedVectorWeakDivergenceIntegrator!");
  info.geom_info = ceed::GeomFactorInfo::Determinant | ceed::GeomFactorInfo::Adjugate;

  // Set up QFunctions.
  CeedInt trial_ncomp, test_ncomp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_ncomp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_ncomp));
  MFEM_VERIFY(
      trial_ncomp == test_ncomp && trial_ncomp == 1,
      "MixedVectorWeakDivergenceIntegrator requires test and trial spaces with a single "
      "component!");
  switch (10 * geom_data->space_dim + geom_data->dim)
  {
    case 22:
      info.apply_qf = assemble_qdata ? f_build_hcurl_22 : f_apply_hcurl_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurl_22_loc : f_apply_hcurl_22_loc);
      break;
    case 33:
      info.apply_qf = assemble_qdata ? f_build_hcurl_33 : f_apply_hcurl_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurl_33_loc : f_apply_hcurl_33_loc);
      break;
    case 21:
      info.apply_qf = assemble_qdata ? f_build_hcurl_21 : f_apply_hcurl_21;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurl_21_loc : f_apply_hcurl_21_loc);
      break;
    case 32:
      info.apply_qf = assemble_qdata ? f_build_hcurl_32 : f_apply_hcurl_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_qdata ? f_build_hcurl_32_loc : f_apply_hcurl_32_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << geom_data->dim << ", " << geom_data->space_dim
                 << ") for MixedVectorWeakDivergenceIntegrator!");
  }
  info.trial_ops = ceed::EvalMode::Interp;
  info.test_ops = ceed::EvalMode::Grad;

  // Set up the coefficient and assemble.
  auto ctx = [&]()
  {
    switch (geom_data->space_dim)
    {
      case 2:
        return ceed::PopulateCoefficientContext<2>(Q, -1.0);
      case 3:
        return ceed::PopulateCoefficientContext<3>(Q, -1.0);
    }
    return std::vector<CeedIntScalar>();
  }();
  ceed::AssembleCeedOperator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar),
                             geom_data, ceed, trial_restr, test_restr, trial_basis,
                             test_basis, op);
}

}  // namespace palace
