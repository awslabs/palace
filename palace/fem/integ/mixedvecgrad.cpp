// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/hcurl_qf.h"
#include "fem/qfunctions/hcurlhdiv_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

using namespace ceed;

void MixedVectorGradientIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                             CeedElemRestriction test_restr,
                                             CeedBasis trial_basis, CeedBasis test_basis,
                                             CeedVector geom_data,
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
  MFEM_VERIFY(trial_num_comp == test_num_comp && trial_num_comp == 1,
              "MixedVectorGradientIntegrator requires test and trial spaces with a single "
              "component!");
  switch (10 * space_dim + dim)
  {
    case 22:
      if (test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = assemble_q_data ? f_build_hcurl_22 : f_apply_hcurl_22;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_hcurl_22_loc : f_apply_hcurl_22_loc);
      }
      else if (test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = assemble_q_data ? f_build_hcurlhdiv_22 : f_apply_hcurlhdiv_22;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_hcurlhdiv_22_loc : f_apply_hcurlhdiv_22_loc);
      }
      else
      {
        MFEM_ABORT("Invalid test map type for MixedVectorGradientIntegrator!");
      }
      break;
    case 33:
      if (test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = assemble_q_data ? f_build_hcurl_33 : f_apply_hcurl_33;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_hcurl_33_loc : f_apply_hcurl_33_loc);
      }
      else if (test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = assemble_q_data ? f_build_hcurlhdiv_33 : f_apply_hcurlhdiv_33;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_hcurlhdiv_33_loc : f_apply_hcurlhdiv_33_loc);
      }
      else
      {
        MFEM_ABORT("Invalid test map type for MixedVectorGradientIntegrator!");
      }
      break;
    case 21:
      if (test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = assemble_q_data ? f_build_hcurl_21 : f_apply_hcurl_21;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_hcurl_21_loc : f_apply_hcurl_21_loc);
      }
      else if (test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = assemble_q_data ? f_build_hcurlhdiv_21 : f_apply_hcurlhdiv_21;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_hcurlhdiv_21_loc : f_apply_hcurlhdiv_21_loc);
      }
      else
      {
        MFEM_ABORT("Invalid test map type for MixedVectorGradientIntegrator!");
      }
      break;
    case 32:
      if (test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = assemble_q_data ? f_build_hcurl_32 : f_apply_hcurl_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_hcurl_32_loc : f_apply_hcurl_32_loc);
      }
      else if (test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = assemble_q_data ? f_build_hcurlhdiv_32 : f_apply_hcurlhdiv_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(
            assemble_q_data ? f_build_hcurlhdiv_32_loc : f_apply_hcurlhdiv_32_loc);
      }
      else
      {
        MFEM_ABORT("Invalid test map type for MixedVectorGradientIntegrator!");
      }
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << dim << ", " << space_dim << ") for MixedVectorGradientIntegrator!");
  }
  info.trial_ops = EvalMode::Grad;
  info.test_ops = EvalMode::Interp;

  // Set up the coefficient and assemble.
  auto ctx = PopulateCoefficientContext(space_dim, Q, transpose);
  AssembleCeedOperator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed,
                       trial_restr, test_restr, trial_basis, test_basis, geom_data,
                       geom_data_restr, op);
}

void MixedVectorWeakDivergenceIntegrator::Assemble(
    Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
    CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
    CeedElemRestriction geom_data_restr, CeedOperator *op) const
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
      "MixedVectorWeakDivergenceIntegrator requires test and trial spaces with a single "
      "component!");
  switch (10 * space_dim + dim)
  {
    case 22:
      info.apply_qf = assemble_q_data ? f_build_hcurl_22 : f_apply_hcurl_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_q_data ? f_build_hcurl_22_loc : f_apply_hcurl_22_loc);
      break;
    case 33:
      info.apply_qf = assemble_q_data ? f_build_hcurl_33 : f_apply_hcurl_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_q_data ? f_build_hcurl_33_loc : f_apply_hcurl_33_loc);
      break;
    case 21:
      info.apply_qf = assemble_q_data ? f_build_hcurl_21 : f_apply_hcurl_21;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_q_data ? f_build_hcurl_21_loc : f_apply_hcurl_21_loc);
      break;
    case 32:
      info.apply_qf = assemble_q_data ? f_build_hcurl_32 : f_apply_hcurl_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(
          assemble_q_data ? f_build_hcurl_32_loc : f_apply_hcurl_32_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << dim << ", " << space_dim
                 << ") for MixedVectorWeakDivergenceIntegrator!");
  }
  info.trial_ops = EvalMode::Interp;
  info.test_ops = EvalMode::Grad;

  // Set up the coefficient and assemble.
  auto ctx = PopulateCoefficientContext(space_dim, Q, transpose, -1.0);
  AssembleCeedOperator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed,
                       trial_restr, test_restr, trial_basis, test_basis, geom_data,
                       geom_data_restr, op);
}

}  // namespace palace
