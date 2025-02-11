// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/hcurlhdiv_qf.h"
#include "fem/qfunctions/hdiv_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

using namespace ceed;

void MixedVectorCurlIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
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
  MFEM_VERIFY(dim == 3 && space_dim == 3,
              "MixedVectorCurlIntegrator is only available in 3D!");
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_num_comp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_num_comp));
  MFEM_VERIFY(
      trial_num_comp == test_num_comp && trial_num_comp == 1,
      "MixedVectorCurlIntegrator requires test and trial spaces with a single component!");
  if (test_map_type == mfem::FiniteElement::H_DIV)
  {
    info.apply_qf = assemble_q_data ? f_build_hdiv_33 : f_apply_hdiv_33;
    info.apply_qf_path = PalaceQFunctionRelativePath(assemble_q_data ? f_build_hdiv_33_loc
                                                                     : f_apply_hdiv_33_loc);
  }
  else if (test_map_type == mfem::FiniteElement::H_CURL)
  {
    info.apply_qf = assemble_q_data ? f_build_hdivhcurl_33 : f_apply_hdivhcurl_33;
    info.apply_qf_path = PalaceQFunctionRelativePath(
        assemble_q_data ? f_build_hdivhcurl_33_loc : f_apply_hdivhcurl_33_loc);
  }
  else
  {
    MFEM_ABORT("Invalid trial/test element map type for MixedVectorCurlIntegrator!");
  }
  info.trial_ops = EvalMode::Curl;
  info.test_ops = EvalMode::Interp;

  // Set up the coefficient and assemble.
  auto ctx = PopulateCoefficientContext(space_dim, Q, transpose);
  AssembleCeedOperator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed,
                       trial_restr, test_restr, trial_basis, test_basis, geom_data,
                       geom_data_restr, op);
}

void MixedVectorWeakCurlIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
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
  MFEM_VERIFY(dim == 3 && space_dim == 3,
              "MixedVectorWeakCurlIntegrator is only available in 3D!");
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_num_comp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_num_comp));
  MFEM_VERIFY(trial_num_comp == test_num_comp && trial_num_comp == 1,
              "MixedVectorWeakCurlIntegrator requires test and trial spaces with a single "
              "component!");
  if (trial_map_type == mfem::FiniteElement::H_DIV)
  {
    info.apply_qf = assemble_q_data ? f_build_hdiv_33 : f_apply_hdiv_33;
    info.apply_qf_path = PalaceQFunctionRelativePath(assemble_q_data ? f_build_hdiv_33_loc
                                                                     : f_apply_hdiv_33_loc);
  }
  else if (trial_map_type == mfem::FiniteElement::H_CURL)
  {
    info.apply_qf = assemble_q_data ? f_build_hcurlhdiv_33 : f_apply_hcurlhdiv_33;
    info.apply_qf_path = PalaceQFunctionRelativePath(
        assemble_q_data ? f_build_hcurlhdiv_33_loc : f_apply_hcurlhdiv_33_loc);
  }
  else
  {
    MFEM_ABORT("Invalid trial/test element map type for MixedVectorWeakCurlIntegrator!");
  }
  info.trial_ops = EvalMode::Interp;
  info.test_ops = EvalMode::Curl;

  // Set up the coefficient and assemble.
  auto ctx = PopulateCoefficientContext(space_dim, Q, transpose, -1.0);
  AssembleCeedOperator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed,
                       trial_restr, test_restr, trial_basis, test_basis, geom_data,
                       geom_data_restr, op);
}

}  // namespace palace
