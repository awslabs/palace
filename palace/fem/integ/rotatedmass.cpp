// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/hcurlrot_qf.h"

namespace palace
{

using namespace ceed;

void RotatedVectorFEMassIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
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
  MFEM_VERIFY(dim == 2 && space_dim == 3,
              "RotatedVectorFEMassIntegrator is only available on boundaries in 3D!");
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_num_comp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_num_comp));
  MFEM_VERIFY(trial_num_comp == test_num_comp && trial_num_comp == 1,
              "RotatedVectorFEMassIntegrator requires test and trial spaces with a single "
              "component!");
  info.apply_qf = assemble_q_data ? f_build_hcurlrot_32 : f_apply_hcurlrot_32;
  info.apply_qf_path = PalaceQFunctionRelativePath(
      assemble_q_data ? f_build_hcurlrot_32_loc : f_apply_hcurlrot_32_loc);

  // The gradient of an H1 function maps like an H(curl) function.
  if (trial_map_type == mfem::FiniteElement::H_CURL &&
      test_map_type == mfem::FiniteElement::H_CURL)
  {
    info.trial_ops = EvalMode::Interp;
    info.test_ops = EvalMode::Interp;
  }
  else if (trial_map_type == mfem::FiniteElement::VALUE &&
           test_map_type == mfem::FiniteElement::VALUE)
  {
    info.trial_ops = EvalMode::Grad;
    info.test_ops = EvalMode::Grad;
  }
  else
  {
    MFEM_ABORT("Invalid trial/test element map type for RotatedVectorFEMassIntegrator!");
  }

  // Set up the coefficient and assemble.
  auto ctx = PopulateCoefficientContext(space_dim, Q, transpose);
  AssembleCeedOperator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed,
                       trial_restr, test_restr, trial_basis, test_basis, geom_data,
                       geom_data_restr, op);
}

}  // namespace palace
