// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/integrator.hpp"
#include "fem/libceed/utils.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/hcurlhdiv_qf.h"
#include "fem/qfunctions/hdiv_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

namespace
{

struct MixedVectorCurlIntegratorInfo : public ceed::IntegratorInfo
{
  bool ctx;  // XX TODO WIP COEFFICIENTS
};

}  // namespace

void MixedVectorCurlIntegrator::Assemble(const ceed::CeedGeomFactorData &geom_data,
                                         Ceed ceed, CeedElemRestriction trial_restr,
                                         CeedElemRestriction test_restr,
                                         CeedBasis trial_basis, CeedBasis test_basis,
                                         CeedOperator *op)
{
  MFEM_VERIFY(geom_data->dim == 3 && geom_data->space_dim == 3,
              "MixedVectorCurlIntegrator is only availble in 3D!");
  MixedVectorCurlIntegratorInfo info;

  // Set up geometry factor quadrature data.
  MFEM_VERIFY(geom_data->wdetJ_vec && geom_data->wdetJ_restr && geom_data->J_vec &&
                  geom_data->J_restr,
              "Missing geometry factor quadrature data for MixedVectorCurlIntegrator!");
  info.geom_info = ceed::GeomFactorInfo::Determinant | ceed::GeomFactorInfo::Jacobian;
  if (test_map_type == mfem::FiniteElement::H_CURL)
  {
    MFEM_VERIFY(geom_data->adjJt_vec && geom_data->adjJt_restr,
                "Missing geometry factor quadrature data for MixedVectorCurlIntegrator!");
    info.geom_info |= ceed::GeomFactorInfo::Adjugate;
  }

  // Set up QFunctions.
  CeedInt trial_ncomp, test_ncomp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_ncomp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_ncomp));
  MFEM_VERIFY(
      trial_ncomp == test_ncomp && trial_ncomp == 1,
      "MixedVectorCurlIntegrator requires test and trial spaces with a single component!");
  if (test_map_type == mfem::FiniteElement::H_DIV)
  {
    info.apply_qf = f_apply_hdiv_33;
    info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdiv_33_loc);
  }
  else if (test_map_type == mfem::FiniteElement::H_CURL)
  {
    info.apply_qf = f_apply_hdivhcurl_33;
    info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdivhcurl_33_loc);
  }
  else
  {
    MFEM_ABORT("Invalid trial/test element map type for MixedVectorCurlIntegrator!");
  }
  info.trial_ops = ceed::EvalMode::Curl;
  info.test_ops = ceed::EvalMode::Interp;

  ceed::AssembleCeedOperator(info, geom_data, ceed, trial_restr, test_restr, trial_basis,
                             test_basis, op);
}

void MixedVectorWeakCurlIntegrator::Assemble(const ceed::CeedGeomFactorData &geom_data,
                                             Ceed ceed, CeedElemRestriction trial_restr,
                                             CeedElemRestriction test_restr,
                                             CeedBasis trial_basis, CeedBasis test_basis,
                                             CeedOperator *op)
{
  MFEM_VERIFY(geom_data->dim == 3 && geom_data->space_dim == 3,
              "MixedVectorWeakCurlIntegrator is only availble in 3D!");
  MixedVectorCurlIntegratorInfo info;

  // Set up geometry factor quadrature data.
  MFEM_VERIFY(geom_data->wdetJ_vec && geom_data->wdetJ_restr && geom_data->J_vec &&
                  geom_data->J_restr,
              "Missing geometry factor quadrature data for MixedVectorWeakCurlIntegrator!");
  info.geom_info = ceed::GeomFactorInfo::Determinant | ceed::GeomFactorInfo::Jacobian;
  if (trial_map_type == mfem::FiniteElement::H_CURL)
  {
    MFEM_VERIFY(
        geom_data->adjJt_vec && geom_data->adjJt_restr,
        "Missing geometry factor quadrature data for MixedVectorWeakCurlIntegrator!");
    info.geom_info |= ceed::GeomFactorInfo::Adjugate;
  }

  // Set up QFunctions.
  CeedInt trial_ncomp, test_ncomp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_ncomp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_ncomp));
  MFEM_VERIFY(trial_ncomp == test_ncomp && trial_ncomp == 1,
              "MixedVectorWeakCurlIntegrator requires test and trial spaces with a single "
              "component!");
  if (trial_map_type == mfem::FiniteElement::H_DIV)
  {
    info.apply_qf = f_apply_hdiv_33;
    info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdiv_33_loc);
  }
  else if (trial_map_type == mfem::FiniteElement::H_CURL)
  {
    info.apply_qf = f_apply_hcurlhdiv_33;
    info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlhdiv_33_loc);
  }
  else
  {
    MFEM_ABORT("Invalid trial/test element map type for MixedVectorWeakCurlIntegrator!");
  }
  info.trial_ops = ceed::EvalMode::Interp;
  info.test_ops = ceed::EvalMode::Curl;

  ceed::AssembleCeedOperator(info, geom_data, ceed, trial_restr, test_restr, trial_basis,
                             test_basis, op);
}

}  // namespace palace
