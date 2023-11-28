// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/integrator.hpp"
#include "fem/libceed/utils.hpp"

#include "fem/qfunctions/hcurl_qf.h"
#include "fem/qfunctions/hcurlhdiv_qf.h"
#include "fem/qfunctions/hdiv_qf.h"

namespace palace
{

namespace
{

struct VectorFEMassIntegratorInfo : public ceed::IntegratorInfo
{
  bool ctx;  // XX TODO WIP COEFFICIENTS
};

}  // namespace

void VectorFEMassIntegrator::Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                                      CeedElemRestriction trial_restr,
                                      CeedElemRestriction test_restr, CeedBasis trial_basis,
                                      CeedBasis test_basis, CeedOperator *op)
{
  VectorFEMassIntegratorInfo info;

  // Set up geometry factor quadrature data.
  MFEM_VERIFY(geom_data->wdetJ_vec && geom_data->wdetJ_restr,
              "Missing geometry factor quadrature data for VectorFEMassIntegrator!");
  info.geom_info = ceed::GeomFactorInfo::Determinant;
  if (trial_map_type == mfem::FiniteElement::H_CURL ||
      test_map_type == mfem::FiniteElement::H_CURL)
  {
    MFEM_VERIFY(geom_data->adjJt_vec && geom_data->adjJt_restr,
                "Missing geometry factor quadrature data for VectorFEMassIntegrator!");
    info.geom_info |= ceed::GeomFactorInfo::Adjugate;
  }
  if (trial_map_type == mfem::FiniteElement::H_DIV ||
      test_map_type == mfem::FiniteElement::H_DIV)
  {
    MFEM_VERIFY(geom_data->J_vec && geom_data->J_restr,
                "Missing geometry factor quadrature data for VectorFEMassIntegrator!");
    info.geom_info |= ceed::GeomFactorInfo::Jacobian;
  }

  // Set up QFunctions.
  CeedInt trial_ncomp, test_ncomp;
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(trial_basis, &trial_ncomp));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(test_basis, &test_ncomp));
  MFEM_VERIFY(
      trial_ncomp == test_ncomp && trial_ncomp == 1,
      "VectorFEMassIntegrator requires test and trial spaces with a single component!");
  switch (10 * geom_data->space_dim + geom_data->dim)
  {
    case 22:
      if (trial_map_type == mfem::FiniteElement::H_CURL &&
          test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = f_apply_hcurl_22;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurl_22_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_DIV &&
               test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = f_apply_hdiv_22;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdiv_22_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_CURL &&
               test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = f_apply_hcurlhdiv_22;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlhdiv_22_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_DIV &&
               test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = f_apply_hdivhcurl_22;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdivhcurl_22_loc);
      }
      else
      {
        MFEM_ABORT("Invalid trial/test element map type for VectorFEMassIntegrator!");
      }
      break;
    case 33:
      if (trial_map_type == mfem::FiniteElement::H_CURL &&
          test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = f_apply_hcurl_33;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurl_33_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_DIV &&
               test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = f_apply_hdiv_33;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdiv_33_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_CURL &&
               test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = f_apply_hcurlhdiv_33;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlhdiv_33_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_DIV &&
               test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = f_apply_hdivhcurl_33;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdivhcurl_33_loc);
      }
      else
      {
        MFEM_ABORT("Invalid trial/test element map type for VectorFEMassIntegrator!");
      }
      break;
    case 21:
      if (trial_map_type == mfem::FiniteElement::H_CURL &&
          test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = f_apply_hcurl_21;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurl_21_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_DIV &&
               test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = f_apply_hdiv_21;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdiv_21_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_CURL &&
               test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = f_apply_hcurlhdiv_21;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlhdiv_21_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_DIV &&
               test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = f_apply_hdivhcurl_21;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdivhcurl_21_loc);
      }
      else
      {
        MFEM_ABORT("Invalid trial/test element map type for VectorFEMassIntegrator!");
      }
      break;
    case 32:
      if (trial_map_type == mfem::FiniteElement::H_CURL &&
          test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = f_apply_hcurl_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurl_32_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_DIV &&
               test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = f_apply_hdiv_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdiv_32_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_CURL &&
               test_map_type == mfem::FiniteElement::H_DIV)
      {
        info.apply_qf = f_apply_hcurlhdiv_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlhdiv_32_loc);
      }
      else if (trial_map_type == mfem::FiniteElement::H_DIV &&
               test_map_type == mfem::FiniteElement::H_CURL)
      {
        info.apply_qf = f_apply_hdivhcurl_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdivhcurl_32_loc);
      }
      else
      {
        MFEM_ABORT("Invalid trial/test element map type for VectorFEMassIntegrator!");
      }
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << geom_data->dim << ", " << geom_data->space_dim
                 << ") for VectorFEMassIntegrator!");
  }
  info.trial_ops = ceed::EvalMode::Interp;
  info.test_ops = ceed::EvalMode::Interp;

  ceed::AssembleCeedOperator(info, geom_data, ceed, trial_restr, test_restr, trial_basis,
                             test_basis, op);
}

}  // namespace palace
