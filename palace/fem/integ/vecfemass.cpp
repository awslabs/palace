// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include <vector>
#include <mfem.hpp>
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/hcurl_qf.h"
#include "fem/qfunctions/hcurlhdiv_qf.h"
#include "fem/qfunctions/hdiv_qf.h"

namespace palace
{

struct VectorFEMassIntegratorInfo : public ceed::IntegratorInfo
{
  VectorFEMassContext ctx;
};

namespace
{

VectorFEMassIntegratorInfo
InitializeIntegratorInfo(const mfem::ParFiniteElementSpace &trial_fespace,
                         const mfem::ParFiniteElementSpace &test_fespace,
                         const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                         bool use_bdr, mfem::Coefficient *Q, mfem::VectorCoefficient *VQ,
                         mfem::MatrixCoefficient *MQ,
                         std::vector<ceed::QuadratureCoefficient> &coeff)
{
  MFEM_VERIFY(trial_fespace.GetVDim() == 1 && test_fespace.GetVDim() == 1,
              "libCEED interface for VectorFEMassIntegrator does not support vdim > 1!");

  VectorFEMassIntegratorInfo info = {{0}};

  mfem::ParMesh &mesh = *trial_fespace.GetParMesh();
  info.ctx.dim = mesh.Dimension() - use_bdr;
  info.ctx.space_dim = mesh.SpaceDimension();

  int trial_map_type = trial_fespace.FEColl()->GetMapType(info.ctx.dim);
  int test_map_type = test_fespace.FEColl()->GetMapType(info.ctx.dim);
  MFEM_VERIFY((trial_map_type == mfem::FiniteElement::H_CURL ||
               trial_map_type == mfem::FiniteElement::H_DIV) &&
                  (test_map_type == mfem::FiniteElement::H_CURL ||
                   test_map_type == mfem::FiniteElement::H_DIV),
              "VectorFEMassIntegrator requires H(div) or H(curl) FE spaces!");

  info.trial_op = ceed::EvalMode::Interp;
  info.test_op = ceed::EvalMode::Interp;
  if (trial_map_type != test_map_type)
  {
    // Quadrature data is nonsymmetric in this case.
    info.qdata_size = info.ctx.dim * info.ctx.dim;
    info.ctx.sym = false;
  }
  else
  {
    info.qdata_size = (info.ctx.dim * (info.ctx.dim + 1)) / 2;
    info.ctx.sym = true;
  }

  mfem::ConstantCoefficient *const_coeff = dynamic_cast<mfem::ConstantCoefficient *>(Q);
  if (const_coeff || !(Q || VQ || MQ))
  {
    info.ctx.coeff = const_coeff ? const_coeff->constant : 1.0;

    if (trial_map_type == mfem::FiniteElement::H_CURL &&
        test_map_type == mfem::FiniteElement::H_CURL)
    {
      info.build_qf = f_build_hcurl_const_scalar;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hcurl_const_scalar_loc);
    }
    else if (trial_map_type == mfem::FiniteElement::H_DIV &&
             test_map_type == mfem::FiniteElement::H_DIV)
    {
      info.build_qf = f_build_hdiv_const_scalar;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hdiv_const_scalar_loc);
    }
    else if (trial_map_type == mfem::FiniteElement::H_CURL &&
             test_map_type == mfem::FiniteElement::H_DIV)
    {
      info.build_qf = f_build_hcurlhdiv_const_scalar;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hcurlhdiv_const_scalar_loc);
    }
    else  // trial_map_type == mfem::FiniteElement::H_DIV && test_map_type ==
          // mfem::FiniteElement::H_CURL
    {
      info.build_qf = f_build_hdivhcurl_const_scalar;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hdivhcurl_const_scalar_loc);
    }
  }
  else if (Q)
  {
    ceed::InitCoefficient(*Q, mesh, ir, indices, use_bdr, coeff.emplace_back());

    if (trial_map_type == mfem::FiniteElement::H_CURL &&
        test_map_type == mfem::FiniteElement::H_CURL)
    {
      info.build_qf = f_build_hcurl_quad_scalar;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hcurl_quad_scalar_loc);
    }
    else if (trial_map_type == mfem::FiniteElement::H_DIV &&
             test_map_type == mfem::FiniteElement::H_DIV)
    {
      info.build_qf = f_build_hdiv_quad_scalar;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hdiv_quad_scalar_loc);
    }
    else if (trial_map_type == mfem::FiniteElement::H_CURL &&
             test_map_type == mfem::FiniteElement::H_DIV)
    {
      info.build_qf = f_build_hcurlhdiv_quad_scalar;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hcurlhdiv_quad_scalar_loc);
    }
    else  // trial_map_type == mfem::FiniteElement::H_DIV && test_map_type ==
          // mfem::FiniteElement::H_CURL
    {
      info.build_qf = f_build_hdivhcurl_quad_scalar;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hdivhcurl_quad_scalar_loc);
    }
  }
  else if (VQ)
  {
    MFEM_VERIFY(VQ->GetVDim() == info.ctx.space_dim,
                "Invalid vector coefficient dimension for VectorFEMassIntegrator!");
    ceed::InitCoefficient(*VQ, mesh, ir, indices, use_bdr, coeff.emplace_back());

    if (trial_map_type == mfem::FiniteElement::H_CURL &&
        test_map_type == mfem::FiniteElement::H_CURL)
    {
      info.build_qf = f_build_hcurl_quad_vector;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hcurl_quad_vector_loc);
    }
    else if (trial_map_type == mfem::FiniteElement::H_DIV &&
             test_map_type == mfem::FiniteElement::H_DIV)
    {
      info.build_qf = f_build_hdiv_quad_vector;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hdiv_quad_vector_loc);
    }
    else if (trial_map_type == mfem::FiniteElement::H_CURL &&
             test_map_type == mfem::FiniteElement::H_DIV)
    {
      info.build_qf = f_build_hcurlhdiv_quad_vector;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hcurlhdiv_quad_vector_loc);
    }
    else  // trial_map_type == mfem::FiniteElement::H_DIV && test_map_type ==
          // mfem::FiniteElement::H_CURL
    {
      info.build_qf = f_build_hdivhcurl_quad_vector;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hdivhcurl_quad_vector_loc);
    }
  }
  else if (MQ)
  {
    MFEM_VERIFY(MQ->GetVDim() == info.ctx.space_dim,
                "Invalid matrix coefficient dimension for VectorFEMassIntegrator!");
    ceed::InitCoefficient(*MQ, mesh, ir, indices, use_bdr, coeff.emplace_back());

    if (trial_map_type == mfem::FiniteElement::H_CURL &&
        test_map_type == mfem::FiniteElement::H_CURL)
    {
      info.build_qf = f_build_hcurl_quad_matrix;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hcurl_quad_matrix_loc);
    }
    else if (trial_map_type == mfem::FiniteElement::H_DIV &&
             test_map_type == mfem::FiniteElement::H_DIV)
    {
      info.build_qf = f_build_hdiv_quad_matrix;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hdiv_quad_matrix_loc);
    }
    else if (trial_map_type == mfem::FiniteElement::H_CURL &&
             test_map_type == mfem::FiniteElement::H_DIV)
    {
      info.build_qf = f_build_hcurlhdiv_quad_matrix;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hcurlhdiv_quad_matrix_loc);
    }
    else  // trial_map_type == mfem::FiniteElement::H_DIV && test_map_type ==
          // mfem::FiniteElement::H_CURL
    {
      info.build_qf = f_build_hdivhcurl_quad_matrix;
      info.build_qf_path = PalaceQFunctionRelativePath(f_build_hdivhcurl_quad_matrix_loc);
    }
  }

  info.apply_qf = f_apply_vecfemass;
  info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_vecfemass_loc);

  return info;
}

}  // namespace

void VectorFEMassIntegrator::Assemble(const mfem::ParFiniteElementSpace &trial_fespace,
                                      const mfem::ParFiniteElementSpace &test_fespace,
                                      const mfem::IntegrationRule &ir,
                                      const std::vector<int> &indices, Ceed ceed,
                                      CeedOperator *op, CeedOperator *op_t)
{
  constexpr bool use_bdr = false;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, test_fespace, ir, indices,
                                             use_bdr, Q, VQ, MQ, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

void VectorFEMassIntegrator::AssembleBoundary(
    const mfem::ParFiniteElementSpace &trial_fespace,
    const mfem::ParFiniteElementSpace &test_fespace, const mfem::IntegrationRule &ir,
    const std::vector<int> &indices, Ceed ceed, CeedOperator *op, CeedOperator *op_t)
{
  constexpr bool use_bdr = true;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, test_fespace, ir, indices,
                                             use_bdr, Q, VQ, MQ, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

}  // namespace palace
