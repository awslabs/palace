// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include <vector>
#include <mfem.hpp>
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/grad_qf.h"

namespace palace
{

struct GradientIntegratorInfo : public ceed::IntegratorInfo
{
  GradContext ctx;
};

namespace
{

GradientIntegratorInfo
InitializeIntegratorInfo(const mfem::ParFiniteElementSpace &trial_fespace,
                         const mfem::ParFiniteElementSpace &test_fespace,
                         const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                         bool use_bdr, mfem::Coefficient *Q, mfem::VectorCoefficient *VQ,
                         mfem::MatrixCoefficient *MQ,
                         std::vector<ceed::QuadratureCoefficient> &coeff)
{
  GradientIntegratorInfo info = {{0}};

  mfem::ParMesh &mesh = *trial_fespace.GetParMesh();
  info.ctx.dim = mesh.Dimension() - use_bdr;
  info.ctx.space_dim = mesh.SpaceDimension();
  MFEM_VERIFY(trial_fespace.GetVDim() == 1 && test_fespace.GetVDim() == info.ctx.space_dim,
              "libCEED interface for GradientIntegrator requires trial space vdim == 1 and "
              "test space vdim == space dimension!");

  info.trial_op = ceed::EvalMode::Grad;
  info.test_op = ceed::EvalMode::Interp;
  info.qdata_size = info.ctx.space_dim * info.ctx.dim;

  mfem::ConstantCoefficient *const_coeff = dynamic_cast<mfem::ConstantCoefficient *>(Q);
  if (const_coeff || !(Q || VQ || MQ))
  {
    info.ctx.coeff = const_coeff ? const_coeff->constant : 1.0;

    info.build_qf = f_build_grad_const_scalar;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_grad_const_scalar_loc);
  }
  else if (Q)
  {
    ceed::InitCoefficient(*Q, mesh, ir, indices, use_bdr, coeff.emplace_back());

    info.build_qf = f_build_grad_quad_scalar;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_grad_quad_scalar_loc);
  }
  else if (VQ)
  {
    MFEM_VERIFY(VQ->GetVDim() == info.ctx.space_dim,
                "Invalid vector coefficient dimension for GradientIntegrator integrator!");
    ceed::InitCoefficient(*VQ, mesh, ir, indices, use_bdr, coeff.emplace_back());

    info.build_qf = f_build_grad_quad_vector;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_grad_quad_vector_loc);
  }
  else if (MQ)
  {
    MFEM_VERIFY(MQ->GetVDim() == info.ctx.space_dim,
                "Invalid matrix coefficient dimension for GradientIntegrator integrator!");
    ceed::InitCoefficient(*MQ, mesh, ir, indices, use_bdr, coeff.emplace_back());

    info.build_qf = f_build_grad_quad_matrix;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_grad_quad_matrix_loc);
  }

  info.apply_qf = f_apply_grad;
  info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_grad_loc);

  return info;
}

}  // namespace

void GradientIntegrator::Assemble(const mfem::ParFiniteElementSpace &trial_fespace,
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

void GradientIntegrator::AssembleBoundary(const mfem::ParFiniteElementSpace &trial_fespace,
                                          const mfem::ParFiniteElementSpace &test_fespace,
                                          const mfem::IntegrationRule &ir,
                                          const std::vector<int> &indices, Ceed ceed,
                                          CeedOperator *op, CeedOperator *op_t)
{
  constexpr bool use_bdr = true;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, test_fespace, ir, indices,
                                             use_bdr, Q, VQ, MQ, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

}  // namespace palace
