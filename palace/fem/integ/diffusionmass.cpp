// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include <vector>
#include <mfem.hpp>
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/diffusionmass_qf.h"

namespace palace
{

struct DiffusionMassIntegratorInfo : public ceed::IntegratorInfo
{
  DiffusionMassContext ctx;
};

namespace
{

DiffusionMassIntegratorInfo
InitializeIntegratorInfo(const mfem::FiniteElementSpace &fespace,
                         const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                         bool use_bdr, mfem::Coefficient *Qd, mfem::VectorCoefficient *VQd,
                         mfem::MatrixCoefficient *MQd, mfem::Coefficient *Qm,
                         std::vector<ceed::QuadratureCoefficient> &coeff)
{
  MFEM_VERIFY(fespace.GetVDim() == 1,
              "libCEED interface for DiffusionMassIntegrator does not support vdim > 1!");

  DiffusionMassIntegratorInfo info = {0};

  mfem::Mesh &mesh = *fespace.GetMesh();
  info.ctx.dim = mesh.Dimension() - use_bdr;
  info.ctx.space_dim = mesh.SpaceDimension();

  info.trial_op = ceed::EvalMode::InterpAndGrad;
  info.test_op = ceed::EvalMode::InterpAndGrad;
  info.qdata_size = (info.ctx.dim * (info.ctx.dim + 1)) / 2 + 1;

  MFEM_VERIFY((Qd || VQd || MQd) && Qm, "libCEED DiffusionMassIntegrator requires both a "
                                        "diffusion and a mass integrator coefficient!");
  if (Qd)
  {
    coeff.emplace_back();
    ceed::InitCoefficient(*Qd, mesh, ir, indices, use_bdr, coeff.back());

    info.build_qf = f_build_diff_mass_quad_scalar;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_diff_mass_quad_scalar_loc);
  }
  else if (VQd)
  {
    MFEM_VERIFY(VQd->GetVDim() == info.ctx.space_dim,
                "Invalid vector coefficient dimension for DiffusionMassIntegrator!");
    coeff.emplace_back();
    ceed::InitCoefficient(*VQd, mesh, ir, indices, use_bdr, coeff.back());

    info.build_qf = f_build_diff_mass_quad_vector;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_diff_mass_quad_vector_loc);
  }
  else if (MQd)
  {
    MFEM_VERIFY(MQd->GetVDim() == info.ctx.space_dim,
                "Invalid matrix coefficient dimension for DiffusionMassIntegrator!");
    coeff.emplace_back();
    ceed::InitCoefficient(*MQd, mesh, ir, indices, use_bdr, coeff.back());

    info.build_qf = f_build_diff_mass_quad_matrix;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_diff_mass_quad_matrix_loc);
  }
  coeff.emplace_back();
  ceed::InitCoefficient(*Qm, mesh, ir, indices, use_bdr, coeff.back());

  info.apply_qf = f_apply_diff_mass;
  info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_diff_mass_loc);

  return info;
}

}  // namespace

void DiffusionMassIntegrator::Assemble(const mfem::FiniteElementSpace &trial_fespace,
                                       const mfem::FiniteElementSpace &test_fespace,
                                       const mfem::IntegrationRule &ir,
                                       const std::vector<int> &indices, Ceed ceed,
                                       CeedOperator *op, CeedOperator *op_t)
{
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "DiffusionMassIntegrator requires the same test and trial spaces!");
  constexpr bool use_bdr = false;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, ir, indices, use_bdr, Qd, VQd,
                                             MQd, Qm, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

void DiffusionMassIntegrator::AssembleBoundary(
    const mfem::FiniteElementSpace &trial_fespace,
    const mfem::FiniteElementSpace &test_fespace, const mfem::IntegrationRule &ir,
    const std::vector<int> &indices, Ceed ceed, CeedOperator *op, CeedOperator *op_t)
{
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "DiffusionMassIntegrator requires the same test and trial spaces!");
  constexpr bool use_bdr = true;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, ir, indices, use_bdr, Qd, VQd,
                                             MQd, Qm, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

}  // namespace palace
