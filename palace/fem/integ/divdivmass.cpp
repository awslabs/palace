// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include <vector>
#include <mfem.hpp>
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/divdivmass_qf.h"

namespace palace
{

struct DivDivMassIntegratorInfo : public ceed::IntegratorInfo
{
  DivDivMassContext ctx;
};

namespace
{

DivDivMassIntegratorInfo
InitializeIntegratorInfo(const mfem::FiniteElementSpace &fespace,
                         const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                         bool use_bdr, mfem::Coefficient *Qd, mfem::Coefficient *Qm,
                         mfem::VectorCoefficient *VQm, mfem::MatrixCoefficient *MQm,
                         std::vector<ceed::QuadratureCoefficient> &coeff)
{
  MFEM_VERIFY(fespace.GetVDim() == 1,
              "libCEED interface for DivDivMassIntegrator does not support vdim > 1!");

  DivDivMassIntegratorInfo info = {0};

  mfem::Mesh &mesh = *fespace.GetMesh();
  info.ctx.dim = mesh.Dimension() - use_bdr;
  info.ctx.space_dim = mesh.SpaceDimension();

  info.trial_op = ceed::EvalMode::InterpAndDiv;
  info.test_op = ceed::EvalMode::InterpAndDiv;
  info.qdata_size = 1 + (info.ctx.dim * (info.ctx.dim + 1)) / 2;

  MFEM_VERIFY(Qd && (Qm || VQm || MQm), "libCEED DivDivMassIntegrator requires both a "
                                        "div-div and a mass integrator coefficient!");
  coeff.emplace_back();
  ceed::InitCoefficient(*Qd, mesh, ir, indices, use_bdr, coeff.back());
  if (Qm)
  {
    coeff.emplace_back();
    ceed::InitCoefficient(*Qm, mesh, ir, indices, use_bdr, coeff.back());

    info.build_qf = f_build_divdiv_mass_quad_scalar;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_divdiv_mass_quad_scalar_loc);
  }
  else if (VQm)
  {
    MFEM_VERIFY(VQm->GetVDim() == info.ctx.space_dim,
                "Invalid vector coefficient dimension for DivDivMassIntegrator!");
    coeff.emplace_back();
    ceed::InitCoefficient(*VQm, mesh, ir, indices, use_bdr, coeff.back());

    info.build_qf = f_build_divdiv_mass_quad_vector;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_divdiv_mass_quad_vector_loc);
  }
  else if (MQm)
  {
    MFEM_VERIFY(MQm->GetVDim() == info.ctx.space_dim,
                "Invalid matrix coefficient dimension for DivDivMassIntegrator!");
    coeff.emplace_back();
    ceed::InitCoefficient(*MQm, mesh, ir, indices, use_bdr, coeff.back());

    info.build_qf = f_build_divdiv_mass_quad_matrix;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_divdiv_mass_quad_matrix_loc);
  }

  info.apply_qf = f_apply_divdiv_mass;
  info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_divdiv_mass_loc);

  return info;
}

}  // namespace

void DivDivMassIntegrator::Assemble(const mfem::FiniteElementSpace &trial_fespace,
                                    const mfem::FiniteElementSpace &test_fespace,
                                    const mfem::IntegrationRule &ir,
                                    const std::vector<int> &indices, Ceed ceed,
                                    CeedOperator *op, CeedOperator *op_t)
{
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "DivDivMassIntegrator requires the same test and trial spaces!");
  constexpr bool use_bdr = false;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, ir, indices, use_bdr, Qd, Qm,
                                             VQm, MQm, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

void DivDivMassIntegrator::AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                                            const mfem::FiniteElementSpace &test_fespace,
                                            const mfem::IntegrationRule &ir,
                                            const std::vector<int> &indices, Ceed ceed,
                                            CeedOperator *op, CeedOperator *op_t)
{
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "DivDivMassIntegrator requires the same test and trial spaces!");
  constexpr bool use_bdr = true;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, ir, indices, use_bdr, Qd, Qm,
                                             VQm, MQm, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

}  // namespace palace
