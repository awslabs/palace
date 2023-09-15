// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include <vector>
#include <mfem.hpp>
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/divdiv_qf.h"

namespace palace
{

struct DivDivIntegratorInfo : public ceed::IntegratorInfo
{
  DivDivContext ctx;
};

namespace
{

DivDivIntegratorInfo
InitializeIntegratorInfo(const mfem::FiniteElementSpace &fespace,
                         const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                         bool use_bdr, mfem::Coefficient *Q,
                         std::vector<ceed::QuadratureCoefficient> &coeff)
{
  MFEM_VERIFY(fespace.GetVDim() == 1,
              "libCEED interface for DivDivIntegrator does not support vdim > 1!");

  DivDivIntegratorInfo info = {0};

  mfem::Mesh &mesh = *fespace.GetMesh();
  info.ctx.dim = mesh.Dimension() - use_bdr;
  info.ctx.space_dim = mesh.SpaceDimension();

  info.trial_op = ceed::EvalMode::Div;
  info.test_op = ceed::EvalMode::Div;
  info.qdata_size = 1;

  mfem::ConstantCoefficient *const_coeff = dynamic_cast<mfem::ConstantCoefficient *>(Q);
  if (const_coeff || !Q)
  {
    info.ctx.coeff = const_coeff ? const_coeff->constant : 1.0;

    info.build_qf = f_build_divdiv_const;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_divdiv_const_loc);
  }
  else if (Q)
  {
    coeff.emplace_back();
    ceed::InitCoefficient(*Q, mesh, ir, indices, use_bdr, coeff.back());

    info.build_qf = f_build_divdiv_quad;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_divdiv_quad_loc);
  }

  info.apply_qf = f_apply_divdiv;
  info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_divdiv_loc);

  return info;
}

}  // namespace

void DivDivIntegrator::Assemble(const mfem::FiniteElementSpace &trial_fespace,
                                const mfem::FiniteElementSpace &test_fespace,
                                const mfem::IntegrationRule &ir,
                                const std::vector<int> &indices, Ceed ceed,
                                CeedOperator *op, CeedOperator *op_t)
{
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "DivDivIntegrator requires the same test and trial spaces!");
  constexpr bool use_bdr = false;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, ir, indices, use_bdr, Q, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

void DivDivIntegrator::AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                                        const mfem::FiniteElementSpace &test_fespace,
                                        const mfem::IntegrationRule &ir,
                                        const std::vector<int> &indices, Ceed ceed,
                                        CeedOperator *op, CeedOperator *op_t)
{
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "DivDivIntegrator requires the same test and trial spaces!");
  constexpr bool use_bdr = true;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, ir, indices, use_bdr, Q, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

}  // namespace palace
