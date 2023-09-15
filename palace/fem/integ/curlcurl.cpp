// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include <vector>
#include <mfem.hpp>
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/curlcurl_qf.h"

namespace palace
{

struct CurlCurlIntegratorInfo : public ceed::IntegratorInfo
{
  CurlCurlContext ctx;
};

namespace
{

CurlCurlIntegratorInfo
InitializeIntegratorInfo(const mfem::FiniteElementSpace &fespace,
                         const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                         bool use_bdr, mfem::Coefficient *Q, mfem::VectorCoefficient *VQ,
                         mfem::MatrixCoefficient *MQ,
                         std::vector<ceed::QuadratureCoefficient> &coeff)
{
  MFEM_VERIFY(fespace.GetVDim() == 1,
              "libCEED interface for CurlCurlIntegrator does not support vdim > 1!");

  CurlCurlIntegratorInfo info = {0};

  mfem::Mesh &mesh = *fespace.GetMesh();
  info.ctx.dim = mesh.Dimension() - use_bdr;
  info.ctx.space_dim = mesh.SpaceDimension();
  info.ctx.curl_dim = (info.ctx.dim < 3) ? 1 : info.ctx.dim;

  info.trial_op = ceed::EvalMode::Curl;
  info.test_op = ceed::EvalMode::Curl;
  info.qdata_size = (info.ctx.curl_dim * (info.ctx.curl_dim + 1)) / 2;

  mfem::ConstantCoefficient *const_coeff = dynamic_cast<mfem::ConstantCoefficient *>(Q);
  if (const_coeff || !(Q || VQ || MQ))
  {
    info.ctx.coeff = const_coeff ? const_coeff->constant : 1.0;

    info.build_qf = f_build_curlcurl_const_scalar;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_curlcurl_const_scalar_loc);
  }
  else if (Q)
  {
    coeff.emplace_back();
    ceed::InitCoefficient(*Q, mesh, ir, indices, use_bdr, coeff.back());

    info.build_qf = f_build_curlcurl_quad_scalar;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_curlcurl_quad_scalar_loc);
  }
  else if (VQ)
  {
    MFEM_VERIFY(VQ->GetVDim() == info.ctx.curl_dim,
                "Invalid vector coefficient dimension for CurlCurlIntegrator!");
    coeff.emplace_back();
    ceed::InitCoefficient(*VQ, mesh, ir, indices, use_bdr, coeff.back());

    info.build_qf = f_build_curlcurl_quad_vector;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_curlcurl_quad_vector_loc);
  }
  else if (MQ)
  {
    MFEM_VERIFY(MQ->GetVDim() == info.ctx.curl_dim,
                "Invalid matrix coefficient dimension for CurlCurlIntegrator!");
    coeff.emplace_back();
    ceed::InitCoefficient(*MQ, mesh, ir, indices, use_bdr, coeff.back());

    info.build_qf = f_build_curlcurl_quad_matrix;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_curlcurl_quad_matrix_loc);
  }

  info.apply_qf = f_apply_curlcurl;
  info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_curlcurl_loc);

  return info;
}

}  // namespace

void CurlCurlIntegrator::Assemble(const mfem::FiniteElementSpace &trial_fespace,
                                  const mfem::FiniteElementSpace &test_fespace,
                                  const mfem::IntegrationRule &ir,
                                  const std::vector<int> &indices, Ceed ceed,
                                  CeedOperator *op, CeedOperator *op_t)
{
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "CurlCurlIntegrator requires the same test and trial spaces!");
  constexpr bool use_bdr = false;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info =
      InitializeIntegratorInfo(trial_fespace, ir, indices, use_bdr, Q, VQ, MQ, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

void CurlCurlIntegrator::AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                                          const mfem::FiniteElementSpace &test_fespace,
                                          const mfem::IntegrationRule &ir,
                                          const std::vector<int> &indices, Ceed ceed,
                                          CeedOperator *op, CeedOperator *op_t)
{
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "CurlCurlIntegrator requires the same test and trial spaces!");
  constexpr bool use_bdr = true;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info =
      InitializeIntegratorInfo(trial_fespace, ir, indices, use_bdr, Q, VQ, MQ, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

}  // namespace palace
