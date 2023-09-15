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

struct MixedVectorCurlIntegratorInfo : public ceed::IntegratorInfo
{
  CurlCurlContext ctx;
};

namespace
{

MixedVectorCurlIntegratorInfo
InitializeIntegratorInfo(const mfem::FiniteElementSpace &trial_fespace,
                         const mfem::FiniteElementSpace &test_fespace,
                         const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                         bool use_bdr, mfem::Coefficient *Q, mfem::VectorCoefficient *VQ,
                         mfem::MatrixCoefficient *MQ,
                         std::vector<ceed::QuadratureCoefficient> &coeff)
{
  MFEM_VERIFY(
      trial_fespace.GetVDim() == 1 && test_fespace.GetVDim() == 1,
      "libCEED interface for MixedVectorCurlIntegrator/MixedVectorWeakCurlIntegrator does "
      "not support vdim > 1!");

  MixedVectorCurlIntegratorInfo info = {0};

  mfem::Mesh &mesh = *trial_fespace.GetMesh();
  info.ctx.dim = mesh.Dimension() - use_bdr;
  info.ctx.space_dim = mesh.SpaceDimension();
  info.ctx.curl_dim = (info.ctx.dim < 3) ? 1 : info.ctx.dim;

  int trial_map_type = trial_fespace.FEColl()->GetMapType(info.ctx.dim);
  int trial_deriv_map_type = trial_fespace.FEColl()->GetDerivMapType(info.ctx.dim);
  int test_map_type = test_fespace.FEColl()->GetMapType(info.ctx.dim);
  int test_deriv_map_type = test_fespace.FEColl()->GetDerivMapType(info.ctx.dim);
  MFEM_VERIFY(
      (trial_map_type == mfem::FiniteElement::H_DIV &&
       test_deriv_map_type == mfem::FiniteElement::H_DIV) ||
          (trial_deriv_map_type == mfem::FiniteElement::H_DIV &&
           test_map_type == mfem::FiniteElement::H_DIV),
      "libCEED interface for MixedVectorCurlIntegrator/MixedVectorWeakCurlIntegrator "
      "requires mixed H(curl) and H(div) FE spaces!");

  info.trial_op = (trial_map_type == mfem::FiniteElement::H_DIV) ? ceed::EvalMode::Interp
                                                                 : ceed::EvalMode::Curl;
  info.test_op = (test_map_type == mfem::FiniteElement::H_DIV) ? ceed::EvalMode::Interp
                                                               : ceed::EvalMode::Curl;
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
                "Invalid vector coefficient dimension for "
                "MixedVectorCurlIntegrator/MixedVectorWeakCurlIntegrator integrator!");
    coeff.emplace_back();
    ceed::InitCoefficient(*VQ, mesh, ir, indices, use_bdr, coeff.back());

    info.build_qf = f_build_curlcurl_quad_vector;
    info.build_qf_path = PalaceQFunctionRelativePath(f_build_curlcurl_quad_vector_loc);
  }
  else if (MQ)
  {
    MFEM_VERIFY(MQ->GetVDim() == info.ctx.curl_dim,
                "Invalid matrix coefficient dimension for "
                "MixedVectorCurlIntegrator/MixedVectorWeakCurlIntegrator integrator!");
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

void MixedVectorCurlIntegrator::Assemble(const mfem::FiniteElementSpace &trial_fespace,
                                         const mfem::FiniteElementSpace &test_fespace,
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

void MixedVectorCurlIntegrator::AssembleBoundary(
    const mfem::FiniteElementSpace &trial_fespace,
    const mfem::FiniteElementSpace &test_fespace, const mfem::IntegrationRule &ir,
    const std::vector<int> &indices, Ceed ceed, CeedOperator *op, CeedOperator *op_t)
{
  constexpr bool use_bdr = true;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, test_fespace, ir, indices,
                                             use_bdr, Q, VQ, MQ, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

void MixedVectorWeakCurlIntegrator::Assemble(const mfem::FiniteElementSpace &trial_fespace,
                                             const mfem::FiniteElementSpace &test_fespace,
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

void MixedVectorWeakCurlIntegrator::AssembleBoundary(
    const mfem::FiniteElementSpace &trial_fespace,
    const mfem::FiniteElementSpace &test_fespace, const mfem::IntegrationRule &ir,
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
