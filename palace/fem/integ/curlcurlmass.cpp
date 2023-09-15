// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include <vector>
#include <mfem.hpp>
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/curlcurlmass_qf.h"

namespace palace
{

struct CurlCurlMassIntegratorInfo : public ceed::IntegratorInfo
{
  CurlCurlMassContext ctx;
};

namespace
{

CurlCurlMassIntegratorInfo
InitializeIntegratorInfo(const mfem::FiniteElementSpace &fespace,
                         const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                         bool use_bdr, mfem::Coefficient *Qc, mfem::VectorCoefficient *VQc,
                         mfem::MatrixCoefficient *MQc, mfem::Coefficient *Qm,
                         mfem::VectorCoefficient *VQm, mfem::MatrixCoefficient *MQm,
                         std::vector<ceed::QuadratureCoefficient> &coeff)
{
  MFEM_VERIFY(fespace.GetVDim() == 1,
              "libCEED interface for CurlCurlMassIntegrator does not support vdim > 1!");

  CurlCurlMassIntegratorInfo info = {0};

  mfem::Mesh &mesh = *fespace.GetMesh();
  info.ctx.dim = mesh.Dimension() - use_bdr;
  info.ctx.space_dim = mesh.SpaceDimension();
  info.ctx.curl_dim = (info.ctx.dim < 3) ? 1 : info.ctx.dim;

  info.trial_op = ceed::EvalMode::InterpAndCurl;
  info.test_op = ceed::EvalMode::InterpAndCurl;
  info.qdata_size = (info.ctx.curl_dim * (info.ctx.curl_dim + 1)) / 2 +
                    (info.ctx.dim * (info.ctx.dim + 1)) / 2;

  MFEM_VERIFY((Qc || VQc || MQc) && (Qm || VQm || MQm),
              "libCEED CurlCurlMassIntegrator requires both a "
              "curl-curl and a mass integrator coefficient!");
  if (Qc)
  {
    coeff.emplace_back();
    ceed::InitCoefficient(*Qc, mesh, ir, indices, use_bdr, coeff.back());

    if (Qm)
    {
      coeff.emplace_back();
      ceed::InitCoefficient(*Qm, mesh, ir, indices, use_bdr, coeff.back());

      info.build_qf = f_build_curlcurl_mass_quad_scalar_scalar;
      info.build_qf_path =
          PalaceQFunctionRelativePath(f_build_curlcurl_mass_quad_scalar_scalar_loc);
    }
    else if (VQm)
    {
      MFEM_VERIFY(VQm->GetVDim() == info.ctx.space_dim,
                  "Invalid vector coefficient dimension for CurlCurlMassIntegrator!");
      coeff.emplace_back();
      ceed::InitCoefficient(*VQm, mesh, ir, indices, use_bdr, coeff.back());

      info.build_qf = f_build_curlcurl_mass_quad_scalar_vector;
      info.build_qf_path =
          PalaceQFunctionRelativePath(f_build_curlcurl_mass_quad_scalar_vector_loc);
    }
    else if (MQm)
    {
      MFEM_VERIFY(MQm->GetVDim() == info.ctx.space_dim,
                  "Invalid matrix coefficient dimension for CurlCurlMassIntegrator!");
      coeff.emplace_back();
      ceed::InitCoefficient(*MQm, mesh, ir, indices, use_bdr, coeff.back());

      info.build_qf = f_build_curlcurl_mass_quad_scalar_matrix;
      info.build_qf_path =
          PalaceQFunctionRelativePath(f_build_curlcurl_mass_quad_scalar_matrix_loc);
    }
  }
  else if (VQc)
  {
    MFEM_VERIFY(VQc->GetVDim() == info.ctx.curl_dim,
                "Invalid vector coefficient dimension for CurlCurlMassIntegrator!");
    coeff.emplace_back();
    ceed::InitCoefficient(*VQc, mesh, ir, indices, use_bdr, coeff.back());

    if (Qm)
    {
      coeff.emplace_back();
      ceed::InitCoefficient(*Qm, mesh, ir, indices, use_bdr, coeff.back());

      info.build_qf = f_build_curlcurl_mass_quad_vector_scalar;
      info.build_qf_path =
          PalaceQFunctionRelativePath(f_build_curlcurl_mass_quad_vector_scalar_loc);
    }
    else if (VQm)
    {
      MFEM_VERIFY(VQm->GetVDim() == info.ctx.space_dim,
                  "Invalid vector coefficient dimension for CurlCurlMassIntegrator!");
      coeff.emplace_back();
      ceed::InitCoefficient(*VQm, mesh, ir, indices, use_bdr, coeff.back());

      info.build_qf = f_build_curlcurl_mass_quad_vector_vector;
      info.build_qf_path =
          PalaceQFunctionRelativePath(f_build_curlcurl_mass_quad_vector_vector_loc);
    }
    else if (MQm)
    {
      MFEM_VERIFY(MQm->GetVDim() == info.ctx.space_dim,
                  "Invalid matrix coefficient dimension for CurlCurlMassIntegrator!");
      coeff.emplace_back();
      ceed::InitCoefficient(*MQm, mesh, ir, indices, use_bdr, coeff.back());

      info.build_qf = f_build_curlcurl_mass_quad_vector_matrix;
      info.build_qf_path =
          PalaceQFunctionRelativePath(f_build_curlcurl_mass_quad_vector_matrix_loc);
    }
  }
  else if (MQc)
  {
    MFEM_VERIFY(MQc->GetVDim() == info.ctx.curl_dim,
                "Invalid matrix coefficient dimension for CurlCurlMassIntegrator!");
    coeff.emplace_back();
    ceed::InitCoefficient(*MQc, mesh, ir, indices, use_bdr, coeff.back());

    if (Qm)
    {
      coeff.emplace_back();
      ceed::InitCoefficient(*Qm, mesh, ir, indices, use_bdr, coeff.back());

      info.build_qf = f_build_curlcurl_mass_quad_matrix_scalar;
      info.build_qf_path =
          PalaceQFunctionRelativePath(f_build_curlcurl_mass_quad_matrix_scalar_loc);
    }
    else if (VQm)
    {
      MFEM_VERIFY(VQm->GetVDim() == info.ctx.space_dim,
                  "Invalid vector coefficient dimension for CurlCurlMassIntegrator!");
      coeff.emplace_back();
      ceed::InitCoefficient(*VQm, mesh, ir, indices, use_bdr, coeff.back());

      info.build_qf = f_build_curlcurl_mass_quad_matrix_vector;
      info.build_qf_path =
          PalaceQFunctionRelativePath(f_build_curlcurl_mass_quad_matrix_vector_loc);
    }
    else if (MQm)
    {
      MFEM_VERIFY(MQm->GetVDim() == info.ctx.space_dim,
                  "Invalid matrix coefficient dimension for CurlCurlMassIntegrator!");
      coeff.emplace_back();
      ceed::InitCoefficient(*MQm, mesh, ir, indices, use_bdr, coeff.back());

      info.build_qf = f_build_curlcurl_mass_quad_matrix_matrix;
      info.build_qf_path =
          PalaceQFunctionRelativePath(f_build_curlcurl_mass_quad_matrix_matrix_loc);
    }
  }

  info.apply_qf = f_apply_curlcurl_mass;
  info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_curlcurl_mass_loc);

  return info;
}

}  // namespace

void CurlCurlMassIntegrator::Assemble(const mfem::FiniteElementSpace &trial_fespace,
                                      const mfem::FiniteElementSpace &test_fespace,
                                      const mfem::IntegrationRule &ir,
                                      const std::vector<int> &indices, Ceed ceed,
                                      CeedOperator *op, CeedOperator *op_t)
{
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "CurlCurlMassIntegrator requires the same test and trial spaces!");
  constexpr bool use_bdr = false;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, ir, indices, use_bdr, Qc, VQc,
                                             MQc, Qm, VQm, MQm, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

void CurlCurlMassIntegrator::AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                                              const mfem::FiniteElementSpace &test_fespace,
                                              const mfem::IntegrationRule &ir,
                                              const std::vector<int> &indices, Ceed ceed,
                                              CeedOperator *op, CeedOperator *op_t)
{
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "CurlCurlMassIntegrator requires the same test and trial spaces!");
  constexpr bool use_bdr = true;
  std::vector<ceed::QuadratureCoefficient> coeff;
  const auto info = InitializeIntegratorInfo(trial_fespace, ir, indices, use_bdr, Qc, VQc,
                                             MQc, Qm, VQm, MQm, coeff);
  ceed::AssembleCeedOperator(info, trial_fespace, test_fespace, ir, indices, use_bdr, coeff,
                             ceed, op, op_t);
}

}  // namespace palace
