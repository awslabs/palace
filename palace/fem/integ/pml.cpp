// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/integrator.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/33/hcurl_pml_33_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

using namespace ceed;

namespace
{

// All PML integrators share everything except the QFunctions and the evaluation modes.
// When quadrature data assembly is requested and a build QFunction is available, the PML
// tensors are evaluated once at assembly and the operator is applied using the generic
// quadrature data QFunctions.
void AssemblePMLOperator(CeedQFunctionUser apply_qf, const char *apply_qf_loc,
                         CeedQFunctionUser build_qf, const char *build_qf_loc,
                         unsigned int ops, bool assemble_q_data, const void *ctx,
                         std::size_t ctx_size, Ceed ceed, CeedElemRestriction trial_restr,
                         CeedElemRestriction test_restr, CeedBasis trial_basis,
                         CeedBasis test_basis, CeedVector geom_data,
                         CeedElemRestriction geom_data_restr, CeedOperator *op)
{
  CeedInt dim, space_dim;
  PalaceCeedCall(ceed, CeedBasisGetDimension(trial_basis, &dim));
  PalaceCeedCall(ceed, CeedGeometryDataGetSpaceDimension(geom_data_restr, dim, &space_dim));
  MFEM_VERIFY(space_dim == 3 && dim == 3,
              "PML integrators are only available for 3D elements (got (space_dim, dim) = ("
                  << space_dim << ", " << dim << "))!");
  MFEM_VERIFY(CeedGeometryDataHasCoordinates(geom_data_restr, dim),
              "PML integrators require the quadrature point coordinates in the geometry "
              "factor data (see BilinearForm::AddDomainIntegratorOnAttributes)!");
  CeedQFunctionInfo info;
  info.assemble_q_data = assemble_q_data && build_qf;
  info.apply_qf = info.assemble_q_data ? build_qf : apply_qf;
  info.apply_qf_path =
      PalaceQFunctionRelativePath(info.assemble_q_data ? build_qf_loc : apply_qf_loc);
  info.trial_ops = ops;
  info.test_ops = ops;
  AssembleCeedOperator(info, const_cast<void *>(ctx), ctx_size, ceed, trial_restr,
                       test_restr, trial_basis, test_basis, geom_data, geom_data_restr, op);
}

}  // namespace

void CurlCurlPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                     CeedElemRestriction test_restr, CeedBasis trial_basis,
                                     CeedBasis test_basis, CeedVector geom_data,
                                     CeedElemRestriction geom_data_restr,
                                     CeedOperator *op) const
{
  AssemblePMLOperator(f_apply_hcurl_pml_curl_33, f_apply_hcurl_pml_curl_33_loc,
                      f_build_hcurl_pml_curl_33, f_build_hcurl_pml_curl_33_loc,
                      EvalMode::Curl, assemble_q_data, ctx, ctx_size, ceed, trial_restr,
                      test_restr, trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void VectorFEMassPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                         CeedElemRestriction test_restr,
                                         CeedBasis trial_basis, CeedBasis test_basis,
                                         CeedVector geom_data,
                                         CeedElemRestriction geom_data_restr,
                                         CeedOperator *op) const
{
  AssemblePMLOperator(f_apply_hcurl_pml_mass_33, f_apply_hcurl_pml_mass_33_loc,
                      f_build_hcurl_pml_mass_33, f_build_hcurl_pml_mass_33_loc,
                      EvalMode::Interp, assemble_q_data, ctx, ctx_size, ceed, trial_restr,
                      test_restr, trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void CurlCurlMassPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                         CeedElemRestriction test_restr,
                                         CeedBasis trial_basis, CeedBasis test_basis,
                                         CeedVector geom_data,
                                         CeedElemRestriction geom_data_restr,
                                         CeedOperator *op) const
{
  AssemblePMLOperator(f_apply_hcurl_pml_curlmass_33, f_apply_hcurl_pml_curlmass_33_loc,
                      f_build_hcurl_pml_curlmass_33, f_build_hcurl_pml_curlmass_33_loc,
                      EvalMode::Interp | EvalMode::Curl, assemble_q_data, ctx, ctx_size,
                      ceed, trial_restr, test_restr, trial_basis, test_basis, geom_data,
                      geom_data_restr, op);
}

void DiffusionPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                      CeedElemRestriction test_restr, CeedBasis trial_basis,
                                      CeedBasis test_basis, CeedVector geom_data,
                                      CeedElemRestriction geom_data_restr,
                                      CeedOperator *op) const
{
  // H1 gradients pull back from the reference element with adj(J)ᵀ / |J| like H(curl)
  // fields, so the mass QFunction computes the diffusion form for the gradient evaluation
  // mode.
  AssemblePMLOperator(f_apply_hcurl_pml_mass_33, f_apply_hcurl_pml_mass_33_loc,
                      f_build_hcurl_pml_mass_33, f_build_hcurl_pml_mass_33_loc,
                      EvalMode::Grad, assemble_q_data, ctx, ctx_size, ceed, trial_restr,
                      test_restr, trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void FloquetMassPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                        CeedElemRestriction test_restr,
                                        CeedBasis trial_basis, CeedBasis test_basis,
                                        CeedVector geom_data,
                                        CeedElemRestriction geom_data_restr,
                                        CeedOperator *op) const
{
  AssemblePMLOperator(
      f_apply_hcurl_pml_floquet_mass_33, f_apply_hcurl_pml_floquet_mass_33_loc,
      f_build_hcurl_pml_floquet_mass_33, f_build_hcurl_pml_floquet_mass_33_loc,
      EvalMode::Interp, assemble_q_data, ctx, ctx_size, ceed, trial_restr, test_restr,
      trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void FloquetDiffusionPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                             CeedElemRestriction test_restr,
                                             CeedBasis trial_basis, CeedBasis test_basis,
                                             CeedVector geom_data,
                                             CeedElemRestriction geom_data_restr,
                                             CeedOperator *op) const
{
  AssemblePMLOperator(
      f_apply_hcurl_pml_floquet_mass_33, f_apply_hcurl_pml_floquet_mass_33_loc,
      f_build_hcurl_pml_floquet_mass_33, f_build_hcurl_pml_floquet_mass_33_loc,
      EvalMode::Grad, assemble_q_data, ctx, ctx_size, ceed, trial_restr, test_restr,
      trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void FloquetCrossPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                         CeedElemRestriction test_restr,
                                         CeedBasis trial_basis, CeedBasis test_basis,
                                         CeedVector geom_data,
                                         CeedElemRestriction geom_data_restr,
                                         CeedOperator *op) const
{
  // The cross term couples the two active fields, so it has no quadrature data form for the
  // generic block-diagonal f_apply_33 QFunction.
  AssemblePMLOperator(
      f_apply_hcurl_pml_floquet_cross_33, f_apply_hcurl_pml_floquet_cross_33_loc, nullptr,
      nullptr, EvalMode::Interp | EvalMode::Curl, false, ctx, ctx_size, ceed, trial_restr,
      test_restr, trial_basis, test_basis, geom_data, geom_data_restr, op);
}

}  // namespace palace
