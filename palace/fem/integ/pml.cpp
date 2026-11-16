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

// All PML integrators share everything except the QFunctions and the evaluation modes. The
// PML tensors are evaluated once at assembly by the build QFunction and stored, and the
// operator is applied using the generic quadrature data QFunctions, f_apply_33_cross for
// the Floquet cross term (cross).
void AssemblePMLOperator(CeedQFunctionUser build_qf, const char *build_qf_loc,
                         unsigned int ops, const void *ctx, std::size_t ctx_size, Ceed ceed,
                         CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                         CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                         CeedElemRestriction geom_data_restr, CeedOperator *op,
                         bool cross = false)
{
  // The geometry factor data of the PML attributes has the coordinates of the quadrature
  // points (see BilinearForm::AddDomainIntegratorOnAttributes).
  CeedInt dim, geom_data_size;
  PalaceCeedCall(ceed, CeedBasisGetDimension(trial_basis, &dim));
  PalaceCeedCall(ceed,
                 CeedElemRestrictionGetNumComponents(geom_data_restr, &geom_data_size));
  MFEM_VERIFY(dim == 3, "PML integrators are only available for 3D elements!");
  MFEM_VERIFY(
      geom_data_size == CeedGeometryDataSize(3, 3, true),
      "PML integrators require the geometry factor data of the PML attributes, with "
      "the coordinates of the quadrature points!");
  CeedQFunctionInfo info;
  info.assemble_q_data = true;
  info.q_data_cross = cross;
  info.apply_qf = build_qf;
  info.apply_qf_path = PalaceQFunctionRelativePath(build_qf_loc);
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
  AssemblePMLOperator(f_build_hcurl_pml_curl_33, f_build_hcurl_pml_curl_33_loc,
                      EvalMode::Curl, ctx, ctx_size, ceed, trial_restr, test_restr,
                      trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void VectorFEMassPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                         CeedElemRestriction test_restr,
                                         CeedBasis trial_basis, CeedBasis test_basis,
                                         CeedVector geom_data,
                                         CeedElemRestriction geom_data_restr,
                                         CeedOperator *op) const
{
  AssemblePMLOperator(f_build_hcurl_pml_mass_33, f_build_hcurl_pml_mass_33_loc,
                      EvalMode::Interp, ctx, ctx_size, ceed, trial_restr, test_restr,
                      trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void CurlCurlMassPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                         CeedElemRestriction test_restr,
                                         CeedBasis trial_basis, CeedBasis test_basis,
                                         CeedVector geom_data,
                                         CeedElemRestriction geom_data_restr,
                                         CeedOperator *op) const
{
  AssemblePMLOperator(f_build_hcurl_pml_curlmass_33, f_build_hcurl_pml_curlmass_33_loc,
                      EvalMode::Interp | EvalMode::Curl, ctx, ctx_size, ceed, trial_restr,
                      test_restr, trial_basis, test_basis, geom_data, geom_data_restr, op);
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
  AssemblePMLOperator(f_build_hcurl_pml_mass_33, f_build_hcurl_pml_mass_33_loc,
                      EvalMode::Grad, ctx, ctx_size, ceed, trial_restr, test_restr,
                      trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void FloquetMassPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                        CeedElemRestriction test_restr,
                                        CeedBasis trial_basis, CeedBasis test_basis,
                                        CeedVector geom_data,
                                        CeedElemRestriction geom_data_restr,
                                        CeedOperator *op) const
{
  AssemblePMLOperator(f_build_hcurl_pml_floquet_mass_33,
                      f_build_hcurl_pml_floquet_mass_33_loc, EvalMode::Interp, ctx,
                      ctx_size, ceed, trial_restr, test_restr, trial_basis, test_basis,
                      geom_data, geom_data_restr, op);
}

void FloquetDiffusionPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                             CeedElemRestriction test_restr,
                                             CeedBasis trial_basis, CeedBasis test_basis,
                                             CeedVector geom_data,
                                             CeedElemRestriction geom_data_restr,
                                             CeedOperator *op) const
{
  AssemblePMLOperator(f_build_hcurl_pml_floquet_mass_33,
                      f_build_hcurl_pml_floquet_mass_33_loc, EvalMode::Grad, ctx, ctx_size,
                      ceed, trial_restr, test_restr, trial_basis, test_basis, geom_data,
                      geom_data_restr, op);
}

void FloquetCrossPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                         CeedElemRestriction test_restr,
                                         CeedBasis trial_basis, CeedBasis test_basis,
                                         CeedVector geom_data,
                                         CeedElemRestriction geom_data_restr,
                                         CeedOperator *op) const
{
  // The cross term couples the field to the curl and the curl to the field, so its
  // quadrature data is applied by f_apply_33_cross rather than f_apply_33.
  AssemblePMLOperator(
      f_build_hcurl_pml_floquet_cross_33, f_build_hcurl_pml_floquet_cross_33_loc,
      EvalMode::Interp | EvalMode::Curl, ctx, ctx_size, ceed, trial_restr, test_restr,
      trial_basis, test_basis, geom_data, geom_data_restr, op, true);
}

}  // namespace palace
