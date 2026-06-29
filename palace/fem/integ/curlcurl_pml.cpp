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

// Build the libCEED operator for a 3D PML integrator. All PML paths share everything
// except the QFunction pointer and the eval mode, so they factor into a common helper.
void AssemblePML33(CeedQFunctionUser apply_qf, const char *apply_qf_path,
                   unsigned int eval_mode, const void *ctx, std::size_t ctx_size,
                   Ceed ceed, CeedElemRestriction trial_restr,
                   CeedElemRestriction test_restr, CeedBasis trial_basis,
                   CeedBasis test_basis, CeedVector geom_data,
                   CeedElemRestriction geom_data_restr, CeedOperator *op)
{
  CeedInt dim, space_dim;
  PalaceCeedCall(ceed, CeedBasisGetDimension(trial_basis, &dim));
  PalaceCeedCall(ceed, CeedGeometryDataGetSpaceDimension(geom_data_restr, dim, &space_dim));
  MFEM_VERIFY(space_dim == 3 && dim == 3,
              "Cartesian PML integrators currently only support 3D elements "
              "(got (space_dim, dim) = ("
                  << space_dim << ", " << dim << ")). 2D is a v2 item.");
  CeedQFunctionInfo info;
  info.assemble_q_data = false;  // σ(x) varies per QP — pre-assembly is not supported.
  info.apply_qf = apply_qf;
  info.apply_qf_path = PalaceQFunctionRelativePath(apply_qf_path);
  info.trial_ops = eval_mode;
  info.test_ops = eval_mode;
  AssembleCeedOperator(info, const_cast<void *>(ctx), ctx_size, ceed, trial_restr,
                       test_restr, trial_basis, test_basis, geom_data, geom_data_restr, op);
}

// Select the μ̃⁻¹ or ε̃ QFunction (real or imaginary part) and its loc string.
struct QFPair
{
  CeedQFunctionUser qf;
  const char *loc;
};
QFPair MuInvQF(PMLTensorPart part)
{
  switch (part)
  {
    case PMLTensorPart::Re:
      return {f_apply_hcurl_pml_muinv_re_33, f_apply_hcurl_pml_muinv_re_33_loc};
    case PMLTensorPart::Im:
      return {f_apply_hcurl_pml_muinv_im_33, f_apply_hcurl_pml_muinv_im_33_loc};
  }
  MFEM_ABORT("Unreachable PMLTensorPart");
}
QFPair EpsQF(PMLTensorPart part)
{
  switch (part)
  {
    case PMLTensorPart::Re:
      return {f_apply_hcurl_pml_eps_re_33, f_apply_hcurl_pml_eps_re_33_loc};
    case PMLTensorPart::Im:
      return {f_apply_hcurl_pml_eps_im_33, f_apply_hcurl_pml_eps_im_33_loc};
  }
  MFEM_ABORT("Unreachable PMLTensorPart");
}
QFPair FloquetMassQF(PMLTensorPart part)
{
  switch (part)
  {
    case PMLTensorPart::Re:
      return {f_apply_hcurl_pml_floquet_mass_re_33,
              f_apply_hcurl_pml_floquet_mass_re_33_loc};
    case PMLTensorPart::Im:
      return {f_apply_hcurl_pml_floquet_mass_im_33,
              f_apply_hcurl_pml_floquet_mass_im_33_loc};
  }
  MFEM_ABORT("Unreachable PMLTensorPart");
}
QFPair FloquetCrossQF(PMLTensorPart part)
{
  switch (part)
  {
    case PMLTensorPart::Re:
      return {f_apply_hcurl_pml_floquet_cross_re_33,
              f_apply_hcurl_pml_floquet_cross_re_33_loc};
    case PMLTensorPart::Im:
      return {f_apply_hcurl_pml_floquet_cross_im_33,
              f_apply_hcurl_pml_floquet_cross_im_33_loc};
  }
  MFEM_ABORT("Unreachable PMLTensorPart");
}

}  // namespace

void CurlCurlPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                     CeedElemRestriction test_restr, CeedBasis trial_basis,
                                     CeedBasis test_basis, CeedVector geom_data,
                                     CeedElemRestriction geom_data_restr,
                                     CeedOperator *op) const
{
  const auto qf = MuInvQF(part);
  AssemblePML33(qf.qf, qf.loc, EvalMode::Curl, ctx, ctx_size, ceed, trial_restr, test_restr,
                trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void VectorFEMassPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                         CeedElemRestriction test_restr,
                                         CeedBasis trial_basis, CeedBasis test_basis,
                                         CeedVector geom_data,
                                         CeedElemRestriction geom_data_restr,
                                         CeedOperator *op) const
{
  const auto qf = EpsQF(part);
  AssemblePML33(qf.qf, qf.loc, EvalMode::Interp, ctx, ctx_size, ceed, trial_restr,
                test_restr, trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void FloquetMassPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                        CeedElemRestriction test_restr,
                                        CeedBasis trial_basis, CeedBasis test_basis,
                                        CeedVector geom_data,
                                        CeedElemRestriction geom_data_restr,
                                        CeedOperator *op) const
{
  const auto qf = FloquetMassQF(part);
  AssemblePML33(qf.qf, qf.loc, EvalMode::Interp, ctx, ctx_size, ceed, trial_restr,
                test_restr, trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void FloquetCrossPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                         CeedElemRestriction test_restr,
                                         CeedBasis trial_basis, CeedBasis test_basis,
                                         CeedVector geom_data,
                                         CeedElemRestriction geom_data_restr,
                                         CeedOperator *op) const
{
  const auto qf = FloquetCrossQF(part);
  constexpr unsigned int ops = EvalMode::Interp | EvalMode::Curl;
  AssemblePML33(qf.qf, qf.loc, ops, ctx, ctx_size, ceed, trial_restr, test_restr,
                trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void DiffusionPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                      CeedElemRestriction test_restr, CeedBasis trial_basis,
                                      CeedBasis test_basis, CeedVector geom_data,
                                      CeedElemRestriction geom_data_restr,
                                      CeedOperator *op) const
{
  // H1 gradient transforms via adj(J)^T just like H(curl) fields, so the eps QFunction
  // (which applies adj(J)^T·ε̃·adj(J)^T·u with weight w|J|) computes the diffusion form
  // (ε̃ ∇φ, ∇ψ) when the basis eval mode is Grad instead of Interp.
  const auto qf = EpsQF(part);
  AssemblePML33(qf.qf, qf.loc, EvalMode::Grad, ctx, ctx_size, ceed, trial_restr, test_restr,
                trial_basis, test_basis, geom_data, geom_data_restr, op);
}

}  // namespace palace
