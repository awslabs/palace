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

// Build the libCEED operator for a 3D PML integrator. Both CurlCurl (μ̃⁻¹) and
// VectorFEMass (ε̃) paths share everything except the QFunction pointer and the eval
// mode, so they factor into a common helper.
void AssemblePML33(CeedQFunctionUser apply_qf, const char *apply_qf_path,
                   EvalMode eval_mode, const void *ctx, std::size_t ctx_size, Ceed ceed,
                   CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                   CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
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

}  // namespace

void CurlCurlPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                     CeedElemRestriction test_restr, CeedBasis trial_basis,
                                     CeedBasis test_basis, CeedVector geom_data,
                                     CeedElemRestriction geom_data_restr,
                                     CeedOperator *op) const
{
  AssemblePML33(imag_part ? f_apply_hcurl_pml_muinv_im_33 : f_apply_hcurl_pml_muinv_re_33,
                imag_part ? f_apply_hcurl_pml_muinv_im_33_loc
                          : f_apply_hcurl_pml_muinv_re_33_loc,
                EvalMode::Curl, ctx, ctx_size, ceed, trial_restr, test_restr, trial_basis,
                test_basis, geom_data, geom_data_restr, op);
}

void VectorFEMassPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                         CeedElemRestriction test_restr,
                                         CeedBasis trial_basis, CeedBasis test_basis,
                                         CeedVector geom_data,
                                         CeedElemRestriction geom_data_restr,
                                         CeedOperator *op) const
{
  AssemblePML33(imag_part ? f_apply_hcurl_pml_eps_im_33 : f_apply_hcurl_pml_eps_re_33,
                imag_part ? f_apply_hcurl_pml_eps_im_33_loc
                          : f_apply_hcurl_pml_eps_re_33_loc,
                EvalMode::Interp, ctx, ctx_size, ceed, trial_restr, test_restr, trial_basis,
                test_basis, geom_data, geom_data_restr, op);
}

}  // namespace palace
