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

void CurlCurlPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                     CeedElemRestriction test_restr, CeedBasis trial_basis,
                                     CeedBasis test_basis, CeedVector geom_data,
                                     CeedElemRestriction geom_data_restr,
                                     CeedOperator *op) const
{
  CeedQFunctionInfo info;
  info.assemble_q_data = false;  // PML path doesn't support pre-assembled quadrature
                                 // data since σ(x) varies per QP.

  CeedInt dim, space_dim;
  PalaceCeedCall(ceed, CeedBasisGetDimension(trial_basis, &dim));
  PalaceCeedCall(ceed, CeedGeometryDataGetSpaceDimension(geom_data_restr, dim, &space_dim));
  MFEM_VERIFY(space_dim == 3 && dim == 3,
              "CurlCurlPMLIntegrator currently only supports 3D elements "
              "(got (space_dim, dim) = ("
                  << space_dim << ", " << dim << ")). Cartesian PML in 2D is a v2 item.");
  info.apply_qf = imag_part ? f_apply_hcurl_pml_im_33 : f_apply_hcurl_pml_re_33;
  info.apply_qf_path = PalaceQFunctionRelativePath(imag_part ? f_apply_hcurl_pml_im_33_loc
                                                             : f_apply_hcurl_pml_re_33_loc);

  info.trial_ops = EvalMode::Curl;
  info.test_ops = EvalMode::Curl;

  AssembleCeedOperator(info, const_cast<void *>(ctx), ctx_size, ceed, trial_restr,
                       test_restr, trial_basis, test_basis, geom_data, geom_data_restr, op);
}

void VectorFEMassPMLIntegrator::Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                                         CeedElemRestriction test_restr,
                                         CeedBasis trial_basis, CeedBasis test_basis,
                                         CeedVector geom_data,
                                         CeedElemRestriction geom_data_restr,
                                         CeedOperator *op) const
{
  CeedQFunctionInfo info;
  info.assemble_q_data = false;

  CeedInt dim, space_dim;
  PalaceCeedCall(ceed, CeedBasisGetDimension(trial_basis, &dim));
  PalaceCeedCall(ceed, CeedGeometryDataGetSpaceDimension(geom_data_restr, dim, &space_dim));
  MFEM_VERIFY(space_dim == 3 && dim == 3,
              "VectorFEMassPMLIntegrator currently only supports 3D elements.");
  info.apply_qf = imag_part ? f_apply_hcurl_pml_im_33 : f_apply_hcurl_pml_re_33;
  info.apply_qf_path = PalaceQFunctionRelativePath(imag_part ? f_apply_hcurl_pml_im_33_loc
                                                             : f_apply_hcurl_pml_re_33_loc);

  info.trial_ops = EvalMode::Interp;
  info.test_ops = EvalMode::Interp;

  AssembleCeedOperator(info, const_cast<void *>(ctx), ctx_size, ceed, trial_restr,
                       test_restr, trial_basis, test_basis, geom_data, geom_data_restr, op);
}

}  // namespace palace
