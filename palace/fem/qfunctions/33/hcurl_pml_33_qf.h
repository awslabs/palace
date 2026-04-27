// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_PML_33_QF_H
#define PALACE_LIBCEED_HCURL_PML_33_QF_H

#include "../coeff/coeff_qf.h"
#include "../coeff/pml_qf.h"
#include "utils_33_qf.h"

// Apply a(u, v) = (μ̃⁻¹ curl u, curl v) for a Cartesian PML region in 3D. The diagonal
// PML tensor μ̃⁻¹ is evaluated per quadrature point from the cached physical coordinate
// x (in the geom_data tail) and the PML profile context. `SELECT_IM` at compile time
// chooses between the real and imaginary parts of μ̃⁻¹ for the two branches of the
// complex stiffness operator.
//
// geom_data layout: {attr, w|J|, adj(J)^T/|J|, x}, size (2 + 3*3 + 3) = 14 components.
//   in[0] + 0       = attr
//   in[0] + Q       = w * |J|
//   in[0] + 2*Q     = adj(J)^T / |J|  (9 components, column-major)
//   in[0] + 11*Q    = x  (3 components: x, y, z)

template <bool SELECT_IM>
CEED_QFUNCTION_HELPER int f_apply_hcurl_pml_33_impl(void *__restrict__ ctx, CeedInt Q,
                                                    const CeedScalar *const *in,
                                                    CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q,
                   *xqp = in[0] + 11 * Q, *u = in[1];
  CeedScalar *v = out[0];
  const CeedIntScalar *pml_ctx = (const CeedIntScalar *)ctx;
  const CeedInt num_attr = PMLNumProfiles(pml_ctx);

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar adjJt_loc[9];
    MatUnpack33(adjJt + i, Q, adjJt_loc);

    CeedScalar coeff[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    const CeedInt a = (CeedInt)attr[i];
    const CeedInt pidx = (a >= 1 && a <= num_attr) ? PMLAttrToProfile(pml_ctx, a) : -1;
    if (pidx >= 0)
    {
      const CeedIntScalar *region = PMLRegion(pml_ctx, num_attr, pidx);
      const CeedScalar x_loc[3] = {xqp[i + Q * 0], xqp[i + Q * 1], xqp[i + Q * 2]};
      CeedScalar mi_re[3], mi_im[3], e_re[3], e_im[3];
      PMLEvalStretchTensors(region, x_loc, mi_re, mi_im, e_re, e_im);
      const CeedScalar *diag = SELECT_IM ? mi_im : mi_re;
      // Column-major diagonal: coeff[0]=xx, coeff[4]=yy, coeff[8]=zz.
      coeff[0] = diag[0];
      coeff[4] = diag[1];
      coeff[8] = diag[2];
    }
    // Non-PML attributes leave coeff at zero ⇒ no contribution (as intended: this
    // integrator is added alongside the bulk curl-curl and only "fills in" PML cells).

    CeedScalar v_loc[3];
    MultAtBCx33(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);
    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurl_pml_re_33)(void *__restrict__ ctx, CeedInt Q,
                                        const CeedScalar *const *in, CeedScalar *const *out)
{
  return f_apply_hcurl_pml_33_impl<false>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_im_33)(void *__restrict__ ctx, CeedInt Q,
                                        const CeedScalar *const *in, CeedScalar *const *out)
{
  return f_apply_hcurl_pml_33_impl<true>(ctx, Q, in, out);
}

#endif  // PALACE_LIBCEED_HCURL_PML_33_QF_H
