// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_PML_33_QF_H
#define PALACE_LIBCEED_HCURL_PML_33_QF_H

#include "../coeff/coeff_qf.h"
#include "../coeff/pml_qf.h"
#include "utils_33_qf.h"

// PML QFunctions for 3D Cartesian UPML. There are two families:
//
//   - MU path (curl-curl): active input is the curl of u, which lives in H(div). The
//     transform from reference to physical for H(div) values with weight |J| is J / |J|
//     (contravariant). We compute v_phys = J · (coeff · (J · u)) · w|J| / |J| = ...
//     effectively `v = J^T · coeff · J · u` multiplied by w, up to the |J| weights.
//     Mirrors f_apply_hdiv_33 from hdiv_33_qf.h.
//
//   - EPS path (mass on H(curl)): active input is u in H(curl), pulled back by
//     adj(J)^T / |J|. Mirrors f_apply_hcurl_33 from hcurl_33_qf.h.
//
// In both cases, `coeff` is a diagonal per-QP matrix built from the PML stretch tensor
// μ̃⁻¹ (MU path) or ε̃ (EPS path), selected real/imag by template parameter SELECT_IM.
//
// geom_data layout: {attr, w|J|, adj(J)^T/|J|, x}, size (2 + 3*3 + 3) = 14 components.
//   in[0] + 0       = attr
//   in[0] + Q       = w * |J|
//   in[0] + 2*Q     = adj(J)^T / |J|  (9 components, column-major)
//   in[0] + 11*Q    = x  (3 components: x, y, z)

// Read μ̃⁻¹ or ε̃ diagonal for attr a at physical x, apply PML scale, return in coeff[9]
// column-major. Non-PML attrs (or scale 0) leave coeff at zero.
CEED_QFUNCTION_HELPER void BuildPMLCoeff33(const CeedIntScalar *pml_ctx, CeedInt num_attr,
                                           CeedInt a, const CeedScalar x[3],
                                           bool select_eps, bool select_im,
                                           CeedScalar scale, CeedScalar coeff[9])
{
  coeff[0] = coeff[1] = coeff[2] = coeff[3] = coeff[4] = coeff[5] = coeff[6] = coeff[7] =
      coeff[8] = 0.0;
  const CeedInt pidx = (a >= 1 && a <= num_attr) ? PMLAttrToProfile(pml_ctx, a) : -1;
  if (pidx < 0)
  {
    return;
  }
  const CeedIntScalar *region = PMLRegion(pml_ctx, num_attr, pidx);
  CeedScalar mi_re[3], mi_im[3], e_re[3], e_im[3];
  PMLEvalStretchTensors(region, x, mi_re, mi_im, e_re, e_im);
  const CeedScalar *diag;
  if (select_eps)
  {
    diag = select_im ? e_im : e_re;
  }
  else
  {
    diag = select_im ? mi_im : mi_re;
  }
  coeff[0] = scale * diag[0];
  coeff[4] = scale * diag[1];
  coeff[8] = scale * diag[2];
}

// Curl-curl PML (uses J for the transform, like f_apply_hdiv_33).
template <bool SELECT_IM>
CEED_QFUNCTION_HELPER int
f_apply_hcurl_pml_curlcurl_33_impl(void *__restrict__ ctx, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q,
                   *xqp = in[0] + 11 * Q, *u = in[1];
  CeedScalar *v = out[0];
  const CeedIntScalar *pml_ctx = (const CeedIntScalar *)ctx;
  const CeedScalar scale = PMLScale(pml_ctx);
  const CeedInt num_attr = PMLNumProfiles(pml_ctx);

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar adjJt_loc[9], J_loc[9], coeff[9];
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    AdjJt33(adjJt_loc, J_loc);

    const CeedScalar x_loc[3] = {xqp[i + Q * 0], xqp[i + Q * 1], xqp[i + Q * 2]};
    BuildPMLCoeff33(pml_ctx, num_attr, (CeedInt)attr[i], x_loc, /*select_eps=*/false,
                    SELECT_IM, scale, coeff);

    CeedScalar v_loc[3];
    MultAtBCx33(J_loc, coeff, J_loc, u_loc, v_loc);
    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
  }
  return 0;
}

// H(curl)-mass PML (uses adj(J)^T for the transform, like f_apply_hcurl_33).
template <bool SELECT_IM>
CEED_QFUNCTION_HELPER int f_apply_hcurl_pml_mass_33_impl(void *__restrict__ ctx, CeedInt Q,
                                                         const CeedScalar *const *in,
                                                         CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q,
                   *xqp = in[0] + 11 * Q, *u = in[1];
  CeedScalar *v = out[0];
  const CeedIntScalar *pml_ctx = (const CeedIntScalar *)ctx;
  const CeedScalar scale = PMLScale(pml_ctx);
  const CeedInt num_attr = PMLNumProfiles(pml_ctx);

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar adjJt_loc[9], coeff[9];
    MatUnpack33(adjJt + i, Q, adjJt_loc);

    const CeedScalar x_loc[3] = {xqp[i + Q * 0], xqp[i + Q * 1], xqp[i + Q * 2]};
    BuildPMLCoeff33(pml_ctx, num_attr, (CeedInt)attr[i], x_loc, /*select_eps=*/true,
                    SELECT_IM, scale, coeff);

    CeedScalar v_loc[3];
    MultAtBCx33(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);
    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurl_pml_muinv_re_33)(void *__restrict__ ctx, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  return f_apply_hcurl_pml_curlcurl_33_impl<false>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_muinv_im_33)(void *__restrict__ ctx, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  return f_apply_hcurl_pml_curlcurl_33_impl<true>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_eps_re_33)(void *__restrict__ ctx, CeedInt Q,
                                            const CeedScalar *const *in,
                                            CeedScalar *const *out)
{
  return f_apply_hcurl_pml_mass_33_impl<false>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_eps_im_33)(void *__restrict__ ctx, CeedInt Q,
                                            const CeedScalar *const *in,
                                            CeedScalar *const *out)
{
  return f_apply_hcurl_pml_mass_33_impl<true>(ctx, Q, in, out);
}

#endif  // PALACE_LIBCEED_HCURL_PML_33_QF_H
