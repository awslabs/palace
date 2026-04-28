// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_PML_33_QF_H
#define PALACE_LIBCEED_HCURL_PML_33_QF_H

#include "../coeff/coeff_qf.h"
#include "../coeff/pml_qf.h"
#include "utils_33_qf.h"

// PML QFunctions for 3D Cartesian UPML. Four entry points: two tensors × two parts.
//
//   - muinv_{re,im}: active input is the curl of u, which lives in H(div). The transform
//     from reference to physical for H(div) is J / |J| (contravariant), so the QFunction
//     computes v = J^T · coeff · J · u (up to |J| weights). Mirrors f_apply_hdiv_33.
//
//   - eps_{re,im}: active input is u in H(curl), pulled back by adj(J)^T / |J| (covariant),
//     so v = adj(J)^T · coeff · adj(J) · u. Mirrors f_apply_hcurl_33.
//
// In all four, `coeff` is a diagonal per-QP 3×3 matrix built from the PML stretch tensor
// μ̃⁻¹ (curl-curl) or ε̃ (mass), real or imaginary part, scaled by the context's scalar.
//
// geom_data layout: {attr, w|J|, adj(J)^T/|J|, x}, size (2 + 3*3 + 3) = 14 components.
//   in[0] + 0       = attr
//   in[0] + Q       = w * |J|
//   in[0] + 2*Q     = adj(J)^T / |J|  (9 components, column-major)
//   in[0] + 11*Q    = x               (3 components)

enum
{
  PML_TENSOR_MUINV = 0,
  PML_TENSOR_EPS = 1
};

// Build the diagonal PML coefficient matrix for attr `a` at physical point `x`. Writes
// 9 entries in column-major order; non-PML attrs or an out-of-range pidx leave coeff at
// zero. TENSOR selects μ̃⁻¹ vs ε̃; IMAG selects the real vs imaginary part.
template <int TENSOR, bool IMAG>
CEED_QFUNCTION_HELPER void BuildPMLDiag33(const CeedIntScalar *pml_ctx, CeedInt num_attr,
                                          CeedInt a, const CeedScalar x[3],
                                          CeedScalar scale, CeedScalar coeff[9])
{
  coeff[0] = coeff[1] = coeff[2] = coeff[3] = coeff[4] = coeff[5] = coeff[6] = coeff[7] =
      coeff[8] = 0.0;
  const CeedInt pidx = (a >= 1 && a <= num_attr) ? PMLAttrToProfile(pml_ctx, a) : -1;
  if (pidx < 0)
  {
    return;
  }
  CeedScalar mi_re[3], mi_im[3], e_re[3], e_im[3];
  PMLEvalStretchTensors(PMLRegion(pml_ctx, num_attr, pidx), x, mi_re, mi_im, e_re, e_im);
  const CeedScalar *diag;
  if (TENSOR == PML_TENSOR_EPS)
  {
    diag = IMAG ? e_im : e_re;
  }
  else
  {
    diag = IMAG ? mi_im : mi_re;
  }
  coeff[0] = scale * diag[0];
  coeff[4] = scale * diag[1];
  coeff[8] = scale * diag[2];
}

// Apply the PML operator v = T^T · coeff · T · u at each QP, with T = J (curl-curl) or
// T = adj(J)^T (mass). Unified across the four entry points below.
template <int TENSOR, bool IMAG, bool USE_J>
CEED_QFUNCTION_HELPER int f_apply_hcurl_pml_33_impl(void *__restrict__ ctx, CeedInt Q,
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
    CeedScalar adjJt_loc[9], J_loc[9], coeff[9];
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    if (USE_J)
    {
      AdjJt33(adjJt_loc, J_loc);
    }
    const CeedScalar *T_loc = USE_J ? J_loc : adjJt_loc;

    const CeedScalar x_loc[3] = {xqp[i + Q * 0], xqp[i + Q * 1], xqp[i + Q * 2]};
    BuildPMLDiag33<TENSOR, IMAG>(pml_ctx, num_attr, (CeedInt)attr[i], x_loc, scale, coeff);

    CeedScalar v_loc[3];
    MultAtBCx33(T_loc, coeff, T_loc, u_loc, v_loc);
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
  return f_apply_hcurl_pml_33_impl<PML_TENSOR_MUINV, false, true>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_muinv_im_33)(void *__restrict__ ctx, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  return f_apply_hcurl_pml_33_impl<PML_TENSOR_MUINV, true, true>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_eps_re_33)(void *__restrict__ ctx, CeedInt Q,
                                            const CeedScalar *const *in,
                                            CeedScalar *const *out)
{
  return f_apply_hcurl_pml_33_impl<PML_TENSOR_EPS, false, false>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_eps_im_33)(void *__restrict__ ctx, CeedInt Q,
                                            const CeedScalar *const *in,
                                            CeedScalar *const *out)
{
  return f_apply_hcurl_pml_33_impl<PML_TENSOR_EPS, true, false>(ctx, Q, in, out);
}

#endif  // PALACE_LIBCEED_HCURL_PML_33_QF_H
