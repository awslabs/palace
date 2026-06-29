// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_PML_33_QF_H
#define PALACE_LIBCEED_HCURL_PML_33_QF_H

#include "../coeff/coeff_qf.h"
#include "../coeff/pml_qf.h"
#include "utils_33_qf.h"

// PML QFunctions for 3D Cartesian UPML. The basic entry points are:
//
//   muinv_{re,im} : active input is the curl of u, which lives in H(div). The transform
//                   from reference to physical for H(div) is J / |J| (contravariant), so
//                   the QFunction computes v = J^T · coeff · J · u (up to |J| weights).
//                   Mirrors f_apply_hdiv_33.
//
//   eps_{re,im}   : active input is u in H(curl), pulled back by adj(J)^T / |J|
//                   (covariant). Mirrors f_apply_hcurl_33.
//
// Floquet variants add the missing quasi-periodic terms inside PML regions:
//
//   floquet_mass_{re,im}  : Kᵀ μ̃⁻¹ K, active input/output is u in H(curl).
//   floquet_cross_{re,im} : (Kᵀ μ̃⁻¹ curl u, v) − (μ̃⁻¹ K u, curl v),
//                           active input/output uses both interp and curl.
//
// In all variants, `coeff` is a diagonal per-QP 3×3 matrix; scale is applied per region.
// The Floquet variants prepend the 3×3 cross-product matrix K to the normal PML context.
//
// geom_data layout: {attr, w|J|, adj(J)^T/|J|, x}, size (2 + 3*3 + 3) = 14 components.

enum
{
  PML_TENSOR_MUINV = 0,
  PML_TENSOR_EPS = 1
};

enum
{
  PML_PART_RE = 0,
  PML_PART_IM = 1
};

// Build the diagonal PML coefficient matrix for attr `a` at physical point `x`. Writes
// 9 entries in column-major order; non-PML attrs or an out-of-range pidx leave coeff at
// zero. TENSOR selects μ̃⁻¹ vs ε̃; PART selects the real or imaginary part.
template <int TENSOR, int PART>
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
  const CeedScalar *diag = (TENSOR == PML_TENSOR_EPS)
                               ? (PART == PML_PART_RE ? e_re : e_im)
                               : (PART == PML_PART_RE ? mi_re : mi_im);
  coeff[0] = scale * diag[0];
  coeff[4] = scale * diag[1];
  coeff[8] = scale * diag[2];
}

CEED_QFUNCTION_HELPER void UnpackFloquetCross33(const CeedIntScalar *ctx,
                                                CeedScalar K[9])
{
  K[0] = ctx[0].second;
  K[1] = ctx[1].second;
  K[2] = ctx[2].second;
  K[3] = ctx[3].second;
  K[4] = ctx[4].second;
  K[5] = ctx[5].second;
  K[6] = ctx[6].second;
  K[7] = ctx[7].second;
  K[8] = ctx[8].second;
}

CEED_QFUNCTION_HELPER void MultAtx33(const CeedScalar A[9], const CeedScalar x[3],
                                     CeedScalar y[3])
{
  y[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
  y[1] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
  y[2] = A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
}

// Apply the PML operator v = T^T · coeff · T · u at each QP, with T = J (curl-curl) or
// T = adj(J)^T (mass/diffusion). Unified across the four entry points below.
template <int TENSOR, int PART, bool USE_J>
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
    BuildPMLDiag33<TENSOR, PART>(pml_ctx, num_attr, (CeedInt)attr[i], x_loc, scale, coeff);

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
  return f_apply_hcurl_pml_33_impl<PML_TENSOR_MUINV, PML_PART_RE, true>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_muinv_im_33)(void *__restrict__ ctx, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  return f_apply_hcurl_pml_33_impl<PML_TENSOR_MUINV, PML_PART_IM, true>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_eps_re_33)(void *__restrict__ ctx, CeedInt Q,
                                            const CeedScalar *const *in,
                                            CeedScalar *const *out)
{
  return f_apply_hcurl_pml_33_impl<PML_TENSOR_EPS, PML_PART_RE, false>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_eps_im_33)(void *__restrict__ ctx, CeedInt Q,
                                            const CeedScalar *const *in,
                                            CeedScalar *const *out)
{
  return f_apply_hcurl_pml_33_impl<PML_TENSOR_EPS, PML_PART_IM, false>(ctx, Q, in, out);
}

template <int PART>
CEED_QFUNCTION_HELPER int f_apply_hcurl_pml_floquet_mass_33_impl(
    void *__restrict__ ctx, CeedInt Q, const CeedScalar *const *in,
    CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q,
                   *xqp = in[0] + 11 * Q, *u = in[1];
  CeedScalar *v = out[0];
  const CeedIntScalar *floquet_ctx = (const CeedIntScalar *)ctx;
  const CeedIntScalar *pml_ctx = floquet_ctx + 9;
  const CeedScalar scale = PMLScale(pml_ctx);
  const CeedInt num_attr = PMLNumProfiles(pml_ctx);
  CeedScalar K[9];
  UnpackFloquetCross33(floquet_ctx, K);

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar adjJt_loc[9], coeff[9], kTcoeffk[9], v_loc[3];
    MatUnpack33(adjJt + i, Q, adjJt_loc);

    const CeedScalar x_loc[3] = {xqp[i + Q * 0], xqp[i + Q * 1], xqp[i + Q * 2]};
    BuildPMLDiag33<PML_TENSOR_MUINV, PART>(pml_ctx, num_attr, (CeedInt)attr[i], x_loc,
                                           scale, coeff);

    MultAtBC33(K, coeff, K, kTcoeffk);
    MultAtBCx33(adjJt_loc, kTcoeffk, adjJt_loc, u_loc, v_loc);
    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurl_pml_floquet_mass_re_33)(void *__restrict__ ctx, CeedInt Q,
                                                     const CeedScalar *const *in,
                                                     CeedScalar *const *out)
{
  return f_apply_hcurl_pml_floquet_mass_33_impl<PML_PART_RE>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_floquet_mass_im_33)(void *__restrict__ ctx, CeedInt Q,
                                                     const CeedScalar *const *in,
                                                     CeedScalar *const *out)
{
  return f_apply_hcurl_pml_floquet_mass_33_impl<PML_PART_IM>(ctx, Q, in, out);
}

template <int PART>
CEED_QFUNCTION_HELPER int f_apply_hcurl_pml_floquet_cross_33_impl(
    void *__restrict__ ctx, CeedInt Q, const CeedScalar *const *in,
    CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q,
                   *xqp = in[0] + 11 * Q, *u = in[1], *curl_u = in[2];
  CeedScalar *v = out[0], *curl_v = out[1];
  const CeedIntScalar *floquet_ctx = (const CeedIntScalar *)ctx;
  const CeedIntScalar *pml_ctx = floquet_ctx + 9;
  const CeedScalar scale = PMLScale(pml_ctx);
  const CeedInt num_attr = PMLNumProfiles(pml_ctx);
  CeedScalar K[9];
  UnpackFloquetCross33(floquet_ctx, K);

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    const CeedScalar curl_u_loc[3] = {curl_u[i + Q * 0], curl_u[i + Q * 1],
                                      curl_u[i + Q * 2]};
    CeedScalar adjJt_loc[9], J_loc[9], coeff[9];
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    AdjJt33(adjJt_loc, J_loc);

    const CeedScalar x_loc[3] = {xqp[i + Q * 0], xqp[i + Q * 1], xqp[i + Q * 2]};
    BuildPMLDiag33<PML_TENSOR_MUINV, PART>(pml_ctx, num_attr, (CeedInt)attr[i], x_loc,
                                           scale, coeff);

    // Curl-to-field half: (K^T μ̃⁻¹ curl u, v), matching the transposed
    // MixedVectorCurlIntegrator used by the non-PML Floquet path.
    CeedScalar curl_phys[3], mu_curl[3], kt_mu_curl[3], v_loc[3];
    MultBx33(J_loc, curl_u_loc, curl_phys);
    MultBx33(coeff, curl_phys, mu_curl);
    MultAtx33(K, mu_curl, kt_mu_curl);
    MultAtx33(adjJt_loc, kt_mu_curl, v_loc);

    // Field-to-curl half: −(μ̃⁻¹ K u, curl v), matching MixedVectorWeakCurlIntegrator.
    CeedScalar u_phys[3], ku_phys[3], mu_ku[3], curl_v_loc[3];
    MultBx33(adjJt_loc, u_loc, u_phys);
    MultBx33(K, u_phys, ku_phys);
    MultBx33(coeff, ku_phys, mu_ku);
    mu_ku[0] = -mu_ku[0];
    mu_ku[1] = -mu_ku[1];
    mu_ku[2] = -mu_ku[2];
    MultAtx33(J_loc, mu_ku, curl_v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
    curl_v[i + Q * 0] = wdetJ[i] * curl_v_loc[0];
    curl_v[i + Q * 1] = wdetJ[i] * curl_v_loc[1];
    curl_v[i + Q * 2] = wdetJ[i] * curl_v_loc[2];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurl_pml_floquet_cross_re_33)(void *__restrict__ ctx, CeedInt Q,
                                                      const CeedScalar *const *in,
                                                      CeedScalar *const *out)
{
  return f_apply_hcurl_pml_floquet_cross_33_impl<PML_PART_RE>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_floquet_cross_im_33)(void *__restrict__ ctx, CeedInt Q,
                                                      const CeedScalar *const *in,
                                                      CeedScalar *const *out)
{
  return f_apply_hcurl_pml_floquet_cross_33_impl<PML_PART_IM>(ctx, Q, in, out);
}

#endif  // PALACE_LIBCEED_HCURL_PML_33_QF_H
