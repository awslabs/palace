// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SURF_32_QF_H
#define PALACE_LIBCEED_SURF_32_QF_H

#include "../33/utils_33_qf.h"
#include "../coeff/coeff_3_qf.h"
#include "utils_32_qf.h"

// QFunctions for surface output functionals (integrals of functions of fields over
// boundary elements of a 3D mesh). The element geometry is computed on the fly from
// mesh node gradients (as in geom_qf.h) instead of stored geometry factor data:
//  - "qw": quadrature weights (CEED_EVAL_WEIGHT).
//  - "grad_x_f": boundary element mesh node Jacobians J_f, shape [3x2, Q]. The surface
//    measure is qw * detJ_f and the surface normal is the normalized cross product of
//    the columns of J_f (matching mfem::CalcOrtho orientation: outward for exterior
//    boundaries, out of element 1 in general; the legacy boundary element orientation
//    flip is applied via the context normal sign).
//  - "attr_k" / "grad_x_k": attached volume element attributes (for material lookup)
//    and mesh node Jacobians J_v at the mapped face quadrature points, for the Piola
//    transformations H(curl): E = adj(J_v)ᵀ/detJ_v u, H(div): B = J_v/detJ_v u, for
//    evaluation side k = 1 (and 2 for two-sided interior boundary evaluation).
//  - "x": coordinates at the face quadrature points (surface flux only).
//  - "u_k": field inputs (reference components at the mapped face quadrature points).
// The output is the integrand times the surface measure, summed over quadrature points
// by the all-ones output element basis, yielding one value per boundary element.

// Surface measure qw * detJ_f and unit normal from the raw boundary element Jacobian.
CEED_QFUNCTION_HELPER CeedScalar SurfMeasure32(const CeedScalar J_f[6], CeedScalar n[3])
{
  n[0] = J_f[1] * J_f[5] - J_f[2] * J_f[4];
  n[1] = J_f[2] * J_f[3] - J_f[0] * J_f[5];
  n[2] = J_f[0] * J_f[4] - J_f[1] * J_f[3];
  const CeedScalar detJ = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  n[0] /= detJ;
  n[1] /= detJ;
  n[2] /= detJ;
  return detJ;
}

// H(curl) field at a point from the raw volume Jacobian: E = adj(J_v)ᵀ/detJ_v u.
CEED_QFUNCTION_HELPER void SurfHcurlField32(CeedInt i, CeedInt Q, const CeedScalar *J_v,
                                            const CeedScalar *u, CeedScalar E[3])
{
  const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
  CeedScalar J_loc[9], adjJt_loc[9];
  MatUnpack33(J_v + i, Q, J_loc);
  const CeedScalar detJ = AdjJt33<true>(J_loc, adjJt_loc);
  MultAx33(adjJt_loc, u_loc, E);
  E[0] /= detJ;
  E[1] /= detJ;
  E[2] /= detJ;
}

// H(div) field at a point from the raw volume Jacobian: B = J_v/detJ_v u.
CEED_QFUNCTION_HELPER void SurfHdivField32(CeedInt i, CeedInt Q, const CeedScalar *J_v,
                                           const CeedScalar *u, CeedScalar B[3])
{
  const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
  CeedScalar J_loc[9], adjJt_loc[9];
  MatUnpack33(J_v + i, Q, J_loc);
  const CeedScalar detJ = AdjJt33<true>(J_loc, adjJt_loc);
  MultAx33(J_loc, u_loc, B);
  B[0] /= detJ;
  B[1] /= detJ;
  B[2] /= detJ;
}

// Two-sided H(curl) average: E = 1/2 (E_1 + E_2).
CEED_QFUNCTION_HELPER void
SurfHcurlField2Avg32(CeedInt i, CeedInt Q, const CeedScalar *J_v1, const CeedScalar *J_v2,
                     const CeedScalar *u_1, const CeedScalar *u_2, CeedScalar E[3])
{
  CeedScalar E_2[3];
  SurfHcurlField32(i, Q, J_v1, u_1, E);
  SurfHcurlField32(i, Q, J_v2, u_2, E_2);
  E[0] = 0.5 * (E[0] + E_2[0]);
  E[1] = 0.5 * (E[1] + E_2[1]);
  E[2] = 0.5 * (E[2] + E_2[2]);
}

// Surface area: v = qw * detJ_f. Inputs: qw, grad_x_f.
CEED_QFUNCTION(f_integ_surf_area_32)(void *, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *qw = in[0], *J_f = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    v[i] = qw[i] * SurfMeasure32(J_f_loc, n);
  }
  return 0;
}

// Squared L2 norm of an H(curl) field: v = qw * detJ_f * |E|². Inputs: qw, grad_x_f,
// attr_1, grad_x_1, u_1.
CEED_QFUNCTION(f_integ_surf_hcurl_norm2_32)(void *, CeedInt Q, const CeedScalar *const *in,
                                            CeedScalar *const *out)
{
  const CeedScalar *qw = in[0], *J_f = in[1], *J_v = in[3], *u = in[4];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], E[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    SurfHcurlField32(i, Q, J_v, u, E);
    v[i] = qw[i] * SurfMeasure32(J_f_loc, n) * (E[0] * E[0] + E[1] * E[1] + E[2] * E[2]);
  }
  return 0;
}

// Interface dielectric energy participation integrands following
// InterfaceDielectricCoefficient (fem/coefficient.hpp):
//   DEFAULT: v = scale0 * |E|²,           scale0 = 0.5 * t * eps
//   MA:      v = scale0 * |n.E|²,         scale0 = 0.5 * t / eps
//   MS:      v = scale0 * |n.(eps_S E)|², scale0 = 0.5 * t / eps, eps_S from the
//            material table looked up with the side 1 volume attribute
//   SA:      v = scale0 * |E_t|² + scale1 * |E_n|², scale0 = 0.5 * t * eps,
//            scale1 = 0.5 * t / eps
// all times the surface measure. The context is a CeedIntScalar array: [0].first =
// epr_type (0 = DEFAULT, 1 = MA, 2 = MS, 3 = SA; selects the integrand at runtime, in its
// own slot since CeedIntScalar is a union), [1].second = scale0, [2].second = scale1,
// followed by an (optional, MS only) material property coefficient context. The shared
// kernel passes ctx + 1 to the per-type helpers. The "_1" variants evaluate from a single
// element; the "_2" variants average the fields from both sides of an interior
// boundary. The normal direction sign is irrelevant (quadratic in n . E).

CEED_QFUNCTION_HELPER CeedScalar SurfEprDefault(const CeedIntScalar *ctx,
                                                const CeedScalar E[3])
{
  return ctx[0].second * (E[0] * E[0] + E[1] * E[1] + E[2] * E[2]);
}

CEED_QFUNCTION_HELPER CeedScalar SurfEprMA(const CeedIntScalar *ctx, const CeedScalar n[3],
                                           const CeedScalar E[3])
{
  const CeedScalar En = n[0] * E[0] + n[1] * E[1] + n[2] * E[2];
  return ctx[0].second * En * En;
}

CEED_QFUNCTION_HELPER CeedScalar SurfEprMS(const CeedIntScalar *ctx, CeedInt attr,
                                           const CeedScalar n[3], const CeedScalar E[3])
{
  CeedScalar eps[9], D[3];
  CoeffUnpack3(ctx + 2, attr, eps);
  MultAx33(eps, E, D);
  const CeedScalar Dn = n[0] * D[0] + n[1] * D[1] + n[2] * D[2];
  return ctx[0].second * Dn * Dn;
}

CEED_QFUNCTION_HELPER CeedScalar SurfEprSA(const CeedIntScalar *ctx, const CeedScalar n[3],
                                           const CeedScalar E[3])
{
  const CeedScalar En = n[0] * E[0] + n[1] * E[1] + n[2] * E[2];
  const CeedScalar Et2 = E[0] * E[0] + E[1] * E[1] + E[2] * E[2] - En * En;
  return ctx[0].second * Et2 + ctx[1].second * En * En;
}

// Single-sided / two-sided interface dielectric. The four interface types (DEFAULT,
// MA, MS, SA) share one compiled kernel: the type is selected at runtime from
// ctx[0].first (0 = DEFAULT, 1 = MA, 2 = MS, 3 = SA). Because the type is uniform across
// the operator this is a divergence-free branch over the existing per-type helpers, so a
// single kernel replaces the per-type specializations at no runtime cost. The type sits
// in its own slot (CeedIntScalar is a union, so .first and .second of one slot alias);
// the per-type helpers see the rest of the context via ctx + 1 (scale0, scale1, ...).
CEED_QFUNCTION_HELPER CeedScalar SurfEpr(const CeedIntScalar *ctx, CeedInt attr,
                                         const CeedScalar n[3], const CeedScalar E[3])
{
  switch (ctx[0].first)
  {
    case 0:
      return SurfEprDefault(ctx + 1, E);
    case 1:
      return SurfEprMA(ctx + 1, n, E);
    case 2:
      return SurfEprMS(ctx + 1, attr, n, E);
    default:  // 3: SA
      return SurfEprSA(ctx + 1, n, E);
  }
}

// Single-sided interface dielectric. Inputs: qw, grad_x_f, attr_1, grad_x_1, u_1.
CEED_QFUNCTION(f_integ_surf_epr_1_32)(void *__restrict__ ctx_, CeedInt Q,
                                      const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *qw = in[0], *J_f = in[1], *attr = in[2], *J_v = in[3], *u = in[4];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], E[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    const CeedScalar wdetJ = qw[i] * SurfMeasure32(J_f_loc, n);
    SurfHcurlField32(i, Q, J_v, u, E);
    v[i] = wdetJ * SurfEpr(ctx, (CeedInt)attr[i], n, E);
  }
  return 0;
}

// Two-sided (averaged) interface dielectric. Inputs: qw, grad_x_f, attr_1, grad_x_1,
// attr_2, grad_x_2, u_1, u_2.
CEED_QFUNCTION(f_integ_surf_epr_2_32)(void *__restrict__ ctx_, CeedInt Q,
                                      const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *qw = in[0], *J_f = in[1], *attr = in[2], *J_v1 = in[3], *J_v2 = in[5],
                   *u_1 = in[6], *u_2 = in[7];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], E[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    const CeedScalar wdetJ = qw[i] * SurfMeasure32(J_f_loc, n);
    SurfHcurlField2Avg32(i, Q, J_v1, J_v2, u_1, u_2, E);
    v[i] = wdetJ * SurfEpr(ctx, (CeedInt)attr[i], n, E);
  }
  return 0;
}

// Surface flux integrands following BdrSurfaceFluxCoefficient (fem/coefficient.hpp):
//   ELECTRIC: V = eps E, MAGNETIC: V = B, POWER: V = E x (mu⁻¹ B),
// with v = V . n times the surface measure. The context is a CeedIntScalar array:
// [0].second = normal sign (legacy boundary element orientation flip), [1].first =
// two_sided flag, [2..4].second = x0, followed by the material property context
// (ELECTRIC: permittivity, POWER: inverse permeability). Two-sided contributions add
// with opposite normals (V_1 - V_2); otherwise interior values average and the flux
// sign is chosen per point so the normal points away from x0. The "_2" variants take
// both sides. Inputs: qw, grad_x_f, attr_1, grad_x_1, [attr_2, grad_x_2,] x, u_k...

CEED_QFUNCTION_HELPER CeedScalar SurfFluxFinalize(const CeedIntScalar *ctx, CeedInt i,
                                                  CeedInt Q, CeedScalar n[3],
                                                  const CeedScalar *x,
                                                  const CeedScalar V[3])
{
  const CeedScalar s_n = ctx[0].second;
  n[0] *= s_n;
  n[1] *= s_n;
  n[2] *= s_n;
  CeedScalar flux = V[0] * n[0] + V[1] * n[1] + V[2] * n[2];
  if (!ctx[1].first)
  {
    // Orient outward from the surface with the given center.
    const CeedScalar dx[3] = {x[i + Q * 0] - ctx[2].second, x[i + Q * 1] - ctx[3].second,
                              x[i + Q * 2] - ctx[4].second};
    if (dx[0] * n[0] + dx[1] * n[1] + dx[2] * n[2] < 0.0)
    {
      flux = -flux;
    }
  }
  return flux;
}

CEED_QFUNCTION_HELPER void SurfFluxCombine(const CeedIntScalar *ctx,
                                           const CeedScalar V_2[3], CeedScalar V[3])
{
  if (ctx[1].first)
  {
    // Two-sided: add the contributions from opposite sides with opposite normals.
    V[0] -= V_2[0];
    V[1] -= V_2[1];
    V[2] -= V_2[2];
  }
  else
  {
    // Average the values from both sides.
    V[0] = 0.5 * (V[0] + V_2[0]);
    V[1] = 0.5 * (V[1] + V_2[1]);
    V[2] = 0.5 * (V[2] + V_2[2]);
  }
}

CEED_QFUNCTION_HELPER void SurfFluxElectric(const CeedIntScalar *ctx, CeedInt i, CeedInt Q,
                                            CeedInt attr, const CeedScalar *J_v,
                                            const CeedScalar *u, CeedScalar V[3])
{
  CeedScalar E[3], eps[9];
  SurfHcurlField32(i, Q, J_v, u, E);
  CoeffUnpack3(ctx + 5, attr, eps);
  MultAx33(eps, E, V);
}

CEED_QFUNCTION_HELPER void SurfFluxPoynting(const CeedIntScalar *ctx, CeedInt i, CeedInt Q,
                                            CeedInt attr, const CeedScalar *J_v,
                                            const CeedScalar *u_e, const CeedScalar *u_b,
                                            CeedScalar V[3])
{
  CeedScalar E[3], B[3], invmu[9], H[3];
  SurfHcurlField32(i, Q, J_v, u_e, E);
  SurfHdivField32(i, Q, J_v, u_b, B);
  CoeffUnpack3(ctx + 5, attr, invmu);
  MultAx33(invmu, B, H);
  V[0] = E[1] * H[2] - E[2] * H[1];
  V[1] = E[2] * H[0] - E[0] * H[2];
  V[2] = E[0] * H[1] - E[1] * H[0];
}

CEED_QFUNCTION(f_integ_surf_flux_e_1_32)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *qw = in[0], *J_f = in[1], *attr = in[2], *J_v = in[3], *x = in[4],
                   *u = in[5];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], V[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    const CeedScalar wdetJ = qw[i] * SurfMeasure32(J_f_loc, n);
    SurfFluxElectric(ctx, i, Q, (CeedInt)attr[i], J_v, u, V);
    v[i] = wdetJ * SurfFluxFinalize(ctx, i, Q, n, x, V);
  }
  return 0;
}

CEED_QFUNCTION(f_integ_surf_flux_e_2_32)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *qw = in[0], *J_f = in[1], *attr_1 = in[2], *J_v1 = in[3],
                   *attr_2 = in[4], *J_v2 = in[5], *x = in[6], *u_1 = in[7], *u_2 = in[8];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], V[3], V_2[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    const CeedScalar wdetJ = qw[i] * SurfMeasure32(J_f_loc, n);
    SurfFluxElectric(ctx, i, Q, (CeedInt)attr_1[i], J_v1, u_1, V);
    SurfFluxElectric(ctx, i, Q, (CeedInt)attr_2[i], J_v2, u_2, V_2);
    SurfFluxCombine(ctx, V_2, V);
    v[i] = wdetJ * SurfFluxFinalize(ctx, i, Q, n, x, V);
  }
  return 0;
}

CEED_QFUNCTION(f_integ_surf_flux_m_1_32)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *qw = in[0], *J_f = in[1], *J_v = in[3], *x = in[4], *u = in[5];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], B[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    const CeedScalar wdetJ = qw[i] * SurfMeasure32(J_f_loc, n);
    SurfHdivField32(i, Q, J_v, u, B);
    v[i] = wdetJ * SurfFluxFinalize(ctx, i, Q, n, x, B);
  }
  return 0;
}

CEED_QFUNCTION(f_integ_surf_flux_m_2_32)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *qw = in[0], *J_f = in[1], *J_v1 = in[3], *J_v2 = in[5], *x = in[6],
                   *u_1 = in[7], *u_2 = in[8];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], B[3], B_2[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    const CeedScalar wdetJ = qw[i] * SurfMeasure32(J_f_loc, n);
    SurfHdivField32(i, Q, J_v1, u_1, B);
    SurfHdivField32(i, Q, J_v2, u_2, B_2);
    SurfFluxCombine(ctx, B_2, B);
    v[i] = wdetJ * SurfFluxFinalize(ctx, i, Q, n, x, B);
  }
  return 0;
}

CEED_QFUNCTION(f_integ_surf_flux_p_1_32)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *qw = in[0], *J_f = in[1], *attr = in[2], *J_v = in[3], *x = in[4],
                   *u_e = in[5], *u_b = in[6];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], S[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    const CeedScalar wdetJ = qw[i] * SurfMeasure32(J_f_loc, n);
    SurfFluxPoynting(ctx, i, Q, (CeedInt)attr[i], J_v, u_e, u_b, S);
    v[i] = wdetJ * SurfFluxFinalize(ctx, i, Q, n, x, S);
  }
  return 0;
}

CEED_QFUNCTION(f_integ_surf_flux_p_2_32)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *qw = in[0], *J_f = in[1], *attr_1 = in[2], *J_v1 = in[3],
                   *attr_2 = in[4], *J_v2 = in[5], *x = in[6], *u_e_1 = in[7],
                   *u_b_1 = in[8], *u_e_2 = in[9], *u_b_2 = in[10];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], S[3], S_2[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    const CeedScalar wdetJ = qw[i] * SurfMeasure32(J_f_loc, n);
    SurfFluxPoynting(ctx, i, Q, (CeedInt)attr_1[i], J_v1, u_e_1, u_b_1, S);
    SurfFluxPoynting(ctx, i, Q, (CeedInt)attr_2[i], J_v2, u_e_2, u_b_2, S_2);
    SurfFluxCombine(ctx, S_2, S);
    v[i] = wdetJ * SurfFluxFinalize(ctx, i, Q, n, x, S);
  }
  return 0;
}

// Stratton-Chu far-field integrand following AddStrattonChuIntegrandAtElement
// (models/strattonchu.cpp): for each observation direction r0,
//   I = (ik/4pi) [n x E - r0 x (n x ZH)] e^{ik r0.r'} dS,  ZH = c0 B, k = omega/c0,
// with complex omega and fields. The context is a CeedIntScalar array: [0].second =
// normal sign, [1].second = omega_re, [2].second = omega_im, [3].first = number of
// directions N, [4..4+3N).second = directions, followed by the (isotropic) light speed
// material property context. Output: 6N components per element (Re I, Im I per
// direction, 3 each). Inputs: qw, grad_x_f, attr_1, grad_x_1, x, u_1..u_4 =
// E_re, E_im, B_re, B_im (external boundaries only).
CEED_QFUNCTION(f_integ_surf_farfield_32)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *qw = in[0], *J_f = in[1], *attr = in[2], *J_v = in[3], *x = in[4],
                   *u_er = in[5], *u_ei = in[6], *u_br = in[7], *u_bi = in[8];
  CeedScalar *v = out[0];
  const CeedScalar s_n = ctx[0].second, omega_re = ctx[1].second, omega_im = ctx[2].second;
  const CeedInt N = ctx[3].first;
  const CeedIntScalar *dirs = ctx + 4, *mat = ctx + 4 + 3 * N;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], E_r[3], E_i[3], B_r[3], B_i[3], c0_mat[9], ZH_r[3],
        ZH_i[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    const CeedScalar wdetJ = qw[i] * SurfMeasure32(J_f_loc, n);
    n[0] *= s_n;
    n[1] *= s_n;
    n[2] *= s_n;
    SurfHcurlField32(i, Q, J_v, u_er, E_r);
    SurfHcurlField32(i, Q, J_v, u_ei, E_i);
    SurfHdivField32(i, Q, J_v, u_br, B_r);
    SurfHdivField32(i, Q, J_v, u_bi, B_i);
    CoeffUnpack3(mat, (CeedInt)attr[i], c0_mat);
    MultAx33(c0_mat, B_r, ZH_r);
    MultAx33(c0_mat, B_i, ZH_i);
    const CeedScalar c0 = c0_mat[0];  // Isotropic (verified at setup)
    const CeedScalar k_re = omega_re / c0, k_im = omega_im / c0;
    const CeedScalar pre_re = -wdetJ * k_im / 12.566370614359172;
    const CeedScalar pre_im = wdetJ * k_re / 12.566370614359172;

    // n x E and n x ZH (real and imaginary parts).
    const CeedScalar nxEr[3] = {n[1] * E_r[2] - n[2] * E_r[1],
                                n[2] * E_r[0] - n[0] * E_r[2],
                                n[0] * E_r[1] - n[1] * E_r[0]};
    const CeedScalar nxEi[3] = {n[1] * E_i[2] - n[2] * E_i[1],
                                n[2] * E_i[0] - n[0] * E_i[2],
                                n[0] * E_i[1] - n[1] * E_i[0]};
    const CeedScalar nxZHr[3] = {n[1] * ZH_r[2] - n[2] * ZH_r[1],
                                 n[2] * ZH_r[0] - n[0] * ZH_r[2],
                                 n[0] * ZH_r[1] - n[1] * ZH_r[0]};
    const CeedScalar nxZHi[3] = {n[1] * ZH_i[2] - n[2] * ZH_i[1],
                                 n[2] * ZH_i[0] - n[0] * ZH_i[2],
                                 n[0] * ZH_i[1] - n[1] * ZH_i[0]};
    const CeedScalar xp[3] = {x[i + Q * 0], x[i + Q * 1], x[i + Q * 2]};

    for (CeedInt d = 0; d < N; d++)
    {
      const CeedScalar r0 = dirs[3 * d].second, r1 = dirs[3 * d + 1].second,
                       r2 = dirs[3 * d + 2].second;
      const CeedScalar dot = r0 * xp[0] + r1 * xp[1] + r2 * xp[2];
      const CeedScalar amp = exp(-k_im * dot);
      const CeedScalar cph = cos(k_re * dot), sph = sin(k_re * dot);
      const CeedScalar w_re = amp * (pre_re * cph - pre_im * sph);
      const CeedScalar w_im = amp * (pre_re * sph + pre_im * cph);

      // A + iB = n x E - r0 x (n x ZH).
      const CeedScalar A[3] = {nxEr[0] - (r1 * nxZHr[2] - r2 * nxZHr[1]),
                               nxEr[1] - (r2 * nxZHr[0] - r0 * nxZHr[2]),
                               nxEr[2] - (r0 * nxZHr[1] - r1 * nxZHr[0])};
      const CeedScalar Bv[3] = {nxEi[0] - (r1 * nxZHi[2] - r2 * nxZHi[1]),
                                nxEi[1] - (r2 * nxZHi[0] - r0 * nxZHi[2]),
                                nxEi[2] - (r0 * nxZHi[1] - r1 * nxZHi[0])};
      for (CeedInt c = 0; c < 3; c++)
      {
        v[i + Q * (6 * d + c)] = A[c] * w_re - Bv[c] * w_im;
        v[i + Q * (6 * d + 3 + c)] = A[c] * w_im + Bv[c] * w_re;
      }
    }
  }
  return 0;
}

// Combined Piola transform selected at runtime by piola: H(curl) (0, E = adj(J)^T u /
// detJ) or H(div) (1, B = J u / detJ). Lets the field and energy boundary-viz kernels
// share one compiled kernel across the ND and RT field spaces.
CEED_QFUNCTION_HELPER void SurfField32(CeedInt piola, CeedInt i, CeedInt Q,
                                       const CeedScalar *J_v, const CeedScalar *u,
                                       CeedScalar V[3])
{
  if (piola)
  {
    SurfHdivField32(i, Q, J_v, u, V);
  }
  else
  {
    SurfHcurlField32(i, Q, J_v, u, V);
  }
}

// Pointwise boundary field values (no quadrature weighting; for visualization output at
// the boundary element lattice points), following BdrFieldVectorCoefficient: the field
// from the attached volume element, averaged over both sides for interior boundaries.
// The H(curl)/H(div) space is selected at runtime from ctx[0].first, so the E (ND) and
// B (RT) boundary fields share one kernel. Inputs ("_1"): attr_1, grad_x_1, u_1; ("_2"):
// attr_1, grad_x_1, attr_2, grad_x_2, u_1, u_2. Output: 3 components per point.
CEED_QFUNCTION(f_eval_bdr_field_1_32)(void *__restrict__ ctx_, CeedInt Q,
                                      const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_v = in[1], *u = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar V[3];
    SurfField32(ctx[0].first, i, Q, J_v, u, V);
    v[i + Q * 0] = ctx[1].second * V[0];
    v[i + Q * 1] = ctx[1].second * V[1];
    v[i + Q * 2] = ctx[1].second * V[2];
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_field_2_32)(void *__restrict__ ctx_, CeedInt Q,
                                      const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_v1 = in[1], *J_v2 = in[3], *u_1 = in[4], *u_2 = in[5];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar V[3], V_2[3];
    SurfField32(ctx[0].first, i, Q, J_v1, u_1, V);
    SurfField32(ctx[0].first, i, Q, J_v2, u_2, V_2);
    v[i + Q * 0] = ctx[1].second * 0.5 * (V[0] + V_2[0]);
    v[i + Q * 1] = ctx[1].second * 0.5 * (V[1] + V_2[1]);
    v[i + Q * 2] = ctx[1].second * 0.5 * (V[2] + V_2[2]);
  }
  return 0;
}

// Pointwise boundary surface charge, surface current, and energy density values at the
// visualization lattice points, following BdrSurfaceFluxCoefficient<ELECTRIC>
// (two-sided), BdrSurfaceCurrentVectorCoefficient, and EnergyDensityCoefficient
// (boundary branch: per-side energies averaged). Context: [0].second = normal sign,
// [1].second = scaling, material table at +2. Inputs ("_1"): grad_x_f, attr_1,
// grad_x_1, u_1; ("_2"): grad_x_f, attr_1, grad_x_1, attr_2, grad_x_2, u_1, u_2.
CEED_QFUNCTION(f_eval_bdr_flux_q_1_32)(void *__restrict__ ctx_, CeedInt Q,
                                       const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_f = in[0], *attr = in[1], *J_v = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], E[3], eps[9], D[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    SurfMeasure32(J_f_loc, n);
    SurfHcurlField32(i, Q, J_v, u, E);
    CoeffUnpack3(ctx + 2, (CeedInt)attr[i], eps);
    MultAx33(eps, E, D);
    v[i] = ctx[0].second * ctx[1].second * (D[0] * n[0] + D[1] * n[1] + D[2] * n[2]);
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_flux_q_2_32)(void *__restrict__ ctx_, CeedInt Q,
                                       const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_f = in[0], *attr_1 = in[1], *J_v1 = in[2], *attr_2 = in[3],
                   *J_v2 = in[4], *u_1 = in[5], *u_2 = in[6];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], E[3], eps[9], D[3], D_2[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    SurfMeasure32(J_f_loc, n);
    SurfHcurlField32(i, Q, J_v1, u_1, E);
    CoeffUnpack3(ctx + 2, (CeedInt)attr_1[i], eps);
    MultAx33(eps, E, D);
    SurfHcurlField32(i, Q, J_v2, u_2, E);
    CoeffUnpack3(ctx + 2, (CeedInt)attr_2[i], eps);
    MultAx33(eps, E, D_2);
    // Two-sided: contributions from opposite sides add with opposite normals.
    v[i] = ctx[0].second * ctx[1].second *
           ((D[0] - D_2[0]) * n[0] + (D[1] - D_2[1]) * n[1] + (D[2] - D_2[2]) * n[2]);
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_current_j_1_32)(void *__restrict__ ctx_, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_f = in[0], *attr = in[1], *J_v = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], B[3], invmu[9], H[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    SurfMeasure32(J_f_loc, n);
    SurfHdivField32(i, Q, J_v, u, B);
    CoeffUnpack3(ctx + 2, (CeedInt)attr[i], invmu);
    MultAx33(invmu, B, H);
    const CeedScalar s = ctx[0].second * ctx[1].second;
    v[i + Q * 0] = s * (n[1] * H[2] - n[2] * H[1]);
    v[i + Q * 1] = s * (n[2] * H[0] - n[0] * H[2]);
    v[i + Q * 2] = s * (n[0] * H[1] - n[1] * H[0]);
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_current_j_2_32)(void *__restrict__ ctx_, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_f = in[0], *attr_1 = in[1], *J_v1 = in[2], *attr_2 = in[3],
                   *J_v2 = in[4], *u_1 = in[5], *u_2 = in[6];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[6], n[3], B[3], invmu[9], H[3], H_2[3];
    MatUnpack32(J_f + i, Q, J_f_loc);
    SurfMeasure32(J_f_loc, n);
    SurfHdivField32(i, Q, J_v1, u_1, B);
    CoeffUnpack3(ctx + 2, (CeedInt)attr_1[i], invmu);
    MultAx33(invmu, B, H);
    SurfHdivField32(i, Q, J_v2, u_2, B);
    CoeffUnpack3(ctx + 2, (CeedInt)attr_2[i], invmu);
    MultAx33(invmu, B, H_2);
    H[0] -= H_2[0];
    H[1] -= H_2[1];
    H[2] -= H_2[2];
    const CeedScalar s = ctx[0].second * ctx[1].second;
    v[i + Q * 0] = s * (n[1] * H[2] - n[2] * H[1]);
    v[i + Q * 1] = s * (n[2] * H[0] - n[0] * H[2]);
    v[i + Q * 2] = s * (n[0] * H[1] - n[1] * H[0]);
  }
  return 0;
}

// Boundary Poynting vector: per-side S = scale * E x (mu^-1 B), averaged over both
// sides for interior boundaries. This matches PoyntingVectorCoefficient's boundary
// branch: full vector, no n-dot, no 1/2 time-average factor. Context: [0].second =
// scaling, material table at +1. Inputs ("_1"): grad_x_f, attr_1, grad_x_1, u_e_1,
// u_b_1; ("_2"): grad_x_f, attr_1, grad_x_1, attr_2, grad_x_2, u_e_1, u_b_1,
// u_e_2, u_b_2.
CEED_QFUNCTION_HELPER void SurfPoynting32(const CeedIntScalar *ctx, CeedInt i, CeedInt Q,
                                          CeedInt attr, const CeedScalar *J_v,
                                          const CeedScalar *u_e, const CeedScalar *u_b,
                                          CeedScalar S[3])
{
  CeedScalar E[3], B[3], invmu[9], H[3];
  SurfHcurlField32(i, Q, J_v, u_e, E);
  SurfHdivField32(i, Q, J_v, u_b, B);
  CoeffUnpack3(ctx + 1, attr, invmu);
  MultAx33(invmu, B, H);
  S[0] = E[1] * H[2] - E[2] * H[1];
  S[1] = E[2] * H[0] - E[0] * H[2];
  S[2] = E[0] * H[1] - E[1] * H[0];
}

CEED_QFUNCTION(f_eval_bdr_poynting_1_32)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[1], *J_v = in[2], *u_e = in[3], *u_b = in[4];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar S[3];
    SurfPoynting32(ctx, i, Q, (CeedInt)attr[i], J_v, u_e, u_b, S);
    const CeedScalar s = ctx[0].second;
    v[i + Q * 0] = s * S[0];
    v[i + Q * 1] = s * S[1];
    v[i + Q * 2] = s * S[2];
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_poynting_2_32)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr_1 = in[1], *J_v1 = in[2], *attr_2 = in[3], *J_v2 = in[4],
                   *u_e_1 = in[5], *u_b_1 = in[6], *u_e_2 = in[7], *u_b_2 = in[8];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar S[3], S_2[3];
    SurfPoynting32(ctx, i, Q, (CeedInt)attr_1[i], J_v1, u_e_1, u_b_1, S);
    SurfPoynting32(ctx, i, Q, (CeedInt)attr_2[i], J_v2, u_e_2, u_b_2, S_2);
    const CeedScalar s = 0.5 * ctx[0].second;
    v[i + Q * 0] = s * (S[0] + S_2[0]);
    v[i + Q * 1] = s * (S[1] + S_2[1]);
    v[i + Q * 2] = s * (S[2] + S_2[2]);
  }
  return 0;
}

// Boundary energy densities: per-side 1/2 (mat F).F with the side attribute, averaged
// over both sides for interior boundaries. The H(curl)/H(div) space is selected at
// runtime from ctx[0].first; ctx[1].second = scaling, material table at +2. Inputs
// ("_1"): grad_x_f, attr_1, grad_x_1, u_1; ("_2"): grad_x_f, attr_1, grad_x_1, attr_2,
// grad_x_2, u_1, u_2.
CEED_QFUNCTION(f_eval_bdr_energy_1_32)(void *__restrict__ ctx_, CeedInt Q,
                                       const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[1], *J_v = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar F[3], mat[9], G[3];
    SurfField32(ctx[0].first, i, Q, J_v, u, F);
    CoeffUnpack3(ctx + 2, (CeedInt)attr[i], mat);
    MultAx33(mat, F, G);
    v[i] = 0.5 * ctx[1].second * (G[0] * F[0] + G[1] * F[1] + G[2] * F[2]);
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_energy_2_32)(void *__restrict__ ctx_, CeedInt Q,
                                       const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr_1 = in[1], *J_v1 = in[2], *attr_2 = in[3], *J_v2 = in[4],
                   *u_1 = in[5], *u_2 = in[6];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar F[3], mat[9], G[3];
    SurfField32(ctx[0].first, i, Q, J_v1, u_1, F);
    CoeffUnpack3(ctx + 2, (CeedInt)attr_1[i], mat);
    MultAx33(mat, F, G);
    CeedScalar U = G[0] * F[0] + G[1] * F[1] + G[2] * F[2];
    SurfField32(ctx[0].first, i, Q, J_v2, u_2, F);
    CoeffUnpack3(ctx + 2, (CeedInt)attr_2[i], mat);
    MultAx33(mat, F, G);
    U = 0.5 * (U + G[0] * F[0] + G[1] * F[1] + G[2] * F[2]);
    v[i] = 0.5 * ctx[1].second * U;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_SURF_32_QF_H
