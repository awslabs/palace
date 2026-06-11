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
// all times the surface measure. The context is a CeedIntScalar array: [0].second =
// scale0, [1].second = scale1, followed by an (optional, MS only) material property
// coefficient context. The "_1" variants evaluate the field from a single volume
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

// Single-sided interface dielectric. Inputs: qw, grad_x_f, attr_1, grad_x_1, u_1.
#define PALACE_SURF_EPR_1_QF(name, integrand)                                            \
  CEED_QFUNCTION(name)(void *__restrict__ ctx_, CeedInt Q, const CeedScalar *const *in,  \
                       CeedScalar *const *out)                                           \
  {                                                                                      \
    const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;                              \
    const CeedScalar *qw = in[0], *J_f = in[1], *attr = in[2], *J_v = in[3], *u = in[4]; \
    CeedScalar *v = out[0];                                                              \
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)                                       \
    {                                                                                    \
      CeedScalar J_f_loc[6], n[3], E[3];                                                 \
      MatUnpack32(J_f + i, Q, J_f_loc);                                                  \
      const CeedScalar wdetJ = qw[i] * SurfMeasure32(J_f_loc, n);                        \
      SurfHcurlField32(i, Q, J_v, u, E);                                                 \
      v[i] = wdetJ * (integrand);                                                        \
    }                                                                                    \
    return 0;                                                                            \
  }

// Two-sided (averaged) interface dielectric. Inputs: qw, grad_x_f, attr_1, grad_x_1,
// attr_2, grad_x_2, u_1, u_2.
#define PALACE_SURF_EPR_2_QF(name, integrand)                                           \
  CEED_QFUNCTION(name)(void *__restrict__ ctx_, CeedInt Q, const CeedScalar *const *in, \
                       CeedScalar *const *out)                                          \
  {                                                                                     \
    const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;                             \
    const CeedScalar *qw = in[0], *J_f = in[1], *attr = in[2], *J_v1 = in[3],           \
                     *J_v2 = in[5], *u_1 = in[6], *u_2 = in[7];                         \
    CeedScalar *v = out[0];                                                             \
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)                                      \
    {                                                                                   \
      CeedScalar J_f_loc[6], n[3], E[3];                                                \
      MatUnpack32(J_f + i, Q, J_f_loc);                                                 \
      const CeedScalar wdetJ = qw[i] * SurfMeasure32(J_f_loc, n);                       \
      SurfHcurlField2Avg32(i, Q, J_v1, J_v2, u_1, u_2, E);                              \
      v[i] = wdetJ * (integrand);                                                       \
    }                                                                                   \
    return 0;                                                                           \
  }

PALACE_SURF_EPR_1_QF(f_integ_surf_epr_def_1_32, SurfEprDefault(ctx, E))
PALACE_SURF_EPR_2_QF(f_integ_surf_epr_def_2_32, SurfEprDefault(ctx, E))
PALACE_SURF_EPR_1_QF(f_integ_surf_epr_ma_1_32, SurfEprMA(ctx, n, E))
PALACE_SURF_EPR_2_QF(f_integ_surf_epr_ma_2_32, SurfEprMA(ctx, n, E))
PALACE_SURF_EPR_1_QF(f_integ_surf_epr_ms_1_32, SurfEprMS(ctx, (CeedInt)attr[i], n, E))
PALACE_SURF_EPR_2_QF(f_integ_surf_epr_ms_2_32, SurfEprMS(ctx, (CeedInt)attr[i], n, E))
PALACE_SURF_EPR_1_QF(f_integ_surf_epr_sa_1_32, SurfEprSA(ctx, n, E))
PALACE_SURF_EPR_2_QF(f_integ_surf_epr_sa_2_32, SurfEprSA(ctx, n, E))

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

// Pointwise boundary field values (no quadrature weighting; for visualization output at
// the boundary element lattice points), following BdrFieldVectorCoefficient: the field
// from the attached volume element, averaged over both sides for interior boundaries.
// Inputs ("_1"): attr_1, grad_x_1, u_1; ("_2"): attr_1, grad_x_1, attr_2, grad_x_2,
// u_1, u_2. Output: 3 components per point.
#define PALACE_SURF_BDR_FIELD_QF(name1, name2, field_helper)                   \
  CEED_QFUNCTION(name1)(void *, CeedInt Q, const CeedScalar *const *in,        \
                        CeedScalar *const *out)                                \
  {                                                                            \
    const CeedScalar *J_v = in[1], *u = in[2];                                 \
    CeedScalar *v = out[0];                                                    \
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)                             \
    {                                                                          \
      CeedScalar V[3];                                                         \
      field_helper(i, Q, J_v, u, V);                                           \
      v[i + Q * 0] = V[0];                                                     \
      v[i + Q * 1] = V[1];                                                     \
      v[i + Q * 2] = V[2];                                                     \
    }                                                                          \
    return 0;                                                                  \
  }                                                                            \
  CEED_QFUNCTION(name2)(void *, CeedInt Q, const CeedScalar *const *in,        \
                        CeedScalar *const *out)                                \
  {                                                                            \
    const CeedScalar *J_v1 = in[1], *J_v2 = in[3], *u_1 = in[4], *u_2 = in[5]; \
    CeedScalar *v = out[0];                                                    \
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)                             \
    {                                                                          \
      CeedScalar V[3], V_2[3];                                                 \
      field_helper(i, Q, J_v1, u_1, V);                                        \
      field_helper(i, Q, J_v2, u_2, V_2);                                      \
      v[i + Q * 0] = 0.5 * (V[0] + V_2[0]);                                    \
      v[i + Q * 1] = 0.5 * (V[1] + V_2[1]);                                    \
      v[i + Q * 2] = 0.5 * (V[2] + V_2[2]);                                    \
    }                                                                          \
    return 0;                                                                  \
  }

PALACE_SURF_BDR_FIELD_QF(f_eval_bdr_hcurl_1_32, f_eval_bdr_hcurl_2_32, SurfHcurlField32)
PALACE_SURF_BDR_FIELD_QF(f_eval_bdr_hdiv_1_32, f_eval_bdr_hdiv_2_32, SurfHdivField32)

#endif  // PALACE_LIBCEED_SURF_32_QF_H
