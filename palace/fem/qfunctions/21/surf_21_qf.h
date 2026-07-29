// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SURF_21_QF_H
#define PALACE_LIBCEED_SURF_21_QF_H

#include "../22/utils_22_qf.h"
#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_2_qf.h"
#include "utils_21_qf.h"

// QFunctions for line output functionals (integrals of functions of fields over boundary
// elements of a 2D mesh). This mirrors the 3D surface functional kernels in surf_32_qf.h,
// but the boundary element is a segment embedded in 2D and H(curl) fields are 2-vectors.

// Line measure qw * |J_f| and unit normal from the raw boundary element Jacobian. This
// matches mfem::CalcOrtho for 2x1 Jacobians: n = (dy, -dx), then normalized.
CEED_QFUNCTION_HELPER CeedScalar SurfMeasure21(const CeedScalar J_f[2], CeedScalar n[2])
{
  n[0] = J_f[1];
  n[1] = -J_f[0];
  const CeedScalar detJ = sqrt(n[0] * n[0] + n[1] * n[1]);
  n[0] /= detJ;
  n[1] /= detJ;
  return detJ;
}

// H(curl) field at a point from the raw volume Jacobian: E = adj(J_v)^T/detJ_v u.
CEED_QFUNCTION_HELPER void SurfHcurlField21(CeedInt i, CeedInt Q, const CeedScalar *J_v,
                                            const CeedScalar *u, CeedScalar E[2])
{
  const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
  CeedScalar J_loc[4], adjJt_loc[4];
  MatUnpack22(J_v + i, Q, J_loc);
  const CeedScalar detJ = AdjJt22<true>(J_loc, adjJt_loc);
  MultBx22(adjJt_loc, u_loc, E);
  E[0] /= detJ;
  E[1] /= detJ;
}

// Two-sided H(curl) average: E = 1/2 (E_1 + E_2).
CEED_QFUNCTION_HELPER void
SurfHcurlField2Avg21(CeedInt i, CeedInt Q, const CeedScalar *J_v1, const CeedScalar *J_v2,
                     const CeedScalar *u_1, const CeedScalar *u_2, CeedScalar E[2])
{
  CeedScalar E_2[2];
  SurfHcurlField21(i, Q, J_v1, u_1, E);
  SurfHcurlField21(i, Q, J_v2, u_2, E_2);
  E[0] = 0.5 * (E[0] + E_2[0]);
  E[1] = 0.5 * (E[1] + E_2[1]);
}

// Runtime-selected 2D field value for boundary visualization. H(curl) fields are full
// in-plane physical vectors. The H(div) branch is not used for Palace's scalar 2D B
// boundary output, but keeps the shared enum path well-defined.
CEED_QFUNCTION_HELPER void SurfField21(CeedInt piola, CeedInt i, CeedInt Q,
                                       const CeedScalar *J_v, const CeedScalar *u,
                                       CeedScalar V[2])
{
  if (piola)
  {
    V[0] = IntegralMap22(i, Q, J_v, u);
    V[1] = 0.0;
  }
  else
  {
    SurfHcurlField21(i, Q, J_v, u, V);
  }
}

// Pointwise boundary vector field values (no quadrature weighting; for visualization
// output at the boundary element lattice points), following BdrFieldVectorCoefficient:
// the full field from the attached volume element, averaged over both sides for interior
// boundaries. Inputs match the 3D field kernels. Output: 2 components per point.
CEED_QFUNCTION(f_eval_bdr_field_1_21)(void *__restrict__ ctx_, CeedInt Q,
                                      const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_v = in[0], *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar V[2];
    SurfField21(ctx[0].first, i, Q, J_v, u, V);
    v[i + Q * 0] = ctx[1].second * V[0];
    v[i + Q * 1] = ctx[1].second * V[1];
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_field_2_21)(void *__restrict__ ctx_, CeedInt Q,
                                      const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_v1 = in[0], *J_v2 = in[1], *u_1 = in[2], *u_2 = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar V[2], V_2[2];
    SurfField21(ctx[0].first, i, Q, J_v1, u_1, V);
    SurfField21(ctx[0].first, i, Q, J_v2, u_2, V_2);
    v[i + Q * 0] = ctx[1].second * 0.5 * (V[0] + V_2[0]);
    v[i + Q * 1] = ctx[1].second * 0.5 * (V[1] + V_2[1]);
  }
  return 0;
}

// Pointwise scalar H1 boundary field values. Inputs: grad_x_1, u_1; or grad_x_1,
// grad_x_2, u_1, u_2. The geometry inputs keep the operator shape shared with other
// boundary field kernels; the scalar VALUE basis already evaluates the physical value.
CEED_QFUNCTION(f_eval_bdr_field_h1_1_21)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = ctx[1].second * u[i];
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_field_h1_2_21)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *u_1 = in[2], *u_2 = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = ctx[1].second * 0.5 * (u_1[i] + u_2[i]);
  }
  return 0;
}

// Line measure: v = qw * |J_f|. Inputs: qw, grad_x_f.
CEED_QFUNCTION(f_integ_surf_area_21)(void *, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *qw = in[0], *J_f = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[2], n[2];
    MatUnpack21(J_f + i, Q, J_f_loc);
    v[i] = qw[i] * SurfMeasure21(J_f_loc, n);
  }
  return 0;
}

// Squared L2 norm of an H(curl) field: v = qw * |J_f| * |E|^2. Inputs: qw,
// grad_x_f, attr_1, grad_x_1, u_1.
CEED_QFUNCTION(f_integ_surf_hcurl_norm2_21)(void *, CeedInt Q, const CeedScalar *const *in,
                                            CeedScalar *const *out)
{
  const CeedScalar *qw = in[0], *J_f = in[1], *J_v = in[3], *u = in[4];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[2], n[2], E[2];
    MatUnpack21(J_f + i, Q, J_f_loc);
    SurfHcurlField21(i, Q, J_v, u, E);
    v[i] = qw[i] * SurfMeasure21(J_f_loc, n) * (E[0] * E[0] + E[1] * E[1]);
  }
  return 0;
}

CEED_QFUNCTION_HELPER CeedScalar SurfEprDefault21(const CeedIntScalar *ctx,
                                                  const CeedScalar E[2])
{
  return ctx[0].second * (E[0] * E[0] + E[1] * E[1]);
}

CEED_QFUNCTION_HELPER CeedScalar SurfEprMA21(const CeedIntScalar *ctx,
                                             const CeedScalar n[2], const CeedScalar E[2])
{
  const CeedScalar En = n[0] * E[0] + n[1] * E[1];
  return ctx[0].second * En * En;
}

CEED_QFUNCTION_HELPER CeedScalar SurfEprMS21(const CeedIntScalar *ctx, CeedInt attr,
                                             const CeedScalar n[2], const CeedScalar E[2])
{
  CeedScalar eps[4], D[2];
  CoeffUnpack2(ctx + 2, attr, eps);
  MultBx22(eps, E, D);
  const CeedScalar Dn = n[0] * D[0] + n[1] * D[1];
  return ctx[0].second * Dn * Dn;
}

CEED_QFUNCTION_HELPER CeedScalar SurfEprSA21(const CeedIntScalar *ctx,
                                             const CeedScalar n[2], const CeedScalar E[2])
{
  const CeedScalar En = n[0] * E[0] + n[1] * E[1];
  const CeedScalar Et2 = E[0] * E[0] + E[1] * E[1] - En * En;
  return ctx[0].second * Et2 + ctx[1].second * En * En;
}

CEED_QFUNCTION_HELPER CeedScalar SurfEpr21(const CeedIntScalar *ctx, CeedInt attr,
                                           const CeedScalar n[2], const CeedScalar E[2])
{
  switch (ctx[0].first)
  {
    case 0:
      return SurfEprDefault21(ctx + 1, E);
    case 1:
      return SurfEprMA21(ctx + 1, n, E);
    case 2:
      return SurfEprMS21(ctx + 1, attr, n, E);
    default:  // 3: SA
      return SurfEprSA21(ctx + 1, n, E);
  }
}

// Single-sided interface dielectric. Inputs: qw, grad_x_f, attr_1, grad_x_1, u_1.
CEED_QFUNCTION(f_integ_surf_epr_1_21)(void *__restrict__ ctx_, CeedInt Q,
                                      const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *qw = in[0], *J_f = in[1], *attr = in[2], *J_v = in[3], *u = in[4];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[2], n[2], E[2];
    MatUnpack21(J_f + i, Q, J_f_loc);
    const CeedScalar wdetJ = qw[i] * SurfMeasure21(J_f_loc, n);
    SurfHcurlField21(i, Q, J_v, u, E);
    v[i] = wdetJ * SurfEpr21(ctx, (CeedInt)attr[i], n, E);
  }
  return 0;
}

// Two-sided (averaged) interface dielectric. Inputs: qw, grad_x_f, attr_1, grad_x_1,
// attr_2, grad_x_2, u_1, u_2.
CEED_QFUNCTION(f_integ_surf_epr_2_21)(void *__restrict__ ctx_, CeedInt Q,
                                      const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *qw = in[0], *J_f = in[1], *attr = in[2], *J_v1 = in[3], *J_v2 = in[5],
                   *u_1 = in[6], *u_2 = in[7];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[2], n[2], E[2];
    MatUnpack21(J_f + i, Q, J_f_loc);
    const CeedScalar wdetJ = qw[i] * SurfMeasure21(J_f_loc, n);
    SurfHcurlField2Avg21(i, Q, J_v1, J_v2, u_1, u_2, E);
    v[i] = wdetJ * SurfEpr21(ctx, (CeedInt)attr[i], n, E);
  }
  return 0;
}

// Pointwise boundary surface charge, surface current, energy density, and Poynting
// values at visualization lattice points for 2D line boundaries. These mirror the 3D
// boundary visualization kernels in surf_32_qf.h using Palace's 2D field convention:
// E is an in-plane H(curl) vector and B is the scalar out-of-plane magnetic flux B_z.
// Vector outputs have two in-plane components to preserve MFEM's 2D ParaView field
// shape. Context layouts match the 3D kernels, except material tables for magnetic
// quantities use the scalar 1x1 out-of-plane inverse permeability.
CEED_QFUNCTION(f_eval_bdr_flux_q_1_21)(void *__restrict__ ctx_, CeedInt Q,
                                       const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_f = in[0], *attr = in[1], *J_v = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[2], n[2], E[2], eps[4], D[2];
    MatUnpack21(J_f + i, Q, J_f_loc);
    SurfMeasure21(J_f_loc, n);
    SurfHcurlField21(i, Q, J_v, u, E);
    CoeffUnpack2(ctx + 2, (CeedInt)attr[i], eps);
    MultBx22(eps, E, D);
    v[i] = ctx[0].second * ctx[1].second * (D[0] * n[0] + D[1] * n[1]);
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_flux_q_2_21)(void *__restrict__ ctx_, CeedInt Q,
                                       const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_f = in[0], *attr_1 = in[1], *J_v1 = in[2], *attr_2 = in[3],
                   *J_v2 = in[4], *u_1 = in[5], *u_2 = in[6];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[2], n[2], E[2], eps[4], D[2], D_2[2];
    MatUnpack21(J_f + i, Q, J_f_loc);
    SurfMeasure21(J_f_loc, n);
    SurfHcurlField21(i, Q, J_v1, u_1, E);
    CoeffUnpack2(ctx + 2, (CeedInt)attr_1[i], eps);
    MultBx22(eps, E, D);
    SurfHcurlField21(i, Q, J_v2, u_2, E);
    CoeffUnpack2(ctx + 2, (CeedInt)attr_2[i], eps);
    MultBx22(eps, E, D_2);
    // Two-sided: contributions from opposite sides add with opposite normals.
    v[i] =
        ctx[0].second * ctx[1].second * ((D[0] - D_2[0]) * n[0] + (D[1] - D_2[1]) * n[1]);
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_current_j_1_21)(void *__restrict__ ctx_, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_f = in[0], *attr = in[1], *J_v = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[2], n[2];
    MatUnpack21(J_f + i, Q, J_f_loc);
    SurfMeasure21(J_f_loc, n);
    const CeedScalar Bz = IntegralMap22(i, Q, J_v, u);
    const CeedScalar Hz = CoeffUnpack1(ctx + 2, (CeedInt)attr[i]) * Bz;
    const CeedScalar s = ctx[0].second * ctx[1].second;
    v[i + Q * 0] = s * Hz * n[1];
    v[i + Q * 1] = -s * Hz * n[0];
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_current_j_2_21)(void *__restrict__ ctx_, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *J_f = in[0], *attr_1 = in[1], *J_v1 = in[2], *attr_2 = in[3],
                   *J_v2 = in[4], *u_1 = in[5], *u_2 = in[6];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_f_loc[2], n[2];
    MatUnpack21(J_f + i, Q, J_f_loc);
    SurfMeasure21(J_f_loc, n);
    const CeedScalar Bz_1 = IntegralMap22(i, Q, J_v1, u_1);
    const CeedScalar Bz_2 = IntegralMap22(i, Q, J_v2, u_2);
    const CeedScalar Hz_1 = CoeffUnpack1(ctx + 2, (CeedInt)attr_1[i]) * Bz_1;
    const CeedScalar Hz_2 = CoeffUnpack1(ctx + 2, (CeedInt)attr_2[i]) * Bz_2;
    const CeedScalar Hz = Hz_1 - Hz_2;
    const CeedScalar s = ctx[0].second * ctx[1].second;
    v[i + Q * 0] = s * Hz * n[1];
    v[i + Q * 1] = -s * Hz * n[0];
  }
  return 0;
}

CEED_QFUNCTION_HELPER void SurfPoynting21(const CeedIntScalar *ctx, CeedInt i, CeedInt Q,
                                          CeedInt attr, const CeedScalar *J_v,
                                          const CeedScalar *u_e, const CeedScalar *u_b,
                                          CeedScalar S[2])
{
  CeedScalar E[2];
  SurfHcurlField21(i, Q, J_v, u_e, E);
  const CeedScalar Bz = IntegralMap22(i, Q, J_v, u_b);
  const CeedScalar Hz = CoeffUnpack1(ctx + 1, attr) * Bz;
  S[0] = E[1] * Hz;
  S[1] = -E[0] * Hz;
}

CEED_QFUNCTION(f_eval_bdr_poynting_1_21)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *J_v = in[1], *u_e = in[2], *u_b = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar S[2];
    SurfPoynting21(ctx, i, Q, (CeedInt)attr[i], J_v, u_e, u_b, S);
    const CeedScalar s = ctx[0].second;
    v[i + Q * 0] = s * S[0];
    v[i + Q * 1] = s * S[1];
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_poynting_2_21)(void *__restrict__ ctx_, CeedInt Q,
                                         const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr_1 = in[0], *J_v1 = in[1], *attr_2 = in[2], *J_v2 = in[3],
                   *u_e_1 = in[4], *u_b_1 = in[5], *u_e_2 = in[6], *u_b_2 = in[7];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar S[2], S_2[2];
    SurfPoynting21(ctx, i, Q, (CeedInt)attr_1[i], J_v1, u_e_1, u_b_1, S);
    SurfPoynting21(ctx, i, Q, (CeedInt)attr_2[i], J_v2, u_e_2, u_b_2, S_2);
    const CeedScalar s = 0.5 * ctx[0].second;
    v[i + Q * 0] = s * (S[0] + S_2[0]);
    v[i + Q * 1] = s * (S[1] + S_2[1]);
  }
  return 0;
}

CEED_QFUNCTION_HELPER CeedScalar SurfEnergy21(const CeedIntScalar *ctx, CeedInt i,
                                              CeedInt Q, CeedInt attr,
                                              const CeedScalar *J_v, const CeedScalar *u)
{
  if (ctx[0].first)
  {
    const CeedScalar Bz = IntegralMap22(i, Q, J_v, u);
    return CoeffUnpack1(ctx + 2, attr) * Bz * Bz;
  }
  CeedScalar E[2], eps[4], D[2];
  SurfHcurlField21(i, Q, J_v, u, E);
  CoeffUnpack2(ctx + 2, attr, eps);
  MultBx22(eps, E, D);
  return D[0] * E[0] + D[1] * E[1];
}

CEED_QFUNCTION(f_eval_bdr_energy_1_21)(void *__restrict__ ctx_, CeedInt Q,
                                       const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *J_v = in[1], *u = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = 0.5 * ctx[1].second * SurfEnergy21(ctx, i, Q, (CeedInt)attr[i], J_v, u);
  }
  return 0;
}

CEED_QFUNCTION(f_eval_bdr_energy_2_21)(void *__restrict__ ctx_, CeedInt Q,
                                       const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr_1 = in[0], *J_v1 = in[1], *attr_2 = in[2], *J_v2 = in[3],
                   *u_1 = in[4], *u_2 = in[5];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar U_1 = SurfEnergy21(ctx, i, Q, (CeedInt)attr_1[i], J_v1, u_1);
    const CeedScalar U_2 = SurfEnergy21(ctx, i, Q, (CeedInt)attr_2[i], J_v2, u_2);
    v[i] = 0.25 * ctx[1].second * (U_1 + U_2);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_SURF_21_QF_H
