// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SURF_21_QF_H
#define PALACE_LIBCEED_SURF_21_QF_H

#include "../22/utils_22_qf.h"
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

#endif  // PALACE_LIBCEED_SURF_21_QF_H
