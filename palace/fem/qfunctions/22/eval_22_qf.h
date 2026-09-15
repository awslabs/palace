// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_EVAL_22_QF_H
#define PALACE_LIBCEED_EVAL_22_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_2_qf.h"
#include "utils_22_qf.h"

// QFunctions for pointwise evaluation of derived field quantities at arbitrary points
// of 2D volume elements. No quadrature weighting is applied. H(curl):
// E = adj(J)^T / detJ u. The 2D magnetic flux density B is scalar-valued (out of
// plane), with H_z = mu^{-1}_{zz} B_z supplied through a scalar coefficient context.

// Electric energy density: v = 0.5 * scale * (eps_r E) . E.
CEED_QFUNCTION(f_eval_energy_e_22)(void *__restrict__ ctx_, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *J = in[1], *u = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar J_loc[4], adjJt_loc[4], E[2], eps[4], D[2];
    MatUnpack22(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt22<true>(J_loc, adjJt_loc);
    MultBx22(adjJt_loc, u_loc, E);
    CoeffUnpack2(ctx + 1, (CeedInt)attr[i], eps);
    MultBx22(eps, E, D);
    v[i] = 0.5 * ctx[0].second * (D[0] * E[0] + D[1] * E[1]) / (detJ * detJ);
  }
  return 0;
}

// Magnetic energy density: v = 0.5 * scale * mu^{-1}_{zz} B_z^2.
CEED_QFUNCTION(f_eval_energy_m_22)(void *__restrict__ ctx_, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *J = in[1], *u = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar invmu_zz = CoeffUnpack1(ctx + 1, (CeedInt)attr[i]);
    const CeedScalar B_z = IntegralMap22(i, Q, J, u);
    v[i] = 0.5 * ctx[0].second * invmu_zz * B_z * B_z;
  }
  return 0;
}

// Poynting vector: v = scale * E x (H_z zhat) = scale * H_z * (E_y, -E_x).
CEED_QFUNCTION(f_eval_poynting_22)(void *__restrict__ ctx_, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *J = in[1], *u_e = in[2], *u_b = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar ue_loc[2] = {u_e[i + Q * 0], u_e[i + Q * 1]};
    CeedScalar J_loc[4], adjJt_loc[4], E[2];
    MatUnpack22(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt22<true>(J_loc, adjJt_loc);
    MultBx22(adjJt_loc, ue_loc, E);
    const CeedScalar B_z = IntegralMap22(i, Q, J, u_b);
    const CeedScalar H_z = CoeffUnpack1(ctx + 1, (CeedInt)attr[i]) * B_z;
    const CeedScalar s = ctx[0].second * H_z / detJ;
    v[i + Q * 0] = s * E[1];
    v[i + Q * 1] = -s * E[0];
  }
  return 0;
}

// Boundary-mode normal Poynting density on a 2D cross-section, matching
// ModeSnCoefficient's sign convention: v = scale * Re{-Ex Hy* + Ey Hx*}, with
// Ht = mu^{-1}_{zz} Bt. This qfunction is applied separately to real and imaginary
// parts and accumulated.
CEED_QFUNCTION(f_eval_mode_sn_22)(void *__restrict__ ctx_, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *J = in[1], *u_e = in[2], *u_b = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar ue_loc[2] = {u_e[i + Q * 0], u_e[i + Q * 1]};
    const CeedScalar ub_loc[2] = {u_b[i + Q * 0], u_b[i + Q * 1]};
    CeedScalar J_loc[4], adjJt_loc[4], E[2], B[2];
    MatUnpack22(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt22<true>(J_loc, adjJt_loc);
    MultBx22(adjJt_loc, ue_loc, E);
    MultBx22(adjJt_loc, ub_loc, B);
    const CeedScalar invmu = CoeffUnpack1(ctx + 1, (CeedInt)attr[i]);
    const CeedScalar s = ctx[0].second * invmu / (detJ * detJ);
    v[i] = s * (-E[0] * B[1] + E[1] * B[0]);
  }
  return 0;
}

// Pointwise H(curl) field value: v = adj(J)^T/detJ u. Inputs: grad_x, u.
CEED_QFUNCTION(f_eval_probe_hcurl_22)(void *, CeedInt Q, const CeedScalar *const *in,
                                      CeedScalar *const *out)
{
  const CeedScalar *J = in[0], *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar J_loc[4], adjJt_loc[4], E[2];
    MatUnpack22(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt22<true>(J_loc, adjJt_loc);
    MultBx22(adjJt_loc, u_loc, E);
    v[i + Q * 0] = E[0] / detJ;
    v[i + Q * 1] = E[1] / detJ;
  }
  return 0;
}

// Pointwise scalar L2 field value. Inputs: grad_x, u. The geometry input is unused but
// kept so this probe has the same field list as the vector-valued domain probes.
CEED_QFUNCTION(f_eval_probe_l2_22)(void *, CeedInt Q, const CeedScalar *const *in,
                                   CeedScalar *const *out)
{
  const CeedScalar *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = u[i];
  }
  return 0;
}

// Pointwise scalar INTEGRAL field value: v = u / |detJ|. Inputs: grad_x, u.
CEED_QFUNCTION(f_eval_probe_integral_22)(void *, CeedInt Q, const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  const CeedScalar *J = in[0], *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = IntegralMap22(i, Q, J, u);
  }
  return 0;
}

// Pointwise H(div) field value: v = J/detJ u. Inputs: grad_x, u.
CEED_QFUNCTION(f_eval_probe_hdiv_22)(void *, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *J = in[0], *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar J_loc[4], B[2];
    MatUnpack22(J + i, Q, J_loc);
    const CeedScalar detJ = DetJ22(J_loc);
    MultBx22(J_loc, u_loc, B);
    v[i + Q * 0] = B[0] / detJ;
    v[i + Q * 1] = B[1] / detJ;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_EVAL_22_QF_H
