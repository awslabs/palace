// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_EVAL_33_QF_H
#define PALACE_LIBCEED_EVAL_33_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

// QFunctions for pointwise evaluation of derived field quantities at arbitrary points
// of 3D volume elements (for example the nodal points of an interpolatory output space
// for visualization). No quadrature weighting is applied; the output is the pointwise
// value. Contributions of repeated applications accumulate at the
// output L-vector (CeedOperatorApplyAdd adds the restriction transpose into the output
// vector), which implements the real/imaginary part sums of the quadratic quantities
// (energy densities, time-averaged Poynting vector). The QFunction itself writes (does
// not accumulate) its output.
//
// Inputs:
//  - in[0]: volume geometry data, shape [11, Q]: {attr, w * detJ (unused),
//           adj(J)^T / detJ (3x3, column-major)} for the Piola transformations and
//           attribute-based material property lookup.
//  - in[1], in[2]: field inputs (reference components at the evaluation points).
// The context is a CeedIntScalar array: [0].second = scale, followed by a material
// property coefficient context.

// Electric energy density: v += 0.5 * scale * (eps_r E) . E with E = adj(J)^T/detJ u.
CEED_QFUNCTION(f_eval_energy_e_33)(void *__restrict__ ctx_, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *adjJt = in[0] + 2 * Q, *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar adjJt_loc[9], E[3], eps[9], D[3];
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    MultAx33(adjJt_loc, u_loc, E);
    CoeffUnpack3(ctx + 1, (CeedInt)attr[i], eps);
    MultAx33(eps, E, D);
    v[i] = 0.5 * ctx[0].second * (D[0] * E[0] + D[1] * E[1] + D[2] * E[2]);
  }
  return 0;
}

// Magnetic energy density: v += 0.5 * scale * (mu^-1 B) . B with B = J/detJ u.
CEED_QFUNCTION(f_eval_energy_m_33)(void *__restrict__ ctx_, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *adjJt = in[0] + 2 * Q, *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar adjJt_loc[9], J_loc[9], B[3], invmu[9], H[3];
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    AdjJt33(adjJt_loc, J_loc);
    MultAx33(J_loc, u_loc, B);
    CoeffUnpack3(ctx + 1, (CeedInt)attr[i], invmu);
    MultAx33(invmu, B, H);
    v[i] = 0.5 * ctx[0].second * (H[0] * B[0] + H[1] * B[1] + H[2] * B[2]);
  }
  return 0;
}

// Poynting vector: v += scale * E x (mu^-1 B) with E = adj(J)^T/detJ u_e (H(curl)) and
// B = J/detJ u_b (H(div)).
CEED_QFUNCTION(f_eval_poynting_33)(void *__restrict__ ctx_, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *adjJt = in[0] + 2 * Q, *u_e = in[1], *u_b = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar ue_loc[3] = {u_e[i + Q * 0], u_e[i + Q * 1], u_e[i + Q * 2]};
    const CeedScalar ub_loc[3] = {u_b[i + Q * 0], u_b[i + Q * 1], u_b[i + Q * 2]};
    CeedScalar adjJt_loc[9], J_loc[9], E[3], B[3], invmu[9], H[3];
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    MultAx33(adjJt_loc, ue_loc, E);
    AdjJt33(adjJt_loc, J_loc);
    MultAx33(J_loc, ub_loc, B);
    CoeffUnpack3(ctx + 1, (CeedInt)attr[i], invmu);
    MultAx33(invmu, B, H);
    const CeedScalar s = ctx[0].second;
    v[i + Q * 0] = s * (E[1] * H[2] - E[2] * H[1]);
    v[i + Q * 1] = s * (E[2] * H[0] - E[0] * H[2]);
    v[i + Q * 2] = s * (E[0] * H[1] - E[1] * H[0]);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_EVAL_33_QF_H
