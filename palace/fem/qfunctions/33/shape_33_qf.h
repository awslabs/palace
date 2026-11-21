// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SHAPE_33_QF_H
#define PALACE_LIBCEED_SHAPE_33_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

// Derivative of (Q u, w) for covariant fields (H(curl) values or gradients of H1 fields):
// S = (w ⋅ Qu) I - w (Qu)ᵀ - u (Qᵀw)ᵀ.
CEED_QFUNCTION(f_shape_hcurl_33)(void *__restrict__ ctx, CeedInt Q,
                                 const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *w = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    const CeedScalar w_loc[3] = {w[i + Q * 0], w[i + Q * 1], w[i + Q * 2]};
    CeedScalar coeff[9], adjJt_loc[9], pu[3], pw[3], qpu[3], qtpw[3];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    MultAx33(adjJt_loc, u_loc, pu);
    MultAx33(adjJt_loc, w_loc, pw);
    MultAx33(coeff, pu, qpu);
    qtpw[0] = coeff[0] * pw[0] + coeff[1] * pw[1] + coeff[2] * pw[2];
    qtpw[1] = coeff[3] * pw[0] + coeff[4] * pw[1] + coeff[5] * pw[2];
    qtpw[2] = coeff[6] * pw[0] + coeff[7] * pw[1] + coeff[8] * pw[2];
    const CeedScalar s = pw[0] * qpu[0] + pw[1] * qpu[1] + pw[2] * qpu[2];

    for (CeedInt c = 0; c < 3; c++)
    {
      const CeedScalar *col = adjJt_loc + 3 * c;
      const CeedScalar a = col[0] * qpu[0] + col[1] * qpu[1] + col[2] * qpu[2];
      const CeedScalar b = col[0] * qtpw[0] + col[1] * qtpw[1] + col[2] * qtpw[2];
      for (CeedInt r = 0; r < 3; r++)
      {
        v[i + Q * (r + 3 * c)] = wdetJ[i] * (s * col[r] - pw[r] * a - pu[r] * b);
      }
    }
  }
  return 0;
}

// Derivative of (Q u, w) for contravariant fields (H(div) values or curls of H(curl)
// fields): S = (Qu) wᵀ + (Qᵀw) uᵀ - (w ⋅ Qu) I.
CEED_QFUNCTION(f_shape_hdiv_33)(void *__restrict__ ctx, CeedInt Q,
                                const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *w = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    const CeedScalar w_loc[3] = {w[i + Q * 0], w[i + Q * 1], w[i + Q * 2]};
    CeedScalar coeff[9], adjJt_loc[9], J_loc[9], pu[3], pw[3], qpu[3], qtpw[3];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    const CeedScalar inv_detJ = AdjJt33<true>(adjJt_loc, J_loc);  // J_loc = J / det(J)
    MultAx33(J_loc, u_loc, pu);
    MultAx33(J_loc, w_loc, pw);
    MultAx33(coeff, pu, qpu);
    qtpw[0] = coeff[0] * pw[0] + coeff[1] * pw[1] + coeff[2] * pw[2];
    qtpw[1] = coeff[3] * pw[0] + coeff[4] * pw[1] + coeff[5] * pw[2];
    qtpw[2] = coeff[6] * pw[0] + coeff[7] * pw[1] + coeff[8] * pw[2];
    const CeedScalar s = pw[0] * qpu[0] + pw[1] * qpu[1] + pw[2] * qpu[2];

    for (CeedInt c = 0; c < 3; c++)
    {
      for (CeedInt r = 0; r < 3; r++)
      {
        v[i + Q * (r + 3 * c)] =
            wdetJ[i] * (inv_detJ * (qpu[r] * w_loc[c] + qtpw[r] * u_loc[c]) -
                        s * adjJt_loc[r + 3 * c]);
      }
    }
  }
  return 0;
}

// Derivative of (Q u, w): S = Q u w I.
CEED_QFUNCTION(f_shape_mass_33)(void *__restrict__ ctx, CeedInt Q,
                                const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *w = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar s =
        CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]) * wdetJ[i] * u[i] * w[i];
    for (CeedInt k = 0; k < 9; k++)
    {
      v[i + Q * k] = s * adjJt[i + Q * k];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_SHAPE_33_QF_H
