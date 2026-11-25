// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SHAPE_22_QF_H
#define PALACE_LIBCEED_SHAPE_22_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_2_qf.h"
#include "utils_22_qf.h"

// Derivative of (Q u, w) for covariant fields (H(curl) values or gradients of H1 fields):
// S = (w ⋅ Qu) I - w (Qu)ᵀ - u (Qᵀw)ᵀ.
CEED_QFUNCTION(f_shape_hcurl_22)(void *__restrict__ ctx, CeedInt Q,
                                 const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *w = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    const CeedScalar w_loc[2] = {w[i + Q * 0], w[i + Q * 1]};
    CeedScalar coeff[4], adjJt_loc[4], pu[2], pw[2], qpu[2], qtpw[2];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    MultBx22(adjJt_loc, u_loc, pu);
    MultBx22(adjJt_loc, w_loc, pw);
    MultBx22(coeff, pu, qpu);
    qtpw[0] = coeff[0] * pw[0] + coeff[1] * pw[1];
    qtpw[1] = coeff[2] * pw[0] + coeff[3] * pw[1];
    const CeedScalar s = pw[0] * qpu[0] + pw[1] * qpu[1];

    for (CeedInt c = 0; c < 2; c++)
    {
      const CeedScalar *col = adjJt_loc + 2 * c;
      const CeedScalar a = col[0] * qpu[0] + col[1] * qpu[1];
      const CeedScalar b = col[0] * qtpw[0] + col[1] * qtpw[1];
      for (CeedInt r = 0; r < 2; r++)
      {
        v[i + Q * (r + 2 * c)] = wdetJ[i] * (s * col[r] - pw[r] * a - pu[r] * b);
      }
    }
  }
  return 0;
}

// Derivative of (Q u, w) for scalar fields mapped by 1 / |J| (curls of H(curl) fields):
// S = -(w Q u) I. in[1] is the quadrature weight, and u and w are in[2] and in[3].
CEED_QFUNCTION(f_shape_l2_22)(void *__restrict__ ctx, CeedInt Q,
                              const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1],
                   *u = in[2], *w = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar s = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]) *
                         qw[i] * qw[i] / wdetJ[i] * u[i] * w[i];
    for (CeedInt k = 0; k < 4; k++)
    {
      v[i + Q * k] = -s * adjJt[i + Q * k];
    }
  }
  return 0;
}

// Derivative of (Q u, w): S = Q u w I.
CEED_QFUNCTION(f_shape_mass_22)(void *__restrict__ ctx, CeedInt Q,
                                const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *w = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar s =
        CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]) * wdetJ[i] * u[i] * w[i];
    for (CeedInt k = 0; k < 4; k++)
    {
      v[i + Q * k] = s * adjJt[i + Q * k];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_SHAPE_22_QF_H
