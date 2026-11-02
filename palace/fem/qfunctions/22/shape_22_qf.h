// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SHAPE_22_QF_H
#define PALACE_LIBCEED_SHAPE_22_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_2_qf.h"
#include "utils_22_qf.h"

// Derivative of (Q grad u, grad w): S = (∇w ⋅ Q∇u) I - ∇w (Q∇u)ᵀ - ∇u (Qᵀ∇w)ᵀ.
CEED_QFUNCTION(f_shape_diffusion_22)(void *__restrict__ ctx, CeedInt Q,
                                     const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *gu = in[1],
                   *gw = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar gu_loc[2] = {gu[i + Q * 0], gu[i + Q * 1]};
    const CeedScalar gw_loc[2] = {gw[i + Q * 0], gw[i + Q * 1]};
    CeedScalar coeff[4], adjJt_loc[4], du[2], dw[2], qdu[2], qtdw[2];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    MultBx22(adjJt_loc, gu_loc, du);
    MultBx22(adjJt_loc, gw_loc, dw);
    MultBx22(coeff, du, qdu);
    qtdw[0] = coeff[0] * dw[0] + coeff[1] * dw[1];
    qtdw[1] = coeff[2] * dw[0] + coeff[3] * dw[1];
    const CeedScalar s = dw[0] * qdu[0] + dw[1] * qdu[1];

    for (CeedInt c = 0; c < 2; c++)
    {
      const CeedScalar *col = adjJt_loc + 2 * c;
      const CeedScalar a = col[0] * qdu[0] + col[1] * qdu[1];
      const CeedScalar b = col[0] * qtdw[0] + col[1] * qtdw[1];
      for (CeedInt r = 0; r < 2; r++)
      {
        v[i + Q * (r + 2 * c)] = wdetJ[i] * (s * col[r] - dw[r] * a - du[r] * b);
      }
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
