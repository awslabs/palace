// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SHAPE_33_QF_H
#define PALACE_LIBCEED_SHAPE_33_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

// Derivative of (Q grad u, grad w): S = (∇w ⋅ Q∇u) I - ∇w (Q∇u)ᵀ - ∇u (Qᵀ∇w)ᵀ.
CEED_QFUNCTION(f_shape_diffusion_33)(void *__restrict__ ctx, CeedInt Q,
                                     const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *gu = in[1],
                   *gw = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar gu_loc[3] = {gu[i + Q * 0], gu[i + Q * 1], gu[i + Q * 2]};
    const CeedScalar gw_loc[3] = {gw[i + Q * 0], gw[i + Q * 1], gw[i + Q * 2]};
    CeedScalar coeff[9], adjJt_loc[9], du[3], dw[3], qdu[3], qtdw[3];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    MultAx33(adjJt_loc, gu_loc, du);
    MultAx33(adjJt_loc, gw_loc, dw);
    MultAx33(coeff, du, qdu);
    qtdw[0] = coeff[0] * dw[0] + coeff[1] * dw[1] + coeff[2] * dw[2];
    qtdw[1] = coeff[3] * dw[0] + coeff[4] * dw[1] + coeff[5] * dw[2];
    qtdw[2] = coeff[6] * dw[0] + coeff[7] * dw[1] + coeff[8] * dw[2];
    const CeedScalar s = dw[0] * qdu[0] + dw[1] * qdu[1] + dw[2] * qdu[2];

    for (CeedInt c = 0; c < 3; c++)
    {
      const CeedScalar *col = adjJt_loc + 3 * c;
      const CeedScalar a = col[0] * qdu[0] + col[1] * qdu[1] + col[2] * qdu[2];
      const CeedScalar b = col[0] * qtdw[0] + col[1] * qtdw[1] + col[2] * qtdw[2];
      for (CeedInt r = 0; r < 3; r++)
      {
        v[i + Q * (r + 3 * c)] = wdetJ[i] * (s * col[r] - dw[r] * a - du[r] * b);
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
