// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SHAPE_32_QF_H
#define PALACE_LIBCEED_SHAPE_32_QF_H

#include "../coeff/coeff_1_qf.h"

// Derivative of (Q u, w) on a surface: S = Q u w P, with P the tangential projection, so
// the output is Q u w times the weighted |J| J (JᵀJ)⁻¹.
CEED_QFUNCTION(f_shape_mass_32)(void *__restrict__ ctx, CeedInt Q,
                                const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *w = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar s =
        CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]) * wdetJ[i] * u[i] * w[i];
    for (CeedInt k = 0; k < 6; k++)
    {
      v[i + Q * k] = s * adjJt[i + Q * k];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_SHAPE_32_QF_H
