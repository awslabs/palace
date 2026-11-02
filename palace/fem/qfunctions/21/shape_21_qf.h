// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SHAPE_21_QF_H
#define PALACE_LIBCEED_SHAPE_21_QF_H

#include "../coeff/coeff_1_qf.h"

// Derivative of (Q u, w) on a curve: S = Q u w P, with P the tangential projection, so
// the output is Q u w times the weighted |J| J (JᵀJ)⁻¹.
CEED_QFUNCTION(f_shape_mass_21)(void *__restrict__ ctx, CeedInt Q,
                                const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *w = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar s =
        CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]) * wdetJ[i] * u[i] * w[i];
    v[i + Q * 0] = s * adjJt[i + Q * 0];
    v[i + Q * 1] = s * adjJt[i + Q * 1];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_SHAPE_21_QF_H
