// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_2_QF_H
#define PALACE_LIBCEED_L2_2_QF_H

#include "../coeff/coeff_1_qf.h"

CEED_QFUNCTION(f_apply_l2_2)(void *__restrict__ ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *qw = in[1], *u = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);
    const CeedScalar w = qw[i] * qw[i] / wdetJ[i];

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = w * coeff * (u0 + u1);
    v[i + Q * 1] = w * coeff * (u0 + u1);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_2_QF_H
