// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_H1_BUILD_3_QF_H
#define PALACE_LIBCEED_H1_BUILD_3_QF_H

#include "../coeff/coeff_1_qf.h"

CEED_QFUNCTION(f_build_h1_3)(void *__restrict__ ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q;
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

    qd[i + Q * 0] = wdetJ[i] * coeff;
    qd[i + Q * 1] = wdetJ[i] * coeff;
    qd[i + Q * 2] = wdetJ[i] * coeff;
    qd[i + Q * 3] = wdetJ[i] * coeff;
    qd[i + Q * 4] = wdetJ[i] * coeff;
    qd[i + Q * 5] = wdetJ[i] * coeff;
    qd[i + Q * 6] = wdetJ[i] * coeff;
    qd[i + Q * 7] = wdetJ[i] * coeff;
    qd[i + Q * 8] = wdetJ[i] * coeff;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_H1_BUILD_3_QF_H
