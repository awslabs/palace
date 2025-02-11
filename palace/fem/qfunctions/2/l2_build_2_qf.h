// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_BUILD_2_QF_H
#define PALACE_LIBCEED_L2_BUILD_2_QF_H

#include "../coeff/coeff_2_qf.h"

CEED_QFUNCTION(f_build_l2_2)(void *__restrict__ ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *qw = in[1];
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[4];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    const CeedScalar w = qw[i] * qw[i] / wdetJ[i];

    qd[i + Q * 0] = w * coeff[0];
    qd[i + Q * 1] = w * coeff[1];
    qd[i + Q * 2] = w * coeff[2];
    qd[i + Q * 3] = w * coeff[3];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_BUILD_2_QF_H
