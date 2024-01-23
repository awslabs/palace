// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_H1_BUILD_1_QF_H
#define PALACE_LIBCEED_H1_BUILD_1_QF_H

#include "../coeff/coeff_1_qf.h"

CEED_QFUNCTION(f_build_h1_1)(void *__restrict__ ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q;
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

    qd[i] = coeff * wdetJ[i];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_H1_BUILD_1_QF_H
