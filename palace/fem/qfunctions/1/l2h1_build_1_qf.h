// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2H1_BUILD_1_QF_H
#define PALACE_LIBCEED_L2H1_BUILD_1_QF_H

#include "../coeff/coeff_1_qf.h"

// Quadrature data of f_apply_l2h1_1.
CEED_QFUNCTION(f_build_l2h1_1)(void *__restrict__ ctx, CeedInt Q,
                               const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *qw = in[1];
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

    qd[i] = coeff * qw[i];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2H1_BUILD_1_QF_H
