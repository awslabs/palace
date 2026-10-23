// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2H1_ERROR_QF_H
#define PALACE_LIBCEED_L2H1_ERROR_QF_H

#include "coeff/coeff_1_qf.h"

// Scalar error estimator qfunction for 2D curl flux estimation.
// Computes ||c1 * u1 - c2 * u2||^2 where u1 and u2 are scalar fields
// (L2 with u1 = 1 / det(J) ̂u1, and H1) and c1, c2 are scalar coefficients.
// Used for the B-field error in 2D where B = curl E is a scalar.
CEED_QFUNCTION(f_apply_l2h1_error)(void *__restrict__ ctx, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *qw = in[1], *u1 = in[2], *u2 = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff1, coeff2;
    CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i], &coeff1);
    CoeffUnpack1(CoeffPairSecond<1>((const CeedIntScalar *)ctx), (CeedInt)attr[i], &coeff2);
    const CeedScalar diff = coeff1 * qw[i] * u1[i] / wdetJ[i] - coeff2 * u2[i];
    v[i] = wdetJ[i] * diff * diff;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2H1_ERROR_QF_H
