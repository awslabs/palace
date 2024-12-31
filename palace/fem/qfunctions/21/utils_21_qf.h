// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_UTILS_21_QF_H
#define PALACE_LIBCEED_UTILS_21_QF_H

#include <math.h>

CEED_QFUNCTION_HELPER CeedScalar DetJ21(const CeedScalar J[2])
{
  // J: 0
  //    1
  return sqrt(J[0] * J[0] + J[1] * J[1]);
}

template <bool ComputeDet = false>
CEED_QFUNCTION_HELPER CeedScalar AdjJt21(const CeedScalar J[2], CeedScalar adjJt[2])
{
  // Compute adj(J)^T / det(J) and store the result.
  // J: 0   adj(J): 1/sqrt(J^T J) J^T
  //    1
  const CeedScalar d = sqrt(J[0] * J[0] + J[1] * J[1]);
  adjJt[0] = J[0] / d;
  adjJt[1] = J[1] / d;
  return ComputeDet ? d : 0.0;
}

CEED_QFUNCTION_HELPER void MatUnpack21(const CeedScalar *A, const CeedInt A_stride,
                                       CeedScalar A_loc[2])
{
  A_loc[0] = A[A_stride * 0];
  A_loc[1] = A[A_stride * 1];
}

CEED_QFUNCTION_HELPER void MultAtBCx21(const CeedScalar A[2], const CeedScalar B[4],
                                       const CeedScalar C[2], const CeedScalar x[1],
                                       CeedScalar y[1])
{
  // A: 0   B: 0 2   C: 0
  //    1      1 3      1
  CeedScalar z[2], t;

  y[0] = C[0] * x[0];
  t = C[1] * x[0];

  z[0] = B[0] * y[0] + B[2] * t;
  z[1] = B[1] * y[0] + B[3] * t;

  y[0] = A[0] * z[0] + A[1] * z[1];
}

CEED_QFUNCTION_HELPER void MultBAx21(const CeedScalar A[2], const CeedScalar B[4],
                                     const CeedScalar x[1], CeedScalar y[2])
{
  // A: 0   B: 0 2
  //    1      1 3
  CeedScalar z[2];

  z[0] = A[0] * x[0];
  z[1] = A[1] * x[0];

  y[0] = B[0] * z[0] + B[2] * z[1];
  y[1] = B[1] * z[0] + B[3] * z[1];
}

CEED_QFUNCTION_HELPER void MultAtBA21(const CeedScalar A[2], const CeedScalar B[4],
                                      CeedScalar C[1])
{
  // A: 0   B: 0 2   C: 0
  //    1      1 3

  // First compute entries of R = B A.
  const CeedScalar R11 = B[0] * A[0] + B[2] * A[1];
  const CeedScalar R21 = B[1] * A[0] + B[3] * A[1];

  C[0] = A[0] * R11 + A[1] * R21;
}

CEED_QFUNCTION_HELPER void MultAtBC21(const CeedScalar A[2], const CeedScalar B[4],
                                      const CeedScalar C[2], CeedScalar D[1])
{
  // A, C: 0   B: 0 2   D: 0
  //       1      1 3

  // First compute entries of R = B C.
  const CeedScalar R11 = B[0] * C[0] + B[2] * C[1];
  const CeedScalar R21 = B[1] * C[0] + B[3] * C[1];

  D[0] = A[0] * R11 + A[1] * R21;
}

CEED_QFUNCTION_HELPER void MultBA21(const CeedScalar A[2], const CeedScalar B[4],
                                    CeedScalar C[2])
{
  // A: 0   B: 0 2   C: 0
  //    1      1 3      1
  C[0] = B[0] * A[0] + B[2] * A[1];
  C[1] = B[1] * A[0] + B[3] * A[1];
}

#endif  // PALACE_LIBCEED_UTILS_21_QF_H
