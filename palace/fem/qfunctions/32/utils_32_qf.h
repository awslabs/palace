// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_UTILS_32_QF_H
#define PALACE_LIBCEED_UTILS_32_QF_H

#ifndef CEED_RUNNING_JIT_PASS
#include <math.h>
#endif

CEED_QFUNCTION_HELPER CeedScalar DetJ32(const CeedScalar J[6])
{
  // J: 0 3
  //    1 4
  //    2 5
  const CeedScalar E = J[0] * J[0] + J[1] * J[1] + J[2] * J[2];
  const CeedScalar G = J[3] * J[3] + J[4] * J[4] + J[5] * J[5];
  const CeedScalar F = J[0] * J[3] + J[1] * J[4] + J[2] * J[5];
  return sqrt(E * G - F * F);
}

template <bool ComputeDet = false>
CEED_QFUNCTION_HELPER CeedScalar AdjJt32(const CeedScalar J[6], CeedScalar adjJt[6])
{
  // Compute adj(J)^T / det(J) and store the result.
  // J: 0 3
  //    1 4
  //    2 5
  const CeedScalar E = J[0] * J[0] + J[1] * J[1] + J[2] * J[2];
  const CeedScalar G = J[3] * J[3] + J[4] * J[4] + J[5] * J[5];
  const CeedScalar F = J[0] * J[3] + J[1] * J[4] + J[2] * J[5];
  const CeedScalar d = sqrt(E * G - F * F);
  adjJt[0] = (G * J[0] - F * J[3]) / d;
  adjJt[1] = (G * J[1] - F * J[4]) / d;
  adjJt[2] = (G * J[2] - F * J[5]) / d;
  adjJt[3] = (E * J[3] - F * J[0]) / d;
  adjJt[4] = (E * J[4] - F * J[1]) / d;
  adjJt[5] = (E * J[5] - F * J[2]) / d;
  return ComputeDet ? d : 0.0;
}

CEED_QFUNCTION_HELPER void MatUnpack32(const CeedScalar *A, const CeedInt A_stride,
                                       CeedScalar A_loc[6])
{
  A_loc[0] = A[A_stride * 0];
  A_loc[1] = A[A_stride * 1];
  A_loc[2] = A[A_stride * 2];
  A_loc[3] = A[A_stride * 3];
  A_loc[4] = A[A_stride * 4];
  A_loc[5] = A[A_stride * 5];
}

CEED_QFUNCTION_HELPER void MultAtBCx32(const CeedScalar A[6], const CeedScalar alpha,
                                       const CeedScalar C[6], const CeedScalar x[2],
                                       CeedScalar y[2])
{
  // A, C: 0 3
  //       1 4
  //       2 5
  CeedScalar R[4];  // R = At * C
  R[0] = A[0] * C[0] + A[1] * C[2] + A[2] * C[4];
  R[1] = A[3] * C[0] + A[4] * C[2] + A[5] * C[4];
  R[2] = A[0] * C[1] + A[1] * C[3] + A[2] * C[5];
  R[3] = A[3] * C[1] + A[4] * C[3] + A[5] * C[5];

  // y = R * x
  y[0] = alpha * (R[0] * x[0] + R[2] * x[1]);
  y[1] = alpha * (R[1] * x[0] + R[3] * x[1]);
}

CEED_QFUNCTION_HELPER void MultBAx32(const CeedScalar A[6], const CeedScalar alpha,
                                     const CeedScalar x[2], CeedScalar y[3])
{
  // A: 0 3
  //    1 4
  //    2 5

  y[0] = alpha * (A[0] * x[0] + A[3] * x[1]);
  y[1] = alpha * (A[1] * x[0] + A[4] * x[1]);
  y[2] = alpha * (A[2] * x[0] + A[5] * x[1]);
}

CEED_QFUNCTION_HELPER void MultAtBA32(const CeedScalar A[6], const CeedScalar alpha,
                                      CeedScalar C[4])
{
  // A: 0 3   C: 0 2
  //    1 4      1 3
  //    2 5

  // This is a symmetric matrix so maybe we can optimize it?
  C[0] = alpha * (A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
  C[1] = alpha * (A[3] * A[0] + A[4] * A[1] + A[5] * A[2]);
  C[2] = alpha * (A[0] * A[3] + A[1] * A[4] + A[2] * A[5]);
  C[3] = alpha * (A[3] * A[3] + A[4] * A[4] + A[5] * A[5]);
}

CEED_QFUNCTION_HELPER void MultAtBC32(const CeedScalar A[6], const CeedScalar alpha,
                                      const CeedScalar C[6], CeedScalar D[4])
{
  // A, C: 0 3   D: 0 2
  //       1 4      1 3
  //       2 5

  D[0] = alpha * (A[0] * C[0] + A[1] * C[1] + A[2] * C[2]);
  D[1] = alpha * (A[3] * C[0] + A[4] * C[1] + A[5] * C[2]);
  D[2] = alpha * (A[0] * C[3] + A[1] * C[4] + A[2] * C[5]);
  D[3] = alpha * (A[3] * C[3] + A[4] * C[4] + A[5] * C[5]);
}

CEED_QFUNCTION_HELPER void MultBA32(const CeedScalar A[6], const CeedScalar alpha,
                                    CeedScalar C[6])
{
  // A, C: 0 3   B: 0 3 6
  //       1 4      1 4 7
  //       2 5      2 5 8
  C[0] = alpha * A[0];
  C[1] = alpha * A[1];
  C[2] = alpha * A[2];
  C[3] = alpha * A[3];
  C[4] = alpha * A[4];
  C[5] = alpha * A[5];
}

#endif  // PALACE_LIBCEED_UTILS_32_QF_H
