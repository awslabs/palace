// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_UTILS_32_QF_H
#define PALACE_LIBCEED_UTILS_32_QF_H

#include <math.h>

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

CEED_QFUNCTION_HELPER void MultAtBCx32(const CeedScalar A[6], const CeedScalar B[9],
                                       const CeedScalar C[6], const CeedScalar x[2],
                                       CeedScalar y[2])
{
  // A, C: 0 3   B: 0 3 6
  //       1 4      1 4 7
  //       2 5      2 5 8
  CeedScalar z[3], t;

  y[0] = C[0] * x[0] + C[3] * x[1];
  y[1] = C[1] * x[0] + C[4] * x[1];
  t = C[2] * x[0] + C[5] * x[1];

  z[0] = B[0] * y[0] + B[3] * y[1] + B[6] * t;
  z[1] = B[1] * y[0] + B[4] * y[1] + B[7] * t;
  z[2] = B[2] * y[0] + B[5] * y[1] + B[8] * t;

  y[0] = A[0] * z[0] + A[1] * z[1] + A[2] * z[2];
  y[1] = A[3] * z[0] + A[4] * z[1] + A[5] * z[2];
}

CEED_QFUNCTION_HELPER void MultBAx32(const CeedScalar A[6], const CeedScalar B[9],
                                     const CeedScalar x[2], CeedScalar y[3])
{
  // A: 0 3   B: 0 3 6
  //    1 4      1 4 7
  //    2 5      2 5 8
  CeedScalar z[3];

  z[0] = A[0] * x[0] + A[3] * x[1];
  z[1] = A[1] * x[0] + A[4] * x[1];
  z[2] = A[2] * x[0] + A[5] * x[1];

  y[0] = B[0] * z[0] + B[3] * z[1] + B[6] * z[2];
  y[1] = B[1] * z[0] + B[4] * z[1] + B[7] * z[2];
  y[2] = B[2] * z[0] + B[5] * z[1] + B[8] * z[2];
}

CEED_QFUNCTION_HELPER void MultAtBA32(const CeedScalar A[6], const CeedScalar B[9],
                                      CeedScalar C[4])
{
  // A: 0 3   B: 0 3 6   C: 0 2
  //    1 4      1 4 7      1 3
  //    2 5      2 5 8

  // First compute entries of R = B A.
  const CeedScalar R11 = B[0] * A[0] + B[3] * A[1] + B[6] * A[2];
  const CeedScalar R21 = B[1] * A[0] + B[4] * A[1] + B[7] * A[2];
  const CeedScalar R31 = B[2] * A[0] + B[5] * A[1] + B[8] * A[2];
  const CeedScalar R12 = B[0] * A[3] + B[3] * A[4] + B[6] * A[5];
  const CeedScalar R22 = B[1] * A[3] + B[4] * A[4] + B[7] * A[5];
  const CeedScalar R32 = B[2] * A[3] + B[5] * A[4] + B[8] * A[5];

  C[0] = A[0] * R11 + A[1] * R21 + A[2] * R31;
  C[1] = A[3] * R11 + A[4] * R21 + A[5] * R31;
  C[2] = A[0] * R12 + A[1] * R22 + A[2] * R32;
  C[3] = A[3] * R12 + A[4] * R22 + A[5] * R32;
}

CEED_QFUNCTION_HELPER void MultAtBC32(const CeedScalar A[6], const CeedScalar B[9],
                                      const CeedScalar C[6], CeedScalar D[4])
{
  // A, C: 0 3   B: 0 3 6   D: 0 2
  //       1 4      1 4 7      1 3
  //       2 5      2 5 8

  // First compute entries of R = B C.
  const CeedScalar R11 = B[0] * C[0] + B[3] * C[1] + B[6] * C[2];
  const CeedScalar R21 = B[1] * C[0] + B[4] * C[1] + B[7] * C[2];
  const CeedScalar R31 = B[2] * C[0] + B[5] * C[1] + B[8] * C[2];
  const CeedScalar R12 = B[0] * C[3] + B[3] * C[4] + B[6] * C[5];
  const CeedScalar R22 = B[1] * C[3] + B[4] * C[4] + B[7] * C[5];
  const CeedScalar R32 = B[2] * C[3] + B[5] * C[4] + B[8] * C[5];

  D[0] = A[0] * R11 + A[1] * R21 + A[2] * R31;
  D[1] = A[3] * R11 + A[4] * R21 + A[5] * R31;
  D[2] = A[0] * R12 + A[1] * R22 + A[2] * R32;
  D[3] = A[3] * R12 + A[4] * R22 + A[5] * R32;
}

CEED_QFUNCTION_HELPER void MultBA32(const CeedScalar A[6], const CeedScalar B[9],
                                    CeedScalar C[6])
{
  // A, C: 0 3   B: 0 3 6
  //       1 4      1 4 7
  //       2 5      2 5 8
  C[0] = B[0] * A[0] + B[3] * A[1] + B[6] * A[2];
  C[1] = B[1] * A[0] + B[4] * A[1] + B[7] * A[2];
  C[2] = B[2] * A[0] + B[5] * A[1] + B[8] * A[2];
  C[3] = B[0] * A[3] + B[3] * A[4] + B[6] * A[5];
  C[4] = B[1] * A[3] + B[4] * A[4] + B[7] * A[5];
  C[5] = B[2] * A[3] + B[5] * A[4] + B[8] * A[5];
}

#endif  // PALACE_LIBCEED_UTILS_32_QF_H
