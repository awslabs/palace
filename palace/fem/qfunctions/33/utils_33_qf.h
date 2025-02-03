// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_UTILS_33_QF_H
#define PALACE_LIBCEED_UTILS_33_QF_H

#include <math.h>

CEED_QFUNCTION_HELPER CeedScalar DetJ33(const CeedScalar J[9])
{
  // J: 0 3 6
  //    1 4 7
  //    2 5 8
  return J[0] * (J[4] * J[8] - J[5] * J[7]) - J[1] * (J[3] * J[8] - J[5] * J[6]) +
         J[2] * (J[3] * J[7] - J[4] * J[6]);
}

template <bool ComputeDet = false>
CEED_QFUNCTION_HELPER CeedScalar AdjJt33(const CeedScalar J[9], CeedScalar adjJt[9])
{
  // Compute adj(J)^T / det(J) and store the result.
  // J: 0 3 6
  //    1 4 7
  //    2 5 8
  adjJt[0] = J[4] * J[8] - J[7] * J[5];
  adjJt[3] = J[7] * J[2] - J[1] * J[8];
  adjJt[6] = J[1] * J[5] - J[4] * J[2];
  adjJt[1] = J[6] * J[5] - J[3] * J[8];
  adjJt[4] = J[0] * J[8] - J[6] * J[2];
  adjJt[7] = J[3] * J[2] - J[0] * J[5];
  adjJt[2] = J[3] * J[7] - J[6] * J[4];
  adjJt[5] = J[6] * J[1] - J[0] * J[7];
  adjJt[8] = J[0] * J[4] - J[3] * J[1];
  return ComputeDet ? (J[0] * adjJt[0] + J[1] * adjJt[1] + J[2] * adjJt[2]) : 0.0;
}

CEED_QFUNCTION_HELPER void MatUnpack33(const CeedScalar *A, const CeedInt A_stride,
                                       CeedScalar A_loc[9])
{
  A_loc[0] = A[A_stride * 0];
  A_loc[1] = A[A_stride * 1];
  A_loc[2] = A[A_stride * 2];
  A_loc[3] = A[A_stride * 3];
  A_loc[4] = A[A_stride * 4];
  A_loc[5] = A[A_stride * 5];
  A_loc[6] = A[A_stride * 6];
  A_loc[7] = A[A_stride * 7];
  A_loc[8] = A[A_stride * 8];
}

CEED_QFUNCTION_HELPER void MultBx33(const CeedScalar B[9], const CeedScalar x[3],
                                    CeedScalar y[3])
{
  // B: 0 3 6
  //    1 4 7
  //    2 5 8
  y[0] = B[0] * x[0] + B[3] * x[1] + B[6] * x[2];
  y[1] = B[1] * x[0] + B[4] * x[1] + B[7] * x[2];
  y[2] = B[2] * x[0] + B[5] * x[1] + B[8] * x[2];
}

CEED_QFUNCTION_HELPER void MultAtBCx33(const CeedScalar A[9], const CeedScalar B[9],
                                       const CeedScalar C[9], const CeedScalar x[3],
                                       CeedScalar y[3])
{
  // A, B, C: 0 3 6
  //          1 4 7
  //          2 5 8
  CeedScalar z[3];

  y[0] = C[0] * x[0] + C[3] * x[1] + C[6] * x[2];
  y[1] = C[1] * x[0] + C[4] * x[1] + C[7] * x[2];
  y[2] = C[2] * x[0] + C[5] * x[1] + C[8] * x[2];

  z[0] = B[0] * y[0] + B[3] * y[1] + B[6] * y[2];
  z[1] = B[1] * y[0] + B[4] * y[1] + B[7] * y[2];
  z[2] = B[2] * y[0] + B[5] * y[1] + B[8] * y[2];

  y[0] = A[0] * z[0] + A[1] * z[1] + A[2] * z[2];
  y[1] = A[3] * z[0] + A[4] * z[1] + A[5] * z[2];
  y[2] = A[6] * z[0] + A[7] * z[1] + A[8] * z[2];
}

CEED_QFUNCTION_HELPER void MultBAx33(const CeedScalar A[9], const CeedScalar B[9],
                                     const CeedScalar x[3], CeedScalar y[3])
{
  // A, B: 0 3 6
  //       1 4 7
  //       2 5 8
  CeedScalar z[3];

  z[0] = A[0] * x[0] + A[3] * x[1] + A[6] * x[2];
  z[1] = A[1] * x[0] + A[4] * x[1] + A[7] * x[2];
  z[2] = A[2] * x[0] + A[5] * x[1] + A[8] * x[2];

  y[0] = B[0] * z[0] + B[3] * z[1] + B[6] * z[2];
  y[1] = B[1] * z[0] + B[4] * z[1] + B[7] * z[2];
  y[2] = B[2] * z[0] + B[5] * z[1] + B[8] * z[2];
}

CEED_QFUNCTION_HELPER void MultAtBA33(const CeedScalar A[9], const CeedScalar B[9],
                                      CeedScalar C[9])
{
  // A, B, C: 0 3 6
  //          1 4 7
  //          2 5 8

  // First compute entries of R = B A.
  const CeedScalar R11 = B[0] * A[0] + B[3] * A[1] + B[6] * A[2];
  const CeedScalar R21 = B[1] * A[0] + B[4] * A[1] + B[7] * A[2];
  const CeedScalar R31 = B[2] * A[0] + B[5] * A[1] + B[8] * A[2];
  const CeedScalar R12 = B[0] * A[3] + B[3] * A[4] + B[6] * A[5];
  const CeedScalar R22 = B[1] * A[3] + B[4] * A[4] + B[7] * A[5];
  const CeedScalar R32 = B[2] * A[3] + B[5] * A[4] + B[8] * A[5];
  const CeedScalar R13 = B[0] * A[6] + B[3] * A[7] + B[6] * A[8];
  const CeedScalar R23 = B[1] * A[6] + B[4] * A[7] + B[7] * A[8];
  const CeedScalar R33 = B[2] * A[6] + B[5] * A[7] + B[8] * A[8];

  C[0] = A[0] * R11 + A[1] * R21 + A[2] * R31;
  C[1] = A[3] * R11 + A[4] * R21 + A[5] * R31;
  C[2] = A[6] * R11 + A[7] * R21 + A[8] * R31;
  C[3] = A[0] * R12 + A[1] * R22 + A[2] * R32;
  C[4] = A[3] * R12 + A[4] * R22 + A[5] * R32;
  C[5] = A[6] * R12 + A[7] * R22 + A[8] * R32;
  C[6] = A[0] * R13 + A[1] * R23 + A[2] * R33;
  C[7] = A[3] * R13 + A[4] * R23 + A[5] * R33;
  C[8] = A[6] * R13 + A[7] * R23 + A[8] * R33;
}

CEED_QFUNCTION_HELPER void MultAtBC33(const CeedScalar A[9], const CeedScalar B[9],
                                      const CeedScalar C[9], CeedScalar D[9])
{
  // A, B, C, D: 0 3 6
  //             1 4 7
  //             2 5 8

  // First compute entries of R = B C.
  const CeedScalar R11 = B[0] * C[0] + B[3] * C[1] + B[6] * C[2];
  const CeedScalar R21 = B[1] * C[0] + B[4] * C[1] + B[7] * C[2];
  const CeedScalar R31 = B[2] * C[0] + B[5] * C[1] + B[8] * C[2];
  const CeedScalar R12 = B[0] * C[3] + B[3] * C[4] + B[6] * C[5];
  const CeedScalar R22 = B[1] * C[3] + B[4] * C[4] + B[7] * C[5];
  const CeedScalar R32 = B[2] * C[3] + B[5] * C[4] + B[8] * C[5];
  const CeedScalar R13 = B[0] * C[6] + B[3] * C[7] + B[6] * C[8];
  const CeedScalar R23 = B[1] * C[6] + B[4] * C[7] + B[7] * C[8];
  const CeedScalar R33 = B[2] * C[6] + B[5] * C[7] + B[8] * C[8];

  D[0] = A[0] * R11 + A[1] * R21 + A[2] * R31;
  D[1] = A[3] * R11 + A[4] * R21 + A[5] * R31;
  D[2] = A[6] * R11 + A[7] * R21 + A[8] * R31;
  D[3] = A[0] * R12 + A[1] * R22 + A[2] * R32;
  D[4] = A[3] * R12 + A[4] * R22 + A[5] * R32;
  D[5] = A[6] * R12 + A[7] * R22 + A[8] * R32;
  D[6] = A[0] * R13 + A[1] * R23 + A[2] * R33;
  D[7] = A[3] * R13 + A[4] * R23 + A[5] * R33;
  D[8] = A[6] * R13 + A[7] * R23 + A[8] * R33;
}

CEED_QFUNCTION_HELPER void MultBA33(const CeedScalar A[9], const CeedScalar B[9],
                                    CeedScalar C[9])
{
  // A, B, C: 0 3 6
  //          1 4 7
  //          2 5 8
  C[0] = B[0] * A[0] + B[3] * A[1] + B[6] * A[2];
  C[1] = B[1] * A[0] + B[4] * A[1] + B[7] * A[2];
  C[2] = B[2] * A[0] + B[5] * A[1] + B[8] * A[2];
  C[3] = B[0] * A[3] + B[3] * A[4] + B[6] * A[5];
  C[4] = B[1] * A[3] + B[4] * A[4] + B[7] * A[5];
  C[5] = B[2] * A[3] + B[5] * A[4] + B[8] * A[5];
  C[6] = B[0] * A[6] + B[3] * A[7] + B[6] * A[8];
  C[7] = B[1] * A[6] + B[4] * A[7] + B[7] * A[8];
  C[8] = B[2] * A[6] + B[5] * A[7] + B[8] * A[8];
}

#endif  // PALACE_LIBCEED_UTILS_33_QF_H
