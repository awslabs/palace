// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_UTILS_22_QF_H
#define PALACE_LIBCEED_UTILS_22_QF_H

#include <math.h>

CEED_QFUNCTION_HELPER CeedScalar DetJ22(const CeedScalar J[4])
{
  // J: 0 2
  //    1 3
  return J[0] * J[3] - J[1] * J[2];
}

template <bool ComputeDet = false>
CEED_QFUNCTION_HELPER CeedScalar AdjJt22(const CeedScalar J[4], CeedScalar adjJt[4])
{
  // Compute adj(J)^T / det(J) and store the result.
  // J: 0 2   adj(J):  J22 -J12
  //    1 3           -J21  J11
  adjJt[0] = J[3];
  adjJt[1] = -J[2];
  adjJt[2] = -J[1];
  adjJt[3] = J[0];
  return ComputeDet ? (J[0] * J[3] - J[1] * J[2]) : 0.0;
}

CEED_QFUNCTION_HELPER void MatUnpack22(const CeedScalar *A, const CeedInt A_stride,
                                       CeedScalar A_loc[4])
{
  A_loc[0] = A[A_stride * 0];
  A_loc[1] = A[A_stride * 1];
  A_loc[2] = A[A_stride * 2];
  A_loc[3] = A[A_stride * 3];
}

CEED_QFUNCTION_HELPER void MultBx22(const CeedScalar B[4], const CeedScalar x[2],
                                    CeedScalar y[2])
{
  // B: 0 2
  //    1 3
  y[0] = B[0] * x[0] + B[2] * x[1];
  y[1] = B[1] * x[0] + B[3] * x[1];
}

CEED_QFUNCTION_HELPER void MultAtBCx22(const CeedScalar A[4], const CeedScalar B[4],
                                       const CeedScalar C[4], const CeedScalar x[2],
                                       CeedScalar y[2])
{
  // A, B, C: 0 2
  //          1 3
  CeedScalar z[2];

  y[0] = C[0] * x[0] + C[2] * x[1];
  y[1] = C[1] * x[0] + C[3] * x[1];

  z[0] = B[0] * y[0] + B[2] * y[1];
  z[1] = B[1] * y[0] + B[3] * y[1];

  y[0] = A[0] * z[0] + A[1] * z[1];
  y[1] = A[2] * z[0] + A[3] * z[1];
}

CEED_QFUNCTION_HELPER void MultBAx22(const CeedScalar A[4], const CeedScalar B[4],
                                     const CeedScalar x[2], CeedScalar y[2])
{
  // A, B: 0 2
  //       1 3
  CeedScalar z[2];

  z[0] = A[0] * x[0] + A[2] * x[1];
  z[1] = A[1] * x[0] + A[3] * x[1];

  y[0] = B[0] * z[0] + B[2] * z[1];
  y[1] = B[1] * z[0] + B[3] * z[1];
}

CEED_QFUNCTION_HELPER void MultAtBA22(const CeedScalar A[4], const CeedScalar B[4],
                                      CeedScalar C[4])
{
  // A, B, C: 0 2
  //          1 3

  // First compute entries of R = B A.
  const CeedScalar R11 = B[0] * A[0] + B[2] * A[1];
  const CeedScalar R21 = B[1] * A[0] + B[3] * A[1];
  const CeedScalar R12 = B[0] * A[2] + B[2] * A[3];
  const CeedScalar R22 = B[1] * A[2] + B[3] * A[3];

  C[0] = A[0] * R11 + A[1] * R21;
  C[1] = A[2] * R11 + A[3] * R21;
  C[2] = A[0] * R12 + A[1] * R22;
  C[3] = A[2] * R12 + A[3] * R22;
}

CEED_QFUNCTION_HELPER void MultAtBC22(const CeedScalar A[4], const CeedScalar B[4],
                                      const CeedScalar C[4], CeedScalar D[4])
{
  // A, B, C, D: 0 2
  //             1 3

  // First compute entries of R = B C.
  const CeedScalar R11 = B[0] * C[0] + B[2] * C[1];
  const CeedScalar R21 = B[1] * C[0] + B[3] * C[1];
  const CeedScalar R12 = B[0] * C[2] + B[2] * C[3];
  const CeedScalar R22 = B[1] * C[2] + B[3] * C[3];

  D[0] = A[0] * R11 + A[1] * R21;
  D[1] = A[2] * R11 + A[3] * R21;
  D[2] = A[0] * R12 + A[1] * R22;
  D[3] = A[2] * R12 + A[3] * R22;
}

CEED_QFUNCTION_HELPER void MultBA22(const CeedScalar A[4], const CeedScalar B[4],
                                    CeedScalar C[4])
{
  // A, B, C: 0 2
  //          1 3
  C[0] = B[0] * A[0] + B[2] * A[1];
  C[1] = B[1] * A[0] + B[3] * A[1];
  C[2] = B[0] * A[2] + B[2] * A[3];
  C[3] = B[1] * A[2] + B[3] * A[3];
}

#endif  // PALACE_LIBCEED_UTILS_22_QF_H
