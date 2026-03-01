// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_UTILS_31_QF_H
#define PALACE_LIBCEED_UTILS_31_QF_H

#ifndef CEED_RUNNING_JIT_PASS
#include <math.h>
#endif

CEED_QFUNCTION_HELPER CeedScalar DetJ31(const CeedScalar J[3])
{
  // J: 0
  //    1
  //    2
  return sqrt(J[0] * J[0] + J[1] * J[1] + J[2] * J[2]);
}

template <bool ComputeDet = false>
CEED_QFUNCTION_HELPER CeedScalar AdjJt31(const CeedScalar J[3], CeedScalar adjJt[3])
{
  // Compute adj(J)^T / det(J) and store the result.
  // J: 0   adj(J): 1/sqrt(J^T J) J^T
  //    1
  //    2
  const CeedScalar d = sqrt(J[0] * J[0] + J[1] * J[1] + J[2] * J[2]);
  adjJt[0] = J[0] / d;
  adjJt[1] = J[1] / d;
  adjJt[2] = J[2] / d;
  return ComputeDet ? d : 0.0;
}

CEED_QFUNCTION_HELPER void MatUnpack31(const CeedScalar *A, const CeedInt A_stride,
                                       CeedScalar A_loc[3])
{
  A_loc[0] = A[A_stride * 0];
  A_loc[1] = A[A_stride * 1];
  A_loc[2] = A[A_stride * 2];
}

CEED_QFUNCTION_HELPER void MultAtBCx31(const CeedScalar A[3], const CeedScalar B[9],
                                       const CeedScalar C[3], const CeedScalar x[1],
                                       CeedScalar y[1])
{
  // A: 0   B: 0 3 6   C: 0
  //    1      1 4 7      1
  //    2      2 5 8      2
  CeedScalar z[3], t;

  y[0] = C[0] * x[0];
  t = C[1] * x[0];
  const CeedScalar s = C[2] * x[0];

  z[0] = B[0] * y[0] + B[3] * t + B[6] * s;
  z[1] = B[1] * y[0] + B[4] * t + B[7] * s;
  z[2] = B[2] * y[0] + B[5] * t + B[8] * s;

  y[0] = A[0] * z[0] + A[1] * z[1] + A[2] * z[2];
}

CEED_QFUNCTION_HELPER void MultBAx31(const CeedScalar A[3], const CeedScalar B[9],
                                     const CeedScalar x[1], CeedScalar y[3])
{
  // A: 0   B: 0 3 6
  //    1      1 4 7
  //    2      2 5 8
  CeedScalar z[3];

  z[0] = A[0] * x[0];
  z[1] = A[1] * x[0];
  z[2] = A[2] * x[0];

  y[0] = B[0] * z[0] + B[3] * z[1] + B[6] * z[2];
  y[1] = B[1] * z[0] + B[4] * z[1] + B[7] * z[2];
  y[2] = B[2] * z[0] + B[5] * z[1] + B[8] * z[2];
}

CEED_QFUNCTION_HELPER void MultAtBA31(const CeedScalar A[3], const CeedScalar B[9],
                                      CeedScalar C[1])
{
  // A: 0   B: 0 3 6   C: 0
  //    1      1 4 7
  //    2      2 5 8

  // First compute entries of R = B A.
  const CeedScalar R11 = B[0] * A[0] + B[3] * A[1] + B[6] * A[2];
  const CeedScalar R21 = B[1] * A[0] + B[4] * A[1] + B[7] * A[2];
  const CeedScalar R31 = B[2] * A[0] + B[5] * A[1] + B[8] * A[2];

  C[0] = A[0] * R11 + A[1] * R21 + A[2] * R31;
}

CEED_QFUNCTION_HELPER void MultAtBC31(const CeedScalar A[3], const CeedScalar B[9],
                                      const CeedScalar C[3], CeedScalar D[1])
{
  // A, C: 0   B: 0 3 6   D: 0
  //       1      1 4 7
  //       2      2 5 8

  // First compute entries of R = B C.
  const CeedScalar R11 = B[0] * C[0] + B[3] * C[1] + B[6] * C[2];
  const CeedScalar R21 = B[1] * C[0] + B[4] * C[1] + B[7] * C[2];
  const CeedScalar R31 = B[2] * C[0] + B[5] * C[1] + B[8] * C[2];

  D[0] = A[0] * R11 + A[1] * R21 + A[2] * R31;
}

CEED_QFUNCTION_HELPER void MultBA31(const CeedScalar A[3], const CeedScalar B[9],
                                    CeedScalar C[3])
{
  // A: 0   B: 0 3 6   C: 0
  //    1      1 4 7      1
  //    2      2 5 8      2
  C[0] = B[0] * A[0] + B[3] * A[1] + B[6] * A[2];
  C[1] = B[1] * A[0] + B[4] * A[1] + B[7] * A[2];
  C[2] = B[2] * A[0] + B[5] * A[1] + B[8] * A[2];
}

#endif  // PALACE_LIBCEED_UTILS_31_QF_H
