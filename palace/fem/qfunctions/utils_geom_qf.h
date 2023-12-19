// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_UTILS_GEOM_QF_H
#define PALACE_LIBCEED_UTILS_GEOM_QF_H

#include <math.h>

CEED_QFUNCTION_HELPER CeedScalar DetJ22(const CeedScalar J[4])
{
  // J: 0 2
  //    1 3
  return J[0] * J[3] - J[1] * J[2];
}

CEED_QFUNCTION_HELPER CeedScalar DetJ33(const CeedScalar J[9])
{
  // J: 0 3 6
  //    1 4 7
  //    2 5 8
  return J[0] * (J[4] * J[8] - J[5] * J[7]) - J[1] * (J[3] * J[8] - J[5] * J[6]) +
         J[2] * (J[3] * J[7] - J[4] * J[6]);
}

CEED_QFUNCTION_HELPER CeedScalar DetJ21(const CeedScalar J[2])
{
  // J: 0
  //    1
  return sqrt(J[0] * J[0] + J[1] * J[1]);
}

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

template <bool ComputeDet>
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

template <bool ComputeDet>
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

template <bool ComputeDet>
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

template <bool ComputeDet>
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

#endif  // PALACE_LIBCEED_UTILS_QF_H
