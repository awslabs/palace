// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_APPLY_3_QF_H
#define PALACE_LIBCEED_APPLY_3_QF_H

CEED_QFUNCTION(f_apply_3)(void *, CeedInt Q, const CeedScalar *const *in,
                          CeedScalar *const *out)
{
  const CeedScalar *__restrict__ qd = in[0], *__restrict__ u = in[1];
  CeedScalar *__restrict__ v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    const CeedScalar u2 = u[i + Q * 2];
    v[i + Q * 0] = qd[i + Q * 0] * u0 + qd[i + Q * 3] * u1 + qd[i + Q * 6] * u2;
    v[i + Q * 1] = qd[i + Q * 1] * u0 + qd[i + Q * 4] * u1 + qd[i + Q * 7] * u2;
    v[i + Q * 2] = qd[i + Q * 2] * u0 + qd[i + Q * 5] * u1 + qd[i + Q * 8] * u2;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_APPLY_3_QF_H
