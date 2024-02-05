// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_APPLY_1_QF_H
#define PALACE_LIBCEED_APPLY_1_QF_H

CEED_QFUNCTION(f_apply_1)(void *, CeedInt Q, const CeedScalar *const *in,
                          CeedScalar *const *out)
{
  const CeedScalar *__restrict__ qd = in[0], *__restrict__ u = in[1];
  CeedScalar *__restrict__ v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = qd[i] * u[i];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_APPLY_1_QF_H
