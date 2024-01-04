// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_APPLY_QF_H
#define PALACE_LIBCEED_APPLY_QF_H

// libCEED QFunctions for application of a generic operator with assembled quadrature data.
// in[0] is (symmetric) quadrature data, shape [ncomp=vdim*(vdim+1)/2, Q]
// in[1] is active vector, shape [ncomp=vdim, Q]
// out[0] is active vector, shape [ncomp=vdim, Q]

// For pairwise apply functions, the inputs and outputs come in pairs and the quadrature
// data is arranged to be applied with the first vdim*(vdim+1)/2 components for the first
// input/output and the remainder for the second.

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

CEED_QFUNCTION(f_apply_2)(void *, CeedInt Q, const CeedScalar *const *in,
                          CeedScalar *const *out)
{
  const CeedScalar *__restrict__ qd = in[0], *__restrict__ u = in[1];
  CeedScalar *__restrict__ v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = qd[i + Q * 0] * u0 + qd[i + Q * 1] * u1;
    v[i + Q * 1] = qd[i + Q * 1] * u0 + qd[i + Q * 2] * u1;
  }
  return 0;
}

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
    v[i + Q * 0] = qd[i + Q * 0] * u0 + qd[i + Q * 1] * u1 + qd[i + Q * 2] * u2;
    v[i + Q * 1] = qd[i + Q * 1] * u0 + qd[i + Q * 3] * u1 + qd[i + Q * 4] * u2;
    v[i + Q * 2] = qd[i + Q * 2] * u0 + qd[i + Q * 4] * u1 + qd[i + Q * 5] * u2;
  }
  return 0;
}

CEED_QFUNCTION(f_apply_22)(void *, CeedInt Q, const CeedScalar *const *in,
                           CeedScalar *const *out)
{
  const CeedScalar *__restrict__ qd1 = in[0], *__restrict__ qd2 = in[0] + 3 * Q,
                                 *__restrict__ u1 = in[1], *__restrict__ u2 = in[2];
  CeedScalar *__restrict__ v1 = out[0], *__restrict__ v2 = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u10 = u1[i + Q * 0];
    const CeedScalar u11 = u1[i + Q * 1];
    v1[i + Q * 0] = qd1[i + Q * 0] * u10 + qd1[i + Q * 1] * u11;
    v1[i + Q * 1] = qd1[i + Q * 1] * u10 + qd1[i + Q * 2] * u11;

    const CeedScalar u20 = u2[i + Q * 0];
    const CeedScalar u21 = u2[i + Q * 1];
    v2[i + Q * 0] = qd2[i + Q * 0] * u20 + qd2[i + Q * 1] * u21;
    v2[i + Q * 1] = qd2[i + Q * 1] * u20 + qd2[i + Q * 2] * u21;
  }
  return 0;
}

CEED_QFUNCTION(f_apply_33)(void *, CeedInt Q, const CeedScalar *const *in,
                           CeedScalar *const *out)
{
  const CeedScalar *__restrict__ qd1 = in[0], *__restrict__ qd2 = in[0] + 6 * Q,
                                 *__restrict__ u1 = in[1], *__restrict__ u2 = in[2];
  CeedScalar *__restrict__ v1 = out[0], *__restrict__ v2 = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u10 = u1[i + Q * 0];
    const CeedScalar u11 = u1[i + Q * 1];
    const CeedScalar u12 = u1[i + Q * 2];
    v1[i + Q * 0] = qd1[i + Q * 0] * u10 + qd1[i + Q * 1] * u11 + qd1[i + Q * 2] * u12;
    v1[i + Q * 1] = qd1[i + Q * 1] * u10 + qd1[i + Q * 3] * u11 + qd1[i + Q * 4] * u12;
    v1[i + Q * 2] = qd1[i + Q * 2] * u10 + qd1[i + Q * 4] * u11 + qd1[i + Q * 5] * u12;

    const CeedScalar u20 = u2[i + Q * 0];
    const CeedScalar u21 = u2[i + Q * 1];
    const CeedScalar u22 = u2[i + Q * 2];
    v2[i + Q * 0] = qd2[i + Q * 0] * u20 + qd2[i + Q * 1] * u21 + qd2[i + Q * 2] * u22;
    v2[i + Q * 1] = qd2[i + Q * 1] * u20 + qd2[i + Q * 3] * u21 + qd2[i + Q * 4] * u22;
    v2[i + Q * 2] = qd2[i + Q * 2] * u20 + qd2[i + Q * 4] * u21 + qd2[i + Q * 5] * u22;
  }
  return 0;
}

CEED_QFUNCTION(f_apply_12)(void *, CeedInt Q, const CeedScalar *const *in,
                           CeedScalar *const *out)
{
  const CeedScalar *__restrict__ qd1 = in[0], *__restrict__ qd2 = in[0] + Q,
                                 *__restrict__ u1 = in[1], *__restrict__ u2 = in[2];
  CeedScalar *__restrict__ v1 = out[0], *__restrict__ v2 = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v1[i] = qd1[i] * u1[i];

    const CeedScalar u20 = u2[i + Q * 0];
    const CeedScalar u21 = u2[i + Q * 1];
    v2[i + Q * 0] = qd2[i + Q * 0] * u20 + qd2[i + Q * 1] * u21;
    v2[i + Q * 1] = qd2[i + Q * 1] * u20 + qd2[i + Q * 2] * u21;
  }
  return 0;
}

CEED_QFUNCTION(f_apply_13)(void *, CeedInt Q, const CeedScalar *const *in,
                           CeedScalar *const *out)
{
  const CeedScalar *__restrict__ qd1 = in[0], *__restrict__ qd2 = in[0] + Q,
                                 *__restrict__ u1 = in[1], *__restrict__ u2 = in[2];
  CeedScalar *__restrict__ v1 = out[0], *__restrict__ v2 = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v1[i] = qd1[i] * u1[i];

    const CeedScalar u20 = u2[i + Q * 0];
    const CeedScalar u21 = u2[i + Q * 1];
    const CeedScalar u22 = u2[i + Q * 2];
    v2[i + Q * 0] = qd2[i + Q * 0] * u20 + qd2[i + Q * 1] * u21 + qd2[i + Q * 2] * u22;
    v2[i + Q * 1] = qd2[i + Q * 1] * u20 + qd2[i + Q * 3] * u21 + qd2[i + Q * 4] * u22;
    v2[i + Q * 2] = qd2[i + Q * 2] * u20 + qd2[i + Q * 4] * u21 + qd2[i + Q * 5] * u22;
  }
  return 0;
}

CEED_QFUNCTION(f_apply_21)(void *, CeedInt Q, const CeedScalar *const *in,
                           CeedScalar *const *out)
{
  const CeedScalar *__restrict__ qd1 = in[0], *__restrict__ qd2 = in[0] + 3 * Q,
                                 *__restrict__ u1 = in[1], *__restrict__ u2 = in[2];
  CeedScalar *__restrict__ v1 = out[0], *__restrict__ v2 = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u10 = u1[i + Q * 0];
    const CeedScalar u11 = u1[i + Q * 1];
    v1[i + Q * 0] = qd1[i + Q * 0] * u10 + qd1[i + Q * 1] * u11;
    v1[i + Q * 1] = qd1[i + Q * 1] * u10 + qd1[i + Q * 2] * u11;

    v2[i] = qd2[i] * u2[i];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_31)(void *, CeedInt Q, const CeedScalar *const *in,
                           CeedScalar *const *out)
{
  const CeedScalar *__restrict__ qd1 = in[0], *__restrict__ qd2 = in[0] + 6 * Q,
                                 *__restrict__ u1 = in[1], *__restrict__ u2 = in[2];
  CeedScalar *__restrict__ v1 = out[0], *__restrict__ v2 = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u10 = u1[i + Q * 0];
    const CeedScalar u11 = u1[i + Q * 1];
    const CeedScalar u12 = u1[i + Q * 2];
    v1[i + Q * 0] = qd1[i + Q * 0] * u10 + qd1[i + Q * 1] * u11 + qd1[i + Q * 2] * u12;
    v1[i + Q * 1] = qd1[i + Q * 1] * u10 + qd1[i + Q * 3] * u11 + qd1[i + Q * 4] * u12;
    v1[i + Q * 2] = qd1[i + Q * 2] * u10 + qd1[i + Q * 4] * u11 + qd1[i + Q * 5] * u12;

    v2[i] = qd2[i] * u2[i];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_APPLY_QF_H
