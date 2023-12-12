// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_QF_H
#define PALACE_LIBCEED_L2_QF_H

#include "coeff_qf.h"

// libCEED QFunctions for L2 operators (Piola transformation u = 1 / det(J) ̂u).
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is quadrature weights, shape [Q]
// in[2] is element attribute, shape [Q]
// in[3] is active vector, shape [ncomp=vdim, Q]
// out[0] is active vector, shape [ncomp=vdim, Q]

CEED_QFUNCTION(f_apply_l2_1)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *qw = in[1], *attr = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

    v[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * u[i];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2_2)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *qw = in[1], *attr = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[3];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    const CeedScalar w = qw[i] * qw[i] / wdetJ[i];

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = w * (coeff[0] * u0 + coeff[1] * u1);
    v[i + Q * 1] = w * (coeff[1] * u0 + coeff[2] * u1);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2_3)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *qw = in[1], *attr = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[6];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    const CeedScalar w = qw[i] * qw[i] / wdetJ[i];

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    const CeedScalar u2 = u[i + Q * 2];
    v[i + Q * 0] = w * (coeff[0] * u0 + coeff[1] * u1 + coeff[2] * u2);
    v[i + Q * 1] = w * (coeff[1] * u0 + coeff[3] * u1 + coeff[4] * u2);
    v[i + Q * 2] = w * (coeff[2] * u0 + coeff[4] * u1 + coeff[5] * u2);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_QF_H
