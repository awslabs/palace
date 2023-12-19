// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_H1_BUILD_QF_H
#define PALACE_LIBCEED_H1_BUILD_QF_H

#include "coeff_qf.h"

// Build functions replace active vector output with quadrature point data, stored as a
// symmetric matrix, and remove active vector input.

CEED_QFUNCTION(f_build_h1_1)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q;
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

    qd[i] = coeff * wdetJ[i];
  }
  return 0;
}

CEED_QFUNCTION(f_build_h1_2)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q;
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[3];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);

    qd[i + Q * 0] = wdetJ[i] * coeff[0];
    qd[i + Q * 1] = wdetJ[i] * coeff[1];
    qd[i + Q * 2] = wdetJ[i] * coeff[2];
  }
  return 0;
}

CEED_QFUNCTION(f_build_h1_3)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q;
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[6];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);

    qd[i + Q * 0] = wdetJ[i] * coeff[0];
    qd[i + Q * 1] = wdetJ[i] * coeff[1];
    qd[i + Q * 2] = wdetJ[i] * coeff[2];
    qd[i + Q * 3] = wdetJ[i] * coeff[3];
    qd[i + Q * 4] = wdetJ[i] * coeff[4];
    qd[i + Q * 5] = wdetJ[i] * coeff[5];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_H1_BUILD_QF_H
