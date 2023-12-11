// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_BUILD_QF_H
#define PALACE_LIBCEED_L2_BUILD_QF_H

#include "types_qf.h"

// Build functions replace active vector output with quadrature point data, stored as a
// symmetric matrix, and remove active vector input.

CEED_QFUNCTION(f_build_l2_1)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *qw = in[1], *attr = in[2];
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack((const MatCoeffContext1 *)ctx, (CeedInt)attr[i]);

    qd[i] = coeff * qw[i] * qw[i] / wdetJ[i];
  }
  return 0;
}

CEED_QFUNCTION(f_build_l2_2)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *qw = in[1], *attr = in[2];
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[3];
    CoeffUnpack((const MatCoeffContext2 *)ctx, (CeedInt)attr[i], coeff);
    const CeedScalar w = qw[i] * qw[i] / wdetJ[i];

    qd[i + Q * 0] = w * coeff[0];
    qd[i + Q * 1] = w * coeff[1];
    qd[i + Q * 2] = w * coeff[2];
  }
  return 0;
}

CEED_QFUNCTION(f_build_l2_3)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *qw = in[1], *attr = in[2];
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[6];
    CoeffUnpack((const MatCoeffContext3 *)ctx, (CeedInt)attr[i], coeff);
    const CeedScalar w = qw[i] * qw[i] / wdetJ[i];

    qd[i + Q * 0] = w * coeff[0];
    qd[i + Q * 1] = w * coeff[1];
    qd[i + Q * 2] = w * coeff[2];
    qd[i + Q * 3] = w * coeff[3];
    qd[i + Q * 4] = w * coeff[4];
    qd[i + Q * 5] = w * coeff[5];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_BUILD_QF_H
