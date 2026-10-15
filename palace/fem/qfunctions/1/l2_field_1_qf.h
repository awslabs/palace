// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_FIELD_1_QF_H
#define PALACE_LIBCEED_L2_FIELD_1_QF_H

#include "../coeff/coeff_1_qf.h"

// L2 (2D curl-curl) operator with a spatially-varying scalar material coefficient sampled
// at the quadrature points. Inputs are in[0] geometry data, in[1] quadrature weights,
// in[2] the active vector, and in[3] the coefficient field, shape [1, Q].

CEED_QFUNCTION(f_apply_l2_field_1)(void *__restrict__ ctx, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *qw = in[1], *u = in[2],
                   *coeff_data = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

    v[i] = (coeff * coeff_data[i] * qw[i] * qw[i] / wdetJ[i]) * u[i];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_FIELD_1_QF_H
