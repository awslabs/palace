// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_H1_FIELD_1_QF_H
#define PALACE_LIBCEED_H1_FIELD_1_QF_H

#include "../coeff/coeff_1_qf.h"

// H1 mass operator with a spatially-varying scalar material coefficient sampled at the
// quadrature points. in[2] is the coefficient field, shape [1, Q]. The field multiplies
// the (possibly identity) attribute-based coefficient from ctx, so the caller can either
// weight the material coefficient or replace it by passing an identity coefficient.

CEED_QFUNCTION(f_apply_h1_field_1)(void *__restrict__ ctx, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *u = in[1], *coeff_data = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

    v[i] = coeff * coeff_data[i] * wdetJ[i] * u[i];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_H1_FIELD_1_QF_H
