// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2H1_1_QF_H
#define PALACE_LIBCEED_L2H1_1_QF_H

#include "../coeff/coeff_1_qf.h"

// Mixed mass with an INTEGRAL trial space (u = û / |J|) and a VALUE test space (v = v̂):
// the |J| of the quadrature measure cancels, leaving the quadrature weight.
CEED_QFUNCTION(f_apply_l2h1_1)(void *__restrict__ ctx, CeedInt Q,
                               const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *qw = in[1], *u = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

    v[i] = (coeff * qw[i]) * u[i];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2H1_1_QF_H
