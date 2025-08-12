// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_H1D_32_QF_H
#define PALACE_LIBCEED_HCURL_H1D_32_QF_H

#include "../coeff/coeff_1_qf.h"
#include "utils_32_qf.h"

CEED_QFUNCTION(f_apply_hcurlh1d_32)(void *__restrict__ ctx, CeedInt Q,
                                    const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar adjJt_loc[6], v_loc[3];
    const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);
    MatUnpack32(adjJt + i, Q, adjJt_loc);
    MultBAx32(adjJt_loc, coeff, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_H1D_32_QF_H
