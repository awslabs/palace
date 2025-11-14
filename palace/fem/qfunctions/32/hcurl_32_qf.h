// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_32_QF_H
#define PALACE_LIBCEED_HCURL_32_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_32_qf.h"

CEED_QFUNCTION(f_apply_hcurl_32)(void *__restrict__ ctx, CeedInt Q,
                                 const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *qdata = in[0], *u = in[1];
  CeedScalar *v = out[0];

  const CeedInt stride = 2 + 6; // attr, w * |J|, (adjJt / |J|) colwise
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    const CeedScalar* adjJt_loc = qdata_i + 2;
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar coeff[9], v_loc[2];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0], coeff);
    MultAtBCx32(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

    v[i + Q * 0] = qdata_i[1] * v_loc[0];
    v[i + Q * 1] = qdata_i[1] * v_loc[1];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_32_QF_H
