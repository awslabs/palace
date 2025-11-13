// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_H1D_21_QF_H
#define PALACE_LIBCEED_HCURL_H1D_21_QF_H

#include "../coeff/coeff_2_qf.h"
#include "utils_21_qf.h"

CEED_QFUNCTION(f_apply_hcurlh1d_21)(void *__restrict__ ctx, CeedInt Q,
                                    const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *qdata = in[0], *u = in[1];
  CeedScalar *v = out[0];

  const CeedInt stride = 2 + 2; // [attr, w * |J|, adjJt_11, adjJt_12]_i, lumped by quad point
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    const CeedScalar* adjJt_loc = qdata_i + 2;
    const CeedScalar u_loc[1] = {u[i + Q * 0]};
    CeedScalar coeff[4], adjJt_loc[2], v_loc[2];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0], coeff);
    MultBAx21(adjJt_loc, coeff, u_loc, v_loc);

    v[i + Q * 0] = qdata_i[1] * v_loc[0];
    v[i + Q * 1] = qdata_i[1] * v_loc[1];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_H1D_21_QF_H
