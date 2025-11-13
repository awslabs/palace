// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_H1D_BUILD_21_QF_H
#define PALACE_LIBCEED_HCURL_H1D_BUILD_21_QF_H

#include "../coeff/coeff_2_qf.h"
#include "utils_21_qf.h"

CEED_QFUNCTION(f_build_hcurlh1d_21)(void *__restrict__ ctx, CeedInt Q,
                                    const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *qdata = in[0];
  CeedScalar *qd = out[0];

  const CeedInt stride = 2 + 2; // [attr, w * |J|, adjJt_11, adjJt_12]_i, lumped by quad point
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    const CeedScalar* adjJt_loc = qdata_i + 2;
    CeedScalar coeff[4], qd_loc[2];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0], coeff);
    MultBA21(adjJt_loc, coeff, qd_loc);

    qd[i + Q * 0] = qdata_i[1] * qd_loc[0];
    qd[i + Q * 1] = qdata_i[1] * qd_loc[1];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_H1D_BUILD_21_QF_H
