// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_H1D_BUILD_32_QF_H
#define PALACE_LIBCEED_HCURL_H1D_BUILD_32_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_32_qf.h"

CEED_QFUNCTION(f_build_hcurlh1d_32)(void *__restrict__ ctx, CeedInt Q,
                                    const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q;
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[9], adjJt_loc[6], qd_loc[6];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack32(adjJt + i, Q, adjJt_loc);
    MultBA32(adjJt_loc, coeff, qd_loc);

    qd[i + Q * 0] = wdetJ[i] * qd_loc[0];
    qd[i + Q * 1] = wdetJ[i] * qd_loc[1];
    qd[i + Q * 2] = wdetJ[i] * qd_loc[2];
    qd[i + Q * 3] = wdetJ[i] * qd_loc[3];
    qd[i + Q * 4] = wdetJ[i] * qd_loc[4];
    qd[i + Q * 5] = wdetJ[i] * qd_loc[5];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_H1D_BUILD_32_QF_H
