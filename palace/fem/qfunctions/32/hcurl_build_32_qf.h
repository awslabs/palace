// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_BUILD_32_QF_H
#define PALACE_LIBCEED_HCURL_BUILD_32_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_32_qf.h"

CEED_QFUNCTION(f_build_hcurl_32)(void *__restrict__ ctx, CeedInt Q,
                                 const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *qdata = in[0];
  CeedScalar *qd = out[0];

  const CeedInt stride = 2 + 6; // attr, w * |J|, (adjJt / |J|) colwise
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    const CeedScalar* adjJt_loc = qdata_i + 2;
    CeedScalar coeff[9], qd_loc[4];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0], coeff);
    MultAtBA32(adjJt_loc, coeff, qd_loc);

    qd[i + Q * 0] = qdata_i[1] * qd_loc[0];
    qd[i + Q * 1] = qdata_i[1] * qd_loc[1];
    qd[i + Q * 2] = qdata_i[1] * qd_loc[2];
    qd[i + Q * 3] = qdata_i[1] * qd_loc[3];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_BUILD_32_QF_H
