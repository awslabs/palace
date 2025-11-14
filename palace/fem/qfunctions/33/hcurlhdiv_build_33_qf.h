// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_HDIV_BUILD_33_QF_H
#define PALACE_LIBCEED_HCURL_HDIV_BUILD_33_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

CEED_QFUNCTION(f_build_hcurlhdiv_33)(void *__restrict__ ctx, CeedInt Q,
                                     const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *qdata = in[0];
  CeedScalar *qd = out[0];

  const CeedInt stride = 2 + 9; // attr, w * |J|, (adjJt / |J|) colwise
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    const CeedScalar* adjJt_loc = qdata_i + 2;
    CeedScalar coeff[9], J_loc[9], qd_loc[9];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0], coeff);
    AdjJt33(adjJt_loc, J_loc);
    MultAtBC33(J_loc, coeff, adjJt_loc, qd_loc);

    qd[i + Q * 0] = qdata_i[1] * qd_loc[0];
    qd[i + Q * 1] = qdata_i[1] * qd_loc[1];
    qd[i + Q * 2] = qdata_i[1] * qd_loc[2];
    qd[i + Q * 3] = qdata_i[1] * qd_loc[3];
    qd[i + Q * 4] = qdata_i[1] * qd_loc[4];
    qd[i + Q * 5] = qdata_i[1] * qd_loc[5];
    qd[i + Q * 6] = qdata_i[1] * qd_loc[6];
    qd[i + Q * 7] = qdata_i[1] * qd_loc[7];
    qd[i + Q * 8] = qdata_i[1] * qd_loc[8];
  }
  return 0;
}

CEED_QFUNCTION(f_build_hdivhcurl_33)(void *__restrict__ ctx, CeedInt Q,
                                     const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *qdata = in[0];
  CeedScalar *qd = out[0];

  const CeedInt stride = 2 + 9; // attr, w * |J|, (adjJt / |J|) colwise
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    const CeedScalar* adjJt_loc = qdata_i + 2;
    CeedScalar coeff[9], J_loc[9], qd_loc[9];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0], coeff);
    AdjJt33(adjJt_loc, J_loc);
    MultAtBC33(adjJt_loc, coeff, J_loc, qd_loc);

    qd[i + Q * 0] = qdata_i[1] * qd_loc[0];
    qd[i + Q * 1] = qdata_i[1] * qd_loc[1];
    qd[i + Q * 2] = qdata_i[1] * qd_loc[2];
    qd[i + Q * 3] = qdata_i[1] * qd_loc[3];
    qd[i + Q * 4] = qdata_i[1] * qd_loc[4];
    qd[i + Q * 5] = qdata_i[1] * qd_loc[5];
    qd[i + Q * 6] = qdata_i[1] * qd_loc[6];
    qd[i + Q * 7] = qdata_i[1] * qd_loc[7];
    qd[i + Q * 8] = qdata_i[1] * qd_loc[8];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_HDIV_BUILD_33_QF_H
