// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURLHDIV_ERROR_33_QF_H
#define PALACE_LIBCEED_HCURLHDIV_ERROR_33_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

CEED_QFUNCTION(f_apply_hcurlhdiv_error_33)(void *__restrict__ ctx, CeedInt Q,
                                           const CeedScalar *const *in,
                                           CeedScalar *const *out)
{
  const CeedScalar *qdata = in[0],*u1 = in[1], *u2 = in[2];
  CeedScalar *v = out[0];

  const CeedInt stride = 2 + 9; // attr, w * |J|, (adjJt / |J|) colwise
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    const CeedScalar* adjJt_loc = qdata_i + 2;
    CeedScalar v1_loc[3], v2_loc[3];
    {
      const CeedScalar u1_loc[3] = {u1[i + Q * 0], u1[i + Q * 1], u1[i + Q * 2]};
      CeedScalar coeff[9];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0], coeff);
      MultBAx33(adjJt_loc, coeff, u1_loc, v1_loc);
    }
    {
      const CeedScalar u2_loc[3] = {u2[i + Q * 0], u2[i + Q * 1], u2[i + Q * 2]};
      CeedScalar coeff[9], J_loc[9];
      CoeffUnpack3(CoeffPairSecond<3>((const CeedIntScalar *)ctx), (CeedInt)qdata_i[0], coeff);
      AdjJt33(adjJt_loc, J_loc);
      MultBAx33(J_loc, coeff, u2_loc, v2_loc);
    }
    v2_loc[0] -= v1_loc[0];
    v2_loc[1] -= v1_loc[1];
    v2_loc[2] -= v1_loc[2];
    v[i] =
        qdata_i[1] * (v2_loc[0] * v2_loc[0] + v2_loc[1] * v2_loc[1] + v2_loc[2] * v2_loc[2]);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivhcurl_error_33)(void *__restrict__ ctx, CeedInt Q,
                                           const CeedScalar *const *in,
                                           CeedScalar *const *out)
{
  const CeedScalar *qdata = in[0], *u1 = in[1], *u2 = in[2];
  CeedScalar *v = out[0];

  const CeedInt stride = 2 + 9; // attr, w * |J|, (adjJt / |J|) colwise
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    const CeedScalar* adjJt_loc = qdata_i + 2;
    CeedScalar v1_loc[3], v2_loc[3];
    {
      const CeedScalar u1_loc[3] = {u1[i + Q * 0], u1[i + Q * 1], u1[i + Q * 2]};
      CeedScalar coeff[9], J_loc[9];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0], coeff);
      AdjJt33(adjJt_loc, J_loc);
      MultBAx33(J_loc, coeff, u1_loc, v1_loc);
    }
    {
      const CeedScalar u2_loc[3] = {u2[i + Q * 0], u2[i + Q * 1], u2[i + Q * 2]};
      CeedScalar coeff[9];
      CoeffUnpack3(CoeffPairSecond<3>((const CeedIntScalar *)ctx), (CeedInt)qdata_i[0], coeff);
      MultBAx33(adjJt_loc, coeff, u2_loc, v2_loc);
    }
    v2_loc[0] -= v1_loc[0];
    v2_loc[1] -= v1_loc[1];
    v2_loc[2] -= v1_loc[2];
    v[i] =
        qdata_i[1] * (v2_loc[0] * v2_loc[0] + v2_loc[1] * v2_loc[1] + v2_loc[2] * v2_loc[2]);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURLHDIV_ERROR_33_QF_H
