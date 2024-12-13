// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURLH1D_ERROR_33_QF_H
#define PALACE_LIBCEED_HCURLH1D_ERROR_33_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

CEED_QFUNCTION(f_apply_hcurlh1d_error_33)(void *__restrict__ ctx, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u1 = in[1],
                   *u2 = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar v1_loc[3], v2_loc[3];
    {
      const CeedScalar u1_loc[3] = {u1[i + Q * 0], u1[i + Q * 1], u1[i + Q * 2]};
      CeedScalar coeff[9], adjJt_loc[9];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack33(adjJt + i, Q, adjJt_loc);
      MultBAx33(adjJt_loc, coeff, u1_loc, v1_loc);
    }
    {
      const CeedScalar u2_loc[3] = {u2[i + Q * 0], u2[i + Q * 1], u2[i + Q * 2]};
      CeedScalar coeff[9];
      CoeffUnpack3(CoeffPairSecond<3>((const CeedIntScalar *)ctx), (CeedInt)attr[i], coeff);
      MultBx33(coeff, u2_loc, v2_loc);
    }
    v2_loc[0] -= v1_loc[0];
    v2_loc[1] -= v1_loc[1];
    v2_loc[2] -= v1_loc[2];
    v[i] =
        wdetJ[i] * (v2_loc[0] * v2_loc[0] + v2_loc[1] * v2_loc[1] + v2_loc[2] * v2_loc[2]);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_h1dhcurl_error_33)(void *__restrict__ ctx, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u1 = in[1],
                   *u2 = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar v1_loc[3], v2_loc[3];
    {
      const CeedScalar u1_loc[3] = {u1[i + Q * 0], u1[i + Q * 1], u1[i + Q * 2]};
      CeedScalar coeff[9];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MultBx33(coeff, u1_loc, v1_loc);
    }
    {
      const CeedScalar u2_loc[3] = {u2[i + Q * 0], u2[i + Q * 1], u2[i + Q * 2]};
      CeedScalar coeff[9], adjJt_loc[9];
      CoeffUnpack3(CoeffPairSecond<3>((const CeedIntScalar *)ctx), (CeedInt)attr[i], coeff);
      MatUnpack33(adjJt + i, Q, adjJt_loc);
      MultBAx33(adjJt_loc, coeff, u2_loc, v2_loc);
    }
    v2_loc[0] -= v1_loc[0];
    v2_loc[1] -= v1_loc[1];
    v2_loc[2] -= v1_loc[2];
    v[i] =
        wdetJ[i] * (v2_loc[0] * v2_loc[0] + v2_loc[1] * v2_loc[1] + v2_loc[2] * v2_loc[2]);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURLH1D_ERROR_33_QF_H
