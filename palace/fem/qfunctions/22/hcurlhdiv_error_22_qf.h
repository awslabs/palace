// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURLHDIV_ERROR_22_QF_H
#define PALACE_LIBCEED_HCURLHDIV_ERROR_22_QF_H

#include "../coeff/coeff_2_qf.h"
#include "utils_22_qf.h"

CEED_QFUNCTION(f_apply_hcurlhdiv_error_22)(void *__restrict__ ctx, CeedInt Q,
                                           const CeedScalar *const *in,
                                           CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u1 = in[1],
                   *u2 = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar adjJt_loc[4], v1_loc[2], v2_loc[2];
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    {
      const CeedScalar u1_loc[2] = {u1[i + Q * 0], u1[i + Q * 1]};
      CeedScalar coeff[4];
      CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MultBAx22(adjJt_loc, coeff, u1_loc, v1_loc);
    }
    {
      const CeedScalar u2_loc[2] = {u2[i + Q * 0], u2[i + Q * 1]};
      CeedScalar coeff[4], J_loc[4];
      CoeffUnpack2(CoeffPairSecond<2>((const CeedIntScalar *)ctx), (CeedInt)attr[i], coeff);
      AdjJt22(adjJt_loc, J_loc);
      MultBAx22(J_loc, coeff, u2_loc, v2_loc);
    }
    v2_loc[0] -= v1_loc[0];
    v2_loc[1] -= v1_loc[1];
    v[i] = wdetJ[i] * (v2_loc[0] * v2_loc[0] + v2_loc[1] * v2_loc[1]);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivhcurl_error_22)(void *__restrict__ ctx, CeedInt Q,
                                           const CeedScalar *const *in,
                                           CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u1 = in[1],
                   *u2 = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar adjJt_loc[4], v1_loc[2], v2_loc[2];
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    {
      const CeedScalar u1_loc[2] = {u1[i + Q * 0], u1[i + Q * 1]};
      CeedScalar coeff[4], J_loc[4];
      CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      AdjJt22(adjJt_loc, J_loc);
      MultBAx22(J_loc, coeff, u1_loc, v1_loc);
    }
    {
      const CeedScalar u2_loc[2] = {u2[i + Q * 0], u2[i + Q * 1]};
      CeedScalar coeff[4];
      CoeffUnpack2(CoeffPairSecond<2>((const CeedIntScalar *)ctx), (CeedInt)attr[i], coeff);
      MultBAx22(adjJt_loc, coeff, u2_loc, v2_loc);
    }
    v2_loc[0] -= v1_loc[0];
    v2_loc[1] -= v1_loc[1];
    v[i] = wdetJ[i] * (v2_loc[0] * v2_loc[0] + v2_loc[1] * v2_loc[1]);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURLHDIV_ERROR_22_QF_H
