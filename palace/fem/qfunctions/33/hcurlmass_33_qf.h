// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_MASS_33_QF_H
#define PALACE_LIBCEED_HCURL_MASS_33_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

CEED_QFUNCTION(f_apply_hcurlmass_33)(void *__restrict__ ctx, CeedInt Q,
                                     const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *gradu = in[2];
  CeedScalar *__restrict__ v = out[0], *__restrict__ gradv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

      v[i] = coeff * wdetJ[i] * u[i];
    }
    {
      const CeedScalar u_loc[3] = {gradu[i + Q * 0], gradu[i + Q * 1], gradu[i + Q * 2]};
      CeedScalar coeff[9], adjJt_loc[9], v_loc[3];
      CoeffUnpack3(CoeffPairSecond<1>((const CeedIntScalar *)ctx), (CeedInt)attr[i], coeff);
      MatUnpack33(adjJt + i, Q, adjJt_loc);
      MultAtBCx33(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

      gradv[i + Q * 0] = wdetJ[i] * v_loc[0];
      gradv[i + Q * 1] = wdetJ[i] * v_loc[1];
      gradv[i + Q * 2] = wdetJ[i] * v_loc[2];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_MASS_33_QF_H
