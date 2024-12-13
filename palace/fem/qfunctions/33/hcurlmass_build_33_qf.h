// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_MASS_BUILD_33_QF_H
#define PALACE_LIBCEED_HCURL_MASS_BUILD_33_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

CEED_QFUNCTION(f_build_hcurlmass_33)(void *__restrict__ ctx, CeedInt Q,
                                     const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q;
  CeedScalar *__restrict__ qd1 = out[0], *__restrict__ qd2 = out[0] + Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

      qd1[i + Q * 0] = coeff * wdetJ[i];
    }
    {
      CeedScalar coeff[9], adjJt_loc[9], qd_loc[9];
      CoeffUnpack3(CoeffPairSecond<1>((const CeedIntScalar *)ctx), (CeedInt)attr[i], coeff);
      MatUnpack33(adjJt + i, Q, adjJt_loc);
      MultAtBA33(adjJt_loc, coeff, qd_loc);

      qd2[i + Q * 0] = wdetJ[i] * qd_loc[0];
      qd2[i + Q * 1] = wdetJ[i] * qd_loc[1];
      qd2[i + Q * 2] = wdetJ[i] * qd_loc[2];
      qd2[i + Q * 3] = wdetJ[i] * qd_loc[3];
      qd2[i + Q * 4] = wdetJ[i] * qd_loc[4];
      qd2[i + Q * 5] = wdetJ[i] * qd_loc[5];
      qd2[i + Q * 6] = wdetJ[i] * qd_loc[6];
      qd2[i + Q * 7] = wdetJ[i] * qd_loc[7];
      qd2[i + Q * 8] = wdetJ[i] * qd_loc[8];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_MASS_BUILD_33_QF_H
