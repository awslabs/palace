// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_MASS_BUILD_21_QF_H
#define PALACE_LIBCEED_HCURL_MASS_BUILD_21_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_2_qf.h"
#include "utils_21_qf.h"

CEED_QFUNCTION(f_build_hcurlmass_21)(void *__restrict__ ctx, CeedInt Q,
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
      CeedScalar adjJt_loc[2], qd_loc[1];
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond<1>((const CeedIntScalar *)ctx), (CeedInt)attr[i]);
      MatUnpack21(adjJt + i, Q, adjJt_loc);
      MultAtBA21(adjJt_loc, coeff, qd_loc);

      qd2[i + Q * 0] = wdetJ[i] * qd_loc[0];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_MASS_BUILD_21_QF_H
