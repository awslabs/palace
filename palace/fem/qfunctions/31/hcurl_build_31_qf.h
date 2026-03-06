// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_BUILD_31_QF_H
#define PALACE_LIBCEED_HCURL_BUILD_31_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_31_qf.h"

CEED_QFUNCTION(f_build_hcurl_31)(void *__restrict__ ctx, CeedInt Q,
                                 const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q;
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[9], adjJt_loc[3], qd_loc[1];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack31(adjJt + i, Q, adjJt_loc);
    MultAtBA31(adjJt_loc, coeff, qd_loc);

    qd[i + Q * 0] = wdetJ[i] * qd_loc[0];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_BUILD_31_QF_H
