// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HDIV_BUILD_22_QF_H
#define PALACE_LIBCEED_HDIV_BUILD_22_QF_H

#include "../coeff/coeff_2_qf.h"
#include "utils_22_qf.h"

CEED_QFUNCTION(f_build_hdiv_22)(void *__restrict__ ctx, CeedInt Q,
                                const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q;
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[4], adjJt_loc[4], J_loc[4], qd_loc[4];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    AdjJt22(adjJt_loc, J_loc);
    MultAtBA22(J_loc, coeff, qd_loc);

    qd[i + Q * 0] = wdetJ[i] * qd_loc[0];
    qd[i + Q * 1] = wdetJ[i] * qd_loc[1];
    qd[i + Q * 2] = wdetJ[i] * qd_loc[2];
    qd[i + Q * 3] = wdetJ[i] * qd_loc[3];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HDIV_BUILD_22_QF_H
