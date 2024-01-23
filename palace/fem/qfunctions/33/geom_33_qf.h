// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_GEOM_33_QF_H
#define PALACE_LIBCEED_GEOM_33_QF_H

#include "utils_33_qf.h"

CEED_QFUNCTION(f_build_geom_factor_33)(void *, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  const CeedScalar *qw = in[0], *J = in[1];
  CeedScalar *attr = out[0], *wdetJ = out[0] + Q, *adjJt = out[0] + 2 * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_loc[9], adjJt_loc[9];
    MatUnpack33(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt33<true>(J_loc, adjJt_loc);

    attr[i] = 0;
    wdetJ[i] = qw[i] * detJ;
    adjJt[i + Q * 0] = adjJt_loc[0] / detJ;
    adjJt[i + Q * 1] = adjJt_loc[1] / detJ;
    adjJt[i + Q * 2] = adjJt_loc[2] / detJ;
    adjJt[i + Q * 3] = adjJt_loc[3] / detJ;
    adjJt[i + Q * 4] = adjJt_loc[4] / detJ;
    adjJt[i + Q * 5] = adjJt_loc[5] / detJ;
    adjJt[i + Q * 6] = adjJt_loc[6] / detJ;
    adjJt[i + Q * 7] = adjJt_loc[7] / detJ;
    adjJt[i + Q * 8] = adjJt_loc[8] / detJ;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_GEOM_33_QF_H
