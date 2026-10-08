// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_GEOM_21_QF_H
#define PALACE_LIBCEED_GEOM_21_QF_H

#include "utils_21_qf.h"

CEED_QFUNCTION(f_build_geom_factor_21)(void *, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *qw = in[1], *J = in[2];
  CeedScalar *qd_attr = out[0], *qd_wdetJ = out[0] + Q, *qd_adjJt = out[0] + 2 * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_loc[2], adjJt_loc[2];
    MatUnpack21(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt21<true>(J_loc, adjJt_loc);

    qd_attr[i] = attr[i];
    qd_wdetJ[i] = qw[i] * detJ;
    qd_adjJt[i + Q * 0] = adjJt_loc[0] / detJ;
    qd_adjJt[i + Q * 1] = adjJt_loc[1] / detJ;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_GEOM_21_QF_H
