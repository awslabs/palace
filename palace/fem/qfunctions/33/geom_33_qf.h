// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_GEOM_33_QF_H
#define PALACE_LIBCEED_GEOM_33_QF_H

#include "utils_33_qf.h"

CEED_QFUNCTION(f_build_geom_factor_33)(void *, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *qw = in[1], *J = in[2];
  // CeedScalar *qd_attr = out[0], *qd_wdetJ = out[0] + Q, *qd_adjJt = out[0] + 2 * Q;
  CeedScalar *qdata = out[0];

  const CeedInt stride = 2 + 9; // attr, w * |J|, (adjJt / |J|) colwise
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_loc[9], adjJt_loc[9];
    MatUnpack33(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt33<true>(J_loc, adjJt_loc);

    // qd_attr[i] = attr[i];
    // qd_wdetJ[i] = qw[i] * detJ;
    // qd_adjJt[i + Q * 0] = adjJt_loc[0] / detJ;
    // qd_adjJt[i + Q * 1] = adjJt_loc[1] / detJ;
    // qd_adjJt[i + Q * 2] = adjJt_loc[2] / detJ;
    // qd_adjJt[i + Q * 3] = adjJt_loc[3] / detJ;
    // qd_adjJt[i + Q * 4] = adjJt_loc[4] / detJ;
    // qd_adjJt[i + Q * 5] = adjJt_loc[5] / detJ;
    // qd_adjJt[i + Q * 6] = adjJt_loc[6] / detJ;
    // qd_adjJt[i + Q * 7] = adjJt_loc[7] / detJ;
    // qd_adjJt[i + Q * 8] = adjJt_loc[8] / detJ;

    qdata[i * stride +  0] = attr[i]; // attribute
    qdata[i * stride +  1] = qw[i] * detJ; // w * |J|
    qdata[i * stride +  2] = adjJt_loc[0] / detJ; // (adjJt / |J|)_11
    qdata[i * stride +  3] = adjJt_loc[1] / detJ; // (adjJt / |J|)_12
    qdata[i * stride +  4] = adjJt_loc[2] / detJ; // (adjJt / |J|)_13
    qdata[i * stride +  5] = adjJt_loc[3] / detJ; // (adjJt / |J|)_21
    qdata[i * stride +  6] = adjJt_loc[4] / detJ; // (adjJt / |J|)_22
    qdata[i * stride +  7] = adjJt_loc[5] / detJ; // (adjJt / |J|)_23
    qdata[i * stride +  8] = adjJt_loc[6] / detJ; // (adjJt / |J|)_31
    qdata[i * stride +  9] = adjJt_loc[7] / detJ; // (adjJt / |J|)_32
    qdata[i * stride + 10] = adjJt_loc[8] / detJ; // (adjJt / |J|)_33
  }
  return 0;
}

#endif  // PALACE_LIBCEED_GEOM_33_QF_H
