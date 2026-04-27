// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_GEOM_32_QF_H
#define PALACE_LIBCEED_GEOM_32_QF_H

#include "utils_32_qf.h"

// Quadrature data layout: {attr, w|J|, adj(J)^T/|J|, x} with x the physical
// coordinate at the QP (needed for spatially-varying coefficients like PML).
// For space_dim=3, dim=2 elements (surfaces embedded in 3D), x has 3 components.
CEED_QFUNCTION(f_build_geom_factor_32)(void *, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *qw = in[1], *J = in[2], *x = in[3];
  CeedScalar *qd_attr = out[0], *qd_wdetJ = out[0] + Q, *qd_adjJt = out[0] + 2 * Q,
             *qd_x = out[0] + (2 + 6) * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_loc[6], adjJt_loc[6];
    MatUnpack32(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt32<true>(J_loc, adjJt_loc);

    qd_attr[i] = attr[i];
    qd_wdetJ[i] = qw[i] * detJ;
    qd_adjJt[i + Q * 0] = adjJt_loc[0] / detJ;
    qd_adjJt[i + Q * 1] = adjJt_loc[1] / detJ;
    qd_adjJt[i + Q * 2] = adjJt_loc[2] / detJ;
    qd_adjJt[i + Q * 3] = adjJt_loc[3] / detJ;
    qd_adjJt[i + Q * 4] = adjJt_loc[4] / detJ;
    qd_adjJt[i + Q * 5] = adjJt_loc[5] / detJ;
    qd_x[i + Q * 0] = x[i + Q * 0];
    qd_x[i + Q * 1] = x[i + Q * 1];
    qd_x[i + Q * 2] = x[i + Q * 2];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_GEOM_32_QF_H
