// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_GEOM_33_QF_H
#define PALACE_LIBCEED_GEOM_33_QF_H

#include "utils_33_qf.h"

// Geometry factor quadrature data {attr, w |J|, adj(J)ᵀ / |J|}, optionally followed by the
// physical coordinates x of the quadrature points (used by QFunctions with spatially
// varying coefficients, such as the PML material tensors).
template <bool COORDS>
CEED_QFUNCTION_HELPER int
f_build_geom_factor_33_impl(CeedInt Q, const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *qw = in[1], *J = in[2], *x = in[COORDS ? 3 : 2];
  CeedScalar *qd_attr = out[0], *qd_wdetJ = out[0] + Q, *qd_adjJt = out[0] + 2 * Q,
             *qd_x = out[0] + 11 * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_loc[9], adjJt_loc[9];
    MatUnpack33(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt33<true>(J_loc, adjJt_loc);

    qd_attr[i] = attr[i];
    qd_wdetJ[i] = qw[i] * detJ;
    qd_adjJt[i + Q * 0] = adjJt_loc[0] / detJ;
    qd_adjJt[i + Q * 1] = adjJt_loc[1] / detJ;
    qd_adjJt[i + Q * 2] = adjJt_loc[2] / detJ;
    qd_adjJt[i + Q * 3] = adjJt_loc[3] / detJ;
    qd_adjJt[i + Q * 4] = adjJt_loc[4] / detJ;
    qd_adjJt[i + Q * 5] = adjJt_loc[5] / detJ;
    qd_adjJt[i + Q * 6] = adjJt_loc[6] / detJ;
    qd_adjJt[i + Q * 7] = adjJt_loc[7] / detJ;
    qd_adjJt[i + Q * 8] = adjJt_loc[8] / detJ;
    if (COORDS)
    {
      qd_x[i + Q * 0] = x[i + Q * 0];
      qd_x[i + Q * 1] = x[i + Q * 1];
      qd_x[i + Q * 2] = x[i + Q * 2];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_geom_factor_33)(void *, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  return f_build_geom_factor_33_impl<false>(Q, in, out);
}

CEED_QFUNCTION(f_build_geom_factor_coords_33)(void *, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  return f_build_geom_factor_33_impl<true>(Q, in, out);
}

#endif  // PALACE_LIBCEED_GEOM_33_QF_H
