// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_GEOM_AXISYMMETRIC_21_QF_H
#define PALACE_LIBCEED_GEOM_AXISYMMETRIC_21_QF_H

#include "utils_21_qf.h"

// Axisymmetric (r, z) variant of f_build_geom_factor_21 (boundary curves of a
// two-dimensional mesh): the quadrature weight carries the revolution measure 2 pi x.
// in[3] is the mesh coordinates at quadrature points, shape [ncomp=space_dim, Q]

CEED_QFUNCTION(f_build_geom_factor_axisymmetric_21)(void *, CeedInt Q,
                                                    const CeedScalar *const *in,
                                                    CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *qw = in[1], *J = in[2], *x = in[3];
  CeedScalar *qd_attr = out[0], *qd_wdetJ = out[0] + Q, *qd_adjJt = out[0] + 2 * Q;
  const CeedScalar two_pi = 6.28318530717958647692528676655900577;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_loc[2], adjJt_loc[2];
    MatUnpack21(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt21<true>(J_loc, adjJt_loc);

    qd_attr[i] = attr[i];
    qd_wdetJ[i] = qw[i] * detJ * two_pi * x[i];
    qd_adjJt[i + Q * 0] = adjJt_loc[0] / detJ;
    qd_adjJt[i + Q * 1] = adjJt_loc[1] / detJ;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_GEOM_AXISYMMETRIC_21_QF_H
