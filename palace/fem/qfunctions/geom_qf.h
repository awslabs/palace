// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_GEOM_QF_H
#define PALACE_LIBCEED_GEOM_QF_H

#include "utils_geom_qf.h"
#include "utils_qf.h"

// libCEED QFunction for building geometry factors for integration and transformations.
// At every quadrature point, compute qw * det(J) and adj(J)^T / |J| and store the result.
// in[0] is quadrature weights, shape [Q]
// in[1] is Jacobians, shape [qcomp=dim, ncomp=space_dim, Q]
// out[0] is quadrature data, stored as {attribute, Jacobian determinant, (transpose)
//        adjugate Jacobian} quadrature data, shape [ncomp=2+space_dim*dim, Q]

CEED_QFUNCTION(f_build_geom_factor_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  const CeedScalar *qw = in[0], *J = in[1];
  CeedScalar *wdetJ = out[0] + Q, *adjJt = out[0] + 2 * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_loc[4], adjJt_loc[4];
    MatUnpack22(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt22<true>(J_loc, adjJt_loc);

    wdetJ[i] = qw[i] * detJ;
    adjJt[i + Q * 0] = adjJt_loc[0] / detJ;
    adjJt[i + Q * 1] = adjJt_loc[1] / detJ;
    adjJt[i + Q * 2] = adjJt_loc[2] / detJ;
    adjJt[i + Q * 3] = adjJt_loc[3] / detJ;
  }
  return 0;
}

CEED_QFUNCTION(f_build_geom_factor_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  const CeedScalar *qw = in[0], *J = in[1];
  CeedScalar *wdetJ = out[0] + Q, *adjJt = out[0] + 2 * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_loc[9], adjJt_loc[9];
    MatUnpack33(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt33<true>(J_loc, adjJt_loc);

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

CEED_QFUNCTION(f_build_geom_factor_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  const CeedScalar *qw = in[0], *J = in[1];
  CeedScalar *wdetJ = out[0] + Q, *adjJt = out[0] + 2 * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_loc[2], adjJt_loc[2];
    MatUnpack21(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt21<true>(J_loc, adjJt_loc);

    wdetJ[i] = qw[i] * detJ;
    adjJt[i + Q * 0] = adjJt_loc[0] / detJ;
    adjJt[i + Q * 1] = adjJt_loc[1] / detJ;
  }
  return 0;
}

CEED_QFUNCTION(f_build_geom_factor_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  const CeedScalar *qw = in[0], *J = in[1];
  CeedScalar *wdetJ = out[0] + Q, *adjJt = out[0] + 2 * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar J_loc[6], adjJt_loc[6];
    MatUnpack32(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt32<true>(J_loc, adjJt_loc);

    wdetJ[i] = qw[i] * detJ;
    adjJt[i + Q * 0] = adjJt_loc[0] / detJ;
    adjJt[i + Q * 1] = adjJt_loc[1] / detJ;
    adjJt[i + Q * 2] = adjJt_loc[2] / detJ;
    adjJt[i + Q * 3] = adjJt_loc[3] / detJ;
    adjJt[i + Q * 4] = adjJt_loc[4] / detJ;
    adjJt[i + Q * 5] = adjJt_loc[5] / detJ;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_GEOM_QF_H
