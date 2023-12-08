// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HDIV_MASS_QF_H
#define PALACE_LIBCEED_HDIV_MASS_QF_H

#include "types_qf.h"
#include "utils_geom_qf.h"
#include "utils_qf.h"

// libCEED QFunctions for H(div) + H(curl) mass operators in 3D (Piola transformations u =
// J / det(J) ̂u and u = adj(J)^T / det(J) ̂u).
// Note: J / det(J) = adj(adj(J)^T / det(J))^T
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is transpose adjugate Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is element attribute, shape [Q]
// in[3] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[4] is active vector curl, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[1] is active vector curl, shape [qcomp=dim, ncomp=1, Q]

// In 2D, this actually uses the L2 Piola transformation on the curl (u = 1 / det(J) ̂u) and
// the curl is has qcomp=1.
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is transpose adjugate Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is quadrature weights, shape [Q]
// in[3] is element attribute, shape [Q]
// in[4] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[5] is active vector curl, shape [ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[1] is active vector curl, shape [ncomp=1, Q]

CEED_QFUNCTION(f_apply_hdivmass_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                    CeedScalar *const *out)
{
  const MatCoeffPairContext21 *bc = (const MatCoeffPairContext21 *)ctx;
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *qw = in[2], *attr = in[3], *u = in[4],
                   *curlu = in[5];
  CeedScalar *v = out[0], *curlv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar coeff[3], adjJt_loc[4], v_loc[2];
    CoeffUnpack(&bc->first, (CeedInt)attr[i], coeff);
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    MultAtBCx22(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack(&bc->second, (CeedInt)attr[i]);

    curlv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * curlu[i];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivmass_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                    CeedScalar *const *out)
{
  const MatCoeffPairContext33 *bc = (const MatCoeffPairContext33 *)ctx;
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3],
                   *curlu = in[4];
  CeedScalar *v = out[0], *curlv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar coeff[6], adjJt_loc[9], v_loc[3];
    CoeffUnpack(&bc->first, (CeedInt)attr[i], coeff);
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    MultAtBCx33(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {curlu[i + Q * 0], curlu[i + Q * 1], curlu[i + Q * 2]};
    CeedScalar coeff[6], adjJt_loc[9], J_loc[9], v_loc[3];
    CoeffUnpack(&bc->second, (CeedInt)attr[i], coeff);
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    AdjJt33<false>(adjJt_loc, J_loc);
    MultAtBCx33(J_loc, coeff, J_loc, u_loc, v_loc);

    curlv[i + Q * 0] = wdetJ[i] * v_loc[0];
    curlv[i + Q * 1] = wdetJ[i] * v_loc[1];
    curlv[i + Q * 2] = wdetJ[i] * v_loc[2];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivmass_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                    CeedScalar *const *out)
{
  const MatCoeffPairContext31 *bc = (const MatCoeffPairContext31 *)ctx;
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *qw = in[2], *attr = in[3], *u = in[4],
                   *curlu = in[5];
  CeedScalar *v = out[0], *curlv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar coeff[6], adjJt_loc[6], v_loc[3];
    CoeffUnpack(&bc->first, (CeedInt)attr[i], coeff);
    MatUnpack32(adjJt + i, Q, adjJt_loc);
    MultAtBCx32(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack(&bc->second, (CeedInt)attr[i]);

    curlv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * curlu[i];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_CURLCURL_MASS_QF_H
