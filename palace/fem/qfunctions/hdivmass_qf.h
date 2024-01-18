// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HDIV_MASS_QF_H
#define PALACE_LIBCEED_HDIV_MASS_QF_H

#include "coeff_qf.h"
#include "utils_geom_qf.h"
#include "utils_qf.h"

// libCEED QFunctions for H(div) + H(curl) mass operators in 3D (Piola transformations u =
// J / det(J) ̂u and u = adj(J)^T / det(J) ̂u).
// Note: J / det(J) = adj(adj(J)^T / det(J))^T
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[2] is active vector curl, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[1] is active vector curl, shape [qcomp=dim, ncomp=1, Q]

// In 2D, this actually uses the L2 Piola transformation on the curl (u = 1 / det(J) ̂u) and
// the curl is has qcomp=1.
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is quadrature weights, shape [Q]
// in[2] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[3] is active vector curl, shape [ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[1] is active vector curl, shape [ncomp=1, Q]

CEED_QFUNCTION(f_apply_hdivmass_22)(void *__restrict__ ctx, CeedInt Q,
                                    const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1],
                   *u = in[2], *curlu = in[3];
  CeedScalar *__restrict__ v = out[0], *__restrict__ curlv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
      CeedScalar coeff[3], adjJt_loc[4], v_loc[2];
      CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack22(adjJt + i, Q, adjJt_loc);
      MultAtBCx22(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

      v[i + Q * 0] = wdetJ[i] * v_loc[0];
      v[i + Q * 1] = wdetJ[i] * v_loc[1];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond2((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      curlv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * curlu[i];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivmass_33)(void *__restrict__ ctx, CeedInt Q,
                                    const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *curlu = in[2];
  CeedScalar *__restrict__ v = out[0], *__restrict__ curlv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
      CeedScalar coeff[6], adjJt_loc[9], v_loc[3];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack33(adjJt + i, Q, adjJt_loc);
      MultAtBCx33(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

      v[i + Q * 0] = wdetJ[i] * v_loc[0];
      v[i + Q * 1] = wdetJ[i] * v_loc[1];
      v[i + Q * 2] = wdetJ[i] * v_loc[2];
    }
    {
      const CeedScalar u_loc[3] = {curlu[i + Q * 0], curlu[i + Q * 1], curlu[i + Q * 2]};
      CeedScalar coeff[6], adjJt_loc[9], J_loc[9], v_loc[3];
      CoeffUnpack3(CoeffPairSecond3((const CeedIntScalar *)ctx), (CeedInt)attr[i], coeff);
      MatUnpack33(adjJt + i, Q, adjJt_loc);
      AdjJt33<false>(adjJt_loc, J_loc);
      MultAtBCx33(J_loc, coeff, J_loc, u_loc, v_loc);

      curlv[i + Q * 0] = wdetJ[i] * v_loc[0];
      curlv[i + Q * 1] = wdetJ[i] * v_loc[1];
      curlv[i + Q * 2] = wdetJ[i] * v_loc[2];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivmass_32)(void *__restrict__ ctx, CeedInt Q,
                                    const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1],
                   *u = in[2], *curlu = in[3];
  CeedScalar *__restrict__ v = out[0], *__restrict__ curlv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
      CeedScalar coeff[6], adjJt_loc[6], v_loc[3];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack32(adjJt + i, Q, adjJt_loc);
      MultAtBCx32(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

      v[i + Q * 0] = wdetJ[i] * v_loc[0];
      v[i + Q * 1] = wdetJ[i] * v_loc[1];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond3((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      curlv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * curlu[i];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_CURLCURL_MASS_QF_H
