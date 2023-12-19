// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_MASS_QF_H
#define PALACE_LIBCEED_L2_MASS_QF_H

#include "coeff_qf.h"
#include "utils_geom_qf.h"
#include "utils_qf.h"

// libCEED QFunctions for L2 + H(div) mass operators (Piola transformations u = 1 / det(J) ̂u
// and u = J / det(J) ̂u).
// Note: J / det(J) = adj(adj(J)^T / det(J))^T
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is quadrature weights, shape [Q]
// in[2] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[3] is active vector divergence, shape [ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[1] is active vector divergence, shape [ncomp=1, Q]

CEED_QFUNCTION(f_apply_l2mass_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1],
                   *u = in[2], *divu = in[3];
  CeedScalar *__restrict__ v = out[0], *__restrict__ divv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
      CeedScalar coeff[3], adjJt_loc[4], J_loc[4], v_loc[2];
      CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack22(adjJt + i, Q, adjJt_loc);
      AdjJt22<false>(adjJt_loc, J_loc);
      MultAtBCx22(J_loc, coeff, J_loc, u_loc, v_loc);

      v[i + Q * 0] = wdetJ[i] * v_loc[0];
      v[i + Q * 1] = wdetJ[i] * v_loc[1];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond2((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2mass_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1],
                   *u = in[2], *divu = in[3];
  CeedScalar *__restrict__ v = out[0], *__restrict__ divv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
      CeedScalar coeff[6], adjJt_loc[9], J_loc[9], v_loc[3];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack33(adjJt + i, Q, adjJt_loc);
      AdjJt33<false>(adjJt_loc, J_loc);
      MultAtBCx33(J_loc, coeff, J_loc, u_loc, v_loc);

      v[i + Q * 0] = wdetJ[i] * v_loc[0];
      v[i + Q * 1] = wdetJ[i] * v_loc[1];
      v[i + Q * 2] = wdetJ[i] * v_loc[2];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond3((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2mass_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1],
                   *u = in[2], *divu = in[3];
  CeedScalar *__restrict__ v = out[0], *__restrict__ divv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar u_loc[1] = {u[i + Q * 0]};
      CeedScalar coeff[3], adjJt_loc[2], J_loc[2], v_loc[2];
      CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack21(adjJt + i, Q, adjJt_loc);
      AdjJt21<false>(adjJt_loc, J_loc);
      MultAtBCx21(J_loc, coeff, J_loc, u_loc, v_loc);

      v[i + Q * 0] = wdetJ[i] * v_loc[0];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond2((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2mass_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1],
                   *u = in[2], *divu = in[3];
  CeedScalar *__restrict__ v = out[0], *__restrict__ divv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
      CeedScalar coeff[6], adjJt_loc[6], J_loc[6], v_loc[3];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack32(adjJt + i, Q, adjJt_loc);
      AdjJt32<false>(adjJt_loc, J_loc);
      MultAtBCx32(J_loc, coeff, J_loc, u_loc, v_loc);

      v[i + Q * 0] = wdetJ[i] * v_loc[0];
      v[i + Q * 1] = wdetJ[i] * v_loc[1];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond3((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_MASS_QF_H
