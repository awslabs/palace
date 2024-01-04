// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_MASS_BUILD_QF_H
#define PALACE_LIBCEED_L2_MASS_BUILD_QF_H

#include "coeff_qf.h"
#include "utils_geom_qf.h"
#include "utils_qf.h"

// Build functions replace active vector output with quadrature point data, stored as a
// symmetric matrix, and remove active vector input.

CEED_QFUNCTION(f_build_l2mass_22)(void *__restrict__ ctx, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1];
  CeedScalar *qd1 = out[0], *qd2 = out[0] + 3 * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      CeedScalar coeff[3], adjJt_loc[4], J_loc[4], qd_loc[3];
      CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack22(adjJt + i, Q, adjJt_loc);
      AdjJt22<false>(adjJt_loc, J_loc);
      MultAtBA22(J_loc, coeff, qd_loc);

      qd1[i + Q * 0] = wdetJ[i] * qd_loc[0];
      qd1[i + Q * 1] = wdetJ[i] * qd_loc[1];
      qd1[i + Q * 2] = wdetJ[i] * qd_loc[2];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond2((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      qd2[i] = coeff * qw[i] * qw[i] / wdetJ[i];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_l2mass_33)(void *__restrict__ ctx, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1];
  CeedScalar *qd1 = out[0], *qd2 = out[0] + 6 * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      CeedScalar coeff[6], adjJt_loc[9], J_loc[9], qd_loc[6];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack33(adjJt + i, Q, adjJt_loc);
      AdjJt33<false>(adjJt_loc, J_loc);
      MultAtBA33(J_loc, coeff, qd_loc);

      qd1[i + Q * 0] = wdetJ[i] * qd_loc[0];
      qd1[i + Q * 1] = wdetJ[i] * qd_loc[1];
      qd1[i + Q * 2] = wdetJ[i] * qd_loc[2];
      qd1[i + Q * 3] = wdetJ[i] * qd_loc[3];
      qd1[i + Q * 4] = wdetJ[i] * qd_loc[4];
      qd1[i + Q * 5] = wdetJ[i] * qd_loc[5];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond3((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      qd2[i] = coeff * qw[i] * qw[i] / wdetJ[i];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_l2mass_21)(void *__restrict__ ctx, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1];
  CeedScalar *qd1 = out[0], *qd2 = out[0] + Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      CeedScalar coeff[3], adjJt_loc[2], J_loc[2], qd_loc[1];
      CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack21(adjJt + i, Q, adjJt_loc);
      AdjJt21<false>(adjJt_loc, J_loc);
      MultAtBA21(J_loc, coeff, qd_loc);

      qd1[i + Q * 0] = wdetJ[i] * qd_loc[0];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond2((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      qd2[i] = coeff * qw[i] * qw[i] / wdetJ[i];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_l2mass_32)(void *__restrict__ ctx, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1];
  CeedScalar *qd1 = out[0], *qd2 = out[0] + 3 * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      CeedScalar coeff[6], adjJt_loc[6], J_loc[6], qd_loc[3];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack32(adjJt + i, Q, adjJt_loc);
      AdjJt32<false>(adjJt_loc, J_loc);
      MultAtBA32(J_loc, coeff, qd_loc);

      qd1[i + Q * 0] = wdetJ[i] * qd_loc[0];
      qd1[i + Q * 1] = wdetJ[i] * qd_loc[1];
      qd1[i + Q * 2] = wdetJ[i] * qd_loc[2];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond3((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      qd2[i] = coeff * qw[i] * qw[i] / wdetJ[i];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_MASS_BUILD_QF_H
