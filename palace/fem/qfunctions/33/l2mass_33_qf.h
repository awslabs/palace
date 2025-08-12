// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_MASS_33_QF_H
#define PALACE_LIBCEED_L2_MASS_33_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

CEED_QFUNCTION(f_apply_l2mass_33)(void *__restrict__ ctx, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1],
                   *u = in[2], *divu = in[3];
  CeedScalar *__restrict__ v = out[0], *__restrict__ divv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
      CeedScalar adjJt_loc[9], J_loc[9], v_loc[3];
      const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);
      MatUnpack33(adjJt + i, Q, adjJt_loc);
      AdjJt33(adjJt_loc, J_loc);
      MultAtBCx33(J_loc, coeff, J_loc, u_loc, v_loc);

      v[i + Q * 0] = wdetJ[i] * v_loc[0];
      v[i + Q * 1] = wdetJ[i] * v_loc[1];
      v[i + Q * 2] = wdetJ[i] * v_loc[2];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond<1>((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_MASS_33_QF_H
