// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_MASS_31_QF_H
#define PALACE_LIBCEED_L2_MASS_31_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_3_qf.h"
#include "utils_31_qf.h"

CEED_QFUNCTION(f_apply_l2mass_31)(void *__restrict__ ctx, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1],
                   *u = in[2], *divu = in[3];
  CeedScalar *__restrict__ v = out[0], *__restrict__ divv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar u_loc[1] = {u[i + Q * 0]};
      CeedScalar coeff[9], adjJt_loc[3], J_loc[3], v_loc[1];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack31(adjJt + i, Q, adjJt_loc);
      AdjJt31(adjJt_loc, J_loc);
      MultAtBCx31(J_loc, coeff, J_loc, u_loc, v_loc);

      v[i + Q * 0] = wdetJ[i] * v_loc[0];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond<3>((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_MASS_31_QF_H
