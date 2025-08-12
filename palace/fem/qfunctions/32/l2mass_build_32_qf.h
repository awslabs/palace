// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_MASS_BUILD_32_QF_H
#define PALACE_LIBCEED_L2_MASS_BUILD_32_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_3_qf.h"
#include "utils_32_qf.h"

CEED_QFUNCTION(f_build_l2mass_32)(void *__restrict__ ctx, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1];
  CeedScalar *qd1 = out[0], *qd2 = out[0] + 4 * Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      CeedScalar adjJt_loc[6], J_loc[6], qd_loc[4];
      const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);
      MatUnpack32(adjJt + i, Q, adjJt_loc);
      AdjJt32(adjJt_loc, J_loc);
      MultAtBA32(J_loc, coeff, qd_loc);

      qd1[i + Q * 0] = wdetJ[i] * qd_loc[0];
      qd1[i + Q * 1] = wdetJ[i] * qd_loc[1];
      qd1[i + Q * 2] = wdetJ[i] * qd_loc[2];
      qd1[i + Q * 3] = wdetJ[i] * qd_loc[3];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond<1>((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      qd2[i] = coeff * qw[i] * qw[i] / wdetJ[i];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_MASS_BUILD_32_QF_H
