// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_MASS_BUILD_21_QF_H
#define PALACE_LIBCEED_L2_MASS_BUILD_21_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_2_qf.h"
#include "utils_21_qf.h"

CEED_QFUNCTION(f_build_l2mass_21)(void *__restrict__ ctx, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1];
  CeedScalar *qd1 = out[0], *qd2 = out[0] + Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      CeedScalar coeff[4], adjJt_loc[2], J_loc[2], qd_loc[1];
      CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack21(adjJt + i, Q, adjJt_loc);
      AdjJt21(adjJt_loc, J_loc);
      MultAtBA21(J_loc, coeff, qd_loc);

      qd1[i + Q * 0] = wdetJ[i] * qd_loc[0];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond<2>((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      qd2[i] = coeff * qw[i] * qw[i] / wdetJ[i];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_MASS_BUILD_21_QF_H
