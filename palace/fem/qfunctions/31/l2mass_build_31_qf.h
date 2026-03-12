// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_MASS_BUILD_31_QF_H
#define PALACE_LIBCEED_L2_MASS_BUILD_31_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_3_qf.h"
#include "utils_31_qf.h"

CEED_QFUNCTION(f_build_l2mass_31)(void *__restrict__ ctx, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1];
  CeedScalar *qd1 = out[0], *qd2 = out[0] + Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      CeedScalar coeff[9], adjJt_loc[3], J_loc[3], qd_loc[1];
      CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
      MatUnpack31(adjJt + i, Q, adjJt_loc);
      AdjJt31(adjJt_loc, J_loc);
      MultAtBA31(J_loc, coeff, qd_loc);

      qd1[i + Q * 0] = wdetJ[i] * qd_loc[0];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond<3>((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      qd2[i] = coeff * qw[i] * qw[i] / wdetJ[i];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_MASS_BUILD_31_QF_H
