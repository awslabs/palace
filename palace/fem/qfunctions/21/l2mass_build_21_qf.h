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
  const CeedScalar *qdata, *qw = in[1];
  CeedScalar *qd1 = out[0], *qd2 = out[0] + Q;

  const CeedInt stride = 2 + 2; // [attr, w * |J|, adjJt_11, adjJt_12]_i, lumped by quad point
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    {
      const CeedScalar* adjJt_loc = qdata_i + 2;
      CeedScalar coeff[4], J_loc[2], qd_loc[1];
      CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0], coeff);
      AdjJt21(adjJt_loc, J_loc);
      MultAtBA21(J_loc, coeff, qd_loc);

      qd1[i + Q * 0] = qdata_i[1] * qd_loc[0];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond<2>((const CeedIntScalar *)ctx), (CeedInt)qdata_i[0]);

      qd2[i] = coeff * qw[i] * qw[i] / qdata_i[1];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_MASS_BUILD_21_QF_H
