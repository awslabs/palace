// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_MASS_21_QF_H
#define PALACE_LIBCEED_L2_MASS_21_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_2_qf.h"
#include "utils_21_qf.h"

CEED_QFUNCTION(f_apply_l2mass_21)(void *__restrict__ ctx, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *qdata = in[0], *qw = in[1], *u = in[2], *divu = in[3];
  CeedScalar *__restrict__ v = out[0], *__restrict__ divv = out[1];

  const CeedInt stride = 2 + 2; // [attr, w * |J|, adjJt_11, adjJt_12]_i, lumped by quad point
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    {
      const CeedScalar* adjJt_loc = qdata_i + 2;
      const CeedScalar u_loc[1] = {u[i + Q * 0]};
      CeedScalar coeff[4], J_loc[2], v_loc[1];
      CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0], coeff);
      AdjJt21(adjJt_loc, J_loc);
      MultAtBCx21(J_loc, coeff, J_loc, u_loc, v_loc);

      v[i + Q * 0] = qdata_i[1] * v_loc[0];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond<2>((const CeedIntScalar *)ctx), (CeedInt)qdata_i[0]);

      divv[i] = (coeff * qw[i] * qw[i] / qdata_i[1]) * divu[i];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_MASS_21_QF_H
