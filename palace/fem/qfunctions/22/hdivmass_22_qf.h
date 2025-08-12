// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HDIV_MASS_22_QF_H
#define PALACE_LIBCEED_HDIV_MASS_22_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_2_qf.h"
#include "utils_22_qf.h"

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
      CeedScalar adjJt_loc[4], v_loc[2];
      const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);
      MatUnpack22(adjJt + i, Q, adjJt_loc);
      MultAtBCx22(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

      v[i + Q * 0] = wdetJ[i] * v_loc[0];
      v[i + Q * 1] = wdetJ[i] * v_loc[1];
    }
    {
      const CeedScalar coeff =
          CoeffUnpack1(CoeffPairSecond<1>((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

      curlv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * curlu[i];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HDIV_MASS_22_QF_H
