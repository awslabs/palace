// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HDIV_21_QF_H
#define PALACE_LIBCEED_HDIV_21_QF_H

#include "../coeff/coeff_2_qf.h"
#include "utils_21_qf.h"

CEED_QFUNCTION(f_apply_hdiv_21)(void *__restrict__ ctx, CeedInt Q,
                                const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[1] = {u[i + Q * 0]};
    CeedScalar coeff[4], adjJt_loc[2], J_loc[2], v_loc[1];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack21(adjJt + i, Q, adjJt_loc);
    AdjJt21(adjJt_loc, J_loc);
    MultAtBCx21(J_loc, coeff, J_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HDIV_21_QF_H
