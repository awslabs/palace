// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURLHDIV_COMPLEX_ERROR_22_QF_H
#define PALACE_LIBCEED_HCURLHDIV_COMPLEX_ERROR_22_QF_H

#include "../coeff/coeff_2_qf.h"
#include "utils_22_qf.h"

// 2D version of f_apply_hcurlhdiv_complex_error_33.
CEED_QFUNCTION(f_apply_hcurlhdiv_complex_error_22)(void *__restrict__ ctx, CeedInt Q,
                                                   const CeedScalar *const *in,
                                                   CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u1r = in[1],
                   *u1i = in[2], *u2r = in[3], *u2i = in[4];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedIntScalar *ctx_r = (const CeedIntScalar *)ctx;
    const CeedIntScalar *ctx_i = CoeffPairSecond<2>(ctx_r);
    const CeedIntScalar *ctx_w = CoeffPairSecond<2>(ctx_i);
    CeedScalar adjJt_loc[4], J_loc[4], eps_r[4], eps_i[4], winv[4];
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    AdjJt22(adjJt_loc, J_loc);
    CoeffUnpack2(ctx_r, (CeedInt)attr[i], eps_r);
    CoeffUnpack2(ctx_i, (CeedInt)attr[i], eps_i);
    CoeffUnpack2(ctx_w, (CeedInt)attr[i], winv);

    // Physical E (covariant Piola) and D (contravariant Piola).
    CeedScalar Er[2], Ei[2], Dr[2], Di[2];
    {
      const CeedScalar ur[2] = {u1r[i + Q * 0], u1r[i + Q * 1]};
      const CeedScalar ui[2] = {u1i[i + Q * 0], u1i[i + Q * 1]};
      MultBx22(adjJt_loc, ur, Er);
      MultBx22(adjJt_loc, ui, Ei);
    }
    {
      const CeedScalar ur[2] = {u2r[i + Q * 0], u2r[i + Q * 1]};
      const CeedScalar ui[2] = {u2i[i + Q * 0], u2i[i + Q * 1]};
      MultBx22(J_loc, ur, Dr);
      MultBx22(J_loc, ui, Di);
    }

    // e = (ε' + i ε'') E - D.
    CeedScalar a[2], b[2], er[2], ei[2];
    MultBx22(eps_r, Er, a);
    MultBx22(eps_i, Ei, b);
    for (CeedInt d = 0; d < 2; d++)
    {
      er[d] = a[d] - b[d] - Dr[d];
    }
    MultBx22(eps_r, Ei, a);
    MultBx22(eps_i, Er, b);
    for (CeedInt d = 0; d < 2; d++)
    {
      ei[d] = a[d] + b[d] - Di[d];
    }
    MultBx22(winv, er, a);
    MultBx22(winv, ei, b);
    v[i] = wdetJ[i] * (er[0] * a[0] + er[1] * a[1] + ei[0] * b[0] + ei[1] * b[1]);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURLHDIV_COMPLEX_ERROR_22_QF_H
