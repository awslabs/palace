// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURLHDIV_COMPLEX_ERROR_33_QF_H
#define PALACE_LIBCEED_HCURLHDIV_COMPLEX_ERROR_33_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

// Complex constitutive flux error (D_h - D)ᴴ W⁻¹ (D_h - D) with D_h = (ε' + i ε'') E for a
// complex H(curl) field E and its recovered H(div) flux D. The context holds three
// coefficient tables: ε', ε'', and W⁻¹. in[1..2] are Re{E}, Im{E} and in[3..4] are Re{D},
// Im{D} (reference values).
CEED_QFUNCTION(f_apply_hcurlhdiv_complex_error_33)(void *__restrict__ ctx, CeedInt Q,
                                                   const CeedScalar *const *in,
                                                   CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u1r = in[1],
                   *u1i = in[2], *u2r = in[3], *u2i = in[4];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedIntScalar *ctx_r = (const CeedIntScalar *)ctx;
    const CeedIntScalar *ctx_i = CoeffPairSecond<3>(ctx_r);
    const CeedIntScalar *ctx_w = CoeffPairSecond<3>(ctx_i);
    CeedScalar adjJt_loc[9], J_loc[9], eps_r[9], eps_i[9], winv[9];
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    AdjJt33(adjJt_loc, J_loc);
    CoeffUnpack3(ctx_r, (CeedInt)attr[i], eps_r);
    CoeffUnpack3(ctx_i, (CeedInt)attr[i], eps_i);
    CoeffUnpack3(ctx_w, (CeedInt)attr[i], winv);

    // Physical E (covariant Piola) and D (contravariant Piola).
    CeedScalar Er[3], Ei[3], Dr[3], Di[3];
    {
      const CeedScalar ur[3] = {u1r[i + Q * 0], u1r[i + Q * 1], u1r[i + Q * 2]};
      const CeedScalar ui[3] = {u1i[i + Q * 0], u1i[i + Q * 1], u1i[i + Q * 2]};
      MultBx33(adjJt_loc, ur, Er);
      MultBx33(adjJt_loc, ui, Ei);
    }
    {
      const CeedScalar ur[3] = {u2r[i + Q * 0], u2r[i + Q * 1], u2r[i + Q * 2]};
      const CeedScalar ui[3] = {u2i[i + Q * 0], u2i[i + Q * 1], u2i[i + Q * 2]};
      MultBx33(J_loc, ur, Dr);
      MultBx33(J_loc, ui, Di);
    }

    // e = (ε' + i ε'') E - D.
    CeedScalar a[3], b[3], er[3], ei[3];
    MultBx33(eps_r, Er, a);
    MultBx33(eps_i, Ei, b);
    for (CeedInt d = 0; d < 3; d++)
    {
      er[d] = a[d] - b[d] - Dr[d];
    }
    MultBx33(eps_r, Ei, a);
    MultBx33(eps_i, Er, b);
    for (CeedInt d = 0; d < 3; d++)
    {
      ei[d] = a[d] + b[d] - Di[d];
    }
    MultBx33(winv, er, a);
    MultBx33(winv, ei, b);
    v[i] = wdetJ[i] * (er[0] * a[0] + er[1] * a[1] + er[2] * a[2] + ei[0] * b[0] +
                       ei[1] * b[1] + ei[2] * b[2]);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURLHDIV_COMPLEX_ERROR_33_QF_H
