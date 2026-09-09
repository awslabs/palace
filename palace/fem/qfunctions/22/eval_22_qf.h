// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_EVAL_22_QF_H
#define PALACE_LIBCEED_EVAL_22_QF_H

#include "utils_22_qf.h"

// Pointwise probes at arbitrary reference points of 2D volume elements. No quadrature
// weighting is applied.

// Pointwise H(curl) field value: v = adj(J)^T/detJ u. Inputs: grad_x, u.
CEED_QFUNCTION(f_eval_probe_hcurl_22)(void *, CeedInt Q, const CeedScalar *const *in,
                                      CeedScalar *const *out)
{
  const CeedScalar *J = in[0], *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar J_loc[4], adjJt_loc[4], E[2];
    MatUnpack22(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt22<true>(J_loc, adjJt_loc);
    MultBx22(adjJt_loc, u_loc, E);
    v[i + Q * 0] = E[0] / detJ;
    v[i + Q * 1] = E[1] / detJ;
  }
  return 0;
}

// Pointwise H(div) field value: v = J/detJ u. Inputs: grad_x, u.
CEED_QFUNCTION(f_eval_probe_hdiv_22)(void *, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *J = in[0], *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar J_loc[4], B[2];
    MatUnpack22(J + i, Q, J_loc);
    const CeedScalar detJ = DetJ22(J_loc);
    MultBx22(J_loc, u_loc, B);
    v[i + Q * 0] = B[0] / detJ;
    v[i + Q * 1] = B[1] / detJ;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_EVAL_22_QF_H
