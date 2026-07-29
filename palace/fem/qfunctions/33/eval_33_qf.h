// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_EVAL_33_QF_H
#define PALACE_LIBCEED_EVAL_33_QF_H

#include "utils_33_qf.h"

// QFunctions for pointwise evaluation at arbitrary points of 3D volume elements.

CEED_QFUNCTION(f_eval_probe_hcurl_33)(void *, CeedInt Q, const CeedScalar *const *in,
                                      CeedScalar *const *out)
{
  const CeedScalar *J = in[0], *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar J_loc[9], adjJt_loc[9], E[3];
    MatUnpack33(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt33<true>(J_loc, adjJt_loc);
    MultAx33(adjJt_loc, u_loc, E);
    v[i + Q * 0] = E[0] / detJ;
    v[i + Q * 1] = E[1] / detJ;
    v[i + Q * 2] = E[2] / detJ;
  }
  return 0;
}

// Pointwise scalar H1/L2 field value. Inputs: grad_x, u. The geometry input is unused
// but kept so this probe has the same field list as the vector-valued probes.
CEED_QFUNCTION(f_eval_probe_l2_33)(void *, CeedInt Q, const CeedScalar *const *in,
                                   CeedScalar *const *out)
{
  const CeedScalar *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = u[i];
  }
  return 0;
}

// Pointwise H(div) field value: v = J/detJ u. Inputs: grad_x, u.
CEED_QFUNCTION(f_eval_probe_hdiv_33)(void *, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *J = in[0], *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar J_loc[9], adjJt_loc[9], B[3];
    MatUnpack33(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt33<true>(J_loc, adjJt_loc);
    MultAx33(J_loc, u_loc, B);
    v[i + Q * 0] = B[0] / detJ;
    v[i + Q * 1] = B[1] / detJ;
    v[i + Q * 2] = B[2] / detJ;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_EVAL_33_QF_H
