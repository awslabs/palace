// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HDIV_X_33_QF_H
#define PALACE_LIBCEED_HDIV_X_33_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

// The H(div) operator of f_apply_hdiv_33 (u = J / det(J) û, v = w / det(J) Jᵀ C J û) with
// the geometry factors computed at the quadrature points from the mesh nodes.
// in[0] is the element attribute, shape [Q]
// in[1] is the quadrature weight, shape [Q]
// in[2] is the Jacobian of the mesh map (the gradient of the nodes), shape [9, Q]
// in[3] is active vector, shape [qcomp=3, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=3, ncomp=1, Q]
CEED_QFUNCTION(f_apply_hdiv_x_33)(void *__restrict__ ctx, CeedInt Q,
                                  const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *qw = in[1], *J = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar coeff[9], J_loc[9], v_loc[3];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack33(J + i, Q, J_loc);
    const CeedScalar w = qw[i] / DetJ33(J_loc);
    MultAtBCx33(J_loc, coeff, J_loc, u_loc, v_loc);

    v[i + Q * 0] = w * v_loc[0];
    v[i + Q * 1] = w * v_loc[1];
    v[i + Q * 2] = w * v_loc[2];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HDIV_X_33_QF_H
