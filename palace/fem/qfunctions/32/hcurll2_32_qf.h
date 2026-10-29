// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_L2_32_QF_H
#define PALACE_LIBCEED_HCURL_L2_32_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_32_qf.h"

// The curl of an H(curl) function on a 2D manifold in 3D is the scalar 1 / det(J) curl ̂u,
// the component of the curl of its tangential trace along the element normal
// n = (J₁ x J₂) / det(J). It is used here as the vector (n · curl u) n = (a₁ x a₂) curl ̂u,
// where a₁, a₂ are the columns of adj(J)^T / det(J), so that a 3 x 3 matrix coefficient can
// act on it (a₁ x a₂ = n / det(J)).

CEED_QFUNCTION(f_apply_hcurll2_32)(void *__restrict__ ctx, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar coeff[9], adjJt_loc[6], n_loc[3], w_loc[3];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack32(adjJt + i, Q, adjJt_loc);
    CrossCols32(adjJt_loc, n_loc);
    MultBAx32(adjJt_loc, coeff, u_loc, w_loc);

    v[i] = wdetJ[i] * (n_loc[0] * w_loc[0] + n_loc[1] * w_loc[1] + n_loc[2] * w_loc[2]);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2hcurl_32)(void *__restrict__ ctx, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[9], adjJt_loc[6], n_loc[3], w_loc[3], v_loc[2];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack32(adjJt + i, Q, adjJt_loc);
    CrossCols32(adjJt_loc, n_loc);
    n_loc[0] *= u[i];
    n_loc[1] *= u[i];
    n_loc[2] *= u[i];
    MultBx32(coeff, n_loc, w_loc);
    MultAtx32(adjJt_loc, w_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_L2_32_QF_H
