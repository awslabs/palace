// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_H1D_QF_H
#define PALACE_LIBCEED_HCURL_H1D_QF_H

#include "types_qf.h"
#include "utils_qf.h"

// libCEED QFunctions for mixed H(curl)-(H1)ᵈ operators (Piola transformation u =
// adj(J)^T / det(J) ̂u and u = ̂u)
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is transpose adjugate Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[3] is element attribute, shape [1]
// out[0] is active vector, shape [ncomp=space_dim, Q]

CEED_QFUNCTION(f_apply_hcurlh1d_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                    CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2];
  CeedScalar *v = out[0];

  MatCoeffContext2 *bc = (MatCoeffContext2 *)ctx;
  const CeedInt attr = (CeedInt)*in[3];
  const CeedScalar *coeff = bc->mat_coeff[bc->attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[4];
    MultBA22(adjJt + i, Q, coeff, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[2] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[3] * u1);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlh1d_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                    CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2];
  CeedScalar *v = out[0];

  MatCoeffContext3 *bc = (MatCoeffContext3 *)ctx;
  const CeedInt attr = (CeedInt)*in[3];
  const CeedScalar *coeff = bc->mat_coeff[bc->attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[9];
    MultBA33(adjJt + i, Q, coeff, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    const CeedScalar u2 = u[i + Q * 2];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[3] * u1 + qd[6] * u2);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[4] * u1 + qd[7] * u2);
    v[i + Q * 2] = wdetJ[i] * (qd[2] * u0 + qd[5] * u1 + qd[8] * u2);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlh1d_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                    CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2];
  CeedScalar *v = out[0];

  MatCoeffContext2 *bc = (MatCoeffContext2 *)ctx;
  const CeedInt attr = (CeedInt)*in[3];
  const CeedScalar *coeff = bc->mat_coeff[bc->attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[2];
    MultBA21(adjJt + i, Q, coeff, qd);

    const CeedScalar u0 = u[i + Q * 0];
    v[i + Q * 0] = wdetJ[i] * qd[0] * u0;
    v[i + Q * 1] = wdetJ[i] * qd[1] * u0;
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlh1d_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                    CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2];
  CeedScalar *v = out[0];

  MatCoeffContext3 *bc = (MatCoeffContext3 *)ctx;
  const CeedInt attr = (CeedInt)*in[3];
  const CeedScalar *coeff = bc->mat_coeff[bc->attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[6];
    MultBA33(adjJt + i, Q, coeff, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[3] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[4] * u1);
    v[i + Q * 2] = wdetJ[i] * (qd[2] * u0 + qd[5] * u1);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_H1D_QF_H
