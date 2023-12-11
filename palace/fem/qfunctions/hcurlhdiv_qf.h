// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_HDIV_QF_H
#define PALACE_LIBCEED_HCURL_HDIV_QF_H

#include "types_qf.h"
#include "utils_geom_qf.h"
#include "utils_qf.h"

// libCEED QFunctions for mixed H(curl)-H(div) operators (Piola transformations u =
// adj(J)^T / det(J) ̂u and u = J / det(J) ̂u).
// Note: J / det(J) = adj(adj(J)^T / det(J))^T
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is transpose adjugate Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is element attribute, shape [Q]
// in[3] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]

CEED_QFUNCTION(f_apply_hcurlhdiv_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar coeff[3], adjJt_loc[4], J_loc[4], v_loc[2];
    CoeffUnpack((const MatCoeffContext2 *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    AdjJt22<false>(adjJt_loc, J_loc);
    MultAtBCx22(J_loc, coeff, adjJt_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlhdiv_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar coeff[6], adjJt_loc[9], J_loc[9], v_loc[3];
    CoeffUnpack((const MatCoeffContext3 *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    AdjJt33<false>(adjJt_loc, J_loc);
    MultAtBCx33(J_loc, coeff, adjJt_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlhdiv_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[1] = {u[i + Q * 0]};
    CeedScalar coeff[3], adjJt_loc[2], J_loc[2], v_loc[2];
    CoeffUnpack((const MatCoeffContext2 *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack21(adjJt + i, Q, adjJt_loc);
    AdjJt21<false>(adjJt_loc, J_loc);
    MultAtBCx21(J_loc, coeff, adjJt_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlhdiv_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar coeff[6], adjJt_loc[6], J_loc[6], v_loc[3];
    CoeffUnpack((const MatCoeffContext3 *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack32(adjJt + i, Q, adjJt_loc);
    AdjJt32<false>(adjJt_loc, J_loc);
    MultAtBCx32(J_loc, coeff, adjJt_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivhcurl_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar coeff[3], adjJt_loc[4], J_loc[4], v_loc[2];
    CoeffUnpack((const MatCoeffContext2 *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    AdjJt22<false>(adjJt_loc, J_loc);
    MultAtBCx22(adjJt_loc, coeff, J_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivhcurl_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar coeff[6], adjJt_loc[9], J_loc[9], v_loc[3];
    CoeffUnpack((const MatCoeffContext3 *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    AdjJt33<false>(adjJt_loc, J_loc);
    MultAtBCx33(adjJt_loc, coeff, J_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivhcurl_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[1] = {u[i + Q * 0]};
    CeedScalar coeff[3], adjJt_loc[2], J_loc[2], v_loc[2];
    CoeffUnpack((const MatCoeffContext2 *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack21(adjJt + i, Q, adjJt_loc);
    AdjJt21<false>(adjJt_loc, J_loc);
    MultAtBCx21(adjJt_loc, coeff, J_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivhcurl_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    CeedScalar coeff[6], adjJt_loc[6], J_loc[6], v_loc[3];
    CoeffUnpack((const MatCoeffContext3 *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack32(adjJt + i, Q, adjJt_loc);
    AdjJt32<false>(adjJt_loc, J_loc);
    MultAtBCx32(adjJt_loc, coeff, J_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_HDIV_QF_H
