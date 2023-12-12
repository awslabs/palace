// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_MASS_QF_H
#define PALACE_LIBCEED_HCURL_MASS_QF_H

#include "coeff_qf.h"
#include "utils_qf.h"

// libCEED QFunctions for H(curl) + H1 mass operators (Piola transformation u =
// adj(J)^T / det(J) ̂u and u = ̂u).
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is transpose adjugate Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is element attribute, shape [Q]
// in[3] is active vector, shape [ncomp=1, Q]
// in[4] is active vector gradient, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [ncomp=1, Q]
// out[1] is active vector gradient, shape [qcomp=dim, ncomp=1, Q]

CEED_QFUNCTION(f_apply_hcurlmass_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3],
                   *gradu = in[4];
  CeedScalar *v = out[0], *gradv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff =
        CoeffUnpack1(CoeffPairSecond2((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

    v[i] = coeff * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {gradu[i + Q * 0], gradu[i + Q * 1]};
    CeedScalar coeff[3], adjJt_loc[4], v_loc[2];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    MultAtBCx22(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

    gradv[i + Q * 0] = wdetJ[i] * v_loc[0];
    gradv[i + Q * 1] = wdetJ[i] * v_loc[1];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlmass_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3],
                   *gradu = in[4];
  CeedScalar *v = out[0], *gradv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff =
        CoeffUnpack1(CoeffPairSecond3((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

    v[i] = coeff * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {gradu[i + Q * 0], gradu[i + Q * 1], gradu[i + Q * 2]};
    CeedScalar coeff[6], adjJt_loc[9], v_loc[3];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    MultAtBCx33(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

    gradv[i + Q * 0] = wdetJ[i] * v_loc[0];
    gradv[i + Q * 1] = wdetJ[i] * v_loc[1];
    gradv[i + Q * 2] = wdetJ[i] * v_loc[2];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlmass_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3],
                   *gradu = in[4];
  CeedScalar *v = out[0], *gradv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff =
        CoeffUnpack1(CoeffPairSecond2((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

    v[i] = coeff * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[1] = {gradu[i + Q * 0]};
    CeedScalar coeff[3], adjJt_loc[2], v_loc[2];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack21(adjJt + i, Q, adjJt_loc);
    MultAtBCx21(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

    gradv[i + Q * 0] = wdetJ[i] * v_loc[0];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlmass_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *attr = in[2], *u = in[3],
                   *gradu = in[4];
  CeedScalar *v = out[0], *gradv = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff =
        CoeffUnpack1(CoeffPairSecond3((const CeedIntScalar *)ctx), (CeedInt)attr[i]);

    v[i] = coeff * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {gradu[i + Q * 0], gradu[i + Q * 1]};
    CeedScalar coeff[6], adjJt_loc[6], v_loc[3];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack32(adjJt + i, Q, adjJt_loc);
    MultAtBCx32(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

    gradv[i + Q * 0] = wdetJ[i] * v_loc[0];
    gradv[i + Q * 1] = wdetJ[i] * v_loc[1];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_MASS_QF_H
