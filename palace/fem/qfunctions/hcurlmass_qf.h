// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_MASS_QF_H
#define PALACE_LIBCEED_HCURL_MASS_QF_H

#include "types_qf.h"
#include "utils_qf.h"

// libCEED QFunctions for H(curl) + H1 mass operators (Piola transformation u =
// adj(J)^T / det(J) ̂u and u = ̂u).
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is transpose adjugate Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is active vector, shape [ncomp=1, Q]
// in[3] is active vector gradient, shape [qcomp=dim, ncomp=1, Q]
// in[4] is element attribute, shape [1]
// out[0] is active vector, shape [ncomp=1, Q]
// out[1] is active vector gradient, shape [qcomp=dim, ncomp=1, Q]

CEED_QFUNCTION(f_apply_hcurlmass_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2], *gradu = in[3];
  CeedScalar *v = out[0], *gradv = out[1];

  const CeedScalar *attr = in[4];
  MatCoeffPairContext21 *bc = (MatCoeffPairContext21 *)ctx;
  // const CeedInt attr = (CeedInt)*in[4];
  // const CeedScalar *coeff = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  // const CeedScalar coeff_mass = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {

    // XXX TODO TESTING
    const CeedScalar coeff_mass = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[(CeedInt)attr[i]]];

    v[i] = coeff_mass * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[3];

    // XXX TODO TESTING
    const CeedScalar *coeff = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[(CeedInt)attr[i]]];

    MultAtBA22(adjJt + i, Q, coeff, qd);

    const CeedScalar gradu0 = gradu[i + Q * 0];
    const CeedScalar gradu1 = gradu[i + Q * 1];
    gradv[i + Q * 0] = wdetJ[i] * (qd[0] * gradu0 + qd[1] * gradu1);
    gradv[i + Q * 1] = wdetJ[i] * (qd[1] * gradu0 + qd[2] * gradu1);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlmass_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2], *gradu = in[3];
  CeedScalar *v = out[0], *gradv = out[1];

  const CeedScalar *attr = in[4];
  MatCoeffPairContext31 *bc = (MatCoeffPairContext31 *)ctx;
  // const CeedInt attr = (CeedInt)*in[4];
  // const CeedScalar *coeff = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  // const CeedScalar coeff_mass = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {

    // XXX TODO TESTING
    const CeedScalar coeff_mass = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[(CeedInt)attr[i]]];

    v[i] = coeff_mass * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[6];

    // XXX TODO TESTING
    const CeedScalar *coeff = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[(CeedInt)attr[i]]];

    MultAtBA33(adjJt + i, Q, coeff, qd);

    const CeedScalar gradu0 = gradu[i + Q * 0];
    const CeedScalar gradu1 = gradu[i + Q * 1];
    const CeedScalar gradu2 = gradu[i + Q * 2];
    gradv[i + Q * 0] = wdetJ[i] * (qd[0] * gradu0 + qd[1] * gradu1 + qd[2] * gradu2);
    gradv[i + Q * 1] = wdetJ[i] * (qd[1] * gradu0 + qd[3] * gradu1 + qd[4] * gradu2);
    gradv[i + Q * 2] = wdetJ[i] * (qd[2] * gradu0 + qd[4] * gradu1 + qd[5] * gradu2);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlmass_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2], *gradu = in[3];
  CeedScalar *v = out[0], *gradv = out[1];

  const CeedScalar *attr = in[4];
  MatCoeffPairContext21 *bc = (MatCoeffPairContext21 *)ctx;
  // const CeedInt attr = (CeedInt)*in[4];
  // const CeedScalar *coeff = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  // const CeedScalar coeff_mass = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {

    // XXX TODO TESTING
    const CeedScalar coeff_mass = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[(CeedInt)attr[i]]];

    v[i] = coeff_mass * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[1];

    // XXX TODO TESTING
    const CeedScalar *coeff = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[(CeedInt)attr[i]]];

    MultAtBA21(adjJt + i, Q, coeff, qd);

    const CeedScalar gradu0 = gradu[i + Q * 0];
    gradv[i + Q * 0] = wdetJ[i] * qd[0] * gradu0;
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlmass_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2], *gradu = in[3];
  CeedScalar *v = out[0], *gradv = out[1];

  const CeedScalar *attr = in[4];
  MatCoeffPairContext31 *bc = (MatCoeffPairContext31 *)ctx;
  // const CeedInt attr = (CeedInt)*in[4];
  // const CeedScalar *coeff = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  // const CeedScalar coeff_mass = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {

    // XXX TODO TESTING
    const CeedScalar coeff_mass = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[(CeedInt)attr[i]]];

    v[i] = coeff_mass * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[3];

    // XXX TODO TESTING
    const CeedScalar *coeff = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[(CeedInt)attr[i]]];

    MultAtBA32(adjJt + i, Q, coeff, qd);

    const CeedScalar gradu0 = gradu[i + Q * 0];
    const CeedScalar gradu1 = gradu[i + Q * 1];
    gradv[i + Q * 0] = wdetJ[i] * (qd[0] * gradu0 + qd[1] * gradu1);
    gradv[i + Q * 1] = wdetJ[i] * (qd[1] * gradu0 + qd[2] * gradu1);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_MASS_QF_H
