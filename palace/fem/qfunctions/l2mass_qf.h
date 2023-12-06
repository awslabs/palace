// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_MASS_QF_H
#define PALACE_LIBCEED_L2_MASS_QF_H

#include "types_qf.h"
#include "utils_qf.h"

// libCEED QFunctions for L2 + H(div) mass operators (Piola transformations u = 1 / det(J) ̂u
// and u = J / det(J) ̂u).
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is quadrature weights, shape [Q]
// in[3] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[4] is active vector divergence, shape [ncomp=1, Q]
// in[5] is element attribute, shape [1]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[1] is active vector divergence, shape [ncomp=1, Q]

CEED_QFUNCTION(f_apply_l2mass_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *qw = in[2], *u = in[3], *divu = in[4];
  CeedScalar *v = out[0], *divv = out[1];

  const CeedScalar *attr = in[5];
  MatCoeffPairContext21 *bc = (MatCoeffPairContext21 *)ctx;
  // const CeedInt attr = (CeedInt)*in[5];
  // const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  // const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[3];

    // XXX TODO TESTING
    const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[(CeedInt)attr[i]]];

    MultAtBA22(J + i, Q, coeff_mass, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[1] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[2] * u1);
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {

    // XXX TODO TESTING
    const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[(CeedInt)attr[i]]];

    divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2mass_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *qw = in[2], *u = in[3], *divu = in[4];
  CeedScalar *v = out[0], *divv = out[1];

  const CeedScalar *attr = in[5];
  MatCoeffPairContext31 *bc = (MatCoeffPairContext31 *)ctx;
  // const CeedInt attr = (CeedInt)*in[5];
  // const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  // const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[6];

    // XXX TODO TESTING
    const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[(CeedInt)attr[i]]];

    MultAtBA33(J + i, Q, coeff_mass, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    const CeedScalar u2 = u[i + Q * 2];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[1] * u1 + qd[2] * u2);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[3] * u1 + qd[4] * u2);
    v[i + Q * 2] = wdetJ[i] * (qd[2] * u0 + qd[4] * u1 + qd[5] * u2);
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {

    // XXX TODO TESTING
    const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[(CeedInt)attr[i]]];

    divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2mass_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *qw = in[2], *u = in[3], *divu = in[4];
  CeedScalar *v = out[0], *divv = out[1];

  const CeedScalar *attr = in[5];
  MatCoeffPairContext21 *bc = (MatCoeffPairContext21 *)ctx;
  // const CeedInt attr = (CeedInt)*in[5];
  // const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  // const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[1];

    // XXX TODO TESTING
    const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[(CeedInt)attr[i]]];

    MultAtBA21(J + i, Q, coeff_mass, qd);

    const CeedScalar u0 = u[i + Q * 0];
    v[i + Q * 0] = wdetJ[i] * qd[0] * u0;
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {

    // XXX TODO TESTING
    const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[(CeedInt)attr[i]]];

    divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2mass_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *qw = in[2], *u = in[3], *divu = in[4];
  CeedScalar *v = out[0], *divv = out[1];

  const CeedScalar *attr = in[5];
  MatCoeffPairContext31 *bc = (MatCoeffPairContext31 *)ctx;
  // const CeedInt attr = (CeedInt)*in[5];
  // const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  // const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[3];

    // XXX TODO TESTING
    const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[(CeedInt)attr[i]]];

    MultAtBA32(J + i, Q, coeff_mass, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[1] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[2] * u1);
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {

    // XXX TODO TESTING
    const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[(CeedInt)attr[i]]];

    divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_MASS_QF_H
