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

  MatCoeffPairContext21 *bc = (MatCoeffPairContext21 *)ctx;
  const CeedInt attr = (CeedInt)*in[5];
  const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[3];
    MultAtBA22(J + i, Q, coeff_mass, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[1] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[2] * u1);
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2mass_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *qw = in[2], *u = in[3], *divu = in[4];
  CeedScalar *v = out[0], *divv = out[1];

  MatCoeffPairContext31 *bc = (MatCoeffPairContext31 *)ctx;
  const CeedInt attr = (CeedInt)*in[5];
  const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[6];
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
    divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2mass_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *qw = in[2], *u = in[3], *divu = in[4];
  CeedScalar *v = out[0], *divv = out[1];

  MatCoeffPairContext21 *bc = (MatCoeffPairContext21 *)ctx;
  const CeedInt attr = (CeedInt)*in[5];
  const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[1];
    MultAtBA21(J + i, Q, coeff_mass, qd);

    const CeedScalar u0 = u[i + Q * 0];
    v[i + Q * 0] = wdetJ[i] * qd[0] * u0;
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_l2mass_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *qw = in[2], *u = in[3], *divu = in[4];
  CeedScalar *v = out[0], *divv = out[1];

  MatCoeffPairContext31 *bc = (MatCoeffPairContext31 *)ctx;
  const CeedInt attr = (CeedInt)*in[5];
  const CeedScalar *coeff_mass = bc->ctx1.mat_coeff[bc->ctx1.attr_mat[attr]];
  const CeedScalar coeff = *bc->ctx2.mat_coeff[bc->ctx2.attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[3];
    MultAtBA32(J + i, Q, coeff_mass, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[1] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[2] * u1);
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    divv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * divu[i];
  }
  return 0;
}

// XX TODO REMOVE AND COEFFICIENTS

// struct DivDivMassContext
// {
//   CeedInt dim, space_dim;
// };

// // libCEED QFunction for building quadrature data for a div-div + mass operator with a
// // scalar coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_divdiv_mass_quad_scalar)(void *ctx, CeedInt Q,
//                                                 const CeedScalar *const *in,
//                                                 CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw * c / det(J) and qw / det(J) J^T C J and store
//   // the result.
//   // in[0] is div-div coefficients with shape [ncomp=1, Q]
//   // in[1] is mass coefficients with shape [ncomp=1, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   DivDivMassContext *bc = (DivDivMassContext *)ctx;
//   const CeedScalar *cd = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdd = out[0], *qdm = out[0] + Q;
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 11:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / J[i];
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * J[i];
//       }
//       break;
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ21(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ21(J + i, Q, cm + i, Q, 1, qw[i], Q, qdm + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ22(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ22(J + i, Q, cm + i, Q, 1, qw[i], Q, qdm + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ32(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ32(J + i, Q, cm + i, Q, 1, qw[i], Q, qdm + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ33(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cm + i, Q, 1, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a div-div + mass operator with a
// // vector coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_divdiv_mass_quad_vector)(void *ctx, CeedInt Q,
//                                                 const CeedScalar *const *in,
//                                                 CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw * c / det(J) and qw / det(J) J^T C J and store
//   // the result.
//   // in[0] is div-div coefficients with shape [ncomp=1, Q]
//   // in[1] is mass coefficients with shape [ncomp=space_dim, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   DivDivMassContext *bc = (DivDivMassContext *)ctx;
//   const CeedScalar *cd = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdd = out[0], *qdm = out[0] + Q;
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ21(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ21(J + i, Q, cm + i, Q, 2, qw[i], Q, qdm + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ22(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ22(J + i, Q, cm + i, Q, 2, qw[i], Q, qdm + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ32(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ32(J + i, Q, cm + i, Q, 3, qw[i], Q, qdm + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ33(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cm + i, Q, 3, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a div-div + mass operator with a
// // matrix coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_divdiv_mass_quad_matrix)(void *ctx, CeedInt Q,
//                                                 const CeedScalar *const *in,
//                                                 CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw * c / det(J) and qw / det(J) J^T C J and store
//   // the result.
//   // in[0] is div-div coefficients with shape [ncomp=1, Q]
//   // in[1] is mass coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   DivDivMassContext *bc = (DivDivMassContext *)ctx;
//   const CeedScalar *cd = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdd = out[0], *qdm = out[0] + Q;
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ21(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ21(J + i, Q, cm + i, Q, 3, qw[i], Q, qdm + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ22(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ22(J + i, Q, cm + i, Q, 3, qw[i], Q, qdm + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ32(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ32(J + i, Q, cm + i, Q, 6, qw[i], Q, qdm + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / DetJ33(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cm + i, Q, 6, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for applying a div-div + mass operator.
// CEED_QFUNCTION(f_apply_divdiv_mass)(void *ctx, CeedInt Q, const CeedScalar *const *in,
//                                     CeedScalar *const *out)
// {
//   // in[0], out[0] have shape [dim, ncomp=1, Q]
//   // in[1], out[1] have shape [ncomp=1, Q]
//   DivDivMassContext *bc = (DivDivMassContext *)ctx;
//   const CeedScalar *u = in[0], *ud = in[1], *qdd = in[2], *qdm = in[2] + Q;
//   CeedScalar *v = out[0], *vd = out[1];
//   CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//   {
//     vd[i] = qdd[i] * ud[i];
//   }
//   switch (bc->dim)
//   {
//     case 1:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         v[i] = qdm[i] * u[i];
//       }
//       break;
//     case 2:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar u0 = u[i + Q * 0];
//         const CeedScalar u1 = u[i + Q * 1];
//         v[i + Q * 0] = qdm[i + Q * 0] * u0 + qdm[i + Q * 1] * u1;
//         v[i + Q * 1] = qdm[i + Q * 1] * u0 + qdm[i + Q * 2] * u1;
//       }
//       break;
//     case 3:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar u0 = u[i + Q * 0];
//         const CeedScalar u1 = u[i + Q * 1];
//         const CeedScalar u2 = u[i + Q * 2];
//         v[i + Q * 0] = qdm[i + Q * 0] * u0 + qdm[i + Q * 1] * u1 + qdm[i + Q * 2] * u2;
//         v[i + Q * 1] = qdm[i + Q * 1] * u0 + qdm[i + Q * 3] * u1 + qdm[i + Q * 4] * u2;
//         v[i + Q * 2] = qdm[i + Q * 2] * u0 + qdm[i + Q * 4] * u1 + qdm[i + Q * 5] * u2;
//       }
//       break;
//   }
//   return 0;
// }

#endif  // PALACE_LIBCEED_L2_MASS_QF_H
