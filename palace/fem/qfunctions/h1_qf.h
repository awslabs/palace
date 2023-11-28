// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_H1_QF_H
#define PALACE_LIBCEED_H1_QF_H

#include "types_qf.h"

// libCEED QFunctions for H1 operators (Piola transformation u = ̂u).
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is active vector, shape [ncomp=vdim, Q]
// in[2] is element attribute, shape [1]
// out[0] is active vector, shape [ncomp=vdim, Q]

CEED_QFUNCTION(f_apply_h1_1)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *u = in[1];
  CeedScalar *v = out[0];

  MatCoeffContext1 *bc = (MatCoeffContext1 *)ctx;
  const CeedInt attr = (CeedInt)*in[2];
  const CeedScalar coeff = *bc->mat_coeff[bc->attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = coeff * wdetJ[i] * u[i];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_h1_2)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *u = in[1];
  CeedScalar *v = out[0];

  MatCoeffContext2 *bc = (MatCoeffContext2 *)ctx;
  const CeedInt attr = (CeedInt)*in[2];
  const CeedScalar *coeff = bc->mat_coeff[bc->attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (coeff[0] * u0 + coeff[1] * u1);
    v[i + Q * 1] = wdetJ[i] * (coeff[1] * u0 + coeff[2] * u1);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_h1_3)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *u = in[1];
  CeedScalar *v = out[0];

  MatCoeffContext2 *bc = (MatCoeffContext2 *)ctx;
  const CeedInt attr = (CeedInt)*in[2];
  const CeedScalar *coeff = bc->mat_coeff[bc->attr_mat[attr]];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    const CeedScalar u2 = u[i + Q * 2];
    v[i + Q * 0] = wdetJ[i] * (coeff[0] * u0 + coeff[1] * u1 + coeff[2] * u2);
    v[i + Q * 1] = wdetJ[i] * (coeff[1] * u0 + coeff[3] * u1 + coeff[4] * u2);
    v[i + Q * 2] = wdetJ[i] * (coeff[2] * u0 + coeff[4] * u1 + coeff[5] * u2);
  }
  return 0;
}

// XX TODO REMOVE AND ADD COEFFICIENT SUPPORT ABOVE

// #include "utils_qf.h"

// struct MassContext
// {
//   CeedInt dim, space_dim, vdim;
//   CeedScalar coeff;
// };

// // libCEED QFunction for building quadrature data for a mass operator with a scalar
// // coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_mass_quad_scalar)(void *ctx, CeedInt Q, const CeedScalar *const
// *in,
//                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute and store qw * c * det(J).
//   // in[0] is coefficients, size (Q)
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   MassContext *bc = (MassContext *)ctx;
//   const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
//   CeedScalar *qd = out[0];
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 11:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qd[i] = qw[i] * c[i] * J[i];
//       }
//       break;
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qd[i] = qw[i] * c[i] * DetJ21(J + i, Q);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qd[i] = qw[i] * c[i] * DetJ22(J + i, Q);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qd[i] = qw[i] * c[i] * DetJ32(J + i, Q);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qd[i] = qw[i] * c[i] * DetJ33(J + i, Q);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a mass operator with a vector
// // coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_mass_quad_vector)(void *ctx, CeedInt Q, const CeedScalar *const
// *in,
//                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute and store qw * det(J) C.
//   // in[0] is coefficients with shape [ncomp=vdim, Q]
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   MassContext *bc = (MassContext *)ctx;
//   const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
//   CeedScalar *qd = out[0];
//   switch (100 * bc->space_dim + 10 * bc->dim + bc->vdim)
//   {
//     case 212:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar wdetJi = qw[i] * DetJ21(J + i, Q);
//         CeedPragmaSIMD for (CeedInt d = 0; d < 2; d++)
//         {
//           qd[i + Q * d] = wdetJi * c[i + Q * d];
//         }
//       }
//       break;
//     case 222:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar wdetJi = qw[i] * DetJ22(J + i, Q);
//         CeedPragmaSIMD for (CeedInt d = 0; d < 2; d++)
//         {
//           qd[i + Q * d] = wdetJi * c[i + Q * d];
//         }
//       }
//       break;
//     case 323:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar wdetJi = qw[i] * DetJ32(J + i, Q);
//         CeedPragmaSIMD for (CeedInt d = 0; d < 3; d++)
//         {
//           qd[i + Q * d] = wdetJi * c[i + Q * d];
//         }
//       }
//       break;
//     case 333:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar wdetJi = qw[i] * DetJ33(J + i, Q);
//         CeedPragmaSIMD for (CeedInt d = 0; d < 3; d++)
//         {
//           qd[i + Q * d] = wdetJi * c[i + Q * d];
//         }
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a mass operator with a matrix
// // coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_mass_quad_matrix)(void *ctx, CeedInt Q, const CeedScalar *const
// *in,
//                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute and store qw * det(J) C.
//   // in[0] is coefficients with shape [ncomp=vdim*(vdim+1)/2, Q]
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   MassContext *bc = (MassContext *)ctx;
//   const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
//   CeedScalar *qd = out[0];
//   switch (100 * bc->space_dim + 10 * bc->dim + bc->vdim)
//   {
//     case 212:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar wdetJi = qw[i] * DetJ21(J + i, Q);
//         CeedPragmaSIMD for (CeedInt d = 0; d < 3; d++)
//         {
//           qd[i + Q * d] = wdetJi * c[i + Q * d];
//         }
//       }
//       break;
//     case 222:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar wdetJi = qw[i] * DetJ22(J + i, Q);
//         CeedPragmaSIMD for (CeedInt d = 0; d < 3; d++)
//         {
//           qd[i + Q * d] = wdetJi * c[i + Q * d];
//         }
//       }
//       break;
//     case 323:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar wdetJi = qw[i] * DetJ32(J + i, Q);
//         CeedPragmaSIMD for (CeedInt d = 0; d < 6; d++)
//         {
//           qd[i + Q * d] = wdetJi * c[i + Q * d];
//         }
//       }
//       break;
//     case 333:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar wdetJi = qw[i] * DetJ33(J + i, Q);
//         CeedPragmaSIMD for (CeedInt d = 0; d < 6; d++)
//         {
//           qd[i + Q * d] = wdetJi * c[i + Q * d];
//         }
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for applying a mass operator with a scalar coefficient.
// CEED_QFUNCTION(f_apply_mass_scalar)(void *ctx, CeedInt Q, const CeedScalar *const *in,
//                                     CeedScalar *const *out)
// {
//   // in[0], out[0] have shape [ncomp=vdim, Q]
//   MassContext *bc = (MassContext *)ctx;
//   const CeedScalar *u = in[0], *qd = in[1];
//   CeedScalar *v = out[0];
//   switch (bc->vdim)
//   {
//     case 1:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         v[i] = qd[i] * u[i];
//       }
//       break;
//     case 2:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar qdi = qd[i];
//         CeedPragmaSIMD for (CeedInt d = 0; d < 2; d++)
//         {
//           v[i + Q * d] = qdi * u[i + Q * d];
//         }
//       }
//       break;
//     case 3:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar qdi = qd[i];
//         CeedPragmaSIMD for (CeedInt d = 0; d < 3; d++)
//         {
//           v[i + Q * d] = qdi * u[i + Q * d];
//         }
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for applying a mass operator with a vector coefficient.
// CEED_QFUNCTION(f_apply_mass_vector)(void *ctx, CeedInt Q, const CeedScalar *const *in,
//                                     CeedScalar *const *out)
// {
//   // in[0], out[0] have shape [ncomp=vdim, Q]
//   MassContext *bc = (MassContext *)ctx;
//   const CeedScalar *u = in[0], *qd = in[1];
//   CeedScalar *v = out[0];
//   switch (bc->vdim)
//   {
//     case 2:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         CeedPragmaSIMD for (CeedInt d = 0; d < 2; d++)
//         {
//           v[i + Q * d] = qd[i + Q * d] * u[i + Q * d];
//         }
//       }
//       break;
//     case 3:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         CeedPragmaSIMD for (CeedInt d = 0; d < 3; d++)
//         {
//           v[i + Q * d] = qd[i + Q * d] * u[i + Q * d];
//         }
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for applying a mass operator with a matrix coefficient.
// CEED_QFUNCTION(f_apply_mass_matrix)(void *ctx, CeedInt Q, const CeedScalar *const *in,
//                                     CeedScalar *const *out)
// {
//   // in[0], out[0] have shape [ncomp=vdim, Q]
//   MassContext *bc = (MassContext *)ctx;
//   const CeedScalar *u = in[0], *qd = in[1];
//   CeedScalar *v = out[0];
//   switch (bc->vdim)
//   {
//     case 2:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar u0 = u[i + Q * 0];
//         const CeedScalar u1 = u[i + Q * 1];
//         v[i + Q * 0] = qd[i + Q * 0] * u0 + qd[i + Q * 1] * u1;
//         v[i + Q * 1] = qd[i + Q * 1] * u0 + qd[i + Q * 2] * u1;
//       }
//       break;
//     case 3:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar u0 = u[i + Q * 0];
//         const CeedScalar u1 = u[i + Q * 1];
//         const CeedScalar u2 = u[i + Q * 2];
//         v[i + Q * 0] = qd[i + Q * 0] * u0 + qd[i + Q * 1] * u1 + qd[i + Q * 2] * u2;
//         v[i + Q * 1] = qd[i + Q * 1] * u0 + qd[i + Q * 3] * u1 + qd[i + Q * 4] * u2;
//         v[i + Q * 2] = qd[i + Q * 2] * u0 + qd[i + Q * 4] * u1 + qd[i + Q * 5] * u2;
//       }
//       break;
//   }
//   return 0;
// }

#endif  // PALACE_LIBCEED_H1_QF_H
