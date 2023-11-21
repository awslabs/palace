// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_H1D_QF_H
#define PALACE_LIBCEED_HCURL_H1D_QF_H

#include "utils_qf.h"

// libCEED QFunctions for mixed H(curl)-(H1)ᵈ operators (Piola transformation u =
// adj(J)^T / det(J) ̂u and u = ̂u)
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is transpose adjugate Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [ncomp=space_dim, Q]

CEED_QFUNCTION(f_apply_hcurlh1d_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                    CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2];
  CeedScalar *v = out[0];
  const CeedScalar coeff[3] = {1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
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
  const CeedScalar coeff[6] = {1.0, 0.0, 0.0,
                               1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
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
  const CeedScalar coeff[3] = {1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
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
  const CeedScalar coeff[6] = {1.0, 0.0, 0.0,
                               1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
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

// XX TODO REMOVE AND ADD COEFFICIENT SUPPORT ABOVE

// struct GradContext
// {
//   CeedInt dim, space_dim;
//   CeedScalar coeff;
// };

// // libCEED QFunction for building quadrature data for a gradient operator with a scalar
// // constant coefficient.
// CEED_QFUNCTION(f_build_grad_const_scalar)(void *ctx, CeedInt Q, const CeedScalar *const
// *in,
//                                           CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw C adj(J)^T and store the result.
//   // in[0] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[1] is quadrature weights, size (Q)
//   GradContext *bc = (GradContext *)ctx;
//   const CeedScalar coeff = bc->coeff;
//   const CeedScalar *J = in[0], *qw = in[1];
//   CeedScalar *qd = out[0];
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 11:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qd[i] = qw[i] * coeff;
//       }
//       break;
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt21(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt22(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt32(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt33(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a gradient operator with a scalar
// // coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_grad_quad_scalar)(void *ctx, CeedInt Q, const CeedScalar *const
// *in,
//                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw C adj(J)^T and store the result.
//   // in[0] is coefficients, size (Q)
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   GradContext *bc = (GradContext *)ctx;
//   const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
//   CeedScalar *qd = out[0];
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 11:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qd[i] = qw[i] * c[i];
//       }
//       break;
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt21(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt22(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt32(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt33(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a gradient operator with a vector
// // coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_grad_quad_vector)(void *ctx, CeedInt Q, const CeedScalar *const
// *in,
//                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw C adj(J)^T and store the result.
//   // in[0] is coefficients with shape [ncomp=vdim, Q]
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   GradContext *bc = (GradContext *)ctx;
//   const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
//   CeedScalar *qd = out[0];
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt21(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt22(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt32(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt33(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a gradient operator with a matrix
// // coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_grad_quad_matrix)(void *ctx, CeedInt Q, const CeedScalar *const
// *in,
//                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw C adj(J)^T and store the result.
//   // in[0] is coefficients with shape [ncomp=vdim*(vdim+1)/2, Q]
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   GradContext *bc = (GradContext *)ctx;
//   const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
//   CeedScalar *qd = out[0];
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt21(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt22(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt32(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultCAdjJt33(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for applying a gradient operator.
// CEED_QFUNCTION(f_apply_grad)(void *ctx, CeedInt Q, const CeedScalar *const *in,
//                              CeedScalar *const *out)
// {
//   // in[0] has shape [dim, ncomp=1, Q]
//   // out[0] has shape [ncomp=space_dim, Q]
//   GradContext *bc = (GradContext *)ctx;
//   const CeedScalar *ug = in[0], *qd = in[1];
//   CeedScalar *v = out[0];
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 11:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         v[i] = qd[i] * ug[i];
//       }
//       break;
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar ug0 = ug[i + Q * 0];
//         v[i + Q * 0] = qd[i + Q * 0] * ug0;
//         v[i + Q * 1] = qd[i + Q * 1] * ug0;
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar ug0 = ug[i + Q * 0];
//         const CeedScalar ug1 = ug[i + Q * 1];
//         v[i + Q * 0] = qd[i + Q * 0] * ug0 + qd[i + Q * 2] * ug1;
//         v[i + Q * 1] = qd[i + Q * 1] * ug0 + qd[i + Q * 3] * ug1;
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar ug0 = ug[i + Q * 0];
//         const CeedScalar ug1 = ug[i + Q * 1];
//         v[i + Q * 0] = qd[i + Q * 0] * ug0 + qd[i + Q * 3] * ug1;
//         v[i + Q * 1] = qd[i + Q * 1] * ug0 + qd[i + Q * 4] * ug1;
//         v[i + Q * 2] = qd[i + Q * 2] * ug0 + qd[i + Q * 5] * ug1;
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar ug0 = ug[i + Q * 0];
//         const CeedScalar ug1 = ug[i + Q * 1];
//         const CeedScalar ug2 = ug[i + Q * 2];
//         v[i + Q * 0] = qd[i + Q * 0] * ug0 + qd[i + Q * 3] * ug1 + qd[i + Q * 6] * ug2;
//         v[i + Q * 1] = qd[i + Q * 1] * ug0 + qd[i + Q * 4] * ug1 + qd[i + Q * 7] * ug2;
//         v[i + Q * 2] = qd[i + Q * 2] * ug0 + qd[i + Q * 5] * ug1 + qd[i + Q * 8] * ug2;
//       }
//       break;
//   }
//   return 0;
// }

#endif  // PALACE_LIBCEED_HCURL_H1D_QF_H
