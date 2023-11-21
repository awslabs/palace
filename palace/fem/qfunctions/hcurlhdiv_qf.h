// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_HDIV_QF_H
#define PALACE_LIBCEED_HCURL_HDIV_QF_H

#include "utils_qf.h"

// libCEED QFunctions for mixed H(curl)-H(div) operators (Piola transformations u =
// adj(J)^T / det(J) ̂u and u = J / det(J) ̂u).
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is transpose adjugate Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[3] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]

CEED_QFUNCTION(f_apply_hcurlhdiv_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *adjJt = in[2], *u = in[3];
  CeedScalar *v = out[0];
  const CeedScalar coeff[3] = {1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[4];
    MultAtBC22(J + i, Q, coeff, adjJt + i, Q, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[2] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[3] * u1);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlhdiv_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *adjJt = in[2], *u = in[3];
  CeedScalar *v = out[0];
  const CeedScalar coeff[6] = {1.0, 0.0, 0.0,
                               1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[9];
    MultAtBC33(J + i, Q, coeff, adjJt + i, Q, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    const CeedScalar u2 = u[i + Q * 2];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[3] * u1 + qd[6] * u2);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[4] * u1 + qd[7] * u2);
    v[i + Q * 2] = wdetJ[i] * (qd[2] * u0 + qd[5] * u1 + qd[8] * u2);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlhdiv_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *adjJt = in[2], *u = in[3];
  CeedScalar *v = out[0];
  const CeedScalar coeff[3] = {1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[1];
    MultAtBC21(J + i, Q, coeff, adjJt + i, Q, qd);

    const CeedScalar u0 = u[i + Q * 0];
    v[i + Q * 0] = wdetJ[i] * qd[0] * u0;
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlhdiv_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *adjJt = in[2], *u = in[3];
  CeedScalar *v = out[0];
  const CeedScalar coeff[6] = {1.0, 0.0, 0.0,
                               1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[4];
    MultAtBC32(J + i, Q, coeff, adjJt + i, Q, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[2] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[3] * u1);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivhcurl_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *adjJt = in[2], *u = in[3];
  CeedScalar *v = out[0];
  const CeedScalar coeff[3] = {1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[4];
    MultAtBC22(adjJt + i, Q, coeff, J + i, Q, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[2] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[3] * u1);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivhcurl_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *adjJt = in[2], *u = in[3];
  CeedScalar *v = out[0];
  const CeedScalar coeff[6] = {1.0, 0.0, 0.0,
                               1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[9];
    MultAtBC33(adjJt + i, Q, coeff, J + i, Q, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    const CeedScalar u2 = u[i + Q * 2];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[3] * u1 + qd[6] * u2);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[4] * u1 + qd[7] * u2);
    v[i + Q * 2] = wdetJ[i] * (qd[2] * u0 + qd[5] * u1 + qd[8] * u2);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivhcurl_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *adjJt = in[2], *u = in[3];
  CeedScalar *v = out[0];
  const CeedScalar coeff[3] = {1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[1];
    MultAtBC21(adjJt + i, Q, coeff, J + i, Q, qd);

    const CeedScalar u0 = u[i + Q * 0];
    v[i + Q * 0] = wdetJ[i] * qd[0] * u0;
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivhcurl_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *J = in[1], *adjJt = in[2], *u = in[3];
  CeedScalar *v = out[0];
  const CeedScalar coeff[6] = {1.0, 0.0, 0.0,
                               1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[4];
    MultAtBC32(adjJt + i, Q, coeff, J + i, Q, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[2] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[3] * u1);
  }
  return 0;
}

// XX TODO REMOVE AND ADD COEFFICIENT SUPPORT ABOVE

// // libCEED QFunction for building quadrature data for a mixed H(curl)-H(div) mass
// operator
// // with a scalar constant coefficient.
// CEED_QFUNCTION(f_build_hcurlhdiv_const_scalar)(void *ctx, CeedInt Q,
//                                                const CeedScalar *const *in,
//                                                CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) adj(J) C J and store the
//   // result.
//   // in[0] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[1] is quadrature weights, size (Q)
//   VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
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
//         MultJtCAdjJt21(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt22(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt32(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt33(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a mixed H(curl)-H(div) mass
// operator
// // with a scalar coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_hcurlhdiv_quad_scalar)(void *ctx, CeedInt Q,
//                                               const CeedScalar *const *in,
//                                               CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) adj(J) C J and store the
//   // result.
//   // in[0] is coefficients with shape [ncomp=1, Q]
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
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
//         MultJtCAdjJt21(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt22(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt32(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt33(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a mixed H(curl)-H(div) mass
// operator
// // with a vector coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_hcurlhdiv_quad_vector)(void *ctx, CeedInt Q,
//                                               const CeedScalar *const *in,
//                                               CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) adj(J) C J and store the
//   // result.
//   // in[0] is coefficients with shape [ncomp=space_dim, Q]
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
//   const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
//   CeedScalar *qd = out[0];
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt21(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt22(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt32(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt33(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a mixed H(curl)-H(div) mass
// operator
// // with a matrix coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_hcurlhdiv_quad_matrix)(void *ctx, CeedInt Q,
//                                               const CeedScalar *const *in,
//                                               CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) adj(J) C J and store the
//   // result.
//   // in[0] is coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
//   const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
//   CeedScalar *qd = out[0];
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt21(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt22(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt32(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt33(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a mixed H(div)-H(curl) mass
// operator
// // with a scalar constant coefficient.
// CEED_QFUNCTION(f_build_hdivhcurl_const_scalar)(void *ctx, CeedInt Q,
//                                                const CeedScalar *const *in,
//                                                CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) adj(J) C J and store the
//   // result.
//   // in[0] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[1] is quadrature weights, size (Q)
//   VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
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
//         MultJtCAdjJt21<true>(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt22<true>(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt32<true>(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt33<true>(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a mixed H(div)-H(curl) mass
// operator
// // with a scalar coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_hdivhcurl_quad_scalar)(void *ctx, CeedInt Q,
//                                               const CeedScalar *const *in,
//                                               CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) adj(J) C J and store the
//   // result.
//   // in[0] is coefficients with shape [ncomp=1, Q]
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
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
//         MultJtCAdjJt21<true>(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt22<true>(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt32<true>(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt33<true>(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a mixed H(div)-H(curl) mass
// operator
// // with a vector coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_hdivhcurl_quad_vector)(void *ctx, CeedInt Q,
//                                               const CeedScalar *const *in,
//                                               CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) adj(J) C J and store the
//   // result.
//   // in[0] is coefficients with shape [ncomp=space_dim, Q]
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
//   const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
//   CeedScalar *qd = out[0];
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt21<true>(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt22<true>(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt32<true>(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt33<true>(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a mixed H(div)-H(curl) mass
// operator
// // with a matrix coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_hdivhcurl_quad_matrix)(void *ctx, CeedInt Q,
//                                               const CeedScalar *const *in,
//                                               CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) adj(J) C J and store the
//   // result.
//   // in[0] is coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
//   // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[2] is quadrature weights, size (Q)
//   VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
//   const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
//   CeedScalar *qd = out[0];
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt21<true>(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt22<true>(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt32<true>(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCAdjJt33<true>(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
//       }
//       break;
//   }
//   return 0;
// }

#endif  // PALACE_LIBCEED_HCURL_HDIV_QF_H
