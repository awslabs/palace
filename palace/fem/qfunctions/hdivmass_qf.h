// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HDIV_MASS_QF_H
#define PALACE_LIBCEED_HDIV_MASS_QF_H

#include "utils_qf.h"

// libCEED QFunctions for H(div) + H(curl) mass operators in 3D (Piola transformations u =
// J / det(J) ̂u and u = adj(J)^T / det(J) ̂u).
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is transpose adjugate Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[3] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[4] is active vector curl, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[1] is active vector curl, shape [qcomp=dim, ncomp=1, Q]

// In 2D, this actually uses the L2 Piola transformation on the curl (u = 1 / det(J) ̂u) and
// the curl is has qcomp=1.
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is transpose adjugate Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is quadrature weights, shape [Q]
// in[3] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[4] is active vector curl, shape [ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[1] is active vector curl, shape [ncomp=1, Q]

CEED_QFUNCTION(f_apply_hdivmass_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                    CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *qw = in[2], *u = in[3], *curlu = in[4];
  CeedScalar *v = out[0], *curlv = out[1];
  const CeedScalar coeff_mass[3] = {1.0, 0.0, 1.0};
  const CeedScalar coeff = 1.0;  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[3];
    MultAtBA22(adjJt + i, Q, coeff_mass, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[1] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[3] * u1);
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    curlv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * curlu[i];
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivmass_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                    CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *J = in[2], *u = in[3], *curlu = in[4];
  CeedScalar *v = out[0], *curlv = out[1];
  const CeedScalar coeff_mass[6] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  const CeedScalar coeff[6] = {1.0, 0.0, 0.0,
                               1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[6];
    MultAtBA33(adjJt + i, Q, coeff_mass, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    const CeedScalar u2 = u[i + Q * 2];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[1] * u1 + qd[2] * u2);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[3] * u1 + qd[4] * u2);
    v[i + Q * 2] = wdetJ[i] * (qd[2] * u0 + qd[4] * u1 + qd[5] * u2);
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[6];
    MultAtBA33(J + i, Q, coeff, qd);

    const CeedScalar curlu0 = curlu[i + Q * 0];
    const CeedScalar curlu1 = curlu[i + Q * 1];
    const CeedScalar curlu2 = curlu[i + Q * 2];
    curlv[i + Q * 0] = wdetJ[i] * (qd[0] * curlu0 + qd[1] * curlu1 + qd[2] * curlu2);
    curlv[i + Q * 1] = wdetJ[i] * (qd[1] * curlu0 + qd[3] * curlu1 + qd[4] * curlu2);
    curlv[i + Q * 2] = wdetJ[i] * (qd[2] * curlu0 + qd[4] * curlu1 + qd[5] * curlu2);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hdivmass_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                    CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *qw = in[2], *u = in[3], *curlu = in[4];
  CeedScalar *v = out[0], *curlv = out[1];
  const CeedScalar coeff_mass[6] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  const CeedScalar coeff = 1.0;  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[3];
    MultAtBA32(adjJt + i, Q, coeff_mass, qd);

    const CeedScalar u0 = u[i + Q * 0];
    const CeedScalar u1 = u[i + Q * 1];
    v[i + Q * 0] = wdetJ[i] * (qd[0] * u0 + qd[1] * u1);
    v[i + Q * 1] = wdetJ[i] * (qd[1] * u0 + qd[2] * u1);
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    curlv[i] = (coeff * qw[i] * qw[i] / wdetJ[i]) * curlu[i];
  }
  return 0;
}

// XX TODO REMOVE AND COEFFICIENTS

// struct CurlCurlMassContext
// {
//   CeedInt dim, space_dim, curl_dim;
// };

// // libCEED QFunction for building quadrature data for a curl-curl + mass operator with
// // scalar coefficients evaluated at quadrature points.
// CEED_QFUNCTION(f_build_curlcurl_mass_quad_scalar_scalar)(void *ctx, CeedInt Q,
//                                                          const CeedScalar *const *in,
//                                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) J^T C J (3D) or qw * c / det(J) (2D)
//   and
//   // qw / det(J) adj(J) C adj(J)^T and store the result.
//   // in[0] is curl-curl coefficients with shape [ncomp=1, Q]
//   // in[1] is mass coefficients with shape [ncomp=1, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   CurlCurlMassContext *bc = (CurlCurlMassContext *)ctx;
//   const CeedScalar *cc = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdc = out[0], *qdm = out[0] + Q * bc->curl_dim * (bc->curl_dim + 1) / 2;
//   switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
//   {
//     case 221:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdc[i] = qw[i] * cc[i] / DetJ22(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt22(J + i, Q, cm + i, Q, 1, qw[i], Q, qdm + i);
//       }
//       break;
//     case 321:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdc[i] = qw[i] * cc[i] / DetJ32(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt32(J + i, Q, cm + i, Q, 1, qw[i], Q, qdm + i);
//       }
//       break;
//     case 333:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cc + i, Q, 1, qw[i], Q, qdc + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cm + i, Q, 1, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a curl-curl + mass operator with
// // scalar and vector coefficients evaluated at quadrature points.
// CEED_QFUNCTION(f_build_curlcurl_mass_quad_scalar_vector)(void *ctx, CeedInt Q,
//                                                          const CeedScalar *const *in,
//                                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) J^T C J (3D) or qw * c / det(J) (2D)
//   and
//   // qw / det(J) adj(J) C adj(J)^T and store the result.
//   // in[0] is curl-curl coefficients with shape [ncomp=1, Q]
//   // in[1] is mass coefficients with shape [ncomp=space_dim, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   CurlCurlMassContext *bc = (CurlCurlMassContext *)ctx;
//   const CeedScalar *cc = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdc = out[0], *qdm = out[0] + Q * bc->curl_dim * (bc->curl_dim + 1) / 2;
//   switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
//   {
//     case 221:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdc[i] = qw[i] * cc[i] / DetJ22(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt22(J + i, Q, cm + i, Q, 2, qw[i], Q, qdm + i);
//       }
//       break;
//     case 321:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdc[i] = qw[i] * cc[i] / DetJ32(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt32(J + i, Q, cm + i, Q, 3, qw[i], Q, qdm + i);
//       }
//       break;
//     case 333:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cc + i, Q, 1, qw[i], Q, qdc + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cm + i, Q, 3, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a curl-curl + mass operator with
// // scalar and matrix coefficients evaluated at quadrature points.
// CEED_QFUNCTION(f_build_curlcurl_mass_quad_scalar_matrix)(void *ctx, CeedInt Q,
//                                                          const CeedScalar *const *in,
//                                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) J^T C J (3D) or qw * c / det(J) (2D)
//   and
//   // qw / det(J) adj(J) C adj(J)^T and store the result.
//   // in[0] is curl-curl coefficients with shape [ncomp=1, Q]
//   // in[1] is mass coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   CurlCurlMassContext *bc = (CurlCurlMassContext *)ctx;
//   const CeedScalar *cc = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdc = out[0], *qdm = out[0] + Q * bc->curl_dim * (bc->curl_dim + 1) / 2;
//   switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
//   {
//     case 221:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdc[i] = qw[i] * cc[i] / DetJ22(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt22(J + i, Q, cm + i, Q, 3, qw[i], Q, qdm + i);
//       }
//       break;
//     case 321:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdc[i] = qw[i] * cc[i] / DetJ32(J + i, Q);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt32(J + i, Q, cm + i, Q, 6, qw[i], Q, qdm + i);
//       }
//       break;
//     case 333:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cc + i, Q, 1, qw[i], Q, qdc + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cm + i, Q, 6, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a curl-curl + mass operator with
// // vector and scalar coefficients evaluated at quadrature points.
// CEED_QFUNCTION(f_build_curlcurl_mass_quad_vector_scalar)(void *ctx, CeedInt Q,
//                                                          const CeedScalar *const *in,
//                                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) J^T C J (3D) or qw * c / det(J) (2D)
//   and
//   // qw / det(J) adj(J) C adj(J)^T and store the result.
//   // in[0] is curl-curl coefficients with shape [ncomp=space_dim, Q]
//   // in[1] is mass coefficients with shape [ncomp=1, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   CurlCurlMassContext *bc = (CurlCurlMassContext *)ctx;
//   const CeedScalar *cc = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdc = out[0], *qdm = out[0] + Q * bc->curl_dim * (bc->curl_dim + 1) / 2;
//   switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
//   {
//     case 333:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cc + i, Q, 3, qw[i], Q, qdc + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cm + i, Q, 1, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a curl-curl + mass operator with
// // vector coefficients evaluated at quadrature points.
// CEED_QFUNCTION(f_build_curlcurl_mass_quad_vector_vector)(void *ctx, CeedInt Q,
//                                                          const CeedScalar *const *in,
//                                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) J^T C J (3D) or qw * c / det(J) (2D)
//   and
//   // qw / det(J) adj(J) C adj(J)^T and store the result.
//   // in[0] is curl-curl coefficients with shape [ncomp=space_dim, Q]
//   // in[1] is mass coefficients with shape [ncomp=space_dim, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   CurlCurlMassContext *bc = (CurlCurlMassContext *)ctx;
//   const CeedScalar *cc = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdc = out[0], *qdm = out[0] + Q * bc->curl_dim * (bc->curl_dim + 1) / 2;
//   switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
//   {
//     case 333:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cc + i, Q, 3, qw[i], Q, qdc + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cm + i, Q, 3, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a curl-curl + mass operator with
// // vector and matrix coefficients evaluated at quadrature points.
// CEED_QFUNCTION(f_build_curlcurl_mass_quad_vector_matrix)(void *ctx, CeedInt Q,
//                                                          const CeedScalar *const *in,
//                                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) J^T C J (3D) or qw * c / det(J) (2D)
//   and
//   // qw / det(J) adj(J) C adj(J)^T and store the result.
//   // in[0] is curl-curl coefficients with shape [ncomp=space_dim, Q]
//   // in[1] is mass coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   CurlCurlMassContext *bc = (CurlCurlMassContext *)ctx;
//   const CeedScalar *cc = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdc = out[0], *qdm = out[0] + Q * bc->curl_dim * (bc->curl_dim + 1) / 2;
//   switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
//   {
//     case 333:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cc + i, Q, 3, qw[i], Q, qdc + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cm + i, Q, 6, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a curl-curl + mass operator with
// // matrix and scalar coefficients evaluated at quadrature points.
// CEED_QFUNCTION(f_build_curlcurl_mass_quad_matrix_scalar)(void *ctx, CeedInt Q,
//                                                          const CeedScalar *const *in,
//                                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) J^T C J (3D) or qw * c / det(J) (2D)
//   and
//   // qw / det(J) adj(J) C adj(J)^T and store the result.
//   // in[0] is curl-curl coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
//   // in[1] is mass coefficients with shape [ncomp=1, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   CurlCurlMassContext *bc = (CurlCurlMassContext *)ctx;
//   const CeedScalar *cc = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdc = out[0], *qdm = out[0] + Q * bc->curl_dim * (bc->curl_dim + 1) / 2;
//   switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
//   {
//     case 333:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cc + i, Q, 6, qw[i], Q, qdc + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cm + i, Q, 1, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a curl-curl + mass operator with
// // matrix and vector coefficients evaluated at quadrature points.
// CEED_QFUNCTION(f_build_curlcurl_mass_quad_matrix_vector)(void *ctx, CeedInt Q,
//                                                          const CeedScalar *const *in,
//                                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) J^T C J (3D) or qw * c / det(J) (2D)
//   and
//   // qw / det(J) adj(J) C adj(J)^T and store the result.
//   // in[0] is curl-curl coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
//   // in[1] is mass coefficients with shape [ncomp=space_dim, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   CurlCurlMassContext *bc = (CurlCurlMassContext *)ctx;
//   const CeedScalar *cc = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdc = out[0], *qdm = out[0] + Q * bc->curl_dim * (bc->curl_dim + 1) / 2;
//   switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
//   {
//     case 333:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cc + i, Q, 6, qw[i], Q, qdc + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cm + i, Q, 3, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a curl-curl + mass operator with
// // matrix coefficients evaluated at quadrature points.
// CEED_QFUNCTION(f_build_curlcurl_mass_quad_matrix_matrix)(void *ctx, CeedInt Q,
//                                                          const CeedScalar *const *in,
//                                                          CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) J^T C J (3D) or qw * c / det(J) (2D)
//   and
//   // qw / det(J) adj(J) C adj(J)^T and store the result.
//   // in[0] is curl-curl coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
//   // in[1] is mass coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   CurlCurlMassContext *bc = (CurlCurlMassContext *)ctx;
//   const CeedScalar *cc = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdc = out[0], *qdm = out[0] + Q * bc->curl_dim * (bc->curl_dim + 1) / 2;
//   switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
//   {
//     case 333:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultJtCJ33(J + i, Q, cc + i, Q, 6, qw[i], Q, qdc + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cm + i, Q, 6, qw[i], Q, qdm + i);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for applying a curl-curl + mass operator.
// CEED_QFUNCTION(f_apply_curlcurl_mass)(void *ctx, CeedInt Q, const CeedScalar *const *in,
//                                       CeedScalar *const *out)
// {
//   // in[0], out[0] have shape [dim, ncomp=1, Q]
//   // in[1], out[1] have shape [curl_dim, ncomp=1, Q]
//   CurlCurlMassContext *bc = (CurlCurlMassContext *)ctx;
//   const CeedScalar *u = in[0], *uc = in[1], *qdc = in[2],
//                    *qdm = in[2] + Q * bc->curl_dim * (bc->curl_dim + 1) / 2;
//   CeedScalar *v = out[0], *vc = out[1];
//   switch (10 * bc->dim + bc->curl_dim)
//   {
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         vc[i] = qdc[i] * uc[i];
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar u0 = u[i + Q * 0];
//         const CeedScalar u1 = u[i + Q * 1];
//         v[i + Q * 0] = qdm[i + Q * 0] * u0 + qdm[i + Q * 1] * u1;
//         v[i + Q * 1] = qdm[i + Q * 1] * u0 + qdm[i + Q * 2] * u1;
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar uc0 = uc[i + Q * 0];
//         const CeedScalar uc1 = uc[i + Q * 1];
//         const CeedScalar uc2 = uc[i + Q * 2];
//         vc[i + Q * 0] = qdc[i + Q * 0] * uc0 + qdc[i + Q * 1] * uc1 + qdc[i + Q * 2] *
//         uc2; vc[i + Q * 1] = qdc[i + Q * 1] * uc0 + qdc[i + Q * 3] * uc1 + qdc[i + Q * 4]
//         * uc2; vc[i + Q * 2] = qdc[i + Q * 2] * uc0 + qdc[i + Q * 4] * uc1 + qdc[i + Q *
//         5] * uc2;
//       }
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

#endif  // PALACE_LIBCEED_CURLCURL_MASS_QF_H
