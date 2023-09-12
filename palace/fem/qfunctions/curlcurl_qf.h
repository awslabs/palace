// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_CURLCURL_QF_H
#define PALACE_LIBCEED_CURLCURL_QF_H

#include "utils_qf.h"

struct CurlCurlContext
{
  CeedInt dim, space_dim, curl_dim;
  CeedScalar coeff;
};

// libCEED QFunction for building quadrature data for a curl-curl operator with a scalar
// constant coefficient.
CEED_QFUNCTION(f_build_curlcurl_const_scalar)(void *ctx, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) J^T C J and store the symmetric part of
  // the result. In 2D, compute and store qw * c / det(J).
  // in[0] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[1] is quadrature weights, size (Q)
  CurlCurlContext *bc = (CurlCurlContext *)ctx;
  const CeedScalar coeff = bc->coeff;
  const CeedScalar *J = in[0], *qw = in[1];
  CeedScalar *qd = out[0];
  switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
  {
    case 221:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        qd[i] = qw[i] * coeff / DetJ22(J + i, Q);
      }
      break;
    case 321:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        qd[i] = qw[i] * coeff / DetJ32(J + i, Q);
      }
      break;
    case 333:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ33(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for a curl-curl operator with a scalar
// coefficient evaluated at quadrature points.
CEED_QFUNCTION(f_build_curlcurl_quad_scalar)(void *ctx, CeedInt Q,
                                             const CeedScalar *const *in,
                                             CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) J^T C J and store the symmetric part of
  // the result. In 2D, compute and store qw * c / det(J).
  // in[0] is coefficients with shape [ncomp=1, Q]
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  CurlCurlContext *bc = (CurlCurlContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
  {
    case 221:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        qd[i] = qw[i] * c[i] / DetJ22(J + i, Q);
      }
      break;
    case 321:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        qd[i] = qw[i] * c[i] / DetJ32(J + i, Q);
      }
      break;
    case 333:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ33(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for a curl-curl operator with a vector
// coefficient evaluated at quadrature points.
CEED_QFUNCTION(f_build_curlcurl_quad_vector)(void *ctx, CeedInt Q,
                                             const CeedScalar *const *in,
                                             CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) J^T C J and store the symmetric part of
  // the result. In 2D, compute and store qw * c / det(J).
  // in[0] is coefficients with shape [ncomp=space_dim, Q]
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  CurlCurlContext *bc = (CurlCurlContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
  {
    case 333:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ33(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for a curl-curl operator
// with a matrix coefficient evaluated at quadrature points
CEED_QFUNCTION(f_build_curlcurl_quad_matrix)(void *ctx, CeedInt Q,
                                             const CeedScalar *const *in,
                                             CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) J^T C J and store the symmetric part of
  // the result. In 2D, compute and store qw * c / det(J).
  // in[0] is coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  CurlCurlContext *bc = (CurlCurlContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (100 * bc->space_dim + 10 * bc->dim + bc->curl_dim)
  {
    case 333:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ33(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for applying a curl-curl operator.
CEED_QFUNCTION(f_apply_curlcurl)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                 CeedScalar *const *out)
{
  // in[0], out[0] have shape [curl_dim, ncomp=1, Q]
  CurlCurlContext *bc = (CurlCurlContext *)ctx;
  const CeedScalar *uc = in[0], *qd = in[1];
  CeedScalar *vc = out[0];
  switch (10 * bc->dim + bc->curl_dim)
  {
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        vc[i] = qd[i] * uc[i];
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        const CeedScalar uc0 = uc[i + Q * 0];
        const CeedScalar uc1 = uc[i + Q * 1];
        const CeedScalar uc2 = uc[i + Q * 2];
        vc[i + Q * 0] = qd[i + Q * 0] * uc0 + qd[i + Q * 1] * uc1 + qd[i + Q * 2] * uc2;
        vc[i + Q * 1] = qd[i + Q * 1] * uc0 + qd[i + Q * 3] * uc1 + qd[i + Q * 4] * uc2;
        vc[i + Q * 2] = qd[i + Q * 2] * uc0 + qd[i + Q * 4] * uc1 + qd[i + Q * 5] * uc2;
      }
      break;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_CURLCURL_QF_H
