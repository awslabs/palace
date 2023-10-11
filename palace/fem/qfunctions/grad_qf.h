// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_GRAD_QF_H
#define PALACE_LIBCEED_GRAD_QF_H

#include "utils_qf.h"

struct GradContext
{
  CeedInt dim, space_dim;
  CeedScalar coeff;
};

// libCEED QFunction for building quadrature data for a gradient operator with a scalar
// constant coefficient.
CEED_QFUNCTION(f_build_grad_const_scalar)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  // At every quadrature point, compute qw C adj(J)^T and store the result.
  // in[0] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[1] is quadrature weights, size (Q)
  GradContext *bc = (GradContext *)ctx;
  const CeedScalar coeff = bc->coeff;
  const CeedScalar *J = in[0], *qw = in[1];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 11:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        qd[i] = qw[i] * coeff;
      }
      break;
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt21(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt22(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt32(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt33(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for a gradient operator with a scalar
// coefficient evaluated at quadrature points.
CEED_QFUNCTION(f_build_grad_quad_scalar)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  // At every quadrature point, compute qw C adj(J)^T and store the result.
  // in[0] is coefficients, size (Q)
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  GradContext *bc = (GradContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 11:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        qd[i] = qw[i] * c[i];
      }
      break;
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt21(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt22(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt32(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt33(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for a gradient operator with a vector
// coefficient evaluated at quadrature points.
CEED_QFUNCTION(f_build_grad_quad_vector)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  // At every quadrature point, compute qw C adj(J)^T and store the result.
  // in[0] is coefficients with shape [ncomp=vdim, Q]
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  GradContext *bc = (GradContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt21(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt22(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt32(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt33(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for a gradient operator with a matrix
// coefficient evaluated at quadrature points.
CEED_QFUNCTION(f_build_grad_quad_matrix)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                         CeedScalar *const *out)
{
  // At every quadrature point, compute qw C adj(J)^T and store the result.
  // in[0] is coefficients with shape [ncomp=vdim*(vdim+1)/2, Q]
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  GradContext *bc = (GradContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt21(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt22(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt32(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultCAdjJt33(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for applying a gradient operator.
CEED_QFUNCTION(f_apply_grad)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                             CeedScalar *const *out)
{
  // in[0] has shape [dim, ncomp=1, Q]
  // out[0] has shape [ncomp=space_dim, Q]
  GradContext *bc = (GradContext *)ctx;
  const CeedScalar *ug = in[0], *qd = in[1];
  CeedScalar *v = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 11:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        v[i] = qd[i] * ug[i];
      }
      break;
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        const CeedScalar ug0 = ug[i + Q * 0];
        v[i + Q * 0] = qd[i + Q * 0] * ug0;
        v[i + Q * 1] = qd[i + Q * 1] * ug0;
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        const CeedScalar ug0 = ug[i + Q * 0];
        const CeedScalar ug1 = ug[i + Q * 1];
        v[i + Q * 0] = qd[i + Q * 0] * ug0 + qd[i + Q * 2] * ug1;
        v[i + Q * 1] = qd[i + Q * 1] * ug0 + qd[i + Q * 3] * ug1;
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        const CeedScalar ug0 = ug[i + Q * 0];
        const CeedScalar ug1 = ug[i + Q * 1];
        v[i + Q * 0] = qd[i + Q * 0] * ug0 + qd[i + Q * 3] * ug1;
        v[i + Q * 1] = qd[i + Q * 1] * ug0 + qd[i + Q * 4] * ug1;
        v[i + Q * 2] = qd[i + Q * 2] * ug0 + qd[i + Q * 5] * ug1;
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        const CeedScalar ug0 = ug[i + Q * 0];
        const CeedScalar ug1 = ug[i + Q * 1];
        const CeedScalar ug2 = ug[i + Q * 2];
        v[i + Q * 0] = qd[i + Q * 0] * ug0 + qd[i + Q * 3] * ug1 + qd[i + Q * 6] * ug2;
        v[i + Q * 1] = qd[i + Q * 1] * ug0 + qd[i + Q * 4] * ug1 + qd[i + Q * 7] * ug2;
        v[i + Q * 2] = qd[i + Q * 2] * ug0 + qd[i + Q * 5] * ug1 + qd[i + Q * 8] * ug2;
      }
      break;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_GRAD_QF_H
