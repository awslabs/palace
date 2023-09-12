// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_VECFEMASS_QF_H
#define PALACE_LIBCEED_VECFEMASS_QF_H

#include "utils_qf.h"

struct VectorFEMassContext
{
  CeedInt dim, space_dim;
  CeedScalar coeff;
};

// libCEED QFunction for building quadrature data for an H(div) mass operator with a scalar
// constant coefficient.
CEED_QFUNCTION(f_build_hdivmass_const_scalar)(void *ctx, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) J^T C J and store the symmetric part of
  // the result.
  // in[0] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[1] is quadrature weights, size (Q)
  VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
  const CeedScalar coeff = bc->coeff;
  const CeedScalar *J = in[0], *qw = in[1];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 11:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        qd[i] = qw[i] * coeff * J[i];
      }
      break;
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ21(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ22(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ32(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ33(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for an H(div) mass operator with a scalar
// coefficient evaluated at quadrature points.
CEED_QFUNCTION(f_build_hdivmass_quad_scalar)(void *ctx, CeedInt Q,
                                             const CeedScalar *const *in,
                                             CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) J^T C J and store the symmetric part of
  // the result.
  // in[0] is coefficients with shape [ncomp=1, Q]
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 11:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        qd[i] = qw[i] * c[i] * J[i];
      }
      break;
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ21(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ22(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ32(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ33(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for an H(div) mass operator with a vector
// coefficient evaluated at quadrature points.
CEED_QFUNCTION(f_build_hdivmass_quad_vector)(void *ctx, CeedInt Q,
                                             const CeedScalar *const *in,
                                             CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) J^T C J and store the symmetric part of
  // the result.
  // in[0] is coefficients with shape [ncomp=space_dim, Q]
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ21(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ22(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ32(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ33(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for an H(div) mass operator with a matrix
// coefficient evaluated at quadrature points.
CEED_QFUNCTION(f_build_hdivmass_quad_matrix)(void *ctx, CeedInt Q,
                                             const CeedScalar *const *in,
                                             CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) J^T C J and store the symmetric part of
  // the result.
  // in[0] is coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ21(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ22(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ32(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultJtCJ33(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for an H(curl) mass operator with a scalar
// constant coefficient.
CEED_QFUNCTION(f_build_hcurlmass_const_scalar)(void *ctx, CeedInt Q,
                                               const CeedScalar *const *in,
                                               CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) adj(J) C adj(J)^T and store the
  // symmetric part of the result. in[0] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[1] is quadrature weights, size (Q)
  VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
  const CeedScalar coeff = bc->coeff;
  const CeedScalar *J = in[0], *qw = in[1];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 11:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        qd[i] = qw[i] * coeff / J[i];
      }
      break;
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt21(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt22(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt32(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt33(J + i, Q, &coeff, 1, 1, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for an H(curl) mass operator with a scalar
// coefficient evaluated at quadrature points.
CEED_QFUNCTION(f_build_hcurlmass_quad_scalar)(void *ctx, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) adj(J) C adj(J)^T and store the
  // symmetric part of the result.
  // in[0] is coefficients with shape [ncomp=1, Q]
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 11:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        qd[i] = qw[i] * c[i] / J[i];
      }
      break;
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt21(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt22(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt32(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt33(J + i, Q, c + i, Q, 1, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for an H(curl) mass operator with a vector
// coefficient evaluated at quadrature points.
CEED_QFUNCTION(f_build_hcurlmass_quad_vector)(void *ctx, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) adj(J) C adj(J)^T and store the
  // symmetric part of the result.
  // in[0] is coefficients with shape [ncomp=space_dim, Q]
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt21(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt22(J + i, Q, c + i, Q, 2, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt32(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt33(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for building quadrature data for an H(curl) mass operator with a matrix
// coefficient evaluated at quadrature points.
CEED_QFUNCTION(f_build_hcurlmass_quad_matrix)(void *ctx, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  // At every quadrature point, compute qw / det(J) adj(J) C adj(J)^T and store the
  // symmetric part of the result.
  // in[0] is coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
  // in[1] is Jacobians with shape [dim, ncomp=space_dim, Q]
  // in[2] is quadrature weights, size (Q)
  VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
  const CeedScalar *c = in[0], *J = in[1], *qw = in[2];
  CeedScalar *qd = out[0];
  switch (10 * bc->space_dim + bc->dim)
  {
    case 21:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt21(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
    case 22:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt22(J + i, Q, c + i, Q, 3, qw[i], Q, qd + i);
      }
      break;
    case 32:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt32(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
      }
      break;
    case 33:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        MultAdjJCAdjJt33(J + i, Q, c + i, Q, 6, qw[i], Q, qd + i);
      }
      break;
  }
  return 0;
}

// libCEED QFunction for applying a vector FE mass operator.
CEED_QFUNCTION(f_apply_vecfemass)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  // in[0], out[0] have shape [dim, ncomp=1, Q]
  VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
  const CeedScalar *u = in[0], *qd = in[1];
  CeedScalar *v = out[0];
  switch (bc->dim)
  {
    case 1:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        v[i] = qd[i] * u[i];
      }
      break;
    case 2:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        const CeedScalar u0 = u[i + Q * 0];
        const CeedScalar u1 = u[i + Q * 1];
        v[i + Q * 0] = qd[i + Q * 0] * u0 + qd[i + Q * 1] * u1;
        v[i + Q * 1] = qd[i + Q * 1] * u0 + qd[i + Q * 2] * u1;
      }
      break;
    case 3:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        const CeedScalar u0 = u[i + Q * 0];
        const CeedScalar u1 = u[i + Q * 1];
        const CeedScalar u2 = u[i + Q * 2];
        v[i + Q * 0] = qd[i + Q * 0] * u0 + qd[i + Q * 1] * u1 + qd[i + Q * 2] * u2;
        v[i + Q * 1] = qd[i + Q * 1] * u0 + qd[i + Q * 3] * u1 + qd[i + Q * 4] * u2;
        v[i + Q * 2] = qd[i + Q * 2] * u0 + qd[i + Q * 4] * u1 + qd[i + Q * 5] * u2;
      }
      break;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_VECFEMASS_QF_H
