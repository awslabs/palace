// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_MESHSENS_33_QF_H
#define PALACE_LIBCEED_MESHSENS_33_QF_H

// Mesh sensitivity QFunction for 3D: computes df/dJ at each quadrature point, where f is
// the CurlCurlMass integrand (stiffness - lambda * mass) for an H(curl) eigenvalue problem.
//
// The differentiation backend is selected at compile time:
//   PALACE_WITH_ENZYME: single reverse-mode pass via __enzyme_autodiff
//   Otherwise: 9 forward-mode passes using a scalar Dual number type
//
// Both produce exact analytic derivatives — no finite differences.

#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

// ============================================================================
// Forward-mode AD: Dual number type and templated utilities
// ============================================================================

struct MeshSensDual
{
  CeedScalar val;
  CeedScalar dot;
};

CEED_QFUNCTION_HELPER MeshSensDual MakeDual(CeedScalar v, CeedScalar d)
{
  return {v, d};
}

CEED_QFUNCTION_HELPER MeshSensDual MakeDual(CeedScalar v)
{
  return {v, 0.0};
}

CEED_QFUNCTION_HELPER MeshSensDual operator+(MeshSensDual a, MeshSensDual b)
{
  return {a.val + b.val, a.dot + b.dot};
}

CEED_QFUNCTION_HELPER MeshSensDual operator-(MeshSensDual a, MeshSensDual b)
{
  return {a.val - b.val, a.dot - b.dot};
}

CEED_QFUNCTION_HELPER MeshSensDual operator*(MeshSensDual a, MeshSensDual b)
{
  return {a.val * b.val, a.val * b.dot + a.dot * b.val};
}

CEED_QFUNCTION_HELPER MeshSensDual operator/(MeshSensDual a, MeshSensDual b)
{
  return {a.val / b.val, (a.dot * b.val - a.val * b.dot) / (b.val * b.val)};
}

// Promote CeedScalar to Dual for mixed operations.
CEED_QFUNCTION_HELPER MeshSensDual operator*(CeedScalar s, MeshSensDual a)
{
  return {s * a.val, s * a.dot};
}

CEED_QFUNCTION_HELPER MeshSensDual operator*(MeshSensDual a, CeedScalar s)
{
  return {a.val * s, a.dot * s};
}

CEED_QFUNCTION_HELPER MeshSensDual operator+(CeedScalar s, MeshSensDual a)
{
  return {s + a.val, a.dot};
}

CEED_QFUNCTION_HELPER MeshSensDual operator+(MeshSensDual a, CeedScalar s)
{
  return {a.val + s, a.dot};
}

CEED_QFUNCTION_HELPER MeshSensDual operator-(CeedScalar s, MeshSensDual a)
{
  return {s - a.val, -a.dot};
}

CEED_QFUNCTION_HELPER MeshSensDual operator-(MeshSensDual a, CeedScalar s)
{
  return {a.val - s, a.dot};
}

CEED_QFUNCTION_HELPER MeshSensDual operator/(MeshSensDual a, CeedScalar s)
{
  return {a.val / s, a.dot / s};
}

CEED_QFUNCTION_HELPER MeshSensDual operator/(CeedScalar s, MeshSensDual a)
{
  return {s / a.val, -s * a.dot / (a.val * a.val)};
}

// ============================================================================
// Templated versions of Palace QFunction utility functions.
// These mirror the implementations in utils_33_qf.h but work with both
// CeedScalar and MeshSensDual.
// ============================================================================

template <typename T>
CEED_QFUNCTION_HELPER T MeshSensDetJ33(const T J[9])
{
  return J[0] * (J[4] * J[8] - J[5] * J[7]) - J[1] * (J[3] * J[8] - J[5] * J[6]) +
         J[2] * (J[3] * J[7] - J[4] * J[6]);
}

template <typename T>
CEED_QFUNCTION_HELPER T MeshSensAdjJt33(const T J[9], T adjJt[9])
{
  // Compute adj(J)^T and return det(J).
  adjJt[0] = J[4] * J[8] - J[7] * J[5];
  adjJt[3] = J[7] * J[2] - J[1] * J[8];
  adjJt[6] = J[1] * J[5] - J[4] * J[2];
  adjJt[1] = J[6] * J[5] - J[3] * J[8];
  adjJt[4] = J[0] * J[8] - J[6] * J[2];
  adjJt[7] = J[3] * J[2] - J[0] * J[5];
  adjJt[2] = J[3] * J[7] - J[6] * J[4];
  adjJt[5] = J[6] * J[1] - J[0] * J[7];
  adjJt[8] = J[0] * J[4] - J[3] * J[1];
  return J[0] * adjJt[0] + J[1] * adjJt[1] + J[2] * adjJt[2];
}

// y = A^T * B * C * x, where A, B, C are 3x3 (column-major) and x, y are 3-vectors.
// B is a CeedScalar matrix (material coefficient, not differentiated).
template <typename T>
CEED_QFUNCTION_HELPER void MeshSensMultAtBCx33(const T A[9], const CeedScalar B[9],
                                                const T C[9], const CeedScalar x[3], T y[3])
{
  // Step 1: z = C * x
  T z[3];
  z[0] = C[0] * x[0] + C[3] * x[1] + C[6] * x[2];
  z[1] = C[1] * x[0] + C[4] * x[1] + C[7] * x[2];
  z[2] = C[2] * x[0] + C[5] * x[1] + C[8] * x[2];

  // Step 2: w = B * z
  T w[3];
  w[0] = B[0] * z[0] + B[3] * z[1] + B[6] * z[2];
  w[1] = B[1] * z[0] + B[4] * z[1] + B[7] * z[2];
  w[2] = B[2] * z[0] + B[5] * z[1] + B[8] * z[2];

  // Step 3: y = A^T * w
  y[0] = A[0] * w[0] + A[1] * w[1] + A[2] * w[2];
  y[1] = A[3] * w[0] + A[4] * w[1] + A[5] * w[2];
  y[2] = A[6] * w[0] + A[7] * w[1] + A[8] * w[2];
}

// ============================================================================
// Integrands: templated on scalar type T.
// These match the math in f_apply_hdivmass_33 exactly.
// ============================================================================

// Stiffness (curl-curl) integrand at one quadrature point:
//   f_K = w * detJ * curl_u^T * J_eff^T * mu_inv * J_eff * curl_u
// where J_eff = adj(adjJt_norm)^T = J/detJ (the H(div) Piola map).
//
// Matches hdivmass_33_qf.h stiffness block: AdjJt33(adjJt_loc, J_loc) then
// MultAtBCx33(J_loc, coeff, J_loc, curlu, curlv) scaled by wdetJ.
template <typename T>
CEED_QFUNCTION_HELPER T MeshSensStiffnessIntegrand(const T *J, const CeedScalar *curl_u,
                                                    const CeedScalar *mu_inv, CeedScalar w)
{
  T adjJt[9];
  const T detJ = MeshSensAdjJt33(J, adjJt);

  // Normalize: adjJt_norm = adjJt / detJ (as stored in Palace geometry data).
  T adjJt_norm[9];
  for (int k = 0; k < 9; k++)
  {
    adjJt_norm[k] = adjJt[k] / detJ;
  }

  // Recover J/detJ via adjugate of adjJt_norm: adj(adjJt/detJ)^T = J/detJ.
  T J_eff[9];
  MeshSensAdjJt33(adjJt_norm, J_eff);

  // v = J_eff^T * mu_inv * J_eff * curl_u
  T v[3];
  MeshSensMultAtBCx33(J_eff, mu_inv, J_eff, curl_u, v);

  // f = w * detJ * dot(curl_u, v)
  return w * detJ * (curl_u[0] * v[0] + curl_u[1] * v[1] + curl_u[2] * v[2]);
}

// Mass integrand at one quadrature point (H(curl) Piola form):
//   f_M = w * detJ * u^T * adjJt_norm^T * eps * adjJt_norm * u
// where adjJt_norm = adj(J)^T / detJ.
//
// Matches hdivmass_33_qf.h mass block: MultAtBCx33(adjJt_loc, coeff, adjJt_loc, u, v)
// scaled by wdetJ.
template <typename T>
CEED_QFUNCTION_HELPER T MeshSensMassIntegrand(const T *J, const CeedScalar *u,
                                               const CeedScalar *eps, CeedScalar w)
{
  T adjJt[9];
  const T detJ = MeshSensAdjJt33(J, adjJt);

  // Normalize: adjJt_norm = adjJt / detJ.
  T adjJt_norm[9];
  for (int k = 0; k < 9; k++)
  {
    adjJt_norm[k] = adjJt[k] / detJ;
  }

  // v = adjJt_norm^T * eps * adjJt_norm * u
  T v[3];
  MeshSensMultAtBCx33(adjJt_norm, eps, adjJt_norm, u, v);

  // f = w * detJ * dot(u, v)
  return w * detJ * (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]);
}

// ============================================================================
// Differentiation: compute df/dJ[9] at one quadrature point.
// Enzyme (reverse-mode) or dual-number (forward-mode) selected at compile time.
// ============================================================================

#ifdef PALACE_WITH_ENZYME

// Enzyme intrinsic declarations.
extern int enzyme_dup;
extern int enzyme_const;
extern double __enzyme_autodiff(void *, ...);

// Enzyme wrappers: thin shims calling the CeedScalar-specialized integrands.
static double MeshSensStiffnessForEnzyme(double *J, const double *curl_u,
                                          const double *mu_inv, double w)
{
  return MeshSensStiffnessIntegrand<CeedScalar>(J, curl_u, mu_inv, w);
}

static double MeshSensMassForEnzyme(double *J, const double *u, const double *eps, double w)
{
  return MeshSensMassIntegrand<CeedScalar>(J, u, eps, w);
}

CEED_QFUNCTION_HELPER void MeshSensDifferentiate(const CeedScalar J_in[9],
                                                  const CeedScalar *curl_u,
                                                  const CeedScalar *mu_inv,
                                                  const CeedScalar *u,
                                                  const CeedScalar *eps, CeedScalar w,
                                                  CeedScalar lambda, CeedScalar df_dJ[9])
{
  CeedScalar J[9], dK[9] = {0}, dM[9] = {0};
  for (int k = 0; k < 9; k++)
  {
    J[k] = J_in[k];
  }

  __enzyme_autodiff((void *)MeshSensStiffnessForEnzyme, enzyme_dup, J, dK, enzyme_const,
                    curl_u, enzyme_const, mu_inv, enzyme_const, w);

  for (int k = 0; k < 9; k++)
  {
    J[k] = J_in[k];
  }

  __enzyme_autodiff((void *)MeshSensMassForEnzyme, enzyme_dup, J, dM, enzyme_const, u,
                    enzyme_const, eps, enzyme_const, w);

  for (int k = 0; k < 9; k++)
  {
    df_dJ[k] = dK[k] - lambda * dM[k];
  }
}

#else  // !PALACE_WITH_ENZYME

// Forward-mode AD via dual numbers: 9 passes, seeding one J component per pass.
CEED_QFUNCTION_HELPER void MeshSensDifferentiate(const CeedScalar J_in[9],
                                                  const CeedScalar *curl_u,
                                                  const CeedScalar *mu_inv,
                                                  const CeedScalar *u,
                                                  const CeedScalar *eps, CeedScalar w,
                                                  CeedScalar lambda, CeedScalar df_dJ[9])
{
  for (int k = 0; k < 9; k++)
  {
    MeshSensDual J_d[9];
    for (int j = 0; j < 9; j++)
    {
      J_d[j] = MakeDual(J_in[j], (j == k) ? 1.0 : 0.0);
    }
    MeshSensDual f_K = MeshSensStiffnessIntegrand(J_d, curl_u, mu_inv, w);
    MeshSensDual f_M = MeshSensMassIntegrand(J_d, u, eps, w);
    df_dJ[k] = f_K.dot - lambda * f_M.dot;
  }
}

#endif  // PALACE_WITH_ENZYME

// ============================================================================
// CeedQFunction: mesh sensitivity for combined stiffness - lambda*mass.
// ============================================================================

// Context layout: same as CurlCurlMass coefficient pair context, followed by lambda.
// [coeff_pair_ctx..., lambda(CeedScalar)]
//
// Inputs:
//   in[0]: element attributes, shape [Q]              (CEED_EVAL_INTERP, passive)
//   in[1]: quadrature weights, shape [Q]              (CEED_EVAL_WEIGHT, passive)
//   in[2]: Jacobian J, shape [9*Q]                    (CEED_EVAL_GRAD, active input)
//   in[3]: eigenvector interp u, shape [3*Q]          (CEED_EVAL_INTERP, passive)
//   in[4]: eigenvector curl curl_u, shape [3*Q]       (CEED_EVAL_CURL, passive)
// Outputs:
//   out[0]: df/dJ sensitivity, shape [9*Q]            (CEED_EVAL_GRAD, active output)

CEED_QFUNCTION(f_apply_meshsens_33)(void *__restrict__ ctx, CeedInt Q,
                                    const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *qw = in[1], *J = in[2], *u = in[3], *curlu = in[4];
  CeedScalar *__restrict__ df_dJ = out[0];

  // The context holds the coefficient pair (mass first, stiffness second) followed by
  // lambda as a CeedScalar at the end.
  const CeedIntScalar *coeff_ctx = (const CeedIntScalar *)ctx;

  // Find lambda at the end of the pair context. Each half has its own size
  // (2 + NumAttr + 9*NumMat), so walk through both to find the appended lambda.
  const CeedIntScalar *stiff_ctx = CoeffPairSecond<3>(coeff_ctx);
  const CeedInt stiff_size = 2 + NumAttr(stiff_ctx) + 9 * NumMat(stiff_ctx);
  const CeedScalar lambda = stiff_ctx[stiff_size].second;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    // Unpack Jacobian at this quadrature point (strided, stride = Q).
    CeedScalar J_loc[9];
    MatUnpack33(J + i, Q, J_loc);

    // Unpack eigenvector values at this quadrature point.
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    const CeedScalar curlu_loc[3] = {curlu[i + Q * 0], curlu[i + Q * 1], curlu[i + Q * 2]};

    // Look up material coefficients from element attribute.
    CeedScalar mu_inv[9], eps[9];
    CoeffUnpack3(CoeffPairSecond<3>(coeff_ctx), (CeedInt)attr[i], mu_inv);
    CoeffUnpack3(coeff_ctx, (CeedInt)attr[i], eps);

    // Differentiate: df/dJ = d(f_K)/dJ - lambda * d(f_M)/dJ.
    CeedScalar df_dJ_loc[9];
    MeshSensDifferentiate(J_loc, curlu_loc, mu_inv, u_loc, eps, qw[i], lambda, df_dJ_loc);

    // Store to strided output.
    df_dJ[i + Q * 0] = df_dJ_loc[0];
    df_dJ[i + Q * 1] = df_dJ_loc[1];
    df_dJ[i + Q * 2] = df_dJ_loc[2];
    df_dJ[i + Q * 3] = df_dJ_loc[3];
    df_dJ[i + Q * 4] = df_dJ_loc[4];
    df_dJ[i + Q * 5] = df_dJ_loc[5];
    df_dJ[i + Q * 6] = df_dJ_loc[6];
    df_dJ[i + Q * 7] = df_dJ_loc[7];
    df_dJ[i + Q * 8] = df_dJ_loc[8];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_MESHSENS_33_QF_H
