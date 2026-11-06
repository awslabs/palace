// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SHAPE_32_QF_H
#define PALACE_LIBCEED_SHAPE_32_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_3_qf.h"
#include "utils_32_qf.h"

// Derivative of (Q u, w) on a surface: S = Q u w P, with P the tangential projection, so
// the output is Q u w times the weighted |J| J (JᵀJ)⁻¹.
CEED_QFUNCTION(f_shape_mass_32)(void *__restrict__ ctx, CeedInt Q,
                                const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *w = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar s =
        CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]) * wdetJ[i] * u[i] * w[i];
    for (CeedInt k = 0; k < 6; k++)
    {
      v[i + Q * k] = s * adjJt[i + Q * k];
    }
  }
  return 0;
}

// Unit normal of a surface element from the columns of adj(J)ᵀ / |J|.
CEED_QFUNCTION_HELPER void UnitNormal32(const CeedScalar A[6], CeedScalar n[3])
{
  n[0] = A[1] * A[5] - A[2] * A[4];
  n[1] = A[2] * A[3] - A[0] * A[5];
  n[2] = A[0] * A[4] - A[1] * A[3];
  const CeedScalar inv_norm = 1.0 / sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  n[0] *= inv_norm;
  n[1] *= inv_norm;
  n[2] *= inv_norm;
}

// Derivative of (Q u, w) for tangential traces of H(curl) fields on a surface, with unit
// normal n and P = I - n nᵀ:
// S = (w ⋅ Qu) P - w (PQu)ᵀ - u (PQᵀw)ᵀ + n (n ⋅ Qu) wᵀ + n (n ⋅ Qᵀw) uᵀ.
CEED_QFUNCTION(f_shape_hcurl_32)(void *__restrict__ ctx, CeedInt Q,
                                 const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *w = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[9], A[6], pu[3], pw[3], qpu[3], qtpw[3], n[3];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack32(adjJt + i, Q, A);
    for (CeedInt r = 0; r < 3; r++)
    {
      pu[r] = A[r] * u[i + Q * 0] + A[r + 3] * u[i + Q * 1];
      pw[r] = A[r] * w[i + Q * 0] + A[r + 3] * w[i + Q * 1];
    }
    for (CeedInt r = 0; r < 3; r++)
    {
      qpu[r] = coeff[r] * pu[0] + coeff[r + 3] * pu[1] + coeff[r + 6] * pu[2];
      qtpw[r] = coeff[3 * r] * pw[0] + coeff[3 * r + 1] * pw[1] + coeff[3 * r + 2] * pw[2];
    }
    UnitNormal32(A, n);
    const CeedScalar s = pw[0] * qpu[0] + pw[1] * qpu[1] + pw[2] * qpu[2];
    const CeedScalar n_qpu = n[0] * qpu[0] + n[1] * qpu[1] + n[2] * qpu[2];
    const CeedScalar n_qtpw = n[0] * qtpw[0] + n[1] * qtpw[1] + n[2] * qtpw[2];

    for (CeedInt c = 0; c < 2; c++)
    {
      const CeedScalar *col = A + 3 * c;
      const CeedScalar a = col[0] * qpu[0] + col[1] * qpu[1] + col[2] * qpu[2];
      const CeedScalar b = col[0] * qtpw[0] + col[1] * qtpw[1] + col[2] * qtpw[2];
      const CeedScalar gu = col[0] * pu[0] + col[1] * pu[1] + col[2] * pu[2];
      const CeedScalar gw = col[0] * pw[0] + col[1] * pw[1] + col[2] * pw[2];
      for (CeedInt r = 0; r < 3; r++)
      {
        v[i + Q * (r + 3 * c)] = wdetJ[i] * (s * col[r] - pw[r] * a - pu[r] * b +
                                             n[r] * (n_qpu * gw + n_qtpw * gu));
      }
    }
  }
  return 0;
}

// Derivative of (Q u, w) for scalar fields mapped by 1 / |J| (curls of H(curl) fields):
// S = -(w Q u) P. in[1] is the quadrature weight, and u and w are in[2] and in[3].
CEED_QFUNCTION(f_shape_l2_32)(void *__restrict__ ctx, CeedInt Q,
                              const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *qw = in[1],
                   *u = in[2], *w = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar s = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]) *
                         qw[i] * qw[i] / wdetJ[i] * u[i] * w[i];
    for (CeedInt k = 0; k < 6; k++)
    {
      v[i + Q * k] = -s * adjJt[i + Q * k];
    }
  }
  return 0;
}

// Derivative of the linear form (f, w) with f fixed, for tangential traces of H(curl)
// fields on a surface: S = (f ⋅ w) P - w (Pf)ᵀ + n (n ⋅ f) wᵀ. The vector f is the first
// column of the coefficient, and w is in[1].
CEED_QFUNCTION(f_shape_hcurl_lf_32)(void *__restrict__ ctx, CeedInt Q,
                                    const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *w = in[1];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[9], A[6], pw[3], n[3];
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack32(adjJt + i, Q, A);
    for (CeedInt r = 0; r < 3; r++)
    {
      pw[r] = A[r] * w[i + Q * 0] + A[r + 3] * w[i + Q * 1];
    }
    UnitNormal32(A, n);
    const CeedScalar *f = coeff;
    const CeedScalar s = f[0] * pw[0] + f[1] * pw[1] + f[2] * pw[2];
    const CeedScalar n_f = n[0] * f[0] + n[1] * f[1] + n[2] * f[2];

    for (CeedInt c = 0; c < 2; c++)
    {
      const CeedScalar *col = A + 3 * c;
      const CeedScalar a = col[0] * f[0] + col[1] * f[1] + col[2] * f[2];
      const CeedScalar gw = col[0] * pw[0] + col[1] * pw[1] + col[2] * pw[2];
      for (CeedInt r = 0; r < 3; r++)
      {
        v[i + Q * (r + 3 * c)] = wdetJ[i] * (s * col[r] - pw[r] * a + n[r] * n_f * gw);
      }
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_SHAPE_32_QF_H
