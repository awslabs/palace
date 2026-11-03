// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SHAPE_21_QF_H
#define PALACE_LIBCEED_SHAPE_21_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_2_qf.h"
#include "utils_21_qf.h"

// Derivative of (Q u, w) on a curve: S = Q u w P, with P the tangential projection, so
// the output is Q u w times the weighted |J| J (JᵀJ)⁻¹.
CEED_QFUNCTION(f_shape_mass_21)(void *__restrict__ ctx, CeedInt Q,
                                const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *w = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar s =
        CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]) * wdetJ[i] * u[i] * w[i];
    v[i + Q * 0] = s * adjJt[i + Q * 0];
    v[i + Q * 1] = s * adjJt[i + Q * 1];
  }
  return 0;
}

// Derivative of (Q u, w) for tangential traces of H(curl) fields on a curve, with unit
// normal n and P = I - n nᵀ:
// S = (w ⋅ Qu) P - w (PQu)ᵀ - u (PQᵀw)ᵀ + n (n ⋅ Qu) wᵀ + n (n ⋅ Qᵀw) uᵀ.
CEED_QFUNCTION(f_shape_hcurl_21)(void *__restrict__ ctx, CeedInt Q,
                                 const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *w = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[4], A[2];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    MatUnpack21(adjJt + i, Q, A);
    const CeedScalar pu[2] = {A[0] * u[i], A[1] * u[i]};
    const CeedScalar pw[2] = {A[0] * w[i], A[1] * w[i]};
    const CeedScalar qpu[2] = {coeff[0] * pu[0] + coeff[2] * pu[1],
                               coeff[1] * pu[0] + coeff[3] * pu[1]};
    const CeedScalar qtpw[2] = {coeff[0] * pw[0] + coeff[1] * pw[1],
                                coeff[2] * pw[0] + coeff[3] * pw[1]};
    const CeedScalar inv_norm = 1.0 / sqrt(A[0] * A[0] + A[1] * A[1]);
    const CeedScalar n[2] = {-A[1] * inv_norm, A[0] * inv_norm};
    const CeedScalar s = pw[0] * qpu[0] + pw[1] * qpu[1];
    const CeedScalar n_qpu = n[0] * qpu[0] + n[1] * qpu[1];
    const CeedScalar n_qtpw = n[0] * qtpw[0] + n[1] * qtpw[1];
    const CeedScalar a = A[0] * qpu[0] + A[1] * qpu[1];
    const CeedScalar b = A[0] * qtpw[0] + A[1] * qtpw[1];
    const CeedScalar gu = A[0] * pu[0] + A[1] * pu[1];
    const CeedScalar gw = A[0] * pw[0] + A[1] * pw[1];
    for (CeedInt r = 0; r < 2; r++)
    {
      v[i + Q * r] =
          wdetJ[i] * (s * A[r] - pw[r] * a - pu[r] * b + n[r] * (n_qpu * gw + n_qtpw * gu));
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_SHAPE_21_QF_H
