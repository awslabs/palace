// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COMPLEX_APPLY_QF_H
#define PALACE_LIBCEED_COMPLEX_APPLY_QF_H

// H(curl)/H(div) fields have layout [physical component][real/imag component][Q].
// Retain separate passive QData fields, including their original restrictions, so
// coefficient updates and backend-specific QData layouts remain valid.
CEED_QFUNCTION(f_apply_complex_3)(void *, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const CeedScalar *__restrict__ qr = in[0], *__restrict__ qi = in[1],
                                 *__restrict__ u = in[2];
  CeedScalar *__restrict__ v = out[0];
  CeedPragmaSIMD for (CeedInt q = 0; q < Q; q++)
  {
    for (CeedInt r = 0; r < 3; r++)
    {
      CeedScalar vr = 0.0, vi = 0.0;
      for (CeedInt c = 0; c < 3; c++)
      {
        const CeedScalar ar = qr[q + (r + 3 * c) * Q];
        const CeedScalar ai = qi[q + (r + 3 * c) * Q];
        const CeedScalar ur = u[q + (2 * c) * Q];
        const CeedScalar ui = u[q + (2 * c + 1) * Q];
        vr += ar * ur - ai * ui;
        vi += ai * ur + ar * ui;
      }
      v[q + (2 * r) * Q] = vr;
      v[q + (2 * r + 1) * Q] = vi;
    }
  }
  return 0;
}

// Real f_apply_33: interpolation tensor first, curl tensor second.
// Imaginary f_apply_3: interpolation tensor only.
CEED_QFUNCTION(f_apply_complex_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                   CeedScalar *const *out)
{
  f_apply_complex_3(ctx, Q, in, out);
  const CeedScalar *__restrict__ qd = in[0] + 9 * Q, *__restrict__ u = in[3];
  CeedScalar *__restrict__ v = out[1];
  CeedPragmaSIMD for (CeedInt q = 0; q < Q; q++)
  {
    for (CeedInt r = 0; r < 3; r++)
    {
      for (CeedInt c = 0; c < 2; c++)
      {
        v[q + (2 * r + c) * Q] = qd[q + (r + 0) * Q] * u[q + (0 + c) * Q] +
                                 qd[q + (r + 3) * Q] * u[q + (2 + c) * Q] +
                                 qd[q + (r + 6) * Q] * u[q + (4 + c) * Q];
      }
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_COMPLEX_APPLY_QF_H
