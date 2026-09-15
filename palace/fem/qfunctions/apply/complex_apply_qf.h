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

// Hermitian transpose of f_apply_complex_3. The forward kernel applies the complex
// 3x3 tensor A with A[r,c] = qr[r+3c] + i qi[r+3c]. The adjoint applies A^H, whose
// entry (r,c) is conj(A[c,r]): the real part transposes the index (qr[c+3r]) and the
// imaginary part transposes AND negates (-qi[c+3r]). Reusing the forward action would
// be correct only for a symmetric tensor; the index transpose and sign flip are
// required for the general (nonsymmetric) case.
CEED_QFUNCTION(f_apply_complex_3_transpose)(void *, CeedInt Q, const CeedScalar *const *in,
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
        // A^H[r,c] = conj(A[c,r]); A[c,r] has index c + 3*r.
        const CeedScalar ar = qr[q + (c + 3 * r) * Q];
        const CeedScalar ai = -qi[q + (c + 3 * r) * Q];
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

// Hermitian transpose of f_apply_complex_33. The curl term uses a real tensor qd that
// acts independently on the real and imaginary blocks, so its adjoint is the plain
// transpose (index k + 3*r instead of r + 3*k) with no conjugation. The complex mass
// block is handled by f_apply_complex_3_transpose.
CEED_QFUNCTION(f_apply_complex_33_transpose)(void *ctx, CeedInt Q,
                                             const CeedScalar *const *in,
                                             CeedScalar *const *out)
{
  f_apply_complex_3_transpose(ctx, Q, in, out);
  const CeedScalar *__restrict__ qd = in[0] + 9 * Q, *__restrict__ u = in[3];
  CeedScalar *__restrict__ v = out[1];
  CeedPragmaSIMD for (CeedInt q = 0; q < Q; q++)
  {
    for (CeedInt r = 0; r < 3; r++)
    {
      for (CeedInt c = 0; c < 2; c++)
      {
        // Transposed real curl tensor: qd^T[r,k] = qd[k + 3*r].
        v[q + (2 * r + c) * Q] = qd[q + (3 * r + 0) * Q] * u[q + (0 + c) * Q] +
                                 qd[q + (3 * r + 1) * Q] * u[q + (2 + c) * Q] +
                                 qd[q + (3 * r + 2) * Q] * u[q + (4 + c) * Q];
      }
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_COMPLEX_APPLY_QF_H
