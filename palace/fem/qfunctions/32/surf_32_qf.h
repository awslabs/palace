// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SURF_32_QF_H
#define PALACE_LIBCEED_SURF_32_QF_H

#include "../33/utils_33_qf.h"
#include "utils_32_qf.h"

// QFunctions for surface output functionals (integrals of functions of fields over
// boundary elements of a 3D mesh).
//
// Inputs are ordered as follows (geometry data first, then fields):
//  - in[0]: face geometry data for the boundary element, shape [8, Q]:
//           {attr, w * detJ_face, adj(J_face)^T / detJ_face (3x2, column-major)}.
//           J_face is the 3x2 Jacobian of the boundary element mapping. The surface
//           normal is computed as the cross product of the columns of
//           adj(J_face)^T / detJ_face, which is parallel (and identically oriented) to
//           the cross product of the columns of J_face. This matches the orientation
//           convention of mfem::CalcOrtho for boundary elements (outward normal for
//           exterior boundaries, pointing out of element 1 in general).
//  - in[1]: volume geometry data for the element attached to the boundary element,
//           evaluated at the (mapped) face quadrature points, shape [11, Q]:
//           {attr, w * detJ_vol (unused), adj(J_vol)^T / detJ_vol (3x3, column-major)}.
//           Used for the Piola transformations of fields from the volume element's
//           reference space: H(curl): E_phys = adj(J)^T/detJ * u_ref, H(div):
//           B_phys = J/detJ * u_ref with J/detJ recovered via AdjJt33 applied to the
//           stored adj(J)^T/detJ.
//  - in[2], in[3], ...: field inputs evaluated at face quadrature points (reference
//           components from the volume element basis tabulation).
//
// The output is the integrand value multiplied by the face quadrature weight, summed
// over quadrature points by the all-ones output element basis, yielding one (or more)
// value(s) per boundary element.

// Compute y = A x for a 3x3 matrix A stored column-major.
CEED_QFUNCTION_HELPER void MultAx33(const CeedScalar A[9], const CeedScalar x[3],
                                    CeedScalar y[3])
{
  y[0] = A[0] * x[0] + A[3] * x[1] + A[6] * x[2];
  y[1] = A[1] * x[0] + A[4] * x[1] + A[7] * x[2];
  y[2] = A[2] * x[0] + A[5] * x[1] + A[8] * x[2];
}

// Compute the (non-unit) surface normal as the cross product of the columns of the
// stored adj(J_face)^T / detJ_face (3x2, column-major), and normalize it.
CEED_QFUNCTION_HELPER void SurfaceNormal32(const CeedScalar adjJt_face[6], CeedScalar n[3])
{
  n[0] = adjJt_face[1] * adjJt_face[5] - adjJt_face[2] * adjJt_face[4];
  n[1] = adjJt_face[2] * adjJt_face[3] - adjJt_face[0] * adjJt_face[5];
  n[2] = adjJt_face[0] * adjJt_face[4] - adjJt_face[1] * adjJt_face[3];
  const CeedScalar norm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  n[0] /= norm;
  n[1] /= norm;
  n[2] /= norm;
}

// Surface area: v = w * detJ_face. No field inputs.
CEED_QFUNCTION(f_integ_surf_area_32)(void *, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0] + Q;
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = wdetJ[i];
  }
  return 0;
}

// Squared L2 norm of an H(curl) field over the surface: v = w * detJ_face * |E|^2 with
// E = adj(J_vol)^T / detJ_vol * u_ref.
CEED_QFUNCTION(f_integ_surf_hcurl_norm2_32)(void *, CeedInt Q, const CeedScalar *const *in,
                                            CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0] + Q, *adjJt_vol = in[1] + 2 * Q, *u = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar adjJt_loc[9], E[3];
    MatUnpack33(adjJt_vol + i, Q, adjJt_loc);
    MultAx33(adjJt_loc, u_loc, E);

    v[i] = wdetJ[i] * (E[0] * E[0] + E[1] * E[1] + E[2] * E[2]);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_SURF_32_QF_H
