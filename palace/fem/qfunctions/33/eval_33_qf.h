// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_EVAL_33_QF_H
#define PALACE_LIBCEED_EVAL_33_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

// QFunctions for pointwise evaluation of derived field quantities at arbitrary points
// of 3D volume elements (for example the nodal points of an interpolatory output space
// for visualization). No quadrature weighting is applied; the element geometry is
// computed on the fly from the mesh nodes gradient (as in geom_qf.h) instead of stored
// geometry factor data. Contributions of repeated applications accumulate at the output
// L-vector (CeedOperatorApplyAdd), which implements the real/imaginary part sums of the
// quadratic quantities.
//
// Inputs: in[0] is element attributes, shape [Q]
//         in[1] is mesh node Jacobians, shape [qcomp=dim, ncomp=space_dim, Q]
//         in[2], in[3] are field inputs (reference components at the points)
// The context is a CeedIntScalar array: [0].second = scale, followed by a material
// property coefficient context. H(curl): E = adj(J)ᵀ/detJ u, H(div): B = J/detJ u.

// Electric energy density: v = 0.5 * scale * (eps_r E) . E.
CEED_QFUNCTION(f_eval_energy_e_33)(void *__restrict__ ctx_, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *J = in[1], *u = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar J_loc[9], adjJt_loc[9], E[3], eps[9], D[3];
    MatUnpack33(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt33<true>(J_loc, adjJt_loc);
    MultAx33(adjJt_loc, u_loc, E);
    CoeffUnpack3(ctx + 1, (CeedInt)attr[i], eps);
    MultAx33(eps, E, D);
    v[i] = 0.5 * ctx[0].second * (D[0] * E[0] + D[1] * E[1] + D[2] * E[2]) / (detJ * detJ);
  }
  return 0;
}

// Magnetic energy density: v = 0.5 * scale * (mu⁻¹ B) . B.
CEED_QFUNCTION(f_eval_energy_m_33)(void *__restrict__ ctx_, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *J = in[1], *u = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar J_loc[9], adjJt_loc[9], B[3], invmu[9], H[3];
    MatUnpack33(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt33<true>(J_loc, adjJt_loc);
    MultAx33(J_loc, u_loc, B);
    CoeffUnpack3(ctx + 1, (CeedInt)attr[i], invmu);
    MultAx33(invmu, B, H);
    v[i] = 0.5 * ctx[0].second * (H[0] * B[0] + H[1] * B[1] + H[2] * B[2]) / (detJ * detJ);
  }
  return 0;
}

// Poynting vector: v = scale * E x (mu⁻¹ B).
CEED_QFUNCTION(f_eval_poynting_33)(void *__restrict__ ctx_, CeedInt Q,
                                   const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedIntScalar *ctx = (const CeedIntScalar *)ctx_;
  const CeedScalar *attr = in[0], *J = in[1], *u_e = in[2], *u_b = in[3];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar ue_loc[3] = {u_e[i + Q * 0], u_e[i + Q * 1], u_e[i + Q * 2]};
    const CeedScalar ub_loc[3] = {u_b[i + Q * 0], u_b[i + Q * 1], u_b[i + Q * 2]};
    CeedScalar J_loc[9], adjJt_loc[9], E[3], B[3], invmu[9], H[3];
    MatUnpack33(J + i, Q, J_loc);
    const CeedScalar detJ = AdjJt33<true>(J_loc, adjJt_loc);
    MultAx33(adjJt_loc, ue_loc, E);
    MultAx33(J_loc, ub_loc, B);
    CoeffUnpack3(ctx + 1, (CeedInt)attr[i], invmu);
    MultAx33(invmu, B, H);
    const CeedScalar s = ctx[0].second / (detJ * detJ);
    v[i + Q * 0] = s * (E[1] * H[2] - E[2] * H[1]);
    v[i + Q * 1] = s * (E[2] * H[0] - E[0] * H[2]);
    v[i + Q * 2] = s * (E[0] * H[1] - E[1] * H[0]);
  }
  return 0;
}

#endif  // PALACE_LIBCEED_EVAL_33_QF_H
