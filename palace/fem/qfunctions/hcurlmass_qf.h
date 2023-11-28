// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_MASS_QF_H
#define PALACE_LIBCEED_HCURL_MASS_QF_H

#include "utils_qf.h"

// libCEED QFunctions for H(curl) + H1 mass operators (Piola transformation u =
// adj(J)^T / det(J) ̂u and u = ̂u).
// in[0] is Jacobian determinant quadrature data, shape [Q]
// in[1] is transpose adjugate Jacobian quadrature data, shape [ncomp=space_dim*dim, Q]
// in[2] is active vector, shape [ncomp=1, Q]
// in[3] is active vector gradient, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [ncomp=1, Q]
// out[1] is active vector gradient, shape [qcomp=dim, ncomp=1, Q]

CEED_QFUNCTION(f_apply_hcurlmass_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2], *gradu = in[3];
  CeedScalar *v = out[0], *gradv = out[1];
  const CeedScalar coeff_mass = 1.0;
  const CeedScalar coeff[3] = {1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = coeff_mass * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[3];
    MultAtBA22(adjJt + i, Q, coeff, qd);

    const CeedScalar gradu0 = gradu[i + Q * 0];
    const CeedScalar gradu1 = gradu[i + Q * 1];
    gradv[i + Q * 0] = wdetJ[i] * (qd[0] * gradu0 + qd[1] * gradu1);
    gradv[i + Q * 1] = wdetJ[i] * (qd[1] * gradu0 + qd[2] * gradu1);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlmass_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2], *gradu = in[3];
  CeedScalar *v = out[0], *gradv = out[1];
  const CeedScalar coeff_mass = 1.0;
  const CeedScalar coeff[6] = {1.0, 0.0, 0.0,
                               1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = coeff_mass * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[6];
    MultAtBA33(adjJt + i, Q, coeff, qd);

    const CeedScalar gradu0 = gradu[i + Q * 0];
    const CeedScalar gradu1 = gradu[i + Q * 1];
    const CeedScalar gradu2 = gradu[i + Q * 2];
    gradv[i + Q * 0] = wdetJ[i] * (qd[0] * gradu0 + qd[1] * gradu1 + qd[2] * gradu2);
    gradv[i + Q * 1] = wdetJ[i] * (qd[1] * gradu0 + qd[3] * gradu1 + qd[4] * gradu2);
    gradv[i + Q * 2] = wdetJ[i] * (qd[2] * gradu0 + qd[4] * gradu1 + qd[5] * gradu2);
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlmass_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2], *gradu = in[3];
  CeedScalar *v = out[0], *gradv = out[1];
  const CeedScalar coeff_mass = 1.0;
  const CeedScalar coeff[3] = {1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = coeff_mass * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[1];
    MultAtBA21(adjJt + i, Q, coeff, qd);

    const CeedScalar gradu0 = gradu[i + Q * 0];
    gradv[i + Q * 0] = wdetJ[i] * qd[0] * gradu0;
  }
  return 0;
}

CEED_QFUNCTION(f_apply_hcurlmass_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *u = in[2], *gradu = in[3];
  CeedScalar *v = out[0], *gradv = out[1];
  const CeedScalar coeff_mass = 1.0;
  const CeedScalar coeff[6] = {1.0, 0.0, 0.0,
                               1.0, 0.0, 1.0};  // XX TODO NON-IDENTITY COEFFICIENTS
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v[i] = coeff_mass * wdetJ[i] * u[i];
  }
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar qd[3];
    MultAtBA32(adjJt + i, Q, coeff, qd);

    const CeedScalar gradu0 = gradu[i + Q * 0];
    const CeedScalar gradu1 = gradu[i + Q * 1];
    gradv[i + Q * 0] = wdetJ[i] * (qd[0] * gradu0 + qd[1] * gradu1);
    gradv[i + Q * 1] = wdetJ[i] * (qd[1] * gradu0 + qd[2] * gradu1);
  }
  return 0;
}

// XX TODO REMOVE AND COEFFICIENTS

// struct DiffusionMassContext
// {
//   CeedInt dim, space_dim;
// };

// // libCEED QFunction for building quadrature data for a diffusion + mass operator with a
// // scalar coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_diff_mass_quad_scalar)(void *ctx, CeedInt Q,
//                                               const CeedScalar *const *in,
//                                               CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) adj(J) C adj(J)^T and qw * c * det(J)
//   // and store the result.
//   // in[0] is diffusion coefficients with shape [ncomp=1, Q]
//   // in[1] is mass coefficients with shape [ncomp=1, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   DiffusionMassContext *bc = (DiffusionMassContext *)ctx;
//   const CeedScalar *cd = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdd = out[0], *qdm = out[0] + Q * bc->dim * (bc->dim + 1) / 2;
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 11:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdd[i] = qw[i] * cd[i] / J[i];
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * J[i];
//       }
//       break;
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt21(J + i, Q, cd + i, Q, 1, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ21(J + i, Q);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt22(J + i, Q, cd + i, Q, 1, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ22(J + i, Q);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt32(J + i, Q, cd + i, Q, 1, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ32(J + i, Q);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cd + i, Q, 1, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ33(J + i, Q);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a diffusion + mass operator with a
// // vector coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_diff_mass_quad_vector)(void *ctx, CeedInt Q,
//                                               const CeedScalar *const *in,
//                                               CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) adj(J) C adj(J)^T and qw * c * det(J)
//   // and store the result.
//   // in[0] is diffusion coefficients with shape [ncomp=space_dim, Q]
//   // in[1] is mass coefficients with shape [ncomp=1, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   DiffusionMassContext *bc = (DiffusionMassContext *)ctx;
//   const CeedScalar *cd = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdd = out[0], *qdm = out[0] + Q * bc->dim * (bc->dim + 1) / 2;
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt21(J + i, Q, cd + i, Q, 2, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ21(J + i, Q);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt22(J + i, Q, cd + i, Q, 2, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ22(J + i, Q);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt32(J + i, Q, cd + i, Q, 3, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ32(J + i, Q);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cd + i, Q, 3, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ33(J + i, Q);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for building quadrature data for a diffusion + mass operator with a
// // matrix coefficient evaluated at quadrature points.
// CEED_QFUNCTION(f_build_diff_mass_quad_matrix)(void *ctx, CeedInt Q,
//                                               const CeedScalar *const *in,
//                                               CeedScalar *const *out)
// {
//   // At every quadrature point, compute qw / det(J) adj(J) C adj(J)^T and qw * c * det(J)
//   // and store the result.
//   // in[0] is diffusion coefficients with shape [ncomp=space_dim*(space_dim+1)/2, Q]
//   // in[1] is mass coefficients with shape [ncomp=1, Q]
//   // in[2] is Jacobians with shape [dim, ncomp=space_dim, Q]
//   // in[3] is quadrature weights, size (Q)
//   DiffusionMassContext *bc = (DiffusionMassContext *)ctx;
//   const CeedScalar *cd = in[0], *cm = in[1], *J = in[2], *qw = in[3];
//   CeedScalar *qdd = out[0], *qdm = out[0] + Q * bc->dim * (bc->dim + 1) / 2;
//   switch (10 * bc->space_dim + bc->dim)
//   {
//     case 21:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt21(J + i, Q, cd + i, Q, 3, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ21(J + i, Q);
//       }
//       break;
//     case 22:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt22(J + i, Q, cd + i, Q, 3, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ22(J + i, Q);
//       }
//       break;
//     case 32:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt32(J + i, Q, cd + i, Q, 6, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ32(J + i, Q);
//       }
//       break;
//     case 33:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         MultAdjJCAdjJt33(J + i, Q, cd + i, Q, 6, qw[i], Q, qdd + i);
//       }
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         qdm[i] = qw[i] * cm[i] * DetJ33(J + i, Q);
//       }
//       break;
//   }
//   return 0;
// }

// // libCEED QFunction for applying a diffusion + mass operator.
// CEED_QFUNCTION(f_apply_diff_mass)(void *ctx, CeedInt Q, const CeedScalar *const *in,
//                                   CeedScalar *const *out)
// {
//   // in[0], out[0] have shape [ncomp=1, Q]
//   // in[1], out[1] have shape [dim, ncomp=1, Q]
//   DiffusionMassContext *bc = (DiffusionMassContext *)ctx;
//   const CeedScalar *u = in[0], *ug = in[1], *qdd = in[2],
//                    *qdm = in[2] + Q * bc->dim * (bc->dim + 1) / 2;
//   CeedScalar *v = out[0], *vg = out[1];
//   switch (bc->dim)
//   {
//     case 1:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         vg[i] = qdd[i] * ug[i];
//       }
//       break;
//     case 2:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar ug0 = ug[i + Q * 0];
//         const CeedScalar ug1 = ug[i + Q * 1];
//         vg[i + Q * 0] = qdd[i + Q * 0] * ug0 + qdd[i + Q * 1] * ug1;
//         vg[i + Q * 1] = qdd[i + Q * 1] * ug0 + qdd[i + Q * 2] * ug1;
//       }
//       break;
//     case 3:
//       CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//       {
//         const CeedScalar ug0 = ug[i + Q * 0];
//         const CeedScalar ug1 = ug[i + Q * 1];
//         const CeedScalar ug2 = ug[i + Q * 2];
//         vg[i + Q * 0] = qdd[i + Q * 0] * ug0 + qdd[i + Q * 1] * ug1 + qdd[i + Q * 2] *
//         ug2; vg[i + Q * 1] = qdd[i + Q * 1] * ug0 + qdd[i + Q * 3] * ug1 + qdd[i + Q * 4]
//         * ug2; vg[i + Q * 2] = qdd[i + Q * 2] * ug0 + qdd[i + Q * 4] * ug1 + qdd[i + Q *
//         5] * ug2;
//       }
//       break;
//   }
//   CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
//   {
//     v[i] = qdm[i] * u[i];
//   }
//   return 0;
// }

#endif  // PALACE_LIBCEED_HCURL_MASS_QF_H
