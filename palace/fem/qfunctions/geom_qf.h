// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_GEOM_QF_H
#define PALACE_LIBCEED_GEOM_QF_H

#include <math.h>

struct GeomFactorContext
{
  bool compute_wdetJ, compute_adjJt, compute_J;
};

CEED_QFUNCTION_HELPER CeedScalar DetJ22(const CeedScalar *J, const CeedInt J_stride)
{
  // J: 0 2
  //    1 3
  return J[J_stride * 0] * J[J_stride * 3] - J[J_stride * 1] * J[J_stride * 2];
}

CEED_QFUNCTION_HELPER CeedScalar DetJ33(const CeedScalar *J, const CeedInt J_stride)
{
  // J: 0 3 6
  //    1 4 7
  //    2 5 8
  return J[J_stride * 0] *
             (J[J_stride * 4] * J[J_stride * 8] - J[J_stride * 5] * J[J_stride * 7]) -
         J[J_stride * 1] *
             (J[J_stride * 3] * J[J_stride * 8] - J[J_stride * 5] * J[J_stride * 6]) +
         J[J_stride * 2] *
             (J[J_stride * 3] * J[J_stride * 7] - J[J_stride * 4] * J[J_stride * 6]);
}

CEED_QFUNCTION_HELPER CeedScalar DetJ21(const CeedScalar *J, const CeedInt J_stride)
{
  // J: 0
  //    1
  return sqrt(J[J_stride * 0] * J[J_stride * 0] + J[J_stride * 1] * J[J_stride * 1]);
}

CEED_QFUNCTION_HELPER CeedScalar DetJ32(const CeedScalar *J, const CeedInt J_stride)
{
  // J: 0 3
  //    1 4
  //    2 5
  const CeedScalar E = J[J_stride * 0] * J[J_stride * 0] +
                       J[J_stride * 1] * J[J_stride * 1] +
                       J[J_stride * 2] * J[J_stride * 2];
  const CeedScalar G = J[J_stride * 3] * J[J_stride * 3] +
                       J[J_stride * 4] * J[J_stride * 4] +
                       J[J_stride * 5] * J[J_stride * 5];
  const CeedScalar F = J[J_stride * 0] * J[J_stride * 3] +
                       J[J_stride * 1] * J[J_stride * 4] +
                       J[J_stride * 2] * J[J_stride * 5];
  return sqrt(E * G - F * F);
}

CEED_QFUNCTION_HELPER void AdjJt22(const CeedScalar *J, const CeedInt J_stride,
                                   const CeedInt qd_stride, CeedScalar *qd)
{
  // Compute adj(J)^T / det(J) and store the result.
  // J: 0 2   adj(J):  J22 -J12   qd: 0 2
  //    1 3           -J21  J11       1 3
  const CeedScalar J11 = J[J_stride * 0];
  const CeedScalar J21 = J[J_stride * 1];
  const CeedScalar J12 = J[J_stride * 2];
  const CeedScalar J22 = J[J_stride * 3];
  const CeedScalar detJ = J11 * J22 - J21 * J12;
  qd[qd_stride * 0] = J22 / detJ;
  qd[qd_stride * 1] = -J12 / detJ;
  qd[qd_stride * 2] = -J21 / detJ;
  qd[qd_stride * 3] = J11 / detJ;
}

CEED_QFUNCTION_HELPER void AdjJt33(const CeedScalar *J, const CeedInt J_stride,
                                   const CeedInt qd_stride, CeedScalar *qd)
{
  // Compute adj(J)^T / det(J) and store the result.
  // J: 0 3 6   qd: 0 3 6
  //    1 4 7       1 4 7
  //    2 5 8       2 5 8
  const CeedScalar J11 = J[J_stride * 0];
  const CeedScalar J21 = J[J_stride * 1];
  const CeedScalar J31 = J[J_stride * 2];
  const CeedScalar J12 = J[J_stride * 3];
  const CeedScalar J22 = J[J_stride * 4];
  const CeedScalar J32 = J[J_stride * 5];
  const CeedScalar J13 = J[J_stride * 6];
  const CeedScalar J23 = J[J_stride * 7];
  const CeedScalar J33 = J[J_stride * 8];
  const CeedScalar A11 = J22 * J33 - J23 * J32;
  const CeedScalar A21 = J23 * J31 - J21 * J33;
  const CeedScalar A31 = J21 * J32 - J22 * J31;
  const CeedScalar A12 = J13 * J32 - J12 * J33;
  const CeedScalar A22 = J11 * J33 - J13 * J31;
  const CeedScalar A32 = J12 * J31 - J11 * J32;
  const CeedScalar A13 = J12 * J23 - J13 * J22;
  const CeedScalar A23 = J13 * J21 - J11 * J23;
  const CeedScalar A33 = J11 * J22 - J12 * J21;
  const CeedScalar detJ = J11 * A11 + J21 * A12 + J31 * A13;
  qd[qd_stride * 0] = A11 / detJ;
  qd[qd_stride * 1] = A12 / detJ;
  qd[qd_stride * 2] = A13 / detJ;
  qd[qd_stride * 3] = A21 / detJ;
  qd[qd_stride * 4] = A22 / detJ;
  qd[qd_stride * 5] = A23 / detJ;
  qd[qd_stride * 6] = A31 / detJ;
  qd[qd_stride * 7] = A32 / detJ;
  qd[qd_stride * 8] = A33 / detJ;
}

CEED_QFUNCTION_HELPER void AdjJt21(const CeedScalar *J, const CeedInt J_stride,
                                   const CeedInt qd_stride, CeedScalar *qd)
{
  // Compute adj(J)^T / det(J) and store the result.
  // J: 0   adj(J): 1/sqrt(J^T J) J^T   qd: 0
  //    1                                   1
  const CeedScalar J11 = J[J_stride * 0];
  const CeedScalar J21 = J[J_stride * 1];
  const CeedScalar d = J11 * J11 + J21 * J21;
  qd[qd_stride * 0] = J11 / d;
  qd[qd_stride * 1] = J21 / d;
}

CEED_QFUNCTION_HELPER void AdjJt32(const CeedScalar *J, const CeedInt J_stride,
                                   const CeedInt qd_stride, CeedScalar *qd)
{
  // Compute adj(J)^T / det(J) and store the result.
  // J: 0 3   qd: 0 3
  //    1 4       1 4
  //    2 5       2 5
  const CeedScalar J11 = J[J_stride * 0];
  const CeedScalar J21 = J[J_stride * 1];
  const CeedScalar J31 = J[J_stride * 2];
  const CeedScalar J12 = J[J_stride * 3];
  const CeedScalar J22 = J[J_stride * 4];
  const CeedScalar J32 = J[J_stride * 5];
  const CeedScalar E = J11 * J11 + J21 * J21 + J31 * J31;
  const CeedScalar G = J12 * J12 + J22 * J22 + J32 * J32;
  const CeedScalar F = J11 * J12 + J21 * J22 + J31 * J32;
  const CeedScalar A11 = G * J11 - F * J12;
  const CeedScalar A12 = G * J21 - F * J22;
  const CeedScalar A13 = G * J31 - F * J32;
  const CeedScalar A21 = E * J12 - F * J11;
  const CeedScalar A22 = E * J22 - F * J21;
  const CeedScalar A23 = E * J32 - F * J31;
  const CeedScalar d = E * G - F * F;
  qd[qd_stride * 0] = A11 / d;
  qd[qd_stride * 1] = A12 / d;
  qd[qd_stride * 2] = A13 / d;
  qd[qd_stride * 3] = A21 / d;
  qd[qd_stride * 4] = A22 / d;
  qd[qd_stride * 5] = A23 / d;
}

CEED_QFUNCTION_HELPER void J22(const CeedScalar *J, const CeedInt J_stride,
                               const CeedInt qd_stride, CeedScalar *qd)
{
  // Compute J / det(J) and store the result.
  // J: 0 2   qd: 0 2
  //    1 3       1 3
  const CeedScalar J11 = J[J_stride * 0];
  const CeedScalar J21 = J[J_stride * 1];
  const CeedScalar J12 = J[J_stride * 2];
  const CeedScalar J22 = J[J_stride * 3];
  const CeedScalar detJ = J11 * J22 - J21 * J12;
  qd[qd_stride * 0] = J11 / detJ;
  qd[qd_stride * 1] = J21 / detJ;
  qd[qd_stride * 2] = J12 / detJ;
  qd[qd_stride * 3] = J22 / detJ;
}

CEED_QFUNCTION_HELPER void J33(const CeedScalar *J, const CeedInt J_stride,
                               const CeedInt qd_stride, CeedScalar *qd)
{
  // Compute J / det(J) and store the result.
  // J: 0 3 6   qd: 0 3 6
  //    1 4 7       1 4 7
  //    2 5 8       2 5 8
  const CeedScalar J11 = J[J_stride * 0];
  const CeedScalar J21 = J[J_stride * 1];
  const CeedScalar J31 = J[J_stride * 2];
  const CeedScalar J12 = J[J_stride * 3];
  const CeedScalar J22 = J[J_stride * 4];
  const CeedScalar J32 = J[J_stride * 5];
  const CeedScalar J13 = J[J_stride * 6];
  const CeedScalar J23 = J[J_stride * 7];
  const CeedScalar J33 = J[J_stride * 8];
  const CeedScalar detJ = J11 * (J22 * J33 - J23 * J32) + J21 * (J13 * J32 - J12 * J33) +
                          J31 * (J12 * J23 - J13 * J22);
  qd[qd_stride * 0] = J11 / detJ;
  qd[qd_stride * 1] = J21 / detJ;
  qd[qd_stride * 2] = J31 / detJ;
  qd[qd_stride * 3] = J12 / detJ;
  qd[qd_stride * 4] = J22 / detJ;
  qd[qd_stride * 5] = J32 / detJ;
  qd[qd_stride * 6] = J13 / detJ;
  qd[qd_stride * 7] = J23 / detJ;
  qd[qd_stride * 8] = J33 / detJ;
}

CEED_QFUNCTION_HELPER void J21(const CeedScalar *J, const CeedInt J_stride,
                               const CeedInt qd_stride, CeedScalar *qd)
{
  // Compute J / det(J) and store the result.
  // J: 0   qd: 0
  //    1       1
  const CeedScalar J11 = J[J_stride * 0];
  const CeedScalar J21 = J[J_stride * 1];
  const CeedScalar d = sqrt(J11 * J11 + J21 * J21);
  qd[qd_stride * 0] = J11 / d;
  qd[qd_stride * 1] = J21 / d;
}

CEED_QFUNCTION_HELPER void J32(const CeedScalar *J, const CeedInt J_stride,
                               const CeedInt qd_stride, CeedScalar *qd)
{
  // Compute J / det(J) and store the result.
  // J: 0 3   qd: 0 3
  //    1 4       1 4
  //    2 5       2 5
  const CeedScalar J11 = J[J_stride * 0];
  const CeedScalar J21 = J[J_stride * 1];
  const CeedScalar J31 = J[J_stride * 2];
  const CeedScalar J12 = J[J_stride * 3];
  const CeedScalar J22 = J[J_stride * 4];
  const CeedScalar J32 = J[J_stride * 5];
  const CeedScalar E = J11 * J11 + J21 * J21 + J31 * J31;
  const CeedScalar G = J12 * J12 + J22 * J22 + J32 * J32;
  const CeedScalar F = J11 * J12 + J21 * J22 + J31 * J32;
  const CeedScalar d = sqrt(E * G - F * F);
  qd[qd_stride * 0] = J11 / d;
  qd[qd_stride * 1] = J21 / d;
  qd[qd_stride * 2] = J31 / d;
  qd[qd_stride * 3] = J12 / d;
  qd[qd_stride * 4] = J22 / d;
  qd[qd_stride * 5] = J32 / d;
}

// libCEED QFunction for building geometry factors for integration and transformations.
// At every quadrature point, compute qw * det(J), adj(J)^T / |J|, and J / |J|, and store
// the result. in[0] is Jacobians, shape [qcomp=dim, ncomp=space_dim, Q] in[1] is quadrature
// weights, shape [Q] out[0] is Jacobian determinant quadrature data, shape [Q] out[1] is
// (transpose) adjugate Jacobian data, shape [ncomp=space_dim*dim, Q] out[2] is Jacobian
// quadrature data, shape [ncomp=space_dim*dim, Q]

CEED_QFUNCTION(f_build_geom_factor_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  GeomFactorContext *bc = (GeomFactorContext *)ctx;
  const CeedScalar *J = in[0], *qw = in[1];
  if (bc->compute_wdetJ)
  {
    CeedScalar *qd = out[0];
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      qd[i] = qw[i] * DetJ22(J + i, Q);
    }
  }
  if (bc->compute_adjJt)
  {
    CeedScalar *qd = bc->compute_wdetJ ? out[1] : out[0];
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      AdjJt22(J + i, Q, Q, qd + i);
    }
  }
  if (bc->compute_J)
  {
    CeedScalar *qd = bc->compute_J ? (bc->compute_wdetJ ? out[2] : out[1])
                                   : (bc->compute_wdetJ ? out[1] : out[0]);
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      J22(J + i, Q, Q, qd + i);
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_geom_factor_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  GeomFactorContext *bc = (GeomFactorContext *)ctx;
  const CeedScalar *J = in[0], *qw = in[1];
  if (bc->compute_wdetJ)
  {
    CeedScalar *qd = out[0];
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      qd[i] = qw[i] * DetJ33(J + i, Q);
    }
  }
  if (bc->compute_adjJt)
  {
    CeedScalar *qd = bc->compute_wdetJ ? out[1] : out[0];
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      AdjJt33(J + i, Q, Q, qd + i);
    }
  }
  if (bc->compute_J)
  {
    CeedScalar *qd = bc->compute_J ? (bc->compute_wdetJ ? out[2] : out[1])
                                   : (bc->compute_wdetJ ? out[1] : out[0]);
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      J33(J + i, Q, Q, qd + i);
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_geom_factor_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  GeomFactorContext *bc = (GeomFactorContext *)ctx;
  const CeedScalar *J = in[0], *qw = in[1];
  if (bc->compute_wdetJ)
  {
    CeedScalar *qd = out[0];
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      qd[i] = qw[i] * DetJ21(J + i, Q);
    }
  }
  if (bc->compute_adjJt)
  {
    CeedScalar *qd = bc->compute_wdetJ ? out[1] : out[0];
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      AdjJt21(J + i, Q, Q, qd + i);
    }
  }
  if (bc->compute_J)
  {
    CeedScalar *qd = bc->compute_J ? (bc->compute_wdetJ ? out[2] : out[1])
                                   : (bc->compute_wdetJ ? out[1] : out[0]);
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      J21(J + i, Q, Q, qd + i);
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_geom_factor_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                       CeedScalar *const *out)
{
  GeomFactorContext *bc = (GeomFactorContext *)ctx;
  const CeedScalar *J = in[0], *qw = in[1];
  if (bc->compute_wdetJ)
  {
    CeedScalar *qd = out[0];
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      qd[i] = qw[i] * DetJ32(J + i, Q);
    }
  }
  if (bc->compute_adjJt)
  {
    CeedScalar *qd = bc->compute_wdetJ ? out[1] : out[0];
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      AdjJt32(J + i, Q, Q, qd + i);
    }
  }
  if (bc->compute_J)
  {
    CeedScalar *qd = bc->compute_J ? (bc->compute_wdetJ ? out[2] : out[1])
                                   : (bc->compute_wdetJ ? out[1] : out[0]);
    CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
    {
      J32(J + i, Q, Q, qd + i);
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_GEOM_QF_H
