// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_MASS_BUILD_QF_H
#define PALACE_LIBCEED_L2_MASS_BUILD_QF_H

#include "types_qf.h"
#include "utils_geom_qf.h"
#include "utils_qf.h"

// Build functions replace active vector output with quadrature point data, stored as a
// symmetric matrix, and remove active vector input.

CEED_QFUNCTION(f_build_l2mass_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const MatCoeffPairContext21 *bc = (const MatCoeffPairContext21 *)ctx;
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *qw = in[2], *attr = in[3];
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[3], adjJt_loc[4], J_loc[4], qd_loc[3];
    CoeffUnpack(&bc->first, (CeedInt)attr[i], coeff);
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    AdjJt22<false>(adjJt_loc, J_loc);
    MultAtBA22(J_loc, coeff, qd_loc);

    qd[i + Q * 0] = wdetJ[i] * qd_loc[0];
    qd[i + Q * 1] = wdetJ[i] * qd_loc[1];
    qd[i + Q * 2] = wdetJ[i] * qd_loc[2];
  }
  qd += 3 * Q;
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack(&bc->second, (CeedInt)attr[i]);

    qd[i] = coeff * qw[i] * qw[i] / wdetJ[i];
  }
  return 0;
}

CEED_QFUNCTION(f_build_l2mass_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const MatCoeffPairContext31 *bc = (const MatCoeffPairContext31 *)ctx;
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *qw = in[2], *attr = in[3];
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[6], adjJt_loc[9], J_loc[9], qd_loc[6];
    CoeffUnpack(&bc->first, (CeedInt)attr[i], coeff);
    MatUnpack33(adjJt + i, Q, adjJt_loc);
    AdjJt33<false>(adjJt_loc, J_loc);
    MultAtBA33(J_loc, coeff, qd_loc);

    qd[i + Q * 0] = wdetJ[i] * qd_loc[0];
    qd[i + Q * 1] = wdetJ[i] * qd_loc[1];
    qd[i + Q * 2] = wdetJ[i] * qd_loc[2];
    qd[i + Q * 3] = wdetJ[i] * qd_loc[3];
    qd[i + Q * 4] = wdetJ[i] * qd_loc[4];
    qd[i + Q * 5] = wdetJ[i] * qd_loc[5];
  }
  qd += 6 * Q;
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack(&bc->second, (CeedInt)attr[i]);

    qd[i] = coeff * qw[i] * qw[i] / wdetJ[i];
  }
  return 0;
}

CEED_QFUNCTION(f_build_l2mass_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const MatCoeffPairContext21 *bc = (const MatCoeffPairContext21 *)ctx;
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *qw = in[2], *attr = in[3];
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[3], adjJt_loc[2], J_loc[2], qd_loc[1];
    CoeffUnpack(&bc->first, (CeedInt)attr[i], coeff);
    MatUnpack21(adjJt + i, Q, adjJt_loc);
    AdjJt21<false>(adjJt_loc, J_loc);
    MultAtBA21(J_loc, coeff, qd_loc);

    qd[i + Q * 0] = wdetJ[i] * qd_loc[0];
  }
  qd += Q;
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack(&bc->second, (CeedInt)attr[i]);

    qd[i] = coeff * qw[i] * qw[i] / wdetJ[i];
  }
  return 0;
}

CEED_QFUNCTION(f_build_l2mass_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  const MatCoeffPairContext31 *bc = (const MatCoeffPairContext31 *)ctx;
  const CeedScalar *wdetJ = in[0], *adjJt = in[1], *qw = in[2], *attr = in[3];
  CeedScalar *qd = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar coeff[6], adjJt_loc[6], J_loc[6], qd_loc[3];
    CoeffUnpack(&bc->first, (CeedInt)attr[i], coeff);
    MatUnpack32(adjJt + i, Q, adjJt_loc);
    AdjJt32<false>(adjJt_loc, J_loc);
    MultAtBA32(J_loc, coeff, qd_loc);

    qd[i + Q * 0] = wdetJ[i] * qd_loc[0];
    qd[i + Q * 1] = wdetJ[i] * qd_loc[1];
    qd[i + Q * 2] = wdetJ[i] * qd_loc[2];
  }
  qd += 3 * Q;
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar coeff = CoeffUnpack(&bc->second, (CeedInt)attr[i]);

    qd[i] = coeff * qw[i] * qw[i] / wdetJ[i];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_L2_MASS_BUILD_QF_H
