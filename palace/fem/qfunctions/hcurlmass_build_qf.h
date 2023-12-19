// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_MASS_BUILD_QF_H
#define PALACE_LIBCEED_HCURL_MASS_BUILD_QF_H

#include "coeff_qf.h"
#include "utils_qf.h"

// Build functions replace active vector output with quadrature point data, stored as a
// symmetric matrix, and remove active vector input.

CEED_QFUNCTION(f_build_hcurlmass_22)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q;
  CeedScalar *__restrict__ qd1 = out[0], *__restrict__ qd2 = out[0] + Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

      qd1[i + Q * 0] = coeff * wdetJ[i];
    }
    {
      CeedScalar coeff[3], adjJt_loc[4], qd_loc[3];
      CoeffUnpack2(CoeffPairSecond1((const CeedIntScalar *)ctx), (CeedInt)attr[i], coeff);
      MatUnpack22(adjJt + i, Q, adjJt_loc);
      MultAtBA22(adjJt_loc, coeff, qd_loc);

      qd2[i + Q * 0] = wdetJ[i] * qd_loc[0];
      qd2[i + Q * 1] = wdetJ[i] * qd_loc[1];
      qd2[i + Q * 2] = wdetJ[i] * qd_loc[2];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_hcurlmass_33)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q;
  CeedScalar *__restrict__ qd1 = out[0], *__restrict__ qd2 = out[0] + Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

      qd1[i + Q * 0] = coeff * wdetJ[i];
    }
    {
      CeedScalar coeff[6], adjJt_loc[9], qd_loc[6];
      CoeffUnpack3(CoeffPairSecond1((const CeedIntScalar *)ctx), (CeedInt)attr[i], coeff);
      MatUnpack33(adjJt + i, Q, adjJt_loc);
      MultAtBA33(adjJt_loc, coeff, qd_loc);

      qd2[i + Q * 0] = wdetJ[i] * qd_loc[0];
      qd2[i + Q * 1] = wdetJ[i] * qd_loc[1];
      qd2[i + Q * 2] = wdetJ[i] * qd_loc[2];
      qd2[i + Q * 3] = wdetJ[i] * qd_loc[3];
      qd2[i + Q * 4] = wdetJ[i] * qd_loc[4];
      qd2[i + Q * 5] = wdetJ[i] * qd_loc[5];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_hcurlmass_21)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q;
  CeedScalar *__restrict__ qd1 = out[0], *__restrict__ qd2 = out[0] + Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

      qd1[i + Q * 0] = coeff * wdetJ[i];
    }
    {
      CeedScalar coeff[3], adjJt_loc[2], qd_loc[1];
      CoeffUnpack2(CoeffPairSecond1((const CeedIntScalar *)ctx), (CeedInt)attr[i], coeff);
      MatUnpack21(adjJt + i, Q, adjJt_loc);
      MultAtBA21(adjJt_loc, coeff, qd_loc);

      qd2[i + Q * 0] = wdetJ[i] * qd_loc[0];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_hcurlmass_32)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                     CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q;
  CeedScalar *__restrict__ qd1 = out[0], *__restrict__ qd2 = out[0] + Q;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    {
      const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)attr[i]);

      qd1[i + Q * 0] = coeff * wdetJ[i];
    }
    {
      CeedScalar coeff[6], adjJt_loc[6], qd_loc[3];
      CoeffUnpack3(CoeffPairSecond1((const CeedIntScalar *)ctx), (CeedInt)attr[i], coeff);
      MatUnpack32(adjJt + i, Q, adjJt_loc);
      MultAtBA32(adjJt_loc, coeff, qd_loc);

      qd2[i + Q * 0] = wdetJ[i] * qd_loc[0];
      qd2[i + Q * 1] = wdetJ[i] * qd_loc[1];
      qd2[i + Q * 2] = wdetJ[i] * qd_loc[2];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_MASS_BUILD_QF_H
