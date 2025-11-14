// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_MASS_33_QF_H
#define PALACE_LIBCEED_HCURL_MASS_33_QF_H

#include "../coeff/coeff_1_qf.h"
#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

CEED_QFUNCTION(f_apply_hcurlmass_33)(void *__restrict__ ctx, CeedInt Q,
                                     const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *qdata = in[0], *u = in[1], *gradu = in[2];
  CeedScalar *__restrict__ v = out[0], *__restrict__ gradv = out[1];

  const CeedInt stride = 2 + 9; // attr, w * |J|, (adjJt / |J|) colwise
  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    {
      const CeedScalar coeff = CoeffUnpack1((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0]);

      v[i] = coeff * qdata_i[1] * u[i];
    }
    const CeedScalar* adjJt_loc = qdata_i + 2;
    {
      const CeedScalar u_loc[3] = {gradu[i + Q * 0], gradu[i + Q * 1], gradu[i + Q * 2]};
      CeedScalar coeff[9], v_loc[3];
      CoeffUnpack3(CoeffPairSecond<1>((const CeedIntScalar *)ctx), (CeedInt)qdata_i[0], coeff);
      MultAtBCx33(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

      gradv[i + Q * 0] = qdata_i[1] * v_loc[0];
      gradv[i + Q * 1] = qdata_i[1] * v_loc[1];
      gradv[i + Q * 2] = qdata_i[1] * v_loc[2];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_MASS_33_QF_H
