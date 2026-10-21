// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_FIELD_22_QF_H
#define PALACE_LIBCEED_HCURL_FIELD_22_QF_H

#include "../coeff/coeff_2_qf.h"
#include "utils_22_qf.h"

// H(curl) mass / H1 gradient operator with a spatially-varying scalar material coefficient
// sampled at the quadrature points. in[2] is the coefficient field, shape [1, Q]; the
// coefficient from ctx (typically an identity/sign) is scaled by it.

CEED_QFUNCTION(f_apply_hcurl_field_22)(void *__restrict__ ctx, CeedInt Q,
                                       const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *adjJt = in[0] + 2 * Q, *u = in[1],
                   *coeff_data = in[2];
  CeedScalar *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar u_loc[2] = {u[i + Q * 0], u[i + Q * 1]};
    const CeedScalar c = coeff_data[i];
    CeedScalar coeff[4], adjJt_loc[4], v_loc[2];
    CoeffUnpack2((const CeedIntScalar *)ctx, (CeedInt)attr[i], coeff);
    coeff[0] *= c;
    coeff[1] *= c;
    coeff[2] *= c;
    coeff[3] *= c;
    MatUnpack22(adjJt + i, Q, adjJt_loc);
    MultAtBCx22(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_FIELD_22_QF_H
