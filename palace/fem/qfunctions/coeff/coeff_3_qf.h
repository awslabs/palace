// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFF_3_QF_H
#define PALACE_LIBCEED_COEFF_3_QF_H

#include "coeff_qf.h"

CEED_QFUNCTION_HELPER void CoeffUnpack3(const CeedIntScalar *ctx, const CeedInt attr,
                                        CeedScalar coeff[9])
{
  const CeedInt k = (NumAttr(ctx) > 0) ? AttrMat(ctx)[attr - 1].first : 0;
  const CeedIntScalar *mat_coeff = MatCoeff(ctx);
  coeff[0] = mat_coeff[9 * k + 0].second;
  coeff[1] = mat_coeff[9 * k + 1].second;
  coeff[2] = mat_coeff[9 * k + 2].second;
  coeff[3] = mat_coeff[9 * k + 3].second;
  coeff[4] = mat_coeff[9 * k + 4].second;
  coeff[5] = mat_coeff[9 * k + 5].second;
  coeff[6] = mat_coeff[9 * k + 6].second;
  coeff[7] = mat_coeff[9 * k + 7].second;
  coeff[8] = mat_coeff[9 * k + 8].second;
}

#endif  // PALACE_LIBCEED_COEFF_3_QF_H
