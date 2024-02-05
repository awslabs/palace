// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFF_3_QF_H
#define PALACE_LIBCEED_COEFF_3_QF_H

#include "coeff_qf.h"

CEED_QFUNCTION_HELPER void CoeffUnpack3(const CeedIntScalar *ctx, const CeedInt attr,
                                        CeedScalar coeff[6])
{
  const CeedInt k = AttrMat(ctx)[attr - 1].first;
  const CeedIntScalar *mat_coeff = MatCoeff(ctx);
  coeff[0] = mat_coeff[6 * k + 0].second;
  coeff[1] = mat_coeff[6 * k + 1].second;
  coeff[2] = mat_coeff[6 * k + 2].second;
  coeff[3] = mat_coeff[6 * k + 3].second;
  coeff[4] = mat_coeff[6 * k + 4].second;
  coeff[5] = mat_coeff[6 * k + 5].second;
}

#endif  // PALACE_LIBCEED_COEFF_3_QF_H
