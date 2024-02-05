// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFF_2_QF_H
#define PALACE_LIBCEED_COEFF_2_QF_H

#include "coeff_qf.h"

CEED_QFUNCTION_HELPER void CoeffUnpack2(const CeedIntScalar *ctx, const CeedInt attr,
                                        CeedScalar coeff[3])
{
  const CeedInt k = AttrMat(ctx)[attr - 1].first;
  const CeedIntScalar *mat_coeff = MatCoeff(ctx);
  coeff[0] = mat_coeff[3 * k + 0].second;
  coeff[1] = mat_coeff[3 * k + 1].second;
  coeff[2] = mat_coeff[3 * k + 2].second;
}

#endif  // PALACE_LIBCEED_COEFF_2_QF_H
