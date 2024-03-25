// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFF_1_QF_H
#define PALACE_LIBCEED_COEFF_1_QF_H

#include "coeff_qf.h"

CEED_QFUNCTION_HELPER CeedScalar CoeffUnpack1(const CeedIntScalar *ctx, const CeedInt attr)
{
  const CeedInt k = (NumAttr(ctx) > 0) ? AttrMat(ctx)[attr - 1].first : 0;
  return MatCoeff(ctx)[k].second;
}

CEED_QFUNCTION_HELPER void CoeffUnpack1(const CeedIntScalar *ctx, const CeedInt attr,
                                        CeedScalar coeff[1])
{
  coeff[0] = CoeffUnpack1(ctx, attr);
}

#endif  // PALACE_LIBCEED_COEFF_1_QF_H
