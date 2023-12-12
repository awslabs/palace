// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFF_QF_H
#define PALACE_LIBCEED_COEFF_QF_H

union CeedIntScalar
{
  CeedInt first;
  CeedScalar second;
};

// The first entry of ctx is the number of (0-based) attributes, followed by the entries of
// the attribute to material index array  (these are also 0-based).
// The next entry is the number of material property coefficients, followed by the
// coefficients.
// Pair coefficients are two coefficient contexts arranged contiguously in memory.

CEED_QFUNCTION_HELPER const CeedIntScalar *AttrMat(const CeedIntScalar *ctx)
{
  return ctx + 1;
}

CEED_QFUNCTION_HELPER const CeedIntScalar *MatCoeff(const CeedIntScalar *ctx)
{
  const CeedInt num_attr = ctx[0].first;
  return ctx + 2 + num_attr;
}

CEED_QFUNCTION_HELPER CeedScalar CoeffUnpack1(const CeedIntScalar *ctx, const CeedInt attr)
{
  const CeedInt k = AttrMat(ctx)[attr].first;
  return MatCoeff(ctx)[k].second;
}

CEED_QFUNCTION_HELPER void CoeffUnpack1(const CeedIntScalar *ctx, const CeedInt attr,
                                        CeedScalar coeff[1])
{
  coeff[0] = CoeffUnpack1(ctx, attr);
}

CEED_QFUNCTION_HELPER void CoeffUnpack2(const CeedIntScalar *ctx, const CeedInt attr,
                                        CeedScalar coeff[3])
{
  const CeedInt k = AttrMat(ctx)[attr].first;
  const CeedIntScalar *mat_coeff = MatCoeff(ctx);
  coeff[0] = mat_coeff[3 * k + 0].second;
  coeff[1] = mat_coeff[3 * k + 1].second;
  coeff[2] = mat_coeff[3 * k + 2].second;
}

CEED_QFUNCTION_HELPER void CoeffUnpack3(const CeedIntScalar *ctx, const CeedInt attr,
                                        CeedScalar coeff[6])
{
  const CeedInt k = AttrMat(ctx)[attr].first;
  const CeedIntScalar *mat_coeff = MatCoeff(ctx);
  coeff[0] = mat_coeff[6 * k + 0].second;
  coeff[1] = mat_coeff[6 * k + 1].second;
  coeff[2] = mat_coeff[6 * k + 2].second;
  coeff[3] = mat_coeff[6 * k + 3].second;
  coeff[4] = mat_coeff[6 * k + 4].second;
  coeff[5] = mat_coeff[6 * k + 5].second;
}

CEED_QFUNCTION_HELPER const CeedIntScalar *CoeffPairSecond1(const CeedIntScalar *ctx)
{
  const CeedInt num_attr = ctx[0].first;
  const CeedInt num_mat = ctx[1 + num_attr].first;
  return ctx + 2 + num_attr + num_mat;
}

CEED_QFUNCTION_HELPER const CeedIntScalar *CoeffPairSecond2(const CeedIntScalar *ctx)
{
  const CeedInt num_attr = ctx[0].first;
  const CeedInt num_mat = ctx[1 + num_attr].first;
  return ctx + 2 + num_attr + 3 * num_mat;
}

CEED_QFUNCTION_HELPER const CeedIntScalar *CoeffPairSecond3(const CeedIntScalar *ctx)
{
  const CeedInt num_attr = ctx[0].first;
  const CeedInt num_mat = ctx[1 + num_attr].first;
  return ctx + 2 + num_attr + 6 * num_mat;
}

#endif  // PALACE_LIBCEED_COEFF_QF_H
