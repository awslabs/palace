// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_TYPES_QF_H
#define PALACE_LIBCEED_TYPES_QF_H

template <int MAX_ATTR, int MAX_NUM_MAT, int DIM>
struct MatCoeffContext
{
  CeedInt attr_mat[MAX_ATTR];
  CeedScalar mat_coeff[MAX_NUM_MAT * DIM];
  constexpr int MaxAttr() { return MAX_ATTR; }
  constexpr int MaxNumMat() { return MAX_NUM_MAT; }
  constexpr int Dim() { return DIM; }
};

template <int MAX_ATTR, int MAX_NUM_MAT, int DIM_FIRST, int DIM_SECOND>
struct MatCoeffPairContext
{
  MatCoeffContext<MAX_ATTR, MAX_NUM_MAT, DIM_FIRST> first;
  MatCoeffContext<MAX_ATTR, MAX_NUM_MAT, DIM_SECOND> second;
};

// XX TODO RUNTIME SCALABLE? MAYBE WITH TYPEDEFS?

#define MAT_COEFF_MAX_ATTR 32
#define MAT_COEFF_MAX_NUM_MAT (8 + 1)

template <int DIM>
using MatCoeffContextN = MatCoeffContext<MAT_COEFF_MAX_ATTR, MAT_COEFF_MAX_NUM_MAT, DIM>;

using MatCoeffContext1 = MatCoeffContextN<1>;
using MatCoeffContext2 = MatCoeffContextN<3>;
using MatCoeffContext3 = MatCoeffContextN<6>;

template <int DIM_FIRST, int DIM_SECOND>
using MatCoeffPairContextMN =
    MatCoeffPairContext<MAT_COEFF_MAX_ATTR, MAT_COEFF_MAX_NUM_MAT, DIM_FIRST, DIM_SECOND>;

using MatCoeffPairContext11 = MatCoeffPairContextMN<1, 1>;
using MatCoeffPairContext22 = MatCoeffPairContextMN<3, 3>;
using MatCoeffPairContext33 = MatCoeffPairContextMN<6, 6>;
using MatCoeffPairContext21 = MatCoeffPairContextMN<3, 1>;
using MatCoeffPairContext31 = MatCoeffPairContextMN<6, 1>;

CEED_QFUNCTION_HELPER CeedScalar CoeffUnpack(const MatCoeffContext1 *ctx, const CeedInt attr)
{
  const CeedInt k = ctx->attr_mat[attr];
  return ctx->mat_coeff[6 * k + 0];
}

CEED_QFUNCTION_HELPER void CoeffUnpack(const MatCoeffContext1 *ctx, const CeedInt attr, CeedScalar coeff[1])
{
  coeff[0] = CoeffUnpack(ctx, attr);
}

CEED_QFUNCTION_HELPER void CoeffUnpack(const MatCoeffContext2 *ctx, const CeedInt attr, CeedScalar coeff[3])
{
  const CeedInt k = ctx->attr_mat[attr];
  coeff[0] = ctx->mat_coeff[6 * k + 0];
  coeff[1] = ctx->mat_coeff[6 * k + 1];
  coeff[2] = ctx->mat_coeff[6 * k + 2];
}

CEED_QFUNCTION_HELPER void CoeffUnpack(const MatCoeffContext3 *ctx, const CeedInt attr, CeedScalar coeff[6])
{
  const CeedInt k = ctx->attr_mat[attr];
  coeff[0] = ctx->mat_coeff[6 * k + 0];
  coeff[1] = ctx->mat_coeff[6 * k + 1];
  coeff[2] = ctx->mat_coeff[6 * k + 2];
  coeff[3] = ctx->mat_coeff[6 * k + 3];
  coeff[4] = ctx->mat_coeff[6 * k + 4];
  coeff[5] = ctx->mat_coeff[6 * k + 5];
}

#endif  // PALACE_LIBCEED_TYPES_QF_H
