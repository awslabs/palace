// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFF_QF_H
#define PALACE_LIBCEED_COEFF_QF_H

union CeedIntScalar
{
  CeedInt first;
  CeedScalar second;
};

// The first entry of ctx is the number of (1-based) attributes, followed by the entries of
// the attribute to material index array (these are 0-based).
// The next entry is the number of material property coefficients, followed by the
// coefficients.
// Pair coefficients are two coefficient contexts arranged contiguously in memory.

CEED_QFUNCTION_HELPER CeedInt NumAttr(const CeedIntScalar *ctx)
{
  return ctx[0].first;
}

CEED_QFUNCTION_HELPER CeedInt NumMat(const CeedIntScalar *ctx)
{
  return ctx[1 + NumAttr(ctx)].first;
}

CEED_QFUNCTION_HELPER const CeedIntScalar *AttrMat(const CeedIntScalar *ctx)
{
  return ctx + 1;
}

CEED_QFUNCTION_HELPER const CeedIntScalar *MatCoeff(const CeedIntScalar *ctx)
{
  return ctx + 2 + NumAttr(ctx);
}

template <int DIM>
CEED_QFUNCTION_HELPER const CeedIntScalar *CoeffPairSecond(const CeedIntScalar *ctx)
{
  return ctx + 2 + NumAttr(ctx) + (DIM * DIM) * NumMat(ctx);
}

#endif  // PALACE_LIBCEED_COEFF_QF_H
