// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_TYPES_QF_H
#define PALACE_LIBCEED_TYPES_QF_H

struct GeomFactorContext
{
  bool compute_wdetJ, compute_adjJt, compute_J;
};

template <int MAX_ATTR, int NUM_MAT, int NUM_COEFF_COMP>
struct MatCoeffContext
{
  CeedInt attr_mat[MAX_ATTR];
  CeedScalar mat_coeff[NUM_MAT][NUM_COEFF_COMP];
  constexpr int MaxAttr() { return MAX_ATTR; }
  constexpr int NumMat() { return NUM_MAT; }
  constexpr int NumCoeffComp() { return NUM_COEFF_COMP; }
  CeedScalar *begin() { return &mat_coeff[0][0]; }
  CeedScalar *end() { return &mat_coeff[NUM_MAT - 1][NUM_COEFF_COMP - 1] + 1; }
};

template <int MAX_ATTR, int NUM_MAT, int NUM_COEFF_COMP1, int NUM_COEFF_COMP2>
struct MatCoeffPairContext
{
  MatCoeffContext<MAX_ATTR, NUM_MAT, NUM_COEFF_COMP1> ctx1;
  MatCoeffContext<MAX_ATTR, NUM_MAT, NUM_COEFF_COMP2> ctx2;
};

#define MAX_ATTR 32
#define NUM_MAT 8

using MatCoeffContext1 = MatCoeffContext<MAX_ATTR, NUM_MAT, 1>;
using MatCoeffContext2 = MatCoeffContext<MAX_ATTR, NUM_MAT, 3>;
using MatCoeffContext3 = MatCoeffContext<MAX_ATTR, NUM_MAT, 6>;

using MatCoeffPairContext11 = MatCoeffPairContext<MAX_ATTR, NUM_MAT, 1, 1>;
using MatCoeffPairContext22 = MatCoeffPairContext<MAX_ATTR, NUM_MAT, 3, 3>;
using MatCoeffPairContext33 = MatCoeffPairContext<MAX_ATTR, NUM_MAT, 6, 6>;
using MatCoeffPairContext21 = MatCoeffPairContext<MAX_ATTR, NUM_MAT, 3, 1>;
using MatCoeffPairContext31 = MatCoeffPairContext<MAX_ATTR, NUM_MAT, 6, 1>;

#endif  // PALACE_LIBCEED_TYPES_QF_H
