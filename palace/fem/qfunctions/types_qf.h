// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_TYPES_QF_H
#define PALACE_LIBCEED_TYPES_QF_H

struct GeomFactorContext
{
  bool compute_wdetJ, compute_adjJt, compute_J;
};

template <int MAX_ATTR, int MAX_NUM_MAT, int NUM_COEFF_COMP>
struct MatCoeffContext
{
  CeedInt attr_mat[MAX_ATTR];
  CeedScalar mat_coeff[MAX_NUM_MAT][NUM_COEFF_COMP];
  constexpr int MaxAttr() { return MAX_ATTR; }
  constexpr int MaxNumMat() { return MAX_NUM_MAT; }
  constexpr int NumCoeffComp() { return NUM_COEFF_COMP; }
};

template <int MAX_ATTR, int MAX_NUM_MAT, int NUM_COEFF_COMP1, int NUM_COEFF_COMP2>
struct MatCoeffPairContext
{
  MatCoeffContext<MAX_ATTR, MAX_NUM_MAT, NUM_COEFF_COMP1> ctx1;
  MatCoeffContext<MAX_ATTR, MAX_NUM_MAT, NUM_COEFF_COMP2> ctx2;
};

// XX TODO RUNTIME SCALABLE? MAYBE WITH TYPEDEFS?

#define MAT_COEFF_MAX_ATTR 32
#define MAT_COEFF_MAX_NUM_MAT (8 + 1)

template <int NUM_COEFF_COMP>
using MatCoeffContextN =
    MatCoeffContext<MAT_COEFF_MAX_ATTR, MAT_COEFF_MAX_NUM_MAT, NUM_COEFF_COMP>;

using MatCoeffContext1 = MatCoeffContextN<1>;
using MatCoeffContext2 = MatCoeffContextN<3>;
using MatCoeffContext3 = MatCoeffContextN<6>;

template <int NUM_COEFF_COMP1, int NUM_COEFF_COMP2>
using MatCoeffPairContextMN = MatCoeffPairContext<MAT_COEFF_MAX_ATTR, MAT_COEFF_MAX_NUM_MAT,
                                                  NUM_COEFF_COMP1, NUM_COEFF_COMP2>;

using MatCoeffPairContext11 = MatCoeffPairContextMN<1, 1>;
using MatCoeffPairContext22 = MatCoeffPairContextMN<3, 3>;
using MatCoeffPairContext33 = MatCoeffPairContextMN<6, 6>;
using MatCoeffPairContext21 = MatCoeffPairContextMN<3, 1>;
using MatCoeffPairContext31 = MatCoeffPairContextMN<6, 1>;

#endif  // PALACE_LIBCEED_TYPES_QF_H
