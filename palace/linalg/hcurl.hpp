// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_HCURL_HPP
#define PALACE_LINALG_HCURL_HPP

#include <memory>
#include <vector>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace mfem
{

template <typename T>
class Array;

}  // namespace mfem

namespace palace
{

class FiniteElementSpaceHierarchy;
class MaterialOperator;

//
// This solver implements a solver for the operator K + M in a Nedelec space.
//
template <typename VecType>
class WeightedHCurlNormSolver
{
  using OperType = typename std::conditional<std::is_same<VecType, ComplexVector>::value,
                                             ComplexOperator, Operator>::type;

private:
  // H(curl) norm operator A = K + M and its projection Gᵀ A G.
  std::unique_ptr<OperType> A;

  // Linear solver for the linear system A y = x;
  std::unique_ptr<BaseKspSolver<OperType>> ksp;

public:
  WeightedHCurlNormSolver(const MaterialOperator &mat_op,
                          FiniteElementSpaceHierarchy &nd_fespaces,
                          FiniteElementSpaceHierarchy &h1_fespaces,
                          const std::vector<mfem::Array<int>> &nd_dbc_tdof_lists,
                          const std::vector<mfem::Array<int>> &h1_dbc_tdof_lists,
                          double tol, int max_it, int print);

  const OperType &GetOperator() { return *A; }

  void Mult(const VecType &x, VecType &y) const { ksp->Mult(x, y); }
};

}  // namespace palace

#endif  // PALACE_LINALG_HCURL_HPP
