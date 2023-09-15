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
class ParFiniteElementSpaceHierarchy;

}  // namespace mfem

namespace palace
{

class MaterialOperator;

//
// This solver implements a solver for the operator K + M in a Nedelec space.
//
class WeightedHCurlNormSolver
{
private:
  // H(curl) norm operator A = K + M and its projection Gᵀ A G.
  std::unique_ptr<Operator> A;

  // Linear solver for the linear system A y = x;
  std::unique_ptr<KspSolver> ksp;

public:
  WeightedHCurlNormSolver(const MaterialOperator &mat_op,
                          const mfem::ParFiniteElementSpaceHierarchy &nd_fespaces,
                          const mfem::ParFiniteElementSpaceHierarchy &h1_fespaces,
                          const std::vector<mfem::Array<int>> &nd_dbc_tdof_lists,
                          const std::vector<mfem::Array<int>> &h1_dbc_tdof_lists,
                          double tol, int max_it, int print, int pa_order_threshold,
                          bool pa_discrete_interp);

  const Operator &GetOperator() { return *A; }

  void Mult(const Vector &x, Vector &y) const { ksp->Mult(x, y); }

  void Mult(const ComplexVector &x, ComplexVector &y)
  {
    Mult(x.Real(), y.Real());
    Mult(x.Imag(), y.Imag());
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_HCURL_HPP
