// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_DIV_FREE_HPP
#define PALACE_LINALG_DIV_FREE_HPP

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
class FiniteElementSpace;
class MaterialOperator;

//
// This solver implements a projection onto a divergence-free space satisfying Gᵀ M x = 0,
// where G represents the discrete gradient matrix with columns spanning the nullspace of
// the curl-curl operator.
//
template <typename VecType>
class DivFreeSolver
{
  using OperType = typename std::conditional<std::is_same<VecType, ComplexVector>::value,
                                             ComplexOperator, Operator>::type;

private:
  // Operators for the divergence-free projection.
  std::unique_ptr<OperType> M;
  std::unique_ptr<Operator> WeakDiv;
  const Operator *Grad;
  const mfem::Array<int> *bdr_tdof_list_M;

  // Optional storage for homogeneous Dirichlet boundary condition on a single true dof,
  // used when the input array of H1 boundary dofs is empty to prevent the Poisson operator
  // from being singular.
  std::vector<mfem::Array<int>> aux_tdof_lists;

  // Linear solver for the projected linear system (Gᵀ M G) y = x.
  std::unique_ptr<BaseKspSolver<OperType>> ksp;

  // Workspace objects for solver application.
  mutable VecType psi, rhs;

public:
  DivFreeSolver(const MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
                FiniteElementSpaceHierarchy &h1_fespaces,
                const std::vector<mfem::Array<int>> &h1_bdr_tdof_lists, double tol,
                int max_it, int print);

  // Given a vector of Nedelec dofs for an arbitrary vector field, compute the Nedelec dofs
  // of the irrotational portion of this vector field. The resulting vector will satisfy
  // ∇ x y = 0.
  void Mult(VecType &y) const;

  void Mult(const VecType &x, VecType &y) const
  {
    y = x;
    Mult(y);
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_DIV_FREE_HPP
