// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_REGION_SOLVER_HPP
#define PALACE_MODELS_REGION_SOLVER_HPP

#include <functional>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/blockcyclic.hpp"
#include "linalg/schursolver.hpp"
#include "utils/labels.hpp"

namespace palace
{

//
// The direct solver of a substructured region with the environment's S_E on Γ. The
// region's system is factored for its Schur complement S_R onto Γ, on the grid of Γ (by
// MUMPS) or on rank 0 and then scattered to it, and S_R + S_E by a dense factorization on
// the grid (BlockCyclicFactorization: Cholesky's for a real system, which must be symmetric
// positive definite, else LU). A solve condenses onto Γ, solves there and expands, which is
// exact but for rounding.
//
template <typename T>
class RegionSolver
{
public:
  using VecType = typename SchurFactorization<T>::VecType;

  // The system from the lower-triangle COO entries of its n_loc local rows (of n_glob, in
  // rank order) and the local index of each in the vectors of the solves (rows), as
  // MakeSchurFactorization takes them, with gamma, its (0-based) rows of Γ in the order of
  // the grid (replicated), factored by the sparse direct solver type. add_interface adds
  // the local part of S_E on the grid to that of a matrix. Collective.
  RegionSolver(const BlockCyclicGrid &grid, LinearSolver type,
               SymbolicFactorization ordering, HYPRE_BigInt n_glob, int n_loc,
               typename SchurFactorization<T>::Coo &&A,
               const std::vector<HYPRE_BigInt> &gamma, std::vector<int> rows,
               const std::function<void(std::vector<T> &)> &add_interface);

  // Solves for a batch, with 0 on the entries of the vectors outside the rows of the
  // system. Y may alias X. Collective.
  void Mult(const std::vector<const VecType *> &X, const std::vector<VecType *> &Y) const;

private:
  const BlockCyclicGrid &grid;
  std::unique_ptr<SchurFactorization<T>> factorization;
  // S_R + S_E on the grid, factored: in MUMPS's Schur buffer, or in dense.
  std::vector<T> dense;
  std::unique_ptr<BlockCyclicFactorization<T>> factors;
};

}  // namespace palace

#endif  // PALACE_MODELS_REGION_SOLVER_HPP
