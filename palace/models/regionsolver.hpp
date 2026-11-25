// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_REGION_SOLVER_HPP
#define PALACE_MODELS_REGION_SOLVER_HPP

#include <complex>
#include <functional>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/blockcyclic.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/schursolver.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

class EnvironmentOperator;
class IoData;
class SideSystem;
class SpaceOperator;

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
  // the local part of S_E on the grid to that of a matrix. With refactor, the entries are
  // kept for Refactor. Collective.
  RegionSolver(const BlockCyclicGrid &grid, LinearSolver type,
               SymbolicFactorization ordering, HYPRE_BigInt n_glob, int n_loc,
               typename SchurFactorization<T>::Coo &&A,
               const std::vector<HYPRE_BigInt> &gamma, std::vector<int> rows,
               const std::function<void(std::vector<T> &)> &add_interface,
               bool refactor = false);

  // The entries of the next Refactor, in the order of the input (with refactor).
  std::vector<T> &Values() { return factorization->Values(); }

  // Factor the entries in Values(), and S_R + S_E with add_interface. Collective.
  void Refactor(const std::function<void(std::vector<T> &)> &add_interface);

  // Solves for a batch, with 0 on the entries of the vectors outside the rows of the
  // system. Y may alias X. Collective.
  void Mult(const std::vector<const VecType *> &X, const std::vector<VecType *> &Y) const;

private:
  const BlockCyclicGrid &grid;
  std::unique_ptr<SchurFactorization<T>> factorization;
  // S_R + S_E on the grid, factored: in MUMPS's Schur buffer, or in dense.
  std::vector<T> dense;
  std::unique_ptr<BlockCyclicFactorization<T>> factors;

  // Factor S_R + S_E, from the last factorization of the system and add_interface.
  void FactorInterface(const std::function<void(std::vector<T> &)> &add_interface);
};

//
// The direct solver of the region of a driven substructuring with the environment's S_E on
// Γ (see SpaceOperator::SetEnvironment): a RegionSolver of the preconditioner matrix of the
// region, which leaves S_E out, on its free DOFs (the Dirichlet and pinned ones are
// identity rows), with S_E(ω) on the grid of the environment (EnvironmentOperator::Grid),
// refactored at each frequency. It is exact but for rounding, except near a
// resonance of the region with Γ held fixed, where the region's interior, factored alone,
// is nearly singular: it is used within FGMRES with the system operator (see
// MakeDrivenKspSolver), which corrects that in a few iterations.
//
class DrivenRegionSolver : public Solver<ComplexOperator>
{
public:
  DrivenRegionSolver(const IoData &iodata, SpaceOperator &space_op,
                     EnvironmentOperator &env);
  ~DrivenRegionSolver() override;

  void SetOperator(const ComplexOperator &op) override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;

private:
  SpaceOperator &space_op;
  EnvironmentOperator &env;
  const LinearSolver type;
  const SymbolicFactorization ordering;

  // The system: the true DOFs which are not eliminated (Dirichlet or pinned), whose solves
  // keep the identity rows of the eliminated ones, and its solver.
  std::unique_ptr<SideSystem> side;
  std::vector<int> eliminated;
  std::unique_ptr<RegionSolver<std::complex<double>>> region;
};

// The linear solver of a driven sweep: that of the configuration or, for a substructured
// region, FGMRES with the DrivenRegionSolver (and the tolerances of the configuration).
std::unique_ptr<ComplexKspSolver> MakeDrivenKspSolver(const IoData &iodata,
                                                      SpaceOperator &space_op);

}  // namespace palace

#endif  // PALACE_MODELS_REGION_SOLVER_HPP
