// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_KSP_HPP
#define PALACE_LINALG_KSP_HPP

#include <memory>
#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include "linalg/iterative.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "utils/labels.hpp"

namespace palace
{

class FiniteElementSpaceHierarchy;
class IoData;

namespace fem::singular
{

struct ParallelFeaturePatches;

}  // namespace fem::singular

namespace config
{

struct LinearSolverData;

}  // namespace config

// Derive the preconditioner matrix symmetry (used by sparse direct solvers) from the
// problem type, shifted-preconditioner flag, and periodicity. This is a derived runtime
// property, not a configuration field, so it is computed here rather than stored on
// LinearSolverData.
MatrixSymmetry GetPreconditionerMatrixSymmetry(const IoData &iodata);

//
// Linear solver class composing an iterative solver and preconditioner object.
//
template <typename OperType>
class BaseKspSolver
{
  static_assert(std::is_same<OperType, Operator>::value ||
                    std::is_same<OperType, ComplexOperator>::value,
                "Solver can only be defined for OperType = Operator or ComplexOperator!");

  using VecType = typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                            ComplexVector, Vector>::type;

protected:
  // The actual solver and preconditioner objects.
  std::unique_ptr<IterativeSolver<OperType>> ksp;
  std::unique_ptr<Solver<OperType>> pc;

  // Counters for number of calls to Mult method for linear solves, cumulative number of
  // iterations, and failed solves.
  mutable int ksp_mult, ksp_mult_it, ksp_mult_failures;

  // Enable timer contribution for Timer::KSP_PRECONDITIONER.
  bool use_timer;

public:
  BaseKspSolver(const config::LinearSolverData &linear, MatrixSymmetry pc_mat_sym,
                int verbose, FiniteElementSpaceHierarchy &fespaces,
                FiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  BaseKspSolver(const IoData &iodata, FiniteElementSpaceHierarchy &fespaces,
                FiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  BaseKspSolver(std::unique_ptr<IterativeSolver<OperType>> &&ksp,
                std::unique_ptr<Solver<OperType>> &&pc);

  int NumTotalMult() const { return ksp_mult; }
  int NumTotalMultIterations() const { return ksp_mult_it; }
  int NumFailedMult() const { return ksp_mult_failures; }

  // Forward tolerance access to the underlying iterative solver.
  double GetRelTol() const { return ksp->GetRelTol(); }
  double GetAbsTol() const { return ksp->GetAbsTol(); }
  void SetRelTol(double tol) { ksp->SetRelTol(tol); }
  void SetAbsTol(double tol) { ksp->SetAbsTol(tol); }

  void EnableTimer() { use_timer = true; }

  void SetOperators(const OperType &op, const OperType &pc_op);

  void Mult(const VecType &x, VecType &y) const;
};

using KspSolver = BaseKspSolver<Operator>;
using ComplexKspSolver = BaseKspSolver<ComplexOperator>;

// Configure the correctness-first electrostatic singular-element
// preconditioner: symmetric standard-space AMG correction and an additive sum
// of exact overlapping straight-feature corrections. The full system remains
// in its original physical coefficient coordinates.
std::unique_ptr<KspSolver>
MakeSingularPatchKspSolver(const IoData &iodata, FiniteElementSpaceHierarchy &fespaces,
                           const Operator &full_operator, const Operator &standard_operator,
                           const fem::singular::ParallelFeaturePatches &feature_patches);
bool UsesSingularPatchKspSolver(const IoData &iodata);

// Configure correctness-first solvers for assembled combined standard-plus-singular
// systems. The exact SuperLU preconditioner acts on the complete operator and therefore
// does not rely on standard MFEM finite element space dimensions or auxiliary-space
// interpolation operators.
std::unique_ptr<KspSolver> MakeSingularDirectKspSolver(const IoData &iodata, MPI_Comm comm);
std::unique_ptr<ComplexKspSolver> MakeSingularComplexKspSolver(const IoData &iodata,
                                                               MPI_Comm comm);

// Construct a coarse solver for one combined standard-plus-singular block. AMS acts on
// the standard H(curl) subspace and is completed by a symmetric enrichment correction;
// the other solver types act on the full combined matrix directly.
std::unique_ptr<Solver<ComplexOperator>>
MakeSingularComplexCoarseSolver(const config::LinearSolverData &linear, LinearSolver type,
                                MatrixSymmetry matrix_symmetry, int verbose, MPI_Comm comm,
                                int standard_size,
                                FiniteElementSpaceHierarchy *primary_fespaces = nullptr,
                                FiniteElementSpaceHierarchy *auxiliary_fespaces = nullptr);

// Configure the combined standard-plus-singular Maxwell solve. The outer Krylov
// method remains user-configured. For multiple polynomial levels, the
// preconditioner is geometric multigrid with combined transfer
// diag(P_standard, I_enrichment) and a complete p=1 coarse operator.
std::unique_ptr<ComplexKspSolver> MakeSingularComplexKspSolver(
    const IoData &iodata, FiniteElementSpaceHierarchy &nd_fespaces,
    FiniteElementSpaceHierarchy &h1_fespaces,
    const std::vector<const Operator *> &combined_prolongations,
    const std::vector<const Operator *> &combined_gradients,
    const std::vector<mfem::Array<int>> &combined_nd_essential_tdofs,
    const std::vector<mfem::Array<int>> &combined_h1_essential_tdofs);

}  // namespace palace

#endif  // PALACE_LINALG_KSP_HPP
