// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_REGION_SOLVER_HPP
#define PALACE_MODELS_REGION_SOLVER_HPP

#include <complex>
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
class SpaceOperator;

//
// The direct solver of the region of a driven substructuring with the environment's S_E on
// Γ (see SpaceOperator::SetEnvironment). The preconditioner matrix of the region, which
// leaves S_E out, is factored on its free DOFs (the Dirichlet and pinned ones are identity
// rows) for its Schur complement S_R onto Γ, and the dense S_R + S_E(ω), 2D block cyclic on
// the grid of the environment (EnvironmentOperator::Grid), by an LU factorization. A solve
// condenses onto Γ, solves there and expands. It is exact but for rounding, except near a
// resonance of the region with Γ held fixed, where the region's interior, factored alone,
// is nearly singular: it is used within FGMRES with the system operator (see
// MakeDrivenKspSolver), which corrects that in a few iterations.
//
class RegionSolver : public Solver<ComplexOperator>
{
public:
  RegionSolver(const IoData &iodata, SpaceOperator &space_op, EnvironmentOperator &env);
  ~RegionSolver() override;

  void SetOperator(const ComplexOperator &op) override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;

private:
  SpaceOperator &space_op;
  EnvironmentOperator &env;
  const LinearSolver type;
  const SymbolicFactorization ordering;

  // The system: the free true DOFs of this process (rows) and the eliminated ones, and the
  // pattern of the local rows of the matrix, with the positions of the system's entries.
  std::vector<int> rows, eliminated, row_ptr, jcn, kept;
  std::unique_ptr<SchurFactorization<std::complex<double>>> factorization;

  // S_R + S_E on the grid and its LU factors: in MUMPS's Schur buffer, or in dense for a
  // Schur complement on rank 0.
  std::vector<std::complex<double>> dense, *C = nullptr;
  std::unique_ptr<BlockCyclicLU> lu;

  // The system, and its first factorization, from the parts of the matrix.
  void Setup(const mfem::HypreParMatrix *Ar, const mfem::HypreParMatrix *Ai);
};

// The linear solver of a driven sweep: that of the configuration or, for a substructured
// region, FGMRES with the RegionSolver (and the tolerances of the configuration).
std::unique_ptr<ComplexKspSolver> MakeDrivenKspSolver(const IoData &iodata,
                                                      SpaceOperator &space_op);

}  // namespace palace

#endif  // PALACE_MODELS_REGION_SOLVER_HPP
