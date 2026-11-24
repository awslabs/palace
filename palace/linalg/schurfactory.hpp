// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SCHUR_FACTORY_HPP
#define PALACE_LINALG_SCHUR_FACTORY_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/schursolver.hpp"
#include "utils/labels.hpp"

namespace palace
{

// The block low-rank (BLR) tolerance of MUMPS's factorizations, at round-off: faster than
// full rank at no loss of accuracy.
inline constexpr double kBlrTol = 1.0e-14;

// What a sparse direct solver factors: a system, or a matrix for its Schur complement.
enum class Factorization
{
  SYSTEM,
  SCHUR_COMPLEMENT
};

// The sparse direct solver of a factorization: the one requested if Palace is built with it
// (STRUMPACK, MUMPS or SuperLU_DIST), else the first it is built with, in the order
// STRUMPACK, MUMPS, SuperLU_DIST for a system, and MUMPS, STRUMPACK, SuperLU_DIST for a
// Schur complement (MUMPS's partial factorization, then STRUMPACK's, else one solve per
// Schur variable).
LinearSolver DirectSolverFor(LinearSolver requested, Factorization kind);

// The name of a sparse direct solver.
const char *DirectSolverName(LinearSolver type);

// The factorization for its Schur complement of a symmetric matrix (real, or complex
// symmetric) from its lower-triangle COO entries, as the constructors of MumpsSchurSolverT
// and DirectSchurSolverT take them, by the direct solver of type (see DirectSolverFor):
// spd if A is real positive definite, and with refactor if it is factored again with new
// values. Collective.
template <typename T>
std::unique_ptr<SchurFactorization<T>>
MakeSchurFactorization(LinearSolver type, SymbolicFactorization ordering, MPI_Comm comm,
                       HYPRE_BigInt n_glob, int n_loc,
                       typename SchurFactorization<T>::Coo &&coo,
                       const std::vector<HYPRE_BigInt> &schur_vars, std::vector<int> rows,
                       bool spd, bool refactor);

}  // namespace palace

#endif  // PALACE_LINALG_SCHUR_FACTORY_HPP
