// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SCHUR_FACTORY_HPP
#define PALACE_LINALG_SCHUR_FACTORY_HPP

#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "linalg/schursolver.hpp"
#include "utils/labels.hpp"

namespace palace
{

class BlockCyclicGrid;

// The block low-rank (BLR) tolerance of MUMPS's factorizations, at round-off: faster than
// full rank at no loss of accuracy.
inline constexpr double kBlrTol = 1.0e-14;

// The factorizations for a Schur complement: the condensation of an environment onto its
// interface, and the Schur complement of a region, which RegionSolver adds S_E to.
enum class SchurKind
{
  ENVIRONMENT,
  REGION
};

// The largest Schur complement of a region, and the fewest processes, for which STRUMPACK
// gives it by default (see SchurSolverFor).
inline constexpr HYPRE_BigInt kRegionStrumpackSchur = 10000;
inline constexpr int kRegionStrumpackProcs = 32;

// The sparse direct solver of a factorization for a Schur complement: the one requested if
// Palace is built with it (STRUMPACK, MUMPS or SuperLU_DIST), else the first it is built
// with of MUMPS, STRUMPACK and SuperLU_DIST. The Schur complement of a region of at most
// kRegionStrumpackSchur variables on num_procs >= kRegionStrumpackProcs processes takes
// STRUMPACK first, with its Schur complement interface: STRUMPACK gathers it on one
// process, and MUMPS distributes it but with a workspace of a few tenths of a GB on each
// process, which from about 30 processes outweighs it. On fewer, MUMPS is faster and
// smaller.
LinearSolver SchurSolverFor(LinearSolver requested, SchurKind kind,
                            HYPRE_BigInt n_schur = 0, int num_procs = 0);

// The name of a sparse direct solver, and of its factorization for a Schur complement:
// partial (MUMPS's, or STRUMPACK's with its Schur complement interface), or with one solve
// per Schur variable.
const char *DirectSolverName(LinearSolver type);
std::string SchurFactorizationName(LinearSolver type);

// The factorization for its Schur complement of a symmetric matrix (real, or complex
// symmetric) from its lower-triangle COO entries, as the constructors of MumpsSchurSolverT
// and DirectSchurSolverT take them, by the direct solver of type (see SchurSolverFor):
// spd if A is real positive definite, and with refactor if it is factored again with new
// values. With a grid, MUMPS's Schur complement is distributed on it (SchurBlocks).
// Collective.
template <typename T>
std::unique_ptr<SchurFactorization<T>>
MakeSchurFactorization(LinearSolver type, SymbolicFactorization ordering, MPI_Comm comm,
                       HYPRE_BigInt n_glob, int n_loc,
                       typename SchurFactorization<T>::Coo &&coo,
                       const std::vector<HYPRE_BigInt> &schur_vars, std::vector<int> rows,
                       bool spd, bool refactor, const BlockCyclicGrid *grid = nullptr);

}  // namespace palace

#endif  // PALACE_LINALG_SCHUR_FACTORY_HPP
