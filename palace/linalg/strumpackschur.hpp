// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_STRUMPACK_SCHUR_HPP
#define PALACE_LINALG_STRUMPACK_SCHUR_HPP

#include <mfem.hpp>

#if defined(PALACE_WITH_STRUMPACK_SCHUR)

#include <vector>
#include <StrumpackSparseSolverMPIDist.hpp>
#include "linalg/schursolver.hpp"
#include "utils/labels.hpp"

namespace palace
{

//
// STRUMPACK's partial factorization with a Schur complement: the reordering places the
// Schur variables last and the factorization leaves them unfactored, so that it gives
// S = A_SS - A_SI A_II^-1 A_IS on rank 0 (which sends their rows to the ranks holding
// them), and its factors of A_II the internal solves.
//
class StrumpackSchurSolver : public SchurSolver<double>
{
public:
  // schur_vars: the global (0-based) rows of the Schur variables, in the order of their
  // rows and columns in S (replicated). ordering is the fill-reducing ordering ("Default":
  // STRUMPACK's); the Schur complement needs a sequential one, so a parallel ordering is
  // replaced by its sequential counterpart. Factors A (collective).
  StrumpackSchurSolver(const mfem::HypreParMatrix &A,
                       const std::vector<HYPRE_BigInt> &schur_vars,
                       SymbolicFactorization ordering = SymbolicFactorization::DEFAULT);

  std::vector<double> ReleaseSchurRows() override;

  void SolveInternal(const std::vector<const mfem::Vector *> &X,
                     const std::vector<mfem::Vector *> &Y) override;

  void Condense(const std::vector<const mfem::Vector *> &X,
                std::vector<double> &g) override;

private:
  MPI_Comm comm;
  int n_loc, n_schur;
  strumpack::SparseSolverMPIDist<double, int> solver;
  SchurRowGather rows;
  std::vector<double> schur;
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_STRUMPACK_SCHUR_HPP
