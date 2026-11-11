// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_STRUMPACK_SCHUR_HPP
#define PALACE_LINALG_STRUMPACK_SCHUR_HPP

#include <mfem.hpp>

#if defined(PALACE_WITH_STRUMPACK_SCHUR)

#include <vector>
#include <StrumpackSparseSolverMPIDist.hpp>
#include "linalg/schursolver.hpp"

namespace palace
{

//
// STRUMPACK's partial factorization with a Schur complement: the reordering places the
// Schur variables last and the factorization leaves them unfactored, so that it gives
// S = A_SS - A_SI A_II^-1 A_IS on rank 0, and its factors of A_II the internal solves.
//
class StrumpackSchurSolver : public SchurSolver<double>
{
public:
  // schur_vars: the global (0-based) rows of the Schur variables, in the order of their
  // rows and columns in S (replicated). Factors A (collective).
  StrumpackSchurSolver(const mfem::HypreParMatrix &A,
                       const std::vector<HYPRE_BigInt> &schur_vars);

  std::vector<double> ReleaseSchur() override;

  void SolveInternal(const std::vector<const mfem::Vector *> &X,
                     const std::vector<mfem::Vector *> &Y) override;

private:
  int n_loc;
  strumpack::SparseSolverMPIDist<double, int> solver;
  std::vector<double> schur;
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_STRUMPACK_SCHUR_HPP
