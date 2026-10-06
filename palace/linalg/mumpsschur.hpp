// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_MUMPS_SCHUR_HPP
#define PALACE_LINALG_MUMPS_SCHUR_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_MUMPS)

#include <vector>
#include <dmumps_c.h>

namespace palace
{

// MUMPS partial factorization with a Schur complement: factors the symmetric parent-space
// operator A (rows/cols outside the eliminated subsystem set to identity) with the given
// Schur variables left unfactored, returning the dense Schur complement on rank 0 -- for
// the environment, S_E = A_GG - A_GE A_EE^-1 A_EG from one partial factorization instead of
// |Gamma| back-solves. The same factorization then solves the internal problem (A_EE, the
// Schur variables fixed at 0), which serves every other environment solve.
class MumpsSchurSolver
{
public:
  // schur_vars: global (0-based) true-DOF indices of the Schur variables, in the order the
  // Schur rows/columns should appear (replicated on all ranks); empty for a plain
  // factorization of A. blr_tol > 0 enables a block low-rank (BLR) factorization with that
  // relative accuracy (0: exact). With serial, rank 0 factors alone (for a small system, it
  // avoids MUMPS's per-rank workspace). With refactor, the input entries are kept so that
  // Refactor can factor new values with the same pattern, reusing the analysis.
  MumpsSchurSolver(const mfem::HypreParMatrix &A,
                   const std::vector<HYPRE_BigInt> &schur_vars, double blr_tol = 0.0,
                   bool serial = false, bool refactor = false);

  ~MumpsSchurSolver();

  // Factor A with new values and the sparsity pattern of the analysis (collective).
  void Refactor(const mfem::HypreParMatrix &A);

  // Dense Schur complement on rank 0 (n_schur x n_schur, column-major, symmetric).
  const std::vector<double> &Schur() const { return schur; }

  // Internal solves for a batch of parent true-DOF vectors (collective): y_k solves A_11 on
  // the internal (non-Schur) variables, with 0 on the Schur variables. y_k may alias x_k.
  void SolveInternal(const std::vector<const mfem::Vector *> &X,
                     const std::vector<mfem::Vector *> &Y);

private:
  void Factor();
  void Check(const char *phase) const;

  MPI_Comm comm;
  int rank = 0;
  bool serial = false, refactor = false;
  bool active = true;  // this rank takes part in the factorization
  HYPRE_BigInt n_glob;
  int n_loc, n_schur;
  std::vector<int> row_cnt, row_disp;
  std::vector<MUMPS_INT> irn, jcn, listvar;
  std::vector<double> val, schur, rhs;
  double blr_tol = 0.0, scale = 1.0;
  DMUMPS_STRUC_C id{};
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_MUMPS_SCHUR_HPP
