// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_DIRECT_SCHUR_HPP
#define PALACE_LINALG_DIRECT_SCHUR_HPP

#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/schursolver.hpp"
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

class DirectSchurFactorization;

//
// The Schur complement S = A_SS - A_SI A_II^-1 A_IS of a complex symmetric sparse matrix A
// onto some of its variables (the Schur variables), with a sparse direct solver other than
// MUMPS: A_II is factored by STRUMPACK or by SuperLU_DIST in its real equivalent form
// [Re A_II, Im A_II; Im A_II, -Re A_II], and S is formed from one solve with A_II per Schur
// variable, in batches (MumpsSchurSolverT's partial factorization gives S without them).
// The vectors of the solves are those of the rows, 0 outside them.
//
class DirectSchurSolver : public SchurFactorization<std::complex<double>>
{
public:
  using T = std::complex<double>;

  // From the entries of an n_glob x n_glob matrix whose rows are numbered by rank (n_loc
  // local rows), each rank passing some of them, in any rows, the Schur variables (global
  // 0-based indices, replicated, in their order in S; none for a plain factorization), and
  // rows: the local index of each local row in the vectors of the solves (-1 for a row of
  // none). Factors the matrix with type (STRUMPACK or SUPERLU), collective.
  DirectSchurSolver(MPI_Comm comm, HYPRE_BigInt n_glob, int n_loc, Coo &&A,
                    const std::vector<HYPRE_BigInt> &schur_vars, std::vector<int> rows,
                    LinearSolver type);
  ~DirectSchurSolver() override;

  std::vector<T> &Values() override { return val; }

  void Refactor() override;

  std::vector<T> ReleaseSchur() override;

  void SolveInternal(const std::vector<const ComplexVector *> &X,
                     const std::vector<ComplexVector *> &Y) override;

  void Reduce(const std::vector<const ComplexVector *> &B, std::vector<T> &red) override;

  void Expand(const std::vector<T> &u, const std::vector<ComplexVector *> &X) override;

private:
  MPI_Comm comm;
  int rank = 0, n_loc = 0, n_schur = 0;
  std::vector<int> rows;
  std::vector<T> val;

  // The local rows: their position among the Schur variables, or their local index in
  // A_II (local rows of A_II numbered from int_off), -1 for the other kind.
  std::vector<int> schur_pos, int_row;
  int n_int = 0;
  HYPRE_BigInt int_off = 0, n_int_glob = 0;

  // Where the entries go (their transposes too, A being symmetric): the input entries in
  // the order sent to each rank, and, received, the slot of each entry in the local values
  // of A_II, A_IS (CSR over the local rows of A_II) or A_SS (dense, on rank 0).
  std::vector<int> send_src, send_cnt, send_dsp, recv_cnt, recv_dsp;
  std::vector<int> recv_kind, recv_slot;
  std::vector<int> ii_ptr, is_ptr, is_col;
  std::vector<HYPRE_BigInt> ii_col;
  std::vector<T> ii_val, is_val, ss_val;

  std::unique_ptr<DirectSchurFactorization> fact;
  std::vector<T> schur, y_last;
  int reduced = 0;

  void Distribute();
  void Factor();
  // x = A_II^-1 b for nb local column-major blocks of A_II's local rows (collective).
  void SolveInterior(const std::vector<T> &b, std::vector<T> &x, int nb) const;
  // y += sign A_SI x for nb local blocks x, y on rank 0 (n_schur x nb, column-major)
  // (collective).
  void AddSchurProduct(const std::vector<T> &x, int nb, double sign, T *y) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_DIRECT_SCHUR_HPP
