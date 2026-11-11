// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_DIRECT_SCHUR_HPP
#define PALACE_LINALG_DIRECT_SCHUR_HPP

#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/vector.hpp"

namespace palace
{

class DirectSchurFactorization;
class StrumpackSchur;  // PROBE

//
// The Schur complement S = A_SS - A_SI A_II^-1 A_IS of a complex symmetric sparse matrix A
// onto some of its variables (the Schur variables), with a sparse direct solver other than
// MUMPS: A_II is factored by STRUMPACK or, without it, by SuperLU_DIST in its real
// equivalent form [Re A_II, Im A_II; Im A_II, -Re A_II], and S is formed from one solve
// with A_II per Schur variable, in batches. It has the interface of MumpsSchurSolverT,
// whose partial factorization gives S without the solves.
//
class DirectSchurSolver
{
public:
  using T = std::complex<double>;

  // Lower-triangle COO entries (1-based global indices) of the matrix, as for
  // MumpsSchurSolverT: each rank passes some of them, in any rows.
  struct Coo
  {
    std::vector<int> irn, jcn;
    std::vector<T> val;
  };

  // From the entries of an n_glob x n_glob matrix whose rows are numbered by rank (n_loc
  // local rows), the Schur variables (global 0-based indices, replicated, in their order in
  // S; none for a plain factorization), and rows: the local index of each local row in the
  // vectors of the solves (-1 for a row of none). Factors the matrix (collective).
  DirectSchurSolver(MPI_Comm comm, HYPRE_BigInt n_glob, int n_loc, Coo &&A,
                    const std::vector<HYPRE_BigInt> &schur_vars, std::vector<int> rows);
  ~DirectSchurSolver();

  // The entries of the next Refactor, in the order of the input (to be overwritten).
  std::vector<T> &Values() { return val; }

  // Factor the entries in Values() (collective).
  void Refactor();

  // The dense Schur complement on rank 0 (n_schur x n_schur, column-major, symmetric up to
  // rounding), handed over (until the next factorization).
  std::vector<T> ReleaseSchur();

  // Internal solves for a batch of vectors of the rows (collective): y_k solves A_II on the
  // internal (non-Schur) variables, 0 on the Schur variables and outside the rows. y_k may
  // alias x_k.
  void SolveInternal(const std::vector<const ComplexVector *> &X,
                     const std::vector<ComplexVector *> &Y);

  // Condensation of a batch of right-hand sides b (vectors of the rows) onto the Schur
  // variables: red = b_S - A_SI A_II^-1 b_I on rank 0 (n_schur x |B|, column-major). The
  // next call must be Expand of the same batch (collective).
  void Reduce(const std::vector<const ComplexVector *> &B, std::vector<T> &red);

  // Expansion of the batch of the last Reduce from the Schur part of its solution, u_S on
  // rank 0 (n_schur x |X|, column-major): x_I = A_II^-1 (b_I - A_IS u_S) and x_S = u_S, 0
  // outside the rows (collective).
  void Expand(const std::vector<T> &u, const std::vector<ComplexVector *> &X);

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
  std::unique_ptr<StrumpackSchur> nat;  // PROBE: STRUMPACK's Schur complement interface
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
