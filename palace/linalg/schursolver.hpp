// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SCHUR_SOLVER_HPP
#define PALACE_LINALG_SCHUR_SOLVER_HPP

#include <functional>
#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include "linalg/vector.hpp"

namespace palace
{

//
// A sparse symmetric matrix A = [A_II, A_IS; A_SI, A_SS] factored for its Schur complement
// S = A_SS - A_SI A_II^-1 A_IS onto some of its variables (the Schur variables), and for
// solves with its internal block A_II. A partial factorization (MumpsSchurSolver,
// StrumpackSchurSolver) gives S directly; BackSolveSchurSolver forms it from one solve
// with A_II per Schur variable.
//
template <typename T>
class SchurSolver
{
public:
  using VecType =
      std::conditional_t<std::is_same_v<T, double>, mfem::Vector, ComplexVector>;

  virtual ~SchurSolver() = default;

  // The dense Schur complement on rank 0 (n_schur x n_schur, column-major, symmetric up to
  // rounding), handed over (collective).
  virtual std::vector<T> ReleaseSchur() = 0;

  // Internal solves for a batch of vectors (collective): y_k solves A_II on the internal
  // variables, with 0 on the Schur variables. y_k may alias x_k.
  virtual void SolveInternal(const std::vector<const VecType *> &X,
                             const std::vector<VecType *> &Y) = 0;
};

//
// A SchurSolver from the entries of the matrix, which it factors again with new values of
// the same pattern, and which condenses right-hand sides onto the Schur variables and
// expands their solutions.
//
template <typename T>
class SchurFactorization : public SchurSolver<T>
{
public:
  using typename SchurSolver<T>::VecType;

  // Lower-triangle COO entries (1-based global indices) of the matrix.
  struct Coo
  {
    std::vector<int> irn, jcn;
    std::vector<T> val;
  };

  // The entries of the next Refactor, in the order of the input (to be overwritten).
  virtual std::vector<T> &Values() = 0;

  // Factor the entries in Values() (collective).
  virtual void Refactor() = 0;

  // Condensation of a batch of right-hand sides onto the Schur variables (collective):
  // red = b_S - A_SI A_II^-1 b_I on rank 0 (n_schur x |B|, column-major). The next call
  // must be Expand of the same batch.
  virtual void Reduce(const std::vector<const VecType *> &B, std::vector<T> &red) = 0;

  // Expansion of the batch of the last Reduce from the Schur part of its solution, u_S on
  // rank 0 (n_schur x |X|, column-major): x_I = A_II^-1 (b_I - A_IS u_S) and x_S = u_S
  // (collective).
  virtual void Expand(const std::vector<T> &u, const std::vector<VecType *> &X) = 0;
};

//
// The Schur complement from solves with A_II, one per Schur variable, in batches, by a
// given solver of the internal problem (a factorization of A_II, or an iterative solver),
// whose variables are the internal ones: the other rows of A are not used.
//
class BackSolveSchurSolver : public SchurSolver<double>
{
public:
  // Solves A_II y_k = x_k on the internal variables for a batch, with y_k = 0 elsewhere
  // (x_k there is not used). y_k may alias x_k.
  using InternalSolver = std::function<void(const std::vector<const mfem::Vector *> &,
                                            const std::vector<mfem::Vector *> &)>;

  // A is referenced (not copied). schur_vars: the global (0-based) rows of the Schur
  // variables, in the order of their rows and columns in S (replicated). block: the number
  // of solves in a batch.
  BackSolveSchurSolver(const mfem::HypreParMatrix &A,
                       const std::vector<HYPRE_BigInt> &schur_vars, InternalSolver solve,
                       int block);

  std::vector<double> ReleaseSchur() override;

  void SolveInternal(const std::vector<const mfem::Vector *> &X,
                     const std::vector<mfem::Vector *> &Y) override;

private:
  const mfem::HypreParMatrix &A;
  InternalSolver solve;
  int n_schur, block;
  // The local rows of the Schur variables and their positions among them, and on rank 0,
  // the positions of those of every rank (in rank order) with their counts and offsets.
  std::vector<int> schur_row, schur_pos, all_pos, pos_cnt, pos_dsp;
};

}  // namespace palace

#endif  // PALACE_LINALG_SCHUR_SOLVER_HPP
