// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SCHUR_SOLVER_HPP
#define PALACE_LINALG_SCHUR_SOLVER_HPP

#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

//
// A sparse symmetric matrix A = [A_II, A_IS; A_SI, A_SS] factored for its Schur complement
// S = A_SS - A_SI A_II^-1 A_IS onto some of its variables (the Schur variables), and for
// solves with its internal block A_II. A partial factorization (MumpsSchurSolverT, and
// DirectSchurSolverT with STRUMPACK's Schur complement interface) gives S directly;
// DirectSchurSolverT otherwise forms it from one solve with A_II per Schur variable.
//
template <typename T>
class SchurSolver
{
public:
  using VecType =
      std::conditional_t<std::is_same_v<T, double>, mfem::Vector, ComplexVector>;

  virtual ~SchurSolver() = default;

  // The rows of the dense Schur complement (symmetric up to rounding) of the Schur
  // variables this rank holds, in the order of their positions among them
  // (|SchurRowGather::Rows()| x n_schur, row-major), handed over (collective).
  virtual std::vector<T> ReleaseSchurRows() = 0;

  // Internal solves for a batch of vectors (collective): y_k solves A_II on the internal
  // variables, with 0 on the Schur variables, and on the entries of the vectors outside the
  // rows of A. y_k may alias x_k.
  virtual void SolveInternal(const std::vector<const VecType *> &X,
                             const std::vector<VecType *> &Y) = 0;

  // Condensations of a batch of vectors onto the Schur variables (collective): g = x_S -
  // A_SI A_II^-1 x_I on rank 0 (n_schur x |X|, column-major). A partial factorization
  // condenses by its forward elimination only, half of an internal solve.
  virtual void Condense(const std::vector<const VecType *> &X, std::vector<T> &g) = 0;
};

//
// A row of w entries (double or std::complex<double>), so that the counts of a row exchange
// are rows (n w entries can overflow an int).
//
template <typename T>
class RowTypeT
{
public:
  explicit RowTypeT(int w)
  {
    MPI_Type_contiguous(w, mpi::DataType<T>(), &type);
    MPI_Type_commit(&type);
  }
  ~RowTypeT() { MPI_Type_free(&type); }
  RowTypeT(const RowTypeT &) = delete;
  RowTypeT &operator=(const RowTypeT &) = delete;
  MPI_Datatype type;
};

using RowType = RowTypeT<double>;

//
// The rows of the Schur variables that each rank holds, to gather their entries of a batch
// of vectors on rank 0 in the order of the Schur variables, and to move the rows of a Schur
// complement between rank 0 and the ranks holding them.
//
class SchurRowGather
{
public:
  // schur_vars: the global (0-based) rows of the Schur variables (replicated); row0 and
  // n_loc: the rows of this rank. Collective.
  SchurRowGather(MPI_Comm comm, HYPRE_BigInt row0, int n_loc,
                 const std::vector<HYPRE_BigInt> &schur_vars);

  // With the rows numbered by rank, n_loc on this one. Collective.
  SchurRowGather(MPI_Comm comm, int n_loc, const std::vector<HYPRE_BigInt> &schur_vars);

  // The local rows of the Schur variables on this rank, and their positions among them.
  const std::vector<int> &Rows() const { return row; }
  const std::vector<int> &Positions() const { return pos; }

  // Gather loc (the entries in Rows() of nb vectors, |Rows()| x nb, column-major) on rank 0
  // into g (n_schur x nb, column-major; not used on the other ranks). Collective.
  void Gather(const std::vector<double> &loc, int nb, double *g) const;

  // The matrix on rank 0 (n_schur x n_schur, column-major; empty on the other ranks) from
  // its rows on the ranks holding them (|Rows()| x n_schur, row-major), which it receives
  // one rank at a time. Collective.
  template <typename T>
  std::vector<T> GatherSchur(const std::vector<T> &rows) const;

  // The rows on each rank (|Rows()| x n_schur, row-major) of a matrix on rank 0 (n_schur x
  // n_schur, column-major), which it sends one rank at a time. Collective.
  template <typename T>
  std::vector<T> ScatterSchur(const std::vector<T> &full) const;

private:
  MPI_Comm comm;
  int n_schur;
  // On rank 0, the positions of the Schur rows of every rank (in rank order), with their
  // counts and offsets.
  std::vector<int> row, pos, all_pos, pos_cnt, pos_dsp;
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

  // The dense Schur complement on rank 0 (n_schur x n_schur, column-major, symmetric up to
  // rounding), handed over instead of its rows (collective).
  virtual std::vector<T> ReleaseSchur() = 0;

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

}  // namespace palace

#endif  // PALACE_LINALG_SCHUR_SOLVER_HPP
