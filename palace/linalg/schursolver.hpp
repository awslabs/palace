// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SCHUR_SOLVER_HPP
#define PALACE_LINALG_SCHUR_SOLVER_HPP

#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include "linalg/blockcyclic.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// A sparse symmetric matrix A = [A_II, A_IS; A_SI, A_SS] factored for its Schur complement
// S = A_SS - A_SI A_II^-1 A_IS onto some of its variables (the Schur variables), and for
// solves with its internal block A_II, from the entries of the matrix, which it factors
// again with new values of the same pattern, and which condenses right-hand sides onto the
// Schur variables and expands their solutions. A partial factorization (MumpsSchurSolverT,
// and DirectSchurSolverT with STRUMPACK's Schur complement interface) gives S directly;
// DirectSchurSolverT otherwise forms it from one solve with A_II per Schur variable.
//
template <typename T>
class SchurFactorization
{
public:
  using VecType =
      std::conditional_t<std::is_same_v<T, double>, mfem::Vector, ComplexVector>;

  virtual ~SchurFactorization() = default;

  // Lower-triangle COO entries (1-based global indices) of the matrix.
  struct Coo
  {
    std::vector<int> irn, jcn;
    std::vector<T> val;
  };

  // The entries of the next Refactor, in the order of the input (to be overwritten).
  virtual std::vector<T> &Values() = 0;

  // The dense Schur complement on rank 0 (n_schur x n_schur, column-major, symmetric up to
  // rounding), handed over (collective).
  virtual std::vector<T> ReleaseSchur() = 0;

  // The local part on a grid of the Schur complement, handed over (collective): by default,
  // that of rank 0 scattered to the grid.
  virtual std::vector<T> ReleaseSchurBlocks(const BlockCyclicGrid &grid)
  {
    return grid.ScatterFromRoot(ReleaseSchur());
  }

  // Factor the entries in Values() (collective).
  virtual void Refactor() = 0;

  // Internal solves for a batch of vectors (collective): y_k solves A_II on the internal
  // variables, with 0 on the Schur variables, and on the entries of the vectors outside the
  // rows of A. y_k may alias x_k.
  virtual void SolveInternal(const std::vector<const VecType *> &X,
                             const std::vector<VecType *> &Y) = 0;

  // Condensations of a batch of vectors onto the Schur variables (collective): g = x_S -
  // A_SI A_II^-1 x_I on rank 0 (n_schur x |X|, column-major). A partial factorization
  // condenses by its forward elimination only, half of an internal solve.
  virtual void Condense(const std::vector<const VecType *> &X, std::vector<T> &g) = 0;

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
