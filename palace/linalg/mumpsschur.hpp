// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_MUMPS_SCHUR_HPP
#define PALACE_LINALG_MUMPS_SCHUR_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_MUMPS)

#include <complex>
#include <type_traits>
#include <vector>
#include <dmumps_c.h>
#include <zmumps_c.h>
#include "linalg/vector.hpp"

namespace palace
{

// Lower-triangle COO pattern (1-based global indices, one entry per position, columns
// ascending in each row) of the local rows of the union of the given symmetric matrices
// (same row distribution; null ones skipped) and the diagonal of the given local rows. The
// entries of each matrix in pattern order are returned in vals (empty for a null one), and
// the first pattern entry of each local row in row_ptr (n_loc + 1).
void LowerTrianglePattern(const std::vector<const mfem::HypreParMatrix *> &A,
                          const mfem::Array<int> &diag_rows, std::vector<MUMPS_INT> &irn,
                          std::vector<MUMPS_INT> &jcn,
                          std::vector<std::vector<double>> &vals,
                          std::vector<int> &row_ptr);

// val += a X on a pattern from LowerTrianglePattern, for a matrix X with entries inside the
// pattern (or zero).
void AddToPattern(const mfem::HypreParMatrix &X, double a, const std::vector<int> &row_ptr,
                  const std::vector<MUMPS_INT> &jcn, std::vector<double> &val);
void AddToPattern(const mfem::HypreParMatrix &X, std::complex<double> a,
                  const std::vector<int> &row_ptr, const std::vector<MUMPS_INT> &jcn,
                  std::vector<std::complex<double>> &val);

// MUMPS partial factorization with a Schur complement: factors the symmetric parent-space
// operator A (rows/cols outside the eliminated subsystem set to identity) with the given
// Schur variables left unfactored, returning the dense Schur complement on rank 0 -- for
// the environment, S_E = A_GG - A_GE A_EE^-1 A_EG from one partial factorization instead of
// |Gamma| back-solves. The same factorization then solves the internal problem (A_EE, the
// Schur variables fixed at 0), which serves every other environment solve. T = double
// (DMUMPS) for a real symmetric A, std::complex<double> (ZMUMPS) for a complex symmetric
// A = Ar + i Ai.
template <typename T>
class MumpsSchurSolverT
{
  static constexpr bool kComplex = !std::is_same_v<T, double>;
  using Struc = std::conditional_t<kComplex, ZMUMPS_STRUC_C, DMUMPS_STRUC_C>;

public:
  using VecType = std::conditional_t<kComplex, ComplexVector, mfem::Vector>;

  // Local lower-triangle COO entries (1-based global indices) of the local rows of an
  // n_glob x n_glob symmetric matrix.
  struct Coo
  {
    std::vector<MUMPS_INT> irn, jcn;
    std::vector<T> val;
  };

  // schur_vars: global (0-based) true-DOF indices of the Schur variables, in the order the
  // Schur rows/columns should appear (replicated on all ranks); empty for a plain
  // factorization of A. blr_tol > 0 enables a block low-rank (BLR) factorization with that
  // relative accuracy (0: exact). With serial, rank 0 factors alone (for a small system, it
  // avoids MUMPS's per-rank workspace). With refactor, the entries are kept so that
  // Refactor can factor new values with the same pattern, reusing the analysis. spd: A is
  // positive definite, so the BLR factorization skips numerical pivoting.
  MumpsSchurSolverT(const mfem::HypreParMatrix &A,
                    const std::vector<HYPRE_BigInt> &schur_vars, double blr_tol = 0.0,
                    bool serial = false, bool refactor = false, bool spd = true);

  // From the local lower-triangle entries of A (n_loc local rows, in rank order).
  MumpsSchurSolverT(MPI_Comm comm, HYPRE_BigInt n_glob, int n_loc, Coo &&A,
                    const std::vector<HYPRE_BigInt> &schur_vars, double blr_tol = 0.0,
                    bool serial = false, bool refactor = false, bool spd = true);

  ~MumpsSchurSolverT();

  // With refactor: the entries of the next Refactor in the order of the input (to be
  // overwritten), and their columns.
  std::vector<T> &Values() { return val; }
  const std::vector<MUMPS_INT> &Columns() const { return jcn; }

  // Factor the entries in Values() with the sparsity pattern of the analysis (collective).
  void Refactor();

  // Dense Schur complement on rank 0 (n_schur x n_schur, column-major, symmetric).
  const std::vector<T> &Schur() const { return schur; }

  // Internal solves for a batch of parent true-DOF vectors (collective): y_k solves A_11 on
  // the internal (non-Schur) variables, with 0 on the Schur variables. y_k may alias x_k.
  void SolveInternal(const std::vector<const VecType *> &X,
                     const std::vector<VecType *> &Y);

  // Condensation of a batch of right-hand sides onto the Schur variables (collective):
  // red = b_S - A_SI A_II^-1 b_I on rank 0 (n_schur x |B|, column-major). The next call
  // must be Expand of the same batch.
  void Reduce(const std::vector<const VecType *> &B, std::vector<T> &red);

  // Expansion of the batch of the last Reduce from the Schur part of its solution, u_S on
  // rank 0 (n_schur x |X|, column-major): x_I = A_II^-1 (b_I - A_IS u_S), with u_S on the
  // Schur variables (collective).
  void Expand(const std::vector<T> &u, const std::vector<VecType *> &X);

private:
  void Init(const std::vector<HYPRE_BigInt> &schur_vars);
  void Factor();
  void Check(const char *phase) const;
  void Call();
  void Gather(const std::vector<const VecType *> &X);
  void Scatter(const std::vector<VecType *> &Y);

  // MUMPS's entry type (std::complex<double> has the layout of ZMUMPS_COMPLEX).
  static auto *Entries(T *v)
  {
    if constexpr (kComplex)
    {
      return reinterpret_cast<ZMUMPS_COMPLEX *>(v);
    }
    else
    {
      return v;
    }
  }

  MPI_Comm comm;
  int rank = 0;
  bool serial = false, refactor = false, spd = true, factored = false;
  bool active = true;  // this rank takes part in the factorization
  HYPRE_BigInt n_glob;
  int n_loc, n_schur;
  std::vector<int> row_cnt, row_disp;
  std::vector<MUMPS_INT> irn, jcn, listvar;
  std::vector<T> val, schur, rhs, redrhs;
  int reduced = 0;  // right-hand sides of the last Reduce, pending Expand
  double blr_tol = 0.0, scale = 1.0;
  Struc id{};
};

using MumpsSchurSolver = MumpsSchurSolverT<double>;
using ComplexMumpsSchurSolver = MumpsSchurSolverT<std::complex<double>>;

}  // namespace palace

#endif

#endif  // PALACE_LINALG_MUMPS_SCHUR_HPP
