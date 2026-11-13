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
#include "linalg/lowertriangle.hpp"
#include "linalg/schursolver.hpp"
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

// MUMPS partial factorization with a Schur complement: factors the symmetric parent-space
// operator A (rows/cols outside the eliminated subsystem set to identity), or the
// subsystem alone, with the given Schur variables left unfactored, returning the dense
// Schur complement on rank 0 -- for the environment, S_E = A_GG - A_GE A_EE^-1 A_EG from
// one partial factorization instead of |Gamma| back-solves. The same factorization then
// solves the internal problem (A_EE, the Schur variables fixed at 0), which serves every
// other environment solve. T = double (DMUMPS) for a real symmetric A, std::complex<double>
// (ZMUMPS) for a complex symmetric A = Ar + i Ai. The right-hand sides and solutions of the
// solves stay distributed by rows; only the condensed (Schur) parts are on rank 0.
template <typename T>
class MumpsSchurSolverT : public SchurFactorization<T>
{
  static constexpr bool kComplex = !std::is_same_v<T, double>;
  using Struc = std::conditional_t<kComplex, ZMUMPS_STRUC_C, DMUMPS_STRUC_C>;

public:
  using typename SchurFactorization<T>::VecType;
  using Coo = typename SchurFactorization<T>::Coo;

  // schur_vars: global (0-based) true-DOF indices of the Schur variables, in the order the
  // Schur rows/columns should appear (replicated on all ranks); empty for a plain
  // factorization of A. blr_tol > 0 enables a block low-rank (BLR) factorization with that
  // relative accuracy (0: exact). With procs > 0, about that many ranks factor (every
  // ceil(size / procs)-th, rank 0 among them), the others sending them their entries: for a
  // small system, it avoids MUMPS's per-rank workspace on every rank (0: all ranks). With
  // refactor, the entries are kept so that Refactor can factor new values with the same
  // pattern, reusing the analysis. spd: A is positive definite, so it is factored without
  // numerical pivoting. ordering is the fill-reducing ordering ("Default": METIS); MUMPS
  // computes a Schur complement only after a sequential analysis, so with Schur variables a
  // parallel ordering is replaced by its sequential counterpart.
  MumpsSchurSolverT(const mfem::HypreParMatrix &A,
                    const std::vector<HYPRE_BigInt> &schur_vars, double blr_tol = 0.0,
                    int procs = 0, bool refactor = false, bool spd = true,
                    SymbolicFactorization ordering = SymbolicFactorization::DEFAULT);

  // From the local lower-triangle entries of A (n_loc local rows, in rank order). With
  // rows, the vectors of the solves are of a larger (parent) space: rows[i] is the local
  // index of local row i in them, or -1 for a row of none (the other entries of a solution
  // are 0).
  MumpsSchurSolverT(MPI_Comm comm, HYPRE_BigInt n_glob, int n_loc, Coo &&A,
                    const std::vector<HYPRE_BigInt> &schur_vars, double blr_tol = 0.0,
                    int procs = 0, bool refactor = false, bool spd = true,
                    std::vector<int> rows = {},
                    SymbolicFactorization ordering = SymbolicFactorization::DEFAULT);

  ~MumpsSchurSolverT() override;

  // With refactor (the factorization releases the entries otherwise).
  std::vector<T> &Values() override { return val; }

  // With the sparsity pattern and the ordering of the analysis.
  void Refactor() override;

  std::vector<T> ReleaseSchur() override;

  void SolveInternal(const std::vector<const VecType *> &X,
                     const std::vector<VecType *> &Y) override;

  void Reduce(const std::vector<const VecType *> &B, std::vector<T> &red) override;

  void Expand(const std::vector<T> &u, const std::vector<VecType *> &X) override;

private:
  void Init(const std::vector<HYPRE_BigInt> &schur_vars);
  void Factor();
  void Check(const char *phase) const;
  void GatherEntries();
  // The largest diagonal entry of the local values (diag_pos), over the ranks: the scale of
  // the BLR tolerance.
  double DiagonalScale() const;
  void Call();
  // The right-hand sides distributed by rows, room for the distributed solution, and the
  // solution on the local rows.
  void SetRhs(const std::vector<const VecType *> &X);
  void SetSolution(int nb);
  void ScatterSolution(const std::vector<VecType *> &Y);

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
  int rank = 0, procs = 0;
  bool refactor = false, spd = true, factored = false;
  SymbolicFactorization ordering = SymbolicFactorization::DEFAULT;
  bool active = true;  // this rank takes part in the factorization
  // With fewer ranks factoring: the ranks factoring, and the group of each (one factoring
  // rank and the ranks sending it their entries, with the counts and offsets of their
  // entries on the factoring rank).
  MPI_Comm sub = MPI_COMM_NULL, grp = MPI_COMM_NULL;
  std::vector<int> grp_cnt, grp_disp, grp_row_cnt, grp_row_disp;
  HYPRE_BigInt n_glob;
  int n_loc, n_schur;
  std::vector<int> row_cnt, row_disp, rows;
  int parent = 0;  // the vectors of the solves are of the parent space of rows
  std::vector<MUMPS_INT> irn, jcn, listvar, irhs_loc, isol_loc;
  std::vector<T> val, grp_val, schur, redrhs, rhs_loc, sol_loc;
  int reduced = 0;  // right-hand sides of the last Reduce, pending Expand
  double blr_tol = 0.0, scale = 1.0;
  std::vector<std::size_t> diag_pos;  // the diagonal entries of the local values (BLR)
  Struc id{};
};

using MumpsSchurSolver = MumpsSchurSolverT<double>;
using ComplexMumpsSchurSolver = MumpsSchurSolverT<std::complex<double>>;

}  // namespace palace

#endif

#endif  // PALACE_LINALG_MUMPS_SCHUR_HPP
