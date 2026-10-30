// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SUPERLU_HPP
#define PALACE_LINALG_SUPERLU_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_SUPERLU)

#include <memory>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/iodata.hpp"

namespace palace
{

//
// A wrapper for the SuperLU_DIST direct solver package.
//
// Every operator is factored from scratch: ReorderingReuse is not honored, so each
// SetOperator repeats the ordering and symbolic factorization. Reusing the row permutation
// of a previous matrix (SamePattern_SameRowPerm) gave an exactly zero pivot for the
// nonsingular coarse matrices of a Quasi-Newton eigenvalue solve with a moving shift, and
// reusing only the ordering (SamePattern) failed on them too.
//
// The native solver is bound to an operator only in the first solve after SetOperator:
// MFEM's solver may only be rebound after it has factored its operator, and destroying one
// bound to an unfactored operator is unsafe. A failed factorization or solve is final,
// and every operator must have the layout (local rows, first row) of the first.
//
class SuperLUSolver : public mfem::Solver
{
private:
  MPI_Comm comm;
  // The operator of the last SetOperator until the next solve, and the operator the native
  // solver is bound to. Declared before solver, whose destructor reads its bound operator.
  mutable std::unique_ptr<mfem::SuperLURowLocMatrix> pending, bound;
  mutable mfem::SuperLUSolver solver;
  mutable bool failed = false;
  mutable int nrhs = 0;  // The number of right-hand sides of every solve (MFEM's limit)
  HYPRE_BigInt layout[3] = {-1, -1, -1};  // Local rows, first row, global rows

  // Bind a pending operator, then run a native solve; a failure in either is final. The
  // native solve rejects another number of right-hand sides only after a binding has
  // released the previous factors, so that is checked first.
  template <typename F>
  void Solve(int num_rhs, F &&f) const
  {
    MFEM_VERIFY(!failed, "SuperLUSolver cannot be used after a failed factorization!");
    MFEM_VERIFY(num_rhs > 0 && (nrhs == 0 || nrhs == num_rhs),
                "SuperLUSolver requires the same positive number of right-hand sides for "
                "every solve!");
    nrhs = num_rhs;
    try
    {
      if (pending)
      {
        // The first binding, or MFEM's rebinding of a factored solver, which refactors.
        solver.SetOperator(*pending);
        bound = std::move(pending);
      }
      f();
    }
    catch (...)
    {
      failed = true;
      throw;
    }
  }

public:
  // The reorder_reuse argument is not used (see above).
  SuperLUSolver(MPI_Comm comm, SymbolicFactorization reorder, bool use_3d,
                bool reorder_reuse, int print);
  SuperLUSolver(const IoData &iodata, MPI_Comm comm, int print)
    : SuperLUSolver(comm, iodata.solver.linear.sym_factorization,
                    iodata.solver.linear.superlu_3d, iodata.solver.linear.reorder_reuse,
                    print)
  {
  }

  // The native solver, for its options only (it is configured in the constructor, and never
  // replaced). Its operator and solves belong to this class.
  mfem::SuperLUSolver &GetSolver() { return solver; }

  void SetOperator(const Operator &op) override;

  void Mult(const Vector &x, Vector &y) const override
  {
    Solve(1, [&] { solver.Mult(x, y); });
  }
  void ArrayMult(const mfem::Array<const Vector *> &X,
                 mfem::Array<Vector *> &Y) const override
  {
    Solve(X.Size(), [&] { solver.ArrayMult(X, Y); });
  }
  void MultTranspose(const Vector &x, Vector &y) const override
  {
    Solve(1, [&] { solver.MultTranspose(x, y); });
  }
  void ArrayMultTranspose(const mfem::Array<const Vector *> &X,
                          mfem::Array<Vector *> &Y) const override
  {
    Solve(X.Size(), [&] { solver.ArrayMultTranspose(X, Y); });
  }

  // Has no effect: every operator is factored from scratch (see above).
  void SetReorderReuse(bool) {}
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_SUPERLU_HPP
