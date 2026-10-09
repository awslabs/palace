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
class SuperLUSolver : public mfem::Solver
{
private:
  MPI_Comm comm;
  std::unique_ptr<mfem::SuperLURowLocMatrix> A;
  mfem::SuperLUSolver solver;
  bool reorder_reuse, robust_pivoting, refine;
  int print;

  // Set the matrix to factor (the factorization is performed on the first solve).
  void SetMatrix(const mfem::HypreParMatrix &hA);

  // Relative residual of a solve with the factorization, for a random solution.
  double FactorizationResidual(const mfem::HypreParMatrix &hA) const;

public:
  // SuperLU_DIST uses static pivoting, for which the elimination is not stable in general
  // for indefinite matrices. For matrices which must be solved to full accuracy (such as
  // the PML subdomain matrices), robust_pivoting checks the accuracy of each factorization
  // and, if it is inaccurate, factors again with tiny pivot replacement and iterative
  // refinement (then kept for all subsequent factorizations). It also disables the
  // symmetric pattern option, which assumes a structurally symmetric matrix after the row
  // permutation. Solves then have a single right-hand side.
  SuperLUSolver(MPI_Comm comm, SymbolicFactorization reorder, bool use_3d,
                bool reorder_reuse, int print, bool robust_pivoting = false);
  SuperLUSolver(const IoData &iodata, MPI_Comm comm, int print)
    : SuperLUSolver(comm, iodata.solver.linear.sym_factorization,
                    iodata.solver.linear.superlu_3d, iodata.solver.linear.reorder_reuse,
                    print)
  {
  }

  mfem::SuperLUSolver &GetSolver() { return solver; }

  void SetOperator(const Operator &op) override;

  void Mult(const Vector &x, Vector &y) const override { solver.Mult(x, y); }
  void ArrayMult(const mfem::Array<const Vector *> &X,
                 mfem::Array<Vector *> &Y) const override
  {
    solver.ArrayMult(X, Y);
  }
  void MultTranspose(const Vector &x, Vector &y) const override
  {
    solver.MultTranspose(x, y);
  }
  void ArrayMultTranspose(const mfem::Array<const Vector *> &X,
                          mfem::Array<Vector *> &Y) const override
  {
    solver.ArrayMultTranspose(X, Y);
  }

  void SetReorderReuse(bool reorder_reuse_) { reorder_reuse = reorder_reuse_; }
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_SUPERLU_HPP
