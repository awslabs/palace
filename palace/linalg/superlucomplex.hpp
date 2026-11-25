// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SUPERLU_COMPLEX_HPP
#define PALACE_LINALG_SUPERLU_COMPLEX_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_SUPERLU)

#include <memory>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "utils/labels.hpp"

namespace palace
{

//
// SuperLU_DIST sparse direct solver for complex-valued operators in complex arithmetic
// (pzgssvx, or pzgssvx3d with a 3D process grid). Its L and U factors need up to about 40%
// less memory than those of the equivalent real-valued system of twice the size used by
// MfemWrapperSolver (the saving grows with the problem size). The process grid is chosen
// as for SuperLUSolver, and small entries are dropped as described for ComplexLocalCSR.
// With robust_pivoting, the accuracy of each factorization is checked as for SuperLUSolver,
// and the matrix is factored again with tiny pivot replacement and iterative refinement if
// it is inaccurate (for a 2D process grid).
//
class SuperLUComplexSolver : public Solver<ComplexOperator>
{
private:
  MPI_Comm comm;
  bool drop_small_entries, reorder_reuse, robust_pivoting, refine;
  int print;

  // SuperLU_DIST data structures.
  struct Data;
  std::unique_ptr<Data> data;

  // Relative residual of a solve with the factorization, for a random solution.
  double FactorizationResidual(const ComplexOperator &op) const;

public:
  SuperLUComplexSolver(MPI_Comm comm, SymbolicFactorization reorder, bool use_3d,
                       bool drop_small_entries, bool reorder_reuse, int print,
                       bool robust_pivoting = false);
  ~SuperLUComplexSolver() override;

  void SetOperator(const ComplexOperator &op) override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_SUPERLU_COMPLEX_HPP
