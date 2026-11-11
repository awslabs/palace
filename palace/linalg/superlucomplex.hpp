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
// (pzgssvx), which halves the memory of the factorization compared to the equivalent
// real-valued system of twice the size used by MfemWrapperSolver. Only for 2D process
// grids. With robust_pivoting, the accuracy of each factorization is checked as in
// SuperLUSolver, and the matrix factored again with tiny pivot replacement and iterative
// refinement if it is inaccurate.
//
class SuperLUComplexSolver : public Solver<ComplexOperator>
{
private:
  MPI_Comm comm;
  SymbolicFactorization reorder;
  bool reorder_reuse, robust_pivoting, refine;
  int print;

  // SuperLU_DIST data structures.
  struct Data;
  std::unique_ptr<Data> data;

  // Factor the current matrix (with a solve for a random right-hand side, also used to
  // check the accuracy of the factorization). Returns the relative residual.
  double Factor(const ComplexOperator &op);

public:
  SuperLUComplexSolver(MPI_Comm comm, SymbolicFactorization reorder, bool reorder_reuse,
                       int print, bool robust_pivoting = false);
  ~SuperLUComplexSolver() override;

  void SetOperator(const ComplexOperator &op) override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_SUPERLU_COMPLEX_HPP
