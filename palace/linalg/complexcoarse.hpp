// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_COMPLEX_COARSE_HPP
#define PALACE_LINALG_COMPLEX_COARSE_HPP

#include <functional>
#include <memory>
#include <mpi.h>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Sparse direct coarse solver for a complex-valued operator Ar + i Ai which is factored in
// complex arithmetic (Solver.Linear.ComplexCoarseSolve), unless its imaginary part is zero
// (see HasZeroImag), as for a model without losses: then Ar is factored in real arithmetic,
// which needs half of the memory, and the same factors solve the real and imaginary parts
// of the right-hand side.
//
// The choice is made for each operator, so it follows the operators of a frequency sweep.
// Only the solver in use exists: it is built when it is first needed, and the other one is
// destroyed, with its factors, before the matrix is factored.
//
class ComplexCoarseSolver : public Solver<ComplexOperator>
{
public:
  using SolverFactory = std::function<std::unique_ptr<Solver<ComplexOperator>>()>;

private:
  MPI_Comm comm;

  // Solver for a real-valued matrix (an MfemWrapperSolver, which factors the real part of
  // an operator without imaginary part) and solver in complex arithmetic.
  SolverFactory make_real_pc, make_complex_pc;
  std::unique_ptr<Solver<ComplexOperator>> real_pc, complex_pc;

  // Whether to print which matrix is factored, when the solver in use is built.
  bool print;

public:
  ComplexCoarseSolver(MPI_Comm comm, SolverFactory make_real_pc,
                      SolverFactory make_complex_pc, bool print);

  // The solver in use: at most one of the two is not null.
  Solver<ComplexOperator> *GetRealSolver() { return real_pc.get(); }
  Solver<ComplexOperator> *GetComplexSolver() { return complex_pc.get(); }

  void SetInitialGuess(bool guess) override;

  void SetOperator(const ComplexOperator &op) override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_COMPLEX_COARSE_HPP
