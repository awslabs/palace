// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_STRUMPACK_COMPLEX_HPP
#define PALACE_LINALG_STRUMPACK_COMPLEX_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_STRUMPACK)

#include <complex>
#include <memory>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "utils/labels.hpp"

namespace strumpack
{

template <typename scalar_t, typename integer_t>
class SparseSolverMPIDist;
template <typename factor_t, typename refine_t, typename integer_t>
class SparseSolverMixedPrecisionMPIDist;

}  // namespace strumpack

namespace palace
{

//
// STRUMPACK sparse direct solver for complex-valued operators in complex arithmetic, which
// halves the memory of the factorization compared to the equivalent real-valued system of
// twice the size used by the other wrappers of sparse direct solvers (see
// MfemWrapperSolver). Optionally, the factorization is computed in single precision, with
// steps of iterative refinement in double precision.
//
class StrumpackComplexSolver : public Solver<ComplexOperator>
{
private:
  MPI_Comm comm;
  SymbolicFactorization reorder;
  bool single;
  int refine_its;
  int print;
  std::unique_ptr<strumpack::SparseSolverMPIDist<std::complex<double>, int>> solver_d;
  std::unique_ptr<strumpack::SparseSolverMixedPrecisionMPIDist<std::complex<float>,
                                                               std::complex<double>, int>>
      solver_s;
  mutable std::vector<std::complex<double>> b, x;

  // Sparsity pattern of the last factored matrix, for the reuse of the reordering.
  std::vector<int> I_prev, J_prev;

public:
  // With single, refine_its steps of iterative refinement in double precision follow the
  // single precision solve.
  StrumpackComplexSolver(MPI_Comm comm, SymbolicFactorization reorder, bool single,
                         int refine_its, int print);
  ~StrumpackComplexSolver() override;

  void SetOperator(const ComplexOperator &op) override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_STRUMPACK_COMPLEX_HPP
