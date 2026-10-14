// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_MUMPS_COMPLEX_HPP
#define PALACE_LINALG_MUMPS_COMPLEX_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_MUMPS)

#include <memory>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "utils/labels.hpp"

namespace palace
{

//
// MUMPS sparse direct solver for complex-valued operators in complex arithmetic (ZMUMPS),
// which halves the memory of the factorization compared to the equivalent real-valued
// system of twice the size used by MfemWrapperSolver. Symmetric matrices are factored as
// complex symmetric (not Hermitian) matrices.
//
class MumpsComplexSolver : public Solver<ComplexOperator>
{
private:
  MPI_Comm comm;
  MatrixSymmetry sym;
  SymbolicFactorization reorder;
  double blr_tol;
  bool reorder_reuse;
  int print;

  // MUMPS data structures.
  struct Data;
  std::unique_ptr<Data> data;

public:
  MumpsComplexSolver(MPI_Comm comm, MatrixSymmetry sym, SymbolicFactorization reorder,
                     double blr_tol, bool reorder_reuse, int print);
  ~MumpsComplexSolver() override;

  void SetOperator(const ComplexOperator &op) override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_MUMPS_COMPLEX_HPP
