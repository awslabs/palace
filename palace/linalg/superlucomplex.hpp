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
// (pzgssvx). Its L and U factors need up to about 40% less memory than those of the
// equivalent real-valued system of twice the size used by MfemWrapperSolver (the saving
// grows with the problem size). Only for 2D process grids (SuperLUSolver is used with a 3D
// process grid).
//
class SuperLUComplexSolver : public Solver<ComplexOperator>
{
private:
  MPI_Comm comm;
  bool reorder_reuse;

  // SuperLU_DIST data structures.
  struct Data;
  std::unique_ptr<Data> data;

public:
  SuperLUComplexSolver(MPI_Comm comm, SymbolicFactorization reorder, bool reorder_reuse,
                       int print);
  ~SuperLUComplexSolver() override;

  void SetOperator(const ComplexOperator &op) override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_SUPERLU_COMPLEX_HPP
