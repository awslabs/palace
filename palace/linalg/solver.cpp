// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "solver.hpp"

namespace palace
{

template <>
void WrapperSolver<Operator>::SetOperator(const Operator &op)
{
  pc->SetOperator(op);
}

template <>
void WrapperSolver<ComplexOperator>::SetOperator(const ComplexOperator &op)
{
  MFEM_VERIFY(op.IsReal() && op.HasReal(),
              "WrapperSolver::SetOperator requires an operator which is purely real for "
              "mfem::Solver!");
  pc->SetOperator(*op.Real());
}

template <>
void WrapperSolver<Operator>::Mult(const Vector &x, Vector &y) const
{
  pc->Mult(x, y);
}

template <>
void WrapperSolver<ComplexOperator>::Mult(const ComplexVector &x, ComplexVector &y) const
{
  mfem::Array<const Vector *> X(2);
  mfem::Array<Vector *> Y(2);
  X[0] = &x.Real();
  X[1] = &x.Imag();
  Y[0] = &y.Real();
  Y[1] = &y.Imag();
  pc->ArrayMult(X, Y);
}

}  // namespace palace
