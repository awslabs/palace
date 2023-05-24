// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_JACOBI_SMOOTHER_HPP
#define PALACE_LINALG_JACOBI_SMOOTHER_HPP

#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Simple Jacobi smoother using the diagonal vector from Operator::AssembleDiagonal(),
// which allows for (approximate) diagonal construction for matrix-free operators.
//
class JacobiSmoother : public mfem::Solver
{
private:
  // Inverse diagonal scaling of the operator.
  Vector dinv;

public:
  JacobiSmoother() : mfem::Solver() {}

  void SetOperator(const Operator &op) override;

  void Mult(const Vector &x, Vector &y) const override;

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

}  // namespace palace

#endif  // PALACE_LINALG_JACOBI_SMOOTHER_HPP
