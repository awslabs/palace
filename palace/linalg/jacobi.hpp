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
// Simple Jacobi smoother using a provided diagonal vector, usually the output of
// AssembleDiagonal() which allows for (approximatE) diagonal construction for matrix-free
// operators.
//
class JacobiSmoother : public mfem::Solver
{
private:
  // Inverse diagonal scaling of the operator.
  Vector dinv;

public:
  JacobiSmoother() : mfem::Solver() {}

  void SetOperator(const Operator &op) override
  {
    MFEM_ABORT("JacobiSmoother requires a ParOperator operator!");
  }
  void SetOperator(const ParOperator &op);

  void Mult(const Vector &x, Vector &y) const override;

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

}  // namespace palace

#endif  // PALACE_LINALG_JACOBI_SMOOTHER_HPP
