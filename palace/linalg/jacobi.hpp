// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_JACOBI_SMOOTHER_HPP
#define PALACE_LINALG_JACOBI_SMOOTHER_HPP

#include <mfem.hpp>

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
  mfem::Vector dinv;

public:
  JacobiSmoother(const mfem::Vector &diag);

  void SetOperator(const mfem::Operator &op) override {}

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override;

  void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const override { Mult(x, y); }
};

}  // namespace palace

#endif  // PALACE_LINALG_JACOBI_SMOOTHER_HPP
