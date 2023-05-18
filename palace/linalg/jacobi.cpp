// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "jacobi.hpp"

#include <general/forall.hpp>

namespace palace
{

void JacobiSmoother::SetOperator(const Operator &op)
{
  height = op.Height();
  width = op.Width();
  dinv.SetSize(height);
  op.AssembleDiagonal(dinv);
  dinv.Reciprocal();
}

void JacobiSmoother::Mult(const Vector &x, Vector &y) const
{
  MFEM_ASSERT(!iterative_mode,
              "JacobiSmoother is not implemented for iterative_mode = true!");
  const int N = height;
  const auto *DI = dinv.Read();
  const auto *X = x.Read();
  auto *Y = y.Write();
  mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { Y[i] = DI[i] * X[i]; });
}

}  // namespace palace
