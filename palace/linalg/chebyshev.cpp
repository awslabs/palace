// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "chebyshev.hpp"

#include <vector>
#include <general/forall.hpp>

namespace palace
{

ChebyshevSmoother::ChebyshevSmoother(int smooth_it, int poly_order)
  : mfem::Solver(), pc_it(smooth_it), order(poly_order), A(nullptr)
{
}

void ChebyshevSmoother::SetOperator(const Operator &op)
{
  const auto *PtAP = dynamic_cast<const ParOperator *>(&op);
  MFEM_VERIFY(PtAP, "ChebyshevSmoother requires a ParOperator operator!");
  A = PtAP;

  height = A->Height();
  width = A->Width();
  r.SetSize(height);
  d.SetSize(height);
  dinv.SetSize(height);
  A->AssembleDiagonal(dinv);
  dinv.Reciprocal();

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. See
  // mfem::OperatorChebyshevSmoother or Adams et al., Parallel multigrid smoothing:
  // polynomial versus Gauss-Seidel, JCP (2003).
  DiagonalOperator Dinv(dinv);
  SymmetricProductOperator DinvA(Dinv, *A);
  lambda_max = 1.1 * linalg::SpectralNorm(A->GetComm(), DinvA, false);
}

void ChebyshevSmoother::Mult(const Vector &x, Vector &y) const
{
  // Apply smoother: y = y + p(A) (x - A y) .
  for (int it = 0; it < pc_it; it++)
  {
    if (iterative_mode || it > 0)
    {
      A->Mult(y, r);
      subtract(x, r, r);
    }
    else
    {
      r = x;
      y = 0.0;
    }

    // 4th-kind Chebyshev smoother, from Phillips and Fischer or Lottes (with k -> k + 1
    // shift due to 1-based indexing).
    {
      const int N = height;
      const auto *DI = dinv.Read();
      const auto *R = r.Read();
      auto *D = d.ReadWrite();
      mfem::forall(N, [=] MFEM_HOST_DEVICE(int i)
                   { D[i] = 4.0 / (3.0 * lambda_max) * DI[i] * R[i]; });
    }
    for (int k = 1; k < order; k++)
    {
      A->AddMult(d, r, -1.0);
      {
        const int N = height;
        const double sd = (2.0 * k - 1.0) / (2.0 * k + 3.0);
        const double sr = (8.0 * k + 4.0) / ((2.0 * k + 3.0) * lambda_max);
        const auto *DI = dinv.Read();
        const auto *R = r.Read();
        auto *Y = y.ReadWrite();
        auto *D = d.ReadWrite();
        mfem::forall(N,
                     [=] MFEM_HOST_DEVICE(int i)
                     {
                       Y[i] += D[i];
                       D[i] = sd * D[i] + sr * DI[i] * R[i];
                     });
      }
    }
    y += d;
  }
}

}  // namespace palace
