// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "chebyshev.hpp"

#include <vector>
#include <mfem/general/forall.hpp>

namespace palace
{

ChebyshevSmoother::ChebyshevSmoother(int smooth_it, int poly_order)
  : mfem::Solver(), pc_it(smooth_it), order(poly_order), A(nullptr)
{
}

void ChebyshevSmoother::SetOperator(const ParOperator &op)
{
  A = &op;

  height = A->Height();
  width = A->Width();
  r.SetSize(height);
  d.SetSize(height);
  dinv.SetSize(height);
  A->AssembleDiagonal(dinv);
  // dinv.Reciprocal();    //XX TODO NEED MFEM PATCH

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. See
  // mfem::OperatorChebyshevSmoother or Adams et al., Parallel multigrid smoothing:
  // polynomial versus Gauss-Seidel, JCP (2003).
  DiagonalOperator Dinv(dinv);
  SymmetricProductOperator DinvA(Dinv, *A);
  lambda_max = 1.1 * linalg::SpectralNorm(A->GetComm(), DinvA, false);
}

void ChebyshevSmoother::Mult(const mfem::Vector &x, mfem::Vector &y) const
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

// void ChebyshevSmoother::ArrayMult(const mfem::Array<const mfem::Vector *> &X,
//                                   mfem::Array<mfem::Vector *> &Y) const
// {
//   // Initialize.
//   const int nrhs = X.Size();
//   const int N = height;
//   mfem::Array<mfem::Vector *> R(nrhs), D(nrhs);
//   std::vector<mfem::Vector> rrefs(nrhs), drefs(nrhs);
//   if (nrhs * N != r.Size())
//   {
//     r.SetSize(nrhs * N);
//     d.SetSize(nrhs * N);
//   }
//   for (int j = 0; j < nrhs; j++)
//   {
//     rrefs[j].MakeRef(r, j * N, N);
//     drefs[j].MakeRef(d, j * N, N);
//     R[j] = &rrefs[j];
//     D[j] = &drefs[j];
//   }

//   // Apply smoother: y = y + p(A) (x - A y) .
//   for (int it = 0; it < pc_it; it++)
//   {
//     if (iterative_mode || it > 0)
//     {
//       A->ArrayMult(Y, R);
//       for (int j = 0; j < nrhs; j++)
//       {
//         subtract(*X[j], *R[j], *R[j]);
//       }
//     }
//     else
//     {
//       for (int j = 0; j < nrhs; j++)
//       {
//         *R[j] = *X[j];
//         *Y[j] = 0.0;
//       }
//     }

//     // 4th-kind Chebyshev smoother
//     {
//       const auto *DI = dinv.Read();
//       for (int j = 0; j < nrhs; j++)
//       {
//         const auto *RR = R[j]->Read();
//         auto *DD = D[j]->ReadWrite();

//         mfem::forall(N, [=] MFEM_HOST_DEVICE(int i)
//                      { DD[i] = 4.0 / (3.0 * lambda_max) * DI[i] * RR[i]; });
//       }
//     }
//     for (int k = 1; k < order; k++)
//     {
//       A->ArrayAddMult(D, R, -1.0);
//       {
//         // From Phillips and Fischer or Lottes (with k -> k + 1 shift due to 1-based
//         // indexing)
//         const double sd = (2.0 * k - 1.0) / (2.0 * k + 3.0);
//         const double sr = (8.0 * k + 4.0) / ((2.0 * k + 3.0) * lambda_max);
//         const auto *DI = dinv.Read();
//         for (int j = 0; j < nrhs; j++)
//         {
//           const auto *RR = R[j]->Read();
//           auto *YY = Y[j]->ReadWrite();
//           auto *DD = D[j]->ReadWrite();
//           mfem::forall(N,
//                        [=] MFEM_HOST_DEVICE(int i)
//                        {
//                          YY[i] += DD[i];
//                          DD[i] = sd * DD[i] + sr * DI[i] * RR[i];
//                        });
//         }
//       }
//     }
//     for (int j = 0; j < nrhs; j++)
//     {
//       *Y[j] += *D[j];
//     }
//   }
// }

}  // namespace palace
