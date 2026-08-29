// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "densematrix.hpp"

#include <functional>
#include <limits>
#include <mfem.hpp>
#include <mfem/linalg/kernels.hpp>

namespace palace
{

namespace
{

// Compute matrix functions for symmetric real-valued 1x1, 2x2, or 3x3 matrices. Returns
// the matrix U * f(Λ) * U' for input U * Λ * U', evaluated from the eigendecomposition
// computed by mfem::DenseMatrix::CalcEigenvalues.
mfem::DenseMatrix MatrixFunction(const mfem::DenseMatrix &M,
                                 const std::function<double(const double &)> &functor)
{
  MFEM_ASSERT(M.Height() == M.Width(),
              "MatrixFunction only available for square matrices!");
  const auto N = M.Height();
  MFEM_VERIFY(1 <= N && N <= 3,
              "MatrixFunction only supports 1x1, 2x2, or 3x3 matrices, N: " << N << "!");
  constexpr auto tol = 10.0 * std::numeric_limits<double>::epsilon();
  for (int i = 0; i < N; i++)
  {
    for (int j = i + 1; j < N; j++)
    {
      MFEM_VERIFY(std::abs(M(i, j) - M(j, i)) < tol,
                  "MatrixFunction only available for symmetric matrices ("
                      << M(i, j) << " != " << M(j, i) << ")!");
    }
  }
  mfem::DenseMatrix Mout(N, N);
  Mout = 0.0;
  if (N == 1)
  {
    Mout(0, 0) = functor(M(0, 0));
    return Mout;
  }
  // N is 2 or 3: evaluate f(M) = Σᵢ functor(λᵢ) vᵢ vᵢᵀ from the eigendecomposition returned
  // by MFEM's robust symmetric eigensolver. Unlike closed-form expressions, it handles
  // every nonzero pattern and returns an orthonormal eigenbasis even for degenerate spectra
  // (a rotated uniaxial material tensor is a common degenerate input). The eigenvectors are
  // stored consecutively, N doubles each; wrap each as a non-owning Vector (no allocation).
  double lambda[3], vec[9];
  M.CalcEigenvalues(lambda, vec);
  for (int i = 0; i < N; i++)
  {
    const mfem::Vector v(vec + N * i, N);
    AddMult_a_VVt(functor(lambda[i]) / (v * v), v, Mout);
  }
  return Mout;
}

}  // namespace

namespace linalg
{

mfem::DenseMatrix MatrixSqrt(const mfem::DenseMatrix &M)
{
  return MatrixFunction(M, [](auto s) { return std::sqrt(s); });
}

mfem::DenseTensor MatrixSqrt(const mfem::DenseTensor &T)
{
  mfem::DenseTensor S(T);
  mfem::DenseMatrix buffS, buffT;
  for (int k = 0; k < T.SizeK(); k++)
  {
    S(k, buffS).Set(1.0, MatrixSqrt(T(k, buffT)));
  }
  return S;
}

mfem::DenseMatrix MatrixPow(const mfem::DenseMatrix &M, double p)
{
  return MatrixFunction(M, [p](auto s) { return std::pow(s, p); });
}

mfem::DenseTensor MatrixPow(const mfem::DenseTensor &T, double p)
{
  mfem::DenseTensor S(T);
  mfem::DenseMatrix buffS, buffT;
  for (int k = 0; k < T.SizeK(); k++)
  {
    S(k, buffS).Set(1.0, MatrixPow(T(k, buffT), p));
  }
  return S;
}

double SingularValueMax(const mfem::DenseMatrix &M)
{
  MFEM_ASSERT(
      M.Height() == M.Width() && M.Height() > 0 && M.Height() <= 3,
      "Matrix singular values only available for square matrices of dimension <= 3!");
  const int N = M.Height();
  if (N == 1)
  {
    return M(0, 0);
  }
  else if (N == 2)
  {
    return mfem::kernels::CalcSingularvalue<2>(M.Data(), 0);
  }
  else
  {
    return mfem::kernels::CalcSingularvalue<3>(M.Data(), 0);
  }
}

double SingularValueMin(const mfem::DenseMatrix &M)
{
  MFEM_ASSERT(
      M.Height() == M.Width() && M.Height() > 0 && M.Height() <= 3,
      "Matrix singular values only available for square matrices of dimension <= 3!");
  const int N = M.Height();
  if (N == 1)
  {
    return M(0, 0);
  }
  else if (N == 2)
  {
    return mfem::kernels::CalcSingularvalue<2>(M.Data(), 1);
  }
  else
  {
    return mfem::kernels::CalcSingularvalue<3>(M.Data(), 2);
  }
}

mfem::DenseTensor Mult(const mfem::DenseTensor &A, const mfem::DenseTensor &B)
{
  MFEM_VERIFY(A.SizeK() == B.SizeK(),
              "Size mismatch for product of two DenseTensor objects!");
  mfem::DenseTensor C(A.SizeI(), B.SizeJ(), A.SizeK());
  mfem::DenseMatrix buffA, buffB, buffC;
  for (int k = 0; k < C.SizeK(); k++)
  {
    Mult(A(k, buffA), B(k, buffB), C(k, buffC));
  }
  return C;
}

}  // namespace linalg

}  // namespace palace
