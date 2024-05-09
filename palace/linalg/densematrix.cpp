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

// Compute matrix functions for symmetric real-valued 2x2 or 3x3 matrices. Returns the
// matrix U * f(Λ) * U' for input U * Λ * U'.
// Reference: Deledalle et al., Closed-form expressions of the eigen decomposition of 2x2
//            and 3x3 Hermitian matrices, HAL hal-01501221 (2017).
mfem::DenseMatrix MatrixFunction(const mfem::DenseMatrix &M,
                                 const std::function<double(const double &)> &functor)
{
  MFEM_ASSERT(M.Height() == M.Width(),
              "MatrixFunction only available for square matrices!");
  const auto N = M.Height();
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
  if (N == 2)
  {
    MFEM_ABORT("2x2 MatrixFunction is not implemented yet!");
  }
  else if (N == 3)
  {
    // Need to specialize based on the number of zeros and their locations.
    const auto &a = M(0, 0), &b = M(1, 1), &c = M(2, 2);
    const auto &d = M(0, 1), &e = M(1, 2), &f = M(0, 2);
    const bool d_non_zero = std::abs(d) > tol;
    const bool e_non_zero = std::abs(e) > tol;
    const bool f_non_zero = std::abs(f) > tol;
    if (!d_non_zero && !e_non_zero && !f_non_zero)
    {
      // a 0 0
      // 0 b 0
      // 0 0 c
      for (int i = 0; i < 3; i++)
      {
        Mout(i, i) = functor(M(i, i));
      }
      return Mout;
    }
    if (d_non_zero && !e_non_zero && !f_non_zero)
    {
      // a d 0
      // d b 0
      // 0 0 c
      const double disc = std::sqrt(a * a - 2.0 * a * b + b * b + 4.0 * d * d);
      const double lambda1 = c;
      const double lambda2 = (a + b - disc) / 2.0;
      const double lambda3 = (a + b + disc) / 2.0;
      const mfem::Vector v1{{0.0, 0.0, 1.0}};
      const mfem::Vector v2{{-(-a + b + disc) / (2.0 * d), 1.0, 0.0}};
      const mfem::Vector v3{{-(-a + b - disc) / (2.0 * d), 1.0, 0.0}};
      AddMult_a_VVt(functor(lambda1), v1, Mout);
      AddMult_a_VVt(functor(lambda2), v2, Mout);
      AddMult_a_VVt(functor(lambda3), v3, Mout);
      return Mout;
    }
    if (!d_non_zero && e_non_zero && !f_non_zero)
    {
      // a 0 0
      // 0 b e
      // 0 e c
      const double disc = std::sqrt(b * b - 2.0 * b * c + c * c + 4.0 * e * e);
      const double lambda1 = a;
      const double lambda2 = 0.5 * (b + c - disc);
      const double lambda3 = 0.5 * (b + c + disc);
      const mfem::Vector v1{{1.0, 0.0, 0.0}};
      const mfem::Vector v2{{0.0, -(-b + c + disc) / (2.0 * e), 1.0}};
      const mfem::Vector v3{{0.0, -(-b + c - disc) / (2.0 * e), 1.0}};
      AddMult_a_VVt(functor(lambda1), v1, Mout);
      AddMult_a_VVt(functor(lambda2), v2, Mout);
      AddMult_a_VVt(functor(lambda3), v3, Mout);
      return Mout;
    }
    if (!d_non_zero && !e_non_zero && f_non_zero)
    {
      // a 0 f
      // 0 b 0
      // f 0 c
      const double disc = std::sqrt(a * a - 2.0 * a * c + c * c + 4.0 * f * f);
      const double lambda1 = b;
      const double lambda2 = 0.5 * (a + c - disc);
      const double lambda3 = 0.5 * (a + c + disc);
      const mfem::Vector v1{{0.0, 1.0, 0.0}};
      const mfem::Vector v2{{-(-a + c + disc) / (2.0 * f), 0.0, 1.0}};
      const mfem::Vector v3{{-(-a + c - disc) / (2.0 * f), 0.0, 1.0}};
      AddMult_a_VVt(functor(lambda1), v1, Mout);
      AddMult_a_VVt(functor(lambda2), v2, Mout);
      AddMult_a_VVt(functor(lambda3), v3, Mout);
      return Mout;
    }
    if ((!d_non_zero && e_non_zero && f_non_zero) ||
        (d_non_zero && !e_non_zero && f_non_zero) ||
        (d_non_zero && e_non_zero && !f_non_zero))
    {
      MFEM_ABORT("This nonzero pattern is not currently supported for MatrixFunction!");
    }
    // General case for all nonzero:
    // a d f
    // d b e
    // f e c
    const double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d, e2 = e * e, f2 = f * f;
    const double a2mbmc = 2.0 * a - b - c;
    const double b2mamc = 2.0 * b - a - c;
    const double c2mamb = 2.0 * c - a - b;
    const double x1 = a2 + b2 + c2 - a * b - b * c + 3.0 * (d2 + e2 + f2);
    const double x2 = -(a2mbmc * b2mamc * c2mamb) +
                      9.0 * (c2mamb * d2 + b2mamc * f2 + a2mbmc * e2) - 54.0 * d * e * f;
    const double phi = std::atan2(std::sqrt(4.0 * x1 * x1 * x1 - x2 * x2), x2);
    const double lambda1 = (a + b + c - 2.0 * std::sqrt(x1) * std::cos(phi / 3.0)) / 3.0;
    const double lambda2 =
        (a + b + c + 2.0 * std::sqrt(x1) * std::cos((phi - M_PI) / 3.0)) / 3.0;
    const double lambda3 =
        (a + b + c + 2.0 * std::sqrt(x1) * std::cos((phi + M_PI) / 3.0)) / 3.0;

    auto SafeDivide = [&](double x, double y)
    {
      if (std::abs(x) <= tol)
      {
        return 0.0;
      }
      if (std::abs(x) >= tol && std::abs(y) <= tol)
      {
        MFEM_ABORT("Logic error: Zero denominator with nonzero numerator!");
        return 0.0;
      }
      return x / y;
    };
    const double m1 = SafeDivide(d * (c - lambda1) - e * f, f * (b - lambda1) - d * e);
    const double m2 = SafeDivide(d * (c - lambda2) - e * f, f * (b - lambda2) - d * e);
    const double m3 = SafeDivide(d * (c - lambda3) - e * f, f * (b - lambda3) - d * e);
    const double l1mcmem1 = lambda1 - c - e * m1;
    const double l2mcmem2 = lambda2 - c - e * m2;
    const double l3mcmem3 = lambda3 - c - e * m3;
    const double n1 = 1.0 + m1 * m1 + SafeDivide(std::pow(l1mcmem1, 2), f2);
    const double n2 = 1.0 + m2 * m2 + SafeDivide(std::pow(l2mcmem2, 2), f2);
    const double n3 = 1.0 + m3 * m3 + SafeDivide(std::pow(l3mcmem3, 2), f2);
    const double tlambda1 = functor(lambda1) / n1;
    const double tlambda2 = functor(lambda2) / n2;
    const double tlambda3 = functor(lambda3) / n3;

    const double at = (tlambda1 * l1mcmem1 * l1mcmem1 + tlambda2 * l2mcmem2 * l2mcmem2 +
                       tlambda3 * l3mcmem3 * l3mcmem3) /
                      f2;
    const double bt = tlambda1 * m1 * m1 + tlambda2 * m2 * m2 + tlambda3 * m3 * m3;
    const double ct = tlambda1 + tlambda2 + tlambda3;
    const double dt =
        (tlambda1 * m1 * l1mcmem1 + tlambda2 * m2 * l2mcmem2 + tlambda3 * m3 * l3mcmem3) /
        f;
    const double et = tlambda1 * m1 + tlambda2 * m2 + tlambda3 * m3;
    const double ft = (tlambda1 * l1mcmem1 + tlambda2 * l2mcmem2 + tlambda3 * l3mcmem3) / f;
    Mout(0, 0) = at;
    Mout(0, 1) = dt;
    Mout(0, 2) = ft;
    Mout(1, 0) = dt;
    Mout(1, 1) = bt;
    Mout(1, 2) = et;
    Mout(2, 0) = ft;
    Mout(2, 1) = et;
    Mout(2, 2) = ct;
    return Mout;
  }
  else
  {
    MFEM_ABORT("MatrixFunction only supports 2x2 or 3x3 matrices!");
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
  for (int k = 0; k < T.SizeK(); k++)
  {
    S(k) = MatrixSqrt(T(k));
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
  for (int k = 0; k < T.SizeK(); k++)
  {
    S(k) = MatrixPow(T(k), p);
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
  for (int k = 0; k < C.SizeK(); k++)
  {
    Mult(A(k), B(k), C(k));
  }
  return C;
}

}  // namespace linalg

}  // namespace palace
