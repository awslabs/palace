// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "operator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <mfem/general/forall.hpp>
#include "utils/communication.hpp"

namespace palace
{

const Operator *ComplexOperator::Real() const
{
  MFEM_ABORT("Real() is not implemented for base class ComplexOperator!");
  return nullptr;
}

const Operator *ComplexOperator::Imag() const
{
  MFEM_ABORT("Imag() is not implemented for base class ComplexOperator!");
  return nullptr;
}

void ComplexOperator::AssembleDiagonal(ComplexVector &diag) const
{
  MFEM_ABORT("Base class ComplexOperator does not implement AssembleDiagonal!");
}

void ComplexOperator::MultTranspose(const ComplexVector &x, ComplexVector &y) const
{
  MFEM_ABORT("Base class ComplexOperator does not implement MultTranspose!");
}

void ComplexOperator::MultHermitianTranspose(const ComplexVector &x, ComplexVector &y) const
{
  MFEM_ABORT("Base class ComplexOperator does not implement MultHermitianTranspose!");
}

void ComplexOperator::AddMult(const ComplexVector &x, ComplexVector &y,
                              const std::complex<double> a) const
{
  MFEM_ABORT("Base class ComplexOperator does not implement AddMult!");
}

void ComplexOperator::AddMultTranspose(const ComplexVector &x, ComplexVector &y,
                                       const std::complex<double> a) const
{
  MFEM_ABORT("Base class ComplexOperator does not implement AddMultTranspose!");
}

void ComplexOperator::AddMultHermitianTranspose(const ComplexVector &x, ComplexVector &y,
                                                const std::complex<double> a) const
{
  MFEM_ABORT("Base class ComplexOperator does not implement AddMultHermitianTranspose!");
}

ComplexWrapperOperator::ComplexWrapperOperator(std::unique_ptr<Operator> &&dAr,
                                               std::unique_ptr<Operator> &&dAi,
                                               const Operator *pAr, const Operator *pAi)
  : ComplexOperator(), data_Ar(std::move(dAr)), data_Ai(std::move(dAi)),
    Ar((data_Ar != nullptr) ? data_Ar.get() : pAr),
    Ai((data_Ai != nullptr) ? data_Ai.get() : pAi)
{
  MFEM_VERIFY(Ar || Ai, "Cannot construct ComplexWrapperOperator from an empty matrix!");
  MFEM_VERIFY((!Ar || !Ai) || (Ar->Height() == Ai->Height() && Ar->Width() == Ai->Width()),
              "Mismatch in dimension of real and imaginary matrix parts!");
  tx.UseDevice(true);
  ty.UseDevice(true);
  height = Ar ? Ar->Height() : Ai->Height();
  width = Ar ? Ar->Width() : Ai->Width();
}

ComplexWrapperOperator::ComplexWrapperOperator(std::unique_ptr<Operator> &&Ar,
                                               std::unique_ptr<Operator> &&Ai)
  : ComplexWrapperOperator(std::move(Ar), std::move(Ai), nullptr, nullptr)
{
}

ComplexWrapperOperator::ComplexWrapperOperator(const Operator *Ar, const Operator *Ai)
  : ComplexWrapperOperator(nullptr, nullptr, Ar, Ai)
{
}

void ComplexWrapperOperator::AssembleDiagonal(ComplexVector &diag) const
{
  diag = 0.0;
  if (Ar)
  {
    Ar->AssembleDiagonal(diag.Real());
  }
  if (Ai)
  {
    Ai->AssembleDiagonal(diag.Imag());
  }
}

void ComplexWrapperOperator::Mult(const ComplexVector &x, ComplexVector &y) const
{
  constexpr bool zero_real = false;
  constexpr bool zero_imag = false;
  const Vector &xr = x.Real();
  const Vector &xi = x.Imag();
  Vector &yr = y.Real();
  Vector &yi = y.Imag();
  if (Ai)
  {
    if (!zero_imag)
    {
      Ai->Mult(xi, yr);
      yr *= -1.0;
    }
    if (!zero_real)
    {
      Ai->Mult(xr, yi);
    }
  }
  else
  {
    yr = 0.0;
    yi = 0.0;
  }
  if (Ar)
  {
    if (!zero_real)
    {
      Ar->AddMult(xr, yr);
    }
    if (!zero_imag)
    {
      Ar->AddMult(xi, yi);
    }
  }
}

void ComplexWrapperOperator::MultTranspose(const ComplexVector &x, ComplexVector &y) const
{
  constexpr bool zero_real = false;
  constexpr bool zero_imag = false;
  const Vector &xr = x.Real();
  const Vector &xi = x.Imag();
  Vector &yr = y.Real();
  Vector &yi = y.Imag();
  if (Ai)
  {
    if (!zero_imag)
    {
      Ai->MultTranspose(xi, yr);
      yr *= -1.0;
    }
    if (!zero_real)
    {
      Ai->MultTranspose(xr, yi);
    }
  }
  else
  {
    yr = 0.0;
    yi = 0.0;
  }
  if (Ar)
  {
    if (!zero_real)
    {
      Ar->AddMultTranspose(xr, yr);
    }
    if (!zero_imag)
    {
      Ar->AddMultTranspose(xi, yi);
    }
  }
}

void ComplexWrapperOperator::MultHermitianTranspose(const ComplexVector &x,
                                                    ComplexVector &y) const
{
  constexpr bool zero_real = false;
  constexpr bool zero_imag = false;
  const Vector &xr = x.Real();
  const Vector &xi = x.Imag();
  Vector &yr = y.Real();
  Vector &yi = y.Imag();
  if (Ai)
  {
    if (!zero_imag)
    {
      Ai->MultTranspose(xi, yr);
    }
    if (!zero_real)
    {
      Ai->MultTranspose(xr, yi);
      yi *= -1.0;
    }
  }
  else
  {
    yr = 0.0;
    yi = 0.0;
  }
  if (Ar)
  {
    if (!zero_real)
    {
      Ar->AddMultTranspose(xr, yr);
    }
    if (!zero_imag)
    {
      Ar->AddMultTranspose(xi, yi);
    }
  }
}

void ComplexWrapperOperator::AddMult(const ComplexVector &x, ComplexVector &y,
                                     const std::complex<double> a) const
{
  constexpr bool zero_real = false;
  constexpr bool zero_imag = false;
  const Vector &xr = x.Real();
  const Vector &xi = x.Imag();
  Vector &yr = y.Real();
  Vector &yi = y.Imag();
  if (a.real() != 0.0 && a.imag() != 0.0)
  {
    ty.SetSize(height);
    Mult(x, ty);
    y.AXPY(a, ty);
  }
  else if (a.real() != 0.0)
  {
    if (Ar)
    {
      if (!zero_real)
      {
        Ar->AddMult(xr, yr, a.real());
      }
      if (!zero_imag)
      {
        Ar->AddMult(xi, yi, a.real());
      }
    }
    if (Ai)
    {
      if (!zero_imag)
      {
        Ai->AddMult(xi, yr, -a.real());
      }
      if (!zero_real)
      {
        Ai->AddMult(xr, yi, a.real());
      }
    }
  }
  else if (a.imag() != 0.0)
  {
    if (Ar)
    {
      if (!zero_real)
      {
        Ar->AddMult(xr, yi, a.imag());
      }
      if (!zero_imag)
      {
        Ar->AddMult(xi, yr, -a.imag());
      }
    }
    if (Ai)
    {
      if (!zero_imag)
      {
        Ai->AddMult(xi, yi, -a.imag());
      }
      if (!zero_real)
      {
        Ai->AddMult(xr, yr, -a.imag());
      }
    }
  }
}

void ComplexWrapperOperator::AddMultTranspose(const ComplexVector &x, ComplexVector &y,
                                              const std::complex<double> a) const
{
  constexpr bool zero_real = false;
  constexpr bool zero_imag = false;
  const Vector &xr = x.Real();
  const Vector &xi = x.Imag();
  Vector &yr = y.Real();
  Vector &yi = y.Imag();
  if (a.real() != 0.0 && a.imag() != 0.0)
  {
    tx.SetSize(width);
    MultTranspose(x, tx);
    y.AXPY(a, tx);
  }
  else if (a.real() != 0.0)
  {
    if (Ar)
    {
      if (!zero_real)
      {
        Ar->AddMultTranspose(xr, yr, a.real());
      }
      if (!zero_imag)
      {
        Ar->AddMultTranspose(xi, yi, a.real());
      }
    }
    if (Ai)
    {
      if (!zero_imag)
      {
        Ai->AddMultTranspose(xi, yr, -a.real());
      }
      if (!zero_real)
      {
        Ai->AddMultTranspose(xr, yi, a.real());
      }
    }
  }
  else if (a.imag() != 0.0)
  {
    if (Ar)
    {
      if (!zero_real)
      {
        Ar->AddMultTranspose(xr, yi, a.imag());
      }
      if (!zero_imag)
      {
        Ar->AddMultTranspose(xi, yr, -a.imag());
      }
    }
    if (Ai)
    {
      if (!zero_imag)
      {
        Ai->AddMultTranspose(xi, yi, -a.imag());
      }
      if (!zero_real)
      {
        Ai->AddMultTranspose(xr, yr, -a.imag());
      }
    }
  }
}

void ComplexWrapperOperator::AddMultHermitianTranspose(const ComplexVector &x,
                                                       ComplexVector &y,
                                                       const std::complex<double> a) const
{
  constexpr bool zero_real = false;
  constexpr bool zero_imag = false;
  const Vector &xr = x.Real();
  const Vector &xi = x.Imag();
  Vector &yr = y.Real();
  Vector &yi = y.Imag();
  if (a.real() != 0.0 && a.imag() != 0.0)
  {
    tx.SetSize(width);
    MultHermitianTranspose(x, tx);
    y.AXPY(a, tx);
  }
  else if (a.real() != 0.0)
  {
    if (Ar)
    {
      if (!zero_real)
      {
        Ar->AddMultTranspose(xr, yr, a.real());
      }
      if (!zero_imag)
      {
        Ar->AddMultTranspose(xi, yi, a.real());
      }
    }
    if (Ai)
    {
      if (!zero_imag)
      {
        Ai->AddMultTranspose(xi, yr, a.real());
      }
      if (!zero_real)
      {
        Ai->AddMultTranspose(xr, yi, -a.real());
      }
    }
  }
  else if (a.imag() != 0.0)
  {
    if (Ar)
    {
      if (!zero_real)
      {
        Ar->AddMultTranspose(xr, yi, a.imag());
      }
      if (!zero_imag)
      {
        Ar->AddMultTranspose(xi, yr, -a.imag());
      }
    }
    if (Ai)
    {
      if (!zero_imag)
      {
        Ai->AddMultTranspose(xi, yi, a.imag());
      }
      if (!zero_real)
      {
        Ai->AddMultTranspose(xr, yr, a.imag());
      }
    }
  }
}

template <>
void BaseDiagonalOperator<Operator>::Mult(const Vector &x, Vector &y) const
{
  const bool use_dev = x.UseDevice() || y.UseDevice();
  const int N = this->height;
  const auto *D = d.Read(use_dev);
  const auto *X = x.Read(use_dev);
  auto *Y = y.Write(use_dev);
  mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE(int i) { Y[i] = D[i] * X[i]; });
}

template <>
void BaseDiagonalOperator<ComplexOperator>::Mult(const ComplexVector &x,
                                                 ComplexVector &y) const
{
  const bool use_dev = x.UseDevice() || y.UseDevice();
  const int N = this->height;
  const auto *DR = d.Real().Read(use_dev);
  const auto *DI = d.Imag().Read(use_dev);
  const auto *XR = x.Real().Read(use_dev);
  const auto *XI = x.Imag().Read(use_dev);
  auto *YR = y.Real().Write(use_dev);
  auto *YI = y.Imag().Write(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        YR[i] = DR[i] * XR[i] - DI[i] * XI[i];
                        YI[i] = DI[i] * XR[i] + DR[i] * XI[i];
                      });
}

template <>
void BaseDiagonalOperator<Operator>::AddMult(const Vector &x, Vector &y,
                                             const double a) const
{
  const bool use_dev = x.UseDevice() || y.UseDevice();
  const int N = this->height;
  const auto *D = d.Read(use_dev);
  const auto *X = x.Read(use_dev);
  auto *Y = y.Write(use_dev);
  mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE(int i) { Y[i] += a * D[i] * X[i]; });
}

template <>
void BaseDiagonalOperator<ComplexOperator>::AddMult(const ComplexVector &x,
                                                    ComplexVector &y,
                                                    const std::complex<double> a) const
{
  const bool use_dev = x.UseDevice() || y.UseDevice();
  const int N = this->height;
  const double ar = a.real();
  const double ai = a.imag();
  const auto *DR = d.Real().Read(use_dev);
  const auto *DI = d.Imag().Read(use_dev);
  const auto *XR = x.Real().Read(use_dev);
  const auto *XI = x.Imag().Read(use_dev);
  auto *YR = y.Real().Write(use_dev);
  auto *YI = y.Imag().Write(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        const auto tr = DR[i] * XR[i] - DI[i] * XI[i];
                        const auto ti = DI[i] * XR[i] + DR[i] * XI[i];
                        YR[i] += ar * tr - ai * ti;
                        YI[i] += ai * tr + ar * ti;
                      });
}

template <>
void DiagonalOperatorHelper<BaseDiagonalOperator<ComplexOperator>,
                            ComplexOperator>::MultHermitianTranspose(const ComplexVector &x,
                                                                     ComplexVector &y) const
{
  const ComplexVector &d =
      static_cast<const BaseDiagonalOperator<ComplexOperator> *>(this)->d;
  const bool use_dev = x.UseDevice() || y.UseDevice();
  const int N = this->height;
  const auto *DR = d.Real().Read(use_dev);
  const auto *DI = d.Imag().Read(use_dev);
  const auto *XR = x.Real().Read(use_dev);
  const auto *XI = x.Imag().Read(use_dev);
  auto *YR = y.Real().Write(use_dev);
  auto *YI = y.Imag().Write(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        YR[i] = DR[i] * XR[i] + DI[i] * XI[i];
                        YI[i] = -DI[i] * XR[i] + DR[i] * XI[i];
                      });
}

template <>
void DiagonalOperatorHelper<BaseDiagonalOperator<ComplexOperator>, ComplexOperator>::
    AddMultHermitianTranspose(const ComplexVector &x, ComplexVector &y,
                              const std::complex<double> a) const
{
  const ComplexVector &d =
      static_cast<const BaseDiagonalOperator<ComplexOperator> *>(this)->d;
  const bool use_dev = x.UseDevice() || y.UseDevice();
  const int N = this->height;
  const double ar = a.real();
  const double ai = a.imag();
  const auto *DR = d.Real().Read(use_dev);
  const auto *DI = d.Imag().Read(use_dev);
  const auto *XR = x.Real().Read(use_dev);
  const auto *XI = x.Imag().Read(use_dev);
  auto *YR = y.Real().Write(use_dev);
  auto *YI = y.Imag().Write(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        const auto tr = DR[i] * XR[i] + DI[i] * XI[i];
                        const auto ti = -DI[i] * XR[i] + DR[i] * XI[i];
                        YR[i] += ar * tr - ai * ti;
                        YI[i] += ai * tr + ar * ti;
                      });
}

namespace linalg
{

template <>
double Norml2(MPI_Comm comm, const Vector &x, const Operator &B, Vector &Bx)
{
  B.Mult(x, Bx);
  double dot = Dot(comm, Bx, x);
  MFEM_ASSERT(dot > 0.0,
              "Non-positive vector norm in normalization (dot = " << dot << ")!");
  return std::sqrt(dot);
}

template <>
double Norml2(MPI_Comm comm, const ComplexVector &x, const Operator &B, ComplexVector &Bx)
{
  // For SPD B, xᴴ B x is real.
  B.Mult(x.Real(), Bx.Real());
  B.Mult(x.Imag(), Bx.Imag());
  std::complex<double> dot = Dot(comm, Bx, x);
  MFEM_ASSERT(dot.real() > 0.0 && std::abs(dot.imag()) < 1.0e-9 * dot.real(),
              "Non-positive vector norm in normalization (dot = " << dot << ")!");
  return std::sqrt(dot.real());
}

std::complex<double> Dot(MPI_Comm comm, const ComplexVector &x, const Operator &A,
                         const ComplexVector &y)
{
  ComplexVector Ax(A.Height());
  Ax.UseDevice(true);
  A.Mult(x.Real(), Ax.Real());
  A.Mult(x.Imag(), Ax.Imag());
  return Dot(comm, Ax, y);
}

std::complex<double> Dot(MPI_Comm comm, const ComplexVector &x, const ComplexOperator &A,
                         const ComplexVector &y)
{
  ComplexVector Ax(A.Height());
  Ax.UseDevice(true);
  A.Mult(x, Ax);
  return Dot(comm, Ax, y);
}

namespace
{

template <typename VecType, typename OperType>
void MultAdjoint(const OperType &A, const VecType &x, VecType &y)
{
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    A.MultHermitianTranspose(x, y);
  }
  else
  {
    A.MultTranspose(x, y);
  }
}

// Smallest and largest eigenvalues of the symmetric tridiagonal matrix with diagonal a and
// off-diagonal b, by bisection with Sturm sequence counts.
std::pair<double, double> TridiagonalExtremeEigenvalues(const std::vector<double> &a,
                                                        const std::vector<double> &b)
{
  const int n = static_cast<int>(a.size());
  double lo = a[0], hi = a[0], bmax = 0.0;
  for (int i = 0; i < n; i++)
  {
    const double r =
        ((i > 0) ? std::abs(b[i - 1]) : 0.0) + ((i < n - 1) ? std::abs(b[i]) : 0.0);
    lo = std::min(lo, a[i] - r);
    hi = std::max(hi, a[i] + r);
    bmax = std::max(bmax, (i < n - 1) ? std::abs(b[i]) : 0.0);
  }
  if (!std::isfinite(lo) || !std::isfinite(hi))
  {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    return {nan, nan};
  }
  const double pivmin = std::numeric_limits<double>::min() * std::max(1.0, bmax * bmax);
  auto Bisect = [&](int k)
  {
    // Smallest x with at least k eigenvalues less than x.
    double l = lo, h = hi;
    for (int it = 0; it < 2100; it++)
    {
      const double x = 0.5 * (l + h);
      if (x <= l || x >= h)
      {
        break;
      }
      int count = 0;
      double q = 1.0;
      for (int i = 0; i < n; i++)
      {
        q = a[i] - x - ((i > 0) ? b[i - 1] * b[i - 1] / q : 0.0);
        if (std::abs(q) < pivmin)
        {
          q = -pivmin;
        }
        count += (q < 0.0);
      }
      (count >= k ? h : l) = x;
    }
    return 0.5 * (l + h);
  };
  return {Bisect(1), Bisect(n)};
}

// Eigenvalues and eigenvectors (the columns of the column-major Q) of the symmetric
// tridiagonal matrix with diagonal a and off-diagonal b, by cyclic Jacobi rotations (the
// matrix is small).
void TridiagonalEigensystem(const std::vector<double> &a, const std::vector<double> &b,
                            std::vector<double> &lambda, std::vector<double> &Q)
{
  const int n = static_cast<int>(a.size());
  std::vector<double> T(n * n, 0.0);
  Q.assign(n * n, 0.0);
  for (int i = 0; i < n; i++)
  {
    T[i + n * i] = a[i];
    Q[i + n * i] = 1.0;
    if (i < n - 1)
    {
      T[i + 1 + n * i] = T[i + n * (i + 1)] = b[i];
    }
  }
  for (int sweep = 0; sweep < 50; sweep++)
  {
    double off = 0.0, diag = 0.0;
    for (int j = 0; j < n; j++)
    {
      diag += T[j + n * j] * T[j + n * j];
      for (int i = 0; i < j; i++)
      {
        off += T[i + n * j] * T[i + n * j];
      }
    }
    if (off <= std::numeric_limits<double>::epsilon() *
                   std::numeric_limits<double>::epsilon() * diag)
    {
      break;
    }
    for (int p = 0; p < n - 1; p++)
    {
      for (int q = p + 1; q < n; q++)
      {
        // Rotation J in the (p, q) plane such that (Jᵀ T J)_pq = 0.
        const double tpq = T[p + n * q];
        if (tpq == 0.0)
        {
          continue;
        }
        const double theta = (T[q + n * q] - T[p + n * p]) / (2.0 * tpq);
        const double t =
            std::copysign(1.0, theta) / (std::abs(theta) + std::hypot(theta, 1.0));
        const double c = 1.0 / std::hypot(t, 1.0), s = t * c;
        for (int k = 0; k < n; k++)
        {
          const double tkp = T[k + n * p], tkq = T[k + n * q];
          T[k + n * p] = c * tkp - s * tkq;
          T[k + n * q] = s * tkp + c * tkq;
        }
        for (int k = 0; k < n; k++)
        {
          const double tpk = T[p + n * k], tqk = T[q + n * k];
          T[p + n * k] = c * tpk - s * tqk;
          T[q + n * k] = s * tpk + c * tqk;
        }
        for (int k = 0; k < n; k++)
        {
          const double qkp = Q[k + n * p], qkq = Q[k + n * q];
          Q[k + n * p] = c * qkp - s * qkq;
          Q[k + n * q] = s * qkp + c * qkq;
        }
      }
    }
  }
  lambda.resize(n);
  for (int i = 0; i < n; i++)
  {
    lambda[i] = T[i + n * i];
  }
}

template <typename VecType, typename OperType>
double LanczosSpectralNorm(MPI_Comm comm, const OperType &A, bool herm, VecType *x0,
                           double tol, int max_it)
{
  // Lanczos iteration for the eigenvalue of largest magnitude of B = A, or B = Aᴴ A when A
  // is not Hermitian. The extreme Ritz values converge from inside the spectrum, and much
  // faster than the power iteration. Without reorthogonalization, the basis loses
  // orthogonality only along converged Ritz vectors, which does not change the extreme
  // Ritz values. The basis restarts from the current Ritz vector when it is full.
  MFEM_VERIFY(A.Height() == A.Width(), "Spectral norm requires a square matrix!");
  constexpr int max_basis = 16;
  const int n = A.Height();
  std::vector<VecType> V;
  V.reserve(max_basis);
  auto AddVector = [&]() -> VecType &
  {
    V.emplace_back(n);
    V.back().UseDevice(true);
    return V.back();
  };
  VecType w(n), t;
  w.UseDevice(true);
  if (!herm)
  {
    t.SetSize(n);
    t.UseDevice(true);
  }
  {
    VecType &v = AddVector();
    if (x0 && x0->Size() == n)
    {
      v = *x0;
    }
    else
    {
      SetRandom(comm, v);
    }
    Normalize(comm, v);
  }

  // Iterate until the Ritz value of largest magnitude has converged. Both ends of the Ritz
  // spectrum move outward, so for an indefinite operator the other end also has to converge
  // when it could still become the largest in magnitude. After a restart, the first Ritz
  // value is the Rayleigh quotient of the restart vector, so convergence is only checked
  // from the second iteration.
  std::vector<double> alpha, beta;
  double l = 0.0, l_min = 0.0, l_max = 0.0, l_min0 = 0.0, l_max0 = 0.0;
  bool test_min = false, test_max = true;
  int it = 0, k = 0;
  bool converged = false;
  auto RitzVector = [&](VecType &y)
  {
    // The Ritz vector of the eigenvalue of largest magnitude, plus the one at the other end
    // of the spectrum when that end is tested, so that a restart keeps both.
    if (k == 0)
    {
      y = V[0];
      return;
    }
    std::vector<double> lambda, Q;
    TridiagonalEigensystem(alpha, beta, lambda, Q);
    const int i_min = std::min_element(lambda.begin(), lambda.end()) - lambda.begin();
    const int i_max = std::max_element(lambda.begin(), lambda.end()) - lambda.begin();
    std::vector<double> s(k, 0.0);
    for (int j : {i_min, i_max})
    {
      if ((j == i_min && test_min) || (j == i_max && test_max && j != i_min))
      {
        for (int i = 0; i < k; i++)
        {
          s[i] += Q[i + k * j];
        }
      }
    }
    y = 0.0;
    for (int i = 0; i < k; i++)
    {
      AXPY(s[i], V[i], y);
    }
    Normalize(comm, y);
  };
  while (it < max_it)
  {
    const VecType &v = V[k];
    if (herm)
    {
      A.Mult(v, w);
    }
    else
    {
      A.Mult(v, t);
      MultAdjoint(A, t, w);
    }
    alpha.push_back(std::real(Dot(comm, w, v)));
    AXPY(-alpha.back(), v, w);
    if (k > 0)
    {
      AXPY(-beta.back(), V[k - 1], w);
    }
    const double b = Norml2(comm, w);
    it++;
    k++;

    std::tie(l_min, l_max) = TridiagonalExtremeEigenvalues(alpha, beta);
    l = std::max(std::abs(l_min), std::abs(l_max));
    if (!std::isfinite(l) || !std::isfinite(b))
    {
      l = std::numeric_limits<double>::quiet_NaN();
      break;
    }
    test_min = (-l_min >= l_max) || (l_max > 0.0 && l_min < -0.5 * l);
    test_max = (l_max >= -l_min) || (l_min < 0.0 && l_max > 0.5 * l);
    if ((k > 1 && (!test_min || std::abs(l_min - l_min0) <= tol * l) &&
         (!test_max || std::abs(l_max - l_max0) <= tol * l)) ||
        b <= 1.0e-12 * l)
    {
      converged = true;
      break;
    }
    l_min0 = l_min;
    l_max0 = l_max;
    if (k < max_basis)
    {
      beta.push_back(b);
      VecType &u = (k < static_cast<int>(V.size())) ? V[k] : AddVector();
      u = w;
      u *= 1.0 / b;
    }
    else
    {
      RitzVector(w);
      V[0] = w;
      alpha.clear();
      beta.clear();
      k = 0;
    }
  }
  if (!converged)
  {
    Mpi::Warning(
        comm, "Lanczos iteration did not converge in {:d} iterations, lambda = {:.3e}!\n",
        it, l);
  }
  if (x0)
  {
    x0->SetSize(n);
    x0->UseDevice(true);
    RitzVector(*x0);
  }
  return herm ? l : std::sqrt(l);
}

}  // namespace

double SpectralNorm(MPI_Comm comm, const Operator &A, bool sym, Vector *x0, double tol,
                    int max_it)
{
  return LanczosSpectralNorm(comm, A, sym, x0, tol, max_it);
}

double SpectralNorm(MPI_Comm comm, const ComplexOperator &A, bool herm, ComplexVector *x0,
                    double tol, int max_it)
{
  return LanczosSpectralNorm(comm, A, herm, x0, tol, max_it);
}

}  // namespace linalg

}  // namespace palace
