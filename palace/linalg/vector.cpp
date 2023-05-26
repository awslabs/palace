// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "vector.hpp"

#include <cstdint>
#include <random>
#include <general/forall.hpp>

namespace palace
{

ComplexVector::ComplexVector(int n) : x(2 * n)
{
  xr.MakeRef(x, 0, n);
  xi.MakeRef(x, n, n);
}

ComplexVector::ComplexVector(const ComplexVector &y) : x(2 * y.Size())
{
  xr.MakeRef(x, 0, y.Size());
  xi.MakeRef(x, y.Size(), y.Size());
  Set(y.Real(), y.Imag());
}

ComplexVector::ComplexVector(const Vector &yr, const Vector &yi) : x(2 * yr.Size())
{
  MFEM_VERIFY(yr.Size() == yi.Size(),
              "Mismatch in dimension of real and imaginary matrix parts in ComplexVector!");
  xr.MakeRef(x, 0, yr.Size());
  xi.MakeRef(x, yr.Size(), yr.Size());
  Set(yr, yi);
}

ComplexVector::ComplexVector(const std::complex<double> *py, int n) : x(2 * n)
{
  xr.MakeRef(x, 0, n);
  xi.MakeRef(x, n, n);
  Set(py, n);
}

void ComplexVector::SetSize(int n)
{
  x.SetSize(2 * n);
  xr.MakeRef(x, 0, n);
  xi.MakeRef(x, n, n);
}

ComplexVector &ComplexVector::operator=(const ComplexVector &y)
{
  Set(y.Real(), y.Imag());
  return *this;
}

void ComplexVector::Set(const Vector &yr, const Vector &yi)
{
  MFEM_VERIFY(yr.Size() == yi.Size() && yr.Size() == Size(),
              "Mismatch in dimension of real and imaginary matrix parts in ComplexVector!");
  Real() = yr;
  Imag() = yi;
}

void ComplexVector::Set(const std::complex<double> *py, int n)
{
  MFEM_VERIFY(n == Size(),
              "Mismatch in dimension for array of std::complex<double> in ComplexVector!");
  Vector y(reinterpret_cast<double *>(const_cast<std::complex<double> *>(py)), 2 * n);
  const int N = n;
  const auto *Y = y.Read();
  auto *XR = Real().Write();
  auto *XI = Imag().Write();
  mfem::forall(N,
               [=] MFEM_HOST_DEVICE(int i)
               {
                 XR[i] = Y[2 * i];
                 XI[i] = Y[2 * i + 1];
               });
}

void ComplexVector::Get(std::complex<double> *py, int n) const
{
  MFEM_VERIFY(n == Size(),
              "Mismatch in dimension for array of std::complex<double> in ComplexVector!");
  Vector y(reinterpret_cast<double *>(py), 2 * n);
  const int N = n;
  const auto *XR = Real().Read();
  const auto *XI = Imag().Read();
  auto *Y = y.Write();
  mfem::forall(N,
               [=] MFEM_HOST_DEVICE(int i)
               {
                 Y[2 * i] = XR[i];
                 Y[2 * i + 1] = XI[i];
               });
  y.HostReadWrite();
}

void ComplexVector::Conj()
{
  Imag() *= -1.0;
}

ComplexVector &ComplexVector::operator=(std::complex<double> s)
{
  Real() = s.real();
  Imag() = s.imag();
  return *this;
}

ComplexVector &ComplexVector::operator*=(std::complex<double> s)
{
  const double sr = s.real();
  const double si = s.imag();
  if (si == 0.0)
  {
    Real() *= sr;
    Imag() *= sr;
  }
  else
  {
    const int N = Size();
    auto *XR = Real().ReadWrite();
    auto *XI = Imag().ReadWrite();
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   const double t = si * XR[i] + sr * XI[i];
                   XR[i] = sr * XR[i] - si * XI[i];
                   XI[i] = t;
                 });
  }
  return *this;
}

void ComplexVector::Reciprocal(bool abs)
{
  const int N = Size();
  auto *XR = Real().ReadWrite();
  auto *XI = Imag().ReadWrite();
  if (abs)
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   const double t = 1.0 / std::sqrt(XR[i] * XR[i] + XI[i] * XI[i]);
                   XR[i] = t;
                   XI[i] = 0.0;
                 });
  }
  else
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   const std::complex<double> t = 1.0 / std::complex<double>(XR[i], XI[i]);
                   XR[i] = t.real();
                   XI[i] = t.imag();
                 });
  }
}

void ComplexVector::SetSubVector(const mfem::Array<int> &rows, std::complex<double> s)
{
  const int N = rows.Size();
  const double sr = s.real();
  const double si = s.imag();
  const auto *idx = rows.Read();
  auto *XR = Real().ReadWrite();
  auto *XI = Imag().ReadWrite();
  mfem::forall(N,
               [=] MFEM_HOST_DEVICE(int i)
               {
                 const int id = idx[i];
                 XR[id] = sr;
               });
  mfem::forall(N,
               [=] MFEM_HOST_DEVICE(int i)
               {
                 const int id = idx[i];
                 XI[id] = si;
               });
}

std::complex<double> ComplexVector::Dot(const ComplexVector &y) const
{
  return {(Real() * y.Real()) + (Imag() * y.Imag()),
          (Imag() * y.Real()) - (Real() * y.Imag())};
}

std::complex<double> ComplexVector::TransposeDot(const ComplexVector &y) const
{
  return {(Real() * y.Real()) - (Imag() * y.Imag()),
          (Imag() * y.Real()) + (Real() * y.Imag())};
}

void ComplexVector::AXPY(std::complex<double> alpha, const ComplexVector &x)
{
  const int N = Size();
  const double ar = alpha.real();
  const double ai = alpha.imag();
  const auto *XR = x.Real().Read();
  const auto *XI = x.Imag().Read();
  auto *YR = Real().ReadWrite();
  auto *YI = Imag().ReadWrite();
  if (ai == 0.0)
  {
    mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { YR[i] += ar * XR[i]; });
    mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { YI[i] += ar * XI[i]; });
  }
  else
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   YR[i] += ar * XR[i] - ai * XI[i];
                   YI[i] += ai * XR[i] + ar * XI[i];
                 });
  }
}

void ComplexVector::AXPBY(std::complex<double> alpha, const ComplexVector &x,
                          std::complex<double> beta)
{
  const int N = Size();
  const double ar = alpha.real();
  const double ai = alpha.imag();
  const auto *XR = x.Real().Read();
  const auto *XI = x.Imag().Read();
  auto *YR = Real().ReadWrite();
  auto *YI = Imag().ReadWrite();
  if (beta == 0.0)
  {
    if (ai == 0.0)
    {
      mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { YR[i] = ar * XR[i]; });
      mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { YI[i] = ar * XI[i]; });
    }
    else
    {
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     YR[i] = ar * XR[i] - ai * XI[i];
                     YI[i] = ai * XR[i] + ar * XI[i];
                   });
    }
  }
  else
  {
    const double br = beta.real();
    const double bi = beta.imag();
    if (ai == 0.0 && bi == 0.0)
    {
      mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { YR[i] = ar * XR[i] + br * YR[i]; });
      mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { YI[i] = ar * XI[i] + br * YI[i]; });
    }
    else
    {
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     const double t = bi * YR[i] + br * YI[i];
                     YR[i] = ar * XR[i] - ai * XI[i] + br * YR[i] - bi * YI[i];
                     YI[i] = ai * XR[i] + ar * XI[i] + t;
                   });
    }
  }
}

void ComplexVector::AXPBYPCZ(std::complex<double> alpha, const ComplexVector &x,
                             std::complex<double> beta, const ComplexVector &y,
                             std::complex<double> gamma)
{
  const int N = Size();
  const double ar = alpha.real();
  const double ai = alpha.imag();
  const double br = beta.real();
  const double bi = beta.imag();
  const auto *XR = x.Real().Read();
  const auto *XI = x.Imag().Read();
  const auto *YR = y.Real().Read();
  const auto *YI = y.Imag().Read();
  auto *ZR = Real().Write();
  auto *ZI = Imag().Write();
  if (gamma == 0.0)
  {
    if (ai == 0.0 && bi == 0.0)
    {
      mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { ZR[i] = ar * XR[i] + br * YR[i]; });
      mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { ZI[i] = ar * XI[i] + br * YI[i]; });
    }
    else
    {
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     ZR[i] = ar * XR[i] - ai * XI[i] + br * YR[i] - bi * YI[i];
                     ZI[i] = ai * XR[i] + ar * XI[i] + bi * YR[i] + br * YI[i];
                   });
    }
  }
  else
  {
    const double gr = gamma.real();
    const double gi = gamma.imag();
    if (ai == 0.0 && bi == 0.0 && gi == 0.0)
    {
      mfem::forall(N, [=] MFEM_HOST_DEVICE(int i)
                   { ZR[i] = ar * XR[i] + br * YR[i] + gr * ZR[i]; });
      mfem::forall(N, [=] MFEM_HOST_DEVICE(int i)
                   { ZI[i] = ar * XI[i] + br * YI[i] + gr * ZI[i]; });
    }
    else
    {
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     const double t = gi * ZR[i] + gr * ZI[i];
                     ZR[i] = ar * XR[i] - ai * XI[i] + br * YR[i] - bi * YI[i] +
                             gr * ZR[i] - gi * ZI[i];
                     ZI[i] = ai * XR[i] + ar * XI[i] + bi * YR[i] + br * YI[i] + t;
                   });
    }
  }
}

namespace linalg
{

template <>
void SetRandom(MPI_Comm comm, Vector &x, int seed)
{
  if (seed == 0)
  {
    std::vector<std::uint32_t> seeds(1);
    std::seed_seq seed_gen{Mpi::Rank(comm)};
    seed_gen.generate(seeds.begin(), seeds.end());
    seed = static_cast<int>(seeds[0]);
  }
  x.Randomize(seed);
}

template <>
void SetRandomSign(MPI_Comm comm, Vector &x, int seed)
{
  SetRandom(comm, x, seed);
  const int N = x.Size();
  auto *X = x.ReadWrite();
  mfem::forall(N, [=] MFEM_HOST_DEVICE(int i)
               { X[i] = (X[i] > 0.0) ? 1.0 : ((X[i] < 0.0) ? -1.0 : 0.0); });
}

template <>
void SetRandom(MPI_Comm comm, ComplexVector &x, int seed)
{
  if (seed == 0)
  {
    std::vector<std::uint32_t> seeds(2);
    std::seed_seq seed_gen{2 * Mpi::Rank(comm), 2 * Mpi::Rank(comm) + 1};
    seed_gen.generate(seeds.begin(), seeds.end());
    SetRandom(comm, x.Real(), static_cast<int>(seeds[0]));
    SetRandom(comm, x.Imag(), static_cast<int>(seeds[1]));
  }
  else
  {
    SetRandom(comm, x.Real(), seed);
    SetRandom(comm, x.Imag(), seed);
  }
}

template <>
void SetRandomSign(MPI_Comm comm, ComplexVector &x, int seed)
{
  SetRandom(comm, x, seed);
  const int N = x.Size();
  auto *XR = x.Real().ReadWrite();
  auto *XI = x.Imag().ReadWrite();
  mfem::forall(N, [=] MFEM_HOST_DEVICE(int i)
               { XR[i] = (XR[i] > 0.0) ? 1.0 : ((XR[i] < 0.0) ? -1.0 : 0.0); });
  mfem::forall(N, [=] MFEM_HOST_DEVICE(int i)
               { XI[i] = (XI[i] > 0.0) ? 1.0 : ((XI[i] < 0.0) ? -1.0 : 0.0); });
}

template <>
double Normalize(MPI_Comm comm, Vector &x, const Operator &B, Vector &Bx)
{
  B.Mult(x, Bx);
  double dot = Dot(comm, x, Bx);
  MFEM_ASSERT(dot > 0.0,
              "Non-positive vector norm in normalization (dot = " << dot << ")!");
  double norm = std::sqrt(dot);
  x *= 1.0 / norm;
  return norm;
}

template <>
double Normalize(MPI_Comm comm, ComplexVector &x, const Operator &B, ComplexVector &Bx)
{
  // For SPD B, xᴴ B x is real.
  B.Mult(x.Real(), Bx.Real());
  B.Mult(x.Imag(), Bx.Imag());
  std::complex<double> dot = Dot(comm, x, Bx);
  MFEM_ASSERT(dot.real() > 0.0 && std::abs(dot.imag()) < 1.0e-9 * dot.real(),
              "Non-positive vector norm in normalization (dot = " << dot << ")!");
  double norm = std::sqrt(dot.real());
  x *= 1.0 / norm;
  return norm;
}

template <>
void AXPY(double alpha, const Vector &x, Vector &y)
{
  if (alpha == 1.0)
  {
    y += x;
  }
  else
  {
    y.Add(alpha, x);
  }
}

template <>
void AXPY(double alpha, const ComplexVector &x, ComplexVector &y)
{
  y.AXPY(alpha, x);
}

template <>
void AXPY(std::complex<double> alpha, const ComplexVector &x, ComplexVector &y)
{
  y.AXPY(alpha, x);
}

template <>
void AXPBY(double alpha, const Vector &x, double beta, Vector &y)
{
  add(alpha, x, beta, y, y);
}

template <>
void AXPBY(std::complex<double> alpha, const ComplexVector &x, std::complex<double> beta,
           ComplexVector &y)
{
  y.AXPBY(alpha, x, beta);
}

template <>
void AXPBY(double alpha, const ComplexVector &x, double beta, ComplexVector &y)
{
  y.AXPBY(alpha, x, beta);
}

template <>
void AXPBYPCZ(double alpha, const Vector &x, double beta, const Vector &y, double gamma,
              Vector &z)
{
  if (gamma == 0.0)
  {
    add(alpha, x, beta, y, z);
  }
  else
  {
    AXPBY(alpha, x, gamma, z);
    z.Add(beta, y);
  }
}

template <>
void AXPBYPCZ(std::complex<double> alpha, const ComplexVector &x, std::complex<double> beta,
              const ComplexVector &y, std::complex<double> gamma, ComplexVector &z)
{
  z.AXPBYPCZ(alpha, x, beta, y, gamma);
}

template <>
void AXPBYPCZ(double alpha, const ComplexVector &x, double beta, const ComplexVector &y,
              double gamma, ComplexVector &z)
{
  z.AXPBYPCZ(alpha, x, beta, y, gamma);
}

}  // namespace linalg

}  // namespace palace
