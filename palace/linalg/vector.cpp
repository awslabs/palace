// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "vector.hpp"

#include <cstdint>
#include <random>
#include <mfem/general/forall.hpp>
#include "linalg/hypre.hpp"
#include "utils/omp.hpp"

namespace palace
{

ComplexVector::ComplexVector(int size) : xr(size), xi(size) {}

ComplexVector::ComplexVector(const ComplexVector &y) : ComplexVector(y.Size())
{
  UseDevice(y.UseDevice());
  Set(y);
}

ComplexVector::ComplexVector(const Vector &yr, const Vector &yi) : ComplexVector(yr.Size())
{
  MFEM_ASSERT(yr.Size() == yi.Size(),
              "Mismatch in dimension of real and imaginary parts in ComplexVector!");
  UseDevice(yr.UseDevice() || yi.UseDevice());
  Set(yr, yi);
}

ComplexVector::ComplexVector(const std::complex<double> *py, int size, bool on_dev)
  : ComplexVector(size)
{
  Set(py, size, on_dev);
}

ComplexVector::ComplexVector(Vector &y, int offset, int size)
{
  MakeRef(y, offset, size);
}

void ComplexVector::UseDevice(bool use_dev)
{
  xr.UseDevice(use_dev);
  xi.UseDevice(use_dev);
}

void ComplexVector::SetSize(int size)
{
  xr.SetSize(size);
  xi.SetSize(size);
}

void ComplexVector::MakeRef(Vector &y, int offset, int size)
{
  MFEM_ASSERT(y.Size() >= offset + 2 * size,
              "Insufficient storage for ComplexVector alias reference of the given size!");
  y.ReadWrite();  // Ensure memory is allocated on device before aliasing
  xr.MakeRef(y, offset, size);
  xi.MakeRef(y, offset + size, size);
}

void ComplexVector::Set(const ComplexVector &y)
{
  MFEM_ASSERT(y.Size() == Size(),
              "Mismatch in dimension of provided parts in ComplexVector!");
  Real() = y.Real();
  Imag() = y.Imag();
}

void ComplexVector::Set(const Vector &yr, const Vector &yi)
{
  MFEM_ASSERT(yr.Size() == yi.Size() && yr.Size() == Size(),
              "Mismatch in dimension of real and imaginary parts in ComplexVector!");
  Real() = yr;
  Imag() = yi;
}

void ComplexVector::Set(const std::complex<double> *py, int size, bool on_dev)
{
  MFEM_ASSERT(size == Size(),
              "Mismatch in dimension for array of std::complex<double> in ComplexVector!");
  auto SetImpl = [this](const double *Y, const int N, bool use_dev)
  {
    auto *XR = Real().Write(use_dev);
    auto *XI = Imag().Write(use_dev);
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          XR[i] = Y[2 * i];
                          XI[i] = Y[2 * i + 1];
                        });
  };
  const bool use_dev = UseDevice();
  if (((!use_dev || !mfem::Device::Allows(mfem::Backend::DEVICE_MASK)) && !on_dev) ||
      (use_dev && mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && on_dev))
  {
    // No copy (host pointer and not using device, or device pointer and using device).
    SetImpl(reinterpret_cast<const double *>(py), size, use_dev);
  }
  else if (!on_dev)
  {
    // Need copy from host to device (host pointer but using device).
    Vector y(2 * size);
    y.UseDevice(true);
    {
      auto *Y = y.HostWrite();
      PalacePragmaOmp(parallel for schedule(static))
      for (int i = 0; i < size; i++)
      {
        Y[2 * i] = py[i].real();
        Y[2 * i + 1] = py[i].imag();
      }
    }
    SetImpl(y.Read(use_dev), size, use_dev);
  }
  else
  {
    MFEM_ABORT("ComplexVector::Set using a device pointer is not implemented when MFEM is "
               "not configured to use the device!");
  }
}

void ComplexVector::Get(std::complex<double> *py, int size, bool on_dev) const
{
  MFEM_ASSERT(size == Size(),
              "Mismatch in dimension for array of std::complex<double> in ComplexVector!");
  auto GetImpl = [this](double *Y, const int N, bool use_dev)
  {
    const auto *XR = Real().Read(use_dev);
    const auto *XI = Imag().Read(use_dev);
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          Y[2 * i] = XR[i];
                          Y[2 * i + 1] = XI[i];
                        });
  };
  const bool use_dev = UseDevice();
  if (((!use_dev || !mfem::Device::Allows(mfem::Backend::DEVICE_MASK)) && !on_dev) ||
      (use_dev && mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && on_dev))
  {
    // No copy (host pointer and not using device, or device pointer and using device).
    GetImpl(reinterpret_cast<double *>(py), size, use_dev);
  }
  else if (!on_dev)
  {
    // Need copy from device to host (host pointer but using device).
    const auto *XR = Real().HostRead();
    const auto *XI = Imag().HostRead();
    PalacePragmaOmp(parallel for schedule(static))
    for (int i = 0; i < size; i++)
    {
      py[i].real(XR[i]);
      py[i].imag(XI[i]);
    }
  }
  else
  {
    MFEM_ABORT("ComplexVector::Get using a device pointer is not implemented when MFEM is "
               "not configured to use the device!");
  }
}

ComplexVector &ComplexVector::operator=(std::complex<double> s)
{
  Real() = s.real();
  Imag() = s.imag();
  return *this;
}

void ComplexVector::SetBlocks(const std::vector<const ComplexVector *> &y,
                              const std::vector<std::complex<double>> &s)
{
  MFEM_ASSERT(s.empty() || y.size() == s.size(),
              "Mismatch in dimension of vector blocks and scaling coefficients!");
  auto *XR = Real().Write();
  auto *XI = Imag().Write();
  int offset = 0;
  for (std::size_t b = 0; b < y.size(); b++)
  {
    MFEM_VERIFY(y[b] && ((b < y.size() - 1 && offset + y[b]->Size() < Size()) ||
                         (b == y.size() - 1 && offset + y[b]->Size() == Size())),
                "Mistmatch between sum of block dimensions and parent vector dimension!");
    const double sr = s.empty() ? 1.0 : s[b].real();
    const double si = s.empty() ? 0.0 : s[b].imag();
    const bool use_dev = UseDevice() || y[b]->UseDevice();
    const int N = y[b]->Size();
    const auto *YR = y[b]->Real().Read();
    const auto *YI = y[b]->Imag().Read();
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          XR[i] = sr * YR[i] - si * YI[i];
                          XI[i] = si * YR[i] + sr * YI[i];
                        });
    XR += N;
    XI += N;
    offset += N;
  }
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
    const bool use_dev = UseDevice();
    const int N = Size();
    auto *XR = Real().ReadWrite(use_dev);
    auto *XI = Imag().ReadWrite(use_dev);
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          const auto t = si * XR[i] + sr * XI[i];
                          XR[i] = sr * XR[i] - si * XI[i];
                          XI[i] = t;
                        });
  }
  return *this;
}

void ComplexVector::Conj()
{
  Imag() *= -1.0;
}

void ComplexVector::Abs()
{
  const bool use_dev = UseDevice();
  const int N = Size();
  auto *XR = Real().ReadWrite(use_dev);
  auto *XI = Imag().ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        XR[i] = std::sqrt(XR[i] * XR[i] + XI[i] * XI[i]);
                        XI[i] = 0.0;
                      });
}

void ComplexVector::Reciprocal()
{
  const bool use_dev = UseDevice();
  const int N = Size();
  auto *XR = Real().ReadWrite(use_dev);
  auto *XI = Imag().ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        const auto s = 1.0 / (XR[i] * XR[i] + XI[i] * XI[i]);
                        XR[i] *= s;
                        XI[i] *= -s;
                      });
}

std::complex<double> ComplexVector::Dot(const ComplexVector &y) const
{
  return {(Real() * y.Real()) + (Imag() * y.Imag()),
          (this == &y) ? 0.0 : ((Imag() * y.Real()) - (Real() * y.Imag()))};
}

std::complex<double> ComplexVector::TransposeDot(const ComplexVector &y) const
{
  return {(Real() * y.Real()) - (Imag() * y.Imag()),
          (this == &y) ? (2.0 * (Imag() * y.Real()))
                       : ((Imag() * y.Real()) + (Real() * y.Imag()))};
}

void ComplexVector::AXPY(std::complex<double> alpha, const ComplexVector &x)
{
  AXPY(alpha, x.Real(), x.Imag(), Real(), Imag());
}

void ComplexVector::AXPY(std::complex<double> alpha, const Vector &xr, const Vector &xi,
                         Vector &yr, Vector &yi)
{
  const bool use_dev = yr.UseDevice() || xr.UseDevice();
  const int N = yr.Size();
  const double ar = alpha.real();
  const double ai = alpha.imag();
  const auto *XR = xr.Read(use_dev);
  const auto *XI = xi.Read(use_dev);
  auto *YR = yr.ReadWrite(use_dev);
  auto *YI = yi.ReadWrite(use_dev);
  if (ai == 0.0)
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          YR[i] += ar * XR[i];
                          YI[i] += ar * XI[i];
                        });
  }
  else
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          const auto t = ai * XR[i] + ar * XI[i];
                          YR[i] += ar * XR[i] - ai * XI[i];
                          YI[i] += t;
                        });
  }
}

void ComplexVector::AXPBY(std::complex<double> alpha, const ComplexVector &x,
                          std::complex<double> beta)
{
  AXPBY(alpha, x.Real(), x.Imag(), beta, Real(), Imag());
}

void ComplexVector::AXPBY(std::complex<double> alpha, const Vector &xr, const Vector &xi,
                          std::complex<double> beta, Vector &yr, Vector &yi)
{
  const bool use_dev = yr.UseDevice() || xr.UseDevice();
  const int N = yr.Size();
  const double ar = alpha.real();
  const double ai = alpha.imag();
  const auto *XR = xr.Read(use_dev);
  const auto *XI = xi.Read(use_dev);
  if (beta == 0.0)
  {
    auto *YR = yr.Write(use_dev);
    auto *YI = yi.Write(use_dev);
    if (ai == 0.0)
    {
      mfem::forall_switch(use_dev, N,
                          [=] MFEM_HOST_DEVICE(int i)
                          {
                            YR[i] = ar * XR[i];
                            YI[i] = ar * XI[i];
                          });
    }
    else
    {
      mfem::forall_switch(use_dev, N,
                          [=] MFEM_HOST_DEVICE(int i)
                          {
                            const auto t = ai * XR[i] + ar * XI[i];
                            YR[i] = ar * XR[i] - ai * XI[i];
                            YI[i] = t;
                          });
    }
  }
  else
  {
    const double br = beta.real();
    const double bi = beta.imag();
    auto *YR = yr.ReadWrite(use_dev);
    auto *YI = yi.ReadWrite(use_dev);
    if (ai == 0.0 && bi == 0.0)
    {
      mfem::forall_switch(use_dev, N,
                          [=] MFEM_HOST_DEVICE(int i)
                          {
                            YR[i] = ar * XR[i] + br * YR[i];
                            YI[i] = ar * XI[i] + br * YI[i];
                          });
    }
    else
    {
      mfem::forall_switch(use_dev, N,
                          [=] MFEM_HOST_DEVICE(int i)
                          {
                            const auto t =
                                ai * XR[i] + ar * XI[i] + bi * YR[i] + br * YI[i];
                            YR[i] = ar * XR[i] - ai * XI[i] + br * YR[i] - bi * YI[i];
                            YI[i] = t;
                          });
    }
  }
}

void ComplexVector::AXPBYPCZ(std::complex<double> alpha, const ComplexVector &x,
                             std::complex<double> beta, const ComplexVector &y,
                             std::complex<double> gamma)
{
  AXPBYPCZ(alpha, x.Real(), x.Imag(), beta, y.Real(), y.Imag(), gamma, Real(), Imag());
}

void ComplexVector::AXPBYPCZ(std::complex<double> alpha, const Vector &xr, const Vector &xi,
                             std::complex<double> beta, const Vector &yr, const Vector &yi,
                             std::complex<double> gamma, Vector &zr, Vector &zi)
{
  const bool use_dev = zr.UseDevice() || xr.UseDevice() || yr.UseDevice();
  const int N = zr.Size();
  const double ar = alpha.real();
  const double ai = alpha.imag();
  const double br = beta.real();
  const double bi = beta.imag();
  const auto *XR = xr.Read(use_dev);
  const auto *XI = xi.Read(use_dev);
  const auto *YR = yr.Read(use_dev);
  const auto *YI = yi.Read(use_dev);
  if (gamma == 0.0)
  {
    auto *ZR = zr.Write(use_dev);
    auto *ZI = zi.Write(use_dev);
    if (ai == 0.0 && bi == 0.0)
    {
      mfem::forall_switch(use_dev, N,
                          [=] MFEM_HOST_DEVICE(int i)
                          {
                            ZR[i] = ar * XR[i] + br * YR[i];
                            ZI[i] = ar * XI[i] + br * YI[i];
                          });
    }
    else
    {
      mfem::forall_switch(use_dev, N,
                          [=] MFEM_HOST_DEVICE(int i)
                          {
                            const auto t =
                                ai * XR[i] + ar * XI[i] + bi * YR[i] + br * YI[i];
                            ZR[i] = ar * XR[i] - ai * XI[i] + br * YR[i] - bi * YI[i];
                            ZI[i] = t;
                          });
    }
  }
  else
  {
    const double gr = gamma.real();
    const double gi = gamma.imag();
    auto *ZR = zr.ReadWrite(use_dev);
    auto *ZI = zi.ReadWrite(use_dev);
    if (ai == 0.0 && bi == 0.0 && gi == 0.0)
    {
      mfem::forall_switch(use_dev, N,
                          [=] MFEM_HOST_DEVICE(int i)
                          {
                            ZR[i] = ar * XR[i] + br * YR[i] + gr * ZR[i];
                            ZI[i] = ar * XI[i] + br * YI[i] + gr * ZI[i];
                          });
    }
    else
    {
      mfem::forall_switch(use_dev, N,
                          [=] MFEM_HOST_DEVICE(int i)
                          {
                            const auto t = ai * XR[i] + ar * XI[i] + bi * YR[i] +
                                           br * YI[i] + gi * ZR[i] + gr * ZI[i];
                            ZR[i] = ar * XR[i] - ai * XI[i] + br * YR[i] - bi * YI[i] +
                                    gr * ZR[i] - gi * ZI[i];
                            ZI[i] = t;
                          });
    }
  }
}

namespace linalg
{

template <>
void SetSubVector(Vector &x, const mfem::Array<int> &rows, double s)
{
  const bool use_dev = x.UseDevice();
  const int N = rows.Size();
  const double sr = s;
  const auto *idx = rows.Read(use_dev);
  auto *X = x.ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        const auto id = idx[i];
                        X[id] = sr;
                      });
}

template <>
void SetSubVector(ComplexVector &x, const mfem::Array<int> &rows, double s)
{
  const bool use_dev = x.UseDevice();
  const int N = rows.Size();
  const double sr = s;
  const auto *idx = rows.Read(use_dev);
  auto *XR = x.Real().ReadWrite(use_dev);
  auto *XI = x.Imag().ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        const int id = idx[i];
                        XR[id] = sr;
                        XI[id] = 0.0;
                      });
}

template <>
void SetSubVector(Vector &x, const mfem::Array<int> &rows, const Vector &y)
{
  const bool use_dev = x.UseDevice();
  const int N = rows.Size();
  const auto *idx = rows.Read(use_dev);
  const auto *Y = y.Read(use_dev);
  auto *X = x.ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        const int id = idx[i];
                        X[id] = Y[id];
                      });
}

template <>
void SetSubVector(ComplexVector &x, const mfem::Array<int> &rows, const ComplexVector &y)
{
  const bool use_dev = x.UseDevice();
  const int N = rows.Size();
  const auto *idx = rows.Read(use_dev);
  const auto *YR = y.Real().Read(use_dev);
  const auto *YI = y.Imag().Read(use_dev);
  auto *XR = x.Real().ReadWrite(use_dev);
  auto *XI = x.Imag().ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        const int id = idx[i];
                        XR[id] = YR[id];
                        XI[id] = YI[id];
                      });
}

template <>
void SetSubVector(Vector &x, int start, const Vector &y)
{
  const bool use_dev = x.UseDevice();
  const int N = y.Size();
  MFEM_ASSERT(start >= 0 && start + N <= x.Size(), "Invalid range for SetSubVector!");
  const auto *Y = y.Read(use_dev);
  auto *X = x.ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        const int id = start + i;
                        X[id] = Y[i];
                      });
}

template <>
void SetSubVector(ComplexVector &x, int start, const ComplexVector &y)
{
  const bool use_dev = x.UseDevice();
  const int N = y.Size();
  MFEM_ASSERT(start >= 0 && start + N <= x.Size(), "Invalid range for SetSubVector!");
  const auto *YR = y.Real().Read(use_dev);
  const auto *YI = y.Imag().Read(use_dev);
  auto *XR = x.Real().ReadWrite(use_dev);
  auto *XI = x.Imag().ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        const int id = start + i;
                        XR[id] = YR[i];
                        XI[id] = YI[i];
                      });
}

template <>
void SetSubVector(Vector &x, int start, int end, double s)
{
  const bool use_dev = x.UseDevice();
  MFEM_ASSERT(start >= 0 && end <= x.Size() && start <= end,
              "Invalid range for SetSubVector!");
  const int N = end - start;
  const double sr = s;
  auto *X = x.ReadWrite(use_dev) + start;
  mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE(int i) { X[i] = sr; });
}

template <>
void SetSubVector(ComplexVector &x, int start, int end, double s)
{
  const bool use_dev = x.UseDevice();
  MFEM_ASSERT(start >= 0 && end <= x.Size() && start <= end,
              "Invalid range for SetSubVector!");
  const int N = end - start;
  const double sr = s;
  auto *XR = x.Real().ReadWrite(use_dev) + start;
  auto *XI = x.Imag().ReadWrite(use_dev) + start;
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        XR[i] = sr;
                        XI[i] = 0.0;
                      });
}

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
  x.Randomize(seed);  // On host always
}

template <>
void SetRandomReal(MPI_Comm comm, Vector &x, int seed)
{
  SetRandom(comm, x, seed);
}

template <>
void SetRandomSign(MPI_Comm comm, Vector &x, int seed)
{
  SetRandom(comm, x, seed);
  const bool use_dev = x.UseDevice();
  const int N = x.Size();
  auto *X = x.ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE(int i)
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
void SetRandomReal(MPI_Comm comm, ComplexVector &x, int seed)
{
  SetRandom(comm, x.Real(), seed);
  x.Imag() = 0.0;
}

template <>
void SetRandomSign(MPI_Comm comm, ComplexVector &x, int seed)
{
  SetRandom(comm, x, seed);
  const bool use_dev = x.UseDevice();
  const int N = x.Size();
  auto *XR = x.Real().ReadWrite(use_dev);
  auto *XI = x.Imag().ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        XR[i] = (XR[i] > 0.0) ? 1.0 : ((XR[i] < 0.0) ? -1.0 : 0.0);
                        XI[i] = (XI[i] > 0.0) ? 1.0 : ((XI[i] < 0.0) ? -1.0 : 0.0);
                      });
}

double LocalDot(const Vector &x, const Vector &y)
{
  static hypre::HypreVector X, Y;
  MFEM_ASSERT(x.Size() == y.Size(), "Size mismatch for vector inner product!");
  X.Update(x);
  Y.Update(y);
  return hypre_SeqVectorInnerProd(X, Y);
}

std::complex<double> LocalDot(const ComplexVector &x, const ComplexVector &y)
{
  if (&x == &y)
  {
    return {LocalDot(x.Real(), y.Real()) + LocalDot(x.Imag(), y.Imag()), 0.0};
  }
  else
  {
    return {LocalDot(x.Real(), y.Real()) + LocalDot(x.Imag(), y.Imag()),
            LocalDot(x.Imag(), y.Real()) - LocalDot(x.Real(), y.Imag())};
  }
}

double LocalSum(const Vector &x)
{
  static hypre::HypreVector X;
  X.Update(x);
  return hypre_SeqVectorSumElts(X);
}

std::complex<double> LocalSum(const ComplexVector &x)
{
  return {LocalSum(x.Real()), LocalSum(x.Imag())};
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

void Sqrt(Vector &x, double s)
{
  const bool use_dev = x.UseDevice();
  const int N = x.Size();
  auto *X = x.ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i) { X[i] = std::sqrt(X[i] * s); });
}

}  // namespace linalg

}  // namespace palace
