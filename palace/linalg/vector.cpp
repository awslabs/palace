// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "vector.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <utility>
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
    MFEM_VERIFY(use_dev, "ComplexVector::Set got a device pointer, but the vector is not "
                         "set to use the device (UseDevice(true))!");
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
    MFEM_VERIFY(use_dev, "ComplexVector::Get got a device pointer, but the vector is not "
                         "set to use the device (UseDevice(true))!");
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

namespace
{

// The sums of HostSums over [0, N). With several threads, each sums a contiguous block and
// the blocks are added in order, so the result only depends on the number of threads.
template <int M, typename Function>
std::array<double, M> ThreadedHostSums(int N, Function &&f)
{
  const int nt = utils::GetMaxThreads();
  if (nt == 1)
  {
    return HostSums<M>(0, N, f);
  }
  std::vector<std::array<double, M>> blocks(nt, std::array<double, M>{});
  PalacePragmaOmp(parallel)
  {
    const std::int64_t t = utils::GetThreadNum(), T = utils::GetNumActiveThreads();
    blocks[t] =
        HostSums<M>(static_cast<int>(N * t / T), static_cast<int>(N * (t + 1) / T), f);
  }
  std::array<double, M> sum = {};
  for (const auto &block : blocks)
  {
    for (int m = 0; m < M; m++)
    {
      sum[m] += block[m];
    }
  }
  return sum;
}

}  // namespace

double LocalDot(const Vector &x, const Vector &y)
{
  MFEM_ASSERT(x.Size() == y.Size(), "Size mismatch for vector inner product!");
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    static hypre::HypreVector X, Y;
    X.Update(x);
    Y.Update(y);
    return hypre_SeqVectorInnerProd(X, Y);
  }
  const auto *X = x.HostRead(), *Y = y.HostRead();
  return ThreadedHostSums<1>(x.Size(), [=](int i) { return std::array{Y[i] * X[i]}; })[0];
}

std::complex<double> LocalDot(const ComplexVector &x, const ComplexVector &y)
{
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
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

  // On the host, a single pass over the data instead of one per real inner product.
  MFEM_ASSERT(x.Size() == y.Size(), "Size mismatch for vector inner product!");
  const auto *XR = x.Real().HostRead(), *XI = x.Imag().HostRead();
  const auto *YR = y.Real().HostRead(), *YI = y.Imag().HostRead();
  if (&x == &y)
  {
    return ThreadedHostSums<1>(x.Size(), [=](int i)
                               { return std::array{XR[i] * XR[i] + XI[i] * XI[i]}; })[0];
  }
  const auto sum = ThreadedHostSums<2>(
      x.Size(), [=](int i)
      { return std::array{YR[i] * XR[i] + YI[i] * XI[i], YR[i] * XI[i] - YI[i] * XR[i]}; });
  return {sum[0], sum[1]};
}

// We implement LocalSum using Hypre instead of using MFEM's Sum because it is
// more efficient on GPUs. TODO: Verify this
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

template <>
void AXPYAndAYPX(double alpha, Vector &x, Vector &y, double beta, const Vector &z)
{
  MFEM_ASSERT(y.Size() == x.Size() && z.Size() == x.Size(),
              "Size mismatch for vector updates!");
  const bool use_dev = x.UseDevice() || y.UseDevice() || z.UseDevice();
  const int N = x.Size();
  const auto *Z = z.Read(use_dev);
  auto *X = x.ReadWrite(use_dev);
  auto *Y = y.ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        const auto t = X[i];
                        Y[i] += alpha * t;
                        X[i] = Z[i] + beta * t;
                      });
}

template <>
void AXPYAndAYPX(std::complex<double> alpha, ComplexVector &x, ComplexVector &y,
                 std::complex<double> beta, const ComplexVector &z)
{
  MFEM_ASSERT(y.Size() == x.Size() && z.Size() == x.Size(),
              "Size mismatch for vector updates!");
  const bool use_dev = x.UseDevice() || y.UseDevice() || z.UseDevice();
  const int N = x.Size();
  const double ar = alpha.real(), ai = alpha.imag(), br = beta.real(), bi = beta.imag();
  const auto *ZR = z.Real().Read(use_dev);
  const auto *ZI = z.Imag().Read(use_dev);
  auto *XR = x.Real().ReadWrite(use_dev);
  auto *XI = x.Imag().ReadWrite(use_dev);
  auto *YR = y.Real().ReadWrite(use_dev);
  auto *YI = y.Imag().ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        const auto tr = XR[i], ti = XI[i];
                        YR[i] += ar * tr - ai * ti;
                        YI[i] += ai * tr + ar * ti;
                        XR[i] = ZR[i] + br * tr - bi * ti;
                        XI[i] = ZI[i] + bi * tr + br * ti;
                      });
}

namespace
{

// y = Σⱼ c[j] X[j] (plus y if add), for the real and, if YI is not null, the imaginary parts
// of the coefficients. The zero parts are skipped, and the terms are added in order.
void LinearCombination(std::span<const Vector> X, const double *CR, const double *CI,
                       Vector &yr, Vector *yi, bool add)
{
  const int n = static_cast<int>(X.size()), N = yr.Size();
  for (const auto &Xj : X)
  {
    MFEM_VERIFY(Xj.Size() == N, "Size mismatch for linear combination!");
  }
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    // One pass over y with all terms, each thread reading one entry of every X[j].
    mfem::Array<const double *> P(n);
    mfem::Array<double> C(2 * n);
    for (int j = 0; j < n; j++)
    {
      P[j] = X[j].Read();
      C[j] = CR[j];
      C[n + j] = CI ? CI[j] : 0.0;
    }
    const auto *XP = P.Read();
    const auto *CC = C.Read();
    auto *YR = add ? yr.ReadWrite() : yr.Write();
    auto *YI = yi ? (add ? yi->ReadWrite() : yi->Write()) : nullptr;
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   double sr = add ? YR[i] : 0.0, si = (add && YI) ? YI[i] : 0.0;
                   for (int j = 0; j < n; j++)
                   {
                     const double x = XP[j][i];
                     if (CC[j] != 0.0)
                     {
                       sr += CC[j] * x;
                     }
                     if (CC[n + j] != 0.0)
                     {
                       si += CC[n + j] * x;
                     }
                   }
                   YR[i] = sr;
                   if (YI)
                   {
                     YI[i] = si;
                   }
                 });
    return;
  }

  // On the host, add all terms to one block of y at a time, which stays in cache, so each
  // X[j] is read once and y is written once. The blocks are smaller when needed to give
  // every thread one.
  const int nt = utils::GetMaxThreads();
  const int B = std::clamp((N + nt - 1) / nt, 64, 1024);
  std::vector<const double *> XP(n);
  for (int j = 0; j < n; j++)
  {
    XP[j] = X[j].HostRead();
  }
  auto *YR = add ? yr.HostReadWrite() : yr.HostWrite();
  auto *YI = yi ? (add ? yi->HostReadWrite() : yi->HostWrite()) : nullptr;
  PalacePragmaOmp(parallel for schedule(static) if (N > B))
  for (int b = 0; b < N; b += B)
  {
    const int e = std::min(b + B, N);
    if (!add)
    {
      std::fill(YR + b, YR + e, 0.0);
      if (YI)
      {
        std::fill(YI + b, YI + e, 0.0);
      }
    }
    for (int j = 0; j < n; j++)
    {
      const double *Xj = XP[j];
      for (auto [Y, cj] : {std::pair{YR, CR[j]}, std::pair{YI, CI ? CI[j] : 0.0}})
      {
        if (Y && cj != 0.0)
        {
          for (int i = b; i < e; i++)
          {
            Y[i] += cj * Xj[i];
          }
        }
      }
    }
  }
}

}  // namespace

void LinearCombination(std::span<const Vector> X, std::span<const std::complex<double>> c,
                       ComplexVector &y, bool add)
{
  MFEM_VERIFY(c.size() == X.size(), "Size mismatch for linear combination coefficients!");
  const auto *C = reinterpret_cast<const double *>(c.data());
  std::vector<double> CR(c.size()), CI(c.size());
  for (std::size_t j = 0; j < c.size(); j++)
  {
    CR[j] = C[2 * j];
    CI[j] = C[2 * j + 1];
  }
  LinearCombination(X, CR.data(), CI.data(), y.Real(), &y.Imag(), add);
}

void LinearCombination(std::span<const Vector> X, std::span<const double> c, Vector &y,
                       bool add)
{
  MFEM_VERIFY(c.size() == X.size(), "Size mismatch for linear combination coefficients!");
  LinearCombination(X, c.data(), nullptr, y, nullptr, add);
}

void LocalDots(std::span<const Vector> X, std::span<const Vector *const> Y, double *H)
{
  MFEM_VERIFY(static_cast<long long>(X.size()) * static_cast<long long>(Y.size()) <=
                  std::numeric_limits<int>::max(),
              "Too many local inner products!");
  const int m = static_cast<int>(X.size()), k = static_cast<int>(Y.size()), mk = m * k;
  if (mk == 0)
  {
    return;
  }
  const int N = X[0].Size();
  mfem::Array<const double *> P(m + k);
  for (int i = 0; i < m + k; i++)
  {
    const Vector &x = (i < m) ? X[i] : *Y[i - m];
    MFEM_VERIFY(x.Size() == N, "Size mismatch for local inner products!");
    P[i] = x.Read();
  }

  // Partial sums over nb slices of the entries of each pair, then added in slice order, so
  // the result does not depend on the number of threads. On the host a slice is a block of
  // contiguous entries, with all pairs of a block next to each other and the Y index varying
  // fastest, so that each X is read once from memory while a block of every Y stays in
  // cache. On the device a slice is every nb-th entry, so that neighboring threads read
  // neighboring entries, with enough slices to fill the device. The partial sums take at
  // most max_partial doubles.
  constexpr int B = 1024, max_partial = 1 << 22;
  const bool dev = mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
  const int nb = std::max(1, std::min({dev ? std::max(N, 1) : (N + B - 1) / B,
                                       dev ? 8192 : std::numeric_limits<int>::max(),
                                       max_partial / mk}));
  const int len = (N + nb - 1) / nb;  // Entries per slice
  Vector partial(nb * mk, mfem::Device::GetDeviceMemoryType()),
      h(mk, mfem::Device::GetDeviceMemoryType());
  partial.UseDevice(true);
  h.UseDevice(true);
  const auto *XY = P.Read();
  auto *S = partial.Write();
  mfem::forall(nb * mk,
               [=] MFEM_HOST_DEVICE(int t)
               {
                 const int b = dev ? t % nb : t / mk, q = dev ? t / nb : t % mk;
                 const int i = dev ? q % m : q / k, j = dev ? q / m : q % k;
                 const double *x = XY[i], *y = XY[m + j];
                 // Entries l = start + r stride of the slice, r < count.
                 const int start = dev ? b : b * len, stride = dev ? nb : 1;
                 const int rest = N - start;
                 const int count =
                     dev ? (rest + nb - 1) / nb : (rest < len ? (rest > 0 ? rest : 0) : len);
                 constexpr int K = 8;
                 double s[K] = {};
                 int r = 0;
                 for (; r + K <= count; r += K)
                 {
                   for (int l = 0; l < K; l++)
                   {
                     const int e = start + (r + l) * stride;
                     s[l] += x[e] * y[e];
                   }
                 }
                 for (int l = 0; r < count; r++, l++)
                 {
                   const int e = start + r * stride;
                   s[l] += x[e] * y[e];
                 }
                 double sum = 0.0;
                 for (int l = 0; l < K; l++)
                 {
                   sum += s[l];
                 }
                 S[b * mk + i + j * m] = sum;
               });
  const auto *PS = partial.Read();
  auto *D = h.Write();
  mfem::forall(mk,
               [=] MFEM_HOST_DEVICE(int p)
               {
                 double sum = 0.0;
                 for (int b = 0; b < nb; b++)
                 {
                   sum += PS[p + b * mk];
                 }
                 D[p] = sum;
               });
  std::copy_n(h.HostRead(), mk, H);
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
