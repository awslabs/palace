// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "vector.hpp"

#include <algorithm>
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
                "Mismatch between sum of block dimensions and parent vector dimension!");
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

namespace
{

// Chunk length of the fused Gram-Schmidt kernels: the chunks of x and y (four arrays of
// 64 KiB for complex vectors) stay in the L2 cache while the basis vectors stream through.
constexpr int GS_CHUNK = 8192;

// Call f(t, begin, end) for contiguous parts of [0, n), one per thread.
template <typename Function>
void ForThreadParts(int n, Function &&f)
{
  if (utils::GetMaxThreads() == 1)
  {
    f(0, 0, n);
    return;
  }
  PalacePragmaOmp(parallel)
  {
    const std::int64_t t = utils::GetThreadNum(), T = utils::GetNumActiveThreads();
    f(static_cast<int>(t), static_cast<int>(n * t / T), static_cast<int>(n * (t + 1) / T));
  }
}

// Add V[j]ᵀ x to hx[j] and V[j]ᵀ y to hy[j] over [begin, end), and the sums of x², xy, y²
// to e[0], e[1], e[2] if e is not null.
template <bool Y>
void RealMultiDotPart(int begin, int end, int m, const double *const *V, const double *X,
                      const double *Yv, double *hx, double *hy, double *e)
{
  constexpr int K = 4;
  for (int c0 = begin; c0 < end; c0 += GS_CHUNK)
  {
    const int len = std::min(GS_CHUNK, end - c0);
    const double *__restrict__ x = X + c0;
    const double *__restrict__ y = Y ? Yv + c0 : nullptr;
    for (int j = 0; j < m; j++)
    {
      const double *__restrict__ v = V[j] + c0;
      double s0[K] = {}, s1[K] = {};
      int i = 0;
      for (; i + K <= len; i += K)
      {
        for (int k = 0; k < K; k++)
        {
          s0[k] += v[i + k] * x[i + k];
          if constexpr (Y)
          {
            s1[k] += v[i + k] * y[i + k];
          }
        }
      }
      for (int k = 0; i < len; i++, k++)
      {
        s0[k] += v[i] * x[i];
        if constexpr (Y)
        {
          s1[k] += v[i] * y[i];
        }
      }
      hx[j] += (s0[0] + s0[1]) + (s0[2] + s0[3]);
      if constexpr (Y)
      {
        hy[j] += (s1[0] + s1[1]) + (s1[2] + s1[3]);
      }
    }
    if (e)
    {
      double t0 = 0.0, t1 = 0.0, t2 = 0.0;
      for (int i = 0; i < len; i++)
      {
        t0 += x[i] * x[i];
        if constexpr (Y)
        {
          t1 += x[i] * y[i];
          t2 += y[i] * y[i];
        }
      }
      e[0] += t0;
      e[1] += t1;
      e[2] += t2;
    }
  }
}

// Complex version: hx, hy hold (real, imaginary) pairs, and e = (xᴴ x, Re xᴴ y, Im xᴴ y,
// yᴴ y).
template <bool Y>
void ComplexMultiDotPart(int begin, int end, int m, const double *const *VR,
                         const double *const *VI, const double *XR, const double *XI,
                         const double *YR, const double *YI, double *hx, double *hy,
                         double *e)
{
  constexpr int K = 4;
  for (int c0 = begin; c0 < end; c0 += GS_CHUNK)
  {
    const int len = std::min(GS_CHUNK, end - c0);
    const double *__restrict__ xr = XR + c0;
    const double *__restrict__ xi = XI + c0;
    const double *__restrict__ yr = Y ? YR + c0 : nullptr;
    const double *__restrict__ yi = Y ? YI + c0 : nullptr;
    for (int j = 0; j < m; j++)
    {
      const double *__restrict__ vr = VR[j] + c0;
      const double *__restrict__ vi = VI[j] + c0;
      double s[8][K] = {};
      auto Add = [&](int i, int k)
      {
        const double a = vr[i], b = vi[i];
        s[0][k] += a * xr[i];
        s[1][k] += b * xi[i];
        s[2][k] += a * xi[i];
        s[3][k] += b * xr[i];
        if constexpr (Y)
        {
          s[4][k] += a * yr[i];
          s[5][k] += b * yi[i];
          s[6][k] += a * yi[i];
          s[7][k] += b * yr[i];
        }
      };
      int i = 0;
      for (; i + K <= len; i += K)
      {
        for (int k = 0; k < K; k++)
        {
          Add(i + k, k);
        }
      }
      for (int k = 0; i < len; i++, k++)
      {
        Add(i, k);
      }
      double r[8];
      for (int l = 0; l < 8; l++)
      {
        r[l] = (s[l][0] + s[l][1]) + (s[l][2] + s[l][3]);
      }
      hx[2 * j] += r[0] + r[1];
      hx[2 * j + 1] += r[2] - r[3];
      if constexpr (Y)
      {
        hy[2 * j] += r[4] + r[5];
        hy[2 * j + 1] += r[6] - r[7];
      }
    }
    if (e)
    {
      double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;
      for (int i = 0; i < len; i++)
      {
        t0 += xr[i] * xr[i] + xi[i] * xi[i];
        if constexpr (Y)
        {
          t1 += xr[i] * yr[i] + xi[i] * yi[i];
          t2 += xr[i] * yi[i] - xi[i] * yr[i];
          t3 += yr[i] * yr[i] + yi[i] * yi[i];
        }
      }
      e[0] += t0;
      e[1] += t1;
      e[2] += t2;
      e[3] += t3;
    }
  }
}

// x = sx (x - Σⱼ hx[j] V[j]) (if HX, else x = sx x) and y = sy (y - Σⱼ hy[j] V[j] - gy x)
// (if Y) over [begin, end); returns the sum of the squares of the updated y (or x).
template <bool HX, bool Y>
double RealMultiSubtractPart(int begin, int end, int m, const double *const *V, double *X,
                             const double *hx, double sx, double *Yv, const double *hy,
                             double gy, double sy)
{
  double nrm2 = 0.0;
  for (int c0 = begin; c0 < end; c0 += GS_CHUNK)
  {
    const int len = std::min(GS_CHUNK, end - c0);
    double *__restrict__ x = X + c0;
    double *__restrict__ y = Y ? Yv + c0 : nullptr;
    if constexpr (HX || Y)
    {
      for (int j = 0; j < m; j++)
      {
        const double *__restrict__ v = V[j] + c0;
        const double a = HX ? hx[j] : 0.0, b = Y ? hy[j] : 0.0;
        for (int i = 0; i < len; i++)
        {
          if constexpr (HX)
          {
            x[i] -= a * v[i];
          }
          if constexpr (Y)
          {
            y[i] -= b * v[i];
          }
        }
      }
    }
    double t = 0.0;
    for (int i = 0; i < len; i++)
    {
      const double q = sx * x[i];
      x[i] = q;
      if constexpr (Y)
      {
        const double z = sy * (y[i] - gy * q);
        y[i] = z;
        t += z * z;
      }
      else
      {
        t += q * q;
      }
    }
    nrm2 += t;
  }
  return nrm2;
}

template <bool HX, bool Y>
double ComplexMultiSubtractPart(int begin, int end, int m, const double *const *VR,
                                const double *const *VI, double *XR, double *XI,
                                const double *hx, double sx, double *YR, double *YI,
                                const double *hy, double gr, double gi, double sy)
{
  double nrm2 = 0.0;
  for (int c0 = begin; c0 < end; c0 += GS_CHUNK)
  {
    const int len = std::min(GS_CHUNK, end - c0);
    double *__restrict__ xr = XR + c0;
    double *__restrict__ xi = XI + c0;
    double *__restrict__ yr = Y ? YR + c0 : nullptr;
    double *__restrict__ yi = Y ? YI + c0 : nullptr;
    if constexpr (HX || Y)
    {
      for (int j = 0; j < m; j++)
      {
        const double *__restrict__ vr = VR[j] + c0;
        const double *__restrict__ vi = VI[j] + c0;
        const double ar = HX ? hx[2 * j] : 0.0, ai = HX ? hx[2 * j + 1] : 0.0;
        const double br = Y ? hy[2 * j] : 0.0, bi = Y ? hy[2 * j + 1] : 0.0;
        for (int i = 0; i < len; i++)
        {
          const double a = vr[i], b = vi[i];
          if constexpr (HX)
          {
            xr[i] -= ar * a - ai * b;
            xi[i] -= ar * b + ai * a;
          }
          if constexpr (Y)
          {
            yr[i] -= br * a - bi * b;
            yi[i] -= br * b + bi * a;
          }
        }
      }
    }
    double t = 0.0;
    for (int i = 0; i < len; i++)
    {
      const double q0 = sx * xr[i], q1 = sx * xi[i];
      xr[i] = q0;
      xi[i] = q1;
      if constexpr (Y)
      {
        const double z0 = sy * (yr[i] - (gr * q0 - gi * q1));
        const double z1 = sy * (yi[i] - (gr * q1 + gi * q0));
        yr[i] = z0;
        yi[i] = z1;
        t += z0 * z0 + z1 * z1;
      }
      else
      {
        t += q0 * q0 + q1 * q1;
      }
    }
    nrm2 += t;
  }
  return nrm2;
}

// Run part(t, begin, end, out) on each thread part of [0, n) with a private zeroed output
// array of size s, and add the outputs in thread order into out.
template <typename Function>
void ThreadedPartSums(int n, int s, double *out, Function &&part)
{
  std::fill(out, out + s, 0.0);
  const int nt = utils::GetMaxThreads();
  if (nt == 1)
  {
    part(0, n, out);
    return;
  }
  std::vector<double> parts(static_cast<std::size_t>(nt) * s, 0.0);
  ForThreadParts(n, [&](int t, int begin, int end)
                 { part(begin, end, parts.data() + static_cast<std::size_t>(t) * s); });
  for (int t = 0; t < nt; t++)
  {
    for (int l = 0; l < s; l++)
    {
      out[l] += parts[static_cast<std::size_t>(t) * s + l];
    }
  }
}

}  // namespace

void LocalMultiDot(const std::vector<Vector> &V, int m, const Vector &x, const Vector *y,
                   double *hx, double *hy, double *extra)
{
  MFEM_ASSERT(m >= 0 && m <= static_cast<int>(V.size()),
              "Out of bounds number of basis vectors for LocalMultiDot!");
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    for (int j = 0; j < m; j++)
    {
      hx[j] = LocalDot(x, V[j]);
      if (y)
      {
        hy[j] = LocalDot(*y, V[j]);
      }
    }
    if (extra)
    {
      extra[0] = LocalDot(x, x);
      if (y)
      {
        extra[1] = LocalDot(*y, x);
        extra[2] = LocalDot(*y, *y);
      }
    }
    return;
  }
  const int n = x.Size();
  std::vector<const double *> P(m);
  for (int j = 0; j < m; j++)
  {
    MFEM_ASSERT(V[j].Size() == n, "Size mismatch for LocalMultiDot!");
    P[j] = V[j].HostRead();
  }
  const double *X = x.HostRead(), *Y = y ? y->HostRead() : nullptr;
  std::vector<double> out(2 * m + 3);
  ThreadedPartSums(n, 2 * m + 3, out.data(),
                   [&](int begin, int end, double *o)
                   {
                     double *e = extra ? o + 2 * m : nullptr;
                     if (y)
                     {
                       RealMultiDotPart<true>(begin, end, m, P.data(), X, Y, o, o + m, e);
                     }
                     else
                     {
                       RealMultiDotPart<false>(begin, end, m, P.data(), X, nullptr, o,
                                               nullptr, e);
                     }
                   });
  std::copy(out.begin(), out.begin() + m, hx);
  if (y)
  {
    std::copy(out.begin() + m, out.begin() + 2 * m, hy);
  }
  if (extra)
  {
    extra[0] = out[2 * m];
    if (y)
    {
      extra[1] = out[2 * m + 1];
      extra[2] = out[2 * m + 2];
    }
  }
}

void LocalMultiDot(const std::vector<ComplexVector> &V, int m, const ComplexVector &x,
                   const ComplexVector *y, std::complex<double> *hx,
                   std::complex<double> *hy, std::complex<double> *extra)
{
  MFEM_ASSERT(m >= 0 && m <= static_cast<int>(V.size()),
              "Out of bounds number of basis vectors for LocalMultiDot!");
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    for (int j = 0; j < m; j++)
    {
      hx[j] = LocalDot(x, V[j]);
      if (y)
      {
        hy[j] = LocalDot(*y, V[j]);
      }
    }
    if (extra)
    {
      extra[0] = LocalDot(x, x);
      if (y)
      {
        extra[1] = LocalDot(*y, x);
        extra[2] = LocalDot(*y, *y);
      }
    }
    return;
  }
  const int n = x.Size();
  std::vector<const double *> PR(m), PI(m);
  for (int j = 0; j < m; j++)
  {
    MFEM_ASSERT(V[j].Size() == n, "Size mismatch for LocalMultiDot!");
    PR[j] = V[j].Real().HostRead();
    PI[j] = V[j].Imag().HostRead();
  }
  const double *XR = x.Real().HostRead(), *XI = x.Imag().HostRead();
  const double *YR = y ? y->Real().HostRead() : nullptr;
  const double *YI = y ? y->Imag().HostRead() : nullptr;
  std::vector<double> out(4 * m + 4);
  ThreadedPartSums(n, 4 * m + 4, out.data(),
                   [&](int begin, int end, double *o)
                   {
                     double *e = extra ? o + 4 * m : nullptr;
                     if (y)
                     {
                       ComplexMultiDotPart<true>(begin, end, m, PR.data(), PI.data(), XR,
                                                 XI, YR, YI, o, o + 2 * m, e);
                     }
                     else
                     {
                       ComplexMultiDotPart<false>(begin, end, m, PR.data(), PI.data(), XR,
                                                  XI, nullptr, nullptr, o, nullptr, e);
                     }
                   });
  for (int j = 0; j < m; j++)
  {
    hx[j] = {out[2 * j], out[2 * j + 1]};
    if (y)
    {
      hy[j] = {out[2 * m + 2 * j], out[2 * m + 2 * j + 1]};
    }
  }
  if (extra)
  {
    extra[0] = out[4 * m];
    if (y)
    {
      extra[1] = {out[4 * m + 1], out[4 * m + 2]};
      extra[2] = out[4 * m + 3];
    }
  }
}

double MultiSubtract(const std::vector<Vector> &V, int m, Vector &x, const double *hx,
                     double sx, Vector *y, const double *hy, double gy, double sy)
{
  MFEM_ASSERT(m >= 0 && m <= static_cast<int>(V.size()),
              "Out of bounds number of basis vectors for MultiSubtract!");
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    for (int j = 0; hx && j < m; j++)
    {
      x.Add(-hx[j], V[j]);
    }
    x *= sx;
    if (!y)
    {
      return LocalDot(x, x);
    }
    for (int j = 0; j < m; j++)
    {
      y->Add(-hy[j], V[j]);
    }
    y->Add(-gy, x);
    *y *= sy;
    return LocalDot(*y, *y);
  }
  const int n = x.Size();
  std::vector<const double *> P(m);
  for (int j = 0; j < m; j++)
  {
    MFEM_ASSERT(V[j].Size() == n, "Size mismatch for MultiSubtract!");
    P[j] = V[j].HostRead();
  }
  double *X = x.HostReadWrite(), *Y = y ? y->HostReadWrite() : nullptr;
  double nrm2;
  ThreadedPartSums(n, 1, &nrm2,
                   [&](int begin, int end, double *o)
                   {
                     if (hx && y)
                     {
                       *o = RealMultiSubtractPart<true, true>(begin, end, m, P.data(), X,
                                                              hx, sx, Y, hy, gy, sy);
                     }
                     else if (hx)
                     {
                       *o = RealMultiSubtractPart<true, false>(
                           begin, end, m, P.data(), X, hx, sx, nullptr, nullptr, 0.0, 1.0);
                     }
                     else if (y)
                     {
                       *o = RealMultiSubtractPart<false, true>(begin, end, m, P.data(), X,
                                                               nullptr, sx, Y, hy, gy, sy);
                     }
                     else
                     {
                       *o = RealMultiSubtractPart<false, false>(begin, end, m, P.data(), X,
                                                                nullptr, sx, nullptr,
                                                                nullptr, 0.0, 1.0);
                     }
                   });
  return nrm2;
}

double MultiSubtract(const std::vector<ComplexVector> &V, int m, ComplexVector &x,
                     const std::complex<double> *hx, double sx, ComplexVector *y,
                     const std::complex<double> *hy, std::complex<double> gy, double sy)
{
  MFEM_ASSERT(m >= 0 && m <= static_cast<int>(V.size()),
              "Out of bounds number of basis vectors for MultiSubtract!");
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    for (int j = 0; hx && j < m; j++)
    {
      x.Add(-hx[j], V[j]);
    }
    x *= sx;
    if (!y)
    {
      return LocalDot(x, x).real();
    }
    for (int j = 0; j < m; j++)
    {
      y->Add(-hy[j], V[j]);
    }
    y->Add(-gy, x);
    *y *= sy;
    return LocalDot(*y, *y).real();
  }
  const int n = x.Size();
  std::vector<const double *> PR(m), PI(m);
  for (int j = 0; j < m; j++)
  {
    MFEM_ASSERT(V[j].Size() == n, "Size mismatch for MultiSubtract!");
    PR[j] = V[j].Real().HostRead();
    PI[j] = V[j].Imag().HostRead();
  }
  double *XR = x.Real().HostReadWrite(), *XI = x.Imag().HostReadWrite();
  double *YR = y ? y->Real().HostReadWrite() : nullptr;
  double *YI = y ? y->Imag().HostReadWrite() : nullptr;
  const double *H = reinterpret_cast<const double *>(hx);
  const double *G = reinterpret_cast<const double *>(hy);
  const double gr = gy.real(), gi = gy.imag();
  double nrm2;
  ThreadedPartSums(
      n, 1, &nrm2,
      [&](int begin, int end, double *o)
      {
        if (hx && y)
        {
          *o = ComplexMultiSubtractPart<true, true>(begin, end, m, PR.data(), PI.data(), XR,
                                                    XI, H, sx, YR, YI, G, gr, gi, sy);
        }
        else if (hx)
        {
          *o = ComplexMultiSubtractPart<true, false>(begin, end, m, PR.data(), PI.data(),
                                                     XR, XI, H, sx, nullptr, nullptr,
                                                     nullptr, 0.0, 0.0, 1.0);
        }
        else if (y)
        {
          *o =
              ComplexMultiSubtractPart<false, true>(begin, end, m, PR.data(), PI.data(), XR,
                                                    XI, nullptr, sx, YR, YI, G, gr, gi, sy);
        }
        else
        {
          *o = ComplexMultiSubtractPart<false, false>(begin, end, m, PR.data(), PI.data(),
                                                      XR, XI, nullptr, sx, nullptr, nullptr,
                                                      nullptr, 0.0, 0.0, 1.0);
        }
      });
  return nrm2;
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
