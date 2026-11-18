// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "vector.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <random>
#include <mfem/general/forall.hpp>
#include <mfem/general/reducers.hpp>
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

// The reducer of K sums for mfem::reduce.
template <int K>
struct SumsReducer
{
  struct value_type
  {
    double v[K];
  };
  static MFEM_HOST_DEVICE void Join(value_type &a, const value_type &b)
  {
    for (int k = 0; k < K; k++)
    {
      a.v[k] += b.v[k];
    }
  }
  static MFEM_HOST_DEVICE void SetInitialValue(value_type &a) { a = {}; }
};

// Whether the sums over vectors with the given UseDevice() run on the device (Read(true)),
// including MFEM's debug device, so that the data are read where they were last written.
bool UseDeviceSums(bool use_dev)
{
  return use_dev && mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
}

// The K sums over i in [0, N) of the terms which f(i, s) adds to s[0], ..., s[K - 1]: one
// reduction kernel and one copy to the host on the device (if dev), ThreadedHostSums on the
// host.
template <int K, typename Function>
std::array<double, K> Sums(int N, bool dev, Function &&f)
{
  if (dev)
  {
    using Value = typename SumsReducer<K>::value_type;
    static mfem::Array<Value> workspace;
    Value sum = {};
    mfem::reduce(
        N, sum, [=] MFEM_HOST_DEVICE(int i, Value &r) { f(i, r.v); }, SumsReducer<K>{},
        true, workspace);
    return std::to_array(sum.v);
  }
  return ThreadedHostSums<K>(N,
                             [&](int i)
                             {
                               std::array<double, K> t = {};
                               f(i, t.data());
                               return t;
                             });
}

}  // namespace

double LocalDot(const Vector &x, const Vector &y)
{
  MFEM_ASSERT(x.Size() == y.Size(), "Size mismatch for vector inner product!");
  const bool dev = UseDeviceSums(x.UseDevice() || y.UseDevice());
  const auto *X = x.Read(dev), *Y = y.Read(dev);
  return Sums<1>(x.Size(), dev,
                 [=] MFEM_HOST_DEVICE(int i, double *s) { s[0] += Y[i] * X[i]; })[0];
}

std::complex<double> LocalDot(const ComplexVector &x, const ComplexVector &y)
{
  // A single pass over the data instead of one per real inner product.
  MFEM_ASSERT(x.Size() == y.Size(), "Size mismatch for vector inner product!");
  const bool dev = UseDeviceSums(x.UseDevice() || y.UseDevice());
  const auto *XR = x.Real().Read(dev), *XI = x.Imag().Read(dev);
  if (&x == &y)
  {
    return Sums<1>(x.Size(), dev, [=] MFEM_HOST_DEVICE(int i, double *s)
                   { s[0] += XR[i] * XR[i] + XI[i] * XI[i]; })[0];
  }
  const auto *YR = y.Real().Read(dev), *YI = y.Imag().Read(dev);
  const auto sum = Sums<2>(x.Size(), dev,
                           [=] MFEM_HOST_DEVICE(int i, double *s)
                           {
                             s[0] += YR[i] * XR[i] + YI[i] * XI[i];
                             s[1] += YR[i] * XI[i] - YI[i] * XR[i];
                           });
  return {sum[0], sum[1]};
}

double LocalSum(const Vector &x)
{
  const bool dev = UseDeviceSums(x.UseDevice());
  const auto *X = x.Read(dev);
  return Sums<1>(x.Size(), dev,
                 [=] MFEM_HOST_DEVICE(int i, double *s) { s[0] += X[i]; })[0];
}

std::complex<double> LocalSum(const ComplexVector &x)
{
  const bool dev = UseDeviceSums(x.UseDevice());
  const auto *XR = x.Real().Read(dev), *XI = x.Imag().Read(dev);
  const auto sum = Sums<2>(x.Size(), dev,
                           [=] MFEM_HOST_DEVICE(int i, double *s)
                           {
                             s[0] += XR[i];
                             s[1] += XI[i];
                           });
  return {sum[0], sum[1]};
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

template <int K>
void LinearCombination(std::span<const Vector> X, const double *C,
                       std::span<Vector *const> Y, const double *s, bool add)
{
  const int n = static_cast<int>(X.size()), N = Y[0]->Size();
  const bool dev = mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
  // An X[j] may alias an output. Without add, read the X[j] first: after Write() of an
  // output, Read() of the same memory would not copy its host data to the device. With add,
  // read the outputs first: on MFEM's debug device, reading an alias protects its host
  // memory, from which ReadWrite() of the output would then copy.
  std::vector<const double *> XP(n);
  double *YP[K], S[K];
  auto outputs = [&]()
  {
    for (int k = 0; k < K; k++)
    {
      MFEM_VERIFY(Y[k]->Size() == N, "Size mismatch for linear combination!");
      YP[k] = add ? Y[k]->ReadWrite(dev) : Y[k]->Write(dev);
      S[k] = s ? s[k] : 1.0;
    }
  };
  if (add)
  {
    outputs();
  }
  for (int j = 0; j < n; j++)
  {
    MFEM_VERIFY(X[j].Size() == N, "Size mismatch for linear combination!");
    XP[j] = X[j].Read(dev);
  }
  if (!add)
  {
    outputs();
  }
  // Terms with only zero coefficients are skipped, so that their vectors may hold anything:
  // under a floating-point model which assumes finite values, a zero coefficient alone does
  // not keep a NaN out of the sum.
  std::vector<int> terms;
  terms.reserve(n);
  for (int j = 0; j < n; j++)
  {
    if (std::any_of(C + K * j, C + K * (j + 1), [](double c) { return c != 0.0; }))
    {
      terms.push_back(j);
    }
  }
  const int nz = static_cast<int>(terms.size());
  if (dev)
  {
    // One pass over the outputs for each chunk of terms, each thread reading one entry of
    // every X[j] of the chunk. The pointers and coefficients of a chunk are kernel
    // arguments (under 4 KB), so that nothing is copied to the device and the host does not
    // wait for it. The terms which alias an output come first, so that they are read before
    // any output is written.
    constexpr int chunk = 3840 / (8 * (K + 1));
    std::vector<int> order(terms);
    const auto others =
        std::stable_partition(order.begin(), order.end(), [&](int j)
                              { return XP[j] && std::find(YP, YP + K, XP[j]) != YP + K; });
    MFEM_VERIFY(others - order.begin() <= chunk,
                "Too many terms of a linear combination alias an output!");
    for (int j0 = 0; j0 == 0 || j0 < nz; j0 += chunk)
    {
      const double *XC[chunk] = {};
      double CC[chunk][K] = {};
      const int nc = std::min(chunk, nz - j0);
      for (int j = 0; j < nc; j++)
      {
        XC[j] = XP[order[j0 + j]];
        std::copy_n(C + K * order[j0 + j], K, CC[j]);
      }
      const bool acc = add || j0 > 0, last = j0 + chunk >= nz;
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     double y[K];
                     for (int k = 0; k < K; k++)
                     {
                       y[k] = acc ? YP[k][i] : 0.0;
                     }
                     for (int j = 0; j < nc; j++)
                     {
                       const double x = XC[j][i];
                       for (int k = 0; k < K; k++)
                       {
                         if (CC[j][k] != 0.0)
                         {
                           y[k] += CC[j][k] * x;
                         }
                       }
                     }
                     for (int k = 0; k < K; k++)
                     {
                       YP[k][i] = last ? S[k] * y[k] : y[k];
                     }
                   });
    }
    return;
  }

  // On the host, add all terms to one block of the outputs at a time, kept in cache, so
  // that each X[j] is read once and each output is written once. The blocks are smaller
  // when needed to give every thread one.
  constexpr int Bmax = 1024;
  const int nt = utils::GetMaxThreads();
  const int B = std::clamp((N + nt - 1) / nt, 64, Bmax);
  PalacePragmaOmp(parallel for schedule(static) if (N > B))
  for (int b = 0; b < N; b += B)
  {
    const int len = std::min(B, N - b);
    double y[K][Bmax];
    for (int k = 0; k < K; k++)
    {
      for (int i = 0; i < len; i++)
      {
        y[k][i] = add ? YP[k][b + i] : 0.0;
      }
    }
    for (int j : terms)
    {
      const double *Xj = XP[j] + b;
      for (int k = 0; k < K; k++)
      {
        if (const double c = C[k + K * j]; c != 0.0)
        {
          for (int i = 0; i < len; i++)
          {
            y[k][i] += c * Xj[i];
          }
        }
      }
    }
    for (int k = 0; k < K; k++)
    {
      for (int i = 0; i < len; i++)
      {
        YP[k][b + i] = S[k] * y[k][i];
      }
    }
  }
}

}  // namespace

void LinearCombination(std::span<const Vector> X, const double *C,
                       std::span<Vector *const> Y, const double *s, bool add)
{
  switch (Y.size())
  {
    case 1:
      return LinearCombination<1>(X, C, Y, s, add);
    case 2:
      return LinearCombination<2>(X, C, Y, s, add);
    case 4:
      return LinearCombination<4>(X, C, Y, s, add);
    default:
      MFEM_ABORT("Unsupported number of outputs for linear combinations!");
  }
}

void LinearCombination(std::span<const Vector> X, std::span<const std::complex<double>> c,
                       ComplexVector &y, bool add)
{
  MFEM_VERIFY(c.size() == X.size(), "Size mismatch for linear combination coefficients!");
  LinearCombination<2>(X, reinterpret_cast<const double *>(c.data()),
                       std::array{&y.Real(), &y.Imag()}, nullptr, add);
}

void LinearCombination(std::span<const Vector> X, std::span<const double> c, Vector &y,
                       bool add)
{
  MFEM_VERIFY(c.size() == X.size(), "Size mismatch for linear combination coefficients!");
  LinearCombination<1>(X, c.data(), std::array{&y}, nullptr, add);
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
  const bool dev = mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
  std::vector<const double *> P(m + k);
  for (int i = 0; i < m + k; i++)
  {
    const Vector &x = (i < m) ? X[i] : *Y[i - m];
    MFEM_VERIFY(x.Size() == N, "Size mismatch for local inner products!");
    P[i] = x.Read(dev);
  }

  // Partial sums over nb slices of the entries of each pair, then added in a fixed order,
  // so the result does not depend on the number of threads. On the host a slice is a block
  // of contiguous entries, with all pairs of a block next to each other and the Y index
  // varying fastest, so that each X is read once from memory while a block of every Y stays
  // in cache. On the device a slice is every nb-th entry, so that neighboring threads read
  // neighboring entries, with enough slices to fill the device, and a thread computes the
  // products of a slice of one X with all Y of a launch, so that each X is read once from
  // memory. The vectors of a launch are kernel arguments (under 4 KB), so that nothing is
  // copied to the device. The partial sums take at most max_partial doubles.
  constexpr int B = 1024, max_partial = 1 << 22;
  const int nb = std::max(
      1, std::min({dev ? std::max(N, 1) : (N + B - 1) / B,
                   dev ? 8192 : std::numeric_limits<int>::max(), max_partial / mk}));
  Vector partial(nb * mk, mfem::Device::GetDeviceMemoryType());
  partial.UseDevice(true);
  auto *S = partial.Write();
  if (dev)
  {
    Vector h(mk, mfem::Device::GetDeviceMemoryType());
    h.UseDevice(true);
    constexpr int nx = 240, ny = 16;
    for (int i0 = 0; i0 < m; i0 += nx)
    {
      for (int j0 = 0; j0 < k; j0 += ny)
      {
        const int mc = std::min(nx, m - i0), kc = std::min(ny, k - j0);
        const double *XC[nx] = {}, *YC[ny] = {};
        std::copy_n(P.begin() + i0, mc, XC);
        std::copy_n(P.begin() + m + j0, kc, YC);
        mfem::forall(nb * mc,
                     [=] MFEM_HOST_DEVICE(int t)
                     {
                       const int b = t % nb, i = t / nb;
                       double s[ny] = {};
                       for (int e = b; e < N; e += nb)
                       {
                         const double x = XC[i][e];
                         for (int j = 0; j < ny; j++)
                         {
                           if (j < kc)
                           {
                             s[j] += x * YC[j][e];
                           }
                         }
                       }
                       for (int j = 0; j < ny; j++)
                       {
                         if (j < kc)
                         {
                           S[b * mk + (i0 + i) + (j0 + j) * m] = s[j];
                         }
                       }
                     });
      }
    }
    // One block of threads adds the partial sums of each pair: with one thread per pair,
    // this would take most of the time for few pairs.
    const auto *PS = partial.Read();
    auto *D = h.Write();
    constexpr int bs = 256;
    mfem::forall_2D(mk, bs, 1,
                    [=] MFEM_HOST_DEVICE(int p)
                    {
                      MFEM_SHARED double sm[bs];
                      MFEM_FOREACH_THREAD(t, x, bs)
                      {
                        double sum = 0.0;
                        for (int b = t; b < nb; b += bs)
                        {
                          sum += PS[p + b * mk];
                        }
                        sm[t] = sum;
                      }
                      for (int w = bs / 2; w > 0; w /= 2)
                      {
                        MFEM_SYNC_THREAD;
                        MFEM_FOREACH_THREAD(t, x, w)
                        {
                          sm[t] += sm[t + w];
                        }
                      }
                      MFEM_FOREACH_THREAD(t, x, 1)
                      {
                        D[p] = sm[0];
                      }
                    });
    std::copy_n(h.HostRead(), mk, H);
    return;
  }

  const int len = (N + nb - 1) / nb;  // Entries per slice
  PalacePragmaOmp(parallel for schedule(static) if (N > B))
  for (int t = 0; t < nb * mk; t++)
  {
    // The products of the entries of slice b.
    const int b = t / mk, q = t % mk, i = q / k, j = q % k;
    const double *x = P[i], *y = P[m + j];
    S[b * mk + i + j * m] = HostSums<1>(b * len, std::min(N, (b + 1) * len),
                                        [=](int e) { return std::array{x[e] * y[e]}; })[0];
  }
  std::fill_n(H, mk, 0.0);
  for (int b = 0; b < nb; b++)
  {
    for (int p = 0; p < mk; p++)
    {
      H[p] += S[p + b * mk];
    }
  }
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
