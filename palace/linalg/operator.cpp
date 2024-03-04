// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "operator.hpp"

#include <mfem/general/forall.hpp>
#include "linalg/slepc.hpp"
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

SumOperator::SumOperator(const Operator &op, double a) : Operator(op.Height(), op.Width())
{
  AddOperator(op, a);
  z.UseDevice(true);
}

void SumOperator::AddOperator(const Operator &op, double a)
{
  MFEM_VERIFY(op.Height() == height && op.Width() == width,
              "Invalid Operator dimensions for SumOperator!");
  ops.emplace_back(&op, a);
}

void SumOperator::Mult(const Vector &x, Vector &y) const
{
  if (ops.size() == 1)
  {
    ops.front().first->Mult(x, y);
    if (ops.front().second != 1.0)
    {
      y *= ops.front().second;
    }
    return;
  }
  y = 0.0;
  AddMult(x, y);
}

void SumOperator::MultTranspose(const Vector &x, Vector &y) const
{
  if (ops.size() == 1)
  {
    ops.front().first->MultTranspose(x, y);
    if (ops.front().second != 1.0)
    {
      y *= ops.front().second;
    }
    return;
  }
  y = 0.0;
  AddMultTranspose(x, y);
}

void SumOperator::AddMult(const Vector &x, Vector &y, const double a) const
{
  z.SetSize(y.Size());
  for (const auto &[op, c] : ops)
  {
    op->Mult(x, z);
    y.Add(a * c, z);
  }
}

void SumOperator::AddMultTranspose(const Vector &x, Vector &y, const double a) const
{
  z.SetSize(y.Size());
  for (const auto &[op, c] : ops)
  {
    op->MultTranspose(x, z);
    y.Add(a * c, z);
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
                        YI[i] += ai * ti + ar * ti;
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
                        YI[i] += ai * ti + ar * ti;
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

double SpectralNorm(MPI_Comm comm, const Operator &A, bool sym, double tol, int max_it)
{
  ComplexWrapperOperator Ar(const_cast<Operator *>(&A), nullptr);  // Non-owning constructor
  return SpectralNorm(comm, Ar, sym, tol, max_it);
}

double SpectralNorm(MPI_Comm comm, const ComplexOperator &A, bool herm, double tol,
                    int max_it)
{
  // XX TODO: Use ARPACK or SLEPc for this when configured.
#if defined(PALACE_WITH_SLEPC)
  return slepc::GetMaxSingularValue(comm, A, herm, tol, max_it);
#else
  // Power iteration loop: ||A||₂² = λₙ(Aᴴ A).
  int it = 0;
  double res = 0.0;
  double l = 0.0, l0 = 0.0;
  ComplexVector u(A.Height()), v(A.Height());
  u.UseDevice(true);
  v.UseDevice(true);
  SetRandom(comm, u);
  Normalize(comm, u);
  while (it < max_it)
  {
    A.Mult(u, v);
    if (herm)
    {
      u = v;
    }
    else
    {
      A.MultHermitianTranspose(v, u);
    }
    l = Normalize(comm, u);
    if (it > 0)
    {
      res = std::abs(l - l0) / l0;
      if (res < tol)
      {
        break;
      }
    }
    l0 = l;
    it++;
  }
  if (it >= max_it)
  {
    Mpi::Warning(comm,
                 "Power iteration did not converge in {:d} iterations, res = {:.3e}, "
                 "lambda = {:.3e}!\n",
                 it, res, l);
  }
  return herm ? l : std::sqrt(l);
#endif
}

}  // namespace linalg

}  // namespace palace
