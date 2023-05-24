// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "operator.hpp"

#include <general/forall.hpp>
#include "linalg/slepc.hpp"
#include "utils/communication.hpp"

namespace palace
{

bool ComplexOperator::IsReal() const
{
  MFEM_ABORT("IsReal() is not implemented for base class ComplexOperator!");
  return false;
}

bool ComplexOperator::IsImag() const
{
  MFEM_ABORT("IsImag() is not implemented for base class ComplexOperator!");
  return false;
}

bool ComplexOperator::HasReal() const
{
  MFEM_ABORT("HasReal() is not implemented for base class ComplexOperator!");
  return false;
}

bool ComplexOperator::HasImag() const
{
  MFEM_ABORT("HasImag() is not implemented for base class ComplexOperator!");
  return false;
}

const Operator *ComplexOperator::Real() const
{
  MFEM_ABORT("Real() is not implemented for base class ComplexOperator!");
  return nullptr;
}

Operator *ComplexOperator::Real()
{
  MFEM_ABORT("Real() is not implemented for base class ComplexOperator!");
  return nullptr;
}

const Operator *ComplexOperator::Imag() const
{
  MFEM_ABORT("Imag() is not implemented for base class ComplexOperator!");
  return nullptr;
}

Operator *ComplexOperator::Imag()
{
  MFEM_ABORT("Imag() is not implemented for base class ComplexOperator!");
  return nullptr;
}

void ComplexOperator::MultTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                    Vector &yi, bool zero_real, bool zero_imag) const
{
  MFEM_ABORT("Base class ComplexOperator does not implement MultTranspose!");
}

void ComplexOperator::MultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                             Vector &yi, bool zero_real,
                                             bool zero_imag) const
{
  MFEM_ABORT("Base class ComplexOperator does not implement MultHermitianTranspose!");
}

void ComplexOperator::AddMult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                              const std::complex<double> a, bool zero_real,
                              bool zero_imag) const
{
  MFEM_ABORT("Base class ComplexOperator does not implement AddMult!");
}

void ComplexOperator::AddMultTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                       Vector &yi, const std::complex<double> a,
                                       bool zero_real, bool zero_imag) const
{
  MFEM_ABORT("Base class ComplexOperator does not implement AddMultTranspose!");
}

void ComplexOperator::AddMultHermitianTranspose(const Vector &xr, const Vector &xi,
                                                Vector &yr, Vector &yi,
                                                const std::complex<double> a,
                                                bool zero_real, bool zero_imag) const
{
  MFEM_ABORT("Base class ComplexOperator does not implement AddMultHermitianTranspose!");
}

ComplexWrapperOperator::ComplexWrapperOperator(std::unique_ptr<Operator> &&data_Ar,
                                               std::unique_ptr<Operator> &&data_Ai,
                                               Operator *Ar, Operator *Ai)
  : ComplexOperator(Ar ? Ar->Height() : (Ai ? Ai->Height() : 0),
                    Ar ? Ar->Width() : (Ai ? Ai->Width() : 0)),
    data_Ar(std::move(data_Ar)), data_Ai(std::move(data_Ai)),
    Ar(this->data_Ar ? this->data_Ar.get() : Ar),
    Ai(this->data_Ai ? this->data_Ai.get() : Ai)
{
  MFEM_VERIFY(Ar || Ai, "Cannot construct ComplexWrapperOperator from an empty matrix!");
  MFEM_VERIFY((!Ar || !Ai) || (Ar->Height() == Ai->Height() && Ar->Width() == Ai->Width()),
              "Mismatch in dimension of real and imaginary matrix parts!");
}

ComplexWrapperOperator::ComplexWrapperOperator(std::unique_ptr<Operator> &&Ar,
                                               std::unique_ptr<Operator> &&Ai)
  : ComplexWrapperOperator(std::move(Ar), std::move(Ai), nullptr, nullptr)
{
}

ComplexWrapperOperator::ComplexWrapperOperator(Operator *Ar, Operator *Ai)
  : ComplexWrapperOperator(nullptr, nullptr, Ar, Ai)
{
}

void ComplexWrapperOperator::Mult(const Vector &xr, const Vector &xi, Vector &yr,
                                  Vector &yi, bool zero_real, bool zero_imag) const
{
  if (Ar)
  {
    if (!zero_real)
    {
      Ar->Mult(xr, yr);
    }
    if (!zero_imag)
    {
      Ar->Mult(xi, yi);
    }
  }
  else
  {
    yr = 0.0;
    yi = 0.0;
  }
  if (Ai)
  {
    if (!zero_imag)
    {
      Ai->AddMult(xi, yr, -1.0);
    }
    if (!zero_real)
    {
      Ai->AddMult(xr, yi, 1.0);
    }
  }
}

void ComplexWrapperOperator::MultTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                           Vector &yi, bool zero_real, bool zero_imag) const
{
  if (Ar)
  {
    if (!zero_real)
    {
      Ar->MultTranspose(xr, yr);
    }
    if (!zero_imag)
    {
      Ar->MultTranspose(xi, yi);
    }
  }
  else
  {
    yr = 0.0;
    yi = 0.0;
  }
  if (Ai)
  {
    if (!zero_imag)
    {
      Ai->AddMultTranspose(xi, yr, -1.0);
    }
    if (!zero_real)
    {
      Ai->AddMultTranspose(xr, yi, 1.0);
    }
  }
}

void ComplexWrapperOperator::MultHermitianTranspose(const Vector &xr, const Vector &xi,
                                                    Vector &yr, Vector &yi, bool zero_real,
                                                    bool zero_imag) const
{
  if (Ar)
  {
    if (!zero_real)
    {
      Ar->MultTranspose(xr, yr);
    }
    if (!zero_imag)
    {
      Ar->MultTranspose(xi, yi);
    }
  }
  else
  {
    yr = 0.0;
    yi = 0.0;
  }
  if (Ai)
  {
    if (!zero_imag)
    {
      Ai->AddMultTranspose(xi, yr, 1.0);
    }
    if (!zero_real)
    {
      Ai->AddMultTranspose(xr, yi, -1.0);
    }
  }
}

void ComplexWrapperOperator::AddMult(const Vector &xr, const Vector &xi, Vector &yr,
                                     Vector &yi, const std::complex<double> a,
                                     bool zero_real, bool zero_imag) const
{
  if (a.real() != 0.0 && a.imag() != 0.0)
  {
    ty.SetSize(height);
    Mult(xr, xi, ty.Real(), ty.Imag(), zero_real, zero_imag);
    const int N = height;
    const double ar = a.real();
    const double ai = a.imag();
    const auto *TYR = ty.Real().Read();
    const auto *TYI = ty.Imag().Read();
    auto *YR = yr.ReadWrite();
    auto *YI = yi.ReadWrite();
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   YR[i] += ar * TYR[i] - ai * TYI[i];
                   YI[i] += ai * TYR[i] + ar * TYI[i];
                 });
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

void ComplexWrapperOperator::AddMultTranspose(const Vector &xr, const Vector &xi,
                                              Vector &yr, Vector &yi,
                                              const std::complex<double> a, bool zero_real,
                                              bool zero_imag) const
{
  if (a.real() != 0.0 && a.imag() != 0.0)
  {
    tx.SetSize(width);
    MultTranspose(xr, xi, tx.Real(), tx.Imag(), zero_real, zero_imag);
    const int N = width;
    const double ar = a.real();
    const double ai = a.imag();
    const auto *TXR = tx.Real().Read();
    const auto *TXI = tx.Imag().Read();
    auto *YR = yr.ReadWrite();
    auto *YI = yi.ReadWrite();
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   YR[i] += ar * TXR[i] - ai * TXI[i];
                   YI[i] += ai * TXR[i] + ar * TXI[i];
                 });
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

void ComplexWrapperOperator::AddMultHermitianTranspose(const Vector &xr, const Vector &xi,
                                                       Vector &yr, Vector &yi,
                                                       const std::complex<double> a,
                                                       bool zero_real, bool zero_imag) const
{
  if (a.real() != 0.0 && a.imag() != 0.0)
  {
    tx.SetSize(width);
    MultHermitianTranspose(xr, xi, tx.Real(), tx.Imag(), zero_real, zero_imag);
    const int N = width;
    const double ar = a.real();
    const double ai = a.imag();
    const auto *TXR = tx.Real().Read();
    const auto *TXI = tx.Imag().Read();
    auto *YR = yr.ReadWrite();
    auto *YI = yi.ReadWrite();
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   YR[i] += ar * TXR[i] - ai * TXI[i];
                   YI[i] += ai * TXR[i] + ar * TXI[i];
                 });
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

SumOperator::SumOperator(const Operator &op, double c) : Operator(op.Height(), op.Width())
{
  AddOperator(op, c);
}

void SumOperator::AddOperator(const Operator &op, double c)
{
  MFEM_VERIFY(op.Height() == height && op.Width() == width,
              "Invalid Operator dimensions for SumOperator!");
  ops.emplace_back(&op, c);
}

void SumOperator::Mult(const Vector &x, Vector &y) const
{
  if (ops.size() == 1 && ops[0].second == 1.0)
  {
    return ops[0].first->Mult(x, y);
  }
  y = 0.0;
  AddMult(x, y);
}

void SumOperator::MultTranspose(const Vector &x, Vector &y) const
{
  if (ops.size() == 1 && ops[0].second == 1.0)
  {
    return ops[0].first->MultTranspose(x, y);
  }
  y = 0.0;
  AddMultTranspose(x, y);
}

void SumOperator::AddMult(const Vector &x, Vector &y, const double a) const
{
  for (const auto &[op, c] : ops)
  {
    op->AddMult(x, y, a * c);
  }
}

void SumOperator::AddMultTranspose(const Vector &x, Vector &y, const double a) const
{
  for (const auto &[op, c] : ops)
  {
    op->AddMultTranspose(x, y, a * c);
  }
}

ComplexSumOperator::ComplexSumOperator(const ComplexOperator &op, std::complex<double> c)
  : ComplexOperator(op.Height(), op.Width())
{
  AddOperator(op, c);
}

void ComplexSumOperator::AddOperator(const ComplexOperator &op, std::complex<double> c)
{
  MFEM_VERIFY(op.Height() == height && op.Width() == width,
              "Invalid Operator dimensions for ComplexSumOperator!");
  ops.emplace_back(&op, c);
}

bool ComplexSumOperator::IsReal() const
{
  for (const auto &[op, c] : ops)
  {
    if (!op->IsReal())
    {
      return false;
    }
  }
  return true;
}

bool ComplexSumOperator::IsImag() const
{
  for (const auto &[op, c] : ops)
  {
    if (!op->IsImag())
    {
      return false;
    }
  }
  return true;
}

void ComplexSumOperator::Mult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                              bool zero_real, bool zero_imag) const
{
  if (ops.Size() == 1 && ops[0].second == 1.0)
  {
    return ops[0].first->Mult(xr, xi, yr, yi, zero_real, zero_imag);
  }
  yr = 0.0;
  yi = 0.0;
  AddMult(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
}

void ComplexSumOperator::MultTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                       Vector &yi, bool zero_real, bool zero_imag) const
{
  if (ops.Size() == 1 && ops[0].second == 1.0)
  {
    return ops[0].first->MultTranspose(xr, xi, yr, yi, zero_real, zero_imag);
  }
  yr = 0.0;
  yi = 0.0;
  AddMultTranspose(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
}

void ComplexSumOperator::MultHermitianTranspose(const Vector &xr, const Vector &xi,
                                                Vector &yr, Vector &yi, bool zero_real,
                                                bool zero_imag) const
{
  if (ops.Size() == 1 && ops[0].second == 1.0)
  {
    return ops[0].first->MultHermitianTranspose(xr, xi, yr, yi, zero_real, zero_imag);
  }
  yr = 0.0;
  yi = 0.0;
  AddMultHermitianTranspose(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
}

void ComplexSumOperator::AddMult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                                 const std::complex<double> a, bool zero_real,
                                 bool zero_imag) const
{
  for (const auto &[op, c] : ops)
  {
    op->AddMult(xr, xi, yr, yi, a * c, zero_real, zero_imag);
  }
}

void ComplexSumOperator::AddMultTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                          Vector &yi, const std::complex<double> a,
                                          bool zero_real, bool zero_imag) const
{
  for (const auto &[op, c] : ops)
  {
    op->AddMultTranspose(xr, xi, yr, yi, a * c, zero_real, zero_imag);
  }
}

void ComplexSumOperator::AddMultHermitianTranspose(const Vector &xr, const Vector &xi,
                                                   Vector &yr, Vector &yi,
                                                   const std::complex<double> a,
                                                   bool zero_real, bool zero_imag) const
{
  for (const auto &[op, c] : ops)
  {
    op->AddMultTranspose(xr, xi, yr, yi, a * c, zero_real, zero_imag);
  }
}

template <>
void DiagonalOperator<Operator>::Mult(const Vector &x, Vector &y) const
{
  const int N = height;
  const auto *D = d.Read();
  const auto *X = x.Read();
  auto *Y = y.Write();
  mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { Y[i] = D[i] * X[i]; });
}

template <>
void DiagonalOperator<ComplexOperator>::Mult(const ComplexVector &x, ComplexVector &y) const
{
  const int N = height;
  const auto *DR = d.Real().Read();
  const auto *DI = d.Imag().Read();
  const auto *XR = x.Real().Read();
  const auto *XI = x.Imag().Read();
  auto *YR = y.Real().Write();
  auto *YI = y.Imag().Write();
  mfem::forall(N,
               [=] MFEM_HOST_DEVICE(int i)
               {
                 YR[i] = DR[i] * XR[i] - DI[i] * XI[i];
                 YI[i] = DI[i] * XR[i] + DR[i] * XI[i];
               });
}

template <>
void DiagonalOperator<ComplexOperator>::MultHermitianTranspose(const ComplexVector &x,
                                                               ComplexVector &y) const
{
  const int N = height;
  const auto *DR = d.Real().Read();
  const auto *DI = d.Imag().Read();
  const auto *XR = x.Real().Read();
  const auto *XI = x.Imag().Read();
  auto *YR = y.Real().Write();
  auto *YI = y.Imag().Write();
  mfem::forall(N,
               [=] MFEM_HOST_DEVICE(int i)
               {
                 YR[i] = DR[i] * XR[i] + DI[i] * XI[i];
                 YI[i] = -DI[i] * XR[i] + DR[i] * XI[i];
               });
}

namespace linalg
{

double SpectralNorm(MPI_Comm comm, const Operator &A, bool sym, double tol, int max_it)
{
  ComplexWrapperOperator Ar(&A, nullptr);  // Non-owning constructor
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
  double l, l0 = 0.0;
  ComplexVector u(A.Height()), v(A.Height());
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
