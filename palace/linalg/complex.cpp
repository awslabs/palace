// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "complex.hpp"

#include <general/forall.hpp>

namespace palace
{

ComplexVector::ComplexVector(int n) : Vector(2 * n)
{
  xr_.MakeRef(*this, 0, n);
  xi_.MakeRef(*this, n, n);
}

ComplexVector::ComplexVector(const ComplexVector &x) : Vector(2 * x.Size())
{
  xr_.MakeRef(*this, 0, x.Size());
  xi_.MakeRef(*this, x.Size(), x.Size());
  Set(x.Real(), x.Imag());
}

ComplexVector::ComplexVector(const Vector &xr, const Vector &xi) : Vector(2 * xr.Size())
{
  MFEM_VERIFY(xr.Size() == xi.Size(),
              "Mismatch in dimension of real and imaginary matrix parts!");
  xr_.MakeRef(*this, 0, xr.Size());
  xi_.MakeRef(*this, xr.Size(), xr.Size());
  Set(xr, xi);
}

int ComplexVector::Size() const
{
  return Vector::Size() / 2;
}

void ComplexVector::SetSize(int n)
{
  Vector::SetSize(2 * n);
  xr_.MakeRef(*this, 0, n);
  xi_.MakeRef(*this, n, n);
}

ComplexVector &ComplexVector::operator=(const ComplexVector &y)
{
  Set(y.Real(), y.Imag());
  return *this;
}

void ComplexVector::Set(const Vector &yr, const Vector &yi)
{
  MFEM_VERIFY(yr.Size() == yi.Size() && yr.Size() == Size(),
              "Mismatch in dimension of real and imaginary matrix parts!");
  Real() = yr;
  Imag() = yi;
  RestoreReal();
  RestoreImag();
}

void ComplexVector::Conj()
{
  Imag() *= -1.0;
  RestoreImag();
}

ComplexVector &ComplexVector::operator*=(std::complex<double> s)
{
  if (s.imag() != 0.0)
  {
    const int N = Size();
    const double sr = s.real();
    const double si = s.imag();
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
  else if (s.real() != 0.0)
  {
    Real() *= s.real();
    Imag() *= s.real();
  }
  RestoreReal();
  RestoreImag();
  return *this;
}

std::complex<double> ComplexVector::Dot(const ComplexVector &y) const
{
  return {mfem::InnerProduct(Real(), y.Real()) + mfem::InnerProduct(Imag(), y.Imag()),
          mfem::InnerProduct(Imag(), y.Real()) - mfem::InnerProduct(Real(), y.Imag())};
}

std::complex<double> ComplexVector::TransposeDot(const ComplexVector &y) const
{
  return {mfem::InnerProduct(Real(), y.Real()) - mfem::InnerProduct(Imag(), y.Imag()),
          mfem::InnerProduct(Imag(), y.Real()) + mfem::InnerProduct(Real(), y.Imag())};
}

void ComplexVector::AXPY(std::complex<double> alpha, const ComplexVector &y)
{
  const int N = Size();
  const double ar = alpha.real();
  const double ai = alpha.imag();
  auto *XR = Real().ReadWrite();
  auto *XI = Imag().ReadWrite();
  const auto *YR = y.Real().Read();
  const auto *YI = y.Imag().Read();
  mfem::forall(N,
               [=] MFEM_HOST_DEVICE(int i)
               {
                 XR[i] += ar * YR[i] - ai * YI[i];
                 XI[i] += ai * YR[i] + ar * YI[i];
               });
  RestoreReal();
  RestoreImag();
}

void ComplexVector::AXPBY(std::complex<double> alpha, const ComplexVector &y,
                          std::complex<double> beta)
{
  const int N = Size();
  const double ar = alpha.real();
  const double ai = alpha.imag();
  const double br = beta.real();
  const double bi = beta.imag();
  auto *XR = Real().ReadWrite();
  auto *XI = Imag().ReadWrite();
  const auto *YR = y.Real().Read();
  const auto *YI = y.Imag().Read();
  mfem::forall(N,
               [=] MFEM_HOST_DEVICE(int i)
               {
                 const double t = bi * XR[i] + br * XI[i];
                 XR[i] = ar * YR[i] - ai * YI[i] + br * XR[i] - bi * XI[i];
                 XI[i] = ai * YR[i] + ar * YI[i] + t;
               });
  RestoreReal();
  RestoreImag();
}

void ComplexVector::AXPBYPCZ(std::complex<double> alpha, const ComplexVector &y,
                             std::complex<double> beta, const ComplexVector &z,
                             std::complex<double> gamma)
{
  const int N = Size();
  const double ar = alpha.real();
  const double ai = alpha.imag();
  const double br = beta.real();
  const double bi = beta.imag();
  const double gr = gamma.real();
  const double gi = gamma.imag();
  auto *XR = Real().ReadWrite();
  auto *XI = Imag().ReadWrite();
  const auto *YR = y.Real().Read();
  const auto *YI = y.Imag().Read();
  const auto *ZR = z.Real().Read();
  const auto *ZI = z.Imag().Read();
  mfem::forall(N,
               [=] MFEM_HOST_DEVICE(int i)
               {
                 const double t = gi * XR[i] + gr * XI[i];
                 XR[i] = ar * YR[i] - ai * YI[i] + br * ZR[i] - bi * ZI[i] + gr * XR[i] -
                         gi * XI[i];
                 XI[i] = ai * YR[i] + ar * YI[i] + bi * ZR[i] + br * ZI[i] + t;
               });
  RestoreReal();
  RestoreImag();
}

void ComplexOperator::Mult(const Vector &x, Vector &y) const
{
  MFEM_ASSERT(x.Size() == 2 * width && y.Size() == 2 * height,
              "Incompatible dimensions for ComplexOperator::Mult!");
  Vector xr, xi, yr, yi;
  xr.MakeRef(const_cast<Vector &>(x), 0, width / 2);
  xi.MakeRef(const_cast<Vector &>(x), width / 2, width / 2);
  yr.MakeRef(y, 0, height / 2);
  yi.MakeRef(y, height / 2, height / 2);
  Mult(xr, xi, yr, yi);
  yr.SyncAliasMemory(y);
  yi.SyncAliasMemory(y);
}

void ComplexOperator::MultTranspose(const Vector &x, Vector &y) const
{
  MFEM_ASSERT(x.Size() == 2 * height && y.Size() == 2 * width,
              "Incompatible dimensions for ComplexOperator::MultTranspose!");
  Vector xr, xi, yr, yi;
  xr.MakeRef(const_cast<Vector &>(x), 0, height / 2);
  xi.MakeRef(const_cast<Vector &>(x), height / 2, height / 2);
  yr.MakeRef(y, 0, width / 2);
  yi.MakeRef(y, width / 2, width / 2);
  MultTranspose(xr, xi, yr, yi);
  yr.SyncAliasMemory(y);
  yi.SyncAliasMemory(y);
}

void ComplexOperator::MultHermitianTranspose(const Vector &x, Vector &y) const
{
  MFEM_ASSERT(x.Size() == 2 * height && y.Size() == 2 * width,
              "Incompatible dimensions for ComplexOperator::MultHermitianTranspose!");
  Vector xr, xi, yr, yi;
  xr.MakeRef(const_cast<Vector &>(x), 0, height / 2);
  xi.MakeRef(const_cast<Vector &>(x), height / 2, height / 2);
  yr.MakeRef(y, 0, width / 2);
  yi.MakeRef(y, width / 2, width / 2);
  MultHermitianTranspose(xr, xi, yr, yi);
  yr.SyncAliasMemory(y);
  yi.SyncAliasMemory(y);
}

void ComplexOperator::AddMult(const Vector &x, Vector &y, const double a) const
{
  MFEM_ASSERT(x.Size() == 2 * width && y.Size() == 2 * height,
              "Incompatible dimensions for ComplexOperator::AddMult!");
  Vector xr, xi, yr, yi;
  xr.MakeRef(const_cast<Vector &>(x), 0, width / 2);
  xi.MakeRef(const_cast<Vector &>(x), width / 2, width / 2);
  yr.MakeRef(y, 0, height / 2);
  yi.MakeRef(y, height / 2, height / 2);
  AddMult(xr, xi, yr, yi, a);
  yr.SyncAliasMemory(y);
  yi.SyncAliasMemory(y);
}

void ComplexOperator::AddMultTranspose(const Vector &x, Vector &y, const double a) const
{
  MFEM_ASSERT(x.Size() == 2 * height && y.Size() == 2 * width,
              "Incompatible dimensions for ComplexOperator::AddMultTranspose!");
  Vector xr, xi, yr, yi;
  xr.MakeRef(const_cast<Vector &>(x), 0, height / 2);
  xi.MakeRef(const_cast<Vector &>(x), height / 2, height / 2);
  yr.MakeRef(y, 0, width / 2);
  yi.MakeRef(y, width / 2, width / 2);
  AddMultTranspose(xr, xi, yr, yi, a);
  yr.SyncAliasMemory(y);
  yi.SyncAliasMemory(y);
}

void ComplexOperator::AddMultHermitianTranspose(const Vector &x, Vector &y,
                                                const double a) const
{
  MFEM_ASSERT(x.Size() == 2 * height && y.Size() == 2 * width,
              "Incompatible dimensions for ComplexOperator::AddMultHermitianTranspose!");
  Vector xr, xi, yr, yi;
  xr.MakeRef(const_cast<Vector &>(x), 0, height / 2);
  xi.MakeRef(const_cast<Vector &>(x), height / 2, height / 2);
  yr.MakeRef(y, 0, width / 2);
  yi.MakeRef(y, width / 2, width / 2);
  AddMultHermitianTranspose(xr, xi, yr, yi, a);
  yr.SyncAliasMemory(y);
  yi.SyncAliasMemory(y);
}

ComplexParOperator::ComplexParOperator(std::unique_ptr<ComplexOperator> &&A,
                                       const mfem::ParFiniteElementSpace &trial_fespace,
                                       const mfem::ParFiniteElementSpace &test_fespace,
                                       bool test_restrict)
  : ComplexOperator(test_fespace.GetTrueVSize(), trial_fespace.GetTrueVSize()),
    A_(std::move(A)), trial_fespace_(trial_fespace), test_fespace_(test_fespace),
    use_R_(test_restrict), trial_dbc_tdof_list_(nullptr), test_dbc_tdof_list_(nullptr),
    diag_policy_(DiagonalPolicy::DIAG_ONE)
{
  MFEM_VERIFY(A_, "Cannot construct ComplexParOperator from an empty matrix!");
  lxr_.SetSize(A_->Width());
  lxi_.SetSize(A_->Width());
  lyr_.SetSize(A_->Height());
  lyi_.SetSize(A_->Height());
  txr_.SetSize(width);
  txi_.SetSize(width);
  if (height != width)
  {
    tyr_.SetSize(height);
    tyi_.SetSize(height);
  }
  else
  {
    tyr_.MakeRef(txr_, 0, height);
    tyi_.MakeRef(txi_, 0, height);
  }
}

void ComplexParOperator::AddMult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                                 const std::complex<double> a, bool zero_real,
                                 bool zero_imag) const
{
  MFEM_ASSERT(xr.Size() == width && xi.Size() == width && yr.Size() == height &&
                  yi.Size() == height,
              "Incompatible dimensions for ComplexParOperator::AddMult!");
  if (trial_dbc_tdof_list_)
  {
    txr_ = xr;
    txi_ = xi;
    txr_.SetSubVector(*trial_dbc_tdof_list_, 0.0);
    txi_.SetSubVector(*trial_dbc_tdof_list_, 0.0);
  }
  if (!zero_real)
  {
    trial_fespace_.GetProlongationMatrix()->Mult(trial_dbc_tdof_list_ ? txr_ : xr, lxr_);
  }
  if (!zero_imag)
  {
    trial_fespace_.GetProlongationMatrix()->Mult(trial_dbc_tdof_list_ ? txi_ : xi, lxi_);
  }

  // Apply the operator on the L-vector.
  lyr_ = 0.0;
  lyi_ = 0.0;
  A_->AddMult(lxr_, lxi_, lyr_, lyi_, a, zero_real, zero_imag);

  if (test_dbc_tdof_list_)
  {
    if (!use_R_)
    {
      test_fespace_.GetProlongationMatrix()->MultTranspose(lyr_, tyr_);
      test_fespace_.GetProlongationMatrix()->MultTranspose(lyi_, tyi_);
    }
    else
    {
      test_fespace_.GetRestrictionMatrix()->Mult(lyr_, tyr_);
      test_fespace_.GetRestrictionMatrix()->Mult(lyi_, tyi_);
    }
    {
      const int N = test_dbc_tdof_list_->Size();
      auto idx = test_dbc_tdof_list_->Read();
      auto TYR = tyr_.ReadWrite();
      auto TYI = tyi_.ReadWrite();
      if (diag_policy_ == DiagonalPolicy::DIAG_ONE && height == width)
      {
        const auto *XR = xr.Read();
        const auto *XI = xi.Read();
        mfem::forall(N,
                     [=] MFEM_HOST_DEVICE(int i)
                     {
                       const int id = idx[i];
                       TYR[id] = XR[id];
                       TYI[id] = XI[id];
                     });
      }
      else if (diag_policy_ == DiagonalPolicy::DIAG_ZERO || height != width)
      {
        mfem::forall(N,
                     [=] MFEM_HOST_DEVICE(int i)
                     {
                       const int id = idx[i];
                       TYR[id] = TYI[id] = 0.0;
                     });
      }
      else
      {
        MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
      }
    }
    yr += tyr_;
    yi += tyi_;
  }
  else
  {
    if (!use_R_)
    {
      test_fespace_.GetProlongationMatrix()->MultTranspose(lyr_, yr);
      test_fespace_.GetProlongationMatrix()->MultTranspose(lyi_, yi);
    }
    else
    {
      test_fespace_.GetRestrictionMatrix()->Mult(lyr_, yr);
      test_fespace_.GetRestrictionMatrix()->Mult(lyi_, yi);
    }
  }
}

void ComplexParOperator::AddMultTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                          Vector &yi, const std::complex<double> a,
                                          bool zero_real, bool zero_imag) const
{
  MFEM_ASSERT(xr.Size() == height && xi.Size() == height && yr.Size() == width &&
                  yi.Size() == width,
              "Incompatible dimensions for ComplexParOperator::AddMultTranspose!");
  if (test_dbc_tdof_list_)
  {
    tyr_ = xr;
    tyi_ = xi;
    tyr_.SetSubVector(*test_dbc_tdof_list_, 0.0);
    tyi_.SetSubVector(*test_dbc_tdof_list_, 0.0);
  }
  if (!use_R_)
  {
    if (!zero_real)
    {
      test_fespace_.GetProlongationMatrix()->Mult(test_dbc_tdof_list_ ? tyr_ : xr, lyr_);
    }
    if (!zero_imag)
    {
      test_fespace_.GetProlongationMatrix()->Mult(test_dbc_tdof_list_ ? tyi_ : xi, lyi_);
    }
  }
  else
  {
    if (!zero_real)
    {
      test_fespace_.GetRestrictionMatrix()->MultTranspose(test_dbc_tdof_list_ ? tyr_ : xr,
                                                          lyr_);
    }
    if (!zero_imag)
    {
      test_fespace_.GetRestrictionMatrix()->MultTranspose(test_dbc_tdof_list_ ? tyi_ : xi,
                                                          lyi_);
    }
  }

  // Apply the operator on the L-vector.
  lxr_ = 0.0;
  lxi_ = 0.0;
  A_->AddMultTranspose(lyr_, lyi_, lxr_, lxi_, a, zero_real, zero_imag);

  if (trial_dbc_tdof_list_)
  {
    trial_fespace_.GetProlongationMatrix()->MultTranspose(lxr_, txr_);
    trial_fespace_.GetProlongationMatrix()->MultTranspose(lxi_, txi_);
    {
      const int N = trial_dbc_tdof_list_->Size();
      auto idx = trial_dbc_tdof_list_->Read();
      auto TXR = txr_.ReadWrite();
      auto TXI = txi_.ReadWrite();
      if (diag_policy_ == DiagonalPolicy::DIAG_ONE && height == width)
      {
        const auto *XR = xr.Read();
        const auto *XI = xi.Read();
        mfem::forall(N,
                     [=] MFEM_HOST_DEVICE(int i)
                     {
                       const int id = idx[i];
                       TXR[id] = XR[id];
                       TXI[id] = XI[id];
                     });
      }
      else if (diag_policy_ == DiagonalPolicy::DIAG_ZERO || height != width)
      {
        mfem::forall(N,
                     [=] MFEM_HOST_DEVICE(int i)
                     {
                       const int id = idx[i];
                       TXR[id] = TXI[id] = 0.0;
                     });
      }
      else
      {
        MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
      }
    }
    yr += txr_;
    yi += txi_;
  }
  else
  {
    trial_fespace_.GetProlongationMatrix()->AddMultTranspose(lxr_, yr);
    trial_fespace_.GetProlongationMatrix()->AddMultTranspose(lxi_, yi);
  }
}

void ComplexParOperator::AddMultHermitianTranspose(const Vector &xr, const Vector &xi,
                                                   Vector &yr, Vector &yi,
                                                   const std::complex<double> a,
                                                   bool zero_real, bool zero_imag) const
{
  MFEM_ASSERT(xr.Size() == height && xi.Size() == height && yr.Size() == width &&
                  yi.Size() == width,
              "Incompatible dimensions for ComplexParOperator::AddMultHermitianTranspose!");
  if (test_dbc_tdof_list_)
  {
    tyr_ = xr;
    tyi_ = xi;
    tyr_.SetSubVector(*test_dbc_tdof_list_, 0.0);
    tyi_.SetSubVector(*test_dbc_tdof_list_, 0.0);
  }
  if (!use_R_)
  {
    if (!zero_real)
    {
      test_fespace_.GetProlongationMatrix()->Mult(test_dbc_tdof_list_ ? tyr_ : xr, lyr_);
    }
    if (!zero_imag)
    {
      test_fespace_.GetProlongationMatrix()->Mult(test_dbc_tdof_list_ ? tyi_ : xi, lyi_);
    }
  }
  else
  {
    if (!zero_real)
    {
      test_fespace_.GetRestrictionMatrix()->MultTranspose(test_dbc_tdof_list_ ? tyr_ : xr,
                                                          lyr_);
    }
    if (!zero_imag)
    {
      test_fespace_.GetRestrictionMatrix()->MultTranspose(test_dbc_tdof_list_ ? tyi_ : xi,
                                                          lyi_);
    }
  }

  // Apply the operator on the L-vector.
  lxr_ = 0.0;
  lxi_ = 0.0;
  A_->AddMultHermitianTranspose(lyr_, lyi_, lxr_, lxi_, a, zero_real, zero_imag);

  if (trial_dbc_tdof_list_)
  {
    trial_fespace_.GetProlongationMatrix()->MultTranspose(lxr_, txr_);
    trial_fespace_.GetProlongationMatrix()->MultTranspose(lxi_, txi_);
    {
      const int N = trial_dbc_tdof_list_->Size();
      auto idx = trial_dbc_tdof_list_->Read();
      auto TXR = txr_.ReadWrite();
      auto TXI = txi_.ReadWrite();
      if (diag_policy_ == DiagonalPolicy::DIAG_ONE && height == width)
      {
        const auto *XR = xr.Read();
        const auto *XI = xi.Read();
        mfem::forall(N,
                     [=] MFEM_HOST_DEVICE(int i)
                     {
                       const int id = idx[i];
                       TXR[id] = XR[id];
                       TXI[id] = XI[id];
                     });
      }
      else if (diag_policy_ == DiagonalPolicy::DIAG_ZERO || height != width)
      {
        mfem::forall(N,
                     [=] MFEM_HOST_DEVICE(int i)
                     {
                       const int id = idx[i];
                       TXR[id] = TXI[id] = 0.0;
                     });
      }
      else
      {
        MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
      }
    }
    yr += txr_;
    yi += txi_;
  }
  else
  {
    trial_fespace_.GetProlongationMatrix()->AddMultTranspose(lxr_, yr);
    trial_fespace_.GetProlongationMatrix()->AddMultTranspose(lxi_, yi);
  }
}

ComplexWrapperOperator::ComplexWrapperOperator(std::unique_ptr<Operator> &&Ar,
                                               std::unique_ptr<Operator> &&Ai)
  : ComplexOperator(Ar ? Ar->Height() : Ai->Height(), Ar ? Ar->Width() : Ai->Width()),
    Ar_(std::move(Ar)), Ai_(std::move(Ai))
{
  MFEM_VERIFY(Ar_ || Ai_, "Cannot construct ComplexWrapperOperator from an empty matrix!");
  MFEM_VERIFY((!Ar_ || !Ai_) ||
                  (Ar_->Height() == Ai_->Height() && Ar_->Width() == Ai_->Width()),
              "Mismatch in dimension of real and imaginary matrix parts!");
  txr_.SetSize(width);
  txi_.SetSize(width);
  if (height != width)
  {
    tyr_.SetSize(height);
    tyi_.SetSize(height);
  }
  else
  {
    tyr_.MakeRef(txr_, 0, height);
    tyi_.MakeRef(txi_, 0, height);
  }
}

void ComplexWrapperOperator::Mult(const Vector &xr, const Vector &xi, Vector &yr,
                                  Vector &yi, bool zero_real, bool zero_imag) const
{
  if (Ar_)
  {
    if (!zero_real)
    {
      Ar_->Mult(xr, yr);
    }
    if (!zero_imag)
    {
      Ar_->Mult(xi, yi);
    }
  }
  else
  {
    yr = 0.0;
    yi = 0.0;
  }
  if (Ai_)
  {
    if (!zero_imag)
    {
      Ai_->AddMult(xi, yr, -1.0);
    }
    if (!zero_real)
    {
      Ai_->AddMult(xr, yi, 1.0);
    }
  }
}

void ComplexWrapperOperator::MultTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                           Vector &yi, bool zero_real, bool zero_imag) const
{
  if (Ar_)
  {
    if (!zero_real)
    {
      Ar_->MultTranspose(xr, yr);
    }
    if (!zero_imag)
    {
      Ar_->MultTranspose(xi, yi);
    }
  }
  else
  {
    yr = 0.0;
    yi = 0.0;
  }
  if (Ai_)
  {
    if (!zero_imag)
    {
      Ai_->AddMultTranspose(xi, yr, -1.0);
    }
    if (!zero_real)
    {
      Ai_->AddMultTranspose(xr, yi, 1.0);
    }
  }
}

void ComplexWrapperOperator::MultHermitianTranspose(const Vector &xr, const Vector &xi,
                                                    Vector &yr, Vector &yi, bool zero_real,
                                                    bool zero_imag) const
{
  if (Ar_)
  {
    if (!zero_real)
    {
      Ar_->MultTranspose(xr, yr);
    }
    if (!zero_imag)
    {
      Ar_->MultTranspose(xi, yi);
    }
  }
  else
  {
    yr = 0.0;
    yi = 0.0;
  }
  if (Ai_)
  {
    if (!zero_imag)
    {
      Ai_->AddMultTranspose(xi, yr, 1.0);
    }
    if (!zero_real)
    {
      Ai_->AddMultTranspose(xr, yi, -1.0);
    }
  }
}

void ComplexWrapperOperator::AddMult(const Vector &xr, const Vector &xi, Vector &yr,
                                     Vector &yi, const std::complex<double> a,
                                     bool zero_real, bool zero_imag) const
{
  if (a.real() != 0.0 && a.imag() != 0.0)
  {
    Mult(xr, xi, tyr_, tyi_, zero_real, zero_imag);
    const int N = height;
    const double ar = a.real();
    const double ai = a.imag();
    auto *YR = yr.ReadWrite();
    auto *YI = yi.ReadWrite();
    const auto *TYR = tyr_.Read();
    const auto *TYI = tyi_.Read();
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   YR[i] += ar * TYR[i] - ai * TYI[i];
                   YI[i] += ai * TYR[i] + ar * TYI[i];
                 });
  }
  else if (a.real() != 0.0)
  {
    if (Ar_)
    {
      if (!zero_real)
      {
        Ar_->AddMult(xr, yr, a.real());
      }
      if (!zero_imag)
      {
        Ar_->AddMult(xi, yi, a.real());
      }
    }
    if (Ai_)
    {
      if (!zero_imag)
      {
        Ai_->AddMult(xi, yr, -a.real());
      }
      if (!zero_real)
      {
        Ai_->AddMult(xr, yi, a.real());
      }
    }
  }
  else if (a.imag() != 0.0)
  {
    if (Ar_)
    {
      if (!zero_real)
      {
        Ar_->AddMult(xr, yi, a.imag());
      }
      if (!zero_imag)
      {
        Ar_->AddMult(xi, yr, -a.imag());
      }
    }
    if (Ai_)
    {
      if (!zero_imag)
      {
        Ai_->AddMult(xi, yi, -a.imag());
      }
      if (!zero_real)
      {
        Ai_->AddMult(xr, yr, -a.imag());
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
    MultTranspose(xr, xi, txr_, txi_, zero_real, zero_imag);
    const int N = width;
    const double ar = a.real();
    const double ai = a.imag();
    auto *YR = yr.ReadWrite();
    auto *YI = yi.ReadWrite();
    const auto *TXR = txr_.Read();
    const auto *TXI = txi_.Read();
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   YR[i] += ar * TXR[i] - ai * TXI[i];
                   YI[i] += ai * TXR[i] + ar * TXI[i];
                 });
  }
  else if (a.real() != 0.0)
  {
    if (Ar_)
    {
      if (!zero_real)
      {
        Ar_->AddMultTranspose(xr, yr, a.real());
      }
      if (!zero_imag)
      {
        Ar_->AddMultTranspose(xi, yi, a.real());
      }
    }
    if (Ai_)
    {
      if (!zero_imag)
      {
        Ai_->AddMultTranspose(xi, yr, -a.real());
      }
      if (!zero_real)
      {
        Ai_->AddMultTranspose(xr, yi, a.real());
      }
    }
  }
  else if (a.imag() != 0.0)
  {
    if (Ar_)
    {
      if (!zero_real)
      {
        Ar_->AddMultTranspose(xr, yi, a.imag());
      }
      if (!zero_imag)
      {
        Ar_->AddMultTranspose(xi, yr, -a.imag());
      }
    }
    if (Ai_)
    {
      if (!zero_imag)
      {
        Ai_->AddMultTranspose(xi, yi, -a.imag());
      }
      if (!zero_real)
      {
        Ai_->AddMultTranspose(xr, yr, -a.imag());
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
    MultHermitianTranspose(xr, xi, txr_, txi_, zero_real, zero_imag);
    const int N = width;
    const double ar = a.real();
    const double ai = a.imag();
    auto *YR = yr.ReadWrite();
    auto *YI = yi.ReadWrite();
    const auto *TXR = txr_.Read();
    const auto *TXI = txi_.Read();
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   YR[i] += ar * TXR[i] - ai * TXI[i];
                   YI[i] += ai * TXR[i] + ar * TXI[i];
                 });
  }
  else if (a.real() != 0.0)
  {
    if (Ar_)
    {
      if (!zero_real)
      {
        Ar_->AddMultTranspose(xr, yr, a.real());
      }
      if (!zero_imag)
      {
        Ar_->AddMultTranspose(xi, yi, a.real());
      }
    }
    if (Ai_)
    {
      if (!zero_imag)
      {
        Ai_->AddMultTranspose(xi, yr, a.real());
      }
      if (!zero_real)
      {
        Ai_->AddMultTranspose(xr, yi, -a.real());
      }
    }
  }
  else if (a.imag() != 0.0)
  {
    if (Ar_)
    {
      if (!zero_real)
      {
        Ar_->AddMultTranspose(xr, yi, a.imag());
      }
      if (!zero_imag)
      {
        Ar_->AddMultTranspose(xi, yr, -a.imag());
      }
    }
    if (Ai_)
    {
      if (!zero_imag)
      {
        Ai_->AddMultTranspose(xi, yi, a.imag());
      }
      if (!zero_real)
      {
        Ai_->AddMultTranspose(xr, yr, a.imag());
      }
    }
  }
}

}  // namespace palace
