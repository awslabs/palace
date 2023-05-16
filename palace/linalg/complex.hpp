// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_COMPLEX_HPP
#define PALACE_LINALG_COMPLEX_HPP

#include <complex>
#include <memory>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Vector and operator classes for complex-valued linear algebra.
//

// A complex-valued vector represented as two real vectors, one for each component. The
// value returned by the vector size is twice the actual complex-valued size.
class ComplexVector : public Vector
{
private:
  Vector xr_, xi_;

public:
  // Create a vector with the given size. The provided size should be the real-valued size,
  // twice the actual complex-valued size, in order to agree with ComplexOperator::Height().
  ComplexVector(int n = 0);

  // Copy constructor.
  ComplexVector(const ComplexVector &x);

  // Copy constructor from separately provided real and imaginary parts.
  ComplexVector(const Vector &xr, const Vector &xi);

  // Copy constructor from an array of complex values. The size provided should be the
  // length of the array x, which is half the resulting real-valued vector size.
  ComplexVector(const std::complex<double> *px, int n);

  // Set the size of the vector. The provided size should be the real-valued size, twice the
  // actual complex-valued size, in order to agree with ComplexOperator::Height(). See the
  // notes for Vector::SetSize for behavior in the cases where n is less than or greater
  // than Size() or Capacity().
  void SetSize(int n);

  // Get const access to the real and imaginary vector parts. Assumes that these are
  // synchronized following a Sync() call.
  const Vector &Real() const { return xr_; }
  const Vector &Imag() const { return xi_; }

  // Get access to the real and imaginary vector parts with required synchronization with
  // the underlying storage.
  Vector &Real()
  {
    xr_.SyncMemory(*this);
    return xr_;
  }
  Vector &Imag()
  {
    xi_.SyncMemory(*this);
    return xi_;
  }
  void RestoreReal() { xr_.SyncAliasMemory(*this); }
  void RestoreImag() { xi_.SyncAliasMemory(*this); }

  // Copy assignment operator. This should probably not be used to modify the size of the
  // vector.
  ComplexVector &operator=(const ComplexVector &y);

  // Copy assignment from separately provided real and imaginary parts.
  void Set(const Vector &yr, const Vector &yi);

  // Copy assignment from an array of complex values. The size provided should be the length
  // of the array x, which is half the real-valued vector size.
  void Set(const std::complex<double> *py, int n);

  // Copy the vector into an array of complex values. The size provided should be the length
  // of the array y, which is half the real-valued vector size.
  void Get(std::complex<double> *py, int n) const;

  // Replace entries with complex conjugate.
  void Conj();

  // Set all entries equal to s.
  ComplexVector &operator=(std::complex<double> s);

  // Scale all entries by s.
  ComplexVector &operator*=(std::complex<double> s);

  // Vector dot product (yᴴ x) or indefinite dot product (yᵀ x) for complex vectors.
  std::complex<double> Dot(const ComplexVector &y) const;
  std::complex<double> TransposeDot(const ComplexVector &y) const;

  // In-place addition x += alpha * y.
  void AXPY(std::complex<double> alpha, const ComplexVector &y);

  // In-place addition x = alpha * y + beta * x.
  void AXPBY(std::complex<double> alpha, const ComplexVector &y, std::complex<double> beta);

  // In-place addition x = alpha * y + beta * z + gamma * x.
  void AXPBYPCZ(std::complex<double> alpha, const ComplexVector &y,
                std::complex<double> beta, const ComplexVector &z,
                std::complex<double> gamma);

  // Update the memory location of the real and imaginary parts to match the underlying
  // storage, or vice versa.
  void Sync()
  {
    xr_.SyncMemory(*this);
    xi_.SyncMemory(*this);
  }
  void SyncAlias()
  {
    xr_.SyncAliasMemory(*this);
    xi_.SyncAliasMemory(*this);
  }
};

// Abstract base class for complex-valued operators. The values returned by the operator
// height and width are twice the actual complex-valued size.
class ComplexOperator : public Operator
{
public:
  ComplexOperator(int s) : Operator(2 * s) {}
  ComplexOperator(int h, int w) : Operator(2 * h, 2 * w) {}

  // Test whether or not the operator is purely real or imaginary.
  virtual bool IsReal() const = 0;
  virtual bool IsImag() const = 0;

  // Get access to the real and imaginary operator parts.
  virtual const Operator &Real() const
  {
    MFEM_ABORT("Real() is not implemented for base class ComplexOperator!");
    return *this;
  }
  virtual Operator &Real()
  {
    MFEM_ABORT("Real() is not implemented for base class ComplexOperator!");
    return *this;
  }
  virtual const Operator &Imag() const
  {
    MFEM_ABORT("Imag() is not implemented for base class ComplexOperator!");
    return *this;
  }
  virtual Operator &Imag()
  {
    MFEM_ABORT("Imag() is not implemented for base class ComplexOperator!");
    return *this;
  }

  void Mult(const Vector &x, Vector &y) const override;

  void Mult(const ComplexVector &x, ComplexVector &y) const
  {
    Mult(x.Real(), x.Imag(), y.Real(), y.Imag());
  }

  virtual void Mult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                    bool zero_real = false, bool zero_imag = false) const = 0;

  void MultTranspose(const Vector &x, Vector &y) const override;

  void MultTranspose(const ComplexVector &x, ComplexVector &y) const
  {
    MultTranspose(x.Real(), x.Imag(), y.Real(), y.Imag());
  }

  virtual void MultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                             bool zero_real = false, bool zero_imag = false) const = 0;

  void MultHermitianTranspose(const Vector &x, Vector &y) const;

  void MultHermitianTranspose(const ComplexVector &x, ComplexVector &y) const
  {
    MultHermitianTranspose(x.Real(), x.Imag(), y.Real(), y.Imag());
  }

  virtual void MultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                      Vector &yi, bool zero_real = false,
                                      bool zero_imag = false) const = 0;

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override;

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const
  {
    AddMult(x.Real(), x.Imag(), y.Real(), y.Imag(), a);
  }

  virtual void AddMult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                       const std::complex<double> a = 1.0, bool zero_real = false,
                       bool zero_imag = false) const = 0;

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override;

  void AddMultTranspose(const ComplexVector &x, ComplexVector &y,
                        const std::complex<double> a = 1.0) const
  {
    AddMultTranspose(x.Real(), x.Imag(), y.Real(), y.Imag(), a);
  }

  virtual void AddMultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                                const std::complex<double> a = 1.0, bool zero_real = false,
                                bool zero_imag = false) const = 0;

  void AddMultHermitianTranspose(const Vector &x, Vector &y, const double a = 1.0) const;

  void AddMultHermitianTranspose(const ComplexVector &x, ComplexVector &y,
                                 const std::complex<double> a = 1.0) const
  {
    AddMultHermitianTranspose(x.Real(), x.Imag(), y.Real(), y.Imag(), a);
  }

  virtual void AddMultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                         Vector &yi, const std::complex<double> a = 1.0,
                                         bool zero_real = false,
                                         bool zero_imag = false) const = 0;
};

// A parallel complex-valued operator represented by RᵀAP for complex-valued A, constructed
// through the actions of Rᵀ, A, and P with possible eliminated essential BC.
class ComplexParOperator : public ComplexOperator
{
private:
  std::unique_ptr<ComplexOperator> A_;
  const mfem::ParFiniteElementSpace &trial_fespace_, &test_fespace_;
  const bool use_R_;

  // Lists of constrained essential boundary true dofs for elimination.
  mutable const mfem::Array<int> *trial_dbc_tdof_list_, *test_dbc_tdof_list_;

  // Diagonal policy for constrained true dofs.
  DiagonalPolicy diag_policy_;

  // Temporary storage for operator application.
  mutable Vector lxr_, lxi_, lyr_, lyi_, txr_, txi_, tyr_, tyi_;

public:
  // Construct the complex-valued parallel operator, inheriting ownership of the local
  // operator.
  ComplexParOperator(std::unique_ptr<ComplexOperator> &&A,
                     const mfem::ParFiniteElementSpace &trial_fespace,
                     const mfem::ParFiniteElementSpace &test_fespace,
                     bool test_restrict = false);

  // Get access to the underlying local (L-vector) operator.
  const ComplexOperator &LocalOperator() const
  {
    MFEM_VERIFY(A_, "No local matrix available for ComplexParOperator::LocalOperator!");
    return *A_;
  }

  // Set essential boundary condition true dofs for square operators.
  void SetEssentialTrueDofs(const mfem::Array<int> &dbc_tdof_list,
                            DiagonalPolicy diag_policy)
  {
    MFEM_VERIFY(height == width, "Set essential true dofs for both test and trial spaces "
                                 "for rectangular ComplexParOperator!");
    trial_dbc_tdof_list_ = &dbc_tdof_list;
    test_dbc_tdof_list_ = &dbc_tdof_list;
    diag_policy_ = diag_policy;
  }

  // Set essential boundary condition true dofs for rectangular operators.
  void SetEssentialTrueDofs(const mfem::Array<int> *trial_dbc_tdof_list,
                            const mfem::Array<int> *test_dbc_tdof_list,
                            DiagonalPolicy diag_policy)
  {
    MFEM_VERIFY(diag_policy == DiagonalPolicy::DIAG_ZERO,
                "Essential boundary condition true dof elimination for rectangular "
                "ParOperator only supports DiagonalPolicy::DIAG_ZERO!");
    trial_dbc_tdof_list_ = trial_dbc_tdof_list;
    test_dbc_tdof_list_ = test_dbc_tdof_list;
    diag_policy_ = diag_policy;
  }

  // Get the essential boundary condition true dofs associated with the operator. May be
  // nullptr.
  const mfem::Array<int> *GetEssentialTrueDofs() const
  {
    MFEM_VERIFY(trial_dbc_tdof_list_ == test_dbc_tdof_list_ && height == width,
                "GetEssentialTrueDofs should only be used for square ComplexParOperator!");
    return trial_dbc_tdof_list_;
  }

  // Set the diagonal policy for the operator.
  void SetDiagonalPolicy(DiagonalPolicy diag_policy) { diag_policy_ = diag_policy; }

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return trial_fespace_.GetComm(); }

  bool IsReal() const override { return A_->IsReal(); }
  bool IsImag() const override { return A_->IsImag(); }

  const Operator &Real() const override { return A_->Real(); }
  Operator &Real() override { return A_->Real(); }
  const Operator &Imag() const override { return A_->Imag(); }
  Operator &Imag() override { return A_->Imag(); }

  using ComplexOperator::AddMult;
  using ComplexOperator::AddMultHermitianTranspose;
  using ComplexOperator::AddMultTranspose;
  using ComplexOperator::Mult;
  using ComplexOperator::MultHermitianTranspose;
  using ComplexOperator::MultTranspose;

  void Mult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
            bool zero_real = false, bool zero_imag = false) const override
  {
    yr = 0.0;
    yi = 0.0;
    AddMult(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
  }

  void MultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                     bool zero_real = false, bool zero_imag = false) const override
  {
    yr = 0.0;
    yi = 0.0;
    AddMultTranspose(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
  }

  void MultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                              bool zero_real = false, bool zero_imag = false) const override
  {
    yr = 0.0;
    yi = 0.0;
    AddMultHermitianTranspose(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
  }

  void AddMult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
               const std::complex<double> a = 1.0, bool zero_real = false,
               bool zero_imag = false) const override;

  void AddMultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                        const std::complex<double> a = 1.0, bool zero_real = false,
                        bool zero_imag = false) const override;

  void AddMultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                                 const std::complex<double> a = 1.0, bool zero_real = false,
                                 bool zero_imag = false) const override;
};

// A complex-valued operator represented using a block 2x2 equivalent-real formulation.
class ComplexWrapperOperator : public ComplexOperator
{
private:
  std::unique_ptr<Operator> Ar_, Ai_;

  // Temporary storage for operator application.
  mutable Vector txr_, txi_, tyr_, tyi_;

public:
  // Construct a complex operator which inherits ownershipt of the input real and imaginary
  // parts.
  ComplexWrapperOperator(std::unique_ptr<Operator> &&Ar, std::unique_ptr<Operator> &&Ai);

  bool IsReal() const override { return Ai_ == nullptr; }
  bool IsImag() const override { return Ar_ == nullptr; }

  const Operator &Real() const override { return *Ar_; }
  Operator &Real() override { return *Ar_; }
  const Operator &Imag() const override { return *Ai_; }
  Operator &Imag() override { return *Ai_; }

  using ComplexOperator::AddMult;
  using ComplexOperator::AddMultHermitianTranspose;
  using ComplexOperator::AddMultTranspose;
  using ComplexOperator::Mult;
  using ComplexOperator::MultHermitianTranspose;
  using ComplexOperator::MultTranspose;

  void Mult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
            bool zero_real = false, bool zero_imag = false) const override;

  void MultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                     bool zero_real = false, bool zero_imag = false) const override;

  void MultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                              bool zero_real = false,
                              bool zero_imag = false) const override;

  void AddMult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
               const std::complex<double> a = 1.0, bool zero_real = false,
               bool zero_imag = false) const override;

  void AddMultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                        const std::complex<double> a = 1.0, bool zero_real = false,
                        bool zero_imag = false) const override;

  void AddMultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                                 const std::complex<double> a = 1.0, bool zero_real = false,
                                 bool zero_imag = false) const override;
};

// Wrap a sequence of operators of the same dimensions and optional coefficients.
class ComplexSumOperator : public ComplexOperator
{
private:
  std::vector<std::pair<const ComplexOperator *, std::complex<double>>> ops_;

public:
  ComplexSumOperator(int s) : ComplexOperator(s) {}
  ComplexSumOperator(int h, int w) : ComplexOperator(h, w) {}
  ComplexSumOperator(const ComplexOperator &op, std::complex<double> c = 1.0)
    : ComplexOperator(op.Height(), op.Width())
  {
    AddOperator(op, c);
  }

  void AddOperator(const ComplexOperator &op, std::complex<double> c = 1.0)
  {
    MFEM_VERIFY(op.Height() == height && op.Width() == width,
                "Invalid Operator dimensions for ComplexSumOperator!");
    ops_.emplace_back(&op, c);
  }

  bool IsReal() const override
  {
    for (const auto &[op, c] : ops_)
    {
      if (!op->IsReal())
      {
        return false;
      }
    }
    return true;
  }

  bool IsImag() const override
  {
    for (const auto &[op, c] : ops_)
    {
      if (!op->IsImag())
      {
        return false;
      }
    }
    return true;
  }

  using ComplexOperator::AddMult;
  using ComplexOperator::AddMultHermitianTranspose;
  using ComplexOperator::AddMultTranspose;
  using ComplexOperator::Mult;
  using ComplexOperator::MultHermitianTranspose;
  using ComplexOperator::MultTranspose;

  void Mult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
            bool zero_real = false, bool zero_imag = false) const override
  {
    yr = 0.0;
    yi = 0.0;
    AddMult(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
  }

  void MultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                     bool zero_real = false, bool zero_imag = false) const override
  {
    yr = 0.0;
    yi = 0.0;
    AddMultTranspose(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
  }

  void MultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                              bool zero_real = false, bool zero_imag = false) const override
  {
    yr = 0.0;
    yi = 0.0;
    AddMultHermitianTranspose(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
  }

  void AddMult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
               const std::complex<double> a = 1.0, bool zero_real = false,
               bool zero_imag = false) const override
  {
    for (const auto &[op, c] : ops_)
    {
      op->AddMult(xr, xi, yr, yi, a * c, zero_real, zero_imag);
    }
  }

  void AddMultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                        const std::complex<double> a = 1.0, bool zero_real = false,
                        bool zero_imag = false) const override
  {
    for (const auto &[op, c] : ops_)
    {
      op->AddMultTranspose(xr, xi, yr, yi, a * c, zero_real, zero_imag);
    }
  }

  void AddMultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                                 const std::complex<double> a = 1.0, bool zero_real = false,
                                 bool zero_imag = false) const override
  {
    for (const auto &[op, c] : ops_)
    {
      op->AddMultTranspose(xr, xi, yr, yi, a * c, zero_real, zero_imag);
    }
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_COMPLEX_HPP
