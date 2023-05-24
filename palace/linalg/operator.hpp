// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_OPERATOR_HPP
#define PALACE_LINALG_OPERATOR_HPP

#include <complex>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>
#include "linalg/vector.hpp"

namespace palace
{

//
// Functionality extending mfem::Operator from MFEM.
//

// Abstract base class for complex-valued operators.
class ComplexOperator
{
protected:
  // The size of the complex-valued operator.
  int height, width;

public:
  ComplexOperator(int s) : height(s), width(s) {}
  ComplexOperator(int h, int w) : height(h), width(w) {}

  // Get the height (size of output) of the operator.
  int Height() const { return height; }

  // Get the width (size of input) of the operator.
  int Width() const { return width; }

  // Test whether or not the operator is purely real or imaginary.
  virtual bool IsReal() const;
  virtual bool IsImag() const;

  // Test whether or not we can access the real and imaginary operator parts.
  virtual bool HasReal() const;
  virtual bool HasImag() const;

  // Get access to the real and imaginary operator parts.
  virtual const Operator *Real() const;
  virtual Operator *Real();
  virtual const Operator *Imag() const;
  virtual Operator *Imag();

  virtual void Mult(const ComplexVector &x, ComplexVector &y) const
  {
    Mult(x.Real(), x.Imag(), y.Real(), y.Imag());
  }

  virtual void Mult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                    bool zero_real = false, bool zero_imag = false) const = 0;

  virtual void MultTranspose(const ComplexVector &x, ComplexVector &y) const
  {
    MultTranspose(x.Real(), x.Imag(), y.Real(), y.Imag());
  }

  virtual void MultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                             bool zero_real = false, bool zero_imag = false) const;

  virtual void MultHermitianTranspose(const ComplexVector &x, ComplexVector &y) const
  {
    MultHermitianTranspose(x.Real(), x.Imag(), y.Real(), y.Imag());
  }

  virtual void MultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                      Vector &yi, bool zero_real = false,
                                      bool zero_imag = false) const;

  virtual void AddMult(const ComplexVector &x, ComplexVector &y,
                       const std::complex<double> a = 1.0) const
  {
    AddMult(x.Real(), x.Imag(), y.Real(), y.Imag(), a);
  }

  virtual void AddMult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                       const std::complex<double> a = 1.0, bool zero_real = false,
                       bool zero_imag = false) const;

  virtual void AddMultTranspose(const ComplexVector &x, ComplexVector &y,
                                const std::complex<double> a = 1.0) const
  {
    AddMultTranspose(x.Real(), x.Imag(), y.Real(), y.Imag(), a);
  }

  virtual void AddMultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                                const std::complex<double> a = 1.0, bool zero_real = false,
                                bool zero_imag = false) const;

  virtual void AddMultHermitianTranspose(const ComplexVector &x, ComplexVector &y,
                                         const std::complex<double> a = 1.0) const
  {
    AddMultHermitianTranspose(x.Real(), x.Imag(), y.Real(), y.Imag(), a);
  }

  virtual void AddMultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                         Vector &yi, const std::complex<double> a = 1.0,
                                         bool zero_real = false,
                                         bool zero_imag = false) const;
};

// A complex-valued operator represented using a block 2x2 equivalent-real formulation.
class ComplexWrapperOperator : public ComplexOperator
{
private:
  // Storage and access for real and imaginary parts of the operator.
  std::unique_ptr<Operator> data_Ar, data_Ai;
  Operator *Ar, *Ai;

  // Temporary storage for operator application.
  mutable ComplexVector tx, ty;

  ComplexWrapperOperator(std::unique_ptr<Operator> &&data_Ar,
                         std::unique_ptr<Operator> &&data_Ai, Operator *Ar, Operator *Ai);

public:
  // Construct a complex operator which inherits ownership of the input real and imaginary
  // parts.
  ComplexWrapperOperator(std::unique_ptr<Operator> &&Ar, std::unique_ptr<Operator> &&Ai);

  // Non-owning constructor.
  ComplexWrapperOperator(Operator *Ar, Operator *Ai);

  bool IsReal() const override { return Ai == nullptr; }
  bool IsImag() const override { return Ar == nullptr; }
  bool HasReal() const override { return Ar != nullptr; }
  bool HasImag() const override { return Ai != nullptr; }
  const Operator *Real() const override { return Ar; }
  Operator *Real() override { return Ar; }
  const Operator *Imag() const override { return Ai; }
  Operator *Imag() override { return Ai; }

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
class SumOperator : public Operator
{
private:
  std::vector<std::pair<const Operator *, double>> ops;

public:
  SumOperator(int s) : Operator(s) {}
  SumOperator(int h, int w) : Operator(h, w) {}
  SumOperator(const Operator &op, double c = 1.0);

  void AddOperator(const Operator &op, double c = 1.0);

  void Mult(const Vector &x, Vector &y) const override;

  void MultTranspose(const Vector &x, Vector &y) const override;

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override;

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override;
};

// Wrap a sequence of operators of the same dimensions and optional coefficients.
class ComplexSumOperator : public ComplexOperator
{
private:
  std::vector<std::pair<const ComplexOperator *, std::complex<double>>> ops;

public:
  ComplexSumOperator(int s) : ComplexOperator(s) {}
  ComplexSumOperator(int h, int w) : ComplexOperator(h, w) {}
  ComplexSumOperator(const ComplexOperator &op, std::complex<double> c = 1.0);

  void AddOperator(const ComplexOperator &op, std::complex<double> c = 1.0);

  bool IsReal() const override;
  bool IsImag() const override;

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

// Wraps two operators such that: (AB)ᵀ = BᵀAᵀ and, for complex symmetric operators, the
// Hermitian transpose operation is (AB)ᴴ = BᴴAᴴ.
template <typename OperType = Operator>
class ProductOperator : public OperType
{
private:
  typedef typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                    ComplexVector, Vector>::type VecType;

  const OperType &A, &B;
  mutable VecType z;

public:
  ProductOperator(const OperType &A, const OperType &B)
    : OperType(A.Height(), B.Width()), A(A), B(B), z(B.Height())
  {
  }

  void Mult(const VecType &x, VecType &y) const override
  {
    B.Mult(x, z);
    A.Mult(z, y);
  }

  void MultTranspose(const VecType &x, VecType &y) const override
  {
    A.MultTranspose(x, z);
    B.MultTranspose(z, y);
  }

  template <typename T = OperType,
            typename = std::enable_if_t<std::is_same<T, ComplexOperator>::value>>
  void MultHermitianTranspose(const VecType &x, VecType &y) const override
  {
    A.MultHermitianTranspose(x, z);
    B.MultHermitianTranspose(z, y);
  }
};

using ComplexProductOperator = ProductOperator<ComplexOperator>;

// Applies the simple (symmetric) operator: diag(d).
template <typename OperType = Operator>
class DiagonalOperator : public OperType
{
private:
  typedef typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                    ComplexVector, Vector>::type VecType;

  const VecType &d;

public:
  DiagonalOperator(const VecType &d) : OperType(d.Size()), d(d) {}

  void Mult(const VecType &x, VecType &y) const override;

  void MultTranspose(const VecType &x, VecType &y) const override { Mult(x, y); }

  template <typename T = OperType,
            typename = std::enable_if_t<std::is_same<T, ComplexOperator>::value>>
  void MultHermitianTranspose(const VecType &x, VecType &y) const override;
};

using ComplexDiagonalOperator = DiagonalOperator<ComplexOperator>;

// A container for a sequence of operators corresponding to a multigrid hierarchy.
// Optionally includes operators for the auxiliary space at each level as well. The
// Operators are stored from coarsest to finest level. The height and width of this operator
// are never set.
template <typename OperType = Operator>
class MultigridOperator : public OperType
{
private:
  typedef typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                    ComplexVector, Vector>::type VecType;

  std::vector<std::unique_ptr<OperType>> ops, aux_ops;

public:
  MultigridOperator(int l) : OperType(0)
  {
    ops.reserve(l);
    aux_ops.reserve(l);
  }

  void AddOperator(std::unique_ptr<OperType> &&op)
  {
    ops.push_back(std::move(op));
    height = ops.back()->Height();
    width = ops.back()->Width();
  }

  void AddAuxiliaryOperator(std::unique_ptr<OperType> &&aux_op)
  {
    aux_ops.push_back(std::move(aux_op));
  }

  bool HasAuxiliaryOperators() const { return !aux_ops.empty(); }

  int GetNumLevels() const { return static_cast<int>(ops.size()); }
  int GetNumAuxiliaryLevels() const { return static_cast<int>(aux_ops.size()); }

  const OperType &GetFinestOperator() const { return *ops.back(); }
  const OperType &GetFinestAuxiliaryOperator() const { return *aux_ops.back(); }

  const OperType &GetOperatorAtLevel(int l) const { return *ops[l]; }
  const OperType &GetAuxiliaryOperatorAtLevel(int l) const { return *aux_ops[l]; }

  void Mult(const VecType &x, VecType &y) const override { GetFinestOperator().Mult(x, y); }
};

using ComplexMultigridOperator = MultigridOperator<ComplexOperator>;

namespace linalg
{

// Estimate operator 2-norm (spectral norm) using power iteration. Assumes the operator is
// not symmetric or Hermitian unless specified.
double SpectralNorm(MPI_Comm comm, const Operator &A, bool sym = false, double tol = 1.0e-4,
                    int max_it = 100);
double SpectralNorm(MPI_Comm comm, const ComplexOperator &A, bool herm = false,
                    double tol = 1.0e-4, int max_it = 100);

}  // namespace linalg

}  // namespace palace

namespace mfem
{

// A symmetric bilinear form operator which replaces *MultTranspose with *Mult.
class SymmetricBilinearForm : public BilinearForm
{
public:
  using BilinearForm::BilinearForm;

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
  void AddMultTranspose(const Vector &x, Vector &y, double c = 1.0) const override
  {
    AddMult(x, y, c);
  }
};

}  // namespace mfem

#endif  // PALACE_LINALG_OPERATOR_HPP
