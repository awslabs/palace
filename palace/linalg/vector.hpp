// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_VECTOR_HPP
#define PALACE_LINALG_VECTOR_HPP

#include <complex>
#include <mfem.hpp>
#include "utils/communication.hpp"

namespace palace
{

using Vector = mfem::Vector;

//
// Functionality extending mfem::Vector from MFEM, including basic functions for parallel
// vectors distributed across MPI processes.
//

// A complex-valued vector represented as two real vectors, one for each component.
class ComplexVector
{
private:
  Vector xr_, xi_;

public:
  // Create a vector with the given size.
  ComplexVector(int n = 0);

  // Copy constructor.
  ComplexVector(const ComplexVector &y);

  // Copy constructor from separately provided real and imaginary parts.
  ComplexVector(const Vector &yr, const Vector &yi);

  // Copy constructor from an array of complex values.
  ComplexVector(const std::complex<double> *py, int n, bool on_dev);

  // Flag for runtime execution on the mfem::Device. See the documentation for mfem::Vector.
  void UseDevice(bool use_dev)
  {
    xr_.UseDevice(use_dev);
    xi_.UseDevice(use_dev);
  }
  bool UseDevice() const { return xr_.UseDevice(); }

  // Return the size of the vector.
  int Size() const { return xr_.Size(); }

  // Set the size of the vector. See the notes for Vector::SetSize for behavior in the cases
  // where n is less than or greater than Size() or Capacity().
  void SetSize(int n);

  // Get access to the real and imaginary vector parts.
  const Vector &Real() const { return xr_; }
  Vector &Real() { return xr_; }
  const Vector &Imag() const { return xi_; }
  Vector &Imag() { return xi_; }

  // Set from a ComplexVector, without resizing.
  ComplexVector &operator=(const ComplexVector &y) { return Set(y); }
  ComplexVector &Set(const ComplexVector &y)
  {
    Set(y.Real(), y.Imag());
    return *this;
  }

  // Set from separately provided real and imaginary parts, without resizing.
  void Set(const Vector &yr, const Vector &yi);

  // Set from an array of complex values, without resizing.
  void Set(const std::complex<double> *py, int n, bool on_dev);

  // Copy the vector into an array of complex values.
  void Get(std::complex<double> *py, int n, bool on_dev) const;

  // Set all entries equal to s.
  ComplexVector &operator=(std::complex<double> s);
  ComplexVector &operator=(double s)
  {
    *this = std::complex<double>(s, 0.0);
    return *this;
  }

  // Scale all entries by s.
  ComplexVector &operator*=(std::complex<double> s);

  // Replace entries with their complex conjugate.
  void Conj();

  // Replace entries with their absolute value.
  void Abs();

  // Set all entries to their reciprocal.
  void Reciprocal();

  // Vector dot product (yᴴ x) or indefinite dot product (yᵀ x) for complex vectors.
  std::complex<double> Dot(const ComplexVector &y) const;
  std::complex<double> TransposeDot(const ComplexVector &y) const;
  std::complex<double> operator*(const ComplexVector &y) const { return Dot(y); }

  // In-place addition (*this) += alpha * x.
  void AXPY(std::complex<double> alpha, const ComplexVector &x);
  void Add(std::complex<double> alpha, const ComplexVector &x) { AXPY(alpha, x); }
  ComplexVector &operator+=(const ComplexVector &x)
  {
    AXPY(1.0, x);
    return *this;
  }

  // In-place addition (*this) = alpha * x + beta * (*this).
  void AXPBY(std::complex<double> alpha, const ComplexVector &x, std::complex<double> beta);

  // In-place addition (*this) = alpha * x + beta * y + gamma * (*this).
  void AXPBYPCZ(std::complex<double> alpha, const ComplexVector &x,
                std::complex<double> beta, const ComplexVector &y,
                std::complex<double> gamma);

  static void AXPY(std::complex<double> alpha, const Vector &xr, const Vector &xi,
                   Vector &yr, Vector &yi);

  static void AXPBY(std::complex<double> alpha, const Vector &xr, const Vector &xi,
                    std::complex<double> beta, Vector &yr, Vector &yi);

  static void AXPBYPCZ(std::complex<double> alpha, const Vector &xr, const Vector &xi,
                       std::complex<double> beta, const Vector &yr, const Vector &yi,
                       std::complex<double> gamma, Vector &zr, Vector &zi);
};

namespace linalg
{

// Returns the global vector size.
template <typename VecType>
inline HYPRE_BigInt GlobalSize(MPI_Comm comm, const VecType &x)
{
  HYPRE_BigInt N = x.Size();
  Mpi::GlobalSum(1, &N, comm);
  return N;
}

// Returns the global vector size for two vectors.
template <typename VecType1, typename VecType2>
inline std::pair<HYPRE_BigInt, HYPRE_BigInt> GlobalSize2(MPI_Comm comm, const VecType1 &x1,
                                                         const VecType2 &x2)
{
  HYPRE_BigInt N[2] = {x1.Size(), x2.Size()};
  Mpi::GlobalSum(2, N, comm);
  return {N[0], N[1]};
}

// Sets all entries of the vector corresponding to the given indices to the given (real)
// value or vector of values.
template <typename VecType>
void SetSubVector(VecType &x, const mfem::Array<int> &rows, double s);
template <typename VecType>
void SetSubVector(VecType &x, const mfem::Array<int> &rows, const VecType &y);

// Sets all entries in the range [start, end) to  the given value.
template <typename VecType>
void SetSubVector(VecType &x, int start, int end, double s);

// Sets all entries of the vector to random numbers sampled from the [-1, 1] or [-1 - 1i,
// 1 + 1i] for complex-valued vectors.
template <typename VecType>
void SetRandom(MPI_Comm comm, VecType &x, int seed = 0);
template <typename VecType>
void SetRandomReal(MPI_Comm comm, VecType &x, int seed = 0);
template <typename VecType>
void SetRandomSign(MPI_Comm comm, VecType &x, int seed = 0);

// Calculate the inner product yᴴ x or yᵀ x.
template <typename VecType>
inline auto Dot(MPI_Comm comm, const VecType &x, const VecType &y)
{
  auto dot = x * y;
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

// Calculate the vector 2-norm.
template <typename VecType>
inline double Norml2(MPI_Comm comm, const VecType &x)
{
  return std::sqrt(std::abs(Dot(comm, x, x)));
}

// Normalize the vector, possibly with respect to an SPD matrix B.
template <typename VecType>
inline double Normalize(MPI_Comm comm, VecType &x)
{
  double norm = Norml2(comm, x);
  MFEM_ASSERT(norm > 0.0, "Zero vector norm in normalization!");
  x *= 1.0 / norm;
  return norm;
}

// Addition y += alpha * x.
template <typename VecType, typename ScalarType>
void AXPY(ScalarType alpha, const VecType &x, VecType &y);

// Addition y = alpha * x + beta * y.
template <typename VecType, typename ScalarType>
void AXPBY(ScalarType alpha, const VecType &x, ScalarType beta, VecType &y);

// Addition z = alpha * x + beta * y + gamma * z.
template <typename VecType, typename ScalarType>
void AXPBYPCZ(ScalarType alpha, const VecType &x, ScalarType beta, const VecType &y,
              ScalarType gamma, VecType &z);

}  // namespace linalg

}  // namespace palace

#endif  // PALACE_LINALG_VECTOR_HPP
