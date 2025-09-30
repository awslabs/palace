// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_VECTOR_HPP
#define PALACE_LINALG_VECTOR_HPP

#include <complex>
#include <vector>
#include <mfem.hpp>
#include "utils/communication.hpp"

namespace palace
{

//
// Functionality extending mfem::Vector from MFEM, including basic functions for parallel
// vectors distributed across MPI processes.
//

using Vector = mfem::Vector;

// A complex-valued vector represented as two real vectors, one for each component.
class ComplexVector
{
private:
  Vector xr, xi;

public:
  // Create a vector with the given size.
  ComplexVector(int size = 0);

  // Copy constructor.
  ComplexVector(const ComplexVector &y);

  // Copy constructor from separately provided real and imaginary parts.
  ComplexVector(const Vector &yr, const Vector &yi);

  // Copy constructor from an array of complex values.
  ComplexVector(const std::complex<double> *py, int size, bool on_dev);

  // Create a vector referencing the memory of another vector, at the given base offset and
  // size.
  ComplexVector(Vector &y, int offset, int size);

  // Flag for runtime execution on the mfem::Device. See the documentation for mfem::Vector.
  void UseDevice(bool use_dev);
  bool UseDevice() const { return xr.UseDevice(); }

  // Return the size of the vector.
  int Size() const { return xr.Size(); }

  // Set the size of the vector. See the notes for Vector::SetSize for behavior in the cases
  // where the new size is less than or greater than Size() or Capacity().
  void SetSize(int size);

  // Set this vector to reference the memory of another vector, at the given base offset and
  // size.
  void MakeRef(Vector &y, int offset, int size);

  // Get access to the real and imaginary vector parts.
  const Vector &Real() const { return xr; }
  Vector &Real() { return xr; }
  const Vector &Imag() const { return xi; }
  Vector &Imag() { return xi; }

  // Set from a ComplexVector, without resizing.
  void Set(const ComplexVector &y);
  ComplexVector &operator=(const ComplexVector &y)
  {
    Set(y);
    return *this;
  }

  // Set from separately provided real and imaginary parts, without resizing.
  void Set(const Vector &yr, const Vector &yi);

  // Set from an array of complex values, without resizing.
  void Set(const std::complex<double> *py, int size, bool on_dev);

  // Copy the vector into an array of complex values.
  void Get(std::complex<double> *py, int size, bool on_dev) const;

  // Set all entries equal to s.
  ComplexVector &operator=(std::complex<double> s);
  ComplexVector &operator=(double s)
  {
    *this = std::complex<double>(s, 0.0);
    return *this;
  }

  // Set the vector from an array of blocks and coefficients, without resizing.
  void SetBlocks(const std::vector<const ComplexVector *> &y,
                 const std::vector<std::complex<double>> &s);

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
  void Subtract(std::complex<double> alpha, const ComplexVector &x) { AXPY(-alpha, x); }
  ComplexVector &operator+=(const ComplexVector &x)
  {
    AXPY(1.0, x);
    return *this;
  }
  ComplexVector &operator-=(const ComplexVector &x)
  {
    AXPY(-1.0, x);
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

// A stack-allocated vector with compile-time fixed size.
//
// StaticVector provides a Vector interface backed by stack memory instead of
// heap allocation. The size N is fixed at compile time, making it suitable for
// small vectors where performance and avoiding dynamic allocation are
// important.
//
// Template parameters:
// - N: The fixed size of the vector (number of elements)
//
// Notes:
// - Inherits from mfem::Vector, so can be used anywhere Vector is expected.
// - Memory is automatically managed (no new/delete needed).
// - Faster than dynamic Vector for small sizes due to stack allocation.
//
// Example usage:
//
// StaticVector<3> vec;  // 3D vector on stack
// vec[0] = 1.0;
// vec[1] = 2.0;
// vec[2] = 3.0;
//
// vec.Sum();
//
// You can also create StaticComplexVectors:
//
// StaticVector<3> vec_real, vec_imag;
// ComplexVector complex_vec(vec_real, vec_imag);
template <int N>
class StaticVector : public Vector
{
private:
  double buff[N];

public:
  StaticVector() : Vector() { SetDataAndSize(buff, N); }

  ~StaticVector()
  {
    MFEM_ASSERT(GetData() == buff,
                "Buffer of StaticVector changed. This indicates a possible bug.");
    MFEM_ASSERT(Size() == N, "Size of StaticVector changed. This indicates a possible bug.")
  }

  using Vector::operator=;  // Extend the implicitly defined assignment operators
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

// Sets contiguous entries from start to the given vector.
template <typename VecType>
void SetSubVector(VecType &x, int start, const VecType &y);

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

// Calculate the local inner product yᴴ x or yᵀ x.
double LocalDot(const Vector &x, const Vector &y);
std::complex<double> LocalDot(const ComplexVector &x, const ComplexVector &y);

// Calculate the parallel inner product yᴴ x or yᵀ x.
template <typename VecType>
inline auto Dot(MPI_Comm comm, const VecType &x, const VecType &y)
{
  auto dot = LocalDot(x, y);
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

// Calculate the vector 2-norm.
template <typename VecType>
inline auto Norml2(MPI_Comm comm, const VecType &x)
{
  return std::sqrt(std::abs(Dot(comm, x, x)));
}

// Normalize the vector, possibly with respect to an SPD matrix B.
template <typename VecType>
inline auto Normalize(MPI_Comm comm, VecType &x)
{
  auto norm = Norml2(comm, x);
  MFEM_ASSERT(norm > 0.0, "Zero vector norm in normalization!");
  x *= 1.0 / norm;
  return norm;
}

// Calculate the local sum of all elements in the vector.
double LocalSum(const Vector &x);
std::complex<double> LocalSum(const ComplexVector &x);

// Calculate the sum of all elements in the vector.
template <typename VecType>
inline auto Sum(MPI_Comm comm, const VecType &x)
{
  auto sum = LocalSum(x);
  Mpi::GlobalSum(1, &sum, comm);
  return sum;
}

// Calculate the mean of all elements in the vector.
template <typename VecType>
inline auto Mean(MPI_Comm comm, const VecType &x)
{
  using ScalarType = typename std::conditional<std::is_same<VecType, ComplexVector>::value,
                                               std::complex<double>, double>::type;
  ScalarType sum[2] = {LocalSum(x), ScalarType(x.Size())};
  Mpi::GlobalSum(2, sum, comm);
  return sum[0] / sum[1];
}

// Normalize a complex vector so its mean is on the positive real axis.
// Returns the original mean phase.
inline double NormalizePhase(MPI_Comm comm, ComplexVector &x)
{
  std::complex<double> mean = Mean(comm, x);
  x *= std::conj(mean) / std::abs(mean);
  return std::atan2(mean.imag(), mean.real());
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

// Compute element-wise square root, optionally with scaling (multiplied before the square
// root).
void Sqrt(Vector &x, double s = 1.0);

// Compute the 3D Cartesian product between A and B and store the result in C.
// If add is true, accumulate the result to C instead of overwriting its
// content.
template <typename VecTypeA, typename VecTypeB, typename VecTypeC>
void Cross3(const VecTypeA &A, const VecTypeB &B, VecTypeC &C, bool add = false)
{
  if (add)
  {
    C[0] += A[1] * B[2] - A[2] * B[1];
    C[1] += A[2] * B[0] - A[0] * B[2];
    C[2] += A[0] * B[1] - A[1] * B[0];
  }
  else
  {
    C[0] = A[1] * B[2] - A[2] * B[1];
    C[1] = A[2] * B[0] - A[0] * B[2];
    C[2] = A[0] * B[1] - A[1] * B[0];
  }
}

}  // namespace linalg

}  // namespace palace

#endif  // PALACE_LINALG_VECTOR_HPP
