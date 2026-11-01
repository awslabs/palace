// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_FUSED_VECTOR_HPP
#define PALACE_LINALG_FUSED_VECTOR_HPP

#include <cmath>
#include <complex>
#include <vector>
#include "linalg/vector.hpp"

namespace palace::linalg::fused
{

// Fused multiply-add a b + c when the hardware supports it (as compilers contract the
// vector updates of the unfused code), and unfused otherwise (std::fma is slow without
// hardware support).
inline double Fma(double a, double b, double c)
{
#if defined(FP_FAST_FMA)
  return std::fma(a, b, c);
#else
  return a * b + c;
#endif
}

//
// Fused host kernels for the Krylov solvers. Each kernel combines vector updates and local
// inner products which are otherwise separate passes over memory, while performing exactly
// the same floating-point operations as the separate operations (mfem::Vector::Add,
// ComplexVector::AXPY, Jacobi scaling, LocalDot), so results are unchanged. Callers must
// check UseHost() and otherwise fall back to the separate (device-capable) operations.
//

// Whether vector data lives on the host, so that the fused kernels can be used.
bool UseHost();

// Local inner product yᴴ x, same as LocalDot.
std::complex<double> LocalDot(const ComplexVector &x, const ComplexVector &y);

// Computes y += a x (as y.Add(a, x)) and returns the local inner product zᴴ y with the
// updated y, or yᴴ y when z is nullptr.
double AddDot(double a, const Vector &x, Vector &y, const Vector *z);
std::complex<double> AddDot(std::complex<double> a, const ComplexVector &x,
                            ComplexVector &y, const ComplexVector *z);

// Conjugate gradient update: x += a p and r -= a z. Then, if dinv is not nullptr, the
// Jacobi preconditioner is applied, z = dinv ∘ r (dinv real-valued), and the local inner
// product rᴴ z is returned (otherwise returns zero).
double CgUpdate(double a, const Vector &p, Vector &x, Vector &r, Vector &z,
                const Vector *dinv);
std::complex<double> CgUpdate(std::complex<double> a, const ComplexVector &p,
                              ComplexVector &x, ComplexVector &r, ComplexVector &z,
                              const Vector *dinv);

// For each row i in rows, y[i] = Σ_k a[k] x[j[k]] (or y[i] += Σ when add is true) over the
// entries [ia[i], ia[i + 1]) of a CSR matrix, accumulated in order from zero with rounded
// products, as hypre's (vectorized) CSR matrix-vector product on the host.
void CsrRowSums(const std::vector<int> &rows, const int *ia, const int *ja, const double *a,
                const double *x, double *y, bool add);

// Same as CsrRowSums, for two vectors at once (each with the same operations).
void CsrRowSums(const std::vector<int> &rows, const int *ia, const int *ja, const double *a,
                const double *xr, const double *xi, double *yr, double *yi, bool add);

// Same as CsrRowSums, for n vectors at once (each with the same operations), stored with
// strides ldx and ldy.
void CsrRowSums(const std::vector<int> &rows, const int *ia, const int *ja, const double *a,
                int n, const double *x, int ldx, double *y, int ldy, bool add);

}  // namespace palace::linalg::fused

#endif  // PALACE_LINALG_FUSED_VECTOR_HPP
