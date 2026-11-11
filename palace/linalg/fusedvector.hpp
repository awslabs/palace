// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_FUSED_VECTOR_HPP
#define PALACE_LINALG_FUSED_VECTOR_HPP

#include <cmath>
#include <cstdint>
#include <vector>

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

// For each row i in rows, y[i] = Σ_k a[k] x[j[k]] (or y[i] += Σ when add is true) over the
// entries [ia[i], ia[i + 1]) of a CSR matrix, accumulated in order from zero with rounded
// products, as hypre's (vectorized) CSR matrix-vector product on the host.
void CsrRowSums(const std::vector<int> &rows, const int *ia, const int *ja, const double *a,
                const double *x, double *y, bool add);

// Same as CsrRowSums, for two vectors at once (each with the same operations).
void CsrRowSums(const std::vector<int> &rows, const int *ia, const int *ja, const double *a,
                const double *xr, const double *xi, double *yr, double *yi, bool add);

// Same as CsrRowSums, with the result of row rows[r] in y[out[r]].
void CsrRowSums(const std::vector<int> &rows, const std::vector<int> &out, const int *ia,
                const int *ja, const double *a, const double *x, double *y, bool add);
void CsrRowSums(const std::vector<int> &rows, const std::vector<int> &out, const int *ia,
                const int *ja, const double *a, const double *xr, const double *xi,
                double *yr, double *yi, bool add);

// Same as CsrRowSums (with add false) for a compact CSR block whose row r, with the entries
// [ia[r], ia[r + 1]), is row rows[r] of y, with the values a[k] or table[v[k]].
void CompactCsrRowSums(const std::vector<int> &rows, const std::vector<int> &ia,
                       const std::vector<int> &ja, const double *a, const double *x,
                       double *y);
void CompactCsrRowSums(const std::vector<int> &rows, const std::vector<int> &ia,
                       const std::vector<int> &ja, const double *a, const double *xr,
                       const double *xi, double *yr, double *yi);
void CompactCsrRowSums(const std::vector<int> &rows, const std::vector<int> &ia,
                       const std::vector<int> &ja, const std::uint16_t *v,
                       const double *table, const double *x, double *y);
void CompactCsrRowSums(const std::vector<int> &rows, const std::vector<int> &ia,
                       const std::vector<int> &ja, const std::uint16_t *v,
                       const double *table, const double *xr, const double *xi, double *yr,
                       double *yi);

}  // namespace palace::linalg::fused

#endif  // PALACE_LINALG_FUSED_VECTOR_HPP
