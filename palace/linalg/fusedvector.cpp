// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fusedvector.hpp"

// This file is compiled with floating-point contraction disabled (see CMakeLists.txt):
// every product is rounded, as in hypre's host CSR matrix-vector product.

namespace palace::linalg::fused
{

void CsrRowSums(const std::vector<int> &rows, const int *ia, const int *ja, const double *a,
                const double *x, double *y, bool add)
{
  // The rows have few entries: a scalar loop is faster than a vectorized in-order reduction
  // (same result).
  for (const int i : rows)
  {
    double t = 0.0;
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#elif defined(__GNUC__) && (__GNUC__ >= 14)
#pragma GCC novector
#endif
    for (int k = ia[i]; k < ia[i + 1]; k++)
    {
      t += a[k] * x[ja[k]];
    }
    y[i] = add ? y[i] + t : t;
  }
}

void CsrRowSums(const std::vector<int> &rows, const int *ia, const int *ja, const double *a,
                const double *xr, const double *xi, double *yr, double *yi, bool add)
{
  for (const int i : rows)
  {
    double tr = 0.0, ti = 0.0;
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#elif defined(__GNUC__) && (__GNUC__ >= 14)
#pragma GCC novector
#endif
    for (int k = ia[i]; k < ia[i + 1]; k++)
    {
      tr += a[k] * xr[ja[k]];
      ti += a[k] * xi[ja[k]];
    }
    yr[i] = add ? yr[i] + tr : tr;
    yi[i] = add ? yi[i] + ti : ti;
  }
}

void CsrRowSums(const std::vector<int> &rows, const std::vector<int> &out, const int *ia,
                const int *ja, const double *a, const double *x, double *y, bool add)
{
  for (std::size_t r = 0; r < rows.size(); r++)
  {
    const int i = rows[r];
    double t = 0.0;
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#elif defined(__GNUC__) && (__GNUC__ >= 14)
#pragma GCC novector
#endif
    for (int k = ia[i]; k < ia[i + 1]; k++)
    {
      t += a[k] * x[ja[k]];
    }
    y[out[r]] = add ? y[out[r]] + t : t;
  }
}

void CsrRowSums(const std::vector<int> &rows, const std::vector<int> &out, const int *ia,
                const int *ja, const double *a, const double *xr, const double *xi,
                double *yr, double *yi, bool add)
{
  for (std::size_t r = 0; r < rows.size(); r++)
  {
    const int i = rows[r];
    double tr = 0.0, ti = 0.0;
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#elif defined(__GNUC__) && (__GNUC__ >= 14)
#pragma GCC novector
#endif
    for (int k = ia[i]; k < ia[i + 1]; k++)
    {
      tr += a[k] * xr[ja[k]];
      ti += a[k] * xi[ja[k]];
    }
    yr[out[r]] = add ? yr[out[r]] + tr : tr;
    yi[out[r]] = add ? yi[out[r]] + ti : ti;
  }
}

namespace
{

template <typename Values>
void CompactRowSums(const std::vector<int> &rows, const std::vector<int> &ia,
                    const std::vector<int> &ja, const Values &a, const double *x, double *y)
{
  for (std::size_t r = 0; r < rows.size(); r++)
  {
    double t = 0.0;
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#elif defined(__GNUC__) && (__GNUC__ >= 14)
#pragma GCC novector
#endif
    for (int k = ia[r]; k < ia[r + 1]; k++)
    {
      t += a(k) * x[ja[k]];
    }
    y[rows[r]] = t;
  }
}

template <typename Values>
void CompactRowSums(const std::vector<int> &rows, const std::vector<int> &ia,
                    const std::vector<int> &ja, const Values &a, const double *xr,
                    const double *xi, double *yr, double *yi)
{
  for (std::size_t r = 0; r < rows.size(); r++)
  {
    double tr = 0.0, ti = 0.0;
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#elif defined(__GNUC__) && (__GNUC__ >= 14)
#pragma GCC novector
#endif
    for (int k = ia[r]; k < ia[r + 1]; k++)
    {
      const double ak = a(k);
      tr += ak * xr[ja[k]];
      ti += ak * xi[ja[k]];
    }
    yr[rows[r]] = tr;
    yi[rows[r]] = ti;
  }
}

}  // namespace

void CompactCsrRowSums(const std::vector<int> &rows, const std::vector<int> &ia,
                       const std::vector<int> &ja, const double *a, const double *x,
                       double *y)
{
  CompactRowSums(rows, ia, ja, [a](int k) { return a[k]; }, x, y);
}

void CompactCsrRowSums(const std::vector<int> &rows, const std::vector<int> &ia,
                       const std::vector<int> &ja, const double *a, const double *xr,
                       const double *xi, double *yr, double *yi)
{
  CompactRowSums(rows, ia, ja, [a](int k) { return a[k]; }, xr, xi, yr, yi);
}

void CompactCsrRowSums(const std::vector<int> &rows, const std::vector<int> &ia,
                       const std::vector<int> &ja, const std::uint16_t *v,
                       const double *table, const double *x, double *y)
{
  CompactRowSums(rows, ia, ja, [v, table](int k) { return table[v[k]]; }, x, y);
}

void CompactCsrRowSums(const std::vector<int> &rows, const std::vector<int> &ia,
                       const std::vector<int> &ja, const std::uint16_t *v,
                       const double *table, const double *xr, const double *xi, double *yr,
                       double *yi)
{
  CompactRowSums(rows, ia, ja, [v, table](int k) { return table[v[k]]; }, xr, xi, yr, yi);
}

}  // namespace palace::linalg::fused
