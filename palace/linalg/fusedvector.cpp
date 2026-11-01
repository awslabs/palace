// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fusedvector.hpp"

#include <cmath>

// This file is compiled with floating-point contraction disabled (see CMakeLists.txt):
// every fused multiply-add is explicit (Fma) and every other product is rounded, matching
// the code generated for the separate operations (mfem::Vector::Add and ComplexVector::AXPY
// contract to FMA, hypre_SeqVectorInnerProd accumulates rounded products in index order).

namespace palace::linalg::fused
{

namespace
{

// y += a x for real vectors, as mfem::Vector::Add (which skips a = 0).
inline void AddReal(double a, double x, double &y)
{
  y = Fma(a, x, y);
}

// y += a x for complex a with nonzero imaginary part, as ComplexVector::AXPY. Its
// vectorized loop and scalar remainder loop contract the products differently.
inline void AddComplexPair(double ar, double ai, double xr, double xi, double &yr,
                           double &yi)
{
  const double t = Fma(ai, xr, ar * xi);
  const double u = Fma(-ai, xi, ar * xr);
  yr = yr + u;
  yi = yi + t;
}

inline void AddComplexTail(double ar, double ai, double xr, double xi, double &yr,
                           double &yi)
{
  const double t = Fma(ai, xr, ar * xi);
  const double u = Fma(ar, xr, -(ai * xi));
  yr = yr + u;
  yi = yi + t;
}

// Number of leading elements processed by the vectorized loop of ComplexVector::AXPY (two
// elements per iteration, entered for more than two elements).
inline int PairCount(int N)
{
  return (N > 2) ? (N / 2) * 2 : 0;
}

}  // namespace

bool UseHost()
{
  return !mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
}

std::complex<double> LocalDot(const ComplexVector &x, const ComplexVector &y)
{
  const int N = x.Size();
  const double *XR = x.Real().HostRead(), *XI = x.Imag().HostRead();
  if (&x == &y)
  {
    double rr = 0.0, ii = 0.0;
    for (int i = 0; i < N; i++)
    {
      rr += XR[i] * XR[i];
      ii += XI[i] * XI[i];
    }
    return {rr + ii, 0.0};
  }
  const double *YR = y.Real().HostRead(), *YI = y.Imag().HostRead();
  double rr = 0.0, ii = 0.0, ir = 0.0, ri = 0.0;
  for (int i = 0; i < N; i++)
  {
    rr += YR[i] * XR[i];
    ii += YI[i] * XI[i];
    ir += YR[i] * XI[i];
    ri += YI[i] * XR[i];
  }
  return {rr + ii, ir - ri};
}

double AddDot(double a, const Vector &x, Vector &y, const Vector *z)
{
  const int N = y.Size();
  const double *X = x.HostRead();
  double *Y = y.HostReadWrite();
  const double *Z = z ? z->HostRead() : Y;
  double dot = 0.0;
  if (a != 0.0)
  {
    for (int i = 0; i < N; i++)
    {
      AddReal(a, X[i], Y[i]);
      dot += Z[i] * Y[i];
    }
  }
  else
  {
    for (int i = 0; i < N; i++)
    {
      dot += Z[i] * Y[i];
    }
  }
  return dot;
}

std::complex<double> AddDot(std::complex<double> a, const ComplexVector &x,
                            ComplexVector &y, const ComplexVector *z)
{
  const int N = y.Size();
  const double ar = a.real(), ai = a.imag();
  const double *XR = x.Real().HostRead(), *XI = x.Imag().HostRead();
  double *YR = y.Real().HostReadWrite(), *YI = y.Imag().HostReadWrite();
  if (!z)
  {
    // Updated yᴴ y as LocalDot(y, y): the two parts are accumulated separately.
    double rr = 0.0, ii = 0.0;
    auto Accumulate = [&](int i)
    {
      rr += YR[i] * YR[i];
      ii += YI[i] * YI[i];
    };
    if (ai == 0.0)
    {
      for (int i = 0; i < N; i++)
      {
        AddReal(ar, XR[i], YR[i]);
        AddReal(ar, XI[i], YI[i]);
        Accumulate(i);
      }
    }
    else
    {
      const int N2 = PairCount(N);
      for (int i = 0; i < N2; i++)
      {
        AddComplexPair(ar, ai, XR[i], XI[i], YR[i], YI[i]);
        Accumulate(i);
      }
      for (int i = N2; i < N; i++)
      {
        AddComplexTail(ar, ai, XR[i], XI[i], YR[i], YI[i]);
        Accumulate(i);
      }
    }
    return {rr + ii, 0.0};
  }
  const double *ZR = z->Real().HostRead(), *ZI = z->Imag().HostRead();
  double rr = 0.0, ii = 0.0, ir = 0.0, ri = 0.0;
  auto Accumulate = [&](int i)
  {
    rr += ZR[i] * YR[i];
    ii += ZI[i] * YI[i];
    ir += ZR[i] * YI[i];
    ri += ZI[i] * YR[i];
  };
  if (ai == 0.0)
  {
    for (int i = 0; i < N; i++)
    {
      AddReal(ar, XR[i], YR[i]);
      AddReal(ar, XI[i], YI[i]);
      Accumulate(i);
    }
  }
  else
  {
    const int N2 = PairCount(N);
    for (int i = 0; i < N2; i++)
    {
      AddComplexPair(ar, ai, XR[i], XI[i], YR[i], YI[i]);
      Accumulate(i);
    }
    for (int i = N2; i < N; i++)
    {
      AddComplexTail(ar, ai, XR[i], XI[i], YR[i], YI[i]);
      Accumulate(i);
    }
  }
  return {rr + ii, ir - ri};
}

double CgUpdate(double a, const Vector &p, Vector &x, Vector &r, Vector &z,
                const Vector *dinv)
{
  const int N = x.Size();
  const double *P = p.HostRead();
  double *X = x.HostReadWrite(), *R = r.HostReadWrite();
  double *Z = dinv ? z.HostReadWrite() : const_cast<double *>(z.HostRead());
  const double *D = dinv ? dinv->HostRead() : nullptr;
  if (a == 0.0)
  {
    // mfem::Vector::Add skips the update for a zero coefficient.
    double dot = 0.0;
    if (D)
    {
      for (int i = 0; i < N; i++)
      {
        Z[i] = D[i] * R[i];
        dot += R[i] * Z[i];
      }
    }
    return dot;
  }
  if (!D)
  {
    for (int i = 0; i < N; i++)
    {
      AddReal(a, P[i], X[i]);
      AddReal(-a, Z[i], R[i]);
    }
    return 0.0;
  }
  double dot = 0.0;
  for (int i = 0; i < N; i++)
  {
    AddReal(a, P[i], X[i]);
    AddReal(-a, Z[i], R[i]);
    Z[i] = D[i] * R[i];
    dot += R[i] * Z[i];
  }
  return dot;
}

std::complex<double> CgUpdate(std::complex<double> a, const ComplexVector &p,
                              ComplexVector &x, ComplexVector &r, ComplexVector &z,
                              const Vector *dinv)
{
  const int N = x.Size();
  const double ar = a.real(), ai = a.imag();
  const double *PR = p.Real().HostRead(), *PI = p.Imag().HostRead();
  double *XR = x.Real().HostReadWrite(), *XI = x.Imag().HostReadWrite();
  double *RR = r.Real().HostReadWrite(), *RI = r.Imag().HostReadWrite();
  double *ZR = dinv ? z.Real().HostReadWrite() : const_cast<double *>(z.Real().HostRead());
  double *ZI = dinv ? z.Imag().HostReadWrite() : const_cast<double *>(z.Imag().HostRead());
  const double *D = dinv ? dinv->HostRead() : nullptr;
  double rr = 0.0, ii = 0.0, ir = 0.0, ri = 0.0;
  auto Precondition = [&](int i)
  {
    ZR[i] = D[i] * RR[i];
    ZI[i] = D[i] * RI[i];
    rr += RR[i] * ZR[i];
    ii += RI[i] * ZI[i];
    ir += RR[i] * ZI[i];
    ri += RI[i] * ZR[i];
  };
  if (ai == 0.0)
  {
    for (int i = 0; i < N; i++)
    {
      AddReal(ar, PR[i], XR[i]);
      AddReal(ar, PI[i], XI[i]);
      AddReal(-ar, ZR[i], RR[i]);
      AddReal(-ar, ZI[i], RI[i]);
      if (D)
      {
        Precondition(i);
      }
    }
  }
  else
  {
    const int N2 = PairCount(N);
    for (int i = 0; i < N2; i++)
    {
      AddComplexPair(ar, ai, PR[i], PI[i], XR[i], XI[i]);
      AddComplexPair(-ar, -ai, ZR[i], ZI[i], RR[i], RI[i]);
      if (D)
      {
        Precondition(i);
      }
    }
    for (int i = N2; i < N; i++)
    {
      AddComplexTail(ar, ai, PR[i], PI[i], XR[i], XI[i]);
      AddComplexTail(-ar, -ai, ZR[i], ZI[i], RR[i], RI[i]);
      if (D)
      {
        Precondition(i);
      }
    }
  }
  return {rr + ii, ir - ri};
}

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

void CsrRowSums(const std::vector<int> &rows, const int *ia, const int *ja, const double *a,
                int n, const double *x, int ldx, double *y, int ldy, bool add)
{
  constexpr int max_n = 32;
  if (n > max_n)
  {
    for (int c = 0; c < n; c++)
    {
      CsrRowSums(rows, ia, ja, a, x + std::size_t(c) * ldx, y + std::size_t(c) * ldy, add);
    }
    return;
  }
  double t[max_n];
  for (const int i : rows)
  {
    for (int c = 0; c < n; c++)
    {
      t[c] = 0.0;
    }
    for (int k = ia[i]; k < ia[i + 1]; k++)
    {
      const double ak = a[k];
      const double *xk = x + ja[k];
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#elif defined(__GNUC__) && (__GNUC__ >= 14)
#pragma GCC novector
#endif
      for (int c = 0; c < n; c++)
      {
        t[c] += ak * xk[std::size_t(c) * ldx];
      }
    }
    for (int c = 0; c < n; c++)
    {
      double &yi = y[std::size_t(c) * ldy + i];
      yi = add ? yi + t[c] : t[c];
    }
  }
}

}  // namespace palace::linalg::fused
