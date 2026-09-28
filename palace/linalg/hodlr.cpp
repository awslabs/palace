// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hodlr.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace palace
{

namespace
{

// Symmetric eigendecomposition by cyclic Jacobi (no LAPACK dependency): A (n x n,
// row-major, destroyed), eigenvectors -> V (row-major, columns are eigenvectors),
// eigenvalues -> lam.
void SymEig(std::vector<double> &A, int n, std::vector<double> &V, std::vector<double> &lam)
{
  V.assign(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; i++)
  {
    V[static_cast<std::size_t>(i) * n + i] = 1.0;
  }
  auto a = [&](int i, int j) -> double & { return A[static_cast<std::size_t>(i) * n + j]; };
  auto v = [&](int i, int j) -> double & { return V[static_cast<std::size_t>(i) * n + j]; };
  for (int sweep = 0; sweep < 100; sweep++)
  {
    double off = 0.0;
    for (int p = 0; p < n; p++)
    {
      for (int q = p + 1; q < n; q++)
      {
        off += a(p, q) * a(p, q);
      }
    }
    if (off < 1e-30)
    {
      break;
    }
    for (int p = 0; p < n; p++)
    {
      for (int q = p + 1; q < n; q++)
      {
        const double apq = a(p, q);
        if (std::abs(apq) < 1e-300)
        {
          continue;
        }
        const double phi = 0.5 * std::atan2(2.0 * apq, a(q, q) - a(p, p));
        const double c = std::cos(phi), s = std::sin(phi);
        for (int k = 0; k < n; k++)
        {
          const double kp = a(k, p), kq = a(k, q);
          a(k, p) = c * kp - s * kq;
          a(k, q) = s * kp + c * kq;
        }
        for (int k = 0; k < n; k++)
        {
          const double pk = a(p, k), qk = a(q, k);
          a(p, k) = c * pk - s * qk;
          a(q, k) = s * pk + c * qk;
        }
        for (int k = 0; k < n; k++)
        {
          const double kp = v(k, p), kq = v(k, q);
          v(k, p) = c * kp - s * kq;
          v(k, q) = s * kp + c * kq;
        }
      }
    }
  }
  lam.assign(n, 0.0);
  for (int i = 0; i < n; i++)
  {
    lam[i] = a(i, i);
  }
}

}  // namespace

long long Hodlr::Storage() const
{
  long long s = 0;
  for (const auto &lf : leaves)
  {
    s += static_cast<long long>(lf.m) * lf.m;
  }
  for (const auto &b : blocks)
  {
    s += static_cast<long long>(b.r) * (b.m1 + b.m2);
  }
  return s;
}

void Hodlr::Mult(const double *xp, double *yp) const
{
  std::fill(yp, yp + n, 0.0);
  for (const auto &lf : leaves)
  {
    for (int i = 0; i < lf.m; i++)
    {
      const double *row = &lf.D[static_cast<std::size_t>(i) * lf.m];
      double acc = 0.0;
      for (int j = 0; j < lf.m; j++)
      {
        acc += row[j] * xp[lf.s + j];
      }
      yp[lf.s + i] += acc;
    }
  }
  std::vector<double> t;
  for (const auto &b : blocks)
  {
    t.assign(b.r, 0.0);  // t = V^T x2
    for (int j = 0; j < b.m2; j++)
    {
      const double *vj = &b.V[static_cast<std::size_t>(j) * b.r];
      const double xj = xp[b.s2 + j];
      for (int k = 0; k < b.r; k++)
      {
        t[k] += vj[k] * xj;
      }
    }
    for (int i = 0; i < b.m1; i++)  // y1 += U t
    {
      const double *ui = &b.U[static_cast<std::size_t>(i) * b.r];
      double acc = 0.0;
      for (int k = 0; k < b.r; k++)
      {
        acc += ui[k] * t[k];
      }
      yp[b.s1 + i] += acc;
    }
    t.assign(b.r, 0.0);  // t = U^T x1
    for (int i = 0; i < b.m1; i++)
    {
      const double *ui = &b.U[static_cast<std::size_t>(i) * b.r];
      const double xi = xp[b.s1 + i];
      for (int k = 0; k < b.r; k++)
      {
        t[k] += ui[k] * xi;
      }
    }
    for (int j = 0; j < b.m2; j++)  // y2 += V t
    {
      const double *vj = &b.V[static_cast<std::size_t>(j) * b.r];
      double acc = 0.0;
      for (int k = 0; k < b.r; k++)
      {
        acc += vj[k] * t[k];
      }
      yp[b.s2 + j] += acc;
    }
  }
}

std::vector<double> Hodlr::Serialize() const
{
  std::vector<double> b;
  b.push_back(n);
  b.push_back(static_cast<double>(leaves.size()));
  b.push_back(static_cast<double>(blocks.size()));
  for (int p : perm)
  {
    b.push_back(p);
  }
  for (const auto &lf : leaves)
  {
    b.push_back(lf.s);
    b.push_back(lf.m);
    b.insert(b.end(), lf.D.begin(), lf.D.end());
  }
  for (const auto &bl : blocks)
  {
    b.push_back(bl.s1);
    b.push_back(bl.m1);
    b.push_back(bl.s2);
    b.push_back(bl.m2);
    b.push_back(bl.r);
    b.insert(b.end(), bl.U.begin(), bl.U.end());
    b.insert(b.end(), bl.V.begin(), bl.V.end());
  }
  return b;
}

Hodlr Hodlr::Deserialize(const std::vector<double> &b)
{
  Hodlr h;
  std::size_t p = 0;
  h.n = static_cast<int>(b[p++]);
  const int nleaf = static_cast<int>(b[p++]);
  const int nblk = static_cast<int>(b[p++]);
  h.perm.resize(h.n);
  for (int i = 0; i < h.n; i++)
  {
    h.perm[i] = static_cast<int>(b[p++]);
  }
  for (int i = 0; i < nleaf; i++)
  {
    Leaf lf;
    lf.s = static_cast<int>(b[p++]);
    lf.m = static_cast<int>(b[p++]);
    lf.D.assign(b.begin() + p, b.begin() + p + static_cast<std::size_t>(lf.m) * lf.m);
    p += static_cast<std::size_t>(lf.m) * lf.m;
    h.leaves.push_back(std::move(lf));
  }
  for (int i = 0; i < nblk; i++)
  {
    Block bl;
    bl.s1 = static_cast<int>(b[p++]);
    bl.m1 = static_cast<int>(b[p++]);
    bl.s2 = static_cast<int>(b[p++]);
    bl.m2 = static_cast<int>(b[p++]);
    bl.r = static_cast<int>(b[p++]);
    bl.U.assign(b.begin() + p, b.begin() + p + static_cast<std::size_t>(bl.m1) * bl.r);
    p += static_cast<std::size_t>(bl.m1) * bl.r;
    bl.V.assign(b.begin() + p, b.begin() + p + static_cast<std::size_t>(bl.m2) * bl.r);
    p += static_cast<std::size_t>(bl.m2) * bl.r;
    h.blocks.push_back(std::move(bl));
  }
  return h;
}

void BuildHodlr(const std::vector<double> &S, int n, const std::vector<int> &idx, int base,
                const std::vector<double> &coords, double tol, int leaf, Hodlr &h)
{
  const int m = static_cast<int>(idx.size());
  if (m <= leaf)
  {
    Hodlr::Leaf lf;
    lf.s = base;
    lf.m = m;
    lf.D.assign(static_cast<std::size_t>(m) * m, 0.0);
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < m; j++)
      {
        lf.D[static_cast<std::size_t>(i) * m + j] =
            S[static_cast<std::size_t>(idx[i]) * n + idx[j]];
      }
      h.perm[base + i] = idx[i];
    }
    h.leaves.push_back(std::move(lf));
    return;
  }
  int axis = 0;
  double best = -1.0;
  for (int d = 0; d < 3; d++)
  {
    double lo = 1e300, hi = -1e300;
    for (int i : idx)
    {
      const double x = coords[static_cast<std::size_t>(i) * 3 + d];
      lo = std::min(lo, x);
      hi = std::max(hi, x);
    }
    if (hi - lo > best)
    {
      best = hi - lo;
      axis = d;
    }
  }
  std::vector<int> s(idx);
  std::sort(s.begin(), s.end(),
            [&](int i, int j)
            {
              return coords[static_cast<std::size_t>(i) * 3 + axis] <
                     coords[static_cast<std::size_t>(j) * 3 + axis];
            });
  std::vector<int> I1(s.begin(), s.begin() + m / 2), I2(s.begin() + m / 2, s.end());
  const int m1 = static_cast<int>(I1.size()), m2 = static_cast<int>(I2.size());
  // Recurse first so the final permuted order within each child range is fixed; the
  // off-diagonal factors below are then indexed consistently with `perm` (the recursion
  // reorders I1/I2 internally).
  BuildHodlr(S, n, I1, base, coords, tol, leaf, h);
  BuildHodlr(S, n, I2, base + m1, coords, tol, leaf, h);
  // Off-diagonal block B = S[rows, cols] with rows/cols the final permuted indices.
  const int *rows = &h.perm[base], *cols = &h.perm[base + m1];
  // C = B^T B (m2 x m2); singular values of B are sqrt(eig(C)).
  std::vector<double> C(static_cast<std::size_t>(m2) * m2, 0.0);
  for (int a = 0; a < m2; a++)
  {
    for (int b = a; b < m2; b++)
    {
      double acc = 0.0;
      for (int k = 0; k < m1; k++)
      {
        acc += S[static_cast<std::size_t>(rows[k]) * n + cols[a]] *
               S[static_cast<std::size_t>(rows[k]) * n + cols[b]];
      }
      C[static_cast<std::size_t>(a) * m2 + b] = acc;
      C[static_cast<std::size_t>(b) * m2 + a] = acc;
    }
  }
  std::vector<double> Vc, lam;
  SymEig(C, m2, Vc, lam);
  double lmax = 0.0;
  for (double l : lam)
  {
    lmax = std::max(lmax, l);
  }
  std::vector<int> keep;
  for (int k = 0; k < m2; k++)
  {
    if (lam[k] > tol * tol * lmax)
    {
      keep.push_back(k);
    }
  }
  const int rB = static_cast<int>(keep.size());
  // B ~= (B V_keep) V_keep^T; store U := B V_keep (m1 x rB), V := V_keep (m2 x rB).
  Hodlr::Block bl;
  bl.s1 = base;
  bl.m1 = m1;
  bl.s2 = base + m1;
  bl.m2 = m2;
  bl.r = rB;
  bl.U.assign(static_cast<std::size_t>(m1) * rB, 0.0);
  bl.V.assign(static_cast<std::size_t>(m2) * rB, 0.0);
  for (int i = 0; i < m1; i++)
  {
    for (int kk = 0; kk < rB; kk++)
    {
      double acc = 0.0;
      for (int l = 0; l < m2; l++)
      {
        acc += S[static_cast<std::size_t>(rows[i]) * n + cols[l]] *
               Vc[static_cast<std::size_t>(l) * m2 + keep[kk]];
      }
      bl.U[static_cast<std::size_t>(i) * rB + kk] = acc;
    }
  }
  for (int j = 0; j < m2; j++)
  {
    for (int kk = 0; kk < rB; kk++)
    {
      bl.V[static_cast<std::size_t>(j) * rB + kk] =
          Vc[static_cast<std::size_t>(j) * m2 + keep[kk]];
    }
  }
  h.blocks.push_back(std::move(bl));
}

}  // namespace palace
