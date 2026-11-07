// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hodlr.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <random>
#include <Eigen/Dense>

namespace palace
{

namespace
{

// The low-rank factors B ~= U V^T of a block with the singular values of B above
// tol * sigma_1(B) (U holding them, V orthonormal): an orthonormal basis Q of the range of
// B, grown by blocks of Gaussian probes until the residual of the probes bounds that of B
// below the truncation (with probability 1 - 10^-kProbes; Halko, Martinsson and Tropp,
// SIAM Review 53, 2011) or reaches the rounding of the probes, past which they add nothing
// to the range, and the truncated SVD of Q^T B.
void LowRankFactors(const Eigen::MatrixXd &B, double tol, Eigen::MatrixXd &U,
                    Eigen::MatrixXd &V)
{
  constexpr int kProbes = 16;
  const Eigen::Index m1 = B.rows(), m2 = B.cols(), rmax = std::min(m1, m2);
  std::mt19937 gen(1);  // the factors are computed once, on one process
  std::normal_distribution<double> normal;
  Eigen::MatrixXd Q(m1, 0);
  double floor = -1.0;  // the residual of a probe below which Q is complete
  while (Q.cols() < rmax)
  {
    const Eigen::Index k = std::min<Eigen::Index>(kProbes, rmax - Q.cols());
    // The probes are drawn once each (an expression could be evaluated more than once).
    const Eigen::MatrixXd W =
        Eigen::MatrixXd::NullaryExpr(m2, k, [&normal, &gen]() { return normal(gen); });
    Eigen::MatrixXd Y = B * W;
    const double rounding =
        100.0 * std::numeric_limits<double>::epsilon() * Y.colwise().norm().maxCoeff();
    for (int pass = 0; pass < 2; pass++)  // twice is enough
    {
      Y -= Q * (Q.transpose() * Y);
    }
    if (Y.colwise().norm().maxCoeff() <= std::max(floor, rounding))
    {
      break;
    }
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(Y);
    Q.conservativeResize(Eigen::NoChange, Q.cols() + k);
    Q.rightCols(k) = qr.householderQ() * Eigen::MatrixXd::Identity(m1, k);
    if (floor < 0.0)
    {
      // ||(I - Q Q^T) B|| <= 10 sqrt(2/pi) max_i ||(I - Q Q^T) B w_i||, and
      // sigma_1(Q^T B) <= sigma_1(B).
      const double s1 =
          Eigen::BDCSVD<Eigen::MatrixXd>(Q.transpose() * B).singularValues()(0);
      floor = tol * s1 / (10.0 * std::sqrt(2.0 / std::numbers::pi));
    }
  }
  // Q spans the range found, but its last columns may have lost orthogonality.
  Q = Eigen::HouseholderQR<Eigen::MatrixXd>(Q).householderQ() *
      Eigen::MatrixXd::Identity(m1, Q.cols());
  Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(
      Q.transpose() * B);
  const Eigen::VectorXd &s = svd.singularValues();
  Eigen::Index r = 0;
  while (r < s.size() && s(r) > tol * s(0))
  {
    r++;
  }
  U = Q * (svd.matrixU().leftCols(r) * s.head(r).asDiagonal());
  V = svd.matrixV().leftCols(r);
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
  Eigen::MatrixXd B(m1, m2), U, V;
  for (int i = 0; i < m1; i++)
  {
    for (int j = 0; j < m2; j++)
    {
      B(i, j) = S[static_cast<std::size_t>(rows[i]) * n + cols[j]];
    }
  }
  LowRankFactors(B, tol, U, V);
  Hodlr::Block bl;
  bl.s1 = base;
  bl.m1 = m1;
  bl.s2 = base + m1;
  bl.m2 = m2;
  bl.r = static_cast<int>(U.cols());
  bl.U.resize(static_cast<std::size_t>(m1) * bl.r);
  bl.V.resize(static_cast<std::size_t>(m2) * bl.r);
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
      bl.U.data(), m1, bl.r) = U;
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
      bl.V.data(), m2, bl.r) = V;
  h.blocks.push_back(std::move(bl));
}

}  // namespace palace
