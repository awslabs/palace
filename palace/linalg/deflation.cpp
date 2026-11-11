// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "deflation.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <mfem.hpp>
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Inner products h[i] = W[i]ᴴ x (one pass over the vectors, one global reduction).
void BlockDot(MPI_Comm comm, const std::vector<ComplexVector> &W, const ComplexVector &x,
              std::vector<std::complex<double>> &h)
{
  h.resize(W.size());
  if (W.empty())
  {
    return;
  }
  linalg::LocalMultiDot(W, static_cast<int>(W.size()), x, nullptr, h.data(), nullptr);
  Mpi::GlobalSum(static_cast<int>(h.size()), h.data(), comm);
}

// y += Σᵢ c[i] W[i] for i < m (one pass over the vectors).
void BlockAdd(const std::vector<ComplexVector> &W, int m, const std::complex<double> *c,
              ComplexVector &y)
{
  std::vector<std::complex<double>> mc(m);
  for (int i = 0; i < m; i++)
  {
    mc[i] = -c[i];
  }
  linalg::MultiSubtract(W, m, y, mc.data(), 1.0);
}

// Single precision basis kernels (host): each vector stores its n real parts, then its n
// imaginary parts. The vectors are streamed once per call, by chunks of the double
// precision operand, which stays in cache.
constexpr int chunk = 1024;

// Local inner products h[j] = W[j]ᴴ x for j < k.
void LocalDotSingle(const std::vector<std::vector<float>> &W, int k, int n,
                    const ComplexVector &x, std::complex<double> *h)
{
  const double *XR = x.Real().HostRead(), *XI = x.Imag().HostRead();
  std::vector<double> sr(k, 0.0), si(k, 0.0);
  for (int i0 = 0; i0 < n; i0 += chunk)
  {
    const int i1 = std::min(n, i0 + chunk);
    for (int j = 0; j < k; j++)
    {
      const float *WR = W[j].data(), *WI = WR + n;
      double ar = 0.0, ai = 0.0;
      for (int i = i0; i < i1; i++)
      {
        const double wr = WR[i], wi = WI[i];
        ar += wr * XR[i] + wi * XI[i];
        ai += wr * XI[i] - wi * XR[i];
      }
      sr[j] += ar;
      si[j] += ai;
    }
  }
  for (int j = 0; j < k; j++)
  {
    h[j] = {sr[j], si[j]};
  }
}

// y += s Σⱼ c[j] W[j] for j < k.
void AddSingle(const std::vector<std::vector<float>> &W, int k, int n,
               const std::complex<double> *c, double s, ComplexVector &y)
{
  double *YR = y.Real().HostReadWrite(), *YI = y.Imag().HostReadWrite();
  for (int i0 = 0; i0 < n; i0 += chunk)
  {
    const int i1 = std::min(n, i0 + chunk);
    for (int j = 0; j < k; j++)
    {
      const float *WR = W[j].data(), *WI = WR + n;
      const double cr = s * c[j].real(), ci = s * c[j].imag();
      for (int i = i0; i < i1; i++)
      {
        const double wr = WR[i], wi = WI[i];
        YR[i] += cr * wr - ci * wi;
        YI[i] += cr * wi + ci * wr;
      }
    }
  }
}

// W[j] = Σₗ V[l] Q(l, j) for j < k, rounded to single precision: one pass over V, by
// chunks of rows accumulated in double precision.
void CombineSingle(const std::vector<ComplexVector> &V, int m, const Eigen::MatrixXcd &Q,
                   int k, int n, std::vector<std::vector<float>> &W)
{
  W.assign(k, std::vector<float>(2 * static_cast<std::size_t>(n)));
  std::vector<const double *> VR(m), VI(m);
  for (int l = 0; l < m; l++)
  {
    VR[l] = V[l].Real().HostRead();
    VI[l] = V[l].Imag().HostRead();
  }
  std::vector<double> acc(2 * static_cast<std::size_t>(k) * chunk);
  for (int i0 = 0; i0 < n; i0 += chunk)
  {
    const int len = std::min(n, i0 + chunk) - i0;
    std::fill(acc.begin(), acc.end(), 0.0);
    for (int l = 0; l < m; l++)
    {
      const double *vr = VR[l] + i0, *vi = VI[l] + i0;
      for (int j = 0; j < k; j++)
      {
        const double qr = Q(l, j).real(), qi = Q(l, j).imag();
        double *ar = acc.data() + 2 * static_cast<std::size_t>(j) * chunk, *ai = ar + chunk;
        for (int i = 0; i < len; i++)
        {
          ar[i] += qr * vr[i] - qi * vi[i];
          ai[i] += qr * vi[i] + qi * vr[i];
        }
      }
    }
    for (int j = 0; j < k; j++)
    {
      const double *ar = acc.data() + 2 * static_cast<std::size_t>(j) * chunk,
                   *ai = ar + chunk;
      float *wr = W[j].data() + i0, *wi = W[j].data() + n + i0;
      for (int i = 0; i < len; i++)
      {
        wr[i] = static_cast<float>(ar[i]);
        wi[i] = static_cast<float>(ai[i]);
      }
    }
  }
}

// Local Gram matrix E(i, j) = W[i]ᴴ X[j] for i, j < k.
void LocalGramSingle(const std::vector<std::vector<float>> &W,
                     const std::vector<std::vector<float>> &X, int k, int n,
                     Eigen::MatrixXcd &E)
{
  E.setZero(k, k);
  for (int i0 = 0; i0 < n; i0 += chunk)
  {
    const int i1 = std::min(n, i0 + chunk);
    for (int j = 0; j < k; j++)
    {
      const float *XR = X[j].data(), *XI = XR + n;
      for (int l = 0; l < k; l++)
      {
        const float *WR = W[l].data(), *WI = WR + n;
        double sr = 0.0, si = 0.0;
        for (int i = i0; i < i1; i++)
        {
          const double wr = WR[i], wi = WI[i], xr = XR[i], xi = XI[i];
          sr += wr * xr + wi * xi;
          si += wr * xi - wi * xr;
        }
        E(l, j) += std::complex<double>(sr, si);
      }
    }
  }
}

std::vector<float> ToSingle(const ComplexVector &x)
{
  const int n = x.Size();
  std::vector<float> v(2 * static_cast<std::size_t>(n));
  const double *XR = x.Real().HostRead(), *XI = x.Imag().HostRead();
  for (int i = 0; i < n; i++)
  {
    v[i] = static_cast<float>(XR[i]);
    v[n + i] = static_cast<float>(XI[i]);
  }
  return v;
}

void FromSingle(const std::vector<float> &v, ComplexVector &x)
{
  const int n = x.Size();
  double *XR = x.Real().HostWrite(), *XI = x.Imag().HostWrite();
  for (int i = 0; i < n; i++)
  {
    XR[i] = v[i];
    XI[i] = v[n + i];
  }
}

bool UseSinglePrecision()
{
  return !mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
}

// Iterations of GMRES for the (small, dense) matrix M started from the first unit vector to
// reach the relative residual tol, at most maxit (maxit + 1 if not reached).
int ModelGmresIterations(const Eigen::MatrixXcd &M, double tol, int maxit)
{
  const int m = static_cast<int>(M.rows());
  maxit = std::min(maxit, m);
  Eigen::MatrixXcd V = Eigen::MatrixXcd::Zero(m, maxit + 1);
  V(0, 0) = 1.0;
  std::vector<std::complex<double>> cs(maxit), sn(maxit);
  std::complex<double> g = 1.0;
  for (int j = 0; j < maxit; j++)
  {
    Eigen::VectorXcd w = M * V.col(j);
    Eigen::VectorXcd h = Eigen::VectorXcd::Zero(j + 2);
    for (int pass = 0; pass < 2; pass++)
    {
      for (int i = 0; i <= j; i++)
      {
        const std::complex<double> d = V.col(i).dot(w);
        h(i) += d;
        w -= d * V.col(i);
      }
    }
    h(j + 1) = w.norm();
    for (int i = 0; i < j; i++)
    {
      const std::complex<double> a = std::conj(cs[i]) * h(i) + std::conj(sn[i]) * h(i + 1);
      h(i + 1) = -sn[i] * h(i) + cs[i] * h(i + 1);
      h(i) = a;
    }
    const double d = std::sqrt(std::norm(h(j)) + std::norm(h(j + 1)));
    if (d == 0.0)
    {
      return j + 1;
    }
    cs[j] = h(j) / d;
    sn[j] = h(j + 1) / d;
    g = -sn[j] * g;
    if (std::abs(g) <= tol)
    {
      return j + 1;
    }
    if (std::abs(h(j + 1)) == 0.0)
    {
      return j + 1;
    }
    V.col(j + 1) = w / std::abs(h(j + 1));
  }
  return maxit + 1;
}

}  // namespace

void DeflationAnalysis::Compute(MPI_Comm comm, const std::complex<double> *H, int ld,
                                int m_, int kmax, const std::vector<int> &k_cand,
                                double rel_tol, double theta_max)
{
  m = m_;
  kmax = std::min(kmax, m);
  // The candidates are chosen on the root once the harmonic Ritz values are known (at most
  // the given ones plus the number of outliers).
  constexpr int max_cand = 16;
  std::vector<int> cand(max_cand, 0), its(max_cand + 1, m + 1);
  int nk = 0;
  num_outliers = 0;
  theta.assign(kmax, 0.0);
  Y = Eigen::MatrixXcd::Zero(m, kmax);
  if (kmax > 0 && Mpi::Root(comm))
  {
    // Harmonic Ritz pairs (θ, y) of B A with respect to zero: the eigenpairs of
    // G = Hₘ + |h|² Hₘ⁻ᴴ eₘ eₘᵀ, with h = H(m, m - 1), whose vectors V y have residuals
    // orthogonal to B A V.
    Eigen::MatrixXcd Hm(m, m);
    for (int j = 0; j < m; j++)
    {
      for (int i = 0; i < m; i++)
      {
        Hm(i, j) = H[static_cast<std::size_t>(j) * ld + i];
      }
    }
    const std::complex<double> h = H[static_cast<std::size_t>(m - 1) * ld + m];
    Eigen::VectorXcd em = Eigen::VectorXcd::Zero(m);
    em(m - 1) = 1.0;
    Eigen::MatrixXcd G = Hm;
    const Eigen::VectorXcd f = Hm.adjoint().fullPivLu().solve(em);
    if (f.allFinite())
    {
      G += std::norm(h) * f * em.transpose();
    }
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(G, true);
    const Eigen::VectorXcd &lambda = es.eigenvalues();
    const Eigen::MatrixXcd &S = es.eigenvectors();
    std::vector<int> idx(m);
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&](int a, int b)
                     { return std::abs(lambda(a)) < std::abs(lambda(b)); });
    for (int j = 0; j < kmax; j++)
    {
      theta[j] = lambda(idx[j]);
      Y.col(j) = S.col(idx[j]);
    }
    if (!Y.allFinite())
    {
      Y.setZero();
    }
    int limit = kmax;
    if (theta_max > 0.0)
    {
      num_outliers = 0;
      while (num_outliers < kmax && std::abs(theta[num_outliers]) < theta_max)
      {
        num_outliers++;
      }
      limit = num_outliers;
    }
    for (int k : k_cand)
    {
      if (k > 0 && k <= limit && nk < max_cand - 1)
      {
        cand[nk++] = k;
      }
    }
    if (theta_max > 0.0 && limit > 0 && (nk == 0 || cand[nk - 1] != limit))
    {
      cand[nk++] = limit;
    }
    if (theta_max > 0.0 && limit > 0)
    {
      // The larger sets, of the values below 2 and 4 times the bound.
      for (double f : {2.0, 4.0})
      {
        int count = limit;
        while (count < kmax && std::abs(theta[count]) < f * theta_max)
        {
          count++;
        }
        if (count > cand[nk - 1] && nk < max_cand)
        {
          cand[nk++] = count;
        }
      }
    }

    // Model iteration counts: GMRES for G, and for S Θ' S⁻¹ with the k values closest to
    // zero of Θ set to one, when the eigenvector basis reproduces G (otherwise the
    // deflated counts stay unknown, m + 1).
    its[0] = ModelGmresIterations(G, rel_tol, m);
    Eigen::PartialPivLU<Eigen::MatrixXcd> slu(S);
    const Eigen::MatrixXcd Sinv = slu.inverse();
    const double rec = (S * lambda.asDiagonal() * Sinv - G).norm() / G.norm();
    if (Sinv.allFinite() && rec < 1.0e-6)
    {
      for (int c = 0; c < nk; c++)
      {
        Eigen::VectorXcd l = lambda;
        for (int j = 0; j < cand[c]; j++)
        {
          l(idx[j]) = 1.0;
        }
        const Eigen::MatrixXcd Gk = S * l.asDiagonal() * Sinv;
        its[c + 1] = ModelGmresIterations(Gk, rel_tol, m);
      }
    }
  }
  MPI_Bcast(Y.data(), 2 * m * kmax, MPI_DOUBLE, 0, comm);
  MPI_Bcast(theta.data(), 2 * kmax, MPI_DOUBLE, 0, comm);
  int head[2] = {nk, num_outliers};
  MPI_Bcast(head, 2, MPI_INT, 0, comm);
  nk = head[0];
  num_outliers = head[1];
  MPI_Bcast(cand.data(), max_cand, MPI_INT, 0, comm);
  MPI_Bcast(its.data(), max_cand + 1, MPI_INT, 0, comm);
  its_plain = its[0];
  k_candidates.assign(cand.begin(), cand.begin() + nk);
  its_deflated.assign(its.begin() + 1, its.begin() + 1 + nk);
}

int DeflationAnalysis::Count(double bound) const
{
  int count = 0;
  while (count < static_cast<int>(theta.size()) && std::abs(theta[count]) < bound)
  {
    count++;
  }
  return count;
}

DeflatedPreconditioner::DeflatedPreconditioner(MPI_Comm comm)
  : Solver<ComplexOperator>(), comm(comm), single(UseSinglePrecision())
{
}

int DeflatedPreconditioner::Build(const ComplexOperator &A_,
                                  const Solver<ComplexOperator> &B_,
                                  const std::vector<ComplexVector> &V,
                                  const DeflationAnalysis &analysis, int k)
{
  Clear();
  A = &A_;
  B = &B_;
  const int m = analysis.m;
  k = std::min(k, static_cast<int>(analysis.Y.cols()));
  if (k <= 0 || m <= 0)
  {
    return 0;
  }

  // Orthonormal coefficients Q of the harmonic Ritz vectors in the basis V (two passes of
  // modified Gram-Schmidt on the root, without the vectors which are numerically dependent
  // on the previous ones): W = V Q is orthonormal as V is.
  Eigen::MatrixXcd Q(m, k);
  prefix.assign(k + 1, 0);
  if (Mpi::Root(comm))
  {
    int kk = 0;
    for (int j = 0; j < k; j++)
    {
      Eigen::VectorXcd q = analysis.Y.col(j);
      const double norm0 = q.norm();
      for (int pass = 0; pass < 2; pass++)
      {
        for (int i = 0; i < kk; i++)
        {
          q -= Q.col(i).dot(q) * Q.col(i);
        }
      }
      const double norm = q.norm();
      if (norm > 1.0e-8 * norm0 && norm > 0.0)
      {
        Q.col(kk++) = q / norm;
      }
      prefix[j + 1] = kk;
    }
  }
  MPI_Bcast(prefix.data(), k + 1, MPI_INT, 0, comm);
  const int kk = prefix[k];
  if (kk == 0)
  {
    Clear();
    return 0;
  }
  MPI_Bcast(Q.data(), 2 * m * kk, MPI_DOUBLE, 0, comm);

  // W = V Q, on the host in one pass over V and rounded to single precision.
  n = V[0].Size();
  if (single)
  {
    CombineSingle(V, m, Q, kk, n, W32);
  }
  else
  {
    for (int j = 0; j < kk; j++)
    {
      ComplexVector w(n);
      w.UseDevice(true);
      w = 0.0;
      BlockAdd(V, m, Q.col(j).data(), w);
      W.push_back(std::move(w));
    }
  }
  t.SetSize(n);
  t.UseDevice(true);
  return kk;
}

int DeflatedPreconditioner::Activate(int k)
{
  k = std::clamp(k, 0, Reserve());
  const int done = single ? static_cast<int>(AW32.size()) : static_cast<int>(AW.size());
  if (k > done)
  {
    // Images of the new vectors (on the host, of the rounded vectors), and E for all the
    // vectors activated so far.
    for (int j = done; j < k; j++)
    {
      ComplexVector aw(n);
      aw.UseDevice(true);
      if (single)
      {
        ComplexVector w(n);
        w.UseDevice(true);
        FromSingle(W32[j], w);
        A->Mult(w, aw);
        AW32.push_back(ToSingle(aw));
      }
      else
      {
        A->Mult(W[j], aw);
        AW.push_back(std::move(aw));
      }
    }
    if (single)
    {
      LocalGramSingle(W32, AW32, k, n, E);
    }
    else
    {
      E.resize(k, k);
      std::vector<std::complex<double>> h(k);
      for (int j = 0; j < k; j++)
      {
        linalg::LocalMultiDot(W, k, AW[j], nullptr, h.data(), nullptr);
        for (int i = 0; i < k; i++)
        {
          E(i, j) = h[i];
        }
      }
    }
    Mpi::GlobalSum(k * k, E.data(), comm);
  }
  if (k == 0)
  {
    active = 0;
    return 0;
  }
  const Eigen::MatrixXcd Ek = E.topLeftCorner(k, k);
  Eigen::FullPivLU<Eigen::MatrixXcd> check(Ek);
  if (Ek.allFinite() && check.isInvertible())
  {
    lu.compute(Ek);
    active = k;
  }
  return active;
}

void DeflatedPreconditioner::Clear()
{
  W.clear();
  AW.clear();
  W32.clear();
  AW32.clear();
  prefix.clear();
  active = 0;
  E.resize(0, 0);
}

std::size_t DeflatedPreconditioner::GetBytes() const
{
  std::size_t bytes = 0;
  for (const auto *X : {&W, &AW})
  {
    for (const auto &v : *X)
    {
      bytes += 2 * static_cast<std::size_t>(v.Size()) * sizeof(double);
    }
  }
  for (const auto *X : {&W32, &AW32})
  {
    for (const auto &v : *X)
    {
      bytes += v.size() * sizeof(float);
    }
  }
  return bytes;
}

void DeflatedPreconditioner::Coefficients(const ComplexVector &x) const
{
  const int k = Size();
  c.resize(k);
  if (single)
  {
    LocalDotSingle(W32, k, n, x, c.data());
  }
  else
  {
    linalg::LocalMultiDot(W, k, x, nullptr, c.data(), nullptr);
  }
  Mpi::GlobalSum(k, c.data(), comm);
  Eigen::Map<Eigen::VectorXcd> cv(c.data(), k);
  cv = lu.solve(Eigen::VectorXcd(cv));
}

void DeflatedPreconditioner::AddBasis(const std::complex<double> *coef,
                                      ComplexVector &y) const
{
  if (single)
  {
    AddSingle(W32, Size(), n, coef, 1.0, y);
  }
  else
  {
    BlockAdd(W, Size(), coef, y);
  }
}

void DeflatedPreconditioner::Mult(const ComplexVector &x, ComplexVector &y) const
{
  MFEM_ASSERT(B, "Missing preconditioner for DeflatedPreconditioner!");
  const int k = Size();
  if (k == 0)
  {
    B->Mult(x, y);
    return;
  }
  Coefficients(x);
  t = x;
  if (single)
  {
    AddSingle(AW32, k, n, c.data(), -1.0, t);
  }
  else
  {
    linalg::MultiSubtract(AW, k, t, c.data(), 1.0);
  }
  B->Mult(t, y);
  AddBasis(c.data(), y);
}

std::pair<double, double> MeasureDeflationCost(MPI_Comm comm,
                                               const std::vector<ComplexVector> &V, int m)
{
  // Enough vectors not to fit in cache, as the deflation vectors do not.
  const int k = std::min(32, m);
  if (k <= 0)
  {
    return {0.0, 0.0};
  }
  const int n = V[0].Size();
  ComplexVector x(n), y(n);
  x.UseDevice(true);
  y.UseDevice(true);
  x = V[0];
  y = 0.0;
  std::vector<std::complex<double>> h(k, 0.0);
  double best = mfem::infinity(), best_reduce = mfem::infinity();
  if (UseSinglePrecision())
  {
    std::vector<std::vector<float>> W(k);
    for (int j = 0; j < k; j++)
    {
      W[j] = ToSingle(V[j]);
    }
    for (int rep = 0; rep < 3; rep++)
    {
      const double t0 = MPI_Wtime();
      LocalDotSingle(W, k, n, x, h.data());
      AddSingle(W, k, n, h.data(), -1.0e-30, y);
      AddSingle(W, k, n, h.data(), 1.0e-30, y);
      best = std::min(best, MPI_Wtime() - t0);
    }
  }
  else
  {
    std::vector<ComplexVector> W(V.begin(), V.begin() + k);
    for (int rep = 0; rep < 3; rep++)
    {
      const double t0 = MPI_Wtime();
      linalg::LocalMultiDot(W, k, x, nullptr, h.data(), nullptr);
      BlockAdd(W, k, h.data(), y);
      BlockAdd(W, k, h.data(), y);
      best = std::min(best, MPI_Wtime() - t0);
    }
  }
  std::vector<std::complex<double>> r(32, 0.0);
  for (int rep = 0; rep < 3; rep++)
  {
    MPI_Barrier(comm);
    const double t0 = MPI_Wtime();
    Mpi::GlobalSum(static_cast<int>(r.size()), r.data(), comm);
    best_reduce = std::min(best_reduce, MPI_Wtime() - t0);
  }
  double t[2] = {best / k, best_reduce};
  Mpi::GlobalMax(2, t, comm);
  return {t[0], t[1]};
}

void DeflationLevels::Reset(const std::vector<int> &sizes_)
{
  sizes.clear();
  for (int k : sizes_)
  {
    if (k > 0 && (sizes.empty() || k > sizes.back()))
    {
      sizes.push_back(k);
    }
  }
  level = 0;
  top = static_cast<int>(sizes.size());
  trying = false;
}

DeflationLevels::Action DeflationLevels::Grow()
{
  if (level + 1 < top)
  {
    level++;
    trying = true;
    return Action::GROW;
  }
  trying = false;
  return Action::KEEP;
}

DeflationLevels::Action DeflationLevels::Verify(bool gain)
{
  trying = false;
  if (gain)
  {
    return Grow();
  }
  if (level > 0)
  {
    top = level--;
    return Action::BACK;
  }
  return Action::STOP;
}

DeflationLevels::Action DeflationLevels::Compare(bool gain)
{
  if (!trying || level == 0)
  {
    return Verify(gain);
  }
  if (gain)
  {
    return Grow();
  }
  trying = false;
  top = level--;
  return Action::BACK;
}

}  // namespace palace
