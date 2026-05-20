// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_AAA_HPP
#define PALACE_UTILS_AAA_HPP

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>
#include <Eigen/Dense>

// LAPACK ZGGEV — used to extract poles from the AAA generalised eigenproblem.
// Declared at file scope (not inside a function or namespace) so it is plain C
// linkage; the symbol is provided by the LAPACK library on the link line.
extern "C" void zggev_(char *, char *, int *, std::complex<double> *, int *,
                       std::complex<double> *, int *, std::complex<double> *,
                       std::complex<double> *, std::complex<double> *, int *,
                       std::complex<double> *, int *, std::complex<double> *, int *,
                       double *, int *);

namespace palace::utils
{

// AAA rational approximation (Nakatsukasa, Sète, Trefethen, 2018). Given samples
// {(zₖ, Fₖ)}, the algorithm greedily selects support points and weights so that
// the barycentric form
//
//   r(z) = Σⱼ wⱼ Fⱼ / (z − zⱼ)  /  Σⱼ wⱼ / (z − zⱼ)
//
// matches F at the support points exactly and minimises the max residual on the
// non-support samples. Iteration stops when the residual is below `tol·‖F‖∞`
// or `m_max` support points have been chosen.
//
// Use case in Palace: rational fit of the wave-port modal propagation residual
// δkₙ(ω) = kₙ(ω) − (α₀ + α₁ω + α₂ω²) for circuit synthesis, augmenting L⁻¹/R⁻¹/C
// with one aux state per support point so the synthesised matrices remain
// quadratic in ω while reproducing the exact dispersion (cf. prom-waveport.md).
struct AAAResult
{
  Eigen::VectorXcd zj;          // support points (selected from samples)
  Eigen::VectorXcd fj;          // function values at support points (= F at zj)
  Eigen::VectorXcd wj;          // barycentric weights (real-only-on-real-axis when
                                // F is real-on-real-axis, but typed complex)
  Eigen::VectorXd err_history;  // max-residual history at each iteration
  bool converged = false;
};

// Pole–residue–asymptote form of the AAA rational interpolant:
//
//   r(z) = d + Σⱼ residues(j) / (z − poles(j))
//
// This is the standard partial-fraction expansion. With d = r(∞) explicitly
// separated, residues are well-conditioned (in contrast to absorbing d into the
// numerator, which makes residues blow up). Used to build the augmented L⁻¹/R⁻¹/C
// state-space realisation: each pole becomes one auxiliary state per port,
// residues become the cross-coupling magnitudes between aux state and port mode.
struct AAAPoleResidue
{
  Eigen::VectorXcd poles;            // length m
  Eigen::VectorXcd residues;         // length m
  std::complex<double> d{0.0, 0.0};  // asymptotic value r(∞)
};

// Evaluate the AAA barycentric interpolant at a (set of) point(s). Numerically
// stable form: at z exactly at a support point, return the corresponding fj.
inline std::complex<double> EvaluateAAA(const AAAResult &r, std::complex<double> z)
{
  const auto m = static_cast<std::size_t>(r.zj.size());
  if (m == 0)
  {
    return 0.0;
  }
  std::complex<double> num = 0.0, den = 0.0;
  for (std::size_t j = 0; j < m; j++)
  {
    auto diff = z - r.zj(j);
    if (std::abs(diff) < 1.0e-14 * std::max(std::abs(r.zj(j)), 1.0))
    {
      return r.fj(j);
    }
    auto inv = 1.0 / diff;
    num += r.wj(j) * r.fj(j) * inv;
    den += r.wj(j) * inv;
  }
  return num / den;
}

// Run AAA. The samples z must be distinct; F has the same size. The algorithm
// is real-on-real symmetric for real-on-real-axis F (weights and values are real
// up to floating point round-off) but the storage type is std::complex<double>
// so the same routine handles general complex-valued F (e.g. evanescent kₙ).
inline AAAResult RunAAA(const Eigen::VectorXcd &z, const Eigen::VectorXcd &F, double tol,
                        std::size_t m_max)
{
  AAAResult r;
  const auto M = static_cast<std::size_t>(z.size());
  if (M == 0)
  {
    return r;
  }
  const double F_norm = F.cwiseAbs().maxCoeff();
  if (F_norm == 0.0)
  {
    // F is zero — return an empty interpolant; EvaluateAAA returns 0.
    r.converged = true;
    return r;
  }
  const double abs_tol = tol * F_norm;

  std::vector<bool> is_support(M, false);
  // Greedy iteration.
  Eigen::VectorXcd R = Eigen::VectorXcd::Constant(M, F.mean());
  std::vector<std::size_t> support_indices;
  support_indices.reserve(m_max);
  Eigen::VectorXcd zj_buf(m_max), fj_buf(m_max);
  Eigen::VectorXcd wj_buf;
  r.err_history.resize(0);
  for (std::size_t m = 0; m < m_max; m++)
  {
    // Pick the index in {non-support} that maximises |F − R|.
    std::size_t j_pick = 0;
    double max_err = -1.0;
    for (std::size_t i = 0; i < M; i++)
    {
      if (is_support[i])
      {
        continue;
      }
      double e = std::abs(F(i) - R(i));
      if (e > max_err)
      {
        max_err = e;
        j_pick = i;
      }
    }
    is_support[j_pick] = true;
    support_indices.push_back(j_pick);
    zj_buf(static_cast<long>(m)) = z(j_pick);
    fj_buf(static_cast<long>(m)) = F(j_pick);

    // Build the Loewner matrix on the non-support rows.
    const std::size_t M_J = M - support_indices.size();
    if (M_J == 0)
    {
      // All samples are now support points; trivial perfect fit.
      wj_buf = Eigen::VectorXcd::Ones(m + 1);
      r.err_history.conservativeResize(static_cast<long>(m + 1));
      r.err_history(static_cast<long>(m)) = 0.0;
      break;
    }
    Eigen::MatrixXcd A(M_J, m + 1);
    {
      std::size_t row = 0;
      for (std::size_t i = 0; i < M; i++)
      {
        if (is_support[i])
        {
          continue;
        }
        for (std::size_t k = 0; k <= m; k++)
        {
          auto inv = 1.0 / (z(i) - zj_buf(static_cast<long>(k)));
          A(row, k) = (F(i) - fj_buf(static_cast<long>(k))) * inv;
        }
        row++;
      }
    }
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A, Eigen::ComputeFullV);
    wj_buf = svd.matrixV().col(static_cast<long>(m));

    // Update R on all samples (support points: r = F by construction).
    double cur_err = 0.0;
    for (std::size_t i = 0; i < M; i++)
    {
      if (is_support[i])
      {
        R(i) = F(i);
        continue;
      }
      std::complex<double> num = 0.0, den = 0.0;
      for (std::size_t k = 0; k <= m; k++)
      {
        auto inv = 1.0 / (z(i) - zj_buf(static_cast<long>(k)));
        num += wj_buf(static_cast<long>(k)) * fj_buf(static_cast<long>(k)) * inv;
        den += wj_buf(static_cast<long>(k)) * inv;
      }
      R(i) = num / den;
      cur_err = std::max(cur_err, std::abs(F(i) - R(i)));
    }
    r.err_history.conservativeResize(static_cast<long>(m + 1));
    r.err_history(static_cast<long>(m)) = cur_err;
    if (cur_err <= abs_tol)
    {
      r.converged = true;
      break;
    }
  }

  const auto m_used = static_cast<long>(support_indices.size());
  r.zj = zj_buf.head(m_used);
  r.fj = fj_buf.head(m_used);
  r.wj = wj_buf;
  return r;
}

// Extract poles, residues, and asymptote from an AAA result.
//
// Asymptote d = r(∞) = (Σ wⱼ fⱼ) / (Σ wⱼ).
// Poles = roots of the denominator D(z) = Σ wⱼ / (z − zⱼ): solve the (m+1)×(m+1)
// generalised eigenproblem (Nakatsukasa et al. 2018, sec 2.7)
//
//        ⎡ 0   wᵀ      ⎤        ⎡ 0   0    ⎤
//   A = ⎢            ⎥,   B = ⎢          ⎥
//        ⎣ 1  diag(zⱼ) ⎦        ⎣ 0   I    ⎦
//
// and discard the spurious eigenvalue at infinity (B is rank-deficient by 1).
// Residues by formula rₖ = N(pₖ) / D'(pₖ) with the d-term subtracted:
//   rₖ = (Σⱼ wⱼ (fⱼ − d) / (pₖ − zⱼ)) · (1 / (− Σⱼ wⱼ / (pₖ − zⱼ)²))
// (the "− d" subtraction is what keeps residues numerically well-conditioned;
// see Nakatsukasa eq 2.21 with the rational shifted to be strictly proper).
inline AAAPoleResidue AAAToPoleResidue(const AAAResult &r)
{
  AAAPoleResidue out;
  const auto m = r.zj.size();
  if (m == 0)
  {
    return out;
  }
  std::complex<double> sum_w = r.wj.sum();
  std::complex<double> sum_wf = (r.wj.array() * r.fj.array()).sum();
  out.d = sum_wf / sum_w;

  // Generalised eigenproblem for poles.
  Eigen::MatrixXcd A_geig = Eigen::MatrixXcd::Zero(m + 1, m + 1);
  Eigen::MatrixXcd B_geig = Eigen::MatrixXcd::Zero(m + 1, m + 1);
  A_geig.row(0).tail(m) = r.wj.transpose();
  for (long j = 0; j < m; j++)
  {
    A_geig(1 + j, 0) = 1.0;
    A_geig(1 + j, 1 + j) = r.zj(j);
    B_geig(1 + j, 1 + j) = 1.0;
  }
  char jobvl = 'N', jobvr = 'N';
  int n = static_cast<int>(m + 1);
  int lwork = 4 * n;
  std::vector<std::complex<double>> alpha(n), beta(n), work(lwork);
  std::vector<double> rwork(8 * n);
  int info = 0;
  zggev_(&jobvl, &jobvr, &n, A_geig.data(), &n, B_geig.data(), &n, alpha.data(),
         beta.data(), nullptr, &n, nullptr, &n, work.data(), &lwork, rwork.data(), &info);
  // Collect finite eigenvalues (drop the single eigenvalue with beta = 0,
  // corresponding to the rank deficiency in B).
  std::vector<std::complex<double>> finite_poles;
  finite_poles.reserve(m);
  for (int j = 0; j < n; j++)
  {
    if (std::abs(beta[j]) > 0.0)
    {
      auto val = alpha[j] / beta[j];
      if (std::isfinite(val.real()) && std::isfinite(val.imag()))
      {
        finite_poles.push_back(val);
      }
    }
  }
  out.poles = Eigen::Map<Eigen::VectorXcd>(finite_poles.data(), finite_poles.size());

  // Residues with d separated out.
  out.residues.resize(out.poles.size());
  for (long k = 0; k < out.poles.size(); k++)
  {
    std::complex<double> pk = out.poles(k);
    std::complex<double> num = 0.0, den = 0.0;
    for (long j = 0; j < m; j++)
    {
      auto inv = 1.0 / (pk - r.zj(j));
      num += r.wj(j) * (r.fj(j) - out.d) * inv;
      den -= r.wj(j) * inv * inv;
    }
    out.residues(k) = num / den;
  }
  return out;
}

// Evaluate the pole–residue partial-fraction form. Used in tests.
inline std::complex<double> EvaluatePoleResidue(const AAAPoleResidue &pr,
                                                std::complex<double> z)
{
  std::complex<double> sum = pr.d;
  for (long k = 0; k < pr.poles.size(); k++)
  {
    sum += pr.residues(k) / (z - pr.poles(k));
  }
  return sum;
}

}  // namespace palace::utils

#endif  // PALACE_UTILS_AAA_HPP
