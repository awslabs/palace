// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_AAA_HPP
#define PALACE_UTILS_AAA_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>

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

// Polynomial-plus-pole–residue form of the AAA rational interpolant:
//
//   x = (z − polynomial_center) / polynomial_scale,
//   r(z) = Σₖ polynomial(k) xᵏ + Σⱼ residues(j) / (z − poles(j)).
//
// A constant asymptote is the common proper-rational case, but a barycentric denominator
// can lose one or more leading coefficients. The resulting interpolant is improper and
// has a nonconstant polynomial quotient, which must be retained when converting to partial
// fractions. Keeping that quotient in the normalized coordinate prevents catastrophic
// cancellation when the support points have a large common translation. Used to build the
// augmented L⁻¹/R⁻¹/C state-space realisation: polynomial terms are folded into the base
// pencil, while each pole becomes an auxiliary state.
struct AAAPoleResidue
{
  Eigen::VectorXcd poles;                        // finite poles in physical z
  Eigen::VectorXcd residues;                     // one residue per finite pole
  Eigen::VectorXcd polynomial;                   // quotient in x, low-degree first
  std::complex<double> polynomial_center = 0.0;  // affine origin for x
  double polynomial_scale = 1.0;                 // positive affine scale for x
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
#if defined(MFEM_DEBUG)
  // Duplicate sample points produce 1/0 = NaN inside the Loewner build; fail loudly.
  for (std::size_t i = 0; i < M; i++)
  {
    for (std::size_t j = i + 1; j < M; j++)
    {
      MFEM_VERIFY(z(static_cast<long>(i)) != z(static_cast<long>(j)),
                  "AAA sample points must be distinct!");
    }
  }
#endif
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
  r.err_history.resize(static_cast<long>(m_max));
  long n_hist = 0;
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
      r.err_history(n_hist++) = 0.0;
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
    r.err_history(n_hist++) = cur_err;
    if (cur_err <= abs_tol)
    {
      r.converged = true;
      break;
    }
  }
  r.err_history.conservativeResize(n_hist);

  const auto m_used = static_cast<long>(support_indices.size());
  r.zj = zj_buf.head(m_used);
  r.fj = fj_buf.head(m_used);
  r.wj = wj_buf;
  return r;
}

// Extract the polynomial quotient, poles, and residues from an AAA result.
//
// Poles = roots of the denominator D(z) = Σ wⱼ / (z − zⱼ): solve the (m+1)×(m+1)
// generalised eigenproblem (Nakatsukasa et al. 2018, sec 2.7)
//
//        ⎡ 0   wᵀ      ⎤        ⎡ 0   0    ⎤
//   A = ⎢            ⎥,   B = ⎢          ⎥
//        ⎣ 1  diag(xⱼ) ⎦        ⎣ 0   I    ⎦
//
// in the affine-normalized coordinate x = (z-center)/scale, then map finite eigenvalues
// back to z. The normalization is essential when the support points have a large common
// translation. Discard the eigenvalues at infinity (B is rank-deficient by at least 1).
// Residues use rₖ = N(pₖ) / D'(pₖ), with the polynomial quotient evaluated at pₖ
// subtracted from every fⱼ. This leaves the numerator unchanged at a denominator root but
// keeps the proper remainder numerically well-conditioned.
inline AAAPoleResidue AAAToPoleResidue(const AAAResult &r)
{
  AAAPoleResidue out;
  const auto m = r.zj.size();
  if (m == 0)
  {
    return out;
  }
  // Clear the barycentric denominators in the dimensionless coordinate
  // x = (z-center)/scale, where the support points are O(1):
  //
  //   N(x) = Σⱼ wⱼ (fⱼ-f_center) ∏ₗ₍ₗ≠ⱼ₎ (x-xₗ),
  //   D(x) = Σⱼ wⱼ              ∏ₗ₍ₗ≠ⱼ₎ (x-xₗ).
  //
  // Centering both z and f avoids cancellation from a large translation. Their leading
  // coefficients are Σw(f-f_center) and Σw, respectively. In particular, Σw ≈ 0 lowers
  // deg(D); it does NOT make the interpolant strictly proper. Coefficients below are stored
  // low-degree first.
  const std::complex<double> center = r.zj.mean();
  const std::complex<double> f_center = r.fj.mean();
  double support_scale = 0.0;
  Eigen::VectorXcd xj(m);
  for (long j = 0; j < m; j++)
  {
    support_scale = std::max(support_scale, std::abs(r.zj(j) - center));
  }
  if (support_scale == 0.0)
  {
    support_scale = 1.0;
  }
  for (long j = 0; j < m; j++)
  {
    xj(j) = (r.zj(j) - center) / support_scale;
  }
  out.polynomial_center = center;
  out.polynomial_scale = support_scale;

  std::vector<std::complex<double>> numerator(static_cast<std::size_t>(m), 0.0);
  std::vector<std::complex<double>> denominator(static_cast<std::size_t>(m), 0.0);
  for (long j = 0; j < m; j++)
  {
    std::vector<std::complex<double>> factor{1.0};
    for (long l = 0; l < m; l++)
    {
      if (l == j)
      {
        continue;
      }
      std::vector<std::complex<double>> next(factor.size() + 1, 0.0);
      for (std::size_t k = 0; k < factor.size(); k++)
      {
        next[k] -= xj(l) * factor[k];
        next[k + 1] += factor[k];
      }
      factor = std::move(next);
    }
    for (std::size_t k = 0; k < factor.size(); k++)
    {
      numerator[k] += r.wj(j) * (r.fj(j) - f_center) * factor[k];
      denominator[k] += r.wj(j) * factor[k];
    }
  }
  // In powers of dimensionless x, coefficients within each polynomial have comparable
  // units, so one relative trimming criterion remains meaningful across each vector.
  auto trim_leading_zeros = [m](std::vector<std::complex<double>> &p)
  {
    double scale = 0.0;
    for (const auto &c : p)
    {
      scale = std::max(scale, std::abs(c));
    }
    if (scale == 0.0)
    {
      p.clear();
      return;
    }
    const double tol =
        64.0 * std::numeric_limits<double>::epsilon() * std::max<long>(m, 1) * scale;
    while (!p.empty() && std::abs(p.back()) <= tol)
    {
      p.pop_back();
    }
  };
  trim_leading_zeros(numerator);
  trim_leading_zeros(denominator);
  MFEM_VERIFY(!denominator.empty(), "AAA interpolant has an identically zero denominator!");

  std::vector<std::complex<double>> quotient;
  if (numerator.size() >= denominator.size())
  {
    auto remainder = numerator;
    quotient.assign(numerator.size() - denominator.size() + 1, 0.0);
    const long degree_den = static_cast<long>(denominator.size()) - 1;
    for (long k = static_cast<long>(quotient.size()) - 1; k >= 0; k--)
    {
      const auto q = remainder[static_cast<std::size_t>(degree_den + k)] /
                     denominator[static_cast<std::size_t>(degree_den)];
      quotient[static_cast<std::size_t>(k)] = q;
      for (long j = 0; j <= degree_den; j++)
      {
        remainder[static_cast<std::size_t>(j + k)] -=
            q * denominator[static_cast<std::size_t>(j)];
      }
    }
    trim_leading_zeros(quotient);
  }

  // Restore the function centering while coefficients still share the same units. Use the
  // scale from before the addition so a proper rational whose shifted quotient is
  // -f_center is recognized as having no polynomial quotient despite roundoff.
  double quotient_scale = std::abs(f_center);
  for (const auto &c : quotient)
  {
    quotient_scale = std::max(quotient_scale, std::abs(c));
  }
  if (quotient.empty())
  {
    quotient.push_back(f_center);
  }
  else
  {
    quotient[0] += f_center;
  }
  const double quotient_tol =
      64.0 * std::numeric_limits<double>::epsilon() * std::max<long>(m, 1) * quotient_scale;
  while (!quotient.empty() && std::abs(quotient.back()) <= quotient_tol)
  {
    quotient.pop_back();
  }

  // Retain q in the normalized coordinate. Expanding q((z-center)/scale) into powers of
  // physical z would immediately reintroduce the large, mutually cancelling coefficients
  // that the affine normalization above was chosen to avoid.
  out.polynomial.resize(static_cast<long>(quotient.size()));
  for (long k = 0; k < out.polynomial.size(); k++)
  {
    out.polynomial(k) = quotient[static_cast<std::size_t>(k)];
  }

  // Solve for poles in the same normalized coordinate used to determine the polynomial
  // degrees above. Forming this pencil with the physical support points can turn exact
  // eigenvalues at infinity into bogus finite poles when all z_j share a large translation.
  Eigen::MatrixXcd A_geig = Eigen::MatrixXcd::Zero(m + 1, m + 1);
  Eigen::MatrixXcd B_geig = Eigen::MatrixXcd::Zero(m + 1, m + 1);
  A_geig.row(0).tail(m) = r.wj.transpose();
  for (long j = 0; j < m; j++)
  {
    A_geig(1 + j, 0) = 1.0;
    A_geig(1 + j, 1 + j) = xj(j);
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
  MFEM_VERIFY(info == 0, "zggev failed with info = " << info);
  // Collect the finite eigenvalues. Their count is known structurally: the poles are the
  // roots of the cleared denominator, so exactly deg(D) of the m+1 pencil eigenvalues are
  // finite and the rest lie at infinity. ZGGEV usually signals infinity with beta = 0
  // exactly, but QZ roundoff (BLAS/architecture dependent) can instead produce a small
  // nonzero beta — a bogus "finite" pole at ~1/eps — and no fixed |beta| threshold is
  // robust when the denominator degree drops (many eigenvalues at infinity, e.g. a purely
  // polynomial interpolant). Instead keep exactly deg(D) eigenvalues, selecting those
  // farthest from infinity in the chordal metric |beta| / sqrt(|alpha|² + |beta|²).
  const long n_finite = static_cast<long>(denominator.size()) - 1;
  std::vector<int> order(static_cast<std::size_t>(n));
  for (int j = 0; j < n; j++)
  {
    order[static_cast<std::size_t>(j)] = j;
  }
  auto chordal = [&](int j)
  {
    const double a = std::abs(alpha[static_cast<std::size_t>(j)]);
    const double b = std::abs(beta[static_cast<std::size_t>(j)]);
    return b / std::hypot(a, b);
  };
  std::sort(order.begin(), order.end(),
            [&](int a, int b) { return chordal(a) > chordal(b); });
  std::vector<std::complex<double>> finite_x_poles;
  std::vector<std::complex<double>> finite_poles;
  finite_x_poles.reserve(static_cast<std::size_t>(n_finite));
  finite_poles.reserve(static_cast<std::size_t>(n_finite));
  for (long k = 0; k < n_finite; k++)
  {
    const int j = order[static_cast<std::size_t>(k)];
    MFEM_VERIFY(std::abs(beta[static_cast<std::size_t>(j)]) > 0.0,
                "AAA pole extraction found fewer finite eigenvalues than the denominator "
                "degree implies!");
    const auto x_pole =
        alpha[static_cast<std::size_t>(j)] / beta[static_cast<std::size_t>(j)];
    const auto pole = center + support_scale * x_pole;
    MFEM_VERIFY(std::isfinite(pole.real()) && std::isfinite(pole.imag()),
                "AAA pole extraction produced a non-finite pole!");
    finite_x_poles.push_back(x_pole);
    finite_poles.push_back(pole);
  }
  out.poles = Eigen::Map<Eigen::VectorXcd>(finite_poles.data(), finite_poles.size());

  // Residues of the strictly proper remainder. Subtracting the polynomial quotient at the
  // pole leaves the residue unchanged because the barycentric denominator vanishes there,
  // while avoiding cancellation against a potentially large polynomial component. Evaluate
  // both the quotient and barycentric sums at the normalized pole xₖ. Converting the ratio
  // from x back to z contributes one factor of support_scale to the physical-z residue.
  out.residues.resize(out.poles.size());
  for (long k = 0; k < out.poles.size(); k++)
  {
    const std::complex<double> x_pole = finite_x_poles[static_cast<std::size_t>(k)];
    std::complex<double> polynomial_at_pole = 0.0;
    for (long j = out.polynomial.size(); j-- > 0;)
    {
      polynomial_at_pole = polynomial_at_pole * x_pole + out.polynomial(j);
    }
    std::complex<double> num = 0.0, den = 0.0;
    for (long j = 0; j < m; j++)
    {
      auto inv = 1.0 / (x_pole - xj(j));
      num += r.wj(j) * (r.fj(j) - polynomial_at_pole) * inv;
      den -= r.wj(j) * inv * inv;
    }
    out.residues(k) = support_scale * num / den;
  }
  return out;
}

// Evaluate the polynomial-plus-pole–residue partial-fraction form. Used in tests.
inline std::complex<double> EvaluatePoleResidue(const AAAPoleResidue &pr,
                                                std::complex<double> z)
{
  const std::complex<double> x = (z - pr.polynomial_center) / pr.polynomial_scale;
  std::complex<double> sum = 0.0;
  for (long k = pr.polynomial.size(); k-- > 0;)
  {
    sum = sum * x + pr.polynomial(k);
  }
  for (long k = 0; k < pr.poles.size(); k++)
  {
    sum += pr.residues(k) / (z - pr.poles(k));
  }
  return sum;
}

}  // namespace palace::utils

#endif  // PALACE_UTILS_AAA_HPP
