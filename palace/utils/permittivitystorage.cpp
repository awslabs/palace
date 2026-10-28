// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Available storage of the RAW (PoleResidue) block of a material, see PreparedStorage in
// permittivitymodel.hpp.
//
// Scaled variables: u = s/Ω, λ = p/Ω, ρ = r/(κΩ), with Ω and κ powers of two so that the
// scaling is exact. With q = λρ, the polarization current is
//
//   Y(u) = u χ_raw(Ωu)/κ = D + Σⱼ qⱼ/(u − λⱼ) + Σ_lossless 2q u/(u² + β²),
//
// D = Σ ρ over every pole (exactly, from the input residues), including the zero pole,
// whose integrator is unobservable in the current and is not a state. Exact lossless poles
// (Re p = 0) are oscillators with a fixed energy, and do not enter F(u) = Y(u) + Y(−u). For
// the stable part, with z = u², βⱼ = λⱼ²,
//
//   F(u) = 2ρ₀ + z Σⱼ 2ρⱼ/(z − βⱼ) = 2D + Σⱼ 2λⱼqⱼ/(z − βⱼ).
//
// The outer spectral factor G(u) = W + L(uI − A)⁻¹B, G(u)G(−u) = F(u), has the stable poles
// λⱼ, zeros μ = −√ζ in the closed left half-plane for the zeros ζ of N(z) = F(u)∏(z − βⱼ),
// and W = √(2D) (zero for D = 0, with a zero of G at infinity). For ρ₀ = 0, N(0) = 0
// exactly and G(0) = 0 is imposed analytically. The available storage of Y is ½xᵀXx with X
// ≥ 0 the solution of AᵀX + XA = −LᵀL, which then also satisfies XB = Cᵀ − LᵀW and W² = 2D:
// the minimal solution of the positive-real equations, because G is outer. Any other
// storage, for example that of a non-outer factor with a zero reflected into the right
// half-plane, is larger. Without division by D, a Riccati equation, or an artificial
// feedthrough.

#include "permittivitymodel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <Eigen/Dense>
#include <fmt/format.h>
#include <mfem.hpp>

namespace palace::permittivity
{

namespace
{

using Complex = std::complex<double>;
constexpr double eps = std::numeric_limits<double>::epsilon();

// The correctly rounded sum of the values, independent of their order (Shewchuk's
// nonoverlapping partials, rounded as Python's math.fsum). Repeated-pole aggregation with
// it gives the same result for any order or split of the residues of a pole.
double ExactSum(const std::vector<double> &xs)
{
  std::vector<double> partials;
  for (double x : xs)
  {
    std::size_t i = 0;
    for (std::size_t k = 0; k < partials.size(); k++)
    {
      double y = partials[k];
      if (std::abs(x) < std::abs(y))
      {
        std::swap(x, y);
      }
      const double hi = x + y, lo = y - (hi - x);
      if (lo != 0.0)
      {
        partials[i++] = lo;
      }
      x = hi;
    }
    partials.resize(i);
    partials.push_back(x);
  }
  std::size_t n = partials.size();
  if (n == 0)
  {
    return 0.0;
  }
  double hi = partials[--n], lo = 0.0;
  while (n > 0)
  {
    const double x = hi, y = partials[--n];
    hi = x + y;
    lo = y - (hi - x);
    if (lo != 0.0)
    {
      break;
    }
  }
  // Round half-even correction: the remaining partials have the sign of lo.
  if (n > 0 && ((lo < 0.0 && partials[n - 1] < 0.0) || (lo > 0.0 && partials[n - 1] > 0.0)))
  {
    const double y = 2.0 * lo, x = hi + y;
    if (y == x - hi)
    {
      hi = x;
    }
  }
  return hi;
}

// A power of two within a factor of two of x > 0.
double PowerOfTwo(double x)
{
  return std::ldexp(1.0, std::ilogb(x));
}

// A stable pole of the scaled RAW block, Re λ < 0, Im λ ≥ 0, with its conjugate if
// Im λ > 0, and its scaled residue ρ.
struct Mode
{
  Complex lambda, rho;
  bool Pair() const { return lambda.imag() > 0.0; }
};

int States(const std::vector<Mode> &modes)
{
  int n = 0;
  for (const auto &mode : modes)
  {
    n += mode.Pair() ? 2 : 1;
  }
  return n;
}

// Real value of the modal sum Σⱼ (vⱼ + conj for a pair), exactly rounded, and the sum of
// the magnitudes of its terms.
std::pair<double, double> ModalSum(const std::vector<Mode> &modes,
                                   const std::vector<Complex> &v)
{
  std::vector<double> parts;
  double scale = 0.0;
  for (std::size_t j = 0; j < modes.size(); j++)
  {
    const double w = modes[j].Pair() ? 2.0 : 1.0;
    parts.push_back(w * v[j].real());
    scale += w * std::abs(v[j]);
  }
  return {ExactSum(parts), scale};
}

// The term c/(z − β) of a real pole, or c/(z − β) + c̄/(z − β̄) of a conjugate pair, and
// its first two derivatives. A pair is evaluated in its real form (n₁z + n₀)/((z − β)(z −
// β̄)), n₁ = 2 Re c, n₀ = −2 Re(c β̄): for a weakly damped pair the two complex terms cancel
// to a value proportional to the damping, which the real form keeps to relative accuracy
// (exactly so for Re c = 0, the residue of a typed Lorentz).
std::array<Complex, 3> PairTerm(Complex c, Complex beta, Complex z, bool pair)
{
  if (!pair)
  {
    const Complex t = 1.0 / (z - beta);
    return {c * t, -c * t * t, 2.0 * c * t * t * t};
  }
  const double n1 = 2.0 * c.real(), n0 = -2.0 * (c * std::conj(beta)).real();
  const Complex Q = (z - beta) * (z - std::conj(beta)), dQ = 2.0 * (z - beta.real());
  const Complex T = (n1 * z + n0) / Q, dT = (n1 - T * dQ) / Q,
                ddT = (-2.0 * dT * dQ - 2.0 * T) / Q;
  return {T, dT, ddT};
}

// The finite zeros, with multiplicity, of the real rational function
//
//   R(z) = d + Σⱼ cⱼ/(z − βⱼ)  (conjugate terms implied for a pair),
//
// given that there are exactly count of them. They are the finite eigenvalues of the
// pencil [[A − zI, B], [C, d]] of a real realization, from Eigen's real QZ algorithm with
// homogeneous eigenvalues (α, β). For d = 0 one infinite eigenvalue is removed exactly by
// an orthogonal compression to (UᵀAV, UᵀV), U and V orthonormal bases of B⊥ and Cᵀ⊥, so
// that a relative degree of one (the common case) leaves no infinite eigenvalue at all.
// Remaining infinite eigenvalues are those with the smallest |β|/|(α, β)|, never a large
// finite zero in their place. Each kept eigenvalue's normwise backward error is computed
// from the smallest singular value of β K − α E.
struct Zeros
{
  // Simple zeros: real, or complex in the upper half-plane with the conjugate implied.
  std::vector<Complex> simple;
  // Double zeros ζ = −b² on the negative real axis (G zeros ±ib), as b², and the boundary
  // frequencies b of the computed axis zeros, simple or double.
  std::vector<double> axis, frequencies;
  double backward_error = 0.0, separation = std::numeric_limits<double>::infinity();
  std::string failure;
};

Zeros RationalZeros(const std::vector<Mode> &modes, const std::vector<Complex> &beta,
                    const std::vector<Complex> &c, double d, int count)
{
  Zeros out;
  if (count == 0)
  {
    return out;
  }
  const int n = States(modes);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
  Eigen::VectorXd B = Eigen::VectorXd::Zero(n), C = Eigen::VectorXd::Zero(n);
  for (std::size_t j = 0, k = 0; j < modes.size(); j++)
  {
    // Balance the input and output of each block, so that |B| = |C|.
    const double sigma = std::sqrt(std::abs(c[j]));
    A(k, k) = beta[j].real();
    B(k) = sigma;
    if (modes[j].Pair())
    {
      A(k, k + 1) = -beta[j].imag();
      A(k + 1, k) = beta[j].imag();
      A(k + 1, k + 1) = beta[j].real();
      C(k) = 2.0 * c[j].real() / sigma;
      C(k + 1) = -2.0 * c[j].imag() / sigma;
      k += 2;
    }
    else
    {
      C(k) = c[j].real() / sigma;
      k += 1;
    }
  }
  Eigen::MatrixXd K, E;
  if (d != 0.0)
  {
    K.setZero(n + 1, n + 1);
    E.setZero(n + 1, n + 1);
    K.topLeftCorner(n, n) = A;
    K.topRightCorner(n, 1) = B;
    K.bottomLeftCorner(1, n) = C.transpose();
    K(n, n) = d;
    E.topLeftCorner(n, n).setIdentity();
  }
  else
  {
    // Zeros z with R(z) = 0: x ∈ ker C with (A − zI)x ∈ span B, so x = Vy and
    // Uᵀ(A − zI)Vy = 0.
    const Eigen::MatrixXd QB = Eigen::HouseholderQR<Eigen::MatrixXd>(B).householderQ();
    const Eigen::MatrixXd QC = Eigen::HouseholderQR<Eigen::MatrixXd>(C).householderQ();
    const auto U = QB.rightCols(n - 1), V = QC.rightCols(n - 1);
    K = U.transpose() * A * V;
    E = U.transpose() * V;
  }
  Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> qz(K, E, false);
  if (qz.info() != Eigen::Success)
  {
    out.failure = "the QZ iteration for the spectral zeros did not converge";
    return out;
  }
  const Eigen::VectorXcd alpha = qz.alphas();
  const Eigen::VectorXd betas = qz.betas();
  const int size = static_cast<int>(alpha.size());
  if (count > size)
  {
    out.failure = "the spectral-zero pencil is smaller than the number of zeros";
    return out;
  }
  std::vector<double> finiteness(size);
  for (int i = 0; i < size; i++)
  {
    finiteness[i] = std::abs(betas(i)) / std::hypot(std::abs(alpha(i)), std::abs(betas(i)));
  }
  std::vector<int> order(size);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(),
                   [&](int a, int b) { return finiteness[a] > finiteness[b]; });
  if (count < size)
  {
    out.separation = (finiteness[order[count]] > 0.0)
                         ? finiteness[order[count - 1]] / finiteness[order[count]]
                         : std::numeric_limits<double>::infinity();
  }
  const double norm_K = K.norm(), norm_E = E.norm();
  std::vector<Complex> raw;
  int lower = 0;
  for (int k = 0; k < count; k++)
  {
    const int i = order[k];
    if (!(betas(i) != 0.0) || !std::isfinite(alpha(i).real()) ||
        !std::isfinite(alpha(i).imag()))
    {
      out.failure = "a spectral zero is not finite";
      return out;
    }
    const double h = std::hypot(std::abs(alpha(i)), std::abs(betas(i)));
    const Complex a = alpha(i) / h, b = betas(i) / h;
    const Eigen::MatrixXcd M = b * K.cast<Complex>() - a * E.cast<Complex>();
    const double sigma = Eigen::JacobiSVD<Eigen::MatrixXcd>(M).singularValues().minCoeff();
    out.backward_error =
        std::max(out.backward_error, sigma / (std::abs(b) * norm_K + std::abs(a) * norm_E));
    const Complex z = alpha(i) / betas(i);
    (z.imag() < 0.0) ? void(lower++) : raw.push_back(z);
  }
  if (static_cast<int>(raw.size()) + lower != count ||
      std::count_if(raw.begin(), raw.end(), [](Complex z) { return z.imag() > 0.0; }) !=
          lower)
  {
    out.failure = "the finite spectral zeros are not closed under conjugation";
    return out;
  }

  // A zero on the imaginary axis (ζ < 0) has even multiplicity for a nonnegative spectral
  // density, and the QZ algorithm splits a double zero into a pair ζ ≈ −b² ± δ, real or
  // conjugate, with δ ~ √ε. Both factor F to rounding, but −√ζ for the pair would move the
  // zeros of G off the axis by O(δ), and the available storage with them: it is only
  // √-conditioned at a boundary zero, and the certificate (sensitive to F, so to δ²)
  // cannot see the difference. A pair is therefore one double zero when that perturbs its
  // quadratic (z − ζ₁)(z − ζ₂) by no more than the rounding of the zero computation, and
  // otherwise two simple boundary zeros, where F changes sign.
  const double cluster = 64.0 * (States(modes) + 1) * eps;
  auto Clustered = [cluster](Complex z1, Complex z2)
  { return std::norm(0.5 * (z1 - z2)) <= cluster * std::abs(z1 * z2); };
  std::vector<double> negative, centers;
  for (const auto z : raw)
  {
    if (z.imag() > 0.0 && z.real() < 0.0 && Clustered(z, std::conj(z)))
    {
      centers.push_back(z.real());
    }
    else if (z.imag() == 0.0 && z.real() < 0.0)
    {
      negative.push_back(z.real());
    }
    else
    {
      out.simple.push_back(z);
    }
  }
  std::sort(negative.begin(), negative.end());
  for (std::size_t i = 0; i < negative.size(); i += 2)
  {
    if (i + 1 == negative.size() || !Clustered(negative[i], negative[i + 1]))
    {
      for (double z : negative)
      {
        out.frequencies.push_back(std::sqrt(-z));
      }
      out.failure = "an imaginary-axis spectral zero is simple";
      return out;
    }
    centers.push_back(0.5 * (negative[i] + negative[i + 1]));
  }

  // QZ is backward stable for the pencil, so a zero has an ABSOLUTE error ~ ε‖K‖, while
  // z = u² squares the range of rates: a slow zero would lose most of its digits. Newton
  // steps on R (on R' for a double zero, whose center is a simple zero of R') recover the
  // relative accuracy. A step is taken only while it decreases the residual, and only
  // within the zero's class: a real zero stays real with its sign, a complex one in the
  // upper half-plane, a double boundary zero negative.
  auto Evaluate = [&](Complex z)
  {
    std::array<Complex, 3> r = {d, 0.0, 0.0};
    for (std::size_t j = 0; j < modes.size(); j++)
    {
      const auto t = PairTerm(c[j], beta[j], z, modes[j].Pair());
      for (int i = 0; i < 3; i++)
      {
        r[i] += t[i];
      }
    }
    return r;
  };
  auto Refine = [&](Complex z, int derivative, auto &&valid)
  {
    auto r = Evaluate(z);
    for (int it = 0; it < 16; it++)
    {
      if (r[derivative] == 0.0 || r[derivative + 1] == 0.0)
      {
        break;
      }
      Complex step = r[derivative] / r[derivative + 1];
      if (z.imag() == 0.0)
      {
        step = step.real();
      }
      const Complex next = z - step;
      const auto rn = Evaluate(next);
      if (!valid(next) || !(std::abs(rn[derivative]) < std::abs(r[derivative])))
      {
        break;
      }
      z = next;
      r = rn;
      if (std::abs(step) <= 4.0 * eps * std::abs(z))
      {
        break;
      }
    }
    return z;
  };
  for (auto &z : out.simple)
  {
    const Complex z0 = z;
    z = Refine(z0, 0,
               [z0](Complex x)
               {
                 return (z0.imag() > 0.0)
                            ? x.imag() > 0.0
                            : (x.imag() == 0.0 && (x.real() < 0.0) == (z0.real() < 0.0));
               });
  }
  for (double z : centers)
  {
    out.axis.push_back(-Refine(z, 1, [](Complex x) { return x.real() < 0.0; }).real());
  }
  for (double b2 : out.axis)
  {
    out.frequencies.push_back(std::sqrt(b2));
  }
  return out;
}

// Blockwise backward error max |R_IJ|/|S_IJ| over the state blocks I, J (a real pole, or a
// conjugate pair whose two real coordinates are one mode: the error must not depend on the
// arbitrary orientation of x = (Re ξ, Im ξ) in its plane), with max norms, and R = 0
// required where S = 0. A one-column (or one-row) residual has a single block there. The
// blocks are given by their first index, ending with the dimension.
double Blockwise(const Eigen::MatrixXd &R, const Eigen::MatrixXd &S,
                 const std::vector<int> &blocks)
{
  const std::vector<int> single = {0, 1};
  const auto &rows = (R.rows() == 1) ? single : blocks;
  const auto &cols = (R.cols() == 1) ? single : blocks;
  double err = 0.0;
  for (std::size_t I = 0; I + 1 < rows.size(); I++)
  {
    for (std::size_t J = 0; J + 1 < cols.size(); J++)
    {
      const int m = rows[I + 1] - rows[I], n = cols[J + 1] - cols[J];
      const double r = R.block(rows[I], cols[J], m, n).cwiseAbs().maxCoeff(),
                   s = S.block(rows[I], cols[J], m, n).cwiseAbs().maxCoeff();
      if (!std::isfinite(r) || !std::isfinite(s))
      {
        return std::numeric_limits<double>::infinity();
      }
      if (r > 0.0)
      {
        err = std::max(err, (s > 0.0) ? r / s : std::numeric_limits<double>::infinity());
      }
    }
  }
  return err;
}

std::string Sci(double x)
{
  return fmt::format("{:.3e}", x);
}

}  // namespace

PreparedStorage::PreparedStorage(const std::vector<config::PermittivityTermData> &terms)
{
  std::vector<config::PermittivityPoleData> raw;
  for (const auto &term : terms)
  {
    if (const auto *t = std::get_if<config::PermittivityPoleData>(&term))
    {
      raw.push_back(*t);
    }
    else
    {
      physical.push_back(term);
    }
  }
  has_raw = !raw.empty();
  if (has_raw)
  {
    PrepareRaw(raw);
  }
  else
  {
    reason = "no PoleResidue terms";
  }
}

void PreparedStorage::PrepareRaw(const std::vector<config::PermittivityPoleData> &input)
{
  auto Finish = [this](RawStorageStatus s, std::string why)
  {
    status = s;
    reason = std::move(why);
  };

  // Canonical RAW block: sort by exact pole and sum the residues of identical poles
  // exactly. Nearby poles are not merged, exactly canceled poles are removed.
  auto raw = input;
  for (const auto &t : raw)
  {
    MFEM_VERIFY(std::isfinite(t.pole.real()) && std::isfinite(t.pole.imag()) &&
                    std::isfinite(t.residue.real()) && std::isfinite(t.residue.imag()) &&
                    t.pole.real() <= 0.0 && t.pole.imag() >= 0.0 &&
                    (t.pole.imag() > 0.0 || t.residue.imag() == 0.0),
                "Invalid material permittivity pole-residue term (see CheckTerms)!");
  }
  auto less = [](Complex a, Complex b)
  { return a.real() < b.real() || (a.real() == b.real() && a.imag() < b.imag()); };
  std::sort(raw.begin(), raw.end(),
            [&](const auto &a, const auto &b) { return less(a.pole, b.pole); });
  std::vector<std::pair<Complex, Complex>> poles;
  for (std::size_t i = 0; i < raw.size();)
  {
    std::vector<double> re, im;
    std::size_t j = i;
    for (; j < raw.size() && raw[j].pole == raw[i].pole; j++)
    {
      re.push_back(raw[j].residue.real());
      im.push_back(raw[j].residue.imag());
    }
    const Complex residue(ExactSum(re), ExactSum(im));
    if (residue != 0.0)
    {
      poles.emplace_back(raw[i].pole, residue);
    }
    i = j;
  }

  // Exact power-of-two rate and amplitude scales.
  double pole_max = 0.0, residue_max = 0.0;
  for (const auto &[p, r] : poles)
  {
    pole_max = std::max(pole_max, std::abs(p));
  }
  const double rate = (pole_max > 0.0) ? PowerOfTwo(pole_max) : 1.0;
  for (const auto &[p, r] : poles)
  {
    residue_max = std::max(residue_max, std::abs(r) / rate);
  }
  const double amplitude = (residue_max > 0.0) ? PowerOfTwo(residue_max) : 1.0;
  certificate.rate_scale = rate;
  certificate.amplitude_scale = amplitude;

  // The zero pole contributes only to D. An exact lossless pole iβ with ρ = iy is the
  // oscillator Y = 2q u/(u² + β²), q = −βy, positive-real only for q > 0 (Re ρ = 0).
  double rho0 = 0.0;
  std::vector<Mode> modes;
  // D is the exact sum of the INPUT residues (2 Re r for a pair): aggregation rounds, for
  // example a shared zero pole of several Drude terms, and D = 0 is structural for them.
  std::vector<double> feedthrough;
  for (const auto &t : input)
  {
    feedthrough.push_back(((t.pole.imag() > 0.0) ? 2.0 : 1.0) * t.residue.real() /
                          (rate * amplitude));
  }
  for (const auto &[p, r] : poles)
  {
    const Complex lambda = p / rate, rho = r / (rate * amplitude);
    if (lambda == 0.0)
    {
      rho0 = rho.real();
    }
    else if (lambda.real() == 0.0)
    {
      const double beta = lambda.imag(), q = -beta * rho.imag();
      if (rho.real() != 0.0)
      {
        const double rel = std::abs(rho.real()) / std::abs(rho);
        if (rel > tolerance)
        {
          return Finish(RawStorageStatus::NonpassiveRaw,
                        fmt::format("the lossless pole {}i has a residue with a nonzero "
                                    "real part (relative {}): Re Y(iω) changes sign at "
                                    "the resonance",
                                    Sci(p.imag()), Sci(rel)));
        }
        return Finish(RawStorageStatus::NumericalUnresolved,
                      fmt::format("the lossless pole {}i has a residue with a real part "
                                  "within the resolution (relative {})",
                                  Sci(p.imag()), Sci(rel)));
      }
      if (q < 0.0)
      {
        return Finish(RawStorageStatus::NonpassiveRaw,
                      fmt::format("the lossless pole {}i has a negative current residue: "
                                  "its resonance stores negative energy",
                                  Sci(p.imag())));
      }
      lossless.emplace_back(beta, 2.0 * q);
    }
    else
    {
      modes.push_back({lambda, rho});
    }
  }
  D = ExactSum(feedthrough);
  certificate.feedthrough = D;
  certificate.lossless = static_cast<int>(lossless.size());
  if (rho0 < 0.0)
  {
    return Finish(RawStorageStatus::NonpassiveRaw,
                  fmt::format("the zero pole has a negative residue {}: a negative "
                              "conductance, Re Y(0) < 0",
                              Sci(rho0 * rate * amplitude)));
  }
  if (D < 0.0)
  {
    double scale = 0.0;
    for (double x : feedthrough)
    {
      scale += std::abs(x);
    }
    if (-D > tolerance * scale)
    {
      return Finish(RawStorageStatus::NonpassiveRaw,
                    fmt::format("the sum of the residues is negative, {} (relative {}): "
                                "Re Y(iω) → D < 0 at high frequency",
                                Sci(D * amplitude), Sci(-D / scale)));
    }
    return Finish(RawStorageStatus::NumericalUnresolved,
                  fmt::format("the sum of the residues is negative within the resolution "
                              "(relative {})",
                              Sci(-D / scale)));
  }

  const int n = States(modes);
  certificate.states = n;
  for (const auto &mode : modes)
  {
    blocks.push_back({mode.lambda.real(), mode.lambda.imag()});
  }

  // Search for a frequency of demonstrated negative Re Y(iv), for a failed certificate:
  // NonpassiveRaw with it, NumericalUnresolved without.
  auto Unresolved = [&](const std::string &why, const std::vector<double> &extra = {})
  {
    blocks.clear();
    std::vector<double> grid = extra;
    double lo = std::numeric_limits<double>::infinity(), hi = 0.0;
    for (const auto &mode : modes)
    {
      lo = std::min(lo, std::abs(mode.lambda));
      hi = std::max(hi, std::abs(mode.lambda));
      if (mode.Pair())
      {
        const double h = std::max(-mode.lambda.real(), 1.0e-6 * mode.lambda.imag());
        for (int k = -64; k <= 64; k++)
        {
          grid.push_back(mode.lambda.imag() + 0.25 * k * h);
        }
      }
    }
    const int per_decade = 64;
    const int count = static_cast<int>(std::ceil(per_decade * (std::log10(hi / lo) + 6.0)));
    for (int k = 0; k <= count; k++)
    {
      grid.push_back(1.0e-3 * lo * std::pow(10.0, static_cast<double>(k) / per_decade));
    }
    std::sort(grid.begin(), grid.end());
    for (double v : grid)
    {
      if (!(v > 0.0) || !std::isfinite(v))
      {
        continue;
      }
      // Re Y(iv) = D + Σ Re of each pole's term, a pair's in its real form
      // Re((n₀ + in₁v)/Q), Q = (iv − λ)(iv − λ̄), n₁ = 2 Re q, n₀ = −2 Re(qλ̄): the two
      // complex terms of a weakly damped pair cancel to a value proportional to its
      // damping, which this keeps to the rounding scale of its own two products.
      std::vector<double> parts = {D};
      double scale = std::abs(D);
      for (const auto &mode : modes)
      {
        const Complex q = mode.lambda * mode.rho;
        if (mode.Pair())
        {
          const double n1 = 2.0 * q.real(), n0 = -2.0 * (q * std::conj(mode.lambda)).real();
          const Complex Q = Complex(-mode.lambda.real(), v - mode.lambda.imag()) *
                            Complex(-mode.lambda.real(), v + mode.lambda.imag());
          const double a = n0 * Q.real() / std::norm(Q),
                       b = n1 * v * Q.imag() / std::norm(Q);
          parts.push_back(a);
          parts.push_back(b);
          scale += std::abs(a) + std::abs(b);
        }
        else
        {
          const double part = -q.real() * mode.lambda.real() /
                              (v * v + mode.lambda.real() * mode.lambda.real());
          parts.push_back(part);
          scale += std::abs(part);
        }
      }
      const double re_y = ExactSum(parts);
      if (re_y < -tolerance * scale)
      {
        return Finish(RawStorageStatus::NonpassiveRaw,
                      fmt::format("{}; the RAW block has gain: Re Y(iω) = {} (relative {}) "
                                  "at ω = {}",
                                  why, Sci(re_y * amplitude * rate), Sci(-re_y / scale),
                                  Sci(v * rate)));
      }
    }
    Finish(RawStorageStatus::NumericalUnresolved,
           why + "; no frequency of negative Re Y(iω) was found");
  };

  if (n == 0)
  {
    // Pure feedthrough (conductance) and lossless oscillators: no other storage.
    W = std::sqrt(2.0 * D);
    certificate.feedthrough_residual =
        (D > 0.0) ? std::abs(W * W - 2.0 * D) / (2.0 * D) : 0.0;
    return Finish(RawStorageStatus::CheckedAvailable,
                  "checked: no stable states, the storage is that of the lossless "
                  "oscillators");
  }

  // Zeros ζ of N(z), with their count from the leading coefficient of N.
  const bool dc_zero = (rho0 == 0.0);
  std::vector<Complex> beta, c;
  for (const auto &mode : modes)
  {
    beta.push_back(mode.lambda * mode.lambda);
    c.push_back(dc_zero ? 2.0 * mode.rho : 2.0 * beta.back() * mode.rho);
  }
  const double d = dc_zero ? 0.0 : 2.0 * D;
  int count = n;
  double lead = d;
  if (d == 0.0)
  {
    // The leading coefficient of H(z)R(z) is the first nonzero moment Σ cⱼβⱼᵏ, of degree
    // n − k − 1. A moment below the rounding of its terms is taken as zero; a wrong
    // decision fails the certificate rather than passing it.
    int degree = 0;
    std::vector<Complex> moment_terms = c;
    for (int k = 0; k < n && degree == 0; k++)
    {
      const auto [moment, scale] = ModalSum(modes, moment_terms);
      if (std::abs(moment) > 16.0 * (n + 1) * eps * scale)
      {
        degree = k + 1;
        lead = moment;
      }
      for (std::size_t j = 0; j < modes.size(); j++)
      {
        moment_terms[j] *= beta[j];
      }
    }
    if (degree == 0)
    {
      return Unresolved("the spectral density has no resolved leading coefficient");
    }
    count = n - degree;
  }
  const int m = count + (dc_zero ? 1 : 0);
  certificate.spectral_zeros = m;
  const double k2 = ((n - m) % 2) ? -lead : lead;
  if (!(k2 > 0.0))
  {
    return Unresolved("the spectral density is negative at high frequency");
  }
  const auto zeros = RationalZeros(modes, beta, c, d, count);
  certificate.zero_backward_error = zeros.backward_error;
  certificate.infinite_separation = zeros.separation;
  // Boundary frequencies and their geometric midpoints for the search for gain.
  std::vector<double> frequencies = zeros.frequencies;
  std::sort(frequencies.begin(), frequencies.end());
  for (std::size_t i = 1, size = frequencies.size(); i < size; i++)
  {
    frequencies.push_back(std::sqrt(frequencies[i - 1] * frequencies[i]));
  }
  if (!zeros.failure.empty())
  {
    return Unresolved(zeros.failure, frequencies);
  }

  // Outer zeros μ = −√ζ in the closed left half-plane, closed under conjugation, and ±ib
  // for each double zero on the imaginary axis.
  std::vector<Complex> mu;
  for (const auto zeta : zeros.simple)
  {
    if (zeta.imag() > 0.0)
    {
      mu.push_back(-std::sqrt(zeta));
      mu.push_back(std::conj(mu.back()));
    }
    else
    {
      mu.push_back(-std::sqrt(zeta.real()));
    }
  }
  for (const double b2 : zeros.axis)
  {
    mu.emplace_back(0.0, std::sqrt(b2));
    mu.emplace_back(0.0, -std::sqrt(b2));
  }
  if (dc_zero)
  {
    mu.emplace_back(0.0, 0.0);
  }
  MFEM_ASSERT(static_cast<int>(mu.size()) == m, "Unexpected spectral zero count!");

  // Residues l of G = M/a at each pole, M(u) = k∏(u − μ), a(u) = ∏(u − λ), as interleaved
  // products of ratios, with the conjugate of the upper pole's for its conjugate.
  const double k = std::sqrt(k2);
  W = (m == n) ? k : 0.0;
  std::vector<Complex> lambda, l, q;
  for (const auto &mode : modes)
  {
    lambda.push_back(mode.lambda);
    q.push_back(mode.lambda * mode.rho);
    if (mode.Pair())
    {
      lambda.push_back(std::conj(mode.lambda));
      q.push_back(std::conj(q.back()));
    }
  }
  for (int j = 0; j < n; j++)
  {
    if (j > 0 && lambda[j] == std::conj(lambda[j - 1]) && lambda[j].imag() < 0.0)
    {
      l.push_back(std::conj(l.back()));
      continue;
    }
    std::vector<Complex> den;
    for (int i = 0; i < n; i++)
    {
      if (i != j)
      {
        den.push_back(lambda[j] - lambda[i]);
      }
    }
    Complex v = k;
    for (int t = 0; t < std::max(m, n - 1); t++)
    {
      if (t < m)
      {
        v *= lambda[j] - mu[t];
      }
      if (t < n - 1)
      {
        v /= den[t];
      }
    }
    l.push_back(v);
  }

  // Complex symmetric modal solution P = −lᵢlⱼ/(λᵢ + λⱼ) of ΛP + PΛ = −lᵀl, and real
  // coordinates ξ = Tx: X = TᵀPT, L = lT, C = qT.
  Eigen::MatrixXcd P(n, n), T = Eigen::MatrixXcd::Zero(n, n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      P(i, j) = -l[i] * l[j] / (lambda[i] + lambda[j]);
    }
  }
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
  Eigen::VectorXd Bv = Eigen::VectorXd::Zero(n);
  std::vector<int> starts;
  for (std::size_t b = 0, i = 0; b < modes.size(); b++)
  {
    starts.push_back(static_cast<int>(i));
    A(i, i) = modes[b].lambda.real();
    Bv(i) = 1.0;
    T(i, i) = 1.0;
    if (modes[b].Pair())
    {
      A(i, i + 1) = -modes[b].lambda.imag();
      A(i + 1, i) = modes[b].lambda.imag();
      A(i + 1, i + 1) = modes[b].lambda.real();
      T(i, i + 1) = Complex(0.0, 1.0);
      T(i + 1, i) = 1.0;
      T(i + 1, i + 1) = Complex(0.0, -1.0);
      i += 2;
    }
    else
    {
      i += 1;
    }
  }
  starts.push_back(n);
  const Eigen::RowVectorXcd lv = Eigen::Map<const Eigen::RowVectorXcd>(l.data(), n),
                            qv = Eigen::Map<const Eigen::RowVectorXcd>(q.data(), n);
  const Eigen::MatrixXcd Xc = T.transpose() * P * T;
  const Eigen::RowVectorXcd Lc = lv * T, Cc = qv * T;
  {
    const Eigen::MatrixXd S = T.cwiseAbs().transpose() * P.cwiseAbs() * T.cwiseAbs();
    const Eigen::RowVectorXd SL = lv.cwiseAbs() * T.cwiseAbs(),
                             SC = qv.cwiseAbs() * T.cwiseAbs();
    certificate.reality_error =
        std::max({Blockwise(Xc.imag(), S, starts), Blockwise(Lc.imag(), SL, starts),
                  Blockwise(Cc.imag(), SC, starts)});
  }
  if (!(certificate.reality_error <= tolerance))
  {
    return Unresolved(fmt::format("the real-coordinate storage has an imaginary part "
                                  "(relative {})",
                                  Sci(certificate.reality_error)),
                      frequencies);
  }
  const Eigen::MatrixXd X0 = 0.5 * (Xc.real() + Xc.real().transpose());
  const Eigen::RowVectorXd Lr = Lc.real(), Cr = Cc.real();

  // Factor X = HHᵀ of the diagonally equilibrated matrix (unit diagonal) by Cholesky with
  // diagonal pivoting. For nearly coincident poles X is nearly singular in exact
  // arithmetic, and rounding can leave a tiny negative Schur complement. The factorization
  // stops when the largest remaining pivot is within its rounding, 64nε: the remainder,
  // bounded by that in magnitude, is dropped, a bounded backward perturbation of X. A
  // remainder beyond it, or a negative pivot, leaves the block unresolved (no eigenvalue
  // clipping). Every check below is of the ACTUAL X = HHᵀ.
  const Eigen::VectorXd diag = X0.diagonal();
  if (!(diag.minCoeff() > 0.0) || !diag.allFinite())
  {
    return Unresolved("the storage matrix has a nonpositive diagonal", frequencies);
  }
  const Eigen::VectorXd s = diag.cwiseSqrt();
  Eigen::MatrixXd Sc = s.cwiseInverse().asDiagonal() * X0 * s.cwiseInverse().asDiagonal();
  Eigen::MatrixXd Hp = Eigen::MatrixXd::Zero(n, n);
  std::vector<int> perm(n);
  std::iota(perm.begin(), perm.end(), 0);
  const double floor = 64.0 * n * eps;
  int rank = 0;
  certificate.min_pivot = 1.0;
  for (; rank < n; rank++)
  {
    Eigen::Index p;
    const double pivot = Sc.diagonal().tail(n - rank).maxCoeff(&p);
    if (!(pivot > floor))
    {
      break;
    }
    p += rank;
    Sc.row(rank).swap(Sc.row(p));
    Sc.col(rank).swap(Sc.col(p));
    Hp.row(rank).swap(Hp.row(p));
    std::swap(perm[rank], perm[p]);
    const double h = std::sqrt(pivot);
    const int rest = n - rank - 1;
    Hp(rank, rank) = h;
    Hp.col(rank).tail(rest) = Sc.col(rank).tail(rest) / h;
    Sc.bottomRightCorner(rest, rest).noalias() -=
        Hp.col(rank).tail(rest) * Hp.col(rank).tail(rest).transpose();
    certificate.min_pivot = std::min(certificate.min_pivot, pivot);
  }
  certificate.zeroed_pivots = n - rank;
  if (rank < n)
  {
    const double remainder = Sc.bottomRightCorner(n - rank, n - rank).cwiseAbs().maxCoeff();
    if (!(remainder <= floor))
    {
      return Unresolved(fmt::format("the storage matrix is indefinite beyond rounding "
                                    "(Schur complement {})",
                                    Sci(remainder)),
                        frequencies);
    }
    certificate.min_pivot = 0.0;
  }
  Eigen::MatrixXd H(n, n);
  for (int i = 0; i < n; i++)
  {
    H.row(perm[i]) = s(perm[i]) * Hp.row(i);
  }

  // Certificate of the ACTUAL X = HHᵀ.
  const Eigen::MatrixXd X = H * H.transpose();
  const Eigen::MatrixXd SX = H.cwiseAbs() * H.cwiseAbs().transpose();
  certificate.lyapunov_residual =
      Blockwise(A.transpose() * X + X * A + Lr.transpose() * Lr,
                A.cwiseAbs().transpose() * SX + SX * A.cwiseAbs() +
                    Lr.cwiseAbs().transpose() * Lr.cwiseAbs(),
                starts);
  certificate.cross_residual = Blockwise(X * Bv - Cr.transpose() + Lr.transpose() * W,
                                         SX * Bv.cwiseAbs() + Cr.cwiseAbs().transpose() +
                                             Lr.cwiseAbs().transpose() * std::abs(W),
                                         starts);
  certificate.feedthrough_residual =
      (D > 0.0) ? std::abs(W * W - 2.0 * D) / (2.0 * D)
                : ((W == 0.0) ? 0.0 : std::numeric_limits<double>::infinity());
  // F at probes from the real form of each pole's term in z = u² (as for the zeros), so
  // that the check is relative to F itself even where weakly damped terms cancel; G from
  // its poles.
  for (const Complex u : {Complex(0.0, 0.37), Complex(0.0, 1.3), Complex(0.0, 4.1),
                          Complex(0.6, 0.0), Complex(0.25, 0.8)})
  {
    Complex F = 2.0 * D, Gp = W, Gm = W;
    double scale = 2.0 * D, gp = std::abs(W), gm = std::abs(W);
    for (std::size_t j = 0; j < modes.size(); j++)
    {
      const Complex f = 2.0 * beta[j] * modes[j].rho;
      const Complex t = PairTerm(f, beta[j], u * u, modes[j].Pair())[0];
      F += t;
      scale += std::abs(t);
    }
    for (int j = 0; j < n; j++)
    {
      const Complex p = l[j] / (u - lambda[j]), m = l[j] / (-u - lambda[j]);
      Gp += p;
      Gm += m;
      gp += std::abs(p);
      gm += std::abs(m);
    }
    certificate.factor_residual =
        std::max(certificate.factor_residual, std::abs(F - Gp * Gm) / (scale + gp * gm));
  }
  const double worst =
      std::max({certificate.zero_backward_error, certificate.lyapunov_residual,
                certificate.cross_residual, certificate.feedthrough_residual,
                certificate.factor_residual});
  if (!(worst <= tolerance))
  {
    return Unresolved(
        fmt::format("the storage certificate failed: residuals (zero {}, Lyapunov {}, "
                    "cross {}, feedthrough {}, factor {}), equilibrated pivot {}",
                    Sci(certificate.zero_backward_error),
                    Sci(certificate.lyapunov_residual), Sci(certificate.cross_residual),
                    Sci(certificate.feedthrough_residual), Sci(certificate.factor_residual),
                    Sci(certificate.min_pivot)),
        frequencies);
  }

  factor.assign(n * n, 0.0);
  C.assign(Cr.data(), Cr.data() + n);
  L.assign(Lr.data(), Lr.data() + n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      factor[j * n + i] = H(i, j);  // Hᵀ, row-major
    }
  }
  Finish(RawStorageStatus::CheckedAvailable,
         fmt::format("checked: largest residual {}, equilibrated pivot {}{}", Sci(worst),
                     Sci(certificate.min_pivot),
                     certificate.zeroed_pivots
                         ? fmt::format(", {} trailing pivots within rounding dropped",
                                       certificate.zeroed_pivots)
                         : std::string()));
}

std::optional<double> PreparedStorage::RawStoredEnergySusceptibility(double omega) const
{
  if (status != RawStorageStatus::CheckedAvailable)
  {
    return std::nullopt;
  }
  if (!has_raw)
  {
    return 0.0;
  }
  const double v = omega / certificate.rate_scale;
  const Complex s(0.0, v);
  const int n = certificate.states;
  std::vector<Complex> z(n);
  for (std::size_t b = 0, i = 0; b < blocks.size(); b++)
  {
    // (iv − A)⁻¹B for each block, B = 1 or (1, 0), with (iv − a)² + w² in product form,
    // accurate near a weakly damped resonance.
    const double a = blocks[b].re, w = blocks[b].im;
    if (w > 0.0)
    {
      const Complex den = Complex(-a, v - w) * Complex(-a, v + w);
      z[i] = (s - a) / den;
      z[i + 1] = w / den;
      i += 2;
    }
    else
    {
      z[i] = 1.0 / (s - a);
      i += 1;
    }
  }
  double energy = 0.0;
  for (int j = 0; j < n; j++)
  {
    Complex y = 0.0;
    for (int i = 0; i < n; i++)
    {
      y += factor[j * n + i] * z[i];
    }
    energy += std::norm(y);
  }
  for (const auto &[beta, strength] : lossless)
  {
    const double den = (beta - v) * (beta + v);
    energy += strength * (v * v + beta * beta) / (den * den);
  }
  return certificate.amplitude_scale * energy;
}

std::optional<double> PreparedStorage::StoredEnergySusceptibility(double omega) const
{
  const auto raw = RawStoredEnergySusceptibility(omega);
  if (!raw)
  {
    return std::nullopt;
  }
  return permittivity::StoredEnergySusceptibility(physical, omega) + *raw;
}

RawStorageRealization PreparedStorage::Realization() const
{
  RawStorageRealization r;
  r.n = (status == RawStorageStatus::CheckedAvailable) ? certificate.states : 0;
  const int n = r.n;
  r.A.assign(n * n, 0.0);
  r.B.assign(n, 0.0);
  for (std::size_t b = 0, i = 0; i < static_cast<std::size_t>(n); b++)
  {
    r.A[i * n + i] = blocks[b].re;
    r.B[i] = 1.0;
    if (blocks[b].im > 0.0)
    {
      r.A[i * n + i + 1] = -blocks[b].im;
      r.A[(i + 1) * n + i] = blocks[b].im;
      r.A[(i + 1) * n + i + 1] = blocks[b].re;
      i += 2;
    }
    else
    {
      i += 1;
    }
  }
  r.C = C;
  r.L = L;
  r.X.assign(n * n, 0.0);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      for (int k = 0; k < n; k++)
      {
        r.X[i * n + j] += factor[k * n + i] * factor[k * n + j];
      }
    }
  }
  r.D = D;
  r.W = W;
  r.lossless = lossless;
  return r;
}

}  // namespace palace::permittivity
