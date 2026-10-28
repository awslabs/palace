// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_PERMITTIVITY_MODEL_HPP
#define PALACE_UTILS_PERMITTIVITY_MODEL_HPP

#include <complex>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include "utils/configfile.hpp"

namespace palace::permittivity
{

//
// Frequency-dependent relative permittivity ε(s) = ε∞ + Σₖ χₖ(s), s = iω, of a material
// with additive terms (config::PermittivityTermData). Rates may be in any consistent unit.
//
// Each term is the polarization P = ε₀χE of an auxiliary differential equation driven by E,
// for which the power delivered to the polarization splits exactly as E·∂P/∂t = dW/dt + D,
// with W the stored energy and D ≥ 0 the dissipated power:
//
//   Relaxation  ẋ + a₀x = E, P = ε₀b₀x:  W = ½ε₀a₀b₀x², D = ε₀b₀ẋ².
//   Oscillator  ẍ + a₁ẋ + a₀x = E, P = ε₀(b₀x + b₁ẋ):
//               W = ε₀(½(b₀ + a₁b₁)ẋ² + a₀b₁xẋ + ½a₀b₀x²), D = ε₀(b₁ẍ² + (a₁b₀ − a₀b₁)ẋ²).
//
// Debye is a relaxation and Djordjevic-Sarkar a positive continuum of them; Lorentz and
// Drude (a₀ = 0) are oscillators with b₁ = 0. For these, W is the only quadratic form of
// the state for which the identity holds for every E, so the stored energy is that of the
// physical model, and positive parameters make W and D nonnegative. A PoleResidue term is a
// relaxation for a real pole and an oscillator for a conjugate pair. With b₁ ≠ 0 the stored
// energy of a passive term is not unique: W + ½δẋ², with |δ| ≤ 2√(b₁(a₁b₀ − a₀b₁)), keeps D
// nonnegative for the same χ. The W above is the midpoint of this interval.
//

// Pole-residue form of the terms other than Djordjevic-Sarkar, for the solver: real poles
// and upper-half-plane poles, whose conjugates are implied, including exact zero poles.
std::vector<config::PermittivityPoleData>
PoleTerms(const std::vector<config::PermittivityTermData> &terms);

// The Djordjevic-Sarkar terms.
std::vector<config::DjordjevicSarkarData>
DjordjevicSarkarTerms(const std::vector<config::PermittivityTermData> &terms);

// Susceptibility χ(s) = Σₖ χₖ(s) at complex s.
std::complex<double> Susceptibility(const std::vector<config::PermittivityTermData> &terms,
                                    std::complex<double> s);

// χ(0), infinite for free carriers (Drude, or a PoleResidue pole at zero).
double StaticSusceptibility(const std::vector<config::PermittivityTermData> &terms);

// The stored electric energy of the polarization at real ω > 0, as a susceptibility: the
// time-averaged stored energy of a field E e^{iωt} in a material with these terms is that
// of a constant permittivity ε∞ + StoredEnergySusceptibility(terms, ω). For a relaxation
// this is Re χ; an oscillator adds the kinetic energy of its motion.
double StoredEnergySusceptibility(const std::vector<config::PermittivityTermData> &terms,
                                  double omega);

// The dissipated power at real ω > 0, as a susceptibility: −Im χ(iω), so that the
// time-averaged dissipated power of the field is that of a loss ω ε₀ (−Im χ) |E|²/2.
double DissipationSusceptibility(const std::vector<config::PermittivityTermData> &terms,
                                 double omega);

//
// Stored electric energy of a material's polarization, prepared once from its terms.
//
// Typed terms (Drude, Debye, Lorentz) and Djordjevic-Sarkar terms keep the physical
// energies of their models, above. All PoleResidue terms of the material together form ONE
// rational block, the RAW block χ_raw(s). Its stored energy is, by convention, the
// AVAILABLE STORAGE of its polarization current J = ∂P/∂t, with admittance Y(s) = s
// χ_raw(s): the largest energy that the polarization can return to the field, the smallest
// stored energy of any passive realization of Y. It depends only on χ_raw(s), not on how it
// is split into terms, their order, or repeated poles. It is NOT a microscopic energy: a
// sum of typed terms and a RAW block with the same χ need not store the same energy. For
// instance χ = 1/(s + 1) + 1/(s + 2) at ω = 1 stores 0.9 as two typed Debye terms and
// 0.8675444679663243 as a RAW block. For a single Drude, Debye, or Lorentz section the two
// agree.
//
// The storage is evaluated as κ‖Hᵀ(iω/Ω − A)⁻¹B‖² (plus the fixed energies of exact
// lossless resonances), for a minimal realization (A, B, C, D) of Y(Ωu)/(κΩ) and X = HHᵀ
// computed from the outer spectral factor of Y(u) + Y(−u) by the Lyapunov equation of its
// residues. What holds when:
//
// * Unconditionally, every finite evaluated value is a squared norm, so nonnegative.
// * If the computed identities (the positive-real or Kalman-Yakubovich-Popov equations, and
//   the outer allocation of the spectral zeros) held exactly, X would be a passive storage
//   of Y and the minimal one, the available storage.
// * CheckedAvailable says that this outer-factor storage approximation passed the
//   numerical checks of RawStorageCertificate, each at most PreparedStorage::tolerance.
//   They are scaled residuals of the identities, not one perturbation of (A, B, C, D) that
//   makes them all exact; they do not bound the error of the energy itself, nor prove
//   passivity or minimality for the given χ_raw. The available storage is ill-conditioned
//   near a spectral zero on the imaginary axis (a double zero, whose computed split pair is
//   projected onto the axis) and for a nearly singular X (nearly coincident poles, whose
//   trailing pivots within rounding are dropped), where an energy error well above the
//   residuals is possible.
//
// The realization, factor, and checks are computed once; each frequency costs O(n²) for n
// states. Ordinary conductivity and ε∞ are not part of the RAW block.
//
// Status of the RAW block. The whole material may be passive even when the RAW block alone
// is not (it is compensated by typed or Djordjevic-Sarkar terms): NonpassiveRaw is not a
// statement of gain of the material, and NumericalUnresolved is not a statement of gain at
// all. In both cases the RAW block has no stored energy under this convention.
enum class RawStorageStatus
{
  // No PoleResidue terms, or the outer-factor storage of the RAW block passed its numerical
  // checks.
  CheckedAvailable,

  // The RAW block alone is demonstrably not passive, beyond the resolution
  // PreparedStorage::tolerance: a negative high-frequency current D, a negative zero-pole
  // conductance, an invalid lossless residue, or a sampled frequency of negative Re Y(iω).
  NonpassiveRaw,

  // The checks failed in double precision (for example for nearly coincident poles), or the
  // block is outside the range of exact power-of-two scaling or overflows when aggregated,
  // and the RAW block was not demonstrated to be nonpassive either.
  NumericalUnresolved
};

// Diagnostics of the preparation of the RAW block, in the scaled variables u = s/Ω and
// χ/κ. The matrix residuals are blockwise: for each pair of state blocks (a real pole, or a
// conjugate pair as one mode), the largest residual entry relative to the largest entry of
// its rounding scale (|A|ᵀ|H||H|ᵀ + |H||H|ᵀ|A| + |L|ᵀ|L| and similar), so that a large fast
// state cannot hide an inaccurate slow state. They are numerical evidence that the
// identities hold to rounding, not a backward error of the data or an interval proof.
struct RawStorageCertificate
{
  // Stable real states n, exact lossless (Re p = 0) oscillators, and finite zeros of the
  // spectral factor G(u), counted with multiplicity, including any zero at u = 0.
  int states = 0, lossless = 0, spectral_zeros = 0;

  // Rate and amplitude scales Ω and κ (powers of two), and the scaled feedthrough D of
  // the current Y(u) → D as u → ∞, including the zero-pole residue.
  double rate_scale = 1.0, amplitude_scale = 1.0, feedthrough = 0.0;

  // Largest normwise backward error of a computed spectral zero as an eigenvalue of its
  // generalized eigenvalue pencil, and the ratio of the smallest |β|/|(α, β)| of a kept
  // (finite) eigenvalue to the largest of a discarded (infinite) one.
  double zero_backward_error = 0.0, infinite_separation = 0.0;

  // Imaginary part of the real-coordinate storage matrix and factor rows, relative.
  double reality_error = 0.0;

  // Smallest pivot of the diagonally pivoted Cholesky factorization of the diagonally
  // equilibrated storage matrix (one is perfectly conditioned), and the number of trailing
  // pivots within the rounding of the factorization, 64nε, whose remainder (bounded by
  // that) was dropped: a bounded perturbation of X, after which every residual is of the
  // actual factor, but whose result is then not the exact minimum. Positive definite when
  // zeroed_pivots = 0 (min_pivot is then positive).
  double min_pivot = 0.0;
  int zeroed_pivots = 0;

  // AᵀX + XA + LᵀL = 0, XB − Cᵀ + LᵀW = 0, and W² − 2D = 0 for the actual X = HHᵀ, and
  // F(u) − G(u)G(−u) at fixed probe points u.
  double lyapunov_residual = 0.0, cross_residual = 0.0, feedthrough_residual = 0.0,
         factor_residual = 0.0;
};

// Realization of the scaled RAW block, for verification only: χ_raw(s) = κ Y(s/Ω)/(s/Ω)
// with Y(u) = D + C(uI − A)⁻¹B on the stable states, plus the exact lossless oscillators,
// and storage X = HHᵀ with spectral factor G(u) = W + L(uI − A)⁻¹B. Matrices are row-major
// n×n (A, X) or length n (B, C, L).
struct RawStorageRealization
{
  int n = 0;
  std::vector<double> A, B, C, L, X;
  double D = 0.0, W = 0.0;
  // (β, 2q) of each exact lossless oscillator Y = 2q u/(u² + β²).
  std::vector<std::pair<double, double>> lossless;
};

class PreparedStorage
{
public:
  // The resolution of the RAW checks: the largest accepted scaled residual of
  // RawStorageCertificate, and the smallest violation, relative to the rounding scale of
  // its terms, that demonstrates nonpassivity. It is the tolerance FindDefects uses for the
  // sign of the dissipated power.
  static constexpr double tolerance = 1.0e-10;

  explicit PreparedStorage(const std::vector<config::PermittivityTermData> &terms);

  bool HasRaw() const { return has_raw; }
  RawStorageStatus Status() const { return status; }

  // Why the RAW block has its status, deterministic for the given terms.
  const std::string &Reason() const { return reason; }
  const RawStorageCertificate &Certificate() const { return certificate; }

  // StoredEnergySusceptibility of the typed and Djordjevic-Sarkar terms plus the
  // outer-factor storage of the RAW block, at finite real ω ≥ 0. Empty unless Status() is
  // CheckedAvailable. Infinite at an exact lossless resonance. The RAW part is otherwise
  // finite (it may underflow to zero at extreme ω); the typed formulas are those of
  // permittivity::StoredEnergySusceptibility.
  std::optional<double> StoredEnergySusceptibility(double omega) const;

  // The outer-factor storage of the RAW block alone at finite real ω ≥ 0, empty unless
  // Status() is CheckedAvailable (zero without PoleResidue terms). Finite except at an
  // exact lossless resonance, where it is infinite.
  std::optional<double> RawStoredEnergySusceptibility(double omega) const;

  RawStorageRealization Realization() const;

private:
  // Typed and Djordjevic-Sarkar terms.
  std::vector<config::PermittivityTermData> physical;

  bool has_raw = false;
  RawStorageStatus status = RawStorageStatus::CheckedAvailable;
  std::string reason;
  RawStorageCertificate certificate;

  // Stable blocks in canonical order: a real pole (im = 0) is one state, an
  // upper-half-plane pole with its conjugate two states, x = (Re ξ, Im ξ).
  struct Block
  {
    double re, im;
  };
  std::vector<Block> blocks;

  // Hᵀ (row-major n×n), C, L, D, W, and lossless oscillators.
  std::vector<double> factor, C, L;
  double D = 0.0, W = 0.0;
  std::vector<std::pair<double, double>> lossless;

  void PrepareRaw(const std::vector<config::PermittivityPoleData> &raw);
};

// Verify the parameters of the terms. A typed term must be passive on its own: a negative
// strength has gain at every frequency. ε∞ must be positive for a positive field energy.
void CheckTerms(double eps_inf, const std::vector<config::PermittivityTermData> &terms);

// Whether every term is a typed, passive model. Otherwise, use FindDefects.
bool IsPassiveByConstruction(const std::vector<config::PermittivityTermData> &terms);

// Frequency ranges ω ∈ [lower, upper] (upper possibly infinite, lower possibly zero) where
// the model has gain, a negative dissipated power, or a negative stored electric energy
// ε∞ + StoredEnergySusceptibility. Sampled on a logarithmic grid spanning the model's
// rates, refined at every resonance on the scale of its linewidth, with bisected endpoints.
struct FrequencyInterval
{
  double lower, upper;
};
struct ModelDefects
{
  std::vector<FrequencyInterval> gain, negative_energy;
};
ModelDefects FindDefects(double eps_inf,
                         const std::vector<config::PermittivityTermData> &terms);

// Whether ω lies in one of the intervals.
bool Contains(const std::vector<FrequencyInterval> &intervals, double omega);

}  // namespace palace::permittivity

#endif  // PALACE_UTILS_PERMITTIVITY_MODEL_HPP
