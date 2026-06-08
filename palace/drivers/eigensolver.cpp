// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <complex>
#include <cstdlib>
#include <map>
#include <string_view>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/arpack.hpp"
#include "linalg/divfree.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/ksp.hpp"
#include "linalg/nleps.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/slepc.hpp"
#include "linalg/vector.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

// Caches the ω-independent boundary mass per wave port (i·M^(p) form, see
// SpaceOperator::GetWavePortBoundaryMassMatrix), and exposes builders for
//   A2_wp(ω)   = Σ_p kₙ,p(ω) · Mwp_p   (ω-dependent only via scalar kₙ,p)
//   dA2_wp/dω = Σ_p (dkₙ,p/dω)(ω) · Mwp_p
// kₙ,p(ω) is sampled from real-ω 2D port EVPs as a *complex* scalar (β + iα where α
// captures wall-loss attenuation from impedance / loss-tan in the port cross-section)
// and stored verbatim — both real and imaginary parts feed the BC stamping and the kₙ²
// fit. The closed-form analytic continuation kₙ²(λ) = α₂ − γ_c·λ² uses *complex*
// coefficients so the imaginary part of kₙ at the fit endpoints is preserved off the
// imaginary axis.
//
// Regime of validity (closed-form continuation): exact for ports whose kₙ²(ω) is
// quadratic in ω with constant complex coefficients — covers PEC and impedance-walled
// lossless / quasi-lossless ports where the wall loss enters as a frequency-independent
// imaginary surface impedance. NOT valid for skin-depth (finite-conductivity wall)
// ports where Z_s ∝ √ω introduces non-quadratic structure, or for strongly dispersive
// media. The midpoint validator inside FitKnSq compares the complex fit to a third
// real-ω probe and warns if the relative error exceeds 1%.
// Pre-assembling Mwp_p once and re-summing scalar coefficients each call avoids
// re-running the bilinear form assembly on every NLEPS Newton iterate. The
// derivative uses central FD on the cheap (cross-section EVP) kₙ,p evaluation.
class WavePortFactor
{
public:
  WavePortFactor(SpaceOperator &space_op, double finite_diff_rel = 1.0e-4)
    : space_op_(space_op), fd_rel_(finite_diff_rel)
  {
    for (const auto &[port_idx, port_data] : space_op_.GetWavePortOp())
    {
      auto Mp = space_op_.GetWavePortBoundaryMassMatrix<ComplexOperator>(
          port_idx, Operator::DIAG_ZERO);
      if (Mp)
      {
        port_indices_.push_back(port_idx);
        Mwp_p_.push_back(std::move(Mp));
      }
    }
  }

  bool empty() const { return Mwp_p_.empty(); }

  // Enable the EXACT complex-frequency cross-section solve for the complex-λ build path
  // (BuildComplex / BuildDComplex). When set, kₙ(λ) is obtained by re-solving the 2D
  // modal EVP at the genuinely complex frequency ω = -i·λ (WavePortData::SolveKnExact)
  // instead of evaluating the closed-form kₙ²(ω) polynomial fit. This is the physically
  // exact analytic continuation and eliminates the fit error that otherwise floors
  // Newton convergence for near-cutoff non-TEM ports. NLEIGS (split form) cannot use
  // this — it needs the closed-form FN scalar — so it keeps the fit via Alpha/Gamma.
  void SetUseExactSolve(bool exact) { use_exact_ = exact; }
  bool UsesExactSolve() const { return use_exact_; }

  // Build Σ_p kₙ,p(ω) · Mwp_p as a fresh sum operator. Each call evaluates kₙ at the
  // requested ω (per-port cross-section EVP, cached internally by WavePortData).
  std::unique_ptr<ComplexOperator> Build(double omega) const
  {
    if (Mwp_p_.empty())
    {
      return {};
    }
    return BuildScaled(MakeKnCoefficients(omega));
  }

  // Exact complex-ω kₙ for the local port `k`, via the cross-section modal EVP solved
  // at ω = -i·λ. The Jacobian path needs a derivative too: dkₙ/dλ is approximated by a
  // central finite difference of SolveKnExact in λ (two extra EVP solves), which is
  // cheap (2D cross-section) and avoids assuming any closed form.
  std::complex<double> KnExact(int local_port, std::complex<double> lambda) const
  {
    // λ = i·ω  ⇒  ω = -i·λ = λ / i.
    const std::complex<double> omega = lambda / std::complex<double>(0.0, 1.0);
    return space_op_.GetWavePortOp().GetWavePortKnExact(port_indices_[local_port], omega);
  }

  // Build Σ_p kₙ,p(λ) · Mwp_p at COMPLEX λ using the analytical kₙ² fit (so kₙ
  // is evaluated as the analytic continuation √(α − γλ²) rather than the real-ω
  // restriction kₙ(|Im λ|)). Used in production by the complex-λ analytic-continuation
  // path (WavePortBCEvaluation = Complex). Caller must have called FitKnSq() first.
  // Optional kn_anchors[k] (one per port) keeps the propagating-sheet sign across
  // line search: if Re(kn · conj(anchor_k)) < 0 the sign of kn is flipped.
  std::unique_ptr<ComplexOperator>
  BuildComplex(std::complex<double> lam,
               const std::vector<std::complex<double>> *kn_anchors = nullptr) const
  {
    if (Mwp_p_.empty())
    {
      return {};
    }
    std::vector<std::complex<double>> coeffs(Mwp_p_.size());
    if (use_exact_)
    {
      // Exact analytic continuation: solve the cross-section EVP at ω = -i·λ per port.
      for (std::size_t k = 0; k < Mwp_p_.size(); k++)
      {
        coeffs[k] = KnExact(static_cast<int>(k), lam);
      }
      return BuildScaled(coeffs);
    }
    MFEM_VERIFY(kn2_coeffs_.size() == Mwp_p_.size(),
                "WavePortFactor::BuildComplex requires FitKnSq() first.");
    for (std::size_t k = 0; k < Mwp_p_.size(); k++)
    {
      const std::complex<double> anchor =
          (kn_anchors && kn_anchors->size() == Mwp_p_.size())
              ? (*kn_anchors)[k]
              : std::complex<double>{0.0, 0.0};
      coeffs[k] = KnComplex(static_cast<int>(k), lam, anchor);
    }
    return BuildScaled(coeffs);
  }

  // Build the analytic Newton-Jacobian wave-port term dA2/dλ at COMPLEX λ via the
  // closed-form derivative of kₙ(λ) = √P(-λ²) where P(s) = Σ_{j=0..D} c_j·s^j is
  // a polynomial in s = ω² fit to kₙ²(ω). The chain rule gives
  //     d kₙ / d λ = -λ · P'(-λ²) / kₙ(λ),  P'(s) = Σ_{j=1..D} j·c_j·s^{j-1}.
  // For the default D=1 case (kₙ² = c₀ + c₁·ω²) this reduces to -c₁·λ/kₙ, the
  // existing low-order formula. Used in production when WavePortBCEvaluation =
  // Complex; the polynomial-pencil path uses BuildJacobianTerm(double ω) instead.
  std::unique_ptr<ComplexOperator>
  BuildDComplex(std::complex<double> lam,
                const std::vector<std::complex<double>> *kn_anchors = nullptr) const
  {
    if (Mwp_p_.empty())
    {
      return {};
    }
    std::vector<std::complex<double>> coeffs(Mwp_p_.size());
    if (use_exact_)
    {
      // Exact path: dkₙ/dλ via central finite difference of the cross-section EVP solve
      // in λ. Two extra (cheap 2D) EVP solves per port. Step scales with |λ| for
      // conditioning; the recovered kₙ is smooth in λ away from cutoff so a modest
      // relative step is accurate.
      for (std::size_t k = 0; k < Mwp_p_.size(); k++)
      {
        const std::complex<double> h =
            fd_rel_ * std::max(std::abs(lam), 1.0) * std::complex<double>(0.0, 1.0);
        const std::complex<double> knp = KnExact(static_cast<int>(k), lam + h);
        const std::complex<double> knm = KnExact(static_cast<int>(k), lam - h);
        coeffs[k] = (knp - knm) / (2.0 * h);
      }
      return BuildScaled(coeffs);
    }
    MFEM_VERIFY(kn2_coeffs_.size() == Mwp_p_.size(),
                "WavePortFactor::BuildDComplex requires FitKnSq() first.");
    for (std::size_t k = 0; k < Mwp_p_.size(); k++)
    {
      const std::complex<double> anchor =
          (kn_anchors && kn_anchors->size() == Mwp_p_.size())
              ? (*kn_anchors)[k]
              : std::complex<double>{0.0, 0.0};
      const std::complex<double> kn = KnComplex(static_cast<int>(k), lam, anchor);
      // P'(s) at s = -λ² evaluated via Horner: dP_ds = Σ_{j≥1} j·c_j·s^{j-1}.
      const std::complex<double> s = -lam * lam;
      const auto &c = kn2_coeffs_[k];
      const std::size_t D = c.size() - 1;
      std::complex<double> dP_ds{0.0, 0.0};
      // Horner of dP/ds: dP/ds = c_1 + 2·c_2·s + 3·c_3·s² + ...
      for (std::size_t j = D; j >= 1; j--)
      {
        dP_ds = dP_ds * s + static_cast<double>(j) * c[j];
        if (j == 1)
        {
          break;  // unsigned underflow guard
        }
      }
      coeffs[k] = -lam * dP_ds / kn;
    }
    return BuildScaled(coeffs);
  }

  // Build the Newton-Jacobian wave-port term, dA2/dλ = -i·dA2/dω, evaluated at
  // λ = i·ω for real positive ω. With Mwp_p stored as i·M_real_p, BuildScaled
  // multiplying by a complex coefficient c = β + iα places (-α, β)·M_real on the
  // (real, imag) slots — so the chain dkₙ/dω·(-i factor) needs to live in c. We
  // therefore set coeffs[k] = -i · dkₙ/dω with kₙ a complex scalar (preserves wall
  // loss). Central FD on the cross-section EVP via the complex kₙ accessor; step is
  // fd_rel_·max(|ω|, 1).
  std::unique_ptr<ComplexOperator> BuildJacobianTerm(double omega) const
  {
    if (Mwp_p_.empty())
    {
      return {};
    }
    const double h = fd_rel_ * std::max(std::abs(omega), 1.0);
    std::vector<std::complex<double>> coeffs(Mwp_p_.size());
    for (std::size_t k = 0; k < port_indices_.size(); k++)
    {
      const std::complex<double> kp =
          space_op_.GetWavePortOp().GetWavePortKnComplex(port_indices_[k], omega + h);
      const std::complex<double> km =
          space_op_.GetWavePortOp().GetWavePortKnComplex(port_indices_[k], omega - h);
      const std::complex<double> dkn_domega = (kp - km) / (2.0 * h);
      // -i scale on a (i·M_real)-shaped operator gives M_real with the conjugated
      // factor (matching the kₙ → -i·dkₙ/dω convention used elsewhere).
      coeffs[k] = std::complex<double>{0.0, -1.0} * dkn_domega;
    }
    return BuildScaled(coeffs);
  }

  // Per-port kₙ,p(ω) (real propagating part). Backward-compat scalar accessor used
  // only by the per-port fit-residual diagnostic; production paths prefer
  // KnComplex / GetWavePortKnComplex which preserve wall-loss attenuation.
  double Kn(int local_port, double omega) const
  {
    return space_op_.GetWavePortOp().GetWavePortKn(port_indices_[local_port], omega);
  }

  std::size_t NumPorts() const { return port_indices_.size(); }

  // Per-port accessors used by the SLEPc NLEIGS split-form path: expose the
  // ω-independent boundary mass M_wp_p (constant matrix) and the kₙ² polynomial fit
  // coefficients (c_0, c_1, ..., c_D ∈ ℂ where kₙ²(ω) = Σ c_j·ω^{2j}) so the solver
  // can build per-port FN_SQRT(FN_RATIONAL) scalars without re-deriving them.
  // FitKnSq() must have been called first.
  const ComplexOperator *GetMatrix(std::size_t local_port) const
  {
    return Mwp_p_[local_port].get();
  }
  // Returns the per-port kₙ² polynomial coefficients in ω², lowest-degree first
  // (c_0 + c_1·ω² + c_2·ω⁴ + ...). Length is fit_order_ + 1.
  const std::vector<std::complex<double>> &Kn2Coeffs(std::size_t local_port) const
  {
    return kn2_coeffs_[local_port];
  }
  // Backward-compat scalar accessors. Alpha = c_0, Gamma = c_1. Only valid when
  // fit_order_ == 1 (the default); higher-order callers must use Kn2Coeffs.
  std::complex<double> Alpha(std::size_t local_port) const
  {
    return kn2_coeffs_[local_port][0];
  }
  std::complex<double> Gamma(std::size_t local_port) const
  {
    MFEM_VERIFY(kn2_coeffs_[local_port].size() >= 2,
                "WavePortFactor::Gamma: fit order must be >= 1.");
    return kn2_coeffs_[local_port][1];
  }
  std::size_t FitOrder() const
  {
    return kn2_coeffs_.empty() ? 0 : kn2_coeffs_[0].size() - 1;
  }

  // Fit kₙ²(ω) = Σ_{j=0..D} c_j·ω^{2j} per port with COMPLEX coefficients c_j ∈ ℂ.
  // D = order argument (default 1 = the historical "α + γ·ω²" two-coefficient fit,
  // exact for textbook closed/open waveguides and impedance-walled ports). Higher D
  // captures sub-quadratic deviations that arise near non-TEM port cutoffs, where
  // modal coupling and edge effects make kₙ²(ω) deviate from quadratic by several
  // percent — sufficient to keep Newton refinement of HYBRID-Complex stuck at the
  // fit error and prevent rational-interpolation NLEIGS from converging.
  //
  // Sampling: D + 1 equally-spaced real-ω anchors on [omega_min, omega_max], each
  // evaluated via GetWavePortKnComplex (re-running the cross-section EVP, cached
  // by WavePortData). Solve a (D+1)·(D+1) Vandermonde system in ω² for c_0..c_D.
  //
  // Midpoint validator: a halfway probe at ω_mid compares the truth kₙ²(ω_mid)
  // (complex) to the polynomial fit. If the relative error in complex magnitude
  // exceeds 1%, the model is inadequate for this port — typically a skin-depth
  // wall whose Z_s ∝ √ω makes kₙ²(ω) non-polynomial, or a dispersive ε(ω). The
  // analytic continuation is still USED in that case — the warning tells the user
  // the result may be less accurate than expected and suggests bumping fit_order.
  void FitKnSq(double omega_min, double omega_max, int order = 1)
  {
    MFEM_VERIFY(order >= 1, "WavePortFactor::FitKnSq: order must be >= 1.");
    kn2_coeffs_.clear();
    kn2_coeffs_.reserve(Mwp_p_.size());
    const std::size_t N = static_cast<std::size_t>(order) + 1;
    // (D+1) anchor points evenly spaced on [ω_min, ω_max]. ω_mid is reserved for
    // the midpoint validator (out-of-sample probe).
    std::vector<double> w_fit(N);
    for (std::size_t j = 0; j < N; j++)
    {
      w_fit[j] = omega_min + (omega_max - omega_min) * static_cast<double>(j) /
                                 static_cast<double>(N - 1);
    }
    const double w_mid = 0.5 * (omega_min + omega_max);
    // Vandermonde V_{j,k} = (w_fit[j])^{2k}, k = 0..D. Solve V·c = kn² per port.
    Eigen::MatrixXd V(N, N);
    for (std::size_t j = 0; j < N; j++)
    {
      double pw = 1.0;
      const double w2 = w_fit[j] * w_fit[j];
      for (std::size_t k = 0; k < N; k++)
      {
        V(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(k)) = pw;
        pw *= w2;
      }
    }
    auto lu = V.partialPivLu();
    for (std::size_t k = 0; k < port_indices_.size(); k++)
    {
      Eigen::VectorXcd kn2(N);
      for (std::size_t j = 0; j < N; j++)
      {
        const std::complex<double> kn =
            space_op_.GetWavePortOp().GetWavePortKnComplex(port_indices_[k], w_fit[j]);
        kn2(static_cast<Eigen::Index>(j)) = kn * kn;
      }
      Eigen::VectorXcd c = lu.solve(kn2);
      std::vector<std::complex<double>> cvec(N);
      for (std::size_t j = 0; j < N; j++)
      {
        cvec[j] = c(static_cast<Eigen::Index>(j));
      }
      kn2_coeffs_.push_back(std::move(cvec));

      // Midpoint validator on the complex truth.
      const std::complex<double> kn_mid_truth =
          space_op_.GetWavePortOp().GetWavePortKnComplex(port_indices_[k], w_mid);
      const std::complex<double> kn_mid_truth_sq = kn_mid_truth * kn_mid_truth;
      // Horner evaluation of the polynomial at ω_mid² (highest degree first).
      const double s_mid = w_mid * w_mid;
      std::complex<double> kn_mid_fit_sq = kn2_coeffs_.back()[N - 1];
      for (std::size_t j = N - 1; j > 0; j--)
      {
        kn_mid_fit_sq = kn_mid_fit_sq * s_mid + kn2_coeffs_.back()[j - 1];
      }
      const double rel_err = std::abs(kn_mid_fit_sq - kn_mid_truth_sq) /
                             std::max(std::abs(kn_mid_truth_sq), 1.0e-30);
      if (rel_err > 0.01)
      {
        Mpi::Print(" Warning: wave port {:d} kₙ²(ω) order-{:d} polynomial fit error "
                   "{:.2e} at midpoint ω={:.3e} (truth kₙ²={:+.3e}{:+.3e}i, fit "
                   "kₙ²={:+.3e}{:+.3e}i) — analytic-continuation BC may be "
                   "inaccurate; raise WavePortFitOrder above {:d} or switch to "
                   "WavePortBCEvaluation=Real\n",
                   port_indices_[k], order, rel_err, w_mid, kn_mid_truth_sq.real(),
                   kn_mid_truth_sq.imag(), kn_mid_fit_sq.real(), kn_mid_fit_sq.imag(),
                   order);
      }
    }
  }

  // Evaluate kₙ_p analytically at complex λ from the kₙ² polynomial fit. λ = i·ω
  // convention; kₙ²(ω) = Σ c_j·ω^{2j} gives kₙ²(λ) = Σ c_j·(-λ²)^j (since ω² = −λ²).
  // Horner evaluation in s = -λ². Optional anchor argument applies a sign flip when
  // Re(kn · conj(anchor)) < 0, keeping the propagating-sheet branch across line
  // search. With a zero anchor the principal sqrt is returned unchanged.
  std::complex<double> KnComplex(int local_port, std::complex<double> lambda,
                                 std::complex<double> anchor = {0.0, 0.0}) const
  {
    const auto &c = kn2_coeffs_[local_port];
    const std::complex<double> s = -lambda * lambda;
    // Horner: P(s) = (((c_D·s + c_{D-1})·s + ...)·s + c_0).
    std::complex<double> kn2 = c.back();
    for (std::size_t j = c.size() - 1; j > 0; j--)
    {
      kn2 = kn2 * s + c[j - 1];
    }
    std::complex<double> kn = std::sqrt(kn2);
    if (anchor != std::complex<double>{0.0, 0.0} && std::real(kn * std::conj(anchor)) < 0.0)
    {
      kn = -kn;
    }
    return kn;
  }

  // Build the three monomial-basis coefficient operators of the seed pencil
  //   A2_seed(λ) = A2_0 + λ·A2_1 + λ²·A2_2
  // by least-squares-fitting a degree-2 polynomial in λ to kₙ_p(λ) per port at
  // `points` (typically 6 — 3 on the imaginary λ axis and 3 displaced into the
  // half-plane). Each port contributes (a₀,p, a₁,p, a₂,p) M_wp_p to the corresponding
  // coefficient operator via the i factor implicit in M_wp_p (= i·M_real). Caller
  // should have called FitKnSq() first to populate the per-port kₙ² coefficients.
  // Returns three operators; missing wave ports yield nullptr.
  struct SeedPencil
  {
    std::unique_ptr<ComplexOperator> A2_0, A2_1, A2_2;
  };
  SeedPencil BuildSeedPencilOffAxis(const std::vector<std::complex<double>> &points) const
  {
    SeedPencil out;
    if (Mwp_p_.empty())
    {
      return out;
    }
    MFEM_VERIFY(kn2_coeffs_.size() == Mwp_p_.size(),
                "BuildSeedPencilOffAxis requires FitKnSq to have been called first.");
    MFEM_VERIFY(points.size() >= 3,
                "BuildSeedPencilOffAxis requires at least 3 sample points.");
    const std::size_t N = points.size();
    // Per-port LSQ fit of a₀ + a₁·λ + a₂·λ² to kₙ_p at each sample point.
    // With M_wp_p = i·M_real_p, the operator coefficient for the j-th monomial is the
    // sum over ports of a_{j,p} · M_wp_p. The implicit i factor is already in M_wp_p.
    std::vector<std::complex<double>> coeff_a0(Mwp_p_.size());
    std::vector<std::complex<double>> coeff_a1(Mwp_p_.size());
    std::vector<std::complex<double>> coeff_a2(Mwp_p_.size());
    Eigen::MatrixXcd V(N, 3);
    for (std::size_t i = 0; i < N; i++)
    {
      V(i, 0) = 1.0;
      V(i, 1) = points[i];
      V(i, 2) = points[i] * points[i];
    }
    auto qr = V.colPivHouseholderQr();
    for (std::size_t k = 0; k < Mwp_p_.size(); k++)
    {
      Eigen::VectorXcd y(N);
      for (std::size_t i = 0; i < N; i++)
      {
        y(i) = KnComplex(static_cast<int>(k), points[i]);
      }
      Eigen::Vector3cd c = qr.solve(y);
      coeff_a0[k] = c(0);
      coeff_a1[k] = c(1);
      coeff_a2[k] = c(2);
    }
    out.A2_0 = BuildScaled(coeff_a0);
    out.A2_1 = BuildScaled(coeff_a1);
    out.A2_2 = BuildScaled(coeff_a2);
    return out;
  }

private:
  std::vector<std::complex<double>> MakeKnCoefficients(double omega) const
  {
    std::vector<std::complex<double>> coeffs(Mwp_p_.size());
    for (std::size_t k = 0; k < port_indices_.size(); k++)
    {
      coeffs[k] = space_op_.GetWavePortOp().GetWavePortKnComplex(port_indices_[k], omega);
    }
    return coeffs;
  }

  std::unique_ptr<ComplexOperator>
  BuildScaled(const std::vector<std::complex<double>> &coeffs) const
  {
    std::vector<const ComplexParOperator *> ops;
    ops.reserve(Mwp_p_.size());
    for (const auto &Mp : Mwp_p_)
    {
      ops.push_back(static_cast<const ComplexParOperator *>(Mp.get()));
    }
    return BuildParSumOperator(coeffs, ops, /*set_essential=*/true);
  }

  SpaceOperator &space_op_;
  std::vector<int> port_indices_;
  std::vector<std::unique_ptr<ComplexOperator>> Mwp_p_;
  // Complex-valued kₙ² polynomial fit coefficients in ω², lowest-degree first:
  // kₙ²(ω) = Σ_{j=0..D} kn2_coeffs_[port][j]·ω^{2j}. Capturing complex coefficients
  // is what lets impedance-walled ports preserve wall-loss attenuation through the
  // closed-form analytic continuation; allowing D > 1 fits sub-quadratic deviation
  // typical of non-TEM ports near cutoff.
  std::vector<std::vector<std::complex<double>>> kn2_coeffs_;
  double fd_rel_;
  // When true, BuildComplex / BuildDComplex re-solve the cross-section EVP at the exact
  // complex frequency ω = -i·λ instead of using the kₙ²(ω) polynomial fit. Set for the
  // HYBRID / SLP complex-λ path; NLEIGS keeps the fit (closed-form FN required).
  bool use_exact_ = false;
};

// Caches the ω-independent boundary curl-curl matrix M_ff for the 2nd-order farfield
// ABC. The full A2 contribution at a complex eigenparameter λ is
//   A2_ff(λ) = -1/(2λ) · M_ff   (analytic continuation of i·(0.5/ω)·M_ff under ω = -i·λ)
// with derivative dA2_ff/dλ = +1/(2λ²) · M_ff. M_ff is stored on the REAL slot of a
// ComplexOperator wrapper so a single complex coefficient applied via
// BuildParSumOperator splits into the (fbr, fbi) halves naturally.
class FarfieldFactor
{
public:
  FarfieldFactor(SpaceOperator &space_op) : space_op_(space_op)
  {
    M_ff_ = space_op_.GetFarfieldExtraBoundaryMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  }

  bool empty() const { return M_ff_ == nullptr; }

  // Real-ω stamping (legacy). +i·(0.5/ω)·M_ff so the contribution lives on the
  // imaginary slot, matching the original FarfieldBoundaryOperator behavior.
  std::unique_ptr<ComplexOperator> Build(double omega) const
  {
    if (!M_ff_)
    {
      return {};
    }
    return BuildScaled(std::complex<double>{0.0, 0.5 / omega});
  }

  // Complex-λ stamping: A2_ff(λ) = -1/(2λ)·M_ff.
  std::unique_ptr<ComplexOperator> BuildComplex(std::complex<double> lambda) const
  {
    if (!M_ff_)
    {
      return {};
    }
    return BuildScaled(-0.5 / lambda);
  }

  // Real-ω derivative: dA2_ff/dω = -i·(0.5/ω²)·M_ff.
  std::unique_ptr<ComplexOperator> BuildJacobianTerm(double omega) const
  {
    if (!M_ff_)
    {
      return {};
    }
    return BuildScaled(std::complex<double>{0.0, -0.5 / (omega * omega)});
  }

  // Complex-λ derivative: dA2_ff/dλ = +1/(2λ²)·M_ff.
  std::unique_ptr<ComplexOperator> BuildDComplex(std::complex<double> lambda) const
  {
    if (!M_ff_)
    {
      return {};
    }
    return BuildScaled(0.5 / (lambda * lambda));
  }

  const ComplexOperator *GetMatrix() const { return M_ff_.get(); }

private:
  std::unique_ptr<ComplexOperator> BuildScaled(std::complex<double> coeff) const
  {
    // Use the runtime std::vector overload (the size-1 array overload would require
    // an explicit-instantiation of BuildParSumOperator<1, ...> that does not exist).
    std::vector<std::complex<double>> coeffs{coeff};
    std::vector<const ComplexParOperator *> ops{
        static_cast<const ComplexParOperator *>(M_ff_.get())};
    return BuildParSumOperator(coeffs, ops, /*set_essential=*/true);
  }

  SpaceOperator &space_op_;
  std::unique_ptr<ComplexOperator> M_ff_;
};

// Caches the ω-independent boundary mass matrices A_σ_g per surface-conductivity
// attribute group. A2_σ(λ) = Σ_g f_g(λ) · A_σ_g where f_g(λ) = i·ω/Z(ω) at ω = -i·λ
// is supplied by SurfaceConductivityOperator::EvaluateScalar. Real-ω stamping evaluates
// f_g at real ω = |Im λ|; complex-λ stamping evaluates at complex λ directly.
class SurfSigmaFactor
{
public:
  SurfSigmaFactor(SpaceOperator &space_op) : space_op_(space_op)
  {
    const auto &sigma_op = space_op_.GetSurfaceConductivityOp();
    for (std::size_t g = 0; g < sigma_op.NumGroups(); g++)
    {
      auto Ag =
          space_op_.GetSurfSigmaBoundaryMassMatrix<ComplexOperator>(g, Operator::DIAG_ZERO);
      if (Ag)
      {
        group_indices_.push_back(g);
        A_sigma_g_.push_back(std::move(Ag));
      }
    }
  }

  bool empty() const { return A_sigma_g_.empty(); }

  std::unique_ptr<ComplexOperator> Build(double omega) const
  {
    if (A_sigma_g_.empty())
    {
      return {};
    }
    return BuildScaled(MakeScalars(std::complex<double>{omega, 0.0}));
  }

  std::unique_ptr<ComplexOperator> BuildComplex(std::complex<double> lambda) const
  {
    if (A_sigma_g_.empty())
    {
      return {};
    }
    // ω = -i·λ
    const std::complex<double> omega = std::complex<double>{0.0, -1.0} * lambda;
    return BuildScaled(MakeScalars(omega));
  }

  // Finite-difference Jacobian on real-ω axis (analytical d/dλ of
  // i·ω/Z(ω) = i·ω·(σ·δ)/(1+i)·F is moderately involved with the finite-thickness
  // F factor; the FD form is correct for either h=0 or h>0 and cheap since EvaluateScalar
  // is closed-form). Step uses fd_rel_·max(|ω|, 1).
  std::unique_ptr<ComplexOperator> BuildJacobianTerm(double omega) const
  {
    if (A_sigma_g_.empty())
    {
      return {};
    }
    const double h = 1.0e-4 * std::max(std::abs(omega), 1.0);
    auto plus = MakeScalars(std::complex<double>{omega + h, 0.0});
    auto minus = MakeScalars(std::complex<double>{omega - h, 0.0});
    std::vector<std::complex<double>> dcoeffs(plus.size());
    for (std::size_t k = 0; k < plus.size(); k++)
    {
      // dA2/dω; A2 was built with i·ω/Z(ω) per group → BuildJacobianTerm is the central
      // FD of that scalar. Multiply by -i to convert to dA2/dλ via dω/dλ = -i, since
      // λ = i·ω and the upstream NLEPS expects d/dλ but the real-ω path stores it as
      // d/dω; conversion happens at the caller for HYBRID's funcDA2DOmega.
      dcoeffs[k] = (plus[k] - minus[k]) / (2.0 * h);
    }
    return BuildScaled(dcoeffs);
  }

  // Analytic complex-λ derivative via central FD on the analytic continuation
  // EvaluateScalar(complex ω). f(λ) = i·ω/Z(ω) at ω = -i·λ has chain-rule factor
  // dω/dλ = -i. Step in λ uses fd_rel_·max(|λ|, 1); EvaluateScalar at the perturbed
  // complex ω is closed-form so the FD is well-conditioned.
  std::unique_ptr<ComplexOperator> BuildDComplex(std::complex<double> lambda) const
  {
    if (A_sigma_g_.empty())
    {
      return {};
    }
    const double h = 1.0e-4 * std::max(std::abs(lambda), 1.0);
    const std::complex<double> lam_plus = lambda + h;
    const std::complex<double> lam_minus = lambda - h;
    auto plus = MakeScalars(std::complex<double>{0.0, -1.0} * lam_plus);
    auto minus = MakeScalars(std::complex<double>{0.0, -1.0} * lam_minus);
    std::vector<std::complex<double>> dcoeffs(plus.size());
    for (std::size_t k = 0; k < plus.size(); k++)
    {
      dcoeffs[k] = (plus[k] - minus[k]) / (2.0 * h);
    }
    return BuildScaled(dcoeffs);
  }

  std::size_t NumGroups() const { return group_indices_.size(); }
  std::size_t GroupIndex(std::size_t k) const { return group_indices_[k]; }
  const ComplexOperator *GetMatrix(std::size_t k) const { return A_sigma_g_[k].get(); }

private:
  std::vector<std::complex<double>> MakeScalars(std::complex<double> omega) const
  {
    const auto &sigma_op = space_op_.GetSurfaceConductivityOp();
    std::vector<std::complex<double>> scalars(group_indices_.size());
    for (std::size_t k = 0; k < group_indices_.size(); k++)
    {
      scalars[k] = sigma_op.EvaluateScalar(group_indices_[k], omega);
    }
    return scalars;
  }

  std::unique_ptr<ComplexOperator>
  BuildScaled(const std::vector<std::complex<double>> &coeffs) const
  {
    std::vector<const ComplexParOperator *> ops;
    ops.reserve(A_sigma_g_.size());
    for (const auto &Ag : A_sigma_g_)
    {
      ops.push_back(static_cast<const ComplexParOperator *>(Ag.get()));
    }
    return BuildParSumOperator(coeffs, ops, /*set_essential=*/true);
  }

  SpaceOperator &space_op_;
  std::vector<std::size_t> group_indices_;
  std::vector<std::unique_ptr<ComplexOperator>> A_sigma_g_;
};

}  // namespace

std::pair<ErrorIndicator, long long int>
EigenSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct and extract the system matrices defining the eigenvalue problem. The diagonal
  // values for the mass matrix PEC dof shift the Dirichlet eigenvalues out of the
  // computational range. The damping matrix may be nullptr.
  BlockTimer bt0(Timer::CONSTRUCT);
  SpaceOperator space_op(iodata, mesh);
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);

  // Per-BC factors: pre-assemble the ω-independent boundary matrices for each
  // nonlinear BC (wave port, 2nd-order farfield ABC, surface conductivity), so the
  // funcA2 hot-path (NLEPS Newton iterates, SLEPc NEP residual evaluations) only
  // re-evaluates the scalar coefficients per λ — never re-runs an FE assembly.
  // Empty when the corresponding BC isn't present in the model.
  WavePortFactor wp_factor(space_op);
  FarfieldFactor ff_factor(space_op);
  SurfSigmaFactor sg_factor(space_op);

  const bool waveport_complex_bc =
      (iodata.solver.eigenmode.waveport_bc_evaluation == WavePortBCEvaluation::COMPLEX);

  // Exact complex-frequency wave-port solve: when the Complex BC is active, the user
  // requested the exact path, and the solver is NOT NLEIGS (which needs a closed-form
  // FN scalar), drive BuildComplex/BuildDComplex to re-solve the cross-section EVP at
  // ω = -i·λ instead of using the kₙ²(ω) polynomial fit. This eliminates the fit error
  // that floors Newton convergence on near-cutoff non-TEM ports.
  const bool use_exact_complex_kn =
      waveport_complex_bc && iodata.solver.eigenmode.waveport_complex_exact &&
      iodata.solver.eigenmode.nonlinear_type != NonlinearEigenSolver::NLEIGS;
  wp_factor.SetUseExactSolve(use_exact_complex_kn);

  // Real-ω stamping (legacy). For each BC, evaluate at real ω = |Im λ|. Wave-port
  // term is added via the factor's cache; bulk farfield/conductivity could go either
  // through SpaceOperator::GetExtraSystemMatrix (legacy assembly) or through the new
  // factors. Routing through the factors when they are non-empty avoids re-running the
  // boundary assembly on every funcA2 invocation; falls back to the legacy bulk
  // assembly otherwise (e.g. for any BC not yet factored).
  auto funcA2_real = [&space_op, &wp_factor, &ff_factor,
                      &sg_factor](double omega) -> std::unique_ptr<ComplexOperator>
  {
    std::vector<std::complex<double>> coeffs;
    std::vector<const ComplexParOperator *> ops;
    auto bulk = space_op.GetExtraSystemMatrix<ComplexOperator>(
        omega, Operator::DIAG_ZERO, /*include_wave_ports=*/false);
    if (bulk)
    {
      coeffs.push_back(1.0 + 0.0i);
      ops.push_back(static_cast<const ComplexParOperator *>(bulk.get()));
    }
    auto wp = wp_factor.Build(omega);
    if (wp)
    {
      coeffs.push_back(1.0 + 0.0i);
      ops.push_back(static_cast<const ComplexParOperator *>(wp.get()));
    }
    if (ops.empty())
    {
      return {};
    }
    if (ops.size() == 1)
    {
      // Single contribution — return it directly to preserve the existing pointer
      // path (BuildParSumOperator with one operand allocates an unnecessary wrapper).
      return bulk ? std::move(bulk) : std::move(wp);
    }
    // Multiple contributions: the sum references the operands' local matrices, so the
    // operands must outlive it. Transfer their ownership into the returned operator.
    auto sum = BuildParSumOperator(coeffs, ops);
    std::vector<std::unique_ptr<ComplexOperator>> operands;
    operands.push_back(std::move(bulk));
    operands.push_back(std::move(wp));
    sum->TakeOwnership(std::move(operands));
    return sum;
  };

  // Complex-λ stamping (analytic continuation). Always installed for the diagnostic
  // comparison printed at SetInitialGuess; promoted to the production residual /
  // Jacobian path when WavePortBCEvaluation = COMPLEX. Each factor exposes
  // BuildComplex(λ) returning Σ (constant matrix) · (scalar holomorphic f_j(λ))
  // contribution for its BC family. Empty contributions are skipped so the result
  // is null when no nonlinear BC is present.
  auto funcA2_complex = [&wp_factor, &ff_factor, &sg_factor](
                            std::complex<double> lambda) -> std::unique_ptr<ComplexOperator>
  {
    std::vector<std::complex<double>> coeffs;
    std::vector<const ComplexParOperator *> ops;
    auto wp = wp_factor.BuildComplex(lambda);
    auto ff = ff_factor.BuildComplex(lambda);
    auto sg = sg_factor.BuildComplex(lambda);
    if (wp)
    {
      coeffs.push_back(1.0 + 0.0i);
      ops.push_back(static_cast<const ComplexParOperator *>(wp.get()));
    }
    if (ff)
    {
      coeffs.push_back(1.0 + 0.0i);
      ops.push_back(static_cast<const ComplexParOperator *>(ff.get()));
    }
    if (sg)
    {
      coeffs.push_back(1.0 + 0.0i);
      ops.push_back(static_cast<const ComplexParOperator *>(sg.get()));
    }
    if (ops.empty())
    {
      return {};
    }
    if (ops.size() == 1)
    {
      return wp ? std::move(wp) : (ff ? std::move(ff) : std::move(sg));
    }
    // Multiple BC contributions (e.g. wave port + 2nd-order ABC): the sum references
    // the operands' local matrices, so keep the operands alive alongside it.
    auto sum = BuildParSumOperator(coeffs, ops);
    std::vector<std::unique_ptr<ComplexOperator>> operands;
    if (wp)
    {
      operands.push_back(std::move(wp));
    }
    if (ff)
    {
      operands.push_back(std::move(ff));
    }
    if (sg)
    {
      operands.push_back(std::move(sg));
    }
    sum->TakeOwnership(std::move(operands));
    return sum;
  };

  // Top-level funcA2 used for the seed pencil and the legacy callers below. Chooses
  // the real-ω or complex-λ path based on the config knob, with a uniform double-ω
  // signature: the complex path interprets the real ω as λ = i·ω (the imaginary λ
  // axis), which is exactly the trace of the complex-λ operator on real frequencies.
  auto funcA2 = [&waveport_complex_bc, &funcA2_real,
                 &funcA2_complex](double omega) -> std::unique_ptr<ComplexOperator>
  {
    if (waveport_complex_bc)
    {
      return funcA2_complex(std::complex<double>{0.0, omega});
    }
    return funcA2_real(omega);
  };
  // Seed-pencil source with the 2nd-order farfield ABC term removed. The ABC
  // contributes −1/(2λ)·M_ff (a simple pole at λ=0). A polynomial seed cannot fit that
  // pole, so the misfit leaks a fictitious, imaginary, rank-deficient a₂·M_ff into the
  // seed's leading (M) coefficient, which displaces the linearization's spurious roots
  // into the +Re(λ) (growing/unphysical) half-plane. We instead remove the ABC from the
  // interpolated source here and re-add it to the seed as a frozen constant in the
  // K-block (see the HYBRID seed assembly below). On the imaginary λ axis (where the
  // interpolation samples) ff_factor.Build(ω) is bit-identical to the farfield part of
  // funcA2(ω), so the subtraction cancels it cleanly in both REAL and COMPLEX modes.
  auto funcA2_seed = [&funcA2, &ff_factor](double omega) -> std::unique_ptr<ComplexOperator>
  {
    auto full = funcA2(omega);
    if (ff_factor.empty() || !full)
    {
      return full;
    }
    auto ff = ff_factor.Build(omega);  // i·(0.5/ω)·M_ff — the farfield part of funcA2(ω)
    if (!ff)
    {
      return full;
    }
    std::vector<std::complex<double>> coeffs{1.0 + 0.0i, -1.0 + 0.0i};
    std::vector<const ComplexParOperator *> ops{
        static_cast<const ComplexParOperator *>(full.get()),
        static_cast<const ComplexParOperator *>(ff.get())};
    auto sum = BuildParSumOperator(coeffs, ops);
    std::vector<std::unique_ptr<ComplexOperator>> operands;
    operands.push_back(std::move(full));
    operands.push_back(std::move(ff));
    sum->TakeOwnership(std::move(operands));
    return sum;
  };
  auto funcP = [&space_op](std::complex<double> a0, std::complex<double> a1,
                           std::complex<double> a2,
                           std::complex<double> omega) -> std::unique_ptr<ComplexOperator>
  { return space_op.GetPreconditionerMatrix<ComplexOperator>(a0, a1, a2, omega); };
  const double target = iodata.solver.eigenmode.target;
  // Populate the wave-port kₙ² fit ahead of any funcA2(target) call. Required when
  // funcA2 routes to funcA2_complex (WavePortBCEvaluation = COMPLEX) since
  // BuildComplex relies on kn2_coeffs_; harmless in REAL mode. Order is the
  // user-configurable polynomial degree in ω² (default 1 = original 2-coefficient
  // fit; raise to capture sub-quadratic kₙ²(ω) deviation near non-TEM cutoffs).
  if (!wp_factor.empty())
  {
    wp_factor.FitKnSq(target, iodata.solver.eigenmode.target_upper,
                      iodata.solver.eigenmode.waveport_fit_order);
  }
  auto A2 = funcA2(target);
  bool has_A2 = (A2 != nullptr);

  // Extend K, C, M operators with interpolated A2 operator.
  // K' = K + A2_0, C' = C + A2_1, M' = M + A2_2
  std::unique_ptr<ComplexOperator> Kp, Cp, Mp;
  std::unique_ptr<Interpolation> interp_op;
  std::unique_ptr<ComplexOperator> A2_0, A2_1, A2_2;
  NonlinearEigenSolver nonlinear_type = iodata.solver.eigenmode.nonlinear_type;
  if (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
  {
    const double target_max = iodata.solver.eigenmode.target_upper;

    // Seed-pencil strategies for HYBRID. Options selectable via env var:
    //   PALACE_NLEPS_SEED unset / "polynomial":
    //     3-point Newton interpolation of A2(ω) in λ. Polynomial extrapolation off the
    //     imaginary λ axis can introduce spurious eigenvalues (eigenvalues of
    //     T_polynomial that don't correspond to T_nonlinear modes).
    //   PALACE_NLEPS_SEED = "off-axis":
    //     6-point LSQ fit (3 on imag axis + 3 mirror points) in λ. Same polynomial
    //     pencil structure; better conditioning for typical Q ~ 5 modes.
    //   PALACE_NLEPS_SEED = "fixed-omega":
    //     Constant A2(ω*=target) ("frozen-frequency" QEP). T_re(λ; ω*) = K + λC +
    //     λ²M + i·kₙ(ω*)·M_p is a strict QEP with no λ-dependence in A2 — matches
    //     the wave-port BC definition exactly at ω* = |Im λ*|. Specialization of
    //     Ruhe's Method of Successive Linear Problems (SIAM JNA 1973); see also
    //     COMSOL "linearization point" iteration (Wu et al., CPC 284, 2023) and
    //     SLAC Omega3P self-consistent iteration (Liao-Bai-Lee-Ko 2010).
    //   PALACE_NLEPS_SEED = "frozen-abc" (DEFAULT when a 2nd-order farfield ABC is
    //   present):
    //     Per-term freeze of ONLY the farfield ABC term. The ABC contributes
    //     −1/(2λ)·M_ff (a simple pole at λ=0); a polynomial seed cannot fit it, and the
    //     misfit injects a fictitious, imaginary, rank-deficient a₂·M_ff into the seed's
    //     leading (M) coefficient — which is what displaces the linearization's spurious
    //     roots into the +Re(λ) half-plane (confirmed empirically: no-ABC and 1st-order
    //     ABC seeds are clean, 2nd-order ABC seeds are not). Freezing the ABC at λ*=i·target
    //     turns it into the bounded constant +i/(2·target)·M_ff, which we place in the
    //     K-block (A2_0) only. The leading coefficient stays the clean bulk mass M, so the
    //     displacement mechanism is removed; the wave-port quadratic fit (the dominant,
    //     well-conditioned part) is untouched in A2_1/A2_2. First-order accurate across
    //     [target, target_upper] — same regime as fixed-omega — but only a seed; the
    //     Quasi-Newton step refines against the exact T(λ) regardless.
    enum class SeedStrategy
    {
      Polynomial,
      OffAxis,
      FixedOmega,
      FrozenABC
    };
    // Default to FrozenABC when a 2nd-order farfield ABC is present (the case the freeze
    // is designed for); otherwise the polynomial seed (with no ABC term) is already clean.
    SeedStrategy seed_strategy =
        ff_factor.empty() ? SeedStrategy::Polynomial : SeedStrategy::FrozenABC;
    if (const char *env = std::getenv("PALACE_NLEPS_SEED"))
    {
      std::string_view sv{env};
      if (sv == "polynomial")
        seed_strategy = SeedStrategy::Polynomial;
      else if (sv == "off-axis")
        seed_strategy = SeedStrategy::OffAxis;
      else if (sv == "fixed-omega")
        seed_strategy = SeedStrategy::FixedOmega;
      else if (sv == "frozen-abc")
        seed_strategy = SeedStrategy::FrozenABC;
    }
    Mpi::Print("\n NLEPS HYBRID seed strategy: {}\n",
               seed_strategy == SeedStrategy::Polynomial ? "polynomial (3-pt on-axis)"
               : seed_strategy == SeedStrategy::OffAxis  ? "off-axis (6-pt LSQ)"
               : seed_strategy == SeedStrategy::FixedOmega
                   ? "fixed-omega (constant A2)"
                   : "frozen-abc (frozen 2nd-order ABC, quadratic wave-port fit)");

    // Detect wave-port-only A2 (used by both off-axis and fixed-omega paths).
    auto bulk_probe = space_op.GetExtraSystemMatrix<ComplexOperator>(
        target, Operator::DIAG_ZERO, /*include_wave_ports=*/false);
    const bool wave_port_only = !bulk_probe && !wp_factor.empty();

    if (seed_strategy == SeedStrategy::FixedOmega)
    {
      // T_re(λ; ω*) = K + λC + λ²M + i·kₙ(ω*)·M_p with kₙ frozen at ω* = target.
      // The wave-port term is constant in λ — put it in the K block (α₀ = full
      // value, α₁ = α₂ = 0).
      A2_0 = funcA2(target);
    }
    else if (seed_strategy == SeedStrategy::OffAxis && wave_port_only)
    {
      wp_factor.FitKnSq(target, target_max, iodata.solver.eigenmode.waveport_fit_order);
      const double displacement = 0.1 * target;
      std::vector<std::complex<double>> pts;
      pts.reserve(6);
      for (int j = 0; j < 3; j++)
      {
        const double w = target + j * (target_max - target) / 2.0;
        pts.push_back(std::complex<double>(0.0, w));
      }
      for (int j = 0; j < 3; j++)
      {
        const double w = target + j * (target_max - target) / 2.0;
        pts.push_back(std::complex<double>(-displacement, w));
      }
      auto seed = wp_factor.BuildSeedPencilOffAxis(pts);
      A2_0 = std::move(seed.A2_0);
      A2_1 = std::move(seed.A2_1);
      A2_2 = std::move(seed.A2_2);
    }
    else if (seed_strategy == SeedStrategy::FrozenABC)
    {
      // Per-term frozen ABC. Interpolate the seed source WITHOUT the farfield term
      // (funcA2_seed = funcA2 − ff_factor.Build), so the wave-port √ keeps its clean
      // quadratic fit distributed through A2_0/A2_1/A2_2 and no fictitious M_ff leaks
      // into A2_2. Then add the ABC back as a frozen constant in the K-block only.
      interp_op = std::make_unique<NewtonInterpolationOperator>(funcA2_seed, A2->Width());
      interp_op->Interpolate(1i * target, 1i * target_max);
      auto A2_0_wp = interp_op->GetInterpolationOperator(0);
      A2_1 = interp_op->GetInterpolationOperator(1);
      A2_2 = interp_op->GetInterpolationOperator(2);
      // Frozen ABC constant: ff_factor.BuildComplex(i·target) = −1/(2·i·target)·M_ff
      //                                                       = +i/(2·target)·M_ff.
      auto A2_0_ff = ff_factor.BuildComplex(1i * target);
      if (A2_0_ff)
      {
        // A2_0 = (wave-port λ⁰ term) + (frozen ABC constant). Both are ComplexParOperator;
        // keep them alive alongside the sum (BuildParSumOperator holds non-owning refs).
        std::vector<std::complex<double>> coeffs{1.0 + 0.0i, 1.0 + 0.0i};
        std::vector<const ComplexParOperator *> ops{
            static_cast<const ComplexParOperator *>(A2_0_wp.get()),
            static_cast<const ComplexParOperator *>(A2_0_ff.get())};
        auto sum = BuildParSumOperator(coeffs, ops);
        std::vector<std::unique_ptr<ComplexOperator>> operands;
        operands.push_back(std::move(A2_0_wp));
        operands.push_back(std::move(A2_0_ff));
        sum->TakeOwnership(std::move(operands));
        A2_0 = std::move(sum);
      }
      else
      {
        A2_0 = std::move(A2_0_wp);
      }
    }
    else
    {
      // Polynomial 3-point Newton interpolation. Distributes the FULL A2(λ) — including
      // any farfield ABC — through K, C, M. Clean when no ABC is present; for the ABC
      // case prefer FrozenABC (the default when ff_factor is non-empty).
      interp_op = std::make_unique<NewtonInterpolationOperator>(funcA2, A2->Width());
      interp_op->Interpolate(1i * target, 1i * target_max);
      A2_0 = interp_op->GetInterpolationOperator(0);
      A2_1 = interp_op->GetInterpolationOperator(1);
      A2_2 = interp_op->GetInterpolationOperator(2);
    }
    // Assemble the seed pencil. Polynomial paths distribute the wave-port
    // contribution through K, C, M (A2_1 and A2_2 always non-null). Fixed-ω
    // puts everything in K (A2_1 = A2_2 = nullptr); Cp and Mp wrap C and M.
    // When fixed-ω is used and input C is null (no damping), Cp is a 0·K
    // wrapper to give the downstream PEP solver a valid (zero) damping
    // operator with the right FE space, since SetOperators(K,C,M) is the
    // only QEP path supported.
    Kp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {K.get(), A2_0.get()});
    if (A2_1)
    {
      Cp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {C.get(), A2_1.get()});
    }
    else if (C)
    {
      Cp = BuildParSumOperator({1.0 + 0i}, {C.get()});
    }
    else
    {
      Cp = BuildParSumOperator({0.0 + 0i}, {K.get()});
    }
    if (A2_2)
    {
      Mp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {M.get(), A2_2.get()});
    }
    else
    {
      Mp = BuildParSumOperator({1.0 + 0i}, {M.get()});
    }

    // Per-port kₙ fit-residual diagnostic. The HYBRID seed pencil is built from a 3-point
    // monomial-basis Newton interpolation of A2(ω) on [target, target_max]; for wave-port
    // dispersion that's a degree-2 polynomial fit on kₙ,p(ω). Sample on a dense grid,
    // compare to the truth, emit the per-port relative residual. A large residual
    // (>~1e-2) flags that the linear seed will be poor and the Quasi-Newton refinement
    // will need many iterations / aggressive Armijo backtracking from each seed.
    if (!wp_factor.empty())
    {
      // Reproduce the 3-point Newton interpolation: the seed pencil's kₙ,p is the
      // unique quadratic in ω matching the true kₙ,p at three evenly spaced points.
      // The fit is exact at those three anchors, so the fit error is only observable
      // OFF the anchors. A few off-anchor probes (the two anchor-midpoints, plus a
      // point just beyond the band) capture the residual just as well as a dense grid
      // — and each probe is a full cross-section EVP solve (GetWavePortKn →
      // WavePortData::Initialize), so a dense grid (formerly 500 points × NumPorts) was
      // a large, pure-diagnostic cost. Keep this cheap.
      std::vector<double> w_fit(3);
      for (int j = 0; j < 3; j++)
      {
        w_fit[j] = target + j * (target_max - target) / 2.0;
      }
      // Off-anchor probe points for the residual estimate. The fit is exact at the
      // anchors (fractions 0, 0.5, 1.0), so probes must avoid those. Use the two
      // anchor-midpoints (0.25, 0.75) where interpolation error peaks.
      const std::vector<double> w_probe = {
          target + 0.25 * (target_max - target),  // between anchors 0 and 1
          target + 0.75 * (target_max - target)   // between anchors 1 and 2
      };
      Eigen::Matrix3d V_fit;
      for (int j = 0; j < 3; j++)
      {
        V_fit(j, 0) = 1.0;
        V_fit(j, 1) = w_fit[j];
        V_fit(j, 2) = w_fit[j] * w_fit[j];
      }
      auto qr = V_fit.partialPivLu();
      for (std::size_t p = 0; p < wp_factor.NumPorts(); p++)
      {
        Eigen::Vector3d kn_fit;
        for (int j = 0; j < 3; j++)
        {
          kn_fit(j) = wp_factor.Kn(static_cast<int>(p), w_fit[j]);
        }
        const Eigen::Vector3d c = qr.solve(kn_fit);
        double max_abs = 0.0, max_rel_num = 0.0;
        for (double w : w_probe)
        {
          const double truth = wp_factor.Kn(static_cast<int>(p), w);
          const double poly = c(0) + c(1) * w + c(2) * w * w;
          max_abs = std::max(max_abs, std::abs(truth));
          max_rel_num = std::max(max_rel_num, std::abs(poly - truth));
        }
        const double rel_err = (max_abs > 0.0) ? max_rel_num / max_abs : 0.0;
        const double f_target =
            iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target) / (2.0 * M_PI);
        const double f_target_max =
            iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target_max) /
            (2.0 * M_PI);
        if (rel_err > 1.0e-2)
        {
          Mpi::Warning("Wave port {:d}: HYBRID linear seed will be poor (order-2 kₙ fit "
                       "residual {:.3e} on [{:.3e}, {:.3e}] GHz exceeds 1%). Quasi-Newton "
                       "refinement may take many iterations or fail to converge.\n",
                       p, rel_err, f_target, f_target_max);
        }
        else
        {
          Mpi::Print(" Wave port {:d} (NLEPS HYBRID seed): order-2 kₙ fit residual {:.3e} "
                     "on [{:.3e}, {:.3e}] GHz\n",
                     p, rel_err, f_target, f_target_max);
        }
      }
    }
  }

  const auto &Curl = space_op.GetCurlMatrix();
  SaveMetadata(space_op.GetNDSpaces());

  // Configure objects for postprocessing.
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);
  ComplexVector E(Curl.Width()), B(Curl.Height());
  E.UseDevice(true);
  B.UseDevice(true);

  // Define and configure the eigensolver to solve the eigenvalue problem:
  //         (K + λ C + λ² M) u = 0    or    K u = -λ² M u
  // with λ = iω. In general, the system matrices are complex and symmetric.
  std::unique_ptr<EigenvalueSolver> eigen;
  const EigenSolverBackend type = iodata.solver.eigenmode.type;
#if !defined(PALACE_WITH_ARPACK) && !defined(PALACE_WITH_SLEPC)
#error "Eigenmode solver requires building with ARPACK or SLEPc!"
#endif
#if !defined(PALACE_WITH_SLEPC)
  if (nonlinear_type == NonlinearEigenSolver::SLP)
  {
    Mpi::Warning("SLP nonlinear eigensolver not available without SLEPc, using Hybrid!\n");
  }
  if (nonlinear_type == NonlinearEigenSolver::NLEIGS)
  {
    Mpi::Warning("NLEIGS nonlinear eigensolver not available with ARPACK, using Hybrid!\n");
  }
  nonlinear_type = NonlinearEigenSolver::HYBRID;
#endif
  if (type == EigenSolverBackend::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    Mpi::Print("\nConfiguring ARPACK eigenvalue solver:\n");
    if (C || has_A2)
    {
      eigen = std::make_unique<arpack::ArpackPEPSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
    }
    else
    {
      eigen = std::make_unique<arpack::ArpackEPSSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
    }
#endif
  }
  else  // EigenSolverBackend::SLEPC
  {
#if defined(PALACE_WITH_SLEPC)
    Mpi::Print("\nConfiguring SLEPc eigenvalue solver:\n");
    std::unique_ptr<slepc::SlepcEigenvalueSolver> slepc;
    if (nonlinear_type == NonlinearEigenSolver::SLP)
    {
      slepc = std::make_unique<slepc::SlepcNEPSolver>(space_op.GetComm(),
                                                      iodata.problem.verbose);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::SLP);
      slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GENERAL);
    }
    else if (nonlinear_type == NonlinearEigenSolver::NLEIGS)
    {
      MFEM_VERIFY(waveport_complex_bc,
                  "NLEIGS requires WavePortBCEvaluation = \"Complex\": the nonlinear "
                  "operator T(λ) must be holomorphic for rational interpolation. "
                  "Either set the BC to Complex or pick a different NonlinearType.");
      MFEM_VERIFY(sg_factor.empty(),
                  "NLEIGS split-form does not yet support surface-conductivity BCs "
                  "(skin-depth Z(ω) ∝ √ω is not natively expressible as a SLEPc FN). "
                  "Use HYBRID or SLP for surf-σ models, or remove the conductivity "
                  "boundaries to use NLEIGS.");
      auto nleigs = std::make_unique<slepc::SlepcNEPNLEIGSSolver>(space_op.GetComm(),
                                                                  iodata.problem.verbose);
      nleigs->SetType(slepc::SlepcEigenvalueSolver::Type::NLEIGS);
      nleigs->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GENERAL);
      nleigs->SetNLEIGSFullBasis(iodata.solver.eigenmode.nleigs_full_basis);
      nleigs->SetNLEIGSInterpolation(iodata.solver.eigenmode.nleigs_interp_tol,
                                     iodata.solver.eigenmode.nleigs_interp_deg);
      nleigs->SetTargetUpper(iodata.solver.eigenmode.target_upper);
      nleigs->SetNLEIGSRKShifts(iodata.solver.eigenmode.nleigs_rk_shifts);
      if (!iodata.solver.eigenmode.nleigs_region.empty())
      {
        nleigs->SetNLEIGSRegion(iodata.solver.eigenmode.nleigs_region);
      }
      if (!iodata.solver.eigenmode.nleigs_singularities.empty())
      {
        // Flat array of (Re, Im) pairs → vector<complex>.
        const auto &flat = iodata.solver.eigenmode.nleigs_singularities;
        std::vector<std::complex<double>> xi(flat.size() / 2);
        for (std::size_t i = 0; i < xi.size(); i++)
        {
          xi[i] = {flat[2 * i], flat[2 * i + 1]};
        }
        nleigs->SetNLEIGSSingularities(xi);
      }
      else if (!wp_factor.empty() &&
               iodata.solver.eigenmode.nleigs_singularities_per_cut > 0)
      {
        // Auto-derive ξ from the per-port wave-port kₙ branch cuts. The split-form
        // term kₙ_p(λ) = √(α_p − γ_p·λ²) is non-analytic on the branch cut connecting
        // its branch points λ = ±√(α_p/γ_p). For a propagating (above-cutoff) port
        // α_p/γ_p < 0, so the branch points are pure-imaginary at ±i·t_bp with
        // t_bp = √(−α_p/γ_p) = ω_cutoff (nondim); the cut is the imaginary-axis
        // segment λ = i·t, |t| < t_bp. SLEPc NLEIGS approximates T(λ) with a rational
        // interpolant whose POLES are chosen (Leja–Bagby) from this ξ set; placing
        // poles on the cut, clustered toward the branch point (the singularity
        // nearest the eigenvalue band), lets the interpolant resolve the √. Units are
        // raw λ = i·ω_nd (NEP runs unscaled, γ = δ = 1) — same as the FN coefficients
        // and the target. Only the upper-half branch point matters (eigenvalues sit
        // just above i·t_bp). Keep the count small (default 1 = the branch point
        // itself): too many poles just below the RG lower bound (= target) pull
        // spurious "converged" modes onto the boundary.
        const int n_per_cut = iodata.solver.eigenmode.nleigs_singularities_per_cut;
        std::vector<std::complex<double>> xi;
        for (std::size_t p = 0; p < wp_factor.NumPorts(); p++)
        {
          const std::complex<double> bp =
              std::sqrt(wp_factor.Alpha(p) / wp_factor.Gamma(p));  // branch point ±bp
          const double t_bp = std::abs(bp.imag());
          if (t_bp < 1.0e-3)
          {
            continue;  // TEM-like port: branch at origin, no usable imaginary cut
          }
          // Geometric clustering toward t_bp on the upper cut, excluding the origin.
          // j = 0 → branch point itself (strongest singularity).
          for (int j = 0; j < n_per_cut; j++)
          {
            const double t = t_bp * std::pow(0.5, j);
            if (t < 0.5)
            {
              break;  // stay clear of the SLEPc near-zero singularity guard
            }
            xi.emplace_back(0.0, t);
          }
        }
        // If no branch-cut singularities were derived (e.g. quasi-TEM ports whose
        // cutoff ≈ 0, so the branch point sits at the origin and is unusable), we must
        // still install SOMETHING when a farfield 2nd-order ABC is active: that BC
        // contributes a genuine pole of T(λ) at λ = 0, and SLEPc's AAA auto-discovery
        // would find it and then abort in Leja–Bagby ("singularity nearly zero"). A far
        // sentinel singularity above the band gives Leja–Bagby a valid pole set and
        // keeps it from sampling the origin pole. The λ=0 pole itself is harmless to
        // the interpolant because the eigenvalue band is bounded away from 0.
        if (xi.empty() && !ff_factor.empty())
        {
          const double target_im = target;  // nondim Im(λ) of the target (λ = i·ω)
          const double upper_im = (iodata.solver.eigenmode.target_upper > 0.0)
                                      ? iodata.solver.eigenmode.target_upper
                                      : 3.0 * target_im;
          xi.emplace_back(0.0, 100.0 * upper_im);
        }
        if (!xi.empty())
        {
          nleigs->SetNLEIGSSingularities(xi);
          Mpi::Print(" NLEIGS: auto-derived {:d} singularit{} (wave-port branch cuts / "
                     "farfield-pole sentinel)\n",
                     xi.size(), xi.size() == 1 ? "y" : "ies");
        }
      }
      // Install the same preconditioner builder SLP uses; NLEIGS' inner KSPs
      // (one per rational interpolation shift) wrap it in a PCSHELL.
      nleigs->SetPreconditionerBuilder(funcP);
      slepc = std::move(nleigs);
    }
    else
    {
      if (C || has_A2)
      {
        if (!iodata.solver.eigenmode.pep_linear)
        {
          slepc = std::make_unique<slepc::SlepcPEPSolver>(space_op.GetComm(),
                                                          iodata.problem.verbose);
          slepc->SetType(slepc::SlepcEigenvalueSolver::Type::TOAR);
        }
        else
        {
          slepc = std::make_unique<slepc::SlepcPEPLinearSolver>(space_op.GetComm(),
                                                                iodata.problem.verbose);
          slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
        }
      }
      else
      {
        slepc = std::make_unique<slepc::SlepcEPSSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
        slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
      }
      slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    }
    slepc->SetOrthogonalization(iodata.solver.linear.gs_orthog == Orthogonalization::MGS,
                                iodata.solver.linear.gs_orthog == Orthogonalization::CGS2);
    eigen = std::move(slepc);
#endif
  }
  EigenvalueSolver::ScaleType scale = iodata.solver.eigenmode.scale
                                          ? EigenvalueSolver::ScaleType::NORM_2
                                          : EigenvalueSolver::ScaleType::NONE;
  // Lifetime extension for NLEIGS' zero-C placeholder. Declared at function scope so
  // it outlives the eigen->Solve() call below (the solver holds a non-owning pointer
  // into it via SlepcNEPNLEIGSSolver::opC_ref).
  std::unique_ptr<ComplexOperator> C_zero_for_nleigs;
  if (nonlinear_type == NonlinearEigenSolver::SLP)
  {
    eigen->SetOperators(*K, *C, *M, EigenvalueSolver::ScaleType::NONE);
    eigen->SetExtraSystemMatrix(funcA2);
    eigen->SetExtraSystemMatrixComplex(funcA2_complex);
    eigen->SetUseComplexA2(waveport_complex_bc);
    eigen->SetPreconditionerUpdate(funcP);
  }
  else if (nonlinear_type == NonlinearEigenSolver::NLEIGS)
  {
#if defined(PALACE_WITH_SLEPC)
    auto *nleigs = dynamic_cast<slepc::SlepcNEPNLEIGSSolver *>(eigen.get());
    MFEM_VERIFY(nleigs, "NLEIGS path requires SLEPc backend.");
    // NLEIGS split form requires a non-null C. If Palace has no damping (no
    // 1st-order farfield, no surface conductivity, no surface impedance, no lumped
    // resistive port — typical for the adapter test case), C is a null unique_ptr.
    // Build a 0·K placeholder so the polynomial λC slot has a valid operator.
    // Mirrors the approach used by the HYBRID seed-pencil path when C is null.
    const ComplexOperator *C_for_split = C.get();
    if (!C_for_split)
    {
      std::vector<std::complex<double>> coeffs{0.0 + 0.0i};
      std::vector<const ComplexParOperator *> ops{
          static_cast<const ComplexParOperator *>(K.get())};
      C_zero_for_nleigs = BuildParSumOperator(coeffs, ops);
      C_for_split = C_zero_for_nleigs.get();
    }
    // SetOperators registers the polynomial pencil terms (K, λC, λ²M) via
    // AddSplitFormTerm. Subsequent AddSplitFormTerm calls append the BC factor terms
    // for wave ports and (optionally) the 2nd-order ABC.
    nleigs->SetOperators(*K, *C_for_split, *M, EigenvalueSolver::ScaleType::NONE);
    // Per-port wave-port term: T_wp_p(λ) = M_wp_p · kn_p(λ) with
    // kn_p(λ) = √(Σ_{j=0..D} c_j · (-λ²)^j). Build the FN as
    // FN_COMBINE_COMPOSE(FN_RATIONAL([a_{2D}, 0, a_{2D-2}, ..., 0, a_0]), FN_SQRT)
    // with a_{2j} = c_j · (-1)^j. SLEPc stores coefficients highest-degree-first.
    // Default fit order D=1 reduces to the textbook two-term [-γ, 0, α].
    for (std::size_t p = 0; p < wp_factor.NumPorts(); p++)
    {
      const auto &c = wp_factor.Kn2Coeffs(p);
      const std::size_t D = c.size() - 1;      // poly degree in ω²
      const std::size_t poly_len = 2 * D + 1;  // poly degree in λ is 2D
      std::vector<PetscScalar> pcoeff(poly_len, 0.0 + 0.0i);
      // kn²(λ) = Σ_{k=0..D} c_k·(-λ²)^k = Σ_k c_k·(-1)^k·λ^{2k}, so the coefficient of
      // λ^{2k} is a_{2k} = (-1)^k·c_k. SLEPc FNRationalSetNumerator stores coefficients
      // HIGHEST-degree-first, so a_{2k} (degree 2k) goes at index (2D − 2k). Odd-power
      // slots stay zero. Default D=1 ⇒ pcoeff = [-c_1, 0, c_0] = [-γ, 0, α], the
      // textbook kn² = α − γλ².
      for (std::size_t k = 0; k <= D; k++)
      {
        const std::size_t idx = 2 * (D - k);  // highest-degree-first index of λ^{2k}
        const double sign = (k % 2 == 0) ? 1.0 : -1.0;
        pcoeff[idx] = sign * c[k];
      }
      FN h, sq, fn_kn;
      PalacePetscCall(FNCreate(space_op.GetComm(), &h));
      PalacePetscCall(FNSetType(h, FNRATIONAL));
      PalacePetscCall(
          FNRationalSetNumerator(h, static_cast<PetscInt>(poly_len), pcoeff.data()));
      PalacePetscCall(FNCreate(space_op.GetComm(), &sq));
      PalacePetscCall(FNSetType(sq, FNSQRT));
      PalacePetscCall(FNCreate(space_op.GetComm(), &fn_kn));
      PalacePetscCall(FNSetType(fn_kn, FNCOMBINE));
      PalacePetscCall(FNCombineSetChildren(fn_kn, FN_COMBINE_COMPOSE, h, sq));
      // The constant matrix Mwp_p_ is stored with the curl-curl coefficient on the
      // imag slot of its ComplexOperator (i·M_real). Multiplying by the bare scalar
      // kn(λ) reproduces the BC stamping +i·kn·M_real.
      nleigs->AddSplitFormTerm(wp_factor.GetMatrix(p), fn_kn);
      // h and sq are children of fn_kn — SLEPc takes references; release our
      // explicit references so destruction follows the parent FN.
      PalacePetscCall(FNDestroy(&h));
      PalacePetscCall(FNDestroy(&sq));
    }
    // Farfield 2nd-order ABC term: f_ff(λ) = -0.5/λ. Constant matrix M_ff is on the
    // real slot of its ComplexOperator wrapper, so the FN scalar carries the full
    // complex coefficient. Build as FN_RATIONAL with numerator [-0.5] and
    // denominator [1, 0] (Horner convention: λ).
    if (!ff_factor.empty())
    {
      FN fn_ff;
      PalacePetscCall(FNCreate(space_op.GetComm(), &fn_ff));
      PalacePetscCall(FNSetType(fn_ff, FNRATIONAL));
      PetscScalar num_ff[1] = {-0.5 + 0.0i};
      PetscScalar den_ff[2] = {1.0 + 0.0i, 0.0 + 0.0i};
      PalacePetscCall(FNRationalSetNumerator(fn_ff, 1, num_ff));
      PalacePetscCall(FNRationalSetDenominator(fn_ff, 2, den_ff));
      nleigs->AddSplitFormTerm(ff_factor.GetMatrix(), fn_ff);
    }
#else
    MFEM_ABORT("NLEIGS requires SLEPc.");
#endif
  }
  else
  {
    if (has_A2)
    {
      eigen->SetOperators(*Kp, *Cp, *Mp, scale);
    }
    else if (C)
    {
      eigen->SetOperators(*K, *C, *M, scale);
    }
    else
    {
      eigen->SetOperators(*K, *M, scale);
    }
  }
  eigen->SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
  const double tol = (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
                         ? iodata.solver.eigenmode.linear_tol
                         : iodata.solver.eigenmode.tol;
  eigen->SetTol(tol);
  eigen->SetMaxIter(iodata.solver.eigenmode.max_it);
  Mpi::Print(" Scaling γ = {:.3e}, δ = {:.3e}\n", eigen->GetScalingGamma(),
             eigen->GetScalingDelta());

  // If desired, use an M-inner product for orthogonalizing the eigenvalue subspace. The
  // constructed matrix just references the real SPD part of the mass matrix (no copy is
  // performed). Boundary conditions don't need to be eliminated here.
  std::unique_ptr<Operator> KM;
  if (iodata.solver.eigenmode.mass_orthog)
  {
    Mpi::Print(" Basis uses M-inner product\n");
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);

    // Mpi::Print(" Basis uses (K + M)-inner product\n");
    // KM = space_op.GetInnerProductMatrix(1.0, 1.0, K.get(), M.get());
    // eigen->SetBMat(*KM);
  }

  // Construct a divergence-free projector so the eigenvalue solve is performed in the space
  // orthogonal to the zero eigenvalues of the stiffness matrix.
  std::unique_ptr<DivFreeSolver<ComplexVector>> divfree;
  if (iodata.solver.linear.divfree_max_it > 0 &&
      !space_op.GetMaterialOp().HasWaveVector() &&
      !space_op.GetMaterialOp().HasLondonDepth())
  {
    Mpi::Print(" Configuring divergence-free projection\n");
    constexpr int divfree_verbose = 0;
    divfree = std::make_unique<DivFreeSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetH1Spaces(),
        space_op.GetAuxBdrTDofLists(), iodata.solver.linear.divfree_tol,
        iodata.solver.linear.divfree_max_it, divfree_verbose);
    eigen->SetDivFreeProjector(*divfree);
  }

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
  }

  // Set up the initial space for the eigenvalue solve. Satisfies boundary conditions and is
  // projected appropriately.
  if (iodata.solver.eigenmode.init_v0)
  {
    ComplexVector v0;
    if (iodata.solver.eigenmode.init_v0_const)
    {
      Mpi::Print(" Using constant starting vector\n");
      space_op.GetConstantInitialVector(v0);
    }
    else
    {
      Mpi::Print(" Using random starting vector\n");
      space_op.GetRandomInitialVector(v0);
    }
    if (divfree)
    {
      divfree->Mult(v0);
    }
    eigen->SetInitialSpace(v0);  // Copies the vector

    // Debug
    // const auto &Grad = space_op.GetGradMatrix();
    // ComplexVector r0(Grad->Width());
    // r0.UseDevice(true);
    // Grad.MultTranspose(v0.Real(), r0.Real());
    // Grad.MultTranspose(v0.Imag(), r0.Imag());
    // r0.Print();
  }

  // Configure the shift-and-invert strategy is employed to solve for the eigenvalues
  // closest to the specified target, σ.
  {
    const double f_target =
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target) / (2 * M_PI);
    Mpi::Print(" Shift-and-invert σ = {:.3e} GHz ({:.3e})\n", f_target, target);
  }
  if (C || has_A2 || nonlinear_type == NonlinearEigenSolver::SLP ||
      nonlinear_type == NonlinearEigenSolver::NLEIGS)
  {
    // Search for eigenvalues closest to λ = iσ.
    eigen->SetShiftInvert(1i * target);
    if (type == EigenSolverBackend::ARPACK)
    {
      // ARPACK searches based on eigenvalues of the transformed problem. The eigenvalue
      // 1 / (λ - σ) will be a large-magnitude negative imaginary number for an eigenvalue
      // λ with frequency close to but not below the target σ.
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::SMALLEST_IMAGINARY);
    }
    else if (nonlinear_type == NonlinearEigenSolver::SLP)
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_MAGNITUDE);
    }
    else if (nonlinear_type == NonlinearEigenSolver::NLEIGS)
    {
      // NLEIGS works on a region (RG), not a target ordering — TARGET_IMAGINARY is
      // the natural choice for eigenvalues clustered near the imaginary λ axis.
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_IMAGINARY);
    }
    else
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_IMAGINARY);
    }
  }
  else
  {
    // Linear EVP has eigenvalues μ = -λ² = ω². Search for eigenvalues closest to μ = σ².
    eigen->SetShiftInvert(target * target);
    if (type == EigenSolverBackend::ARPACK)
    {
      // ARPACK searches based on eigenvalues of the transformed problem. 1 / (μ - σ²)
      // will be a large-magnitude positive real number for an eigenvalue μ with frequency
      // close to but below the target σ².
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
    }
    else
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_REAL);
    }
  }

  // Set up the linear solver required for solving systems involving the shifted operator
  // (K - σ² M) or P(iσ) = (K + iσ C - σ² M) during the eigenvalue solve. The
  // preconditioner for complex linear systems is constructed from a real approximation
  // to the complex system matrix.
  auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * target, -target * target + 0.0i,
                                    K.get(), C.get(), M.get(), A2.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
      1.0 + 0.0i, 1i * target, -target * target + 0.0i, target + 0.0i);
  auto ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                                &space_op.GetH1Spaces());
  ksp->SetOperators(*A, *P);
  eigen->SetLinearSolver(*ksp);

  // Initialize structures for storing and reducing the results of error estimation.
  const bool is_2d = (space_op.GetNDSpace().Dimension() < 3);
  std::unique_ptr<TimeDependentFluxErrorEstimator<ComplexVector>> estimator_3d;
  std::unique_ptr<BoundaryModeFluxErrorEstimator<ComplexVector>> estimator_2d;
  if (is_2d)
  {
    estimator_2d = std::make_unique<BoundaryModeFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        space_op.GetCurlSpace(), space_op.GetH1Spaces(), iodata.solver.linear.estimator_tol,
        iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg);
  }
  else
  {
    estimator_3d = std::make_unique<TimeDependentFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
        iodata.solver.linear.estimator_mg);
  }
  auto AddEstimate =
      [&](const ComplexVector &E, const ComplexVector &B, double Et, ErrorIndicator &ind)
  {
    if (is_2d)
    {
      estimator_2d->AddErrorIndicator(E, B, Et, ind);
    }
    else
    {
      estimator_3d->AddErrorIndicator(E, B, Et, ind);
    }
  };
  ErrorIndicator indicator;

  // Eigenvalue problem solve.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\n");
  int num_conv = eigen->Solve();
  {
    std::complex<double> lambda = (num_conv > 0) ? eigen->GetEigenvalue(0) : 0.0;
    Mpi::Print(" Found {:d} converged eigenvalue{}{}\n", num_conv,
               (num_conv > 1) ? "s" : "",
               (num_conv > 0)
                   ? fmt::format(" (first = {:.3e}{:+.3e}i)", lambda.real(), lambda.imag())
                   : "");
  }

  if (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
  {
    Mpi::Print("\n Refining eigenvalues with Quasi-Newton solver\n");
    auto qn = std::make_unique<QuasiNewtonSolver>(space_op.GetComm(), std::move(eigen),
                                                  num_conv, iodata.problem.verbose,
                                                  iodata.solver.eigenmode.refine_nonlinear);
    qn->SetTol(iodata.solver.eigenmode.tol);
    qn->SetMaxIter(iodata.solver.eigenmode.max_it);
    if (C)
    {
      qn->SetOperators(*K, *C, *M, EigenvalueSolver::ScaleType::NONE);
    }
    else
    {
      qn->SetOperators(*K, *M, EigenvalueSolver::ScaleType::NONE);
    }
    // Install the residual builder. funcA2 is the top-level switch (real-ω vs
    // complex-λ) that returns the SAME operator either way when called on the
    // imaginary λ axis — so we install funcA2 (the switched form) as the legacy
    // double-ω callback.
    qn->SetExtraSystemMatrix(funcA2);
    // Install the complex-λ A2 builder unconditionally so the diagnostic at
    // SetInitialGuess can always print abs(T_analytic) for comparison. The flag
    // SetUseComplexA2 (set just below) toggles whether the production residual /
    // Jacobian path uses it; in REAL mode it does not, preserving baselines.
    qn->SetExtraSystemMatrixComplex(funcA2_complex);
    qn->SetUseComplexA2(waveport_complex_bc);

    // Analytic Jacobians. The complex-λ Jacobian is the sum of per-factor
    // BuildDComplex terms; the real-ω Jacobian uses BuildJacobianTerm(double ω).
    auto funcDA2DOmega = [&wp_factor, &ff_factor,
                          &sg_factor](double omega) -> std::unique_ptr<ComplexOperator>
    {
      std::vector<std::complex<double>> coeffs;
      std::vector<const ComplexParOperator *> ops;
      auto wp = wp_factor.BuildJacobianTerm(omega);
      auto ff = ff_factor.BuildJacobianTerm(omega);
      auto sg = sg_factor.BuildJacobianTerm(omega);
      if (wp)
      {
        coeffs.push_back(1.0 + 0.0i);
        ops.push_back(static_cast<const ComplexParOperator *>(wp.get()));
      }
      if (ff)
      {
        coeffs.push_back(1.0 + 0.0i);
        ops.push_back(static_cast<const ComplexParOperator *>(ff.get()));
      }
      if (sg)
      {
        coeffs.push_back(1.0 + 0.0i);
        ops.push_back(static_cast<const ComplexParOperator *>(sg.get()));
      }
      if (ops.empty())
        return {};
      if (ops.size() == 1)
        return wp ? std::move(wp) : (ff ? std::move(ff) : std::move(sg));
      // Multiple Jacobian-term contributions (e.g. wave port + 2nd-order ABC): the sum
      // references the operands' local matrices, so the operands must outlive it.
      auto sum = BuildParSumOperator(coeffs, ops);
      std::vector<std::unique_ptr<ComplexOperator>> operands;
      if (wp)
      {
        operands.push_back(std::move(wp));
      }
      if (ff)
      {
        operands.push_back(std::move(ff));
      }
      if (sg)
      {
        operands.push_back(std::move(sg));
      }
      sum->TakeOwnership(std::move(operands));
      return sum;
    };
    auto funcDA2DLambdaComplex =
        [&wp_factor, &ff_factor,
         &sg_factor](std::complex<double> lambda) -> std::unique_ptr<ComplexOperator>
    {
      std::vector<std::complex<double>> coeffs;
      std::vector<const ComplexParOperator *> ops;
      auto wp = wp_factor.BuildDComplex(lambda);
      auto ff = ff_factor.BuildDComplex(lambda);
      auto sg = sg_factor.BuildDComplex(lambda);
      if (wp)
      {
        coeffs.push_back(1.0 + 0.0i);
        ops.push_back(static_cast<const ComplexParOperator *>(wp.get()));
      }
      if (ff)
      {
        coeffs.push_back(1.0 + 0.0i);
        ops.push_back(static_cast<const ComplexParOperator *>(ff.get()));
      }
      if (sg)
      {
        coeffs.push_back(1.0 + 0.0i);
        ops.push_back(static_cast<const ComplexParOperator *>(sg.get()));
      }
      if (ops.empty())
        return {};
      if (ops.size() == 1)
        return wp ? std::move(wp) : (ff ? std::move(ff) : std::move(sg));
      // Multiple Jacobian-term contributions (e.g. wave port + 2nd-order ABC): the sum
      // references the operands' local matrices, so the operands must outlive it.
      auto sum = BuildParSumOperator(coeffs, ops);
      std::vector<std::unique_ptr<ComplexOperator>> operands;
      if (wp)
      {
        operands.push_back(std::move(wp));
      }
      if (ff)
      {
        operands.push_back(std::move(ff));
      }
      if (sg)
      {
        operands.push_back(std::move(sg));
      }
      sum->TakeOwnership(std::move(operands));
      return sum;
    };
    const bool any_factor_active =
        !wp_factor.empty() || !ff_factor.empty() || !sg_factor.empty();
    if (any_factor_active)
    {
      qn->SetExtraSystemMatrixDerivative(funcDA2DOmega);
      // Install the complex-λ Jacobian unconditionally; it's only consulted by the
      // solver when use_complex_a2 is true, and being installed lets us swap
      // between modes without re-installing.
      qn->SetExtraSystemMatrixDerivativeComplex(funcDA2DLambdaComplex);
    }
    qn->SetPreconditionerUpdate(funcP);
    qn->SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
    qn->SetPreconditionerLag(iodata.solver.eigenmode.preconditioner_lag,
                             iodata.solver.eigenmode.preconditioner_lag_tol);
    qn->SetMaxRestart(iodata.solver.eigenmode.max_restart);
    qn->SetLinearSolver(*ksp);
    qn->SetShiftInvert(1i * target);
    eigen = std::move(qn);

    // Suppress wave port output during nonlinear eigensolver iterations.
    space_op.GetWavePortOp().SetSuppressOutput(true);
    num_conv = eigen->Solve();
    space_op.GetWavePortOp().SetSuppressOutput(false);
  }

  BlockTimer bt2(Timer::POSTPRO);
  SaveMetadata(*ksp);

  // Calculate and record the error indicators, and postprocess the results.
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n");
  if (!KM)
  {
    // Normalize the finalized eigenvectors with respect to mass matrix (unit electric field
    // energy) even if they are not computed to be orthogonal with respect to it.
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);
    eigen->RescaleEigenvectors(num_conv);
  }
  Mpi::Print("\n");

  for (int i = 0; i < num_conv; i++)
  {
    // Get the eigenvalue and relative error.
    std::complex<double> omega = eigen->GetEigenvalue(i);
    double error_bkwd = eigen->GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error_abs = eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);
    if (!C && !has_A2)
    {
      // Linear EVP has eigenvalue μ = -λ² = ω².
      omega = std::sqrt(omega);
    }
    else
    {
      // Quadratic EVP solves for eigenvalue λ = iω.
      omega /= 1i;
    }

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
    eigen->GetEigenvector(i, E);

    linalg::NormalizePhase(space_op.GetComm(), E);

    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs.
      // B = -1/(iω) ∇ x E + 1/ω kp x E.
      floquet_corr->AddMult(E, B, 1.0 / omega);
    }

    auto total_domain_energy =
        post_op.MeasureAndPrintAll(i, E, B, omega, error_abs, error_bkwd, num_conv);

    // Calculate and record the error indicators.
    if (i < iodata.solver.eigenmode.n)
    {
      AddEstimate(E, B, total_domain_energy, indicator);
    }

    // Final write: Different condition than end of loop (i = num_conv - 1).
    if (i == iodata.solver.eigenmode.n - 1)
    {
      post_op.MeasureFinalize(indicator);
    }
  }
  MFEM_VERIFY(num_conv >= iodata.solver.eigenmode.n, "Eigenmode solve only found "
                                                         << num_conv << " modes when "
                                                         << iodata.solver.eigenmode.n
                                                         << " were requested!");
  return {indicator, space_op.GlobalTrueVSize()};
}

}  // namespace palace
