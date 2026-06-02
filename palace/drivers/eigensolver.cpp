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
#include "models/lumpedportoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
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

  // Build Σ_p kₙ,p(λ) · Mwp_p at COMPLEX λ using the analytical kₙ² fit (so kₙ
  // is evaluated as the analytic continuation √(α − γλ²) rather than the real-ω
  // restriction kₙ(|Im λ|)). Diagnostic only — used to compare T_nonlinear (with
  // funcA2 on real ω axis) against T_analytic (with kₙ extended into the complex
  // plane). Caller must have called FitKnSq() first.
  std::unique_ptr<ComplexOperator> BuildComplex(std::complex<double> lam) const
  {
    if (Mwp_p_.empty())
    {
      return {};
    }
    MFEM_VERIFY(kn2_alpha_.size() == Mwp_p_.size(),
                "WavePortFactor::BuildComplex requires FitKnSq() first.");
    std::vector<std::complex<double>> coeffs(Mwp_p_.size());
    for (std::size_t k = 0; k < Mwp_p_.size(); k++)
    {
      coeffs[k] = KnComplex(static_cast<int>(k), lam);
    }
    return BuildScaled(coeffs);
  }

  // Build the Newton-Jacobian wave-port term, dA2/dλ = -i·dA2/dω, evaluated at
  // λ = i·ω for real positive ω. With Mwp_p purely imaginary (= i·M_real_p), this is
  //   -i · i · Σ (dkₙ_p/dω) · M_real_p = Σ (dkₙ_p/dω) · M_real_p
  // i.e. the imaginary part of Mwp_p with a real scalar coefficient (no extra i wrap).
  // We assemble it as a sum on the i*M_real ports with coefficient -i*(dkₙ/dω) so the
  // final operator stores the right combination in real/imag slots. Central FD on the
  // cheap (cross-section EVP) kₙ; step is fd_rel_·max(|ω|,1).
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
      const double kp = space_op_.GetWavePortOp().GetWavePortKn(port_indices_[k],
                                                                omega + h);
      const double km = space_op_.GetWavePortOp().GetWavePortKn(port_indices_[k],
                                                                omega - h);
      const double dkn_domega = (kp - km) / (2.0 * h);
      // -i scale on a (i·M_real)-shaped operator gives M_real with real coefficient.
      coeffs[k] = std::complex<double>(0.0, -dkn_domega);
    }
    return BuildScaled(coeffs);
  }

  // Per-port kₙ,p(ω) (real positive for a propagating mode). Exposed so callers can
  // diagnose the seed-pencil fit residual against the truth.
  double Kn(int local_port, double omega) const
  {
    return space_op_.GetWavePortOp().GetWavePortKn(port_indices_[local_port], omega);
  }

  std::size_t NumPorts() const { return port_indices_.size(); }

  // Fit kₙ²(ω) = γ_p ω² + α_p per port (exact for textbook closed/open waveguides:
  // kₙ² = ω²/c² − ωc²/c²). Two-point fit on real-axis kₙ samples is enough since the
  // model is two-coefficient. Used to evaluate kₙ analytically at complex λ for the
  // 6-point off-axis seed sampling, where the cross-section EVP can't be re-run at
  // complex frequency. After fit, kₙ_p(λ) = √(α_p − γ_p λ²) (where λ = i·ω so
  // λ² = −ω²).
  void FitKnSq(double omega_min, double omega_max)
  {
    kn2_alpha_.clear();
    kn2_gamma_.clear();
    kn2_alpha_.reserve(Mwp_p_.size());
    kn2_gamma_.reserve(Mwp_p_.size());
    const double w1 = omega_min;
    const double w2 = omega_max;
    for (std::size_t k = 0; k < port_indices_.size(); k++)
    {
      const double kn1 = space_op_.GetWavePortOp().GetWavePortKn(port_indices_[k], w1);
      const double kn2 = space_op_.GetWavePortOp().GetWavePortKn(port_indices_[k], w2);
      const double kn1_sq = kn1 * kn1;
      const double kn2_sq = kn2 * kn2;
      // Solve [1, w1²; 1, w2²] [α; γ] = [kn1²; kn2²]
      const double det = w2 * w2 - w1 * w1;
      const double gamma = (kn2_sq - kn1_sq) / det;
      const double alpha = kn1_sq - gamma * w1 * w1;
      kn2_alpha_.push_back(alpha);
      kn2_gamma_.push_back(gamma);
    }
  }

  // Evaluate kₙ_p analytically at complex λ from the kₙ² fit. λ = i·ω convention;
  // kₙ²(ω) = γω² + α gives kₙ²(λ) = α − γλ² (since ω² = −λ²).
  std::complex<double> KnComplex(int local_port, std::complex<double> lambda) const
  {
    return std::sqrt(kn2_alpha_[local_port] - kn2_gamma_[local_port] * lambda * lambda);
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
    MFEM_VERIFY(kn2_alpha_.size() == Mwp_p_.size(),
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
      coeffs[k] = std::complex<double>(
          space_op_.GetWavePortOp().GetWavePortKn(port_indices_[k], omega), 0.0);
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
  std::vector<double> kn2_alpha_, kn2_gamma_;
  double fd_rel_;
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

  // Wave-port factor: pre-assembled per-port boundary mass operators. The wave-port
  // contribution to A2(ω) is Σ_p kₙ,p(ω)·Mwp_p — only the scalar kₙ depends on ω, so we
  // skip rebuilding the boundary mass on every funcA2 call (a hot path for HYBRID
  // Newton refinement). Empty if there are no wave ports.
  WavePortFactor wp_factor(space_op);
  // Check if there are nonlinear terms and, if so, setup interpolation operator.
  auto funcA2 = [&space_op, &wp_factor](double omega) -> std::unique_ptr<ComplexOperator>
  {
    // Bulk contributions (farfield, surface conductivity). Wave-port is excluded — it
    // is added via the factored wp_factor below using the cached Mwp_p.
    auto bulk = space_op.GetExtraSystemMatrix<ComplexOperator>(
        omega, Operator::DIAG_ZERO, /*include_wave_ports=*/false);
    auto wp = wp_factor.Build(omega);
    if (bulk && wp)
    {
      return BuildParSumOperator({1.0 + 0.0i, 1.0 + 0.0i}, {bulk.get(), wp.get()});
    }
    if (bulk)
    {
      return bulk;
    }
    return wp;
  };
  auto funcP = [&space_op](std::complex<double> a0, std::complex<double> a1,
                           std::complex<double> a2,
                           double omega) -> std::unique_ptr<ComplexOperator>
  { return space_op.GetPreconditionerMatrix<ComplexOperator>(a0, a1, a2, omega); };
  const double target = iodata.solver.eigenmode.target;
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

    // Seed-pencil strategies for HYBRID. Three options selectable via env var:
    //   PALACE_NLEPS_SEED unset / "polynomial":
    //     3-point Newton interpolation of A2(ω) in λ — current default. Polynomial
    //     extrapolation off the imaginary λ axis can introduce spurious eigenvalues
    //     (eigenvalues of T_polynomial that don't correspond to T_nonlinear modes).
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
    enum class SeedStrategy { Polynomial, OffAxis, FixedOmega };
    SeedStrategy seed_strategy = SeedStrategy::Polynomial;
    if (const char *env = std::getenv("PALACE_NLEPS_SEED"))
    {
      std::string_view sv{env};
      if (sv == "off-axis") seed_strategy = SeedStrategy::OffAxis;
      else if (sv == "fixed-omega") seed_strategy = SeedStrategy::FixedOmega;
    }
    Mpi::Print(
        "\n NLEPS HYBRID seed strategy: {}\n",
        seed_strategy == SeedStrategy::Polynomial   ? "polynomial (3-pt on-axis)"
        : seed_strategy == SeedStrategy::OffAxis    ? "off-axis (6-pt LSQ)"
                                                    : "fixed-omega (constant A2)");

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
      wp_factor.FitKnSq(target, target_max);
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
    else
    {
      // Polynomial 3-point Newton interpolation (default).
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
      constexpr int n_dense = 500;
      std::vector<double> w_fit(3), w_dense(n_dense);
      for (int j = 0; j < 3; j++)
      {
        w_fit[j] = target + j * (target_max - target) / 2.0;
      }
      for (int i = 0; i < n_dense; i++)
      {
        w_dense[i] = target + (target_max - target) * i / (n_dense - 1.0);
      }
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
        for (int i = 0; i < n_dense; i++)
        {
          const double w = w_dense[i];
          const double truth = wp_factor.Kn(static_cast<int>(p), w);
          const double poly = c(0) + c(1) * w + c(2) * w * w;
          max_abs = std::max(max_abs, std::abs(truth));
          max_rel_num = std::max(max_rel_num, std::abs(poly - truth));
        }
        const double rel_err = (max_abs > 0.0) ? max_rel_num / max_abs : 0.0;
        const double f_target =
            iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target) /
            (2.0 * M_PI);
        const double f_target_max =
            iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target_max) /
            (2.0 * M_PI);
        if (rel_err > 1.0e-2)
        {
          Mpi::Warning(
              "Wave port {:d}: HYBRID linear seed will be poor (order-2 kₙ fit "
              "residual {:.3e} on [{:.3e}, {:.3e}] GHz exceeds 1%). Quasi-Newton "
              "refinement may take many iterations or fail to converge.\n",
              p, rel_err, f_target, f_target_max);
        }
        else
        {
          Mpi::Print(
              " Wave port {:d} (NLEPS HYBRID seed): order-2 kₙ fit residual {:.3e} "
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
  if (nonlinear_type == NonlinearEigenSolver::SLP)
  {
    eigen->SetOperators(*K, *C, *M, EigenvalueSolver::ScaleType::NONE);
    eigen->SetExtraSystemMatrix(funcA2);
    eigen->SetPreconditionerUpdate(funcP);
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
  if (C || has_A2 || nonlinear_type == NonlinearEigenSolver::SLP)
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
      1.0 + 0.0i, 1i * target, -target * target + 0.0i, target);
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
    qn->SetExtraSystemMatrix(funcA2);
    if (!wp_factor.empty())
    {
      // Analytical wave-port Jacobian via Σ_p (dkₙ,p/dω)·Mwp_p. Bulk farfield/surface-σ
      // contributions to dA2/dω are not yet exposed analytically — when present, the
      // wave-port-only derivative is incomplete and the solver falls back to FD instead
      // (handled by NOT installing the derivative callback). For most of the wave-port
      // dispersion benchmarks the bulk part is absent and the analytical path applies.
      auto bulk_probe = space_op.GetExtraSystemMatrix<ComplexOperator>(
          target, Operator::DIAG_ZERO, /*include_wave_ports=*/false);
      if (!bulk_probe)
      {
        qn->SetExtraSystemMatrixDerivative(
            [&wp_factor](double omega) -> std::unique_ptr<ComplexOperator>
            { return wp_factor.BuildJacobianTerm(omega); });
        // Diagnostic: provide the analytic-continuation A2(λ) so SetInitialGuess
        // can compare T_re (kₙ on real axis) vs T_analytic (full complex λ). The
        // off-axis seed path above already calls FitKnSq(); call here too in case
        // wave_port_only branch above wasn't taken.
        wp_factor.FitKnSq(target, iodata.solver.eigenmode.target_upper);
        qn->SetExtraSystemMatrixComplex(
            [&wp_factor](std::complex<double> lam) -> std::unique_ptr<ComplexOperator>
            { return wp_factor.BuildComplex(lam); });
      }
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
