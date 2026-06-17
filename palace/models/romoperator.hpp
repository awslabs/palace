// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_ROM_OPERATOR_HPP
#define PALACE_MODELS_ROM_OPERATOR_HPP

#include <complex>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/units.hpp"

namespace palace
{

class IoData;
class SpaceOperator;

namespace config
{

struct LinearSolverData;

}  // namespace config

// Helper classs to define inner product used in the PROM when doing circuit synthesis. The
// weight matrix is the simple domain mass matrix of the space with the tdof of the port
// zeroed out (using the diag policy materials) and replaced by the boundary port matrix.
//
// This 'hybrid' approach is different from the conventional mass matrix $M_{ab} = (v,
// w)_{\Omega} + \sum_i (v, w)_{\Gamma_i}$ where the tdof of the port contributed both to
// the domain $\Omega$ overlap and the port overlap $\Gamma_i$. Our method basically the
// tdof on the ports fully with $\Gamma_i$. This ensures correct orthogonality of PROM basis
// with the boundary form used when calculating the system response matrix (S-Matrix).
//
// The ParSumOperator defined in linalg/rap.hpp can only have a single diag policy so we
// can't use that approach here.
struct HybridBulkBoundaryOperator
{
  // Need this since SetEssentialTrueDofs makes pointer to these values, even though it
  // looks like it is passing by references.
  mfem::Array<int> port_tdof_list;

  std::unique_ptr<ParOperator> W_inner_product_weight_bulk = {};
  std::unique_ptr<ParOperator> W_inner_product_weight_port = {};

  void validate_operators_zero_bulk_tdof()
  {
    // We don't need to zero of PEC in the actual model, since vectors entering inner
    // product should already have the PEC contribution removed.
    MFEM_VERIFY(W_inner_product_weight_port,
                "Port weight operator must be constructed before validation!");
    if (W_inner_product_weight_bulk)
    {
      W_inner_product_weight_bulk->SetEssentialTrueDofs(port_tdof_list,
                                                        Operator::DIAG_ZERO);
      MFEM_VERIFY((W_inner_product_weight_bulk->NumCols() ==
                   W_inner_product_weight_port->NumCols()) &&
                      (W_inner_product_weight_bulk->NumRows() ==
                       W_inner_product_weight_port->NumRows()),
                  "Mismatch sizes of inner product matrices!");
    }
  }

  HybridBulkBoundaryOperator(const SpaceOperator &space_op,
                             DomainOrthogonalizationWeight domain_orthog_type);

  HybridBulkBoundaryOperator(mfem::Array<int> &&port_tdof_list_,
                             std::unique_ptr<ParOperator> &&W_inner_product_weight_bulk_,
                             std::unique_ptr<ParOperator> &&W_inner_product_weight_port_)
    : port_tdof_list{std::move(port_tdof_list_)},
      W_inner_product_weight_bulk{std::move(W_inner_product_weight_bulk_)},
      W_inner_product_weight_port{std::move(W_inner_product_weight_port_)}
  {
    validate_operators_zero_bulk_tdof();
  }

  auto NumCols() const { return W_inner_product_weight_port->NumCols(); }

  auto NumRows() const { return W_inner_product_weight_port->NumRows(); }

  void Mult(const Vector &x, Vector &y) const
  {
    if (W_inner_product_weight_bulk)
    {
      W_inner_product_weight_bulk->Mult(x, y);
      W_inner_product_weight_port->AddMult(x, y);
    }
    else
    {
      // Identify bulk matrix with port tdof removed
      y.Set(1.0, x);
      linalg::SetSubVector(y, port_tdof_list, 0.0);
      W_inner_product_weight_port->AddMult(x, y);
    }
  }

  double InnerProduct(const Vector &x, const Vector &y, Vector &v_workspace) const
  {
    // TODO(future): Optimize to not use mult if InnerProduct is wrapped using our ceed
    // operator machinery.
    v_workspace.SetSize(x.Size());
    v_workspace.UseDevice(x.UseDevice());
    Mult(x, v_workspace);
    return linalg::LocalDot(y, v_workspace);
  }

  double InnerProduct(MPI_Comm comm, const Vector &x, const Vector &y,
                      Vector &v_workspace) const
  {
    auto dot = InnerProduct(x, y, v_workspace);
    Mpi::GlobalSum(1, &dot, comm);
    return dot;
  }
};

// Class for handling minimal-rational interpolation of solutions in frequency space. Used
// as an error indicator and for selecting the next frequency sample points in PROM
// construction. Each excitation gets a separate MRI, so sample frequencies are not shared.
class MinimalRationalInterpolation
{
private:
  // (Complex-valued) upper-triangular matrix R from orthogonalization of the HDM samples.
  // Minimal rational interpolant (MRI) defined by the vector q of interpolation weights and
  // support points z is used as an error indicator.
  std::vector<ComplexVector> Q;
  std::size_t dim_Q = 0;
  Eigen::MatrixXcd R;
  Eigen::VectorXcd q;
  std::vector<double> z;

public:
  MinimalRationalInterpolation(std::size_t max_size);
  void AddSolutionSample(double omega, const ComplexVector &u, MPI_Comm comm,
                         Orthogonalization orthog_type);
  std::vector<double> FindMaxError(std::size_t N) const;

  const auto &GetSamplePoints() const { return z; }
};

//
// A class handling projection-based reduced order model (PROM) construction and use for
// adaptive fast frequency sweeps.
//
class RomOperator
{
protected:
  // Reference to HDM discretization (not owned).
  // TODO(C++20): Use std::reference_wrapper with incomplete types.
  SpaceOperator &space_op;

  // Used for constructing & reuse of RHS1.
  int excitation_idx_cache = 0;

  // HDM system matrices and excitation:
  // - System matrix is: A(ω) = K + iω C - ω² M + A2(ω).
  // - Excitation / drive: = iω RHS1 + RHS2(ω).
  // - Vector r is internal vector workspace of size RHS
  // - The non-quadratic in ω operators A2(ω) and RHS2(ω) are built on fly in SolveHDM.
  // - Need to recompute RHS1 when excitation index changes (cf excitation_idx_cache).
  std::unique_ptr<ComplexOperator> K, M, C, A2;
  ComplexVector RHS1, RHS2, r;

  // Per-port wave-port boundary mass operators M^(p)_{μ⁻¹} (ω-independent). The
  // wave-port contribution to A(ω) factors as i·Σ_p kₙ,p(ω)·M^(p)_{μ⁻¹}; storing the
  // ω-independent operator once and the per-port reduced projection separately lets the
  // PROM online phase apply the wave-port term without touching any HDM-size object.
  // Keys are wave-port indices. Mwp_p carries the boundary mass in its imaginary slot
  // only (the i in i·Σ_p kₙ,p(ω)·M is baked in at assembly time), so Mwp_p_r is purely
  // imaginary; SolvePROM multiplies by the real scalar kₙ,p(ω) to get the i·kₙ·M
  // contribution to Aᵣ(ω) without an extra factor of i. The synthesis path inspects
  // Mp_r.imag() to recover the real symmetric M_proj used for SVD-based augmentation.
  std::map<int, std::unique_ptr<ComplexOperator>> Mwp_p;
  std::map<int, Eigen::MatrixXcd> Mwp_p_r;
  // True iff GetExtraSystemMatrix has any non-wave-port contributors (currently
  // second-order farfield or surface conductivity). Set in the ctor; controls the slow
  // SolvePROM A2 fallback path. NOTE: these contributors are now ALSO folded into the
  // circuit synthesis (see M_ff_, Asig_g_ below); has_other_A2 still gates the online
  // SolvePROM fallback, which is independent of synthesis.
  bool has_other_A2 = false;
  // One-time self-check state for the factored online A2_other path (SolvePROM). On the
  // first factored online solve we verify that the factored sum (0.5/ω·M_ff_r +
  // Σ_g EvaluateScalar(g,ω)/i·Asig_g_r) reproduces the full HDM projection of A2_other(ω)
  // to round-off. This guards against a future ω-dependent non-wave-port BC being added to
  // GetExtraSystemMatrix without a matching factored operator (which would otherwise be
  // silently dropped). On failure we set other_A2_factored_ok = false and permanently use
  // the slow per-ω HDM fallback.
  mutable bool other_A2_self_checked = false;
  mutable bool other_A2_factored_ok = true;

  // ω-independent boundary operators for the other frequency-dependent BCs, folded into
  // circuit synthesis the same way as the wave ports (each contributes i·f(ω)·M_proj·v to
  // Aᵣ(ω); project M_proj once, fit/inject the scalar f(ω)). Stored with the boundary mass
  // in the IMAGINARY slot (matching Mwp_p) so the synthesis path treats them uniformly.
  //   - Second-order farfield ABC: f(ω) = 0.5/ω, a single exact pole at ω=0 (residue 0.5),
  //     injected analytically (no fit). M_ff_ is the curl-curl boundary matrix.
  //   - Surface conductivity, one entry per attribute group g: f(ω) = ω/Z_g(ω) (the i is
  //     the implicit slot factor; surf_sigma_op.EvaluateScalar returns i·ω/Z, so the scalar
  //     used here is EvaluateScalar/i, generally complex), fit by poly+AAA.
  std::unique_ptr<ComplexOperator> M_ff_;
  Eigen::MatrixXcd M_ff_r;
  std::vector<std::unique_ptr<ComplexOperator>> Asig_g_;
  std::vector<Eigen::MatrixXcd> Asig_g_r;

  // Sweep band [ω_min, ω_max] (nondimensional, rad) captured from iodata at construction
  // time. Used to (a) sample kₙ,p(ω) for the synthesis polynomial fit, and (b) define the
  // dense grid over which the fit residual is evaluated.
  double sweep_omega_min = 0.0;
  double sweep_omega_max = 0.0;
  // Synthesis-side polynomial-fit configuration captured from iodata.
  double waveport_synthesis_tol = 0.0;
  std::size_t waveport_synthesis_order_max = 4;
  // Numerical-rank cutoff for the SVD of M_proj when augmenting in regime 2.
  double waveport_synthesis_rank_tol = 1.0e-6;
  // Force regime (auto/polynomial/augmented), see utils/labels.hpp.
  WavePortSynthesisRegime waveport_synthesis_force = WavePortSynthesisRegime::AUTO;

  // System properties: will be set when calling SetExcitationIndex & SolveHDM.
  bool has_A2 = true;
  bool has_RHS1 = true;
  bool has_RHS2 = true;

  // Per-excitation RHS2 sample cache. The wave-port excitation vector RHS2(ω) =
  // Σ_p∈excited 2ω·(n×Hₙ,p(ω)) does not factor cleanly in ω because the modal field
  // Hₙ,p(ω) deforms with ω. Cache the HDM RHS2 at every offline sample frequency,
  // project it into the basis (with rolling extensions in UpdatePROM), and interpolate
  // the small projected vector RHS2_r(ω) ∈ ℂⁿ in the online phase using barycentric
  // Lagrange. Samples accumulate across all greedy excitations and persist for the
  // online sweep over excitations.
  // Map<excitation_idx, Map<ω⋆, …>>.
  std::map<int, std::map<double, ComplexVector>> RHS2_hdm_samples;
  std::map<int, std::map<double, Eigen::VectorXcd>> RHS2_r_samples;

  // HDM linear system solver and preconditioner.
  std::unique_ptr<ComplexKspSolver> ksp;

  // PROM matrices and vectors. Projected matrices are Mr = Vᴴ M V where V is the reduced
  // order basis defined below.
  Eigen::MatrixXcd Kr, Mr, Cr;  // Extend during UpdatePROM as modes are added
  Eigen::VectorXcd RHS1r;       // Need to recompute drive vector on excitation change.

  // Frequency dependant PROM matrix Ar and RHSr are assembled and used only during
  // SolvePROM. Define them here so memory allocation can be reused in "online" evaluation.
  Eigen::MatrixXcd Ar;
  Eigen::VectorXcd RHSr;

  // PROM reduced-order basis (real-valued).
  std::vector<Vector> V;
  Orthogonalization orthog_type;

  // Weight operator for PROM basis if doing synthesis, in order to have correct
  // orthogonality to port vectors and converge with mesh / order. Default empty, which
  // means the identity weight in the finite element basis.
  std::optional<HybridBulkBoundaryOperator> weight_op_W = {};

  // Label to distinguish port modes from solution projection and to print PROM matrices.
  std::vector<std::string> v_node_label;

  // Upper-triangular matrix R from orthogonalization procedure U = VR. Here U the HDM
  // fields added by `UpdatePROM`.
  Eigen::MatrixXd orth_R;

  // MRIs: one for each excitation index. Only used to pick new frequency sample point.
  std::map<int, MinimalRationalInterpolation> mri;

  // Number of leading PROM basis rows occupied by lumped port modes. Each lumped port
  // flagged with IncludeInSynthesis contributes exactly one real basis vector, added
  // port-first and in order by AddLumpedPortModesForSynthesis. Ports excluded via
  // IncludeInSynthesis = false contribute nothing. This count is the invariant that the
  // port-port block scaling in CalculateNormalizedPROMMatrices and the basis reservation
  // in the constructor both rely on; iterating to the total port count instead would
  // over-run the per-port scaling vector whenever a port is excluded.
  std::size_t NumSynthesisPortModes() const;

  // Number of wave ports included in synthesis (IncludeInSynthesis = true). Each
  // contributes one port mode at the reference frequency; that mode is generally complex,
  // so AddWavePortModesForSynthesis adds up to two real basis vectors (real + imaginary
  // parts) per included port. Excluded ports add nothing. Used to size the basis
  // reservation so it never exceeds the count actually added.
  std::size_t NumSynthesisWavePortModes() const;

  // Helper function that normalizes PROM matrices so that they correspond to proper
  // admittance matrices on the ports. Also does unit conversion to physical (input) units.
  // Define so that 1.0 on port i corresponds to full (un-normalized solution), so you can
  // use Linv, Rinv, C directly.
  //
  // For wave ports whose order-2 polynomial fit residual exceeds the user-set tolerance,
  // the synthesised matrices can be enlarged by auxiliary rows/columns for the selected
  // rational realization. The returned `aux_labels` is the (possibly empty) list of those
  // new labels in the same order as the rows/columns appended to the matrices (i.e. each
  // matrix has size (V.size() + aux_labels.size())).
  struct NormalizedMatrices
  {
    std::unique_ptr<Eigen::MatrixXcd> L_inv;
    std::unique_ptr<Eigen::MatrixXcd> R_inv;  // null if no dissipative contribution
    std::unique_ptr<Eigen::MatrixXcd> C;
    std::vector<std::string> aux_labels;
  };
  NormalizedMatrices CalculateNormalizedPROMMatrices(const Units &units) const;

  // Wave-port dispersion synthesis helpers used by CalculateNormalizedPROMMatrices.
  //
  // The chosen regime determines whether k_n,p(omega) is approximated by an order-2
  // polynomial (Polynomial: clean absorption into Kr/Cr/Mr), augmented with the legacy
  // AAA-fit residual poles (Augmented), or replaced by a passive rational admittance
  // model (PassiveRational: extra aux states extend the pencil).
  enum class WavePortRegime
  {
    Polynomial,
    Augmented,
    PassiveRational,
    DtnRational
  };

  // Per-pole singular-vector aux block for a single wave port in regime 2. One aux
  // state per (pole, kept singular direction); the augmented pencil stores them as
  // appended rows/columns to Kr/Cr/Mr.
  struct WavePortAuxBlock
  {
    int port_idx = 0;
    // Aux-state label prefix. Empty -> "waveport_<port_idx>" (the wave-port default). Set
    // to e.g. "farfield" or "surfsigma_<g>" for the other frequency-dependent BCs so their
    // aux rows are distinguishable in the emitted matrices.
    std::string label;
    std::vector<double> sigmas;           // singular values kept above tolerance
    std::vector<Eigen::VectorXd> u_dirs;  // matching left singular vectors
    std::vector<std::complex<double>> poles;
    std::vector<std::complex<double>> residues;
  };

  // Passive fixed-pole admittance realization for one wave port. The scalar modal
  // admittance is approximated as
  //   y(s) = d + Σ_k r_k / (s + a_k),  a_k > 0, r_k >= 0, d >= 0,
  // and the boundary contribution is s*y(s)*M_proj. Eliminating the aux states realizes
  // each term r_k*s/(s+a_k) = r_k - a_k*r_k/(s+a_k) with stable dissipative poles.
  struct PassiveWavePortAuxBlock
  {
    int port_idx = 0;
    std::string label;
    std::vector<double> sigmas;
    std::vector<Eigen::VectorXd> u_dirs;
    std::vector<double> poles;     // Positive a_k in s+a_k.
    std::vector<double> residues;  // Nonnegative r_k.
  };

  // Result of fitting k_n,p(omega) on the sweep band for a single port. Polynomial regime
  // populates {alpha0, alpha1, alpha2}; Augmented regime additionally populates `aux`
  // and folds the AAA constant `d` into `alpha0`. PassiveRational populates `passive_aux`
  // and `passive_d` instead; the polynomial coefficients are retained only for diagnostics.
  struct WavePortDispersionFit
  {
    int port_idx = 0;
    WavePortRegime regime = WavePortRegime::Polynomial;
    double alpha0 = 0.0;
    double alpha1 = 0.0;
    double alpha2 = 0.0;
    // Complex polynomial coefficients, used by the generalized scalar-dispersion path
    // (FitScalarDispersion, for surface conductivity) where f(ω) is complex. The wave-port
    // path leaves these zero and uses the real alpha0/1/2 above.
    std::complex<double> alpha0c = 0.0;
    std::complex<double> alpha1c = 0.0;
    std::complex<double> alpha2c = 0.0;
    std::optional<WavePortAuxBlock> aux;
    std::optional<PassiveWavePortAuxBlock> passive_aux;
    double passive_d = 0.0;
    // DtnRational regime: structured √ fit of q(s)=i·kₙ via q²=c₀+c₁s+c₂s², realized as
    // q(s)=dtn_poly0 + dtn_poly1·s + Σ R_k/(s−p_k). The complex s-domain pole/residue list is
    // stored in `dtn_aux` (poles p_k, residues R_k, with the per-direction σ_j/u_j SVD data),
    // and dtn_poly0/dtn_poly1 carry the proper-part polynomial. The contribution to Aᵣ(ω) is
    // q(iω)·M_proj. Handles lossless (jω-axis poles, Re(q(iω))=0 by construction) and lossy
    // (LHP poles) uniformly; no spurious loss on a lossless port.
    std::optional<WavePortAuxBlock> dtn_aux;
    std::complex<double> dtn_poly0 = 0.0;
    std::complex<double> dtn_poly1 = 0.0;
    double rel_err_polynomial = 0.0;  // residual of α-only fit on dense grid
    double rel_err_augmented = 0.0;   // residual after AAA augmentation (Augmented only)
    double rel_err_passive = 0.0;     // residual after passive rational fit
    double rel_err_dtn = 0.0;         // residual after DtN structured-√ fit
  };

  // Choose a synthesis regime given the polynomial-fit residual. `meets_tol` is the
  // outcome of comparing rel_err against waveport_synthesis_tol; `force_setting`
  // overrides AUTO when set. Emits a Mpi::Warning if the chosen regime is incompatible
  // with the residual (e.g. forced Polynomial with rel_err > tol).
  WavePortRegime SelectWavePortRegime(int port_idx, double rel_err, bool meets_tol) const;

  // Sample kₙ,ₚ(ω) on the sweep band, fit a quadratic, run AAA on the residual when
  // augmenting, and pack the result. The dense-grid residual is captured for both
  // regimes for diagnostic logging.
  WavePortDispersionFit FitWavePortDispersion(int port_idx,
                                              const Eigen::MatrixXcd &Mp_r) const;

  // Evaluate the fitted real wave-port dispersion model at frequency ω:
  //   kₙ(ω) ≈ α₀ + α₁ω + α₂ω² + Σₖ Re(rₖ / (ω − pₖ))
  // (the AAA pole–residue residual term is present only in the Augmented regime). This is
  // the single source of truth for evaluating a WavePortDispersionFit on the real axis —
  // shared by the dense-grid residual check in FitWavePortDispersion and the online
  // SolvePROM surrogate path so the two never diverge. Only valid for the real-coefficient
  // wave-port fit (alpha0/1/2 + aux poles); the complex scalar-dispersion fit
  // (alpha0c/1c/2c, surface conductivity) is evaluated separately.
  static double EvaluateWavePortKnFit(const WavePortDispersionFit &fit, double omega);

  // Evaluate the scalar multiplying Mp_r in the online reduced system. For the polynomial
  // and legacy augmented regimes this is real-valued k_n(omega). For PassiveRational it is
  // omega*y(i omega), generally complex, so scalar*Mp_r realizes s*y(s)*M_proj on the real
  // frequency axis just like the synthesized augmented pencil.
  static std::complex<double>
  EvaluateWavePortMultiplierFit(const WavePortDispersionFit &fit, double omega);

  // Accumulate the polynomial-fit corrections from one port into the running Kr_corr,
  // Cr_corr, Mr_corr buffers. Folds α₀+d into Im(Kr), -α₁ into Re(Cr), -α₂ into Im(Mr).
  static void ApplyPolynomialFitCorrections(const WavePortDispersionFit &fit,
                                            const Eigen::MatrixXcd &Mp_r,
                                            Eigen::MatrixXcd &Kr_corr,
                                            Eigen::MatrixXcd &Cr_corr,
                                            Eigen::MatrixXcd &Mr_corr);

  // Accumulate a passive modal admittance fit y(s)=d+Σr/(s+a) as
  // s*y(s)*M_proj = s*d*M_proj + Σ r*M_proj - Σ a*r*M_proj/(s+a). The pole terms are
  // realized by PassiveWavePortAuxBlock in BuildAugmentedPencil.
  static void ApplyPassiveRationalFitCorrections(const WavePortDispersionFit &fit,
                                                 const Eigen::MatrixXcd &Mp_r,
                                                 Eigen::MatrixXcd &Kr_corr,
                                                 Eigen::MatrixXcd &Cr_corr);

  // Fit the DtN multiplier q(s)=i·kₙ(ω) of one wave port by the structured square root:
  // fit q²=c₀+c₁s+c₂s² (s=iω), complete the square q=√c₂·u·√(1−t) with u=s+c₁/(2c₂),
  // t=ω_s²/u², fit √(1−t)≈P(t)/Q(t) (real coeffs, column-scaled LSQ), then convert to the
  // s-domain pole-residue form q(s)=poly0+poly1·s+Σ R_k/(s−p_k). Populates fit.dtn_poly0,
  // dtn_poly1 and fit.dtn_aux (complex s-poles/residues + SVD directions of M_proj), and the
  // dense-grid residual fit.rel_err_dtn. Lossless ports give jω-axis poles (Re(q(iω))=0 by
  // construction, no spurious loss); lossy ports give LHP poles carrying the physical loss.
  void FitWavePortDtnRational(WavePortDispersionFit &fit, const Eigen::MatrixXcd &Mp_r,
                              const std::vector<double> &fit_omegas,
                              const Eigen::VectorXd &y_fit,
                              const std::vector<double> &dense_omegas) const;

  // Accumulate the DtN proper-polynomial part: q(s)=poly0+poly1·s+(pole part). The
  // contribution to Aᵣ(ω) is q(iω)·M_proj, so poly0·M_proj folds into Kr and poly1·M_proj
  // into Cr (the poly1·s = poly1·iω → iω·C channel). M_proj = (−i)·Mp_r. The pole part is
  // realized by the DtN branch of BuildAugmentedPencil. (poly0/poly1 are generally complex
  // for a lossy port; purely imaginary-on-jω structure is preserved for lossless.)
  static void ApplyDtnRationalFitCorrections(const WavePortDispersionFit &fit,
                                             const Eigen::MatrixXcd &Mp_r,
                                             Eigen::MatrixXcd &Kr_corr,
                                             Eigen::MatrixXcd &Cr_corr);

  // Generalized scalar-dispersion fit for the other frequency-dependent BCs (surface
  // conductivity). Like FitWavePortDispersion but driven by an arbitrary scalar function
  // f(ω) (here ω/Z_g(ω), generally COMPLEX) instead of the wave-port kₙ(ω). The projected
  // boundary mass `Mp_r` is purely imaginary (= i·M_proj) so the assembled contribution is
  // i·f(ω)·M_proj·v. Produces complex α₀/α₁/α₂ and, when augmenting, an aux block with
  // complex residues. `label` is used in the diagnostic print and aux-state labels.
  WavePortDispersionFit
  FitScalarDispersion(const std::string &label, const Eigen::MatrixXcd &Mp_r,
                      const std::function<std::complex<double>(std::complex<double>)> &f,
                      bool allow_augment) const;

  // Complex-coefficient variant of ApplyPolynomialFitCorrections for the other BCs: folds
  // α₀·Mp_r into Kr, -i·α₁·Mp_r into Cr, -α₂·Mp_r into Mr with COMPLEX α (the wave-port
  // version assumes real α). With Mp_r = i·M_proj this reproduces i·(α₀+α₁ω+α₂ω²)·M_proj·v.
  static void ApplyComplexPolynomialFitCorrections(
      std::complex<double> alpha0, std::complex<double> alpha1, std::complex<double> alpha2,
      const Eigen::MatrixXcd &Mp_r, Eigen::MatrixXcd &Kr_corr, Eigen::MatrixXcd &Cr_corr,
      Eigen::MatrixXcd &Mr_corr);

  // Build an aux block for one (matrix, pole-residue list) contribution, used by the ABC
  // (analytic single pole at ω=0, residue 0.5) and surf-σ (fitted poles). `Mp_r` is the
  // projected purely-imaginary boundary mass; the SVD of its imaginary part gives the
  // coupling directions exactly as in the wave-port augmented path.
  static std::optional<WavePortAuxBlock>
  MakeAuxBlock(int label_idx, const Eigen::MatrixXcd &Mp_r,
               const std::vector<std::complex<double>> &poles,
               const std::vector<std::complex<double>> &residues, double rank_tol);

  // Append aux-state rows/columns for regime-2 wave ports onto an n×n base pencil
  // (Kr_total, Cr_total, Mr_total). Returns the (n+aux)×(n+aux) augmented matrices and
  // appends labels to `aux_labels`. The base matrices are left at size n×n if no aux
  // states are present (no copy).
  struct AugmentedPencil
  {
    Eigen::MatrixXcd Kr;
    Eigen::MatrixXcd Cr;
    Eigen::MatrixXcd Mr;
  };
  static AugmentedPencil
  BuildAugmentedPencil(const Eigen::MatrixXcd &Kr_total, const Eigen::MatrixXcd &Cr_total,
                       const Eigen::MatrixXcd &Mr_total,
                       const std::vector<WavePortAuxBlock> &aux_blocks,
                       const std::vector<PassiveWavePortAuxBlock> &passive_aux_blocks,
                       const std::vector<WavePortAuxBlock> &dtn_aux_blocks,
                       std::vector<std::string> &aux_labels);

public:
  RomOperator(const config::LinearSolverData &linear, int verbose, SpaceOperator &space_op,
              std::size_t max_size_per_excitation);
  RomOperator(const IoData &iodata, SpaceOperator &space_op,
              std::size_t max_size_per_excitation);

  // Return the HDM linear solver.
  const ComplexKspSolver &GetLinearSolver() const { return *ksp; }

  // Return PROM dimension.
  auto GetReducedDimension() const { return V.size(); }

  // Return set of sampled parameter points for basis construction.
  const auto &GetSamplePoints(int excitation_idx) const
  {
    return mri.at(excitation_idx).GetSamplePoints();
  }

  // Set excitation index to build corresponding RHS vector (linear in frequency part).
  void SetExcitationIndex(int excitation_idx);

  // Assemble and solve the HDM at the specified frequency.
  void SolveHDM(int excitation_idx, double omega, ComplexVector &u);

  // Add port fields to PROM for circuit synthesis to connect to the outside world. Requires
  // fields that are orthogonal to boundary overlap and primary fields, which requires
  // custom weight matrix in orthogonality — see weight_op_W.
  void AddLumpedPortModesForSynthesis();

  // As above for wave ports. Adds one basis vector per (port, mode) seeded from the
  // cross-section EVP at `omega_ref`. The reference frequency should be inside the
  // sweep band (typically the band centre); the choice rescales the basis vector but
  // does not change correctness. Mode field is generally complex, so up to two basis
  // vectors are added per wave port (real + imaginary parts that survive the
  // orthogonalisation tolerance).
  void AddWavePortModesForSynthesis(double omega_ref);

  // Add field configuration to the reduced-order basis and update the PROM. Requires a name
  // "node_label". This will be printed in the header of the csv files when printing PROM
  // matrices. It is needed to distinguish port and solution field configuration as well as
  // to reconstruct if field configuration are pure real, imaginary or complex.
  void UpdatePROM(const ComplexVector &u, std::string_view node_label);

  // Add solution u to the minimal-rational interpolation for error estimation. MRI are
  // separated by excitation index.
  void UpdateMRI(int excitation_idx, double omega, const ComplexVector &u);

  // Assemble and solve the PROM at the specified frequency, expanding the solution back
  // into the high-dimensional space.
  void SolvePROM(int excitation_idx, double omega, ComplexVector &u);

  // Compute the location(s) of the maximum error in the range of the previously sampled
  // parameter points.
  std::vector<double> FindMaxError(int excitation_idx, std::size_t N = 1) const
  {
    return mri.at(excitation_idx).FindMaxError(N);
  }

  // Compute eigenvalue estimates for the current PROM system.
  std::vector<std::complex<double>> ComputeEigenvalueEstimates() const;

  // Print PROM matrices to file include in input (SI) units.
  void PrintPROMMatrices(const Units &units, const fs::path &post_dir) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_ROM_OPERATOR_HPP
