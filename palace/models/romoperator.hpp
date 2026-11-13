// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_ROM_OPERATOR_HPP
#define PALACE_MODELS_ROM_OPERATOR_HPP

#include <complex>
#include <limits>
#include <map>
#include <set>
#include <memory>
#include <optional>
#include <span>
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
struct FloquetOrder;

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
// The samples are given by their coefficients in an orthonormal basis, which may grow
// between samples (the coefficients of the later basis vectors are zero for the earlier
// samples), and are replicated on every process.
class MinimalRationalInterpolation
{
private:
  // QR decomposition of the snapshot matrix with columns [c; iω c] for the coefficients c
  // of each sample, with rows interleaved so that a new basis vector only adds rows.
  // Minimal rational interpolant (MRI) defined by the vector q of interpolation weights and
  // support points z is used as an error indicator.
  Eigen::MatrixXcd Q, R;
  Eigen::VectorXcd q;
  std::vector<double> z;

public:
  void AddSolutionSample(double omega, const Eigen::VectorXcd &c);
  // With a communicator, the search for N = 1 is split over its processes (same result).
  std::vector<double> FindMaxError(std::size_t N, MPI_Comm comm = MPI_COMM_NULL) const;

  const auto &GetSamplePoints() const { return z; }
};

// Dense solver for the PROM system A X = B, which is solved redundantly on every process,
// overwriting B with X. LU with partial pivoting is blocked, so it stays cheap at large n
// when all processes factor concurrently, and backward stable unless the pivots grow. It is
// also more accurate than full-pivoting QR near a lossless resonance, where the QR
// truncates the numerical rank and drops the resonant mode. The QR is used, on all
// processes together, when an LU solution is not finite (an exactly singular matrix) or a
// column has a normwise backward error above max(n, 16) ε (pivot growth) or not finite. The
// componentwise backward error, while invariant to scaling, exceeds n ε for LU solutions of
// badly scaled systems without pivot growth. Returns true if the QR was used.
bool SolveReducedSystem(MPI_Comm comm, const Eigen::MatrixXcd &A, Eigen::MatrixXcd &B);

// Projections Vᵀ A V of symmetric operators A onto the real vectors V, replicated on every
// process.
std::vector<Eigen::MatrixXd> ProjectSymmetricOperators(MPI_Comm comm,
                                                       const std::vector<Vector> &V,
                                                       std::span<const Operator *const> A);

//
// A class handling projection-based reduced order model (PROM) construction and use for
// adaptive fast frequency sweeps.
//
class RomOperator
{
protected:
  // Reference to HDM discretization (not owned).
  SpaceOperator &space_op;

  // HDM system matrices and excitation:
  // - System matrix is: A(ω) = K + iω C - ω² M + A2(ω) + F(ω).
  // - Excitation / drive: = iω RHS1 + RHS2(ω).
  // - Vector r is internal vector workspace of size RHS
  // - The non-quadratic in ω operators A2(ω), F(ω), and RHS2(ω) are built on the fly.
  // - RHS1 is real and assembled once per excitation (empty if zero on every process), in
  //   the order of the port excitations, and projected onto the basis as it grows.
  std::unique_ptr<ComplexOperator> K, M, C;
  std::vector<Vector> RHS1;
  std::vector<bool> has_RHS1;  // The same on every process
  ComplexVector r;

  // The sparse frequency-dependent boundary terms of A(ω), each a scalar f(ω) times an
  // ω-independent operator i·M stored with the boundary form M on the imaginary slot, so
  // that the PROM gets f(ω)·Vᵀ(i·M)V without any HDM-size work online and circuit
  // synthesis can fit f(ω) on the sweep band:
  //   - Wave port p: f = kₙ,ₚ(ω) and the µ⁻¹ boundary mass. Inactive ports included in
  //     synthesis are unloaded nodes: kept for their node scaling and reference data, they
  //     do not contribute to A(ω).
  //   - Second-order farfield ABC: f = 0.5/ω, with the negative sign of β inside M (see
  //     SpaceOperator::GetFarfieldBoundaryCurlCurlMatrix).
  //   - Surface conductivity, per attribute group: f = ω/Z(ω).
  //   - Rational surface impedance, per boundary: f = g(iω)/i for the Robin coefficient
  //     g(s) = s·D(s)/N(s).
  //   - Floquet port Robin BC, per port: f = γ₀,ₚ(ω), the specular (0,0) propagation
  //     constant.
  enum class BoundaryTermKind
  {
    WavePort,
    Farfield,
    SurfaceConductivity,
    RationalImpedance,
    FloquetPort
  };
  struct BoundaryTerm
  {
    BoundaryTermKind kind;
    int idx;  // Port, group, or boundary index
    std::unique_ptr<ComplexOperator> op;
    bool symmetric;
    Eigen::MatrixXcd Ar;  // Vᵀ op V
  };
  std::vector<BoundaryTerm> bdr_terms;

  // Return f(ω) of a boundary term, zero if it does not contribute to A(ω).
  std::complex<double> GetBoundaryTermCoefficient(const BoundaryTerm &term,
                                                  double omega) const;

  // Port-space evaluation of the wave-port terms of the PROM. The modal correction
  // Wᵣ(ω) = Σ_k g_k(ω) (Vᵀs_k)(Vᵀs_k)ᵀ, the excitation RHS2ᵣ(ω) = −2iω Vᵀs_full, and the
  // S-parameters only need the pairings Vᵀs of the basis with the per-port modal n×H shape
  // vectors. Since s is supported on the port only, Vᵀs is computed from the basis
  // restricted to the port submesh (V_wp, kept on the host and extended as the basis
  // grows) and the port-space mode forms refreshed by WavePortData::Initialize, without
  // any per-frequency assembly and projection of HDM-size vectors. sV_wp_full caches
  // Vᵀs_full of every port at the frequency of the last SolvePROM.
  std::map<int, std::vector<Vector>> V_wp;
  std::size_t V_wp_dim = 0;
  std::map<int, Eigen::VectorXcd> sV_wp_full;

  // Add the wave-port modal correction Wᵣ(ω) to Ar and cache the pairings of the basis
  // with the mode of every wave port at ω.
  void AddWavePortModalCorrection(double omega);

  // Sweep band [ω_min, ω_max] (nondimensional, rad) captured from iodata at construction
  // time. Used to (a) sample kₙ,p(ω) for the synthesis polynomial fit, and (b) define the
  // dense grid over which the fit residual is evaluated.
  std::vector<double> sweep_omega_samples;
  double sweep_omega_min = 0.0;
  double sweep_omega_max = 0.0;
  // Synthesis-side fit configuration derived from AdaptiveTol. The AAA order and SVD rank
  // cutoff are internal guards so users only configure the adaptive sweep tolerance.
  double waveport_synthesis_tol = 0.0;
  std::size_t waveport_synthesis_order_max = 12;
  double waveport_synthesis_rank_tol = 1.0e-6;

  // HDM linear system solver and preconditioner.
  std::unique_ptr<ComplexKspSolver> ksp;

  // PROM matrices and vectors. Projected matrices are Mr = Vᴴ M V where V is the reduced
  // order basis defined below.
  Eigen::MatrixXcd Kr, Mr, Cr;  // Extend during UpdatePROM as modes are added
  Eigen::MatrixXcd RHS1r;       // Vᴴ RHS1, one column per excitation

  // Projections onto the basis of the Floquet port modes, for the DtN correction
  // F(ω) = Σ_k g_k(ω) v_k v_kᴴ and the excitation of the (0,0) modes (see FloquetPortData).
  // A mode v = Σ_d e_d v_d combines the Cartesian component projections v_d with the
  // polarization e, which rotates with frequency when k_F scales with it, so Vᵀv_d is
  // projected once as the basis grows: three columns per order.
  std::vector<std::pair<int, const FloquetOrder *>> floquet_orders;  // Port, order
  Eigen::MatrixXcd floquet_comp_V;

  // Return Vᵀv of the TE or TM mode of order k at the current frequency.
  Eigen::VectorXcd GetFloquetModeProjection(std::size_t k, bool is_te) const;

  // PROM matrix Ar at the frequency of the last SolvePROM, and the reduced solutions of all
  // excitations, one column per excitation.
  Eigen::MatrixXcd Ar, Y;
  bool warned_qr = false;
  std::set<int> warned_modal_reaction;  // Ports

  // Column of an excitation in RHS1, RHS1r, and Y.
  Eigen::Index GetExcitationColumn(int excitation_idx) const;

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

  // With a weight operator, the Cholesky factor L of the Euclidean Gram matrix of the
  // basis, VᵀV = L Lᵀ, so that ||V c|| = ||Lᵀ c|| for the error estimates. Empty otherwise
  // (V is orthonormal).
  Eigen::MatrixXd gram_L;

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

  // Wave-port dispersion synthesis helpers used by CalculateNormalizedPROMMatrices.
  //
  // The chosen regime determines whether k_n,p(omega) is approximated by an order-2
  // polynomial (clean absorption into Kr/Cr/Mr), or by that polynomial plus AAA-fit
  // residual poles.
  enum class WavePortRegime
  {
    Polynomial,
    Augmented
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
    // Signed weights and matching real directions in the symmetric factorization
    // M_proj = Σ_j weights[j] u_j u_jᵀ. The weights are eigenvalues, not unsigned
    // singular values: modal-correction residue matrices are generally indefinite.
    std::vector<double> weights;
    std::vector<Eigen::VectorXd> u_dirs;
    std::vector<std::complex<double>> poles;
    std::vector<std::complex<double>> residues;
  };

  // Result of fitting k_n,p(omega) (or the scalar dispersion f(ω) of another frequency-
  // dependent BC) on the sweep band. Polynomial regime populates {alpha0c, alpha1c,
  // alpha2c}; Augmented regime additionally populates `aux` and folds the AAA constant `d`
  // into `alpha0c`.
  struct WavePortDispersionFit
  {
    int port_idx = 0;
    WavePortRegime regime = WavePortRegime::Polynomial;
    // Complex polynomial coefficients. For a wave port, k_n,p(ω) is complex for a lossy
    // port mode (its imaginary part is the attenuation of the mode), and real up to
    // roundoff for a lossless one.
    std::complex<double> alpha0c = 0.0;
    std::complex<double> alpha1c = 0.0;
    std::complex<double> alpha2c = 0.0;
    std::optional<WavePortAuxBlock> aux;
    double rel_err_polynomial = 0.0;  // residual of α-only fit on dense grid
    double rel_err_augmented = 0.0;   // residual after AAA augmentation (Augmented only)
  };

  // Helper function that normalizes PROM matrices so that they correspond to proper
  // admittance matrices on the ports. Also does unit conversion to physical (input) units.
  // Define so that 1.0 on port i corresponds to full (un-normalized solution), so you can
  // use Linv, Rinv, C directly.
  //
  // For wave ports whose order-2 polynomial fit residual exceeds AdaptiveTol, the
  // synthesised matrices are enlarged by auxiliary rows/columns for the AAA residual
  // realization. The returned `aux_labels` is the (possibly empty) list of those
  // new labels in the same order as the rows/columns appended to the matrices (i.e. each
  // matrix has size (V.size() + aux_labels.size())). `wave_port_fits` captures the exact
  // per-port fit used in the emitted matrices so the port-reference table can be evaluated
  // consistently with the L/R/C realization.
  struct NormalizedMatrices
  {
    struct PortLoad
    {
      std::string label;
      std::unique_ptr<Eigen::MatrixXcd> L_inv;
      std::unique_ptr<Eigen::MatrixXcd> R_inv;  // null if no dissipative contribution
      std::unique_ptr<Eigen::MatrixXcd> C;
    };

    std::unique_ptr<Eigen::MatrixXcd> L_inv;
    std::unique_ptr<Eigen::MatrixXcd> R_inv;  // null if no dissipative contribution
    std::unique_ptr<Eigen::MatrixXcd> C;
    std::vector<std::string> aux_labels;
    std::vector<PortLoad> port_loads;
    std::vector<WavePortDispersionFit> wave_port_fits;
  };
  NormalizedMatrices CalculateNormalizedPROMMatrices(const Units &units) const;

  // Print the matched reference admittance used to convert synthesized port admittances to
  // Palace/Kurokawa S-parameters. The L/R/C matrices remain loaded; downstream cascade
  // tools should use the per-port rom-portload-* matrices for matrix-level port-load
  // removal and this table for frequency-domain S-parameter normalization.
  void PrintPortReferenceData(const Units &units, const fs::path &post_dir,
                              const NormalizedMatrices &matrices) const;

  // Sample kₙ,ₚ(ω) on the sweep band, fit a quadratic, run AAA on the residual when
  // augmenting, and pack the result. The dense-grid residual is captured for both
  // regimes for diagnostic logging.
  WavePortDispersionFit FitWavePortDispersion(int port_idx,
                                              const Eigen::MatrixXcd &Mp_r) const;

  // Evaluate the fitted wave-port dispersion model at frequency ω:
  //   kₙ(ω) ≈ α₀ + α₁ω + α₂ω² + Σₖ rₖ / (ω − pₖ)
  // (the AAA pole–residue residual term is present only in the Augmented regime), the same
  // expression the synthesized pencil realizes. This is the single source of truth for
  // evaluating a wave-port WavePortDispersionFit.
  static std::complex<double> EvaluateWavePortKnFit(const WavePortDispersionFit &fit,
                                                    std::complex<double> omega);

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

  // Exact-split fit for a rational surface impedance boundary: the Robin coefficient
  // g(s) = s·D(s)/N(s) = P(s) + R(s)/N(s) (long division) has a polynomial part P that is
  // exactly representable in the pencil (deg(P) <= 2 for any passive input), so only the
  // strictly proper remainder R/N is AAA-fit — directly.
  WavePortDispersionFit FitRationalImpedanceDispersion(const std::string &label,
                                                       const Eigen::MatrixXcd &Mp_r,
                                                       int idx) const;

  // Accumulate the polynomial-fit corrections of one frequency-dependent BC (a wave port,
  // or another BC) into the running Kr_corr, Cr_corr, Mr_corr buffers: folds α₀·Mp_r into
  // Kr, -i·α₁·Mp_r into Cr, -α₂·Mp_r into Mr with complex α. With Mp_r = i·M_proj this
  // reproduces i·(α₀+α₁ω+α₂ω²)·M_proj·v in K + iωC - ω²M. Real α (a lossless wave port)
  // fold into Im(Kr), Re(Cr), and Im(Mr).
  static void ApplyComplexPolynomialFitCorrections(
      std::complex<double> alpha0, std::complex<double> alpha1, std::complex<double> alpha2,
      const Eigen::MatrixXcd &Mp_r, Eigen::MatrixXcd &Kr_corr, Eigen::MatrixXcd &Cr_corr,
      Eigen::MatrixXcd &Mr_corr);

  // Build an aux block for one (matrix, pole-residue list) contribution, used by the ABC
  // (analytic single pole at ω=0, residue 0.5), surf-σ, and rational impedance (fitted
  // poles). `label` (required, nonempty) prefixes the aux-state row names in the
  // synthesized matrices. `Mp_r` is the projected purely-imaginary boundary mass; a signed
  // symmetric eigendecomposition of its imaginary part gives the coupling directions
  // exactly, including indefinite modal-correction residue matrices.
  static std::optional<WavePortAuxBlock>
  MakeAuxBlock(std::string label, const Eigen::MatrixXcd &Mp_r,
               const std::vector<std::complex<double>> &poles,
               const std::vector<std::complex<double>> &residues, double rank_tol);
  static bool AddAuxBlockDirections(WavePortAuxBlock &blk, const Eigen::MatrixXcd &Mp_r,
                                    double rank_tol);

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
                       std::vector<std::string> &aux_labels);

public:
  RomOperator(const IoData &iodata, SpaceOperator &space_op,
              std::size_t max_size_per_excitation);

  // Return the HDM linear solver.
  const ComplexKspSolver &GetLinearSolver() const { return *ksp; }

  // Return PROM dimension and the reduced solution of an excitation at the frequency of the
  // last SolvePROM.
  auto GetReducedDimension() const { return V.size(); }
  const auto &GetBasis() const { return V; }
  Eigen::VectorXcd GetReducedSolution(int excitation_idx) const
  {
    return Y.col(GetExcitationColumn(excitation_idx));
  }

  // Return the pairing Vᵀs_full of the basis with the mode of a wave port at the frequency
  // of the last PROM assembly, or nullptr before the first one.
  const Eigen::VectorXcd *GetWavePortModePairing(int port_idx) const
  {
    const auto it = sV_wp_full.find(port_idx);
    return (it != sV_wp_full.end()) ? &it->second : nullptr;
  }

  // Return set of sampled parameter points for basis construction.
  const auto &GetSamplePoints(int excitation_idx) const
  {
    return mri.at(excitation_idx).GetSamplePoints();
  }

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
  // to reconstruct if field configuration are pure real, imaginary or complex. Returns the
  // coefficients c of u in the updated basis, u = V c.
  Eigen::VectorXcd UpdatePROM(const ComplexVector &u, std::string_view node_label);

  // Add the solution with basis coefficients c to the minimal-rational interpolation for
  // error estimation. MRI are separated by excitation index.
  void UpdateMRI(int excitation_idx, double omega, const Eigen::VectorXcd &c);

  // Return ||V c|| for the coefficients c of a vector in the basis (zero-padded to the
  // basis dimension).
  double BasisNorm(const Eigen::VectorXcd &c) const;

  // Assemble and solve the PROM at the specified frequency, for every excitation. The
  // reduced solutions are available from GetReducedSolution(excitation_idx).
  void SolvePROM(double omega);

  // Expand the reduced solution of an excitation into the high-dimensional space.
  void ExpandPROMSolution(int excitation_idx, ComplexVector &u) const;

  // Compute the location(s) of the maximum error in the range of the previously sampled
  // parameter points.
  std::vector<double> FindMaxError(int excitation_idx, std::size_t N = 1) const;

  // Print PROM matrices to file include in input (SI) units.
  void PrintPROMMatrices(const Units &units, const fs::path &post_dir) const;

  // Compute eigenvalue estimates of the synthesized L⁻¹/R⁻¹/C system via the companion
  // linearization (in s = iω) of the quadratic pencil (L⁻¹ + iωR⁻¹ − ω²C)v = 0. Returns
  // every finite root whose Re(f) lies inside [fmin_GHz, fmax_GHz]; spurious roots are left
  // in and flagged by their reported HDM residuals rather than a fixed cutoff. The
  // complex-symmetric pencil yields decaying modes with positive Im{f}, the same sign
  // convention as the eigenmode solver's eig.csv.
  struct EigenvalueEstimate
  {
    double freq_re_GHz;  // Re(f): resonant frequency
    double freq_im_GHz;  // Im(f): half-linewidth (positive = decaying)
    double Q;            // quality factor = |Re(ω)| / (2|Im(ω)|)
    // Eigenvector of the synthesized pencil (basis + aux rows), used to evaluate the
    // eigenpair residual on the HDM below.
    Eigen::VectorXcd eigvec;
    // HDM eigenpair errors, following the eigenmode solver conventions: prolongate the
    // eigenvector to the finite element space and evaluate the full frequency-dependent
    // operator at the complex eigenfrequency.
    double error_abs = -1.0;   // ‖(K + iωC − ω²M + A2(ω)) u‖₂ / ‖u‖₂
    double error_bkwd = -1.0;  // error_abs / (‖K‖₂ + |ω|‖C‖₂ + |ω|²‖M‖₂)
  };
  static std::vector<EigenvalueEstimate>
  ComputeEigenvalueEstimates(const Eigen::MatrixXcd &L_inv, const Eigen::MatrixXcd *R_inv,
                             const Eigen::MatrixXcd &C, double fmin_GHz, double fmax_GHz);

  // Fill error_abs / error_bkwd for each estimate by prolongating the eigenvector's
  // basis block to the HDM space.
  void ComputeEigenvalueEstimateErrors(const Units &units,
                                       std::vector<EigenvalueEstimate> &estimates) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_ROM_OPERATOR_HPP
