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
  // SolvePROM A2 fallback path.
  bool has_other_A2 = false;

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
  // For wave ports whose order-2 polynomial fit residual exceeds the user-set tolerance
  // (regime 2), the synthesised matrices are enlarged by one auxiliary row/column per
  // AAA pole used to fit δkₙ(ω). Each aux row carries a synthetic node label of the
  // form `waveport_<idx>_aux_<k>`. The returned `aux_labels` is the (possibly empty)
  // list of those new labels in the same order as the rows/columns appended to the
  // matrices (i.e. each matrix has size (V.size() + aux_labels.size())).
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
  // The chosen regime determines whether kₙ,ₚ(ω) is approximated by an order-2
  // polynomial (Polynomial: clean absorption into Kr/Cr/Mr) or augmented with AAA-fit
  // residual poles (Augmented: extra aux states extend the pencil).
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
    std::vector<double> sigmas;           // singular values kept above tolerance
    std::vector<Eigen::VectorXd> u_dirs;  // matching left singular vectors
    std::vector<std::complex<double>> poles;
    std::vector<std::complex<double>> residues;
  };

  // Result of fitting kₙ,ₚ(ω) on the sweep band for a single port. Polynomial regime
  // populates {alpha0, alpha1, alpha2}; Augmented regime additionally populates `aux`
  // and folds the AAA constant `d` into `alpha0`.
  struct WavePortDispersionFit
  {
    int port_idx = 0;
    WavePortRegime regime = WavePortRegime::Polynomial;
    double alpha0 = 0.0;
    double alpha1 = 0.0;
    double alpha2 = 0.0;
    std::optional<WavePortAuxBlock> aux;
    double rel_err_polynomial = 0.0;  // residual of α-only fit on dense grid
    double rel_err_augmented = 0.0;   // residual after AAA augmentation (Augmented only)
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

  // Accumulate the polynomial-fit corrections from one port into the running Kr_corr,
  // Cr_corr, Mr_corr buffers. Folds α₀+d into Im(Kr), -α₁ into Re(Cr), -α₂ into Im(Mr).
  static void ApplyPolynomialFitCorrections(const WavePortDispersionFit &fit,
                                            const Eigen::MatrixXcd &Mp_r,
                                            Eigen::MatrixXcd &Kr_corr,
                                            Eigen::MatrixXcd &Cr_corr,
                                            Eigen::MatrixXcd &Mr_corr);

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
  static AugmentedPencil BuildAugmentedPencil(
      const Eigen::MatrixXcd &Kr_total, const Eigen::MatrixXcd &Cr_total,
      const Eigen::MatrixXcd &Mr_total,
      const std::vector<WavePortAuxBlock> &aux_blocks,
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
