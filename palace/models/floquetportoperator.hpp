// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_FLOQUET_PORT_OPERATOR_HPP
#define PALACE_MODELS_FLOQUET_PORT_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "fem/gridfunction.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;

namespace config
{
struct FloquetPortData;
}  // namespace config

// Flags controlling which subsystems use a given Floquet mode.
enum class FloquetModeUse : std::uint8_t
{
  Output = 1,  // S-parameter CSV output (user-requested ±MaxOrder)
  Dtn = 2,     // DtN boundary correction (BZ-centered range)
  Both = 3     // Output | Dtn
};
inline bool HasFlag(FloquetModeUse u, FloquetModeUse flag)
{
  return (static_cast<uint8_t>(u) & static_cast<uint8_t>(flag)) != 0;
}

// Represents a single diffraction order (m, n) with a specific polarization (TE/TM).
struct FloquetMode
{
  int m, n;            // Lattice indices (physical convention)
  bool is_te;          // true = TE (s-pol), false = TM (p-pol)
  FloquetModeUse use;  // Which subsystems use this mode
  mfem::Vector B_mn;   // Transverse wavevector B_mn = m*b1 + n*b2
  mfem::Vector e_pol;  // Polarization unit vector (3D, tangential to port)
  ComplexVector v;     // Fourier projection: v_j = int_Gamma (nxnxN_j).e_p exp(-iB.r) dS
  double gamma_sq;     // gamma_mn^2 = omega^2*mu*eps - |B_mn - k_F|^2 (freq-dependent)

  // Mass-consistent boundary mass matrix: A[j,k] = int_Gamma (N_j.e_pol)(N_k.e_pol) dS.
  // Assembled once (frequency-independent, real). Only populated when MassConsistent mode
  // is active. Stored as HypreParMatrix on the true DOF space.
  std::unique_ptr<mfem::HypreParMatrix> A_mass;
};

// Low-rank complex operator: F*x = sum_k g_k conj(v_k) (v_k^T x).
// Bilinear convention: v^T x (no conjugation in inner product), conj(v) in outer vector.
// Used for the Floquet port boundary matrix.
class LowRankComplexOperator : public ComplexOperator
{
public:
  struct Term
  {
    const ComplexVector *v;  // Projection vector (not owned)
    std::complex<double> g;  // Frequency-dependent scalar weight
  };

private:
  std::vector<Term> terms;
  MPI_Comm comm;
  mutable ComplexVector tmp;

public:
  LowRankComplexOperator(MPI_Comm comm, int n) : ComplexOperator(n), comm(comm), tmp(n)
  {
    tmp.UseDevice(true);
  }

  void AddTerm(const ComplexVector *v, std::complex<double> g) { terms.push_back({v, g}); }

  void Mult(const ComplexVector &x, ComplexVector &y) const override;

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override;
};

// Dense boundary DtN operator: F_bdr is an N_bdr × N_bdr Eigen matrix applied via
// restrict → Allgatherv (on port_comm) → dense matvec → scatter. Cost: O(N_bdr²) per
// application, independent of mode count. Replaces LowRankComplexOperator for performance.
class DenseBoundaryOperator : public ComplexOperator
{
private:
  MPI_Comm port_comm;              // Sub-communicator for port-owning ranks only
  int n_full;                      // Full system size (number of true DOFs)
  mfem::Array<int> bdr_tdof_list;  // Local boundary DOF indices into the full vector
  int n_bdr_local;                 // Number of local boundary DOFs
  int n_bdr_global;                // Total boundary DOFs across port_comm ranks
  std::vector<int> recv_counts;    // Per-rank counts for Allgatherv on port_comm
  std::vector<int> recv_displs;    // Per-rank displacements for Allgatherv
  Eigen::MatrixXcd F;              // Dense DtN matrix, size n_bdr_global × n_bdr_global
  mutable Eigen::VectorXcd x_bdr_local, x_bdr_global, y_bdr_global;

public:
  DenseBoundaryOperator(MPI_Comm port_comm, int n_full,
                        const mfem::Array<int> &bdr_tdof_list, int n_bdr_global);

  // Get the dense boundary matrix (for reuse in the preconditioner Schur complement).
  const Eigen::MatrixXcd &GetDenseMatrix() const { return F; }
  int GetNBdrGlobal() const { return n_bdr_global; }

  // Set the dense matrix directly (e.g., after assembling from mode projections).
  void SetDenseMatrix(const Eigen::MatrixXcd &F_new) { F = F_new; }

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override;
};

// Mass-consistent DtN operator: F*x = sum_k g_k A_k x, where A_k is a pre-assembled
// real boundary mass matrix weighted by the mode polarization direction.
// Replaces the rank-1 LowRankComplexOperator to fix the Cauchy-Schwarz overestimate.
class MassConsistentDtNOperator : public ComplexOperator
{
public:
  struct Term
  {
    mfem::HypreParMatrix *A;  // Pre-assembled real boundary mass (not owned)
    std::complex<double> g;   // Frequency-dependent DtN coefficient
  };

private:
  std::vector<Term> terms;
  mutable Vector tmp;

public:
  MassConsistentDtNOperator(int n) : ComplexOperator(n), tmp(n) { tmp.UseDevice(true); }

  void AddTerm(mfem::HypreParMatrix *A, std::complex<double> g) { terms.push_back({A, g}); }

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override;
};

// Simple solver that applies a scalar times identity: y = scale * x.
// Used as the s-block preconditioner in the auxiliary system (A_SS = -I → P_SS = -I).
class ScaledIdentitySolver : public Solver<ComplexOperator>
{
private:
  std::complex<double> scale;

public:
  ScaledIdentitySolver(int n, std::complex<double> s = -1.0)
    : Solver<ComplexOperator>(false), scale(s)
  {
    this->height = n;
    this->width = n;
  }
  void SetOperator(const ComplexOperator &) override {}
  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    y.Real().Set(scale.real(), x.Real());
    y.Real().Add(-scale.imag(), x.Imag());
    y.Imag().Set(scale.imag(), x.Real());
    y.Imag().Add(scale.real(), x.Imag());
  }
};

// Auxiliary S-parameter DOF operator for the coupled [E, s] system.
// Applies:
//   [A_EE   A_ES] [E]     A_EE E + Σ_i g_full_i × s_i × conj(v_i)
//   [A_SE   A_SS] [s]  =  v_i^T E − s_i
//
// The E block uses the existing system operator (no Robin on Floquet boundaries).
// The s block couples modes to the boundary field through Fourier projections.
// The DtN boundary condition is enforced through the coupling, not through Robin.
//
// MPI convention: s DOFs are "owned" by rank 0. On other ranks, s entries in the
// extended vector are zero. This ensures correct norms and dot products in GMRES.
class FloquetAuxSystemOperator : public ComplexOperator
{
private:
  const ComplexOperator *A_EE;  // Existing system operator for E (not owned)
  MPI_Comm comm;
  int n_E;  // Size of E block (distributed DOFs)
  int n_s;  // Number of auxiliary mode DOFs

  struct AuxMode
  {
    const ComplexVector *v;       // Fourier projection vector (not owned)
    std::complex<double> g_full;  // Full DtN coefficient jλ/(μ|Γ|)
  };
  std::vector<AuxMode> aux_modes;

  mutable ComplexVector E_tmp, y_E_tmp;

public:
  FloquetAuxSystemOperator(MPI_Comm comm, const ComplexOperator *A_EE, int n_E, int n_s);

  void SetAEE(const ComplexOperator *op) { A_EE = op; }
  void ClearModes() { aux_modes.clear(); }
  void AddMode(const ComplexVector *v, std::complex<double> g_full);
  int GetNE() const { return n_E; }
  int GetNS() const { return n_s; }

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override;
};

// Sum of two ComplexOperators: (A + B)*x = A*x + B*x.
class SumComplexOperator : public ComplexOperator
{
private:
  std::unique_ptr<ComplexOperator> op_A, op_B;

public:
  SumComplexOperator(std::unique_ptr<ComplexOperator> &&A,
                     std::unique_ptr<ComplexOperator> &&B)
    : ComplexOperator(A->Height(), A->Width()), op_A(std::move(A)), op_B(std::move(B))
  {
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    op_A->Mult(x, y);
    op_B->AddMult(x, y);
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    op_A->AddMult(x, y, a);
    op_B->AddMult(x, y, a);
  }
};

//
// Data for a single Floquet port boundary.
//
class FloquetPortData
{
public:
  int excitation;
  bool active;

  [[nodiscard]] constexpr bool HasExcitation() const { return excitation != 0; }
  const auto &GetAttrList() const { return attr_list; }

  // Access the enumerated modes (for CSV column setup).
  const auto &GetModes() const { return modes; }

  // Check if mode (m, n, is_te) is part of the incident excitation at this port.
  // For linear polarization (TE/TM), only one matches. For circular (RHC/LHC), both do.
  [[nodiscard]] bool IsIncidentMode(int m, int n, bool is_te) const
  {
    if (m != 0 || n != 0)
    {
      return false;
    }
    return is_te ? (std::abs(inc_alpha_te) > 1e-14) : (std::abs(inc_alpha_tm) > 1e-14);
  }

  // Get the complex polarization coefficient for the incident mode.
  [[nodiscard]] std::complex<double> GetIncidentAlpha(bool is_te) const
  {
    return is_te ? inc_alpha_te : inc_alpha_tm;
  }

  FloquetPortData(const config::FloquetPortData &data, const IoData &iodata,
                  const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &nd_fespace);

  // Update propagation constants for given omega. Cached.
  void Initialize(double omega);

  // Static utility functions (public for testability).
  static void ComputeReciprocalLattice(const mfem::Vector &a1, const mfem::Vector &a2,
                                       mfem::Vector &b1, mfem::Vector &b2);
  static int ComputeBZOffset(const mfem::Vector &kF_unwrapped,
                             const mfem::Vector &kF_wrapped, const mfem::Vector &b,
                             double b_sq);

  // Get the dense boundary DtN operator F(omega) for this port.
  // Returns a DenseBoundaryOperator with the pre-computed N_bdr × N_bdr matrix.
  // Uses g_correction (Robin+correction mode) or g_full (FullDtN mode).
  std::unique_ptr<ComplexOperator> GetBoundaryOperator() const;

  // Compute the DtN correction coefficient g_correction for a single mode at the current
  // frequency. Returns 0 if the mode should be skipped (negligible or near-cutoff cap).
  // Used by both GetBoundaryOperator and the ROM projection.
  std::complex<double> ComputeDtNCorrectionCoeff(const FloquetMode &mode) const;

  // Compute the full DtN coefficient g_full (NOT the correction g_full - g_uniform).
  // Used by FullDtN mode and the preconditioner Schur complement assembly.
  std::complex<double> ComputeDtNFullCoeff(const FloquetMode &mode) const;

  bool UseFullDtN() const { return use_full_dtn; }
  bool UseAuxiliary() const { return use_auxiliary; }
  bool UseMassConsistent() const { return use_mass_consistent; }
  const mfem::Array<int> &GetBdrTDofList() const { return bdr_tdof_list; }
  int GetNBdrGlobal() const { return n_bdr_global; }
  MPI_Comm GetPortComm() const { return port_comm; }
  const std::vector<int> &GetRecvCounts() const { return recv_counts; }
  const std::vector<int> &GetRecvDispls() const { return recv_displs; }

  // S-parameter for all propagating orders at the current frequency.
  // If subtract_incident is true, subtracts the incident field contribution from the
  // driving port's (0,0) modes (total → scattered field conversion).
  // Returns S-parameters in the TE/TM linear basis (always).
  // Circular rotation, if needed, is applied by the caller.
  std::map<std::tuple<int, int, bool>, std::complex<double>>
  GetAllSParameters(const GridFunction &E, bool subtract_incident = false) const;

  // Get the number of propagating orders at current frequency.
  int NumPropagatingOrders() const;

  // Add excitation vector contribution directly to the RHS.
  // The incident field is normalized to inject 1 W (nondimensional) through the port,
  // consistent with lumped and wave port conventions.
  bool AddExcitationVector(double omega, ComplexVector &RHS) const;

  // Returns 1 W: the port mode is normalized such that ∫ (E_inc × H_inc⋆) · n̂ dS = 1.
  double GetExcitationPower() const { return HasExcitation() ? 1.0 : 0.0; }

  // Material properties at the port (nondimensional, from adjacent volume element).
  double GetMuEpsPort() const { return mu_eps_port; }
  double GetMuRPort() const { return mu_r_port; }
  double GetPortArea() const { return port_area; }

private:
  const MaterialOperator &mat_op;
  mfem::Array<int> attr_list;

  // Lattice and reciprocal lattice vectors.
  mfem::Vector a1, a2, b1, b2;

  // Bloch wave vector (BZ-wrapped, from MaterialOperator).
  mfem::Vector k_F;

  // Port geometry.
  mfem::Vector port_normal;
  double port_area;

  // Material properties at the port.
  double mu_eps_port;  // mu_r * eps_r (for propagation constant)
  double mu_r_port;    // mu_r (for DtN coefficient and excitation)

private:
  // Diffraction order limits.
  int max_order_m, max_order_n;

  // All Floquet modes with pre-assembled projection vectors.
  std::vector<FloquetMode> modes;

  // Full dense DtN: use dense boundary matrix instead of Robin+correction.
  bool use_full_dtn = false;

  // Auxiliary S-parameter DOF mode: couple scalar S-parameter unknowns to boundary E
  // through FE weak form. Avoids mass-vs-projection mismatch of rank-1 corrections.
  bool use_auxiliary = false;

  // Mass-consistent DtN: use per-mode boundary mass matrices A_i[j,k] = ∫(N_j·ê)(N_k·ê)dS
  // instead of rank-1 outer products conj(v)v^T. Fixes the Cauchy-Schwarz overestimate.
  bool use_mass_consistent = false;

  // Boundary DOF infrastructure (always computed, used by DenseBoundaryOperator and
  // preconditioner Schur complement).
  MPI_Comm port_comm = MPI_COMM_NULL;  // Sub-communicator for port-owning ranks
  mfem::Array<int> bdr_tdof_list;      // Local boundary true DOF indices
  int n_bdr_global = 0;                // Total boundary DOFs across port_comm
  std::vector<int> recv_counts;        // Per-rank counts for Allgatherv on port_comm
  std::vector<int> recv_displs;        // Per-rank displacements

  // FE space for re-assembling mode vectors when polarization changes with frequency.
  mfem::ParFiniteElementSpace *nd_fespace_ptr = nullptr;

  // Frequency cache.
  double omega0 = 0.0;

  // BZ wrapping offset: when MaterialOperator wraps k_F by subtracting G = bz_m*b1+bz_n*b2,
  // the Fourier projection kernel must be shifted accordingly. Mode labels remain physical.
  // The offset is frequency-dependent when FloquetReferenceFrequency is used.
  int bz_m = 0, bz_n = 0;

  // Unwrapped config k_F (in Lc units, before BZ wrapping and frequency scaling).
  // Stored for computing frequency-dependent BZ offsets.
  mfem::Vector kF_config_Lc;

  // Incident polarization coefficients: E_inc = α_TE ê_TE + α_TM ê_TM.
  // TE: (1,0). TM: (0,1). RHC: (1,j)/√2. LHC: (1,-j)/√2.
  std::complex<double> inc_alpha_te, inc_alpha_tm;

  // MPI communicator.
  MPI_Comm comm;

  // Compute the effective DtN eigenvalue and unit-power normalization for the incident
  // polarization. Used by both AddExcitationVector and GetAllSParameters.
  struct IncidentNormalization
  {
    double gamma_00;      // Propagation constant of the (0,0) mode
    double lambda_te_00;  // TE DtN eigenvalue: γ
    double lambda_tm_00;  // TM DtN eigenvalue: ω²με/γ
    double lambda_eff;    // Weighted: |α_TE|²λ_TE + |α_TM|²λ_TM
    double c_inc;         // Unit-power scale factor: 1/√(λ_eff |Γ| / (2ωμ))
  };
  IncidentNormalization ComputeIncidentNormalization(double omega) const;

  void EnumerateOrders();
  void AssembleFourierProjections(mfem::ParFiniteElementSpace &nd_fespace);
};

//
// Manager class for all Floquet port boundaries.
//
class FloquetPortOperator
{
private:
  const MaterialOperator &mat_op;
  std::map<int, FloquetPortData> ports;
  bool suppress_output = false;

public:
  FloquetPortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                      mfem::ParFiniteElementSpace &nd_fespace);

  const FloquetPortData &GetPort(int idx) const { return ports.at(idx); }
  auto begin() const { return ports.begin(); }
  auto end() const { return ports.end(); }
  auto Size() const { return ports.size(); }
  bool Empty() const { return ports.empty(); }

  void SetSuppressOutput(bool suppress) { suppress_output = suppress; }

  // Initialize all ports for the given frequency.
  void Initialize(double omega);

  // Get the combined low-rank boundary operator F(omega) for all Floquet ports.
  // Returns nullptr if no Floquet ports are configured or all ports use Auxiliary mode.
  std::unique_ptr<ComplexOperator> GetExtraSystemOperator(double omega);

  // Returns true if any active port uses Auxiliary mode.
  bool HasAuxiliaryPorts() const;

  // Total number of auxiliary mode DOFs across all active Auxiliary ports.
  int GetAuxModeCount(double omega);

  // Populate an existing FloquetAuxSystemOperator with modes for the current frequency.
  void PopulateAuxModes(double omega, FloquetAuxSystemOperator &aux_op);

  // Add Robin boundary mass coefficient (iγ₀/μ) for absorption on Floquet port faces.
  // When for_preconditioner is true, always add Robin (needed for preconditioner
  // stability). When false (system matrix), skip Robin for FullDtN ports (the dense
  // operator handles it).
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi,
                                     bool for_preconditioner = false);

  // Add excitation vector contributions for the given excitation index.
  bool AddExcitationVector(int excitation_idx, double omega, ComplexVector &RHS);

  // Collect all boundary attributes from all Floquet ports.
  mfem::Array<int> GetAttrList() const;
};

}  // namespace palace

#endif  // PALACE_MODELS_FLOQUET_PORT_OPERATOR_HPP
