// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP
#define PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP

#include <array>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "fem/substructure.hpp"
#include "linalg/operator.hpp"
#include "linalg/schursolver.hpp"
#include "linalg/vector.hpp"
#include "models/waveportoperator.hpp"
#include "utils/labels.hpp"

namespace palace
{

class IoData;

namespace config
{

struct BoundaryData;
struct DomainData;

}  // namespace config

class PortExcitations;
class SpaceOperator;

//
// Barycentric rational interpolation of a vector-valued function of the frequency from its
// samples x_j = x(ω_j), with one set of weights w for all entries,
//   x̃(ω) = Σ_j a_j(ω) x_j,   a_j(ω) = (w_j / (ω - ω_j)) / Σ_k w_k / (ω - ω_k),
// which interpolates the samples. The weights are those of the minimal rational interpolant
// (MRI) of the snapshots {x_j, iω_j x_j}, as for the error indicator of the adaptive driven
// solver (MinimalRationalInterpolation). The samples are distributed by rows over comm (the
// local rows of each rank), and the weights are replicated.
//
class BarycentricInterpolant
{
public:
  explicit BarycentricInterpolant(MPI_Comm comm = MPI_COMM_SELF) : comm(comm) {}

  // Add a sample (scaled by the caller so that its parts weigh alike), updating the
  // weights.
  void AddSample(double omega, const std::vector<std::complex<double>> &x);

  // The interpolant at ω (the local rows).
  std::vector<std::complex<double>> Evaluate(double omega) const;

  // The interpolant over an orthonormal basis q_i of the samples, x̃(ω) = Σ_i y_i(ω) q_i.
  const std::vector<std::vector<std::complex<double>>> &Basis() const { return Q; }
  std::vector<std::complex<double>> BasisCoefficients(double omega) const;

  // The frequency between the outermost samples where the denominator Σ_j w_j / (ω - ω_j)
  // is smallest (near a pole, or far from the samples): the next sample.
  double FindMaxError() const;

  const std::vector<double> &Samples() const { return z; }
  const std::vector<std::complex<double>> &Weights() const { return w; }

  // Sample j as added (the local rows), x_j = Σ_i R_ij q_i.
  std::vector<std::complex<double>> Sample(std::size_t j) const;

  // The coefficients a_j(ω) of samples at z with weights w.
  static std::vector<std::complex<double>>
  Coefficients(const std::vector<double> &z, const std::vector<std::complex<double>> &w,
               double omega);

private:
  MPI_Comm comm;
  std::vector<double> z;
  std::vector<std::vector<std::complex<double>>> Q;  // orthonormal basis of the samples
  std::vector<std::complex<double>> R;               // samples = Q R (m x m, column-major)
  std::vector<std::complex<double>> w;
};

//
// A model saved by an offline driven substructuring sweep, written and read on rank 0 (the
// header broadcast): the header (the interface signatures, the fingerprints of the
// environment and of its sources, the environment's ports, the boundary attributes of the
// region, the frequencies, the reference length of the nondimensionalization), the record
// of the environment (see EnvironmentRecord) and one record per frequency: S_E(ω), the
// condensation g of the environment's sources (per excitation with environment sources),
// and the condensed functionals of the environment's ports, whose values are V_j = c_jk +
// h_j^T u_Γ for excitation k. An exact model holds the frequencies of a uniform sweep, a
// rational one (with weights) the samples of an adaptive sweep, interpolated at any
// frequency between them (BarycentricInterpolant).
//
struct DrivenSubstructureModel
{
  static constexpr int kVersion = 13;  // of the header
  static constexpr int kSourceFp = 7;  // see DrivenSubstructure::SourceFingerprint
  // Relative tolerance of the comparisons of the environment with that of a saved model:
  // its fingerprint and element moments are sums over elements, which the partition changes
  // by rounding only.
  static constexpr double kEnvironmentTol = 1.0e-10;
  int nG = 0, sig_w = 0;
  double length_scale = 0.0;  // model.Lc of the offline run
  EnvironmentRecord environment;
  std::vector<double> signatures;                        // nG x sig_w
  std::vector<double> env_fp;                            // environment fingerprint
  std::array<std::array<double, 3>, 3> fp_centers = {};  // of its fields
  std::vector<int> excitations;  // excitations with environment sources
  std::vector<double> exc_fp;    // their source fingerprints (kSourceFp each)
  std::vector<int> ports;        // environment lumped ports (index) and wave ports (-index)
  std::vector<int> port_attributes;        // of each port: their number, then them (sorted)
  std::vector<int> region_bdr_attributes;  // see BoundaryAttributesInRegion
  std::vector<double> omega;               // frequencies of the records (nondimensional)
  std::vector<std::complex<double>> weights;  // barycentric weights (rational model)
  // The adaptive sampling of a rational model: whether it met its tolerance tol (false
  // until it ends), the error estimate of its last sample, and the least eigenvalue of
  // Im S_E / ‖S_E‖ between its samples below that at the samples (0 if passive).
  bool converged = true;
  double error = 0.0, tol = 0.0, passivity = 0.0;

  struct Record
  {
    // S (nG x nG), g (nG x excitations), h (nG x ports), c (ports x excitations), all
    // column-major.
    std::vector<std::complex<double>> S, g, h, c;
  };

  // A record as one vector, in the order of the file (S by its lower triangle), and the
  // sizes of its parts.
  std::vector<std::complex<double>> Flatten(const Record &r) const;
  std::array<std::size_t, 4> PartSizes() const;

  // The records of the model in its file (see substructuring::ModelFile), written and read
  // by the calling process (ReadHeader is collective): the interface records of a new file
  // (substructuring::InterfaceRecords: nG, sig_w, the signatures and length_scale), the
  // header and the record of the environment, in a new file or appended (superseding the
  // previous ones), and the records of the frequencies.
  void WriteHeader(const std::string &path, bool create = true) const;
  void AppendRecord(const std::string &path, const Record &r) const;
  void ReadHeader(const std::string &path, MPI_Comm comm);
  std::vector<std::complex<double>> ReadFlatRecord(const std::string &path, int j) const;
  // Number of complete records of frequencies in the file.
  int NumRecords(const std::string &path) const;

  // The positions in a model file of the payloads of its records, in order, each checked to
  // be of RecordBytes() (on every process). Collective.
  std::vector<long long> RecordPositions(const std::string &path, MPI_Comm comm) const;
  // The file cut to the header and its first n records (for a restarted sweep).
  void TruncateRecords(const std::string &path, int n) const;

  // Whether a saved header is that of this model: the same interface (and so partition),
  // environment (its fingerprint and record), sources (to source_tol, see SameSource),
  // ports, boundaries and frequencies.
  bool SameHeader(const DrivenSubstructureModel &saved, double source_tol) const;

  // Whether fingerprints agree, to a relative tolerance: the environment's (two counts,
  // then complex quadratic forms) entry by entry (FingerprintsMatch), a source's (a nonzero
  // flag, then complex pairings with the fingerprint fields) relative to its largest
  // pairing (PairingsMatch), since a pairing can vanish up to a rounding that depends on
  // the partition. The default tol suits the exact sources of lumped ports; those of wave
  // ports follow their modes, to the tolerance of their eigensolver.
  static bool SameEnvironment(const std::vector<double> &a, const std::vector<double> &b);
  static bool SameSource(const double *a, const double *b, double tol = 1.0e-10);

  // The passivity of a condensed environment S_E (n x n, complex symmetric, by its lower
  // triangle in the order of a record): the least eigenvalue of Im S_E, which is positive
  // semidefinite for a passive environment (Im u_Γ^H S_E u_Γ is the power it absorbs),
  // relative to ‖S_E‖_F. A negative value is a violation.
  static double Passivity(const std::complex<double> *lower, int n);

  // Whether S_E is passive to tol: Im S_E + tol ‖S_E‖_F I is positive definite, that is,
  // Passivity >= -tol (up to rounding), tested by a Cholesky factorization.
  static bool Passive(const std::complex<double> *lower, int n, double tol);

private:
  std::size_t RecordBytes() const;
};

//
// Exact per-frequency substructuring of a driven problem: the environment operator
// A_E(ω) = K_E + iω C_E - ω² M_E + A2_E(ω), assembled with Palace's operators restricted
// to the environment, is condensed onto the interface Γ by a partial complex symmetric
// MUMPS factorization, which gives S_E(ω). The region is factored with S_E(ω) added on Γ
// as one system. Online, the environment is not assembled: S_E(ω) and the condensation of
// the environment's sources come from a saved model, and the environment interior of the
// solution is not computed.
//
class DrivenSubstructure
{
public:
  // Online, with the saved model: on the mesh of the region (region_mesh), Γ is the set of
  // the region's boundary DOFs at the points of the saved interface DOFs. The
  // factorizations are by the sparse direct solver given (STRUMPACK, MUMPS or SUPERLU),
  // else: of the environment with its Schur complement onto Γ, by MUMPS (its partial
  // factorization), STRUMPACK or SuperLU_DIST (one solve per interface DOF), of the region
  // with S_E (one system) by STRUMPACK, MUMPS or SuperLU_DIST, in these orders among those
  // Palace is built with, with the fill-reducing ordering given.
  DrivenSubstructure(SpaceOperator &space_op, const std::vector<int> &region_attributes,
                     const std::vector<int> &environment_attributes, bool online = false,
                     const DrivenSubstructureModel *saved = nullptr,
                     bool region_mesh = false, LinearSolver solver = LinearSolver::DEFAULT,
                     SymbolicFactorization ordering = SymbolicFactorization::DEFAULT);
  ~DrivenSubstructure();

  // Condense the environment and factor the region against it at the angular frequency ω
  // (nondimensional), or only condense the environment (without region, for its condensed
  // data). The analyses of the first frequency are reused at the next ones. Online, S_E(ω)
  // is given on rank 0 (|Γ| x |Γ|, column-major, in interface order). Collective.
  void Condense(double omega, bool region = true);
  void Condense(double omega, std::vector<std::complex<double>> &&S_env);

  // Solve the coupled problem at the frequency of the last condensation for a batch of
  // right-hand sides (full true-DOF vectors): the region with the interface loads of the
  // environment sources, then the environment interior. Online, g_env on rank 0
  // (|Γ| x |rhs|, column-major) is the condensation of the environment's sources, which
  // replaces the environment part of the right-hand sides, and the environment interior
  // of u is zero. Collective.
  void Solve(const std::vector<const ComplexVector *> &rhs, std::vector<ComplexVector> &u,
             const std::vector<std::complex<double>> *g_env = nullptr);

  // Condensation onto Γ of vectors supported on the environment interior and Γ with the
  // environment factorization of the last condensation, b_Γ - A_ΓE A_EE^-1 b_E, on rank 0
  // (|Γ| x |b|, column-major), and with x the environment interior solutions A_EE^-1 b_E
  // (true DOFs, 0 elsewhere). Offline only. Collective.
  std::vector<std::complex<double>>
  CondenseEnvironment(const std::vector<const ComplexVector *> &b,
                      std::vector<ComplexVector> *x = nullptr);

  // On rank 0, after Solve: the condensation of the environment's sources (|Γ| x |rhs|,
  // offline) and the interface solution (|Γ| x |rhs|), column-major.
  const std::vector<std::complex<double>> &EnvironmentSourceCondensation() const
  {
    return g_last;
  }
  const std::vector<std::complex<double>> &InterfaceSolution() const { return u_last; }

  // S_E(ω) of the last condensation on rank 0 (|Γ| x |Γ|, column-major, complex
  // symmetric), in interface order (offline), or handed over.
  const std::vector<std::complex<double>> &Schur() const { return S; }
  std::vector<std::complex<double>> TakeSchur() { return std::move(S); }

  // Global true DOFs of the interface in interface order (replicated), and its size.
  const std::vector<HYPRE_BigInt> &InterfaceTrueDofs() const { return gamma_tdofs; }
  int InterfaceSize() const { return static_cast<int>(gamma_tdofs.size()); }

  // Local true-DOF classification: environment interior, interface (excluding the
  // Dirichlet DOFs).
  const std::vector<char> &EnvironmentInterior() const { return is_env_int; }
  const std::vector<char> &Interface() const { return is_gamma; }
  const std::vector<char> &RegionInterior() const { return is_region_int; }

  // Interface index of each local true DOF (-1 if not on Γ).
  std::vector<int> InterfaceIndex() const;

  MPI_Comm GetComm() const;

  // Local true DOFs touched by the boundary elements with the given attributes.
  std::vector<char> BoundaryTrueDofs(const std::vector<int> &bdr_attributes) const;

  // The attributes among bdr_attributes of boundary elements on Γ: each of their DOFs an
  // interface DOF or a Dirichlet one, and one at least an interface DOF. Collective.
  std::vector<int>
  InterfaceBoundaryAttributes(const std::vector<int> &bdr_attributes) const;

  // Write a mesh of the region for online runs, in mesh length units (coordinates times
  // length_scale): the region, and the environment elements on the boundary elements with
  // the given attributes (the environment's ports, whose excitations and outputs need
  // them). Returns the local environment elements written. Collective.
  std::vector<char> WriteRegionMesh(const std::string &path,
                                    const std::vector<int> &port_attributes,
                                    double length_scale) const;

  // The wave ports in the environment (each wave port lies inside one side).
  std::vector<int> EnvironmentWavePorts() const;

  // Fingerprint of the environment, independent of the region, of the partition and of the
  // DOF numbering: the global counts of environment-interior and interface DOFs, and
  // r_k^T A_E(ω) r_k for three fixed fields r_k, applied with the partially assembled
  // environment operators. Collective.
  std::vector<double> EnvironmentFingerprint(double omega) const;

  // The centers of the fields r_k: those of the saved model, else the FingerprintCenters of
  // the environment.
  const std::array<std::array<double, 3>, 3> &FingerprintCenters() const
  {
    return fp_centers;
  }

  // Fingerprint of the environment part b_E of a right-hand side, independent of the
  // partition (which changes the basis of some Nédélec DOFs): whether b_E is nonzero, and
  // r_k^T b_E for the three fixed fields r_k (complex). Collective.
  std::vector<double> SourceFingerprint(const ComplexVector &b) const;

  // Whether a right-hand side has entries on Γ, which a model cannot hold: they would come
  // from the elements of both sides, and from the offline run. Collective.
  bool OnInterface(const ComplexVector &b) const;

private:
  SpaceOperator &space_op;
  LinearSolver solver;
  SymbolicFactorization ordering;
  bool online, region_only = false;
  std::vector<char> is_env_int, is_gamma, is_region_int;
  std::vector<HYPRE_BigInt> gamma_tdofs;
  std::vector<int> gamma_cnt, gamma_disp;  // interface DOFs per rank

  // The operator K + iω C - ω² M + A2(ω) of one side (environment or region), factored on
  // the side's interior and Γ: the lower-triangle pattern and the entries of the
  // frequency-independent parts in pattern order (the assembled operators are not kept).
  struct Side
  {
    std::vector<int> attrs;
    // Local true DOFs pinned in the side's operator: outside its interior and Γ.
    mfem::Array<int> pinned;
    // The factored system, on the side's interior and Γ: its local true DOFs (the local
    // rows), global size, and its indices of the interface DOFs.
    std::vector<int> rows;
    HYPRE_BigInt n_sys = 0;
    std::vector<HYPRE_BigInt> gamma_sys;
    // Lower-triangle pattern of the local rows (row_ptr over the local true DOFs, global
    // columns jcn, 1-based), with the entries of the parts on it.
    std::vector<int> row_ptr, jcn;
    std::vector<std::vector<double>> parts;  // Kr, Ki, Cr, Ci, Mr, Mi (empty if absent)
    bool extra = false;                      // A2(ω) is present
    // The local entries (1-based COO in the system's indices) at the first frequency,
    // until factored.
    std::vector<int> irn_sys, jcn_sys;
    std::vector<std::complex<double>> val;
    std::unique_ptr<SchurFactorization<std::complex<double>>> schur;
    // The rank-one modal terms g s sᵀ of the side's wave ports, bordering the system: an
    // extra unknown y per term after its n_dofs DOFs, on the last rank, with the entries
    // √(gσ) s against the term's DOFs (local true DOFs, in the system) and -σ on its
    // diagonal (eliminating y adds g s sᵀ), σ = border. The entries follow those of the
    // pattern.
    std::vector<int> wave_ports, term_ports;
    std::vector<std::vector<int>> term_dofs;
    HYPRE_BigInt n_dofs = 0;
    double border = 1.0;
  };
  Side env, region;
  bool solvable = false;  // the region is factored against the environment
  // The region with S_E on Γ, factored as one system; its entries are the region's, then
  // (on rank 0) S_E's lower triangle.
  std::unique_ptr<SchurFactorization<std::complex<double>>> coupled;
  std::size_t coupled_region_entries = 0;
  // On rank 0: S_E (offline), the environment's source condensation and the interface
  // solution of the last Solve.
  std::vector<std::complex<double>> S, g_last, u_last;

  // The local true DOFs of the wave ports (each in the wave_ports of its side).
  std::map<int, std::vector<char>> wave_port_dofs;

  // The fixed fields r_k of the fingerprints (true DOFs), projected on first use.
  std::array<std::array<double, 3>, 3> fp_centers;
  mutable std::array<Vector, 3> fp_fields;
  const Vector &FingerprintField(int k) const;

  // The real and imaginary parts of the stiffness (op = 0), damping (1) or mass (2) matrix
  // of a side, assembled on its elements, with its pinned DOFs eliminated (null if zero).
  std::array<std::unique_ptr<mfem::HypreParMatrix>, 2> SideParts(const Side &side,
                                                                 int op) const;

  // Factor the region with S_E (on rank 0) added on Γ as one system, releasing S_E once its
  // entries are in the factored values if release (set_up: Setup of the region was called
  // for its first factorization).
  void FactorCoupled(double omega, std::vector<std::complex<double>> &S_env, bool release,
                     bool set_up);

  // The region and Γ from the coupled factorization, with g_env (rank 0) on Γ.
  void SolveCoupled(const std::vector<const ComplexVector *> &rhs,
                    std::vector<ComplexVector> &u,
                    const std::vector<std::complex<double>> &g_env);

  // The rank-one modal terms of a side's wave ports at ω with their ports, in a fixed
  // order.
  std::vector<std::pair<int, WavePortOperator::ModalCorrectionTerm>>
  ModalTerms(const Side &side, double omega);

  // The assembled parts of A2(ω) of a side, with the pinned DOFs eliminated.
  std::vector<std::unique_ptr<mfem::HypreParMatrix>> ExtraParts(const Side &side,
                                                                double omega);

  // The pattern and entries of a side at the first frequency.
  void Setup(Side &side, double omega);

  // Factor a side (set up at ω), or refactor it at another frequency.
  void Factor(Side &side, double omega);

  // The entries of a side at ω on its pattern, given its A2(ω) and modal terms.
  void Fill(const Side &side, double omega,
            const std::vector<std::unique_ptr<mfem::HypreParMatrix>> &extra,
            const std::vector<std::pair<int, WavePortOperator::ModalCorrectionTerm>> &terms,
            std::vector<std::complex<double>> &val) const;
};

// The lumped ports on the environment side (and its wave ports, as -index), with the true
// DOFs of its lumped ports and the attributes of each (see
// DrivenSubstructureModel::port_attributes), aborting for a lumped port touching Γ (its
// excitation and voltage would depend on both sides). Collective.
std::vector<int> EnvironmentPorts(const IoData &iodata, const DrivenSubstructure &ds,
                                  std::map<int, std::vector<char>> &port_dofs,
                                  std::vector<int> &port_attributes);

// The full S_E (n x n, row-major) of a record's lower triangle (by columns), mapped to the
// current interface as M^-T S_E M^-1 with map, if given (SignatureMap::DualMatrix).
std::vector<std::complex<double>> FullSchur(const std::complex<double> *lower, int n,
                                            const SignatureMap *map);

// The index of x in y (relative tolerance tol), -1 if none.
int FindValue(const std::vector<double> &y, double x, double tol = 1.0e-10);

// The tolerance of the source fingerprint of an excitation (see
// DrivenSubstructureModel::SameSource): the modes of the environment's wave ports, and so
// their sources, are reproducible to their eigensolver tolerance over the spectral gap
// only.
double SourceTolerance(const IoData &iodata, const PortExcitations &port_excitations,
                       const std::vector<int> &env_ports, int excitation_idx);

// The postprocessing of a region whose field is known in the region and on Γ only (from a
// model, online or adaptive): the energies of the domain postprocessing in the region,
// without totals or participation ratios, the surface fluxes on boundaries in the region
// (on the mesh of the region, region_bdr_attributes, those the model records), and the
// probes found in its elements, with a warning (what: "Online" or "Adaptive") naming those
// left out, and the power of the environment's ports if it has some. Collective.
void RegionPostprocessing(const std::vector<int> &region_attributes,
                          const mfem::ParMesh &mesh,
                          const std::vector<int> *region_bdr_attributes, bool env_ports,
                          const char *what, config::DomainData &domains,
                          config::BoundaryData &boundaries);

}  // namespace palace

#endif  // PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP
