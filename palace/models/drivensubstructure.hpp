// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP
#define PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP

#include <array>
#include <complex>
#include <map>
#include <memory>
#include <optional>
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

class BlockCyclicGrid;
class PortExcitations;
class SpaceOperator;

//
// Barycentric rational interpolation of a vector-valued function of the frequency from its
// samples x_j = x(ω_j), with one set of weights w for all entries,
//   x̃(ω) = Σ_j a_j(ω) x_j,   a_j(ω) = (w_j / (ω - ω_j)) / Σ_k w_k / (ω - ω_k),
// which interpolates the samples. The weights are those of the minimal rational interpolant
// (MRI) of the snapshots {x_j, iω_j x_j}, as for the error indicator of the adaptive driven
// solver (MinimalRationalInterpolation). The entries of the samples are distributed over
// comm (each process holds its entries of every sample), and the weights are replicated.
//
class BarycentricInterpolant
{
public:
  explicit BarycentricInterpolant(MPI_Comm comm = MPI_COMM_SELF) : comm(comm) {}

  // Add a sample (scaled by the caller so that its parts weigh alike), updating the
  // weights.
  void AddSample(double omega, const std::vector<std::complex<double>> &x);

  // The interpolant at ω (the entries of this process).
  std::vector<std::complex<double>> Evaluate(double omega) const;

  // The frequency between the outermost samples where the denominator Σ_j w_j / (ω - ω_j)
  // is smallest (near a pole, or far from the samples): the next sample.
  double FindMaxError() const;

  const std::vector<double> &Samples() const { return z; }
  const std::vector<std::complex<double>> &Weights() const { return w; }

  // Sample j as added (the entries of this process), x_j = Σ_i R_ij q_i.
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

  // The coefficients y_i(ω) of the interpolant over the basis q_i of the samples,
  // x̃(ω) = Σ_i y_i(ω) q_i.
  std::vector<std::complex<double>> BasisCoefficients(double omega) const;
};

//
// A model saved by an offline driven substructuring sweep, written and read on rank 0 (the
// header broadcast): the header (the interface signatures, the fingerprints of the
// environment and of its sources, the environment's ports, the frequencies, the reference
// length of the nondimensionalization), the record
// of the environment (see EnvironmentRecord) and one record per frequency: S_E(ω), the
// condensation g of the environment's sources (per excitation with environment sources),
// and the condensed functionals of the environment's ports, whose values are V_j = c_jk +
// h_j^T u_Γ for excitation k. An exact model holds the frequencies of a uniform sweep, a
// rational one (with weights) the samples of an adaptive sweep, interpolated at any
// frequency between them (BarycentricInterpolant).
//
struct DrivenSubstructureModel
{
  static constexpr int kVersion = 14;  // of the model (a new header record leaves it)
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
  std::vector<int> port_attributes;  // of each port: their number, then them (sorted)
  std::vector<double> omega;         // frequencies of the records (nondimensional)
  std::vector<std::complex<double>> weights;  // barycentric weights (rational model)
  // The adaptive sampling of a rational model: whether it met its tolerance tol (false
  // until it ends), the error estimate of its last sample, and the least eigenvalue of
  // Im S_E / ‖S_E‖ between its samples below that at the samples (0 if passive).
  bool converged = true;
  double error = 0.0, tol = 0.0, passivity = 0.0;

  struct Record
  {
    // S (nG x nG, complex symmetric) by the share of this process of its lower triangle
    // (see TriangleShare), and on rank 0 g (nG x excitations), h (nG x ports) and c (ports
    // x excitations), column-major.
    std::vector<std::complex<double>> S, g, h, c;
  };

  // A record as the vector of this process (see EnvironmentOperator::RecordFunction): its
  // share of S_E, then the parts after it, handed over. The sizes of the parts of a record
  // in the order of the file: S by the columns of its lower triangle, then g, h and c.
  static std::vector<std::complex<double>> Flatten(Record &&r);
  std::array<std::size_t, 4> PartSizes() const;

  // The records of the model in its file (see substructuring::ModelFile), written and read
  // by the calling process (AppendRecord and ReadHeader are collective): the interface
  // records of a new file (substructuring::InterfaceRecords: nG, sig_w, the signatures and
  // length_scale), the header and the record of the environment, in a new file or appended
  // (superseding the previous ones), and the records of the frequencies, which rank 0
  // writes, S from the shares of the processes in chunks.
  void WriteHeader(const std::string &path, bool create = true) const;
  void AppendRecord(const std::string &path, const Record &r, MPI_Comm comm) const;
  void ReadHeader(const std::string &path, MPI_Comm comm);
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

  // The passivity of a condensed environment S_E (n x n, complex symmetric), given by the
  // share of this process of its lower triangle (see TriangleShare, over the communicator
  // of the grid): the least eigenvalue of Im S_E, which is positive semidefinite for a
  // passive environment (Im u_Γ^H S_E u_Γ is the power it absorbs), relative to ‖S_E‖_F,
  // computed on the grid. A negative value is a violation. Collective.
  static double Passivity(const std::complex<double> *share, int n,
                          const BlockCyclicGrid &grid);

  // Whether S_E is passive to tol: Im S_E + tol ‖S_E‖_F I is positive definite, that is,
  // Passivity >= -tol (up to rounding), tested by a Cholesky factorization on the grid.
  // Collective.
  static bool Passive(const std::complex<double> *share, int n, double tol,
                      const BlockCyclicGrid &grid);

private:
  std::size_t RecordBytes() const;
};

//
// The environment of a driven substructuring, condensed at a frequency: the environment
// operator A_E(ω) = K_E + iω C_E - ω² M_E + A2_E(ω), assembled with Palace's operators
// restricted to the environment, is condensed onto the interface Γ by a partial complex
// symmetric factorization, which gives S_E(ω) and the condensations of vectors of the
// environment (see CondenseEnvironment). Online, with a saved model, the environment is
// not assembled: the class classifies the DOFs of Γ and of the two sides, and gives the
// fingerprints and the port and boundary queries of the run (see EnvironmentOperator).
//
class DrivenSubstructure
{
public:
  // Offline: the environment is condensed (see Condense), factored with its Schur
  // complement onto Γ by the sparse direct solver given (STRUMPACK, MUMPS or SUPERLU), else
  // by MUMPS, STRUMPACK or SuperLU_DIST (see SchurSolverFor), with the fill-reducing
  // ordering given.
  DrivenSubstructure(SpaceOperator &space_op, const std::vector<int> &region_attributes,
                     const std::vector<int> &environment_attributes,
                     LinearSolver solver = LinearSolver::DEFAULT,
                     SymbolicFactorization ordering = SymbolicFactorization::DEFAULT);

  // Online, against a model of the environment (saved, or sampled in this run): its
  // fingerprint centers, and on the mesh of the region (region_mesh) Γ is the set of the
  // region's boundary DOFs at the points of its interface DOFs. The environment is not
  // condensed.
  DrivenSubstructure(SpaceOperator &space_op, const std::vector<int> &region_attributes,
                     const std::vector<int> &environment_attributes,
                     const DrivenSubstructureModel &model, bool region_mesh);
  ~DrivenSubstructure();

  // Condense the environment at the angular frequency ω (nondimensional): S_E(ω) (see
  // TakeSchur), and its factorization for CondenseEnvironment. The analysis of the first
  // frequency is reused at the next ones. Offline only. Collective.
  void Condense(double omega);

  // Condensation onto Γ of vectors supported on the environment interior and Γ with the
  // environment factorization of the last condensation, b_Γ - A_ΓE A_EE^-1 b_E, on rank 0
  // (|Γ| x |b|, column-major), and with x the environment interior solutions A_EE^-1 b_E
  // (true DOFs, 0 elsewhere). Offline only. Collective.
  std::vector<std::complex<double>>
  CondenseEnvironment(const std::vector<const ComplexVector *> &b,
                      std::vector<ComplexVector> *x = nullptr);

  // S_E(ω) of the last condensation (complex symmetric, in interface order), by the share
  // of this process of its lower triangle (see TriangleShare), handed over.
  std::vector<std::complex<double>> TakeSchur() { return std::move(S); }

  // The interface Γ, in its order, with the classes of the local true DOFs: the
  // environment interior, Γ and the region interior, none of them Dirichlet.
  const SubstructureInterface &GetInterface() const { return iface; }

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
  std::vector<int> region_attrs;
  SubstructureInterface iface;

  DrivenSubstructure(SpaceOperator &space_op, const std::vector<int> &region_attributes,
                     const std::vector<int> &environment_attributes,
                     const DrivenSubstructureModel *model, bool region_mesh,
                     LinearSolver solver, SymbolicFactorization ordering);

  // The operator K + iω C - ω² M + A2(ω) of the environment, factored on its interior and
  // Γ: the lower-triangle pattern and the entries of the
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
  Side env;
  // The share of this process of the lower triangle of S_E of the last condensation, and
  // the grid of the interface, of all the processes, on which the environment is condensed
  // (offline).
  std::vector<std::complex<double>> S;
  std::unique_ptr<BlockCyclicGrid> gamma_grid;

  // The local true DOFs of the wave ports (each in the wave_ports of its side).
  std::map<int, std::vector<char>> wave_port_dofs;

  // The fixed fields r_k of the fingerprints (true DOFs), projected on first use, on every
  // process alike (also one without true DOFs, whose fields are empty).
  std::array<std::array<double, 3>, 3> fp_centers;
  mutable std::array<std::optional<Vector>, 3> fp_fields;
  const Vector &FingerprintField(int k) const;

  // The real and imaginary parts of the stiffness (op = 0), damping (1) or mass (2) matrix
  // of a side, assembled on its elements, with its pinned DOFs eliminated (null if zero).
  std::array<std::unique_ptr<mfem::HypreParMatrix>, 2> SideParts(const Side &side,
                                                                 int op) const;

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

// The lumped ports on the environment side (and its wave ports, as -index), with the
// attributes of each (see DrivenSubstructureModel::port_attributes), aborting for a lumped
// port touching Γ (its excitation and voltage would depend on both sides). Collective.
std::vector<int> EnvironmentPorts(const IoData &iodata, const DrivenSubstructure &ds,
                                  std::vector<int> &port_attributes);

// The index of x in y (relative tolerance tol), -1 if none.
int FindValue(const std::vector<double> &y, double x, double tol = 1.0e-10);

// The tolerance of the source fingerprint of an excitation (see
// DrivenSubstructureModel::SameSource): the modes of the environment's wave ports, and so
// their sources, are reproducible to their eigensolver tolerance over the spectral gap
// only.
double SourceTolerance(const IoData &iodata, const PortExcitations &port_excitations,
                       const std::vector<int> &env_ports, int excitation_idx);

}  // namespace palace

#endif  // PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP
