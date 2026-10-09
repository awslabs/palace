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
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/waveportoperator.hpp"

namespace palace
{

template <typename T>
class MumpsSchurSolverT;
class SpaceOperator;

//
// The saved environment model of a driven substructuring sweep: the interface signatures,
// the fingerprints of the environment and of its sources, the environment's lumped ports,
// the frequencies, and one record per frequency: S_E(ω), the condensation of the
// environment's sources (per excitation with environment sources) and the condensed
// voltage functionals of the environment's lumped ports, V_j = c_jk + h_j^T u_Γ for
// excitation k. Written and read on rank 0 (the header is broadcast).
//
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
// A model saved by an offline driven substructuring sweep: the header (interface
// signatures, fingerprints, environment ports) and one record per frequency. An exact model
// holds the frequencies of a uniform sweep. A rational one (with weights) holds the samples
// of an adaptive sweep, interpolated at any frequency between them
// (BarycentricInterpolant).
//
struct DrivenSubstructureModel
{
  static constexpr int kSourceFp = 7;  // see DrivenSubstructure::SourceFingerprint
  int version = 2, nG = 0, sig_w = 0;
  std::vector<double> signatures;  // nG x sig_w
  std::vector<double> env_fp;      // environment fingerprint
  std::vector<int> excitations;    // excitations with environment sources
  std::vector<double> exc_fp;      // their source fingerprints (kSourceFp each)
  std::vector<int> ports;     // environment lumped ports (index) and wave ports (-index)
  std::vector<double> omega;  // frequencies of the records (nondimensional)
  std::vector<std::complex<double>> weights;  // barycentric weights (rational model)
  int capacity = 0;  // records the header has room for (at least omega.size())

  struct Record
  {
    // S (nG x nG), g (nG x excitations), h (nG x ports), c (ports x excitations), all
    // column-major.
    std::vector<std::complex<double>> S, g, h, c;
  };

  // A record as one vector, in the order of the file (S by its lower triangle), and the
  // sizes of its parts.
  std::vector<std::complex<double>> Flatten(const Record &r) const;
  Record Unflatten(const std::vector<std::complex<double>> &x) const;
  std::array<std::size_t, 4> PartSizes() const;

  // The header, in a new file or rewritten in place (with the same capacity).
  void WriteHeader(const std::string &path, bool create = true) const;
  void AppendRecord(const std::string &path, const Record &r) const;
  void ReadHeader(const std::string &path, MPI_Comm comm);
  std::vector<std::complex<double>> ReadFlatRecord(const std::string &path, int j) const;
  Record ReadRecord(const std::string &path, int j) const
  {
    return Unflatten(ReadFlatRecord(path, j));
  }
  // Number of complete records in the file.
  int NumRecords(const std::string &path) const;

  // Whether fingerprints agree, to a relative tolerance. The environment's (two counts,
  // then complex quadratic forms) entry by entry. A source's (a nonzero flag, then complex
  // pairings with the fingerprint fields, kSourceFp entries) relative to its largest
  // pairing: a pairing can vanish up to rounding, which depends on the partition (for a
  // source in a plane where a fingerprint field is normal to it).
  static bool SameEnvironment(const std::vector<double> &a, const std::vector<double> &b);
  static bool SameSource(const double *a, const double *b);

  // The passivity of a condensed environment S_E (n x n, complex symmetric, by its lower
  // triangle in the order of a record): the least eigenvalue of Im S_E, which is positive
  // semidefinite for a passive environment (Im u_Γ^H S_E u_Γ is the power it absorbs),
  // relative to ‖S_E‖_F. A negative value is a violation.
  static double Passivity(const std::complex<double> *lower, int n);

private:
  std::size_t HeaderBytes() const;
  std::size_t RecordBytes() const;
};

//
// Exact per-frequency substructuring of a driven problem: the environment operator
// A_E(ω) = K_E + iω C_E - ω² M_E + A2_E(ω), assembled with Palace's operators restricted
// to the environment, is condensed onto the interface Γ by a partial complex symmetric
// MUMPS factorization, which gives S_E(ω). The region is condensed onto Γ in the same way,
// and the interface system S_R(ω) + S_E(ω) is factored densely. Online, the environment is
// not assembled: S_E(ω) and the condensation of the environment's sources come from a
// saved model, and the environment interior of the solution is not computed.
//
class DrivenSubstructure
{
public:
  DrivenSubstructure(SpaceOperator &space_op, const std::vector<int> &region_attributes,
                     const std::vector<int> &environment_attributes, bool online = false);
  ~DrivenSubstructure();

  // Condense the environment and factor the region against it at the angular frequency ω
  // (nondimensional). The analyses of the first frequency are reused at the next ones.
  // Online, S_E(ω) is given on rank 0 (|Γ| x |Γ|, column-major, in interface order).
  // Collective.
  void Condense(double omega);
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
  // (|Γ| x |b|, column-major). Offline only. Collective.
  std::vector<std::complex<double>>
  CondenseEnvironment(const std::vector<const ComplexVector *> &b);

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

  // The wave ports in the environment (each wave port lies inside one side).
  std::vector<int> EnvironmentWavePorts() const;

  // A reduced model of the region against a given environment (online): Galerkin projection
  // onto the real span of field snapshots, as in the adaptive driven solver. The
  // environment enters through its reduced terms on Γ, which the caller forms from the
  // interface rows of the basis.
  //
  // Add the real and imaginary parts of the fields u (true DOFs, zero in the environment
  // interior) to the orthonormal basis (collective).
  void AddReducedBasis(const std::vector<ComplexVector> &u);
  int ReducedDimension() const { return static_cast<int>(red_basis.size()); }

  // The interface rows of the basis on rank 0 (|Γ| x r, column-major, in interface order).
  const std::vector<double> &ReducedInterfaceBasis() const { return red_gamma; }

  // Solve the reduced system at ω for the region parts of rhs, with the environment's
  // terms A_env (r x r) and b_env (r x |rhs|), column-major on rank 0, and expand: u = V y,
  // and the interface solution (InterfaceSolution) on rank 0. Collective.
  void SolveReduced(double omega, const std::vector<const ComplexVector *> &rhs,
                    const std::vector<std::complex<double>> &A_env,
                    const std::vector<std::complex<double>> &b_env,
                    std::vector<ComplexVector> &u);

  // Fingerprint of the environment, independent of the region, of the partition and of the
  // DOF numbering: the global counts of environment-interior and interface DOFs, and
  // r_k^T A_E(ω) r_k for three fixed fields r_k, applied with the partially assembled
  // environment operators. Collective.
  std::vector<double> EnvironmentFingerprint(double omega) const;

  // Fingerprint of the environment part b_E of a right-hand side, independent of the
  // partition (which changes the basis of some Nédélec DOFs): whether b_E is nonzero, and
  // r_k^T b_E for the three fixed fields r_k (complex). Collective.
  std::vector<double> SourceFingerprint(const ComplexVector &b) const;

private:
  SpaceOperator &space_op;
  bool online;
  std::vector<char> is_env_int, is_gamma, is_region_int;
  std::vector<HYPRE_BigInt> gamma_tdofs;
  std::vector<int> gamma_cnt, gamma_disp;  // interface DOFs per rank

  // The operator K + iω C - ω² M + A2(ω) of one side (environment or region), with its
  // pinned DOFs, factored by MUMPS: the lower-triangle pattern, the entries of the
  // frequency-independent parts in pattern order (the assembled operators are not kept),
  // and the positions of the unit diagonal of the pinned DOFs.
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
    std::unique_ptr<MumpsSchurSolverT<std::complex<double>>> schur;
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
  // On rank 0: S_E, and the factored interface system S_R + S_E with its pivots; the
  // environment's source condensation and the interface solution of the last Solve.
  std::vector<std::complex<double>> S, T, g_last, u_last;
  std::vector<int> T_piv;

  // The local true DOFs of the wave ports (each in the wave_ports of its side).
  std::map<int, std::vector<char>> wave_port_dofs;

  // The reduced region: its operator parts (Kr, Ki, Cr, Ci, Mr, Mi on the region, with the
  // pinned DOFs eliminated), the basis V, the parts' projections V^T P V (r x r,
  // column-major), and the interface rows of V on rank 0.
  std::vector<std::unique_ptr<mfem::HypreParMatrix>> red_ops;
  std::vector<Vector> red_basis;
  std::vector<std::vector<double>> red_parts;
  std::vector<double> red_gamma;

  // Factor S_R + S_E on rank 0.
  void FactorInterface();

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

}  // namespace palace

#endif  // PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP
