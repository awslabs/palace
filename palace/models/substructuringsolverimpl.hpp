// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SUBSTRUCTURING_SOLVER_IMPL_HPP
#define PALACE_MODELS_SUBSTRUCTURING_SOLVER_IMPL_HPP

#include "models/substructuringsolver.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <set>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fem/substructure.hpp"
#include "fem/substructuremodel.hpp"
#include "linalg/blockcyclic.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/schurfactory.hpp"
#include "linalg/schursolver.hpp"
#include "linalg/vector.hpp"
#include "models/curlcurloperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/regionsolver.hpp"
#include "models/superconductorsheetoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace substructuring
{

// Materialized DtN: applies the precomputed interface operator S_E to a distributed
// interface vector (gathered by an Allreduce): each rank multiplies its own dense row
// block.
class MaterializedDtN : public mfem::Operator
{
public:
  MaterializedDtN(const std::vector<double> &S_rows, int row_off,
                  const std::vector<int> &gamma_global, int nG_global, MPI_Comm comm)
    : mfem::Operator(static_cast<int>(gamma_global.size())), S_rows(S_rows),
      row_off(row_off), gamma_global(gamma_global), nG_global(nG_global), comm(comm),
      xg(nG_global)
  {
    // Each rank owns a contiguous block of the global interface order, in rank order.
    const int nloc = static_cast<int>(std::count_if(
        gamma_global.begin(), gamma_global.end(), [](int g) { return g >= 0; }));
    cnt.resize(Mpi::Size(comm));
    disp.assign(cnt.size(), 0);
    MPI_Allgather(&nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    for (std::size_t r = 1; r < cnt.size(); r++)
    {
      disp[r] = disp[r - 1] + cnt[r - 1];
    }
    xo.resize(nloc);
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    const double *xx = x.HostRead();
    for (int i = 0; i < height; i++)
    {
      if (gamma_global[i] >= 0)
      {
        xo[gamma_global[i] - row_off] = xx[i];
      }
    }
    MPI_Allgatherv(xo.data(), static_cast<int>(xo.size()), MPI_DOUBLE, xg.data(),
                   cnt.data(), disp.data(), MPI_DOUBLE, comm);
    double *yy = y.HostWrite();
    std::fill(yy, yy + height, 0.0);
    for (int i = 0; i < height; i++)
    {
      if (gamma_global[i] >= 0)
      {
        const int r = gamma_global[i] - row_off;
        const double *row = &S_rows[static_cast<std::size_t>(r) * nG_global];
        double s = 0.0;
        for (int j = 0; j < nG_global; j++)
        {
          s += row[j] * xg[j];
        }
        yy[i] = s;
      }
    }
  }

private:
  const std::vector<double> &S_rows;
  int row_off;
  const std::vector<int> &gamma_global;
  int nG_global;
  MPI_Comm comm;
  std::vector<int> cnt, disp;
  mutable std::vector<double> xo, xg;
};

}  // namespace substructuring

// Region/environment operators are those of Palace's regular solve, restricted to the
// domains of each side; the interface is identified in true-DOF space and the environment
// condensed onto it (a distributed dense DtN).
struct SubstructuringSolver::Impl
{
  const IoData &iodata;
  Mesh &mesh;
  mfem::ParMesh &parent;
  bool magnetostatic;
  // The head of a model, from rank 0: its format version (-1: the file cannot be read, 0:
  // not a model; a container of another version reads as another version), order,
  // interface size, signature type and width (0: another version, or an incomplete model),
  // its records, the interface signatures, the position of S_E, the record of the
  // environment, magnetostatic regularization, and sides of the sheet faces on Gamma (0 or
  // empty: not saved). Collective.
  struct SavedHead
  {
    int version = -1, order = 0, nG = 0, sig_type = 0, width = 0;
    substructuring::ModelFile file;
    std::vector<double> sig;
    EnvironmentRecord environment;
    double regularization = 0.0;
    std::vector<double> sheet_sides;
    long long S_pos = 0;
  };
  static SavedHead ReadSavedHead(const std::string &path, MPI_Comm comm);
  // The head of the model of an online run, checked against the run (empty offline).
  // Collective.
  static SavedHead ReadOnlineHead(const IoData &iodata, bool magnetostatic, MPI_Comm comm);
  SavedHead saved;
  // Online on a mesh of the region (written offline next to the model, whose environment
  // elements it holds, which tell it from the full mesh): Γ also holds the region's
  // boundary DOFs at the saved interface DOFs, and the environment is that of the model.
  bool region_only;
  // The operator of the regular electrostatic or magnetostatic solve on the finest mesh,
  // which assembles the operators of the sides (see GetTerm), and its (unbroken)
  // finite element space.
  std::unique_ptr<LaplaceOperator> laplace_op;
  std::unique_ptr<CurlCurlOperator> curlcurl_op;
  FiniteElementSpace &parent_space;
  mfem::ParFiniteElementSpace &parent_fes;
  // Magnetostatics with two-sided superconductor films: the solution space of the operator,
  // broken across the films, on which the operators and true-DOF vectors act (see
  // SolutionSpace()), its maps to the parent true DOFs, and the parent true DOF of each of
  // its true DOFs.
  const FiniteElementSpace *broken = nullptr;
  std::unique_ptr<mfem::HypreParMatrix> first_version, all_versions;
  std::vector<int> version_parent;
  int nt;

  std::unique_ptr<mfem::HypreParMatrix> A_region;
  mfem::Array<int> ra_arr, ea_arr;
  // The grid of Gamma, for the dense interface system of the region solve, and the region
  // solver.
  std::unique_ptr<BlockCyclicGrid> gamma_grid;
  std::unique_ptr<RegionSolver<double>> region_solver;
  // Local true DOFs: region-free (the region solve), of region elements (where the region
  // operators act), and interface.
  mfem::Array<int> rf_idx, rc_idx, gamma_idx;
  // The environment's side (its interior and Gamma) with Gamma as the Schur variables,
  // which gives S_E and every environment interior solve, and how it is factored.
  std::unique_ptr<SchurFactorization<double>> env_schur;
  std::string env_factorization;
  // Block size of the batched multi-RHS environment solves.
  static constexpr int kMaterializeBlock = 32;
  // Relative tolerance of the fingerprint checks of a saved model: the fingerprints are
  // sums over elements, which the partition changes by rounding only.
  static constexpr double kFingerprintTol = 1.0e-9;
  // Region energy operator: the solve operator (electrostatics), the volume terms
  // (magnetostatics).
  std::unique_ptr<mfem::HypreParMatrix> A_region_energy;
  mfem::HypreParMatrix *K_region_e = nullptr;
  std::vector<char> is_gamma, is_env_int, is_region_free;
  mfem::Array<int> dbc_tdofs;
  std::vector<char> env_dbc;  // the Dirichlet true DOFs of environment elements
  std::map<int, std::vector<int>> terminal_tdofs;  // terminal index -> its true DOFs
  std::vector<char> region_dbc;                    // Dirichlet DOFs of region elements
  // The domain attributes of the environment: those of the configuration, or by default the
  // attributes of the mesh not in the region.
  std::vector<int> env_attributes;

  // Interface operator S_E over a global interface enumeration, distributed by rows: each
  // rank stores its own rows [gamma_off, gamma_off + gamma_nloc) (row-major, x nG_global).
  std::vector<int> gamma_global;  // owned parent true DOF -> global interface index, or -1
  int nG_global = 0;
  std::vector<double> S_rows;
  int gamma_off = 0, gamma_nloc = 0;
  std::unique_ptr<substructuring::MaterializedDtN> mat_dtn;

  // Dirichlet-lift modes (electrostatics): the environment condensed onto Gamma + a set of
  // Dirichlet lifts x_k (the terminal unit potentials). With t_k = A_env x_k and
  // w_k = A_EE^-1 t_k|_E: the region-solve interface load g_k = (t_k - A_env w_k)|_Gamma
  // (G_rows, owned interface rows: G_rows[r * K + k]), and the environment's exact
  // condensed energy [u_G; V]^T [S_E, G; G^T, Cmode] [u_G; V] with the replicated Cmode. A
  // partition-invariant fingerprint per mode (LiftFingerprint) guards reuse.
  std::vector<int> mode_ids;
  std::vector<double> mode_fp;  // 2 per mode
  std::vector<double> G_rows, Cmode;
  bool modes_ready = false;
  // Magnetostatic energy interface operator S^K (distributed like S_rows), materialized
  // with a saved model; otherwise the energies use environment solves (electrostatics: S^K
  // = S_E).
  std::vector<double> SK_rows;
  // Rank-independent presence flags (the distributed arrays are empty on ranks that own no
  // interface DOFs, so their emptiness must not gate collectives).
  bool have_sk = false;
  std::unique_ptr<substructuring::MaterializedDtN> mat_dtn_K;
  // Kinetic energy interface operator S^Kin = P^T M_kin^E P of the London sheets and
  // domains of the environment, M_kin^E = M_sheet^E + M_London^E (distributed like S_rows),
  // materialized with S^K when the environment has them (EnvKinetic); the region's kinetic
  // energy is that of its solution.
  std::vector<double> SKin_rows;
  bool have_skin = false;
  std::unique_ptr<substructuring::MaterializedDtN> mat_dtn_Kin;
  bool env_built = false;  // environment interior solver set up

  Impl(const IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh);

  // A term of the operator of the solves on the elements with the given domain attributes
  // and the boundary elements bounding them (all, without), from the operator of the
  // regular solve: the volume terms (curl-curl, with the London mass of the London domains,
  // or Laplace), the superconductor sheet terms (magnetostatics; null without sheets), the
  // London mass alone (magnetostatics; null without London domains), or the regularization
  // eps M (magnetostatics, with attrs). Assembled, or applied matrix-free. Collective.
  enum class Term
  {
    VOLUME,
    SHEETS,
    LONDON,
    MASS
  };
  std::unique_ptr<ParOperator> GetTerm(const std::vector<int> *attrs, Term term,
                                       bool assemble) const;

  // The assembled term (null if there is none), symmetric on a broken space. Collective.
  std::unique_ptr<mfem::HypreParMatrix> AssembleTerm(const std::vector<int> *attrs,
                                                     Term term) const;

  // The region operators: of the solves, of the energies, and the sheet terms. Collective.
  void AssembleRegion();

  // Environment operators (whole-mesh size), assembled on first use. Collective. An online
  // run with a matching saved model does not assemble them (its reuse checks apply the
  // environment operator matrix-free, see EnvApply).
  std::unique_ptr<mfem::HypreParMatrix> A_env, K_env, D_env, M_sheet, M_sheet_env;
  const mfem::HypreParMatrix &Aenv();
  // Energy operator: A_env for electrostatics, the volume terms for magnetostatics.
  const mfem::HypreParMatrix &Kenv();
  // D_env = A_env - K_env = eps M_env (magnetostatics).
  const mfem::HypreParMatrix &Denv();
  // Sheet masses on both sides and on the environment side.
  const mfem::HypreParMatrix &Msheet();
  const mfem::HypreParMatrix &MsheetEnv();
  // London masses of the London domains on both sides and on the environment side (null
  // without them).
  std::unique_ptr<mfem::HypreParMatrix> M_london, M_london_env;
  const mfem::HypreParMatrix *Mlondon();
  const mfem::HypreParMatrix *MlondonEnv();
  // Whether the environment has London sheets or domains, of kinetic energy (collective).
  bool EnvKinetic();
  // y = M_kin^E x, the kinetic load of x in the environment. Collective.
  void KineticEnvMult(const mfem::Vector &x, mfem::Vector &y);

  // Y_k = A_env X_k for true-DOF vectors X_k (M_sheet^E X_k without the volume terms), by
  // the assembled operator if there is one, else matrix-free (set up on first use).
  // Collective.
  std::vector<mfem::Vector> EnvApply(const std::vector<const mfem::Vector *> &X,
                                     bool volume = true) const;
  mutable std::unique_ptr<Operator> env_op, env_sheet_op;

  // The space of the operators and true-DOF vectors: the parent space, broken across the
  // two-sided superconductor films.
  const FiniteElementSpace &SolutionSpace() const;

  // The parent true DOF of a true DOF of SolutionSpace() (the DOF it is a version of).
  int ParentTrueDof(int t) const { return broken ? version_parent[t] : t; }

  // A true-DOF vector of SolutionSpace() on the parent space: its first version, the field
  // of one side on the two-sided films (for the output fields).
  const mfem::Vector &FirstVersion(const mfem::Vector &x, mfem::Vector &y) const;

  // With two-sided superconductor films, the solution space of the operator, broken across
  // them, and its version maps. Collective.
  void SetUpBrokenSpace();

  // Mass regularization eps (1/mu) M / L^2 making the singular magnetostatic curl-curl
  // definite outside the London domains, L the extent of the mesh (the largest side of its
  // bounding box, in units of Lc: 1 for the default Lc), so that it scales as the curl-curl
  // operator whatever the permeability and Lc; mag_eps is eps / L^2. The energies of the
  // sources use the unregularized operator in a form stationary in eps, with an error
  // O(eps^2), but the part of a surface current that its quadrature leaves incompatible
  // with the null space is amplified by 1/eps. The kinetic energies, and with them the
  // energies in the volume, are linear in eps, so London sheets and domains take a much
  // smaller eps. On the mesh of the region, which may lack the sheets of the environment
  // and has another extent, mag_eps is that of the model.
  static constexpr double kMagRegularization = 1.0e-3;
  static constexpr double kMagRegularizationSheets = 1.0e-7;
  // Largest regularization energy of a surface current relative to its magnetic energy
  // (about eps / lambda_1 for one that closes exactly, lambda_1 the smallest nonzero
  // eigenvalue), which bounds the error of its energy.
  static constexpr double kMaxRegularizationFraction = 1.0e-2;
  double mag_eps = kMagRegularization;

  // Whether the surface current J (on the true DOFs of SolutionSpace()) closes: r = G^T J,
  // G the discrete gradient from the H1 solution space, vanishes at the H1 true DOFs off
  // the conductors (the PEC boundaries and superconductor sheets) and sums to zero over
  // each connected conductor, relative to the largest entry of |G|^T |J| (or, for a
  // conductor, the sum of its entries if larger), up to the quadrature of the excitation.
  // This is the compatibility of J with the null space of the curl-curl and sheet
  // operators, the gradients of the functions constant on each conductor, so it does not
  // depend on the regularization. On the mesh of the region, whose currents may close in
  // the environment, neither the environment elements the mesh holds nor the conductors
  // touching them are checked. Collective.
  bool CurrentCloses(const Vector &J);
  static constexpr double kClosureTol = 1.0e-4;
  // |G|, and the conductor of each H1 true DOF: the least global H1 true DOF of its
  // connected conductor, kOffConductor off the conductors, or kUnchecked. Set up on first
  // use. Collective.
  std::unique_ptr<mfem::HypreParMatrix> closure_G_abs;
  std::vector<HYPRE_BigInt> conductor;
  static constexpr HYPRE_BigInt kOffConductor = -1, kUnchecked = -2;
  const Operator &ClosureGradient() const;
  void SetUpClosure();

  // London superconductor sheets (magnetostatics): whether there are any, the local
  // boundary elements of the one-sided sheets on Gamma, and the sides of these. Like any
  // boundary term of the native operators, a sheet face belongs to the side of its neighbor
  // with the lower domain attribute (see Mesh::GetCeedBdrAttributesByNeighbor), so exactly
  // one side; on Gamma, either side, which the mesh of the region keeps the environment
  // neighbors to tell (see WriteRegionMesh), and which an online run checks against the
  // model, since it assembles only the region's: the moments of AttributeMoments of the
  // faces on Gamma that the region holds, of which the area and the first moments are kept
  // by a refinement of the region. Checks the two-sided films against the interface.
  // Collective.
  bool has_sheets = false, has_london = false;
  std::vector<int> gamma_sheet_faces;
  std::array<double, 6> gamma_sheet_sides = {};
  void CheckSheets();

  // Sheet and London masses of the region, for the flux-loop sources and energies.
  std::unique_ptr<mfem::HypreParMatrix> M_sheet_region, M_london_region;

  // Source modes of magnetostatic sources (see ComputeSourceModes): ids, fingerprints of
  // the environment sources, the interface loads g and energy couplings h (Gamma rows), the
  // constants c and the works B_ij = b^E_i^T W_j, and with S^Kin their kinetic parts
  // h^kin_k = Z^T M_kin^E (W_k - a_k) (Gamma rows) and c^kin_ij = (W_i - a_i)^T
  // M_kin^E (W_j - a_j). Only sources with an environment part have modes.
  std::vector<int> src_ids;
  std::vector<double> src_fp, SG_rows, SH_rows, SC, SB, SHkin_rows, SCkin;
  bool src_ready = false;

  // Environment source b^E = M_sheet^E a + J|_E (J on the environment interior; its
  // interface entries load the region), zero on the Dirichlet DOFs, and its fingerprint,
  // independent of the partition (which changes the basis of some Nédélec DOFs): whether
  // b^E is nonzero (fp[0] = 0: no environment part), and r_k^T b^E for two fixed fields
  // r_k. Collective.
  mfem::Vector EnvSource(const Source &src, std::array<double, 3> &fp);

  // Region source: M_sheet^R a + J on the region-free DOFs.
  mfem::Vector RegionSource(const Source &src) const;

  // Source modes of magnetostatic sources (a_k, J_k). With the environment field
  // u^E = Z u_Gamma + W_k (Z the A-harmonic extension, W_k = A_EE^-1 b^E_k), the
  // environment contributes
  //   g_k = (b^E_k - A_env W_k)|_Gamma                   to the region-condensed load, and
  //   u_i,G^T S^K u_j,G + u_i,G^T h_j + h_i^T u_j,G + c_ij   to the energies, with
  //   h_k = Z^T [(A_env - D_env) W_k - b^E_k],
  //   c_ij = W_i^T K_env W_j + (W_i - a_i)^T M_sheet^E (W_j - a_j) - J_i^E^T W_j
  //          - J_j^E^T W_i
  // (the environment part of the energy of SourceEnergyMatrix). A batched environment
  // solve and a condensation per source (two with S^Kin, for the kinetic parts).
  void ComputeSourceModes(const std::vector<int> &ids, const std::vector<Source> &src);

  // Whether the modes cover every source with an environment part (need[j]; same id and
  // fingerprints), with col[j] its mode column (-1: none, zero modes). On the mesh of the
  // region, a source takes the mode of its id whenever the model has one, as do the
  // model_sources. Collective.
  bool SourceModesMatch(const std::vector<int> &ids, const std::vector<Source> &src,
                        std::vector<int> &col, std::vector<char> &need);

  // Append a record of size bytes to the model file, whose payload rank 0 writes with
  // write (see substructuring::ModelFile::Append). Collective.
  void AppendRecord(const std::string &path, std::uint32_t tag, long long size,
                    const std::function<void(std::ostream &)> &write) const;

  void AppendSourceModes(const std::string &path);

  // Global parent true DOFs of the interface, in global interface order (replicated).
  // Collective.
  std::vector<HYPRE_BigInt> InterfaceTrueDofs() const;

  // Whether the interface holds versions of DOFs split by two-sided films crossing Gamma.
  // Collective.
  bool InterfaceSplit() const;

  // Geometric signatures of the interface DOFs (see TrueDofSignatures), replicated, in
  // global interface order. With width 13 (see InterfaceSplit), the versions of a DOF split
  // by a two-sided film crossing Gamma, which share its signature, have the side moment
  // t·Σ_f (c_f - m) × n_f: t the tangent of the DOF (its e_b moments), m the midpoint of
  // its edge, and f the boundary faces of the region (on Gamma or on the boundary of the
  // mesh) at the edge on the side of the version, of center c_f and unit normal n_f out of
  // the region. It is zero for the other DOFs and flips with the DOF orientation, and every
  // mesh of the region has these faces. Collective.
  std::vector<double> InterfaceSignatures(int width) const;

  // Factor the environment's side (its interior and Gamma) for its Schur complement onto
  // Gamma on first use, which then gives every environment interior solve (an online run
  // with a saved model needs it only for environment fields). Collective.
  void EnsureEnv();

  // Recover the environment interior of full fields (region + Dirichlet parts set),
  // batched: u_E = A_EE^-1 (f - A_env u)|_E, with the environment sources f (if any).
  void RecoverEnvInterior(std::vector<mfem::Vector> &u,
                          const std::vector<mfem::Vector> *f = nullptr);

  // E_ij = u_i^T K u_j with the energy operator of the whole domain, applied once per
  // field. Collective.
  mfem::DenseMatrix EnergyGram(const std::vector<mfem::Vector> &u);

  // Unit Dirichlet lift of terminal idx: 1 on its true DOFs.
  mfem::Vector TerminalMode(int idx) const;

  // Whether the energy operator differs from the solve operator (magnetostatics).
  bool EnergyDiffers() const { return magnetostatic; }

  // Fingerprint of a lift's environment side (t = A_env x), independent of the partition
  // (which changes the basis of some Nédélec DOFs): the energies x^T A_env x and
  // (x + r)^T A_env (x + r) for a fixed field r. Collective.
  std::array<double, 2> LiftFingerprint(const mfem::Vector &x, const mfem::Vector &t) const;

  // Modes (see the members) of the lifts x_k from K batched environment solves. Lifts that
  // do not reach the environment interior skip the solve.
  void ComputeModes(const std::vector<int> &ids, const std::vector<mfem::Vector> &x);

  // Whether the current modes cover every lift that reaches the environment (need[j]; same
  // id and fingerprints, or on the mesh of the region an id of the model), with col[j] its
  // mode column (-1: none, zero modes). On the mesh of the region, an excitation the model
  // lacks must have DOFs in it (one without is of the environment). Collective.
  bool ModesMatch(const std::vector<int> &ids, const std::vector<mfem::Vector> &xd,
                  std::vector<int> &col, std::vector<char> &need) const;

  // Per-rank counts/displacements of the rows of an interface-row-distributed array
  // (gamma_nloc each).
  void RowLayout(std::vector<int> &cnt, std::vector<int> &disp) const;

  // Gather an interface-row-distributed array (width w) to rank 0 (full, global row order).
  std::vector<double> GatherRows(const std::vector<double> &rows, int w) const;

  // Scatter a full (rank 0) interface-row array of width w onto the owned rows.
  void ScatterRows(const std::vector<double> &full, int w, std::vector<double> &rows) const;

  // Append an interface operator (symmetric nG x nG, row-major, each rank its rows like
  // S_rows) to the model as the record tag, by its lower triangle (see
  // TriangleSize): column b is row b from its diagonal, so that the rows of
  // a rank give a contiguous part of it, which rank 0 writes in rank order (see
  // ToRoot). Collective.
  void AppendRows(const std::string &path, std::uint32_t tag,
                  const std::vector<double> &rows) const;

  // The rows of this rank of an interface operator of the model at byte offset pos, mapped
  // onto the current interface (S_cur = M^-T S M^-1, see SignatureMap): each rank reads its
  // share of the lower triangle, whose entries go to the ranks holding the rows they add to
  // (see RedistributeTriangle). Collective.
  std::vector<double> ReadRows(const std::string &path, std::streamoff pos,
                               const SignatureMap &map) const;

  // The records of a model file (see substructuring::ModelFile): HEAD [kModelVersion,
  // Solver.Order, nG, signature type], SIGS and LENC (substructuring::InterfaceRecords: the
  // signatures nG x width and Model.Lc of the offline run), SENV [S_E, lower triangle],
  // ENVF [environment fingerprint, 5 doubles], in magnetostatics MEPS [the regularization
  // eps] and SIDS [the sides of the sheet faces on Gamma, 6 doubles], ENVR [the
  // EnvironmentRecord], ELPS [the flux loops of the environment, which the mesh of the
  // region leaves out], and as they are computed SKEN [S^K, lower triangle]
  // (magnetostatics), SKIN [S^Kin, lower triangle] (with London sheets or domains in the
  // environment), MODE [K][ids][fingerprints 2K][G nG x K][Cmode KxK] (electrostatics),
  // SRCM [K][ids][fingerprints 2K][g nG x K][h nG x K][c KxK][B KxK] (then [h^kin nG x
  // K][c^kin KxK] in a model with S^Kin), and FLXD [K][L][excitation ids K][loop ids L]
  // [fluxoids KxL] (flux loops driven together).
  static constexpr int kModelVersion = 11;
  using ModelFile = substructuring::ModelFile;
  static constexpr std::uint32_t kHeadTag = ModelFile::Tag("HEAD"),
                                 kSETag = ModelFile::Tag("SENV"),
                                 kEnvFpTag = ModelFile::Tag("ENVF"),
                                 kEpsTag = ModelFile::Tag("MEPS"),
                                 kSheetSidesTag = ModelFile::Tag("SIDS"),
                                 kEnvRecordTag = ModelFile::Tag("ENVR"),
                                 kSKenTag = ModelFile::Tag("SKEN"),
                                 kSKinTag = ModelFile::Tag("SKIN"),
                                 kModesTag = ModelFile::Tag("MODE"),
                                 kSourceModesTag = ModelFile::Tag("SRCM"),
                                 kEnvLoopsTag = ModelFile::Tag("ELPS"),
                                 kFluxoidsTag = ModelFile::Tag("FLXD");

  // The flux loops of the environment, left out of the mesh of the region: no DOF of their
  // films and holes in the region or on Gamma, and no surface current or aperture on them.
  // Set by WriteRegionMesh, read from the model online (ReadEnvLoops, collective).
  std::vector<int> env_loops;
  static std::vector<int> ReadEnvLoops(const std::string &path, MPI_Comm comm);

  // Online: whether the excitations of the environment are those of the model (always on
  // the mesh of the region), and the sources taken from the model whole (see
  // SubstructuringSolver::TakeModelExcitations).
  bool env_excitations_unchanged = false;
  std::set<int> model_sources;

  // Append the fluxoids that the flux-loop excitations ids impose on the loops loop_ids
  // (see SubstructuringSolver::EnvironmentFluxoids) to the model file. Collective.
  void AppendFluxoids(const std::string &path, const std::vector<int> &ids,
                      const std::vector<int> &loop_ids,
                      const mfem::DenseMatrix &fluxoid) const;

  // The width of the interface signatures of a model by their type: 1 H1 coordinates, 2
  // H(curl) moments, 3 with the side moments of split DOFs (0: not a model).
  static int SignatureFileWidth(int sig_type);

  // Fixed projected polynomial field k of the fingerprints (the same on both sides of the
  // two-sided films), about its center: fixed fractions of the bounding box of the centers
  // of the environment elements, which does not depend on the region. Collective.
  mfem::Vector FingerprintField(int k) const;
  const std::array<std::array<double, 3>, 3> &FingerprintCenters() const;
  mutable std::optional<std::array<std::array<double, 3>, 3>> fp_centers;

  // Fingerprint of the environment operator: global counts of the environment-closure DOFs
  // and of the Dirichlet DOFs of environment elements, and the energies r^T A_env r of the
  // three fixed fields. Independent of the region (A_env has no region contribution), of
  // the partition, the DOF numbering and the H(curl) orientations, but it changes with the
  // environment mesh, materials, order, physics or Dirichlet boundaries. With A_env applied
  // by EnvApply (an online run assembles no A_env). Collective.
  static constexpr int kEnvFpSize = 5;
  std::array<double, kEnvFpSize> EnvironmentFingerprint() const;
  mutable std::optional<std::array<double, kEnvFpSize>> env_fp;

  // The fixed fields 0 and 1 of the fingerprints, cached.
  const mfem::Vector &CachedFingerprintField(int k) const;
  mutable std::array<mfem::Vector, 2> fp_fields;

  // Environment fingerprint of a loaded model.
  std::array<double, kEnvFpSize> saved_env_fp = {};
  bool have_env_fp = false;

  // Append the EnvironmentRecord, with written marking the elements written in the mesh of
  // the region. Collective.
  void AppendEnvRecord(const std::string &path, const std::vector<char> &written);

  void AppendEnvFingerprint(const std::string &path);

  // The true DOFs of the environment elements touching the marked true DOFs, those written
  // in the mesh of the region (see palace::WriteRegionMesh). Collective.
  std::vector<char> WrittenEnvTrueDofs(const std::vector<char> &keep) const;

  // The two-sided films of the mesh of the region (see WriteRegionMesh), whose interface
  // DOFs must have the versions of the full mesh: the films of the region, and the films
  // given whole (their DOFs are kept). If a written environment element touches an
  // interface DOF on a film crossing Gamma, the DOFs of the environment film faces at such
  // DOFs are kept, which writes the environment along every crossing. An interface DOF at
  // the edge of a film of the region is kept, for the environment to join its sides.
  // Returns the other film faces of the environment, written as plain boundaries (the
  // environment DOFs are not solved for). Collective.
  std::vector<char> RegionFilmClosure(std::vector<char> &keep,
                                      const std::set<int> &whole) const;

  // The mesh of the region for online runs: the region, with the environment elements
  // touching the region's Dirichlet DOFs (set on Γ by boundary conditions of the
  // environment) and, in magnetostatics, those of the surface currents and of the flux
  // loops not in env_loops (their films whole), whose excitations need their boundaries,
  // and the two-sided films of RegionFilmClosure. Collective.
  std::vector<char> WriteRegionMesh(const std::string &path);

  // Abort if a loaded model was condensed from a different environment than the current
  // one.
  void CheckEnvFingerprint();

  void AppendModes(const std::string &path);

  // Load the records after S_E (collective). Interface rows are mapped onto the online
  // interface like S_E (dual quantities: x_cur = M^-T x_saved, rows and S^K columns); modes
  // absent here are recomputed on demand.
  void LoadSections(const std::string &path, const SignatureMap &map);

  // Selection matrix E (parent true DOFs x marked DOFs): E(i, k) = 1 for the k-th marked
  // DOF, numbered globally rank by rank in local order (for the interface, as
  // gamma_global). Each rank's marked DOFs are its own contiguous block, so the matrix is
  // local.
  std::unique_ptr<mfem::HypreParMatrix>
  AssembleSelection(const std::vector<char> &mark) const;

  // Factor the region on its first solve: its side (its interior and Gamma, the region-free
  // DOFs) for its Schur complement onto Gamma, with S_E (see RegionSolver). Collective.
  void EnsureRegionSolver();

  // Region solves on the region-free DOFs (rf_idx), in place (right-hand sides in,
  // solutions out), in one batch. Collective.
  void SolveRegion(std::vector<mfem::Vector> &x);

  // Full parent vectors from Dirichlet values xd (on dbc_tdofs) and from region-free values
  // x_rf (and optional Dirichlet values), zero elsewhere.
  mfem::Vector DirichletField(const mfem::Vector &xd) const;
  mfem::Vector RegionField(const mfem::Vector &x_rf,
                           const mfem::Vector *xd = nullptr) const;

  // Batched environment interior solves y_k = A_EE^-1 x_k (inputs/outputs are parent
  // true-DOF vectors; only the environment interior of y_k is kept, the rest is zeroed), by
  // the environment solver. y_k may alias x_k.
  void ApplyAeeInvMulti(const std::vector<const mfem::Vector *> &X,
                        const std::vector<mfem::Vector *> &Y) const;

  // The condensations P^T v_k = v_k|_Gamma - A_Gamma,E A_EE^-1 v_k|_E on the owned
  // interface rows (row-major, width |V|), by the environment solver. Collective.
  void CondenseRows(const std::vector<const mfem::Vector *> &V,
                    std::vector<double> &rows) const;

  // Zero the entries outside the environment interior.
  void MaskToEnvInterior(const std::vector<mfem::Vector *> &Y) const;

  // Number the interface true DOFs globally, in rank order, and size the rows of S_E.
  // Collective.
  void EnumerateInterface();

  // Load S_E and the sections of the saved model onto this interface (the head was read
  // and checked by the constructor). Collective.
  void LoadModel(const std::string &model_path, int sig_w);

  // Materialize S_E (from the partial factorization, or S_E e_c = (A_env e_c)|_Gamma -
  // (A_env A_EE^-1 (A_env e_c)|_E)|_Gamma), S^K = S_E - P^T D P and, with London sheets or
  // domains in the environment, S^Kin = P^T M_kin^E P, for a magnetostatic model to save:
  // by batched multi-RHS environment solves for the lifts P e_c, and the condensations P^T
  // D P e_c and P^T M_kin^E P e_c. Collective.
  void MaterializeWithEnergy();

  // Save S_E, the region mesh and the sections to the model file. Collective.
  void SaveModel(const std::string &model_path, int sig_w);
};

}  // namespace palace

#endif  // PALACE_MODELS_SUBSTRUCTURING_SOLVER_IMPL_HPP
