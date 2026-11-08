// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SUBSTRUCTURING_SOLVER_IMPL_HPP
#define PALACE_MODELS_SUBSTRUCTURING_SOLVER_IMPL_HPP

#include "models/substructuringsolver.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fem/substructure.hpp"
#include "linalg/gmg.hpp"
#include "linalg/hodlr.hpp"
#include "linalg/ksp.hpp"
#include "linalg/mumpsschur.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "linalg/vector.hpp"
#include "models/superconductorsheetoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace substructuring
{

// Sparse direct solver for the environment and region factorizations: MUMPS, SuperLU_DIST
// or STRUMPACK, in this order of preference (null without one).
constexpr const char *kDirectSolverName =
#if defined(MFEM_USE_MUMPS)
    "MUMPS";
#elif defined(MFEM_USE_SUPERLU)
    "SuperLU_DIST";
#elif defined(MFEM_USE_STRUMPACK)
    "STRUMPACK";
#else
    "";
#endif
constexpr bool kHasDirectSolver = (kDirectSolverName[0] != '\0');

#if defined(MFEM_USE_MUMPS)
// MUMPS factorization of a symmetric operator (MumpsSchurSolver without Schur variables:
// silent, and retried with a larger workspace when the estimate is too small), optionally
// on rank 0 alone.
class MumpsDirectSolver : public mfem::Solver
{
public:
  MumpsDirectSolver(bool serial, SymbolicFactorization ordering)
    : serial(serial), ordering(ordering)
  {
  }
  void SetOperator(const mfem::Operator &op) override
  {
    const auto *A = dynamic_cast<const mfem::HypreParMatrix *>(&op);
    MFEM_VERIFY(A, "MumpsDirectSolver requires a HypreParMatrix operator!");
    height = width = A->Height();
    mumps = std::make_unique<MumpsSchurSolver>(*A, std::vector<HYPRE_BigInt>{}, 0.0, serial,
                                               ordering);
  }
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    mumps->SolveInternal({&x}, {&y});
  }
  void ArrayMult(const mfem::Array<const mfem::Vector *> &X,
                 mfem::Array<mfem::Vector *> &Y) const override
  {
    mumps->SolveInternal(std::vector<const mfem::Vector *>(X.begin(), X.end()),
                         std::vector<mfem::Vector *>(Y.begin(), Y.end()));
  }

private:
  bool serial;
  SymbolicFactorization ordering;
  std::unique_ptr<MumpsSchurSolver> mumps;
};
#endif

// With serial, a small system is factored on rank 0 alone (MUMPS only).
inline std::unique_ptr<mfem::Solver> MakeDirectSolver([[maybe_unused]] const IoData &iodata,
                                                      [[maybe_unused]] MPI_Comm comm,
                                                      [[maybe_unused]] bool serial = false)
{
#if defined(MFEM_USE_MUMPS)
  return std::make_unique<MumpsDirectSolver>(serial,
                                             iodata.solver.linear.sym_factorization);
#elif defined(MFEM_USE_SUPERLU)
  return std::make_unique<SuperLUSolver>(iodata, comm, 0);
#elif defined(MFEM_USE_STRUMPACK)
  return std::make_unique<StrumpackSolver>(iodata, comm, 0);
#else
  return nullptr;
#endif
}

// Multi-RHS solves with a direct factor whose RHS count is fixed by its first solve (MFEM's
// SuperLU wrapper does not allow changing it): the block size is set by the first call (or
// up front), and batches are split and zero-padded to it. Inputs and outputs may alias.
class BlockedDirectSolver
{
public:
  BlockedDirectSolver(std::unique_ptr<mfem::Solver> &&solver, int block = 0)
    : solver(std::move(solver)), block(block)
  {
  }

  void SetOperator(const mfem::Operator &op) { solver->SetOperator(op); }

  void Mult(const std::vector<const mfem::Vector *> &X,
            const std::vector<mfem::Vector *> &Y) const
  {
    const int n = static_cast<int>(X.size());
    if (n == 0)
    {
      return;
    }
    const int size = X[0]->Size();
    if (block == 0)
    {
      block = std::min(kMaxBlock, n);
    }
    pad_in.SetSize(size);
    pad_in = 0.0;
    pad_out.SetSize(size);
    mfem::Array<const mfem::Vector *> Xp(block);
    mfem::Array<mfem::Vector *> Yp(block);
    for (int c0 = 0; c0 < n; c0 += block)
    {
      const int nb = std::min(block, n - c0);
      for (int k = 0; k < block; k++)
      {
        Xp[k] = (k < nb) ? X[c0 + k] : &pad_in;
        Yp[k] = (k < nb) ? Y[c0 + k] : &pad_out;
        Yp[k]->SetSize(size);
      }
      if (block == 1)
      {
        pad_out = *Xp[0];  // the solve may not be done in place
        solver->Mult(pad_out, *Yp[0]);
      }
      else
      {
        solver->ArrayMult(Xp, Yp);
      }
    }
  }

  static constexpr int kMaxBlock = 32;

private:
  std::unique_ptr<mfem::Solver> solver;
  mutable int block;
  mutable mfem::Vector pad_in, pad_out;
};

// Region-condensed system operator: (region operator on region-free true DOFs) + implicit
// DtN on the interface, with non-region-free true DOFs pinned to identity.
class RegionCondensedOperator : public mfem::Operator
{
public:
  RegionCondensedOperator(mfem::HypreParMatrix &A_region_free, const mfem::Operator &dtn,
                          const std::vector<char> &is_region_free)
    : mfem::Operator(A_region_free.Height()), A_region_free(A_region_free), dtn(dtn),
      is_region_free(is_region_free), td(A_region_free.Height())
  {
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    A_region_free.Mult(x, y);
    dtn.Mult(x, td);
    const double *t = td.HostRead();
    double *yy = y.HostReadWrite();
    for (int i = 0; i < height; i++)
    {
      if (is_region_free[i])
      {
        yy[i] += t[i];
      }
    }
  }

private:
  mfem::HypreParMatrix &A_region_free;
  const mfem::Operator &dtn;
  const std::vector<char> &is_region_free;
  mutable mfem::Vector td;
};

// Materialized DtN: applies the precomputed interface operator S_E to a distributed
// interface vector (gathered by an Allreduce): each rank multiplies its own dense row
// block, or applies the replicated HODLR form and reads back its rows.
class MaterializedDtN : public mfem::Operator
{
public:
  MaterializedDtN(const std::vector<double> &S_rows, int row_off,
                  const std::vector<int> &gamma_global, int nG_global, MPI_Comm comm,
                  const Hodlr *hodlr = nullptr)
    : mfem::Operator(static_cast<int>(gamma_global.size())), S_rows(S_rows),
      row_off(row_off), gamma_global(gamma_global), nG_global(nG_global), comm(comm),
      hodlr(hodlr), xg(nG_global)
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
    if (hodlr)
    {
      xp.resize(nG_global);
      yp.resize(nG_global);
      yg.resize(nG_global);
    }
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
    if (hodlr)
    {
      for (int k = 0; k < nG_global; k++)
      {
        xp[k] = xg[hodlr->perm[k]];
      }
      hodlr->Mult(xp.data(), yp.data());
      for (int k = 0; k < nG_global; k++)
      {
        yg[hodlr->perm[k]] = yp[k];
      }
      for (int i = 0; i < height; i++)
      {
        if (gamma_global[i] >= 0)
        {
          yy[i] = yg[gamma_global[i]];
        }
      }
      return;
    }
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
  const Hodlr *hodlr;
  std::vector<int> cnt, disp;
  mutable std::vector<double> xo, xg, xp, yp, yg;
};

}  // namespace substructuring

// Region/environment operators are assembled on the parent space with domain-restricted
// coefficients; the interface is identified in true-DOF space and the environment condensed
// onto it (a distributed dense DtN).
struct SubstructuringSolver::Impl
{
  const IoData &iodata;
  Mesh &mesh;
  mfem::ParMesh &parent;
  bool magnetostatic;
  std::unique_ptr<mfem::FiniteElementCollection> fec;
  FiniteElementSpace parent_space;
  mfem::ParFiniteElementSpace &parent_fes;
  // Magnetostatics with two-sided superconductor films: the parent space broken across the
  // films, on which the operators and true-DOF vectors act (see SolutionSpace()), its maps
  // to the parent true DOFs, and the parent true DOF of each of its true DOFs.
  std::unique_ptr<FiniteElementSpace> broken;
  std::unique_ptr<mfem::HypreParMatrix> first_version, all_versions;
  std::vector<int> version_parent;
  int nt;

  std::unique_ptr<mfem::HypreParMatrix> A_region, A_region_free;
  mfem::Array<int> ra_arr, ea_arr;
  // Iterative environment interior solve on the environment submesh (Gamma essential),
  // with geometric multigrid for higher-order H1.
  std::unique_ptr<mfem::ParSubMesh> env_submesh_owned;    // single-level ownership
  std::vector<std::unique_ptr<Mesh>> env_mesh_vec;        // GMG: [0] owns the ParSubMesh
  mfem::ParSubMesh *env_submesh = nullptr;                // raw ptr to the owned submesh
  std::unique_ptr<mfem::ParFiniteElementSpace> env_sfes;  // single-level solve space
  mfem::ParFiniteElementSpace *env_solve_fes = nullptr;   // fespace the solver acts on
  std::vector<std::unique_ptr<mfem::H1_FECollection>> env_fecs;
  std::unique_ptr<FiniteElementSpaceHierarchy> env_hierarchy;
  std::unique_ptr<MultigridOperator> env_mg_op;
  std::vector<mfem::Array<int>> env_dbc_lists;     // per-level essential DOFs
  std::unique_ptr<mfem::HypreParMatrix> env_A_ee;  // single-level submesh operator
  std::unique_ptr<KspSolver> env_ksp;
  // Direct factorization of the region-free block (the region pc, or the exact region solve
  // when region_exact) and its (local) region-free true DOFs.
  std::unique_ptr<substructuring::BlockedDirectSolver> reg_lu;
  bool region_exact = false;
  // Local true DOFs: region-free (the region solve), of region elements (where the region
  // operators act), and interface.
  mfem::Array<int> rf_idx, rc_idx, gamma_idx;
  mutable mfem::Vector reg_r, reg_z;
  std::unique_ptr<mfem::HypreParMatrix> env_A_ee_parent;  // parent-space A_EE (eliminated)
  std::unique_ptr<substructuring::BlockedDirectSolver> env_lu_parent;  // its direct factor
#if defined(MFEM_USE_MUMPS)
  // MUMPS partial factorization of the environment with Gamma as the Schur variables: gives
  // S_E directly and serves every environment interior solve (set up by MaterializeMumps).
  std::unique_ptr<MumpsSchurSolver> env_mumps;
#endif
  // Block size of the batched multi-RHS environment solves.
  static constexpr int kMaterializeBlock = substructuring::BlockedDirectSolver::kMaxBlock;
  // Largest interface for which the region preconditioner factors the dense Gamma block of
  // the exact condensed operator (dense block + its fill ~ 3 |Gamma|^2 doubles in total).
  static constexpr int kDirectCondensedMaxInterface = 10000;
  // The exact dense-Gamma region factor (~|Gamma|^3 once) pays off against CG on
  // A_region_free (~|Gamma|^2 per iteration, per excitation) from about |Gamma| /
  // kDenseInterfacePerExcitation excitations in a batch.
  static constexpr int kDenseInterfacePerExcitation = 1000;
  // Region-free size up to which the region is factored on one rank (MUMPS's per-rank
  // workspace dominates a small factorization spread over many ranks).
  static constexpr long long kSerialRegionFactorDofs = 100000;
  // Conjugate gradients without a direct factorization: the environment interior solves
  // condense S_E column by column, so their error enters S_E itself, while the energies
  // are stationary in the error of the region solves.
  static constexpr double kEnvCgTol = 1.0e-12, kRegionCgTol = 1.0e-10;
  static constexpr int kEnvCgMaxIts = 1000, kRegionCgMaxIts = 2000;
  // Relative tolerance of the fingerprint checks of a saved model: the fingerprints are
  // sums over elements, which the partition changes by rounding only.
  static constexpr double kFingerprintTol = 1.0e-9;
  static constexpr const char *kEnvironmentDiffers =
      "The environment differs from the one the saved substructuring model was condensed "
      "from (environment mesh, materials, boundary conditions, order or problem type "
      "changed): rerun in \"Offline\" mode to condense it again!";
  bool region_direct = false;          // region pc is a direct factor (created lazily)
  bool region_dense_eligible = false;  // ... and a dense S_E is available for the exact one
  mfem::Array<int> non_env_int;        // parent true DOFs that are not environment-interior
  bool env_parent_direct = false;
  mfem::Array<int> env_ess;  // solve-space essential true DOFs (Gamma + env Dirichlet)
  mutable mfem::ParGridFunction env_pgf, env_sgf;  // parent / submesh transfer buffers
  mutable mfem::Vector env_srhs, env_ssol;
  // Region preconditioner by geometric multigrid on the region submesh (higher-order H1;
  // region Dirichlet terminals essential, Gamma free).
  std::vector<std::unique_ptr<Mesh>> reg_mesh_vec;
  mfem::ParSubMesh *reg_submesh = nullptr;
  mfem::ParFiniteElementSpace *reg_solve_fes = nullptr;
  std::vector<std::unique_ptr<mfem::H1_FECollection>> reg_fecs;
  std::unique_ptr<FiniteElementSpaceHierarchy> reg_hierarchy;
  std::unique_ptr<MultigridOperator> reg_mg_op;
  std::vector<mfem::Array<int>> reg_dbc_lists;
  std::unique_ptr<GeometricMultigridSolver<Operator>> reg_gmg;  // one cycle per apply
  mfem::Array<int> reg_ess;
  mutable mfem::ParGridFunction reg_pgf, reg_sgf;
  mutable mfem::Vector reg_srhs, reg_ssol;
  // Region energy operator: the solve operator (electrostatics), the pure curl-curl
  // (magnetostatics).
  std::unique_ptr<mfem::HypreParMatrix> A_region_energy;
  mfem::HypreParMatrix *K_region_e = nullptr;
  std::vector<char> is_gamma, is_env_int, is_region_free;
  mfem::Array<int> dbc_tdofs;
  int env_dbc_local = 0;                           // owned Dirichlet DOFs of env elements
  std::map<int, std::vector<int>> terminal_tdofs;  // terminal index -> its true DOFs
  std::vector<char> region_dbc;                    // Dirichlet DOFs of region elements
  // Online on a mesh of the region (written offline next to the model, whose environment
  // elements it holds): Γ also holds the region's boundary DOFs at the saved interface
  // DOFs, and the environment is that of the model.
  bool region_only = false;

  std::unique_ptr<substructuring::RegionCondensedOperator> region_op;
  std::unique_ptr<KspSolver> region_ksp;  // CG on the region-condensed operator

  // Interface operator S_E over a global interface enumeration, distributed by rows: each
  // rank stores its own rows [gamma_off, gamma_off + gamma_nloc) (row-major, x nG_global).
  std::vector<int> gamma_global;  // owned parent true DOF -> global interface index, or -1
  int nG_global = 0;
  std::vector<double> S_rows;
  int gamma_off = 0, gamma_nloc = 0;
  // Optional compressed (HODLR) form of S_E, replicated, used in place of S_rows.
  std::unique_ptr<Hodlr> hodlr;
  std::unique_ptr<substructuring::MaterializedDtN> mat_dtn;

  // Dirichlet-lift modes: the environment condensed onto Gamma + a set of Dirichlet lifts
  // x_k (e.g. the terminal unit potentials). With t_k = A_env x_k and
  // w_k = A_EE^-1 t_k|_E: the region-solve interface load g_k = (t_k - A_env w_k)|_Gamma
  // (G_rows, owned interface rows: G_rows[r * K + k]). For the energy (operator K_env;
  // equal to A_env for electrostatics, the pure curl-curl for magnetostatics) the
  // environment's exact condensed energy is [u_G; V]^T [S^K, G^K; G^K^T, Cmode] [u_G; V]
  // with S^K = P^T K P = S_E - P^T D P (D = A_env - K_env, P the A-harmonic lifting), G^K
  // (GK_rows; empty when K_env = A_env) and the replicated Cmode. A partition-invariant
  // fingerprint per mode (LiftFingerprint) guards reuse.
  std::vector<int> mode_ids;
  std::vector<double> mode_fp;  // 2 per mode
  std::vector<double> G_rows, GK_rows, Cmode;
  bool modes_ready = false;
  bool modes_saved = false;  // appended to the model file in this run
  // Magnetostatic energy interface operator S^K (distributed like S_rows), materialized
  // with a saved model; otherwise the energies use environment solves (electrostatics: S^K
  // = S_E).
  std::vector<double> SK_rows;
  // Rank-independent presence flags (the distributed arrays are empty on ranks that own no
  // interface DOFs, so their emptiness must not gate collectives).
  bool have_sk = false, have_gk = false;
  std::unique_ptr<Hodlr> hodlr_K;
  std::unique_ptr<substructuring::MaterializedDtN> mat_dtn_K;
  bool env_built = false;  // environment interior solver set up

  Impl(const IoData &iodata, Mesh &mesh);

  // Environment operators (whole-mesh size), assembled on first use. Collective. An online
  // run with a matching saved model does not assemble them (its reuse checks run element by
  // element, see EnvApply).
  std::unique_ptr<mfem::HypreParMatrix> A_env, K_env, D_env, M_sheet, M_sheet_env;
  const mfem::HypreParMatrix &Aenv();
  // Energy operator: A_env for electrostatics, the pure curl-curl for magnetostatics.
  const mfem::HypreParMatrix &Kenv();
  // D_env = A_env - K_env = eps M_env (magnetostatics).
  const mfem::HypreParMatrix &Denv();
  // Sheet masses on both sides and on the environment side.
  const mfem::HypreParMatrix &Msheet();
  const mfem::HypreParMatrix &MsheetEnv();

  // The terms of AssembleParent element by element, without assembling them: f(vdofs, A_e)
  // for each local element with an attribute of coef_by_attr and each sheet face of sheet
  // (if non-null) for which use(vdofs) is true, with the integrators and quadrature of
  // AssembleParent, MFEM's DOF transformation and the DOFs of SolutionSpace().
  void ForEachElement(const std::map<int, mfem::DenseMatrix> &coef_by_attr, bool with_mass,
                      bool with_curl, const std::vector<double> *sheet,
                      const std::function<bool(const mfem::Array<int> &)> &use,
                      const std::function<void(const mfem::Array<int> &,
                                               const mfem::DenseMatrix &)> &f) const;

  // Y_k = A_env X_k for true-DOF vectors X_k, element by element over the elements the X_k
  // touch (volume and/or sheet terms). Collective.
  std::vector<mfem::Vector> EnvApply(const std::vector<const mfem::Vector *> &X,
                                     bool volume = true, bool sheets = true) const;

  // The space of the operators and true-DOF vectors: the parent space, broken across the
  // two-sided superconductor films.
  const FiniteElementSpace &SolutionSpace() const;

  // The parent true DOF of a true DOF of SolutionSpace() (the DOF it is a version of).
  int ParentTrueDof(int t) const { return broken ? version_parent[t] : t; }

  // A true-DOF vector of SolutionSpace() on the parent space: its first version, the field
  // of one side on the two-sided films (for the output fields).
  const mfem::Vector &FirstVersion(const mfem::Vector &x, mfem::Vector &y) const;

  // Break the parent space across the two-sided superconductor films of the mesh, and
  // return the two-port terms of the film faces bounding the region and of those bounding
  // the environment (null without films), assembled as in the native operator. Collective.
  std::array<std::unique_ptr<mfem::HypreParMatrix>, 2> SetUpBrokenSpace();

  // Mass regularization eps M making the singular magnetostatic curl-curl definite; the
  // energies use the unregularized operator and are stationary, so their error is O(eps^2).
  // With London sheets, the flux-loop and surface-current sources are orthogonal to the
  // null space of the curl-curl + sheet operator, so a much smaller eps is safe. On the
  // mesh of the region, which may lack the sheets of the environment, eps is that of the
  // model.
  static constexpr double kMagRegularization = 1.0e-3;
  static constexpr double kMagRegularizationSheets = 1.0e-7;
  // Largest regularization energy of a surface current relative to its magnetic energy
  // (about eps / lambda_1 for one that closes, lambda_1 the smallest nonzero eigenvalue).
  static constexpr double kMaxRegularizationFraction = 1.0e-2;
  double mag_eps = kMagRegularization;

  // London superconductor sheets (magnetostatics) of a side: per local boundary element,
  // 1/L_ksq (with the factor for a cracked sheet, whose two sides each carry half) on the
  // side it belongs to, and the two-port terms of the faces of the two-sided films of the
  // side (null without). A sheet face belongs to the region if it bounds a region element
  // (all its DOFs are then region or interface DOFs), else to the environment: independent
  // of the partition and of the region design, and exact either way for a face on Gamma.
  struct Sheets
  {
    std::vector<double> coef;
    std::unique_ptr<mfem::HypreParMatrix> two_sided;
  };
  Sheets sheet_region, sheet_env, sheet_all;
  bool has_sheets = false;
  void SetUpSheets(const mfem::Array<int> &ra,
                   std::array<std::unique_ptr<mfem::HypreParMatrix>, 2> two_sided);

  // Sheet masses (both sides, and per side), for the flux-loop sources and energies.
  std::unique_ptr<mfem::HypreParMatrix> M_sheet_region;

  // Source modes of magnetostatic sources (see ComputeSourceModes): ids, fingerprints of
  // the environment sources, the interface loads g and energy couplings h (Gamma rows), the
  // constants c and the works B_ij = b^E_i^T W_j. Only sources with an environment part
  // have modes.
  static constexpr int kSourceModesMagic = 0x33435253;  // "SRC3"
  std::vector<int> src_ids;
  std::vector<double> src_fp, SG_rows, SH_rows, SC, SB;
  bool src_ready = false, src_saved = false;

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
  // (the environment part of the energy of SourceEnergyMatrix). Two batched environment
  // solves per source.
  void ComputeSourceModes(const std::vector<int> &ids, const std::vector<Source> &src);

  // Whether the modes cover every source with an environment part (need[j]; same id and
  // fingerprints), with col[j] its mode column (-1: none, zero modes). Collective.
  bool SourceModesMatch(const std::vector<int> &ids, const std::vector<Source> &src,
                        std::vector<int> &col, std::vector<char> &need);

  // Write to the model file from rank 0, anew or appended, and abort on every rank if it
  // cannot be written. Collective.
  void WriteModel(const std::string &path, bool append,
                  const std::function<void(std::ofstream &)> &write) const;

  void AppendSourceModes(const std::string &path);

  std::unique_ptr<mfem::HypreParMatrix>
  AssembleParent(const std::map<int, mfem::DenseMatrix> &coef_by_attr, bool with_mass,
                 bool with_curl = true, const Sheets *sheet = nullptr);

  std::map<int, mfem::DenseMatrix> region_eps, env_eps;

  // Global number of environment closure DOFs that a direct factorization would carry.
  long long EnvDirectSize() const;

#if defined(MFEM_USE_MUMPS)
  // Materialize S_E by one MUMPS partial factorization (Schur complement on Gamma) of A_env
  // with the DOFs outside E + Gamma eliminated; the factor is kept as the environment
  // solver. Fills S_rows. Collective.
  void MaterializeMumps();
#endif

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

  // Set up the environment interior solver on first use (an online run with a saved model
  // needs it only for environment fields).
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

  // Modes (see the members) of the lifts x_k from batched environment solves: K solves,
  // plus K more for the energy correction when it differs. Lifts that do not reach the
  // environment interior skip the solve.
  void ComputeModes(const std::vector<int> &ids, const std::vector<mfem::Vector> &x);

  // Whether the current modes cover the given lifts (same ids + environment fingerprints);
  // col[j] receives the mode column of request j. Collective.
  bool ModesMatch(const std::vector<int> &ids, const std::vector<mfem::Vector> &xd,
                  std::vector<int> &col) const;

  // Per-rank counts/displacements of an interface-row-distributed array (gamma_nloc * w
  // each).
  void RowLayout(int w, std::vector<int> &cnt, std::vector<int> &disp) const;

  // Gather an interface-row-distributed array (width w) to rank 0 (full, global row order).
  std::vector<double> GatherRows(const std::vector<double> &rows, int w) const;

  // Scatter a full (rank 0) interface-row array of width w onto the owned rows.
  void ScatterRows(const std::vector<double> &full, int w, std::vector<double> &rows) const;

  // A model file: [kModelMagic][kModelVersion][Solver.Order][nG][signature type]
  // [signatures nG x width][S_E nG x nG], then the sections [kLengthMagic][Model.Lc of the
  // offline run], [kEnvFpMagic][environment fingerprint, 5 doubles], [kEpsMagic][the
  // magnetostatic regularization eps] (magnetostatics), [kEnvMeshMagic][ElementMoments of
  // the environment and of its elements in the mesh of the region, 12 doubles],
  // [kSKenMagic][S^K nG x nG] (magnetostatics), [kModesMagic][K][ids][fingerprints
  // 2K][has_GK][G nG x K][G^K nG x K if has_GK][Cmode KxK], and
  // [kSourceModesMagic][K][ids][fingerprints 2K][g nG x K][h nG x K][c KxK][B KxK].
  static constexpr int kModelMagic = 0x4c444d53;  // "SMDL"
  static constexpr int kModelVersion = 1;
  static constexpr int kHeadInts = 5;
  static constexpr int kSKenMagic = 0x4e454b53;     // "SKEN"
  static constexpr int kModesMagic = 0x33444f4d;    // "MOD3"
  static constexpr int kEnvFpMagic = 0x33564e45;    // "ENV3"
  static constexpr int kLengthMagic = 0x434e454c;   // "LENC"
  static constexpr int kEpsMagic = 0x5350454d;      // "MEPS"
  static constexpr int kEnvMeshMagic = 0x48534d45;  // "EMSH"

  // The width of the interface signatures of a model by their type: 1 H1 coordinates, 2
  // H(curl) moments, 3 with the side moments of split DOFs (0: not a model).
  static int SignatureFileWidth(int sig_type);

  // The head of a model, from rank 0: its format version (-1: the file cannot be read, 0:
  // not a model), order, interface size, signature type and width (0: another version, or
  // an incomplete model), the interface signatures, the offsets of S_E and of the sections
  // after it, the reference length, environment fingerprint, element moments of the
  // environment and of the mesh of the region, and magnetostatic regularization (0 or
  // empty: not saved). Collective.
  struct SavedHead
  {
    int version = -1, order = 0, nG = 0, sig_type = 0, width = 0;
    std::vector<double> sig, env_fp;
    std::array<double, 6> env_moments = {}, region_moments = {};
    double length = 0.0, regularization = 0.0;
    std::streamoff S_pos = 0, sections_pos = 0;
  };
  static SavedHead ReadSavedHead(const std::string &path, MPI_Comm comm);
  // The head of the model of an online run (read by the constructor).
  SavedHead saved;

  // Fixed projected polynomial field k of the fingerprints (the same on both sides of the
  // two-sided films).
  mfem::Vector FingerprintField(int k) const;

  // Fingerprint of the environment: global counts of the environment-closure DOFs and of
  // the Dirichlet DOFs of environment elements, and the energies r^T A_env r of three fixed
  // fields (projected polynomials). Independent of the region (A_env has no region
  // contribution), of the partition, the DOF numbering and the H(curl) orientations, but it
  // changes with the environment mesh, materials, order, physics or Dirichlet boundaries.
  // Computed element by element (no assembled A_env). Collective.
  std::array<double, 5> EnvironmentFingerprint() const;
  mutable std::optional<std::array<double, 5>> env_fp;

  // The fixed fields 0 and 1 of the fingerprints, cached.
  const mfem::Vector &CachedFingerprintField(int k) const;
  mutable std::array<mfem::Vector, 2> fp_fields;

  // Environment fingerprint of a loaded model.
  std::array<double, 5> saved_env_fp = {0.0, 0.0, 0.0, 0.0, 0.0};
  bool have_env_fp = false;

  void AppendEnvMoments(const std::string &path, const std::vector<char> &written);

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
  // environment) and, in magnetostatics, those of the surface currents and flux loops,
  // whose excitations need their boundaries (the films of the flux loops whole), and the
  // two-sided films of RegionFilmClosure. Collective.
  std::vector<char> WriteRegionMesh(const std::string &path);

  // Abort if a loaded model was condensed from a different environment than the current
  // one.
  void CheckEnvFingerprint();

  void AppendSK(const std::string &path) const;
  void AppendModes(const std::string &path);

  // Load the sections at byte offset `pos` (collective). Interface rows are mapped onto the
  // online interface like S_E (dual quantities: x_cur = M^-T x_saved, rows and S^K
  // columns). Reading stops at the first unknown section; a later mode section supersedes
  // an earlier one, and modes absent here are recomputed on demand.
  void LoadSections(const std::string &path, std::streamoff pos, const SignatureMap &map);

  void BuildEnvSubmeshSolver();

  // Essential submesh DOFs (Gamma and the environment Dirichlet DOFs), from a transfer of
  // the parent marker onto the solve space.
  void ComputeEnvEss();

  // Boundary attributes of a submesh whose true DOFs all carry the parent true-DOF marker
  // (1 or 0), from a transfer of the marker. Decided globally, so identical on every rank
  // (the GMG and single-level paths use different collectives).
  mfem::Array<int> EssentialSubmeshAttributes(mfem::ParSubMesh &submesh,
                                              const mfem::Vector &marker);

  // Higher-order H1 environment Dirichlet solve via geometric p-multigrid on the submesh.
  // Returns false (falling back to the single-level solve) if a hierarchy cannot be built.
  bool BuildEnvGmg();

  std::unique_ptr<mfem::HypreParMatrix> AssembleEnvSubmesh();

  // Selection matrix E (parent true DOFs x marked DOFs): E(i, k) = 1 for the k-th marked
  // DOF, numbered globally rank by rank in local order (for the interface, as
  // gamma_global). Each rank's marked DOFs are its own contiguous block, so the matrix is
  // local.
  std::unique_ptr<mfem::HypreParMatrix>
  AssembleSelection(const std::vector<char> &mark) const;

  // Create the direct region preconditioner on the first solve batch (n excitations): the
  // exact dense-Gamma factor A_region_free + S_E when the batch amortizes it, else
  // A_region_free, restricted to the region-free DOFs. Created once (destroying a
  // never-factored SuperLU object crashes in SuperLU_DIST). Collective.
  void EnsureRegionFactor(int n);

  // Region preconditioner apply with the direct factor: identity off the region-free DOFs
  // (A_region_free is identity there), the factor on them.
  void ApplyRegionFactor(const mfem::Vector &r, mfem::Vector &z) const;

  // Region-condensed solves on the region-free DOFs (rf_idx), in place: right-hand sides
  // in, solutions out. One multi-RHS solve with the exact factor, else CG per excitation.
  // Collective.
  void SolveRegion(std::vector<mfem::Vector> &x) const;

  // Full parent vectors from Dirichlet values xd (on dbc_tdofs) and from region-free values
  // x_rf (and optional Dirichlet values), zero elsewhere.
  mfem::Vector DirichletField(const mfem::Vector &xd) const;
  mfem::Vector RegionField(const mfem::Vector &x_rf,
                           const mfem::Vector *xd = nullptr) const;

  // S_E as a parent-space HypreParMatrix (dense Gamma x Gamma block), for factoring the
  // exact condensed region operator.
  std::unique_ptr<mfem::HypreParMatrix> AssembleDenseInterface() const;

  // Batched environment interior solves y_k = A_EE^-1 x_k (inputs/outputs are parent
  // true-DOF vectors; only the environment interior of y_k is kept, the rest is zeroed).
  // The direct path runs multi-RHS blocks through the single factor; the iterative path
  // loops. y_k may alias x_k.
  void ApplyAeeInvMulti(const std::vector<const mfem::Vector *> &X,
                        const std::vector<mfem::Vector *> &Y) const;

  void ApplyAeeInv(const mfem::Vector &x, mfem::Vector &y) const;

  bool EnvMumps() const;

  // Environment interior solve on the environment submesh (the transfer to and from the
  // parent space around an iterative solve).
  void ApplyAeeInvIterative(const mfem::Vector &x_parent, mfem::Vector &y_parent) const;

  // Region preconditioner by geometric multigrid on the region submesh (order >= 2 H1).
  // Returns false to fall back to the single-level preconditioner.
  bool BuildRegionGmg();

  // Region preconditioner apply (~ A_region_free^-1): identity off the region-free DOFs,
  // the region-submesh GMG solve on them.
  void ApplyRegionGmgPc(const mfem::Vector &r, mfem::Vector &z) const;

  // Number the interface true DOFs globally, in rank order, and size the rows of S_E.
  // Collective.
  void EnumerateInterface();

  // Load S_E and the sections of the saved model onto this interface (the head was read
  // and checked by the constructor). Collective.
  void LoadModel(const std::string &model_path, int sig_w);

  // Materialize S_E e_c = (A_env e_c)|_Gamma - (A_env A_EE^-1 (A_env e_c)|_E)|_Gamma by
  // batched multi-RHS environment solves; with_sk also gives S^K = S_E - P^T D P (one more
  // batched solve per column), for a magnetostatic model to save. Collective.
  void MaterializeByBackSolves(bool with_sk);

  // Save S_E, the region mesh and the sections to the model file. Collective.
  void SaveModel(const std::string &model_path, int sig_w);

  // HODLR compression of S_E (and S^K): its well-separated interface blocks are low rank,
  // although S_E is not. Built on rank 0 from the gathered rows, broadcast, and applied
  // replicated. Collective.
  void CompressInterface(double hodlr_tol);

  // The region-condensed solver: CG on A_region_free + S_E, preconditioned by a direct
  // factor, region-submesh multigrid, or AMS / BoomerAMG on A_region_free (with the exact
  // factor, SolveRegion bypasses CG). Returns its description. Collective.
  std::string SetUpRegionSolver();
};

}  // namespace palace

#endif  // PALACE_MODELS_SUBSTRUCTURING_SOLVER_IMPL_HPP
