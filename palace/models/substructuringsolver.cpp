// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolver.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <vector>
#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "fem/substructure.hpp"
#include "linalg/amg.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/gmg.hpp"
#include "linalg/hodlr.hpp"
#include "linalg/iterative.hpp"
#include "linalg/ksp.hpp"
#include "linalg/mumpsschur.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/solver.hpp"
#include "linalg/superlu.hpp"
#include "models/materialoperator.hpp"
#include "models/superconductorsheetoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

// MFEM integrator using the quadrature order of Palace's own operators, so the condensed
// discretization matches the native one also on curved meshes.
// Coefficient given per (domain or boundary) element, for the integrators on those
// elements.
class ElementCoefficient : public mfem::Coefficient
{
  const std::vector<double> &value;

public:
  explicit ElementCoefficient(const std::vector<double> &value) : value(value) {}
  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &) override
  {
    return value[T.ElementNo];
  }
};

template <typename Integrator>
class NativeQuadrature : public Integrator
{
public:
  using Integrator::Integrator;
  void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &T,
                             mfem::DenseMatrix &M) override
  {
    this->IntRule =
        &mfem::IntRules.Get(el.GetGeomType(), fem::DefaultIntegrationOrder::Get(T));
    Integrator::AssembleElementMatrix(el, T, M);
  }
};

// Piecewise-constant (per element attribute) matrix coefficient for anisotropic materials.
// Attributes absent from the map contribute a zero tensor (e.g. environment attributes when
// assembling the region operator).
class PWMatrixCoefficient : public mfem::MatrixCoefficient
{
  const std::map<int, mfem::DenseMatrix> &mats;

public:
  PWMatrixCoefficient(int dim, const std::map<int, mfem::DenseMatrix> &m)
    : mfem::MatrixCoefficient(dim), mats(m)
  {
  }
  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    auto it = mats.find(T.Attribute);
    if (it != mats.end())
    {
      K = it->second;
    }
    else
    {
      K.SetSize(width);
      K = 0.0;
    }
  }
};

// Implicit environment Dirichlet-to-Neumann action on parent true DOFs:
//   y|_Gamma = A_GG x - A_GE A_EE^-1 A_EG x
// via A_env matvecs and one A_EE solve (A_EE = A_env with non-interior true DOFs
// eliminated).
class ImplicitDtN : public mfem::Operator
{
public:
  using ApplyFn = std::function<void(const mfem::Vector &, mfem::Vector &)>;
  ImplicitDtN(mfem::HypreParMatrix &A_env, ApplyFn Aee_inv,
              const std::vector<char> &is_gamma, const std::vector<char> &is_env_int)
    : mfem::Operator(A_env.Height()), A_env(A_env), Aee_inv(std::move(Aee_inv)),
      is_gamma(is_gamma), is_env_int(is_env_int), t(A_env.Height()), rhs(A_env.Height()),
      ye(A_env.Height()), t2(A_env.Height())
  {
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    A_env.Mult(x, t);
    rhs = 0.0;
    for (int i = 0; i < height; i++)
    {
      if (is_env_int[i])
      {
        rhs(i) = t(i);
      }
    }
    ye = 0.0;
    Aee_inv(rhs, ye);
    for (int i = 0; i < height; i++)
    {
      if (!is_env_int[i])
      {
        ye(i) = 0.0;
      }
    }
    A_env.Mult(ye, t2);
    y = 0.0;
    for (int i = 0; i < height; i++)
    {
      if (is_gamma[i])
      {
        y(i) = t(i) - t2(i);
      }
    }
  }

private:
  mfem::HypreParMatrix &A_env;
  ApplyFn Aee_inv;
  const std::vector<char> &is_gamma, &is_env_int;
  mutable mfem::Vector t, rhs, ye, t2;
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
    for (int i = 0; i < height; i++)
    {
      if (is_region_free[i])
      {
        y(i) += td(i);
      }
    }
  }

private:
  mfem::HypreParMatrix &A_region_free;
  const mfem::Operator &dtn;
  const std::vector<char> &is_region_free;
  mutable mfem::Vector td;
};

// Preconditioner adapter: forwards Solver::Mult to a callable (used to route the region
// preconditioner through the region-submesh geometric multigrid). SetOperator is a no-op
// because the underlying multigrid is built once, up front.
class CallableSolver : public Solver<Operator>
{
public:
  CallableSolver(int h, std::function<void(const mfem::Vector &, mfem::Vector &)> f)
    : Solver<Operator>(false), apply(std::move(f))
  {
    height = width = h;
  }
  void SetOperator(const Operator &) override {}
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override { apply(x, y); }

private:
  std::function<void(const mfem::Vector &, mfem::Vector &)> apply;
};

// (MaterializedDtN is defined after the Hodlr helpers below, since it can apply either the
// dense per-rank row blocks or a compressed hierarchical operator.)

// Materialized DtN: applies the precomputed interface operator S_E to a distributed
// interface vector. The global interface vector is gathered by an Allreduce, then either
// each rank multiplies its own dense row block (S_rows) or -- when a compressed Hodlr is
// supplied -- the replicated hierarchical operator is applied and each rank reads back its
// interface rows.
class MaterializedDtN : public mfem::Operator
{
public:
  MaterializedDtN(const std::vector<double> &S_rows, int row_off,
                  const std::vector<int> &gamma_global, int nG_global, MPI_Comm comm,
                  const Hodlr *hodlr = nullptr)
    : mfem::Operator(static_cast<int>(gamma_global.size())), S_rows(S_rows),
      row_off(row_off), gamma_global(gamma_global), nG_global(nG_global), comm(comm),
      hodlr(hodlr)
  {
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    std::vector<double> xl(nG_global, 0.0), xg(nG_global, 0.0);
    for (int i = 0; i < height; i++)
    {
      if (gamma_global[i] >= 0)
      {
        xl[gamma_global[i]] = x(i);
      }
    }
    MPI_Allreduce(xl.data(), xg.data(), nG_global, MPI_DOUBLE, MPI_SUM, comm);
    y = 0.0;
    if (hodlr)
    {
      // Compressed apply: permute to tree order, apply the replicated Hodlr, un-permute,
      // and read back this rank's interface rows.
      std::vector<double> xp(nG_global), yp(nG_global), yg(nG_global);
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
          y(i) = yg[gamma_global[i]];
        }
      }
      return;
    }
    for (int i = 0; i < height; i++)
    {
      if (gamma_global[i] >= 0)
      {
        const int r = gamma_global[i] - row_off;  // owned row block is contiguous
        const double *row = &S_rows[static_cast<std::size_t>(r) * nG_global];
        double s = 0.0;
        for (int j = 0; j < nG_global; j++)
        {
          s += row[j] * xg[j];
        }
        y(i) = s;
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
};

}  // namespace

// Parallel region-condensed static solve (electrostatic H1 or magnetostatic H(curl)):
// region/environment operators assembled on the parent space with domain-restricted
// material coefficients, interface identified in true-DOF space, environment condensed via
// a distributed DtN. Magnetostatics uses a small curl-curl mass regularization (see
// kMagRegularization).
struct SubstructuringSolver::Impl
{
  const IoData &iodata;
  mfem::ParMesh &parent;
  bool magnetostatic;
  std::unique_ptr<mfem::FiniteElementCollection> fec;
  mfem::ParFiniteElementSpace parent_fes;
  int nt;

  std::unique_ptr<mfem::HypreParMatrix> A_region, A_env, A_region_free;
  // Environment interior solve on the environment submesh (Gamma as an essential boundary):
  // avoids the all-identity ranks of the parent-space eliminated operator and enables
  // geometric multigrid on a standard Dirichlet problem.
  mfem::Array<int> ra_arr, ea_arr;
  std::unique_ptr<mfem::ParSubMesh> env_submesh_owned;    // single-level ownership
  std::vector<std::unique_ptr<Mesh>> env_mesh_vec;        // GMG: [0] owns the ParSubMesh
  mfem::ParSubMesh *env_submesh = nullptr;                // raw ptr to the owned submesh
  std::unique_ptr<mfem::ParFiniteElementSpace> env_sfes;  // single-level solve space
  mfem::ParFiniteElementSpace *env_solve_fes = nullptr;   // fespace the solver acts on
  std::vector<std::unique_ptr<mfem::H1_FECollection>> env_fecs;
  std::unique_ptr<FiniteElementSpaceHierarchy> env_hierarchy;
  std::unique_ptr<MultigridOperator> env_mg_op;
  std::vector<mfem::Array<int>>
      env_dbc_lists;  // per-level essential (ParOperator MakeRefs)
  std::unique_ptr<mfem::HypreParMatrix> env_A_ee;  // single-level submesh operator
  std::unique_ptr<KspSolver> env_ksp;
#if defined(MFEM_USE_SUPERLU)
  std::unique_ptr<SuperLUSolver> reg_lu;  // direct A_region_free factorization (region pc)
  std::unique_ptr<mfem::HypreParMatrix>
      env_A_ee_parent;  // parent-space env A_EE (eliminated)
  std::unique_ptr<SuperLUSolver>
      env_lu_parent;  // parent-space direct env solve (no transfer)
#endif
#if defined(MFEM_USE_MUMPS)
  // MUMPS partial factorization of the environment with Gamma as the Schur variables: gives
  // S_E directly and serves every environment interior solve (set up by MaterializeMumps).
  std::unique_ptr<MumpsSchurSolver> env_mumps;
#endif
  // Block size for the batched multi-RHS S_E materialization (SuperLU triangular solves are
  // per-call-overhead bound for a single RHS; ~4x cheaper per column at 32).
  static constexpr int kMaterializeBlock = 32;
  // Largest interface for which the region preconditioner factors the dense Gamma block of
  // the exact condensed operator (dense block + its fill ~ 3 |Gamma|^2 doubles in total).
  static constexpr int kDirectCondensedMaxInterface = 10000;
  // The exact dense-Gamma region factor costs ~|Gamma|^3 once, CG on A_region_free ~50
  // iterations of ~|Gamma|^2 per excitation: the factor pays off from about |Gamma| /
  // kDenseInterfacePerExcitation excitations in a batch (measured break-even ~4 at |Gamma|
  // = 3750, ~1 at |Gamma| ~ 1000).
  static constexpr int kDenseInterfacePerExcitation = 1000;
  bool region_direct = false;          // region pc is a SuperLU factor (created lazily)
  bool region_dense_eligible = false;  // ... and a dense S_E is available for the exact one
  // RHS count of the direct env factor, fixed by its first solve (MFEM's SuperLU wrapper
  // does not allow changing it): kMaterializeBlock when S_E is materialized, else the size
  // of the first excitation batch. Smaller batches are zero-padded to it.
  mutable int env_block = 0;
  mutable mfem::Vector env_pad_in,
      env_pad_out;               // shared zero input / scratch output (padding)
  mfem::Array<int> non_env_int;  // parent true DOFs that are not environment-interior
  bool env_parent_direct = false;
  mfem::Array<int> env_ess;  // solve-space essential true DOFs (Gamma + env Dirichlet)
  mutable mfem::ParGridFunction env_pgf, env_sgf;  // parent / submesh transfer buffers
  mutable mfem::Vector env_srhs, env_ssol;
  // Region-solve geometric multigrid (Phase C): the region preconditioner runs on the
  // region ParSubMesh (region Dirichlet terminals essential, interface Gamma free),
  // avoiding the all-identity ranks of the parent-space region-free operator.
  std::vector<std::unique_ptr<Mesh>> reg_mesh_vec;
  mfem::ParSubMesh *reg_submesh = nullptr;
  mfem::ParFiniteElementSpace *reg_solve_fes = nullptr;
  std::vector<std::unique_ptr<mfem::H1_FECollection>> reg_fecs;
  std::unique_ptr<FiniteElementSpaceHierarchy> reg_hierarchy;
  std::unique_ptr<MultigridOperator> reg_mg_op;
  std::vector<mfem::Array<int>> reg_dbc_lists;
  std::unique_ptr<KspSolver> reg_gmg_ksp;  // inner GMG solve on the region submesh
  mfem::Array<int> reg_ess;
  bool region_gmg = false;
  mutable mfem::ParGridFunction reg_pgf, reg_sgf;
  mutable mfem::Vector reg_srhs, reg_ssol;
  // Energy (QoI) operators. Electrostatic: alias the solve operators. Magnetostatic: pure
  // curl-curl (no mass), so the magnetic energy / inductance is physical.
  std::unique_ptr<mfem::HypreParMatrix> A_region_energy, A_env_energy;
  mfem::HypreParMatrix *K_region_e = nullptr, *K_env_e = nullptr;
  std::vector<char> is_gamma, is_env_int, is_region_free;
  mfem::Array<int> dbc_tdofs;
  std::map<int, std::vector<int>> terminal_tdofs;  // terminal index -> its true DOFs

  std::unique_ptr<RegionCondensedOperator> region_op;
  std::unique_ptr<KspSolver> region_ksp;  // Palace CG + wrapped AMG/AMS on the region block
  std::unique_ptr<ImplicitDtN> dtn;

  // Materialized (reusable) interface operator: replicated dense S_E + load g_E over a
  // global interface enumeration, computed once so region solves need no environment
  // solves.
  std::vector<int> gamma_global;  // owned parent true DOF -> global interface index, or -1
  int nG_global = 0;
  // Distributed dense interface operator S_E: each rank stores the contiguous block of rows
  // [gamma_off, gamma_off + gamma_nloc) it owns (row-major, gamma_nloc x nG_global),
  // instead of the full nG x nG matrix replicated on every rank. Cuts per-rank storage
  // O(nG^2) -> O(nG^2 / P) and parallelizes the apply.
  std::vector<double> S_rows;
  int gamma_off = 0, gamma_nloc = 0;
  // Optional compressed (HODLR) form of S_E, used in place of S_rows when interface
  // compression is requested (replicated across ranks).
  std::unique_ptr<Hodlr> hodlr;
  std::unique_ptr<MaterializedDtN> mat_dtn;

  // Dirichlet-lift modes: the environment condensed onto Gamma + a set of Dirichlet lifts
  // x_k (terminal unit potentials, magnetostatic flux-loop lifts). With t_k = A_env x_k and
  // w_k = A_EE^-1 t_k|_E: the region-solve interface load g_k = (t_k - A_env w_k)|_Gamma
  // (G_rows, owned interface rows: G_rows[r * K + k]). For the energy (operator K_env;
  // equal to A_env for electrostatics, the pure curl-curl for magnetostatics) the
  // environment's exact condensed energy is [u_G; V]^T [S^K, G^K; G^K^T, Cmode] [u_G; V]
  // with S^K = P^T K P = S_E - P^T D P (D = A_env - K_env, P the A-harmonic lifting), G^K
  // (GK_rows; empty when K_env = A_env) and the replicated Cmode. Each mode carries a
  // partition-invariant fingerprint (x^T A_env x, |(A_env x)|_Gamma|) so a saved mode is
  // reused only for the same environment-side lift.
  std::vector<int> mode_ids;
  std::vector<double> mode_fp;  // 2 per mode
  std::vector<double> G_rows, GK_rows, Cmode;
  bool modes_ready = false;
  bool modes_saved = false;  // appended to the model file in this run
  // Magnetostatic energy interface operator S^K (distributed like S_rows); empty for
  // electrostatics (S^K = S_E) or when not materialized (then the environment path is
  // used).
  std::vector<double> SK_rows;
  // Rank-independent presence flags (the distributed arrays are empty on ranks that own no
  // interface DOFs, so their emptiness must not gate collectives).
  bool have_sk = false, have_gk = false;
  std::unique_ptr<Hodlr> hodlr_K;
  std::unique_ptr<MaterializedDtN> mat_dtn_K;
  std::unique_ptr<mfem::HypreParMatrix>
      D_env;               // A_env - K_env = eps M_env (magnetostatics)
  bool env_built = false;  // environment interior operator factored / set up

  Impl(const IoData &iodata, mfem::ParMesh &parent)
    : iodata(iodata), parent(parent),
      magnetostatic(iodata.problem.type == ProblemType::MAGNETOSTATIC),
      fec(magnetostatic
              ? std::unique_ptr<mfem::FiniteElementCollection>(
                    new mfem::ND_FECollection(iodata.solver.order, parent.Dimension()))
              : std::unique_ptr<mfem::FiniteElementCollection>(
                    new mfem::H1_FECollection(iodata.solver.order, parent.Dimension()))),
      parent_fes(&parent, fec.get()), nt(parent_fes.GetTrueVSize())
  {
    const auto &sub = *iodata.solver.substructuring;
    mfem::Array<int> ra(static_cast<int>(sub.region_attributes.size())),
        ea(static_cast<int>(sub.environment_attributes.size()));
    std::copy(sub.region_attributes.begin(), sub.region_attributes.end(), ra.begin());
    std::copy(sub.environment_attributes.begin(), sub.environment_attributes.end(),
              ea.begin());
    ra_arr = ra;
    ea_arr = ea;

    // Interface / region / environment true-DOF markers.
    mfem::Array<int> rm, em, im;
    MarkInterfaceTrueDofs(parent_fes, ra, ea, rm, em, im);

    // Terminal Dirichlet true DOFs (per terminal). The Dirichlet DOF *set* is the union of
    // all terminals and is fixed across excitations; only the prescribed values change, so
    // the interface/interior partition below is excitation-independent.
    const auto &terminals = iodata.boundaries.terminal;
    mfem::Array<int> dir_mark(nt);
    dir_mark = 0;
    {
      const int maxb = parent.bdr_attributes.Size() ? parent.bdr_attributes.Max() : 0;
      for (const auto &[idx, term] : terminals)
      {
        mfem::Array<int> ess_bdr(maxb), ess;
        ess_bdr = 0;
        for (int a : term.attributes)
        {
          if (a >= 1 && a <= maxb)
          {
            ess_bdr[a - 1] = 1;
          }
        }
        parent_fes.GetEssentialTrueDofs(ess_bdr, ess);
        auto &list = terminal_tdofs[idx];
        for (int i = 0; i < ess.Size(); i++)
        {
          dir_mark[ess[i]] = 1;
          list.push_back(ess[i]);
        }
      }
    }
    // Grounded and perfectly conducting boundaries, fixed at zero as in the native
    // operators. Flux-loop films are London sheets (free DOFs with the kinetic sheet term),
    // not Dirichlet boundaries.
    {
      const int maxb = parent.bdr_attributes.Size() ? parent.bdr_attributes.Max() : 0;
      mfem::Array<int> ess_bdr(maxb), ess;
      ess_bdr = 0;
      auto mark = [&](int a)
      {
        if (a >= 1 && a <= maxb)
        {
          ess_bdr[a - 1] = 1;
        }
      };
      for (int a : iodata.boundaries.pec.attributes)
      {
        mark(a);
      }
      parent_fes.GetEssentialTrueDofs(ess_bdr, ess);
      for (int i = 0; i < ess.Size(); i++)
      {
        dir_mark[ess[i]] = 1;
      }
    }
    for (int i = 0; i < nt; i++)
    {
      if (dir_mark[i])
      {
        dbc_tdofs.Append(i);
      }
    }

    is_gamma.assign(nt, 0);
    is_env_int.assign(nt, 0);
    is_region_free.assign(nt, 0);
    for (int i = 0; i < nt; i++)
    {
      if (dir_mark[i])
      {
        continue;
      }
      const bool r = rm[i], e = em[i];
      if (r && e)
      {
        is_gamma[i] = 1;
      }
      else if (e)
      {
        is_env_int[i] = 1;
      }
      if (r)
      {
        is_region_free[i] = 1;  // region-free includes the interface
      }
    }

    // Domain-restricted scalar permittivity coefficients (zero outside the subdomain).
    const int max_attr = parent.attributes.Size() ? parent.attributes.Max() : 1;
    const int dim = parent.Dimension();
    auto in = [](const mfem::Array<int> &s, int a)
    {
      for (int x : s)
      {
        if (x == a)
        {
          return true;
        }
      }
      return false;
    };
    // Reconstruct the (possibly anisotropic) material tensor per attribute from its eigen-
    // decomposition: M_ij = sum_k s[k] v[k]_i v[k]_j. For magnetostatics the curl-curl
    // coefficient is the inverse permeability, so the reconstructed mu tensor is inverted.
    auto tensor = [dim](const config::SymmetricMatrixData<3> &prop, bool invert)
    {
      mfem::DenseMatrix e(dim);
      e = 0.0;
      for (int k = 0; k < 3; k++)
      {
        for (int i = 0; i < dim; i++)
        {
          for (int j = 0; j < dim; j++)
          {
            e(i, j) += prop.s[k] * prop.v[k][i] * prop.v[k][j];
          }
        }
      }
      if (invert)
      {
        e.Invert();
      }
      return e;
    };
    for (const auto &mat : iodata.domains.materials)
    {
      for (int a : mat.attributes)
      {
        if (a < 1 || a > max_attr)
        {
          continue;
        }
        if (in(ra, a))
        {
          region_eps[a] = tensor(magnetostatic ? mat.mu_r : mat.epsilon_r, magnetostatic);
        }
        else if (in(ea, a))
        {
          env_eps[a] = tensor(magnetostatic ? mat.mu_r : mat.epsilon_r, magnetostatic);
        }
      }
    }
    if (magnetostatic)
    {
      SetUpSheets(ra);
      mag_eps = has_sheets ? kMagRegularizationSheets : kMagRegularization;
    }
    A_region = AssembleParent(region_eps, magnetostatic, true, &sheet_region);
    A_env = AssembleParent(env_eps, magnetostatic, true, &sheet_env);
    if (magnetostatic)
    {
      // Energy operators: the pure curl-curl (the London kinetic energy is added
      // separately, see SheetEnergyMatrix).
      A_region_energy = AssembleParent(region_eps, false);
      A_env_energy = AssembleParent(env_eps, false);
      K_region_e = A_region_energy.get();
      K_env_e = A_env_energy.get();
      if (has_sheets)
      {
        M_sheet = AssembleParent({}, false, false, &sheet_all);
        M_sheet_region = AssembleParent({}, false, false, &sheet_region);
        M_sheet_env = AssembleParent({}, false, false, &sheet_env);
      }
    }
    else
    {
      K_region_e = A_region.get();
      K_env_e = A_env.get();
    }
  }

  // Small mass regularization making the singular magnetostatic curl-curl definite (AMS-
  // convergent). Exact for divergence-free-compatible excitations (flux loops); a
  // gauge-free (pseudo-inverse) treatment for general surface currents is a follow-up.
  static constexpr double kMagRegularization = 1.0e-3;
  // With London sheets: the flux-loop sources M_sheet a are orthogonal to the null space
  // of the curl-curl + sheet operator (gradients of potentials constant on the sheets), so
  // a much smaller regularization is safe, and the O(eps^2) energy error is negligible.
  static constexpr double kMagRegularizationSheets = 1.0e-7;
  double mag_eps = kMagRegularization;

  // Environment/region size (global true DOFs) below which a sparse-direct factorization is
  // used by default for the interior solve; above it, fall back to iterative / GMG. Chosen
  // so typical region-in-chip environments factor comfortably while very large domains do
  // not exhaust memory on the SuperLU factorization.
  static constexpr long long kDirectMaxDofs = 2000000;

  // London superconductor sheets (magnetostatics): per local boundary element, 1/L_ksq
  // (with the factor for a cracked sheet, whose two sides each carry half) on the side it
  // belongs to. A sheet face belongs to the region if it bounds a region element (all its
  // DOFs are then region or interface DOFs), else to the environment: independent of the
  // partition and of the region design, and exact either way for a face on Gamma.
  std::vector<double> sheet_region, sheet_env, sheet_all, sheet_none;
  bool has_sheets = false;
  void SetUpSheets(const mfem::Array<int> &ra)
  {
    const int nbe = parent.GetNBE();
    sheet_region.assign(nbe, 0.0);
    sheet_env.assign(nbe, 0.0);
    sheet_all.assign(nbe, 0.0);
    sheet_none.assign(nbe, 0.0);
    std::map<int, double> coef;
    for (const auto &sc : iodata.boundaries.superconductor)
    {
      const double Ls = (sc.Ls > 0.0) ? sc.Ls
                                      : SuperconductorSheetOperator::KineticSheetInductance(
                                            sc.lambda_L, sc.thickness);
      for (int a : sc.attributes)
      {
        coef[a] = 1.0 / (Ls * (iodata.boundaries.cracked_attributes.count(a) ? 2.0 : 1.0));
      }
    }
    int local = 0, global = 0;
    parent.ExchangeFaceNbrData();
    mfem::FaceElementTransformations FET;
    mfem::IsoparametricTransformation T1, T2;
    for (int be = 0; be < nbe; be++)
    {
      auto it = coef.find(parent.GetBdrAttribute(be));
      if (it == coef.end())
      {
        continue;
      }
      BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(be, parent, FET, T1,
                                                                       T2);
      const bool region = ra.Find(FET.Elem1->Attribute) >= 0 ||
                          (FET.Elem2 && ra.Find(FET.Elem2->Attribute) >= 0);
      (region ? sheet_region : sheet_env)[be] = it->second;
      sheet_all[be] = it->second;
      local = 1;
    }
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, parent.GetComm());
    has_sheets = (global > 0);
  }

  // Sheet masses (both sides, and per side), for the flux-loop sources and energies.
  std::unique_ptr<mfem::HypreParMatrix> M_sheet, M_sheet_region, M_sheet_env;

  // Source modes of London flux states (see ComputeSourceModes): ids, fingerprints of the
  // environment sources, the interface load g, the energy coupling h (Gamma rows) and the
  // constants c.
  static constexpr int kSourceModesMagic = 0x31435253;  // "SRC1"
  std::vector<int> src_ids;
  std::vector<double> src_fp, SG_rows, SH_rows, SC;
  bool src_ready = false, src_saved = false;

  // Environment source b^E = M_sheet^E a (zero on the Dirichlet DOFs) and its fingerprint.
  mfem::Vector EnvSource(const mfem::Vector &a, std::array<double, 2> &fp) const
  {
    mfem::Vector b(nt);
    M_sheet_env->Mult(a, b);
    for (int d = 0; d < dbc_tdofs.Size(); d++)
    {
      b(dbc_tdofs[d]) = 0.0;
    }
    double loc[2] = {b * b, a * b}, glob[2];
    MPI_Allreduce(loc, glob, 2, MPI_DOUBLE, MPI_SUM, parent_fes.GetComm());
    fp = {glob[0], glob[1]};
    return b;
  }

  // Source modes of the London flux states with generators a_k. With the environment field
  // u^E = Z u_Gamma + W_k (Z the A-harmonic extension, W_k = A_EE^-1 b^E_k), the
  // environment contributes
  //   g_k = (b^E_k - A_env W_k)|_Gamma                   to the region-condensed load, and
  //   u_i,G^T S^K u_j,G + u_i,G^T h_j + h_i^T u_j,G + c_ij   to the energies, with
  //   h_k = Z^T [(A_env - D_env) W_k - b^E_k],
  //   c_ij = W_i^T K_env W_j + (W_i - a_i)^T M_sheet^E (W_j - a_j).
  // Two batched environment solves per state.
  void ComputeSourceModes(const std::vector<int> &ids, const std::vector<mfem::Vector> &a)
  {
    EnsureEnv();
    MPI_Comm comm = parent_fes.GetComm();
    const int K = static_cast<int>(ids.size());
    src_ids = ids;
    src_fp.assign(2 * K, 0.0);
    std::vector<mfem::Vector> b(K), W(K, mfem::Vector(nt)), q(K, mfem::Vector(nt)),
        z(K, mfem::Vector(nt));
    std::vector<const mfem::Vector *> X(K);
    std::vector<mfem::Vector *> Y(K);
    for (int k = 0; k < K; k++)
    {
      std::array<double, 2> fp;
      b[k] = EnvSource(a[k], fp);
      src_fp[2 * k] = fp[0];
      src_fp[2 * k + 1] = fp[1];
      W[k] = 0.0;
      X[k] = &b[k];
      Y[k] = &W[k];
    }
    if (K > 0)
    {
      ApplyAeeInvMulti(X, Y);
    }
    mfem::Vector t(nt);
    auto rows = [&](std::vector<double> &out, int k, const mfem::Vector &v)
    {
      for (int i = 0; i < nt; i++)
      {
        if (is_gamma[i])
        {
          out[static_cast<std::size_t>(gamma_global[i] - gamma_off) * K + k] = v(i);
        }
      }
    };
    SG_rows.assign(static_cast<std::size_t>(gamma_nloc) * K, 0.0);
    for (int k = 0; k < K; k++)
    {
      A_env->Mult(W[k], t);
      t.Neg();
      t += b[k];  // g_k = b^E_k|_Gamma - (A_env W_k)|_Gamma
      rows(SG_rows, k, t);
      // q_k = (A_env - D_env) W_k - b^E_k, then z_k = A_EE^-1 q_k|_E.
      A_env->Mult(W[k], q[k]);
      Denv().Mult(W[k], t);
      q[k] -= t;
      q[k] -= b[k];
      z[k] = 0.0;
      X[k] = &q[k];
      Y[k] = &z[k];
    }
    if (K > 0)
    {
      ApplyAeeInvMulti(X, Y);
    }
    SH_rows.assign(static_cast<std::size_t>(gamma_nloc) * K, 0.0);
    for (int k = 0; k < K; k++)
    {
      A_env->Mult(z[k], t);
      t.Neg();
      t += q[k];  // h_k = q_k|_Gamma - (A_env z_k)|_Gamma
      rows(SH_rows, k, t);
    }
    SC.assign(static_cast<std::size_t>(K) * K, 0.0);
    std::vector<mfem::Vector> KW(K, mfem::Vector(nt)), MD(K, mfem::Vector(nt)),
        D(K, mfem::Vector(nt));
    for (int k = 0; k < K; k++)
    {
      K_env_e->Mult(W[k], KW[k]);
      D[k] = W[k];
      D[k] -= a[k];
      M_sheet_env->Mult(D[k], MD[k]);
    }
    for (int i = 0; i < K; i++)
    {
      for (int j = 0; j < K; j++)
      {
        SC[static_cast<std::size_t>(i) * K + j] = (W[i] * KW[j]) + (D[i] * MD[j]);
      }
    }
    if (K > 0)
    {
      MPI_Allreduce(MPI_IN_PLACE, SC.data(), K * K, MPI_DOUBLE, MPI_SUM, comm);
    }
    src_ready = true;
    src_saved = false;
  }

  // Whether the source modes cover the given generators (same ids and environment source
  // fingerprints); col[j] receives the mode column of request j. Collective.
  bool SourceModesMatch(const std::vector<int> &ids, const std::vector<mfem::Vector> &a,
                        std::vector<int> &col) const
  {
    col.assign(ids.size(), -1);
    bool ok = src_ready;
    for (std::size_t j = 0; j < ids.size(); j++)
    {
      for (std::size_t k = 0; k < src_ids.size(); k++)
      {
        if (src_ids[k] == ids[j])
        {
          col[j] = static_cast<int>(k);
        }
      }
      std::array<double, 2> fp;
      (void)EnvSource(a[j], fp);  // collective: always evaluated
      if (col[j] < 0)
      {
        ok = false;
        continue;
      }
      for (int q = 0; q < 2; q++)
      {
        const double ref = src_fp[2 * col[j] + q];
        if (std::abs(fp[q] - ref) >
            1.0e-9 * std::max(std::abs(ref), std::abs(fp[q])) + 1e-300)
        {
          ok = false;
        }
      }
    }
    return ok;
  }

  void AppendSourceModes(const std::string &path)
  {
    const int K = static_cast<int>(src_ids.size());
    const std::vector<double> G_full = GatherRows(SG_rows, K);
    const std::vector<double> H_full = GatherRows(SH_rows, K);
    if (Mpi::Root(parent_fes.GetComm()))
    {
      std::ofstream f(path, std::ios::binary | std::ios::app);
      f.write(reinterpret_cast<const char *>(&kSourceModesMagic), sizeof(int));
      f.write(reinterpret_cast<const char *>(&K), sizeof(int));
      f.write(reinterpret_cast<const char *>(src_ids.data()), sizeof(int) * K);
      f.write(reinterpret_cast<const char *>(src_fp.data()), sizeof(double) * 2 * K);
      f.write(reinterpret_cast<const char *>(G_full.data()),
              sizeof(double) * G_full.size());
      f.write(reinterpret_cast<const char *>(H_full.data()),
              sizeof(double) * H_full.size());
      f.write(reinterpret_cast<const char *>(SC.data()), sizeof(double) * SC.size());
    }
    src_saved = true;
  }

  std::unique_ptr<mfem::HypreParMatrix>
  AssembleParent(const std::map<int, mfem::DenseMatrix> &coef_by_attr, bool with_mass,
                 bool with_curl = true, const std::vector<double> *sheet = nullptr)
  {
    PWMatrixCoefficient coef(parent.Dimension(), coef_by_attr);
    mfem::ParBilinearForm a(&parent_fes);
    ElementCoefficient sheet_coef(sheet ? *sheet : sheet_none);
    if (magnetostatic && sheet && has_sheets)
    {
      // London kinetic sheet term (1/L_ksq) A_t . v_t on the given boundary elements.
      a.AddBoundaryIntegrator(
          new NativeQuadrature<mfem::VectorFEMassIntegrator>(sheet_coef));
    }
    if (magnetostatic)
    {
      if (with_curl)
      {
        a.AddDomainIntegrator(new NativeQuadrature<mfem::CurlCurlIntegrator>(coef));
      }
      const int max_attr = parent.attributes.Size() ? parent.attributes.Max() : 1;
      mfem::Vector mass(max_attr);
      mass = 0.0;
      if (with_mass)
      {
        for (const auto &[a_attr, m] : coef_by_attr)
        {
          mass(a_attr - 1) = mag_eps;
        }
      }
      mfem::PWConstCoefficient mcoef(mass);  // outlives Assemble() below
      if (with_mass)
      {
        a.AddDomainIntegrator(new NativeQuadrature<mfem::VectorFEMassIntegrator>(mcoef));
      }
      a.Assemble();
      a.Finalize();
      return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
    }
    a.AddDomainIntegrator(new NativeQuadrature<mfem::DiffusionIntegrator>(coef));
    a.Assemble();
    a.Finalize();
    return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
  }

  std::map<int, mfem::DenseMatrix> region_eps, env_eps;

  // Environment interior solver: A_EE^-1 (env-interior with Gamma + env Dirichlet held
  // fixed), built on the environment submesh so every rank owns real DOFs (no all-identity
  // ranks).
  // Global number of environment closure DOFs that a direct factorization would carry.
  long long EnvDirectSize() const
  {
    long long env_loc = 0;
    for (int i = 0; i < nt; i++)
    {
      if (is_env_int[i] || is_gamma[i])
      {
        env_loc++;
      }
    }
    long long env_glob = 0;
    MPI_Allreduce(&env_loc, &env_glob, 1, MPI_LONG_LONG, MPI_SUM, parent_fes.GetComm());
    return env_glob;
  }

#if defined(MFEM_USE_MUMPS)
  // Materialize S_E by ONE MUMPS partial factorization (Schur complement on Gamma) of A_env
  // with every DOF outside the environment closure-minus-Dirichlet (E + Gamma) eliminated
  // to identity; the factor is kept as the environment solver. Fills S_rows (owned
  // interface rows). Collective.
  void MaterializeMumps()
  {
    MPI_Comm comm = parent_fes.GetComm();
    mfem::Array<int> other;
    for (int i = 0; i < nt; i++)
    {
      if (!is_env_int[i] && !is_gamma[i])
      {
        other.Append(i);
      }
    }
    mfem::HypreParMatrix A_sch(*A_env);
    {
      std::unique_ptr<mfem::HypreParMatrix> e(A_sch.EliminateRowsCols(other));
    }
    // Schur variables: interface true DOFs in global interface order (replicated).
    const HYPRE_BigInt tstart = parent_fes.GetMyTDofOffset();
    std::vector<HYPRE_BigInt> mine(gamma_nloc), g2t(nG_global);
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        mine[gamma_global[i] - gamma_off] = tstart + i;
      }
    }
    std::vector<int> cnt, disp;
    RowLayout(1, cnt, disp);
    MPI_Allgatherv(mine.data(), gamma_nloc, HYPRE_MPI_BIG_INT, g2t.data(), cnt.data(),
                   disp.data(), HYPRE_MPI_BIG_INT, comm);
    env_mumps = std::make_unique<MumpsSchurSolver>(
        A_sch, g2t, iodata.solver.substructuring->factorization_tol);
    // The Schur is symmetric (column-major == row-major): scatter its rows to their owners.
    ScatterRows(env_mumps->Schur(), nG_global, S_rows);
    env_built = true;
  }
#endif

  // Build (factor) the environment interior solver on first use: an online run that loads
  // S_E and the terminal modes never needs it unless environment fields are requested.
  void EnsureEnv()
  {
    if (!env_built)
    {
      BuildEnvSubmeshSolver();
      env_built = true;
    }
  }

  // Recover the environment interior of full fields (region + Dirichlet parts set),
  // batched: u_E = -A_EE^-1 (A_env u)|_E.
  void RecoverEnvInterior(std::vector<mfem::Vector> &u)
  {
    EnsureEnv();
    const int n = static_cast<int>(u.size());
    std::vector<mfem::Vector> rhs(n, mfem::Vector(nt));
    std::vector<const mfem::Vector *> X(n);
    std::vector<mfem::Vector *> Y(n);
    for (int k = 0; k < n; k++)
    {
      A_env->Mult(u[k], rhs[k]);
      rhs[k].Neg();
      X[k] = &rhs[k];
      Y[k] = &rhs[k];  // solve in place
    }
    if (n > 0)
    {
      ApplyAeeInvMulti(X, Y);
    }
    for (int k = 0; k < n; k++)
    {
      for (int i = 0; i < nt; i++)
      {
        if (is_env_int[i])
        {
          u[k](i) = rhs[k](i);
        }
      }
    }
  }

  // Unit Dirichlet lift of terminal idx: 1 on its true DOFs.
  mfem::Vector TerminalMode(int idx) const
  {
    mfem::Vector x(nt);
    x = 0.0;
    for (int d : terminal_tdofs.at(idx))
    {
      x(d) = 1.0;
    }
    return x;
  }

  // Whether the energy operator differs from the solve operator (magnetostatic curl-curl
  // vs. curl-curl + mass), i.e. the energy needs the D-corrected condensed forms.
  bool EnergyDiffers() const { return magnetostatic; }

  const mfem::HypreParMatrix &Denv()
  {
    if (!D_env)
    {
      D_env = AssembleParent(env_eps, true, false);  // eps M_env, exactly A_env - K_env
    }
    return *D_env;
  }

  // Partition-invariant fingerprint of a lift's environment side: (x^T A_env x,
  // |(A_env x)|_Gamma|_2). Collective.
  std::array<double, 2> LiftFingerprint(const mfem::Vector &x, const mfem::Vector &t) const
  {
    double loc[2] = {x * t, 0.0};
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        loc[1] += t(i) * t(i);
      }
    }
    double glob[2];
    MPI_Allreduce(loc, glob, 2, MPI_DOUBLE, MPI_SUM, parent_fes.GetComm());
    return {glob[0], std::sqrt(glob[1])};
  }

  // Modes (see the members) of the lifts x_k from batched environment solves: K solves,
  // plus K more for the energy correction when it differs. Lifts that do not reach the
  // environment interior skip the solve.
  void ComputeModes(const std::vector<int> &ids, const std::vector<mfem::Vector> &x)
  {
    EnsureEnv();
    MPI_Comm comm = parent_fes.GetComm();
    const int K = static_cast<int>(ids.size());
    mode_ids = ids;
    std::vector<mfem::Vector> t(K, mfem::Vector(nt)), w(K, mfem::Vector(nt));
    std::vector<int> touches(K, 0);
    mode_fp.assign(2 * K, 0.0);
    for (int k = 0; k < K; k++)
    {
      A_env->Mult(x[k], t[k]);
      const auto fp = LiftFingerprint(x[k], t[k]);
      mode_fp[2 * k] = fp[0];
      mode_fp[2 * k + 1] = fp[1];
      for (int i = 0; i < nt && !touches[k]; i++)
      {
        touches[k] = (is_env_int[i] && t[k](i) != 0.0);
      }
    }
    if (K > 0)
    {
      MPI_Allreduce(MPI_IN_PLACE, touches.data(), K, MPI_INT, MPI_MAX, comm);
    }
    auto solve_touching = [&](std::vector<mfem::Vector> &in, std::vector<mfem::Vector> &out)
    {
      std::vector<const mfem::Vector *> X;
      std::vector<mfem::Vector *> Y;
      for (int k = 0; k < K; k++)
      {
        out[k] = 0.0;
        if (touches[k])
        {
          X.push_back(&in[k]);
          Y.push_back(&out[k]);
        }
      }
      if (!X.empty())
      {
        ApplyAeeInvMulti(X, Y);
      }
    };
    solve_touching(t, w);
    // g = (a - A_env b)|_Gamma, C_kl = x_k . a_l - t_k . b_l  for (a, b) = (t, w) [solve]
    // or (v, z) [energy correction].
    mfem::Vector tmp(nt);
    auto couplings = [&](const std::vector<mfem::Vector> &a,
                         const std::vector<mfem::Vector> &b, std::vector<double> &Grows,
                         std::vector<double> &C)
    {
      Grows.assign(static_cast<std::size_t>(gamma_nloc) * K, 0.0);
      for (int k = 0; k < K; k++)
      {
        A_env->Mult(b[k], tmp);
        for (int i = 0; i < nt; i++)
        {
          if (is_gamma[i])
          {
            Grows[static_cast<std::size_t>(gamma_global[i] - gamma_off) * K + k] =
                a[k](i) - tmp(i);
          }
        }
      }
      C.assign(static_cast<std::size_t>(K) * K, 0.0);
      for (int k = 0; k < K; k++)
      {
        for (int l = 0; l < K; l++)
        {
          C[static_cast<std::size_t>(k) * K + l] = (x[k] * a[l]) - (t[k] * b[l]);
        }
      }
      if (K > 0)
      {
        MPI_Allreduce(MPI_IN_PLACE, C.data(), K * K, MPI_DOUBLE, MPI_SUM, comm);
      }
    };
    couplings(t, w, G_rows, Cmode);
    GK_rows.clear();
    have_gk = EnergyDiffers();
    if (have_gk)
    {
      // P^T D P blocks: p_k = x_k - w_k (the A-harmonic lift), v_k = D p_k,
      // z_k = A_EE^-1 v_k|_E; subtract them from the A-blocks.
      std::vector<mfem::Vector> v(K, mfem::Vector(nt)), z(K, mfem::Vector(nt));
      for (int k = 0; k < K; k++)
      {
        tmp = x[k];
        tmp -= w[k];
        Denv().Mult(tmp, v[k]);
        touches[k] = 1;  // v_k is generally nonzero in E
      }
      solve_touching(v, z);
      std::vector<double> GD, CD;
      couplings(v, z, GD, CD);
      GK_rows = G_rows;
      for (std::size_t q = 0; q < GK_rows.size(); q++)
      {
        GK_rows[q] -= GD[q];
      }
      for (std::size_t q = 0; q < Cmode.size(); q++)
      {
        Cmode[q] -= CD[q];
      }
    }
    modes_ready = true;
    modes_saved = false;
  }

  // Whether the current modes cover the given lifts (same ids + environment fingerprints);
  // col[j] receives the mode column of request j. Collective.
  bool ModesMatch(const std::vector<int> &ids, const std::vector<mfem::Vector> &x,
                  std::vector<int> &col) const
  {
    col.assign(ids.size(), -1);
    bool ok = modes_ready;
    mfem::Vector t(nt);
    for (std::size_t j = 0; j < ids.size(); j++)
    {
      for (std::size_t k = 0; k < mode_ids.size(); k++)
      {
        if (mode_ids[k] == ids[j])
        {
          col[j] = static_cast<int>(k);
        }
      }
      A_env->Mult(x[j], t);
      const auto fp = LiftFingerprint(x[j], t);  // collective: always evaluated
      if (col[j] < 0)
      {
        ok = false;
        continue;
      }
      for (int q = 0; q < 2; q++)
      {
        const double ref = mode_fp[2 * col[j] + q];
        if (std::abs(fp[q] - ref) >
            1.0e-9 * std::max(std::abs(ref), std::abs(fp[q])) + 1e-300)
        {
          ok = false;
        }
      }
    }
    return ok;
  }

  // Per-rank counts/displacements of an interface-row-distributed array (gamma_nloc * w
  // each).
  void RowLayout(int w, std::vector<int> &cnt, std::vector<int> &disp) const
  {
    MPI_Comm comm = parent_fes.GetComm();
    const int nranks = Mpi::Size(comm);
    cnt.assign(nranks, 0);
    disp.assign(nranks, 0);
    const int mine = gamma_nloc * w;
    MPI_Allgather(&mine, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    for (int r = 1; r < nranks; r++)
    {
      disp[r] = disp[r - 1] + cnt[r - 1];
    }
  }

  // Gather an interface-row-distributed array (width w) to rank 0 (full, global row order).
  std::vector<double> GatherRows(const std::vector<double> &rows, int w) const
  {
    std::vector<int> cnt, disp;
    RowLayout(w, cnt, disp);
    MPI_Comm comm = parent_fes.GetComm();
    const bool root = (Mpi::Rank(comm) == 0);
    std::vector<double> full;
    if (root)
    {
      full.assign(static_cast<std::size_t>(nG_global) * w, 0.0);
    }
    MPI_Gatherv(rows.data(), gamma_nloc * w, MPI_DOUBLE, root ? full.data() : nullptr,
                cnt.data(), disp.data(), MPI_DOUBLE, 0, comm);
    return full;
  }

  // Scatter a full (rank 0) interface-row array of width w onto the owned rows.
  void ScatterRows(const std::vector<double> &full, int w, std::vector<double> &rows) const
  {
    std::vector<int> cnt, disp;
    RowLayout(w, cnt, disp);
    MPI_Comm comm = parent_fes.GetComm();
    rows.assign(static_cast<std::size_t>(gamma_nloc) * w, 0.0);
    MPI_Scatterv(Mpi::Rank(comm) == 0 ? full.data() : nullptr, cnt.data(), disp.data(),
                 MPI_DOUBLE, rows.data(), gamma_nloc * w, MPI_DOUBLE, 0, comm);
  }

  // Model-file sections following S_E: [kEnvFpMagic][environment fingerprint, 4 doubles],
  // [kSKenMagic][S^K nG x nG] (magnetostatics) and
  // [kModesMagic][K][ids][fingerprints 2K][has_GK][G nG x K][G^K nG x K if has_GK][Cmode
  // KxK].
  static constexpr int kSKenMagic = 0x4e454b53;     // "SKEN"
  static constexpr int kModesMagic = 0x32444f4d;    // "MOD2"
  static constexpr int kEnvFpMagicV1 = 0x46564e45;  // "ENVF" (no Dirichlet count)
  static constexpr int kEnvFpMagic = 0x32564e45;    // "ENV2"
  // Fingerprint of the environment: global counts of the environment-closure DOFs and of
  // its Dirichlet DOFs, and the trace and Frobenius norm of A_env. Independent of the
  // region (A_env has no region contribution), of the partition and of the DOF numbering
  // (and of H(curl) orientation signs), but it changes with the environment mesh,
  // materials, order, physics or Dirichlet boundaries. Collective.
  std::array<double, 4> EnvironmentFingerprint() const
  {
    mfem::SparseMatrix diag;
    A_env->GetDiag(diag);
    double trace = 0.0;
    for (int i = 0; i < diag.Height(); i++)
    {
      trace += diag.Elem(i, i);
    }
    // Environment-closure DOFs: free (interior and interface) and Dirichlet (touched by an
    // environment element: a nonzero diagonal of A_env, which is assembled with a zero
    // coefficient on the region).
    double loc[3] = {trace, 0.0, 0.0};
    for (int i = 0; i < nt; i++)
    {
      loc[1] += (is_env_int[i] || is_gamma[i]) ? 1.0 : 0.0;
    }
    for (int d = 0; d < dbc_tdofs.Size(); d++)
    {
      const int i = dbc_tdofs[d];
      loc[2] += (diag.Elem(i, i) != 0.0) ? 1.0 : 0.0;
    }
    double glob[3];
    MPI_Allreduce(loc, glob, 3, MPI_DOUBLE, MPI_SUM, parent_fes.GetComm());
    return {glob[1], glob[2], glob[0], A_env->FNorm()};
  }

  // Saved-model environment fingerprint (valid when have_env_fp).
  std::array<double, 4> saved_env_fp = {0.0, 0.0, 0.0, 0.0};
  bool have_env_fp = false;

  void AppendEnvFingerprint(const std::string &path) const
  {
    const std::array<double, 4> fp = EnvironmentFingerprint();  // collective
    if (Mpi::Root(parent_fes.GetComm()))
    {
      std::ofstream f(path, std::ios::binary | std::ios::app);
      f.write(reinterpret_cast<const char *>(&kEnvFpMagic), sizeof(int));
      f.write(reinterpret_cast<const char *>(fp.data()), sizeof(double) * fp.size());
    }
  }

  // Abort if a loaded model was condensed from a different environment than the current
  // one.
  void CheckEnvFingerprint() const
  {
    const std::array<double, 4> fp = EnvironmentFingerprint();  // collective
    if (!have_env_fp)
    {
      return;  // model written before the fingerprint existed
    }
    bool ok = (fp[0] == saved_env_fp[0] && fp[1] == saved_env_fp[1]);
    for (int q = 2; q < 4; q++)
    {
      ok = ok && std::abs(fp[q] - saved_env_fp[q]) <=
                     1.0e-10 * std::max(std::abs(fp[q]), std::abs(saved_env_fp[q]));
    }
    MFEM_VERIFY(ok,
                "The environment differs from the one the saved substructuring model was "
                "condensed from (environment mesh, materials, boundary conditions, "
                "order or problem type changed): rerun in \"Offline\" mode to "
                "condense it again!");
  }

  void AppendSK(const std::string &path) const
  {
    const std::vector<double> full = GatherRows(SK_rows, nG_global);
    if (Mpi::Root(parent_fes.GetComm()))
    {
      std::ofstream f(path, std::ios::binary | std::ios::app);
      f.write(reinterpret_cast<const char *>(&kSKenMagic), sizeof(int));
      f.write(reinterpret_cast<const char *>(full.data()), sizeof(double) * full.size());
    }
  }
  void AppendModes(const std::string &path)
  {
    const int K = static_cast<int>(mode_ids.size());
    const int has_gk = have_gk ? 1 : 0;
    const std::vector<double> G_full = GatherRows(G_rows, K);
    const std::vector<double> GK_full =
        has_gk ? GatherRows(GK_rows, K) : std::vector<double>();
    if (Mpi::Root(parent_fes.GetComm()))
    {
      std::ofstream f(path, std::ios::binary | std::ios::app);
      f.write(reinterpret_cast<const char *>(&kModesMagic), sizeof(int));
      f.write(reinterpret_cast<const char *>(&K), sizeof(int));
      f.write(reinterpret_cast<const char *>(mode_ids.data()), sizeof(int) * K);
      f.write(reinterpret_cast<const char *>(mode_fp.data()), sizeof(double) * 2 * K);
      f.write(reinterpret_cast<const char *>(&has_gk), sizeof(int));
      f.write(reinterpret_cast<const char *>(G_full.data()),
              sizeof(double) * G_full.size());
      f.write(reinterpret_cast<const char *>(GK_full.data()),
              sizeof(double) * GK_full.size());
      f.write(reinterpret_cast<const char *>(Cmode.data()), sizeof(double) * Cmode.size());
    }
    modes_saved = true;
  }

  // Load the optional sections at byte offset `pos` (collective). Interface rows are
  // re-ordered onto the online numbering like S_E: row g <- saved row perm[g] times sgn[g]
  // (and S^K columns likewise). Reading stops at the first unknown / missing section; a
  // later mode section supersedes an earlier one. Modes absent here are recomputed on
  // demand.
  void LoadSections(const std::string &path, std::streamoff pos,
                    const std::vector<int> &perm, const std::vector<double> &sgn)
  {
    MPI_Comm comm = parent_fes.GetComm();
    const bool root = (Mpi::Rank(comm) == 0);
    const int nG = nG_global;
    std::ifstream f;
    if (root)
    {
      f.open(path, std::ios::binary);
      f.seekg(pos);
    }
    while (true)
    {
      int magic = 0;
      if (root && !f.read(reinterpret_cast<char *>(&magic), sizeof(int)))
      {
        magic = 0;
      }
      MPI_Bcast(&magic, 1, MPI_INT, 0, comm);
      if (magic == kEnvFpMagic)
      {
        if (root)
        {
          f.read(reinterpret_cast<char *>(saved_env_fp.data()),
                 sizeof(double) * saved_env_fp.size());
        }
        MPI_Bcast(saved_env_fp.data(), 4, MPI_DOUBLE, 0, comm);
        have_env_fp = true;
      }
      else if (magic == kEnvFpMagicV1)
      {
        double legacy[3];  // older fingerprint without the Dirichlet count: not checked
        if (root)
        {
          f.read(reinterpret_cast<char *>(legacy), sizeof(legacy));
        }
      }
      else if (magic == kSKenMagic)
      {
        std::vector<double> full;
        if (root)
        {
          std::vector<double> off(static_cast<std::size_t>(nG) * nG);
          f.read(reinterpret_cast<char *>(off.data()), sizeof(double) * off.size());
          full.assign(off.size(), 0.0);
          for (int a = 0; a < nG; a++)
          {
            for (int b = 0; b < nG; b++)
            {
              full[static_cast<std::size_t>(a) * nG + b] =
                  sgn[a] * sgn[b] * off[static_cast<std::size_t>(perm[a]) * nG + perm[b]];
            }
          }
        }
        ScatterRows(full, nG, SK_rows);
        have_sk = true;
      }
      else if (magic == kSourceModesMagic)
      {
        int K = 0;
        std::vector<int> ids;
        std::vector<double> fp, C, G_full, H_full;
        if (root)
        {
          f.read(reinterpret_cast<char *>(&K), sizeof(int));
        }
        MPI_Bcast(&K, 1, MPI_INT, 0, comm);
        ids.resize(K);
        fp.resize(2 * K);
        C.resize(static_cast<std::size_t>(K) * K);
        if (root)
        {
          f.read(reinterpret_cast<char *>(ids.data()), sizeof(int) * K);
          f.read(reinterpret_cast<char *>(fp.data()), sizeof(double) * 2 * K);
        }
        auto read_rows = [&](std::vector<double> &full)
        {
          if (root)
          {
            std::vector<double> off(static_cast<std::size_t>(nG) * K);
            f.read(reinterpret_cast<char *>(off.data()), sizeof(double) * off.size());
            full.assign(off.size(), 0.0);
            for (int g = 0; g < nG; g++)
            {
              for (int k = 0; k < K; k++)
              {
                full[static_cast<std::size_t>(g) * K + k] =
                    sgn[g] * off[static_cast<std::size_t>(perm[g]) * K + k];
              }
            }
          }
        };
        read_rows(G_full);
        read_rows(H_full);
        if (root)
        {
          f.read(reinterpret_cast<char *>(C.data()), sizeof(double) * C.size());
        }
        MPI_Bcast(ids.data(), K, MPI_INT, 0, comm);
        MPI_Bcast(fp.data(), 2 * K, MPI_DOUBLE, 0, comm);
        MPI_Bcast(C.data(), K * K, MPI_DOUBLE, 0, comm);
        ScatterRows(G_full, K, SG_rows);
        ScatterRows(H_full, K, SH_rows);
        src_ids = ids;
        src_fp = fp;
        SC = C;
        src_ready = true;
        src_saved = true;  // already in the file
      }
      else if (magic == kModesMagic)
      {
        int hdr[2] = {0, 0};  // K, has_gk
        std::vector<int> ids;
        std::vector<double> fp, C, G_full, GK_full;
        if (root)
        {
          f.read(reinterpret_cast<char *>(&hdr[0]), sizeof(int));
          ids.resize(hdr[0]);
          fp.resize(2 * hdr[0]);
          f.read(reinterpret_cast<char *>(ids.data()), sizeof(int) * hdr[0]);
          f.read(reinterpret_cast<char *>(fp.data()), sizeof(double) * 2 * hdr[0]);
          f.read(reinterpret_cast<char *>(&hdr[1]), sizeof(int));
        }
        MPI_Bcast(hdr, 2, MPI_INT, 0, comm);
        const int K = hdr[0];
        ids.resize(K);
        fp.resize(2 * K);
        C.resize(static_cast<std::size_t>(K) * K);
        auto read_rows = [&](std::vector<double> &full)
        {
          if (root)
          {
            std::vector<double> off(static_cast<std::size_t>(nG) * K);
            f.read(reinterpret_cast<char *>(off.data()), sizeof(double) * off.size());
            full.assign(off.size(), 0.0);
            for (int g = 0; g < nG; g++)
            {
              for (int k = 0; k < K; k++)
              {
                full[static_cast<std::size_t>(g) * K + k] =
                    sgn[g] * off[static_cast<std::size_t>(perm[g]) * K + k];
              }
            }
          }
        };
        read_rows(G_full);
        if (hdr[1])
        {
          read_rows(GK_full);
        }
        if (root)
        {
          f.read(reinterpret_cast<char *>(C.data()), sizeof(double) * C.size());
        }
        MPI_Bcast(ids.data(), K, MPI_INT, 0, comm);
        MPI_Bcast(fp.data(), 2 * K, MPI_DOUBLE, 0, comm);
        MPI_Bcast(C.data(), K * K, MPI_DOUBLE, 0, comm);
        ScatterRows(G_full, K, G_rows);
        GK_rows.clear();
        have_gk = (hdr[1] != 0);
        if (have_gk)
        {
          ScatterRows(GK_full, K, GK_rows);
        }
        mode_ids = ids;
        mode_fp = fp;
        Cmode = C;
        modes_ready = true;
        modes_saved = true;  // already in the file
      }
      else
      {
        break;
      }
    }
  }

  void BuildEnvSubmeshSolver()
  {
    const int mg_levels = iodata.solver.linear.mg_max_levels;
#if defined(MFEM_USE_SUPERLU)
    // Direct A_EE (default when the environment fits): factor A_env with the
    // non-env-interior true DOFs eliminated, in PARENT space. SuperLU tolerates the
    // all-identity ranks that rule out a parent-space *iterative* solve, and no ParSubMesh
    // transfer is needed per solve (measured ~40% of the S_E materialization). Order >= 2
    // H(curl) always takes this path: the transfer mishandles higher-order tetrahedral
    // edge/face orientation across a cut.
    const bool direct_fits = (EnvDirectSize() <= kDirectMaxDofs);
    if (direct_fits || (magnetostatic && iodata.solver.order > 1))
    {
      non_env_int.SetSize(0);
      for (int i = 0; i < nt; i++)
      {
        if (!is_env_int[i])
        {
          non_env_int.Append(i);
        }
      }
      env_A_ee_parent = std::make_unique<mfem::HypreParMatrix>(*A_env);
      {
        std::unique_ptr<mfem::HypreParMatrix> e(
            env_A_ee_parent->EliminateRowsCols(non_env_int));
      }
      env_lu_parent = std::make_unique<SuperLUSolver>(iodata, parent_fes.GetComm(), 0);
      env_lu_parent->SetOperator(*env_A_ee_parent);
      env_parent_direct = true;
      return;
    }
#endif
    // Iterative fallback (very large environment, or no SuperLU): environment submesh
    // solve.
    const bool use_gmg = !magnetostatic && iodata.solver.order > 1 && mg_levels > 1;

    // Create the environment submesh. For the GMG path it must be owned by a Palace Mesh
    // (for the CEED attribute maps + the FE-space hierarchy); otherwise a standalone
    // ParSubMesh suffices for the MFEM assembly + transfer.
    if (use_gmg)
    {
      env_mesh_vec.clear();
      env_mesh_vec.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParSubMesh>(
          mfem::ParSubMesh::CreateFromDomain(parent, ea_arr))));
      env_submesh = dynamic_cast<mfem::ParSubMesh *>(&env_mesh_vec[0]->Get());
      env_mesh_vec[0]->RebuildCeedAttributes();
    }
    else
    {
      env_submesh_owned = std::make_unique<mfem::ParSubMesh>(
          mfem::ParSubMesh::CreateFromDomain(parent, ea_arr));
      env_submesh = env_submesh_owned.get();
    }
    env_pgf.SetSpace(&parent_fes);

    if (!(use_gmg && BuildEnvGmg()))
    {
      // Single-level submesh Dirichlet solve (order 1, magnetostatic H(curl), or GMG
      // fallback): standalone FE space + wrapped AMG / AMS.
      env_sfes = std::make_unique<mfem::ParFiniteElementSpace>(env_submesh, fec.get());
      env_solve_fes = env_sfes.get();
      env_sgf.SetSpace(env_solve_fes);
      ComputeEnvEss();
      env_A_ee = AssembleEnvSubmesh();
      {
        std::unique_ptr<mfem::HypreParMatrix> e(env_A_ee->EliminateRowsCols(env_ess));
      }
      MPI_Comm comm = env_solve_fes->GetComm();
      {
        std::unique_ptr<Solver<Operator>> pc;
        if (magnetostatic)
        {
          auto ams = std::make_unique<mfem::HypreAMS>(env_solve_fes);
          ams->SetPrintLevel(0);
          pc = std::make_unique<MfemWrapperSolver<Operator>>(std::move(ams), true, false,
                                                             false);
        }
        else
        {
          pc = std::make_unique<MfemWrapperSolver<Operator>>(
              std::make_unique<BoomerAmgSolver>(1, 1, true, 0), true, false, false);
        }
        auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
        pcg->SetInitialGuess(false);
        pcg->SetRelTol(1.0e-12);
        pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
        pcg->SetMaxIter(1000);
        env_ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(pc));
        env_ksp->SetOperators(*env_A_ee, *env_A_ee);
      }
    }

    env_srhs.SetSize(env_solve_fes->GetTrueVSize());
    env_ssol.SetSize(env_solve_fes->GetTrueVSize());
  }

  // Essential submesh DOFs = the parent non-(environment-interior) DOFs (interface Gamma +
  // environment Dirichlet), mapped onto the solve space by transferring the parent marker.
  // Robust to Palace inserting material-interface boundary elements at Gamma.
  void ComputeEnvEss()
  {
    mfem::Vector t(nt);
    for (int i = 0; i < nt; i++)
    {
      t(i) = is_env_int[i] ? 0.0 : 1.0;
    }
    mfem::ParGridFunction pg(&parent_fes), sg(env_solve_fes);
    pg.SetFromTrueDofs(t);
    sg = 0.0;
    env_submesh->Transfer(pg, sg);
    mfem::Vector st(env_solve_fes->GetTrueVSize());
    sg.GetTrueDofs(st);
    env_ess.SetSize(0);
    for (int i = 0; i < st.Size(); i++)
    {
      if (std::abs(st(i)) > 0.5)  // fabs: ND transfer may flip the marker's sign
      {
        env_ess.Append(i);
      }
    }
  }

  // Higher-order H1 environment Dirichlet solve via geometric p-multigrid on the submesh.
  // Returns false (falling back to the single-level solve) if a hierarchy cannot be built.
  bool BuildEnvGmg()
  {
    const int order = iodata.solver.order;
    const int dim = parent.Dimension();
    const int mg_levels = iodata.solver.linear.mg_max_levels;

    // Identify essential boundary attributes on the submesh (Gamma + env Dirichlet) from a
    // marker transfer onto a scratch order-p space: a boundary attribute is essential iff
    // all of its true DOFs lie in the essential set.
    auto scratch = std::make_unique<mfem::ParFiniteElementSpace>(env_submesh, fec.get());
    std::set<int> ess_set;
    {
      mfem::Vector t(nt);
      for (int i = 0; i < nt; i++)
      {
        t(i) = is_env_int[i] ? 0.0 : 1.0;
      }
      mfem::ParGridFunction pg(&parent_fes), sg(scratch.get());
      pg.SetFromTrueDofs(t);
      sg = 0.0;
      env_submesh->Transfer(pg, sg);
      mfem::Vector st(scratch->GetTrueVSize());
      sg.GetTrueDofs(st);
      for (int i = 0; i < st.Size(); i++)
      {
        if (std::abs(st(i)) > 0.5)
        {
          ess_set.insert(i);
        }
      }
    }
    mfem::Array<int> ess_attr;
    // The essential-attribute decision must be identical on every rank (otherwise ranks
    // diverge between the GMG and single-level paths and deadlock in the collectives
    // below). An attribute is essential iff, globally, it has DOFs and none of them are
    // non-essential.
    MPI_Comm comm = env_submesh->GetComm();
    const int lbmax =
        env_submesh->bdr_attributes.Size() ? env_submesh->bdr_attributes.Max() : 0;
    int bmax = 0;
    MPI_Allreduce(&lbmax, &bmax, 1, MPI_INT, MPI_MAX, comm);
    std::vector<int> all_in(bmax, 1), has_dofs(bmax, 0);
    for (int a = 1; a <= bmax; a++)
    {
      mfem::Array<int> m(bmax);
      m = 0;
      m[a - 1] = 1;
      mfem::Array<int> adofs;
      scratch->GetEssentialTrueDofs(m, adofs);
      if (adofs.Size() > 0)
      {
        has_dofs[a - 1] = 1;
      }
      for (int d : adofs)
      {
        if (!ess_set.count(d))
        {
          all_in[a - 1] = 0;
          break;
        }
      }
    }
    std::vector<int> g_all_in(bmax), g_has(bmax);
    MPI_Allreduce(all_in.data(), g_all_in.data(), bmax, MPI_INT, MPI_LAND, comm);
    MPI_Allreduce(has_dofs.data(), g_has.data(), bmax, MPI_INT, MPI_LOR, comm);
    for (int a = 1; a <= bmax; a++)
    {
      if (g_all_in[a - 1] && g_has[a - 1])
      {
        ess_attr.Append(a);
      }
    }
    if (ess_attr.Size() == 0)
    {
      return false;
    }

    // p-multigrid hierarchy on the submesh, with ess_attr as the Dirichlet boundary.
    env_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
        order, dim, mg_levels, iodata.solver.linear.mg_coarsening, false);
    std::vector<mfem::Array<int>> &dbc_lists = env_dbc_lists;
    dbc_lists.clear();
    env_hierarchy = std::make_unique<FiniteElementSpaceHierarchy>(
        fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
            mg_levels, env_mesh_vec, env_fecs, &ess_attr, &dbc_lists));
    if (env_hierarchy->GetNumLevels() < 2)
    {
      env_hierarchy.reset();
      env_fecs.clear();
      return false;
    }
    env_solve_fes = &env_hierarchy->GetFinestFESpace().Get();

    // Ceed-consistent material coefficient (cf. divfree.cpp): built via a MaterialOperator
    // on the submesh so the attribute-to-material map matches Palace's local CEED numbering
    // (a hand-built coefficient keyed by global attribute assembles the interior to zero).
    MaterialOperator env_mat_op(iodata, *env_mesh_vec[0]);
    MaterialPropertyCoefficient coef(env_mat_op.GetAttributeToMaterial(),
                                     env_mat_op.GetPermittivityReal());
    BilinearForm a(env_hierarchy->GetFinestFESpace());
    a.AddDomainIntegrator<DiffusionIntegrator>(coef);
    auto a_vec = a.Assemble(*env_hierarchy, false);

    const std::size_t nl = env_hierarchy->GetNumLevels();
    env_mg_op = std::make_unique<MultigridOperator>(nl);
    for (std::size_t l = 0; l < nl; l++)
    {
      auto &fes_l = env_hierarchy->GetFESpaceAtLevel(l);
      auto A_l = std::make_unique<ParOperator>(std::move(a_vec[l]), fes_l);
      A_l->SetEssentialTrueDofs(dbc_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
      env_mg_op->AddOperator(std::move(A_l));
    }
    env_ess = dbc_lists.back();  // finest essential, consistent with the operators

    auto amg = std::make_unique<MfemWrapperSolver<Operator>>(
        std::make_unique<BoomerAmgSolver>(1, 1, true, 0));
    amg->SetDropSmallEntries(false);
    auto gmg = std::make_unique<GeometricMultigridSolver<Operator>>(
        iodata, comm, std::move(amg), env_hierarchy->GetProlongationOperators());
    auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
    pcg->SetInitialGuess(false);
    pcg->SetRelTol(1.0e-12);
    pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
    pcg->SetMaxIter(1000);
    env_ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(gmg));
    env_ksp->SetOperators(*env_mg_op, *env_mg_op);

    env_sgf.SetSpace(env_solve_fes);
    return true;
  }

  std::unique_ptr<mfem::HypreParMatrix> AssembleEnvSubmesh()
  {
    MFEM_VERIFY(!has_sheets, "Superconductor sheets in substructuring need the direct "
                             "environment factorization (the environment is too large)!");
    PWMatrixCoefficient coef(parent.Dimension(), env_eps);
    mfem::ParBilinearForm a(env_sfes.get());
    if (magnetostatic)
    {
      a.AddDomainIntegrator(new NativeQuadrature<mfem::CurlCurlIntegrator>(coef));
      const int am = env_submesh->attributes.Size() ? env_submesh->attributes.Max() : 1;
      mfem::Vector mass(am);
      mass = 0.0;
      for (const auto &[attr, t] : env_eps)
      {
        mass(attr - 1) = mag_eps;
      }
      mfem::PWConstCoefficient mcoef(mass);
      a.AddDomainIntegrator(new NativeQuadrature<mfem::VectorFEMassIntegrator>(mcoef));
      a.Assemble();
      a.Finalize();
      return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
    }
    a.AddDomainIntegrator(new NativeQuadrature<mfem::DiffusionIntegrator>(coef));
    a.Assemble();
    a.Finalize();
    return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
  }

  // Apply A_EE^-1 to a parent true-DOF vector (nonzero on the environment interior):
  // transfer to the submesh, solve the Dirichlet problem, transfer back, and restrict to
  // the environment interior.
  // Interface selection E_Gamma (parent true DOFs x global interface index): E(i, g_i) = 1
  // for each interface true DOF. Each rank's interface indices are its own contiguous
  // block, so the matrix is local (diag block only).
  std::unique_ptr<mfem::HypreParMatrix> AssembleInterfaceSelection() const
  {
    MPI_Comm comm = parent_fes.GetComm();
    std::vector<int> I(nt + 1, 0);
    std::vector<HYPRE_BigInt> J;
    std::vector<double> V;
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        J.push_back(gamma_global[i]);
        V.push_back(1.0);
      }
      I[i + 1] = static_cast<int>(J.size());
    }
    // Column (interface) partitioning in the layout HYPRE expects for this build.
    std::vector<HYPRE_BigInt> cols;
    if (HYPRE_AssumedPartitionCheck())
    {
      cols = {gamma_off, gamma_off + gamma_nloc, nG_global};
    }
    else
    {
      const int nranks = Mpi::Size(comm);
      std::vector<int> cnt(nranks);
      MPI_Allgather(&gamma_nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
      cols.assign(nranks + 1, 0);
      for (int r = 0; r < nranks; r++)
      {
        cols[r + 1] = cols[r] + cnt[r];
      }
    }
    return std::make_unique<mfem::HypreParMatrix>(
        comm, nt, parent_fes.GlobalTrueVSize(), static_cast<HYPRE_BigInt>(nG_global),
        I.data(), J.data(), V.data(), parent_fes.GetTrueDofOffsets(), cols.data());
  }

  // Create the direct region preconditioner on the first solve batch (n excitations): the
  // exact dense-Gamma factor A_region_free + S_E when the batch is large enough to amortize
  // it, else A_region_free (CG then needs ~50 iterations). Decided once: a SuperLU object
  // that was never factored must not be destroyed (SuperLU_DIST teardown bug). Collective;
  // n is the same on all ranks.
  void EnsureRegionFactor(int n)
  {
#if defined(MFEM_USE_SUPERLU)
    if (!region_direct || reg_lu)
    {
      return;
    }
    reg_lu = std::make_unique<SuperLUSolver>(iodata, parent_fes.GetComm(), 0);
    if (region_dense_eligible &&
        static_cast<long long>(std::max(n, 1)) * kDenseInterfacePerExcitation >= nG_global)
    {
      std::unique_ptr<mfem::HypreParMatrix> S_mat = AssembleDenseInterface();
      std::unique_ptr<mfem::HypreParMatrix> P(
          mfem::ParAdd(A_region_free.get(), S_mat.get()));
      S_mat.reset();
      reg_lu->SetOperator(*P);  // copies P into SuperLU's distributed format
    }
    else
    {
      reg_lu->SetOperator(*A_region_free);
    }
#endif
  }

  // S_E as a parent-space HypreParMatrix (dense Gamma x Gamma block, zero elsewhere), for
  // factoring the exact condensed region operator. Each rank's S_E rows are its own
  // interface true DOFs; columns map interface index -> global true DOF via a replicated
  // table.
  std::unique_ptr<mfem::HypreParMatrix> AssembleDenseInterface() const
  {
    MPI_Comm comm = parent_fes.GetComm();
    const int nG = nG_global, nloc = gamma_nloc;
    const HYPRE_BigInt tstart = parent_fes.GetMyTDofOffset();
    std::vector<HYPRE_BigInt> mine(nloc), g2t(nG);
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        mine[gamma_global[i] - gamma_off] = tstart + i;
      }
    }
    const int nranks = Mpi::Size(comm);
    std::vector<int> cnt(nranks), disp(nranks);
    MPI_Allgather(&nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    for (int r = 0, acc = 0; r < nranks; r++)
    {
      disp[r] = acc;
      acc += cnt[r];
    }
    MPI_Allgatherv(mine.data(), nloc, HYPRE_MPI_BIG_INT, g2t.data(), cnt.data(),
                   disp.data(), HYPRE_MPI_BIG_INT, comm);
    std::vector<int> I(nt + 1, 0);
    std::vector<HYPRE_BigInt> J;
    std::vector<double> V;
    J.reserve(static_cast<std::size_t>(nloc) * nG);
    V.reserve(static_cast<std::size_t>(nloc) * nG);
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        const double *row =
            &S_rows[static_cast<std::size_t>(gamma_global[i] - gamma_off) * nG];
        for (int j = 0; j < nG; j++)
        {
          J.push_back(g2t[j]);
          V.push_back(row[j]);
        }
      }
      I[i + 1] = static_cast<int>(J.size());
    }
    const HYPRE_BigInt glob = parent_fes.GlobalTrueVSize();
    return std::make_unique<mfem::HypreParMatrix>(comm, nt, glob, glob, I.data(), J.data(),
                                                  V.data(), parent_fes.GetTrueDofOffsets(),
                                                  parent_fes.GetTrueDofOffsets());
  }

  // Batched environment interior solves y_k = A_EE^-1 x_k (inputs/outputs are parent
  // true-DOF vectors; only the environment interior of y_k is kept, the rest is zeroed).
  // The direct path runs multi-RHS blocks through the single factor; the iterative path
  // loops. y_k may alias x_k.
  void ApplyAeeInvMulti(const std::vector<const mfem::Vector *> &X,
                        const std::vector<mfem::Vector *> &Y) const
  {
    MFEM_ASSERT(X.size() == Y.size(), "ApplyAeeInvMulti size mismatch!");
    MFEM_VERIFY(env_built, "Environment solver used before EnsureEnv()!");
#if defined(MFEM_USE_MUMPS)
    if (env_mumps)
    {
      // Internal solve of the Schur factorization: A_EE on the interior (Gamma held at 0,
      // eliminated DOFs are identity), then masked to the interior like the SuperLU path.
      env_mumps->SolveInternal(X, Y);
      for (auto *y : Y)
      {
        for (int i = 0; i < nt; i++)
        {
          if (!is_env_int[i])
          {
            (*y)(i) = 0.0;
          }
        }
      }
      return;
    }
#endif
#if defined(MFEM_USE_SUPERLU)
    if (env_parent_direct)
    {
      if (env_block == 0)
      {
        env_block = std::max(1, std::min(kMaterializeBlock, static_cast<int>(X.size())));
        env_pad_in.SetSize(nt);
        env_pad_in = 0.0;
        env_pad_out.SetSize(nt);
      }
      // The factored operator is block diagonal [A_EE, 0; 0, I] (non-interior rows/cols
      // eliminated), so the interior solution depends only on the interior rhs: the inputs
      // are passed as-is and only the output is masked. The tail block is padded with a
      // shared zero rhs (SuperLU fixes the RHS count at the first solve).
      mfem::Array<const mfem::Vector *> Xp(env_block);
      mfem::Array<mfem::Vector *> Yp(env_block);
      const int n = static_cast<int>(X.size());
      for (int c0 = 0; c0 < n; c0 += env_block)
      {
        const int nb = std::min(env_block, n - c0);
        for (int k = 0; k < env_block; k++)
        {
          Xp[k] = (k < nb) ? X[c0 + k] : &env_pad_in;
          Yp[k] = (k < nb) ? Y[c0 + k] : &env_pad_out;
          Yp[k]->SetSize(nt);
        }
        if (env_block == 1)
        {
          env_lu_parent->Mult(*Xp[0], *Yp[0]);
        }
        else
        {
          env_lu_parent->ArrayMult(Xp, Yp);
        }
        for (int k = 0; k < nb; k++)
        {
          mfem::Vector &y = *Y[c0 + k];
          for (int i = 0; i < nt; i++)
          {
            if (!is_env_int[i])
            {
              y(i) = 0.0;
            }
          }
        }
      }
      return;
    }
#endif
    for (std::size_t k = 0; k < X.size(); k++)
    {
      ApplyAeeInv(*X[k], *Y[k]);
    }
  }

  void ApplyAeeInv(const mfem::Vector &x_parent, mfem::Vector &y_parent) const
  {
#if defined(MFEM_USE_MUMPS)
    if (env_mumps)
    {
      ApplyAeeInvMulti({&x_parent}, {&y_parent});
      return;
    }
#endif
#if defined(MFEM_USE_SUPERLU)
    if (env_parent_direct)
    {
      // Parent-space direct env interior solve (no submesh transfer), via the batched path
      // so the factor's fixed RHS count is respected.
      ApplyAeeInvMulti({&x_parent}, {&y_parent});
      return;
    }
#endif
    env_pgf.SetFromTrueDofs(x_parent);
    env_sgf = 0.0;
    env_submesh->Transfer(env_pgf, env_sgf);
    env_sgf.GetTrueDofs(env_srhs);
    for (int i = 0; i < env_ess.Size(); i++)
    {
      env_srhs(env_ess[i]) = 0.0;
    }
    env_ssol = 0.0;
    env_ksp->Mult(env_srhs, env_ssol);
    env_sgf.SetFromTrueDofs(env_ssol);
    env_pgf = 0.0;
    env_submesh->Transfer(env_sgf, env_pgf);
    y_parent.SetSize(nt);
    env_pgf.GetTrueDofs(y_parent);
    for (int i = 0; i < nt; i++)
    {
      if (!is_env_int[i])
      {
        y_parent(i) = 0.0;
      }
    }
  }

  // Region-solve geometric multigrid preconditioner on the region submesh (order>=2 H1).
  // Returns false to fall back to the single-level wrapped BoomerAMG preconditioner.
  bool BuildRegionGmg()
  {
    const int order = iodata.solver.order;
    const int dim = parent.Dimension();
    const int mg_levels = iodata.solver.linear.mg_max_levels;
    if (magnetostatic || order <= 1 || mg_levels <= 1)
    {
      return false;
    }

    reg_mesh_vec.clear();
    reg_mesh_vec.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParSubMesh>(
        mfem::ParSubMesh::CreateFromDomain(parent, ra_arr))));
    reg_submesh = dynamic_cast<mfem::ParSubMesh *>(&reg_mesh_vec[0]->Get());
    reg_mesh_vec[0]->RebuildCeedAttributes();
    reg_pgf.SetSpace(&parent_fes);

    // Essential submesh DOFs = the region Dirichlet terminals (Gamma stays free), found by
    // transferring the parent Dirichlet marker onto a scratch order-p space.
    auto scratch = std::make_unique<mfem::ParFiniteElementSpace>(reg_submesh, fec.get());
    std::set<int> ess_set;
    {
      mfem::Vector t(nt);
      t = 0.0;
      for (int i = 0; i < dbc_tdofs.Size(); i++)
      {
        t(dbc_tdofs[i]) = 1.0;
      }
      mfem::ParGridFunction pg(&parent_fes), sg(scratch.get());
      pg.SetFromTrueDofs(t);
      sg = 0.0;
      reg_submesh->Transfer(pg, sg);
      mfem::Vector st(scratch->GetTrueVSize());
      sg.GetTrueDofs(st);
      for (int i = 0; i < st.Size(); i++)
      {
        if (std::abs(st(i)) > 0.5)
        {
          ess_set.insert(i);
        }
      }
    }
    // Global essential-attribute decision (identical on every rank, see BuildEnvGmg).
    MPI_Comm comm = reg_submesh->GetComm();
    const int lbmax =
        reg_submesh->bdr_attributes.Size() ? reg_submesh->bdr_attributes.Max() : 0;
    int bmax = 0;
    MPI_Allreduce(&lbmax, &bmax, 1, MPI_INT, MPI_MAX, comm);
    std::vector<int> all_in(bmax, 1), has_dofs(bmax, 0);
    for (int a = 1; a <= bmax; a++)
    {
      mfem::Array<int> m(bmax);
      m = 0;
      m[a - 1] = 1;
      mfem::Array<int> adofs;
      scratch->GetEssentialTrueDofs(m, adofs);
      if (adofs.Size() > 0)
      {
        has_dofs[a - 1] = 1;
      }
      for (int d : adofs)
      {
        if (!ess_set.count(d))
        {
          all_in[a - 1] = 0;
          break;
        }
      }
    }
    std::vector<int> g_all_in(bmax), g_has(bmax);
    MPI_Allreduce(all_in.data(), g_all_in.data(), bmax, MPI_INT, MPI_LAND, comm);
    MPI_Allreduce(has_dofs.data(), g_has.data(), bmax, MPI_INT, MPI_LOR, comm);
    mfem::Array<int> ess_attr;
    for (int a = 1; a <= bmax; a++)
    {
      if (g_all_in[a - 1] && g_has[a - 1])
      {
        ess_attr.Append(a);
      }
    }
    if (ess_attr.Size() == 0)
    {
      return false;  // pure-Neumann region block: keep the single-level preconditioner
    }

    reg_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
        order, dim, mg_levels, iodata.solver.linear.mg_coarsening, false);
    reg_dbc_lists.clear();
    reg_hierarchy = std::make_unique<FiniteElementSpaceHierarchy>(
        fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
            mg_levels, reg_mesh_vec, reg_fecs, &ess_attr, &reg_dbc_lists));
    if (reg_hierarchy->GetNumLevels() < 2)
    {
      reg_hierarchy.reset();
      reg_fecs.clear();
      return false;
    }
    reg_solve_fes = &reg_hierarchy->GetFinestFESpace().Get();

    MaterialOperator reg_mat_op(iodata, *reg_mesh_vec[0]);
    MaterialPropertyCoefficient coef(reg_mat_op.GetAttributeToMaterial(),
                                     reg_mat_op.GetPermittivityReal());
    BilinearForm a(reg_hierarchy->GetFinestFESpace());
    a.AddDomainIntegrator<DiffusionIntegrator>(coef);
    auto a_vec = a.Assemble(*reg_hierarchy, false);

    const std::size_t nl = reg_hierarchy->GetNumLevels();
    reg_mg_op = std::make_unique<MultigridOperator>(nl);
    for (std::size_t l = 0; l < nl; l++)
    {
      auto &fes_l = reg_hierarchy->GetFESpaceAtLevel(l);
      auto A_l = std::make_unique<ParOperator>(std::move(a_vec[l]), fes_l);
      A_l->SetEssentialTrueDofs(reg_dbc_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
      reg_mg_op->AddOperator(std::move(A_l));
    }
    reg_ess = reg_dbc_lists.back();

    auto amg = std::make_unique<MfemWrapperSolver<Operator>>(
        std::make_unique<BoomerAmgSolver>(1, 1, true, 0));
    amg->SetDropSmallEntries(false);
    auto gmg = std::make_unique<GeometricMultigridSolver<Operator>>(
        iodata, comm, std::move(amg), reg_hierarchy->GetProlongationOperators());
    auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
    pcg->SetInitialGuess(false);
    pcg->SetRelTol(1.0e-10);
    pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
    pcg->SetMaxIter(1000);
    reg_gmg_ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(gmg));
    reg_gmg_ksp->SetOperators(*reg_mg_op, *reg_mg_op);

    reg_sgf.SetSpace(reg_solve_fes);
    reg_srhs.SetSize(reg_solve_fes->GetTrueVSize());
    reg_ssol.SetSize(reg_solve_fes->GetTrueVSize());
    region_gmg = true;
    return true;
  }

  // Region preconditioner apply (approximates A_region_free^-1): identity on the
  // non-region-free DOFs (A_region_free is identity there), and the region-submesh GMG
  // solve on the region-free DOFs (transfer parent -> submesh, solve, transfer back,
  // restrict).
  void ApplyRegionGmgPc(const mfem::Vector &r, mfem::Vector &z) const
  {
    z.SetSize(nt);
    for (int i = 0; i < nt; i++)
    {
      z(i) = is_region_free[i] ? 0.0 : r(i);
    }
    reg_pgf.SetFromTrueDofs(r);
    reg_sgf = 0.0;
    reg_submesh->Transfer(reg_pgf, reg_sgf);
    reg_sgf.GetTrueDofs(reg_srhs);
    for (int i = 0; i < reg_ess.Size(); i++)
    {
      reg_srhs(reg_ess[i]) = 0.0;
    }
    reg_ssol = 0.0;
    reg_gmg_ksp->Mult(reg_srhs, reg_ssol);
    reg_sgf.SetFromTrueDofs(reg_ssol);
    reg_pgf = 0.0;
    reg_submesh->Transfer(reg_sgf, reg_pgf);
    mfem::Vector zr(nt);
    reg_pgf.GetTrueDofs(zr);
    for (int i = 0; i < nt; i++)
    {
      if (is_region_free[i])
      {
        z(i) = zr(i);
      }
    }
  }
};

SubstructuringSolver::SubstructuringSolver(const IoData &iodata,
                                           const std::vector<std::unique_ptr<Mesh>> &mesh)
  : impl(std::make_unique<Impl>(iodata, mesh.back()->Get()))
{
  MFEM_VERIFY(iodata.solver.substructuring,
              "SubstructuringSolver requires a Solver.Substructuring configuration!");
}

SubstructuringSolver::~SubstructuringSolver() = default;

void SubstructuringSolver::CondenseEnvironment()
{
  // The environment interior solver is built lazily (EnsureEnv): an online run that loads
  // S_E and the terminal modes computes capacitance without ever factoring the environment.
  Impl *pi = impl.get();
  impl->dtn = std::make_unique<ImplicitDtN>(
      *impl->A_env, [pi](const mfem::Vector &x, mfem::Vector &y) { pi->ApplyAeeInv(x, y); },
      impl->is_gamma, impl->is_env_int);

  // A_region restricted to region-free true DOFs (non-region-free pinned to identity).
  impl->A_region_free = std::make_unique<mfem::HypreParMatrix>(*impl->A_region);
  mfem::Array<int> non_rfree;
  for (int i = 0; i < impl->nt; i++)
  {
    if (!impl->is_region_free[i])
    {
      non_rfree.Append(i);
    }
  }
  {
    std::unique_ptr<mfem::HypreParMatrix> tmp(
        impl->A_region_free->EliminateRowsCols(non_rfree));
  }

  // Materialize the interface operator S_E and load g_E once, so region solves reuse them
  // without any environment solves. Global interface enumeration via MPI_Exscan.
  MPI_Comm comm = impl->parent_fes.GetComm();
  int nloc = 0;
  for (int i = 0; i < impl->nt; i++)
  {
    if (impl->is_gamma[i])
    {
      nloc++;
    }
  }
  int off = 0;
  MPI_Exscan(&nloc, &off, 1, MPI_INT, MPI_SUM, comm);
  impl->nG_global = 0;
  MPI_Allreduce(&nloc, &impl->nG_global, 1, MPI_INT, MPI_SUM, comm);
  impl->gamma_global.assign(impl->nt, -1);
  {
    int c = off;
    for (int i = 0; i < impl->nt; i++)
    {
      if (impl->is_gamma[i])
      {
        impl->gamma_global[i] = c++;
      }
    }
  }

  const int nG = impl->nG_global;
  impl->gamma_off = off;
  impl->gamma_nloc = nloc;
  impl->S_rows.assign(static_cast<std::size_t>(nloc) * nG, 0.0);
  // Per-rank row counts/displacements for gather/scatter of the distributed S_E (row
  // blocks).
  const int nranks = Mpi::Size(comm);
  std::vector<int> row_cnt(nranks), row_disp(nranks);
  {
    std::vector<int> nloc_all(nranks);
    MPI_Allgather(&nloc, 1, MPI_INT, nloc_all.data(), 1, MPI_INT, comm);
    int acc = 0;
    for (int r = 0; r < nranks; r++)
    {
      row_cnt[r] = nloc_all[r] * nG;
      row_disp[r] = acc;
      acc += row_cnt[r];
    }
  }

  // Interface true-DOF geometric signature (replicated), for re-ordering a saved S_E onto
  // the current interface after region re-meshing / re-partitioning (fixed Gamma). H1: DOF
  // coordinates (width 3). H(curl): per-edge moments dof(e_b) (edge vector) and dof(x_a
  // e_b) (edge_b * mid_a), width 12 -- position + signed orientation, matching edge DOFs up
  // to an orientation flip.
  const int sdim = impl->parent.Dimension();
  const int sig_w = impl->magnetostatic ? 12 : 3;
  auto gamma_sig = [&]()
  {
    std::vector<double> loc(static_cast<std::size_t>(nG) * sig_w, 0.0),
        glob(static_cast<std::size_t>(nG) * sig_w, 0.0);
    mfem::ParGridFunction gf(&impl->parent_fes);
    Vector td(impl->nt);
    auto stamp = [&](int slot, mfem::Coefficient *sc, mfem::VectorCoefficient *vc)
    {
      if (sc)
      {
        gf.ProjectCoefficient(*sc);
      }
      else
      {
        gf.ProjectCoefficient(*vc);
      }
      gf.GetTrueDofs(td);
      for (int i = 0; i < impl->nt; i++)
      {
        if (impl->is_gamma[i])
        {
          loc[static_cast<std::size_t>(impl->gamma_global[i]) * sig_w + slot] = td(i);
        }
      }
    };
    if (!impl->magnetostatic)
    {
      for (int d = 0; d < sdim; d++)
      {
        mfem::FunctionCoefficient xc([d](const mfem::Vector &x) { return x(d); });
        stamp(d, &xc, nullptr);
      }
    }
    else
    {
      // Edge (tangential-integral) moments: o_b = dof(e_b) = edge vector; and
      // M[a][b] = dof(x_a e_b) = (edge vector)_b * midpoint_a. Together (12 numbers) these
      // uniquely identify an edge DOF up to an orientation flip (both are linear in the
      // tangent, so a flip negates all of them).
      for (int b = 0; b < 3; b++)
      {
        mfem::Vector e(3);
        e = 0.0;
        e(b) = 1.0;
        mfem::VectorConstantCoefficient ec(e);
        stamp(b, nullptr, &ec);
      }
      for (int a = 0; a < 3; a++)
      {
        for (int b = 0; b < 3; b++)
        {
          mfem::VectorFunctionCoefficient xc(3,
                                             [a, b](const mfem::Vector &x, mfem::Vector &v)
                                             {
                                               v = 0.0;
                                               v(b) = x(a);
                                             });
          stamp(3 + a * 3 + b, nullptr, &xc);
        }
      }
    }
    MPI_Allreduce(loc.data(), glob.data(), nG * sig_w, MPI_DOUBLE, MPI_SUM, comm);
    return glob;
  };

  // Offline/online: the interface enumeration above is cheap and deterministic to rebuild,
  // but materializing S_E costs |Gamma| environment solves. In Online mode with a saved
  // model, load S_E instead; in Offline mode with a path set, materialize and save it.
  const auto &subcfg = *impl->iodata.solver.substructuring;
  const bool online = (subcfg.mode == SubstructuringMode::ONLINE);
  const std::string &model_path = subcfg.save_model;
  const int rank = Mpi::Rank(comm);
  bool loaded = false;
  if (online && !model_path.empty())
  {
    int nG_file = -1;
    if (rank == 0)
    {
      std::ifstream f(model_path, std::ios::binary);
      if (f.good())
      {
        f.read(reinterpret_cast<char *>(&nG_file), sizeof(int));
      }
    }
    MPI_Bcast(&nG_file, 1, MPI_INT, 0, comm);
    MFEM_VERIFY(nG_file == nG,
                "Saved substructuring model interface size ("
                    << nG_file << ") does not match this run (" << nG
                    << "); the interface (Gamma) must be identical between the offline and "
                       "online runs.");
    int sig_type = 0;  // 0: none (identical mesh), 1: H1 coords, 2: H(curl) edge signature
    if (rank == 0)
    {
      std::ifstream f(model_path, std::ios::binary);
      f.seekg(sizeof(int));
      f.read(reinterpret_cast<char *>(&sig_type), sizeof(int));
    }
    MPI_Bcast(&sig_type, 1, MPI_INT, 0, comm);
    const int file_w = (sig_type == 2) ? 12 : (sig_type == 1 ? 3 : 0);
    std::vector<double> saved_sig;
    if (file_w > 0)
    {
      saved_sig.assign(static_cast<std::size_t>(nG) * file_w, 0.0);
      if (rank == 0)
      {
        std::ifstream f(model_path, std::ios::binary);
        f.seekg(2 * sizeof(int));
        f.read(reinterpret_cast<char *>(saved_sig.data()), sizeof(double) * nG * file_w);
      }
      MPI_Bcast(saved_sig.data(), nG * file_w, MPI_DOUBLE, 0, comm);
    }
    // Rank 0 holds the full saved S_E transiently, re-orders it onto the current interface
    // numbering, then scatters contiguous row blocks -- no replicated full matrix persists.
    std::vector<double> S_full;
    if (rank == 0)
    {
      S_full.assign(static_cast<std::size_t>(nG) * nG, 0.0);
      std::ifstream f(model_path, std::ios::binary);
      f.seekg(static_cast<std::streamoff>(2 * sizeof(int) + sizeof(double) * nG * file_w));
      f.read(reinterpret_cast<char *>(S_full.data()), sizeof(double) * nG * nG);
    }
    // Signature match -> (perm, sgn): online interface index g corresponds to saved index
    // perm[g] with orientation sgn[g] (sgn=1 for H1; +/-1 for H(curl)). Computed on all
    // ranks from the replicated signatures. gamma_global stays fresh (contiguous), and S_E
    // is re-indexed to the online order: S_on[i][j] = sgn_i sgn_j S_off[perm_i][perm_j].
    std::vector<int> perm(nG);
    std::iota(perm.begin(), perm.end(), 0);  // identity unless re-ordered by signature
    std::vector<double> sgn(nG, 1.0);
    if (file_w > 0)
    {
      const std::vector<double> cur = gamma_sig();
      const bool signed_match = (sig_type == 2);
      double worst = 0.0;
      for (int g = 0; g < nG; g++)
      {
        int best = 0;
        double bs = 1.0, bd = 1e300;
        for (int s = 0; s < nG; s++)
        {
          double dp = 0.0, dm = 0.0;
          for (int d = 0; d < file_w; d++)
          {
            const double a = cur[static_cast<std::size_t>(g) * file_w + d];
            const double b = saved_sig[static_cast<std::size_t>(s) * file_w + d];
            dp += (a - b) * (a - b);
            if (signed_match)
            {
              dm += (a + b) * (a + b);
            }
          }
          if (dp < bd)
          {
            bd = dp;
            best = s;
            bs = 1.0;
          }
          if (signed_match && dm < bd)
          {
            bd = dm;
            best = s;
            bs = -1.0;
          }
        }
        perm[g] = best;
        sgn[g] = bs;
        worst = std::max(worst, bd);
      }
      double gworst = 0.0;
      MPI_Allreduce(&worst, &gworst, 1, MPI_DOUBLE, MPI_MAX, comm);
      MFEM_VERIFY(std::sqrt(gworst) < 1e-8,
                  "Online interface DOFs do not match the saved model (max mismatch "
                      << std::sqrt(gworst) << "); the interface Gamma must be identical.");
    }
    if (rank == 0 && file_w > 0)
    {
      const std::vector<double> S_off = S_full;
      for (int i = 0; i < nG; i++)
      {
        for (int j = 0; j < nG; j++)
        {
          S_full[static_cast<std::size_t>(i) * nG + j] =
              sgn[i] * sgn[j] * S_off[static_cast<std::size_t>(perm[i]) * nG + perm[j]];
        }
      }
    }
    MPI_Scatterv(rank == 0 ? S_full.data() : nullptr, row_cnt.data(), row_disp.data(),
                 MPI_DOUBLE, impl->S_rows.data(), nloc * nG, MPI_DOUBLE, 0, comm);
    loaded = true;
    // Optional sections after S_E: the magnetostatic energy operator S^K and the lift
    // modes.
    const std::streamoff sections_pos =
        static_cast<std::streamoff>(2 * sizeof(int)) +
        static_cast<std::streamoff>(sizeof(double)) * nG * file_w +
        static_cast<std::streamoff>(sizeof(double)) * nG * nG;
    impl->LoadSections(model_path, sections_pos, perm, sgn);
    impl->CheckEnvFingerprint();
  }
  // MUMPS Schur materialization (one partial factorization, no |Gamma| back-solves) when
  // available and the environment fits a direct factorization. Not used when the
  // magnetostatic energy operator S^K is also materialized (its correction needs the lifted
  // columns, which the batched back-solve path produces anyway).
  bool mumps_done = false;
#if defined(MFEM_USE_MUMPS)
  if (!loaded && nG > 0 && !(impl->EnergyDiffers() && !model_path.empty()) &&
      impl->EnvDirectSize() <= Impl::kDirectMaxDofs)
  {
    impl->MaterializeMumps();
    mumps_done = true;
  }
#endif
  if (!loaded && !mumps_done && impl->iodata.solver.substructuring->factorization_tol > 0.0)
  {
    Mpi::Warning(comm, "Solver.Substructuring.FactorizationTol is only used with the MUMPS "
                       "environment factorization; the environment is factored exactly.\n");
  }
  if (!loaded)
  {
    impl->EnsureEnv();  // no-op after the MUMPS materialization (it is the env solver)
    // Map each owned interface column (global index in [off, off+nloc)) to its local true
    // DOF, so each unit-vector RHS is set in O(1) instead of scanning all true DOFs per
    // column.
    std::vector<int> col_to_dof(nloc, -1);
    for (int i = 0; i < impl->nt; i++)
    {
      if (impl->is_gamma[i])
      {
        col_to_dof[impl->gamma_global[i] - off] = i;
      }
    }
    // S_E e_c = (A_env e_c)|_Gamma - (A_env A_EE^-1 (A_env e_c)|_E)|_Gamma. The results are
    // complete (assembled) true-DOF vectors, so each rank stores its own S_E rows directly
    // -- no interface Allreduce (O(nG^2) communication) is needed.
    auto store_col = [&](int c, const Vector &y)
    {
      for (int r = 0; r < nloc; r++)
      {
        impl->S_rows[static_cast<std::size_t>(r) * nG + c] = y(col_to_dof[r]);
      }
    };
    auto set_unit = [&](Vector &e, int c)
    {
      e = 0.0;
      if (c >= off && c < off + nloc)
      {
        e(col_to_dof[c - off]) = 1.0;
      }
    };
    if (mumps_done)
    {
      // S_E already materialized by the MUMPS Schur factorization.
    }
#if defined(MFEM_USE_SUPERLU)
    else if (impl->env_parent_direct)
    {
      // Batched direct materialization: kMaterializeBlock columns per multi-RHS back-solve
      // through the (single) environment factor. This is its first solve, which fixes the
      // factor's block size at kMaterializeBlock for the online excitation batches too.
      // Thin couplings built once: G = A_env E_Gamma (column c = A_env e_c, supported next
      // to Gamma) and R = E_Gamma^T A_env (the interface rows), so each column costs O(nnz
      // near Gamma) instead of two full A_env matvecs.
      std::unique_ptr<mfem::HypreParMatrix> E = impl->AssembleInterfaceSelection();
      std::unique_ptr<mfem::HypreParMatrix> G(mfem::ParMult(impl->A_env.get(), E.get()));
      std::unique_ptr<mfem::HypreParMatrix> Et(E->Transpose());
      std::unique_ptr<mfem::HypreParMatrix> R(mfem::ParMult(Et.get(), impl->A_env.get()));
      E.reset();
      Et.reset();
      const int B = Impl::kMaterializeBlock;
      std::vector<Vector> t(B, Vector(impl->nt));
      std::vector<const mfem::Vector *> X(B);
      std::vector<mfem::Vector *> Y(B);
      Vector ec(nloc), z(nloc);
      // Magnetostatic model for reuse: also materialize the energy interface operator
      // S^K = S_E - P^T D P (one more batched solve per column), so an online run needs no
      // environment solve for the inductance energies.
      const bool with_sk = impl->EnergyDiffers() && !model_path.empty();
      std::vector<Vector> v(with_sk ? B : 0, Vector(impl->nt));
      if (with_sk)
      {
        impl->SK_rows.assign(static_cast<std::size_t>(nloc) * nG, 0.0);
        impl->have_sk = true;
      }
      const mfem::HypreParMatrix &Dm = with_sk ? impl->Denv() : *impl->A_env;
      std::vector<double> agg(static_cast<std::size_t>(B) * nloc);  // A_GG e_c, own rows
      for (int c0 = 0; c0 < nG; c0 += B)
      {
        const int nb = std::min(B, nG - c0);
        X.resize(nb);
        Y.resize(nb);
        for (int k = 0; k < nb; k++)
        {
          const int c = c0 + k;
          ec = 0.0;
          if (c >= off && c < off + nloc)
          {
            ec(c - off) = 1.0;
          }
          // A_env e_c: its env-interior part is the solve RHS, its Gamma part is A_GG e_c.
          G->Mult(ec, t[k]);
          for (int r = 0; r < nloc; r++)
          {
            agg[static_cast<std::size_t>(k) * nloc + r] = t[k](col_to_dof[r]);
          }
          X[k] = &t[k];
          Y[k] = &t[k];  // solve in place
        }
        impl->ApplyAeeInvMulti(X, Y);
        for (int k = 0; k < nb; k++)
        {
          R->Mult(t[k], z);  // (A_env y)|_Gamma = A_GE A_EE^-1 A_EG e_c, this rank's rows
          const int c = c0 + k;
          for (int r = 0; r < nloc; r++)
          {
            impl->S_rows[static_cast<std::size_t>(r) * nG + c] =
                agg[static_cast<std::size_t>(k) * nloc + r] - z(r);
          }
        }
        if (with_sk)
        {
          // Energy correction P^T D P e_c: p = e_c - y (the A-harmonic lift of e_c; t[k]
          // holds y), v = D p, then S^K col = S col - (v|_Gamma - R A_EE^-1 v|_E).
          for (int k = 0; k < nb; k++)
          {
            const int c = c0 + k;
            t[k].Neg();
            if (c >= off && c < off + nloc)
            {
              t[k](col_to_dof[c - off]) += 1.0;
            }
            Dm.Mult(t[k], v[k]);
            for (int r = 0; r < nloc; r++)
            {
              agg[static_cast<std::size_t>(k) * nloc + r] = v[k](col_to_dof[r]);
            }
            X[k] = &v[k];
            Y[k] = &v[k];
          }
          impl->ApplyAeeInvMulti(X, Y);
          for (int k = 0; k < nb; k++)
          {
            R->Mult(v[k], z);
            const int c = c0 + k;
            for (int r = 0; r < nloc; r++)
            {
              const std::size_t q = static_cast<std::size_t>(r) * nG + c;
              impl->SK_rows[q] =
                  impl->S_rows[q] - (agg[static_cast<std::size_t>(k) * nloc + r] - z(r));
            }
          }
        }
      }
    }
    else
#endif
    {
      Vector e(impl->nt), y(impl->nt);
      for (int c = 0; c < nG; c++)
      {
        set_unit(e, c);
        impl->dtn->Mult(e, y);
        store_col(c, y);
      }
    }
    if (!impl->magnetostatic)
    {
      // Terminal modes: K batched env solves (cheap vs. the |Gamma| columns). Magnetostatic
      // lifts (flux loops) are only known to the driver; their modes are computed (and
      // appended to the model) at the first EnergyMatrix call.
      std::vector<int> ids;
      std::vector<Vector> lifts;
      for (const auto &[idx, dofs] : impl->terminal_tdofs)
      {
        ids.push_back(idx);
        lifts.push_back(impl->TerminalMode(idx));
      }
      impl->ComputeModes(ids, lifts);
    }
    if (!model_path.empty())
    {
      const int sig_type = impl->magnetostatic ? 2 : 1;
      const std::vector<double> sig = gamma_sig();  // collective: all ranks participate
      // Gather the distributed row blocks to rank 0 to write the full matrix to disk.
      std::vector<double> S_full;
      if (rank == 0)
      {
        S_full.assign(static_cast<std::size_t>(nG) * nG, 0.0);
      }
      MPI_Gatherv(impl->S_rows.data(), nloc * nG, MPI_DOUBLE,
                  rank == 0 ? S_full.data() : nullptr, row_cnt.data(), row_disp.data(),
                  MPI_DOUBLE, 0, comm);
      if (rank == 0)
      {
        std::ofstream f(model_path, std::ios::binary);
        f.write(reinterpret_cast<const char *>(&nG), sizeof(int));
        f.write(reinterpret_cast<const char *>(&sig_type), sizeof(int));
        f.write(reinterpret_cast<const char *>(sig.data()), sizeof(double) * nG * sig_w);
        f.write(reinterpret_cast<const char *>(S_full.data()), sizeof(double) * nG * nG);
      }
      impl->AppendEnvFingerprint(model_path);  // collective
      if (impl->have_sk)
      {
        impl->AppendSK(model_path);  // collective
      }
      if (impl->modes_ready)
      {
        impl->AppendModes(model_path);  // collective
      }
    }
  }
  // Optional hierarchical (HODLR) off-diagonal compression of the assembled S_E: the DtN's
  // well-separated interface-block couplings are low-rank, so this compresses S_E storage +
  // apply to a controlled relative tolerance while the near-field / diagonal stays exact
  // (unlike a global low-rank, which fails -- S_E is full rank). Built on rank 0 from the
  // gathered dense S_E, then broadcast (compressed) and applied replicated on every rank.
  {
    const double hodlr_tol = impl->iodata.solver.substructuring->interface_offdiag_tol;
    if (hodlr_tol > 0.0 && nG > 0)
    {
      // Interface DOF coordinates (replicated) for the coordinate-median clustering. Each
      // true DOF is placed at its interpolation point: H1 -> node coordinate, Nedelec ->
      // edge midpoint. Computed per element (GetNodes mapped through the element
      // transformation), reduced onto the owning rank's true DOF, then summed to a
      // replicated array.
      std::vector<double> coords_loc(static_cast<std::size_t>(nG) * 3, 0.0),
          coords(static_cast<std::size_t>(nG) * 3, 0.0);
      std::vector<double> tdof_xyz(static_cast<std::size_t>(impl->nt) * 3, 0.0);
      std::vector<char> have(impl->nt, 0);
      mfem::Array<int> edofs;
      mfem::Vector phys;
      for (int e = 0; e < impl->parent_fes.GetNE(); e++)
      {
        const mfem::FiniteElement *fe = impl->parent_fes.GetFE(e);
        mfem::ElementTransformation *T = impl->parent.GetElementTransformation(e);
        const mfem::IntegrationRule &nodes = fe->GetNodes();
        impl->parent_fes.GetElementDofs(e, edofs);
        for (int j = 0; j < edofs.Size(); j++)
        {
          const int ldof = edofs[j] >= 0 ? edofs[j] : -1 - edofs[j];
          const int t = impl->parent_fes.GetLocalTDofNumber(ldof);
          if (t < 0 || have[t])
          {
            continue;
          }
          T->Transform(nodes.IntPoint(j), phys);
          for (int d = 0; d < phys.Size() && d < 3; d++)
          {
            tdof_xyz[static_cast<std::size_t>(t) * 3 + d] = phys(d);
          }
          have[t] = 1;
        }
      }
      for (int i = 0; i < impl->nt; i++)
      {
        if (impl->is_gamma[i])
        {
          for (int d = 0; d < 3; d++)
          {
            coords_loc[static_cast<std::size_t>(impl->gamma_global[i]) * 3 + d] =
                tdof_xyz[static_cast<std::size_t>(i) * 3 + d];
          }
        }
      }
      MPI_Allreduce(coords_loc.data(), coords.data(), nG * 3, MPI_DOUBLE, MPI_SUM, comm);
      // Gather a dense interface operator to rank 0, build the Hodlr there, broadcast it,
      // and free the dense row blocks (the compressed operator supplants them).
      auto compress = [&](std::vector<double> &rows, const char *name)
      {
        std::vector<double> S_full;
        if (rank == 0)
        {
          S_full.assign(static_cast<std::size_t>(nG) * nG, 0.0);
        }
        MPI_Gatherv(rows.data(), nloc * nG, MPI_DOUBLE, rank == 0 ? S_full.data() : nullptr,
                    row_cnt.data(), row_disp.data(), MPI_DOUBLE, 0, comm);
        std::vector<double> buf;
        int buflen = 0;
        if (rank == 0)
        {
          Hodlr h;
          h.n = nG;
          h.perm.assign(nG, 0);
          std::vector<int> idx(nG);
          std::iota(idx.begin(), idx.end(), 0);
          BuildHodlr(S_full, nG, idx, 0, coords, hodlr_tol, 32, h);
          buf = h.Serialize();
          buflen = static_cast<int>(buf.size());
        }
        MPI_Bcast(&buflen, 1, MPI_INT, 0, comm);
        buf.resize(buflen);
        MPI_Bcast(buf.data(), buflen, MPI_DOUBLE, 0, comm);
        auto h = std::make_unique<Hodlr>(Hodlr::Deserialize(buf));
        const long long ret = h->Storage();
        rows.clear();
        rows.shrink_to_fit();
        Mpi::Print("[HODLR] tol={:.1e}: {} storage {:d}/{:d} = {:.3f} of dense\n",
                   hodlr_tol, name, ret, static_cast<long long>(nG) * nG,
                   static_cast<double>(ret) / (static_cast<double>(nG) * nG));
        return h;
      };
      impl->hodlr = compress(impl->S_rows, "S_E");
      if (impl->have_sk)
      {
        impl->hodlr_K = compress(impl->SK_rows, "S^K");
      }
    }
  }

  // g_E is excitation-dependent; it is computed per excitation in the region solve.
  impl->mat_dtn =
      std::make_unique<MaterializedDtN>(impl->S_rows, impl->gamma_off, impl->gamma_global,
                                        impl->nG_global, comm, impl->hodlr.get());
  if (impl->have_sk)
  {
    impl->mat_dtn_K = std::make_unique<MaterializedDtN>(impl->SK_rows, impl->gamma_off,
                                                        impl->gamma_global, impl->nG_global,
                                                        comm, impl->hodlr_K.get());
  }

  // Region-condensed solver: Palace CG preconditioned by a wrapped AMS (H(curl)) or
  // BoomerAMG (H1) on the region-free block. Built once and reused across excitations (the
  // condensed interface operator S_E is fixed). The region operator carries the same mass
  // regularization as the environment, so it is definite and AMS is used without the
  // singular-problem option.
  {
    std::unique_ptr<Solver<Operator>> pc;
    bool use_direct = false;
#if defined(MFEM_USE_SUPERLU)
    // Default to a direct factorization of A_region_free when the region is small enough to
    // factor (a near-exact preconditioner -> the outer CG converges in a few iterations);
    // fall back to GMG / AMG for very large regions.
    {
      long long reg_loc = 0;
      for (int i = 0; i < impl->nt; i++)
      {
        if (impl->is_region_free[i])
        {
          reg_loc++;
        }
      }
      long long reg_glob = 0;
      MPI_Allreduce(&reg_loc, &reg_glob, 1, MPI_LONG_LONG, MPI_SUM, comm);
      use_direct = (reg_glob <= Impl::kDirectMaxDofs);
    }
#endif
#if defined(MFEM_USE_SUPERLU)
    if (use_direct)
    {
      // Direct factorization of the region preconditioner. When the interface is small
      // enough for a dense Gamma block (and S_E is dense, not HODLR), factor the EXACT
      // condensed operator A_region_free + S_E, so the outer CG converges in one iteration:
      // the non-local DtN cannot be sparsified (truncating it to the interface stencil
      // makes the preconditioner far worse than omitting it), and A_region_free alone needs
      // O(50) iterations. Otherwise fall back to A_region_free (S_E omitted).
      // The factor is created at the first solve batch, once the number of excitations is
      // known (EnsureRegionFactor chooses the exact dense-Gamma factor or A_region_free).
      impl->region_direct = true;
      impl->region_dense_eligible =
          !impl->hodlr && nG > 0 && nG <= Impl::kDirectCondensedMaxInterface;
      pc = std::make_unique<CallableSolver>(impl->A_region_free->Height(),
                                            [pi](const mfem::Vector &r, mfem::Vector &z)
                                            { pi->reg_lu->Mult(r, z); });
    }
    else
#endif
        if (impl->BuildRegionGmg())
    {
      // Order>=2 H1: geometric multigrid on the region submesh, routed through the region
      // preconditioner apply.
      pc = std::make_unique<CallableSolver>(impl->A_region_free->Height(),
                                            [pi](const mfem::Vector &r, mfem::Vector &z)
                                            { pi->ApplyRegionGmgPc(r, z); });
    }
    else if (impl->magnetostatic)
    {
      auto ams = std::make_unique<mfem::HypreAMS>(&impl->parent_fes);
      ams->SetPrintLevel(0);
      pc = std::make_unique<MfemWrapperSolver<Operator>>(
          std::move(ams), /*save_assembled=*/true, /*complex_matrix=*/false,
          /*drop_small_entries=*/false);
    }
    else
    {
      pc = std::make_unique<MfemWrapperSolver<Operator>>(
          std::make_unique<BoomerAmgSolver>(1, 1, true, 0), /*save_assembled=*/true,
          /*complex_matrix=*/false, /*drop_small_entries=*/false);
    }
    auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
    pcg->SetInitialGuess(false);
    pcg->SetRelTol(1.0e-10);
    pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
    pcg->SetMaxIter(2000);
    impl->region_op = std::make_unique<RegionCondensedOperator>(
        *impl->A_region_free, *impl->mat_dtn, impl->is_region_free);
    impl->region_ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(pc));
    impl->region_ksp->SetOperators(*impl->region_op, *impl->A_region_free);
  }
}

namespace
{

// Terminal excitation Dirichlet field: the driven terminal at 1 V, all others grounded.
Vector TerminalDbc(const std::map<int, std::vector<int>> &terminal_tdofs, int drive_idx,
                   int nt)
{
  Vector dbc(nt);
  dbc = 0.0;
  for (const auto &[idx, dofs] : terminal_tdofs)
  {
    const double value = (idx == drive_idx) ? 1.0 : 0.0;
    for (int d : dofs)
    {
      dbc(d) = value;
    }
  }
  return dbc;
}

}  // namespace

Vector SubstructuringSolver::SolveExcitation(int drive_idx)
{
  return SolveExcitations({drive_idx})[0];
}

std::vector<Vector>
SubstructuringSolver::SolveExcitations(const std::vector<int> &drive_terminal_indices)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  std::vector<Vector> dbcs;
  dbcs.reserve(drive_terminal_indices.size());
  for (int idx : drive_terminal_indices)
  {
    dbcs.push_back(TerminalDbc(impl->terminal_tdofs, idx, impl->nt));
  }
  return SolveDirichletBatch(dbcs);
}

Vector SubstructuringSolver::SolveDirichlet(const Vector &dbc_values)
{
  return SolveDirichlets({dbc_values})[0];
}

std::vector<Vector>
SubstructuringSolver::SolveDirichlets(const std::vector<Vector> &dbc_values)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  // Prescribe an arbitrary boundary field on the Dirichlet DOF set (e.g. a flux-loop lift).
  std::vector<Vector> dbcs(dbc_values.size(), Vector(impl->nt));
  for (std::size_t k = 0; k < dbc_values.size(); k++)
  {
    dbcs[k] = 0.0;
    for (int i = 0; i < impl->dbc_tdofs.Size(); i++)
    {
      const int d = impl->dbc_tdofs[i];
      dbcs[k](d) = dbc_values[k](d);
    }
  }
  return SolveDirichletBatch(dbcs);
}

std::vector<Vector>
SubstructuringSolver::SolveDirichletBatch(const std::vector<Vector> &dbcs)
{
  impl->EnsureEnv();  // general Dirichlet data needs environment solves
  const int nt = impl->nt;
  const int n = static_cast<int>(dbcs.size());
  impl->EnsureRegionFactor(n);
  MPI_Comm comm = impl->parent_fes.GetComm();

  // g_E per excitation: the interface response to the Dirichlet data,
  //   g_E = (A_env x)|_Gamma - (A_env A_EE^-1 (A_env x)|_E)|_Gamma.
  // The environment solve is skipped (exactly) when (A_env x)|_E vanishes, i.e. the
  // Dirichlet data does not touch the environment (e.g. a region terminal); the others run
  // as one batched multi-RHS solve.
  std::vector<Vector> t(n, Vector(nt)), gE(n, Vector(nt));
  std::vector<int> touches(n, 0);
  for (int k = 0; k < n; k++)
  {
    impl->A_env->Mult(dbcs[k], t[k]);
    for (int i = 0; i < nt && !touches[k]; i++)
    {
      if (impl->is_env_int[i] && t[k](i) != 0.0)
      {
        touches[k] = 1;
      }
    }
  }
  if (n > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, touches.data(), n, MPI_INT, MPI_MAX, comm);
  }
  {
    std::vector<const mfem::Vector *> X;
    std::vector<mfem::Vector *> Y;
    std::vector<int> ks;
    for (int k = 0; k < n; k++)
    {
      if (touches[k])
      {
        X.push_back(&t[k]);
        Y.push_back(&gE[k]);
        ks.push_back(k);
      }
    }
    if (!X.empty())
    {
      impl->ApplyAeeInvMulti(X, Y);
    }
    Vector t2(nt);
    for (int k = 0; k < n; k++)
    {
      if (touches[k])
      {
        impl->A_env->Mult(gE[k], t2);
        for (int i = 0; i < nt; i++)
        {
          gE[k](i) = t[k](i) - t2(i);
        }
      }
      else
      {
        gE[k] = t[k];
      }
    }
  }

  // Region-condensed solve per excitation using the materialized S_E (no environment solves
  // in the loop). RHS: region Dirichlet elimination + environment load g_E on the
  // interface. Both are complete (assembled) true-DOF vectors, so each rank reads its own
  // entries.
  std::vector<Vector> u(n, Vector(nt));
  {
    Vector b(nt), tr(nt);
    for (int k = 0; k < n; k++)
    {
      impl->A_region->Mult(dbcs[k], tr);
      b = 0.0;
      for (int i = 0; i < nt; i++)
      {
        if (impl->is_region_free[i])
        {
          b(i) -= tr(i);
        }
        if (impl->is_gamma[i])
        {
          b(i) -= gE[k](i);
        }
      }
      u[k] = 0.0;
      impl->region_ksp->Mult(b, u[k]);
      for (int i = 0; i < impl->dbc_tdofs.Size(); i++)
      {
        u[k](impl->dbc_tdofs[i]) = dbcs[k](impl->dbc_tdofs[i]);
      }
    }
  }

  impl->RecoverEnvInterior(u);  // batched: u_E = -A_EE^-1 (A_env u)|_E
  return u;
}

Vector SubstructuringSolver::SolveRegion()
{
  // Default single excitation: drive the lowest-index terminal.
  MFEM_VERIFY(!impl->terminal_tdofs.empty(), "No terminals configured!");
  return SolveExcitation(impl->terminal_tdofs.begin()->first);
}

Vector SubstructuringSolver::SolveSource(const Vector &f)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  impl->EnsureEnv();
  impl->EnsureRegionFactor(1);
  const int nt = impl->nt;

  // Environment interior source response: solve A_EE w = f_E, correction (A_env w)|_Gamma.
  Vector fE(nt), w(nt), Aw(nt);
  fE = 0.0;
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_env_int[i])
    {
      fE(i) = f(i);
    }
  }
  w = 0.0;
  impl->ApplyAeeInv(fE, w);
  impl->A_env->Mult(w, Aw);  // assembled: each rank reads its own interface entries

  // RHS: region/interface source minus the environment source correction on the interface.
  Vector b(nt);
  b = 0.0;
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_region_free[i])
    {
      b(i) = f(i);
    }
  }
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_gamma[i])
    {
      b(i) -= Aw(i);
    }
  }

  Vector u(nt);
  u = 0.0;
  impl->region_ksp->Mult(b, u);

  // Recover environment interior: u_E = A_EE^-1 (f_E - (A_env u)|_E).
  Vector Au(nt), rhs(nt), uE(nt);
  impl->A_env->Mult(u, Au);
  rhs = 0.0;
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_env_int[i])
    {
      rhs(i) = fE(i) - Au(i);
    }
  }
  uE = 0.0;
  impl->ApplyAeeInv(rhs, uE);
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_env_int[i])
    {
      u(i) = uE(i);
    }
  }
  return u;
}

mfem::DenseMatrix SubstructuringSolver::SheetEnergyMatrix(const std::vector<int> &ids,
                                                          const std::vector<Vector> &a,
                                                          std::vector<Vector> *fields,
                                                          int n_fields)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  MFEM_VERIFY(impl->magnetostatic && impl->has_sheets,
              "SheetEnergyMatrix needs London superconductor sheets!");
  MFEM_VERIFY(ids.size() == a.size(), "SheetEnergyMatrix: one id per generator!");
  const int nt = impl->nt;
  const int n = static_cast<int>(a.size());
  const int nf = fields ? std::max(0, std::min(n_fields, n)) : 0;
  MPI_Comm comm = impl->parent_fes.GetComm();

  // Energies of the London flux states, as the native solver:
  //   E_ij = u_i^T K_cc u_j + (u_i - a_i)^T M_sheet (u_j - a_j),
  // with the kinetic term formed from the differences (u - a is small for a stiff sheet,
  // and the expanded form would cancel).
  auto kinetic = [&](const Vector &ui, const Vector &uj, const Vector &ai, const Vector &aj,
                     const mfem::HypreParMatrix &Ms)
  {
    Vector di(ui), dj(uj), Mdj(nt);
    di -= ai;
    dj -= aj;
    Ms.Mult(dj, Mdj);
    return di * Mdj;
  };

  // Without the energy interface operator S^K (no model saved or loaded), or when fields
  // are requested: full solves with the environment.
  if (!impl->mat_dtn_K || nf > 0)
  {
    std::vector<Vector> u(n, Vector(nt));
    Vector b(nt);
    for (int k = 0; k < n; k++)
    {
      impl->M_sheet->Mult(a[k], b);
      for (int d = 0; d < impl->dbc_tdofs.Size(); d++)
      {
        b(impl->dbc_tdofs[d]) = 0.0;
      }
      u[k] = SolveSource(b);
    }
    mfem::DenseMatrix E(n);
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        double loc = kinetic(u[i], u[j], a[i], a[j], *impl->M_sheet), glob = 0.0;
        MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE, MPI_SUM, comm);
        E(i, j) = MutualEnergy(u[i], u[j]) + glob;
      }
    }
    if (fields)
    {
      fields->assign(u.begin(), u.begin() + nf);
    }
    return E;
  }

  // Source modes of these generators: reuse the saved / current ones if they match (id +
  // environment source fingerprint), else compute them (needs the environment once;
  // appended to the model in an offline run that saves one).
  std::vector<int> col;
  if (!impl->SourceModesMatch(ids, a, col))
  {
    impl->ComputeSourceModes(ids, a);
    col.resize(n);
    std::iota(col.begin(), col.end(), 0);
    const auto &subcfg = *impl->iodata.solver.substructuring;
    if (subcfg.mode != SubstructuringMode::ONLINE && !subcfg.save_model.empty())
    {
      impl->AppendSourceModes(subcfg.save_model);
    }
  }
  const int K = static_cast<int>(impl->src_ids.size());
  auto row = [&](int i)
  { return static_cast<std::size_t>(impl->gamma_global[i] - impl->gamma_off) * K; };

  // Region-condensed solve per state (environment interior left at zero): the region source
  // M_sheet^R a_k plus the environment's interface load g_k.
  impl->EnsureRegionFactor(n);
  std::vector<Vector> u(n, Vector(nt)), Ku(n, Vector(nt)), Su(n, Vector(nt));
  {
    Vector b(nt), bR(nt);
    for (int j = 0; j < n; j++)
    {
      impl->M_sheet_region->Mult(a[j], bR);
      b = 0.0;
      for (int i = 0; i < nt; i++)
      {
        if (impl->is_region_free[i])
        {
          b(i) = bR(i);
        }
        if (impl->is_gamma[i])
        {
          b(i) += impl->SG_rows[row(i) + col[j]];
        }
      }
      u[j] = 0.0;
      impl->region_ksp->Mult(b, u[j]);
      for (int d = 0; d < impl->dbc_tdofs.Size(); d++)
      {
        u[j](impl->dbc_tdofs[d]) = 0.0;
      }
      impl->K_region_e->Mult(u[j], Ku[j]);
      impl->mat_dtn_K->Mult(u[j], Su[j]);
    }
  }

  // E_ij = [u_i^T K_R u_j + (u_i - a_i)^T M_sheet^R (u_j - a_j)]  (region)
  //      + u_i,G^T S^K u_j,G + u_i,G^T h_j + h_i^T u_j,G + c_ij  (environment).
  mfem::DenseMatrix E(n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      double e = (u[i] * Ku[j]) + kinetic(u[i], u[j], a[i], a[j], *impl->M_sheet_region);
      for (int r = 0; r < nt; r++)
      {
        if (impl->is_gamma[r])
        {
          e += u[i](r) * Su[j](r) + u[i](r) * impl->SH_rows[row(r) + col[j]] +
               impl->SH_rows[row(r) + col[i]] * u[j](r);
        }
      }
      E(i, j) = e;
    }
  }
  if (n > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, E.GetData(), n * n, MPI_DOUBLE, MPI_SUM, comm);
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      E(i, j) += impl->SC[static_cast<std::size_t>(col[i]) * K + col[j]];
    }
  }
  return E;
}

mfem::DenseMatrix SubstructuringSolver::EnergyMatrix(const std::vector<int> &ids,
                                                     const std::vector<Vector> &lifts,
                                                     std::vector<Vector> *fields,
                                                     int n_fields,
                                                     std::vector<Vector> *region_fields)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  MFEM_VERIFY(ids.size() == lifts.size(), "EnergyMatrix: one id per lift!");
  const int nt = impl->nt;
  const int n = static_cast<int>(lifts.size());
  const int nf = fields ? std::max(0, std::min(n_fields, n)) : 0;
  MPI_Comm comm = impl->parent_fes.GetComm();
  // Dirichlet data restricted to the Dirichlet DOF set (as SolveDirichlets does).
  std::vector<Vector> x(n, Vector(nt));
  for (int k = 0; k < n; k++)
  {
    x[k] = 0.0;
    for (int d = 0; d < impl->dbc_tdofs.Size(); d++)
    {
      const int i = impl->dbc_tdofs[d];
      x[k](i) = lifts[k](i);
    }
  }

  // Magnetostatics without the energy interface operator S^K: environment path (full
  // fields, then the energy u_i^T K u_j).
  if (impl->EnergyDiffers() && !impl->mat_dtn_K)
  {
    std::vector<Vector> u = SolveDirichletBatch(x);
    mfem::DenseMatrix E(n);
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        E(i, j) = MutualEnergy(u[i], u[j]);
      }
    }
    if (fields)
    {
      fields->assign(u.begin(), u.begin() + nf);
    }
    if (region_fields)
    {
      *region_fields = u;
    }
    return E;
  }

  // Modes of these lifts: reuse the saved / current ones if they match (id + environment
  // fingerprint), else compute them (needs the environment once; appended to the model in
  // an offline run that saves one).
  std::vector<int> col;
  if (!impl->ModesMatch(ids, x, col))
  {
    impl->ComputeModes(ids, x);
    col.resize(n);
    std::iota(col.begin(), col.end(), 0);
    const auto &subcfg = *impl->iodata.solver.substructuring;
    if (subcfg.mode != SubstructuringMode::ONLINE && !subcfg.save_model.empty())
    {
      impl->AppendModes(subcfg.save_model);
    }
  }
  const int K = static_cast<int>(impl->mode_ids.size());
  const std::vector<double> &GE = impl->have_gk ? impl->GK_rows : impl->G_rows;
  const MaterializedDtN &SK = impl->mat_dtn_K ? *impl->mat_dtn_K : *impl->mat_dtn;
  auto row = [&](int i)
  { return static_cast<std::size_t>(impl->gamma_global[i] - impl->gamma_off) * K; };

  // Region-condensed solve per lift (environment interior left at zero): the interface load
  // is the precomputed g_k, so no environment solve is needed.
  impl->EnsureRegionFactor(n);
  std::vector<Vector> u(n, Vector(nt)), Ku(n, Vector(nt)), Su(n, Vector(nt));
  {
    Vector b(nt), tr(nt);
    for (int j = 0; j < n; j++)
    {
      impl->A_region->Mult(x[j], tr);
      b = 0.0;
      for (int i = 0; i < nt; i++)
      {
        if (impl->is_region_free[i])
        {
          b(i) -= tr(i);
        }
        if (impl->is_gamma[i])
        {
          b(i) -= impl->G_rows[row(i) + col[j]];
        }
      }
      u[j] = 0.0;
      impl->region_ksp->Mult(b, u[j]);
      for (int d = 0; d < impl->dbc_tdofs.Size(); d++)
      {
        const int i = impl->dbc_tdofs[d];
        u[j](i) = x[j](i);
      }
      impl->K_region_e->Mult(u[j], Ku[j]);
      SK.Mult(u[j], Su[j]);
    }
  }

  // E_ij = u_i^T K_R u_j + [u_G,i; e_i]^T [S^K, G^K; G^K^T, Cmode] [u_G,j; e_j]: the region
  // energy plus the environment's exact condensed energy.
  mfem::DenseMatrix E(n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      double a = u[i] * Ku[j];
      for (int r = 0; r < nt; r++)
      {
        if (impl->is_gamma[r])
        {
          a += u[i](r) * Su[j](r) + u[i](r) * GE[row(r) + col[j]] +
               GE[row(r) + col[i]] * u[j](r);
        }
      }
      E(i, j) = a;
    }
  }
  if (n > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, E.GetData(), n * n, MPI_DOUBLE, MPI_SUM, comm);
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      E(i, j) += impl->Cmode[static_cast<std::size_t>(col[i]) * K + col[j]];
    }
  }

  if (region_fields)
  {
    *region_fields = u;
  }
  // Optional full fields (environment interior recovered on demand).
  if (fields)
  {
    fields->assign(u.begin(), u.begin() + nf);
    if (nf > 0)
    {
      impl->RecoverEnvInterior(*fields);
    }
  }
  return E;
}

mfem::DenseMatrix
SubstructuringSolver::CapacitanceMatrix(const std::vector<int> &terminal_indices,
                                        std::vector<Vector> *fields, int n_fields,
                                        std::vector<Vector> *region_fields)
{
  MFEM_VERIFY(!impl->magnetostatic,
              "CapacitanceMatrix is for electrostatic substructuring problems!");
  std::vector<Vector> lifts;
  for (int idx : terminal_indices)
  {
    lifts.push_back(TerminalLift(idx));
  }
  return EnergyMatrix(terminal_indices, lifts, fields, n_fields, region_fields);
}

ErrorIndicator
SubstructuringSolver::RegionErrorIndicator(const std::vector<Vector> &region_fields,
                                           const mfem::DenseMatrix &E) const
{
  MFEM_VERIFY(!impl->magnetostatic,
              "Region error indicators are only available for electrostatics!");
  MFEM_VERIFY(E.Height() == static_cast<int>(region_fields.size()),
              "RegionErrorIndicator: one energy per field!");
  const auto &iodata = impl->iodata;
  const int order = iodata.solver.order, dim = impl->parent.Dimension();

  // Palace's gradient-flux estimator on the region submesh: its flux recovery sees only the
  // region, so the unknown environment field does not pollute the indicators near the
  // interface.
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromDomain(impl->parent, impl->ra_arr))));
  mesh[0]->RebuildCeedAttributes();
  auto &submesh = static_cast<mfem::ParSubMesh &>(mesh[0]->Get());
  MaterialOperator mat_op(iodata, *mesh[0]);
  mfem::H1_FECollection h1_fec(order, dim);
  mfem::ND_FECollection nd_fec(order, dim);
  FiniteElementSpace h1_fespace(*mesh[0], &h1_fec), nd_fespace(*mesh[0], &nd_fec);
  auto rt_fecs = fem::ConstructFECollections<mfem::RT_FECollection>(
      order - 1, dim, 1, iodata.solver.linear.mg_coarsening, false);
  auto rt_fespaces =
      fem::ConstructFiniteElementSpaceHierarchy<mfem::RT_FECollection>(1, mesh, rt_fecs);
  GradFluxErrorEstimator<Vector> estimator(mat_op, nd_fespace, rt_fespaces,
                                           iodata.solver.linear.estimator_tol,
                                           iodata.solver.linear.estimator_max_it, 0, false);
  const auto &grad = nd_fespace.GetDiscreteInterpolator(h1_fespace);

  ErrorIndicator region;
  mfem::ParGridFunction pgf(&impl->parent_fes), sgf(&h1_fespace.Get());
  Vector v(h1_fespace.GetTrueVSize()), e(nd_fespace.GetTrueVSize());
  for (std::size_t k = 0; k < region_fields.size(); k++)
  {
    pgf.SetFromTrueDofs(region_fields[k]);
    sgf = 0.0;
    submesh.Transfer(pgf, sgf);
    sgf.GetTrueDofs(v);
    e = 0.0;
    grad.AddMult(v, e, -1.0);  // E = -grad(V)
    const double energy = 0.5 * E(static_cast<int>(k), static_cast<int>(k));
    estimator.AddErrorIndicator(e, energy, region);
  }

  // Region element indicators onto the parent mesh (zero on the environment).
  Vector local(impl->parent.GetNE());
  local = 0.0;
  if (region.Local().Size() > 0)
  {
    const auto &parent_id = submesh.GetParentElementIDMap();
    const double *r = region.Local().HostRead();
    for (int i = 0; i < parent_id.Size(); i++)
    {
      local(parent_id[i]) = r[i];
    }
  }
  return ErrorIndicator(std::move(local));
}

Vector SubstructuringSolver::TerminalLift(int terminal_index) const
{
  MFEM_VERIFY(impl->terminal_tdofs.count(terminal_index),
              "Unknown terminal index " << terminal_index << "!");
  return impl->TerminalMode(terminal_index);
}

bool SubstructuringSolver::HasSheets() const
{
  return impl->has_sheets;
}

bool SubstructuringSolver::EnvironmentFactored() const
{
  return impl->env_built;
}

std::vector<int> SubstructuringSolver::TerminalIndices() const
{
  std::vector<int> idx;
  for (const auto &[i, dofs] : impl->terminal_tdofs)
  {
    idx.push_back(i);
  }
  return idx;
}

double SubstructuringSolver::MutualEnergy(const Vector &ui, const Vector &uj) const
{
  Vector t(impl->nt), t2(impl->nt);
  impl->K_region_e->Mult(uj, t);
  impl->K_env_e->Mult(uj, t2);
  t += t2;
  double local = 0.0;
  for (int i = 0; i < impl->nt; i++)
  {
    local += ui(i) * t(i);
  }
  double global = 0.0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, impl->parent_fes.GetComm());
  return global;
}

long long int SubstructuringSolver::RegionGlobalTrueVSize() const
{
  return impl->parent_fes.GlobalTrueVSize();
}

void SubstructuringSolver::WriteParaView(const std::string &dir,
                                         const std::vector<int> &terminals,
                                         const std::vector<Vector> &fields) const
{
  MFEM_VERIFY(terminals.size() == fields.size(),
              "WriteParaView requires one field per terminal!");
  const int order = impl->iodata.solver.order;
  mfem::ParGridFunction phi(&impl->parent_fes);
  mfem::ParaViewDataCollection pv("paraview", &impl->parent);
  pv.SetPrefixPath(dir);
  pv.SetLevelsOfDetail(order);
  pv.SetHighOrderOutput(order > 1);
  pv.SetDataFormat(mfem::VTKFormat::BINARY);
  pv.RegisterField("V", &phi);
  for (std::size_t j = 0; j < fields.size(); j++)
  {
    phi.SetFromTrueDofs(fields[j]);
    pv.SetCycle(static_cast<int>(j));
    pv.SetTime(static_cast<double>(terminals[j]));
    pv.Save();
  }
}

double SubstructuringSolver::ElectrostaticEnergy(const Vector &u) const
{
  Vector t(impl->nt), t2(impl->nt);
  impl->K_region_e->Mult(u, t);
  impl->K_env_e->Mult(u, t2);
  t += t2;
  double local = 0.0;
  for (int i = 0; i < impl->nt; i++)
  {
    local += u(i) * t(i);
  }
  double global = 0.0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, impl->parent_fes.GetComm());
  return 0.5 * global;
}

}  // namespace palace
