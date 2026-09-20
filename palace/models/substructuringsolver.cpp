// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolver.hpp"

#include <fstream>
#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/mesh.hpp"
#include "fem/substructure.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

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
// using the parent-space environment matrix A_env (with A_EE the environment-interior
// block, realized as A_env with all non-interior true DOFs identity-eliminated) and its
// solver. All operations are distributed matvecs plus one distributed A_EE solve, so this
// is parallel by construction.
class ImplicitDtN : public mfem::Operator
{
public:
  ImplicitDtN(mfem::HypreParMatrix &A_env, mfem::Solver &Aee_inv,
              const std::vector<char> &is_gamma, const std::vector<char> &is_env_int)
    : mfem::Operator(A_env.Height()), A_env(A_env), Aee_inv(Aee_inv), is_gamma(is_gamma),
      is_env_int(is_env_int), t(A_env.Height()), rhs(A_env.Height()), ye(A_env.Height()),
      t2(A_env.Height())
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
    Aee_inv.Mult(rhs, ye);
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
  mfem::Solver &Aee_inv;
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

// Materialized DtN: applies a replicated dense interface operator S_E (computed once) to a
// distributed interface vector via a gather (Allreduce over the global interface
// enumeration) and a local dense apply on owned interface rows. Cheap and reusable across
// region solves.
class MaterializedDtN : public mfem::Operator
{
public:
  MaterializedDtN(const mfem::DenseMatrix &S, const std::vector<int> &gamma_global,
                  int nG_global, MPI_Comm comm)
    : mfem::Operator(static_cast<int>(gamma_global.size())), S(S),
      gamma_global(gamma_global), nG_global(nG_global), comm(comm)
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
    for (int i = 0; i < height; i++)
    {
      if (gamma_global[i] >= 0)
      {
        double s = 0.0;
        for (int j = 0; j < nG_global; j++)
        {
          s += S(gamma_global[i], j) * xg[j];
        }
        y(i) = s;
      }
    }
  }

private:
  const mfem::DenseMatrix &S;
  const std::vector<int> &gamma_global;
  int nG_global;
  MPI_Comm comm;
};

}  // namespace

// Parallel region-condensed static solve (electrostatic H1 or magnetostatic H(curl)).
// Region/environment operators are assembled on the parent finite element space with
// domain- restricted material coefficients; the interface is identified in true-DOF space;
// the environment is condensed through an implicit distributed DtN. The magnetostatic
// curl-curl operator carries a subdomain mass regularization (gauge-free formulation is a
// follow-up).
struct SubstructuringSolver::Impl
{
  const IoData &iodata;
  mfem::ParMesh &parent;
  bool magnetostatic;
  std::unique_ptr<mfem::FiniteElementCollection> fec;
  mfem::ParFiniteElementSpace parent_fes;
  int nt;

  std::unique_ptr<mfem::HypreParMatrix> A_region, A_env, A_env_int, A_region_free;
  // Energy (QoI) operators. Electrostatic: alias the solve operators. Magnetostatic: pure
  // curl-curl (no mass), so the magnetic energy / inductance is physical.
  std::unique_ptr<mfem::HypreParMatrix> A_region_energy, A_env_energy;
  mfem::HypreParMatrix *K_region_e = nullptr, *K_env_e = nullptr;
  std::vector<char> is_gamma, is_env_int, is_region_free;
  mfem::Array<int> dbc_tdofs;
  mfem::Vector dbc_values;  // full parent true-DOF vector, prescribed values on Dirichlet
  std::map<int, std::vector<int>> terminal_tdofs;  // terminal index -> its true DOFs

  std::unique_ptr<mfem::HypreSolver> prec_env;
  std::unique_ptr<mfem::HypreSolver> prec_region;
  std::unique_ptr<mfem::HyprePCG> solver_env;
  std::unique_ptr<ImplicitDtN> dtn;

  // Materialized (reusable) interface operator: replicated dense S_E + load g_E over a
  // global interface enumeration, computed once so region solves need no environment
  // solves.
  std::vector<int> gamma_global;  // owned parent true DOF -> global interface index, or -1
  int nG_global = 0;
  mfem::DenseMatrix S_dense;
  std::unique_ptr<MaterializedDtN> mat_dtn;

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

    // Interface / region / environment true-DOF markers.
    mfem::Array<int> rm, em, im;
    MarkInterfaceTrueDofs(parent_fes, ra, ea, rm, em, im);

    // Terminal Dirichlet true DOFs (per terminal). The Dirichlet DOF *set* is the union of
    // all terminals and is fixed across excitations; only the prescribed values change, so
    // the interface/interior partition below is excitation-independent.
    const auto &terminals = iodata.boundaries.terminal;
    mfem::Array<int> dir_mark(nt);
    dir_mark = 0;
    dbc_values.SetSize(nt);
    dbc_values = 0.0;
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
    A_region = AssembleParent(region_eps, magnetostatic);
    A_env = AssembleParent(env_eps, magnetostatic);
    if (magnetostatic)
    {
      A_region_energy = AssembleParent(region_eps, false);  // pure curl-curl
      A_env_energy = AssembleParent(env_eps, false);
      K_region_e = A_region_energy.get();
      K_env_e = A_env_energy.get();
    }
    else
    {
      K_region_e = A_region.get();
      K_env_e = A_env.get();
    }
  }

  // Small mass regularization making the magnetostatic curl-curl solve operator positive
  // definite (so AMS converges without the singular-problem option). Kept small so the
  // magnetic energy, measured with the pure curl-curl operators, is physical to ~1e-2 %.
  // Requires a divergence-free excitation (physical current). Exact (gauge-free, pseudo-
  // inverse DtN) region solves via a submesh AMS are a follow-up refinement.
  static constexpr double kMagRegularization = 1.0e-3;

  std::unique_ptr<mfem::HypreParMatrix>
  AssembleParent(const std::map<int, mfem::DenseMatrix> &coef_by_attr, bool with_mass)
  {
    PWMatrixCoefficient coef(parent.Dimension(), coef_by_attr);
    mfem::ParBilinearForm a(&parent_fes);
    if (magnetostatic)
    {
      a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(coef));
      const int max_attr = parent.attributes.Size() ? parent.attributes.Max() : 1;
      mfem::Vector mass(max_attr);
      mass = 0.0;
      if (with_mass)
      {
        for (const auto &[a_attr, m] : coef_by_attr)
        {
          mass(a_attr - 1) = kMagRegularization;
        }
      }
      mfem::PWConstCoefficient mcoef(mass);  // outlives Assemble() below
      if (with_mass)
      {
        a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(mcoef));
      }
      a.Assemble();
      a.Finalize();
      return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
    }
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(coef));
    a.Assemble();
    a.Finalize();
    return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
  }

  std::map<int, mfem::DenseMatrix> region_eps, env_eps;
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
  // A_EE: the environment operator with all non-(environment-interior) true DOFs
  // identity-eliminated, so its inverse acts only on the environment interior.
  impl->A_env_int = std::make_unique<mfem::HypreParMatrix>(*impl->A_env);
  mfem::Array<int> non_int;
  for (int i = 0; i < impl->nt; i++)
  {
    if (!impl->is_env_int[i])
    {
      non_int.Append(i);
    }
  }
  {
    std::unique_ptr<mfem::HypreParMatrix> tmp(impl->A_env_int->EliminateRowsCols(non_int));
  }
  if (impl->magnetostatic)
  {
    auto ams = std::make_unique<mfem::HypreAMS>(*impl->A_env_int, &impl->parent_fes);
    ams->SetPrintLevel(0);
    impl->prec_env = std::move(ams);
  }
  else
  {
    auto amg = std::make_unique<mfem::HypreBoomerAMG>(*impl->A_env_int);
    amg->SetPrintLevel(0);
    impl->prec_env = std::move(amg);
  }
  impl->solver_env = std::make_unique<mfem::HyprePCG>(*impl->A_env_int);
  impl->solver_env->SetTol(1.0e-13);
  impl->solver_env->SetMaxIter(1000);
  impl->solver_env->SetPrintLevel(0);
  impl->solver_env->SetPreconditioner(*impl->prec_env);
  impl->dtn = std::make_unique<ImplicitDtN>(*impl->A_env, *impl->solver_env, impl->is_gamma,
                                            impl->is_env_int);

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

  // Region-free preconditioner for the region-condensed CG. For magnetostatics the
  // curl-curl block is singular, so an AMS singular-problem preconditioner is essential for
  // CG to converge on the gradient nullspace.
  if (impl->magnetostatic)
  {
    auto ams = std::make_unique<mfem::HypreAMS>(*impl->A_region_free, &impl->parent_fes);
    ams->SetPrintLevel(0);
    impl->prec_region = std::move(ams);
  }
  else
  {
    auto amg = std::make_unique<mfem::HypreBoomerAMG>(*impl->A_region_free);
    amg->SetPrintLevel(0);
    impl->prec_region = std::move(amg);
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
  auto gather_interface = [&](const Vector &y, double *col)
  {
    std::vector<double> loc(nG, 0.0);
    for (int i = 0; i < impl->nt; i++)
    {
      if (impl->is_gamma[i])
      {
        loc[impl->gamma_global[i]] = y(i);
      }
    }
    MPI_Allreduce(loc.data(), col, nG, MPI_DOUBLE, MPI_SUM, comm);
  };
  impl->S_dense.SetSize(nG);
  impl->S_dense = 0.0;

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
    MFEM_VERIFY(nG_file == nG, "Saved substructuring model interface size ("
                                   << nG_file << ") does not match this run (" << nG
                                   << "); the mesh and partition count must be identical.");
    if (rank == 0)
    {
      std::ifstream f(model_path, std::ios::binary);
      f.seekg(sizeof(int));
      f.read(reinterpret_cast<char *>(impl->S_dense.GetData()), sizeof(double) * nG * nG);
    }
    MPI_Bcast(impl->S_dense.GetData(), nG * nG, MPI_DOUBLE, 0, comm);
    loaded = true;
  }
  if (!loaded)
  {
    Vector e(impl->nt), y(impl->nt);
    std::vector<double> col(nG);
    for (int c = 0; c < nG; c++)
    {
      e = 0.0;
      for (int i = 0; i < impl->nt; i++)
      {
        if (impl->is_gamma[i] && impl->gamma_global[i] == c)
        {
          e(i) = 1.0;
        }
      }
      impl->dtn->Mult(e, y);
      gather_interface(y, col.data());
      for (int r = 0; r < nG; r++)
      {
        impl->S_dense(r, c) = col[r];
      }
    }
    if (!model_path.empty() && rank == 0)
    {
      std::ofstream f(model_path, std::ios::binary);
      f.write(reinterpret_cast<const char *>(&nG), sizeof(int));
      f.write(reinterpret_cast<const char *>(impl->S_dense.GetData()),
              sizeof(double) * nG * nG);
    }
  }
  // g_E is excitation-dependent; it is computed per excitation in the region solve.
  impl->mat_dtn = std::make_unique<MaterializedDtN>(impl->S_dense, impl->gamma_global,
                                                    impl->nG_global, comm);
}

Vector SubstructuringSolver::SolveExcitation(int drive_idx)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  const int nt = impl->nt;
  MPI_Comm comm = impl->parent_fes.GetComm();

  // Prescribed terminal values for this excitation: driven terminal at 1 V, others
  // grounded.
  impl->dbc_values = 0.0;
  for (const auto &[idx, dofs] : impl->terminal_tdofs)
  {
    const double value = (idx == drive_idx) ? 1.0 : 0.0;
    for (int d : dofs)
    {
      impl->dbc_values(d) = value;
    }
  }

  // g_E for this excitation: interface response to the (excitation-specific) Dirichlet
  // data, via one environment solve through the implicit DtN, gathered to the global
  // interface.
  std::vector<double> g_glob(impl->nG_global, 0.0);
  {
    Vector gy(nt);
    impl->dtn->Mult(impl->dbc_values, gy);
    std::vector<double> loc(impl->nG_global, 0.0);
    for (int i = 0; i < nt; i++)
    {
      if (impl->is_gamma[i])
      {
        loc[impl->gamma_global[i]] = gy(i);
      }
    }
    MPI_Allreduce(loc.data(), g_glob.data(), impl->nG_global, MPI_DOUBLE, MPI_SUM, comm);
  }

  // RHS: region Dirichlet elimination + environment DtN load g_E.
  Vector b(nt);
  b = 0.0;
  {
    Vector t(nt);
    impl->A_region->Mult(impl->dbc_values, t);
    for (int i = 0; i < nt; i++)
    {
      if (impl->is_region_free[i])
      {
        b(i) -= t(i);
      }
    }
  }
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_gamma[i])
    {
      b(i) -= g_glob[impl->gamma_global[i]];
    }
  }

  // Region-condensed solve using the materialized S_E (no environment solves in the loop).
  RegionCondensedOperator sysop(*impl->A_region_free, *impl->mat_dtn, impl->is_region_free);
  Vector u(nt);
  u = 0.0;
  mfem::CGSolver cg(comm);
  cg.SetOperator(sysop);
  cg.SetPreconditioner(*impl->prec_region);
  cg.SetRelTol(1.0e-10);
  cg.SetMaxIter(2000);
  cg.SetPrintLevel(0);
  cg.Mult(b, u);
  MFEM_VERIFY(cg.GetConverged(), "Region-condensed CG solve did not converge!");
  for (int i = 0; i < impl->dbc_tdofs.Size(); i++)
  {
    u(impl->dbc_tdofs[i]) = impl->dbc_values(impl->dbc_tdofs[i]);
  }

  // Recover the environment interior: u_E = -A_EE^-1 (A_env u)|_E.
  {
    Vector t(nt);
    impl->A_env->Mult(u, t);
    Vector rhs(nt);
    rhs = 0.0;
    for (int i = 0; i < nt; i++)
    {
      if (impl->is_env_int[i])
      {
        rhs(i) = -t(i);
      }
    }
    Vector uE(nt);
    uE = 0.0;
    impl->solver_env->Mult(rhs, uE);
    for (int i = 0; i < nt; i++)
    {
      if (impl->is_env_int[i])
      {
        u(i) = uE(i);
      }
    }
  }
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
  const int nt = impl->nt;
  MPI_Comm comm = impl->parent_fes.GetComm();

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
  impl->solver_env->Mult(fE, w);
  impl->A_env->Mult(w, Aw);
  std::vector<double> src_glob(impl->nG_global, 0.0), loc(impl->nG_global, 0.0);
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_gamma[i])
    {
      loc[impl->gamma_global[i]] = Aw(i);
    }
  }
  MPI_Allreduce(loc.data(), src_glob.data(), impl->nG_global, MPI_DOUBLE, MPI_SUM, comm);

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
      b(i) -= src_glob[impl->gamma_global[i]];
    }
  }

  RegionCondensedOperator sysop(*impl->A_region_free, *impl->mat_dtn, impl->is_region_free);
  Vector u(nt);
  u = 0.0;
  mfem::CGSolver cg(comm);
  cg.SetOperator(sysop);
  cg.SetPreconditioner(*impl->prec_region);
  cg.SetRelTol(1.0e-10);
  cg.SetMaxIter(2000);
  cg.SetPrintLevel(0);
  cg.Mult(b, u);
  MFEM_VERIFY(cg.GetConverged(), "Region-condensed CG solve did not converge!");

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
  impl->solver_env->Mult(rhs, uE);
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_env_int[i])
    {
      u(i) = uE(i);
    }
  }
  return u;
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
