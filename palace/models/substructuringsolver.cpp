// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolver.hpp"

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

// Parallel region-condensed electrostatic solve. Region/environment operators are assembled
// on the parent finite element space with domain-restricted (isotropic scalar) permittivity
// coefficients; the interface is identified in true-DOF space; the environment is condensed
// through an implicit distributed DtN. Reusing LaplaceOperator/KspSolver, tensor materials,
// the capacitance sweep, and reuse-optimized (materialized) DtN are follow-ups.
struct SubstructuringSolver::Impl
{
  const IoData &iodata;
  mfem::ParMesh &parent;
  mfem::H1_FECollection fec;
  mfem::ParFiniteElementSpace parent_fes;
  int nt;

  std::unique_ptr<mfem::HypreParMatrix> A_region, A_env, A_env_int, A_region_free;
  std::vector<char> is_gamma, is_env_int, is_region_free;
  mfem::Array<int> dbc_tdofs;
  mfem::Vector dbc_values;  // full parent true-DOF vector, prescribed values on Dirichlet

  std::unique_ptr<mfem::HypreBoomerAMG> amg_env;
  std::unique_ptr<mfem::HyprePCG> solver_env;
  std::unique_ptr<ImplicitDtN> dtn;

  // Materialized (reusable) interface operator: replicated dense S_E + load g_E over a
  // global interface enumeration, computed once so region solves need no environment
  // solves.
  std::vector<int> gamma_global;  // owned parent true DOF -> global interface index, or -1
  int nG_global = 0;
  mfem::DenseMatrix S_dense;
  mfem::Vector g_dense;
  std::unique_ptr<MaterializedDtN> mat_dtn;

  Impl(const IoData &iodata, mfem::ParMesh &parent)
    : iodata(iodata), parent(parent), fec(iodata.solver.order, parent.Dimension()),
      parent_fes(&parent, &fec), nt(parent_fes.GetTrueVSize())
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

    // Dirichlet terminals (single excitation: lowest index -> 1 V, others grounded).
    const auto &terminals = iodata.boundaries.terminal;
    mfem::Array<int> dir_mark(nt);
    dir_mark = 0;
    dbc_values.SetSize(nt);
    dbc_values = 0.0;
    if (!terminals.empty())
    {
      const int drive = terminals.begin()->first;
      const int maxb = parent.bdr_attributes.Size() ? parent.bdr_attributes.Max() : 0;
      for (const auto &[idx, term] : terminals)
      {
        const double value = (idx == drive) ? 1.0 : 0.0;
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
        for (int i = 0; i < ess.Size(); i++)
        {
          dir_mark[ess[i]] = 1;
          dbc_values(ess[i]) = value;
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
    mfem::Vector er(max_attr), ee(max_attr);
    er = 0.0;
    ee = 0.0;
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
          er(a - 1) = mat.epsilon_r.s[0];
        }
        else if (in(ea, a))
        {
          ee(a - 1) = mat.epsilon_r.s[0];
        }
      }
    }
    A_region = AssembleParent(er);
    A_env = AssembleParent(ee);
  }

  std::unique_ptr<mfem::HypreParMatrix> AssembleParent(const mfem::Vector &eps_by_attr)
  {
    mfem::PWConstCoefficient eps(const_cast<mfem::Vector &>(eps_by_attr));
    mfem::ParBilinearForm a(&parent_fes);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(eps));
    a.Assemble();
    a.Finalize();
    return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
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
  impl->amg_env = std::make_unique<mfem::HypreBoomerAMG>(*impl->A_env_int);
  impl->amg_env->SetPrintLevel(0);
  impl->solver_env = std::make_unique<mfem::HyprePCG>(*impl->A_env_int);
  impl->solver_env->SetTol(1.0e-13);
  impl->solver_env->SetMaxIter(1000);
  impl->solver_env->SetPrintLevel(0);
  impl->solver_env->SetPreconditioner(*impl->amg_env);
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
  // g_E: interface response to the environment Dirichlet data.
  impl->g_dense.SetSize(nG);
  {
    Vector gy(impl->nt);
    impl->dtn->Mult(impl->dbc_values, gy);
    std::vector<double> gc(nG);
    gather_interface(gy, gc.data());
    for (int r = 0; r < nG; r++)
    {
      impl->g_dense(r) = gc[r];
    }
  }
  impl->mat_dtn = std::make_unique<MaterializedDtN>(impl->S_dense, impl->gamma_global,
                                                    impl->nG_global, comm);
}

Vector SubstructuringSolver::SolveRegion()
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before SolveRegion!");
  const int nt = impl->nt;

  // RHS: region Dirichlet elimination + the materialized environment DtN load g_E.
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
      b(i) -= impl->g_dense(impl->gamma_global[i]);
    }
  }

  // Region-condensed solve using the materialized S_E (no environment solves in the loop).
  RegionCondensedOperator sysop(*impl->A_region_free, *impl->mat_dtn, impl->is_region_free);
  Vector u(nt);
  u = 0.0;
  mfem::CGSolver cg(impl->parent_fes.GetComm());
  cg.SetOperator(sysop);
  cg.SetRelTol(1.0e-10);
  cg.SetMaxIter(2000);
  cg.SetPrintLevel(0);
  cg.Mult(b, u);
  MFEM_VERIFY(cg.GetConverged(), "Region-condensed CG solve did not converge!");
  for (int i = 0; i < impl->dbc_tdofs.Size(); i++)
  {
    u(impl->dbc_tdofs[i]) = impl->dbc_values(impl->dbc_tdofs[i]);
  }
  return u;
}

long long int SubstructuringSolver::RegionGlobalTrueVSize() const
{
  return impl->parent_fes.GlobalTrueVSize();
}

}  // namespace palace
