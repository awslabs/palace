// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolver.hpp"

#include <mfem.hpp>
#include "fem/mesh.hpp"
#include "fem/substructure.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"

namespace palace
{

// Implementation note: for this first electrostatic increment the region and environment
// operators are assembled directly with a BilinearForm using a scalar (isotropic)
// permittivity coefficient from the material configuration, and the region is solved with
// CG. Reusing LaplaceOperator + KspSolver/AMG (and tensor materials, postprocessing, the
// full capacitance sweep) is a follow-up; the public interface is unchanged by that.
struct SubstructuringSolver::Impl
{
  const IoData &iodata;
  mfem::ParMesh &parent;
  mfem::H1_FECollection fec;
  mfem::ParFiniteElementSpace parent_fes;
  std::unique_ptr<Substructure> region, environment;
  std::vector<int> owner, gamma_index;
  mfem::Vector eps_by_attr;  // relative permittivity indexed by (attribute - 1)
  std::unique_ptr<DtNBoundaryOperator> dtn;

  Impl(const IoData &iodata, mfem::ParMesh &parent)
    : iodata(iodata), parent(parent), fec(iodata.solver.order, parent.Dimension()),
      parent_fes(&parent, &fec)
  {
    const auto &sub = *iodata.solver.substructuring;
    mfem::Array<int> ra(sub.region_attributes.size()),
        ea(sub.environment_attributes.size());
    std::copy(sub.region_attributes.begin(), sub.region_attributes.end(), ra.begin());
    std::copy(sub.environment_attributes.begin(), sub.environment_attributes.end(),
              ea.begin());
    region = std::make_unique<Substructure>(parent_fes, ra, fec);
    environment = std::make_unique<Substructure>(parent_fes, ea, fec);
    MarkParentDofOwnership({region.get(), environment.get()}, parent_fes.GetVSize(), owner);
    BuildInterfaceIndex(owner, gamma_index);

    // Scalar (isotropic) relative permittivity per attribute from the material config.
    int max_attr = parent.attributes.Size() ? parent.attributes.Max() : 1;
    eps_by_attr.SetSize(max_attr);
    eps_by_attr = 1.0;
    for (const auto &mat : iodata.domains.materials)
    {
      for (int a : mat.attributes)
      {
        if (a >= 1 && a <= max_attr)
        {
          eps_by_attr(a - 1) = mat.epsilon_r.s[0];
        }
      }
    }
  }

  // Assemble a substructure's Laplace (grad eps grad) operator on its submesh FE space.
  std::unique_ptr<mfem::SparseMatrix> AssembleStiffness(Substructure &s)
  {
    mfem::PWConstCoefficient eps(eps_by_attr);
    mfem::BilinearForm a(&s.GetFESpace());
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(eps));
    a.Assemble();
    a.Finalize();
    return std::make_unique<mfem::SparseMatrix>(a.SpMat());
  }

  // Build the Dirichlet terminal marker/values for a substructure, for a single excitation:
  // the terminal with the lowest index is driven to 1, all other terminals grounded (0).
  void BuildTerminals(Substructure &s, std::vector<char> &marker, mfem::Vector &vals)
  {
    auto &fes = s.GetFESpace();
    const int nv = fes.GetVSize();
    marker.assign(nv, 0);
    vals.SetSize(nv);
    vals = 0.0;
    const auto &terminals = iodata.boundaries.terminal;
    if (terminals.empty())
    {
      return;
    }
    const int drive_idx = terminals.begin()->first;
    const int max_bdr =
        s.GetSubMesh().bdr_attributes.Size() ? s.GetSubMesh().bdr_attributes.Max() : 0;
    for (const auto &[idx, term] : terminals)
    {
      const double value = (idx == drive_idx) ? 1.0 : 0.0;
      mfem::Array<int> ess_bdr(max_bdr);
      ess_bdr = 0;
      for (int a : term.attributes)
      {
        if (a >= 1 && a <= max_bdr)
        {
          ess_bdr[a - 1] = 1;
        }
      }
      mfem::Array<int> ess_vdofs;
      fes.GetEssentialVDofs(ess_bdr, ess_vdofs);
      for (int i = 0; i < nv; i++)
      {
        if (ess_vdofs[i])
        {
          marker[i] = 1;
          vals(i) = value;
        }
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
  MFEM_VERIFY(impl->region->ConformingMap() && impl->environment->ConformingMap(),
              "Substructuring requires a conforming interface (order-1, matching mesh)!");
}

SubstructuringSolver::~SubstructuringSolver() = default;

void SubstructuringSolver::CondenseEnvironment()
{
  auto A_env = impl->AssembleStiffness(*impl->environment);
  std::vector<char> marker;
  mfem::Vector vals;
  impl->BuildTerminals(*impl->environment, marker, vals);
  mfem::Vector f_env(A_env->Height());
  f_env = 0.0;
  impl->dtn = std::make_unique<DtNBoundaryOperator>(*impl->environment, impl->gamma_index,
                                                    *A_env, f_env, marker, vals);
}

Vector SubstructuringSolver::SolveRegion()
{
  MFEM_VERIFY(impl->dtn, "CondenseEnvironment must be called before SolveRegion!");
  auto &region = *impl->region;
  auto A_region = impl->AssembleStiffness(region);
  const auto &par = region.GetParentDof();
  const int nr = A_region->Height();

  std::vector<char> rm;
  mfem::Vector rv;
  impl->BuildTerminals(region, rm, rv);

  // Compact region free system: interface (gamma index) + region interior (not Dirichlet,
  // not interface), with region Dirichlet terminals eliminated into the RHS.
  const int nG = impl->dtn->Size();
  const int N = impl->parent_fes.GetVSize();
  std::vector<int> ri(N, -1), rint;
  for (int i = 0; i < nr; i++)
  {
    if (!rm[i] && impl->gamma_index[par[i]] < 0 && ri[par[i]] < 0)
    {
      ri[par[i]] = static_cast<int>(rint.size());
      rint.push_back(par[i]);
    }
  }
  const int ndof = nG + static_cast<int>(rint.size());
  auto idx = [&](int i)
  { return impl->gamma_index[par[i]] >= 0 ? impl->gamma_index[par[i]] : nG + ri[par[i]]; };

  mfem::DenseMatrix K(ndof);
  K = 0.0;
  mfem::Vector b(ndof);
  b = 0.0;
  for (int i = 0; i < nr; i++)
  {
    if (rm[i])
    {
      continue;
    }
    const int *cols = A_region->GetRowColumns(i);
    const double *vals = A_region->GetRowEntries(i);
    for (int k = 0; k < A_region->RowSize(i); k++)
    {
      const int j = cols[k];
      const double v = vals[k];
      if (rm[j])
      {
        b(idx(i)) -= v * rv(j);
        continue;
      }
      K(idx(i), idx(j)) += v;
    }
  }
  for (int a = 0; a < nG; a++)
  {
    b(a) += impl->dtn->Load()(a);
    for (int bb = 0; bb < nG; bb++)
    {
      K(a, bb) += impl->dtn->Schur()(a, bb);
    }
  }

  mfem::DenseMatrix Ki(K);
  Ki.Invert();
  mfem::Vector u(ndof);
  Ki.Mult(b, u);

  // Scatter to region submesh DOFs (Dirichlet DOFs carry their prescribed value).
  Vector u_region(nr);
  u_region = 0.0;
  for (int i = 0; i < nr; i++)
  {
    u_region(i) = rm[i] ? rv(i) : u(idx(i));
  }
  return u_region;
}

long long int SubstructuringSolver::RegionGlobalTrueVSize() const
{
  return impl->region->GetFESpace().GlobalTrueVSize();
}

}  // namespace palace
