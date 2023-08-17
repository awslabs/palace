// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_MULTIGRID_HPP
#define PALACE_FEM_MULTIGRID_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "utils/iodata.hpp"

namespace palace::utils
{

//
// Methods for constructing hierarchies of finite element spaces for multigrid.
//

// Helper function for getting the order of the finite element space underlying a bilinear
// form.
inline auto GetMaxElementOrder(mfem::BilinearForm &a)
{
  return a.FESpace()->GetMaxElementOrder();
}

// Helper function for getting the order of the finite element space underlying a mixed
// bilinear form.
inline auto GetMaxElementOrder(mfem::MixedBilinearForm &a)
{
  return std::max(a.TestFESpace()->GetMaxElementOrder(),
                  a.TrialFESpace()->GetMaxElementOrder());
}

// Assemble a bilinear or mixed bilinear form. If the order is lower than the specified
// threshold, the operator is assembled as a sparse matrix.
template <typename BilinearForm>
inline std::unique_ptr<Operator>
AssembleOperator(std::unique_ptr<BilinearForm> &&a, bool mfem_pa_support,
                 int pa_order_threshold, int skip_zeros = 1)
{
  mfem::AssemblyLevel assembly_level =
      (mfem::DeviceCanUseCeed() ||
       (mfem_pa_support && GetMaxElementOrder(*a) >= pa_order_threshold))
          ? mfem::AssemblyLevel::PARTIAL
          : mfem::AssemblyLevel::LEGACY;
  a->SetAssemblyLevel(assembly_level);
  a->Assemble(skip_zeros);
  a->Finalize(skip_zeros);
  if (assembly_level == mfem::AssemblyLevel::LEGACY ||
      (assembly_level == mfem::AssemblyLevel::PARTIAL &&
       GetMaxElementOrder(*a) < pa_order_threshold &&
       std::is_base_of<mfem::BilinearForm, BilinearForm>::value))
  {
    // libCEED full assembly does not support mixed forms.
#ifdef MFEM_USE_CEED
    mfem::SparseMatrix *spm =
        a->HasExt() ? mfem::ceed::CeedOperatorFullAssemble(*a) : a->LoseMat();
#else
    mfem::SparseMatrix *spm = a->LoseMat();
#endif
    MFEM_VERIFY(spm, "Missing assembled sparse matrix!");
    return std::unique_ptr<Operator>(spm);
  }
  else
  {
    return std::move(a);
  }
}

// Construct sequence of FECollection objects.
template <typename FECollection>
std::vector<std::unique_ptr<FECollection>> inline ConstructFECollections(
    int p, int dim, int mg_max_levels,
    config::LinearSolverData::MultigridCoarsenType mg_coarsen_type, bool mat_lor)
{
  // If the solver will use a LOR preconditioner, we need to construct with a specific basis
  // type.
  constexpr int pmin = (std::is_base_of<mfem::H1_FECollection, FECollection>::value ||
                        std::is_base_of<mfem::ND_FECollection, FECollection>::value)
                           ? 1
                           : 0;
  MFEM_VERIFY(p >= pmin, "FE space order must not be less than " << pmin << "!");
  int b1 = mfem::BasisType::GaussLobatto, b2 = mfem::BasisType::GaussLegendre;
  if (mat_lor)
  {
    b2 = mfem::BasisType::IntegratedGLL;
  }

  // Construct the p-multigrid hierarchy, first finest to coarsest and then reverse the
  // order.
  std::vector<std::unique_ptr<FECollection>> fecs;
  for (int l = 0; l < std::max(1, mg_max_levels); l++)
  {
    if constexpr (std::is_base_of<mfem::ND_FECollection, FECollection>::value ||
                  std::is_base_of<mfem::RT_FECollection, FECollection>::value)
    {
      fecs.push_back(std::make_unique<FECollection>(p, dim, b1, b2));
    }
    else
    {
      fecs.push_back(std::make_unique<FECollection>(p, dim, b1));
      MFEM_CONTRACT_VAR(b2);
    }
    if (p == pmin)
    {
      break;
    }
    switch (mg_coarsen_type)
    {
      case config::LinearSolverData::MultigridCoarsenType::LINEAR:
        p--;
        break;
      case config::LinearSolverData::MultigridCoarsenType::LOGARITHMIC:
        p = (p + pmin) / 2;
        break;
    }
  }
  std::reverse(fecs.begin(), fecs.end());
  return fecs;
}

// Construct a hierarchy of finite element spaces given a sequence of meshes and
// finite element collections. Uses geometric multigrid and p multigrid.
template <typename FECollection>
inline mfem::ParFiniteElementSpaceHierarchy ConstructFiniteElementSpaceHierarchy(
    int mg_max_levels, bool mg_legacy_transfer, int pa_order_threshold,
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    const std::vector<std::unique_ptr<FECollection>> &fecs, int dim = 1,
    int ordering = mfem::Ordering::byNODES)
{
  MFEM_VERIFY(!mesh.empty() && !fecs.empty(),
              "Empty mesh or FE collection for FE space construction!");
  int coarse_mesh_l =
      std::max(0, static_cast<int>(mesh.size() + fecs.size()) - 1 - mg_max_levels);
  auto *fespace = new mfem::ParFiniteElementSpace(mesh[coarse_mesh_l].get(), fecs[0].get(),
                                                  dim, ordering);
  mfem::ParFiniteElementSpaceHierarchy fespaces(mesh[coarse_mesh_l].get(), fespace, false,
                                                true);

  // h-refinement
  for (std::size_t l = coarse_mesh_l + 1; l < mesh.size(); l++)
  {
    fespace = new mfem::ParFiniteElementSpace(mesh[l].get(), fecs[0].get(), dim, ordering);

    auto *P = new ParOperator(
        std::make_unique<mfem::TransferOperator>(fespaces.GetFinestFESpace(), *fespace),
        fespaces.GetFinestFESpace(), *fespace, true);
    fespaces.AddLevel(mesh[l].get(), fespace, P, false, true, true);
  }

  // p-refinement
  for (std::size_t l = 1; l < fecs.size(); l++)
  {
    fespace =
        new mfem::ParFiniteElementSpace(mesh.back().get(), fecs[l].get(), dim, ordering);

    ParOperator *P;
    if (!mg_legacy_transfer && mfem::DeviceCanUseCeed())
    {
      // Partial and full assembly for this operator is only available with libCEED backend.
      auto p = std::make_unique<mfem::DiscreteLinearOperator>(&fespaces.GetFinestFESpace(),
                                                              fespace);
      p->AddDomainInterpolator(new mfem::IdentityInterpolator);
      P = new ParOperator(AssembleOperator(std::move(p), false, pa_order_threshold),
                          fespaces.GetFinestFESpace(), *fespace, true);
    }
    else
    {
      P = new ParOperator(
          std::make_unique<mfem::TransferOperator>(fespaces.GetFinestFESpace(), *fespace),
          fespaces.GetFinestFESpace(), *fespace, true);
    }
    fespaces.AddLevel(mesh.back().get(), fespace, P, false, true, true);
  }
  return fespaces;
}

// Overload for treatment of Dirichlet boundary conditions, extracts the true
// dof vectors for each level within a finite element space hierarchy.
template <typename FECollection, template <class... C> typename Container,
          typename... MeshT>
mfem::ParFiniteElementSpaceHierarchy ConstructFiniteElementSpaceHierarchy(
    int mg_max_levels, bool mg_legacy_transfer, int pa_order_threshold,
    const Container<MeshT...> &mesh, const std::vector<std::unique_ptr<FECollection>> &fecs,
    const mfem::Array<int> &dbc_marker, std::vector<mfem::Array<int>> &dbc_tdof_lists,
    int dim = 1, int ordering = mfem::Ordering::byNODES)
{
  dbc_tdof_lists.clear();
  auto fespaces = ConstructFiniteElementSpaceHierarchy(
      mg_max_levels, mg_legacy_transfer, pa_order_threshold, mesh, fecs, dim, ordering);
  for (int l = 0; l < fespaces.GetNumLevels(); l++)
  {
    fespaces.GetFESpaceAtLevel(l).GetEssentialTrueDofs(dbc_marker,
                                                       dbc_tdof_lists.emplace_back());
  }
  return fespaces;
}

}  // namespace palace::utils

#endif  // PALACE_FEM_MULTIGRID_HPP
