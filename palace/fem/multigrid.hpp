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
// Methods for constructing hierarchies of finite element spaces for geometric multigrid.
//

// Helper function for accessing the finite element space from a bilinear form.
inline auto GetFESpace(mfem::BilinearForm &a)
{
  return *a.FESpace();
}

// Helper function for accessing the finite element space from a mixed bilinear form.
inline auto GetFESpace(mfem::MixedBilinearForm &a)
{
  return *a.TrialFESpace();
}

// Get the default assembly level for a bilinear or mixed bilinear form.
inline auto GetAssemblyLevel(int order, int pa_order_threshold)
{
  return (mfem::DeviceCanUseCeed() || order > pa_order_threshold)
             ? mfem::AssemblyLevel::PARTIAL
             : mfem::AssemblyLevel::LEGACY;
}

// Get the default assembly level for a bilinear or mixed bilinear form which only supports
// partial assembly for the libCEED backend.
inline auto GetAssemblyLevel()
{
  return mfem::DeviceCanUseCeed() ? mfem::AssemblyLevel::PARTIAL
                                  : mfem::AssemblyLevel::LEGACY;
}

// Assembly a bilinear or mixed bilinear form. If the order is lower than the specified
// threshold, the operator is assembled as a sparse matrix.
template <typename BilinearForm>
inline std::unique_ptr<Operator> AssembleOperator(std::unique_ptr<BilinearForm> &&a,
                                                  int pa_order_threshold)
{
  if (a->GetAssemblyLevel() == mfem::AssemblyLevel::PARTIAL &&
      GetFESpace(*a).GetMaxElementOrder() > pa_order_threshold)
  {
    return std::move(a);
  }
  else
  {
#ifdef MFEM_USE_CEED
    mfem::SparseMatrix *spm =
        a->HasExt() ? mfem::ceed::CeedOperatorFullAssemble(*a) : a->LoseMat();
#else
    mfem::SparseMatrix *spm = a->LoseMat();
#endif
    MFEM_VERIFY(spm, "Missing assembled sparse matrix!");
    return std::unique_ptr<Operator>(spm);
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
  MFEM_VERIFY(p >= 1, "FE space order must not be less than 1!");
  int b1 = mfem::BasisType::GaussLobatto, b2 = mfem::BasisType::GaussLegendre;
  if (mat_lor)
  {
    b2 = mfem::BasisType::IntegratedGLL;
  }
  constexpr int pm1 = (std::is_same<FECollection, mfem::H1_FECollection>::value ||
                       std::is_same<FECollection, mfem::ND_FECollection>::value)
                          ? 0
                          : 1;

  // Construct the p-multigrid hierarchy, first finest to coarsest and then reverse the
  // order.
  std::vector<std::unique_ptr<FECollection>> fecs;
  for (int l = 0; l < std::max(1, mg_max_levels); l++)
  {
    if constexpr (std::is_same<FECollection, mfem::ND_FECollection>::value ||
                  std::is_same<FECollection, mfem::RT_FECollection>::value)
    {
      fecs.push_back(std::make_unique<FECollection>(p - pm1, dim, b1, b2));
    }
    else
    {
      fecs.push_back(std::make_unique<FECollection>(p - pm1, dim, b1));
      MFEM_CONTRACT_VAR(b2);
    }
    if (p == 1)
    {
      break;
    }
    switch (mg_coarsen_type)
    {
      case config::LinearSolverData::MultigridCoarsenType::LINEAR:
        p--;
        break;
      case config::LinearSolverData::MultigridCoarsenType::LOGARITHMIC:
        p = (p + 1) / 2;
        break;
      default:
        MFEM_ABORT("Invalid coarsening type for p-multigrid levels!");
        break;
    }
  }
  std::reverse(fecs.begin(), fecs.end());
  return fecs;
}

// Construct a hierarchy of finite element spaces given a sequence of meshes and
// finite element collections. Dirichlet boundary conditions are additionally
// marked.
template <typename FECollection>
inline mfem::ParFiniteElementSpaceHierarchy ConstructFiniteElementSpaceHierarchy(
    int mg_max_levels, bool mg_legacy_transfer, int pa_order_threshold,
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    const std::vector<std::unique_ptr<FECollection>> &fecs,
    const mfem::Array<int> *dbc_marker = nullptr,
    std::vector<mfem::Array<int>> *dbc_tdof_lists = nullptr)
{
  MFEM_VERIFY(!mesh.empty() && !fecs.empty() &&
                  (!dbc_tdof_lists || dbc_tdof_lists->empty()),
              "Empty mesh or FE collection for FE space construction!");
  auto mesh_levels = std::min(mesh.size() - 1, mg_max_levels - fecs.size());
  auto *fespace = new mfem::ParFiniteElementSpace(mesh[mesh.size() - mesh_levels - 1].get(),
                                                  fecs[0].get());
  if (dbc_marker && dbc_tdof_lists)
  {
    fespace->GetEssentialTrueDofs(*dbc_marker, dbc_tdof_lists->emplace_back());
  }
  mfem::ParFiniteElementSpaceHierarchy fespaces(mesh[mesh.size() - mesh_levels - 1].get(),
                                                fespace, false, true);

  // h-refinement
  for (std::size_t l = mesh.size() - mesh_levels; l < mesh.size(); l++)
  {
    fespace = new mfem::ParFiniteElementSpace(mesh[l].get(), fecs[0].get());
    if (dbc_marker && dbc_tdof_lists)
    {
      fespace->GetEssentialTrueDofs(*dbc_marker, dbc_tdof_lists->emplace_back());
    }
    auto *P = new ParOperator(
        std::make_unique<mfem::TransferOperator>(fespaces.GetFinestFESpace(), *fespace),
        fespaces.GetFinestFESpace(), *fespace, true);
    fespaces.AddLevel(mesh[l].get(), fespace, P, false, true, true);
  }

  // p-refinement
  for (std::size_t l = 1; l < fecs.size(); l++)
  {
    fespace = new mfem::ParFiniteElementSpace(mesh.back().get(), fecs[l].get());
    if (dbc_marker && dbc_tdof_lists)
    {
      fespace->GetEssentialTrueDofs(*dbc_marker, dbc_tdof_lists->emplace_back());
    }
    ParOperator *P;
    if (!mg_legacy_transfer)
    {
      // Partial assembly for this operator is only available with libCEED backend.
      auto p = std::make_unique<mfem::DiscreteLinearOperator>(&fespaces.GetFinestFESpace(),
                                                              fespace);
      p->AddDomainInterpolator(new mfem::IdentityInterpolator);
      p->SetAssemblyLevel(GetAssemblyLevel());
      p->Assemble();
      p->Finalize();
      P = new ParOperator(AssembleOperator(std::move(p), pa_order_threshold),
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

}  // namespace palace::utils

#endif  // PALACE_FEM_MULTIGRID_HPP
