// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_MULTIGRID_HPP
#define PALACE_FEM_MULTIGRID_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "utils/iodata.hpp"

namespace palace::fem
{

//
// Methods for constructing hierarchies of finite element spaces for geometric multigrid.
//

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
// finite element collections. Dirichlet boundary conditions are additionally
// marked.
template <typename FECollection>
inline mfem::ParFiniteElementSpaceHierarchy ConstructFiniteElementSpaceHierarchy(
    int mg_max_levels, bool mg_legacy_transfer, int pa_order_threshold,
    bool pa_discrete_interp, const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    const std::vector<std::unique_ptr<FECollection>> &fecs,
    const mfem::Array<int> *dbc_marker = nullptr,
    std::vector<mfem::Array<int>> *dbc_tdof_lists = nullptr)
{
  MFEM_VERIFY(!mesh.empty() && !fecs.empty() &&
                  (!dbc_tdof_lists || dbc_tdof_lists->empty()),
              "Empty mesh or FE collection for FE space construction!");
  int coarse_mesh_l =
      std::max(0, static_cast<int>(mesh.size() + fecs.size()) - 1 - mg_max_levels);
  auto *fespace = new mfem::ParFiniteElementSpace(mesh[coarse_mesh_l].get(), fecs[0].get());
  if (dbc_marker && dbc_tdof_lists)
  {
    fespace->GetEssentialTrueDofs(*dbc_marker, dbc_tdof_lists->emplace_back());
  }
  mfem::ParFiniteElementSpaceHierarchy fespaces(mesh[coarse_mesh_l].get(), fespace, false,
                                                true);

  // h-refinement
  for (std::size_t l = coarse_mesh_l + 1; l < mesh.size(); l++)
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
      constexpr bool skip_zeros_interp = true;
      DiscreteLinearOperator p(fespaces.GetFinestFESpace(), *fespace);
      p.AddDomainInterpolator(std::make_unique<IdentityInterpolator>());
      P = new ParOperator(
          p.Assemble(pa_discrete_interp ? pa_order_threshold : 99, skip_zeros_interp),
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

}  // namespace palace::fem

#endif  // PALACE_FEM_MULTIGRID_HPP
