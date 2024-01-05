// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_MULTIGRID_HPP
#define PALACE_FEM_MULTIGRID_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace::fem
{

//
// Methods for constructing hierarchies of finite element spaces for geometric multigrid.
//

// Construct sequence of FECollection objects.
template <typename FECollection>
inline std::vector<std::unique_ptr<FECollection>>
ConstructFECollections(int p, int dim, int mg_max_levels,
                       config::LinearSolverData::MultigridCoarsenType mg_coarsen_type,
                       bool mat_lor)
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

// Construct a hierarchy of finite element spaces given a sequence of meshes and finite
// element collections. Additionally, Dirichlet boundary conditions are marked.
template <typename FECollection>
inline FiniteElementSpaceHierarchy ConstructFiniteElementSpaceHierarchy(
    int mg_max_levels, const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    const std::vector<std::unique_ptr<FECollection>> &fecs,
    const mfem::Array<int> *dbc_attr = nullptr,
    std::vector<mfem::Array<int>> *dbc_tdof_lists = nullptr)
{
  MFEM_VERIFY(!mesh.empty() && !fecs.empty() &&
                  (!dbc_tdof_lists || dbc_tdof_lists->empty()),
              "Empty mesh or FE collection for FE space construction!");
  int coarse_mesh_l = std::max(0, static_cast<int>(mesh.size() + fecs.size()) - 1 -
                                      std::max(1, mg_max_levels));
  FiniteElementSpaceHierarchy fespaces(
      std::make_unique<FiniteElementSpace>(mesh[coarse_mesh_l].get(), fecs[0].get()));

  mfem::Array<int> dbc_marker;
  if (dbc_attr && dbc_tdof_lists)
  {
    int bdr_attr_max = mesh[coarse_mesh_l]->bdr_attributes.Size()
                           ? mesh[coarse_mesh_l]->bdr_attributes.Max()
                           : 0;
    dbc_marker = mesh::AttrToMarker(bdr_attr_max, *dbc_attr);
    fespaces.GetFinestFESpace().GetEssentialTrueDofs(dbc_marker,
                                                     dbc_tdof_lists->emplace_back());
  }

  // h-refinement
  for (std::size_t l = coarse_mesh_l + 1; l < mesh.size(); l++)
  {
    fespaces.AddLevel(std::make_unique<FiniteElementSpace>(mesh[l].get(), fecs[0].get()));
    if (dbc_attr && dbc_tdof_lists)
    {
      fespaces.GetFinestFESpace().GetEssentialTrueDofs(dbc_marker,
                                                       dbc_tdof_lists->emplace_back());
    }
  }

  // p-refinement
  for (std::size_t l = 1; l < fecs.size(); l++)
  {
    fespaces.AddLevel(
        std::make_unique<FiniteElementSpace>(mesh.back().get(), fecs[l].get()));
    if (dbc_attr && dbc_tdof_lists)
    {
      fespaces.GetFinestFESpace().GetEssentialTrueDofs(dbc_marker,
                                                       dbc_tdof_lists->emplace_back());
    }
  }

  return fespaces;
}

// Similar to ConstructFiniteElementSpaceHierarchy above, but in this case the finite
// element space at each level is an auxiliary space associated with the coresponding level
// of the provided finite element space objects.
template <typename FECollection>
inline AuxiliaryFiniteElementSpaceHierarchy ConstructAuxiliaryFiniteElementSpaceHierarchy(
    FiniteElementSpaceHierarchy &primal_fespaces,
    const std::vector<std::unique_ptr<FECollection>> &fecs,
    const mfem::Array<int> *dbc_attr = nullptr,
    std::vector<mfem::Array<int>> *dbc_tdof_lists = nullptr)
{
  MFEM_VERIFY((primal_fespaces.GetNumLevels() > 0) && !fecs.empty() &&
                  (!dbc_tdof_lists || dbc_tdof_lists->empty()),
              "Empty mesh or FE collection for FE space construction!");
  mfem::ParMesh *mesh = primal_fespaces.GetFESpaceAtLevel(0).GetParMesh();
  AuxiliaryFiniteElementSpaceHierarchy fespaces(
      std::make_unique<AuxiliaryFiniteElementSpace>(primal_fespaces.GetFESpaceAtLevel(0),
                                                    mesh, fecs[0].get()));

  mfem::Array<int> dbc_marker;
  if (dbc_attr && dbc_tdof_lists)
  {
    int bdr_attr_max = mesh->bdr_attributes.Size() ? mesh->bdr_attributes.Max() : 0;
    dbc_marker = mesh::AttrToMarker(bdr_attr_max, *dbc_attr);
    fespaces.GetFinestFESpace().GetEssentialTrueDofs(dbc_marker,
                                                     dbc_tdof_lists->emplace_back());
  }

  // h-refinement
  std::size_t l;
  for (l = 1; l < primal_fespaces.GetNumLevels(); l++)
  {
    if (primal_fespaces.GetFESpaceAtLevel(l).GetParMesh() == mesh)
    {
      break;
    }
    fespaces.AddLevel(std::make_unique<AuxiliaryFiniteElementSpace>(
        primal_fespaces.GetFESpaceAtLevel(l),
        primal_fespaces.GetFESpaceAtLevel(l).GetParMesh(), fecs[0].get()));
    if (dbc_attr && dbc_tdof_lists)
    {
      fespaces.GetFinestFESpace().GetEssentialTrueDofs(dbc_marker,
                                                       dbc_tdof_lists->emplace_back());
    }

    mesh = primal_fespaces.GetFESpaceAtLevel(l).GetParMesh();
  }

  // p-refinement
  const auto l0 = l - 1;
  for (; l < primal_fespaces.GetNumLevels(); l++)
  {
    fespaces.AddLevel(std::make_unique<AuxiliaryFiniteElementSpace>(
        primal_fespaces.GetFESpaceAtLevel(l), mesh, fecs[l - l0].get()));
    if (dbc_attr && dbc_tdof_lists)
    {
      fespaces.GetFinestFESpace().GetEssentialTrueDofs(dbc_marker,
                                                       dbc_tdof_lists->emplace_back());
    }
  }

  return fespaces;
}

}  // namespace palace::fem

#endif  // PALACE_FEM_MULTIGRID_HPP
