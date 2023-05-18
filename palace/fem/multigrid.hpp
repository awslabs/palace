// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_MULTIGRID_HPP
#define PALACE_FEM_MULTIGRID_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"

namespace palace::utils
{

//
// Methods for constructing hierarchies of finite element spaces for geometric multigrid.
//

// Construct sequence of FECollection objects.
template <typename FECollection>
std::vector<std::unique_ptr<FECollection>> ConstructFECollections(bool pc_pmg, bool pc_lor,
                                                                  int p, int dim)
{
  // If the solver will use a LOR preconditioner, we need to construct with a specific basis
  // type.
  constexpr int pmin = (std::is_same<FECollection, mfem::H1_FECollection>::value ||
                        std::is_same<FECollection, mfem::ND_FECollection>::value)
                           ? 1
                           : 0;
  MFEM_VERIFY(p >= pmin, "FE space order must not be less than " << pmin << "!");
  int b1 = mfem::BasisType::GaussLobatto, b2 = mfem::BasisType::GaussLegendre;
  if (pc_lor)
  {
    b2 = mfem::BasisType::IntegratedGLL;
  }
  std::vector<std::unique_ptr<FECollection>> fecs;
  if (pc_pmg)
  {
    fecs.reserve(p);
    for (int o = pmin; o <= p; o++)
    {
      if constexpr (std::is_same<FECollection, mfem::ND_FECollection>::value ||
                    std::is_same<FECollection, mfem::RT_FECollection>::value)
      {
        fecs.push_back(std::make_unique<FECollection>(o, dim, b1, b2));
      }
      else
      {
        fecs.push_back(std::make_unique<FECollection>(o, dim, b1));
      }
    }
  }
  else
  {
    fecs.reserve(1);
    if constexpr (std::is_same<FECollection, mfem::ND_FECollection>::value ||
                  std::is_same<FECollection, mfem::RT_FECollection>::value)
    {
      fecs.push_back(std::make_unique<FECollection>(p, dim, b1, b2));
    }
    else
    {
      fecs.push_back(std::make_unique<FECollection>(p, dim, b1));
    }
  }
  return fecs;
}

// Construct a hierarchy of finite element spaces given a sequence of meshes and
// finite element collections. Dirichlet boundary conditions are additionally
// marked.
template <typename FECollection>
mfem::ParFiniteElementSpaceHierarchy ConstructFiniteElementSpaceHierarchy(
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    const std::vector<std::unique_ptr<FECollection>> &fecs,
    const mfem::Array<int> *dbc_marker = nullptr,
    std::vector<mfem::Array<int>> *dbc_tdof_lists = nullptr)
{
  MFEM_VERIFY(!mesh.empty() && !fecs.empty() &&
                  (!dbc_tdof_lists || dbc_tdof_lists->empty()),
              "Empty mesh or FE collection for FE space construction!");
  auto *fespace = new mfem::ParFiniteElementSpace(mesh[0].get(), fecs[0].get());
  if (dbc_marker && dbc_tdof_lists)
  {
    fespace->GetEssentialTrueDofs(*dbc_marker, dbc_tdof_lists->emplace_back());
  }
  mfem::ParFiniteElementSpaceHierarchy fespaces(mesh[0].get(), fespace, false, true);

  // XX TODO: LibCEED transfer operators!

  // h-refinement
  for (std::size_t l = 1; l < mesh.size(); l++)
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
    auto *P = new ParOperator(
        std::make_unique<mfem::TransferOperator>(fespaces.GetFinestFESpace(), *fespace),
        fespaces.GetFinestFESpace(), *fespace, true);
    fespaces.AddLevel(mesh.back().get(), fespace, P, false, true, true);
  }
  return fespaces;
}

// Construct a single-level finite element space hierarchy from a single mesh and
// finite element collection. Unnecessary to pass the dirichlet boundary
// conditions as they need not be incorporated in any inter-space projectors.
template <typename FECollection>
mfem::ParFiniteElementSpaceHierarchy
ConstructFiniteElementSpaceHierarchy(mfem::ParMesh &mesh, const FECollection &fec,
                                     const mfem::Array<int> *dbc_marker = nullptr,
                                     mfem::Array<int> *dbc_tdof_list = nullptr)
{
  auto *fespace = new mfem::ParFiniteElementSpace(&mesh, &fec);
  if (dbc_marker && dbc_tdof_list)
  {
    fespace->GetEssentialTrueDofs(*dbc_marker, *dbc_tdof_list);
  }
  return mfem::ParFiniteElementSpaceHierarchy(&mesh, fespace, false, true);
}

}  // namespace palace::utils

#endif  // PALACE_FEM_MULTIGRID_HPP
