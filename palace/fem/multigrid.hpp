// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_MULTIGRID_HPP
#define PALACE_FEM_MULTIGRID_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
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
ConstructFECollections(int p, int dim, int mg_max_levels, MultigridCoarsening mg_coarsening,
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
    switch (mg_coarsening)
    {
      case MultigridCoarsening::LINEAR:
        p--;
        break;
      case MultigridCoarsening::LOGARITHMIC:
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
    int mg_max_levels, const std::vector<std::unique_ptr<Mesh>> &mesh,
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
      std::make_unique<FiniteElementSpace>(*mesh[coarse_mesh_l], fecs[0].get()));

  mfem::Array<int> dbc_marker;
  if (dbc_attr && dbc_tdof_lists)
  {
    int bdr_attr_max = mesh[coarse_mesh_l]->Get().bdr_attributes.Size()
                           ? mesh[coarse_mesh_l]->Get().bdr_attributes.Max()
                           : 0;
    dbc_marker = mesh::AttrToMarker(bdr_attr_max, *dbc_attr);
    fespaces.GetFinestFESpace().Get().GetEssentialTrueDofs(dbc_marker,
                                                           dbc_tdof_lists->emplace_back());
  }

  // test copying mesh?
  Mpi::Print("multigrid.hpp L104\n");
std::vector<std::unique_ptr<Mesh>> internal_mesh;
Mpi::Print("multigrid.hpp L106\n");
internal_mesh.reserve(mesh.size());
Mpi::Print("multigrid.hpp L108\n");
for (const auto& m : mesh) {
  Mpi::Print("multigrid.hpp L110\n");
    internal_mesh.push_back(std::make_unique<Mesh>(*m));
}
Mpi::Print("multigrid.hpp L113\n");

  // h-refinement.
  for (std::size_t l = coarse_mesh_l + 1; l < mesh.size(); l++)
  {
    Mpi::Print("l: {}, previous mesh NE: {}, current mesh NE: {}\n", l, mesh[l - 1]->Get().GetGlobalNE(), mesh[l]->Get().GetGlobalNE());
    if (mesh[l - 1]->Get().GetGlobalNE() == mesh[l]->Get().GetGlobalNE())
    {
      // Rebalance
      Mpi::Print("Mesh level {} is a rebalancing, compute transfer operator\n", l);
      auto refine_op = std::make_unique<mfem::TransferOperator>(fespaces.GetFESpaceAtLevel(l - 2), fespaces.GetFESpaceAtLevel(l - 1));

      fespaces.GetFESpaceAtLevel(l - 1).GetParMesh().Rebalance();//? This is bad, modifies the refined but not rebalanced mesh[l - 1]!
      fespaces.GetFESpaceAtLevel(l - 1).GetMesh().Update(); // not sure needed
      fespaces.GetFESpaceAtLevel(l - 1).Update(); // not sure needed
      //auto reb_fe = std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get());
      //auto rebalance_op[l].reset(const_cast<Operator*>(reb_fe.Get().GetUpdateOperator()));
      //auto rebalance_op = const_cast<Operator*>(reb_fe->Get().GetUpdateOperator());
      auto rebalance_op = const_cast<Operator*>(fespaces.GetFESpaceAtLevel(l - 1).Get().GetUpdateOperator());
      std::cout << "rank: " << Mpi::Rank(mesh[l]->GetComm()) << " refine width/height: " << refine_op->Width() << " " << refine_op->Height()
                << " rebalance  width/height: " << rebalance_op->Width() << " " << rebalance_op->Height() << "\n";
      fespaces.UpdateLevel(std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get())); // pass refine and rebalance op too!
      if (dbc_attr && dbc_tdof_lists)
      {
        dbc_tdof_lists->pop_back();
      }
    }
    else
    {
      fespaces.AddLevel(std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get()));
    }
    if (dbc_attr && dbc_tdof_lists)
    {
      fespaces.GetFinestFESpace().Get().GetEssentialTrueDofs(
          dbc_marker, dbc_tdof_lists->emplace_back());
    }
  }

  // p-refinement.
  for (std::size_t l = 1; l < fecs.size(); l++)
  {
    fespaces.AddLevel(std::make_unique<FiniteElementSpace>(*mesh.back(), fecs[l].get()));
    if (dbc_attr && dbc_tdof_lists)
    {
      fespaces.GetFinestFESpace().Get().GetEssentialTrueDofs(
          dbc_marker, dbc_tdof_lists->emplace_back());
    }
  }

  return fespaces;
}

}  // namespace palace::fem

#endif  // PALACE_FEM_MULTIGRID_HPP
