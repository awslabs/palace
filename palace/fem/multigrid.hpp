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
    int mg_max_levels, /*const*/ std::vector<std::unique_ptr<Mesh>> &mesh,
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
  Mpi::Print("In ConstructFiniteElementSpaceHierarchy\n");
/*
  // test copying mesh?
  Mpi::Print("multigrid.hpp L104\n");
//std::vector<Mesh> internal_mesh;
std::vector<std::unique_ptr<Mesh>> internal_mesh;
Mpi::Print("multigrid.hpp L106\n");
internal_mesh.reserve(mesh.size());
Mpi::Print("multigrid.hpp L108\n");
for (const auto& m : mesh)
{
  Mpi::Print("multigrid.hpp L110\n");
  internal_mesh.emplace_back(std::make_unique<Mesh>(*m));
  //internal_mesh.emplace_back(Mesh(*m));
  Mpi::Print("multigrid.hpp L112\n");
}
Mpi::Print("multigrid.hpp L114\n");
*/
  int n_rebal = 0;
  int n_skipped = 0;
  int prev_sequence = mesh[coarse_mesh_l]->Get().GetSequence();
  // h-refinement.
  for (std::size_t l = coarse_mesh_l + 1; l < mesh.size(); l++)
  {
    Mpi::Print("l: {}, previous mesh NE: {}, current mesh NE: {}, sequence: {}, last op is refine: {}, last op is rebalance: {}, \n", l, mesh[l - 1]->Get().GetGlobalNE(), mesh[l]->Get().GetGlobalNE(), mesh[l]->Get().GetSequence(), (mesh[l]->Get().GetLastOperation() == mfem::Mesh::REFINE), (mesh[l]->Get().GetLastOperation() == mfem::Mesh::REBALANCE));
    if (mesh[l]->Get().GetSequence() != prev_sequence + 1)
    {
      Mpi::Print("Skipping it because prev_sequence: {}\n", prev_sequence);
      n_skipped++;
      continue; // duplicate mesh, skipping it!
    }
    else if (mesh[l]->Get().GetLastOperation() == mfem::Mesh::REBALANCE)
    {
      // Rebalance
      Mpi::Print("Mesh level {} is a rebalancing (n_rebal: {}), compute transfer operator\n", l, n_rebal);
      auto refine_op = std::make_unique<mfem::TransferOperator>(fespaces.GetFESpaceAtLevel(l - 2 - n_rebal - n_skipped), fespaces.GetFESpaceAtLevel(l - 1 - n_rebal - n_skipped));
     // Mpi::Print("call mesh[l-1]->Update()\n");
     // mesh[l - 1]->Update();
     // Mpi::Print("mesh[l - 1] sequence: {}, last operation == REFINE: {}, == REBALANCE: {}\n", mesh[l - 1]->Get().GetSequence(), mesh[l - 1]->Get().GetLastOperation() == mfem::Mesh::REFINE, mesh[l - 1]->Get().GetLastOperation() == mfem::Mesh::REBALANCE);
     // //Mpi::Print("create mesh_copy\n");
     // //auto mesh_copy = std::make_unique<Mesh>(false, *mesh[l - 1]);
     // Mpi::Print("create pmesh_copy\n");
     // auto pmesh_copy = std::make_unique<mfem::ParMesh>(*mesh[l - 1]);
     // Mpi::Print("pmesh_copy sequence: {}, last operation == REFINE: {}, == REBALANCE: {}\n", pmesh_copy->GetSequence(), pmesh_copy->GetLastOperation() == mfem::Mesh::REFINE, pmesh_copy->GetLastOperation() == mfem::Mesh::REBALANCE);
      Mpi::Print("call fespace[{}] rebalance()\n", l - 1 - n_rebal - n_skipped);
      fespaces.GetFESpaceAtLevel(l - 1 - n_rebal - n_skipped).GetParMesh().Rebalance();//? This is bad, modifies the refined but not rebalanced mesh[l - 1]!
      Mpi::Print("get rebalance_op\n");
      auto rebalance_op = std::unique_ptr<Operator>(const_cast<Operator*>(fespaces.GetFESpaceAtLevel(l - 1 - n_rebal - n_skipped).Get().GetUpdateOperator()));
     // Mpi::Print("replace modified mesh by its copy\n");
     // //mesh[l - 1] = std::move(mesh_copy);
     // //mesh[l - 1] = std::make_unique<Mesh>(false, std::move(pmesh_copy));
     // mesh[l - 1] = std::make_unique<Mesh>(std::move(pmesh_copy), false);
     // Mpi::Print("mesh[l - 1] sequence: {}, last operation == REFINE: {}, == REBALANCE: {}\n", mesh[l - 1]->Get().GetSequence(), mesh[l - 1]->Get().GetLastOperation() == mfem::Mesh::REFINE, mesh[l - 1]->Get().GetLastOperation() == mfem::Mesh::REBALANCE);
     //Mpi::Print("Removing l-1 element of vector \n");
     // mesh.erase(mesh.begin() + l - 1); // error container overflow?
      //Mpi::Print("call mesh[l-1]->Update() ??\n");
      //mesh[l - 1]->Update();
      std::cout << "rank: " << Mpi::Rank(mesh[l]->GetComm()) << " refine width/height: " << refine_op->Width() << " " << refine_op->Height()
                << " rebalance width/height: " << rebalance_op->Width() << " " << rebalance_op->Height() << "\n";
      fespaces.UpdateLevel(
        std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get()),
        std::move(refine_op),
        std::move(rebalance_op)
      );
      prev_sequence++;
      Mpi::Print("After fespaces.UpdateLevel, num levels: {}\n", fespaces.GetNumLevels());
      n_rebal++;
      if (dbc_attr && dbc_tdof_lists)
      {
        dbc_tdof_lists->pop_back();
      }
    }
    else
    {
      prev_sequence++;
      fespaces.AddLevel(std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get()));
      Mpi::Print("fespaces added level {}: with mesh level {}\n", fespaces.GetNumLevels(), l);
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
