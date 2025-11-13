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
    int mg_max_levels, /*const*/ std::vector<std::unique_ptr<Mesh>> &mesh, int call_id, int amr_max_its,
    const std::vector<std::unique_ptr<FECollection>> &fecs,
    const mfem::Array<int> *dbc_attr = nullptr,
    std::vector<mfem::Array<int>> *dbc_tdof_lists = nullptr)
{
  MFEM_VERIFY(!mesh.empty() && !fecs.empty() &&
                  (!dbc_tdof_lists || dbc_tdof_lists->empty()),
              "Empty mesh or FE collection for FE space construction!");
  int coarse_mesh_l = std::max(0, static_cast<int>(mesh.size() + fecs.size()) - 1 -
                                      std::max(1, mg_max_levels)); // CAREFUL THIS NEEDS TO BE CHANGED SINCE mesh CONTAINS A BUNCH OF DUPLICATES!!!!!
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

  // mesh vector has [orig, ref1 * 2 * max_its + 1, reb1, ref2 * 2 * (max_its - 1) + 1, reb2, ...]
  // Can infer the AMR iteration number from the size of the mesh vector.
  // it = 1 -> 1 + 2 * max_its + 2         = 2 * max_its + 3 = 9
  // it = 2 -> it1 + 2 * (max_its - 1) + 2 = 4 * max_its + 3 = 15
  // it = 3 -> it2 + 2 * (max_its - 2) + 2 = 6 * max_its + 1 = 19
  std::size_t mesh_size = mesh.size();
  std::size_t amr_it = (mesh_size == (2 * amr_max_its + 3)) ? 1 : (mesh_size == (4 * amr_max_its + 3)) ? 2 : (mesh_size == (6 * amr_max_its + 1)) ? 3 : 0;
  if (amr_it == 0) Mpi::Print("\n\n\n COULD NOT DETECT AMR IT!!\n\n");
  else Mpi::Print("AMR It # {}\n", amr_it);
  // AMR it = 3 (assuming max_its = 3)
  // 0 orig fe0
  // 1 ref1 fe1_copy
  // 2 ref1 skip (modified by it=1, call=0)
  // 3 ref1 skip (modified by it=1, call=1)
  // 4 ref1 skip (modified by it=2 call=0)
  // 5 ref1 skip (modified by it=2 call=1)
  // 6 ref1 fe1 (modified by it=3, call=0)
  // 7 ref1 fe1 (modified by it=3, call=1)
  // 8 reb1 fe1 (ref transfer between fe0 and fe1_copy + rebalance modifies #6 for call=0, #7 for call=1)
  // 9 ref2 fe2_copy
  // 10 ref2 skip (modified by it=2 call=0)
  // 11 ref2 skip (modified by it=2 call=1)
  // 12 ref2 fe2 (modified by it=3 call=0)
  // 13 ref2 fe2 (modified by it=3 call=1)
  // 14 reb2 (refine transfer between fe1 and fe2_copy modifies #13 and #14)
  // 15 ref3 fe3_copy
  // 16 ref3 fe3 (modified by it=3 call=0)
  // 17 ref3 fe3 (modified by it=3 call=1)
  // 18 reb3 fe3

  // h-refinement.
  std::size_t fe_idx = 1, copy_idx = 1, copy_increment = 2 * (amr_max_its - (fe_idx - 1)) + 2;
  for (std::size_t l = coarse_mesh_l + 1; l < mesh.size(); l++)
  {
    Mpi::Print("Mesh level {}, ", l);
    if (l == copy_idx + 2 * (amr_it - fe_idx) + call_id + 1)
    {
      fespaces.AddLevel(std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get()));
      Mpi::Print("fespaces added fespace level {}\n", fe_idx);
    }
    else if (l == copy_idx)
    {
      fespaces.CopyLevel(std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get()));
      Mpi::Print("fespaces created copy of fespace level {}\n", fe_idx);
      continue;
    }
    else if (l == copy_idx + copy_increment - 1)
    {
      Mpi::Print("Verifying that mesh level {} is a rebalance: {}, ", l, mesh[l]->Get().GetLastOperation() == mfem::Mesh::REBALANCE);
      Mpi::Print("rebalancing and updating fespace {}\n", fe_idx);
      auto refine_op = std::make_unique<mfem::TransferOperator>(fespaces.GetFESpaceAtLevel(fe_idx - 1), fespaces.GetFESpaceAtLevel(fe_idx, true));
      Vector t1(refine_op->Width()); Vector t2(refine_op->Height());
      t1 = 1.0; t2 = 0.0;
      refine_op->Mult(t1, t2); // refine_op is wrong at second AMR it between levels 1 and 2
      Mpi::Print("before rebalance t2 min/max: {}, {}\n", t2.Min(), t2.Max());
      Mpi::Print("GetNumLevels: {}\n", fespaces.GetNumLevels());
      Mpi::Print("fespaces.GetFESpaceAtLevel({}).GetParMesh().GetLastOperation() == REFINE: {}\n", fe_idx, fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().GetLastOperation() == mfem::Mesh::REFINE);
      fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().Rebalance(); // modifies the mesh associated with that fespace!
      Mpi::Print("mesh {} was modified by rebalance_op\n", copy_idx + 2 * (amr_it - fe_idx) + call_id + 1);
      t1 = 1.0; t2 = 0.0;
      //refine_op->Mult(t1, t2);
      //Mpi::Print("after rebalance t2 min/max: {}, {}\n", t2.Min(), t2.Max());
      Mpi::Print("get rebalance_op\n");
      auto rebalance_op = std::unique_ptr<Operator>(const_cast<Operator*>(fespaces.GetFESpaceAtLevel(fe_idx).Get().GetUpdateOperator()));
      Mpi::Print("SetUpdateOperatorOwner to false\n");
      fespaces.GetFESpaceAtLevel(fe_idx).Get().SetUpdateOperatorOwner(false); // test?
      std::cout << "rank: " << Mpi::Rank(mesh[l]->GetComm()) << " refine width/height: " << refine_op->Width() << " " << refine_op->Height()
                << " rebalance width/height: " << rebalance_op->Width() << " " << rebalance_op->Height() << "\n";
      fespaces.UpdateLevel(
        std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get()),
        std::move(refine_op),
        std::move(rebalance_op)
      );
      fe_idx++;
      copy_idx += copy_increment;
      copy_increment = 2 * (amr_max_its - (fe_idx - 1)) + 2;
      if (dbc_attr && dbc_tdof_lists)
      {
        dbc_tdof_lists->pop_back();
      }
    }
    else
    {
      Mpi::Print("skipping\n");
      continue;
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
