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
    int mg_max_levels, /*const*/ std::vector<std::unique_ptr<Mesh>> &mesh, int call_id,
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

  // n = 6
  // mesh vector has [orig, ref1 x n, reb1, ref2 x n, reb2, ...]
  // first call to ConstructFiniteElementSpaceHierarchy will do
  // 0   orig fe0
  // 1   ref1 fe1
  // 2   ref1 fe1_copy
  // 3   ref1 skip
  // 4.  ref1 skip
  // 5   ref1 skip
  // 6   ref1 skip
  // 7   reb1 fe1 (modifies #1) refine between fe0 and fe1_copy, rebalance
  // 8.  ref2 fe2
  // 9.  ref2 fe2_copy
  // 10  ref2 skip
  // ...
  // 13 ref2 skip
  // 14 reb2 fe2 (modifies #n+2) refine between fe1_copy and fe2_copy
  // 15 ref2 fe3
  // ..

  // second call to ConstructFiniteElementSpaceHierarchy will do
  // 0   orig fe0
  // 1.  ref1 skip
  // 2.  ref1 skip
  // 3.  ref1 fe1
  // 4.  ref1 fe1_copy
  // 5.  ref1 skip
  // 6.  ref1 skip
  // 7.  reb1 fe1 (modifies #3)
  // 8.  ref2 skip
  // 9   ref2 skip
  // 10. ref2 fe2
  // 11  ref2 fe2_copy
  // 12. ref2 skip
  // 13  ref2 skip
  // 14. reb2 fe2 (modifies #10)

  std::size_t n_copies = 6;
  /*
  for (std::size_t c = 0; c < 3; c++)
  {
    Mpi::Print("call # {} to ConstructFiniteElementSpaceHierarchy\n", c);
    std::size_t fespace_idx = 1;
    for (std::size_t l = 1; l < 29; l++)
    {
      Mpi::Print("Mesh level {}, ", l);
      if (l == (fespace_idx - 1) * (n_copies + 1) + 1 + c*2) {Mpi::Print("Creating fespace {}\n", fespace_idx);}
      else if (l == (fespace_idx - 1) * (n_copies + 1) + 2 + c*2) {Mpi::Print("Creating copy of fespace {}\n", fespace_idx);}
      else if (l == fespace_idx * (n_copies + 1)) {Mpi::Print("Rebalancing and updating fespace {}\n", fespace_idx); fespace_idx++;}
      else {Mpi::Print("skipping\n");}
    }
  }
  */
  //int n_rebal = 0;
  //int n_skipped = 0;

  // h-refinement.
  std::size_t fe_idx = 1;
  for (std::size_t l = coarse_mesh_l + 1; l < mesh.size(); l++)
  {
    Mpi::Print("Mesh level {}, ", l);
    if (l == (fe_idx - 1) * (n_copies + 1) + 1 + call_id*2)
    {
      fespaces.AddLevel(std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get()));
      Mpi::Print("fespaces added level {}\n", fe_idx);
    }
    else if (l == (fe_idx - 1) * (n_copies + 1) + 2 + call_id*2)
    {
      fespaces.CopyLevel(std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get()));
      Mpi::Print("fespaces created copy of level {}\n", fe_idx);
      continue;
    }
    else if (l == fe_idx * (n_copies + 1))
    {
      Mpi::Print("Verifying that mesh level {} is a rebalance: {}, ", l, mesh[l]->Get().GetLastOperation() == mfem::Mesh::REBALANCE);
      Mpi::Print("rebalancing and updating fespace {}\n", fe_idx);
      auto refine_op = std::make_unique<mfem::TransferOperator>(fespaces.GetFESpaceAtLevel(fe_idx - 1, true), fespaces.GetFESpaceAtLevel(fe_idx, true));
      Vector t1(refine_op->Width()); Vector t2(refine_op->Height());
      t1 = 1.0; t2 = 0.0;
      refine_op->Mult(t1, t2);
      Mpi::Print("before rebalance t2 min/max: {}, {}\n", t2.Min(), t2.Max());
      //Mpi::Print("call fespace[{}] rebalance()\n", l - 1 - n_rebal - n_skipped);
      fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().Rebalance();
      //^^^is bad, modifies the refined but not rebalanced mesh[l - 1] and probably messes up refine_op which depends on fespaces[l - 1 - n_rebal - n_skipped]
      t1 = 1.0; t2 = 0.0;
      refine_op->Mult(t1, t2);
      Mpi::Print("after rebalance t2 min/max: {}, {}\n", t2.Min(), t2.Max());
      Mpi::Print("get rebalance_op\n");
      auto rebalance_op = std::unique_ptr<Operator>(const_cast<Operator*>(fespaces.GetFESpaceAtLevel(fe_idx).Get().GetUpdateOperator()));

      std::cout << "rank: " << Mpi::Rank(mesh[l]->GetComm()) << " refine width/height: " << refine_op->Width() << " " << refine_op->Height()
                << " rebalance width/height: " << rebalance_op->Width() << " " << rebalance_op->Height() << "\n";
      fespaces.UpdateLevel(
        std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get()),
        std::move(refine_op),
        std::move(rebalance_op)
      );
      fe_idx++;
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


/*
    if ((mesh[l]->Get().GetLastOperation() == mfem::Mesh::REFINE) && (mesh[l - 1]->Get().GetGlobalNE() == mesh[l]->Get().GetGlobalNE()))
    {
      // Duplicate mesh, skip it!
      Mpi::Print("Skipping duplicate mesh\n");
      n_skipped++;
      continue;
    }
    else if (mesh[l]->Get().GetLastOperation() == mfem::Mesh::REBALANCE)
    {
      // Rebalance
      Mpi::Print("Mesh level {} is a rebalancing (n_rebal: {}), compute transfer operator\n", l, n_rebal);
      auto refine_op = std::make_unique<mfem::TransferOperator>(fespaces.GetFESpaceAtLevel(l - 2 - n_rebal - n_skipped), fespaces.GetFESpaceAtLevel(l - 1 - n_rebal - n_skipped));
      Vector t1(refine_op->Width()); Vector t2(refine_op->Height());
      t1 = 1.0; t2 = 0.0;
      refine_op->Mult(t1, t2);
      Mpi::Print("before rebalance t2 min/max: {}, {}\n", t2.Min(), t2.Max());
      Mpi::Print("call fespace[{}] rebalance()\n", l - 1 - n_rebal - n_skipped);
      fespaces.GetFESpaceAtLevel(l - 1 - n_rebal - n_skipped).GetParMesh().Rebalance();
      //^^^is bad, modifies the refined but not rebalanced mesh[l - 1] and probably messes up refine_op which depends on fespaces[l - 1 - n_rebal - n_skipped]
      t1 = 1.0; t2 = 0.0;
      refine_op->Mult(t1, t2);
      Mpi::Print("after rebalance t2 min/max: {}, {}\n", t2.Min(), t2.Max());
      Mpi::Print("get rebalance_op\n");
      auto rebalance_op = std::unique_ptr<Operator>(const_cast<Operator*>(fespaces.GetFESpaceAtLevel(l - 1 - n_rebal - n_skipped).Get().GetUpdateOperator()));

      std::cout << "rank: " << Mpi::Rank(mesh[l]->GetComm()) << " refine width/height: " << refine_op->Width() << " " << refine_op->Height()
                << " rebalance width/height: " << rebalance_op->Width() << " " << rebalance_op->Height() << "\n";
      fespaces.UpdateLevel(
        std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get()),
        std::move(refine_op),
        std::move(rebalance_op)
      );

      Mpi::Print("After fespaces.UpdateLevel, num levels: {}\n", fespaces.GetNumLevels());
      n_rebal++;
      if (dbc_attr && dbc_tdof_lists)
      {
        dbc_tdof_lists->pop_back();
      }
    }
    else
    {
      fespaces.AddLevel(std::make_unique<FiniteElementSpace>(*mesh[l], fecs[0].get()));
      Mpi::Print("fespaces added level {}: with mesh level {}\n", fespaces.GetNumLevels(), l);
    }
*/
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
