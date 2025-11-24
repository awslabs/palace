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
  Mpi::Print("In ConstructFiniteElementSpaceHierarchy with mg_max_levels: {}\n", mg_max_levels);

  // mesh vector has [orig, ref1 * 2 * max_its + 1, reb1, ref2 * 2 * (max_its - 1) + 1, reb2, ...]
  // Can infer the AMR iteration number from the size of the mesh vector.
  // it1 = 1 + 2 * max_its + 2         = 2 * max_its + 3
  // it2 = it1 + 2 * (max_its - 1) + 2 = 4 * max_its + 3
  // it3 = it2 + 2 * (max_its - 2) + 2 = 6 * max_its + 1
  // it4 = it3 + 2 * (max_its - 3) + 2 = 8 * max_its - 3
  // it5 = it4 + 2 * (max_its - 4) + 2 = 10* max_its - 9
  std::size_t mesh_size = mesh.size();
  //std::size_t amr_it = (mesh_size == (2 * amr_max_its + 3)) ? 1 : (mesh_size == (4 * amr_max_its + 3)) ? 2 : (mesh_size == (6 * amr_max_its + 1)) ? 3 : 0;
  std::size_t amr_it = int( ((2.0*amr_max_its + 3.0) - std::sqrt(std::pow(2.0*amr_max_its + 3.0, 2.0) - 4.0*(mesh_size - 1.0))) / 2.0);
  if (amr_it < 0 || amr_it > amr_max_its) Mpi::Print("\n\n\n COULD NOT DETECT AMR IT!!\n\n");
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
      mfem::ParGridFunction g1(&fespaces.GetFESpaceAtLevel(fe_idx - 1).Get()); //g1 = 0.0;//g1 = fespaces.GetFESpaceAtLevel(fe_idx - 1).GetParMesh().GetMyRank() + 1;
      mfem::ParGridFunction g2(&fespaces.GetFESpaceAtLevel(fe_idx, true).Get()); //g2 = 0.0;
      mfem::ParGridFunction g2_2(&fespaces.GetFESpaceAtLevel(fe_idx).Get()); //g2_2 = 0.0;
      //mfem::ParGridFunction g2_comp(&fespaces.GetFESpaceAtLevel(fe_idx, true).Get()); g2_comp = fespaces.GetFESpaceAtLevel(fe_idx, true).GetParMesh().GetMyRank() + 1;
      mfem::VectorFunctionCoefficient vcoef(3, [](const Vector &x, Vector &v) {
         v(0) = 1.0; v(1) = 2.0; v(2) = 3.0;
      });
      mfem::FunctionCoefficient scoef([](const Vector &x) {
         return x(0) + 2 * x(1) + 3 * x(2);
      });
      if (fecs[0]->GetRangeType(3) == mfem::FiniteElement::VECTOR)
      {
        Mpi::Print("Vector FiniteElementSpace!\n");
        g1.ProjectCoefficient(vcoef);
        g1.ExchangeFaceNbrData();
        //std::cout << "rank: " << fespaces.GetFESpaceAtLevel(fe_idx - 1).GetParMesh().GetMyRank() << " g1.Size(): " << g1.Size() << " g1 min/max: " << g1.Min() << " " << g1.Max() << " g1 max error: " << g1.ComputeMaxError(vcoef) << "\n";
        Mpi::Print("g1 max error: {}\n", g1.ComputeMaxError(vcoef));
        refine_op->Mult(g1, g2);
        g2.ExchangeFaceNbrData();
        //std::cout << "rank: " << fespaces.GetFESpaceAtLevel(fe_idx, true).GetParMesh().GetMyRank() << " g2.Size(): " << g2.Size() << " rebal g2 min/max: " << g2.Min() << " " << g2.Max() << " g2 max error: " << g2.ComputeMaxError(vcoef) << "\n";
        Mpi::Print("g2 max error: {}\n", g2.ComputeMaxError(vcoef));
        g2_2.ProjectCoefficient(vcoef);
        Mpi::Print("g2_2 max error before rebalance: {}\n", g2_2.ComputeMaxError(vcoef));
      }
      else
      {
        Mpi::Print("Scalar FiniteElementSpace!\n");
        g1.ProjectCoefficient(scoef);
        g1.ExchangeFaceNbrData();
        //std::cout << "rank: " << fespaces.GetFESpaceAtLevel(fe_idx - 1).GetParMesh().GetMyRank() << " g1.Size(): " << g1.Size() << " g1 min/max: " << g1.Min() << " " << g1.Max() << " g1 max error: " << g1.ComputeMaxError(scoef) << "\n";
        Mpi::Print("g1 max error: {}\n", g1.ComputeMaxError(scoef));
        refine_op->Mult(g1, g2);
        g2.ExchangeFaceNbrData();
        //std::cout << "rank: " << fespaces.GetFESpaceAtLevel(fe_idx, true).GetParMesh().GetMyRank() << " g2.Size(): " << g2.Size() << " rebal g2 min/max: " << g2.Min() << " " << g2.Max() << " g2 max error: " << g2.ComputeMaxError(scoef) << "\n";
        Mpi::Print("g2 max error: {}\n", g2.ComputeMaxError(scoef));
        g2_2.ProjectCoefficient(scoef); Mpi::Print("g2_2 max error before rebalance: {}\n", g2_2.ComputeMaxError(scoef));
      }
      //g2_comp2 -= g2_comp;
      //g2_comp -= g2;
      //std::cout << "rank: " << fespaces.GetFESpaceAtLevel(fe_idx, true).GetParMesh().GetMyRank() << " g2.Size(): " << g2.Size() << " rebal g2 min/max: " << g2.Min() << " " << g2.Max() << " g2_comp min/max: " << g2_comp.Min() << " " << g2_comp.Max() << "\n";

      //std::cout << "rank: " << fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().GetMyRank() << " g2_comp2 min/max: " << g2_comp2.Min() << " " << g2_comp2.Max() << "\n";
      Mpi::Print("GetNumLevels: {}\n", fespaces.GetNumLevels());
      Mpi::Print("fespaces.GetFESpaceAtLevel({}).GetParMesh().GetLastOperation() == REFINE: {}\n", fe_idx, fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().GetLastOperation() == mfem::Mesh::REFINE);
      fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().Rebalance(); // modifies the mesh associated with that fespace!
      Mpi::Print("mesh {} was modified by rebalance_op\n", copy_idx + 2 * (amr_it - fe_idx) + call_id + 1);
      Mpi::Print("get rebalance_op\n");
      auto rebalance_op = std::unique_ptr<Operator>(const_cast<Operator*>(fespaces.GetFESpaceAtLevel(fe_idx).Get().GetUpdateOperator()));
      // Mpi::Print("SetUpdateOperatorOwner to false\n");
      fespaces.GetFESpaceAtLevel(fe_idx).Get().SetUpdateOperatorOwner(false); // test
      //mfem::OperatorPtr T(mfem::Operator::Hypre_ParCSR); //make sure Operator type is correct...
      //fespaces.GetFESpaceAtLevel(fe_idx).Get().Update(); // need to call Update() otherwise T.Ptr is null?!
      //fespaces.GetFESpaceAtLevel(fe_idx).Get().GetUpdateOperator(T);
      //mfem::HypreParMatrix *Thm = T.Is<mfem::HypreParMatrix>();
      //std::cout << "rank " << Mpi::Rank(mesh[l]->GetComm()) << " T.Ptr() != nullptr: " << (T.Ptr() != nullptr) <<  " Thm != nullptr: " << (Thm != nullptr) << "\n";
      //std::unique_ptr<Operator> rebalance_op;
      //T.SetOperatorOwner(false); // should it be before or after assigning T.Ptr() to rebalance_op?
      //rebalance_op.reset(T.Ptr());
      //rebalance_op.reset(Thm);
      std::cout << "rank: " << Mpi::Rank(mesh[l]->GetComm()) << " GetNE: " << mesh[l]->GetNE() << " refine width/height: " << refine_op->Width() << " " << refine_op->Height()
                << " rebalance width/height: " << rebalance_op->Width() << " " << rebalance_op->Height() << "\n";
      mfem::ParGridFunction g3(&fespaces.GetFESpaceAtLevel(fe_idx).Get()); g3 = 0.0;
      //mfem::ParGridFunction g3_comp(&fespaces.GetFESpaceAtLevel(fe_idx).Get()); g3_comp = fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().GetMyRank() + 1;
      //std::cout << "before rebalance_op->Mult rank: " << Mpi::Rank(mesh[l]->GetComm()) << "\n";
      rebalance_op->Mult(g2, g3);
      g3.ExchangeFaceNbrData();
      if (fecs[0]->GetRangeType(3) == mfem::FiniteElement::VECTOR)
      {
        g2_2.Update(); Mpi::Print("g2_2 max error after rebalance and update: {}\n", g2_2.ComputeMaxError(vcoef));
        //std::cout << "rank: " << fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().GetMyRank() << " after rebal g3 min/max: " << g3.Min() << " " << g3.Max() << " g3 max error: " << g3.ComputeMaxError(vcoef) << "\n";
        Mpi::Print("g3 max error: {}\n", g3.ComputeMaxError(vcoef));
      }
      else
      {
        g2_2.Update(); Mpi::Print("g2_2 max edrror after rebalance and update: {}\n", g2_2.ComputeMaxError(scoef));
        //std::cout << "rank: " << fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().GetMyRank() << " after rebal g3 min/max: " << g3.Min() << " " << g3.Max() << " g3 max error: " << g3.ComputeMaxError(scoef) << "\n";
        Mpi::Print("g3 max error: {}\n", g1.ComputeMaxError(scoef));
      }
      //std::cout << "after rebalance_op->Mult rank: " << Mpi::Rank(mesh[l]->GetComm()) << "\n";
      //g3_comp -= g3;
      //std::cout << "rank: " << fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().GetMyRank() << " after rebal g3 min/max: " << g3.Min() << " " << g3.Max() << " g3_comp min/max: " << g3_comp.Min() << " " << g3_comp.Max() << "\n";

      Mpi::Print("Compare rebalanced meshes\n");
      // Compare element counts per rank
      int local_elements1 = mesh[l]->Get().GetNE();
      int local_elements2 = fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().GetNE();
      // Get global element numbering
      mfem::Array<int> global_elem_ids1, global_elem_ids2;
      mesh[l]->Get().GetGlobalElementIndices(global_elem_ids1);
      fespaces.GetFESpaceAtLevel(fe_idx).GetParMesh().GetGlobalElementIndices(global_elem_ids2);
      // Compare on each rank
      bool identical = (local_elements1 == local_elements2);
      if (identical && local_elements1 > 0) {
        for (int i = 0; i < local_elements1; i++) {
          if (global_elem_ids1[i] != global_elem_ids2[i]) {
              identical = false;
              break;
          }
        }
      }
      // Gather results across all ranks
      int local_match = identical ? 1 : 0;
      int global_match;
      MPI_Allreduce(&local_match, &global_match, 1, MPI_INT, MPI_LAND, mesh[l]->GetComm());
      Mpi::Print("Rebalanced meshes partionings are identical: {}\n", global_match);

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
    //fespaces.AddLevel(std::make_unique<FiniteElementSpace>(fespaces.GetFESpaceAtLevel(fe_idx-1).GetMesh(), fecs[l].get()));
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
