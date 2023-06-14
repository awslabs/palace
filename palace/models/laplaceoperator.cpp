// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "laplaceoperator.hpp"

#include "fem/coefficient.hpp"
#include "fem/multigrid.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

namespace
{

mfem::Array<int> SetUpBoundaryProperties(const IoData &iodata, const mfem::ParMesh &mesh)
{
  int bdr_attr_max = mesh.bdr_attributes.Max();
  if (!iodata.boundaries.pec.empty() || !iodata.boundaries.lumpedport.empty())
  {
    // Check that boundary attributes have been specified correctly.
    mfem::Array<int> bdr_attr_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    bool first = true;
    for (auto attr : iodata.boundaries.pec.attributes)
    {
      // MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
      //             "Ground boundary attribute tags must be non-negative and correspond to
      //             " attributes in the mesh!");
      // MFEM_VERIFY(bdr_attr_marker[attr-1],
      //             "Unknown ground boundary attribute " << attr << "!");
      if (attr <= 0 || attr > bdr_attr_marker.Size() || !bdr_attr_marker[attr - 1])
      {
        if (first)
        {
          Mpi::Print("\n");
          first = false;
        }
        Mpi::Warning(
            "Unknown ground boundary attribute {:d}!\nSolver will just ignore it!\n", attr);
      }
    }
    for (const auto &[idx, data] : iodata.boundaries.lumpedport)
    {
      for (const auto &node : data.nodes)
      {
        for (auto attr : node.attributes)
        {
          MFEM_VERIFY(
              attr > 0 && attr <= bdr_attr_max,
              "Terminal boundary attribute tags must be non-negative and correspond to "
              "attributes in the mesh!");
          MFEM_VERIFY(bdr_attr_marker[attr - 1] > 0,
                      "Unknown terminal boundary attribute " << attr << "!");
        }
      }
    }
  }

  // Mark selected boundary attributes from the mesh as essential (Dirichlet).
  mfem::Array<int> dbc_bcs, dbc_marker;
  for (auto attr : iodata.boundaries.pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max)
    {
      continue;  // Can just ignore if wrong
    }
    dbc_bcs.Append(attr);
  }
  for (const auto &[idx, data] : iodata.boundaries.lumpedport)
  {
    for (const auto &node : data.nodes)
    {
      for (auto attr : node.attributes)
      {
        dbc_bcs.Append(attr);
      }
    }
  }
  MFEM_VERIFY(dbc_bcs.Size() > 0,
              "Electrostatic problem is ill-posed without any Dirichlet boundaries!");
  mesh::AttrToMarker(bdr_attr_max, dbc_bcs, dbc_marker);
  return dbc_marker;
}

std::map<int, mfem::Array<int>> ConstructSources(const IoData &iodata)
{
  // Construct mapping from terminal index to list of associated attributes.
  std::map<int, mfem::Array<int>> source_attr_lists;
  for (const auto &[idx, data] : iodata.boundaries.lumpedport)
  {
    mfem::Array<int> &attr_list = source_attr_lists[idx];
    for (const auto &node : data.nodes)
    {
      for (auto attr : node.attributes)
      {
        attr_list.Append(attr);
      }
    }
  }
  return source_attr_lists;
}

}  // namespace

LaplaceOperator::LaplaceOperator(const IoData &iodata,
                                 const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh)
  : assembly_level(iodata.solver.linear.mat_pa ? mfem::AssemblyLevel::PARTIAL
                                               : mfem::AssemblyLevel::LEGACY),
    skip_zeros(0), pc_mg(iodata.solver.linear.pc_mg), print_hdr(true),
    dbc_marker(SetUpBoundaryProperties(iodata, *mesh.back())),
    h1_fecs(utils::ConstructFECollections<mfem::H1_FECollection>(
        pc_mg, false, iodata.solver.order, mesh.back()->Dimension())),
    nd_fec(iodata.solver.order, mesh.back()->Dimension()),
    h1_fespaces(pc_mg ? utils::ConstructFiniteElementSpaceHierarchy(
                            mesh, h1_fecs, &dbc_marker, &dbc_tdof_lists)
                      : utils::ConstructFiniteElementSpaceHierarchy(
                            *mesh.back(), *h1_fecs.back(), &dbc_marker,
                            &dbc_tdof_lists.emplace_back())),
    nd_fespace(mesh.back().get(), &nd_fec), mat_op(iodata, *mesh.back()),
    source_attr_lists(ConstructSources(iodata))
{
  // Print essential BC information.
  if (dbc_marker.Max() > 0)
  {
    Mpi::Print("\nConfiguring Dirichlet BC at attributes:\n");
    utils::PrettyPrintMarker(dbc_marker);
  }
}

std::unique_ptr<Operator> LaplaceOperator::GetStiffnessMatrix()
{
  if (print_hdr)
  {
    Mpi::Print("\nAssembling system matrices, number of global unknowns:\n"
               " H1: {:d}, ND: {:d}\n",
               GetH1Space().GlobalTrueVSize(), GetNDSpace().GlobalTrueVSize());
    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
  auto K = std::make_unique<MultigridOperator>(h1_fespaces.GetNumLevels());
  for (int l = 0; l < h1_fespaces.GetNumLevels(); l++)
  {
    auto &h1_fespace_l = h1_fespaces.GetFESpaceAtLevel(l);
    constexpr MaterialPropertyType MatType = MaterialPropertyType::PERMITTIVITY_REAL;
    MaterialPropertyCoefficient<MatType> epsilon_func(mat_op);
    auto k = std::make_unique<mfem::SymmetricBilinearForm>(&h1_fespace_l);
    k->AddDomainIntegrator(new mfem::MixedGradGradIntegrator(epsilon_func));
    k->SetAssemblyLevel(assembly_level);
    k->Assemble(skip_zeros);
    k->Finalize(skip_zeros);
    if (print_hdr)
    {
      Mpi::Print(" Level {:d}: {:d} unknowns", l, h1_fespace_l.GlobalTrueVSize());
      if (assembly_level == mfem::AssemblyLevel::LEGACY)
      {
        HYPRE_BigInt nnz = k->SpMat().NumNonZeroElems();
        Mpi::GlobalSum(1, &nnz, h1_fespace_l.GetComm());
        Mpi::Print(", {:d} NNZ\n", nnz);
      }
      else
      {
        Mpi::Print("\n");
      }
    }
    auto K_l = std::make_unique<ParOperator>(std::move(k), h1_fespace_l);
    K_l->SetEssentialTrueDofs(dbc_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
    K->AddOperator(std::move(K_l));
  }
  print_hdr = false;
  return K;
}

std::unique_ptr<Operator> LaplaceOperator::GetGradMatrix()
{
  auto grad = std::make_unique<mfem::DiscreteLinearOperator>(&GetH1Space(), &GetNDSpace());
  grad->AddDomainInterpolator(new mfem::GradientInterpolator);
  grad->SetAssemblyLevel(assembly_level);
  grad->Assemble();
  grad->Finalize();
  return std::make_unique<ParOperator>(std::move(grad), GetH1Space(), GetNDSpace(), true);
}

void LaplaceOperator::GetExcitationVector(int idx, const Operator &K, Vector &X,
                                          Vector &RHS)
{
  // Apply the Dirichlet BCs to the solution vector: V = 1 on terminal boundaries with the
  // given index, V = 0 on all ground and other terminal boundaries.
  mfem::ParGridFunction x(&GetH1Space());
  x = 0.0;

  // Get a marker of all boundary attributes with the given source surface index.
  mfem::Array<int> source_marker;
  const mfem::Array<int> &source_list = source_attr_lists[idx];
  mesh::AttrToMarker(dbc_marker.Size(), source_list, source_marker);
  mfem::ConstantCoefficient one(1.0);
  x.ProjectBdrCoefficient(one, source_marker);  // Values are only correct on master

  // Eliminate the essential BC to get the RHS vector.
  X.SetSize(GetH1Space().GetTrueVSize());
  RHS.SetSize(GetH1Space().GetTrueVSize());
  X = 0.0;
  RHS = 0.0;
  x.ParallelProject(X);  // Restrict to the true dofs
  const auto *mg_K = dynamic_cast<const MultigridOperator *>(&K);
  const auto *PtAP_K = mg_K ? dynamic_cast<const ParOperator *>(&mg_K->GetFinestOperator())
                            : dynamic_cast<const ParOperator *>(&K);
  MFEM_VERIFY(PtAP_K, "LaplaceOperator requires ParOperator for RHS elimination!");
  PtAP_K->EliminateRHS(X, RHS);
}

}  // namespace palace
