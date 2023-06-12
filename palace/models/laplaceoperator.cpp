// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "laplaceoperator.hpp"

#include "fem/coefficient.hpp"
#include "fem/multigrid.hpp"
#include "fem/operator.hpp"
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
  : dbc_marker(SetUpBoundaryProperties(iodata, *mesh.back())), skip_zeros(0),
    pc_gmg(iodata.solver.linear.mat_gmg), print_hdr(true),
    h1_fecs(utils::ConstructFECollections<mfem::H1_FECollection>(
        pc_gmg, false, iodata.solver.order, mesh.back()->Dimension())),
    nd_fec(iodata.solver.order, mesh.back()->Dimension()),
    h1_fespaces(
        pc_gmg
            ? utils::ConstructFiniteElementSpaceHierarchy(mesh, h1_fecs, dbc_marker)
            : utils::ConstructFiniteElementSpaceHierarchy(*mesh.back(), *h1_fecs.back())),
    nd_fespace(mesh.back().get(), &nd_fec), mat_op(iodata, *mesh.back()),
    source_attr_lists(ConstructSources(iodata))
{
  // Finalize setup.
  h1_fespaces.GetFinestFESpace().GetEssentialTrueDofs(dbc_marker, dbc_tdof_list);

  // Print essential BC information.
  if (dbc_marker.Max() > 0)
  {
    Mpi::Print("\nConfiguring Dirichlet BC at attributes:\n");
    utils::PrettyPrintMarker(dbc_marker);
  }
}

void LaplaceOperator::PrintHeader()
{
  if (print_hdr)
  {
    Mpi::Print("\nConfiguring system matrices, number of global unknowns: {:d}\n",
               h1_fespaces.GetFinestFESpace().GlobalTrueVSize());
    print_hdr = false;
  }
}

void LaplaceOperator::GetStiffnessMatrix(std::vector<std::unique_ptr<mfem::Operator>> &K,
                                         std::vector<std::unique_ptr<mfem::Operator>> &Ke)
{
  K.clear();
  Ke.clear();
  K.reserve(h1_fespaces.GetNumLevels());
  Ke.reserve(h1_fespaces.GetNumLevels());
  for (int l = 0; l < h1_fespaces.GetNumLevels(); l++)
  {
    auto &h1_fespace_l = h1_fespaces.GetFESpaceAtLevel(l);
    mfem::Array<int> dbc_tdof_list_l;
    h1_fespace_l.GetEssentialTrueDofs(dbc_marker, dbc_tdof_list_l);

    MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_REAL> epsilon_func(
        mat_op);
    mfem::ParBilinearForm k(&h1_fespace_l);
    k.AddDomainIntegrator(new mfem::MixedGradGradIntegrator(epsilon_func));
    // k.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    k.Assemble(skip_zeros);
    k.Finalize(skip_zeros);
    mfem::HypreParMatrix *hK = k.ParallelAssemble();
    mfem::HypreParMatrix *hKe = hK->EliminateRowsCols(dbc_tdof_list_l);
    PrintHeader();
    {
      std::string str = "";
      if (pc_gmg)
      {
        str =
            fmt::format(" (Level {:d}, {:d} unknowns)", l, h1_fespace_l.GlobalTrueVSize());
      }
      Mpi::Print(" K{}: NNZ = {:d}, norm = {:e}\n", str, hK->NNZ(),
                 hypre_ParCSRMatrixFnorm(*hK));
    }
    K.emplace_back(hK);
    Ke.emplace_back(hKe);
  }
}

std::unique_ptr<mfem::Operator> LaplaceOperator::GetNegGradMatrix()
{
  mfem::ParDiscreteLinearOperator grad(&h1_fespaces.GetFinestFESpace(), &nd_fespace);
  grad.AddDomainInterpolator(new mfem::GradientInterpolator);
  // grad.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  grad.Assemble();
  grad.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> NegGrad(grad.ParallelAssemble());
  *NegGrad *= -1.0;
  return NegGrad;
}

void LaplaceOperator::GetExcitationVector(int idx, const mfem::Operator &K,
                                          const mfem::Operator &Ke, mfem::Vector &X,
                                          mfem::Vector &RHS)
{
  // Apply the Dirichlet BCs to the solution vector: V = 1 on terminal boundaries with the
  // given index, V = 0 on all ground and other terminal boundaries.
  mfem::ParGridFunction x(&h1_fespaces.GetFinestFESpace());
  x = 0.0;

  // Get a marker of all boundary attributes with the given source surface index.
  mfem::Array<int> source_marker;
  const mfem::Array<int> &source_list = source_attr_lists[idx];
  mesh::AttrToMarker(dbc_marker.Size(), source_list, source_marker);
  mfem::ConstantCoefficient one(1.0);
  x.ProjectBdrCoefficient(one, source_marker);  // Values are only correct on master

  // Eliminate the essential BC to get the RHS vector.
  X.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
  RHS.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
  X = 0.0;
  RHS = 0.0;
  x.ParallelProject(X);  // Restrict to the true dofs
  dynamic_cast<const mfem::HypreParMatrix &>(K).EliminateBC(
      dynamic_cast<const mfem::HypreParMatrix &>(Ke), dbc_tdof_list, X, RHS);
}

}  // namespace palace
