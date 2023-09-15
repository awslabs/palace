// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "laplaceoperator.hpp"

#include "fem/bilinearform.hpp"
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
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
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
      for (const auto &elem : data.elements)
      {
        for (auto attr : elem.attributes)
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
    for (const auto &elem : data.elements)
    {
      for (auto attr : elem.attributes)
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
    for (const auto &elem : data.elements)
    {
      for (auto attr : elem.attributes)
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
  : pa_order_threshold(iodata.solver.pa_order_threshold),
    pa_discrete_interp(iodata.solver.pa_discrete_interp), skip_zeros(false),
    print_hdr(true), dbc_marker(SetUpBoundaryProperties(iodata, *mesh.back())),
    h1_fecs(fem::ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsen_type, false)),
    nd_fec(iodata.solver.order, mesh.back()->Dimension()),
    h1_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        iodata.solver.linear.mg_max_levels, iodata.solver.linear.mg_legacy_transfer,
        pa_order_threshold, pa_discrete_interp, mesh, h1_fecs, &dbc_marker,
        &dbc_tdof_lists)),
    nd_fespace(mesh.back().get(), &nd_fec), mat_op(iodata, *mesh.back()),
    source_attr_lists(ConstructSources(iodata))
{
  // Print essential BC information.
  if (dbc_marker.Size() && dbc_marker.Max() > 0)
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
               " H1: {:d}, ND: {:d}\n Operator assembly level: {}\n",
               GetH1Space().GlobalTrueVSize(), GetNDSpace().GlobalTrueVSize(),
               GetH1Space().GetMaxElementOrder() > pa_order_threshold ? "Partial" : "Full");
    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
  auto K = std::make_unique<MultigridOperator>(h1_fespaces.GetNumLevels());
  for (int l = 0; l < h1_fespaces.GetNumLevels(); l++)
  {
    // Force coarse level operator to be fully assembled always.
    const auto &h1_fespace_l = h1_fespaces.GetFESpaceAtLevel(l);
    if (print_hdr)
    {
      Mpi::Print(" Level {:d} (p = {:d}): {:d} unknowns", l,
                 h1_fespace_l.GetMaxElementOrder(), h1_fespace_l.GlobalTrueVSize());
    }
    constexpr auto MatType = MaterialPropertyType::PERMITTIVITY_REAL;
    MaterialPropertyCoefficient<MatType> epsilon_func(mat_op);
    BilinearForm k(h1_fespace_l);
    k.AddDomainIntegrator(std::make_unique<DiffusionIntegrator>(epsilon_func));
    auto K_l = std::make_unique<ParOperator>(
        k.Assemble((l > 0) ? pa_order_threshold : 99, skip_zeros), h1_fespace_l);
    if (print_hdr)
    {
      if (const auto *k_spm =
              dynamic_cast<const mfem::SparseMatrix *>(&K_l->LocalOperator()))
      {
        HYPRE_BigInt nnz = k_spm->NumNonZeroElems();
        Mpi::GlobalSum(1, &nnz, h1_fespace_l.GetComm());
        Mpi::Print(", {:d} NNZ\n", nnz);
      }
      else
      {
        Mpi::Print("\n");
      }
    }
    K_l->SetEssentialTrueDofs(dbc_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
    K->AddOperator(std::move(K_l));
  }
  print_hdr = false;
  return K;
}

std::unique_ptr<Operator> LaplaceOperator::GetGradMatrix()
{
  constexpr bool skip_zeros_interp = true;
  DiscreteLinearOperator grad(GetH1Space(), GetNDSpace());
  grad.AddDomainInterpolator(std::make_unique<GradientInterpolator>());
  return std::make_unique<ParOperator>(
      grad.Assemble(pa_discrete_interp ? pa_order_threshold : 99, skip_zeros_interp),
      GetH1Space(), GetNDSpace(), true);
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
