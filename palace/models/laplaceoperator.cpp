// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "laplaceoperator.hpp"

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/multigrid.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

LaplaceOperator::LaplaceOperator(const IoData &iodata,
                                 const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh)
  : print_hdr(true), dbc_attr(SetUpBoundaryProperties(iodata, *mesh.back())),
    h1_fecs(fem::ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsen_type, false)),
    nd_fec(std::make_unique<mfem::ND_FECollection>(iodata.solver.order,
                                                   mesh.back()->Dimension())),
    h1_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        iodata.solver.linear.mg_max_levels, mesh, h1_fecs, &dbc_attr, &dbc_tdof_lists)),
    nd_fespace(h1_fespaces.GetFinestFESpace(), mesh.back().get(), nd_fec.get()),
    mat_op(iodata, *mesh.back()), source_attr_lists(ConstructSources(iodata))
{
  // Finalize setup.
  BilinearForm::pa_order_threshold = iodata.solver.pa_order_threshold;
  fem::DefaultIntegrationOrder::q_order_jac = iodata.solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = iodata.solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = iodata.solver.q_order_extra;

  // Print essential BC information.
  if (dbc_attr.Size())
  {
    Mpi::Print("\nConfiguring Dirichlet BC at attributes:\n");
    std::sort(dbc_attr.begin(), dbc_attr.end());
    utils::PrettyPrint(dbc_attr);
  }
}

mfem::Array<int> LaplaceOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                          const mfem::ParMesh &mesh)
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
  mfem::Array<int> dbc_bcs;
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
  return dbc_bcs;
}

std::map<int, mfem::Array<int>> LaplaceOperator::ConstructSources(const IoData &iodata)
{
  // Construct mapping from terminal index to list of associated attributes.
  std::map<int, mfem::Array<int>> attr_lists;
  for (const auto &[idx, data] : iodata.boundaries.lumpedport)
  {
    mfem::Array<int> &attr_list = attr_lists[idx];
    for (const auto &elem : data.elements)
    {
      for (auto attr : elem.attributes)
      {
        attr_list.Append(attr);
      }
    }
  }
  return attr_lists;
}

namespace
{

void PrintHeader(const mfem::ParFiniteElementSpace &h1_fespace,
                 const mfem::ParFiniteElementSpace &nd_fespace, bool &print_hdr)
{
  if (print_hdr)
  {
    Mpi::Print("\nAssembling system matrices, number of global unknowns:\n"
               " H1 (p = {:d}): {:d}, ND (p = {:d}): {:d}\n Operator "
               "assembly level: {}\n",
               h1_fespace.GetMaxElementOrder(), h1_fespace.GlobalTrueVSize(),
               nd_fespace.GetMaxElementOrder(), nd_fespace.GlobalTrueVSize(),
               nd_fespace.GetMaxElementOrder() > BilinearForm::pa_order_threshold
                   ? "Partial"
                   : "Full");

    auto &mesh = *h1_fespace.GetParMesh();
    const int q_order = fem::DefaultIntegrationOrder::Get(
        *h1_fespace.GetFE(0), *h1_fespace.GetFE(0), *mesh.GetElementTransformation(0));
    Mpi::Print(" Mesh geometries:\n");
    for (auto geom : mesh::CheckElements(mesh).GetGeomTypes())
    {
      const auto *fe = h1_fespace.FEColl()->FiniteElementForGeometry(geom);
      MFEM_VERIFY(fe, "MFEM does not support H1 spaces on geometry = "
                          << mfem::Geometry::Name[geom] << "!");
      Mpi::Print("  {}: P = {:d}, Q = {:d} (quadrature order = {:d})\n",
                 mfem::Geometry::Name[geom], fe->GetDof(),
                 mfem::IntRules.Get(geom, q_order).GetNPoints(), q_order);
    }

    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
}

}  // namespace

std::unique_ptr<Operator> LaplaceOperator::GetStiffnessMatrix()
{
  PrintHeader(GetH1Space(), GetNDSpace(), print_hdr);
  auto K = std::make_unique<MultigridOperator>(GetH1Spaces().GetNumLevels());
  for (std::size_t l = 0; l < GetH1Spaces().GetNumLevels(); l++)
  {
    // Force coarse level operator to be fully assembled always.
    const auto &h1_fespace_l = GetH1Spaces().GetFESpaceAtLevel(l);
    if (print_hdr)
    {
      Mpi::Print(" Level {:d} (p = {:d}): {:d} unknowns", l,
                 h1_fespace_l.GetMaxElementOrder(), h1_fespace_l.GlobalTrueVSize());
    }
    constexpr bool skip_zeros = false;
    MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal());
    BilinearForm k(h1_fespace_l);
    k.AddDomainIntegrator<DiffusionIntegrator>((mfem::MatrixCoefficient &)epsilon_func);
    auto K_l = std::make_unique<ParOperator>(
        (l > 0) ? k.Assemble(skip_zeros) : k.FullAssemble(skip_zeros), h1_fespace_l);
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

void LaplaceOperator::GetExcitationVector(int idx, const Operator &K, Vector &X,
                                          Vector &RHS)
{
  // Apply the Dirichlet BCs to the solution vector: V = 1 on terminal boundaries with the
  // given index, V = 0 on all ground and other terminal boundaries.
  mfem::ParGridFunction x(&GetH1Space());
  x = 0.0;

  // Get a marker of all boundary attributes with the given source surface index.
  const mfem::ParMesh &mesh = GetMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> source_marker = mesh::AttrToMarker(bdr_attr_max, source_attr_lists[idx]);
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
