// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "laplaceoperator.hpp"

#include <set>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

LaplaceOperator::LaplaceOperator(const IoData &iodata,
                                 const std::vector<std::unique_ptr<Mesh>> &mesh)
  : print_hdr(true), dbc_attr(SetUpBoundaryProperties(iodata, *mesh.back())),
    h1_fecs(fem::ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsening, false)),
    nd_fec(std::make_unique<mfem::ND_FECollection>(iodata.solver.order,
                                                   mesh.back()->Dimension())),
    rt_fecs(fem::ConstructFECollections<mfem::RT_FECollection>(
        iodata.solver.order - 1, mesh.back()->Dimension(),
        iodata.solver.linear.estimator_mg ? iodata.solver.linear.mg_max_levels : 1,
        iodata.solver.linear.mg_coarsening, false)),
    h1_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        iodata.solver.linear.mg_max_levels, mesh, h1_fecs, &dbc_attr, &dbc_tdof_lists)),
    nd_fespace(*mesh.back(), nd_fec.get()),
    rt_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::RT_FECollection>(
        iodata.solver.linear.estimator_mg ? iodata.solver.linear.mg_max_levels : 1, mesh,
        rt_fecs)),
    mat_op(iodata, *mesh.back()), source_attr_lists(ConstructSources(iodata))
{
  // Print essential BC information.
  if (dbc_attr.Size())
  {
    Mpi::Print("\nConfiguring Dirichlet BC at attributes:\n");
    utils::PrettyPrint(dbc_attr);
  }
}

mfem::Array<int> LaplaceOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                          const mfem::ParMesh &mesh)
{
  // Check that boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!iodata.boundaries.pec.empty() || !iodata.boundaries.lumpedport.empty())
  {
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    std::set<int> bdr_warn_list;
    for (auto attr : iodata.boundaries.pec.attributes)
    {
      // MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
      //             "Ground boundary attribute tags must be non-negative and correspond to
      //             " attributes in the mesh!");
      // MFEM_VERIFY(bdr_attr_marker[attr - 1],
      //             "Unknown ground boundary attribute " << attr << "!");
      if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
      {
        bdr_warn_list.insert(attr);
      }
    }
    if (!bdr_warn_list.empty())
    {
      Mpi::Print("\n");
      Mpi::Warning("Unknown ground boundary attributes!\nSolver will just ignore them!");
      utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
      Mpi::Print("\n");
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
  dbc_bcs.Reserve(static_cast<int>(iodata.boundaries.pec.attributes.size()) +
                  static_cast<int>(iodata.boundaries.lumpedport.size()));
  for (auto attr : iodata.boundaries.pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
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
    attr_list.Reserve(
        static_cast<int>(data.elements.size()));  // Average one attribute per element
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
                 const mfem::ParFiniteElementSpace &nd_fespace,
                 const mfem::ParFiniteElementSpace &rt_fespace, bool &print_hdr)
{
  if (print_hdr)
  {
    Mpi::Print("\nAssembling system matrices, number of global unknowns:\n"
               " H1 (p = {:d}): {:d}, ND (p = {:d}): {:d}, RT (p = {:d}): {:d}\n Operator "
               "assembly level: {}\n",
               h1_fespace.GetMaxElementOrder(), h1_fespace.GlobalTrueVSize(),
               nd_fespace.GetMaxElementOrder(), nd_fespace.GlobalTrueVSize(),
               rt_fespace.GetMaxElementOrder(), rt_fespace.GlobalTrueVSize(),
               (h1_fespace.GetMaxElementOrder() >= BilinearForm::pa_order_threshold)
                   ? "Partial"
                   : "Full");

    const auto &mesh = *h1_fespace.GetParMesh();
    const auto geom_types = mesh::CheckElements(mesh).GetGeomTypes();
    Mpi::Print(" Mesh geometries:\n");
    for (auto geom : geom_types)
    {
      const auto *fe = nd_fespace.FEColl()->FiniteElementForGeometry(geom);
      MFEM_VERIFY(fe, "MFEM does not support ND spaces on geometry = "
                          << mfem::Geometry::Name[geom] << "!");
      const int q_order = fem::DefaultIntegrationOrder::Get(mesh, geom);
      Mpi::Print("  {}: P = {:d}, Q = {:d} (quadrature order = {:d}){}\n",
                 mfem::Geometry::Name[geom], fe->GetDof(),
                 mfem::IntRules.Get(geom, q_order).GetNPoints(), q_order,
                 (geom == geom_types.back()) ? "" : ",");
    }

    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
}

}  // namespace

std::unique_ptr<Operator> LaplaceOperator::GetStiffnessMatrix()
{
  // When partially assembled, the coarse operators can reuse the fine operator quadrature
  // data if the spaces correspond to the same mesh.
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);

  constexpr bool skip_zeros = false;
  MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityReal());
  BilinearForm k(GetH1Space());
  k.AddDomainIntegrator<DiffusionIntegrator>(epsilon_func);
  // k.AssembleQuadratureData();
  auto k_vec = k.Assemble(GetH1Spaces(), skip_zeros);
  auto K = std::make_unique<MultigridOperator>(GetH1Spaces().GetNumLevels());
  for (std::size_t l = 0; l < GetH1Spaces().GetNumLevels(); l++)
  {
    const auto &h1_fespace_l = GetH1Spaces().GetFESpaceAtLevel(l);
    if (print_hdr)
    {
      Mpi::Print(" Level {:d} (p = {:d}): {:d} unknowns", l,
                 h1_fespace_l.GetMaxElementOrder(), h1_fespace_l.GlobalTrueVSize());
      if (const auto *k_spm = dynamic_cast<const hypre::HypreCSRMatrix *>(k_vec[l].get()))
      {
        HYPRE_BigInt nnz = k_spm->NNZ();
        Mpi::GlobalSum(1, &nnz, h1_fespace_l.GetComm());
        Mpi::Print(", {:d} NNZ\n", nnz);
      }
      else
      {
        Mpi::Print("\n");
      }
    }
    auto K_l = std::make_unique<ParOperator>(std::move(k_vec[l]), h1_fespace_l);
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
  mfem::ParGridFunction x(&GetH1Space().Get());
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
  X.UseDevice(true);
  RHS.UseDevice(true);
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
