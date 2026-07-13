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
#include "utils/tablecsv.hpp"

namespace palace
{

LaplaceOperator::LaplaceOperator(const config::BoundaryData &boundaries,
                                 const config::SolverData &solver,
                                 const std::vector<config::MaterialData> &materials,
                                 ProblemType problem_type,
                                 const std::vector<std::unique_ptr<Mesh>> &mesh)
  : print_hdr(true),
    dbc_attr(SetUpBoundaryProperties(boundaries.pec, boundaries.terminal, *mesh.back())),
    h1_fecs(fem::ConstructFECollections<mfem::H1_FECollection>(
        solver.order, mesh.back()->Dimension(), solver.linear.mg_max_levels,
        solver.linear.mg_coarsening, false)),
    nd_fec(std::make_unique<mfem::ND_FECollection>(solver.order, mesh.back()->Dimension())),
    rt_fecs(fem::ConstructFECollections<mfem::RT_FECollection>(
        solver.order - 1, mesh.back()->Dimension(),
        solver.linear.estimator_mg ? solver.linear.mg_max_levels : 1,
        solver.linear.mg_coarsening, false)),
    h1_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        solver.linear.mg_max_levels, mesh, h1_fecs, &dbc_attr, &dbc_tdof_lists)),
    nd_fespace(*mesh.back(), nd_fec.get()),
    rt_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::RT_FECollection>(
        solver.linear.estimator_mg ? solver.linear.mg_max_levels : 1, mesh, rt_fecs)),
    mat_op(materials, boundaries.periodic, problem_type, *mesh.back()),
    source_attr_lists(ConstructSources(boundaries.terminal))
{
  // Print essential BC information.
  if (dbc_attr.Size())
  {
    Mpi::Print("\nConfiguring Dirichlet BC at attributes:\n");
    utils::PrettyPrint(dbc_attr);
  }
}

LaplaceOperator::LaplaceOperator(const IoData &iodata,
                                 const std::vector<std::unique_ptr<Mesh>> &mesh)
  : LaplaceOperator(iodata.boundaries, iodata.solver, iodata.domains.materials,
                    iodata.problem.type, mesh)
{
}

mfem::Array<int> LaplaceOperator::SetUpBoundaryProperties(
    const config::PecBoundaryData &pec, const std::map<int, config::TerminalData> &terminal,
    const mfem::ParMesh &mesh)
{
  // Check that boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!pec.empty() || !terminal.empty())
  {
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    std::set<int> bdr_warn_list;
    for (auto attr : pec.attributes)
    {
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
    for (const auto &[idx, data] : terminal)
    {
      for (auto attr : data.attributes)
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

  // Mark selected boundary attributes from the mesh as essential (Dirichlet).
  mfem::Array<int> dbc_bcs;
  dbc_bcs.Reserve(static_cast<int>(pec.attributes.size()) +
                  static_cast<int>(terminal.size()));
  for (auto attr : pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
    {
      continue;
    }
    dbc_bcs.Append(attr);
  }
  for (const auto &[idx, data] : terminal)
  {
    for (auto attr : data.attributes)
    {
      dbc_bcs.Append(attr);
    }
  }
  MFEM_VERIFY(dbc_bcs.Size() > 0,
              "Electrostatic problem is ill-posed without any Dirichlet boundaries!");
  return dbc_bcs;
}

std::map<int, mfem::Array<int>>
LaplaceOperator::ConstructSources(const std::map<int, config::TerminalData> &terminal)
{
  // Construct mapping from terminal index to list of associated attributes.
  std::map<int, mfem::Array<int>> attr_lists;
  for (const auto &[idx, data] : terminal)
  {
    mfem::Array<int> &attr_list = attr_lists[idx];
    attr_list.Reserve(static_cast<int>(data.attributes.size()));
    for (auto attr : data.attributes)
    {
      attr_list.Append(attr);
    }
  }
  return attr_lists;
}

namespace
{

const ParOperator &GetFinestParOperator(const Operator &K)
{
  const Operator *K_finest = &K;
  if (const auto *K_mg = dynamic_cast<const MultigridOperator *>(&K))
  {
    K_finest = &K_mg->GetFinestOperator();
  }
  const auto *K_par = dynamic_cast<const ParOperator *>(K_finest);
  MFEM_VERIFY(K_par, "LaplaceOperator requires a ParOperator at the finest level!");
  return *K_par;
}

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
    const auto geom_types = mesh::CheckElements(mesh).GetGeomTypes(mesh.Dimension());
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
  auto K = AssembleStiffnessMatrix(GetH1Spaces(), mat_op, dbc_tdof_lists, print_hdr);
  print_hdr = false;
  return K;
}

std::unique_ptr<Operator> LaplaceOperator::AssembleStiffnessMatrix(
    FiniteElementSpaceHierarchy &h1_fespaces, const MaterialOperator &mat_op,
    const std::vector<mfem::Array<int>> &dbc_tdof_lists, bool print_hdr)
{
  MFEM_VERIFY(h1_fespaces.GetNumLevels() == dbc_tdof_lists.size(),
              "Laplace operator requires one essential true dof list per H1 level!");
  constexpr bool skip_zeros = false;
  MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityReal());
  BilinearForm k(h1_fespaces.GetFinestFESpace());
  k.AddDomainIntegrator<DiffusionIntegrator>(epsilon_func);
  // k.AssembleQuadratureData();
  auto k_vec = k.Assemble(h1_fespaces, skip_zeros);
  auto K = std::make_unique<MultigridOperator>(h1_fespaces.GetNumLevels());
  for (std::size_t l = 0; l < h1_fespaces.GetNumLevels(); l++)
  {
    const auto &h1_fespace_l = h1_fespaces.GetFESpaceAtLevel(l);
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
  return K;
}

std::unique_ptr<Operator>
LaplaceOperator::GetUnconstrainedStiffnessOperator(const Operator &K)
{
  const auto &K_finest = GetFinestParOperator(K);
  return std::make_unique<ParOperator>(K_finest.LocalOperator(),
                                       K_finest.TrialFiniteElementSpace());
}

void LaplaceOperator::GetExcitationVector(int idx, const Operator &K, Vector &X,
                                          Vector &RHS)
{
  GetExcitationVector(GetH1Space(), source_attr_lists.at(idx), K, X, RHS);
}

void LaplaceOperator::GetExcitationVector(FiniteElementSpace &h1_fespace,
                                          const mfem::Array<int> &source_attr,
                                          const Operator &K, Vector &X, Vector &RHS)
{
  // Apply the Dirichlet BCs to the solution vector: V = 1 on terminal boundaries with the
  // given index, V = 0 on all ground and other terminal boundaries.
  mfem::ParGridFunction x(&h1_fespace.Get());
  x = 0.0;

  // Get a marker of all boundary attributes with the given source surface index.
  const mfem::ParMesh &mesh = h1_fespace.GetMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> source_marker = mesh::AttrToMarker(bdr_attr_max, source_attr);
  mfem::ConstantCoefficient one(1.0);
  x.ProjectBdrCoefficient(one, source_marker);  // Values are only correct on master

  // Eliminate the essential BC to get the RHS vector.
  X.SetSize(h1_fespace.GetTrueVSize());
  RHS.SetSize(h1_fespace.GetTrueVSize());
  X.UseDevice(true);
  RHS.UseDevice(true);
  X = 0.0;
  RHS = 0.0;
  x.ParallelProject(X);  // Restrict to the true dofs
  GetFinestParOperator(K).EliminateRHS(X, RHS);
}

void LaplaceOperator::GetExcitationVector(FiniteElementSpace &h1_fespace,
                                          const std::vector<int> &source_attr,
                                          const Operator &K, Vector &X, Vector &RHS)
{
  mfem::Array<int> attr;
  attr.Reserve(static_cast<int>(source_attr.size()));
  for (int a : source_attr)
  {
    attr.Append(a);
  }
  GetExcitationVector(h1_fespace, attr, K, X, RHS);
}

mfem::DenseMatrix LaplaceOperator::ComputeCapacitanceMatrix(MPI_Comm comm,
                                                            const Operator &energy_op,
                                                            const std::vector<Vector> &V)
{
  const int n_term = static_cast<int>(V.size());
  MFEM_VERIFY(n_term > 0,
              "Cannot compute a capacitance matrix without terminal solutions!");
  MFEM_VERIFY(energy_op.Height() == V.front().Size() &&
                  energy_op.Width() == V.front().Size(),
              "Invalid electrostatic energy operator dimensions!");

  mfem::DenseMatrix C(n_term);
  Vector KV(energy_op.Height());
  KV.UseDevice(true);
  for (int i = 0; i < n_term; i++)
  {
    MFEM_VERIFY(V[i].Size() == energy_op.Width(),
                "Invalid terminal solution dimension for capacitance calculation!");
    energy_op.Mult(V[i], KV);
    for (int j = i; j < n_term; j++)
    {
      C(i, j) = C(j, i) = linalg::Dot<Vector>(comm, V[j], KV);
    }
  }
  return C;
}

void LaplaceOperator::WriteTerminalMatrix(MPI_Comm comm, const fs::path &post_dir,
                                          std::string_view file, std::string_view name,
                                          std::string_view unit,
                                          const std::vector<int> &terminal_indices,
                                          const mfem::DenseMatrix &mat, double scale)
{
  MFEM_VERIFY(mat.Height() == static_cast<int>(terminal_indices.size()) &&
                  mat.Width() == static_cast<int>(terminal_indices.size()),
              "Terminal matrix dimensions do not match the number of terminal indices!");
  if (!Mpi::Root(comm))
  {
    return;
  }

  TableWithCSVFile output(post_dir / file);
  output.table.insert(Column("i", "i", 0, 0, 2, ""));
  for (int idx : terminal_indices)
  {
    output.table["i"] << idx;
  }
  for (int j = 0; j < mat.Width(); j++)
  {
    const int idx = terminal_indices[j];
    const auto key = fmt::format("i2{}", idx);
    output.table.insert(key, fmt::format("{}[i][{}] {}", name, idx, unit));
    auto &col = output.table[key];
    for (int i = 0; i < mat.Height(); i++)
    {
      col << mat(i, j) * scale;
    }
  }
  output.WriteFullTableTrunc();
}

}  // namespace palace
