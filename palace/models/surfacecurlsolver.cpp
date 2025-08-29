// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacecurlsolver.hpp"
#include "surfacefluxoperator.hpp"

#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

std::unique_ptr<Vector> SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data,
                                                const IoData &iodata, const Mesh &mesh,
                                                const FiniteElementSpace &nd_fespace,
                                                int flux_loop_idx)
{
  // Use flux loop configuration from SurfaceFluxData

  std::vector<double> flux_values = flux_data.flux_amounts;
  std::vector<int> surface_attrs = flux_data.fluxloop_pec;
  surface_attrs.insert(surface_attrs.end(), flux_data.hole_attributes.begin(),
                       flux_data.hole_attributes.end());

  Vector loop_normal(3);
  for (int i = 0; i < 3; i++)
  {
    loop_normal[i] = flux_data.direction[i];
  }

  // Use parameters
  const mfem::ParFiniteElementSpace *fespace = &nd_fespace.Get();
  int order = iodata.solver.order;

  MPI_Comm comm = mesh.GetComm();
  const auto &pmesh = mesh.Get();

  // First attribute is the metal surface, rest are holes
  int bdr_attr_metal = surface_attrs[0];
  mfem::Array<int> hole_surface_attrs;
  for (int i = 1; i < static_cast<int>(surface_attrs.size()); i++)
  {
    hole_surface_attrs.Append(surface_attrs[i]);
  }
  int num_holes = hole_surface_attrs.Size();

  // Extract hole boundary DOFs
  std::unordered_map<int, mfem::Array<int>> attr_to_elements;
  std::vector<std::unordered_map<int, int>> hole_dof_to_edge_maps(num_holes);
  std::vector<mfem::Array<int>> hole_ess_tdof_lists(num_holes);
  std::vector<mfem::Array<int>> hole_ldof_markers(num_holes);
  std::vector<std::unordered_set<int>> hole_boundary_edge_ldofs(num_holes);

  const_cast<mfem::ParFiniteElementSpace *>(fespace)->GetBoundaryElementsByAttribute(
      hole_surface_attrs, attr_to_elements);
  for (int h = 0; h < num_holes; h++)
  {
    const_cast<mfem::ParFiniteElementSpace *>(fespace)->GetBoundaryEdgeDoFs(
        attr_to_elements[hole_surface_attrs[h]], hole_ess_tdof_lists[h],
        hole_ldof_markers[h], hole_boundary_edge_ldofs[h], &hole_dof_to_edge_maps[h],
        nullptr, nullptr, nullptr);
  }

  // Create submesh
  mfem::Array<int> bdr_attr_marker(1);
  bdr_attr_marker[0] = bdr_attr_metal;
  mfem::ParSubMesh boundary_submesh =
      mfem::ParSubMesh::CreateFromBoundary(pmesh, bdr_attr_marker);

  // Extract submesh boundary edges
  mfem::Array<int> submesh_boundary_edge_ids;
  for (int i = 0; i < boundary_submesh.GetNBE(); i++)
  {
    mfem::Array<int> edges, orientations;
    boundary_submesh.GetBdrElementEdges(i, edges, orientations);
    for (int j = 0; j < edges.Size(); j++)
      submesh_boundary_edge_ids.Append(edges[j]);
  }

  // Map submesh to parent edges
  const mfem::Array<int> &parent_edge_ids = boundary_submesh.GetParentEdgeIDMap();
  std::unordered_map<int, int> submesh_to_parent_bdr_edge_map;
  for (int submesh_edge : submesh_boundary_edge_ids)
    submesh_to_parent_bdr_edge_map[submesh_edge] = parent_edge_ids[submesh_edge];

  // Match hole boundaries
  std::vector<mfem::Array<int>> hole_boundary_edges;
  mesh::MatchBoundaryEdges(pmesh, boundary_submesh, submesh_boundary_edge_ids,
                           submesh_to_parent_bdr_edge_map, hole_dof_to_edge_maps,
                           hole_boundary_edges);
  // Assign boundary attributes and create edge sets
  int current_attr =
      boundary_submesh.bdr_attributes.Size() > 0 ? boundary_submesh.bdr_attributes[0] : 1;
  std::vector<int> hole_boundary_attrs(num_holes);
  std::vector<std::unordered_set<int>> hole_edge_sets(num_holes);
  for (int h = 0; h < num_holes; h++)
  {
    hole_boundary_attrs[h] = current_attr + 1 + h;
    hole_edge_sets[h] = std::unordered_set<int>(hole_boundary_edges[h].begin(),
                                                hole_boundary_edges[h].end());
  }

  for (int i = 0; i < boundary_submesh.GetNBE(); i++)
  {
    mfem::Array<int> edges, orientations;
    boundary_submesh.GetBdrElementEdges(i, edges, orientations);
    for (int h = 0; h < num_holes; h++)
    {
      bool is_hole_boundary = false;
      for (int j = 0; j < edges.Size(); j++)
      {
        if (hole_edge_sets[h].count(edges[j]))
        {
          is_hole_boundary = true;
          break;
        }
      }
      if (is_hole_boundary)
      {
        boundary_submesh.GetBdrElement(i)->SetAttribute(hole_boundary_attrs[h]);
        break;
      }
    }
  }
  boundary_submesh.SetAttributes();

  // Compute hole properties using input flux values
  std::vector<mfem::Array<int>> hole_bdr_markers(num_holes);
  std::vector<double> hole_perimeters(num_holes);
  std::vector<double> hole_field_values(num_holes);
  std::vector<std::unordered_map<int, double>> hole_edge_lengths(num_holes);
  std::vector<std::unordered_map<int, int>> hole_edge_orientations(num_holes);

  mfem::H1_FECollection h1_fec(order, boundary_submesh.Dimension());
  mfem::ParFiniteElementSpace h1_fespace(&boundary_submesh, &h1_fec);
  mfem::ConstantCoefficient one(1.0);

  for (int h = 0; h < num_holes; h++)
  {
    hole_bdr_markers[h].SetSize(boundary_submesh.bdr_attributes.Max());
    hole_bdr_markers[h] = 0;
    hole_bdr_markers[h][hole_boundary_attrs[h] - 1] = 1;

    mfem::ParLinearForm perimeter_form(&h1_fespace);
    perimeter_form.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(one),
                                         hole_bdr_markers[h]);
    perimeter_form.Assemble();

    double local_perimeter = perimeter_form.Sum();
    MPI_Allreduce(&local_perimeter, &hole_perimeters[h], 1, MPI_DOUBLE, MPI_SUM,
                  boundary_submesh.GetComm());

    hole_field_values[h] = flux_values[h] / hole_perimeters[h];

    mesh::ComputeSubmeshBoundaryEdgeOrientations(boundary_submesh, hole_boundary_edges[h],
                                                 loop_normal, hole_edge_orientations[h],
                                                 hole_edge_lengths[h], order);
    for (auto &pair : hole_edge_lengths[h])
      pair.second = -pair.second;
  }

  // Create Nedelec space and solve
  mfem::ND_FECollection nd_fec(order, boundary_submesh.Dimension());
  mfem::ParFiniteElementSpace nd_fespace_submesh(&boundary_submesh, &nd_fec);

  std::vector<std::unordered_map<int, mfem::Array<int>>> hole_edge_to_dofs_maps(num_holes);
  std::vector<mfem::Array<int>> hole_boundary_edge_dofs(num_holes);

  for (int h = 0; h < num_holes; h++)
  {
    for (int i = 0; i < hole_boundary_edges[h].Size(); i++)
    {
      int edge_idx = hole_boundary_edges[h][i];
      mfem::Array<int> edge_dofs;
      nd_fespace_submesh.GetEdgeDofs(edge_idx, edge_dofs);
      hole_edge_to_dofs_maps[h][edge_idx] = edge_dofs;
      for (int j = 0; j < edge_dofs.Size(); j++)
        hole_boundary_edge_dofs[h].Append(edge_dofs[j]);
    }
  }

  mfem::Array<int> combined_inner_bdr_marker(boundary_submesh.bdr_attributes.Max());
  combined_inner_bdr_marker = 0;
  mfem::Array<int> ldof_marker_submesh(nd_fespace_submesh.GetVSize());
  ldof_marker_submesh = 0;
  mfem::ParGridFunction A(&nd_fespace_submesh);
  A = 0.0;

  // Directly apply loop BC by computing the integration of 1D Nedelec elements on bounndary
  // edges
  const mfem::IntegrationRule *ir = &mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * order);
  mfem::ND_SegmentElement nd_seg(order);
  for (int h = 0; h < num_holes; h++)
  {
    combined_inner_bdr_marker[hole_boundary_attrs[h] - 1] = 1;
    double field_strength = hole_field_values[h];

    for (const auto &pair : hole_edge_to_dofs_maps[h])
    {
      int edge = pair.first;
      const mfem::Array<int> &edge_dofs = pair.second;
      int orientation = hole_edge_orientations[h][edge];
      double coeff = -order * field_strength * orientation;

      mfem::IsoparametricTransformation edge_trans;
      boundary_submesh.GetEdgeTransformation(edge, &edge_trans);

      mfem::DenseMatrix vshape(edge_dofs.Size(), 1);

      for (int j = 0; j < edge_dofs.Size(); j++)
      {
        double dof_value = 0.0;

        for (int q = 0; q < ir->GetNPoints(); q++)
        {
          const mfem::IntegrationPoint &ip = ir->IntPoint(q);
          edge_trans.SetIntPoint(&ip);
          nd_seg.CalcVShape(ip, vshape);

          dof_value += ip.weight * edge_trans.Weight() * coeff * vshape(j, 0);
        }

        A(edge_dofs[j]) = dof_value;
      }
    }
    // Mark DOFs for synchronization
    for (int i = 0; i < hole_boundary_edge_dofs[h].Size(); i++)
      ldof_marker_submesh[hole_boundary_edge_dofs[h][i]] = 1;
  }

  // Notify start of 2D surface curl problem solving
  Mpi::Print("\nSolving 2D surface curl problem for flux loop boundary conditions...\n");

  // Create Palace finite element space hierarchy for the submesh
  std::vector<std::unique_ptr<mfem::ND_FECollection>> submesh_nd_fecs;
  submesh_nd_fecs.push_back(
      std::make_unique<mfem::ND_FECollection>(order, boundary_submesh.Dimension()));
  std::vector<std::unique_ptr<mfem::H1_FECollection>> submesh_h1_fecs;
  submesh_h1_fecs.push_back(
      std::make_unique<mfem::H1_FECollection>(order, boundary_submesh.Dimension()));
  std::vector<std::unique_ptr<Mesh>> submesh_vec;
  submesh_vec.push_back(std::make_unique<Mesh>(boundary_submesh, 1));

  FiniteElementSpaceHierarchy submesh_nd_fespaces(
      fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
          1, submesh_vec, submesh_nd_fecs, nullptr, nullptr));
  FiniteElementSpaceHierarchy submesh_h1_fespaces(
      fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
          1, submesh_vec, submesh_h1_fecs, nullptr, nullptr));

  mfem::Array<int> submesh_ess_tdof_list;
  submesh_nd_fespaces.GetFESpaceAtLevel(0).Get().GetEssentialTrueDofs(
      combined_inner_bdr_marker, submesh_ess_tdof_list);

  // First synchronize the marker itself to ensure all processors agree on which DoFs to
  // sync
  mfem::GroupCommunicator *gc =
      submesh_nd_fespaces.GetFESpaceAtLevel(0).Get().ScalarGroupComm();
  mfem::Array<int> global_marker(ldof_marker_submesh);
  gc->Reduce<int>(global_marker.GetData(), mfem::GroupCommunicator::BitOR<int>);
  gc->Bcast(global_marker);

  // Synchronize the edge dof boundary condition across processors with MaxAbs reduction
  mfem::Array<double> values(A.GetData(), A.Size());
  gc->ReduceBegin(values.GetData());
  gc->ReduceMarked<double>(values.GetData(), global_marker, 0,
                           mfem::GroupCommunicator::MaxAbs<double>);
  gc->Bcast(values.GetData());
  delete gc;

  // Create unit coefficient for curl-curl term and the regularization coefficient
  MaterialPropertyCoefficient reg_coeff(boundary_submesh.attributes.Max());
  MaterialPropertyCoefficient unit_coeff(boundary_submesh.attributes.Max());
  for (int attr = 1; attr <= boundary_submesh.attributes.Max(); attr++)
  {
    unit_coeff.AddMaterialProperty(attr, 1.0);
    reg_coeff.AddMaterialProperty(attr, flux_data.regularization);
  }

  // Use Palace's BilinearForm and assemble system matrix
  BilinearForm a(submesh_nd_fespaces.GetFESpaceAtLevel(0));
  a.AddDomainIntegrator<CurlCurlIntegrator>(unit_coeff);
  a.AddDomainIntegrator<VectorFEMassIntegrator>(reg_coeff);

  auto k_vec = a.Assemble(submesh_nd_fespaces, false);
  auto K_op = std::make_unique<ParOperator>(std::move(k_vec[0]),
                                            submesh_nd_fespaces.GetFESpaceAtLevel(0));

  // Compute boundary-interior coupling BEFORE applying boundary conditions
  Vector RHS(submesh_nd_fespaces.GetFESpaceAtLevel(0).GetTrueVSize());
  Vector X(submesh_nd_fespaces.GetFESpaceAtLevel(0).GetTrueVSize());
  Vector boundary_vals(submesh_nd_fespaces.GetFESpaceAtLevel(0).GetTrueVSize());
  RHS.UseDevice(true);
  X.UseDevice(true);
  boundary_vals.UseDevice(true);

  // Get boundary values directly from MFEM GridFunction
  A.GetTrueDofs(boundary_vals);

  // Compute RHS = -K * boundary_values (coupling term)
  K_op->Mult(boundary_vals, RHS);
  RHS *= -1.0;

  // Set boundary values in RHS
  linalg::SetSubVector(RHS, submesh_ess_tdof_list, boundary_vals);

  // Apply boundary conditions to operator
  K_op->SetEssentialTrueDofs(submesh_ess_tdof_list, Operator::DiagonalPolicy::DIAG_ONE);

  // Set up Palace KSP solver
  KspSolver ksp(iodata, submesh_nd_fespaces, &submesh_h1_fespaces);
  ksp.SetOperators(*K_op, *K_op);

  X = 0.0;
  ksp.Mult(RHS, X);

  // Set solution directly in MFEM GridFunction
  A.SetFromTrueDofs(X);

  /*
  // Set up and solve system
  mfem::ParBilinearForm a(&nd_fespace_submesh);

  // Add curl term: ∫ (curl A) · (curl v) dΩ
  mfem::ConstantCoefficient curl_reg(1.0);
  a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(curl_reg));

  // Add small regularization for stability
  mfem::ConstantCoefficient reg_param(1e-6);
  a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(reg_param));

  nd_fespace_submesh.GetEssentialTrueDofs(combined_inner_bdr_marker, submesh_ess_tdof_list);
  a.Assemble();
  mfem::ParLinearForm b(&nd_fespace_submesh);
  b.Assemble();

  mfem::HypreParMatrix A_mat;
  Vector B, X;
  a.FormLinearSystem(submesh_ess_tdof_list, A, b, A_mat, X, B);

  mfem::GMRESSolver gmres(MPI_COMM_WORLD);
  gmres.SetOperator(A_mat);
  gmres.SetRelTol(1e-8);
  gmres.SetMaxIter(1000);
  gmres.SetPrintLevel(0);
  gmres.Mult(B, X);

  a.RecoverFEMSolution(X, b, A);
  */

  // Transfer to parent mesh
  mfem::ParGridFunction A_3d(const_cast<mfem::ParFiniteElementSpace *>(fespace));
  A_3d = 0.0;
  mfem::ParSubMesh::Transfer(A, A_3d);

  // Extract true DOFs and return as Vector
  auto result = std::make_unique<Vector>();
  result->SetSize(fespace->GetTrueVSize());
  result->UseDevice(true);
  A_3d.GetTrueDofs(*result);

  // Clear large temporary objects to free memory
  hole_dof_to_edge_maps.clear();
  hole_ess_tdof_lists.clear();
  hole_ldof_markers.clear();
  hole_boundary_edge_ldofs.clear();
  hole_boundary_edges.clear();
  hole_edge_sets.clear();
  hole_edge_to_dofs_maps.clear();
  hole_boundary_edge_dofs.clear();
  hole_bdr_markers.clear();
  hole_edge_lengths.clear();
  attr_to_elements.clear();
  submesh_to_parent_bdr_edge_map.clear();

  // Clear Palace-specific objects
  submesh_nd_fecs.clear();
  submesh_h1_fecs.clear();
  submesh_vec.clear();

  return result;
}

void VerifyFluxThroughHoles(const mfem::ParGridFunction &B_gf,
                            const std::vector<int> &hole_attributes,
                            const std::vector<double> &target_fluxes, const Mesh &mesh,
                            MPI_Comm comm)
{
  for (int h = 0; h < static_cast<int>(hole_attributes.size()); h++)
  {
    int hole_attr = hole_attributes[h];
    double target_flux = target_fluxes[h];

    // Compute flux Φ = ∫ B·n dS through hole surface
    mfem::Array<int> hole_marker(mesh.Get().bdr_attributes.Max());
    hole_marker = 0;
    hole_marker[hole_attr - 1] = 1;

    mfem::ParLinearForm flux_form(B_gf.ParFESpace());
    flux_form.AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator,
                                    hole_marker);
    flux_form.Assemble();

    double computed_flux = flux_form * B_gf;
    double global_flux;
    MPI_Allreduce(&computed_flux, &global_flux, 1, MPI_DOUBLE, MPI_SUM, comm);

    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0)
    {
      Mpi::Print("Hole attribute {:d}: Target flux = {:.6e}, Computed flux = {:.6e}, Error "
                 "= {:.6e}\n",
                 hole_attr, target_flux, global_flux, std::abs(global_flux - target_flux));
    }
  }
}

void VerifyFluxThroughAllHoles(const mfem::ParGridFunction &B_gf, const IoData &iodata,
                               int current_flux_loop_idx, const Mesh &mesh, MPI_Comm comm)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  if (rank == 0)
  {
    Mpi::Print("FluxLoop {:d} excitation - Flux through all holes:\n",
               current_flux_loop_idx);
  }

  // Compute flux through all holes in all flux loops
  for (const auto &[loop_idx, flux_data] : iodata.boundaries.fluxloop)
  {
    for (int h = 0; h < static_cast<int>(flux_data.hole_attributes.size()); h++)
    {
      int hole_attr = flux_data.hole_attributes[h];
      double target_flux =
          (loop_idx == current_flux_loop_idx) ? flux_data.flux_amounts[h] : 0.0;

      mfem::Array<int> hole_marker(mesh.Get().bdr_attributes.Max());
      hole_marker = 0;
      hole_marker[hole_attr - 1] = 1;

      mfem::ParLinearForm flux_form(B_gf.ParFESpace());
      flux_form.AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator,
                                      hole_marker);
      flux_form.Assemble();

      double computed_flux = flux_form * B_gf;
      double global_flux;
      MPI_Allreduce(&computed_flux, &global_flux, 1, MPI_DOUBLE, MPI_SUM, comm);

      if (rank == 0)
      {
        Mpi::Print(
            "  Loop {:d} Hole {:d}: Target = {:.6e}, Computed = {:.6e}, Error = {:.6e}\n",
            loop_idx, hole_attr, target_flux, global_flux,
            std::abs(global_flux - target_flux));
      }
    }
  }
}

}  // namespace palace