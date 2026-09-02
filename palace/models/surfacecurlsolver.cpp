// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacecurlsolver.hpp"
#include "surfacefluxoperator.hpp"

#include <algorithm>
#include <cstdlib>
#include <unordered_set>
#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

// Build the London drive a_h as a 3D curl-free cohomology generator carrying the hole
// fluxoid ∮_∂hole a_h·dl = Φ: a_h = ∇θ with θ jumping by Φ across a cut half-plane S
// bounded by the vertical flux line L = {x=cx, y=cy, all z} through the hole centroid. It
// is assembled as a lowest-order (Whitney) edge cochain (edge DOF = ±Φ across S, 0
// otherwise) then projected exactly into the order-p ND space (ND_1 ⊂ ND_p). curl a_h is
// supported only on faces pierced by L, which threads the hole opening off the film Σ, so
// a_h|Σ is curl-free with circulation Φ. Its gradient part is absorbable by A → A + ∇χ, so
// the extracted inductance depends only on the cohomology class (Φ) and is gauge-invariant.
//
// PALACE_CUT_DIR selects the cut orientation (default "y": S = {x=cx, y≥cy}; "x": the
// gauge-equivalent S = {y=cy, x≥cx}) for the gauge-invariance gate. Overall sign is
// irrelevant: the downstream London normalization rescales a_h so cᵀa_h = Φ.
Vector BuildCutCohomologyGenerator(const SurfaceFluxData &flux_data,
                                   const mfem::ParFiniteElementSpace &ndp_fespace,
                                   const Mesh &mesh)
{
  auto &pmesh = const_cast<mfem::ParMesh &>(mesh.Get());
  MPI_Comm comm = pmesh.GetComm();
  const int sdim = pmesh.SpaceDimension();
  MFEM_VERIFY(sdim == 3, "London cut cohomology generator requires a 3D mesh!");

  // Total imposed fluxoid Φ.
  double phi = 0.0;
  for (double f : flux_data.flux_amounts)
  {
    phi += f;
  }

  // Hole centroid (cx, cy): average of hole-boundary vertex coordinates. The vertical flux
  // line L = {x=cx, y=cy, all z} passes through the hole opening (empty region), never the
  // film.
  std::unordered_set<int> hole_attrs(flux_data.hole_attributes.begin(),
                                     flux_data.hole_attributes.end());
  double csum[3] = {0.0, 0.0, 0.0};
  double cnt = 0.0;
  mfem::Array<int> bverts;
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    if (!hole_attrs.count(pmesh.GetBdrAttribute(be)))
    {
      continue;
    }
    pmesh.GetBdrElementVertices(be, bverts);
    for (int v : bverts)
    {
      const double *x = pmesh.GetVertex(v);
      csum[0] += x[0];
      csum[1] += x[1];
      csum[2] += x[2];
      cnt += 1.0;
    }
  }
  Mpi::GlobalSum(3, csum, comm);
  Mpi::GlobalSum(1, &cnt, comm);
  MFEM_VERIFY(cnt > 0.0, "No hole boundary elements found for London cut generator!");
  const double cx = csum[0] / cnt, cy = csum[1] / cnt;

  // Cut half-plane orientation (see function comment). Default S = {x=cx, y≥cy};
  // PALACE_CUT_DIR=x selects the gauge-equivalent S = {y=cy, x≥cx} sharing the same flux
  // line L.
  const char *cut_env = std::getenv("PALACE_CUT_DIR");
  const bool cut_x = (cut_env && (cut_env[0] == 'x' || cut_env[0] == 'X'));

  // Lowest-order (Whitney) cut cochain on an ND_1 space over the full 3D mesh. Each rank
  // sets its local edge DOFs from vertex COORDINATES (identical on all ranks for a shared
  // edge) in the canonical GetEdgeVertices orientation (ev0→ev1), so the owner's true-DOF
  // value is correct and GetTrueDofs needs no cross-rank sign reconciliation.
  mfem::ND_FECollection nd1_fec(1, sdim);
  mfem::ParFiniteElementSpace nd1_fespace(&pmesh, &nd1_fec);
  mfem::ParGridFunction ah1(&nd1_fespace);
  ah1.UseDevice(false);
  ah1 = 0.0;
  mfem::Array<int> ev, edofs;
  for (int e = 0; e < pmesh.GetNEdges(); e++)
  {
    pmesh.GetEdgeVertices(e, ev);
    const double *x0 = pmesh.GetVertex(ev[0]);
    const double *x1 = pmesh.GetVertex(ev[1]);
    // d = signed distance to the cut PLANE; s = the in-plane half-plane selector
    // coordinate.
    double d0, d1, s0, s1;
    if (!cut_x)
    {
      d0 = x0[0] - cx;  // plane x=cx
      d1 = x1[0] - cx;
      s0 = x0[1];  // half-plane y≥cy
      s1 = x1[1];
    }
    else
    {
      d0 = x0[1] - cy;  // plane y=cy
      d1 = x1[1] - cy;
      s0 = x0[0];  // half-plane x≥cx
      s1 = x1[0];
    }
    if (d0 * d1 >= 0.0)
    {
      continue;  // edge does not cross the cut plane
    }
    const double t = d0 / (d0 - d1);
    const double scoord = s0 + t * (s1 - s0);
    if ((!cut_x && scoord < cy) || (cut_x && scoord < cx))
    {
      continue;  // crossing point lies off the half-plane
    }
    // Edge ev0→ev1 crosses the cut: θ jumps +Φ toward the +normal (d1>0) side of the plane.
    nd1_fespace.GetEdgeDofs(e, edofs);
    ah1(edofs[0]) = (d1 > 0.0) ? phi : -phi;
  }

  // Project the ND_1 cut field exactly into the order-p ND space. The per-element field is
  // a degree-1 polynomial, integrated exactly by the ND_p DOF functionals, so ND_1 ⊂ ND_p
  // is reproduced and tangential continuity (hence curl-free-on-Σ, circulation Φ) is
  // preserved.
  mfem::VectorGridFunctionCoefficient ah1_coeff(&ah1);
  mfem::ParGridFunction ahp(const_cast<mfem::ParFiniteElementSpace *>(&ndp_fespace));
  ahp.UseDevice(false);
  ahp = 0.0;
  ahp.ProjectCoefficient(ah1_coeff);

  Vector result(ndp_fespace.GetTrueVSize());
  result.UseDevice(true);
  ahp.GetTrueDofs(result);
  return result;
}

}  // namespace

Vector SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const IoData &iodata,
                               const Mesh &mesh, const FiniteElementSpace &nd_fespace,
                               int flux_loop_idx,
                               PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                               bool harmonic_generator)
{
  Vector result;
  SolveSurfaceCurlProblem(flux_data, iodata, mesh, nd_fespace, flux_loop_idx, post_op,
                          result, harmonic_generator);
  return result;
}

void SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const IoData &iodata,
                             const Mesh &mesh, const FiniteElementSpace &nd_fespace,
                             int flux_loop_idx,
                             PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                             Vector &result, bool harmonic_generator)
{
  const mfem::ParFiniteElementSpace *fespace = &nd_fespace.Get();
  int order = iodata.solver.order;

  const auto &pmesh = mesh.Get();

  // London drive: build a_h directly on the 3D ND space as a curl-free cohomology generator
  // via a topological cut cochain (gauge-invariant; see BuildCutCohomologyGenerator). No
  // submesh needed.
  if (harmonic_generator)
  {
    result = BuildCutCohomologyGenerator(flux_data, *fespace, mesh);
    // Populate the post_op A buffer for consistency with the pure-PEC path (used only as
    // scratch).
    auto &A_3d = post_op.GetAGridFunction().Real();
    A_3d.SetFromTrueDofs(result);
    return;
  }

  // Extract metal surface and hole attributes from flux_data
  mfem::Array<int> metal_surface_attrs;
  for (int metal_attr : flux_data.fluxloop_pec)
  {
    metal_surface_attrs.Append(metal_attr);
  }

  // Validate that we have at least one metal surface attribute
  MFEM_VERIFY(metal_surface_attrs.Size() > 0,
              "At least one metal surface attribute must be specified in FluxLoopPEC!");
  mfem::Array<int> hole_surface_attrs;
  for (int hole_attr : flux_data.hole_attributes)
  {
    hole_surface_attrs.Append(hole_attr);
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

  // Create submesh from all metal surface attributes
  mfem::ParSubMesh boundary_submesh =
      mfem::ParSubMesh::CreateFromBoundary(pmesh, metal_surface_attrs);

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

  // Assign boundary attributes and create edge sets.
  // Find the maximum existing boundary attribute across ALL ranks to avoid conflicts.
  int current_attr = 1;
  if (boundary_submesh.bdr_attributes.Size() > 0)
  {
    current_attr = boundary_submesh.bdr_attributes.Max();
  }
  Mpi::GlobalMax(1, &current_attr, boundary_submesh.GetComm());
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
      for (int j = 0; j < edges.Size() && !is_hole_boundary; j++)
      {
        if (hole_edge_sets[h].count(edges[j]))
        {
          is_hole_boundary = true;
          boundary_submesh.GetBdrElement(i)->SetAttribute(hole_boundary_attrs[h]);
        }
      }
      if (is_hole_boundary)
      {
        break;  // Exit hole loop once we've assigned an attribute
      }
    }
  }
  boundary_submesh.SetAttributes();

  // Compute hole properties using input flux values.
  // Use the maximum hole boundary attribute for marker sizing since local
  // bdr_attributes.Max() may be smaller on ranks without hole boundary elements.
  int marker_size = hole_boundary_attrs.back();
  std::vector<mfem::Array<int>> hole_bdr_markers(num_holes);
  std::vector<double> hole_perimeters(num_holes);
  std::vector<double> hole_field_values(num_holes);
  std::vector<std::unordered_map<int, double>> hole_edge_lengths(num_holes);

  mfem::GridFunction *nodes = boundary_submesh.GetNodes();
  const mfem::FiniteElementCollection *h1_fec = nodes->FESpace()->FEColl();
  mfem::ParFiniteElementSpace h1_pfespace(&boundary_submesh, h1_fec);
  mfem::ConstantCoefficient one(1.0);

  // Create loop normal vector from flux data direction
  Vector loop_normal(const_cast<double *>(flux_data.direction.data()), 3);

  for (int h = 0; h < num_holes; h++)
  {
    hole_bdr_markers[h].SetSize(marker_size);
    hole_bdr_markers[h] = 0;
    hole_bdr_markers[h][hole_boundary_attrs[h] - 1] = 1;

    mfem::ParLinearForm perimeter_form(&h1_pfespace);
    perimeter_form.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(one),
                                         hole_bdr_markers[h]);
    perimeter_form.Assemble();
    double local_perimeter = perimeter_form.Sum();
    hole_perimeters[h] = local_perimeter;
    Mpi::GlobalSum(1, &hole_perimeters[h], boundary_submesh.GetComm());
    hole_field_values[h] = flux_data.flux_amounts[h] / hole_perimeters[h];

    mesh::ComputeSubmeshBoundaryEdgeOrientations(boundary_submesh, hole_boundary_edges[h],
                                                 loop_normal, hole_edge_lengths[h], order);
  }

  // Create Nedelec space and solve
  mfem::ND_FECollection nd_fec(order, boundary_submesh.Dimension());
  mfem::ParFiniteElementSpace nd_fespace_submesh(&boundary_submesh, &nd_fec);

  std::vector<std::unordered_map<int, mfem::Array<int>>> hole_edge_to_dofs_maps(num_holes);
  std::vector<mfem::Array<int>> hole_boundary_edge_dofs(num_holes);

  mfem::Array<int> edge_dofs;
  for (int h = 0; h < num_holes; h++)
  {
    for (int i = 0; i < hole_boundary_edges[h].Size(); i++)
    {
      int edge_idx = hole_boundary_edges[h][i];
      nd_fespace_submesh.GetEdgeDofs(edge_idx, edge_dofs);
      hole_edge_to_dofs_maps[h][edge_idx] = edge_dofs;
      hole_boundary_edge_dofs[h].Append(edge_dofs);
    }
  }

  mfem::Array<int> combined_inner_bdr_marker(marker_size);
  combined_inner_bdr_marker = 0;
  mfem::Array<int> ldof_marker_submesh(nd_fespace_submesh.GetVSize());
  ldof_marker_submesh = 0;
  mfem::ParGridFunction A(&nd_fespace_submesh);
  A.UseDevice(false);
  A = 0.0;

  // Directly apply loop BC by computing the integration of 1D Nedelec elements on boundary
  // edges
  for (int h = 0; h < num_holes; h++)
  {
    combined_inner_bdr_marker[hole_boundary_attrs[h] - 1] = 1;
    for (const auto &pair : hole_edge_to_dofs_maps[h])
    {
      int edge = pair.first;
      const mfem::Array<int> &edge_dofs = pair.second;
      double oriented_length = hole_edge_lengths[h][edge];
      for (int j = 0; j < edge_dofs.Size(); j++)
        A(edge_dofs[j]) = hole_field_values[h] * oriented_length;
    }

    // Mark DOFs for synchronization
    for (int i = 0; i < hole_boundary_edge_dofs[h].Size(); i++)
    {
      ldof_marker_submesh[hole_boundary_edge_dofs[h][i]] = 1;
    }
  }

  // Notify start of 2D surface curl problem solving
  Mpi::Print("\nSolving 2D surface curl problem for flux loop boundary conditions...\n");

  // Create Palace finite element space hierarchy for the submesh using P-multigrid
  // Construct P-multigrid FE collections using solver configuration
  int mg_max_levels = iodata.solver.linear.mg_max_levels;
  auto mg_coarsening = iodata.solver.linear.mg_coarsening;
  std::vector<std::unique_ptr<mfem::ND_FECollection>> submesh_nd_fecs =
      fem::ConstructFECollections<mfem::ND_FECollection>(
          order, boundary_submesh.Dimension(), mg_max_levels, mg_coarsening, false);
  std::vector<std::unique_ptr<mfem::H1_FECollection>> submesh_h1_fecs =
      fem::ConstructFECollections<mfem::H1_FECollection>(
          order, boundary_submesh.Dimension(), mg_max_levels, mg_coarsening, false);
  std::vector<std::unique_ptr<Mesh>> submesh_vec;
  submesh_vec.push_back(std::make_unique<Mesh>(boundary_submesh, 1));

  FiniteElementSpaceHierarchy submesh_nd_fespaces(
      fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
          mg_max_levels, submesh_vec, submesh_nd_fecs, nullptr, nullptr));
  FiniteElementSpaceHierarchy submesh_h1_fespaces(
      fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
          mg_max_levels, submesh_vec, submesh_h1_fecs, nullptr, nullptr));

  // Create separate essential DoF lists for each level in the hierarchy
  std::vector<mfem::Array<int>> mg_submesh_ess_tdof_lists(
      submesh_nd_fespaces.GetNumLevels());
  for (std::size_t l = 0; l < submesh_nd_fespaces.GetNumLevels(); l++)
  {
    const auto &nd_fespace_l = submesh_nd_fespaces.GetFESpaceAtLevel(l);
    nd_fespace_l.Get().GetEssentialTrueDofs(combined_inner_bdr_marker,
                                            mg_submesh_ess_tdof_lists[l]);
  }

  // Get finest level essential DoF list for boundary conditions
  auto submesh_ess_tdof_list =
      mg_submesh_ess_tdof_lists[submesh_nd_fespaces.GetNumLevels() - 1];

  // First synchronize the marker itself to ensure all processors agree on which DoFs to
  // sync
  auto gc = std::unique_ptr<mfem::GroupCommunicator>(
      submesh_nd_fespaces.GetFinestFESpace().Get().ScalarGroupComm());
  mfem::Array<int> global_marker(ldof_marker_submesh);
  gc->Reduce<int>(global_marker.GetData(), mfem::GroupCommunicator::BitOR<int>);
  gc->Bcast(global_marker);

  // Synchronize the edge DoF boundary condition across processors with MaxAbs reduction
  mfem::Array<double> values(A.GetData(), A.Size());
  gc->ReduceBegin(values.GetData());
  gc->ReduceMarked<double>(values.GetData(), global_marker, 0,
                           mfem::GroupCommunicator::MaxAbs<double>);
  gc->Bcast(values.GetData());
  A.HostReadWrite();
  A.SetTrueVector();

  // Create unit coefficient for curl-curl term and the regularization coefficient
  MaterialPropertyCoefficient reg_coeff(boundary_submesh.attributes.Max());
  MaterialPropertyCoefficient unit_coeff(boundary_submesh.attributes.Max());
  for (int attr = 1; attr <= boundary_submesh.attributes.Max(); attr++)
  {
    unit_coeff.AddMaterialProperty(attr, 1.0);
    reg_coeff.AddMaterialProperty(attr, flux_data.regularization);
  }

  // Use Palace's BilinearForm with forced full assembly for the submesh solve.
  // The 2D submesh is small and must stay on host (GroupCommunicator and direct DOF
  // indexing are incompatible with device memory).
  BilinearForm a(submesh_nd_fespaces.GetFinestFESpace());
  a.AddDomainIntegrator<CurlCurlIntegrator>(unit_coeff);
  a.AddDomainIntegrator<VectorFEMassIntegrator>(reg_coeff);

  std::vector<std::unique_ptr<Operator>> k_vec;
  k_vec.reserve(submesh_nd_fespaces.GetNumLevels());
  for (std::size_t l = 0; l < submesh_nd_fespaces.GetNumLevels(); l++)
  {
    BilinearForm a_l(submesh_nd_fespaces.GetFESpaceAtLevel(l));
    a_l.AddDomainIntegrator<CurlCurlIntegrator>(unit_coeff);
    a_l.AddDomainIntegrator<VectorFEMassIntegrator>(reg_coeff);
    k_vec.push_back(a_l.FullAssemble(false));
  }
  auto K_op = std::make_unique<MultigridOperator>(submesh_nd_fespaces.GetNumLevels());

  // Add operators for each level using pre-computed essential DoF lists
  for (std::size_t l = 0; l < submesh_nd_fespaces.GetNumLevels(); l++)
  {
    const auto &nd_fespace_l = submesh_nd_fespaces.GetFESpaceAtLevel(l);
    auto K_l = std::make_unique<ParOperator>(std::move(k_vec[l]), nd_fespace_l);
    K_l->SetEssentialTrueDofs(mg_submesh_ess_tdof_lists[l],
                              Operator::DiagonalPolicy::DIAG_ONE);
    K_op->AddOperator(std::move(K_l));
  }

  // Compute boundary-interior coupling before applying boundary conditions.
  // Keep all submesh vectors on host — the 2D solve is small and uses host-side
  // operations (GroupCommunicator, direct DOF indexing) incompatible with device memory.
  Vector RHS(submesh_nd_fespaces.GetFinestFESpace().GetTrueVSize());
  Vector X(submesh_nd_fespaces.GetFinestFESpace().GetTrueVSize());
  Vector boundary_vals(submesh_nd_fespaces.GetFinestFESpace().GetTrueVSize());

  // Get boundary values directly from MFEM GridFunction
  A.GetTrueDofs(boundary_vals);

  // Compute RHS = -K * boundary_values (coupling term)
  K_op->Mult(boundary_vals, RHS);
  RHS *= -1.0;

  // Set boundary values in RHS
  linalg::SetSubVector(RHS, submesh_ess_tdof_list, boundary_vals);

  // Set up Palace KSP solver
  KspSolver ksp(iodata, submesh_nd_fespaces, &submesh_h1_fespaces);
  ksp.SetOperators(*K_op, *K_op);

  // Solve the 2D surface curl problem.
  X = 0.0;
  ksp.Mult(RHS, X);

  // Set solution directly in MFEM GridFunction.
  A.SetFromTrueDofs(X);

  // Transfer to parent mesh
  auto &A_3d = post_op.GetAGridFunction().Real();
  A_3d = 0.0;  // Clear the buffer
  mfem::ParSubMesh::Transfer(A, A_3d);

  // Extract true DOFs and populate result vector
  result.SetSize(fespace->GetTrueVSize());
  result.UseDevice(true);
  A_3d.GetTrueDofs(result);
}

double ComputeFluxThroughSurface(const mfem::ParGridFunction &B_gf,
                                 const std::vector<int> &attributes, const Mesh &mesh,
                                 const MaterialOperator &mat_op,
                                 const mfem::Vector &flux_direction, MPI_Comm comm)
{
  if (attributes.empty())
  {
    return 0.0;
  }

  // Ensure face neighbor data is available for shared boundary evaluation.
  const_cast<mfem::ParGridFunction &>(B_gf).ExchangeFaceNbrData();

  // Create magnetic flux coefficient with direction-based orientation
  BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC> flux_coeff(
      nullptr, &B_gf, mat_op, false, flux_direction,
      BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC>::OrientationMode::DIRECTION_BASED);

  mfem::ParMesh *pmesh = const_cast<mfem::ParMesh *>(&mesh.Get());
  mfem::ParFiniteElementSpace *fes =
      const_cast<mfem::ParFiniteElementSpace *>(B_gf.ParFESpace());

  // Integrate flux_coeff over the selected boundaries
  double local_flux = 0.0;
  int nbdr = pmesh->GetNBE();
  // Precompute integration rules for each geometry/order combination
  std::map<std::pair<mfem::Geometry::Type, int>, const mfem::IntegrationRule *> ir_map;
  for (int be = 0; be < nbdr; ++be)
  {
    int attr = pmesh->GetBdrAttribute(be);
    if (std::find(attributes.begin(), attributes.end(), attr) == attributes.end())
    {
      continue;
    }
    const mfem::FiniteElement *fe = fes->GetBE(be);
    mfem::ElementTransformation *Tr = pmesh->GetBdrElementTransformation(be);
    const int int_order = fem::DefaultIntegrationOrder::Get(*Tr);
    auto geom_order = std::make_pair(fe->GetGeomType(), int_order);
    const mfem::IntegrationRule *ir;
    auto it = ir_map.find(geom_order);
    if (it == ir_map.end())
    {
      ir = &mfem::IntRules.Get(fe->GetGeomType(), int_order);
      ir_map[geom_order] = ir;
    }
    else
    {
      ir = it->second;
    }
    for (int q = 0; q < ir->GetNPoints(); ++q)
    {
      const mfem::IntegrationPoint &ip = ir->IntPoint(q);
      Tr->SetIntPoint(&ip);
      double val = flux_coeff.Eval(*Tr, ip);
      local_flux += val * ip.weight * Tr->Weight();
    }
  }
  double computed_flux = local_flux;
  Mpi::GlobalSum(1, &computed_flux, comm);
  return computed_flux;
}

void VerifyFluxThroughHoles(const mfem::ParGridFunction &B_gf,
                            const std::vector<int> &hole_attributes,
                            const std::vector<double> &target_fluxes, const Mesh &mesh,
                            const MaterialOperator &mat_op,
                            const mfem::Vector &flux_direction, MPI_Comm comm)
{
  for (std::size_t h = 0; h < hole_attributes.size(); h++)
  {
    int hole_attr = hole_attributes[h];
    double target_flux = target_fluxes[h];

    double computed_flux =
        ComputeFluxThroughSurface(B_gf, {hole_attr}, mesh, mat_op, flux_direction, comm);

    if (Mpi::Root(comm))
    {
      Mpi::Print("Hole attribute {:d}: Target flux = {:.6e}, Computed flux = {:.6e}, Error "
                 "= {:.6e}\n",
                 hole_attr, target_flux, computed_flux,
                 std::abs(computed_flux - target_flux));
    }
  }
}

void VerifyFluxThroughAllHoles(const mfem::ParGridFunction &B_gf, const IoData &iodata,
                               int current_flux_loop_idx, const Mesh &mesh,
                               const MaterialOperator &mat_op, MPI_Comm comm)
{
  if (Mpi::Root(comm))
  {
    Mpi::Print("FluxLoop {:d} excitation - Flux through all holes:\n",
               current_flux_loop_idx);
  }

  // Compute flux through all holes in all flux loops
  for (const auto &[loop_idx, flux_data] : iodata.boundaries.fluxloop)
  {
    for (std::size_t h = 0; h < flux_data.hole_attributes.size(); h++)
    {
      int hole_attr = flux_data.hole_attributes[h];
      double target_flux =
          (loop_idx == current_flux_loop_idx) ? flux_data.flux_amounts[h] : 0.0;

      mfem::Vector flux_direction(const_cast<double *>(flux_data.direction.data()), 3);
      double computed_flux =
          ComputeFluxThroughSurface(B_gf, {hole_attr}, mesh, mat_op, flux_direction, comm);

      Mpi::Print(
          "  Loop {:d} Hole {:d}: Target = {:.6e}, Computed = {:.6e}, Error = {:.6e}\n",
          loop_idx, hole_attr, target_flux, computed_flux,
          std::abs(computed_flux - target_flux));
    }
  }
}

}  // namespace palace
