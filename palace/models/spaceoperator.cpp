// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "spaceoperator.hpp"

#include <limits>
#include <set>
#include <string_view>
#include <type_traits>
#include <fmt/format.h>
#include "drivers/singularsolver.hpp"
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "fem/singulardofs.hpp"
#include "fem/singularsystem.hpp"
#include "linalg/hypre.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/ksp.hpp"
#include "linalg/rap.hpp"
#include "models/floquetportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

struct SparseMatrixStatistics
{
  HYPRE_BigInt stored = 0;
  HYPRE_BigInt exact_zeros = 0;
  HYPRE_BigInt small_14 = 0;
  HYPRE_BigInt small_12 = 0;
  HYPRE_BigInt nonfinite = 0;
  double maximum_absolute_value = 0.0;
};

SparseMatrixStatistics GetSparseMatrixStatistics(const mfem::HypreParMatrix &matrix)
{
  mfem::SparseMatrix diagonal, off_diagonal;
  HYPRE_BigInt *off_diagonal_columns = nullptr;
  matrix.GetDiag(diagonal);
  matrix.GetOffd(off_diagonal, off_diagonal_columns);

  SparseMatrixStatistics statistics;
  const auto inspect_maximum = [&](const mfem::SparseMatrix &block)
  {
    const auto *values = block.HostReadData();
    const int number_entries = block.NumNonZeroElems();
    statistics.stored += number_entries;
    for (int entry = 0; entry < number_entries; entry++)
    {
      if (!std::isfinite(values[entry]))
      {
        statistics.nonfinite++;
      }
      else
      {
        statistics.maximum_absolute_value =
            std::max(statistics.maximum_absolute_value, std::abs(values[entry]));
      }
    }
  };
  inspect_maximum(diagonal);
  inspect_maximum(off_diagonal);
  Mpi::GlobalMax(1, &statistics.maximum_absolute_value, matrix.GetComm());

  const double threshold_14 = 1.0e-14 * statistics.maximum_absolute_value;
  const double threshold_12 = 1.0e-12 * statistics.maximum_absolute_value;
  const auto inspect_small = [&](const mfem::SparseMatrix &block)
  {
    const auto *values = block.HostReadData();
    const int number_entries = block.NumNonZeroElems();
    for (int entry = 0; entry < number_entries; entry++)
    {
      const double magnitude = std::abs(values[entry]);
      statistics.exact_zeros += values[entry] == 0.0;
      statistics.small_14 += magnitude <= threshold_14;
      statistics.small_12 += magnitude <= threshold_12;
    }
  };
  inspect_small(diagonal);
  inspect_small(off_diagonal);
  Mpi::GlobalSum(1, &statistics.stored, matrix.GetComm());
  Mpi::GlobalSum(1, &statistics.exact_zeros, matrix.GetComm());
  Mpi::GlobalSum(1, &statistics.small_14, matrix.GetComm());
  Mpi::GlobalSum(1, &statistics.small_12, matrix.GetComm());
  Mpi::GlobalSum(1, &statistics.nonfinite, matrix.GetComm());
  return statistics;
}

void PrintSparseMatrixStatistics(std::string_view label, const mfem::HypreParMatrix &matrix)
{
  const auto statistics = GetSparseMatrixStatistics(matrix);
  Mpi::Print(" Singular sparse audit, {}: {:d} stored, {:d} exact zero, {:d} <= "
             "1e-14 max, {:d} <= 1e-12 max, {:d} nonfinite, max |a| = {:.3e}\n",
             label, statistics.stored, statistics.exact_zeros, statistics.small_14,
             statistics.small_12, statistics.nonfinite, statistics.maximum_absolute_value);
}

template <std::size_t NumVertices>
std::vector<std::array<fem::singular::GlobalVertexId, NumVertices>>
GatherBoundaryEntities(const mfem::ParMesh &mesh,
                       const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
                       const std::set<int> &boundary_attributes,
                       mfem::Geometry::Type geometry)
{
  if (boundary_attributes.empty())
  {
    return {};
  }
  MFEM_VERIFY(source_vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()),
              "Singular boundary classification requires complete source vertex IDs!");

  std::set<std::array<fem::singular::GlobalVertexId, NumVertices>> local_entities;
  for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
  {
    if (boundary_attributes.count(mesh.GetBdrAttribute(boundary)) == 0)
    {
      continue;
    }
    const auto &element = *mesh.GetBdrElement(boundary);
    MFEM_VERIFY(element.GetGeometryType() == geometry &&
                    element.GetNVertices() == static_cast<int>(NumVertices),
                "Singular boundary classification received an unexpected boundary "
                "element geometry!");
    std::array<fem::singular::GlobalVertexId, NumVertices> entity;
    for (std::size_t vertex = 0; vertex < NumVertices; vertex++)
    {
      const int mesh_vertex = element.GetVertices()[vertex];
      MFEM_VERIFY(mesh_vertex >= 0 && mesh_vertex < mesh.GetNV(),
                  "A singular boundary entity has an invalid vertex!");
      entity[vertex] = source_vertex_ids[mesh_vertex];
    }
    std::sort(entity.begin(), entity.end());
    MFEM_VERIFY(entity[0] >= 0 &&
                    std::adjacent_find(entity.begin(), entity.end()) == entity.end(),
                "A singular boundary entity has invalid source vertices!");
    local_entities.insert(entity);
  }

  std::vector<fem::singular::GlobalVertexId> local_data;
  local_data.reserve(NumVertices * local_entities.size());
  for (const auto &entity : local_entities)
  {
    local_data.insert(local_data.end(), entity.begin(), entity.end());
  }
  MFEM_VERIFY(local_data.size() <=
                  static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "Too many local singular boundary entities!");
  const int local_count = static_cast<int>(local_data.size());
  std::vector<int> counts(Mpi::Size(mesh.GetComm()));
  Mpi::Allgather(1, &local_count, counts.data(), mesh.GetComm());
  std::vector<int> offsets(counts.size());
  int global_count = 0;
  for (std::size_t rank = 0; rank < counts.size(); rank++)
  {
    MFEM_VERIFY(counts[rank] >= 0 && counts[rank] % static_cast<int>(NumVertices) == 0 &&
                    counts[rank] <= std::numeric_limits<int>::max() - global_count,
                "Invalid distributed singular boundary entity count!");
    offsets[rank] = global_count;
    global_count += counts[rank];
  }
  std::vector<fem::singular::GlobalVertexId> global_data(global_count);
  Mpi::Allgatherv(local_count, local_data.data(), global_data.data(), counts.data(),
                  offsets.data(), mesh.GetComm());

  std::set<std::array<fem::singular::GlobalVertexId, NumVertices>> global_entities;
  for (int offset = 0; offset < global_count; offset += NumVertices)
  {
    std::array<fem::singular::GlobalVertexId, NumVertices> entity;
    std::copy_n(global_data.begin() + offset, NumVertices, entity.begin());
    MFEM_VERIFY(entity[0] >= 0 && std::is_sorted(entity.begin(), entity.end()) &&
                    std::adjacent_find(entity.begin(), entity.end()) == entity.end(),
                "A gathered singular boundary entity is invalid!");
    global_entities.insert(entity);
  }
  return {global_entities.begin(), global_entities.end()};
}

}  // namespace

SpaceOperator::SpaceOperator(const config::SolverData &solver,
                             const config::DomainData &domains,
                             const config::BoundaryData &boundaries,
                             ProblemType problem_type, const Units &units,
                             const std::vector<std::unique_ptr<Mesh>> &mesh)
  : pc_mat_real(solver.linear.pc_mat_real), pc_mat_shifted(solver.linear.pc_mat_shifted),
    print_hdr(true), print_prec_hdr(true),
    dbc_attr(SetUpBoundaryProperties(boundaries.pec, *mesh.back())),
    nd_fecs(fem::ConstructFECollections<mfem::ND_FECollection>(
        solver.order, mesh.back()->Dimension(), solver.linear.mg_max_levels,
        solver.linear.mg_coarsening, false)),
    h1_fecs(fem::ConstructFECollections<mfem::H1_FECollection>(
        solver.order, mesh.back()->Dimension(), solver.linear.mg_max_levels,
        solver.linear.mg_coarsening, false)),
    rt_fecs(fem::ConstructFECollections<mfem::RT_FECollection>(
        solver.order - 1, mesh.back()->Dimension(),
        solver.linear.estimator_mg ? solver.linear.mg_max_levels : 1,
        solver.linear.mg_coarsening, false)),
    nd_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        solver.linear.mg_max_levels, mesh, nd_fecs, &dbc_attr, &nd_dbc_tdof_lists)),
    h1_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        solver.linear.mg_max_levels, mesh, h1_fecs, &dbc_attr, &h1_dbc_tdof_lists)),
    rt_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::RT_FECollection>(
        solver.linear.estimator_mg ? solver.linear.mg_max_levels : 1, mesh, rt_fecs)),
    mat_op(domains.materials, boundaries.periodic, problem_type, *mesh.back()),
    current_dipole_op(domains.current_dipole, units, *mesh.back()),
    farfield_op(boundaries.farfield, problem_type, mat_op, *mesh.back()),
    surf_sigma_op(boundaries.conductivity, problem_type, units, mat_op, *mesh.back()),
    surf_z_op(boundaries.impedance, boundaries.cracked_attributes, units, mat_op,
              *mesh.back()),
    surf_rz_op(boundaries.rational_impedance, boundaries.cracked_attributes, problem_type,
               units, mat_op, *mesh.back()),
    lumped_port_op(boundaries.lumpedport, units, mat_op, *mesh.back()),
    wave_port_op(boundaries, domains, solver, problem_type, units, mat_op, GetNDSpace(),
                 GetH1Space()),
    floquet_port_op(boundaries.floquetport, boundaries.periodic, problem_type, units,
                    mat_op, GetNDSpace().Get()),
    surf_j_op(boundaries.current, *mesh.back()),
    port_excitation_helper(lumped_port_op, wave_port_op, floquet_port_op, surf_j_op,
                           current_dipole_op)
{
  // In 2D, curl maps H(curl) → L2 (scalar), so we need an L2 FE space for B = curl E.
  // Must use INTEGRAL map type so the discrete interpolator recognizes this as the curl
  // target space.
  if (mesh.back()->Dimension() == 2)
  {
    const int l2_order = solver.order - 1;
    const int dim = mesh.back()->Dimension();
    l2_curl_fecs.push_back(std::make_unique<mfem::L2_FECollection>(
        l2_order, dim, mfem::BasisType::GaussLegendre, mfem::FiniteElement::INTEGRAL));
    l2_curl_fespaces = std::make_unique<FiniteElementSpaceHierarchy>(
        fem::ConstructFiniteElementSpaceHierarchy<mfem::L2_FECollection>(1, mesh,
                                                                         l2_curl_fecs));
  }

  // Finalize setup.
  CheckBoundaryProperties();
  combined_nd_dbc_tdof_lists = nd_dbc_tdof_lists;
  combined_h1_dbc_tdof_lists = h1_dbc_tdof_lists;
  combined_h1_aux_bdr_tdof_lists = aux_bdr_tdof_lists;

  // Print essential BC information.
  if (dbc_attr.Size())
  {
    Mpi::Print("\nConfiguring Dirichlet PEC BC at attributes:\n");
    utils::PrettyPrint(dbc_attr);
  }
}

SpaceOperator::SpaceOperator(
    const IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh,
    const fem::singular::FeatureTopology *singular_features_in,
    const fem::singular::TriangleFeatureTopology *triangle_singular_features_in,
    const std::vector<fem::singular::GlobalVertexId> *source_vertex_ids_in)
  : SpaceOperator(iodata.solver, iodata.domains, iodata.boundaries, iodata.problem.type,
                  iodata.units, mesh)
{
  singular_features = singular_features_in;
  triangle_singular_features = triangle_singular_features_in;
  source_vertex_ids = source_vertex_ids_in;
  if (singular_features || triangle_singular_features)
  {
    SetUpSingularEnrichment(iodata);
    CheckSingularExcitations(iodata.problem.type);
  }

  // Validate excitations after wave port setup is complete.
  CheckExcitations(iodata.problem.type);
}

void SpaceOperator::SetUpSingularEnrichment(const IoData &iodata)
{
  const auto &solver = iodata.solver;
  const auto &boundaries = iodata.boundaries;
  const bool tetrahedral = singular_features != nullptr;
  const bool triangular = triangle_singular_features != nullptr;
  MFEM_VERIFY(tetrahedral != triangular && source_vertex_ids &&
                  source_vertex_ids->size() ==
                      static_cast<std::size_t>(GetMesh().Get().GetNV()) &&
                  solver.singular_elements.Enabled(),
              "Full-wave singular enrichment requires exactly one complete simplex "
              "feature topology and source vertex IDs!");
  MFEM_VERIFY(
      (tetrahedral && GetMesh().Dimension() == 3) ||
          (triangular && GetMesh().Dimension() == 2 && solver.singular_elements.order == 1),
      "Full-wave singular feature topology is inconsistent with the solve mesh!");
  const auto number_levels = GetNDSpaces().GetNumLevels();
  MFEM_VERIFY(number_levels == GetH1Spaces().GetNumLevels(),
              "Full-wave singular ND and H1 hierarchies have inconsistent levels!");
  for (std::size_t level = 0; level < number_levels; level++)
  {
    MFEM_VERIFY(&GetNDSpaces().GetFESpaceAtLevel(level).GetMesh() == &GetMesh() &&
                    &GetH1Spaces().GetFESpaceAtLevel(level).GetMesh() == &GetMesh(),
                "Full-wave singular multigrid currently supports polynomial levels on "
                "one mesh; singular AMR topology transfer remains unsupported!");
  }

  Mpi::Print("\nSetting up singular full-wave enrichment:\n");
  const auto setup_start = Timer::Now();
  auto stage_start = setup_start;
  const auto report_stage = [&](std::string_view stage)
  {
    const auto now = Timer::Now();
    double elapsed = Timer::Duration(now - stage_start).count();
    double minimum = elapsed;
    double maximum = elapsed;
    Mpi::GlobalMin(1, &minimum, GetComm());
    Mpi::GlobalMax(1, &maximum, GetComm());
    Mpi::GlobalSum(1, &elapsed, GetComm());
    Mpi::Print(" Singular setup timing, {} (s): min. {:.3f}, max. {:.3f}, avg. {:.3f}\n",
               stage, minimum, maximum, elapsed / Mpi::Size(GetComm()));
    stage_start = Timer::Now();
  };

  std::set<int> lumped_boundary_attributes;
  const auto collect_lumped_attributes =
      [&lumped_boundary_attributes](const std::map<int, double> &coefficients)
  {
    for (const auto &[attribute, coefficient] : coefficients)
    {
      if (coefficient != 0.0)
      {
        lumped_boundary_attributes.insert(attribute);
      }
    }
  };
  collect_lumped_attributes(lumped_port_op.GetStiffnessBdrCoefficientMap());
  collect_lumped_attributes(lumped_port_op.GetDampingBdrCoefficientMap());
  collect_lumped_attributes(lumped_port_op.GetMassBdrCoefficientMap());
  const auto excluded_boundary_trace_faces =
      tetrahedral
          ? GatherBoundaryEntities<3>(GetMesh().Get(), *source_vertex_ids,
                                      lumped_boundary_attributes, mfem::Geometry::TRIANGLE)
          : std::vector<std::array<fem::singular::GlobalVertexId, 3>>{};
  if (!excluded_boundary_trace_faces.empty())
  {
    Mpi::Print(" Excluding nonintegrable edge-gradient traces on {:d} lumped-port "
               "face{}\n",
               excluded_boundary_trace_faces.size(),
               excluded_boundary_trace_faces.size() == 1 ? "" : "s");
  }

  if (tetrahedral)
  {
    singular_dofs =
        std::make_unique<fem::singular::DofTopology>(fem::singular::BuildLocalDofTopology(
            GetMesh(), *singular_features, *source_vertex_ids,
            solver.singular_elements.order, excluded_boundary_trace_faces));
    singular_numbering = std::make_unique<fem::singular::ParallelDofNumbering>(
        fem::singular::BuildParallelDofNumbering(GetComm(), *singular_dofs));
  }
  else
  {
    triangle_singular_dofs = std::make_unique<fem::singular::TriangleDofTopology>(
        fem::singular::BuildLocalTriangleDofTopology(GetMesh(), *triangle_singular_features,
                                                     *source_vertex_ids,
                                                     solver.singular_elements.order));
    singular_numbering = std::make_unique<fem::singular::ParallelDofNumbering>(
        fem::singular::BuildParallelDofNumbering(GetComm(), *triangle_singular_dofs));
  }
  report_stage("DOF topology and parallel numbering");
  MFEM_VERIFY(singular_numbering->h1.owned_size <= std::numeric_limits<int>::max() &&
                  singular_numbering->nd.owned_size <= std::numeric_limits<int>::max(),
              "Full-wave singular local true-DOF count exceeds integer limits!");
  singular_nd_gradient_true_dofs.SetSize(
      static_cast<int>(singular_numbering->h1.owned_size));
  int gradient_dof = 0;
  for (std::size_t local_h1 = 0; local_h1 < singular_numbering->h1.local_to_true.size();
       local_h1++)
  {
    if (singular_numbering->h1.owner[local_h1] != Mpi::Rank(GetComm()))
    {
      continue;
    }
    const HYPRE_BigInt nd_true = singular_numbering->h1_to_nd_true[local_h1];
    MFEM_VERIFY(nd_true >= singular_numbering->nd.owned_offset &&
                    nd_true < singular_numbering->nd.owned_offset +
                                  singular_numbering->nd.owned_size,
                "Owned singular H1 and gradient ND true DOFs must share one MPI owner!");
    singular_nd_gradient_true_dofs[gradient_dof++] =
        static_cast<int>(nd_true - singular_numbering->nd.owned_offset);
  }
  MFEM_VERIFY(gradient_dof == singular_nd_gradient_true_dofs.Size(),
              "Singular gradient true-DOF list has an inconsistent local size!");
  singular_nd_gradient_true_dofs.Sort();
  singular_nd_gradient_true_dofs.Unique();
  MFEM_VERIFY(gradient_dof == singular_nd_gradient_true_dofs.Size(),
              "Singular gradient true-DOF list contains duplicate entries!");

  long long local_enriched_elements = 0;
  long long local_h1_basis_incidences = 0;
  long long local_nd_basis_incidences = 0;
  const auto count_element_dofs = [&](const auto &topology)
  {
    for (const auto &element : topology.elements)
    {
      if (!element.h1.empty() || !element.nd.empty())
      {
        local_enriched_elements++;
        local_h1_basis_incidences += static_cast<long long>(element.h1.size());
        local_nd_basis_incidences += static_cast<long long>(element.nd.size());
      }
    }
  };
  if (tetrahedral)
  {
    count_element_dofs(*singular_dofs);
  }
  else
  {
    count_element_dofs(*triangle_singular_dofs);
  }
  long long global_enriched_elements = local_enriched_elements;
  long long global_h1_basis_incidences = local_h1_basis_incidences;
  long long global_nd_basis_incidences = local_nd_basis_incidences;
  long long minimum_enriched_elements = local_enriched_elements;
  long long maximum_enriched_elements = local_enriched_elements;
  long long minimum_h1_basis_incidences = local_h1_basis_incidences;
  long long maximum_h1_basis_incidences = local_h1_basis_incidences;
  long long minimum_nd_basis_incidences = local_nd_basis_incidences;
  long long maximum_nd_basis_incidences = local_nd_basis_incidences;
  Mpi::GlobalMin(1, &minimum_enriched_elements, GetComm());
  Mpi::GlobalMax(1, &maximum_enriched_elements, GetComm());
  Mpi::GlobalMin(1, &minimum_h1_basis_incidences, GetComm());
  Mpi::GlobalMax(1, &maximum_h1_basis_incidences, GetComm());
  Mpi::GlobalMin(1, &minimum_nd_basis_incidences, GetComm());
  Mpi::GlobalMax(1, &maximum_nd_basis_incidences, GetComm());
  Mpi::GlobalSum(1, &global_enriched_elements, GetComm());
  Mpi::GlobalSum(1, &global_h1_basis_incidences, GetComm());
  Mpi::GlobalSum(1, &global_nd_basis_incidences, GetComm());
  Mpi::Print(" Enriched elements: {:d} global; basis incidences: {:d} H1, {:d} ND\n",
             global_enriched_elements, global_h1_basis_incidences,
             global_nd_basis_incidences);
  Mpi::Print(" Singular rank workload: enriched elements {:d}-{:d}; H1 incidences "
             "{:d}-{:d}; ND incidences {:d}-{:d}\n",
             minimum_enriched_elements, maximum_enriched_elements,
             minimum_h1_basis_incidences, maximum_h1_basis_incidences,
             minimum_nd_basis_incidences, maximum_nd_basis_incidences);

  singular_material_coefficients.assign(GetMesh().GetNE(), {1.0, 1.0});
  singular_imag_material_coefficients.assign(GetMesh().GetNE(), {0.0, 1.0});
  singular_abs_material_coefficients.assign(GetMesh().GetNE(), {1.0, 1.0});
  auto &materials = singular_material_coefficients;
  auto &imaginary_materials = singular_imag_material_coefficients;
  auto &absolute_materials = singular_abs_material_coefficients;
  for (int element = 0; element < GetMesh().GetNE(); element++)
  {
    const int attribute = GetMesh().Get().GetAttribute(element);
    const bool enriched =
        tetrahedral ? !singular_features->elements[element].nodes.empty() ||
                          !singular_features->elements[element].edges.empty()
                    : !triangle_singular_features->elements[element].nodes.empty();
    if (!enriched)
    {
      continue;
    }
    MFEM_VERIFY(mat_op.IsIsotropic(attribute),
                "Full-wave singular enrichment requires isotropic material in every "
                "enriched simplex! Domain attribute: "
                    << attribute);
    materials[element] = {mat_op.GetPermittivityReal(attribute)(0, 0),
                          mat_op.GetInvPermeability(attribute)(0, 0)};
    imaginary_materials[element] = {mat_op.GetPermittivityImag(attribute)(0, 0), 1.0};
    absolute_materials[element] = {mat_op.GetPermittivityAbs(attribute)(0, 0), 1.0};
  }

  const fem::singular::AdaptiveAssemblyOptions options{
      solver.singular_elements.quadrature_order, solver.singular_elements.abs_tol,
      solver.singular_elements.rel_tol, solver.singular_elements.max_subdivisions};
  const auto constrained_impedance_attributes =
      tetrahedral
          ? GetConstrainedSingularImpedanceAttributes(iodata, *singular_features)
          : GetConstrainedSingularImpedanceAttributes(iodata, *triangle_singular_features);
  const auto filter_impedance_coefficients = [&](std::map<int, double> coefficients)
  {
    for (int attribute : constrained_impedance_attributes)
    {
      coefficients.erase(attribute);
    }
    return coefficients;
  };
  const auto impedance_stiffness_coefficients =
      filter_impedance_coefficients(surf_z_op.GetStiffnessBdrCoefficientMap());
  const auto impedance_damping_coefficients =
      filter_impedance_coefficients(surf_z_op.GetDampingBdrCoefficientMap());
  const auto impedance_mass_coefficients =
      filter_impedance_coefficients(surf_z_op.GetMassBdrCoefficientMap());
  singular_assembly_options = options;
  singular_domain_matrices.resize(number_levels);
  singular_domain_imag_matrices.resize(number_levels);
  singular_domain_abs_matrices.resize(number_levels);
  singular_lumped_stiffness_matrices.resize(number_levels);
  singular_lumped_damping_matrices.resize(number_levels);
  singular_lumped_mass_matrices.resize(number_levels);
  singular_impedance_stiffness_matrices.resize(number_levels);
  singular_impedance_damping_matrices.resize(number_levels);
  singular_impedance_mass_matrices.resize(number_levels);
  singular_gradients.reserve(number_levels);
  report_stage("material and boundary coefficient preparation");

  const std::size_t finest_level = number_levels - 1;
  auto &finest_h1_space = GetH1Spaces().GetFESpaceAtLevel(finest_level);
  auto &finest_nd_space = GetNDSpaces().GetFESpaceAtLevel(finest_level);
  std::vector<std::vector<fem::singular::IsotropicMaterialCoefficients>> material_batches{
      materials, absolute_materials};
  if (mat_op.HasLossTangent())
  {
    material_batches.push_back(imaginary_materials);
  }
  Mpi::Print(" Singular setup finest level {:d}/{:d}: assembling local domain tensors "
             "and sparse blocks for {:d} material batch{}\n",
             finest_level + 1, number_levels, material_batches.size(),
             material_batches.size() == 1 ? "" : "es");
  auto local_matrices = tetrahedral
                            ? fem::singular::AssembleLocalSparseEnrichmentMatricesBatch(
                                  *singular_dofs, finest_h1_space.Get(),
                                  finest_nd_space.Get(), material_batches, options)
                            : fem::singular::AssembleLocalSparseEnrichmentMatricesBatch(
                                  *triangle_singular_dofs, finest_h1_space.Get(),
                                  finest_nd_space.Get(), material_batches, options);
  HYPRE_BigInt reference_entries =
      static_cast<HYPRE_BigInt>(local_matrices.front().affine_reference_table_entries);
  HYPRE_BigInt reference_hits =
      static_cast<HYPRE_BigInt>(local_matrices.front().affine_reference_cache_hits);
  Mpi::GlobalMax(1, &reference_entries, GetComm());
  Mpi::GlobalSum(1, &reference_hits, GetComm());
  Mpi::Print(" Singular affine reference cache: {:d} entries/rank maximum, {:d} "
             "global hits\n",
             reference_entries, reference_hits);
  HYPRE_BigInt affine_contraction_counts[6] = {
      static_cast<HYPRE_BigInt>(local_matrices.front().affine_nd_mass_contraction_count),
      static_cast<HYPRE_BigInt>(local_matrices.front().affine_nd_mass_reintegration_count),
      static_cast<HYPRE_BigInt>(
          local_matrices.front().affine_nd_mass_reintegration_batch_count),
      static_cast<HYPRE_BigInt>(local_matrices.front().affine_nd_curl_contraction_count),
      static_cast<HYPRE_BigInt>(local_matrices.front().affine_nd_curl_reintegration_count),
      static_cast<HYPRE_BigInt>(
          local_matrices.front().affine_nd_curl_reintegration_batch_count)};
  Mpi::GlobalSum(6, affine_contraction_counts, GetComm());
  Mpi::Print(" Singular affine ND coupling: {:d} mass contractions, {:d} mass "
             "reintegrations in {:d} batches; {:d} curl contractions, {:d} curl "
             "reintegrations in {:d} batches\n",
             affine_contraction_counts[0], affine_contraction_counts[1],
             affine_contraction_counts[2], affine_contraction_counts[3],
             affine_contraction_counts[4], affine_contraction_counts[5]);
  const auto report_local_assembly_stage = [&](std::string_view stage, double elapsed)
  {
    double minimum = elapsed;
    double maximum = elapsed;
    Mpi::GlobalMin(1, &minimum, GetComm());
    Mpi::GlobalMax(1, &maximum, GetComm());
    Mpi::GlobalSum(1, &elapsed, GetComm());
    Mpi::Print(" Singular local assembly timing, {} (s): min. {:.3f}, max. {:.3f}, "
               "avg. {:.3f}\n",
               stage, minimum, maximum, elapsed / Mpi::Size(GetComm()));
  };
  report_local_assembly_stage("enrichment tensor evaluation",
                              local_matrices.front().enrichment_evaluation_time);
  report_local_assembly_stage("standard-enrichment coupling evaluation",
                              local_matrices.front().standard_enrichment_evaluation_time);
  report_local_assembly_stage("standard-reference generation",
                              local_matrices.front().standard_reference_generation_time);
  report_local_assembly_stage("coupling setup",
                              local_matrices.front().standard_enrichment_setup_time);
  report_local_assembly_stage("ND coupling contraction",
                              local_matrices.front().nd_coupling_time);
  report_local_assembly_stage("H1 exact-gradient coupling",
                              local_matrices.front().h1_gradient_coupling_time);
  report_local_assembly_stage("H1 mass coupling",
                              local_matrices.front().h1_mass_coupling_time);
  report_local_assembly_stage("material copies and transformations",
                              local_matrices.front().material_transformation_time);
  report_local_assembly_stage("sparse insertion",
                              local_matrices.front().sparse_insertion_time);
  report_local_assembly_stage("sparse finalization",
                              local_matrices.front().sparse_finalization_time);
  report_stage("finest-level local domain assembly");
  MFEM_VERIFY(local_matrices.size() == material_batches.size(),
              "Full-wave singular material batch assembly returned an inconsistent "
              "number of operators!");
  const auto assemble_domain = [&](std::size_t batch)
  {
    return fem::singular::AssembleParallelSparseEnrichmentMatrices(
        local_matrices[batch], *singular_numbering, finest_h1_space.Get(),
        finest_nd_space.Get());
  };
  const auto compact_blocks = [](fem::singular::ParallelSparseOperatorBlocks &blocks)
  {
    if (!blocks.standard_enrichment)
    {
      MFEM_VERIFY(!blocks.enrichment_standard && !blocks.enrichment_enrichment,
                  "Cannot compact incomplete singular operator blocks!");
      return HYPRE_BigInt{0};
    }
    MFEM_VERIFY(blocks.enrichment_standard && blocks.enrichment_enrichment,
                "Cannot compact incomplete singular operator blocks!");
    HYPRE_BigInt removed = 0;
    removed += fem::singular::RemoveExplicitZeros(*blocks.standard_enrichment);
    removed += fem::singular::RemoveExplicitZeros(*blocks.enrichment_standard);
    removed += fem::singular::RemoveExplicitZeros(*blocks.enrichment_enrichment);
    return removed;
  };
  const auto compact_domain =
      [&compact_blocks](fem::singular::ParallelSparseEnrichmentMatrices &matrices)
  {
    return compact_blocks(matrices.h1_diffusion) + compact_blocks(matrices.h1_mass) +
           compact_blocks(matrices.nd_mass) + compact_blocks(matrices.nd_curl_curl);
  };
  singular_domain_matrices[finest_level] = assemble_domain(0);
  singular_domain_abs_matrices[finest_level] = assemble_domain(1);
  HYPRE_BigInt source_zeros_removed =
      compact_domain(singular_domain_matrices[finest_level]) +
      compact_domain(singular_domain_abs_matrices[finest_level]);
  if (mat_op.HasLossTangent())
  {
    singular_domain_imag_matrices[finest_level] = assemble_domain(2);
    source_zeros_removed += compact_domain(singular_domain_imag_matrices[finest_level]);
  }
  singular_nd_element_patch_matrices.clear();
  Mpi::Print(" Singular source compaction removed {:d} exact zeros from finest-level "
             "domain blocks\n",
             source_zeros_removed);
  report_stage("finest-level parallel domain assembly");

  fem::singular::LocalSparseOperatorBlocks local_lumped_stiffness, local_lumped_damping,
      local_lumped_mass, local_impedance_stiffness, local_impedance_damping,
      local_impedance_mass;
  const auto assemble_boundary =
      [&](const std::map<int, double> &coefficients,
          fem::singular::LocalSparseOperatorBlocks &local_boundary)
  {
    fem::singular::ParallelSparseOperatorBlocks result;
    if (coefficients.empty())
    {
      return result;
    }
    local_boundary =
        tetrahedral
            ? fem::singular::AssembleLocalSparseNDBoundaryMassMatrices(
                  *singular_dofs, finest_nd_space.Get(), coefficients, options)
            : fem::singular::AssembleLocalSparseNDBoundaryMassMatrices(
                  *triangle_singular_dofs, finest_nd_space.Get(), coefficients, options);
    return fem::singular::AssembleParallelSparseNDBoundaryMassMatrices(
        local_boundary, *singular_numbering, finest_nd_space.Get());
  };
  singular_lumped_stiffness_matrices[finest_level] = assemble_boundary(
      lumped_port_op.GetStiffnessBdrCoefficientMap(), local_lumped_stiffness);
  singular_lumped_damping_matrices[finest_level] =
      assemble_boundary(lumped_port_op.GetDampingBdrCoefficientMap(), local_lumped_damping);
  singular_lumped_mass_matrices[finest_level] =
      assemble_boundary(lumped_port_op.GetMassBdrCoefficientMap(), local_lumped_mass);
  singular_impedance_stiffness_matrices[finest_level] =
      assemble_boundary(impedance_stiffness_coefficients, local_impedance_stiffness);
  singular_impedance_damping_matrices[finest_level] =
      assemble_boundary(impedance_damping_coefficients, local_impedance_damping);
  singular_impedance_mass_matrices[finest_level] =
      assemble_boundary(impedance_mass_coefficients, local_impedance_mass);
  source_zeros_removed =
      compact_blocks(singular_lumped_stiffness_matrices[finest_level]) +
      compact_blocks(singular_lumped_damping_matrices[finest_level]) +
      compact_blocks(singular_lumped_mass_matrices[finest_level]) +
      compact_blocks(singular_impedance_stiffness_matrices[finest_level]) +
      compact_blocks(singular_impedance_damping_matrices[finest_level]) +
      compact_blocks(singular_impedance_mass_matrices[finest_level]);
  Mpi::Print(" Singular source compaction removed {:d} exact zeros from finest-level "
             "boundary blocks\n",
             source_zeros_removed);
  report_stage("finest-level singular boundary assembly");

  for (std::size_t fine_level = finest_level; fine_level > 0; fine_level--)
  {
    const std::size_t coarse_level = fine_level - 1;
    const auto *h1_prolongation = dynamic_cast<const ParOperator *>(
        &GetH1Spaces().GetProlongationAtLevel(coarse_level));
    const auto *nd_prolongation = dynamic_cast<const ParOperator *>(
        &GetNDSpaces().GetProlongationAtLevel(coarse_level));
    MFEM_VERIFY(h1_prolongation && nd_prolongation,
                "Full-wave singular p-multigrid requires assembled standard H1 and ND "
                "prolongation operators!");
    const auto &h1_parallel = h1_prolongation->ParallelAssemble();
    const auto &nd_parallel = nd_prolongation->ParallelAssemble();
    singular_domain_matrices[coarse_level] =
        fem::singular::RestrictParallelSparseEnrichmentMatrices(
            singular_domain_matrices[fine_level], h1_parallel, nd_parallel);
    singular_domain_abs_matrices[coarse_level] =
        fem::singular::RestrictParallelSparseEnrichmentMatrices(
            singular_domain_abs_matrices[fine_level], h1_parallel, nd_parallel);
    if (mat_op.HasLossTangent())
    {
      singular_domain_imag_matrices[coarse_level] =
          fem::singular::RestrictParallelSparseEnrichmentMatrices(
              singular_domain_imag_matrices[fine_level], h1_parallel, nd_parallel);
    }
    const auto restrict_boundary =
        [&nd_parallel](const fem::singular::ParallelSparseOperatorBlocks &fine)
    { return fem::singular::RestrictParallelSparseOperatorBlocks(fine, nd_parallel); };
    singular_lumped_stiffness_matrices[coarse_level] =
        restrict_boundary(singular_lumped_stiffness_matrices[fine_level]);
    singular_lumped_damping_matrices[coarse_level] =
        restrict_boundary(singular_lumped_damping_matrices[fine_level]);
    singular_lumped_mass_matrices[coarse_level] =
        restrict_boundary(singular_lumped_mass_matrices[fine_level]);
    singular_impedance_stiffness_matrices[coarse_level] =
        restrict_boundary(singular_impedance_stiffness_matrices[fine_level]);
    singular_impedance_damping_matrices[coarse_level] =
        restrict_boundary(singular_impedance_damping_matrices[fine_level]);
    singular_impedance_mass_matrices[coarse_level] =
        restrict_boundary(singular_impedance_mass_matrices[fine_level]);
    compact_domain(singular_domain_matrices[coarse_level]);
    compact_domain(singular_domain_abs_matrices[coarse_level]);
    if (mat_op.HasLossTangent())
    {
      compact_domain(singular_domain_imag_matrices[coarse_level]);
    }
    compact_blocks(singular_lumped_stiffness_matrices[coarse_level]);
    compact_blocks(singular_lumped_damping_matrices[coarse_level]);
    compact_blocks(singular_lumped_mass_matrices[coarse_level]);
    compact_blocks(singular_impedance_stiffness_matrices[coarse_level]);
    compact_blocks(singular_impedance_damping_matrices[coarse_level]);
    compact_blocks(singular_impedance_mass_matrices[coarse_level]);
  }
  const auto discard_error_bounds = [](fem::singular::ParallelSparseOperatorBlocks &blocks)
  {
    blocks.enrichment_enrichment_estimated_absolute_error.reset();
    blocks.standard_enrichment_estimated_absolute_error.reset();
  };
  const auto discard_domain_error_bounds =
      [&discard_error_bounds](fem::singular::ParallelSparseEnrichmentMatrices &matrices)
  {
    discard_error_bounds(matrices.h1_diffusion);
    discard_error_bounds(matrices.h1_mass);
    discard_error_bounds(matrices.nd_mass);
    discard_error_bounds(matrices.nd_curl_curl);
  };
  for (std::size_t level = 0; level < number_levels; level++)
  {
    discard_domain_error_bounds(singular_domain_matrices[level]);
    discard_domain_error_bounds(singular_domain_abs_matrices[level]);
    if (mat_op.HasLossTangent())
    {
      discard_domain_error_bounds(singular_domain_imag_matrices[level]);
    }
    discard_error_bounds(singular_lumped_stiffness_matrices[level]);
    discard_error_bounds(singular_lumped_damping_matrices[level]);
    discard_error_bounds(singular_lumped_mass_matrices[level]);
    discard_error_bounds(singular_impedance_stiffness_matrices[level]);
    discard_error_bounds(singular_impedance_damping_matrices[level]);
    discard_error_bounds(singular_impedance_mass_matrices[level]);
  }
  report_stage("coarse-level sparse restriction");

  for (std::size_t level = 0; level < number_levels; level++)
  {
    auto &h1_space = GetH1Spaces().GetFESpaceAtLevel(level);
    auto &nd_space = GetNDSpaces().GetFESpaceAtLevel(level);
    auto enrichment_gradient =
        fem::singular::BuildParallelEnrichmentGradient(GetComm(), *singular_numbering);
    const auto &standard_gradient_operator = nd_space.GetDiscreteInterpolator(h1_space);
    const auto *standard_gradient =
        dynamic_cast<const ParOperator *>(&standard_gradient_operator);
    MFEM_VERIFY(standard_gradient,
                "Full-wave singular enrichment requires an assembled standard gradient!");
    singular_gradients.push_back(fem::singular::BuildParallelEnrichedGradient(
        standard_gradient->ParallelAssemble(), *enrichment_gradient));
  }
  report_stage("enriched gradient assembly");

  mfem::Array<int> singular_essential_attributes = dbc_attr;
  singular_essential_attributes.Append(
      boundaries.auxpec.attributes.data(),
      static_cast<int>(boundaries.auxpec.attributes.size()));
  for (int attribute : constrained_impedance_attributes)
  {
    singular_essential_attributes.Append(attribute);
  }
  singular_essential_attributes.Sort();
  singular_essential_attributes.Unique();
  if (tetrahedral)
  {
    singular_h1_essential_true_dofs = fem::singular::GetEssentialH1TrueDofs(
        GetComm(), *singular_features, *singular_dofs, *singular_numbering,
        singular_essential_attributes);
    singular_nd_essential_true_dofs = fem::singular::GetEssentialNDTrueDofs(
        GetComm(), *singular_features, *singular_dofs, *singular_numbering,
        singular_essential_attributes);
  }
  else
  {
    singular_h1_essential_true_dofs = fem::singular::GetEssentialTriangleH1TrueDofs(
        GetComm(), *triangle_singular_features, *triangle_singular_dofs,
        *singular_numbering, singular_essential_attributes);
    singular_nd_essential_true_dofs = fem::singular::GetEssentialTriangleNDTrueDofs(
        GetComm(), *triangle_singular_features, *triangle_singular_dofs,
        *singular_numbering, singular_essential_attributes);
  }

  mfem::Array<int> singular_h1_auxiliary_true_dofs;
  const std::set<int> auxiliary_boundary_attributes(aux_bdr_attr.begin(),
                                                    aux_bdr_attr.end());
  if (tetrahedral)
  {
    const auto auxiliary_boundary_faces =
        GatherBoundaryEntities<3>(GetMesh().Get(), *source_vertex_ids,
                                  auxiliary_boundary_attributes, mfem::Geometry::TRIANGLE);
    singular_h1_auxiliary_true_dofs = fem::singular::GetEssentialH1TrueDofsOnFaces(
        GetComm(), auxiliary_boundary_faces, *singular_dofs, *singular_numbering);
  }
  else
  {
    const auto auxiliary_boundary_segments =
        GatherBoundaryEntities<2>(GetMesh().Get(), *source_vertex_ids,
                                  auxiliary_boundary_attributes, mfem::Geometry::SEGMENT);
    singular_h1_auxiliary_true_dofs =
        fem::singular::GetEssentialTriangleH1TrueDofsOnSegments(
            GetComm(), auxiliary_boundary_segments, *triangle_singular_dofs,
            *singular_numbering);
  }

  combined_nd_dbc_tdof_lists = nd_dbc_tdof_lists;
  combined_h1_dbc_tdof_lists = h1_dbc_tdof_lists;
  combined_h1_aux_bdr_tdof_lists = aux_bdr_tdof_lists;
  for (std::size_t level = 0; level < number_levels; level++)
  {
    auto &combined_nd = combined_nd_dbc_tdof_lists[level];
    const int standard_nd_size = GetNDSpaces().GetFESpaceAtLevel(level).GetTrueVSize();
    for (int dof : singular_nd_essential_true_dofs)
    {
      combined_nd.Append(standard_nd_size + dof);
    }
    combined_nd.Sort();

    auto &combined_h1 = combined_h1_dbc_tdof_lists[level];
    const int standard_h1_size = GetH1Spaces().GetFESpaceAtLevel(level).GetTrueVSize();
    for (int dof : singular_h1_essential_true_dofs)
    {
      combined_h1.Append(standard_h1_size + dof);
    }
    combined_h1.Sort();

    auto &combined_h1_auxiliary = combined_h1_aux_bdr_tdof_lists[level];
    for (int dof : singular_h1_auxiliary_true_dofs)
    {
      combined_h1_auxiliary.Append(standard_h1_size + dof);
    }
    combined_h1_auxiliary.Sort();
    MFEM_VERIFY(
        std::adjacent_find(combined_nd.begin(), combined_nd.end()) == combined_nd.end() &&
            std::adjacent_find(combined_h1.begin(), combined_h1.end()) ==
                combined_h1.end() &&
            std::adjacent_find(combined_h1_auxiliary.begin(),
                               combined_h1_auxiliary.end()) == combined_h1_auxiliary.end(),
        "Full-wave singular essential true DOFs are not unique!");
  }
  report_stage("essential DOF construction");

  singular_nd_prolongations.reserve(number_levels > 0 ? number_levels - 1 : 0);
  singular_h1_prolongations.reserve(number_levels > 0 ? number_levels - 1 : 0);
  singular_nd_coordinate_shifts.clear();
  for (std::size_t level = 0; level + 1 < number_levels; level++)
  {
    const auto &standard_nd_operator = GetNDSpaces().GetProlongationAtLevel(level);
    const auto &standard_h1_operator = GetH1Spaces().GetProlongationAtLevel(level);
    const auto *standard_nd_prolongation =
        dynamic_cast<const ParOperator *>(&standard_nd_operator);
    const auto *standard_h1_prolongation =
        dynamic_cast<const ParOperator *>(&standard_h1_operator);
    MFEM_VERIFY(standard_nd_prolongation && standard_h1_prolongation,
                "Full-wave singular p-multigrid requires assembled standard ND and H1 "
                "prolongation operators!");
    singular_h1_prolongations.push_back(fem::singular::BuildParallelEnrichedProlongation(
        standard_h1_prolongation->ParallelAssemble(), singular_numbering->h1));
    singular_nd_prolongations.push_back(fem::singular::BuildParallelEnrichedProlongation(
        standard_nd_prolongation->ParallelAssemble(), singular_numbering->nd));
  }
  report_stage("enriched prolongation assembly");

  double total = Timer::Duration(Timer::Now() - setup_start).count();
  double total_minimum = total;
  double total_maximum = total;
  Mpi::GlobalMin(1, &total_minimum, GetComm());
  Mpi::GlobalMax(1, &total_maximum, GetComm());
  Mpi::GlobalSum(1, &total, GetComm());
  Mpi::Print(" Singular setup timing, total (s): min. {:.3f}, max. {:.3f}, avg. {:.3f}\n",
             total_minimum, total_maximum, total / Mpi::Size(GetComm()));

  Mpi::Print(" Singular full-wave enrichment: {:d} ND + {:d} H1 global true DOFs on "
             "{:d} polynomial level{}\n",
             singular_numbering->nd.global_size, singular_numbering->h1.global_size,
             number_levels, number_levels == 1 ? "" : "s");
}

const Operator &SpaceOperator::GetGradMatrix() const
{
  return !singular_gradients.empty()
             ? static_cast<const Operator &>(*singular_gradients.back())
             : GetNDSpace().GetDiscreteInterpolator(GetH1Space());
}

std::vector<const Operator *> SpaceOperator::GetCombinedNDProlongationOperators() const
{
  std::vector<const Operator *> operators;
  operators.reserve(singular_nd_prolongations.size());
  for (const auto &prolongation : singular_nd_prolongations)
  {
    operators.push_back(prolongation.get());
  }
  return operators;
}

std::vector<const Operator *> SpaceOperator::GetCombinedH1ProlongationOperators() const
{
  std::vector<const Operator *> operators;
  operators.reserve(singular_h1_prolongations.size());
  for (const auto &prolongation : singular_h1_prolongations)
  {
    operators.push_back(prolongation.get());
  }
  return operators;
}

const mfem::HypreParMatrix *SpaceOperator::GetFinestSingularNDCoordinateShift() const
{
  if (singular_nd_coordinate_shifts.empty())
  {
    return nullptr;
  }
  return singular_nd_coordinate_shifts.back().get();
}

std::vector<const Operator *> SpaceOperator::GetCombinedGradientOperators() const
{
  std::vector<const Operator *> operators;
  operators.reserve(singular_gradients.size());
  for (const auto &gradient : singular_gradients)
  {
    operators.push_back(gradient.get());
  }
  return operators;
}

void SpaceOperator::CheckSingularExcitations(ProblemType problem_type) const
{
  if (problem_type != ProblemType::DRIVEN)
  {
    return;
  }
  MFEM_VERIFY(current_dipole_op.Empty(),
              "Driven singular elements do not yet assemble current-dipole enrichment "
              "load vectors!");

  const auto source_attributes = surf_j_op.GetAttrList();
  if (source_attributes.Size() == 0)
  {
    return;
  }
  const int maximum_attribute =
      GetMesh().Get().bdr_attributes.Size() ? GetMesh().Get().bdr_attributes.Max() : 0;
  const auto source_marker = mesh::AttrToMarker(maximum_attribute, source_attributes);
  bool intersects_enrichment = false;
  for (int boundary_element = 0; boundary_element < GetMesh().GetNBE(); boundary_element++)
  {
    const int attribute = GetMesh().Get().GetBdrAttribute(boundary_element);
    if (attribute <= 0 || attribute > source_marker.Size() || !source_marker[attribute - 1])
    {
      continue;
    }
    int element = -1, face = -1;
    GetMesh().Get().GetBdrElementAdjacentElement(boundary_element, element, face);
    if (element < 0)
    {
      continue;
    }
    intersects_enrichment =
        intersects_enrichment ||
        (singular_features && (!singular_features->elements[element].nodes.empty() ||
                               !singular_features->elements[element].edges.empty())) ||
        (triangle_singular_features &&
         !triangle_singular_features->elements[element].nodes.empty());
  }
  Mpi::GlobalOr(1, &intersects_enrichment, GetComm());
  MFEM_VERIFY(!intersects_enrichment,
              "A driven surface-current source touches an enriched element. Its singular "
              "load-vector entries must be implemented before this source is supported!");
}

void SpaceOperator::CheckExcitations(ProblemType problem_type) const
{
  if (problem_type == ProblemType::DRIVEN)
  {
    MFEM_VERIFY(!port_excitation_helper.Empty(),
                "Driven problems must specify at least one excitation!");
  }
  else if (problem_type == ProblemType::EIGENMODE)
  {
    MFEM_VERIFY(port_excitation_helper.Empty(),
                "Eigenmode problems must not specify any excitation!");
  }
  else if (problem_type == ProblemType::TRANSIENT)
  {
    MFEM_VERIFY(
        port_excitation_helper.Size() == 1,
        "Transient problems currently only support a single excitation per simulation!");
  }
}

mfem::Array<int> SpaceOperator::SetUpBoundaryProperties(const config::PecBoundaryData &pec,
                                                        const mfem::ParMesh &mesh)
{
  // Check that boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!pec.empty())
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
      // MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
      //             "PEC boundary attribute tags must be non-negative and correspond to "
      //             "attributes in the mesh!");
      // MFEM_VERIFY(bdr_attr_marker[attr - 1],
      //             "Unknown PEC boundary attribute " << attr << "!");
      if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
      {
        bdr_warn_list.insert(attr);
      }
    }
    if (!bdr_warn_list.empty())
    {
      Mpi::Print("\n");
      Mpi::Warning("Unknown PEC boundary attributes!\nSolver will just ignore them!");
      utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
      Mpi::Print("\n");
    }
  }

  // Mark selected boundary attributes from the mesh as essential (Dirichlet).
  mfem::Array<int> dbc_bcs;
  dbc_bcs.Reserve(static_cast<int>(pec.attributes.size()));
  for (auto attr : pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
    {
      continue;  // Can just ignore if wrong
    }
    dbc_bcs.Append(attr);
  }
  return dbc_bcs;
}

void SpaceOperator::CheckBoundaryProperties()
{
  // Mark selected boundary attributes from the mesh as having some Dirichlet, Neumann, or
  // mixed BC applied.
  const mfem::ParMesh &mesh = GetMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, dbc_attr);
  const auto farfield_marker = mesh::AttrToMarker(bdr_attr_max, farfield_op.GetAttrList());
  const auto surf_sigma_marker =
      mesh::AttrToMarker(bdr_attr_max, surf_sigma_op.GetAttrList());
  const auto surf_z_Rs_marker = mesh::AttrToMarker(bdr_attr_max, surf_z_op.GetRsAttrList());
  const auto surf_z_Ls_marker = mesh::AttrToMarker(bdr_attr_max, surf_z_op.GetLsAttrList());
  const auto surf_rz_marker = mesh::AttrToMarker(bdr_attr_max, surf_rz_op.GetAttrList());
  const auto lumped_port_Rs_marker =
      mesh::AttrToMarker(bdr_attr_max, lumped_port_op.GetRsAttrList());
  const auto lumped_port_Ls_marker =
      mesh::AttrToMarker(bdr_attr_max, lumped_port_op.GetLsAttrList());
  const auto wave_port_marker =
      mesh::AttrToMarker(bdr_attr_max, wave_port_op.GetAttrList());
  const auto floquet_port_marker =
      mesh::AttrToMarker(bdr_attr_max, floquet_port_op.GetAttrList());
  mfem::Array<int> aux_bdr_marker(dbc_marker.Size());
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    aux_bdr_marker[i] = (dbc_marker[i] || farfield_marker[i] || surf_sigma_marker[i] ||
                         surf_z_Rs_marker[i] || surf_rz_marker[i] || surf_z_Ls_marker[i] ||
                         lumped_port_Rs_marker[i] || lumped_port_Ls_marker[i] ||
                         wave_port_marker[i] || floquet_port_marker[i]);
    if (aux_bdr_marker[i])
    {
      aux_bdr_attr.Append(i + 1);
    }
  }
  // aux_bdr_marker = 1;  // Mark all boundaries (including material interfaces
  //                      // added during mesh preprocessing)
  //                      // As tested, this does not eliminate all DC modes!
  for (std::size_t l = 0; l < GetH1Spaces().GetNumLevels(); l++)
  {
    GetH1Spaces().GetFESpaceAtLevel(l).Get().GetEssentialTrueDofs(
        aux_bdr_marker, aux_bdr_tdof_lists.emplace_back());
  }

  // A final check that no boundary attribute is assigned multiple boundary conditions.
  const auto surf_z_marker = mesh::AttrToMarker(bdr_attr_max, surf_z_op.GetAttrList());
  const auto lumped_port_marker =
      mesh::AttrToMarker(bdr_attr_max, lumped_port_op.GetAttrList());
  const auto surf_j_marker = mesh::AttrToMarker(bdr_attr_max, surf_j_op.GetAttrList());
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    MFEM_VERIFY(dbc_marker[i] + farfield_marker[i] + surf_sigma_marker[i] +
                        surf_z_marker[i] + surf_rz_marker[i] + lumped_port_marker[i] +
                        wave_port_marker[i] + floquet_port_marker[i] + surf_j_marker[i] <=
                    1,
                "Boundary attributes should not be specified with multiple BC!");
  }
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
               (nd_fespace.GetMaxElementOrder() >= BilinearForm::pa_order_threshold)
                   ? "Partial"
                   : "Full");

    const auto &mesh = *nd_fespace.GetParMesh();
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
  }
  print_hdr = false;
}

void AddIntegrators(BilinearForm &a, const MaterialPropertyCoefficient *df,
                    const MaterialPropertyCoefficient *f,
                    const MaterialPropertyCoefficient *dfb,
                    const MaterialPropertyCoefficient *fb,
                    const MaterialPropertyCoefficient *fp, bool assemble_q_data = false)
{
  if (df && !df->empty() && f && !f->empty())
  {
    a.AddDomainIntegrator<CurlCurlMassIntegrator>(*df, *f);
  }
  else
  {
    if (df && !df->empty())
    {
      a.AddDomainIntegrator<CurlCurlIntegrator>(*df);
    }
    if (f && !f->empty())
    {
      a.AddDomainIntegrator<VectorFEMassIntegrator>(*f);
    }
  }
  if (dfb && !dfb->empty() && fb && !fb->empty())
  {
    a.AddBoundaryIntegrator<CurlCurlMassIntegrator>(*dfb, *fb);
  }
  else
  {
    if (dfb && !dfb->empty())
    {
      a.AddBoundaryIntegrator<CurlCurlIntegrator>(*dfb);
    }
    if (fb && !fb->empty())
    {
      a.AddBoundaryIntegrator<VectorFEMassIntegrator>(*fb);
    }
  }
  if (fp && !fp->empty())
  {
    a.AddDomainIntegrator<MixedVectorWeakCurlIntegrator>(*fp);
    a.AddDomainIntegrator<MixedVectorCurlIntegrator>(*fp, true);
  }
  if (assemble_q_data)
  {
    a.AssembleQuadratureData();
  }
}

void AddAuxIntegrators(BilinearForm &a, const MaterialPropertyCoefficient *f,
                       const MaterialPropertyCoefficient *fb, bool assemble_q_data = false)
{
  if (f && !f->empty())
  {
    a.AddDomainIntegrator<DiffusionIntegrator>(*f);
  }
  if (fb && !fb->empty())
  {
    a.AddBoundaryIntegrator<DiffusionIntegrator>(*fb);
  }
  if (assemble_q_data)
  {
    a.AssembleQuadratureData();
  }
}

auto AssembleOperator(const FiniteElementSpace &fespace,
                      const MaterialPropertyCoefficient *df,
                      const MaterialPropertyCoefficient *f,
                      const MaterialPropertyCoefficient *dfb,
                      const MaterialPropertyCoefficient *fb,
                      const MaterialPropertyCoefficient *fp, bool skip_zeros = false,
                      bool assemble_q_data = false)
{
  BilinearForm a(fespace);
  AddIntegrators(a, df, f, dfb, fb, fp, assemble_q_data);
  return a.Assemble(skip_zeros);
}

auto AssembleOperators(const FiniteElementSpaceHierarchy &fespaces,
                       const MaterialPropertyCoefficient *df,
                       const MaterialPropertyCoefficient *f,
                       const MaterialPropertyCoefficient *dfb,
                       const MaterialPropertyCoefficient *fb,
                       const MaterialPropertyCoefficient *fp, bool skip_zeros = false,
                       bool assemble_q_data = false, std::size_t l0 = 0)
{
  BilinearForm a(fespaces.GetFinestFESpace());
  AddIntegrators(a, df, f, dfb, fb, fp, assemble_q_data);
  return a.Assemble(fespaces, skip_zeros, l0);
}

auto AssembleAuxOperators(const FiniteElementSpaceHierarchy &fespaces,
                          const MaterialPropertyCoefficient *f,
                          const MaterialPropertyCoefficient *fb, bool skip_zeros = false,
                          bool assemble_q_data = false, std::size_t l0 = 0)
{
  BilinearForm a(fespaces.GetFinestFESpace());
  AddAuxIntegrators(a, f, fb, assemble_q_data);
  return a.Assemble(fespaces, skip_zeros, l0);
}

fem::singular::ParallelSparseOperatorBlocks
AddSingularOperatorBlocks(const fem::singular::ParallelSparseOperatorBlocks &domain,
                          const fem::singular::ParallelSparseOperatorBlocks &boundary)
{
  MFEM_VERIFY(domain.standard_enrichment && domain.enrichment_standard &&
                  domain.enrichment_enrichment,
              "Cannot combine incomplete singular domain operator blocks!");
  const bool has_boundary = boundary.standard_enrichment != nullptr;
  MFEM_VERIFY(!has_boundary ||
                  (boundary.enrichment_standard && boundary.enrichment_enrichment),
              "Cannot combine incomplete singular boundary operator blocks!");
  const auto add = [has_boundary](const mfem::HypreParMatrix &a,
                                  const std::unique_ptr<mfem::HypreParMatrix> &b)
  {
    return has_boundary ? std::unique_ptr<mfem::HypreParMatrix>(mfem::Add(1.0, a, 1.0, *b))
                        : std::make_unique<mfem::HypreParMatrix>(a);
  };

  fem::singular::ParallelSparseOperatorBlocks result;
  result.standard_enrichment =
      add(*domain.standard_enrichment, boundary.standard_enrichment);
  result.enrichment_standard =
      add(*domain.enrichment_standard, boundary.enrichment_standard);
  result.enrichment_enrichment =
      add(*domain.enrichment_enrichment, boundary.enrichment_enrichment);
  if (domain.transformed_enrichment_diagonal &&
      (!has_boundary || boundary.transformed_enrichment_diagonal))
  {
    result.transformed_enrichment_diagonal =
        std::make_unique<Vector>(*domain.transformed_enrichment_diagonal);
    if (has_boundary)
    {
      *result.transformed_enrichment_diagonal += *boundary.transformed_enrichment_diagonal;
    }
  }
  return result;
}

fem::singular::ParallelSparseOperatorBlocks AddScaledSingularOperatorBlocks(
    std::initializer_list<
        std::pair<double, const fem::singular::ParallelSparseOperatorBlocks *>>
        terms)
{
  fem::singular::ParallelSparseOperatorBlocks result;
  const fem::singular::ParallelSparseOperatorBlocks *zero_template = nullptr;
  bool transformed_diagonal_complete = true;
  const auto add = [](std::unique_ptr<mfem::HypreParMatrix> &sum, double coefficient,
                      const mfem::HypreParMatrix &matrix)
  {
    if (!sum)
    {
      sum = std::make_unique<mfem::HypreParMatrix>(matrix);
      *sum *= coefficient;
    }
    else
    {
      sum.reset(mfem::Add(1.0, *sum, coefficient, matrix));
    }
  };
  for (const auto &[coefficient, blocks] : terms)
  {
    if (!blocks || !blocks->standard_enrichment)
    {
      continue;
    }
    MFEM_VERIFY(blocks->enrichment_standard && blocks->enrichment_enrichment,
                "Cannot combine incomplete singular operator blocks!");
    zero_template = zero_template ? zero_template : blocks;
    if (coefficient == 0.0)
    {
      continue;
    }
    add(result.standard_enrichment, coefficient, *blocks->standard_enrichment);
    add(result.enrichment_standard, coefficient, *blocks->enrichment_standard);
    add(result.enrichment_enrichment, coefficient, *blocks->enrichment_enrichment);
    if (blocks->transformed_enrichment_diagonal)
    {
      if (!result.transformed_enrichment_diagonal)
      {
        result.transformed_enrichment_diagonal =
            std::make_unique<Vector>(*blocks->transformed_enrichment_diagonal);
        *result.transformed_enrichment_diagonal *= coefficient;
      }
      else
      {
        result.transformed_enrichment_diagonal->Add(
            coefficient, *blocks->transformed_enrichment_diagonal);
      }
    }
    else
    {
      transformed_diagonal_complete = false;
    }
  }
  if (!result.standard_enrichment && zero_template)
  {
    add(result.standard_enrichment, 0.0, *zero_template->standard_enrichment);
    add(result.enrichment_standard, 0.0, *zero_template->enrichment_standard);
    add(result.enrichment_enrichment, 0.0, *zero_template->enrichment_enrichment);
  }
  MFEM_VERIFY(result.standard_enrichment && result.enrichment_standard &&
                  result.enrichment_enrichment,
              "Singular preconditioner has no compatible enrichment block structure!");
  if (!transformed_diagonal_complete)
  {
    result.transformed_enrichment_diagonal.reset();
  }
  return result;
}

}  // namespace

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetStiffnessMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient df(mat_op.MaxCeedAttribute()), f(mat_op.MaxCeedAttribute()),
      fb(mat_op.MaxCeedBdrAttribute()), fc(mat_op.MaxCeedAttribute());
  AddStiffnessCoefficients(1.0, df, f);
  AddStiffnessBdrCoefficients(1.0, fb);
  if (!mat_op.HasFloquetFrequencyScaling())
  {
    AddRealPeriodicCoefficients(1.0, f);
    AddImagPeriodicCoefficients(1.0, fc);
  }
  int empty[2] = {(df.empty() && f.empty() && fb.empty()), (fc.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> kr, ki;
  if (!empty[0])
  {
    kr = AssembleOperator(GetNDSpace(), &df, &f, nullptr, &fb, nullptr, skip_zeros);
  }
  if (!empty[1])
  {
    ki =
        AssembleOperator(GetNDSpace(), nullptr, nullptr, nullptr, nullptr, &fc, skip_zeros);
  }
  if (HasSingularEnrichment())
  {
    MFEM_VERIFY(kr && !ki, "Full-wave singular stiffness assembly requires a real domain "
                           "curl-curl operator without Floquet terms!");
    auto standard = std::make_unique<ParOperator>(std::move(kr), GetNDSpace());
    standard->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    auto enrichment = AddScaledSingularOperatorBlocks(
        {{1.0, &singular_domain_matrices.back().nd_curl_curl},
         {1.0, &singular_lumped_stiffness_matrices.back()},
         {1.0, &singular_impedance_stiffness_matrices.back()}});
    auto combined = std::make_unique<fem::singular::ParallelHybridEnrichedOperator>(
        std::move(standard), std::move(enrichment), nd_dbc_tdof_lists.back(),
        singular_nd_essential_true_dofs, diag_policy);
    if constexpr (std::is_same<OperType, ComplexOperator>::value)
    {
      return std::make_unique<ComplexWrapperOperator>(std::move(combined), nullptr);
    }
    else
    {
      return combined;
    }
  }
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto K =
        std::make_unique<ComplexParOperator>(std::move(kr), std::move(ki), GetNDSpace());
    K->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return K;
  }
  else
  {
    MFEM_VERIFY(!ki, "Unexpected imaginary part in GetStiffnessMatrix<Operator>!");
    auto K = std::make_unique<ParOperator>(std::move(kr), GetNDSpace());
    K->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return K;
  }
}

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetDampingMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient f(mat_op.MaxCeedAttribute()),
      fb(mat_op.MaxCeedBdrAttribute());
  MaterialPropertyCoefficient fp(mat_op.MaxCeedAttribute());
  AddDampingCoefficients(1.0, f);
  AddDampingBdrCoefficients(1.0, fb);
  if (mat_op.HasFloquetFrequencyScaling())
  {
    AddImagPeriodicCoefficients(1.0, fp);
  }
  int empty = (f.empty() && fb.empty() && fp.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto c = AssembleOperator(GetNDSpace(), nullptr, &f, nullptr, &fb, &fp, skip_zeros);
  if (HasSingularEnrichment())
  {
    MFEM_VERIFY(c && !mat_op.HasConductivity() && !mat_op.HasFloquetFrequencyScaling(),
                "Full-wave singular damping assembly supports only lumped-port and "
                "surface-impedance resistance!");
    auto standard = std::make_unique<ParOperator>(std::move(c), GetNDSpace());
    standard->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    auto enrichment = AddScaledSingularOperatorBlocks(
        {{1.0, &singular_lumped_damping_matrices.back()},
         {1.0, &singular_impedance_damping_matrices.back()}});
    auto combined = std::make_unique<fem::singular::ParallelHybridEnrichedOperator>(
        std::move(standard), std::move(enrichment), nd_dbc_tdof_lists.back(),
        singular_nd_essential_true_dofs, diag_policy);
    if constexpr (std::is_same<OperType, ComplexOperator>::value)
    {
      return std::make_unique<ComplexWrapperOperator>(std::move(combined), nullptr);
    }
    else
    {
      return combined;
    }
  }
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto C = std::make_unique<ComplexParOperator>(std::move(c), nullptr, GetNDSpace());
    C->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return C;
  }
  else
  {
    auto C = std::make_unique<ParOperator>(std::move(c), GetNDSpace());
    C->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return C;
  }
}

template <typename OperType>
std::unique_ptr<OperType> SpaceOperator::GetMassMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient fr(mat_op.MaxCeedAttribute()), fi(mat_op.MaxCeedAttribute()),
      fbr(mat_op.MaxCeedBdrAttribute()), fbi(mat_op.MaxCeedBdrAttribute());
  AddRealMassCoefficients(1.0, fr);
  AddRealMassBdrCoefficients(1.0, fbr);
  if (mat_op.HasFloquetFrequencyScaling())
  {
    AddRealPeriodicCoefficients(-1.0, fr);
  }
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    AddImagMassCoefficients(1.0, fi);
  }
  int empty[2] = {(fr.empty() && fbr.empty()), (fi.empty() && fbi.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> mr, mi;
  if (!empty[0])
  {
    mr = AssembleOperator(GetNDSpace(), nullptr, &fr, nullptr, &fbr, nullptr, skip_zeros);
  }
  if (!empty[1])
  {
    mi = AssembleOperator(GetNDSpace(), nullptr, &fi, nullptr, &fbi, nullptr, skip_zeros);
  }
  if (HasSingularEnrichment())
  {
    MFEM_VERIFY(mr, "Full-wave singular mass assembly requires a positive real "
                    "domain permittivity operator!");
    auto standard_real = std::make_unique<ParOperator>(std::move(mr), GetNDSpace());
    standard_real->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    auto real_enrichment =
        AddScaledSingularOperatorBlocks({{1.0, &singular_domain_matrices.back().nd_mass},
                                         {1.0, &singular_lumped_mass_matrices.back()},
                                         {1.0, &singular_impedance_mass_matrices.back()}});
    auto combined_real = std::make_unique<fem::singular::ParallelHybridEnrichedOperator>(
        std::move(standard_real), std::move(real_enrichment), nd_dbc_tdof_lists.back(),
        singular_nd_essential_true_dofs, diag_policy);
    if constexpr (std::is_same<OperType, ComplexOperator>::value)
    {
      std::unique_ptr<Operator> combined_imaginary;
      if (mi)
      {
        MFEM_VERIFY(singular_domain_imag_matrices.back().nd_mass.enrichment_enrichment,
                    "Full-wave singular imaginary mass blocks were not assembled!");
        auto standard_imaginary =
            std::make_unique<ParOperator>(std::move(mi), GetNDSpace());
        standard_imaginary->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(),
                                                 Operator::DiagonalPolicy::DIAG_ZERO);
        combined_imaginary =
            std::make_unique<fem::singular::ParallelHybridEnrichedOperator>(
                std::move(standard_imaginary), singular_domain_imag_matrices.back().nd_mass,
                nd_dbc_tdof_lists.back(), singular_nd_essential_true_dofs,
                Operator::DiagonalPolicy::DIAG_ZERO);
      }
      return std::make_unique<ComplexWrapperOperator>(std::move(combined_real),
                                                      std::move(combined_imaginary));
    }
    else
    {
      return combined_real;
    }
  }
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto M =
        std::make_unique<ComplexParOperator>(std::move(mr), std::move(mi), GetNDSpace());
    M->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M;
  }
  else
  {
    auto M = std::make_unique<ParOperator>(std::move(mr), GetNDSpace());
    M->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M;
  }
}

std::unique_ptr<Operator>
SpaceOperator::GetBulkStiffnessMatrix(Operator::DiagonalPolicy diag_policy)
{
  MFEM_VERIFY(HasSingularEnrichment(),
              "The combined bulk stiffness matrix requires singular enrichment!");
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient df(mat_op.MaxCeedAttribute()), f(mat_op.MaxCeedAttribute());
  AddStiffnessCoefficients(1.0, df, f);
  constexpr bool skip_zeros = false;
  auto kr = AssembleOperator(GetNDSpace(), &df, &f, nullptr, nullptr, nullptr, skip_zeros);
  MFEM_VERIFY(kr, "Full-wave singular field energy requires a real bulk curl-curl "
                  "operator!");
  auto standard = std::make_unique<ParOperator>(std::move(kr), GetNDSpace());
  standard->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
  return std::make_unique<fem::singular::ParallelHybridEnrichedOperator>(
      std::move(standard), singular_domain_matrices.back().nd_curl_curl,
      nd_dbc_tdof_lists.back(), singular_nd_essential_true_dofs, diag_policy);
}

std::unique_ptr<Operator>
SpaceOperator::GetBulkMassMatrix(Operator::DiagonalPolicy diag_policy)
{
  MFEM_VERIFY(HasSingularEnrichment(),
              "The combined bulk mass matrix requires singular enrichment!");
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient fr(mat_op.MaxCeedAttribute());
  AddRealMassCoefficients(1.0, fr);
  constexpr bool skip_zeros = false;
  auto mr =
      AssembleOperator(GetNDSpace(), nullptr, &fr, nullptr, nullptr, nullptr, skip_zeros);
  MFEM_VERIFY(mr, "Full-wave singular divergence projection requires a positive real "
                  "bulk permittivity operator!");
  auto standard = std::make_unique<ParOperator>(std::move(mr), GetNDSpace());
  standard->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
  return std::make_unique<fem::singular::ParallelHybridEnrichedOperator>(
      std::move(standard), singular_domain_matrices.back().nd_mass,
      nd_dbc_tdof_lists.back(), singular_nd_essential_true_dofs, diag_policy);
}

std::unique_ptr<mfem::HypreParMatrix> SpaceOperator::GetBulkScalarDiffusionMatrix()
{
  MFEM_VERIFY(HasSingularEnrichment(),
              "The combined scalar diffusion matrix requires singular enrichment!");
  MaterialPropertyCoefficient epsilon(mat_op.MaxCeedAttribute());
  AddRealMassCoefficients(1.0, epsilon);
  constexpr bool skip_zeros = false;
  auto standard_levels = AssembleAuxOperators(GetH1Spaces(), &epsilon, nullptr, skip_zeros);
  MFEM_VERIFY(standard_levels.size() == GetH1Spaces().GetNumLevels() &&
                  standard_levels.back(),
              "Full-wave singular divergence projection requires a positive real "
              "scalar diffusion operator!");
  auto standard = ParOperator(std::move(standard_levels.back()), GetH1Space())
                      .StealParallelAssemble(skip_zeros);
  auto combined = fem::singular::BuildParallelEnrichedOperator(
      *standard, singular_domain_matrices.back().h1_diffusion);
  fem::singular::RemoveExplicitZeros(*combined);
  return combined;
}

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetExtraSystemMatrix(double omega, Operator::DiagonalPolicy diag_policy)
{
  return GetExtraSystemMatrix<OperType>(omega, diag_policy, /*include_wave_ports=*/true);
}

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetExtraSystemMatrix(double omega, Operator::DiagonalPolicy diag_policy,
                                    bool include_wave_ports)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient dfbr(mat_op.MaxCeedBdrAttribute()),
      dfbi(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
      fbi(mat_op.MaxCeedBdrAttribute());
  AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi, include_wave_ports);
  int empty[2] = {(dfbr.empty() && fbr.empty()), (dfbi.empty() && fbi.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> ar, ai;
  if (!empty[0])
  {
    ar = AssembleOperator(GetNDSpace(), nullptr, nullptr, &dfbr, &fbr, nullptr, skip_zeros);
  }
  if (!empty[1])
  {
    ai = AssembleOperator(GetNDSpace(), nullptr, nullptr, &dfbi, &fbi, nullptr, skip_zeros);
  }
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto A =
        std::make_unique<ComplexParOperator>(std::move(ar), std::move(ai), GetNDSpace());
    A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return A;
  }
  else
  {
    MFEM_VERIFY(!ai, "Unexpected imaginary part in GetExtraSystemMatrix<Operator>!");
    auto A = std::make_unique<ParOperator>(std::move(ar), GetNDSpace());
    A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return A;
  }
}

std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemMatrix(std::complex<double> omega,
                                    Operator::DiagonalPolicy diag_policy)
{
  // Complex-ω A2(λ) for the eigenmode nonlinear solve: identical assembly to the real-ω
  // overload but the frequency-dependent boundary terms (2nd-order ABC, surface
  // conductivity, rational impedance, numeric wave ports) are evaluated at the genuinely
  // complex frequency (ω = -i·λ). Always returns a ComplexOperator since these terms
  // carry a real-slot contribution at complex ω.
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient dfbr(mat_op.MaxCeedBdrAttribute()),
      dfbi(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
      fbi(mat_op.MaxCeedBdrAttribute());
  AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi);
  int empty[2] = {(dfbr.empty() && fbr.empty()), (dfbi.empty() && fbi.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> ar, ai;
  if (!empty[0])
  {
    ar = AssembleOperator(GetNDSpace(), nullptr, nullptr, &dfbr, &fbr, nullptr, skip_zeros);
  }
  if (!empty[1])
  {
    ai = AssembleOperator(GetNDSpace(), nullptr, nullptr, &dfbi, &fbi, nullptr, skip_zeros);
  }
  auto A = std::make_unique<ComplexParOperator>(std::move(ar), std::move(ai), GetNDSpace());
  A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
  return A;
}

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetWavePortBoundaryMassMatrix(int port_idx,
                                             Operator::DiagonalPolicy diag_policy)
{
  // Per-port μ⁻¹ boundary mass matrix, ω-independent — see
  // WavePortOperator::AddBoundaryMassBdrCoefficients. Pure imaginary part of A2(ω) when
  // the per-ω scalar i·k_n(ω) is reattached.
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
  wave_port_op.AddBoundaryMassBdrCoefficients(port_idx, fb);
  int empty = fb.empty();
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto m =
      AssembleOperator(GetNDSpace(), nullptr, nullptr, nullptr, &fb, nullptr, skip_zeros);
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto M_op = std::make_unique<ComplexParOperator>(nullptr, std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
  else
  {
    auto M_op = std::make_unique<ParOperator>(std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
}

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetFarfieldBoundaryCurlCurlMatrix(Operator::DiagonalPolicy diag_policy,
                                                 bool imag_slot)
{
  // ω-independent boundary curl-curl matrix M_ff for the 2nd-order farfield ABC, with unit
  // coefficient. By default stored on the REAL slot of the resulting ComplexParOperator so
  // that downstream BuildParSumOperator can scale it by an arbitrary complex coefficient
  // (the real-ω path uses i·(0.5/ω); the complex-λ path uses -0.5/λ). With imag_slot=true
  // it is placed on the IMAGINARY slot, matching the wave-port boundary-mass convention
  // (the i in i·f(ω)·M baked in) so circuit synthesis can fold it in uniformly. Returns
  // null if the farfield ABC order < 2 or it contributes no DoFs on any rank (the check is
  // collective, so the null contract is rank-uniform).
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient df(mat_op.MaxCeedBdrAttribute());
  farfield_op.AddExtraSystemBoundaryCurlCurlBdrCoefficients(1.0, df);
  int empty = df.empty();
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto m =
      AssembleOperator(GetNDSpace(), nullptr, nullptr, &df, nullptr, nullptr, skip_zeros);
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto M_op =
        imag_slot
            ? std::make_unique<ComplexParOperator>(nullptr, std::move(m), GetNDSpace())
            : std::make_unique<ComplexParOperator>(std::move(m), nullptr, GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
  else
  {
    MFEM_VERIFY(!imag_slot,
                "imag_slot is only meaningful for the ComplexOperator instantiation of "
                "GetFarfieldBoundaryCurlCurlMatrix!");
    auto M_op = std::make_unique<ParOperator>(std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
}

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetSurfaceConductivityBoundaryMatrix(int group_idx,
                                                    Operator::DiagonalPolicy diag_policy)
{
  // ω-independent unit-coefficient boundary mass A_σ for surface-conductivity group
  // group_idx, on the IMAGINARY slot (matching the wave-port convention). The full A2
  // contribution at frequency ω is (i·ω/Z_g(ω))·A_σ; circuit synthesis factors out A_σ here
  // and applies the scalar via a dispersion fit. Returns null if the group contributes no
  // DoFs on any rank (the check is collective, so the null contract is rank-uniform).
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
  surf_sigma_op.AddBoundaryMassBdrCoefficients(static_cast<std::size_t>(group_idx), fb);
  int empty = fb.empty();
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto m =
      AssembleOperator(GetNDSpace(), nullptr, nullptr, nullptr, &fb, nullptr, skip_zeros);
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto M_op = std::make_unique<ComplexParOperator>(nullptr, std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
  else
  {
    auto M_op = std::make_unique<ParOperator>(std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
}

template <typename OperType>
std::unique_ptr<OperType> SpaceOperator::GetRationalImpedanceBoundaryMassMatrix(
    int idx, Operator::DiagonalPolicy diag_policy, bool imag_slot)
{
  // λ-independent boundary mass matrix M_b for rational impedance boundary idx, with unit
  // coefficient (including crack scaling). By default stored on the REAL slot so that
  // downstream BuildParSumOperator can scale it by the arbitrary complex Robin coefficient
  // g(λ) from GetRationalImpedanceOp().EvalRobinCoefficient(idx, λ) (the NLEPS HYBRID
  // fit-or-freeze seed strategy). With imag_slot=true it is placed on the IMAGINARY slot,
  // matching the wave-port boundary-mass convention (the i in i·f(ω)·M baked in, with
  // f(ω) = g(iω)/i) so circuit synthesis can fold it in uniformly. Returns null if the
  // boundary contributes no DoFs on any rank (the check is collective, so the null contract
  // is rank-uniform).
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
  surf_rz_op.AddUnitBdrCoefficients(idx, fb);
  int empty = fb.empty();
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto m =
      AssembleOperator(GetNDSpace(), nullptr, nullptr, nullptr, &fb, nullptr, skip_zeros);
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto M_op =
        imag_slot
            ? std::make_unique<ComplexParOperator>(nullptr, std::move(m), GetNDSpace())
            : std::make_unique<ComplexParOperator>(std::move(m), nullptr, GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
  else
  {
    MFEM_VERIFY(!imag_slot,
                "imag_slot is only meaningful for the ComplexOperator instantiation of "
                "GetRationalImpedanceBoundaryMassMatrix!");
    auto M_op = std::make_unique<ParOperator>(std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
}

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetFloquetRobinBoundaryMassMatrix(int port_idx,
                                                 Operator::DiagonalPolicy diag_policy)
{
  // ω-independent µ⁻¹ boundary mass for a single Floquet port's Robin BC, on the
  // IMAGINARY slot (matching the wave-port / farfield convention: the i in i·γ₀·M is baked
  // in). The full online contribution is i·γ₀(ω)·M_floquet, with γ₀ the specular (0,0)
  // propagation constant. Returns null if the port contributes no DoFs on any rank (the
  // check is collective, so the null contract is rank-uniform).
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  const auto &port = floquet_port_op.GetPort(port_idx);
  MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
  MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  muinv_func.RestrictCoefficient(mat_op.GetCeedBdrAttributes(port.GetAttrList()));
  fb.AddCoefficient(muinv_func.GetAttributeToMaterial(), muinv_func.GetMaterialProperties(),
                    1.0);
  int empty = fb.empty();
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto m =
      AssembleOperator(GetNDSpace(), nullptr, nullptr, nullptr, &fb, nullptr, skip_zeros);
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto M_op = std::make_unique<ComplexParOperator>(nullptr, std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
  else
  {
    auto M_op = std::make_unique<ParOperator>(std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
}

template <typename OperType, typename ScalarType>
std::unique_ptr<OperType>
SpaceOperator::GetSystemMatrix(ScalarType a0, ScalarType a1, ScalarType a2,
                               const OperType *K, const OperType *C, const OperType *M,
                               const OperType *A2)
{
  if (HasSingularEnrichment())
  {
    if constexpr (std::is_same_v<OperType, ComplexOperator>)
    {
      auto combined = std::make_unique<SumComplexOperator>(GetNDTrueVSize());
      const auto add =
          [&combined](const ComplexOperator *op, std::complex<double> coefficient)
      {
        if (op && coefficient != 0.0)
        {
          combined->AddOperator(*op, coefficient);
          return true;
        }
        return false;
      };
      bool nonempty = add(K, a0);
      nonempty = add(C, a1) || nonempty;
      nonempty = add(M, a2) || nonempty;
      nonempty = add(A2, 1.0) || nonempty;
      MFEM_VERIFY(nonempty, "Full-wave singular system operator is empty!");
      return combined;
    }
    else
    {
      auto combined = std::make_unique<SumOperator>(GetNDTrueVSize());
      const auto add = [&combined](const Operator *op, double coefficient)
      {
        if (op && coefficient != 0.0)
        {
          combined->AddOperator(*op, coefficient);
          return true;
        }
        return false;
      };
      bool nonempty = add(K, a0);
      nonempty = add(C, a1) || nonempty;
      nonempty = add(M, a2) || nonempty;
      nonempty = add(A2, 1.0) || nonempty;
      MFEM_VERIFY(nonempty, "Full-wave singular real system operator is empty!");
      return combined;
    }
  }

  // When A2 is an abstract ComplexOperator (not a sparse ComplexParOperator), it cannot
  // participate in BuildParSumOperator. Build the sparse KCM sum separately and add A2 via
  // a non-owning SumComplexOperator.
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    if (A2 && !dynamic_cast<const ComplexParOperator *>(A2))
    {
      auto A_KCM = BuildParSumOperator({a0, a1, a2}, {K, C, M});
      return std::make_unique<SumComplexOperator>(std::move(A_KCM), *A2);
    }
  }
  return BuildParSumOperator({a0, a1, a2, ScalarType{1}}, {K, C, M, A2});
}

std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemOperator(double omega, Operator::DiagonalPolicy diag_policy)
{
  auto A2 = GetExtraSystemMatrix<ComplexOperator>(omega, diag_policy);
  auto F = floquet_port_op.GetExtraSystemOperator(omega);
  if (A2 && F)
  {
    return std::make_unique<SumComplexOperator>(std::move(A2), std::move(F));
  }
  if (A2)
  {
    return A2;
  }
  return F;
}

std::unique_ptr<Operator> SpaceOperator::GetInnerProductMatrix(double a0, double a2,
                                                               const ComplexOperator *K,
                                                               const ComplexOperator *M)
{
  if (HasSingularEnrichment())
  {
    auto combined = std::make_unique<SumOperator>(GetNDTrueVSize());
    bool nonempty = false;
    if (K && K->Real() && a0 != 0.0)
    {
      combined->AddOperator(*K->Real(), a0);
      nonempty = true;
    }
    if (M && M->Real() && a2 != 0.0)
    {
      combined->AddOperator(*M->Real(), a2);
      nonempty = true;
    }
    MFEM_VERIFY(nonempty, "Full-wave singular inner-product operator is empty!");
    return combined;
  }
  const auto *PtAP_K = (K) ? dynamic_cast<const ComplexParOperator *>(K) : nullptr;
  const auto *PtAP_M = (M) ? dynamic_cast<const ComplexParOperator *>(M) : nullptr;
  return BuildParSumOperator(
      {a0, a2}, {PtAP_K ? PtAP_K->Real() : nullptr, PtAP_M ? PtAP_M->Real() : nullptr});
}

namespace
{

template <typename OperType>
auto BuildLevelParOperator(std::unique_ptr<Operator> &&br, std::unique_ptr<Operator> &&bi,
                           const FiniteElementSpace &fespace);

template <>
auto BuildLevelParOperator<Operator>(std::unique_ptr<Operator> &&br,
                                     std::unique_ptr<Operator> &&bi,
                                     const FiniteElementSpace &fespace)
{
  MFEM_VERIFY(
      !bi,
      "Should not be constructing a real-valued ParOperator with non-zero imaginary part!");
  return std::make_unique<ParOperator>(std::move(br), fespace);
}

template <>
auto BuildLevelParOperator<ComplexOperator>(std::unique_ptr<Operator> &&br,
                                            std::unique_ptr<Operator> &&bi,
                                            const FiniteElementSpace &fespace)
{
  return std::make_unique<ComplexParOperator>(std::move(br), std::move(bi), fespace);
}

// L2-project a boundary VectorCoefficient onto the Nedelec tangent-trace space on the
// port surface: assemble the boundary mass M_bdr on port elements and CG-solve
// M_bdr * e = f. MFEM's ProjectBdrCoefficientTangent is a DOF-functional interpolant,
// not an L2 projection, so it returns a different vector (see PR #684 revert).
//
// TODO(future): Restrict the system to port-supported DOFs. The current full-T-space
// solve with DIAG_ONE on non-port DOFs is correct but pays for a trivially-identity
// block over most unknowns; a port-DOF-only system would be much smaller.
void ProjectBdrCoefficientViaMassSolve(SumVectorCoefficient &fb, const LumpedPortData &data,
                                       const MaterialOperator &mat_op,
                                       FiniteElementSpace &nd_fespace, MPI_Comm comm,
                                       Vector &result)
{
  // Assemble the boundary linear form f = ∫ φ_i · coeff dS (parallel-correct via P^T).
  mfem::LinearForm rhs(&nd_fespace.Get());
  rhs.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb));
  rhs.UseFastAssembly(false);
  rhs.UseDevice(false);
  rhs.Assemble();
  rhs.UseDevice(true);
  nd_fespace.GetProlongationMatrix()->MultTranspose(rhs, result);

  // Assemble boundary mass matrix M_bdr = ∫ φ_i · φ_j dS on port surface.
  MaterialPropertyCoefficient fb_mass(mat_op.MaxCeedBdrAttribute());
  for (const auto &elem : data.elems)
  {
    fb_mass.AddMaterialProperty(mat_op.GetCeedBdrAttributes(elem->GetAttrList()), 1.0);
  }
  BilinearForm m_bdr(nd_fespace);
  if (!fb_mass.empty())
  {
    m_bdr.AddBoundaryIntegrator<VectorFEMassIntegrator>(fb_mass);
  }
  auto M_bdr = std::make_unique<ParOperator>(m_bdr.Assemble(false), nd_fespace);

  // M_bdr is zero off the port, so set non-port DOFs essential with DIAG_ONE to get a
  // full-rank SPD system on the full T-space. non_port_tdof_list must outlive M_bdr.
  mfem::Array<int> attr_list;
  for (const auto &elem : data.elems)
  {
    attr_list.Append(elem->GetAttrList());
  }
  const auto &mesh = nd_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker;
  mesh::AttrToMarker(bdr_attr_max, attr_list, attr_marker);

  mfem::Array<int> port_tdof_list;
  nd_fespace.Get().GetEssentialTrueDofs(attr_marker, port_tdof_list);

  mfem::Array<int> non_port_tdof_list;
  {
    std::vector<bool> is_port(nd_fespace.GetTrueVSize(), false);
    for (int i = 0; i < port_tdof_list.Size(); i++)
    {
      is_port[port_tdof_list[i]] = true;
    }
    for (int i = 0; i < nd_fespace.GetTrueVSize(); i++)
    {
      if (!is_port[i])
      {
        non_port_tdof_list.Append(i);
      }
    }
  }
  M_bdr->SetEssentialTrueDofs(non_port_tdof_list, Operator::DIAG_ONE);

  // CG solve M_bdr * e = f entirely in T-vector space.
  // TODO: Make solver parameters configurable from IoData, or inherit other settings.
  auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(1.0e-14);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(200);
  auto jac = std::make_unique<JacobiSmoother<Operator>>(comm);
  auto ksp = std::make_unique<BaseKspSolver<Operator>>(std::move(pcg), std::move(jac));
  ksp->SetOperators(*M_bdr, *M_bdr);

  Vector sol(nd_fespace.GetTrueVSize());
  sol.UseDevice(true);
  sol = 0.0;
  ksp->Mult(result, sol);
  result = sol;
}

}  // namespace

template <typename A3Type>
void SpaceOperator::AssemblePreconditioner(
    std::complex<double> a0, std::complex<double> a1, std::complex<double> a2, A3Type a3,
    std::vector<std::unique_ptr<Operator>> &br_vec,
    std::vector<std::unique_ptr<Operator>> &br_aux_vec,
    std::vector<std::unique_ptr<Operator>> &bi_vec,
    std::vector<std::unique_ptr<Operator>> &bi_aux_vec)
{
  constexpr bool skip_zeros = false, assemble_q_data = false;
  MaterialPropertyCoefficient dfr(mat_op.MaxCeedAttribute()),
      dfi(mat_op.MaxCeedAttribute()), fr(mat_op.MaxCeedAttribute()),
      fi(mat_op.MaxCeedAttribute()), dfbr(mat_op.MaxCeedBdrAttribute()),
      dfbi(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
      fbi(mat_op.MaxCeedBdrAttribute()), fpi(mat_op.MaxCeedAttribute()),
      fpr(mat_op.MaxCeedAttribute());
  AddStiffnessCoefficients(a0.real(), dfr, fr);
  AddStiffnessCoefficients(a0.imag(), dfi, fi);
  AddStiffnessBdrCoefficients(a0.real(), fbr);
  AddStiffnessBdrCoefficients(a0.imag(), fbi);
  AddDampingCoefficients(a1.real(), fr);
  AddDampingCoefficients(a1.imag(), fi);
  AddDampingBdrCoefficients(a1.real(), fbr);
  AddDampingBdrCoefficients(a1.imag(), fbi);
  AddRealMassCoefficients(pc_mat_shifted ? std::abs(a2.real()) : a2.real(), fr);
  AddRealMassCoefficients(a2.imag(), fi);
  AddRealMassBdrCoefficients(pc_mat_shifted ? std::abs(a2.real()) : a2.real(), fbr);
  AddRealMassBdrCoefficients(a2.imag(), fbi);
  AddImagMassCoefficients(a2.real(), fi);
  AddImagMassCoefficients(-a2.imag(), fr);
  AddExtraSystemBdrCoefficients(a3, dfbr, dfbi, fbr, fbi);
  if (mat_op.HasFloquetFrequencyScaling())
  {
    // k₀-based tensors: cross terms scale with a1, mass term with a2.
    AddImagPeriodicCoefficients(a1.imag(), fpi);
    AddImagPeriodicCoefficients(a1.real(), fpr);
    AddRealPeriodicCoefficients(-a2.real(), fr);
    AddRealPeriodicCoefficients(-a2.imag(), fi);
  }
  else
  {
    AddRealPeriodicCoefficients(a0.real(), fr);
    AddRealPeriodicCoefficients(a0.imag(), fi);
    AddImagPeriodicCoefficients(a0.real(), fpi);
    AddImagPeriodicCoefficients(-a0.imag(), fpr);
  }
  int empty[2] = {
      (dfr.empty() && fr.empty() && dfbr.empty() && fbr.empty() && fpr.empty()),
      (dfi.empty() && fi.empty() && dfbi.empty() && fbi.empty() && fpi.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  if (!empty[0])
  {
    br_vec = AssembleOperators(GetNDSpaces(), &dfr, &fr, &dfbr, &fbr, &fpr, skip_zeros,
                               assemble_q_data);
    br_aux_vec =
        AssembleAuxOperators(GetH1Spaces(), &fr, &fbr, skip_zeros, assemble_q_data);
  }
  if (!empty[1])
  {
    bi_vec = AssembleOperators(GetNDSpaces(), &dfi, &fi, &dfbi, &fbi, &fpi, skip_zeros,
                               assemble_q_data);
    bi_aux_vec =
        AssembleAuxOperators(GetH1Spaces(), &fi, &fbi, skip_zeros, assemble_q_data);
  }
}

template <typename A3Type>
void SpaceOperator::AssemblePreconditioner(
    std::complex<double> a0, std::complex<double> a1, std::complex<double> a2, A3Type a3,
    std::vector<std::unique_ptr<Operator>> &br_vec,
    std::vector<std::unique_ptr<Operator>> &br_aux_vec)
{
  constexpr bool skip_zeros = false, assemble_q_data = false;
  MaterialPropertyCoefficient dfr(mat_op.MaxCeedAttribute()), fr(mat_op.MaxCeedAttribute()),
      dfbr(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute());
  AddStiffnessCoefficients(a0.real(), dfr, fr);
  AddStiffnessBdrCoefficients(a0.real(), fbr);
  AddDampingCoefficients(a1.imag(), fr);
  AddDampingBdrCoefficients(a1.imag(), fbr);
  AddAbsMassCoefficients(pc_mat_shifted ? std::abs(a2.real()) : a2.real(), fr);
  AddRealMassBdrCoefficients(pc_mat_shifted ? std::abs(a2.real()) : a2.real(), fbr);
  AddExtraSystemBdrCoefficients(a3, dfbr, dfbr, fbr, fbr);
  if (mat_op.HasFloquetFrequencyScaling())
  {
    AddRealPeriodicCoefficients(-(pc_mat_shifted ? std::abs(a2.real()) : a2.real()), fr);
  }
  else
  {
    AddRealPeriodicCoefficients(a0.real(), fr);
  }
  int empty = (dfr.empty() && fr.empty() && dfbr.empty() && fbr.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (!empty)
  {
    br_vec = AssembleOperators(GetNDSpaces(), &dfr, &fr, &dfbr, &fbr, nullptr, skip_zeros,
                               assemble_q_data);
    br_aux_vec =
        AssembleAuxOperators(GetH1Spaces(), &fr, &fbr, skip_zeros, assemble_q_data);
  }
}

void SpaceOperator::AssemblePreconditioner(
    double a0, double a1, double a2, double a3,
    std::vector<std::unique_ptr<Operator>> &br_vec,
    std::vector<std::unique_ptr<Operator>> &br_aux_vec)
{
  constexpr bool skip_zeros = false, assemble_q_data = false;
  MaterialPropertyCoefficient dfr(mat_op.MaxCeedAttribute()), fr(mat_op.MaxCeedAttribute()),
      dfbr(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute());
  AddStiffnessCoefficients(a0, dfr, fr);
  AddStiffnessBdrCoefficients(a0, fbr);
  AddDampingCoefficients(a1, fr);
  AddDampingBdrCoefficients(a1, fbr);
  AddAbsMassCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fr);
  AddRealMassBdrCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fbr);
  AddExtraSystemBdrCoefficients(a3, dfbr, dfbr, fbr, fbr);
  if (mat_op.HasFloquetFrequencyScaling())
  {
    AddRealPeriodicCoefficients(-(pc_mat_shifted ? std::abs(a2) : a2), fr);
  }
  else
  {
    AddRealPeriodicCoefficients(a0, fr);
  }
  int empty = (dfr.empty() && fr.empty() && dfbr.empty() && fbr.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (!empty)
  {
    br_vec = AssembleOperators(GetNDSpaces(), &dfr, &fr, &dfbr, &fbr, nullptr, skip_zeros,
                               assemble_q_data);
    br_aux_vec =
        AssembleAuxOperators(GetH1Spaces(), &fr, &fbr, skip_zeros, assemble_q_data);
  }
}

template <typename OperType, typename ScalarType, typename A3Type>
std::unique_ptr<OperType> SpaceOperator::GetPreconditionerMatrix(ScalarType a0,
                                                                 ScalarType a1,
                                                                 ScalarType a2, A3Type a3)
{
  if (HasSingularEnrichment())
  {
    const auto real_part = [](ScalarType value)
    {
      if constexpr (std::is_same_v<ScalarType, std::complex<double>>)
      {
        return value.real();
      }
      else
      {
        return value;
      }
    };
    const auto imaginary_part = [](ScalarType value)
    {
      if constexpr (std::is_same_v<ScalarType, std::complex<double>>)
      {
        return value.imag();
      }
      else
      {
        return 0.0;
      }
    };
    const bool complex_preconditioner =
        std::is_same_v<OperType, ComplexOperator> && !pc_mat_real;
    const double stiffness_coefficient = real_part(a0);
    const double mass_coefficient = [&]()
    {
      const double value = real_part(a2);
      return pc_mat_shifted ? std::abs(value) : value;
    }();
    const double damping_coefficient =
        complex_preconditioner
            ? real_part(a1)
            : (std::is_same_v<ScalarType, std::complex<double>> ? imaginary_part(a1)
                                                                : real_part(a1));

    if (print_prec_hdr)
    {
      Mpi::Print("\nAssembling combined singular multigrid hierarchy:\n");
    }
    const auto number_levels = GetNDSpaces().GetNumLevels();
    MFEM_VERIFY(number_levels == GetH1Spaces().GetNumLevels() &&
                    singular_domain_matrices.size() == number_levels &&
                    singular_domain_abs_matrices.size() == number_levels &&
                    singular_lumped_stiffness_matrices.size() == number_levels &&
                    singular_lumped_damping_matrices.size() == number_levels &&
                    singular_lumped_mass_matrices.size() == number_levels &&
                    singular_impedance_stiffness_matrices.size() == number_levels &&
                    singular_impedance_damping_matrices.size() == number_levels &&
                    singular_impedance_mass_matrices.size() == number_levels &&
                    singular_gradients.size() == number_levels &&
                    combined_nd_dbc_tdof_lists.size() == number_levels &&
                    combined_h1_dbc_tdof_lists.size() == number_levels,
                "Full-wave singular multigrid hierarchy is inconsistent!");

    std::vector<std::unique_ptr<Operator>> standard_operators(number_levels),
        standard_imaginary_operators(number_levels),
        standard_auxiliary_operators(number_levels),
        standard_imaginary_auxiliary_operators(number_levels);
    if constexpr (std::is_same_v<OperType, ComplexOperator>)
    {
      if (complex_preconditioner)
      {
        AssemblePreconditioner(a0, a1, a2, a3, standard_operators,
                               standard_auxiliary_operators, standard_imaginary_operators,
                               standard_imaginary_auxiliary_operators);
      }
      else
      {
        AssemblePreconditioner(a0, a1, a2, a3, standard_operators,
                               standard_auxiliary_operators);
      }
    }
    else
    {
      AssemblePreconditioner(a0, a1, a2, a3, standard_operators,
                             standard_auxiliary_operators);
    }
    MFEM_VERIFY(standard_operators.back(),
                "Full-wave singular preconditioner finest standard level is empty!");
    auto enrichment_gradient =
        fem::singular::BuildParallelEnrichmentGradient(GetComm(), *singular_numbering);
    auto hierarchy = std::make_unique<BaseMultigridOperator<OperType>>(number_levels);
    for (std::size_t level = 0; level < number_levels; level++)
    {
      auto &nd_space = GetNDSpaces().GetFESpaceAtLevel(level);
      auto &h1_space = GetH1Spaces().GetFESpaceAtLevel(level);
      auto real_enrichment =
          complex_preconditioner
              ? AddScaledSingularOperatorBlocks(
                    {{stiffness_coefficient, &singular_domain_matrices[level].nd_curl_curl},
                     {stiffness_coefficient, &singular_lumped_stiffness_matrices[level]},
                     {stiffness_coefficient, &singular_impedance_stiffness_matrices[level]},
                     {damping_coefficient, &singular_lumped_damping_matrices[level]},
                     {damping_coefficient, &singular_impedance_damping_matrices[level]},
                     {mass_coefficient, &singular_domain_matrices[level].nd_mass},
                     {-imaginary_part(a2), &singular_domain_imag_matrices[level].nd_mass},
                     {mass_coefficient, &singular_lumped_mass_matrices[level]},
                     {mass_coefficient, &singular_impedance_mass_matrices[level]}})
              : AddScaledSingularOperatorBlocks(
                    {{stiffness_coefficient, &singular_domain_matrices[level].nd_curl_curl},
                     {stiffness_coefficient, &singular_lumped_stiffness_matrices[level]},
                     {stiffness_coefficient, &singular_impedance_stiffness_matrices[level]},
                     {damping_coefficient, &singular_lumped_damping_matrices[level]},
                     {damping_coefficient, &singular_impedance_damping_matrices[level]},
                     {mass_coefficient, &singular_domain_abs_matrices[level].nd_mass},
                     {mass_coefficient, &singular_lumped_mass_matrices[level]},
                     {mass_coefficient, &singular_impedance_mass_matrices[level]}});
      const bool has_imaginary =
          complex_preconditioner && standard_imaginary_operators[level] != nullptr;
      fem::singular::ParallelSparseOperatorBlocks imaginary_enrichment;
      if (has_imaginary)
      {
        imaginary_enrichment = AddScaledSingularOperatorBlocks(
            {{imaginary_part(a0), &singular_domain_matrices[level].nd_curl_curl},
             {imaginary_part(a0), &singular_lumped_stiffness_matrices[level]},
             {imaginary_part(a0), &singular_impedance_stiffness_matrices[level]},
             {imaginary_part(a1), &singular_lumped_damping_matrices[level]},
             {imaginary_part(a1), &singular_impedance_damping_matrices[level]},
             {imaginary_part(a2), &singular_domain_matrices[level].nd_mass},
             {real_part(a2), &singular_domain_imag_matrices[level].nd_mass},
             {imaginary_part(a2), &singular_lumped_mass_matrices[level]},
             {imaginary_part(a2), &singular_impedance_mass_matrices[level]}});
      }
      const auto *standard_gradient =
          dynamic_cast<const ParOperator *>(&nd_space.GetDiscreteInterpolator(h1_space));
      MFEM_VERIFY(standard_gradient,
                  "Singular multigrid requires an assembled standard gradient!");
      auto real_auxiliary_enrichment =
          fem::singular::ProjectParallelSparseOperatorBlocksToH1(
              real_enrichment, standard_gradient->ParallelAssemble(), *enrichment_gradient);
      fem::singular::ParallelSparseOperatorBlocks imaginary_auxiliary_enrichment;
      if (has_imaginary)
      {
        imaginary_auxiliary_enrichment =
            fem::singular::ProjectParallelSparseOperatorBlocksToH1(
                imaginary_enrichment, standard_gradient->ParallelAssemble(),
                *enrichment_gradient);
      }
      const auto compact_enrichment =
          [](fem::singular::ParallelSparseOperatorBlocks &blocks)
      {
        HYPRE_BigInt removed = 0;
        removed += fem::singular::RemoveExplicitZeros(*blocks.standard_enrichment);
        removed += fem::singular::RemoveExplicitZeros(*blocks.enrichment_standard);
        removed += fem::singular::RemoveExplicitZeros(*blocks.enrichment_enrichment);
        return removed;
      };
      HYPRE_BigInt source_zeros_removed = compact_enrichment(real_enrichment);
      source_zeros_removed += compact_enrichment(real_auxiliary_enrichment);
      if (has_imaginary)
      {
        source_zeros_removed += compact_enrichment(imaginary_enrichment);
        source_zeros_removed += compact_enrichment(imaginary_auxiliary_enrichment);
      }
      const HYPRE_BigInt primary_sparse_nnz = real_enrichment.standard_enrichment->NNZ() +
                                              real_enrichment.enrichment_standard->NNZ() +
                                              real_enrichment.enrichment_enrichment->NNZ();
      const HYPRE_BigInt auxiliary_sparse_nnz =
          real_auxiliary_enrichment.standard_enrichment->NNZ() +
          real_auxiliary_enrichment.enrichment_standard->NNZ() +
          real_auxiliary_enrichment.enrichment_enrichment->NNZ();

      std::unique_ptr<Operator> primary_real, primary_imaginary, auxiliary_real,
          auxiliary_imaginary;
      mfem::HypreParMatrix *assembled_primary_real = nullptr;
      if (level == 0)
      {
        auto standard = ParOperator(std::move(standard_operators[level]), nd_space)
                            .StealParallelAssemble(false);
        source_zeros_removed += fem::singular::RemoveExplicitZeros(*standard);
        auto combined =
            fem::singular::BuildParallelEnrichedOperator(*standard, real_enrichment);
        combined->EliminateBC(combined_nd_dbc_tdof_lists[level], Operator::DIAG_ONE);
        const HYPRE_BigInt removed = fem::singular::RemoveExplicitZeros(*combined);
        assembled_primary_real = combined.get();
        primary_real = std::move(combined);
        if (has_imaginary)
        {
          auto standard_imaginary =
              ParOperator(std::move(standard_imaginary_operators[level]), nd_space)
                  .StealParallelAssemble(false);
          auto combined_imaginary = fem::singular::BuildParallelEnrichedOperator(
              *standard_imaginary, imaginary_enrichment);
          combined_imaginary->EliminateBC(combined_nd_dbc_tdof_lists[level],
                                          Operator::DIAG_ZERO);
          fem::singular::RemoveExplicitZeros(*combined_imaginary);
          primary_imaginary = std::move(combined_imaginary);
        }
        if (print_prec_hdr)
        {
          Mpi::Print(" Singular coarse sparse compaction removed {:d} source and {:d} "
                     "post-assembly explicit zeros\n",
                     source_zeros_removed, removed);
          PrintSparseMatrixStatistics("coarse combined after essential elimination",
                                      *assembled_primary_real);
        }
      }
      else
      {
        auto standard =
            std::make_unique<ParOperator>(std::move(standard_operators[level]), nd_space);
        standard->SetEssentialTrueDofs(nd_dbc_tdof_lists[level], Operator::DIAG_ONE);
        primary_real = std::make_unique<fem::singular::ParallelHybridEnrichedOperator>(
            std::move(standard), std::move(real_enrichment), nd_dbc_tdof_lists[level],
            singular_nd_essential_true_dofs, Operator::DIAG_ONE);
        if (has_imaginary)
        {
          auto standard_imaginary = std::make_unique<ParOperator>(
              std::move(standard_imaginary_operators[level]), nd_space);
          standard_imaginary->SetEssentialTrueDofs(nd_dbc_tdof_lists[level],
                                                   Operator::DIAG_ZERO);
          primary_imaginary =
              std::make_unique<fem::singular::ParallelHybridEnrichedOperator>(
                  std::move(standard_imaginary), std::move(imaginary_enrichment),
                  nd_dbc_tdof_lists[level], singular_nd_essential_true_dofs,
                  Operator::DIAG_ZERO);
        }
      }

      MFEM_VERIFY(standard_auxiliary_operators[level],
                  "Singular multigrid requires a standard H1 auxiliary operator!");
      auto standard_auxiliary = std::make_unique<ParOperator>(
          std::move(standard_auxiliary_operators[level]), h1_space);
      standard_auxiliary->SetEssentialTrueDofs(h1_dbc_tdof_lists[level],
                                               Operator::DIAG_ONE);
      auxiliary_real = std::make_unique<fem::singular::ParallelHybridEnrichedOperator>(
          std::move(standard_auxiliary), std::move(real_auxiliary_enrichment),
          h1_dbc_tdof_lists[level], singular_h1_essential_true_dofs, Operator::DIAG_ONE);
      if (has_imaginary)
      {
        MFEM_VERIFY(standard_imaginary_auxiliary_operators[level],
                    "Complex singular multigrid requires an imaginary standard H1 "
                    "auxiliary operator!");
        auto standard_imaginary_auxiliary = std::make_unique<ParOperator>(
            std::move(standard_imaginary_auxiliary_operators[level]), h1_space);
        standard_imaginary_auxiliary->SetEssentialTrueDofs(h1_dbc_tdof_lists[level],
                                                           Operator::DIAG_ZERO);
        auxiliary_imaginary =
            std::make_unique<fem::singular::ParallelHybridEnrichedOperator>(
                std::move(standard_imaginary_auxiliary),
                std::move(imaginary_auxiliary_enrichment), h1_dbc_tdof_lists[level],
                singular_h1_essential_true_dofs, Operator::DIAG_ZERO);
      }

      if (print_prec_hdr)
      {
        Mpi::Print(" Level {:d} (p = {:d}): {:d} combined ND unknowns, {}\n"
                   " Level {:d} (auxiliary): {:d} combined H1 unknowns, standard partial "
                   "assembly + {:d} sparse enrichment NNZ\n",
                   level, nd_space.GetMaxElementOrder(),
                   nd_space.GlobalTrueVSize() + singular_numbering->nd.global_size,
                   level == 0
                       ? fmt::format("{:d} assembled NNZ", assembled_primary_real->NNZ())
                       : fmt::format("standard partial assembly + {:d} sparse "
                                     "enrichment NNZ",
                                     primary_sparse_nnz),
                   level, h1_space.GlobalTrueVSize() + singular_numbering->h1.global_size,
                   auxiliary_sparse_nnz);
      }
      if constexpr (std::is_same_v<OperType, ComplexOperator>)
      {
        hierarchy->AddOperator(std::make_unique<ComplexWrapperOperator>(
            std::move(primary_real), std::move(primary_imaginary)));
        hierarchy->AddAuxiliaryOperator(std::make_unique<ComplexWrapperOperator>(
            std::move(auxiliary_real), std::move(auxiliary_imaginary)));
      }
      else
      {
        hierarchy->AddOperator(std::move(primary_real));
        hierarchy->AddAuxiliaryOperator(std::move(auxiliary_real));
      }
    }
    print_prec_hdr = false;
    return hierarchy;
  }

  // When partially assembled, the coarse operators can reuse the fine operator quadrature
  // data if the spaces correspond to the same mesh. When appropriate, we build the
  // preconditioner on all levels based on the actual complex-valued system matrix. The
  // coarse operator is always fully assembled.
  if (print_prec_hdr)
  {
    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
  MFEM_VERIFY(GetH1Spaces().GetNumLevels() == GetNDSpaces().GetNumLevels(),
              "Multigrid hierarchy mismatch for auxiliary space preconditioning!");

  const auto n_levels = GetNDSpaces().GetNumLevels();
  std::vector<std::unique_ptr<Operator>> br_vec(n_levels), bi_vec(n_levels),
      br_aux_vec(n_levels), bi_aux_vec(n_levels);
  if (std::is_same<OperType, ComplexOperator>::value && !pc_mat_real)
  {
    AssemblePreconditioner(a0, a1, a2, a3, br_vec, br_aux_vec, bi_vec, bi_aux_vec);
  }
  else
  {
    AssemblePreconditioner(a0, a1, a2, a3, br_vec, br_aux_vec);
  }

  auto B = std::make_unique<BaseMultigridOperator<OperType>>(n_levels);
  for (bool aux : {false, true})
  {
    for (std::size_t l = 0; l < n_levels; l++)
    {
      const auto &fespace_l =
          aux ? GetH1Spaces().GetFESpaceAtLevel(l) : GetNDSpaces().GetFESpaceAtLevel(l);
      const auto &dbc_tdof_lists_l = aux ? h1_dbc_tdof_lists[l] : nd_dbc_tdof_lists[l];
      auto &br_l = aux ? br_aux_vec[l] : br_vec[l];
      auto &bi_l = aux ? bi_aux_vec[l] : bi_vec[l];
      if (print_prec_hdr)
      {
        Mpi::Print(" Level {:d}{} (p = {:d}): {:d} unknowns", l, aux ? " (auxiliary)" : "",
                   fespace_l.GetMaxElementOrder(), fespace_l.GlobalTrueVSize());
        const auto *b_spm = dynamic_cast<const hypre::HypreCSRMatrix *>(br_l.get());
        if (!b_spm)
        {
          b_spm = dynamic_cast<const hypre::HypreCSRMatrix *>(bi_l.get());
        }
        if (b_spm)
        {
          HYPRE_BigInt nnz = b_spm->NNZ();
          Mpi::GlobalSum(1, &nnz, fespace_l.GetComm());
          Mpi::Print(", {:d} NNZ\n", nnz);
        }
        else
        {
          Mpi::Print("\n");
        }
      }
      auto B_l =
          BuildLevelParOperator<OperType>(std::move(br_l), std::move(bi_l), fespace_l);
      B_l->SetEssentialTrueDofs(dbc_tdof_lists_l, Operator::DiagonalPolicy::DIAG_ONE);
      if (aux)
      {
        B->AddAuxiliaryOperator(std::move(B_l));
      }
      else
      {
        B->AddOperator(std::move(B_l));
      }
    }
  }

  print_prec_hdr = false;
  return B;
}

void SpaceOperator::AddStiffnessCoefficients(double coeff, MaterialPropertyCoefficient &df,
                                             MaterialPropertyCoefficient &f)
{
  // Contribution from material permeability. In 2D, curl is scalar so the curl-curl
  // coefficient is scalar (1x1).
  df.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetCurlCurlInvPermeability(),
                    coeff);

  // Contribution for London superconductors.
  if (mat_op.HasLondonDepth())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetInvLondonDepth(), coeff);
  }
}

void SpaceOperator::AddStiffnessBdrCoefficients(double coeff,
                                                MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance and lumped ports (inductance).
  surf_z_op.AddStiffnessBdrCoefficients(coeff, fb);
  lumped_port_op.AddStiffnessBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddDampingCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  // Contribution for domain conductivity.
  if (mat_op.HasConductivity())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetConductivity(), coeff);
  }
}

void SpaceOperator::AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance, lumped ports, and absorbing
  // boundaries (resistance).
  farfield_op.AddDampingBdrCoefficients(coeff, fb);
  surf_z_op.AddDampingBdrCoefficients(coeff, fb);
  lumped_port_op.AddDampingBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddRealMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityReal(), coeff);
}

void SpaceOperator::AddRealMassBdrCoefficients(double coeff,
                                               MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance and lumped ports (capacitance).
  surf_z_op.AddMassBdrCoefficients(coeff, fb);
  lumped_port_op.AddMassBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddImagMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (mat_op.HasLossTangent())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityImag(), coeff);
  }
}

void SpaceOperator::AddAbsMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityAbs(), coeff);
}

void SpaceOperator::AddExtraSystemBdrCoefficients(double omega,
                                                  MaterialPropertyCoefficient &dfbr,
                                                  MaterialPropertyCoefficient &dfbi,
                                                  MaterialPropertyCoefficient &fbr,
                                                  MaterialPropertyCoefficient &fbi,
                                                  bool include_wave_ports)
{
  // Contribution for second-order farfield boundaries and finite conductivity boundaries.
  farfield_op.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi);
  surf_sigma_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
  surf_rz_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);

  // Contribution for numeric wave ports. Skipped when the caller will apply the wave-port
  // contribution separately via per-port operators (see GetWavePortBoundaryMassMatrix).
  if (include_wave_ports)
  {
    wave_port_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
  }

  // Contribution for Floquet ports (Robin BC).
  floquet_port_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
}

void SpaceOperator::AddExtraSystemBdrCoefficients(std::complex<double> omega,
                                                  MaterialPropertyCoefficient &dfbr,
                                                  MaterialPropertyCoefficient &dfbi,
                                                  MaterialPropertyCoefficient &fbr,
                                                  MaterialPropertyCoefficient &fbi)
{
  // Complex-ω overload for the eigenmode nonlinear solve and its matching preconditioner:
  // all frequency-dependent boundary terms (2nd-order farfield ABC, surface conductivity,
  // rational surface impedance, numeric wave ports) are evaluated at the genuinely complex
  // frequency (ω = -i·λ) so the assembled A2(λ) is the exact analytic continuation. For
  // real ω the operators' complex overloads reduce to their double overloads, except the
  // wave-port term, which additionally carries the line attenuation -Im(k_n)·M on the
  // real slot (the real-ω overload intentionally stamps only Re(k_n)).
  farfield_op.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi);
  surf_sigma_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
  surf_rz_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
  wave_port_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
}

void SpaceOperator::AddRealPeriodicCoefficients(double coeff,
                                                MaterialPropertyCoefficient &f)
{
  // Floquet periodicity contributions.
  if (mat_op.HasWaveVector())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetFloquetMass(), coeff);
  }
}

void SpaceOperator::AddImagPeriodicCoefficients(double coeff,
                                                MaterialPropertyCoefficient &f)
{
  // Floquet periodicity contributions.
  if (mat_op.HasWaveVector())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetFloquetCurl(), coeff);
  }
}

bool SpaceOperator::GetExcitationVector(int excitation_idx, Vector &RHS)
{
  // Time domain excitation vector.
  Vector standard_rhs(GetNDSpace().GetTrueVSize());
  standard_rhs.UseDevice(true);
  standard_rhs = 0.0;
  Vector singular_rhs;
  Vector *singular_rhs_ptr = nullptr;
  if (HasSingularEnrichment())
  {
    singular_rhs.SetSize(static_cast<int>(singular_numbering->nd.owned_size));
    singular_rhs.UseDevice(true);
    singular_rhs = 0.0;
    singular_rhs_ptr = &singular_rhs;
  }
  bool nnz = AddExcitationVector1Internal(excitation_idx, standard_rhs, singular_rhs_ptr);

  RHS.SetSize(GetNDTrueVSize());
  RHS.UseDevice(true);
  RHS = 0.0;
  for (int i = 0; i < standard_rhs.Size(); i++)
  {
    RHS[i] = standard_rhs[i];
  }
  for (int i = 0; i < singular_rhs.Size(); i++)
  {
    RHS[standard_rhs.Size() + i] = singular_rhs[i];
  }
  linalg::SetSubVector(RHS, GetCombinedNDDbcTDofList(), 0.0);
  return nnz;
}

bool SpaceOperator::GetExcitationVector(int excitation_idx, double omega,
                                        ComplexVector &RHS)
{
  // Frequency domain excitation vector: RHS = iω RHS1 + RHS2(ω).
  ComplexVector standard_rhs(GetNDSpace().GetTrueVSize());
  standard_rhs.UseDevice(true);
  standard_rhs = 0.0;
  Vector singular_rhs;
  Vector *singular_rhs_ptr = nullptr;
  if (HasSingularEnrichment())
  {
    singular_rhs.SetSize(static_cast<int>(singular_numbering->nd.owned_size));
    singular_rhs.UseDevice(true);
    singular_rhs = 0.0;
    singular_rhs_ptr = &singular_rhs;
  }
  bool nnz1 =
      AddExcitationVector1Internal(excitation_idx, standard_rhs.Real(), singular_rhs_ptr);
  standard_rhs *= 1i * omega;
  bool nnz2 = AddExcitationVector2Internal(excitation_idx, omega, standard_rhs);

  RHS.SetSize(GetNDTrueVSize());
  RHS.UseDevice(true);
  RHS = 0.0;
  for (int i = 0; i < standard_rhs.Size(); i++)
  {
    RHS.Real()[i] = standard_rhs.Real()[i];
    RHS.Imag()[i] = standard_rhs.Imag()[i];
  }
  for (int i = 0; i < singular_rhs.Size(); i++)
  {
    RHS.Imag()[standard_rhs.Size() + i] = omega * singular_rhs[i];
  }
  linalg::SetSubVector(RHS, GetCombinedNDDbcTDofList(), 0.0);
  return nnz1 || nnz2;
}

void SpaceOperator::GetLumpedPortExcitationVectorPrimaryEt(int port_idx,
                                                           ComplexVector &Et_primary)
{
  const auto &data = GetLumpedPortOp().GetPort(port_idx);

  SumVectorCoefficient fb(GetMesh().SpaceDimension());
  for (const auto &elem : data.elems)
  {
    const double Rs = 1.0 * data.GetToSquare(*elem);
    const double Einc = std::sqrt(
        Rs / (elem->GetGeometryWidth() * elem->GetGeometryLength() * data.elems.size()));
    fb.AddCoefficient(elem->GetModeCoefficient(Einc));
  }

  Et_primary.SetSize(GetNDSpace().GetTrueVSize());
  Et_primary.UseDevice(true);
  Et_primary = 0.0;

  // L2-project the port mode onto the ND tangent-trace space; DOF interpolation
  // gives a different vector and breaks the W-norm = 1 invariant Einc relies on.
  ProjectBdrCoefficientViaMassSolve(fb, data, mat_op, GetNDSpace(), GetComm(),
                                    Et_primary.Real());

  linalg::SetSubVector(Et_primary.Real(), GetNDDbcTDofLists().back(), 0.0);
}

void SpaceOperator::GetLumpedPortExcitationVectorPrimaryHtcn(int port_idx,
                                                             ComplexVector &Htcn_primary)
{
  const auto &data = lumped_port_op.GetPort(port_idx);

  SumVectorCoefficient fb(GetMesh().SpaceDimension());
  for (const auto &elem : data.elems)
  {
    const double Rs = 1.0 * data.GetToSquare(*elem);
    const double Hinc = 1.0 / std::sqrt(Rs * elem->GetGeometryWidth() *
                                        elem->GetGeometryLength() * data.elems.size());
    fb.AddCoefficient(elem->GetModeCoefficient(Hinc));
  }

  Htcn_primary.SetSize(GetNDSpace().GetTrueVSize());
  Htcn_primary.UseDevice(true);
  Htcn_primary = 0.0;

  // L2-project the port mode onto the ND tangent-trace space; DOF interpolation
  // gives a different vector and breaks the W-norm = 1 invariant Einc relies on.
  ProjectBdrCoefficientViaMassSolve(fb, data, mat_op, GetNDSpace(), GetComm(),
                                    Htcn_primary.Real());

  linalg::SetSubVector(Htcn_primary.Real(), GetNDDbcTDofLists().back(), 0.0);
}

void SpaceOperator::GetWavePortFieldVectorPrimaryEt(int port_idx, double omega_ref,
                                                    ComplexVector &Et_primary)
{
  // Trigger (or refresh) the cross-section EVP at omega_ref. WavePortData::Initialize
  // caches by omega so calling this is cheap if already initialised.
  wave_port_op.GetWavePortKn(port_idx, omega_ref);
  const auto &data = wave_port_op.GetPort(port_idx);

  // Project the modal tangential E-field onto the parent ND space, restricted to the
  // port boundary attributes. Real and imaginary parts are projected separately, since
  // the full waveport mode is generally complex.
  Et_primary.SetSize(GetNDSpace().GetTrueVSize());
  Et_primary.UseDevice(true);
  Et_primary = 0.0;

  const auto &mesh = GetNDSpace().GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker;
  mesh::AttrToMarker(bdr_attr_max, data.GetAttrList(), attr_marker);

  GridFunction rhs_re(GetNDSpace());
  GridFunction rhs_im(GetNDSpace());
  rhs_re = 0.0;
  rhs_im = 0.0;
  {
    auto fb_re = data.GetModeFieldCoefficientReal(1.0);
    rhs_re.Real().ProjectBdrCoefficientTangent(*fb_re, attr_marker);
  }
  {
    auto fb_im = data.GetModeFieldCoefficientImag(1.0);
    rhs_im.Real().ProjectBdrCoefficientTangent(*fb_im, attr_marker);
  }
  GetNDSpace().GetRestrictionMatrix()->Mult(rhs_re.Real(), Et_primary.Real());
  GetNDSpace().GetRestrictionMatrix()->Mult(rhs_im.Real(), Et_primary.Imag());

  linalg::SetSubVector(Et_primary.Real(), GetNDDbcTDofLists().back(), 0.0);
  linalg::SetSubVector(Et_primary.Imag(), GetNDDbcTDofLists().back(), 0.0);
}

bool SpaceOperator::GetExcitationVector1(int excitation_idx, ComplexVector &RHS1)
{
  // Assemble the frequency domain excitation term with linear frequency dependence
  // (coefficient iω, see GetExcitationVector above, is accounted for later).
  Vector standard_rhs(GetNDSpace().GetTrueVSize());
  standard_rhs.UseDevice(true);
  standard_rhs = 0.0;
  Vector singular_rhs;
  Vector *singular_rhs_ptr = nullptr;
  if (HasSingularEnrichment())
  {
    singular_rhs.SetSize(static_cast<int>(singular_numbering->nd.owned_size));
    singular_rhs.UseDevice(true);
    singular_rhs = 0.0;
    singular_rhs_ptr = &singular_rhs;
  }
  bool nnz1 = AddExcitationVector1Internal(excitation_idx, standard_rhs, singular_rhs_ptr);

  RHS1.SetSize(GetNDTrueVSize());
  RHS1.UseDevice(true);
  RHS1 = 0.0;
  for (int i = 0; i < standard_rhs.Size(); i++)
  {
    RHS1.Real()[i] = standard_rhs[i];
  }
  for (int i = 0; i < singular_rhs.Size(); i++)
  {
    RHS1.Real()[standard_rhs.Size() + i] = singular_rhs[i];
  }
  linalg::SetSubVector(RHS1, GetCombinedNDDbcTDofList(), 0.0);
  return nnz1;
}

bool SpaceOperator::GetExcitationVector2(int excitation_idx, double omega,
                                         ComplexVector &RHS2)
{
  ComplexVector standard_rhs(GetNDSpace().GetTrueVSize());
  standard_rhs.UseDevice(true);
  standard_rhs = 0.0;
  bool nnz2 = AddExcitationVector2Internal(excitation_idx, omega, standard_rhs);

  RHS2.SetSize(GetNDTrueVSize());
  RHS2.UseDevice(true);
  RHS2 = 0.0;
  for (int i = 0; i < standard_rhs.Size(); i++)
  {
    RHS2.Real()[i] = standard_rhs.Real()[i];
    RHS2.Imag()[i] = standard_rhs.Imag()[i];
  }
  linalg::SetSubVector(RHS2, GetCombinedNDDbcTDofList(), 0.0);
  return nnz2;
}

bool SpaceOperator::AddExcitationVector1Internal(int excitation_idx, Vector &RHS1,
                                                 Vector *singular_RHS1)
{
  // Assemble the time domain excitation -g'(t) J or frequency domain excitation -iω J.
  // The g'(t) or iω factors are not accounted for here, they are accounted for in the time
  // integration or frequency sweep later.
  MFEM_VERIFY(RHS1.Size() == GetNDSpace().GetTrueVSize(),
              "Invalid T-vector size for AddExcitationVector1Internal!");
  MFEM_VERIFY(!singular_RHS1 || (HasSingularEnrichment() &&
                                 singular_RHS1->Size() ==
                                     static_cast<int>(singular_numbering->nd.owned_size)),
              "Invalid singular T-vector size for AddExcitationVector1Internal!");

  // Boundary sources
  SumVectorCoefficient fb(GetMesh().SpaceDimension());
  SumVectorCoefficient singular_fb(GetMesh().SpaceDimension());
  lumped_port_op.AddExcitationBdrCoefficients(excitation_idx, fb);
  lumped_port_op.AddExcitationBdrCoefficients(excitation_idx, singular_fb);
  surf_j_op.AddExcitationBdrCoefficients(fb);  // No excitation_idx — currently in all

  // Domain sources (current dipoles) - use integrator-based approach
  bool has_current_dipoles = !current_dipole_op.Empty();

  int empty[2] = {(fb.empty()), (!has_current_dipoles)};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return false;
  }

  mfem::ParLinearForm rhs1(&GetNDSpace().Get());

  // Add boundary integrators
  if (!empty[0])
  {
    rhs1.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb));
  }

  // Add domain integrators for current dipoles
  if (!empty[1])
  {
    current_dipole_op.AddExcitationDomainIntegrators(rhs1);
  }

  rhs1.UseFastAssembly(false);
  rhs1.UseDevice(false);
  rhs1.Assemble();
  rhs1.UseDevice(true);
  GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs1, RHS1);
  int singular_empty = singular_fb.empty();
  Mpi::GlobalMin(1, &singular_empty, GetComm());
  if (singular_RHS1 && !singular_empty)
  {
    if (singular_dofs)
    {
      fem::singular::AssembleParallelNDBoundaryLinearForm(
          *singular_dofs, *singular_numbering, GetNDSpace().Get(), singular_fb,
          singular_assembly_options, *singular_RHS1);
    }
    else
    {
      MFEM_VERIFY(triangle_singular_dofs,
                  "Singular excitation assembly requires simplex ND topology!");
      fem::singular::AssembleParallelNDBoundaryLinearForm(
          *triangle_singular_dofs, *singular_numbering, GetNDSpace().Get(), singular_fb,
          singular_assembly_options, *singular_RHS1);
    }
  }
  return true;
}

bool SpaceOperator::AddExcitationVector2Internal(int excitation_idx, double omega,
                                                 ComplexVector &RHS2)
{
  // Assemble the contribution of wave ports to the frequency domain excitation term at
  // the specified frequency.
  MFEM_VERIFY(RHS2.Size() == GetNDSpace().GetTrueVSize(),
              "Invalid T-vector size for AddExcitationVector2Internal!");
  // Floquet port excitation: directly adds to RHS (not via boundary linear form).
  bool nnz_floquet = floquet_port_op.AddExcitationVector(excitation_idx, omega, RHS2);

  SumVectorCoefficient fbr(GetMesh().SpaceDimension()), fbi(GetMesh().SpaceDimension());
  wave_port_op.AddExcitationBdrCoefficients(excitation_idx, omega, fbr, fbi);
  int empty = (fbr.empty() && fbi.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return nnz_floquet;
  }
  {
    mfem::LinearForm rhs2(&GetNDSpace().Get());
    rhs2.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr));
    rhs2.UseFastAssembly(false);
    rhs2.UseDevice(false);
    rhs2.Assemble();
    rhs2.UseDevice(true);
    GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs2, RHS2.Real());
  }
  {
    mfem::LinearForm rhs2(&GetNDSpace().Get());
    rhs2.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi));
    rhs2.UseFastAssembly(false);
    rhs2.UseDevice(false);
    rhs2.Assemble();
    rhs2.UseDevice(true);
    GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs2, RHS2.Imag());
  }
  return true;
}

std::unique_ptr<Vector>
SpaceOperator::BuildSingularBoundaryFunctional(SumVectorCoefficient &coefficient,
                                               mfem::Array<int> attribute_marker)
{
  MFEM_VERIFY(HasSingularEnrichment(),
              "Singular boundary functional requires singular enrichment!");

  mfem::ParLinearForm standard_form(&GetNDSpace().Get());
  standard_form.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(coefficient),
                                      attribute_marker);
  standard_form.UseFastAssembly(false);
  standard_form.UseDevice(false);
  standard_form.Assemble();
  standard_form.UseDevice(true);
  Vector standard_true_dofs(GetNDSpace().GetTrueVSize());
  standard_true_dofs.UseDevice(true);
  standard_true_dofs = 0.0;
  GetNDSpace().GetProlongationMatrix()->MultTranspose(standard_form, standard_true_dofs);

  Vector enrichment_true_dofs;
  if (singular_dofs)
  {
    fem::singular::AssembleParallelNDBoundaryLinearForm(
        *singular_dofs, *singular_numbering, GetNDSpace().Get(), coefficient,
        singular_assembly_options, enrichment_true_dofs, &attribute_marker);
  }
  else
  {
    MFEM_VERIFY(triangle_singular_dofs,
                "Singular boundary functional requires simplex ND topology!");
    fem::singular::AssembleParallelNDBoundaryLinearForm(
        *triangle_singular_dofs, *singular_numbering, GetNDSpace().Get(), coefficient,
        singular_assembly_options, enrichment_true_dofs, &attribute_marker);
  }

  auto functional = std::make_unique<Vector>(GetNDTrueVSize());
  functional->UseDevice(true);
  *functional = 0.0;
  for (int i = 0; i < standard_true_dofs.Size(); i++)
  {
    (*functional)[i] = standard_true_dofs[i];
  }
  for (int i = 0; i < enrichment_true_dofs.Size(); i++)
  {
    (*functional)[standard_true_dofs.Size() + i] = enrichment_true_dofs[i];
  }
  return functional;
}

std::complex<double>
SpaceOperator::ApplySingularBoundaryFunctional(const Vector &functional,
                                               const ComplexVector &field) const
{
  MFEM_VERIFY(HasSingularEnrichment() && field.Size() == GetNDTrueVSize() &&
                  functional.Size() == field.Size(),
              "Singular boundary functional requires a combined ND true-DOF field!");
  std::complex<double> result = 0.0;
  const auto *values = functional.HostRead();
  const auto *field_real = field.Real().HostRead();
  const auto *field_imaginary = field.Imag().HostRead();
  for (int i = 0; i < functional.Size(); i++)
  {
    result += values[i] * std::complex<double>(field_real[i], field_imaginary[i]);
  }
  Mpi::GlobalSum(1, &result, GetComm());
  return result;
}

void SpaceOperator::CacheSingularLumpedPortFunctionals(bool include_sparameters)
{
  MFEM_VERIFY(HasSingularEnrichment(),
              "Singular lumped-port functionals require singular enrichment!");
  for (const auto &[port_idx, port] : lumped_port_op)
  {
    mfem::Array<int> attributes;
    for (const auto &element : port.elems)
    {
      attributes.Append(element->GetAttrList());
    }
    mfem::Array<int> attribute_marker;
    const auto &mesh = GetMesh().Get();
    const int maximum_attribute =
        mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mesh::AttrToMarker(maximum_attribute, attributes, attribute_marker);
    if (singular_lumped_voltage_functionals.count(port_idx) == 0)
    {
      SumVectorCoefficient coefficient(GetMesh().SpaceDimension());
      port.AddVoltageBdrCoefficients(coefficient);
      singular_lumped_voltage_functionals.emplace(
          port_idx, BuildSingularBoundaryFunctional(coefficient, attribute_marker));
    }
    if (include_sparameters && singular_lumped_sparameter_functionals.count(port_idx) == 0)
    {
      SumVectorCoefficient coefficient(GetMesh().SpaceDimension());
      port.AddSParameterBdrCoefficients(coefficient);
      singular_lumped_sparameter_functionals.emplace(
          port_idx, BuildSingularBoundaryFunctional(coefficient, attribute_marker));
    }
  }
}

std::complex<double> SpaceOperator::GetSingularLumpedPortVoltage(int port_idx,
                                                                 const ComplexVector &field)
{
  if (singular_lumped_voltage_functionals.count(port_idx) == 0)
  {
    CacheSingularLumpedPortFunctionals(false);
  }
  return ApplySingularBoundaryFunctional(*singular_lumped_voltage_functionals.at(port_idx),
                                         field);
}

std::complex<double>
SpaceOperator::GetSingularLumpedPortSParameter(int port_idx, const ComplexVector &field)
{
  if (singular_lumped_sparameter_functionals.count(port_idx) == 0)
  {
    CacheSingularLumpedPortFunctionals(true);
  }
  return ApplySingularBoundaryFunctional(
      *singular_lumped_sparameter_functionals.at(port_idx), field);
}

void SpaceOperator::GetConstantInitialVector(ComplexVector &v)
{
  v.SetSize(GetNDTrueVSize());
  v.UseDevice(true);
  v = 1.0;
  linalg::SetSubVector(v.Real(), GetCombinedNDDbcTDofList(), 0.0);
}

void SpaceOperator::GetRandomInitialVector(ComplexVector &v)
{
  v.SetSize(GetNDTrueVSize());
  v.UseDevice(true);
  linalg::SetRandom(GetNDSpace().GetComm(), v);
  linalg::SetSubVector(v, GetCombinedNDDbcTDofList(), 0.0);
}

template std::unique_ptr<Operator>
    SpaceOperator::GetStiffnessMatrix(Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
    SpaceOperator::GetStiffnessMatrix(Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
    SpaceOperator::GetDampingMatrix(Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
    SpaceOperator::GetDampingMatrix(Operator::DiagonalPolicy);

template std::unique_ptr<Operator> SpaceOperator::GetMassMatrix(Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
    SpaceOperator::GetMassMatrix(Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
SpaceOperator::GetExtraSystemMatrix(double, Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemMatrix(double, Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
SpaceOperator::GetExtraSystemMatrix(double, Operator::DiagonalPolicy, bool);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemMatrix(double, Operator::DiagonalPolicy, bool);

template std::unique_ptr<Operator>
SpaceOperator::GetWavePortBoundaryMassMatrix(int, Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetWavePortBoundaryMassMatrix(int, Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
SpaceOperator::GetFarfieldBoundaryCurlCurlMatrix(Operator::DiagonalPolicy, bool);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetFarfieldBoundaryCurlCurlMatrix(Operator::DiagonalPolicy, bool);

template std::unique_ptr<Operator>
SpaceOperator::GetSurfaceConductivityBoundaryMatrix(int, Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetSurfaceConductivityBoundaryMatrix(int, Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
SpaceOperator::GetRationalImpedanceBoundaryMassMatrix(int, Operator::DiagonalPolicy, bool);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetRationalImpedanceBoundaryMassMatrix(int, Operator::DiagonalPolicy, bool);

template std::unique_ptr<Operator>
SpaceOperator::GetFloquetRobinBoundaryMassMatrix(int, Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetFloquetRobinBoundaryMassMatrix(int, Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
SpaceOperator::GetSystemMatrix<Operator, double>(double, double, double, const Operator *,
                                                 const Operator *, const Operator *,
                                                 const Operator *);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetSystemMatrix<ComplexOperator, std::complex<double>>(
    std::complex<double>, std::complex<double>, std::complex<double>,
    const ComplexOperator *, const ComplexOperator *, const ComplexOperator *,
    const ComplexOperator *);

template std::unique_ptr<Operator>
SpaceOperator::GetPreconditionerMatrix<Operator, double, double>(double, double, double,
                                                                 double);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetPreconditionerMatrix<ComplexOperator, std::complex<double>, double>(
    std::complex<double>, std::complex<double>, std::complex<double>, double);
template std::unique_ptr<ComplexOperator>
    SpaceOperator::GetPreconditionerMatrix<ComplexOperator, std::complex<double>,
                                           std::complex<double>>(std::complex<double>,
                                                                 std::complex<double>,
                                                                 std::complex<double>,
                                                                 std::complex<double>);

}  // namespace palace
