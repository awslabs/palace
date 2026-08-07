// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "laplaceoperator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
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

LaplaceOperator::LaplaceOperator(
    const config::BoundaryData &boundaries, const config::SolverData &solver,
    const std::vector<config::MaterialData> &materials, ProblemType problem_type,
    const std::vector<std::unique_ptr<Mesh>> &mesh,
    const fem::singular::FeatureTopology *singular_features,
    const std::vector<fem::singular::GlobalVertexId> *source_vertex_ids,
    const fem::singular::TriangleFeatureTopology *triangle_singular_features)
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
    source_attr_lists(ConstructSources(boundaries.terminal)),
    singular_features(singular_features),
    triangle_singular_features(triangle_singular_features),
    source_vertex_ids(source_vertex_ids), standard_order(solver.order),
    singular_order(solver.singular_elements.order),
    singular_assembly_options{solver.singular_elements.quadrature_order,
                              solver.singular_elements.abs_tol,
                              solver.singular_elements.rel_tol,
                              solver.singular_elements.max_subdivisions,
                              solver.singular_elements.UsesFixedSubdivision(),
                              solver.singular_elements.subdivisions}
{
  const bool has_tetrahedral_features = singular_features != nullptr;
  const bool has_triangular_features = triangle_singular_features != nullptr;
  MFEM_VERIFY(!(has_tetrahedral_features && has_triangular_features),
              "A Laplace operator cannot use both triangular and tetrahedral singular "
              "features!");
  MFEM_VERIFY((has_tetrahedral_features || has_triangular_features) ==
                  (source_vertex_ids != nullptr),
              "Singular feature topology and source vertex IDs must be supplied together!");
  MFEM_VERIFY(
      !HasSingularEnrichment() ||
          (solver.singular_elements.Enabled() &&
           source_vertex_ids->size() ==
               static_cast<std::size_t>(mesh.back()->Get().GetNV()) &&
           ((!has_tetrahedral_features &&
             triangle_singular_features->elements.size() ==
                 static_cast<std::size_t>(mesh.back()->GetNE()) &&
             mesh.back()->Dimension() == 2 && solver.singular_elements.order == 1) ||
            (!has_triangular_features &&
             singular_features->elements.size() ==
                 static_cast<std::size_t>(mesh.back()->GetNE()) &&
             mesh.back()->Dimension() == 3))),
      "Singular feature topology does not match the electrostatic solve mesh!");
  if (HasSingularEnrichment())
  {
    for (std::size_t level = 0; level < h1_fespaces.GetNumLevels(); level++)
    {
      MFEM_VERIFY(&h1_fespaces.GetFESpaceAtLevel(level).GetMesh() == &GetMesh(),
                  "Singular electrostatic multigrid currently supports polynomial "
                  "levels on one mesh!");
    }
  }

  // Print essential BC information.
  if (dbc_attr.Size())
  {
    Mpi::Print("\nConfiguring Dirichlet BC at attributes:\n");
    utils::PrettyPrint(dbc_attr);
  }
}

LaplaceOperator::LaplaceOperator(
    const IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh,
    const fem::singular::FeatureTopology *singular_features,
    const std::vector<fem::singular::GlobalVertexId> *source_vertex_ids,
    const fem::singular::TriangleFeatureTopology *triangle_singular_features)
  : LaplaceOperator(iodata.boundaries, iodata.solver, iodata.domains.materials,
                    iodata.problem.type, mesh, singular_features, source_vertex_ids,
                    triangle_singular_features)
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

SingularOperatorDiagnostics BuildSingularOperatorDiagnostics(
    const mfem::HypreParMatrix &matrix, int standard_local_size,
    const fem::singular::ParallelDofNumbering &numbering,
    const fem::singular::AdaptiveAssemblyOptions &options, int standard_order,
    int singular_order, std::size_t local_quadrature_leaf_count,
    int local_quadrature_maximum_depth, std::size_t local_duffy_table_entries,
    std::size_t local_duffy_cache_hits)
{
  const auto hypre_max = static_cast<std::size_t>(std::numeric_limits<HYPRE_BigInt>::max());
  if (local_quadrature_leaf_count > hypre_max || local_duffy_table_entries > hypre_max ||
      local_duffy_cache_hits > hypre_max)
  {
    throw std::overflow_error("Singular quadrature diagnostics overflow!");
  }
  HYPRE_BigInt quadrature_leaf_count =
      static_cast<HYPRE_BigInt>(local_quadrature_leaf_count);
  HYPRE_BigInt duffy_table_entries = static_cast<HYPRE_BigInt>(local_duffy_table_entries);
  HYPRE_BigInt duffy_cache_hits = static_cast<HYPRE_BigInt>(local_duffy_cache_hits);
  int quadrature_maximum_depth = local_quadrature_maximum_depth;
  Mpi::GlobalSum(1, &quadrature_leaf_count, matrix.GetComm());
  Mpi::GlobalMax(1, &quadrature_maximum_depth, matrix.GetComm());
  Mpi::GlobalMax(1, &duffy_table_entries, matrix.GetComm());
  Mpi::GlobalSum(1, &duffy_cache_hits, matrix.GetComm());

  mfem::Vector diagonal;
  matrix.GetDiag(diagonal);
  bool valid =
      standard_local_size >= 0 &&
      numbering.h1.owned_size <= std::numeric_limits<int>::max() &&
      diagonal.Size() == standard_local_size + static_cast<int>(numbering.h1.owned_size);
  std::array<double, 2> diagonal_minimum{std::numeric_limits<double>::infinity(),
                                         std::numeric_limits<double>::infinity()};
  std::array<double, 2> diagonal_maximum{0.0, 0.0};
  std::array<HYPRE_BigInt, 2> diagonal_count{0, 0};
  if (valid)
  {
    for (int i = 0; i < diagonal.Size(); i++)
    {
      const double value = diagonal[i];
      if (!std::isfinite(value) || !(value > 0.0))
      {
        valid = false;
        break;
      }
      const int block = (i < standard_local_size) ? 0 : 1;
      diagonal_minimum[block] = std::min(diagonal_minimum[block], value);
      diagonal_maximum[block] = std::max(diagonal_maximum[block], value);
      diagonal_count[block]++;
    }
  }
  Mpi::GlobalAnd(1, &valid, matrix.GetComm());
  if (!valid)
  {
    throw std::runtime_error(
        "Singular stiffness matrix has invalid local dimensions or diagonal entries!");
  }
  Mpi::GlobalMin(2, diagonal_minimum.data(), matrix.GetComm());
  Mpi::GlobalMax(2, diagonal_maximum.data(), matrix.GetComm());
  Mpi::GlobalSum(2, diagonal_count.data(), matrix.GetComm());
  if (diagonal_count[0] != matrix.GetGlobalNumRows() - numbering.h1.global_size ||
      diagonal_count[1] != numbering.h1.global_size ||
      !std::all_of(diagonal_minimum.begin(), diagonal_minimum.end(),
                   [](double value) { return std::isfinite(value) && value > 0.0; }))
  {
    throw std::runtime_error(
        "Singular stiffness diagonal diagnostics have inconsistent global sizes!");
  }

  return {fem::singular::ReferenceIntegral::ConventionVersion,
          standard_order,
          singular_order,
          options.quadrature_order,
          options.fixed_subdivision,
          options.subdivisions,
          options.absolute_tolerance,
          options.relative_tolerance,
          options.maximum_subdivisions,
          quadrature_leaf_count,
          quadrature_maximum_depth,
          duffy_table_entries,
          duffy_cache_hits,
          numbering.h1.global_size,
          numbering.nd.global_size,
          diagonal_minimum[0],
          diagonal_maximum[0],
          diagonal_maximum[0] / diagonal_minimum[0],
          diagonal_minimum[1],
          diagonal_maximum[1],
          diagonal_maximum[1] / diagonal_minimum[1],
          std::max(diagonal_maximum[0], diagonal_maximum[1]) /
              std::min(diagonal_minimum[0], diagonal_minimum[1])};
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
  if (HasSingularEnrichment())
  {
    std::vector<fem::singular::IsotropicMaterialCoefficients> materials(GetMesh().GetNE(),
                                                                        {1.0, 1.0});
    fem::singular::LocalSparseH1EnrichmentMatrices local_enrichment;
    if (singular_features)
    {
      singular_dofs =
          std::make_unique<fem::singular::DofTopology>(fem::singular::BuildLocalDofTopology(
              GetMesh(), *singular_features, *source_vertex_ids, singular_order));
      singular_numbering = std::make_unique<fem::singular::ParallelDofNumbering>(
          fem::singular::BuildParallelDofNumbering(GetComm(), *singular_dofs));
      for (int element = 0; element < GetMesh().GetNE(); element++)
      {
        const auto &incidence = singular_features->elements[element];
        if (incidence.nodes.empty() && incidence.edges.empty())
        {
          continue;
        }
        const int attribute = GetMesh().Get().GetAttribute(element);
        MFEM_VERIFY(mat_op.IsPermittivityIsotropic(attribute),
                    "Singular electrostatic enrichment initially requires isotropic "
                    "permittivity in every enriched tetrahedron! Domain attribute: "
                        << attribute);
        materials[element].electric = mat_op.GetPermittivityReal(attribute)(0, 0);
      }
      local_enrichment = fem::singular::AssembleLocalSparseH1EnrichmentMatrices(
          *singular_dofs, GetH1Space().Get(), materials, singular_assembly_options);
    }
    else
    {
      triangle_singular_dofs = std::make_unique<fem::singular::TriangleDofTopology>(
          fem::singular::BuildLocalTriangleDofTopology(
              GetMesh(), *triangle_singular_features, *source_vertex_ids, singular_order));
      singular_numbering = std::make_unique<fem::singular::ParallelDofNumbering>(
          fem::singular::BuildParallelDofNumbering(GetComm(), *triangle_singular_dofs));
      for (int element = 0; element < GetMesh().GetNE(); element++)
      {
        if (triangle_singular_features->elements[element].nodes.empty())
        {
          continue;
        }
        const int attribute = GetMesh().Get().GetAttribute(element);
        MFEM_VERIFY(mat_op.IsPermittivityIsotropic(attribute),
                    "Triangular singular electrostatic enrichment requires isotropic "
                    "permittivity in every enriched triangle! Domain attribute: "
                        << attribute);
        materials[element].electric = mat_op.GetPermittivityReal(attribute)(0, 0);
      }
      local_enrichment = fem::singular::AssembleLocalSparseH1EnrichmentMatrices(
          *triangle_singular_dofs, GetH1Space().Get(), materials,
          singular_assembly_options);
    }
    auto parallel_enrichment = fem::singular::AssembleParallelSparseH1EnrichmentMatrices(
        local_enrichment, *singular_numbering, GetH1Space().Get());

    const std::size_t number_levels = GetH1Spaces().GetNumLevels();
    std::vector<std::unique_ptr<mfem::HypreParMatrix>> standard_operators(number_levels);
    for (std::size_t level = 0; level < number_levels; level++)
    {
      auto &h1_fespace = GetH1Spaces().GetFESpaceAtLevel(level);
      BilinearForm k_level(h1_fespace);
      k_level.AddDomainIntegrator<DiffusionIntegrator>(epsilon_func);
      standard_operators[level] =
          ParOperator(k_level.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();
    }
    const auto &standard_stiffness = standard_operators.back();
    auto unconstrained = fem::singular::BuildParallelEnrichedOperator(*standard_stiffness,
                                                                      parallel_enrichment);
    singular_unconstrained_stiffness =
        std::make_unique<mfem::HypreParMatrix>(*unconstrained);
    singular_stiffness_enrichment_error =
        std::move(parallel_enrichment.enrichment_enrichment_estimated_absolute_error);
    singular_stiffness_coupling_error =
        std::move(parallel_enrichment.standard_enrichment_estimated_absolute_error);
    singular_diagnostics =
        std::make_unique<SingularOperatorDiagnostics>(BuildSingularOperatorDiagnostics(
            *singular_unconstrained_stiffness, GetH1Space().GetTrueVSize(),
            *singular_numbering, singular_assembly_options, standard_order, singular_order,
            local_enrichment.total_quadrature_leaf_count,
            local_enrichment.maximum_subdivision_depth,
            local_enrichment.duffy_reference_table_entries,
            local_enrichment.duffy_reference_cache_hits));

    auto enrichment_gradient =
        fem::singular::BuildParallelEnrichmentGradient(GetComm(), *singular_numbering);
    const auto &standard_gradient_operator =
        GetNDSpace().GetDiscreteInterpolator(GetH1Space());
    const auto *standard_gradient =
        dynamic_cast<const ParOperator *>(&standard_gradient_operator);
    MFEM_VERIFY(standard_gradient,
                "Singular electrostatics requires an assembled standard gradient!");
    singular_gradient = fem::singular::BuildParallelEnrichedGradient(
        standard_gradient->ParallelAssemble(), *enrichment_gradient);

    mfem::Array<int> enrichment_essential;
    std::vector<std::vector<std::size_t>> feature_membership;
    std::size_t feature_count;
    if (singular_features)
    {
      enrichment_essential = fem::singular::GetEssentialH1TrueDofs(
          GetComm(), *singular_features, *singular_dofs, *singular_numbering);
      feature_membership =
          fem::singular::BuildH1DofFeatureMembership(*singular_features, *singular_dofs);
      feature_count = singular_features->features.size();
    }
    else
    {
      enrichment_essential = fem::singular::GetEssentialTriangleH1TrueDofs(
          GetComm(), *triangle_singular_features, *triangle_singular_dofs,
          *singular_numbering);
      feature_membership = fem::singular::BuildTriangleH1DofFeatureMembership(
          *triangle_singular_features, *triangle_singular_dofs);
      feature_count = triangle_singular_features->vertices.size();
    }

    combined_h1_dbc_tdof_lists = dbc_tdof_lists;
    for (std::size_t level = 0; level < number_levels; level++)
    {
      auto &combined_essential = combined_h1_dbc_tdof_lists[level];
      const int standard_size = GetH1Spaces().GetFESpaceAtLevel(level).GetTrueVSize();
      for (int dof : enrichment_essential)
      {
        combined_essential.Append(standard_size + dof);
      }
      combined_essential.Sort();
      MFEM_VERIFY(std::adjacent_find(combined_essential.begin(),
                                     combined_essential.end()) == combined_essential.end(),
                  "Singular electrostatic essential true DOFs are not unique!");
    }

    singular_h1_prolongations.reserve(number_levels > 0 ? number_levels - 1 : 0);
    for (std::size_t level = 0; level + 1 < number_levels; level++)
    {
      const auto *standard_prolongation =
          dynamic_cast<const ParOperator *>(&GetH1Spaces().GetProlongationAtLevel(level));
      MFEM_VERIFY(standard_prolongation,
                  "Singular electrostatic p-multigrid requires assembled standard H1 "
                  "prolongation operators!");
      singular_h1_prolongations.push_back(fem::singular::BuildParallelEnrichedProlongation(
          standard_prolongation->ParallelAssemble(), singular_numbering->h1));
    }

    auto constrained_blocks = fem::singular::BuildParallelConstrainedOperatorBlocks(
        *standard_stiffness, parallel_enrichment, dbc_tdof_lists.back(),
        enrichment_essential);
    auto dirichlet = fem::singular::BuildParallelDirichletSystem(
        std::move(unconstrained), GetH1Space().GetTrueVSize(), dbc_tdof_lists.back(),
        enrichment_essential);

    std::vector<std::unique_ptr<mfem::HypreParMatrix>> combined_operators(number_levels),
        coupling_operators(number_levels);
    const std::size_t finest_level = number_levels - 1;
    combined_operators[finest_level] = std::move(dirichlet.constrained);
    standard_operators[finest_level] = std::move(constrained_blocks.standard_standard);
    coupling_operators[finest_level] = std::move(constrained_blocks.standard_enrichment);
    // Rediscretize the standard block to preserve local FE sparsity. The singular space is
    // identical on every p-level, so A_ee is invariant and A_se,c = P_s^T A_se,f.
    for (std::size_t fine_level = finest_level; fine_level > 0; fine_level--)
    {
      const std::size_t coarse_level = fine_level - 1;
      std::unique_ptr<mfem::HypreParMatrix> discarded(
          standard_operators[coarse_level]->EliminateRowsCols(
              dbc_tdof_lists[coarse_level]));

      const auto *standard_par_operator = dynamic_cast<const ParOperator *>(
          &GetH1Spaces().GetProlongationAtLevel(coarse_level));
      MFEM_VERIFY(standard_par_operator,
                  "Singular electrostatic p-multigrid requires assembled standard H1 "
                  "prolongation operators!");
      const auto &standard_prolongation = standard_par_operator->ParallelAssemble();
      std::unique_ptr<mfem::HypreParMatrix> standard_restriction(
          standard_prolongation.Transpose());
      coupling_operators[coarse_level].reset(
          mfem::ParMult(standard_restriction.get(), coupling_operators[fine_level].get()));
      MFEM_VERIFY(coupling_operators[coarse_level],
                  "Failed to restrict the singular electrostatic coupling block!");
      coupling_operators[coarse_level]->EliminateRows(dbc_tdof_lists[coarse_level]);
      discarded.reset(
          coupling_operators[coarse_level]->EliminateCols(enrichment_essential));

      fem::singular::ParallelSparseOperatorBlocks coarse_enrichment;
      coarse_enrichment.enrichment_enrichment =
          std::make_unique<mfem::HypreParMatrix>(*constrained_blocks.enrichment_enrichment);
      coarse_enrichment.standard_enrichment = std::move(coupling_operators[coarse_level]);
      coarse_enrichment.enrichment_standard.reset(
          coarse_enrichment.standard_enrichment->Transpose());
      combined_operators[coarse_level] = fem::singular::BuildParallelEnrichedOperator(
          *standard_operators[coarse_level], coarse_enrichment);
      coupling_operators[coarse_level] = std::move(coarse_enrichment.standard_enrichment);
    }

    auto feature_patches = fem::singular::BuildParallelFeaturePatches(
        *combined_operators.front(), *coupling_operators.front(),
        GetH1Spaces().GetFESpaceAtLevel(0).GetTrueVSize(), feature_membership,
        singular_numbering->h1, feature_count);
    singular_constrained_standard_stiffness = std::move(standard_operators.front());
    singular_diagnostics->feature_patch_count = feature_patches.patches.size();
    singular_diagnostics->feature_patch_sum_standard_dofs =
        feature_patches.sum_global_standard_dofs;
    singular_diagnostics->feature_patch_sum_enrichment_dofs =
        feature_patches.sum_global_enrichment_dofs;
    singular_diagnostics->feature_patch_maximum_standard_dofs =
        feature_patches.maximum_global_standard_dofs;
    singular_diagnostics->feature_patch_maximum_enrichment_dofs =
        feature_patches.maximum_global_enrichment_dofs;
    singular_diagnostics->feature_patch_minimum_enrichment_multiplicity =
        feature_patches.minimum_enrichment_multiplicity;
    singular_diagnostics->feature_patch_maximum_enrichment_multiplicity =
        feature_patches.maximum_enrichment_multiplicity;
    singular_feature_patches = std::move(feature_patches);
    singular_eliminated_stiffness = std::move(dirichlet.eliminated);
    singular_essential_true_dofs = std::move(dirichlet.essential_true_dofs);
    MFEM_VERIFY(singular_essential_true_dofs.Size() ==
                        combined_h1_dbc_tdof_lists.back().Size() &&
                    std::equal(singular_essential_true_dofs.begin(),
                               singular_essential_true_dofs.end(),
                               combined_h1_dbc_tdof_lists.back().begin()),
                "Finest singular electrostatic essential true DOFs are inconsistent!");

    if (triangle_singular_features)
    {
      Mpi::Print(" Singular H1 enrichment: {:d} global true DOFs, {:d} graded Duffy "
                 "quadrature point evaluations\n",
                 singular_numbering->h1.global_size,
                 singular_diagnostics->quadrature_leaf_count);
    }
    else if (singular_diagnostics->quadrature_fixed_subdivision)
    {
      Mpi::Print(" Singular H1 enrichment: {:d} global true DOFs, {:d} fixed "
                 "quadrature leaves (uniform depth = {:d})\n",
                 singular_numbering->h1.global_size,
                 singular_diagnostics->quadrature_leaf_count,
                 singular_diagnostics->quadrature_subdivisions);
    }
    else
    {
      Mpi::Print(" Singular H1 enrichment: {:d} global true DOFs, {:d} adaptive "
                 "quadrature leaves (maximum depth = {:d}), {:d} cached Duffy tables "
                 "per rank maximum ({:d} cache hits)\n",
                 singular_numbering->h1.global_size,
                 singular_diagnostics->quadrature_leaf_count,
                 singular_diagnostics->quadrature_maximum_depth,
                 singular_diagnostics->duffy_reference_table_maximum_entries,
                 singular_diagnostics->duffy_reference_cache_hits);
    }
    Mpi::Print(" Singular stiffness diagonal spread: standard = {:.3e}, enrichment = "
               "{:.3e}, combined = {:.3e}\n",
               singular_diagnostics->standard_diagonal_spread,
               singular_diagnostics->enrichment_diagonal_spread,
               singular_diagnostics->combined_diagonal_spread);
    Mpi::Print(
        " Singular feature patches: {:d} patches, maximum {:d} global standard + "
        "{:d} global enrichment true DOFs, enrichment overlap multiplicity {:d}-{:d}\n",
        singular_diagnostics->feature_patch_count,
        singular_diagnostics->feature_patch_maximum_standard_dofs,
        singular_diagnostics->feature_patch_maximum_enrichment_dofs,
        singular_diagnostics->feature_patch_minimum_enrichment_multiplicity,
        singular_diagnostics->feature_patch_maximum_enrichment_multiplicity);
    print_hdr = false;
    if (number_levels == 1)
    {
      return std::move(combined_operators.front());
    }

    auto hierarchy = std::make_unique<MultigridOperator>(number_levels);
    for (std::size_t level = 0; level < number_levels; level++)
    {
      // mfem::HypreParMatrix::NNZ() reports the global nonzero count. Unlike the local
      // HypreCSRMatrix count used below, it must not be summed again across ranks.
      const HYPRE_BigInt nnz = combined_operators[level]->NNZ();
      Mpi::Print(" Level {:d} (p = {:d}): {:d} combined H1 unknowns, {:d} NNZ\n", level,
                 GetH1Spaces().GetFESpaceAtLevel(level).GetMaxElementOrder(),
                 combined_operators[level]->GetGlobalNumRows(), nnz);
      hierarchy->AddOperator(std::move(combined_operators[level]));
    }
    return hierarchy;
  }

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

const Operator &LaplaceOperator::GetGradMatrix() const
{
  if (HasSingularEnrichment())
  {
    MFEM_VERIFY(singular_gradient,
                "Singular gradient requested before stiffness assembly!");
    return *singular_gradient;
  }
  return GetNDSpace().GetDiscreteInterpolator(GetH1Space());
}

const mfem::HypreParMatrix &LaplaceOperator::GetUnconstrainedStiffnessMatrix() const
{
  MFEM_VERIFY(singular_unconstrained_stiffness,
              "Unconstrained stiffness is available only after singular assembly!");
  return *singular_unconstrained_stiffness;
}

const mfem::HypreParMatrix &LaplaceOperator::GetSingularStandardStiffnessMatrix() const
{
  MFEM_VERIFY(singular_constrained_standard_stiffness,
              "Constrained standard stiffness is available only after singular "
              "assembly!");
  return *singular_constrained_standard_stiffness;
}

const fem::singular::ParallelFeaturePatches &
LaplaceOperator::GetSingularFeaturePatches() const
{
  MFEM_VERIFY(!singular_feature_patches.patches.empty(),
              "Feature patches are available only after singular assembly!");
  return singular_feature_patches;
}

std::vector<const Operator *> LaplaceOperator::GetCombinedH1ProlongationOperators() const
{
  std::vector<const Operator *> operators;
  operators.reserve(singular_h1_prolongations.size());
  for (const auto &prolongation : singular_h1_prolongations)
  {
    operators.push_back(prolongation.get());
  }
  return operators;
}

const SingularOperatorDiagnostics &LaplaceOperator::GetSingularDiagnostics() const
{
  MFEM_VERIFY(singular_diagnostics,
              "Singular diagnostics are available only after singular assembly!");
  return *singular_diagnostics;
}

double LaplaceOperator::GetSingularStiffnessEnergyErrorBound(
    const mfem::Vector &combined_true_dofs) const
{
  const int standard_size = GetH1Space().GetTrueVSize();
  const int enrichment_size =
      singular_numbering ? static_cast<int>(singular_numbering->h1.owned_size) : 0;
  bool valid = singular_stiffness_enrichment_error && singular_stiffness_coupling_error &&
               combined_true_dofs.Size() == standard_size + enrichment_size &&
               singular_stiffness_enrichment_error->Height() == enrichment_size &&
               singular_stiffness_enrichment_error->Width() == enrichment_size &&
               singular_stiffness_coupling_error->Height() == standard_size &&
               singular_stiffness_coupling_error->Width() == enrichment_size;
  Mpi::GlobalAnd(1, &valid, GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Singular stiffness error bound received inconsistent assembly data!");
  }

  Vector standard(standard_size), enrichment(enrichment_size);
  for (int i = 0; i < standard_size; i++)
  {
    standard[i] = std::abs(combined_true_dofs[i]);
  }
  for (int i = 0; i < enrichment_size; i++)
  {
    enrichment[i] = std::abs(combined_true_dofs[standard_size + i]);
  }

  Vector enrichment_error(enrichment_size), coupling_error(standard_size);
  singular_stiffness_enrichment_error->Mult(enrichment, enrichment_error);
  singular_stiffness_coupling_error->Mult(enrichment, coupling_error);
  const double bound = linalg::Dot(GetComm(), enrichment, enrichment_error) +
                       2.0 * linalg::Dot(GetComm(), standard, coupling_error);
  if (!std::isfinite(bound) || bound < 0.0)
  {
    throw std::runtime_error(
        "Singular stiffness energy error bound is nonfinite or negative!");
  }
  return bound;
}

std::unique_ptr<fem::singular::EnrichedH1FieldEvaluator>
LaplaceOperator::GetSingularFieldEvaluator()
{
  MFEM_VERIFY(singular_features && singular_dofs && singular_numbering,
              "Tetrahedral singular field evaluation is available only after singular "
              "assembly!");
  return std::make_unique<fem::singular::EnrichedH1FieldEvaluator>(
      *singular_dofs, *singular_numbering, GetH1Space().Get());
}

std::unique_ptr<fem::singular::TriangleEnrichedH1FieldEvaluator>
LaplaceOperator::GetTriangleSingularFieldEvaluator()
{
  MFEM_VERIFY(triangle_singular_features && triangle_singular_dofs && singular_numbering,
              "Triangular singular field evaluation is available only after singular "
              "assembly!");
  return std::make_unique<fem::singular::TriangleEnrichedH1FieldEvaluator>(
      *triangle_singular_dofs, *singular_numbering, GetH1Space().Get());
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
  X.SetSize(K.Width());
  RHS.SetSize(K.Height());
  X.UseDevice(true);
  RHS.UseDevice(true);
  X = 0.0;
  RHS = 0.0;
  if (HasSingularEnrichment())
  {
    HYPRE_BigInt global_essential_true_dofs = singular_essential_true_dofs.Size();
    Mpi::GlobalSum(1, &global_essential_true_dofs, mesh.GetComm());
    MFEM_VERIFY(singular_eliminated_stiffness && global_essential_true_dofs > 0 &&
                    K.Width() >= GetH1Space().GetTrueVSize(),
                "Singular Dirichlet system was not assembled!");
    Vector standard_x(GetH1Space().GetTrueVSize());
    x.ParallelProject(standard_x);
    for (int i = 0; i < standard_x.Size(); i++)
    {
      X[i] = standard_x[i];
    }
    const auto *hierarchy = dynamic_cast<const MultigridOperator *>(&K);
    const auto *hypre_K = dynamic_cast<const mfem::HypreParMatrix *>(
        hierarchy ? &hierarchy->GetFinestOperator() : &K);
    MFEM_VERIFY(hypre_K, "Singular LaplaceOperator requires an assembled Hypre matrix!");
    hypre_K->EliminateBC(*singular_eliminated_stiffness, singular_essential_true_dofs, X,
                         RHS);
    return;
  }

  x.ParallelProject(X);  // Restrict to the true dofs
  const auto *mg_K = dynamic_cast<const MultigridOperator *>(&K);
  const auto *PtAP_K = mg_K ? dynamic_cast<const ParOperator *>(&mg_K->GetFinestOperator())
                            : dynamic_cast<const ParOperator *>(&K);
  MFEM_VERIFY(PtAP_K, "LaplaceOperator requires ParOperator for RHS elimination!");
  PtAP_K->EliminateRHS(X, RHS);
}

}  // namespace palace
