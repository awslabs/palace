// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "spaceoperator.hpp"

#include <limits>
#include <set>
#include <type_traits>
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

namespace palace
{

using namespace std::complex_literals;

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

  if (tetrahedral)
  {
    singular_dofs =
        std::make_unique<fem::singular::DofTopology>(fem::singular::BuildLocalDofTopology(
            GetMesh(), *singular_features, *source_vertex_ids,
            solver.singular_elements.order));
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
  MFEM_VERIFY(singular_numbering->h1.owned_size <= std::numeric_limits<int>::max() &&
                  singular_numbering->nd.owned_size <= std::numeric_limits<int>::max(),
              "Full-wave singular local true-DOF count exceeds integer limits!");

  std::vector<fem::singular::IsotropicMaterialCoefficients> materials(GetMesh().GetNE(),
                                                                      {1.0, 1.0});
  std::vector<fem::singular::IsotropicMaterialCoefficients> imaginary_materials(
      GetMesh().GetNE(), {0.0, 1.0});
  std::vector<fem::singular::IsotropicMaterialCoefficients> absolute_materials(
      GetMesh().GetNE(), {1.0, 1.0});
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
  for (std::size_t level = 0; level < number_levels; level++)
  {
    auto &h1_space = GetH1Spaces().GetFESpaceAtLevel(level);
    auto &nd_space = GetNDSpaces().GetFESpaceAtLevel(level);
    std::vector<std::vector<fem::singular::IsotropicMaterialCoefficients>> material_batches{
        materials, absolute_materials};
    if (mat_op.HasLossTangent())
    {
      material_batches.push_back(imaginary_materials);
    }
    const auto local_matrices =
        tetrahedral
            ? fem::singular::AssembleLocalSparseEnrichmentMatricesBatch(
                  *singular_dofs, h1_space.Get(), nd_space.Get(), material_batches, options)
            : fem::singular::AssembleLocalSparseEnrichmentMatricesBatch(
                  *triangle_singular_dofs, h1_space.Get(), nd_space.Get(), material_batches,
                  options);
    MFEM_VERIFY(local_matrices.size() == material_batches.size(),
                "Full-wave singular material batch assembly returned an inconsistent "
                "number of operators!");
    const auto assemble_domain = [&](std::size_t batch)
    {
      return fem::singular::AssembleParallelSparseEnrichmentMatrices(
          local_matrices[batch], *singular_numbering, h1_space.Get(), nd_space.Get());
    };
    singular_domain_matrices[level] = assemble_domain(0);
    singular_domain_abs_matrices[level] = assemble_domain(1);
    if (mat_op.HasLossTangent())
    {
      singular_domain_imag_matrices[level] = assemble_domain(2);
    }

    const auto assemble_lumped_boundary = [&](const std::map<int, double> &coefficients)
    {
      fem::singular::ParallelSparseOperatorBlocks result;
      if (coefficients.empty())
      {
        return result;
      }
      const auto local_boundary =
          tetrahedral ? fem::singular::AssembleLocalSparseNDBoundaryMassMatrices(
                            *singular_dofs, nd_space.Get(), coefficients, options)
                      : fem::singular::AssembleLocalSparseNDBoundaryMassMatrices(
                            *triangle_singular_dofs, nd_space.Get(), coefficients, options);
      return fem::singular::AssembleParallelSparseNDBoundaryMassMatrices(
          local_boundary, *singular_numbering, nd_space.Get());
    };
    singular_lumped_stiffness_matrices[level] =
        assemble_lumped_boundary(lumped_port_op.GetStiffnessBdrCoefficientMap());
    singular_lumped_damping_matrices[level] =
        assemble_lumped_boundary(lumped_port_op.GetDampingBdrCoefficientMap());
    singular_lumped_mass_matrices[level] =
        assemble_lumped_boundary(lumped_port_op.GetMassBdrCoefficientMap());
    singular_impedance_stiffness_matrices[level] =
        assemble_lumped_boundary(impedance_stiffness_coefficients);
    singular_impedance_damping_matrices[level] =
        assemble_lumped_boundary(impedance_damping_coefficients);
    singular_impedance_mass_matrices[level] =
        assemble_lumped_boundary(impedance_mass_coefficients);

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

  combined_nd_dbc_tdof_lists = nd_dbc_tdof_lists;
  combined_h1_dbc_tdof_lists = h1_dbc_tdof_lists;
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
    MFEM_VERIFY(
        std::adjacent_find(combined_nd.begin(), combined_nd.end()) == combined_nd.end() &&
            std::adjacent_find(combined_h1.begin(), combined_h1.end()) == combined_h1.end(),
        "Full-wave singular essential true DOFs are not unique!");
  }

  singular_nd_prolongations.reserve(number_levels > 0 ? number_levels - 1 : 0);
  singular_h1_prolongations.reserve(number_levels > 0 ? number_levels - 1 : 0);
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
    singular_nd_prolongations.push_back(fem::singular::BuildParallelEnrichedProlongation(
        standard_nd_prolongation->ParallelAssemble(), singular_numbering->nd));
    singular_h1_prolongations.push_back(fem::singular::BuildParallelEnrichedProlongation(
        standard_h1_prolongation->ParallelAssemble(), singular_numbering->h1));
  }

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
  return result;
}

fem::singular::ParallelSparseOperatorBlocks AddScaledSingularOperatorBlocks(
    std::initializer_list<
        std::pair<double, const fem::singular::ParallelSparseOperatorBlocks *>>
        terms)
{
  fem::singular::ParallelSparseOperatorBlocks result;
  const fem::singular::ParallelSparseOperatorBlocks *zero_template = nullptr;
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
    auto standard =
        ParOperator(std::move(kr), GetNDSpace()).StealParallelAssemble(skip_zeros);
    const auto enrichment = AddScaledSingularOperatorBlocks(
        {{1.0, &singular_domain_matrices.back().nd_curl_curl},
         {1.0, &singular_lumped_stiffness_matrices.back()},
         {1.0, &singular_impedance_stiffness_matrices.back()}});
    auto combined = fem::singular::BuildParallelEnrichedOperator(*standard, enrichment);
    combined->EliminateBC(GetCombinedNDDbcTDofList(), diag_policy);
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
    auto standard =
        ParOperator(std::move(c), GetNDSpace()).StealParallelAssemble(skip_zeros);
    const auto enrichment = AddScaledSingularOperatorBlocks(
        {{1.0, &singular_lumped_damping_matrices.back()},
         {1.0, &singular_impedance_damping_matrices.back()}});
    auto combined = fem::singular::BuildParallelEnrichedOperator(*standard, enrichment);
    combined->EliminateBC(GetCombinedNDDbcTDofList(), diag_policy);
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
    auto standard_real =
        ParOperator(std::move(mr), GetNDSpace()).StealParallelAssemble(skip_zeros);
    const auto real_enrichment =
        AddScaledSingularOperatorBlocks({{1.0, &singular_domain_matrices.back().nd_mass},
                                         {1.0, &singular_lumped_mass_matrices.back()},
                                         {1.0, &singular_impedance_mass_matrices.back()}});
    auto combined_real =
        fem::singular::BuildParallelEnrichedOperator(*standard_real, real_enrichment);
    combined_real->EliminateBC(GetCombinedNDDbcTDofList(), diag_policy);
    if constexpr (std::is_same<OperType, ComplexOperator>::value)
    {
      std::unique_ptr<mfem::HypreParMatrix> combined_imaginary;
      if (mi)
      {
        MFEM_VERIFY(singular_domain_imag_matrices.back().nd_mass.enrichment_enrichment,
                    "Full-wave singular imaginary mass blocks were not assembled!");
        auto standard_imaginary =
            ParOperator(std::move(mi), GetNDSpace()).StealParallelAssemble(skip_zeros);
        combined_imaginary = fem::singular::BuildParallelEnrichedOperator(
            *standard_imaginary, singular_domain_imag_matrices.back().nd_mass);
        combined_imaginary->EliminateBC(GetCombinedNDDbcTDofList(),
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
  auto standard =
      ParOperator(std::move(kr), GetNDSpace()).StealParallelAssemble(skip_zeros);
  auto combined = fem::singular::BuildParallelEnrichedOperator(
      *standard, singular_domain_matrices.back().nd_curl_curl);
  combined->EliminateBC(GetCombinedNDDbcTDofList(), diag_policy);
  return combined;
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
  auto standard =
      ParOperator(std::move(mr), GetNDSpace()).StealParallelAssemble(skip_zeros);
  auto combined = fem::singular::BuildParallelEnrichedOperator(
      *standard, singular_domain_matrices.back().nd_mass);
  combined->EliminateBC(GetCombinedNDDbcTDofList(), diag_policy);
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

namespace
{

std::unique_ptr<mfem::HypreParMatrix> AddScaledHypreMatrices(
    const std::vector<std::pair<double, const mfem::HypreParMatrix *>> &terms)
{
  std::unique_ptr<mfem::HypreParMatrix> result;
  for (const auto &[coefficient, matrix] : terms)
  {
    if (coefficient == 0.0 || !matrix)
    {
      continue;
    }
    if (!result)
    {
      result = std::make_unique<mfem::HypreParMatrix>(*matrix);
      *result *= coefficient;
    }
    else
    {
      result.reset(mfem::Add(1.0, *result, coefficient, *matrix));
    }
  }
  return result;
}

const mfem::HypreParMatrix *GetHyprePart(const ComplexOperator *op, bool imaginary)
{
  if (!op)
  {
    return nullptr;
  }
  return dynamic_cast<const mfem::HypreParMatrix *>(imaginary ? op->Imag() : op->Real());
}

}  // namespace

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
      const auto sparse = [](const ComplexOperator *op)
      {
        return !op ||
               ((!op->Real() || dynamic_cast<const mfem::HypreParMatrix *>(op->Real())) &&
                (!op->Imag() || dynamic_cast<const mfem::HypreParMatrix *>(op->Imag())));
      };
      MFEM_VERIFY(sparse(K) && sparse(C) && sparse(M) && sparse(A2),
                  "Full-wave singular system composition requires sparse assembled "
                  "real and imaginary operator parts!");
      const std::array<std::pair<std::complex<double>, const ComplexOperator *>, 4> terms{
          std::pair{a0, K}, std::pair{a1, C}, std::pair{a2, M},
          std::pair{std::complex<double>(1.0, 0.0), A2}};
      std::vector<std::pair<double, const mfem::HypreParMatrix *>> real_terms;
      std::vector<std::pair<double, const mfem::HypreParMatrix *>> imaginary_terms;
      real_terms.reserve(2 * terms.size());
      imaginary_terms.reserve(2 * terms.size());
      for (const auto &[coefficient, op] : terms)
      {
        const auto *real = GetHyprePart(op, false);
        const auto *imaginary = GetHyprePart(op, true);
        real_terms.emplace_back(coefficient.real(), real);
        real_terms.emplace_back(-coefficient.imag(), imaginary);
        imaginary_terms.emplace_back(coefficient.real(), imaginary);
        imaginary_terms.emplace_back(coefficient.imag(), real);
      }
      auto Ar = AddScaledHypreMatrices(real_terms);
      auto Ai = AddScaledHypreMatrices(imaginary_terms);
      MFEM_VERIFY(Ar || Ai, "Full-wave singular system matrix is empty!");
      return std::make_unique<ComplexWrapperOperator>(std::move(Ar), std::move(Ai));
    }
    else
    {
      const auto *Kr = dynamic_cast<const mfem::HypreParMatrix *>(K);
      const auto *Cr = dynamic_cast<const mfem::HypreParMatrix *>(C);
      const auto *Mr = dynamic_cast<const mfem::HypreParMatrix *>(M);
      const auto *A2r = dynamic_cast<const mfem::HypreParMatrix *>(A2);
      auto A = AddScaledHypreMatrices({{a0, Kr}, {a1, Cr}, {a2, Mr}, {1.0, A2r}});
      MFEM_VERIFY(A, "Full-wave singular real system matrix is empty!");
      return A;
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
    const auto *Kr = GetHyprePart(K, false);
    const auto *Mr = GetHyprePart(M, false);
    MFEM_VERIFY((!K || Kr) && (!M || Mr),
                "Full-wave singular inner product requires sparse real matrix parts!");
    return AddScaledHypreMatrices({{a0, Kr}, {a2, Mr}});
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
        standard_imaginary_operators(number_levels), unused_auxiliary(number_levels),
        unused_imaginary_auxiliary(number_levels);
    if constexpr (std::is_same_v<OperType, ComplexOperator>)
    {
      if (complex_preconditioner)
      {
        AssemblePreconditioner(a0, a1, a2, a3, standard_operators, unused_auxiliary,
                               standard_imaginary_operators, unused_imaginary_auxiliary);
      }
      else
      {
        AssemblePreconditioner(a0, a1, a2, a3, standard_operators, unused_auxiliary);
      }
    }
    else
    {
      AssemblePreconditioner(a0, a1, a2, a3, standard_operators, unused_auxiliary);
    }
    MFEM_VERIFY(standard_operators.back(),
                "Full-wave singular preconditioner finest standard level is empty!");
    std::vector<std::unique_ptr<mfem::HypreParMatrix>> combined_operators(number_levels),
        combined_imaginary_operators(number_levels);

    const std::size_t finest_level = number_levels - 1;
    auto &finest_nd_space = GetNDSpaces().GetFESpaceAtLevel(finest_level);
    auto finest_standard =
        ParOperator(std::move(standard_operators[finest_level]), finest_nd_space)
            .StealParallelAssemble(false);
    const auto finest_enrichment =
        complex_preconditioner
            ? AddScaledSingularOperatorBlocks(
                  {{stiffness_coefficient,
                    &singular_domain_matrices[finest_level].nd_curl_curl},
                   {stiffness_coefficient,
                    &singular_lumped_stiffness_matrices[finest_level]},
                   {stiffness_coefficient,
                    &singular_impedance_stiffness_matrices[finest_level]},
                   {damping_coefficient, &singular_lumped_damping_matrices[finest_level]},
                   {damping_coefficient,
                    &singular_impedance_damping_matrices[finest_level]},
                   {mass_coefficient, &singular_domain_matrices[finest_level].nd_mass},
                   {-imaginary_part(a2),
                    &singular_domain_imag_matrices[finest_level].nd_mass},
                   {mass_coefficient, &singular_lumped_mass_matrices[finest_level]},
                   {mass_coefficient, &singular_impedance_mass_matrices[finest_level]}})
            : AddScaledSingularOperatorBlocks(
                  {{stiffness_coefficient,
                    &singular_domain_matrices[finest_level].nd_curl_curl},
                   {stiffness_coefficient,
                    &singular_lumped_stiffness_matrices[finest_level]},
                   {stiffness_coefficient,
                    &singular_impedance_stiffness_matrices[finest_level]},
                   {damping_coefficient, &singular_lumped_damping_matrices[finest_level]},
                   {damping_coefficient,
                    &singular_impedance_damping_matrices[finest_level]},
                   {mass_coefficient, &singular_domain_abs_matrices[finest_level].nd_mass},
                   {mass_coefficient, &singular_lumped_mass_matrices[finest_level]},
                   {mass_coefficient, &singular_impedance_mass_matrices[finest_level]}});
    combined_operators[finest_level] =
        fem::singular::BuildParallelEnrichedOperator(*finest_standard, finest_enrichment);
    combined_operators[finest_level]->EliminateBC(combined_nd_dbc_tdof_lists[finest_level],
                                                  Operator::DIAG_ONE);
    if (complex_preconditioner && standard_imaginary_operators[finest_level])
    {
      auto finest_standard_imaginary =
          ParOperator(std::move(standard_imaginary_operators[finest_level]),
                      finest_nd_space)
              .StealParallelAssemble(false);
      const auto finest_imaginary_enrichment = AddScaledSingularOperatorBlocks(
          {{imaginary_part(a0), &singular_domain_matrices[finest_level].nd_curl_curl},
           {imaginary_part(a0), &singular_lumped_stiffness_matrices[finest_level]},
           {imaginary_part(a0), &singular_impedance_stiffness_matrices[finest_level]},
           {imaginary_part(a1), &singular_lumped_damping_matrices[finest_level]},
           {imaginary_part(a1), &singular_impedance_damping_matrices[finest_level]},
           {imaginary_part(a2), &singular_domain_matrices[finest_level].nd_mass},
           {real_part(a2), &singular_domain_imag_matrices[finest_level].nd_mass},
           {imaginary_part(a2), &singular_lumped_mass_matrices[finest_level]},
           {imaginary_part(a2), &singular_impedance_mass_matrices[finest_level]}});
      combined_imaginary_operators[finest_level] =
          fem::singular::BuildParallelEnrichedOperator(*finest_standard_imaginary,
                                                       finest_imaginary_enrichment);
      combined_imaginary_operators[finest_level]->EliminateBC(
          combined_nd_dbc_tdof_lists[finest_level], Operator::DIAG_ZERO);
    }
    for (std::size_t fine_level = finest_level; fine_level > 0; fine_level--)
    {
      const std::size_t coarse_level = fine_level - 1;
      auto *prolongation = singular_nd_prolongations[coarse_level].get();
      combined_operators[coarse_level].reset(
          mfem::RAP(prolongation, combined_operators[fine_level].get(), prolongation));
      MFEM_VERIFY(combined_operators[coarse_level],
                  "Failed to Galerkin-project the combined singular preconditioner!");
      combined_operators[coarse_level]->EliminateBC(
          combined_nd_dbc_tdof_lists[coarse_level], Operator::DIAG_ONE);
      if (combined_imaginary_operators[fine_level])
      {
        combined_imaginary_operators[coarse_level].reset(mfem::RAP(
            prolongation, combined_imaginary_operators[fine_level].get(), prolongation));
        MFEM_VERIFY(combined_imaginary_operators[coarse_level],
                    "Failed to Galerkin-project the imaginary singular preconditioner!");
        combined_imaginary_operators[coarse_level]->EliminateBC(
            combined_nd_dbc_tdof_lists[coarse_level], Operator::DIAG_ZERO);
      }
    }

    auto hierarchy = std::make_unique<BaseMultigridOperator<OperType>>(number_levels);
    for (std::size_t level = 0; level < number_levels; level++)
    {
      auto &nd_space = GetNDSpaces().GetFESpaceAtLevel(level);
      std::unique_ptr<mfem::HypreParMatrix> auxiliary(
          mfem::RAP(singular_gradients[level].get(), combined_operators[level].get(),
                    singular_gradients[level].get()));
      MFEM_VERIFY(auxiliary,
                  "Failed to project the combined singular preconditioner into H1!");
      auxiliary->EliminateBC(combined_h1_dbc_tdof_lists[level], Operator::DIAG_ONE);
      std::unique_ptr<mfem::HypreParMatrix> imaginary_auxiliary;
      if (combined_imaginary_operators[level])
      {
        imaginary_auxiliary.reset(mfem::RAP(singular_gradients[level].get(),
                                            combined_imaginary_operators[level].get(),
                                            singular_gradients[level].get()));
        MFEM_VERIFY(imaginary_auxiliary,
                    "Failed to project the imaginary singular preconditioner into H1!");
        imaginary_auxiliary->EliminateBC(combined_h1_dbc_tdof_lists[level],
                                         Operator::DIAG_ZERO);
      }
      if (print_prec_hdr)
      {
        HYPRE_BigInt primary_nnz = combined_operators[level]->NNZ();
        HYPRE_BigInt auxiliary_nnz = auxiliary->NNZ();
        Mpi::GlobalSum(1, &primary_nnz, GetComm());
        Mpi::GlobalSum(1, &auxiliary_nnz, GetComm());
        Mpi::Print(" Level {:d} (p = {:d}): {:d} combined ND unknowns, {:d} NNZ\n"
                   " Level {:d} (auxiliary): {:d} combined H1 unknowns, {:d} NNZ\n",
                   level, nd_space.GetMaxElementOrder(),
                   combined_operators[level]->GetGlobalNumRows(), primary_nnz, level,
                   auxiliary->GetGlobalNumRows(), auxiliary_nnz);
      }
      if constexpr (std::is_same_v<OperType, ComplexOperator>)
      {
        hierarchy->AddOperator(std::make_unique<ComplexWrapperOperator>(
            std::move(combined_operators[level]),
            std::move(combined_imaginary_operators[level])));
        hierarchy->AddAuxiliaryOperator(std::make_unique<ComplexWrapperOperator>(
            std::move(auxiliary), std::move(imaginary_auxiliary)));
      }
      else
      {
        hierarchy->AddOperator(std::move(combined_operators[level]));
        hierarchy->AddAuxiliaryOperator(std::move(auxiliary));
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

std::complex<double>
SpaceOperator::GetSingularBoundaryFunctional(SumVectorCoefficient &coefficient,
                                             const ComplexVector &field)
{
  MFEM_VERIFY(HasSingularEnrichment() && field.Size() == GetNDTrueVSize(),
              "Singular boundary functional requires a combined ND true-DOF field!");

  mfem::ParLinearForm standard_form(&GetNDSpace().Get());
  standard_form.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(coefficient));
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
        singular_assembly_options, enrichment_true_dofs);
  }
  else
  {
    MFEM_VERIFY(triangle_singular_dofs,
                "Singular boundary functional requires simplex ND topology!");
    fem::singular::AssembleParallelNDBoundaryLinearForm(
        *triangle_singular_dofs, *singular_numbering, GetNDSpace().Get(), coefficient,
        singular_assembly_options, enrichment_true_dofs);
  }

  std::complex<double> result = 0.0;
  for (int i = 0; i < standard_true_dofs.Size(); i++)
  {
    result +=
        standard_true_dofs[i] * std::complex<double>(field.Real()[i], field.Imag()[i]);
  }
  for (int i = 0; i < enrichment_true_dofs.Size(); i++)
  {
    const int combined_dof = standard_true_dofs.Size() + i;
    result += enrichment_true_dofs[i] *
              std::complex<double>(field.Real()[combined_dof], field.Imag()[combined_dof]);
  }
  Mpi::GlobalSum(1, &result, GetComm());
  return result;
}

std::complex<double> SpaceOperator::GetSingularLumpedPortVoltage(int port_idx,
                                                                 const ComplexVector &field)
{
  SumVectorCoefficient coefficient(GetMesh().SpaceDimension());
  lumped_port_op.GetPort(port_idx).AddVoltageBdrCoefficients(coefficient);
  return GetSingularBoundaryFunctional(coefficient, field);
}

std::complex<double>
SpaceOperator::GetSingularLumpedPortSParameter(int port_idx, const ComplexVector &field)
{
  SumVectorCoefficient coefficient(GetMesh().SpaceDimension());
  lumped_port_op.GetPort(port_idx).AddSParameterBdrCoefficients(coefficient);
  return GetSingularBoundaryFunctional(coefficient, field);
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
