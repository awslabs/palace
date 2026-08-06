// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hierarchicalmaxwellestimator.hpp"

#include <algorithm>
#include <cmath>
#include <map>

#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/singularassembly.hpp"
#include "fem/singularsystem.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"

namespace palace
{
namespace
{

int UnsignedDof(int dof)
{
  return dof >= 0 ? dof : -1 - dof;
}

double DofSign(int dof)
{
  return dof >= 0 ? 1.0 : -1.0;
}

std::map<int, double> AbsoluteBoundaryCoefficients(std::map<int, double> coefficients)
{
  for (auto &[attribute, coefficient] : coefficients)
  {
    (void)attribute;
    coefficient = std::abs(coefficient);
  }
  return coefficients;
}

void SetSignedMatrix(const mfem::DenseMatrix &input, const mfem::Array<int> &standard_dofs,
                     int local_size, mfem::DenseMatrix &output)
{
  const int standard_size = standard_dofs.Size();
  MFEM_VERIFY(input.Height() == local_size && input.Width() == local_size,
              "Hierarchical Maxwell element matrix has an inconsistent size!");
  output.SetSize(local_size);
  for (int row = 0; row < local_size; row++)
  {
    const double row_sign = row < standard_size ? DofSign(standard_dofs[row]) : 1.0;
    for (int column = 0; column < local_size; column++)
    {
      const double column_sign =
          column < standard_size ? DofSign(standard_dofs[column]) : 1.0;
      output(row, column) = row_sign * column_sign * input(row, column);
    }
  }
}

}  // namespace

HierarchicalMaxwellDomainData::HierarchicalMaxwellDomainData(SpaceOperator &space_op_)
  : space_op(&space_op_)
{
  auto &space_op = space_op_;
  MFEM_VERIFY(space_op.HasSingularEnrichment(),
              "Hierarchical Maxwell domain data requires singular enrichment!");
  MFEM_VERIFY(!space_op.GetMaterialOp().HasLondonDepth(),
              "The first hierarchical Maxwell domain adapter does not yet include the "
              "inverse-London-depth mass contribution to K!");
  const int dimension = space_op.GetMesh().Dimension();
  MFEM_VERIFY(dimension == 2 || dimension == 3,
              "Hierarchical Maxwell estimation requires a two- or three-dimensional mesh!");
  const int coarse_order = space_op.GetNDSpace().GetMaxElementOrder();
  fine_nd_collection = std::make_unique<mfem::ND_FECollection>(coarse_order + 1, dimension);
  fine_h1_collection = std::make_unique<mfem::H1_FECollection>(coarse_order + 1, dimension);
  fine_nd_space =
      std::make_unique<FiniteElementSpace>(space_op.GetMesh(), fine_nd_collection.get());
  fine_h1_space =
      std::make_unique<FiniteElementSpace>(space_op.GetMesh(), fine_h1_collection.get());

  auto &coarse_nd = space_op.GetNDSpace().Get();
  MFEM_VERIFY(coarse_nd.GetConformingProlongation() == nullptr &&
                  fine_nd_space->Get().GetConformingProlongation() == nullptr,
              "The hierarchical Maxwell domain adapter requires a conforming mesh without "
              "hanging-node local-to-true constraints!");
  injection = fem::hierarchical::BuildSparsePInjection(space_op.GetMesh().Get(), coarse_nd,
                                                       fine_nd_space->Get());

  const auto &materials = space_op.GetSingularMaterialCoefficients();
  const auto &imaginary_materials = space_op.GetSingularImagMaterialCoefficients();
  const auto &absolute_materials = space_op.GetSingularAbsMaterialCoefficients();
  MFEM_VERIFY(materials.size() == static_cast<std::size_t>(space_op.GetMesh().GetNE()) &&
                  imaginary_materials.size() == materials.size() &&
                  absolute_materials.size() == materials.size(),
              "Hierarchical Maxwell material batches do not cover the local mesh!");
  // This initial domain adapter handles all-isotropic devices (including the target
  // transmon). Anisotropic standard-only elements will use tensor MFEM integrators in the
  // parallel adapter rather than the scalar retained batches.
  for (int element = 0; element < space_op.GetMesh().GetNE(); element++)
  {
    const int attribute = space_op.GetMesh().Get().GetAttribute(element);
    MFEM_VERIFY(space_op.GetMaterialOp().IsIsotropic(attribute),
                "The serial hierarchical Maxwell domain adapter currently requires "
                "isotropic materials on every element! Domain attribute: "
                    << attribute);
  }

  std::vector<fem::singular::LocalSparseEnrichmentMatrices> local;
  if (dimension == 3)
  {
    const auto *topology = space_op.GetSingularDofTopology();
    MFEM_VERIFY(topology, "Missing tetrahedral singular topology!");
    local = fem::singular::AssembleLocalSparseEnrichmentMatricesBatch(
        *topology, fine_h1_space->Get(), fine_nd_space->Get(),
        {materials, imaginary_materials, absolute_materials},
        space_op.GetSingularAssemblyOptions(),
        fem::singular::RetainAllNDElementPatchBatches);
    element_enrichment_guests.resize(topology->elements.size());
    for (std::size_t element = 0; element < topology->elements.size(); element++)
    {
      for (const auto &dof : topology->elements[element].nd)
      {
        element_enrichment_guests[element].push_back(static_cast<int>(dof.dof));
      }
    }
  }
  else
  {
    const auto *topology = space_op.GetTriangleSingularDofTopology();
    MFEM_VERIFY(topology, "Missing triangular singular topology!");
    local = fem::singular::AssembleLocalSparseEnrichmentMatricesBatch(
        *topology, fine_h1_space->Get(), fine_nd_space->Get(),
        {materials, imaginary_materials, absolute_materials},
        space_op.GetSingularAssemblyOptions(),
        fem::singular::RetainAllNDElementPatchBatches);
    element_enrichment_guests.resize(topology->elements.size());
    for (std::size_t element = 0; element < topology->elements.size(); element++)
    {
      for (const auto &dof : topology->elements[element].nd)
      {
        element_enrichment_guests[element].push_back(static_cast<int>(dof.dof));
      }
    }
  }
  MFEM_VERIFY(local.size() == 3,
              "Hierarchical Maxwell domain assembly did not return all material slots!");

  enrichment_size =
      dimension == 3
          ? static_cast<int>(space_op.GetSingularDofTopology()->nd_dofs.size())
          : static_cast<int>(space_op.GetTriangleSingularDofTopology()->nd_dofs.size());
  const int fine_standard_size = fine_nd_space->GetVSize();
  const int coarse_standard_size = coarse_nd.GetVSize();
  MFEM_VERIFY(static_cast<int>(space_op.GetSingularParallelNumbering().nd.local_size) ==
                  enrichment_size,
              "Hierarchical Maxwell local enrichment numbering is inconsistent!");
  enrichment_prolongation = fem::singular::BuildParallelEnrichmentProlongation(
      space_op.GetComm(), space_op.GetSingularParallelNumbering().nd);

  // Essential true DOFs for the p+1 standard space and the order-independent singular
  // space. The serial shared-engine masks are populated only when local and true layouts
  // coincide; parallel lifting uses these true lists directly.
  // GetEssentialTrueDofs consumes a boundary marker, not an attribute list.
  const int boundary_attribute_maximum = space_op.GetMesh().Get().bdr_attributes.Size()
                                             ? space_op.GetMesh().Get().bdr_attributes.Max()
                                             : 0;
  const auto dbc_marker =
      mesh::AttrToMarker(boundary_attribute_maximum, space_op.GetNDDbcAttributes());
  fine_nd_space->Get().GetEssentialTrueDofs(dbc_marker, fine_standard_essential_true_dofs);
  const int coarse_standard_true_size = coarse_nd.GetTrueVSize();
  for (int dof : space_op.GetCombinedNDDbcTDofList())
  {
    if (dof >= coarse_standard_true_size)
    {
      enrichment_essential_true_dofs.Append(dof - coarse_standard_true_size);
    }
  }
  const bool serial_identity_layout =
      Mpi::Size(space_op.GetComm()) == 1 &&
      coarse_nd.GetVSize() == coarse_nd.GetTrueVSize() &&
      fine_nd_space->GetVSize() == fine_nd_space->GetTrueVSize() &&
      static_cast<int>(space_op.GetSingularParallelNumbering().nd.owned_size) ==
          enrichment_size;
  entity_patches_available = serial_identity_layout;
  if (serial_identity_layout)
  {
    coarse_essential.assign(coarse_standard_size + enrichment_size, false);
    for (int dof : space_op.GetCombinedNDDbcTDofList())
    {
      coarse_essential[dof] = true;
    }
    fine_essential.assign(fine_standard_size + enrichment_size, false);
    for (int dof : fine_standard_essential_true_dofs)
    {
      fine_essential[dof] = true;
    }
    for (int dof : enrichment_essential_true_dofs)
    {
      fine_essential[fine_standard_size + dof] = true;
    }
  }

  // Map retained enriched patches in each material batch by element. The assembly retains
  // only elements with singular basis incidence; ordinary elements are completed below.
  std::vector<std::vector<int>> retained_index(
      3, std::vector<int>(space_op.GetMesh().GetNE(), -1));
  for (int batch = 0; batch < 3; batch++)
  {
    for (std::size_t patch = 0; patch < local[batch].nd_element_patches.size(); patch++)
    {
      retained_index[batch][local[batch].nd_element_patches[patch].element] =
          static_cast<int>(patch);
    }
  }

  mfem::CurlCurlIntegrator curl_integrator;
  mfem::VectorFEMassIntegrator mass_integrator;
  mfem::Array<int> standard_dofs;
  elements.resize(space_op.GetMesh().GetNE());
  for (int element = 0; element < space_op.GetMesh().GetNE(); element++)
  {
    auto &data = elements[element];
    data.support_element = element;
    if (retained_index[0][element] >= 0)
    {
      MFEM_VERIFY(retained_index[1][element] >= 0 && retained_index[2][element] >= 0,
                  "Retained hierarchical material slots disagree on enriched elements!");
      const auto &real = local[0].nd_element_patches[retained_index[0][element]];
      const auto &imaginary = local[1].nd_element_patches[retained_index[1][element]];
      const auto &absolute = local[2].nd_element_patches[retained_index[2][element]];
      const int local_size = real.mass.Height();
      standard_dofs = real.standard_dofs;
      data.standard_dofs = real.standard_dofs;
      data.enrichment_dofs = real.enrichment_dofs;
      for (int dof : standard_dofs)
      {
        data.dofs.push_back(UnsignedDof(dof));
      }
      for (int dof : real.enrichment_dofs)
      {
        data.dofs.push_back(fine_standard_size + dof);
      }
      SetSignedMatrix(real.curl_curl, standard_dofs, local_size, data.curl_curl);
      SetSignedMatrix(real.mass, standard_dofs, local_size, data.mass_real);
      SetSignedMatrix(imaginary.mass, standard_dofs, local_size, data.mass_imag);
      SetSignedMatrix(absolute.mass, standard_dofs, local_size, data.mass_abs);
      continue;
    }

    mfem::DofTransformation transformation;
    fine_nd_space->Get().GetElementVDofs(element, standard_dofs, transformation);
    data.standard_dofs = standard_dofs;
    for (int dof : standard_dofs)
    {
      data.dofs.push_back(UnsignedDof(dof));
    }
    auto &T = *space_op.GetMesh().Get().GetElementTransformation(element);
    const auto &fe = *fine_nd_space->Get().GetFE(element);
    mfem::DenseMatrix unweighted_curl, unweighted_mass;
    curl_integrator.AssembleElementMatrix(fe, T, unweighted_curl);
    mass_integrator.AssembleElementMatrix(fe, T, unweighted_mass);
    transformation.TransformDual(unweighted_curl);
    transformation.TransformDual(unweighted_mass);
    const int local_size = standard_dofs.Size();
    const auto set_standard =
        [&](const mfem::DenseMatrix &unweighted, double scale, mfem::DenseMatrix &matrix)
    {
      matrix.SetSize(local_size);
      for (int row = 0; row < local_size; row++)
      {
        for (int column = 0; column < local_size; column++)
        {
          matrix(row, column) = DofSign(standard_dofs[row]) *
                                DofSign(standard_dofs[column]) * scale *
                                unweighted(row, column);
        }
      }
    };
    set_standard(unweighted_curl, materials[element].inverse_magnetic, data.curl_curl);
    set_standard(unweighted_mass, materials[element].electric, data.mass_real);
    set_standard(unweighted_mass, imaginary_materials[element].electric, data.mass_imag);
    set_standard(unweighted_mass, absolute_materials[element].electric, data.mass_abs);
  }

  const auto lumped_stiffness = space_op.GetLumpedPortOp().GetStiffnessBdrCoefficientMap();
  const auto lumped_damping = space_op.GetLumpedPortOp().GetDampingBdrCoefficientMap();
  const auto lumped_mass = space_op.GetLumpedPortOp().GetMassBdrCoefficientMap();
  const auto impedance_stiffness =
      space_op.GetSurfaceImpedanceOp().GetStiffnessBdrCoefficientMap();
  const auto impedance_damping =
      space_op.GetSurfaceImpedanceOp().GetDampingBdrCoefficientMap();
  const auto impedance_mass = space_op.GetSurfaceImpedanceOp().GetMassBdrCoefficientMap();
  const bool has_polynomial_boundary =
      !lumped_stiffness.empty() || !lumped_damping.empty() || !lumped_mass.empty() ||
      !impedance_stiffness.empty() || !impedance_damping.empty() || !impedance_mass.empty();
  unsupported_polynomial_boundary = dimension != 3 && has_polynomial_boundary;
  if (dimension == 3)
  {
    const auto *topology = space_op.GetSingularDofTopology();
    const std::set<int> no_excluded;
    const auto &excluded = space_op.GetSingularConstrainedImpedanceAttributes();
    const auto assemble_boundary =
        [&](const std::map<int, double> &coefficients,
            std::vector<fem::singular::LocalNDBoundaryPatchMatrices> &patches,
            const std::set<int> &excluded_singular_attributes)
    {
      if (!coefficients.empty())
      {
        (void)fem::singular::AssembleLocalSparseNDBoundaryMassMatrices(
            *topology, fine_nd_space->Get(), coefficients,
            space_op.GetSingularAssemblyOptions(), &patches, excluded_singular_attributes);
      }
    };
    assemble_boundary(lumped_stiffness, boundary_stiffness, no_excluded);
    assemble_boundary(impedance_stiffness, boundary_stiffness, excluded);
    assemble_boundary(lumped_damping, boundary_damping, no_excluded);
    assemble_boundary(impedance_damping, boundary_damping, excluded);
    assemble_boundary(lumped_mass, boundary_mass, no_excluded);
    assemble_boundary(impedance_mass, boundary_mass, excluded);
    assemble_boundary(AbsoluteBoundaryCoefficients(lumped_stiffness),
                      boundary_stiffness_abs, no_excluded);
    assemble_boundary(AbsoluteBoundaryCoefficients(impedance_stiffness),
                      boundary_stiffness_abs, excluded);
    assemble_boundary(AbsoluteBoundaryCoefficients(lumped_damping), boundary_damping_abs,
                      no_excluded);
    assemble_boundary(AbsoluteBoundaryCoefficients(impedance_damping), boundary_damping_abs,
                      excluded);
    assemble_boundary(AbsoluteBoundaryCoefficients(lumped_mass), boundary_mass_abs,
                      no_excluded);
    assemble_boundary(AbsoluteBoundaryCoefficients(impedance_mass), boundary_mass_abs,
                      excluded);
  }
}

HierarchicalMaxwellDomainData::~HierarchicalMaxwellDomainData() = default;

int HierarchicalMaxwellDomainData::GetFineStandardSize() const
{
  return fine_nd_space->GetVSize();
}

std::vector<fem::hierarchical::ComplexLocalOperatorContribution>
HierarchicalMaxwellDomainData::BuildComplexDomainContributions(
    std::complex<double> omega) const
{
  const std::complex<double> mass_scale = -omega * omega;
  std::vector<fem::hierarchical::ComplexLocalOperatorContribution> contributions;
  contributions.reserve(elements.size());
  for (const auto &element : elements)
  {
    auto &data = contributions.emplace_back();
    data.support_element = element.support_element;
    data.dofs = element.dofs;
    data.matrix_real = element.curl_curl;
    data.matrix_real.Add(mass_scale.real(), element.mass_real);
    data.matrix_real.Add(-mass_scale.imag(), element.mass_imag);
    data.matrix_imag.SetSize(element.curl_curl.Height());
    data.matrix_imag = 0.0;
    data.matrix_imag.Add(mass_scale.imag(), element.mass_real);
    data.matrix_imag.Add(mass_scale.real(), element.mass_imag);
    data.rhs_real.SetSize(static_cast<int>(data.dofs.size()));
    data.rhs_imag.SetSize(static_cast<int>(data.dofs.size()));
    data.rhs_real = 0.0;
    data.rhs_imag = 0.0;
  }
  return contributions;
}

std::vector<fem::hierarchical::LocalOperatorContribution>
HierarchicalMaxwellDomainData::BuildDomainMetricContributions(
    std::complex<double> omega) const
{
  MFEM_VERIFY(std::abs(omega) > 0.0,
              "Hierarchical Maxwell graph metric requires nonzero frequency so mass "
              "controls the curl-gradient nullspace!");
  const double mass_scale = std::norm(omega);
  std::vector<fem::hierarchical::LocalOperatorContribution> contributions;
  contributions.reserve(elements.size());
  for (const auto &element : elements)
  {
    auto &data = contributions.emplace_back();
    data.support_element = element.support_element;
    data.dofs = element.dofs;
    data.matrix = element.curl_curl;
    data.matrix.Add(mass_scale, element.mass_abs);
    data.rhs.SetSize(static_cast<int>(data.dofs.size()));
    data.rhs = 0.0;
  }
  return contributions;
}

std::vector<fem::hierarchical::ComplexLocalOperatorContribution>
HierarchicalMaxwellDomainData::BuildComplexPolynomialContributions(
    std::complex<double> omega) const
{
  MFEM_VERIFY(!unsupported_polynomial_boundary,
              "Hierarchical Maxwell polynomial boundary contributions are not yet "
              "implemented in two dimensions!");
  auto contributions = BuildComplexDomainContributions(omega);
  const auto append =
      [&](const std::vector<fem::singular::LocalNDBoundaryPatchMatrices> &patches,
          std::complex<double> scale)
  {
    for (const auto &patch : patches)
    {
      auto &data = contributions.emplace_back();
      data.support_element = patch.element;
      data.dofs = patch.dofs;
      data.matrix_real = patch.mass;
      data.matrix_real *= scale.real();
      data.matrix_imag = patch.mass;
      data.matrix_imag *= scale.imag();
      data.rhs_real.SetSize(static_cast<int>(data.dofs.size()));
      data.rhs_imag.SetSize(static_cast<int>(data.dofs.size()));
      data.rhs_real = 0.0;
      data.rhs_imag = 0.0;
    }
  };
  append(boundary_stiffness, 1.0);
  append(boundary_damping, std::complex<double>(0.0, 1.0) * omega);
  append(boundary_mass, -omega * omega);
  return contributions;
}

std::vector<fem::hierarchical::LocalOperatorContribution>
HierarchicalMaxwellDomainData::BuildPolynomialMetricContributions(
    std::complex<double> omega) const
{
  MFEM_VERIFY(!unsupported_polynomial_boundary,
              "Hierarchical Maxwell polynomial boundary contributions are not yet "
              "implemented in two dimensions!");
  auto contributions = BuildDomainMetricContributions(omega);
  const auto append =
      [&](const std::vector<fem::singular::LocalNDBoundaryPatchMatrices> &patches,
          double scale)
  {
    for (const auto &patch : patches)
    {
      auto &data = contributions.emplace_back();
      data.support_element = patch.element;
      data.dofs = patch.dofs;
      data.matrix = patch.mass;
      data.matrix *= scale;
      data.rhs.SetSize(static_cast<int>(data.dofs.size()));
      data.rhs = 0.0;
    }
  };
  append(boundary_stiffness_abs, 1.0);
  append(boundary_damping_abs, std::abs(omega));
  append(boundary_mass_abs, std::norm(omega));
  return contributions;
}

std::vector<fem::singular::LocalNDElementPatchMatrices>
HierarchicalMaxwellDomainData::BuildPolynomialMetricElementPatches(
    std::complex<double> omega) const
{
  const auto contributions = BuildPolynomialMetricContributions(omega);
  std::vector<mfem::DenseMatrix> assembled(elements.size());
  for (std::size_t element = 0; element < elements.size(); element++)
  {
    assembled[element].SetSize(static_cast<int>(elements[element].dofs.size()));
    assembled[element] = 0.0;
  }
  for (const auto &contribution : contributions)
  {
    MFEM_VERIFY(contribution.support_element >= 0 &&
                    contribution.support_element < static_cast<int>(elements.size()),
                "Hierarchical Maxwell metric contribution has invalid support!");
    const auto &element = elements[contribution.support_element];
    auto &matrix = assembled[contribution.support_element];
    for (int row = 0; row < static_cast<int>(contribution.dofs.size()); row++)
    {
      const auto found_row =
          std::find(element.dofs.begin(), element.dofs.end(), contribution.dofs[row]);
      MFEM_VERIFY(found_row != element.dofs.end(),
                  "Hierarchical boundary metric row is absent from its support element!");
      const int element_row = static_cast<int>(found_row - element.dofs.begin());
      for (int column = 0; column < static_cast<int>(contribution.dofs.size()); column++)
      {
        const auto found_column =
            std::find(element.dofs.begin(), element.dofs.end(), contribution.dofs[column]);
        MFEM_VERIFY(found_column != element.dofs.end(),
                    "Hierarchical boundary metric column is absent from its support "
                    "element!");
        const int element_column = static_cast<int>(found_column - element.dofs.begin());
        matrix(element_row, element_column) += contribution.matrix(row, column);
      }
    }
  }

  std::vector<fem::singular::LocalNDElementPatchMatrices> patches;
  patches.reserve(elements.size());
  for (std::size_t element_index = 0; element_index < elements.size(); element_index++)
  {
    const auto &element = elements[element_index];
    auto &patch = patches.emplace_back();
    patch.element = static_cast<int>(element_index);
    patch.standard_dofs = element.standard_dofs;
    patch.enrichment_dofs = element.enrichment_dofs;
    const int local_size = static_cast<int>(element.dofs.size());
    patch.curl_curl.SetSize(local_size);
    for (int row = 0; row < local_size; row++)
    {
      const double row_sign =
          row < patch.standard_dofs.Size() ? DofSign(patch.standard_dofs[row]) : 1.0;
      for (int column = 0; column < local_size; column++)
      {
        const double column_sign = column < patch.standard_dofs.Size()
                                       ? DofSign(patch.standard_dofs[column])
                                       : 1.0;
        // ParallelElementPatchInverse applies these signs while gathering/scattering, so
        // convert the final unsigned matrix back to its local element orientation here.
        patch.curl_curl(row, column) =
            row_sign * column_sign * assembled[element_index](row, column);
      }
    }
    patch.mass.SetSize(local_size);
    patch.mass = 0.0;
    patch.mass_estimated_absolute_error.SetSize(local_size);
    patch.mass_estimated_absolute_error = 0.0;
    patch.curl_curl_estimated_absolute_error.SetSize(local_size);
    patch.curl_curl_estimated_absolute_error = 0.0;
  }
  return patches;
}

HierarchicalMaxwellDomainData::ParallelEstimate
HierarchicalMaxwellDomainData::LiftLocalComplexResidual(
    std::complex<double> omega, fem::hierarchical::ComplexResidual residual,
    PatchShape patch_shape) const
{
  auto &coarse_fespace = space_op->GetNDSpace().Get();
  auto &fine_fespace = fine_nd_space->Get();
  const auto &numbering = space_op->GetSingularParallelNumbering().nd;
  const int fine_standard_true = fine_fespace.GetTrueVSize();
  const int enrichment_owned = static_cast<int>(numbering.owned_size);

  if (patch_shape == PatchShape::ENTITY)
  {
    // Certified edge/face/interior lifting through the shared engine, currently limited
    // to one rank with identity local-to-true layouts. On that layout the local residual
    // is complete, so essential elimination happens directly in local coordinates.
    MFEM_VERIFY(entity_patches_available,
                "Entity-patch hierarchical lifting requires one rank with identity "
                "local-to-true layouts; use the element-patch shape instead!");
    for (int dof = 0; dof < residual.real.Size(); dof++)
    {
      if (fine_essential[dof])
      {
        residual.real(dof) = 0.0;
        residual.imag(dof) = 0.0;
      }
    }
    const auto metric = BuildPolynomialMetricContributions(omega);
    const auto lifting = fem::hierarchical::LiftComplexResidualByPatches(
        space_op->GetMesh().Get(), coarse_fespace, fine_fespace, injection, metric,
        fine_essential, coarse_essential, residual, element_enrichment_guests);
    ParallelEstimate estimate;
    estimate.indicator_energy.SetSize(static_cast<int>(lifting.indicator.size()));
    for (std::size_t element = 0; element < lifting.indicator.size(); element++)
    {
      estimate.indicator_energy(static_cast<int>(element)) = lifting.indicator[element];
    }
    estimate.total_energy = lifting.energy;
    return estimate;
  }

  const auto local_to_true = [&](const Vector &local, Vector &combined_true)
  {
    Vector standard_local(fine_fespace.GetVSize());
    Vector enrichment_local(enrichment_size);
    for (int dof = 0; dof < standard_local.Size(); dof++)
    {
      standard_local(dof) = local(dof);
    }
    for (int dof = 0; dof < enrichment_size; dof++)
    {
      enrichment_local(dof) = local(fine_fespace.GetVSize() + dof);
    }
    Vector standard_true(fine_standard_true), enrichment_true(enrichment_owned);
    fine_fespace.Dof_TrueDof_Matrix()->MultTranspose(standard_local, standard_true);
    enrichment_prolongation->MultTranspose(enrichment_local, enrichment_true);
    combined_true.SetSize(fine_standard_true + enrichment_owned);
    for (int dof = 0; dof < fine_standard_true; dof++)
    {
      combined_true(dof) = standard_true(dof);
    }
    for (int dof = 0; dof < enrichment_owned; dof++)
    {
      combined_true(fine_standard_true + dof) = enrichment_true(dof);
    }
    for (int dof : fine_standard_essential_true_dofs)
    {
      combined_true(dof) = 0.0;
    }
    for (int dof : enrichment_essential_true_dofs)
    {
      combined_true(fine_standard_true + dof) = 0.0;
    }
  };

  Vector residual_true_real, residual_true_imag;
  local_to_true(residual.real, residual_true_real);
  local_to_true(residual.imag, residual_true_imag);
  auto metric_patches = BuildPolynomialMetricElementPatches(omega);
  fem::singular::ParallelElementPatchInverse inverse(
      fine_fespace, numbering, metric_patches, 1.0, 1.0, fine_standard_essential_true_dofs,
      enrichment_essential_true_dofs);
  Vector correction_true_real, correction_true_imag;
  inverse.Mult(residual_true_real, correction_true_real);
  inverse.Mult(residual_true_imag, correction_true_imag);

  const auto true_to_local = [&](const Vector &combined_true, Vector &local)
  {
    Vector standard_true(fine_standard_true), enrichment_true(enrichment_owned);
    for (int dof = 0; dof < fine_standard_true; dof++)
    {
      standard_true(dof) = combined_true(dof);
    }
    for (int dof = 0; dof < enrichment_owned; dof++)
    {
      enrichment_true(dof) = combined_true(fine_standard_true + dof);
    }
    Vector standard_local(fine_fespace.GetVSize()), enrichment_local(enrichment_size);
    fine_fespace.Dof_TrueDof_Matrix()->Mult(standard_true, standard_local);
    enrichment_prolongation->Mult(enrichment_true, enrichment_local);
    local.SetSize(fine_fespace.GetVSize() + enrichment_size);
    for (int dof = 0; dof < standard_local.Size(); dof++)
    {
      local(dof) = standard_local(dof);
    }
    for (int dof = 0; dof < enrichment_size; dof++)
    {
      local(fine_fespace.GetVSize() + dof) = enrichment_local(dof);
    }
  };

  Vector correction_local_real, correction_local_imag;
  true_to_local(correction_true_real, correction_local_real);
  true_to_local(correction_true_imag, correction_local_imag);
  const auto metric = BuildPolynomialMetricContributions(omega);
  ParallelEstimate estimate;
  estimate.indicator_energy.SetSize(space_op->GetMesh().GetNE());
  estimate.indicator_energy = 0.0;
  for (const auto &contribution : metric)
  {
    Vector local_real(static_cast<int>(contribution.dofs.size()));
    Vector local_imag(static_cast<int>(contribution.dofs.size()));
    for (int i = 0; i < local_real.Size(); i++)
    {
      local_real(i) = correction_local_real(contribution.dofs[i]);
      local_imag(i) = correction_local_imag(contribution.dofs[i]);
    }
    Vector action(local_real.Size());
    contribution.matrix.Mult(local_real, action);
    estimate.indicator_energy(contribution.support_element) += local_real * action;
    contribution.matrix.Mult(local_imag, action);
    estimate.indicator_energy(contribution.support_element) += local_imag * action;
  }
  estimate.total_energy = estimate.indicator_energy.Sum();
  Mpi::GlobalSum(1, &estimate.total_energy, space_op->GetComm());
  return estimate;
}

namespace
{

// Inject one component of the coarse combined true field into the uneliminated local
// combined p+1 layout.
void InjectCoarseComponent(mfem::ParFiniteElementSpace &coarse_fespace,
                           mfem::ParFiniteElementSpace &fine_fespace,
                           const mfem::HypreParMatrix &enrichment_prolongation,
                           int enrichment_size, const Vector &coarse_true,
                           Vector &fine_local)
{
  const int coarse_standard_true = coarse_fespace.GetTrueVSize();
  const int enrichment_owned = enrichment_prolongation.Width();
  Vector coarse_standard_true_vector(coarse_standard_true);
  Vector coarse_enrichment_true(enrichment_owned);
  for (int dof = 0; dof < coarse_standard_true; dof++)
  {
    coarse_standard_true_vector(dof) = coarse_true(dof);
  }
  for (int dof = 0; dof < enrichment_owned; dof++)
  {
    coarse_enrichment_true(dof) = coarse_true(coarse_standard_true + dof);
  }
  Vector coarse_standard_local(coarse_fespace.GetVSize());
  coarse_fespace.Dof_TrueDof_Matrix()->Mult(coarse_standard_true_vector,
                                            coarse_standard_local);
  Vector fine_standard_local(fine_fespace.GetVSize());
  mfem::PRefinementTransferOperator transfer(coarse_fespace, fine_fespace);
  transfer.Mult(coarse_standard_local, fine_standard_local);
  Vector fine_enrichment_local(enrichment_size);
  enrichment_prolongation.Mult(coarse_enrichment_true, fine_enrichment_local);
  fine_local.SetSize(fine_fespace.GetVSize() + enrichment_size);
  for (int dof = 0; dof < fine_fespace.GetVSize(); dof++)
  {
    fine_local(dof) = fine_standard_local(dof);
  }
  for (int dof = 0; dof < enrichment_size; dof++)
  {
    fine_local(fine_fespace.GetVSize() + dof) = fine_enrichment_local(dof);
  }
}

}  // namespace

HierarchicalMaxwellDomainData::ParallelEstimate
HierarchicalMaxwellDomainData::EstimatePolynomialEigenResidual(
    std::complex<double> omega, const ComplexVector &coarse_field,
    PatchShape patch_shape) const
{
  auto &coarse_fespace = space_op->GetNDSpace().Get();
  auto &fine_fespace = fine_nd_space->Get();
  const auto &numbering = space_op->GetSingularParallelNumbering().nd;
  const int coarse_standard_true = coarse_fespace.GetTrueVSize();
  const int enrichment_owned = static_cast<int>(numbering.owned_size);
  MFEM_VERIFY(coarse_field.Size() == coarse_standard_true + enrichment_owned,
              "Hierarchical Maxwell estimator received an invalid coarse true field!");
  Vector injected_real, injected_imag;
  InjectCoarseComponent(coarse_fespace, fine_fespace, *enrichment_prolongation,
                        enrichment_size, coarse_field.Real(), injected_real);
  InjectCoarseComponent(coarse_fespace, fine_fespace, *enrichment_prolongation,
                        enrichment_size, coarse_field.Imag(), injected_imag);
  const auto physical = BuildComplexPolynomialContributions(omega);
  const std::vector<bool> no_local_essential(injected_real.Size(), false);
  auto residual = fem::hierarchical::AssembleComplexResidual(
      injected_real.Size(), physical, injected_real, injected_imag, no_local_essential);
  return LiftLocalComplexResidual(omega, std::move(residual), patch_shape);
}

fem::hierarchical::ComplexResidual
HierarchicalMaxwellDomainData::AssembleFineExcitation(int excitation_idx, double omega,
                                                      const ComplexVector &coarse_rhs) const
{
  auto &fine_fespace = fine_nd_space->Get();
  const auto &numbering = space_op->GetSingularParallelNumbering().nd;
  const int coarse_standard_true = space_op->GetNDSpace().GetTrueVSize();
  const int enrichment_owned = static_cast<int>(numbering.owned_size);
  MFEM_VERIFY(coarse_rhs.Size() == coarse_standard_true + enrichment_owned,
              "Hierarchical Maxwell excitation received an invalid coarse right-hand "
              "side!");
  const int local_size = fine_fespace.GetVSize() + enrichment_size;
  fem::hierarchical::ComplexResidual excitation{Vector(local_size), Vector(local_size)};
  excitation.real = 0.0;
  excitation.imag = 0.0;

  // Standard rows: the same boundary excitation coefficients as the coarse assembly,
  // integrated against p+1 test functions. b = i omega RHS1 puts the load on the
  // imaginary slot for the real mode coefficients used by lumped ports.
  SumVectorCoefficient fb(space_op->GetMesh().SpaceDimension());
  space_op->GetLumpedPortOp().AddExcitationBdrCoefficients(excitation_idx, fb);
  space_op->GetSurfaceCurrentOp().AddExcitationBdrCoefficients(fb);
  int empty = fb.empty();
  Mpi::GlobalMin(1, &empty, space_op->GetComm());
  if (!empty)
  {
    mfem::ParLinearForm rhs1(&fine_fespace);
    rhs1.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb));
    rhs1.UseFastAssembly(false);
    rhs1.Assemble();
    for (int dof = 0; dof < fine_fespace.GetVSize(); dof++)
    {
      excitation.imag(dof) = omega * rhs1[dof];
    }
  }

  // Enrichment rows are independent of the standard order: splice the exact coarse
  // combined entries (already i omega scaled and essential-eliminated).
  Vector owned_real(enrichment_owned), owned_imag(enrichment_owned);
  for (int dof = 0; dof < enrichment_owned; dof++)
  {
    owned_real(dof) = coarse_rhs.Real()(coarse_standard_true + dof);
    owned_imag(dof) = coarse_rhs.Imag()(coarse_standard_true + dof);
  }
  Vector local_real(enrichment_size), local_imag(enrichment_size);
  enrichment_prolongation->Mult(owned_real, local_real);
  enrichment_prolongation->Mult(owned_imag, local_imag);
  for (int dof = 0; dof < enrichment_size; dof++)
  {
    excitation.real(fine_fespace.GetVSize() + dof) += local_real(dof);
    excitation.imag(fine_fespace.GetVSize() + dof) += local_imag(dof);
  }
  return excitation;
}

HierarchicalMaxwellDomainData::ParallelEstimate
HierarchicalMaxwellDomainData::EstimatePolynomialDrivenResidual(
    double omega, const ComplexVector &coarse_field, int excitation_idx,
    const ComplexVector &coarse_rhs, PatchShape patch_shape) const
{
  auto &coarse_fespace = space_op->GetNDSpace().Get();
  auto &fine_fespace = fine_nd_space->Get();
  const auto &numbering = space_op->GetSingularParallelNumbering().nd;
  const int coarse_standard_true = coarse_fespace.GetTrueVSize();
  const int enrichment_owned = static_cast<int>(numbering.owned_size);
  MFEM_VERIFY(coarse_field.Size() == coarse_standard_true + enrichment_owned,
              "Hierarchical Maxwell estimator received an invalid coarse true field!");
  Vector injected_real, injected_imag;
  InjectCoarseComponent(coarse_fespace, fine_fespace, *enrichment_prolongation,
                        enrichment_size, coarse_field.Real(), injected_real);
  InjectCoarseComponent(coarse_fespace, fine_fespace, *enrichment_prolongation,
                        enrichment_size, coarse_field.Imag(), injected_imag);
  const auto physical = BuildComplexPolynomialContributions({omega, 0.0});
  const std::vector<bool> no_local_essential(injected_real.Size(), false);
  auto residual = fem::hierarchical::AssembleComplexResidual(
      injected_real.Size(), physical, injected_real, injected_imag, no_local_essential);
  const auto excitation = AssembleFineExcitation(excitation_idx, omega, coarse_rhs);
  residual.real += excitation.real;
  residual.imag += excitation.imag;
  return LiftLocalComplexResidual({omega, 0.0}, std::move(residual), patch_shape);
}

}  // namespace palace
