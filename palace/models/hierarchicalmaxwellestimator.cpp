// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hierarchicalmaxwellestimator.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

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
  // Constrained (hanging-node) prolongation rows are not yet certified for the record
  // congruences, and NCMesh-backed spaces must not be probed through the serial
  // conforming-interpolation path.
  MFEM_VERIFY(space_op.GetMesh().Get().Conforming(),
              "The hierarchical Maxwell estimator requires a conforming mesh; "
              "nonconforming singular AMR keeps the sliced recovery estimator!");
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
  entity_patches_available = true;
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
      // The three material slots share one mass block up to their electric scalars, so
      // store the real-weighted block once with the scalar ratios.
      SetSignedMatrix(real.mass, standard_dofs, local_size, data.mass);
      (void)imaginary;
      (void)absolute;
      data.electric_real = 1.0;
      data.electric_imag =
          imaginary_materials[element].electric / materials[element].electric;
      data.electric_abs =
          absolute_materials[element].electric / materials[element].electric;
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
    set_standard(unweighted_mass, 1.0, data.mass);
    data.electric_real = materials[element].electric;
    data.electric_imag = imaginary_materials[element].electric;
    data.electric_abs = absolute_materials[element].electric;
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

void HierarchicalMaxwellDomainData::LocalToTrue(const Vector &local,
                                                Vector &combined_true) const
{
  auto &fine_fespace = fine_nd_space->Get();
  const auto &numbering = space_op->GetSingularParallelNumbering().nd;
  const int fine_standard_true = fine_fespace.GetTrueVSize();
  const int enrichment_owned = static_cast<int>(numbering.owned_size);
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
}

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
    data.matrix_real.Add(mass_scale.real() * element.electric_real -
                             mass_scale.imag() * element.electric_imag,
                         element.mass);
    data.matrix_imag.SetSize(element.curl_curl.Height());
    data.matrix_imag = 0.0;
    data.matrix_imag.Add(mass_scale.imag() * element.electric_real +
                             mass_scale.real() * element.electric_imag,
                         element.mass);
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
    data.matrix.Add(mass_scale * element.electric_abs, element.mass);
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

std::vector<HierarchicalMaxwellDomainData::TrueElementRecord>
HierarchicalMaxwellDomainData::BuildTrueMetricElementRecords(
    std::complex<double> omega) const
{
  auto &coarse_fespace = space_op->GetNDSpace().Get();
  auto &fine_fespace = fine_nd_space->Get();
  const auto &numbering = space_op->GetSingularParallelNumbering().nd;
  const auto &mesh = space_op->GetMesh().Get();

  // Per-ldof rows of the conforming prolongation. Rows are single +-1 entries for edge
  // DOFs, but shared tetrahedral face DOFs couple through 2x2 orientation blocks across
  // ranks, so records fold local-to-true congruences instead of scalar signs.
  const auto build_ldof_rows = [](mfem::ParFiniteElementSpace &fespace)
  {
    std::vector<std::vector<std::pair<HYPRE_BigInt, double>>> rows(fespace.GetVSize());
    mfem::SparseMatrix diag, offd;
    HYPRE_BigInt *colmap = nullptr;
    fespace.Dof_TrueDof_Matrix()->GetDiag(diag);
    fespace.Dof_TrueDof_Matrix()->GetOffd(offd, colmap);
    const HYPRE_BigInt true_offset = fespace.GetMyTDofOffset();
    for (int ldof = 0; ldof < fespace.GetVSize(); ldof++)
    {
      for (int entry = diag.GetI()[ldof]; entry < diag.GetI()[ldof + 1]; entry++)
      {
        if (std::abs(diag.GetData()[entry]) > 1.0e-14)
        {
          rows[ldof].push_back({true_offset + diag.GetJ()[entry], diag.GetData()[entry]});
        }
      }
      for (int entry = offd.GetI()[ldof]; entry < offd.GetI()[ldof + 1]; entry++)
      {
        if (std::abs(offd.GetData()[entry]) > 1.0e-14)
        {
          rows[ldof].push_back({colmap[offd.GetJ()[entry]], offd.GetData()[entry]});
        }
      }
      MFEM_VERIFY(!rows[ldof].empty() && rows[ldof].size() <= 4,
                  "Conforming true-DOF prolongation has an unexpected row pattern!");
    }
    return rows;
  };
  const auto fine_rows = build_ldof_rows(fine_fespace);
  const auto coarse_rows = build_ldof_rows(coarse_fespace);
  const HYPRE_BigInt fine_standard_global = fine_fespace.GlobalTrueVSize();
  const HYPRE_BigInt coarse_standard_global = coarse_fespace.GlobalTrueVSize();

  // Rank-local essential flags for every local DOF, including non-owned shared copies.
  const int boundary_attribute_maximum =
      mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto dbc_marker =
      mesh::AttrToMarker(boundary_attribute_maximum, space_op->GetNDDbcAttributes());
  mfem::Array<int> fine_essential_vdofs, coarse_essential_vdofs;
  fine_fespace.GetEssentialVDofs(dbc_marker, fine_essential_vdofs);
  coarse_fespace.GetEssentialVDofs(dbc_marker, coarse_essential_vdofs);
  Vector enrichment_essential_owned(static_cast<int>(numbering.owned_size));
  enrichment_essential_owned = 0.0;
  for (int dof : enrichment_essential_true_dofs)
  {
    enrichment_essential_owned(dof) = 1.0;
  }
  Vector enrichment_essential_local(enrichment_size);
  enrichment_prolongation->Mult(enrichment_essential_owned, enrichment_essential_local);

  // Merge domain and facet metric blocks directly (still in the unsigned local
  // convention) without materializing a contribution vector: metric = curl-curl +
  // |omega| damping facets + |omega|^2 (mass facets and absolute domain mass).
  MFEM_VERIFY(std::abs(omega) > 0.0,
              "Hierarchical Maxwell graph metric requires nonzero frequency so mass "
              "controls the curl-gradient nullspace!");
  const double metric_mass_scale = std::norm(omega);
  std::vector<mfem::DenseMatrix> merged(mesh.GetNE());
  std::vector<std::map<int, int>> dof_positions(mesh.GetNE());
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    merged[element] = elements[element].curl_curl;
    merged[element].Add(metric_mass_scale * elements[element].electric_abs,
                        elements[element].mass);
    for (int i = 0; i < static_cast<int>(elements[element].dofs.size()); i++)
    {
      dof_positions[element].emplace(elements[element].dofs[i], i);
    }
  }
  const auto merge_boundary =
      [&](const std::vector<fem::singular::LocalNDBoundaryPatchMatrices> &patches,
          double scale)
  {
    for (const auto &facet : patches)
    {
      const auto &positions = dof_positions[facet.element];
      for (int i = 0; i < static_cast<int>(facet.dofs.size()); i++)
      {
        const auto row = positions.find(facet.dofs[i]);
        MFEM_VERIFY(row != positions.end(),
                    "Hierarchical facet metric DOF is absent from its support element!");
        for (int j = 0; j < static_cast<int>(facet.dofs.size()); j++)
        {
          const auto column = positions.find(facet.dofs[j]);
          MFEM_VERIFY(column != positions.end(),
                      "Hierarchical facet metric DOF is absent from its support element!");
          merged[facet.element](row->second, column->second) += scale * facet.mass(i, j);
        }
      }
    }
  };
  merge_boundary(boundary_stiffness_abs, 1.0);
  merge_boundary(boundary_damping_abs, std::abs(omega));
  merge_boundary(boundary_mass_abs, metric_mass_scale);

  std::vector<TrueElementRecord> records(mesh.GetNE());
  mfem::Array<int> coarse_dofs;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &record = records[element];
    record.element_id = mesh.GetGlobalElementNum(element);
    const auto &domain = elements[element];
    const int local_size = static_cast<int>(domain.dofs.size());

    // Element block of the local-to-true prolongation: standard DOFs may couple pairs of
    // true DOFs through shared-face orientation blocks; enrichment DOFs are identity.
    std::vector<HYPRE_BigInt> true_ids;
    std::vector<char> true_essential;
    std::map<HYPRE_BigInt, int> true_positions;
    const auto true_position = [&](HYPRE_BigInt id, bool essential)
    {
      const auto found = true_positions.find(id);
      if (found != true_positions.end())
      {
        return found->second;
      }
      const int position = static_cast<int>(true_ids.size());
      true_positions.emplace(id, position);
      true_ids.push_back(id);
      true_essential.push_back(essential);
      return position;
    };
    std::vector<std::vector<std::pair<int, double>>> local_to_true(local_size);
    for (int i = 0; i < local_size; i++)
    {
      const int dof = domain.dofs[i];
      if (dof < fine_fespace.GetVSize())
      {
        for (const auto &[id, value] : fine_rows[dof])
        {
          local_to_true[i].push_back(
              {true_position(id, fine_essential_vdofs[dof] != 0), value});
        }
      }
      else
      {
        const int enrichment_dof = dof - fine_fespace.GetVSize();
        local_to_true[i].push_back(
            {true_position(fine_standard_global + numbering.local_to_true[enrichment_dof],
                           std::abs(enrichment_essential_local(enrichment_dof)) > 0.5),
             1.0});
      }
    }
    const int true_size = static_cast<int>(true_ids.size());
    record.fine_dofs = true_ids;
    record.fine_essential = true_essential;

    // Congruence to the true basis: A_true = P_e^T A_local P_e.
    mfem::DenseMatrix prolongation_block(local_size, true_size);
    prolongation_block = 0.0;
    for (int i = 0; i < local_size; i++)
    {
      for (const auto &[position, value] : local_to_true[i])
      {
        prolongation_block(i, position) = value;
      }
    }
    mfem::DenseMatrix scratch(local_size, true_size);
    mfem::Mult(merged[element], prolongation_block, scratch);
    record.metric.SetSize(true_size);
    mfem::MultAtB(prolongation_block, scratch, record.metric);

    // Coarse element block, with its own true union and prolongation block.
    mfem::DofTransformation transformation;
    coarse_fespace.GetElementVDofs(element, coarse_dofs, transformation);
    std::vector<HYPRE_BigInt> coarse_ids;
    std::vector<char> coarse_flags;
    std::map<HYPRE_BigInt, int> coarse_positions;
    const auto coarse_position = [&](HYPRE_BigInt id, bool essential)
    {
      const auto found = coarse_positions.find(id);
      if (found != coarse_positions.end())
      {
        return found->second;
      }
      const int position = static_cast<int>(coarse_ids.size());
      coarse_positions.emplace(id, position);
      coarse_ids.push_back(id);
      coarse_flags.push_back(essential);
      return position;
    };
    const int coarse_local_size = coarse_dofs.Size();
    std::vector<std::vector<std::pair<int, double>>> coarse_to_true(coarse_local_size);
    std::vector<int> coarse_unsigned(coarse_local_size);
    for (int i = 0; i < coarse_local_size; i++)
    {
      const int unsigned_dof = coarse_dofs[i] >= 0 ? coarse_dofs[i] : -1 - coarse_dofs[i];
      coarse_unsigned[i] = unsigned_dof;
      for (const auto &[id, value] : coarse_rows[unsigned_dof])
      {
        coarse_to_true[i].push_back(
            {coarse_position(id, coarse_essential_vdofs[unsigned_dof] != 0), value});
      }
    }
    for (int guest = 0; guest < static_cast<int>(element_enrichment_guests[element].size());
         guest++)
    {
      const int enrichment_dof = element_enrichment_guests[element][guest];
      coarse_position(coarse_standard_global + numbering.local_to_true[enrichment_dof],
                      std::abs(enrichment_essential_local(enrichment_dof)) > 0.5);
    }
    record.coarse_dofs = coarse_ids;
    record.coarse_essential = coarse_flags;

    // Injection block in true coordinates:
    //   I_true = P_c^T T_local pinv(P_f)^T,  pinv(P_f) = (P_f^T P_f)^{-1} P_f^T,
    // exact because the injected coarse function is conforming. Enrichment rows are the
    // identity onto their fine enrichment true DOF.
    mfem::DenseMatrix injection_local(coarse_local_size, local_size);
    injection_local = 0.0;
    for (int i = 0; i < coarse_local_size; i++)
    {
      const auto &column = injection.columns[coarse_unsigned[i]];
      for (std::size_t entry = 0; entry < column.dofs.size(); entry++)
      {
        const auto position = dof_positions[element].find(column.dofs[entry]);
        if (position != dof_positions[element].end())
        {
          injection_local(i, position->second) = column.values[entry];
        }
      }
    }
    mfem::DenseMatrix gram(true_size);
    mfem::MultAtB(prolongation_block, prolongation_block, gram);
    mfem::DenseMatrixInverse gram_inverse(gram);
    mfem::DenseMatrix gram_inverse_matrix;
    gram_inverse.GetInverseMatrix(gram_inverse_matrix);
    mfem::DenseMatrix pinv_transpose(local_size, true_size);
    mfem::Mult(prolongation_block, gram_inverse_matrix, pinv_transpose);
    mfem::DenseMatrix injection_true_standard(coarse_local_size, true_size);
    mfem::Mult(injection_local, pinv_transpose, injection_true_standard);
    record.injection.SetSize(static_cast<int>(coarse_ids.size()), true_size);
    record.injection = 0.0;
    for (int i = 0; i < coarse_local_size; i++)
    {
      for (const auto &[position, value] : coarse_to_true[i])
      {
        for (int j = 0; j < true_size; j++)
        {
          record.injection(position, j) += value * injection_true_standard(i, j);
        }
      }
    }
    for (int guest = 0; guest < static_cast<int>(element_enrichment_guests[element].size());
         guest++)
    {
      const int enrichment_dof = element_enrichment_guests[element][guest];
      const HYPRE_BigInt guest_id =
          coarse_standard_global + numbering.local_to_true[enrichment_dof];
      const HYPRE_BigInt fine_id =
          fine_standard_global + numbering.local_to_true[enrichment_dof];
      const auto fine_position = true_positions.find(fine_id);
      MFEM_VERIFY(fine_position != true_positions.end(),
                  "Element enrichment guest is absent from its own element record!");
      record.injection(coarse_positions.at(guest_id), fine_position->second) = 1.0;
    }
  }
  return records;
}

std::vector<HierarchicalMaxwellDomainData::TrueElementRecord>
HierarchicalMaxwellDomainData::ExchangeHaloRecords(std::vector<TrueElementRecord> records,
                                                   int *ghost_count) const
{
  auto &mesh = space_op->GetMesh().Get();
  auto &fine_fespace = fine_nd_space->Get();
  if (ghost_count)
  {
    *ghost_count = 0;
  }
  if (Mpi::Size(space_op->GetComm()) == 1)
  {
    return records;
  }

  // Interface seed elements: any shared fine standard DOF or shared mesh vertex (the
  // latter covers vertex-anchored enrichment between ranks that share no ND DOF).
  std::vector<char> shared_ldof(fine_fespace.GetVSize(), false);
  const auto &group_ldof = fine_fespace.GroupComm().GroupLDofTable();
  for (int group = 1; group < group_ldof.Size(); group++)
  {
    const int *ldofs = group_ldof.GetRow(group);
    for (int i = 0; i < group_ldof.RowSize(group); i++)
    {
      shared_ldof[ldofs[i]] = true;
    }
  }
  std::vector<char> shared_vertex(mesh.GetNV(), false);
  for (int group = 1; group < mesh.GetNGroups(); group++)
  {
    for (int vertex = 0; vertex < mesh.GroupNVertices(group); vertex++)
    {
      shared_vertex[mesh.GroupVertex(group, vertex)] = true;
    }
  }

  // Three layers of DOF-connectivity outward from the interface. Owned complement modes
  // reach the entity's element ring, coarse guest columns one more ring, and their
  // support elements a third.
  std::map<HYPRE_BigInt, std::vector<int>> dof_elements;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    for (const auto dof : records[element].fine_dofs)
    {
      dof_elements[dof].push_back(element);
    }
  }
  std::vector<int> layer(mesh.GetNE(), 0);
  mfem::Array<int> vdofs, vertices;
  mfem::DofTransformation transformation;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    bool seed = false;
    fine_fespace.GetElementVDofs(element, vdofs, transformation);
    for (int dof : vdofs)
    {
      seed = seed || shared_ldof[dof >= 0 ? dof : -1 - dof];
    }
    mesh.GetElementVertices(element, vertices);
    for (int vertex : vertices)
    {
      seed = seed || shared_vertex[vertex];
    }
    layer[element] = seed ? 1 : 0;
  }
  for (int pass = 1; pass < 3; pass++)
  {
    std::vector<int> next(layer);
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      if (layer[element] != pass)
      {
        continue;
      }
      for (const auto dof : records[element].fine_dofs)
      {
        for (int neighbor : dof_elements[dof])
        {
          if (next[neighbor] == 0)
          {
            next[neighbor] = pass + 1;
          }
        }
      }
    }
    layer = next;
  }

  // Pack the halo once and send it to every mesh-topology neighbor rank.
  std::vector<double> payload;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    if (layer[element] == 0)
    {
      continue;
    }
    const auto &record = records[element];
    payload.push_back(static_cast<double>(record.element_id));
    payload.push_back(static_cast<double>(record.fine_dofs.size()));
    payload.push_back(static_cast<double>(record.coarse_dofs.size()));
    for (const auto dof : record.fine_dofs)
    {
      payload.push_back(static_cast<double>(dof));
    }
    for (const auto dof : record.coarse_dofs)
    {
      payload.push_back(static_cast<double>(dof));
    }
    for (const char flag : record.fine_essential)
    {
      payload.push_back(flag ? 1.0 : 0.0);
    }
    for (const char flag : record.coarse_essential)
    {
      payload.push_back(flag ? 1.0 : 0.0);
    }
    for (int i = 0; i < record.metric.Height(); i++)
    {
      for (int j = 0; j < record.metric.Width(); j++)
      {
        payload.push_back(record.metric(i, j));
      }
    }
    for (int i = 0; i < record.injection.Height(); i++)
    {
      for (int j = 0; j < record.injection.Width(); j++)
      {
        payload.push_back(record.injection(i, j));
      }
    }
  }

  const auto &group_topology = mesh.gtopo;
  std::vector<int> neighbors;
  for (int neighbor = 1; neighbor < group_topology.GetNumNeighbors(); neighbor++)
  {
    neighbors.push_back(group_topology.GetNeighborRank(neighbor));
  }
  std::vector<MPI_Request> requests;
  requests.reserve(2 * neighbors.size());
  const long long send_size = static_cast<long long>(payload.size());
  std::vector<long long> receive_sizes(neighbors.size(), 0);
  for (std::size_t i = 0; i < neighbors.size(); i++)
  {
    requests.emplace_back();
    MPI_Isend(&send_size, 1, MPI_LONG_LONG, neighbors[i], 71, space_op->GetComm(),
              &requests.back());
    requests.emplace_back();
    MPI_Irecv(&receive_sizes[i], 1, MPI_LONG_LONG, neighbors[i], 71, space_op->GetComm(),
              &requests.back());
  }
  MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
  requests.clear();
  std::vector<std::vector<double>> received(neighbors.size());
  for (std::size_t i = 0; i < neighbors.size(); i++)
  {
    received[i].resize(receive_sizes[i]);
    requests.emplace_back();
    MPI_Isend(payload.data(), static_cast<int>(payload.size()), MPI_DOUBLE, neighbors[i],
              72, space_op->GetComm(), &requests.back());
    requests.emplace_back();
    MPI_Irecv(received[i].data(), static_cast<int>(received[i].size()), MPI_DOUBLE,
              neighbors[i], 72, space_op->GetComm(), &requests.back());
  }
  MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);

  std::vector<TrueElementRecord> augmented = std::move(records);
  std::set<HYPRE_BigInt> known;
  for (const auto &record : records)
  {
    known.insert(record.element_id);
  }
  for (const auto &buffer : received)
  {
    std::size_t cursor = 0;
    while (cursor < buffer.size())
    {
      TrueElementRecord record;
      record.element_id = static_cast<HYPRE_BigInt>(buffer[cursor++]);
      const int fine_count = static_cast<int>(buffer[cursor++]);
      const int coarse_count = static_cast<int>(buffer[cursor++]);
      record.fine_dofs.resize(fine_count);
      record.coarse_dofs.resize(coarse_count);
      record.fine_essential.resize(fine_count);
      record.coarse_essential.resize(coarse_count);
      for (int i = 0; i < fine_count; i++)
      {
        record.fine_dofs[i] = static_cast<HYPRE_BigInt>(buffer[cursor++]);
      }
      for (int i = 0; i < coarse_count; i++)
      {
        record.coarse_dofs[i] = static_cast<HYPRE_BigInt>(buffer[cursor++]);
      }
      for (int i = 0; i < fine_count; i++)
      {
        record.fine_essential[i] = buffer[cursor++] > 0.5;
      }
      for (int i = 0; i < coarse_count; i++)
      {
        record.coarse_essential[i] = buffer[cursor++] > 0.5;
      }
      record.metric.SetSize(fine_count);
      for (int i = 0; i < fine_count; i++)
      {
        for (int j = 0; j < fine_count; j++)
        {
          record.metric(i, j) = buffer[cursor++];
        }
      }
      record.injection.SetSize(coarse_count, fine_count);
      for (int i = 0; i < coarse_count; i++)
      {
        for (int j = 0; j < fine_count; j++)
        {
          record.injection(i, j) = buffer[cursor++];
        }
      }
      if (known.insert(record.element_id).second)
      {
        augmented.push_back(std::move(record));
        if (ghost_count)
        {
          (*ghost_count)++;
        }
      }
    }
    MFEM_VERIFY(cursor == buffer.size(),
                "Hierarchical halo record buffer decoded inconsistently!");
  }
  return augmented;
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
    // Certified edge/face/interior lifting, now rank-count-agnostic: reduce the exact
    // local residual to combined true DOFs (which also eliminates essential equations),
    // then lift over the halo-augmented record set.
    Vector residual_true_real, residual_true_imag;
    LocalToTrue(residual.real, residual_true_real);
    LocalToTrue(residual.imag, residual_true_imag);
    return LiftTrueComplexResidualByEntityPatches(omega, residual_true_real,
                                                  residual_true_imag);
  }

  Vector residual_true_real, residual_true_imag;
  LocalToTrue(residual.real, residual_true_real);
  LocalToTrue(residual.imag, residual_true_imag);
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

fem::hierarchical::ComplexResidual
HierarchicalMaxwellDomainData::AssembleStreamingComplexResidual(
    std::complex<double> omega, const Vector &injected_real,
    const Vector &injected_imag) const
{
  const std::complex<double> mass_scale = -omega * omega;
  const std::complex<double> damping_scale = std::complex<double>(0.0, 1.0) * omega;
  const int combined_size = injected_real.Size();
  fem::hierarchical::ComplexResidual residual{Vector(combined_size), Vector(combined_size)};
  residual.real = 0.0;
  residual.imag = 0.0;
  mfem::DenseMatrix block_real, block_imag;
  Vector local_real, local_imag, action;
  const auto apply = [&](const std::vector<int> &dofs, const mfem::DenseMatrix &real,
                         const mfem::DenseMatrix &imag)
  {
    const int local_size = static_cast<int>(dofs.size());
    local_real.SetSize(local_size);
    local_imag.SetSize(local_size);
    action.SetSize(local_size);
    for (int i = 0; i < local_size; i++)
    {
      local_real(i) = injected_real(dofs[i]);
      local_imag(i) = injected_imag(dofs[i]);
    }
    // rr -= Ar xr - Ai xi; ri -= Ai xr + Ar xi.
    real.Mult(local_real, action);
    for (int i = 0; i < local_size; i++)
    {
      residual.real(dofs[i]) -= action(i);
    }
    real.Mult(local_imag, action);
    for (int i = 0; i < local_size; i++)
    {
      residual.imag(dofs[i]) -= action(i);
    }
    imag.Mult(local_real, action);
    for (int i = 0; i < local_size; i++)
    {
      residual.imag(dofs[i]) -= action(i);
    }
    imag.Mult(local_imag, action);
    for (int i = 0; i < local_size; i++)
    {
      residual.real(dofs[i]) += action(i);
    }
  };
  for (const auto &element : elements)
  {
    const int local_size = static_cast<int>(element.dofs.size());
    block_real.SetSize(local_size);
    block_imag.SetSize(local_size);
    block_real = element.curl_curl;
    block_real.Add(mass_scale.real() * element.electric_real -
                       mass_scale.imag() * element.electric_imag,
                   element.mass);
    block_imag = 0.0;
    block_imag.Add(mass_scale.imag() * element.electric_real +
                       mass_scale.real() * element.electric_imag,
                   element.mass);
    apply(element.dofs, block_real, block_imag);
  }
  const auto apply_boundary =
      [&](const std::vector<fem::singular::LocalNDBoundaryPatchMatrices> &patches,
          std::complex<double> scale)
  {
    for (const auto &patch : patches)
    {
      const int local_size = static_cast<int>(patch.dofs.size());
      block_real.SetSize(local_size);
      block_imag.SetSize(local_size);
      block_real = patch.mass;
      block_real *= scale.real();
      block_imag = patch.mass;
      block_imag *= scale.imag();
      apply(patch.dofs, block_real, block_imag);
    }
  };
  apply_boundary(boundary_stiffness, 1.0);
  apply_boundary(boundary_damping, damping_scale);
  apply_boundary(boundary_mass, mass_scale);
  return residual;
}

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
  auto residual = AssembleStreamingComplexResidual(omega, injected_real, injected_imag);
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
  auto residual =
      AssembleStreamingComplexResidual({omega, 0.0}, injected_real, injected_imag);
  const auto excitation = AssembleFineExcitation(excitation_idx, omega, coarse_rhs);
  residual.real += excitation.real;
  residual.imag += excitation.imag;
  return LiftLocalComplexResidual({omega, 0.0}, std::move(residual), patch_shape);
}

HierarchicalMaxwellDomainData::ParallelEstimate
HierarchicalMaxwellDomainData::LiftTrueComplexResidualByEntityPatches(
    std::complex<double> omega, const Vector &residual_true_real,
    const Vector &residual_true_imag) const
{
  auto &coarse_fespace = space_op->GetNDSpace().Get();
  auto &fine_fespace = fine_nd_space->Get();
  auto &mesh = space_op->GetMesh().Get();
  const auto &numbering = space_op->GetSingularParallelNumbering().nd;
  MPI_Comm comm = space_op->GetComm();
  const int fine_standard_true = fine_fespace.GetTrueVSize();
  const int enrichment_owned = static_cast<int>(numbering.owned_size);
  MFEM_VERIFY(residual_true_real.Size() == fine_standard_true + enrichment_owned &&
                  residual_true_imag.Size() == fine_standard_true + enrichment_owned,
              "Entity-patch lifting received an invalid combined true residual!");
  const HYPRE_BigInt fine_standard_global = fine_fespace.GlobalTrueVSize();
  const HYPRE_BigInt fine_true_offset = fine_fespace.GetMyTDofOffset();
  const HYPRE_BigInt enrichment_offset = numbering.owned_offset;

  auto records = ExchangeHaloRecords(BuildTrueMetricElementRecords(omega));
  const int own_elements = mesh.GetNE();
  const int record_count = static_cast<int>(records.size());

  // Dense halo DOF universes: every quantity below is a flat array over these indices,
  // which replaces the previous node-based map storage (whose allocation overhead
  // dominated time and memory at production scale).
  std::vector<HYPRE_BigInt> fine_universe, coarse_universe;
  for (const auto &record : records)
  {
    fine_universe.insert(fine_universe.end(), record.fine_dofs.begin(),
                         record.fine_dofs.end());
    coarse_universe.insert(coarse_universe.end(), record.coarse_dofs.begin(),
                           record.coarse_dofs.end());
  }
  std::sort(fine_universe.begin(), fine_universe.end());
  fine_universe.erase(std::unique(fine_universe.begin(), fine_universe.end()),
                      fine_universe.end());
  std::sort(coarse_universe.begin(), coarse_universe.end());
  coarse_universe.erase(std::unique(coarse_universe.begin(), coarse_universe.end()),
                        coarse_universe.end());
  const int nf = static_cast<int>(fine_universe.size());
  const int nc = static_cast<int>(coarse_universe.size());
  const auto fine_index = [&](HYPRE_BigInt id)
  {
    const auto it = std::lower_bound(fine_universe.begin(), fine_universe.end(), id);
    return (it != fine_universe.end() && *it == id)
               ? static_cast<int>(it - fine_universe.begin())
               : -1;
  };
  const auto coarse_index = [&](HYPRE_BigInt id)
  {
    const auto it = std::lower_bound(coarse_universe.begin(), coarse_universe.end(), id);
    return (it != coarse_universe.end() && *it == id)
               ? static_cast<int>(it - coarse_universe.begin())
               : -1;
  };

  // Record DOF lists in dense indices, essential flags, and DOF->record incidence (CSR).
  std::vector<std::vector<int>> record_fine(record_count), record_coarse(record_count);
  std::vector<char> fine_essential_flag(nf, false), coarse_essential_flag(nc, false);
  std::vector<int> incidence_offsets(nf + 1, 0);
  for (int r = 0; r < record_count; r++)
  {
    record_fine[r].resize(records[r].fine_dofs.size());
    for (std::size_t i = 0; i < records[r].fine_dofs.size(); i++)
    {
      const int f = fine_index(records[r].fine_dofs[i]);
      record_fine[r][i] = f;
      fine_essential_flag[f] = fine_essential_flag[f] || records[r].fine_essential[i];
      incidence_offsets[f + 1]++;
    }
    record_coarse[r].resize(records[r].coarse_dofs.size());
    for (std::size_t i = 0; i < records[r].coarse_dofs.size(); i++)
    {
      const int c = coarse_index(records[r].coarse_dofs[i]);
      record_coarse[r][i] = c;
      coarse_essential_flag[c] = coarse_essential_flag[c] || records[r].coarse_essential[i];
    }
  }
  for (int f = 0; f < nf; f++)
  {
    incidence_offsets[f + 1] += incidence_offsets[f];
  }
  std::vector<int> incidence_records(incidence_offsets[nf]);
  {
    std::vector<int> cursor(incidence_offsets.begin(), incidence_offsets.end() - 1);
    for (int r = 0; r < record_count; r++)
    {
      for (const int f : record_fine[r])
      {
        incidence_records[cursor[f]++] = r;
      }
    }
  }

  // Guest columns as one shared CSR over the coarse universe (entries deduplicated by
  // (coarse, fine); duplicated element contributions are consistent by construction).
  std::vector<int> column_offsets(nc + 1, 0);
  std::vector<int> column_fine;
  std::vector<double> column_value;
  {
    struct Entry
    {
      int coarse;
      int fine;
      double value;
    };
    std::vector<Entry> entries;
    for (int r = 0; r < record_count; r++)
    {
      const auto &record = records[r];
      for (int row = 0; row < record.injection.Height(); row++)
      {
        for (int j = 0; j < record.injection.Width(); j++)
        {
          const double value = record.injection(row, j);
          if (std::abs(value) > 1.0e-13)
          {
            entries.push_back({record_coarse[r][row], record_fine[r][j], value});
          }
        }
      }
    }
    std::sort(entries.begin(), entries.end(), [](const Entry &a, const Entry &b)
              { return a.coarse != b.coarse ? a.coarse < b.coarse : a.fine < b.fine; });
    entries.erase(std::unique(entries.begin(), entries.end(),
                              [](const Entry &a, const Entry &b)
                              { return a.coarse == b.coarse && a.fine == b.fine; }),
                  entries.end());
    column_fine.reserve(entries.size());
    column_value.reserve(entries.size());
    for (const auto &entry : entries)
    {
      column_offsets[entry.coarse + 1]++;
      column_fine.push_back(entry.fine);
      column_value.push_back(entry.value);
    }
    for (int c = 0; c < nc; c++)
    {
      column_offsets[c + 1] += column_offsets[c];
    }
  }
  // Only each record's metric is needed from here on (restricted congruence, sweeps, and
  // indicator energies); release the injection blocks and id/flag arrays.
  for (auto &record : records)
  {
    record.injection.SetSize(0, 0);
    record.fine_dofs.clear();
    record.fine_dofs.shrink_to_fit();
    record.coarse_dofs.clear();
    record.coarse_dofs.shrink_to_fit();
    record.fine_essential.clear();
    record.fine_essential.shrink_to_fit();
    record.coarse_essential.clear();
    record.coarse_essential.shrink_to_fit();
  }

  // Ownership and complete residual values on the halo.
  const auto owns_dof_id = [&](HYPRE_BigInt dof)
  {
    if (dof < fine_standard_global)
    {
      return dof >= fine_true_offset && dof < fine_true_offset + fine_standard_true;
    }
    const HYPRE_BigInt enrichment = dof - fine_standard_global;
    return enrichment >= enrichment_offset &&
           enrichment < enrichment_offset + enrichment_owned;
  };
  std::vector<char> fine_owned(nf, false);
  std::vector<int> fine_true_index(nf, -1);
  for (int f = 0; f < nf; f++)
  {
    const HYPRE_BigInt id = fine_universe[f];
    if (owns_dof_id(id))
    {
      fine_owned[f] = true;
      fine_true_index[f] =
          id < fine_standard_global
              ? static_cast<int>(id - fine_true_offset)
              : fine_standard_true +
                    static_cast<int>(id - fine_standard_global - enrichment_offset);
    }
  }

  // Neighbor merge of flat (id, values...) streams; unknown ids are beyond this rank's
  // halo and provably cannot touch its owned DOFs, so they are skipped.
  const auto &group_topology = mesh.gtopo;
  std::vector<int> neighbors;
  for (int neighbor = 1; neighbor < group_topology.GetNumNeighbors(); neighbor++)
  {
    neighbors.push_back(group_topology.GetNeighborRank(neighbor));
  }
  const auto exchange_sum = [&](const std::vector<HYPRE_BigInt> &universe,
                                const auto &index_of, std::vector<double> &values,
                                int width)
  {
    if (neighbors.empty())
    {
      return;
    }
    std::vector<double> payload;
    const int n = static_cast<int>(universe.size());
    for (int i = 0; i < n; i++)
    {
      bool nonzero = false;
      for (int k = 0; k < width; k++)
      {
        nonzero = nonzero || values[static_cast<std::size_t>(i) * width + k] != 0.0;
      }
      if (nonzero)
      {
        payload.push_back(static_cast<double>(universe[i]));
        for (int k = 0; k < width; k++)
        {
          payload.push_back(values[static_cast<std::size_t>(i) * width + k]);
        }
      }
    }
    std::vector<MPI_Request> requests;
    const long long send_size = static_cast<long long>(payload.size());
    std::vector<long long> receive_sizes(neighbors.size(), 0);
    for (std::size_t i = 0; i < neighbors.size(); i++)
    {
      requests.emplace_back();
      MPI_Isend(&send_size, 1, MPI_LONG_LONG, neighbors[i], 73, comm, &requests.back());
      requests.emplace_back();
      MPI_Irecv(&receive_sizes[i], 1, MPI_LONG_LONG, neighbors[i], 73, comm,
                &requests.back());
    }
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
    requests.clear();
    std::vector<std::vector<double>> received(neighbors.size());
    for (std::size_t i = 0; i < neighbors.size(); i++)
    {
      received[i].resize(receive_sizes[i]);
      requests.emplace_back();
      MPI_Isend(payload.data(), static_cast<int>(payload.size()), MPI_DOUBLE, neighbors[i],
                74, comm, &requests.back());
      requests.emplace_back();
      MPI_Irecv(received[i].data(), static_cast<int>(received[i].size()), MPI_DOUBLE,
                neighbors[i], 74, comm, &requests.back());
    }
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
    for (const auto &buffer : received)
    {
      for (std::size_t cursor = 0; cursor < buffer.size();
           cursor += static_cast<std::size_t>(width) + 1)
      {
        const int i = index_of(static_cast<HYPRE_BigInt>(buffer[cursor]));
        if (i < 0)
        {
          continue;
        }
        for (int k = 0; k < width; k++)
        {
          values[static_cast<std::size_t>(i) * width + k] += buffer[cursor + 1 + k];
        }
      }
    }
  };

  std::vector<double> residual_values(2 * static_cast<std::size_t>(nf), 0.0);
  for (int f = 0; f < nf; f++)
  {
    if (fine_owned[f])
    {
      residual_values[2 * static_cast<std::size_t>(f)] =
          residual_true_real(fine_true_index[f]);
      residual_values[2 * static_cast<std::size_t>(f) + 1] =
          residual_true_imag(fine_true_index[f]);
    }
  }
  exchange_sum(fine_universe, fine_index, residual_values, 2);

  // Reverse ldof -> group map classifies entity ownership.
  const auto build_group_map = [](mfem::ParFiniteElementSpace &fespace)
  {
    std::vector<int> group(fespace.GetVSize(), 0);
    const auto &table = fespace.GroupComm().GroupLDofTable();
    for (int g = 1; g < table.Size(); g++)
    {
      const int *ldofs = table.GetRow(g);
      for (int i = 0; i < table.RowSize(g); i++)
      {
        group[ldofs[i]] = g;
      }
    }
    return group;
  };
  const auto fine_group = build_group_map(fine_fespace);
  const auto build_ldof_rows = [](mfem::ParFiniteElementSpace &fespace)
  {
    std::vector<std::vector<std::pair<HYPRE_BigInt, double>>> rows(fespace.GetVSize());
    mfem::SparseMatrix diag, offd;
    HYPRE_BigInt *colmap = nullptr;
    fespace.Dof_TrueDof_Matrix()->GetDiag(diag);
    fespace.Dof_TrueDof_Matrix()->GetOffd(offd, colmap);
    const HYPRE_BigInt true_offset = fespace.GetMyTDofOffset();
    for (int ldof = 0; ldof < fespace.GetVSize(); ldof++)
    {
      for (int entry = diag.GetI()[ldof]; entry < diag.GetI()[ldof + 1]; entry++)
      {
        if (std::abs(diag.GetData()[entry]) > 1.0e-14)
        {
          rows[ldof].push_back({true_offset + diag.GetJ()[entry], diag.GetData()[entry]});
        }
      }
      for (int entry = offd.GetI()[ldof]; entry < offd.GetI()[ldof + 1]; entry++)
      {
        if (std::abs(offd.GetData()[entry]) > 1.0e-14)
        {
          rows[ldof].push_back({colmap[offd.GetJ()[entry]], offd.GetData()[entry]});
        }
      }
    }
    return rows;
  };
  const auto fine_rows = build_ldof_rows(fine_fespace);
  const auto coarse_rows = build_ldof_rows(coarse_fespace);

  // Patch storage: owned complement columns live in one pooled sparse array; guest
  // columns are referenced by coarse universe index into the shared column CSR instead of
  // being copied per patch.
  struct Patch
  {
    int owned = 0;
    std::vector<int> owned_begin;  // owned + 1 offsets into the pools.
    std::vector<int> guests;       // coarse universe indices.
    // Packed lower-triangular Cholesky factor of the restricted SPD metric; the dense
    // restricted matrix itself is never retained.
    std::vector<double> factor;
    std::vector<int> support;  // record indices, sorted unique.
  };
  const auto packed_cholesky_factor =
      [](mfem::DenseMatrix &matrix, std::vector<double> &factor)
  {
    const int n = matrix.Height();
    factor.assign(static_cast<std::size_t>(n) * (n + 1) / 2, 0.0);
    for (int j = 0; j < n; j++)
    {
      double diagonal = matrix(j, j);
      const std::size_t row_j = static_cast<std::size_t>(j) * (j + 1) / 2;
      for (int k = 0; k < j; k++)
      {
        diagonal -= factor[row_j + k] * factor[row_j + k];
      }
      MFEM_VERIFY(diagonal > 0.0,
                  "Entity-patch restricted metric lost positive definiteness!");
      const double pivot = std::sqrt(diagonal);
      factor[row_j + j] = pivot;
      for (int i = j + 1; i < n; i++)
      {
        const std::size_t row_i = static_cast<std::size_t>(i) * (i + 1) / 2;
        double value = matrix(i, j);
        for (int k = 0; k < j; k++)
        {
          value -= factor[row_i + k] * factor[row_j + k];
        }
        factor[row_i + j] = value / pivot;
      }
    }
  };
  const auto packed_cholesky_solve =
      [](const std::vector<double> &factor, int n, Vector &rhs)
  {
    for (int i = 0; i < n; i++)
    {
      const std::size_t row_i = static_cast<std::size_t>(i) * (i + 1) / 2;
      double value = rhs(i);
      for (int k = 0; k < i; k++)
      {
        value -= factor[row_i + k] * rhs(k);
      }
      rhs(i) = value / factor[row_i + i];
    }
    for (int i = n - 1; i >= 0; i--)
    {
      const std::size_t row_i = static_cast<std::size_t>(i) * (i + 1) / 2;
      double value = rhs(i);
      for (int k = i + 1; k < n; k++)
      {
        value -= factor[static_cast<std::size_t>(k) * (k + 1) / 2 + i] * rhs(k);
      }
      rhs(i) = value / factor[row_i + i];
    }
  };
  std::vector<Patch> patches;
  std::vector<int> pool_dof;
  std::vector<double> pool_value;
  int owned_modes = 0;

  // Generation-stamped scratch arrays for patch-local gathers (no per-patch allocation).
  std::vector<int> fine_stamp(nf, -1);
  std::vector<int> fine_slot(nf, -1);
  int generation = 0;

  const auto add_patch =
      [&](const std::vector<int> &fine_entity, const std::vector<int> &coarse_entity)
  {
    std::vector<int> rows;
    for (const int f : fine_entity)
    {
      if (f >= 0 && !fine_essential_flag[f])
      {
        rows.push_back(f);
      }
    }
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    if (rows.empty())
    {
      return;
    }
    const int row_count = static_cast<int>(rows.size());
    generation++;
    for (int i = 0; i < row_count; i++)
    {
      fine_stamp[rows[i]] = generation;
      fine_slot[rows[i]] = i;
    }
    // Range of the injected coarse entity trace restricted to the entity rows.
    std::vector<mfem::Vector> range;
    for (const int c : coarse_entity)
    {
      if (c < 0 || coarse_essential_flag[c])
      {
        continue;
      }
      mfem::Vector vector(row_count);
      vector = 0.0;
      for (int entry = column_offsets[c]; entry < column_offsets[c + 1]; entry++)
      {
        const int f = column_fine[entry];
        if (fine_stamp[f] == generation)
        {
          vector(fine_slot[f]) = column_value[entry];
        }
      }
      for (int repeat = 0; repeat < 2; repeat++)
      {
        for (const auto &basis : range)
        {
          vector.Add(-mfem::InnerProduct(vector, basis), basis);
        }
      }
      const double norm = vector.Norml2();
      if (norm > 1.0e-12)
      {
        vector /= norm;
        range.push_back(vector);
      }
    }
    // Deterministic complement in sorted-row coordinates.
    std::vector<mfem::Vector> owned_columns;
    for (int direction = 0; direction < row_count; direction++)
    {
      mfem::Vector vector(row_count);
      vector = 0.0;
      vector(direction) = 1.0;
      for (int repeat = 0; repeat < 2; repeat++)
      {
        for (const auto &basis : range)
        {
          vector.Add(-mfem::InnerProduct(vector, basis), basis);
        }
        for (const auto &basis : owned_columns)
        {
          vector.Add(-mfem::InnerProduct(vector, basis), basis);
        }
      }
      const double norm = vector.Norml2();
      if (norm > 1.0e-10)
      {
        vector /= norm;
        owned_columns.push_back(vector);
      }
    }
    if (owned_columns.empty())
    {
      return;
    }
    Patch patch;
    patch.owned = static_cast<int>(owned_columns.size());
    owned_modes += patch.owned;
    patch.owned_begin.push_back(static_cast<int>(pool_dof.size()));
    for (const auto &column : owned_columns)
    {
      for (int i = 0; i < row_count; i++)
      {
        if (std::abs(column(i)) > 1.0e-13)
        {
          pool_dof.push_back(rows[i]);
          pool_value.push_back(column(i));
        }
      }
      patch.owned_begin.push_back(static_cast<int>(pool_dof.size()));
    }

    // Guests: non-essential coarse DOFs of every record incident to the entity.
    {
      std::vector<int> owner_records;
      for (const int f : rows)
      {
        for (int entry = incidence_offsets[f]; entry < incidence_offsets[f + 1]; entry++)
        {
          owner_records.push_back(incidence_records[entry]);
        }
      }
      std::sort(owner_records.begin(), owner_records.end());
      owner_records.erase(std::unique(owner_records.begin(), owner_records.end()),
                          owner_records.end());
      std::vector<int> guests;
      for (const int r : owner_records)
      {
        for (const int c : record_coarse[r])
        {
          if (!coarse_essential_flag[c])
          {
            guests.push_back(c);
          }
        }
      }
      std::sort(guests.begin(), guests.end());
      guests.erase(std::unique(guests.begin(), guests.end()), guests.end());
      patch.guests = std::move(guests);
    }

    // Support union over every basis column DOF, through the incidence CSR.
    {
      std::vector<int> support;
      const auto add_dof_support = [&](int f)
      {
        MFEM_VERIFY(!fine_essential_flag[f],
                    "Entity-patch basis columns must avoid essential fine equations!");
        for (int entry = incidence_offsets[f]; entry < incidence_offsets[f + 1]; entry++)
        {
          support.push_back(incidence_records[entry]);
        }
      };
      for (int entry = patch.owned_begin.front(); entry < patch.owned_begin.back(); entry++)
      {
        add_dof_support(pool_dof[entry]);
      }
      for (const int c : patch.guests)
      {
        for (int entry = column_offsets[c]; entry < column_offsets[c + 1]; entry++)
        {
          add_dof_support(column_fine[entry]);
        }
      }
      std::sort(support.begin(), support.end());
      support.erase(std::unique(support.begin(), support.end()), support.end());
      patch.support = std::move(support);
    }

    // Restricted metric over the support union, assembled into a transient dense block
    // and immediately factored. Column membership per record DOF is gathered through the
    // generation-stamped scratch.
    const int dimension = patch.owned + static_cast<int>(patch.guests.size());
    mfem::DenseMatrix restricted(dimension);
    restricted = 0.0;
    std::vector<std::vector<std::pair<int, double>>> dof_columns_scratch;
    for (const int r : patch.support)
    {
      const auto &fine = record_fine[r];
      const int local_size = static_cast<int>(fine.size());
      dof_columns_scratch.assign(local_size, {});
      bool touched = false;
      generation++;
      for (int i = 0; i < local_size; i++)
      {
        fine_stamp[fine[i]] = generation;
        fine_slot[fine[i]] = i;
      }
      for (int b = 0; b < patch.owned; b++)
      {
        for (int entry = patch.owned_begin[b]; entry < patch.owned_begin[b + 1]; entry++)
        {
          if (fine_stamp[pool_dof[entry]] == generation)
          {
            dof_columns_scratch[fine_slot[pool_dof[entry]]].push_back(
                {b, pool_value[entry]});
            touched = true;
          }
        }
      }
      for (std::size_t g = 0; g < patch.guests.size(); g++)
      {
        const int c = patch.guests[g];
        for (int entry = column_offsets[c]; entry < column_offsets[c + 1]; entry++)
        {
          if (fine_stamp[column_fine[entry]] == generation)
          {
            dof_columns_scratch[fine_slot[column_fine[entry]]].push_back(
                {patch.owned + static_cast<int>(g), column_value[entry]});
            touched = true;
          }
        }
      }
      if (!touched)
      {
        continue;
      }
      const auto &metric = records[r].metric;
      for (int i = 0; i < local_size; i++)
      {
        if (dof_columns_scratch[i].empty())
        {
          continue;
        }
        for (int j = 0; j < local_size; j++)
        {
          const double value = metric(i, j);
          if (value == 0.0 || dof_columns_scratch[j].empty())
          {
            continue;
          }
          for (const auto &[bi, vi] : dof_columns_scratch[i])
          {
            for (const auto &[bj, vj] : dof_columns_scratch[j])
            {
              restricted(bi, bj) += vi * value * vj;
            }
          }
        }
      }
    }
    packed_cholesky_factor(restricted, patch.factor);
    patches.push_back(std::move(patch));
  };

  // Owned entities: shared entities belong to their DOF-group master; local entities to
  // this rank. Entity DOF ids come from the interior accessors through the prolongation
  // rows (which may couple pairs of true DOFs on shared tetrahedral faces).
  mfem::Array<int> fine_entity_dofs, coarse_entity_dofs;
  const auto entity_indices =
      [&](mfem::Array<int> &dofs, const auto &rows, const auto &index_of)
  {
    std::vector<int> indices;
    mfem::FiniteElementSpace::AdjustVDofs(dofs);
    for (int dof : dofs)
    {
      for (const auto &[id, value] : rows[dof])
      {
        (void)value;
        indices.push_back(index_of(id));
      }
    }
    return indices;
  };
  const auto owns_entity = [&](const mfem::Array<int> &dofs)
  {
    for (int dof : dofs)
    {
      const int group = fine_group[dof];
      if (group != 0)
      {
        return static_cast<bool>(
            fine_fespace.GroupComm().GetGroupTopology().IAmMaster(group));
      }
    }
    return true;
  };
  for (int edge = 0; edge < mesh.GetNEdges(); edge++)
  {
    fine_fespace.GetEdgeInteriorDofs(edge, fine_entity_dofs);
    if (fine_entity_dofs.Size() == 0 || !owns_entity(fine_entity_dofs))
    {
      continue;
    }
    coarse_fespace.GetEdgeInteriorDofs(edge, coarse_entity_dofs);
    add_patch(entity_indices(fine_entity_dofs, fine_rows, fine_index),
              entity_indices(coarse_entity_dofs, coarse_rows, coarse_index));
  }
  if (mesh.Dimension() == 3)
  {
    for (int face = 0; face < mesh.GetNFaces(); face++)
    {
      fine_fespace.GetFaceInteriorDofs(face, fine_entity_dofs);
      if (fine_entity_dofs.Size() == 0 || !owns_entity(fine_entity_dofs))
      {
        continue;
      }
      coarse_fespace.GetFaceInteriorDofs(face, coarse_entity_dofs);
      add_patch(entity_indices(fine_entity_dofs, fine_rows, fine_index),
                entity_indices(coarse_entity_dofs, coarse_rows, coarse_index));
    }
  }
  for (int element = 0; element < own_elements; element++)
  {
    fine_fespace.GetElementInteriorDofs(element, fine_entity_dofs);
    if (fine_entity_dofs.Size() == 0)
    {
      continue;
    }
    coarse_fespace.GetElementInteriorDofs(element, coarse_entity_dofs);
    add_patch(entity_indices(fine_entity_dofs, fine_rows, fine_index),
              entity_indices(coarse_entity_dofs, coarse_rows, coarse_index));
  }

  // Global ownership identity: every conforming complement direction is owned exactly
  // once across all ranks.
  {
    long long counts[3] = {owned_modes, fine_standard_true, 0};
    mfem::Array<int> essential_true;
    const int boundary_attribute_maximum =
        mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    const auto marker =
        mesh::AttrToMarker(boundary_attribute_maximum, space_op->GetNDDbcAttributes());
    fine_fespace.GetEssentialTrueDofs(marker, essential_true);
    counts[1] -= essential_true.Size();
    mfem::Array<int> coarse_essential_true;
    coarse_fespace.GetEssentialTrueDofs(marker, coarse_essential_true);
    counts[2] = coarse_fespace.GetTrueVSize() - coarse_essential_true.Size();
    Mpi::GlobalSum(3, counts, comm);
    MFEM_VERIFY(counts[0] == counts[1] - counts[2],
                "Entity-patch ownership must cover every conforming complement "
                "direction exactly once across all ranks!");
  }

  // Bounded undamped additive-Schwarz sweeps applied identically to both components.
  std::vector<double> correction(2 * static_cast<std::size_t>(nf), 0.0);
  std::vector<double> current(2 * static_cast<std::size_t>(nf), 0.0);
  std::vector<double> delta(4 * static_cast<std::size_t>(nf), 0.0);
  std::vector<double> guest_data(4 * static_cast<std::size_t>(nc), 0.0);
  Vector record_real, record_imag, record_action;
  constexpr int sweeps = 4;
  for (int sweep = 0; sweep < sweeps; sweep++)
  {
    std::copy(residual_values.begin(), residual_values.end(), current.begin());
    if (sweep > 0)
    {
      for (int r = 0; r < record_count; r++)
      {
        const auto &fine = record_fine[r];
        const int local_size = static_cast<int>(fine.size());
        record_real.SetSize(local_size);
        record_imag.SetSize(local_size);
        record_action.SetSize(local_size);
        for (int i = 0; i < local_size; i++)
        {
          record_real(i) = correction[2 * static_cast<std::size_t>(fine[i])];
          record_imag(i) = correction[2 * static_cast<std::size_t>(fine[i]) + 1];
        }
        records[r].metric.Mult(record_real, record_action);
        for (int i = 0; i < local_size; i++)
        {
          current[2 * static_cast<std::size_t>(fine[i])] -= record_action(i);
        }
        records[r].metric.Mult(record_imag, record_action);
        for (int i = 0; i < local_size; i++)
        {
          current[2 * static_cast<std::size_t>(fine[i]) + 1] -= record_action(i);
        }
      }
    }
    for (int f = 0; f < nf; f++)
    {
      if (fine_essential_flag[f])
      {
        current[2 * static_cast<std::size_t>(f)] = 0.0;
        current[2 * static_cast<std::size_t>(f) + 1] = 0.0;
      }
    }

    std::fill(delta.begin(), delta.end(), 0.0);
    std::fill(guest_data.begin(), guest_data.end(), 0.0);
    for (std::size_t p = 0; p < patches.size(); p++)
    {
      const auto &patch = patches[p];
      const int dimension = patch.owned + static_cast<int>(patch.guests.size());
      Vector rhs_real(dimension), rhs_imag(dimension);
      for (int b = 0; b < patch.owned; b++)
      {
        double sum_real = 0.0, sum_imag = 0.0;
        for (int entry = patch.owned_begin[b]; entry < patch.owned_begin[b + 1]; entry++)
        {
          sum_real +=
              pool_value[entry] * current[2 * static_cast<std::size_t>(pool_dof[entry])];
          sum_imag += pool_value[entry] *
                      current[2 * static_cast<std::size_t>(pool_dof[entry]) + 1];
        }
        rhs_real(b) = sum_real;
        rhs_imag(b) = sum_imag;
      }
      for (std::size_t g = 0; g < patch.guests.size(); g++)
      {
        const int c = patch.guests[g];
        double sum_real = 0.0, sum_imag = 0.0;
        for (int entry = column_offsets[c]; entry < column_offsets[c + 1]; entry++)
        {
          sum_real += column_value[entry] *
                      current[2 * static_cast<std::size_t>(column_fine[entry])];
          sum_imag += column_value[entry] *
                      current[2 * static_cast<std::size_t>(column_fine[entry]) + 1];
        }
        rhs_real(patch.owned + static_cast<int>(g)) = sum_real;
        rhs_imag(patch.owned + static_cast<int>(g)) = sum_imag;
      }
      Vector &c_real = rhs_real;
      Vector &c_imag = rhs_imag;
      packed_cholesky_solve(patch.factor, dimension, c_real);
      packed_cholesky_solve(patch.factor, dimension, c_imag);
      for (int b = 0; b < patch.owned; b++)
      {
        for (int entry = patch.owned_begin[b]; entry < patch.owned_begin[b + 1]; entry++)
        {
          delta[4 * static_cast<std::size_t>(pool_dof[entry])] +=
              c_real(b) * pool_value[entry];
          delta[4 * static_cast<std::size_t>(pool_dof[entry]) + 1] +=
              c_imag(b) * pool_value[entry];
        }
      }
      for (std::size_t g = 0; g < patch.guests.size(); g++)
      {
        guest_data[4 * static_cast<std::size_t>(patch.guests[g])] +=
            c_real(patch.owned + static_cast<int>(g));
        guest_data[4 * static_cast<std::size_t>(patch.guests[g]) + 1] +=
            c_imag(patch.owned + static_cast<int>(g));
        guest_data[4 * static_cast<std::size_t>(patch.guests[g]) + 2] += 1.0;
      }
    }
    exchange_sum(coarse_universe, coarse_index, guest_data, 4);
    for (int c = 0; c < nc; c++)
    {
      const double count = guest_data[4 * static_cast<std::size_t>(c) + 2];
      if (count <= 0.0)
      {
        continue;
      }
      const double average_real = guest_data[4 * static_cast<std::size_t>(c)] / count;
      const double average_imag = guest_data[4 * static_cast<std::size_t>(c) + 1] / count;
      for (int entry = column_offsets[c]; entry < column_offsets[c + 1]; entry++)
      {
        const int f = column_fine[entry];
        // Guest averages scatter only at owned DOFs so the merge counts each entry once.
        if (!fine_essential_flag[f] && fine_owned[f])
        {
          delta[4 * static_cast<std::size_t>(f) + 2] += average_real * column_value[entry];
          delta[4 * static_cast<std::size_t>(f) + 3] += average_imag * column_value[entry];
        }
      }
    }
    exchange_sum(fine_universe, fine_index, delta, 4);
    for (int f = 0; f < nf; f++)
    {
      correction[2 * static_cast<std::size_t>(f)] +=
          delta[4 * static_cast<std::size_t>(f)] +
          delta[4 * static_cast<std::size_t>(f) + 2];
      correction[2 * static_cast<std::size_t>(f) + 1] +=
          delta[4 * static_cast<std::size_t>(f) + 1] +
          delta[4 * static_cast<std::size_t>(f) + 3];
    }
  }

  // Element indicator energies over this rank's own elements.
  ParallelEstimate estimate;
  estimate.indicator_energy.SetSize(own_elements);
  estimate.indicator_energy = 0.0;
  for (int element = 0; element < own_elements; element++)
  {
    const auto &fine = record_fine[element];
    const int local_size = static_cast<int>(fine.size());
    record_real.SetSize(local_size);
    record_imag.SetSize(local_size);
    record_action.SetSize(local_size);
    for (int i = 0; i < local_size; i++)
    {
      record_real(i) = correction[2 * static_cast<std::size_t>(fine[i])];
      record_imag(i) = correction[2 * static_cast<std::size_t>(fine[i]) + 1];
    }
    records[element].metric.Mult(record_real, record_action);
    estimate.indicator_energy(element) += record_real * record_action;
    records[element].metric.Mult(record_imag, record_action);
    estimate.indicator_energy(element) += record_imag * record_action;
  }
  estimate.total_energy = estimate.indicator_energy.Sum();
  Mpi::GlobalSum(1, &estimate.total_energy, comm);
  return estimate;
}

}  // namespace palace
