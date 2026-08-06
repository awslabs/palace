// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hierarchicalmaxwellestimator.hpp"

#include <algorithm>
#include <cmath>

#include "fem/fespace.hpp"
#include "fem/singularassembly.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"

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

HierarchicalMaxwellDomainData::HierarchicalMaxwellDomainData(SpaceOperator &space_op)
{
  MFEM_VERIFY(space_op.HasSingularEnrichment(),
              "Hierarchical Maxwell domain data requires singular enrichment!");
  MFEM_VERIFY(Mpi::Size(space_op.GetComm()) == 1,
              "The first hierarchical Maxwell domain adapter is serial-only; parallel "
              "shared-entity ownership must be enabled before MPI production use!");
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
  MFEM_VERIFY(coarse_nd.GetVSize() == coarse_nd.GetTrueVSize() &&
                  fine_nd_space->GetVSize() == fine_nd_space->GetTrueVSize(),
              "The first hierarchical Maxwell domain adapter requires a conforming "
              "rank-one mesh with identity local-to-true constraints!");
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

  const int enrichment_size =
      dimension == 3
          ? static_cast<int>(space_op.GetSingularDofTopology()->nd_dofs.size())
          : static_cast<int>(space_op.GetTriangleSingularDofTopology()->nd_dofs.size());
  const int fine_standard_size = fine_nd_space->GetVSize();
  const int coarse_standard_size = coarse_nd.GetVSize();
  MFEM_VERIFY(static_cast<int>(space_op.GetSingularParallelNumbering().nd.owned_size) ==
                  enrichment_size,
              "Serial hierarchical Maxwell enrichment numbering is inconsistent!");

  // Essential masks in the local combined layouts. On one rank, conforming true and local
  // standard DOF numbers agree; singular true identities are order-independent.
  coarse_essential.assign(coarse_standard_size + enrichment_size, false);
  for (int dof : space_op.GetCombinedNDDbcTDofList())
  {
    MFEM_VERIFY(dof >= 0 && dof < static_cast<int>(coarse_essential.size()),
                "Invalid coarse combined essential DOF!");
    coarse_essential[dof] = true;
  }
  fine_essential.assign(fine_standard_size + enrichment_size, false);
  mfem::Array<int> fine_standard_essential;
  fine_nd_space->Get().GetEssentialTrueDofs(space_op.GetNDDbcAttributes(),
                                            fine_standard_essential);
  for (int dof : fine_standard_essential)
  {
    fine_essential[dof] = true;
  }
  for (int dof : space_op.GetCombinedNDDbcTDofList())
  {
    if (dof >= coarse_standard_size)
    {
      fine_essential[fine_standard_size + dof - coarse_standard_size] = true;
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
}

HierarchicalMaxwellDomainData::~HierarchicalMaxwellDomainData() = default;

int HierarchicalMaxwellDomainData::GetFineStandardSize() const
{
  return fine_nd_space->GetVSize();
}

int HierarchicalMaxwellDomainData::GetEnrichmentSize() const
{
  return static_cast<int>(fine_essential.size()) - GetFineStandardSize();
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

}  // namespace palace
