// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularsolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <nlohmann/json.hpp>

#include "linalg/operator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

std::vector<fem::singular::TriangleMaterial>
GetSingularTriangleMaterials(const IoData &iodata)
{
  std::map<int, double> permittivity;
  for (const auto &material : iodata.domains.materials)
  {
    const double epsilon = material.epsilon_r.s[0];
    const bool isotropic =
        std::all_of(material.epsilon_r.s.begin() + 1, material.epsilon_r.s.end(),
                    [epsilon](double value) { return value == epsilon; });
    if (!isotropic)
    {
      // The corner extractor diagnoses an anisotropic material only when it occurs in
      // an enriched fan. Unrelated anisotropic domains must not disable an otherwise
      // isotropic electrostatic corner.
      continue;
    }
    MFEM_VERIFY(std::isfinite(epsilon) && epsilon > 0.0,
                "Triangular singular feature extraction requires positive permittivity "
                "in every isotropic material sector!");
    for (int attribute : material.attributes)
    {
      MFEM_VERIFY(attribute > 0 && permittivity.emplace(attribute, epsilon).second,
                  "Triangular singular feature extraction found an invalid or duplicate "
                  "material domain attribute "
                      << attribute << "!");
    }
  }
  std::vector<fem::singular::TriangleMaterial> result;
  result.reserve(permittivity.size());
  for (const auto &[attribute, epsilon] : permittivity)
  {
    result.push_back({attribute, epsilon});
  }
  return result;
}

namespace
{

std::map<int, double> GetIsotropicLossTangents(const IoData &iodata)
{
  std::map<int, double> result;
  for (const auto &material : iodata.domains.materials)
  {
    const double loss_tangent = material.tandelta.s[0];
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() *
                             std::max(1.0, std::abs(loss_tangent));
    MFEM_VERIFY(std::isfinite(loss_tangent) &&
                    std::all_of(material.tandelta.s.begin() + 1, material.tandelta.s.end(),
                                [loss_tangent, tolerance](double value)
                                {
                                  return std::isfinite(value) &&
                                         std::abs(value - loss_tangent) <= tolerance;
                                }),
                "Singular bulk dielectric loss requires isotropic loss tangent!");
    for (int attribute : material.attributes)
    {
      MFEM_VERIFY(attribute > 0 && result.emplace(attribute, loss_tangent).second,
                  "Singular bulk dielectric loss found an invalid or duplicate material "
                  "domain attribute "
                      << attribute << "!");
    }
  }
  return result;
}

void VerifyCommonWedgeLossTangent(const std::map<int, double> &loss_tangents,
                                  const std::vector<int> &attributes,
                                  const char *description)
{
  MFEM_VERIFY(!attributes.empty(), description << " has no material sectors!");
  const auto first = loss_tangents.find(attributes.front());
  MFEM_VERIFY(first != loss_tangents.end(),
              description << " references an unknown material domain attribute "
                          << attributes.front() << "!");
  const double reference = first->second;
  const double tolerance =
      128.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(reference));
  for (int attribute : attributes)
  {
    const auto current = loss_tangents.find(attribute);
    MFEM_VERIFY(current != loss_tangents.end(),
                description << " references an unknown material domain attribute "
                            << attribute << "!");
    MFEM_VERIFY(
        std::abs(current->second - reference) <= tolerance,
        description
            << " has unequal material loss tangents. The exact transmission-wedge "
               "exponent is then complex, but the current singular basis supports only "
               "real exponents. Use a common loss tangent in all sectors or disable "
               "enrichment at this finite-metal feature!");
  }
}

}  // namespace

void ValidateSingularLossTangents(const IoData &iodata, const mfem::Mesh &mesh,
                                  const fem::singular::FeatureTopology &features)
{
  const auto loss_tangents = GetIsotropicLossTangents(iodata);
  for (std::size_t segment = 0; segment < features.segments.size(); segment++)
  {
    if (features.segments[segment].type !=
        fem::singular::FeatureSegmentType::TRANSMISSION_WEDGE)
    {
      continue;
    }
    std::vector<int> attributes;
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      const auto &incidence = features.elements.at(element);
      if (std::any_of(incidence.edges.begin(), incidence.edges.end(),
                      [segment](const auto &edge) { return edge.segment == segment; }))
      {
        attributes.push_back(mesh.GetAttribute(element));
      }
    }
    std::sort(attributes.begin(), attributes.end());
    attributes.erase(std::unique(attributes.begin(), attributes.end()), attributes.end());
    VerifyCommonWedgeLossTangent(loss_tangents, attributes,
                                 "A three-dimensional singular transmission wedge");
  }
}

void ValidateSingularLossTangents(const IoData &iodata,
                                  const fem::singular::TriangleFeatureTopology &features)
{
  const auto loss_tangents = GetIsotropicLossTangents(iodata);
  for (const auto &vertex : features.vertices)
  {
    if (vertex.type != fem::singular::FeatureVertexType::CORNER)
    {
      continue;
    }
    std::vector<int> attributes;
    attributes.reserve(vertex.sectors.size());
    for (const auto &sector : vertex.sectors)
    {
      attributes.push_back(sector.domain_attribute);
    }
    std::sort(attributes.begin(), attributes.end());
    attributes.erase(std::unique(attributes.begin(), attributes.end()), attributes.end());
    VerifyCommonWedgeLossTangent(loss_tangents, attributes,
                                 "A two-dimensional singular transmission wedge");
  }
}

nlohmann::json GetSingularSurfaceParticipationMetadata(const IoData &iodata)
{
  auto interfaces = nlohmann::json::array();
  for (const auto &[index, data] : iodata.boundaries.postpro.dielectric)
  {
    const double cutoff_mesh_units =
        data.edge_cutoff * iodata.units.GetMeshLengthRelativeScale();
    const double cutoff_meters =
        iodata.units.Dimensionalize<Units::ValueType::LENGTH>(data.edge_cutoff);
    interfaces.push_back(
        {{"Index", index},
         {"Regularization", data.edge_cutoff > 0.0 ? "explicit edge cutoff" : "none"},
         {"EdgeCutoffMeshUnits", cutoff_mesh_units},
         {"EdgeCutoffMeters", cutoff_meters}});
  }
  return {{"OutputFile", "surface-Q.csv"},
          {"Definition", "thin-layer surface integral; a positive edge cutoff excludes the "
                         "corresponding singular endpoint neighborhood"},
          {"Interfaces", std::move(interfaces)},
          {"IncludesExcludedEdgeNeighborhood", false},
          {"ResponseCorrected", false},
          {"UnregularizedNuAtMostOneHalfAccepted", false}};
}

namespace
{

nlohmann::json MakeSingularSurfaceIntegrabilityMetadata(std::vector<double> exponents)
{
  MFEM_VERIFY(!exponents.empty(),
              "Singular surface integrability metadata requires at least one exponent!");
  std::sort(exponents.begin(), exponents.end());
  exponents.erase(std::unique(exponents.begin(), exponents.end()), exponents.end());
  MFEM_VERIFY(std::isfinite(exponents.front()) && exponents.front() > 0.0 &&
                  std::isfinite(exponents.back()) && exponents.back() < 1.0,
              "Singular surface integrability metadata received an invalid exponent!");
  const double minimum_exponent = exponents.front();
  const double maximum_exponent = exponents.back();
  const bool finite = minimum_exponent > 0.5;
  return {{"Criterion", "finite exactly when every active electric exponent has nu > 1/2"},
          {"Exponents", std::move(exponents)},
          {"MinimumExponent", minimum_exponent},
          {"MaximumExponent", maximum_exponent},
          {"UnregularizedFinite", finite},
          {"RequiresPhysicalRegularization", !finite}};
}

}  // namespace

nlohmann::json
GetSingularSurfaceIntegrabilityMetadata(const fem::singular::FeatureTopology &features)
{
  std::vector<double> exponents;
  exponents.reserve(features.features.size());
  for (const auto &feature : features.features)
  {
    exponents.push_back(feature.nu);
  }
  return MakeSingularSurfaceIntegrabilityMetadata(std::move(exponents));
}

nlohmann::json GetSingularSurfaceIntegrabilityMetadata(
    const fem::singular::TriangleFeatureTopology &features)
{
  std::vector<double> exponents;
  exponents.reserve(features.vertices.size());
  for (const auto &vertex : features.vertices)
  {
    exponents.push_back(vertex.nu);
  }
  return MakeSingularSurfaceIntegrabilityMetadata(std::move(exponents));
}

void FullWaveSingularFeatures::Preprocess(const IoData &iodata,
                                          const std::unique_ptr<mfem::Mesh> &serial_mesh,
                                          MPI_Comm comm)
{
  serial_sheet_features = {};
  local_sheet_features = {};
  serial_line_features = {};
  local_line_features = {};
  source_vertex_ids.clear();
  source_element_ids.clear();
  dimension = 0;
  if (!iodata.solver.singular_elements.Enabled())
  {
    return;
  }

  if (Mpi::Root(comm))
  {
    MFEM_VERIFY(serial_mesh,
                "Root rank has no serial mesh for full-wave singular feature extraction!");
    dimension = serial_mesh->Dimension();
    if (dimension == 3)
    {
      serial_sheet_features = fem::singular::ExtractSerialSheetFeatures(
          *serial_mesh, iodata.solver.singular_elements.attributes,
          GetSingularTriangleMaterials(iodata));
      ValidateSingularLossTangents(iodata, *serial_mesh, serial_sheet_features);
    }
    else
    {
      MFEM_VERIFY(dimension == 2 && serial_mesh->SpaceDimension() == 2,
                  "Full-wave singular enrichment requires a 2D triangular or 3D "
                  "tetrahedral mesh!");
      MFEM_VERIFY(iodata.solver.singular_elements.order == 1,
                  "Triangular full-wave singular enrichment currently supports only "
                  "Solver.SingularElements.Order = 1!");
      serial_line_features = fem::singular::ExtractSerialLineFeatures(
          *serial_mesh, iodata.solver.singular_elements.attributes,
          GetSingularTriangleMaterials(iodata));
      ValidateSingularLossTangents(iodata, serial_line_features);
    }
  }
  Mpi::Broadcast(1, &dimension, 0, comm);
  if (dimension == 2)
  {
    fem::singular::BroadcastSerialLineTipFeatures(serial_line_features, comm);
    MFEM_VERIFY(!serial_line_features.Empty(),
                "Full-wave singular extraction produced no PEC line tips or corners!");
  }
  else
  {
    MFEM_VERIFY(dimension == 3, "Full-wave singular enrichment requires a 2D or 3D mesh!");
    fem::singular::BroadcastSerialSheetFeatures(serial_sheet_features, comm);
    MFEM_VERIFY(!serial_sheet_features.Empty(),
                "Full-wave singular extraction produced no PEC sheet-edge features!");
  }
}

void FullWaveSingularFeatures::ProcessPartitionedMesh(
    const IoData &iodata, const mfem::ParMesh &parallel_mesh,
    const mesh::PartitionMetadata &metadata)
{
  MFEM_VERIFY(iodata.solver.singular_elements.Enabled() &&
                  parallel_mesh.Dimension() == dimension,
              "Unexpected or inconsistent full-wave singular partition metadata!");
  if (dimension == 2)
  {
    local_line_features = fem::singular::DistributeSerialLineTipFeatures(
        serial_line_features, parallel_mesh, metadata.source_vertex_ids,
        metadata.source_element_ids);
  }
  else
  {
    local_sheet_features = fem::singular::DistributeSerialSheetFeatures(
        serial_sheet_features, parallel_mesh, metadata.source_vertex_ids,
        metadata.source_element_ids);
  }
  source_vertex_ids = metadata.source_vertex_ids;
  source_element_ids = metadata.source_element_ids;
  MFEM_VERIFY(source_vertex_ids.size() == static_cast<std::size_t>(parallel_mesh.GetNV()) &&
                  source_element_ids.size() ==
                      static_cast<std::size_t>(parallel_mesh.GetNE()),
              "Full-wave singular source-entity metadata is incomplete!");
}

const fem::singular::FeatureTopology *FullWaveSingularFeatures::GetSheetFeatures() const
{
  return dimension == 3 && !local_sheet_features.Empty() ? &local_sheet_features : nullptr;
}

const fem::singular::TriangleFeatureTopology *
FullWaveSingularFeatures::GetLineFeatures() const
{
  return dimension == 2 && !local_line_features.Empty() ? &local_line_features : nullptr;
}

const std::vector<fem::singular::GlobalVertexId> *
FullWaveSingularFeatures::GetSourceVertexIds() const
{
  return source_vertex_ids.empty() ? nullptr : &source_vertex_ids;
}

double SingularComplexQuadraticForm(MPI_Comm comm, const mfem::Operator &op,
                                    const ComplexVector &x)
{
  MFEM_VERIFY(op.Height() == op.Width() && op.Width() == x.Size(),
              "Invalid dimensions for a singular full-wave quadratic form!");
  Vector work(op.Height());
  work.UseDevice(true);
  op.Mult(x.Real(), work);
  double value = linalg::Dot(comm, x.Real(), work);
  op.Mult(x.Imag(), work);
  value += linalg::Dot(comm, x.Imag(), work);
  MFEM_VERIFY(std::isfinite(value),
              "Singular full-wave quadratic form produced a nonfinite value!");
  return value;
}

SingularFullWaveEnergy MeasureSingularFullWaveEnergy(MPI_Comm comm,
                                                     const mfem::Operator &mass,
                                                     const mfem::Operator &stiffness,
                                                     const ComplexVector &electric_field,
                                                     std::complex<double> omega)
{
  const double omega_squared = std::norm(omega);
  MFEM_VERIFY(omega_squared > 0.0 && std::isfinite(omega_squared),
              "Singular full-wave energy requires a finite nonzero frequency!");
  const double electric_form = SingularComplexQuadraticForm(comm, mass, electric_field);
  const double magnetic_form =
      SingularComplexQuadraticForm(comm, stiffness, electric_field) / omega_squared;
  const double scale = std::max({1.0, std::abs(electric_form), std::abs(magnetic_form)});
  constexpr double positivity_tolerance = 256.0 * std::numeric_limits<double>::epsilon();
  MFEM_VERIFY(electric_form >= -positivity_tolerance * scale &&
                  magnetic_form >= -positivity_tolerance * scale,
              "Singular full-wave energy is negative beyond floating-point roundoff!");
  return {0.5 * std::max(0.0, electric_form), 0.5 * std::max(0.0, magnetic_form)};
}

}  // namespace palace
