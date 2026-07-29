// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularsolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <type_traits>
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

void UpdateSingularSourceEntityIds(
    const mfem::ParMesh &mesh,
    std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
    std::vector<fem::singular::GlobalVertexId> &source_element_ids)
{
  using GlobalVertexId = fem::singular::GlobalVertexId;
  if (mesh.Nonconforming())
  {
    source_vertex_ids.assign(mesh.GetNV(), -1);
    const auto &ncmesh = *mesh.ncmesh;
    for (int node = 0; node < ncmesh.GetNumNodes(); node++)
    {
      const auto &record = ncmesh.GetNode(node);
      if (record.HasVertex() && record.vert_index >= 0 && record.vert_index < mesh.GetNV())
      {
        source_vertex_ids[record.vert_index] = node;
      }
    }
    MFEM_VERIFY(std::find(source_vertex_ids.begin(), source_vertex_ids.end(),
                          GlobalVertexId{-1}) == source_vertex_ids.end(),
                "Nonconforming singular mesh contains a vertex without a persistent "
                "NCMesh node ID!");
  }
  else
  {
    const std::size_t old_vertices = source_vertex_ids.size();
    MFEM_VERIFY(old_vertices <= static_cast<std::size_t>(mesh.GetNV()),
                "Conforming refinement removed an existing singular source vertex!");

    GlobalVertexId maximum_id = -1;
    for (GlobalVertexId id : source_vertex_ids)
    {
      maximum_id = std::max(maximum_id, id);
    }
    Mpi::GlobalMax(1, &maximum_id, mesh.GetComm());

    mfem::Array<HYPRE_BigInt> global_vertices;
    mesh.GetGlobalVertexIndices(global_vertices);
    MFEM_VERIFY(global_vertices.Size() == mesh.GetNV(),
                "Conforming refined mesh has incomplete global vertex numbering!");
    source_vertex_ids.resize(mesh.GetNV());
    for (int vertex = static_cast<int>(old_vertices); vertex < mesh.GetNV(); vertex++)
    {
      const HYPRE_BigInt global_vertex = global_vertices[vertex];
      MFEM_VERIFY(global_vertex >= 0 &&
                      global_vertex <=
                          std::numeric_limits<GlobalVertexId>::max() - maximum_id - 1,
                  "Refined singular vertex identity exceeds the supported integer range!");
      source_vertex_ids[vertex] =
          maximum_id + 1 + static_cast<GlobalVertexId>(global_vertex);
    }
  }

  std::set<GlobalVertexId> unique_vertices(source_vertex_ids.begin(),
                                           source_vertex_ids.end());
  MFEM_VERIFY(unique_vertices.size() == source_vertex_ids.size(),
              "Refined singular mesh contains duplicate local source vertex IDs!");

  source_element_ids.resize(mesh.GetNE());
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    source_element_ids[element] = mesh.GetGlobalElementNum(element);
  }
}

namespace
{

using StableFaceKey = std::array<fem::singular::GlobalVertexId, 3>;

template <typename FeatureTopology>
mfem::Array<int> BuildSingularRefinementProtectionImpl(
    const mfem::ParMesh &mesh, const FeatureTopology &features,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids)
{
  MFEM_VERIFY(features.elements.size() == static_cast<std::size_t>(mesh.GetNE()) &&
                  source_vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()),
              "Singular refinement protection received inconsistent feature topology!");
  const int dimension = mesh.Dimension();
  MFEM_VERIFY((dimension == 2 && mesh.SpaceDimension() == 2) ||
                  (dimension == 3 && mesh.SpaceDimension() == 3),
              "Singular refinement protection requires a 2D or 3D simplex mesh!");
  auto ElementFaceKeys = [&](int element)
  {
    std::vector<StableFaceKey> keys;
    const auto &record = *mesh.GetElement(element);
    const int *vertices = record.GetVertices();
    if (dimension == 2)
    {
      constexpr std::array<std::array<int, 2>, 3> faces{
          std::array<int, 2>{0, 1}, {1, 2}, {2, 0}};
      keys.reserve(faces.size());
      for (const auto &face : faces)
      {
        std::array<fem::singular::GlobalVertexId, 2> edge{
            source_vertex_ids[vertices[face[0]]], source_vertex_ids[vertices[face[1]]]};
        std::sort(edge.begin(), edge.end());
        keys.push_back({edge[0], edge[1], -1});
      }
    }
    else
    {
      constexpr std::array<std::array<int, 3>, 4> faces{
          std::array<int, 3>{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
      keys.reserve(faces.size());
      for (const auto &face : faces)
      {
        StableFaceKey key{source_vertex_ids[vertices[face[0]]],
                          source_vertex_ids[vertices[face[1]]],
                          source_vertex_ids[vertices[face[2]]]};
        std::sort(key.begin(), key.end());
        keys.push_back(key);
      }
    }
    return keys;
  };

  std::map<StableFaceKey, bool> conforming_faces;
  if (mesh.Nonconforming())
  {
    mfem::Array<int> vertices;
    for (int face = 0; face < mesh.GetNumFaces(); face++)
    {
      mesh.GetFaceVertices(face, vertices);
      MFEM_VERIFY(vertices.Size() == dimension,
                  "Singular refinement protection found a nonsimplex mesh face!");
      StableFaceKey key{-1, -1, -1};
      for (int vertex = 0; vertex < vertices.Size(); vertex++)
      {
        key[vertex] = source_vertex_ids[vertices[vertex]];
      }
      std::sort(key.begin(), key.begin() + vertices.Size());
      const auto information = mesh.GetFaceInformation(face);
      const bool conforming = information.IsBoundary() || information.IsConforming();
      auto [record, inserted] = conforming_faces.emplace(key, conforming);
      if (!inserted)
      {
        record->second = record->second && conforming;
      }
    }
  }

  std::vector<StableFaceKey> local_faces;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &incidence = features.elements[element];
    const bool enriched = [&]()
    {
      if constexpr (std::is_same_v<FeatureTopology, fem::singular::TriangleFeatureTopology>)
      {
        return !incidence.nodes.empty();
      }
      else
      {
        return !incidence.nodes.empty() || !incidence.edges.empty();
      }
    }();
    if (!enriched)
    {
      continue;
    }
    auto keys = ElementFaceKeys(element);
    if (mesh.Nonconforming())
    {
      for (const auto &key : keys)
      {
        const auto face = conforming_faces.find(key);
        MFEM_VERIFY(
            face != conforming_faces.end() && face->second,
            "A singular-enriched element touches a nonconforming interface. Refinement "
            "through the protected singular patch requires certified parent/child "
            "enrichment constraints!");
      }
    }
    local_faces.insert(local_faces.end(), keys.begin(), keys.end());
  }
  std::sort(local_faces.begin(), local_faces.end());
  local_faces.erase(std::unique(local_faces.begin(), local_faces.end()), local_faces.end());

  const std::size_t local_scalar_size = 3 * local_faces.size();
  MFEM_VERIFY(local_scalar_size <=
                  static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "Local singular refinement protection exceeds MPI integer counts!");
  const int local_count = static_cast<int>(local_scalar_size);
  std::vector<int> counts(Mpi::Size(mesh.GetComm()));
  Mpi::Allgather(1, &local_count, counts.data(), mesh.GetComm());
  std::vector<int> offsets(counts.size());
  std::partial_sum(counts.begin(), counts.end() - 1, offsets.begin() + 1);
  const std::size_t global_scalar_size =
      std::accumulate(counts.begin(), counts.end(), std::size_t{0});
  MFEM_VERIFY(global_scalar_size <=
                  static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "Global singular refinement protection exceeds MPI integer counts!");
  std::vector<fem::singular::GlobalVertexId> local_data(local_scalar_size);
  std::vector<fem::singular::GlobalVertexId> global_data(global_scalar_size);
  for (std::size_t face = 0; face < local_faces.size(); face++)
  {
    std::copy(local_faces[face].begin(), local_faces[face].end(),
              local_data.begin() + 3 * face);
  }
  Mpi::Allgatherv(local_count, local_data.data(), global_data.data(), counts.data(),
                  offsets.data(), mesh.GetComm());
  MFEM_VERIFY(global_data.size() % 3 == 0,
              "Gathered singular refinement protection is malformed!");

  std::set<StableFaceKey> protected_faces;
  for (std::size_t offset = 0; offset < global_data.size(); offset += 3)
  {
    protected_faces.insert(
        {global_data[offset], global_data[offset + 1], global_data[offset + 2]});
  }

  mfem::Array<int> marker(mesh.GetNE());
  marker = 0;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto keys = ElementFaceKeys(element);
    marker[element] =
        std::any_of(keys.begin(), keys.end(), [&protected_faces](const auto &key)
                    { return protected_faces.find(key) != protected_faces.end(); });
  }
  return marker;
}

}  // namespace

mfem::Array<int> BuildSingularRefinementProtection(
    const mfem::ParMesh &mesh, const fem::singular::FeatureTopology &features,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids)
{
  return BuildSingularRefinementProtectionImpl(mesh, features, source_vertex_ids);
}

mfem::Array<int> BuildSingularRefinementProtection(
    const mfem::ParMesh &mesh, const fem::singular::TriangleFeatureTopology &features,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids)
{
  return BuildSingularRefinementProtectionImpl(mesh, features, source_vertex_ids);
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

void FullWaveSingularFeatures::ProcessRefinedMesh(const IoData &iodata,
                                                  const mfem::ParMesh &parallel_mesh)
{
  MFEM_VERIFY(iodata.solver.singular_elements.Enabled() &&
                  parallel_mesh.Dimension() == dimension,
              "Unexpected or inconsistent refined full-wave singular mesh!");
  UpdateSingularSourceEntityIds(parallel_mesh, source_vertex_ids, source_element_ids);
  mesh::PartitionMetadata metadata{source_vertex_ids, source_element_ids};
  ProcessPartitionedMesh(iodata, parallel_mesh, metadata);
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

mesh::PartitionMetadata FullWaveSingularFeatures::GetSourceEntityMetadata() const
{
  return {source_vertex_ids, source_element_ids};
}

mfem::Array<int>
FullWaveSingularFeatures::GetRefinementProtection(const mfem::ParMesh &mesh) const
{
  MFEM_VERIFY(mesh.Dimension() == dimension && !source_vertex_ids.empty(),
              "Full-wave singular refinement protection is not initialized!");
  return dimension == 3 ? BuildSingularRefinementProtection(mesh, local_sheet_features,
                                                            source_vertex_ids)
                        : BuildSingularRefinementProtection(mesh, local_line_features,
                                                            source_vertex_ids);
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
