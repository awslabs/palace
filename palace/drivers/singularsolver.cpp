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

namespace
{

std::set<int> GetPecSingularAttributes(const IoData &iodata)
{
  std::set<int> attributes(iodata.boundaries.pec.attributes.begin(),
                           iodata.boundaries.pec.attributes.end());
  attributes.insert(iodata.boundaries.auxpec.attributes.begin(),
                    iodata.boundaries.auxpec.attributes.end());
  return attributes;
}

std::set<int> GetImpedanceSingularAttributes(const IoData &iodata)
{
  std::set<int> attributes;
  for (const auto &impedance : iodata.boundaries.impedance)
  {
    attributes.insert(impedance.attributes.begin(), impedance.attributes.end());
  }
  return attributes;
}

void VerifyImpedanceExponent(double nu, const char *description)
{
  MFEM_VERIFY(
      std::isfinite(nu) && nu > 0.5,
      description << " has nu <= 1/2. Its free singular tangential trace is not square "
                     "integrable in the surface-impedance Robin form. Use a resolved "
                     "finite-thickness corner with nu > 1/2, remove that attribute from "
                     "Solver.SingularElements, or use PEC enrichment!");
}

}  // namespace

void ValidateSingularImpedanceFeatures(const IoData &iodata,
                                       const fem::singular::FeatureTopology &features)
{
  const auto pec_attributes = GetPecSingularAttributes(iodata);
  const auto impedance_attributes = GetImpedanceSingularAttributes(iodata);
  for (const auto &segment : features.segments)
  {
    bool has_pec = false, has_impedance = false;
    for (int attribute : segment.boundary_attributes)
    {
      has_pec = has_pec || pec_attributes.count(attribute) > 0;
      has_impedance = has_impedance || impedance_attributes.count(attribute) > 0;
    }
    MFEM_VERIFY(!(has_pec && has_impedance),
                "A three-dimensional singular feature edge joins PEC and surface-impedance "
                "faces. The current enrichment cannot impose distinct trace spaces on one "
                "shared feature DOF!");
    if (has_impedance)
    {
      MFEM_VERIFY(segment.feature < features.features.size(),
                  "A singular impedance edge references an invalid straight feature!");
      VerifyImpedanceExponent(features.features[segment.feature].nu,
                              "A three-dimensional singular impedance edge");
    }
  }
}

void ValidateSingularImpedanceFeatures(
    const IoData &iodata, const fem::singular::TriangleFeatureTopology &features)
{
  const auto pec_attributes = GetPecSingularAttributes(iodata);
  const auto impedance_attributes = GetImpedanceSingularAttributes(iodata);
  for (const auto &vertex : features.vertices)
  {
    bool has_pec = false, has_impedance = false;
    for (std::size_t segment_index : vertex.selected_segments)
    {
      MFEM_VERIFY(segment_index < features.selected_segments.size(),
                  "A singular impedance corner references an invalid boundary segment!");
      const int attribute = features.selected_segments[segment_index].boundary_attribute;
      has_pec = has_pec || pec_attributes.count(attribute) > 0;
      has_impedance = has_impedance || impedance_attributes.count(attribute) > 0;
    }
    MFEM_VERIFY(
        !(has_pec && has_impedance),
        "A two-dimensional singular feature corner joins PEC and surface-impedance "
        "segments. The current enrichment cannot impose distinct trace spaces on one "
        "shared feature DOF!");
    if (has_impedance)
    {
      VerifyImpedanceExponent(vertex.nu, "A two-dimensional singular impedance corner");
    }
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

class NonconformingVertexIdentity::Impl
{
private:
  using GlobalVertexId = fem::singular::GlobalVertexId;
  using VertexKey = std::array<GlobalVertexId, 3>;

  struct GatheredKeys
  {
    std::vector<int> counts;
    std::vector<int> offsets;
    std::vector<VertexKey> keys;
  };

  std::map<VertexKey, GlobalVertexId> key_ids;
  std::map<GlobalVertexId, VertexKey> id_keys;
  std::map<int, GlobalVertexId> local_node_ids;
  GlobalVertexId next_id = 0;

  static VertexKey LeafKey(int node) { return {0, static_cast<GlobalVertexId>(node), 0}; }

  static VertexKey ParentKey(GlobalVertexId first, GlobalVertexId second)
  {
    if (second < first)
    {
      std::swap(first, second);
    }
    return {1, first, second};
  }

  static GatheredKeys Gather(MPI_Comm comm, const std::vector<VertexKey> &local_keys)
  {
    MFEM_VERIFY(local_keys.size() <
                    static_cast<std::size_t>(std::numeric_limits<int>::max() / 3),
                "Nonconforming vertex key count exceeds MPI integer counts!");
    GatheredKeys gathered;
    gathered.counts.resize(Mpi::Size(comm));
    const int local_count = static_cast<int>(local_keys.size());
    Mpi::Allgather(1, &local_count, gathered.counts.data(), comm);
    gathered.offsets.resize(gathered.counts.size());
    std::partial_sum(gathered.counts.begin(), gathered.counts.end() - 1,
                     gathered.offsets.begin() + 1);

    std::vector<int> value_counts(gathered.counts), value_offsets(gathered.offsets);
    for (int &count : value_counts)
    {
      count *= 3;
    }
    for (int &offset : value_offsets)
    {
      offset *= 3;
    }
    const int total = gathered.offsets.back() + gathered.counts.back();
    std::vector<GlobalVertexId> values(Mpi::Root(comm) ? 3 * total : 0);
    MPI_Gatherv(local_keys.data(), 3 * local_count, MPI_INT64_T, values.data(),
                value_counts.data(), value_offsets.data(), MPI_INT64_T, 0, comm);
    if (Mpi::Root(comm))
    {
      gathered.keys.resize(total);
      for (int i = 0; i < total; i++)
      {
        gathered.keys[i] = {values[3 * i], values[3 * i + 1], values[3 * i + 2]};
      }
    }
    return gathered;
  }

  void RegisterKeys(MPI_Comm comm, const std::vector<VertexKey> &local_keys,
                    const std::vector<GlobalVertexId> &local_ids)
  {
    MFEM_VERIFY(local_keys.size() == local_ids.size(),
                "Observed nonconforming vertex key and ID counts differ!");
    const auto gathered = Gather(comm, local_keys);
    const int total = gathered.offsets.back() + gathered.counts.back();
    std::vector<GlobalVertexId> ids(Mpi::Root(comm) ? total : 0);
    MPI_Gatherv(local_ids.data(), static_cast<int>(local_ids.size()), MPI_INT64_T,
                ids.data(), gathered.counts.data(), gathered.offsets.data(), MPI_INT64_T, 0,
                comm);
    if (Mpi::Root(comm))
    {
      for (int i = 0; i < total; i++)
      {
        const auto [key, inserted_key] = key_ids.emplace(gathered.keys[i], ids[i]);
        MFEM_VERIFY(inserted_key || key->second == ids[i],
                    "One refined vertex key has inconsistent persistent IDs!");
        const auto [id, inserted_id] = id_keys.emplace(ids[i], gathered.keys[i]);
        MFEM_VERIFY(inserted_id || id->second == gathered.keys[i],
                    "Distinct refined vertex keys share one persistent ID!");
      }
    }
  }

  std::vector<GlobalVertexId> AssignKeys(MPI_Comm comm,
                                         const std::vector<VertexKey> &local_keys)
  {
    const auto gathered = Gather(comm, local_keys);
    const int total = gathered.offsets.back() + gathered.counts.back();
    std::vector<GlobalVertexId> ids(Mpi::Root(comm) ? total : 0);
    if (Mpi::Root(comm))
    {
      std::set<VertexKey> missing;
      for (const auto &key : gathered.keys)
      {
        if (key_ids.find(key) == key_ids.end())
        {
          missing.insert(key);
        }
      }
      for (const auto &key : missing)
      {
        MFEM_VERIFY(next_id < std::numeric_limits<GlobalVertexId>::max(),
                    "Persistent nonconforming vertex identity range is exhausted!");
        const GlobalVertexId id = next_id++;
        key_ids.emplace(key, id);
        const auto [record, inserted] = id_keys.emplace(id, key);
        MFEM_VERIFY(inserted && record->second == key,
                    "Failed to assign a unique nonconforming vertex identity!");
      }
      for (int i = 0; i < total; i++)
      {
        ids[i] = key_ids.at(gathered.keys[i]);
      }
    }
    std::vector<GlobalVertexId> local_ids(local_keys.size());
    MPI_Scatterv(ids.data(), gathered.counts.data(), gathered.offsets.data(), MPI_INT64_T,
                 local_ids.data(), static_cast<int>(local_ids.size()), MPI_INT64_T, 0,
                 comm);
    return local_ids;
  }

public:
  void Clear()
  {
    key_ids.clear();
    id_keys.clear();
    local_node_ids.clear();
    next_id = 0;
  }

  void Observe(const mfem::ParMesh &mesh, const std::vector<GlobalVertexId> &vertex_ids)
  {
    MPI_Comm comm = mesh.GetComm();
    bool valid = vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV());
    Mpi::GlobalAnd(1, &valid, comm);
    MFEM_VERIFY(valid,
                "Observed nonconforming vertex identities have an invalid local size!");

    local_node_ids.clear();
    const auto vertex_nodes = mesh::GetNonconformingVertexNodeIds(mesh);
    std::map<GlobalVertexId, int> local_ids;
    for (int vertex = 0; vertex < mesh.GetNV(); vertex++)
    {
      const auto id = vertex_ids[vertex];
      MFEM_VERIFY(id >= 0, "Observed a negative nonconforming vertex identity!");
      const auto [node, inserted_node] = local_node_ids.emplace(vertex_nodes[vertex], id);
      MFEM_VERIFY(inserted_node || node->second == id,
                  "One local NCMesh node has inconsistent persistent IDs!");
      const auto [identity, inserted_identity] =
          local_ids.emplace(id, vertex_nodes[vertex]);
      MFEM_VERIFY(inserted_identity || identity->second == vertex_nodes[vertex],
                  "Distinct local NCMesh vertices share one persistent ID!");
    }
    GlobalVertexId maximum_id = -1;
    for (const auto &[node, id] : local_node_ids)
    {
      MFEM_CONTRACT_VAR(node);
      maximum_id = std::max(maximum_id, id);
    }
    Mpi::GlobalMax(1, &maximum_id, comm);
    MFEM_VERIFY(maximum_id < std::numeric_limits<GlobalVertexId>::max(),
                "Persistent nonconforming vertex identity range is exhausted!");
    if (Mpi::Root(comm))
    {
      next_id = std::max(next_id, maximum_id + 1);
    }

    std::vector<VertexKey> local_keys;
    std::vector<GlobalVertexId> local_key_ids;
    const auto &ncmesh = *mesh.ncmesh;
    for (const auto &[node, id] : local_node_ids)
    {
      const auto &record = ncmesh.GetNode(node);
      if (record.p1 == record.p2 && record.p1 == node)
      {
        local_keys.push_back(LeafKey(node));
        local_key_ids.push_back(id);
      }
      else
      {
        const auto first = local_node_ids.find(record.p1);
        const auto second = local_node_ids.find(record.p2);
        if (first != local_node_ids.end() && second != local_node_ids.end())
        {
          local_keys.push_back(ParentKey(first->second, second->second));
          local_key_ids.push_back(id);
        }
      }
    }
    RegisterKeys(comm, local_keys, local_key_ids);
  }

  void Update(const mfem::ParMesh &mesh, std::vector<GlobalVertexId> &vertex_ids)
  {
    MPI_Comm comm = mesh.GetComm();
    MFEM_VERIFY(mesh.Nonconforming() && mesh.ncmesh,
                "Nonconforming vertex identity update requires an NCMesh!");
    const auto vertex_nodes = mesh::GetNonconformingVertexNodeIds(mesh);
    std::set<int> unresolved;
    const auto &ncmesh = *mesh.ncmesh;
    const auto Collect = [&](const auto &self, int node) -> void
    {
      if (local_node_ids.find(node) != local_node_ids.end() ||
          unresolved.find(node) != unresolved.end())
      {
        return;
      }
      MFEM_VERIFY(node >= 0 && node < ncmesh.GetNumNodes(),
                  "A newly refined NCMesh vertex has an unavailable parent node!");
      unresolved.insert(node);
      const auto &record = ncmesh.GetNode(node);
      if (record.p1 != record.p2)
      {
        self(self, record.p1);
        self(self, record.p2);
      }
      else
      {
        MFEM_VERIFY(record.p1 == node,
                    "An NCMesh top-level vertex has an invalid parent relation!");
      }
    };
    for (int node : vertex_nodes)
    {
      Collect(Collect, node);
    }

    int global_unresolved = unresolved.size();
    Mpi::GlobalSum(1, &global_unresolved, comm);
    while (global_unresolved > 0)
    {
      std::vector<int> ready_nodes;
      std::vector<VertexKey> ready_keys;
      for (int node : unresolved)
      {
        const auto &record = ncmesh.GetNode(node);
        if (record.p1 == record.p2)
        {
          ready_nodes.push_back(node);
          ready_keys.push_back(LeafKey(node));
          continue;
        }
        const auto first = local_node_ids.find(record.p1);
        const auto second = local_node_ids.find(record.p2);
        if (first != local_node_ids.end() && second != local_node_ids.end())
        {
          ready_nodes.push_back(node);
          ready_keys.push_back(ParentKey(first->second, second->second));
        }
      }
      int global_ready = ready_nodes.size();
      Mpi::GlobalSum(1, &global_ready, comm);
      MFEM_VERIFY(global_ready > 0,
                  "Unable to resolve newly refined NCMesh vertex parent identities!");
      const auto ready_ids = AssignKeys(comm, ready_keys);
      for (std::size_t i = 0; i < ready_nodes.size(); i++)
      {
        const auto [record, inserted] =
            local_node_ids.emplace(ready_nodes[i], ready_ids[i]);
        MFEM_VERIFY(inserted || record->second == ready_ids[i],
                    "A newly refined NCMesh node received inconsistent identities!");
        unresolved.erase(ready_nodes[i]);
      }
      global_unresolved = unresolved.size();
      Mpi::GlobalSum(1, &global_unresolved, comm);
    }

    vertex_ids.resize(mesh.GetNV());
    std::set<GlobalVertexId> unique_ids;
    for (int vertex = 0; vertex < mesh.GetNV(); vertex++)
    {
      vertex_ids[vertex] = local_node_ids.at(vertex_nodes[vertex]);
      unique_ids.insert(vertex_ids[vertex]);
    }
    MFEM_VERIFY(unique_ids.size() == vertex_ids.size(),
                "Distinct local refined vertices share one persistent identity!");
  }
};

NonconformingVertexIdentity::NonconformingVertexIdentity() : impl(std::make_unique<Impl>())
{
}

NonconformingVertexIdentity::~NonconformingVertexIdentity() = default;

void NonconformingVertexIdentity::Clear()
{
  impl->Clear();
}

void NonconformingVertexIdentity::Observe(
    const mfem::ParMesh &mesh, const std::vector<fem::singular::GlobalVertexId> &vertex_ids)
{
  impl->Observe(mesh, vertex_ids);
}

void NonconformingVertexIdentity::Update(
    const mfem::ParMesh &mesh, std::vector<fem::singular::GlobalVertexId> &vertex_ids)
{
  impl->Update(mesh, vertex_ids);
}

void UpdateSingularSourceEntityIds(
    const mfem::ParMesh &mesh,
    std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
    std::vector<fem::singular::GlobalVertexId> &source_element_ids,
    NonconformingVertexIdentity &vertex_identity)
{
  using GlobalVertexId = fem::singular::GlobalVertexId;
  if (mesh.Nonconforming())
  {
    vertex_identity.Update(mesh, source_vertex_ids);
  }
  else
  {
    const std::size_t old_vertices = source_vertex_ids.size();
    MFEM_VERIFY(old_vertices <= static_cast<std::size_t>(mesh.GetNV()),
                "Conforming refinement removed an existing singular source vertex!");

    GlobalVertexId maximum_id = -1;
    for (GlobalVertexId id : source_vertex_ids)
    {
      MFEM_VERIFY(id >= 0,
                  "Conforming singular mesh contains an invalid source vertex ID!");
      maximum_id = std::max(maximum_id, id);
    }
    Mpi::GlobalMax(1, &maximum_id, mesh.GetComm());
    MFEM_VERIFY(maximum_id <= std::numeric_limits<int>::max(),
                "Conforming singular vertex identity exceeds the supported mesh-ID "
                "range!");

    const int new_vertices = mesh.GetNV() - static_cast<int>(old_vertices);
    std::vector<int> counts(Mpi::Size(mesh.GetComm()));
    Mpi::Allgather(1, &new_vertices, counts.data(), mesh.GetComm());
    std::vector<int> offsets(counts.size());
    std::partial_sum(counts.begin(), counts.end() - 1, offsets.begin() + 1);
    const std::size_t gathered_size =
        std::accumulate(counts.begin(), counts.end(), std::size_t{0});
    MFEM_VERIFY(gathered_size <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
                "Conforming refined vertex metadata exceeds MPI integer counts!");

    mfem::Array<HYPRE_BigInt> global_vertices;
    mesh.GetGlobalVertexIndices(global_vertices);
    MFEM_VERIFY(global_vertices.Size() == mesh.GetNV(),
                "Conforming refined mesh has incomplete global vertex numbering!");
    std::vector<GlobalVertexId> local_global_vertices(new_vertices);
    std::vector<double> local_coordinates(3 * new_vertices);
    const int dimension = mesh.SpaceDimension();
    MFEM_VERIFY(dimension == 2 || dimension == 3,
                "Singular refinement requires a two- or three-dimensional mesh!");
    for (int local = 0; local < new_vertices; local++)
    {
      const int vertex = static_cast<int>(old_vertices) + local;
      const HYPRE_BigInt global_vertex = global_vertices[vertex];
      MFEM_VERIFY(global_vertex >= 0 && static_cast<unsigned long long>(global_vertex) <=
                                            static_cast<unsigned long long>(
                                                std::numeric_limits<GlobalVertexId>::max()),
                  "Conforming refined mesh has an invalid global vertex number!");
      local_global_vertices[local] = static_cast<GlobalVertexId>(global_vertex);
      const double *coordinate = mesh.GetVertex(vertex);
      for (int d = 0; d < dimension; d++)
      {
        MFEM_VERIFY(std::isfinite(coordinate[d]),
                    "Conforming refined mesh has a nonfinite vertex coordinate!");
        local_coordinates[3 * local + d] = coordinate[d] == 0.0 ? 0.0 : coordinate[d];
      }
      for (int d = dimension; d < 3; d++)
      {
        local_coordinates[3 * local + d] = 0.0;
      }
    }

    std::vector<int> coordinate_counts(counts), coordinate_offsets(offsets);
    for (int &count : coordinate_counts)
    {
      MFEM_VERIFY(count <= std::numeric_limits<int>::max() / 3,
                  "Conforming refined vertex coordinates exceed MPI integer counts!");
      count *= 3;
    }
    for (int &offset : coordinate_offsets)
    {
      MFEM_VERIFY(offset <= std::numeric_limits<int>::max() / 3,
                  "Conforming refined vertex coordinates exceed MPI integer counts!");
      offset *= 3;
    }

    const bool root = Mpi::Root(mesh.GetComm());
    std::vector<GlobalVertexId> gathered_global_vertices(root ? gathered_size : 0);
    std::vector<double> gathered_coordinates(root ? 3 * gathered_size : 0);
    MPI_Gatherv(local_global_vertices.data(), new_vertices, MPI_INT64_T,
                gathered_global_vertices.data(), counts.data(), offsets.data(), MPI_INT64_T,
                0, mesh.GetComm());
    MPI_Gatherv(local_coordinates.data(), 3 * new_vertices, MPI_DOUBLE,
                gathered_coordinates.data(), coordinate_counts.data(),
                coordinate_offsets.data(), MPI_DOUBLE, 0, mesh.GetComm());

    struct VertexRecord
    {
      std::array<double, 3> coordinate{};
      GlobalVertexId id = -1;
    };
    std::vector<GlobalVertexId> gathered_ids(root ? gathered_size : 0);
    bool valid = true;
    if (root)
    {
      std::map<GlobalVertexId, VertexRecord> records;
      for (std::size_t occurrence = 0; occurrence < gathered_size; occurrence++)
      {
        std::array<double, 3> coordinate{gathered_coordinates[3 * occurrence],
                                         gathered_coordinates[3 * occurrence + 1],
                                         gathered_coordinates[3 * occurrence + 2]};
        const auto [record, inserted] = records.emplace(
            gathered_global_vertices[occurrence], VertexRecord{coordinate, -1});
        if (!inserted && record->second.coordinate != coordinate)
        {
          valid = false;
        }
      }

      std::vector<VertexRecord *> ordered;
      ordered.reserve(records.size());
      for (auto &[global_vertex, record] : records)
      {
        ordered.push_back(&record);
      }
      std::sort(ordered.begin(), ordered.end(),
                [](const VertexRecord *left, const VertexRecord *right)
                { return left->coordinate < right->coordinate; });
      for (std::size_t vertex = 1; vertex < ordered.size(); vertex++)
      {
        if (ordered[vertex - 1]->coordinate == ordered[vertex]->coordinate)
        {
          valid = false;
        }
      }
      if (ordered.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max() - maximum_id))
      {
        valid = false;
      }
      if (valid)
      {
        for (std::size_t vertex = 0; vertex < ordered.size(); vertex++)
        {
          ordered[vertex]->id = maximum_id + 1 + static_cast<GlobalVertexId>(vertex);
        }
        for (std::size_t occurrence = 0; occurrence < gathered_size; occurrence++)
        {
          gathered_ids[occurrence] = records.at(gathered_global_vertices[occurrence]).id;
        }
      }
    }
    MPI_Bcast(&valid, 1, MPI_C_BOOL, 0, mesh.GetComm());
    MFEM_VERIFY(valid,
                "Conforming singular refinement could not assign canonical new vertex "
                "IDs. Distinct refined vertices must not have coincident coordinates!");

    source_vertex_ids.resize(mesh.GetNV());
    MPI_Scatterv(gathered_ids.data(), counts.data(), offsets.data(), MPI_INT64_T,
                 source_vertex_ids.data() + old_vertices, new_vertices, MPI_INT64_T, 0,
                 mesh.GetComm());
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

using GlobalVertexId = fem::singular::GlobalVertexId;
using StableEdgeKey = std::array<GlobalVertexId, 2>;
using StableFaceKey = std::array<GlobalVertexId, 3>;

template <std::size_t N>
std::map<std::array<GlobalVertexId, N>, int>
GatherSelectedBoundaryEntities(const mfem::ParMesh &mesh,
                               const std::vector<int> &boundary_attributes,
                               const std::vector<GlobalVertexId> &source_vertex_ids)
{
  MFEM_VERIFY(source_vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()),
              "Singular boundary reconstruction received incomplete vertex identities!");
  const std::set<int> selected(boundary_attributes.begin(), boundary_attributes.end());
  std::vector<GlobalVertexId> local;
  std::set<int> local_attributes;
  for (int boundary_element = 0; boundary_element < mesh.GetNBE(); boundary_element++)
  {
    const int attribute = mesh.GetBdrAttribute(boundary_element);
    if (selected.find(attribute) == selected.end())
    {
      continue;
    }
    local_attributes.insert(attribute);
    const auto &element = *mesh.GetBdrElement(boundary_element);
    MFEM_VERIFY(element.GetNVertices() == static_cast<int>(N),
                "A selected singular boundary entity has invalid simplex topology!");
    std::array<GlobalVertexId, N> key;
    const int *vertices = element.GetVertices();
    for (std::size_t vertex = 0; vertex < N; vertex++)
    {
      MFEM_VERIFY(vertices[vertex] >= 0 &&
                      vertices[vertex] < static_cast<int>(mesh.GetNV()),
                  "A selected singular boundary entity has an invalid mesh vertex!");
      key[vertex] = source_vertex_ids[vertices[vertex]];
    }
    std::sort(key.begin(), key.end());
    MFEM_VERIFY(key.front() >= 0 && std::adjacent_find(key.begin(), key.end()) == key.end(),
                "A selected singular boundary entity has invalid stable vertex IDs!");
    local.insert(local.end(), key.begin(), key.end());
    local.push_back(attribute);
  }

  constexpr std::size_t width = N + 1;
  MFEM_VERIFY(local.size() % width == 0 &&
                  local.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "Local singular boundary topology exceeds MPI integer counts!");
  const int local_count = static_cast<int>(local.size());
  std::vector<int> counts(Mpi::Size(mesh.GetComm()));
  Mpi::Allgather(1, &local_count, counts.data(), mesh.GetComm());
  std::vector<int> offsets(counts.size());
  std::partial_sum(counts.begin(), counts.end() - 1, offsets.begin() + 1);
  const std::size_t global_count =
      std::accumulate(counts.begin(), counts.end(), std::size_t{0});
  MFEM_VERIFY(global_count % width == 0 &&
                  global_count <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "Global singular boundary topology exceeds MPI integer counts!");
  std::vector<GlobalVertexId> global(global_count);
  Mpi::Allgatherv(local_count, local.data(), global.data(), counts.data(), offsets.data(),
                  mesh.GetComm());

  std::map<std::array<GlobalVertexId, N>, int> result;
  std::set<int> present;
  for (std::size_t offset = 0; offset < global.size(); offset += width)
  {
    std::array<GlobalVertexId, N> key;
    std::copy_n(global.begin() + offset, N, key.begin());
    const auto raw_attribute = global[offset + N];
    MFEM_VERIFY(raw_attribute > 0 && raw_attribute <= static_cast<GlobalVertexId>(
                                                          std::numeric_limits<int>::max()),
                "A gathered singular boundary entity has an invalid attribute!");
    const int attribute = static_cast<int>(raw_attribute);
    const auto [record, inserted] = result.emplace(key, attribute);
    MFEM_VERIFY(inserted || record->second == attribute,
                "One singular boundary entity has inconsistent attributes across MPI "
                "ranks!");
    present.insert(attribute);
  }
  MFEM_VERIFY(present == selected,
              "Adaptive refinement removed a selected singular boundary attribute!");
  return result;
}

int CheckedMeshEntityId(GlobalVertexId id)
{
  MFEM_VERIFY(id >= 0 && id <= static_cast<GlobalVertexId>(std::numeric_limits<int>::max()),
              "Refined singular feature identity exceeds the supported mesh-ID range!");
  return static_cast<int>(id);
}

std::map<GlobalVertexId, std::array<double, 3>>
GatherStableVertexCoordinates(const mfem::ParMesh &mesh,
                              const std::vector<GlobalVertexId> &source_vertex_ids,
                              const std::set<GlobalVertexId> &requested)
{
  const int dimension = mesh.SpaceDimension();
  MFEM_VERIFY(dimension == 2 || dimension == 3,
              "Singular feature coordinates require a two- or three-dimensional mesh!");
  std::vector<GlobalVertexId> ids(requested.begin(), requested.end());
  std::map<GlobalVertexId, std::size_t> index;
  for (std::size_t i = 0; i < ids.size(); i++)
  {
    index.emplace(ids[i], i);
  }
  std::vector<int> counts(ids.size(), 0);
  std::vector<double> minima(3 * ids.size(), std::numeric_limits<double>::infinity());
  std::vector<double> maxima(3 * ids.size(), -std::numeric_limits<double>::infinity());
  std::array<double, 3> coordinate{};
  for (int vertex = 0; vertex < mesh.GetNV(); vertex++)
  {
    const auto requested_vertex = index.find(source_vertex_ids[vertex]);
    if (requested_vertex == index.end())
    {
      continue;
    }
    const auto i = requested_vertex->second;
    counts[i] = 1;
    mesh.GetNode(vertex, coordinate.data());
    for (int d = 0; d < dimension; d++)
    {
      MFEM_VERIFY(std::isfinite(coordinate[d]),
                  "A refined singular feature vertex has a nonfinite coordinate!");
      minima[3 * i + d] = coordinate[d];
      maxima[3 * i + d] = coordinate[d];
    }
    for (int d = dimension; d < 3; d++)
    {
      minima[3 * i + d] = maxima[3 * i + d] = 0.0;
    }
  }
  if (!ids.empty())
  {
    Mpi::GlobalSum(static_cast<int>(counts.size()), counts.data(), mesh.GetComm());
    Mpi::GlobalMin(static_cast<int>(minima.size()), minima.data(), mesh.GetComm());
    Mpi::GlobalMax(static_cast<int>(maxima.size()), maxima.data(), mesh.GetComm());
  }

  std::map<GlobalVertexId, std::array<double, 3>> result;
  for (std::size_t i = 0; i < ids.size(); i++)
  {
    MFEM_VERIFY(counts[i] > 0,
                "A refined singular feature references a missing stable vertex!");
    std::array<double, 3> value{};
    for (int d = 0; d < dimension; d++)
    {
      const double scale =
          std::max({1.0, std::abs(minima[3 * i + d]), std::abs(maxima[3 * i + d])});
      const double tolerance = 4096.0 * std::numeric_limits<double>::epsilon() * scale;
      MFEM_VERIFY(maxima[3 * i + d] - minima[3 * i + d] <= tolerance,
                  "A shared refined singular vertex has inconsistent coordinates!");
      value[d] = 0.5 * (minima[3 * i + d] + maxima[3 * i + d]);
    }
    result.emplace(ids[i], value);
  }
  return result;
}

double Norm(const std::array<double, 3> &value)
{
  return std::sqrt(value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
}

struct ChildFeatureSegment
{
  StableEdgeKey edge;
  std::size_t parent = 0;
  double first = 0.0;
  double second = 0.0;
};

bool IsChildSegment(const StableEdgeKey &candidate, const StableEdgeKey &parent,
                    const std::map<GlobalVertexId, std::array<double, 3>> &coordinates,
                    ChildFeatureSegment &child)
{
  const auto &a = coordinates.at(parent[0]);
  const auto &b = coordinates.at(parent[1]);
  const auto &p = coordinates.at(candidate[0]);
  const auto &q = coordinates.at(candidate[1]);
  std::array<double, 3> direction{}, offset_p{}, offset_q{};
  double length_squared = 0.0;
  for (int d = 0; d < 3; d++)
  {
    direction[d] = b[d] - a[d];
    offset_p[d] = p[d] - a[d];
    offset_q[d] = q[d] - a[d];
    length_squared += direction[d] * direction[d];
  }
  MFEM_VERIFY(length_squared > 0.0 && std::isfinite(length_squared),
              "A parent singular segment has invalid physical length!");
  const double first = (offset_p[0] * direction[0] + offset_p[1] * direction[1] +
                        offset_p[2] * direction[2]) /
                       length_squared;
  const double second = (offset_q[0] * direction[0] + offset_q[1] * direction[1] +
                         offset_q[2] * direction[2]) /
                        length_squared;
  std::array<double, 3> residual_p{}, residual_q{};
  for (int d = 0; d < 3; d++)
  {
    residual_p[d] = offset_p[d] - first * direction[d];
    residual_q[d] = offset_q[d] - second * direction[d];
  }
  const double length = std::sqrt(length_squared);
  const double coordinate_scale = std::max({1.0, Norm(a), Norm(b), Norm(p), Norm(q)});
  const double distance_tolerance =
      1.0e-10 * length + 4096.0 * std::numeric_limits<double>::epsilon() * coordinate_scale;
  const double parameter_tolerance = distance_tolerance / length;
  if (Norm(residual_p) > distance_tolerance || Norm(residual_q) > distance_tolerance ||
      first < -parameter_tolerance || first > 1.0 + parameter_tolerance ||
      second < -parameter_tolerance || second > 1.0 + parameter_tolerance ||
      std::abs(second - first) <= parameter_tolerance)
  {
    return false;
  }
  child.edge = candidate;
  child.first = std::clamp(std::min(first, second), 0.0, 1.0);
  child.second = std::clamp(std::max(first, second), 0.0, 1.0);
  return true;
}

template <typename FeatureTopology>
mfem::Array<int> BuildSingularRefinementProtectionImpl(
    const mfem::ParMesh &mesh, const FeatureTopology &features,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids, bool *conforming,
    mfem::Array<int> *repair)
{
  MFEM_VERIFY(features.elements.size() == static_cast<std::size_t>(mesh.GetNE()) &&
                  source_vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()),
              "Singular refinement protection received inconsistent feature topology!");
  const int dimension = mesh.Dimension();
  MFEM_VERIFY((dimension == 2 && mesh.SpaceDimension() == 2) ||
                  (dimension == 3 && mesh.SpaceDimension() == 3),
              "Singular refinement protection requires a 2D or 3D simplex mesh!");

  mfem::Array<int> local_enriched(mesh.GetNE());
  local_enriched = 0;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &incidence = features.elements[element];
    if constexpr (std::is_same_v<FeatureTopology, fem::singular::TriangleFeatureTopology>)
    {
      local_enriched[element] = !incidence.nodes.empty();
    }
    else
    {
      local_enriched[element] = !incidence.nodes.empty() || !incidence.edges.empty();
    }
  }

  // Exchange one scalar per element so the closure follows MFEM's actual shared and
  // nonconforming master/slave face topology. Exact face-vertex matching is insufficient
  // once one side of a face has been refined.
  auto &parallel_mesh = const_cast<mfem::ParMesh &>(mesh);
  mfem::L2_FECollection collection(0, dimension);
  mfem::ParFiniteElementSpace space(&parallel_mesh, &collection);
  mfem::ParGridFunction enriched(&space);
  enriched = 0.0;
  mfem::Array<int> dofs;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    if (!local_enriched[element])
    {
      continue;
    }
    space.GetElementDofs(element, dofs);
    MFEM_VERIFY(dofs.Size() == 1,
                "Singular refinement metadata must have one value per element!");
    enriched[dofs[0]] = 1.0;
  }
  enriched.ExchangeFaceNbrData();

  std::vector<bool> all_enriched(mesh.GetNE() + parallel_mesh.GetNFaceNeighborElements(),
                                 false);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    all_enriched[element] = local_enriched[element] != 0;
  }
  mfem::Vector value;
  for (int neighbor = 0; neighbor < parallel_mesh.GetNFaceNeighborElements(); neighbor++)
  {
    enriched.GetElementDofValues(mesh.GetNE() + neighbor, value);
    MFEM_VERIFY(value.Size() == 1 &&
                    (std::abs(value[0]) <= 1.0e-12 || std::abs(value[0] - 1.0) <= 1.0e-12),
                "Singular refinement metadata exchange is inconsistent!");
    all_enriched[mesh.GetNE() + neighbor] = value[0] > 0.5;
  }

  std::unique_ptr<mfem::Table> face_elements(parallel_mesh.GetFaceToAllElementTable());
  MFEM_VERIFY(face_elements,
              "Singular refinement closure requires face-to-element topology!");
  mfem::Array<int> marker(mesh.GetNE());
  marker = 0;
  mfem::Array<int> local_repair(mesh.GetNE());
  local_repair = 0;
  bool locally_conforming = true;
  mfem::Array<int> elements;
  for (int face = 0; face < face_elements->Size(); face++)
  {
    face_elements->GetRow(face, elements);
    bool has_local = false;
    bool touches_enrichment = false;
    for (int element : elements)
    {
      if (element < 0)
      {
        // Parallel nonconforming face tables use -1 when a face has no local
        // element on this rank.
        continue;
      }
      MFEM_VERIFY(
          element >= 0 && element < static_cast<int>(all_enriched.size()),
          "Singular refinement face topology has invalid element "
              << element << " on face " << face << " (local elements = " << mesh.GetNE()
              << ", face-neighbor elements = " << parallel_mesh.GetNFaceNeighborElements()
              << ", metadata size = " << all_enriched.size() << ")!");
      if (element < mesh.GetNE())
      {
        has_local = true;
      }
      touches_enrichment = touches_enrichment || all_enriched[element];
    }
    if (!has_local || !touches_enrichment)
    {
      continue;
    }
    const auto information = mesh.GetFaceInformation(face);
    const bool face_conforming = information.IsBoundary() || information.IsConforming();
    locally_conforming = locally_conforming && face_conforming;
    if (!face_conforming)
    {
      for (const auto &side : information.element)
      {
        if (side.location == mfem::Mesh::ElementLocation::Local &&
            side.conformity == mfem::Mesh::ElementConformity::Superset)
        {
          MFEM_VERIFY(side.index >= 0 && side.index < mesh.GetNE(),
                      "Singular refinement repair has an invalid coarse element!");
          local_repair[side.index] = 1;
        }
      }
      if (information.tag == mfem::Mesh::FaceInfoTag::MasterNonconforming)
      {
        const int coarse = information.element[0].index;
        MFEM_VERIFY(coarse >= 0 && coarse < mesh.GetNE(),
                    "Singular refinement repair has an invalid master element!");
        local_repair[coarse] = 1;
      }
    }
    for (int element : elements)
    {
      if (element >= 0 && element < mesh.GetNE())
      {
        marker[element] = 1;
      }
    }
  }
  Mpi::GlobalAnd(1, &locally_conforming, mesh.GetComm());
  if (conforming)
  {
    *conforming = locally_conforming;
  }
  if (repair)
  {
    *repair = std::move(local_repair);
  }
  return marker;
}

}  // namespace

void RebuildRefinedSingularFeatures(
    const mfem::ParMesh &mesh, const std::vector<int> &boundary_attributes,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
    fem::singular::TriangleFeatureTopology &features)
{
  MFEM_VERIFY(mesh.Dimension() == 2 && mesh.SpaceDimension() == 2 &&
                  source_vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()) &&
                  !features.vertices.empty(),
              "Refined triangular singular topology is inconsistent!");
  const auto selected =
      GatherSelectedBoundaryEntities<2>(mesh, boundary_attributes, source_vertex_ids);

  fem::singular::TriangleFeatureTopology updated;
  updated.vertices = features.vertices;
  updated.selected_segments.reserve(selected.size());
  std::map<GlobalVertexId, std::vector<std::size_t>> incident_segments;
  for (const auto &[edge, attribute] : selected)
  {
    const std::size_t segment = updated.selected_segments.size();
    updated.selected_segments.push_back(
        {static_cast<int>(segment),
         static_cast<int>(segment),
         {CheckedMeshEntityId(edge[0]), CheckedMeshEntityId(edge[1])},
         attribute});
    incident_segments[edge[0]].push_back(segment);
    incident_segments[edge[1]].push_back(segment);
  }
  for (auto &vertex : updated.vertices)
  {
    const GlobalVertexId id = vertex.mesh_vertex;
    const auto incident = incident_segments.find(id);
    MFEM_VERIFY(incident != incident_segments.end() && !incident->second.empty(),
                "Adaptive refinement disconnected a singular line tip or corner from "
                "its selected PEC boundary!");
    vertex.selected_segments = incident->second;
  }
  updated.elements.clear();
  features = std::move(updated);
}

void RebuildRefinedSingularFeatures(
    const mfem::ParMesh &mesh, const std::vector<int> &boundary_attributes,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
    fem::singular::FeatureTopology &features)
{
  MFEM_VERIFY(mesh.Dimension() == 3 && mesh.SpaceDimension() == 3 &&
                  source_vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()) &&
                  !features.features.empty() && !features.segments.empty(),
              "Refined tetrahedral singular topology is inconsistent!");
  const auto selected_faces =
      GatherSelectedBoundaryEntities<3>(mesh, boundary_attributes, source_vertex_ids);

  std::set<StableEdgeKey> selected_edges;
  std::set<GlobalVertexId> coordinate_ids;
  for (const auto &[face, attribute] : selected_faces)
  {
    MFEM_CONTRACT_VAR(attribute);
    coordinate_ids.insert(face.begin(), face.end());
    selected_edges.insert({face[0], face[1]});
    selected_edges.insert({face[0], face[2]});
    selected_edges.insert({face[1], face[2]});
  }
  std::vector<StableEdgeKey> parent_edges(features.segments.size());
  std::map<GlobalVertexId, std::vector<std::size_t>> parents_by_endpoint;
  for (std::size_t parent = 0; parent < features.segments.size(); parent++)
  {
    const auto &segment = features.segments[parent];
    parent_edges[parent] = {segment.mesh_vertices[0], segment.mesh_vertices[1]};
    std::sort(parent_edges[parent].begin(), parent_edges[parent].end());
    coordinate_ids.insert(parent_edges[parent].begin(), parent_edges[parent].end());
    parents_by_endpoint[parent_edges[parent][0]].push_back(parent);
    parents_by_endpoint[parent_edges[parent][1]].push_back(parent);
  }
  const auto coordinates =
      GatherStableVertexCoordinates(mesh, source_vertex_ids, coordinate_ids);

  std::vector<ChildFeatureSegment> children;
  children.reserve(features.segments.size() * 2);
  for (const auto &candidate : selected_edges)
  {
    std::set<std::size_t> possible_parents;
    for (GlobalVertexId endpoint : candidate)
    {
      const auto parent = parents_by_endpoint.find(endpoint);
      if (parent != parents_by_endpoint.end())
      {
        possible_parents.insert(parent->second.begin(), parent->second.end());
      }
    }
    ChildFeatureSegment child;
    std::vector<ChildFeatureSegment> matches;
    for (std::size_t parent : possible_parents)
    {
      if (IsChildSegment(candidate, parent_edges[parent], coordinates, child))
      {
        child.parent = parent;
        matches.push_back(child);
      }
    }
    if (matches.empty())
    {
      // A refinement implementation may split one old edge more than once in a single
      // call, leaving an interior child without an old endpoint. Fall back to a complete
      // geometric search in that uncommon case.
      for (std::size_t parent = 0; parent < parent_edges.size(); parent++)
      {
        if (possible_parents.find(parent) != possible_parents.end())
        {
          continue;
        }
        if (IsChildSegment(candidate, parent_edges[parent], coordinates, child))
        {
          child.parent = parent;
          matches.push_back(child);
        }
      }
    }
    MFEM_VERIFY(matches.size() <= 1,
                "A refined selected PEC edge maps to multiple parent singular segments!");
    if (!matches.empty())
    {
      children.push_back(matches.front());
    }
  }
  std::sort(children.begin(), children.end(),
            [](const auto &a, const auto &b) { return a.edge < b.edge; });
  MFEM_VERIFY(std::adjacent_find(children.begin(), children.end(),
                                 [](const auto &a, const auto &b)
                                 { return a.edge == b.edge; }) == children.end(),
              "Refined singular segment reconstruction produced duplicate child edges!");

  std::vector<std::vector<std::pair<double, double>>> parent_intervals(
      features.segments.size());
  for (const auto &child : children)
  {
    parent_intervals[child.parent].emplace_back(child.first, child.second);
  }
  constexpr double coverage_tolerance = 2.0e-9;
  for (std::size_t parent = 0; parent < parent_intervals.size(); parent++)
  {
    auto &intervals = parent_intervals[parent];
    MFEM_VERIFY(!intervals.empty(),
                "Adaptive refinement removed a three-dimensional singular feature "
                "segment!");
    std::sort(intervals.begin(), intervals.end());
    MFEM_VERIFY(intervals.front().first <= coverage_tolerance &&
                    intervals.back().second >= 1.0 - coverage_tolerance,
                "Refined child singular segments do not cover their parent edge!");
    double endpoint = intervals.front().second;
    for (std::size_t interval = 1; interval < intervals.size(); interval++)
    {
      MFEM_VERIFY(
          std::abs(intervals[interval].first - endpoint) <= coverage_tolerance,
          "Refined child singular segments contain a gap or overlap on their parent "
          "edge!");
      endpoint = intervals[interval].second;
    }
  }

  fem::singular::FeatureTopology updated;
  updated.features = features.features;
  for (auto &feature : updated.features)
  {
    feature.segments.clear();
    feature.mesh_vertices.clear();
  }
  updated.sheet_faces.reserve(selected_faces.size());
  for (const auto &[face, attribute] : selected_faces)
  {
    updated.sheet_faces.push_back(
        {static_cast<int>(updated.sheet_faces.size()), face, attribute});
  }

  updated.segments.reserve(children.size());
  std::map<GlobalVertexId, std::vector<std::size_t>> incident_segments;
  for (const auto &child : children)
  {
    const auto &parent = features.segments[child.parent];
    const std::size_t segment = updated.segments.size();
    updated.segments.push_back(
        {static_cast<int>(segment),
         {CheckedMeshEntityId(child.edge[0]), CheckedMeshEntityId(child.edge[1])},
         parent.feature,
         parent.boundary_attributes,
         parent.type});
    auto &feature = updated.features.at(parent.feature);
    feature.segments.push_back(segment);
    feature.mesh_vertices.push_back(CheckedMeshEntityId(child.edge[0]));
    feature.mesh_vertices.push_back(CheckedMeshEntityId(child.edge[1]));
    incident_segments[child.edge[0]].push_back(segment);
    incident_segments[child.edge[1]].push_back(segment);
  }
  for (auto &feature : updated.features)
  {
    std::sort(feature.mesh_vertices.begin(), feature.mesh_vertices.end());
    feature.mesh_vertices.erase(
        std::unique(feature.mesh_vertices.begin(), feature.mesh_vertices.end()),
        feature.mesh_vertices.end());
    MFEM_VERIFY(!feature.segments.empty(),
                "Adaptive refinement removed a complete straight singular feature!");
  }

  std::map<GlobalVertexId, fem::singular::FeatureVertex> old_vertices;
  for (const auto &vertex : features.vertices)
  {
    old_vertices.emplace(vertex.mesh_vertex, vertex);
  }
  updated.vertices.reserve(incident_segments.size());
  for (const auto &[id, segments] : incident_segments)
  {
    fem::singular::FeatureVertex vertex;
    const auto old = old_vertices.find(id);
    if (old != old_vertices.end())
    {
      vertex = old->second;
    }
    else
    {
      MFEM_VERIFY(segments.size() == 2 && updated.segments[segments[0]].feature ==
                                              updated.segments[segments[1]].feature,
                  "Refinement created a nonregular vertex on a straight singular "
                  "feature!");
      vertex.mesh_vertex = CheckedMeshEntityId(id);
      vertex.type = fem::singular::FeatureVertexType::REGULAR;
      vertex.nu = updated.features[updated.segments[segments[0]].feature].nu;
    }
    vertex.id = updated.vertices.size();
    vertex.mesh_vertex = CheckedMeshEntityId(id);
    vertex.segments = segments;
    vertex.features.clear();
    for (std::size_t segment : segments)
    {
      vertex.features.push_back(updated.segments[segment].feature);
      vertex.nu =
          std::min(vertex.nu, updated.features[updated.segments[segment].feature].nu);
    }
    std::sort(vertex.features.begin(), vertex.features.end());
    vertex.features.erase(std::unique(vertex.features.begin(), vertex.features.end()),
                          vertex.features.end());
    updated.vertices.push_back(std::move(vertex));
  }
  updated.elements.clear();
  features = std::move(updated);
}

mfem::Array<int> BuildSingularRefinementProtection(
    const mfem::ParMesh &mesh, const fem::singular::FeatureTopology &features,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids, bool *conforming,
    mfem::Array<int> *repair)
{
  return BuildSingularRefinementProtectionImpl(mesh, features, source_vertex_ids,
                                               conforming, repair);
}

mfem::Array<int> BuildSingularRefinementProtection(
    const mfem::ParMesh &mesh, const fem::singular::TriangleFeatureTopology &features,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids, bool *conforming,
    mfem::Array<int> *repair)
{
  return BuildSingularRefinementProtectionImpl(mesh, features, source_vertex_ids,
                                               conforming, repair);
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
  vertex_identity.Clear();
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
      ValidateSingularImpedanceFeatures(iodata, serial_sheet_features);
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
      ValidateSingularImpedanceFeatures(iodata, serial_line_features);
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
  if (parallel_mesh.Nonconforming())
  {
    vertex_identity.Observe(parallel_mesh, source_vertex_ids);
  }
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
  UpdateSingularSourceEntityIds(parallel_mesh, source_vertex_ids, source_element_ids,
                                vertex_identity);
  if (dimension == 2)
  {
    RebuildRefinedSingularFeatures(parallel_mesh,
                                   iodata.solver.singular_elements.attributes,
                                   source_vertex_ids, serial_line_features);
  }
  else
  {
    RebuildRefinedSingularFeatures(parallel_mesh,
                                   iodata.solver.singular_elements.attributes,
                                   source_vertex_ids, serial_sheet_features);
  }
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

mfem::Array<int> FullWaveSingularFeatures::GetRefinementProtection(
    const mfem::ParMesh &mesh, bool *conforming, mfem::Array<int> *repair) const
{
  MFEM_VERIFY(mesh.Dimension() == dimension && !source_vertex_ids.empty(),
              "Full-wave singular refinement protection is not initialized!");
  return dimension == 3
             ? BuildSingularRefinementProtection(mesh, local_sheet_features,
                                                 source_vertex_ids, conforming, repair)
             : BuildSingularRefinementProtection(mesh, local_line_features,
                                                 source_vertex_ids, conforming, repair);
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
