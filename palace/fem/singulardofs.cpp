// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singulardofs.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "utils/communication.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

namespace
{

struct ProvisionalDof
{
  DofKey key;
  HigherOrderBasis basis;
};

struct ProvisionalTriangleDof
{
  DofKey key;
  TriangleBasis basis;
};

constexpr std::size_t PackedDofKeySize = 21;
using PackedDofKey = std::array<GlobalVertexId, PackedDofKeySize>;

EntityKey MakeEntityKey(std::vector<GlobalVertexId> vertices)
{
  if (vertices.empty() || vertices.size() > 4)
  {
    throw std::runtime_error("A singular DOF has invalid geometric support!");
  }
  std::sort(vertices.begin(), vertices.end());
  if (std::adjacent_find(vertices.begin(), vertices.end()) != vertices.end())
  {
    throw std::runtime_error("A singular DOF entity contains duplicate vertices!");
  }
  EntityKey result;
  result.size = vertices.size();
  std::copy(vertices.begin(), vertices.end(), result.vertices.begin());
  return result;
}

GlobalVertexId GetGlobalVertex(const mfem::Element &element, int local,
                               const std::vector<GlobalVertexId> &serial_vertex_ids)
{
  if (local < 0 || local >= 4)
  {
    throw std::runtime_error("A singular basis has an invalid local vertex index!");
  }
  const int mesh_vertex = element.GetVertices()[local];
  if (mesh_vertex < 0 || mesh_vertex >= static_cast<int>(serial_vertex_ids.size()))
  {
    throw std::runtime_error("A singular basis references an invalid mesh vertex!");
  }
  return serial_vertex_ids[mesh_vertex];
}

GlobalVertexId GetTriangleGlobalVertex(const mfem::Element &element, int local,
                                       const std::vector<GlobalVertexId> &serial_vertex_ids)
{
  if (local < 0 || local >= 3)
  {
    throw std::runtime_error(
        "A triangular singular basis has an invalid local vertex index!");
  }
  const int mesh_vertex = element.GetVertices()[local];
  if (mesh_vertex < 0 || mesh_vertex >= static_cast<int>(serial_vertex_ids.size()))
  {
    throw std::runtime_error(
        "A triangular singular basis references an invalid mesh vertex!");
  }
  return serial_vertex_ids[mesh_vertex];
}

DofKey MakeDofKey(const mfem::Element &element, const HigherOrderBasis &basis,
                  const std::vector<GlobalVertexId> &serial_vertex_ids)
{
  std::vector<std::pair<GlobalVertexId, int>> support;
  for (int i = 0; i < 4; i++)
  {
    if (basis.interpolation_indices[i] > 0)
    {
      support.emplace_back(GetGlobalVertex(element, basis.nodes[i], serial_vertex_ids),
                           basis.interpolation_indices[i]);
    }
  }
  std::sort(support.begin(), support.end());
  if (support.empty() ||
      std::adjacent_find(support.begin(), support.end(), [](const auto &a, const auto &b)
                         { return a.first == b.first; }) != support.end())
  {
    throw std::runtime_error("A singular basis has invalid interpolation-point support!");
  }

  DofKey key;
  key.family = basis.family;
  key.order = basis.order;
  std::vector<GlobalVertexId> support_vertices;
  for (std::size_t i = 0; i < support.size(); i++)
  {
    support_vertices.push_back(support[i].first);
    key.interpolation_weights[i] = support[i].second;
  }
  key.support_entity = MakeEntityKey(std::move(support_vertices));

  switch (basis.family)
  {
    case HigherOrderBasisFamily::NODE_GRADIENT:
      key.singular_entity =
          MakeEntityKey({GetGlobalVertex(element, basis.nodes[0], serial_vertex_ids)});
      key.component_entity =
          MakeEntityKey({GetGlobalVertex(element, basis.nodes[0], serial_vertex_ids),
                         GetGlobalVertex(element, basis.nodes[1], serial_vertex_ids)});
      break;
    case HigherOrderBasisFamily::NODE_ROTATIONAL:
      key.singular_entity =
          MakeEntityKey({GetGlobalVertex(element, basis.nodes[0], serial_vertex_ids)});
      key.component_entity =
          MakeEntityKey({GetGlobalVertex(element, basis.nodes[1], serial_vertex_ids),
                         GetGlobalVertex(element, basis.nodes[2], serial_vertex_ids)});
      break;
    case HigherOrderBasisFamily::EDGE_GRADIENT:
      key.singular_entity =
          MakeEntityKey({GetGlobalVertex(element, basis.nodes[0], serial_vertex_ids),
                         GetGlobalVertex(element, basis.nodes[1], serial_vertex_ids)});
      key.component_entity =
          MakeEntityKey({GetGlobalVertex(element, basis.nodes[0], serial_vertex_ids),
                         GetGlobalVertex(element, basis.nodes[1], serial_vertex_ids),
                         GetGlobalVertex(element, basis.nodes[2], serial_vertex_ids)});
      break;
    case HigherOrderBasisFamily::EDGE_ROTATIONAL:
      key.singular_entity =
          MakeEntityKey({GetGlobalVertex(element, basis.nodes[0], serial_vertex_ids),
                         GetGlobalVertex(element, basis.nodes[1], serial_vertex_ids)});
      key.component_entity =
          MakeEntityKey({GetGlobalVertex(element, basis.nodes[2], serial_vertex_ids),
                         GetGlobalVertex(element, basis.nodes[3], serial_vertex_ids)});
      break;
  }
  return key;
}

void AppendBases(const mfem::Element &element, std::vector<HigherOrderBasis> basis,
                 const std::vector<GlobalVertexId> &serial_vertex_ids,
                 std::vector<ProvisionalDof> &h1, std::vector<ProvisionalDof> &nd)
{
  for (auto &entry : basis)
  {
    ProvisionalDof dof{MakeDofKey(element, entry, serial_vertex_ids), std::move(entry)};
    if (IsGradientFamily(dof.key.family))
    {
      h1.push_back(dof);
    }
    nd.push_back(std::move(dof));
  }
}

void AppendTriangleBasis(const mfem::Element &element, TriangleBasis basis,
                         const std::vector<GlobalVertexId> &serial_vertex_ids,
                         std::vector<ProvisionalTriangleDof> &h1,
                         std::vector<ProvisionalTriangleDof> &nd)
{
  const auto vertex = [&element, &serial_vertex_ids](int local)
  { return GetTriangleGlobalVertex(element, local, serial_vertex_ids); };
  const GlobalVertexId singular = vertex(basis.nodes[0]);

  DofKey key;
  key.family = basis.family;
  key.order = basis.order;
  key.singular_entity = MakeEntityKey({singular});
  if (basis.family == HigherOrderBasisFamily::NODE_GRADIENT)
  {
    const GlobalVertexId radial = vertex(basis.nodes[1]);
    key.support_entity = MakeEntityKey({singular, radial});
    key.component_entity = key.support_entity;
    key.interpolation_weights = {1, 1, 0, 0};
  }
  else if (basis.family == HigherOrderBasisFamily::NODE_ROTATIONAL)
  {
    const GlobalVertexId first = vertex(basis.nodes[1]);
    const GlobalVertexId second = vertex(basis.nodes[2]);
    key.support_entity = MakeEntityKey({singular, first, second});
    key.component_entity = MakeEntityKey({first, second});
    key.interpolation_weights = {1, 1, 1, 0};
  }
  else
  {
    throw std::invalid_argument("The triangular singular basis has an unsupported family!");
  }

  ProvisionalTriangleDof dof{std::move(key), std::move(basis)};
  if (IsGradientFamily(dof.key.family))
  {
    h1.push_back(dof);
  }
  nd.push_back(std::move(dof));
}

template <typename T>
void CheckUniqueLocalKeys(const T &dofs)
{
  std::set<DofKey> keys;
  for (const auto &dof : dofs)
  {
    if (!keys.insert(dof.key).second)
    {
      throw std::runtime_error(
          "Singular basis enumeration produced a duplicate local DOF!");
    }
  }
}

std::array<int, 4> CanonicalNodeTuple(const mfem::Element &element, int singular_vertex,
                                      const std::vector<GlobalVertexId> &serial_vertex_ids)
{
  int singular_local = -1;
  std::array<std::pair<GlobalVertexId, int>, 3> remaining;
  int next = 0;
  for (int local = 0; local < 4; local++)
  {
    const int mesh_vertex = element.GetVertices()[local];
    if (mesh_vertex == singular_vertex)
    {
      singular_local = local;
    }
    else
    {
      remaining[next++] = {serial_vertex_ids[mesh_vertex], local};
    }
  }
  if (singular_local < 0 || next != 3)
  {
    throw std::invalid_argument(
        "A rank-local singular node is not in its incident tetrahedron!");
  }
  std::sort(remaining.begin(), remaining.end());
  return {singular_local, remaining[0].second, remaining[1].second, remaining[2].second};
}

std::array<int, 4> CanonicalEdgeTuple(const mfem::Element &element,
                                      const std::array<GlobalVertexId, 2> &singular_edge,
                                      const std::vector<GlobalVertexId> &serial_vertex_ids)
{
  std::array<std::pair<GlobalVertexId, int>, 2> singular;
  std::array<std::pair<GlobalVertexId, int>, 2> remaining;
  int singular_count = 0, remaining_count = 0;
  for (int local = 0; local < 4; local++)
  {
    const GlobalVertexId vertex = serial_vertex_ids[element.GetVertices()[local]];
    if (vertex == singular_edge[0] || vertex == singular_edge[1])
    {
      singular[singular_count++] = {vertex, local};
    }
    else
    {
      remaining[remaining_count++] = {vertex, local};
    }
  }
  if (singular_count != 2 || remaining_count != 2)
  {
    throw std::invalid_argument(
        "A rank-local singular edge is not in its incident tetrahedron!");
  }
  std::sort(singular.begin(), singular.end());
  std::sort(remaining.begin(), remaining.end());
  return {singular[0].second, singular[1].second, remaining[0].second, remaining[1].second};
}

void ValidateFeatureTopology(const mfem::Mesh &mesh, const FeatureTopology &features,
                             const std::vector<GlobalVertexId> &serial_vertex_ids)
{
  if (mesh.Dimension() != 3 || mesh.SpaceDimension() != 3 ||
      features.elements.size() != static_cast<std::size_t>(mesh.GetNE()) ||
      serial_vertex_ids.size() != static_cast<std::size_t>(mesh.GetNV()))
  {
    throw std::invalid_argument(
        "Singular feature topology does not match the tetrahedral mesh!");
  }
  std::set<GlobalVertexId> unique_vertex_ids;
  for (GlobalVertexId vertex : serial_vertex_ids)
  {
    if (vertex < 0 || !unique_vertex_ids.insert(vertex).second)
    {
      throw std::invalid_argument(
          "Singular DOF enumeration requires unique stable vertex IDs!");
    }
  }
  for (std::size_t i = 0; i < features.vertices.size(); i++)
  {
    const auto &vertex = features.vertices[i];
    if (vertex.id != i || vertex.mesh_vertex < 0 || !std::isfinite(vertex.nu) ||
        vertex.nu <= 0.0 || vertex.nu >= 1.0)
    {
      throw std::invalid_argument("Singular node topology is inconsistent!");
    }
  }
  for (std::size_t i = 0; i < features.features.size(); i++)
  {
    const auto &feature = features.features[i];
    if (feature.id != i || !std::isfinite(feature.nu) || feature.nu <= 0.0 ||
        feature.nu >= 1.0)
    {
      throw std::invalid_argument("Straight singular feature is inconsistent!");
    }
  }
  for (std::size_t i = 0; i < features.segments.size(); i++)
  {
    const auto &segment = features.segments[i];
    if (segment.mesh_edge < 0 || segment.feature >= features.features.size())
    {
      throw std::invalid_argument("Singular edge topology is inconsistent!");
    }
  }

  mfem::Array<int> edge_vertices;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &tetrahedron = *mesh.GetElement(element);
    for (const auto &node : features.elements[element].nodes)
    {
      if (node.vertex >= features.vertices.size() || node.mesh_vertex < 0 ||
          node.mesh_vertex >= mesh.GetNV() ||
          serial_vertex_ids[node.mesh_vertex] !=
              features.vertices[node.vertex].mesh_vertex ||
          node.canonical_nodes !=
              CanonicalNodeTuple(tetrahedron, node.mesh_vertex, serial_vertex_ids))
      {
        throw std::invalid_argument(
            "Singular node incidence is inconsistent with its tetrahedron!");
      }
    }
    for (const auto &edge : features.elements[element].edges)
    {
      if (edge.segment >= features.segments.size() || edge.mesh_edge < 0 ||
          edge.mesh_edge >= mesh.GetNEdges())
      {
        throw std::invalid_argument(
            "Singular edge incidence is inconsistent with its tetrahedron!");
      }
      const auto &segment = features.segments[edge.segment];
      mesh.GetEdgeVertices(edge.mesh_edge, edge_vertices);
      if (edge.feature != segment.feature || edge_vertices.Size() != 2)
      {
        throw std::invalid_argument(
            "Singular edge incidence is inconsistent with its segment!");
      }
      std::array<GlobalVertexId, 2> stable_edge{serial_vertex_ids[edge_vertices[0]],
                                                serial_vertex_ids[edge_vertices[1]]};
      std::sort(stable_edge.begin(), stable_edge.end());
      const std::array<GlobalVertexId, 2> source_edge{segment.mesh_vertices[0],
                                                      segment.mesh_vertices[1]};
      if (stable_edge != source_edge ||
          edge.canonical_nodes !=
              CanonicalEdgeTuple(tetrahedron, stable_edge, serial_vertex_ids))
      {
        throw std::invalid_argument(
            "Singular edge incidence has inconsistent stable vertex IDs!");
      }
    }
  }
}

std::array<int, 3>
CanonicalTriangleNodeTuple(const mfem::Element &element, int singular_vertex,
                           const std::vector<GlobalVertexId> &serial_vertex_ids)
{
  int singular_local = -1;
  std::array<std::pair<GlobalVertexId, int>, 2> remaining;
  int next = 0;
  for (int local = 0; local < 3; local++)
  {
    const int mesh_vertex = element.GetVertices()[local];
    if (mesh_vertex == singular_vertex)
    {
      singular_local = local;
    }
    else
    {
      remaining[next++] = {serial_vertex_ids[mesh_vertex], local};
    }
  }
  if (singular_local < 0 || next != 2)
  {
    throw std::invalid_argument(
        "A rank-local singular tip is not in its incident triangle!");
  }
  std::sort(remaining.begin(), remaining.end());
  return {singular_local, remaining[0].second, remaining[1].second};
}

void ValidateTriangleFeatureTopology(const mfem::Mesh &mesh,
                                     const TriangleFeatureTopology &features,
                                     const std::vector<GlobalVertexId> &serial_vertex_ids)
{
  if (mesh.Dimension() != 2 || mesh.SpaceDimension() != 2 ||
      features.elements.size() != static_cast<std::size_t>(mesh.GetNE()) ||
      serial_vertex_ids.size() != static_cast<std::size_t>(mesh.GetNV()))
  {
    throw std::invalid_argument(
        "Singular line-tip topology does not match the triangular mesh!");
  }
  std::set<GlobalVertexId> unique_vertex_ids;
  for (GlobalVertexId vertex : serial_vertex_ids)
  {
    if (vertex < 0 || !unique_vertex_ids.insert(vertex).second)
    {
      throw std::invalid_argument(
          "Triangular singular DOF enumeration requires unique stable vertex IDs!");
    }
  }
  for (std::size_t i = 0; i < features.vertices.size(); i++)
  {
    const auto &vertex = features.vertices[i];
    if (vertex.id != i || vertex.mesh_vertex < 0 ||
        vertex.mesh_vertex >= static_cast<int>(serial_vertex_ids.size()) ||
        !std::isfinite(vertex.nu) || vertex.nu <= 0.0 || vertex.nu >= 1.0 ||
        (vertex.type != FeatureVertexType::ENDPOINT &&
         vertex.type != FeatureVertexType::CORNER))
    {
      throw std::invalid_argument("Singular line-tip topology is inconsistent!");
    }
  }
  for (const auto &segment : features.selected_segments)
  {
    if (segment.boundary_element < 0 || segment.mesh_edge < 0 ||
        segment.mesh_vertices[0] < 0 ||
        segment.mesh_vertices[0] >= segment.mesh_vertices[1] ||
        segment.boundary_attribute <= 0)
    {
      throw std::invalid_argument(
          "Selected singular line-segment topology is inconsistent!");
    }
  }
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &triangle = *mesh.GetElement(element);
    if (triangle.GetGeometryType() != mfem::Geometry::TRIANGLE)
    {
      throw std::invalid_argument(
          "Triangular singular DOF enumeration requires triangular elements!");
    }
    for (const auto &node : features.elements[element].nodes)
    {
      if (node.vertex >= features.vertices.size() || node.mesh_vertex < 0 ||
          node.mesh_vertex >= mesh.GetNV() ||
          serial_vertex_ids[node.mesh_vertex] !=
              features.vertices[node.vertex].mesh_vertex ||
          node.canonical_nodes !=
              CanonicalTriangleNodeTuple(triangle, node.mesh_vertex, serial_vertex_ids))
      {
        throw std::invalid_argument(
            "Singular tip incidence is inconsistent with its triangle!");
      }
    }
  }
}

PackedDofKey PackDofKey(const DofKey &key)
{
  PackedDofKey packed;
  std::size_t next = 0;
  packed[next++] = static_cast<int>(key.family);
  packed[next++] = key.order;
  const auto pack_entity = [&packed, &next](const EntityKey &entity)
  {
    if (entity.size == 0 || entity.size > entity.vertices.size())
    {
      throw std::invalid_argument(
          "Cannot communicate a singular DOF with an invalid entity!");
    }
    packed[next++] = static_cast<GlobalVertexId>(entity.size);
    for (GlobalVertexId vertex : entity.vertices)
    {
      packed[next++] = vertex;
    }
  };
  pack_entity(key.singular_entity);
  pack_entity(key.support_entity);
  pack_entity(key.component_entity);
  for (int weight : key.interpolation_weights)
  {
    packed[next++] = weight;
  }
  if (next != packed.size())
  {
    throw std::logic_error("Singular DOF key packing has an inconsistent size!");
  }
  return packed;
}

std::uint64_t HashPackedDofKey(const PackedDofKey &key)
{
  // Bytewise FNV-1a gives a decomposition-independent rendezvous rank without
  // depending on the host representation of the fixed-width integer key.
  std::uint64_t hash = 14695981039346656037ULL;
  for (GlobalVertexId value : key)
  {
    const auto word = static_cast<std::uint64_t>(value);
    for (int byte = 0; byte < 8; byte++)
    {
      hash ^= (word >> (8 * byte)) & 0xffULL;
      hash *= 1099511628211ULL;
    }
  }
  return hash;
}

template <std::size_t RecordSize>
using IntegerRecord = std::array<GlobalVertexId, RecordSize>;

template <std::size_t RecordSize>
std::vector<IntegerRecord<RecordSize>> ExchangeIntegerRecords(
    MPI_Comm comm, const std::vector<std::vector<IntegerRecord<RecordSize>>> &send_records)
{
  static_assert(RecordSize > 0);
  const int ranks = Mpi::Size(comm);
  bool valid = send_records.size() == static_cast<std::size_t>(ranks);
  std::vector<int> send_counts(ranks), receive_counts(ranks);
  std::int64_t send_total = 0;
  if (valid)
  {
    for (int rank = 0; rank < ranks; rank++)
    {
      if (send_records[rank].size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max() / RecordSize))
      {
        valid = false;
        break;
      }
      send_counts[rank] = static_cast<int>(send_records[rank].size() * RecordSize);
      send_total += send_counts[rank];
      if (send_total > std::numeric_limits<int>::max())
      {
        valid = false;
        break;
      }
    }
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::overflow_error(
        "Distributed singular DOF exchange exceeds MPI integer counts!");
  }

  MPI_Alltoall(send_counts.data(), 1, MPI_INT, receive_counts.data(), 1, MPI_INT, comm);
  std::vector<int> send_displacements(ranks), receive_displacements(ranks);
  std::int64_t receive_total = 0;
  valid = true;
  for (int rank = 0; rank < ranks; rank++)
  {
    if (receive_counts[rank] < 0 ||
        receive_counts[rank] % static_cast<int>(RecordSize) != 0 ||
        receive_total > std::numeric_limits<int>::max() - receive_counts[rank])
    {
      valid = false;
      break;
    }
    receive_displacements[rank] = static_cast<int>(receive_total);
    receive_total += receive_counts[rank];
  }
  int displacement = 0;
  for (int rank = 0; rank < ranks; rank++)
  {
    send_displacements[rank] = displacement;
    displacement += send_counts[rank];
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::overflow_error(
        "Distributed singular DOF exchange exceeds MPI integer counts!");
  }

  std::vector<GlobalVertexId> send_buffer(static_cast<std::size_t>(send_total));
  std::size_t next = 0;
  for (const auto &destination : send_records)
  {
    for (const auto &record : destination)
    {
      std::copy(record.begin(), record.end(), send_buffer.begin() + next);
      next += RecordSize;
    }
  }
  if (next != send_buffer.size())
  {
    throw std::logic_error(
        "Distributed singular DOF exchange packed inconsistent send dimensions!");
  }

  std::vector<GlobalVertexId> receive_buffer(static_cast<std::size_t>(receive_total));
  MPI_Alltoallv(send_buffer.data(), send_counts.data(), send_displacements.data(),
                mpi::DataType<GlobalVertexId>(), receive_buffer.data(),
                receive_counts.data(), receive_displacements.data(),
                mpi::DataType<GlobalVertexId>(), comm);

  std::vector<IntegerRecord<RecordSize>> result(static_cast<std::size_t>(receive_total) /
                                                RecordSize);
  for (std::size_t record = 0; record < result.size(); record++)
  {
    std::copy(receive_buffer.begin() + record * RecordSize,
              receive_buffer.begin() + (record + 1) * RecordSize, result[record].begin());
  }
  return result;
}

constexpr std::size_t KeyOccurrenceRecordSize = PackedDofKeySize + 2;
using KeyOccurrenceRecord = IntegerRecord<KeyOccurrenceRecordSize>;
using NumberingResponseRecord = IntegerRecord<3>;

PackedDofKey GetRecordKey(const KeyOccurrenceRecord &record)
{
  PackedDofKey key;
  std::copy(record.begin(), record.begin() + PackedDofKeySize, key.begin());
  return key;
}

KeyOccurrenceRecord MakeKeyOccurrenceRecord(const PackedDofKey &key, int rank, int local)
{
  KeyOccurrenceRecord record;
  std::copy(key.begin(), key.end(), record.begin());
  record[PackedDofKeySize] = rank;
  record[PackedDofKeySize + 1] = local;
  return record;
}

TrueDofMap BuildTrueDofMap(MPI_Comm comm, const std::vector<DofKey> &keys)
{
  bool locally_valid =
      std::is_sorted(keys.begin(), keys.end()) &&
      std::adjacent_find(keys.begin(), keys.end()) == keys.end() &&
      keys.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max());
  std::vector<PackedDofKey> packed_local(keys.size());
  PackedDofKey previous{};
  if (locally_valid)
  {
    try
    {
      for (std::size_t i = 0; i < keys.size(); i++)
      {
        const auto packed = PackDofKey(keys[i]);
        if (i > 0 && !(previous < packed))
        {
          locally_valid = false;
          break;
        }
        packed_local[i] = packed;
        previous = packed;
      }
    }
    catch (const std::exception &)
    {
      locally_valid = false;
    }
  }
  Mpi::GlobalAnd(1, &locally_valid, comm);
  if (!locally_valid)
  {
    throw std::invalid_argument(
        "Parallel singular DOF numbering received an invalid local key table!");
  }

  const int ranks = Mpi::Size(comm);
  const int my_rank = Mpi::Rank(comm);
  std::vector<std::vector<KeyOccurrenceRecord>> rendezvous_send(ranks);
  for (std::size_t local = 0; local < packed_local.size(); local++)
  {
    const int rendezvous = static_cast<int>(HashPackedDofKey(packed_local[local]) %
                                            static_cast<std::uint64_t>(ranks));
    rendezvous_send[rendezvous].push_back(
        MakeKeyOccurrenceRecord(packed_local[local], my_rank, static_cast<int>(local)));
  }
  const auto rendezvous_records = ExchangeIntegerRecords(comm, rendezvous_send);

  using Occurrence = std::pair<int, int>;
  std::map<PackedDofKey, std::vector<Occurrence>> rendezvous_groups;
  locally_valid = true;
  for (const auto &record : rendezvous_records)
  {
    const auto key = GetRecordKey(record);
    const GlobalVertexId origin_rank = record[PackedDofKeySize];
    const GlobalVertexId origin_local = record[PackedDofKeySize + 1];
    if (HashPackedDofKey(key) % static_cast<std::uint64_t>(ranks) !=
            static_cast<std::uint64_t>(my_rank) ||
        origin_rank < 0 || origin_rank >= ranks || origin_local < 0 ||
        origin_local > std::numeric_limits<int>::max())
    {
      locally_valid = false;
      break;
    }
    rendezvous_groups[key].emplace_back(static_cast<int>(origin_rank),
                                        static_cast<int>(origin_local));
  }
  std::vector<std::vector<KeyOccurrenceRecord>> owner_send(ranks);
  if (locally_valid)
  {
    for (auto &[key, occurrences] : rendezvous_groups)
    {
      std::sort(occurrences.begin(), occurrences.end());
      if (occurrences.empty() ||
          std::adjacent_find(occurrences.begin(), occurrences.end(),
                             [](const auto &left, const auto &right)
                             { return left.first == right.first; }) != occurrences.end())
      {
        locally_valid = false;
        break;
      }
      const int owner = occurrences.front().first;
      for (const auto &[origin_rank, origin_local] : occurrences)
      {
        owner_send[owner].push_back(
            MakeKeyOccurrenceRecord(key, origin_rank, origin_local));
      }
    }
  }
  Mpi::GlobalAnd(1, &locally_valid, comm);
  if (!locally_valid)
  {
    throw std::invalid_argument(
        "Distributed singular DOF rendezvous received inconsistent key occurrences!");
  }
  const auto owner_records = ExchangeIntegerRecords(comm, owner_send);

  std::map<PackedDofKey, std::vector<Occurrence>> owned_groups;
  locally_valid = true;
  for (const auto &record : owner_records)
  {
    const auto key = GetRecordKey(record);
    const GlobalVertexId origin_rank = record[PackedDofKeySize];
    const GlobalVertexId origin_local = record[PackedDofKeySize + 1];
    if (origin_rank < 0 || origin_rank >= ranks || origin_local < 0 ||
        origin_local > std::numeric_limits<int>::max())
    {
      locally_valid = false;
      break;
    }
    owned_groups[key].emplace_back(static_cast<int>(origin_rank),
                                   static_cast<int>(origin_local));
  }
  if (locally_valid)
  {
    for (auto &[key, occurrences] : owned_groups)
    {
      (void)key;
      std::sort(occurrences.begin(), occurrences.end());
      if (occurrences.empty() || occurrences.front().first != my_rank ||
          std::adjacent_find(occurrences.begin(), occurrences.end(),
                             [](const auto &left, const auto &right)
                             { return left.first == right.first; }) != occurrences.end())
      {
        locally_valid = false;
        break;
      }
    }
  }
  locally_valid = locally_valid &&
                  owned_groups.size() <=
                      static_cast<std::size_t>(std::numeric_limits<HYPRE_BigInt>::max());
  Mpi::GlobalAnd(1, &locally_valid, comm);
  if (!locally_valid)
  {
    throw std::invalid_argument(
        "Distributed singular DOF ownership exchange received inconsistent records!");
  }

  const std::array<HYPRE_BigInt, 2> local_sizes{
      static_cast<HYPRE_BigInt>(keys.size()),
      static_cast<HYPRE_BigInt>(owned_groups.size())};
  std::vector<HYPRE_BigInt> gathered_sizes(2 * ranks);
  Mpi::Allgather(2, local_sizes.data(), gathered_sizes.data(), comm);
  std::vector<HYPRE_BigInt> owned_offsets(ranks), owned_sizes(ranks);
  TrueDofMap result;
  for (int rank = 0; rank < ranks; rank++)
  {
    const HYPRE_BigInt rank_local_size = gathered_sizes[2 * rank];
    const HYPRE_BigInt rank_owned_size = gathered_sizes[2 * rank + 1];
    if (rank_local_size < 0 || rank_owned_size < 0 ||
        rank_local_size >
            std::numeric_limits<HYPRE_BigInt>::max() - result.global_local_size ||
        rank_owned_size > std::numeric_limits<HYPRE_BigInt>::max() - result.global_size)
    {
      locally_valid = false;
      break;
    }
    if (rank < my_rank)
    {
      result.local_offset += rank_local_size;
    }
    result.global_local_size += rank_local_size;
    owned_offsets[rank] = result.global_size;
    owned_sizes[rank] = rank_owned_size;
    result.global_size += rank_owned_size;
  }
  Mpi::GlobalAnd(1, &locally_valid, comm);
  if (!locally_valid)
  {
    throw std::overflow_error("Parallel singular DOF count overflow!");
  }
  result.local_size = local_sizes[0];
  result.owned_offset = owned_offsets[my_rank];
  result.owned_size = owned_sizes[my_rank];

  static_assert(std::numeric_limits<GlobalVertexId>::is_signed);
  static_assert(std::numeric_limits<HYPRE_BigInt>::is_signed);
  static_assert(sizeof(GlobalVertexId) >= sizeof(HYPRE_BigInt));
  std::vector<std::vector<NumberingResponseRecord>> response_send(ranks);
  HYPRE_BigInt owner_position = 0;
  for (const auto &[key, occurrences] : owned_groups)
  {
    (void)key;
    const HYPRE_BigInt true_dof = result.owned_offset + owner_position++;
    for (const auto &[origin_rank, origin_local] : occurrences)
    {
      response_send[origin_rank].push_back(
          {origin_local, my_rank, static_cast<GlobalVertexId>(true_dof)});
    }
  }
  if (owner_position != result.owned_size)
  {
    throw std::logic_error(
        "Distributed singular DOF ownership produced inconsistent local dimensions!");
  }
  const auto responses = ExchangeIntegerRecords(comm, response_send);

  result.owner.resize(keys.size());
  result.local_to_true.resize(keys.size());
  std::vector<bool> assigned(keys.size());
  locally_valid = responses.size() == keys.size();
  if (locally_valid)
  {
    for (const auto &response : responses)
    {
      const GlobalVertexId local = response[0];
      const GlobalVertexId owner = response[1];
      const GlobalVertexId true_dof = response[2];
      if (local < 0 || local >= static_cast<GlobalVertexId>(keys.size()) || owner < 0 ||
          owner >= ranks || owner > my_rank || true_dof < 0 ||
          true_dof >= result.global_size ||
          true_dof < owned_offsets[static_cast<int>(owner)] ||
          true_dof >= owned_offsets[static_cast<int>(owner)] +
                          owned_sizes[static_cast<int>(owner)] ||
          assigned[static_cast<std::size_t>(local)])
      {
        locally_valid = false;
        break;
      }
      assigned[static_cast<std::size_t>(local)] = true;
      result.owner[static_cast<std::size_t>(local)] = static_cast<int>(owner);
      result.local_to_true[static_cast<std::size_t>(local)] =
          static_cast<HYPRE_BigInt>(true_dof);
    }
    locally_valid = locally_valid && std::all_of(assigned.begin(), assigned.end(),
                                                 [](bool value) { return value; });
  }
  Mpi::GlobalAnd(1, &locally_valid, comm);
  if (!locally_valid)
  {
    throw std::runtime_error(
        "Distributed singular DOF numbering received inconsistent owner responses!");
  }
  return result;
}

template <typename Topology>
ParallelDofNumbering BuildParallelDofNumberingImpl(MPI_Comm comm, const Topology &topology)
{
  bool locally_valid = topology.h1_to_nd.size() == topology.h1_dofs.size();
  Mpi::GlobalAnd(1, &locally_valid, comm);
  if (!locally_valid)
  {
    throw std::invalid_argument(
        "Parallel singular DOF numbering requires a complete H1-to-ND map!");
  }

  ParallelDofNumbering result;
  result.h1 = BuildTrueDofMap(comm, topology.h1_dofs);
  result.nd = BuildTrueDofMap(comm, topology.nd_dofs);
  result.h1_to_nd_true.resize(topology.h1_dofs.size());
  locally_valid = true;
  for (std::size_t h1 = 0; h1 < topology.h1_dofs.size(); h1++)
  {
    const std::size_t nd = topology.h1_to_nd[h1];
    if (nd >= topology.nd_dofs.size() || !(topology.h1_dofs[h1] == topology.nd_dofs[nd]) ||
        result.h1.owner[h1] != result.nd.owner[nd])
    {
      locally_valid = false;
      break;
    }
    result.h1_to_nd_true[h1] = result.nd.local_to_true[nd];
  }
  Mpi::GlobalAnd(1, &locally_valid, comm);
  if (!locally_valid)
  {
    throw std::invalid_argument("Parallel singular H1-to-ND topology is inconsistent!");
  }
  return result;
}

template <typename Topology>
mfem::Array<int> GetEssentialH1TrueDofsImpl(MPI_Comm comm,
                                            const std::set<EntityKey> &essential_support,
                                            const Topology &topology,
                                            const ParallelDofNumbering &parallel_numbering)
{
  const auto &numbering = parallel_numbering.h1;
  if (numbering.owner.size() != topology.h1_dofs.size() ||
      numbering.local_to_true.size() != topology.h1_dofs.size() ||
      numbering.owned_offset < 0 || numbering.owned_size < 0 ||
      numbering.owned_size > std::numeric_limits<int>::max())
  {
    throw std::invalid_argument(
        "Essential singular H1 classification received inconsistent DOF numbering!");
  }

  mfem::Array<int> result;
  for (std::size_t local = 0; local < topology.h1_dofs.size(); local++)
  {
    if (essential_support.find(topology.h1_dofs[local].support_entity) ==
            essential_support.end() ||
        numbering.owner[local] != Mpi::Rank(comm))
    {
      continue;
    }
    const HYPRE_BigInt true_dof = numbering.local_to_true[local] - numbering.owned_offset;
    if (true_dof < 0 || true_dof >= numbering.owned_size)
    {
      throw std::invalid_argument(
          "An essential singular H1 DOF is outside its owner partition!");
    }
    result.Append(static_cast<int>(true_dof));
  }
  result.Sort();
  if (std::adjacent_find(result.begin(), result.end()) != result.end())
  {
    throw std::logic_error("Essential singular H1 true DOFs are not unique!");
  }
  return result;
}

template <typename Topology>
mfem::Array<int> GetEssentialNDTrueDofsImpl(MPI_Comm comm,
                                            const std::set<EntityKey> &essential_support,
                                            const Topology &topology,
                                            const ParallelDofNumbering &parallel_numbering)
{
  const auto &numbering = parallel_numbering.nd;
  if (numbering.owner.size() != topology.nd_dofs.size() ||
      numbering.local_to_true.size() != topology.nd_dofs.size() ||
      numbering.owned_offset < 0 || numbering.owned_size < 0 ||
      numbering.owned_size > std::numeric_limits<int>::max())
  {
    throw std::invalid_argument(
        "Essential singular ND classification received inconsistent DOF numbering!");
  }

  mfem::Array<int> result;
  for (std::size_t local = 0; local < topology.nd_dofs.size(); local++)
  {
    if (essential_support.find(topology.nd_dofs[local].support_entity) ==
            essential_support.end() ||
        numbering.owner[local] != Mpi::Rank(comm))
    {
      continue;
    }
    const HYPRE_BigInt true_dof = numbering.local_to_true[local] - numbering.owned_offset;
    if (true_dof < 0 || true_dof >= numbering.owned_size)
    {
      throw std::invalid_argument(
          "An essential singular ND DOF is outside its owner partition!");
    }
    result.Append(static_cast<int>(true_dof));
  }
  result.Sort();
  if (std::adjacent_find(result.begin(), result.end()) != result.end())
  {
    throw std::logic_error("Essential singular ND true DOFs are not unique!");
  }
  return result;
}

std::set<EntityKey>
GetSelectedTriangleSegments(const TriangleFeatureTopology &features,
                            const std::set<int> *essential_boundary_attributes = nullptr)
{
  std::set<EntityKey> selected_segments;
  std::set<EntityKey> unique_segments;
  for (const auto &segment : features.selected_segments)
  {
    const auto key = MakeEntityKey({segment.mesh_vertices[0], segment.mesh_vertices[1]});
    if (segment.boundary_attribute <= 0 || segment.mesh_vertices[0] < 0 ||
        segment.mesh_vertices[0] >= segment.mesh_vertices[1] ||
        !unique_segments.insert(key).second)
    {
      throw std::invalid_argument(
          "Essential triangular classification received invalid selected segments!");
    }
    if (!essential_boundary_attributes ||
        essential_boundary_attributes->count(segment.boundary_attribute) > 0)
    {
      selected_segments.insert(key);
    }
  }
  return selected_segments;
}

std::vector<std::array<GlobalVertexId, 3>>
GetSelectedSheetFaces(const FeatureTopology &features,
                      const std::set<int> *essential_boundary_attributes = nullptr)
{
  std::set<std::array<GlobalVertexId, 3>> unique_faces;
  std::vector<std::array<GlobalVertexId, 3>> selected_faces;
  for (const auto &face : features.sheet_faces)
  {
    if (face.boundary_attribute <= 0 ||
        !std::is_sorted(face.mesh_vertices.begin(), face.mesh_vertices.end()) ||
        std::adjacent_find(face.mesh_vertices.begin(), face.mesh_vertices.end()) !=
            face.mesh_vertices.end() ||
        !unique_faces.insert(face.mesh_vertices).second)
    {
      throw std::invalid_argument(
          "Essential singular classification received invalid sheet faces!");
    }
    if (!essential_boundary_attributes ||
        essential_boundary_attributes->count(face.boundary_attribute) > 0)
    {
      selected_faces.push_back(face.mesh_vertices);
    }
  }
  return selected_faces;
}

template <typename Dof>
bool HasSheetTrace(const Dof &dof,
                   const std::vector<std::array<GlobalVertexId, 3>> &sheet_faces)
{
  const auto &support = dof.support_entity;
  if (support.size == 0 || support.size > 4)
  {
    throw std::invalid_argument(
        "An enriched singular DOF has invalid trace-support topology!");
  }
  return std::any_of(
      sheet_faces.begin(), sheet_faces.end(),
      [&support](const auto &face)
      {
        return support.size <= face.size() &&
               std::all_of(
                   support.vertices.begin(), support.vertices.begin() + support.size,
                   [&face](GlobalVertexId vertex)
                   { return std::binary_search(face.begin(), face.end(), vertex); });
      });
}

}  // namespace

bool EntityKey::operator==(const EntityKey &other) const
{
  return size == other.size && vertices == other.vertices;
}

bool EntityKey::operator<(const EntityKey &other) const
{
  return std::tie(size, vertices) < std::tie(other.size, other.vertices);
}

bool DofKey::operator==(const DofKey &other) const
{
  return family == other.family && order == other.order &&
         singular_entity == other.singular_entity &&
         support_entity == other.support_entity &&
         component_entity == other.component_entity &&
         interpolation_weights == other.interpolation_weights;
}

bool DofKey::operator<(const DofKey &other) const
{
  return std::tie(family, order, singular_entity, support_entity, component_entity,
                  interpolation_weights) <
         std::tie(other.family, other.order, other.singular_entity, other.support_entity,
                  other.component_entity, other.interpolation_weights);
}

bool IsGradientFamily(HigherOrderBasisFamily family)
{
  return family == HigherOrderBasisFamily::NODE_GRADIENT ||
         family == HigherOrderBasisFamily::EDGE_GRADIENT;
}

std::vector<std::vector<std::size_t>>
BuildH1DofFeatureMembership(const FeatureTopology &features, const DofTopology &topology)
{
  std::map<GlobalVertexId, std::vector<std::size_t>> vertex_features;
  for (const auto &vertex : features.vertices)
  {
    std::vector<std::size_t> membership = vertex.features;
    std::sort(membership.begin(), membership.end());
    const bool valid =
        vertex.mesh_vertex >= 0 && !membership.empty() &&
        std::adjacent_find(membership.begin(), membership.end()) == membership.end() &&
        std::all_of(membership.begin(), membership.end(),
                    [&features, &vertex](std::size_t feature)
                    {
                      return feature < features.features.size() &&
                             std::binary_search(
                                 features.features[feature].mesh_vertices.begin(),
                                 features.features[feature].mesh_vertices.end(),
                                 vertex.mesh_vertex);
                    });
    if (!valid || !vertex_features
                       .emplace(static_cast<GlobalVertexId>(vertex.mesh_vertex),
                                std::move(membership))
                       .second)
    {
      throw std::invalid_argument("Singular node-to-feature membership is inconsistent!");
    }
  }

  std::map<std::array<GlobalVertexId, 2>, std::size_t> edge_features;
  for (const auto &segment : features.segments)
  {
    const std::array<GlobalVertexId, 2> edge{
        static_cast<GlobalVertexId>(segment.mesh_vertices[0]),
        static_cast<GlobalVertexId>(segment.mesh_vertices[1])};
    if (segment.feature >= features.features.size() || edge[0] < 0 || edge[0] >= edge[1])
    {
      throw std::invalid_argument("Singular edge-to-feature membership is inconsistent!");
    }
    const auto [iterator, inserted] = edge_features.emplace(edge, segment.feature);
    if (!inserted && iterator->second != segment.feature)
    {
      throw std::invalid_argument(
          "One singular edge belongs to multiple straight features!");
    }
  }

  std::vector<std::vector<std::size_t>> result;
  result.reserve(topology.h1_dofs.size());
  for (const auto &dof : topology.h1_dofs)
  {
    if (!IsGradientFamily(dof.family))
    {
      throw std::invalid_argument("The singular H1 topology contains a rotational basis!");
    }
    if (dof.family == HigherOrderBasisFamily::NODE_GRADIENT)
    {
      if (dof.singular_entity.size != 1)
      {
        throw std::invalid_argument(
            "A node-gradient H1 basis has an invalid singular entity!");
      }
      const auto membership = vertex_features.find(dof.singular_entity.vertices[0]);
      if (membership == vertex_features.end())
      {
        throw std::invalid_argument(
            "A node-gradient H1 basis is absent from the feature topology!");
      }
      result.push_back(membership->second);
    }
    else
    {
      if (dof.singular_entity.size != 2)
      {
        throw std::invalid_argument(
            "An edge-gradient H1 basis has an invalid singular entity!");
      }
      const std::array<GlobalVertexId, 2> edge{dof.singular_entity.vertices[0],
                                               dof.singular_entity.vertices[1]};
      const auto membership = edge_features.find(edge);
      if (membership == edge_features.end())
      {
        throw std::invalid_argument(
            "An edge-gradient H1 basis is absent from the feature topology!");
      }
      result.push_back({membership->second});
    }
  }
  return result;
}

std::vector<std::vector<std::size_t>>
BuildTriangleH1DofFeatureMembership(const TriangleFeatureTopology &features,
                                    const TriangleDofTopology &topology)
{
  std::map<GlobalVertexId, std::size_t> tip_features;
  for (std::size_t feature = 0; feature < features.vertices.size(); feature++)
  {
    const auto &vertex = features.vertices[feature];
    if (vertex.id != feature || vertex.mesh_vertex < 0 ||
        vertex.selected_segments.empty() || vertex.selected_segments.size() > 2 ||
        std::any_of(vertex.selected_segments.begin(), vertex.selected_segments.end(),
                    [&features](std::size_t segment)
                    { return segment >= features.selected_segments.size(); }) ||
        !tip_features.emplace(static_cast<GlobalVertexId>(vertex.mesh_vertex), feature)
             .second)
    {
      throw std::invalid_argument(
          "Triangular singular tip-to-feature membership is inconsistent!");
    }
  }

  std::vector<std::vector<std::size_t>> result;
  result.reserve(topology.h1_dofs.size());
  for (const auto &dof : topology.h1_dofs)
  {
    if (dof.family != HigherOrderBasisFamily::NODE_GRADIENT ||
        dof.singular_entity.size != 1)
    {
      throw std::invalid_argument(
          "The triangular singular H1 topology contains an invalid basis!");
    }
    const auto feature = tip_features.find(dof.singular_entity.vertices[0]);
    if (feature == tip_features.end())
    {
      throw std::invalid_argument(
          "A triangular H1 basis is absent from the line-tip topology!");
    }
    result.push_back({feature->second});
  }
  return result;
}

std::vector<GlobalVertexId>
MapPartitionedSerialVertexIds(const mfem::Mesh &serial_mesh,
                              const mfem::ParMesh &parallel_mesh, const int *partitioning)
{
  if (!partitioning)
  {
    throw std::invalid_argument("Serial vertex-ID mapping requires an element partition!");
  }
  if (serial_mesh.Dimension() != parallel_mesh.Dimension() ||
      serial_mesh.SpaceDimension() != parallel_mesh.SpaceDimension())
  {
    throw std::invalid_argument("Serial and parallel meshes have incompatible dimensions!");
  }

  std::vector<GlobalVertexId> result(parallel_mesh.GetNV(), -1);
  int local_element = 0;
  for (int serial_element = 0; serial_element < serial_mesh.GetNE(); serial_element++)
  {
    const int rank = partitioning[serial_element];
    if (rank < 0 || rank >= parallel_mesh.GetNRanks())
    {
      throw std::invalid_argument("Element partition contains an invalid rank!");
    }
    if (rank != parallel_mesh.GetMyRank())
    {
      continue;
    }
    if (local_element >= parallel_mesh.GetNE())
    {
      throw std::invalid_argument(
          "Parallel mesh element ordering does not match its serial partition!");
    }

    const auto &serial = *serial_mesh.GetElement(serial_element);
    const int parallel_element = local_element++;
    const auto &parallel = *parallel_mesh.GetElement(parallel_element);
    if (serial.GetGeometryType() != parallel.GetGeometryType() ||
        serial.GetNVertices() != parallel.GetNVertices())
    {
      throw std::invalid_argument(
          "Parallel element does not match its serial source element!");
    }
    const auto *serial_reference_vertices =
        mfem::Geometries.GetVertices(serial.GetGeometryType());
    const auto *parallel_reference_vertices =
        mfem::Geometries.GetVertices(parallel.GetGeometryType());
    if (!serial_reference_vertices || !parallel_reference_vertices ||
        serial_reference_vertices->GetNPoints() != serial.GetNVertices() ||
        parallel_reference_vertices->GetNPoints() != parallel.GetNVertices())
    {
      throw std::runtime_error(
          "Unable to obtain reference vertices for a serial-to-parallel element map!");
    }
    mfem::IsoparametricTransformation serial_transformation;
    mfem::IsoparametricTransformation parallel_transformation;
    serial_mesh.GetElementTransformation(serial_element, &serial_transformation);
    parallel_mesh.GetElementTransformation(parallel_element, &parallel_transformation);
    for (int i = 0; i < serial.GetNVertices(); i++)
    {
      const int serial_vertex = serial.GetVertices()[i];
      const int parallel_vertex = parallel.GetVertices()[i];
      if (serial_vertex < 0 || serial_vertex >= serial_mesh.GetNV() ||
          parallel_vertex < 0 || parallel_vertex >= parallel_mesh.GetNV())
      {
        throw std::invalid_argument(
            "Serial-to-parallel element map contains an invalid vertex!");
      }
      auto &mapped = result[parallel_vertex];
      if (mapped >= 0 && mapped != serial_vertex)
      {
        throw std::invalid_argument(
            "A parallel vertex maps to inconsistent serial vertices!");
      }
      mapped = serial_vertex;

      mfem::Vector serial_coordinate(serial_mesh.SpaceDimension());
      mfem::Vector parallel_coordinate(parallel_mesh.SpaceDimension());
      serial_transformation.Transform(serial_reference_vertices->IntPoint(i),
                                      serial_coordinate);
      parallel_transformation.Transform(parallel_reference_vertices->IntPoint(i),
                                        parallel_coordinate);
      for (int d = 0; d < serial_mesh.SpaceDimension(); d++)
      {
        const double scale = 1.0 + std::max(std::abs(serial_coordinate[d]),
                                            std::abs(parallel_coordinate[d]));
        if (std::abs(serial_coordinate[d] - parallel_coordinate[d]) >
            32.0 * std::numeric_limits<double>::epsilon() * scale)
        {
          std::ostringstream message;
          message.precision(17);
          message << "Mapped serial and parallel vertices have different coordinates "
                  << "(serial element " << serial_element << ", parallel element "
                  << parallel_element << ", local vertex " << i << ", component " << d
                  << ": " << serial_coordinate[d] << " != " << parallel_coordinate[d]
                  << ")!";
          throw std::invalid_argument(message.str());
        }
      }
    }
  }
  if (local_element != parallel_mesh.GetNE() ||
      std::find(result.begin(), result.end(), -1) != result.end())
  {
    throw std::invalid_argument(
        "Parallel mesh is not the supplied partition of the serial mesh!");
  }
  return result;
}

DofTopology BuildLocalDofTopology(const mfem::Mesh &mesh, const FeatureTopology &features,
                                  const std::vector<GlobalVertexId> &serial_vertex_ids,
                                  int order)
{
  ValidateFeatureTopology(mesh, features, serial_vertex_ids);
  if (order < 1)
  {
    throw std::invalid_argument("Singular basis order must be positive!");
  }

  std::vector<std::vector<ProvisionalDof>> element_h1(mesh.GetNE());
  std::vector<std::vector<ProvisionalDof>> element_nd(mesh.GetNE());
  std::set<DofKey> h1_keys, nd_keys;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &tetrahedron = *mesh.GetElement(element);
    if (tetrahedron.GetGeometryType() != mfem::Geometry::TETRAHEDRON)
    {
      throw std::invalid_argument(
          "Singular DOF enumeration requires tetrahedral elements!");
    }

    auto &h1 = element_h1[element];
    auto &nd = element_nd[element];
    for (const auto &node : features.elements[element].nodes)
    {
      const double nu = features.vertices[node.vertex].nu;
      AppendBases(tetrahedron,
                  EnumerateHigherOrderNodeGradientBases(node.canonical_nodes, order, nu),
                  serial_vertex_ids, h1, nd);
      AppendBases(tetrahedron,
                  EnumerateHigherOrderNodeRotationalBases(node.canonical_nodes, order, nu),
                  serial_vertex_ids, h1, nd);
    }
    for (const auto &edge : features.elements[element].edges)
    {
      const double nu = features.features[edge.feature].nu;
      AppendBases(tetrahedron,
                  EnumerateHigherOrderEdgeGradientBases(edge.canonical_nodes, order, nu),
                  serial_vertex_ids, h1, nd);
      AppendBases(tetrahedron,
                  EnumerateHigherOrderEdgeRotationalBases(edge.canonical_nodes, order, nu),
                  serial_vertex_ids, h1, nd);
    }
    CheckUniqueLocalKeys(h1);
    CheckUniqueLocalKeys(nd);
    for (const auto &dof : h1)
    {
      h1_keys.insert(dof.key);
    }
    for (const auto &dof : nd)
    {
      nd_keys.insert(dof.key);
    }
  }

  DofTopology result;
  result.h1_dofs.assign(h1_keys.begin(), h1_keys.end());
  result.nd_dofs.assign(nd_keys.begin(), nd_keys.end());
  result.elements.resize(mesh.GetNE());
  std::map<DofKey, std::size_t> h1_numbering, nd_numbering;
  for (std::size_t i = 0; i < result.h1_dofs.size(); i++)
  {
    h1_numbering.emplace(result.h1_dofs[i], i);
  }
  for (std::size_t i = 0; i < result.nd_dofs.size(); i++)
  {
    nd_numbering.emplace(result.nd_dofs[i], i);
  }

  result.h1_to_nd.resize(result.h1_dofs.size());
  for (std::size_t i = 0; i < result.h1_dofs.size(); i++)
  {
    const auto nd = nd_numbering.find(result.h1_dofs[i]);
    if (nd == nd_numbering.end())
    {
      throw std::runtime_error(
          "An enriched H1 gradient basis is absent from the ND space!");
    }
    result.h1_to_nd[i] = nd->second;
  }

  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &element_map = result.elements[element];
    for (auto &dof : element_h1[element])
    {
      element_map.h1.push_back({h1_numbering.at(dof.key), std::move(dof.basis)});
    }
    for (auto &dof : element_nd[element])
    {
      element_map.nd.push_back({nd_numbering.at(dof.key), std::move(dof.basis)});
    }
  }
  return result;
}

DofTopology BuildSerialDofTopology(const mfem::Mesh &mesh, const FeatureTopology &features,
                                   int order)
{
  const auto *parallel_mesh = dynamic_cast<const mfem::ParMesh *>(&mesh);
  if (parallel_mesh && parallel_mesh->GetNRanks() > 1)
  {
    throw std::invalid_argument(
        "Serial singular DOF enumeration cannot use a multi-rank ParMesh!");
  }
  std::vector<GlobalVertexId> serial_vertex_ids(mesh.GetNV());
  for (int vertex = 0; vertex < mesh.GetNV(); vertex++)
  {
    serial_vertex_ids[vertex] = vertex;
  }
  return BuildLocalDofTopology(mesh, features, serial_vertex_ids, order);
}

TriangleDofTopology BuildLocalTriangleDofTopology(
    const mfem::Mesh &mesh, const TriangleFeatureTopology &features,
    const std::vector<GlobalVertexId> &serial_vertex_ids, int order)
{
  ValidateTriangleFeatureTopology(mesh, features, serial_vertex_ids);
  if (order != 1)
  {
    throw std::invalid_argument(
        "The triangular singular basis currently supports only order one!");
  }

  std::vector<std::vector<ProvisionalTriangleDof>> element_h1(mesh.GetNE());
  std::vector<std::vector<ProvisionalTriangleDof>> element_nd(mesh.GetNE());
  std::set<DofKey> h1_keys, nd_keys;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &triangle = *mesh.GetElement(element);
    auto &h1 = element_h1[element];
    auto &nd = element_nd[element];
    for (const auto &node : features.elements[element].nodes)
    {
      const double nu = features.vertices[node.vertex].nu;
      const int i = node.canonical_nodes[0];
      const int j = node.canonical_nodes[1];
      const int k = node.canonical_nodes[2];
      AppendTriangleBasis(
          triangle,
          TriangleBasis{HigherOrderBasisFamily::NODE_GRADIENT, {i, j, k}, order, nu},
          serial_vertex_ids, h1, nd);
      AppendTriangleBasis(
          triangle,
          TriangleBasis{HigherOrderBasisFamily::NODE_GRADIENT, {i, k, j}, order, nu},
          serial_vertex_ids, h1, nd);
      AppendTriangleBasis(
          triangle,
          TriangleBasis{HigherOrderBasisFamily::NODE_ROTATIONAL, {i, j, k}, order, nu},
          serial_vertex_ids, h1, nd);
    }
    CheckUniqueLocalKeys(h1);
    CheckUniqueLocalKeys(nd);
    for (const auto &dof : h1)
    {
      h1_keys.insert(dof.key);
    }
    for (const auto &dof : nd)
    {
      nd_keys.insert(dof.key);
    }
  }

  TriangleDofTopology result;
  result.h1_dofs.assign(h1_keys.begin(), h1_keys.end());
  result.nd_dofs.assign(nd_keys.begin(), nd_keys.end());
  result.elements.resize(mesh.GetNE());
  std::map<DofKey, std::size_t> h1_numbering, nd_numbering;
  for (std::size_t i = 0; i < result.h1_dofs.size(); i++)
  {
    h1_numbering.emplace(result.h1_dofs[i], i);
  }
  for (std::size_t i = 0; i < result.nd_dofs.size(); i++)
  {
    nd_numbering.emplace(result.nd_dofs[i], i);
  }
  result.h1_to_nd.resize(result.h1_dofs.size());
  for (std::size_t i = 0; i < result.h1_dofs.size(); i++)
  {
    const auto nd = nd_numbering.find(result.h1_dofs[i]);
    if (nd == nd_numbering.end())
    {
      throw std::runtime_error(
          "A triangular enriched H1 gradient basis is absent from the ND space!");
    }
    result.h1_to_nd[i] = nd->second;
  }
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &element_map = result.elements[element];
    for (auto &dof : element_h1[element])
    {
      element_map.h1.push_back({h1_numbering.at(dof.key), std::move(dof.basis)});
    }
    for (auto &dof : element_nd[element])
    {
      element_map.nd.push_back({nd_numbering.at(dof.key), std::move(dof.basis)});
    }
  }
  return result;
}

TriangleDofTopology BuildSerialTriangleDofTopology(const mfem::Mesh &mesh,
                                                   const TriangleFeatureTopology &features,
                                                   int order)
{
  const auto *parallel_mesh = dynamic_cast<const mfem::ParMesh *>(&mesh);
  if (parallel_mesh && parallel_mesh->GetNRanks() > 1)
  {
    throw std::invalid_argument(
        "Serial triangular singular DOF enumeration cannot use a multi-rank ParMesh!");
  }
  std::vector<GlobalVertexId> serial_vertex_ids(mesh.GetNV());
  for (int vertex = 0; vertex < mesh.GetNV(); vertex++)
  {
    serial_vertex_ids[vertex] = vertex;
  }
  return BuildLocalTriangleDofTopology(mesh, features, serial_vertex_ids, order);
}

ParallelDofNumbering BuildParallelDofNumbering(MPI_Comm comm, const DofTopology &topology)
{
  return BuildParallelDofNumberingImpl(comm, topology);
}

ParallelDofNumbering BuildParallelDofNumbering(MPI_Comm comm,
                                               const TriangleDofTopology &topology)
{
  return BuildParallelDofNumberingImpl(comm, topology);
}

mfem::Array<int> GetEssentialH1TrueDofs(MPI_Comm comm, const FeatureTopology &features,
                                        const DofTopology &topology,
                                        const ParallelDofNumbering &parallel_numbering)
{
  mfem::Array<int> attributes;
  for (const auto &face : features.sheet_faces)
  {
    attributes.Append(face.boundary_attribute);
  }
  attributes.Sort();
  attributes.Unique();
  return GetEssentialH1TrueDofs(comm, features, topology, parallel_numbering, attributes);
}

mfem::Array<int>
GetEssentialH1TrueDofs(MPI_Comm comm, const FeatureTopology &features,
                       const DofTopology &topology,
                       const ParallelDofNumbering &parallel_numbering,
                       const mfem::Array<int> &essential_boundary_attributes)
{
  const auto &numbering = parallel_numbering.h1;
  if (numbering.owner.size() != topology.h1_dofs.size() ||
      numbering.local_to_true.size() != topology.h1_dofs.size() ||
      numbering.owned_offset < 0 || numbering.owned_size < 0 ||
      numbering.owned_size > std::numeric_limits<int>::max())
  {
    throw std::invalid_argument(
        "Essential singular H1 classification received inconsistent DOF numbering!");
  }

  const std::set<int> attributes(essential_boundary_attributes.begin(),
                                 essential_boundary_attributes.end());
  const auto sheet_faces = GetSelectedSheetFaces(features, &attributes);

  mfem::Array<int> result;
  for (std::size_t local = 0; local < topology.h1_dofs.size(); local++)
  {
    if (!HasSheetTrace(topology.h1_dofs[local], sheet_faces) ||
        numbering.owner[local] != Mpi::Rank(comm))
    {
      continue;
    }
    const HYPRE_BigInt true_dof = numbering.local_to_true[local] - numbering.owned_offset;
    if (true_dof < 0 || true_dof >= numbering.owned_size)
    {
      throw std::invalid_argument(
          "An essential singular H1 DOF is outside its owner partition!");
    }
    result.Append(static_cast<int>(true_dof));
  }
  result.Sort();
  if (std::adjacent_find(result.begin(), result.end()) != result.end())
  {
    throw std::logic_error("Essential singular H1 true DOFs are not unique!");
  }
  return result;
}

mfem::Array<int> GetEssentialNDTrueDofs(MPI_Comm comm, const FeatureTopology &features,
                                        const DofTopology &topology,
                                        const ParallelDofNumbering &parallel_numbering)
{
  mfem::Array<int> attributes;
  for (const auto &face : features.sheet_faces)
  {
    attributes.Append(face.boundary_attribute);
  }
  attributes.Sort();
  attributes.Unique();
  return GetEssentialNDTrueDofs(comm, features, topology, parallel_numbering, attributes);
}

mfem::Array<int>
GetEssentialNDTrueDofs(MPI_Comm comm, const FeatureTopology &features,
                       const DofTopology &topology,
                       const ParallelDofNumbering &parallel_numbering,
                       const mfem::Array<int> &essential_boundary_attributes)
{
  const auto &numbering = parallel_numbering.nd;
  if (numbering.owner.size() != topology.nd_dofs.size() ||
      numbering.local_to_true.size() != topology.nd_dofs.size() ||
      numbering.owned_offset < 0 || numbering.owned_size < 0 ||
      numbering.owned_size > std::numeric_limits<int>::max())
  {
    throw std::invalid_argument(
        "Essential singular ND classification received inconsistent DOF numbering!");
  }
  const std::set<int> attributes(essential_boundary_attributes.begin(),
                                 essential_boundary_attributes.end());
  const auto sheet_faces = GetSelectedSheetFaces(features, &attributes);

  mfem::Array<int> result;
  for (std::size_t local = 0; local < topology.nd_dofs.size(); local++)
  {
    if (!HasSheetTrace(topology.nd_dofs[local], sheet_faces) ||
        numbering.owner[local] != Mpi::Rank(comm))
    {
      continue;
    }
    const HYPRE_BigInt true_dof = numbering.local_to_true[local] - numbering.owned_offset;
    if (true_dof < 0 || true_dof >= numbering.owned_size)
    {
      throw std::invalid_argument(
          "An essential singular ND DOF is outside its owner partition!");
    }
    result.Append(static_cast<int>(true_dof));
  }
  result.Sort();
  if (std::adjacent_find(result.begin(), result.end()) != result.end())
  {
    throw std::logic_error("Essential singular ND true DOFs are not unique!");
  }
  return result;
}

mfem::Array<int>
GetEssentialTriangleH1TrueDofs(MPI_Comm comm, const TriangleFeatureTopology &features,
                               const TriangleDofTopology &topology,
                               const ParallelDofNumbering &parallel_numbering)
{
  const auto selected_segments = GetSelectedTriangleSegments(features);
  return GetEssentialH1TrueDofsImpl(comm, selected_segments, topology, parallel_numbering);
}

mfem::Array<int>
GetEssentialTriangleH1TrueDofs(MPI_Comm comm, const TriangleFeatureTopology &features,
                               const TriangleDofTopology &topology,
                               const ParallelDofNumbering &parallel_numbering,
                               const mfem::Array<int> &essential_boundary_attributes)
{
  const std::set<int> attributes(essential_boundary_attributes.begin(),
                                 essential_boundary_attributes.end());
  const auto selected_segments = GetSelectedTriangleSegments(features, &attributes);
  return GetEssentialH1TrueDofsImpl(comm, selected_segments, topology, parallel_numbering);
}

mfem::Array<int>
GetEssentialTriangleNDTrueDofs(MPI_Comm comm, const TriangleFeatureTopology &features,
                               const TriangleDofTopology &topology,
                               const ParallelDofNumbering &parallel_numbering)
{
  const auto selected_segments = GetSelectedTriangleSegments(features);
  return GetEssentialNDTrueDofsImpl(comm, selected_segments, topology, parallel_numbering);
}

mfem::Array<int>
GetEssentialTriangleNDTrueDofs(MPI_Comm comm, const TriangleFeatureTopology &features,
                               const TriangleDofTopology &topology,
                               const ParallelDofNumbering &parallel_numbering,
                               const mfem::Array<int> &essential_boundary_attributes)
{
  const std::set<int> attributes(essential_boundary_attributes.begin(),
                                 essential_boundary_attributes.end());
  const auto selected_segments = GetSelectedTriangleSegments(features, &attributes);
  return GetEssentialNDTrueDofsImpl(comm, selected_segments, topology, parallel_numbering);
}

}  // namespace singular

}  // namespace fem

}  // namespace palace
