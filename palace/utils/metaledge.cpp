// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "metaledge.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <set>
#include <tuple>
#include <utility>
#include "fem/coefficient.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/diagnostic.hpp"
#include "utils/edgedistance.hpp"
#include "utils/geodata.hpp"

namespace palace
{

namespace
{

bool IsCoincident(const mesh::BoundaryEdgeSegment &segment, const EdgeDistanceTree &tree,
                  double tolerance_squared)
{
  mfem::Vector point(3);
  for (const double t : {0.0, 0.25, 0.5, 0.75, 1.0})
  {
    for (int d = 0; d < 3; d++)
    {
      point[d] = (1.0 - t) * segment.p0[d] + t * segment.p1[d];
    }
    if (tree.DistanceSquared(point) > tolerance_squared)
    {
      return false;
    }
  }
  return true;
}

std::shared_ptr<const EdgeDistanceTree> BuildSupportTree(const mfem::ParMesh &mesh,
                                                         const std::vector<int> &attributes,
                                                         bool perimeter,
                                                         bool exterior_only = false)
{
  auto marker = mesh::BdrAttrToMarker(mesh, attributes, true);
  auto segments = perimeter
                      ? mesh::GetBoundaryEdgeSegments(mesh, marker)
                      : mesh::GetBoundaryElementEdgeSegments(mesh, marker, exterior_only);
  return segments.empty() ? nullptr
                          : std::make_shared<EdgeDistanceTree>(std::move(segments));
}

void SortAndUnique(std::vector<MetalBoundaryCondition> &conditions)
{
  std::sort(conditions.begin(), conditions.end(), [](const auto &a, const auto &b)
            { return std::tie(a.type, a.index) < std::tie(b.type, b.index); });
  conditions.erase(std::unique(conditions.begin(), conditions.end(),
                               [](const auto &a, const auto &b)
                               { return a.type == b.type && a.index == b.index; }),
                   conditions.end());
}

}  // namespace

MetalEdgeGeometry ExtractMetalEdgeGeometry(const mfem::ParMesh &mesh,
                                           const config::BoundaryData &boundaries,
                                           MetalSurfaceExtraction surface)
{
  MFEM_VERIFY(mesh.Dimension() == 3 && mesh.SpaceDimension() == 3,
              "Automatic metal edge extraction requires a three-dimensional mesh!");

  std::map<int, std::vector<MetalBoundaryCondition>> attribute_conditions;
  auto AddCondition =
      [&](const std::vector<int> &attributes, MetalBoundaryConditionType type, int index)
  {
    for (const int attribute : attributes)
    {
      MFEM_VERIFY(attribute > 0, "Metal boundary attributes must be positive!");
      attribute_conditions[attribute].push_back({type, index});
    }
  };

  AddCondition(boundaries.pec.attributes, MetalBoundaryConditionType::PEC, 0);
  AddCondition(boundaries.auxpec.attributes, MetalBoundaryConditionType::PEC, 0);
  for (const auto &[index, terminal] : boundaries.terminal)
  {
    (void)index;
    AddCondition(terminal.attributes, MetalBoundaryConditionType::PEC, 0);
  }
  for (const auto &[index, potential] : boundaries.prescribed_potential)
  {
    (void)index;
    AddCondition(potential.attributes, MetalBoundaryConditionType::PEC, 0);
    AddCondition(potential.terminal_attributes, MetalBoundaryConditionType::PEC, 0);
  }
  for (std::size_t i = 0; i < boundaries.conductivity.size(); i++)
  {
    AddCondition(boundaries.conductivity[i].attributes,
                 MetalBoundaryConditionType::CONDUCTIVITY, static_cast<int>(i));
  }
  for (std::size_t i = 0; i < boundaries.impedance.size(); i++)
  {
    AddCondition(boundaries.impedance[i].attributes, MetalBoundaryConditionType::IMPEDANCE,
                 static_cast<int>(i));
  }
  for (std::size_t i = 0; i < boundaries.rational_impedance.size(); i++)
  {
    AddCondition(boundaries.rational_impedance[i].attributes,
                 MetalBoundaryConditionType::RATIONAL_IMPEDANCE, static_cast<int>(i));
  }

  MetalEdgeGeometry result;
  if (attribute_conditions.empty())
  {
    return result;
  }
  for (auto &[attribute, conditions] : attribute_conditions)
  {
    (void)attribute;
    SortAndUnique(conditions);
    MFEM_VERIFY(conditions.size() == 1,
                "A metal boundary attribute is assigned multiple metal boundary "
                "conditions!");
  }

  std::vector<int> metal_attributes;
  metal_attributes.reserve(attribute_conditions.size());
  for (const auto &[attribute, conditions] : attribute_conditions)
  {
    (void)conditions;
    metal_attributes.push_back(attribute);
  }
  auto metal_marker = mesh::BdrAttrToMarker(mesh, metal_attributes, true);
  const auto perimeter_segments = mesh::GetBoundaryEdgeSegments(mesh, metal_marker);
  if (perimeter_segments.empty())
  {
    return result;
  }

  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
  double extent = 0.0;
  for (int d = 0; d < 3; d++)
  {
    extent = std::max(extent, bbmax[d] - bbmin[d]);
  }
  MFEM_VERIFY(extent > 0.0, "Degenerate mesh geometry for metal edge extraction!");
  const double tolerance_squared = 1.0e-18 * extent * extent;

  std::map<int, std::shared_ptr<const EdgeDistanceTree>> attribute_support;
  for (const auto &[attribute, conditions] : attribute_conditions)
  {
    (void)conditions;
    attribute_support.emplace(attribute,
                              BuildSupportTree(mesh, std::vector<int>{attribute}, false));
  }

  struct InterfaceSupport
  {
    int index;
    InterfaceDielectric type;
    std::shared_ptr<const EdgeDistanceTree> tree;
  };
  std::vector<InterfaceSupport> interface_support;
  std::set<int> interface_attributes;
  std::map<std::vector<int>, std::shared_ptr<const EdgeDistanceTree>>
      interface_support_trees;
  for (const auto &[index, dielectric] : boundaries.postpro.dielectric)
  {
    interface_attributes.insert(dielectric.attributes.begin(), dielectric.attributes.end());
    if (dielectric.type == InterfaceDielectric::DEFAULT)
    {
      continue;
    }
    auto attributes = dielectric.attributes;
    std::sort(attributes.begin(), attributes.end());
    attributes.erase(std::unique(attributes.begin(), attributes.end()), attributes.end());
    auto [tree_it, inserted] = interface_support_trees.try_emplace(attributes);
    if (inserted)
    {
      tree_it->second = BuildSupportTree(mesh, attributes, true);
    }
    auto tree = tree_it->second;
    if (tree)
    {
      interface_support.push_back({index, dielectric.type, std::move(tree)});
    }
  }

  // An exterior, nonmetal boundary surface which is not itself a dielectric interface is
  // a simulation cut surface. Internal lumped ports and sources are deliberately omitted.
  // Where an exterior face edge coincides with the metal perimeter, the metal has been
  // cut by the simulation domain (for example at a wave port or an outer box).
  std::map<int, std::shared_ptr<const EdgeDistanceTree>> truncation_support;
  const int maximum_boundary_attribute = mesh::GetMaxBdrAttribute(mesh);
  std::vector<int> boundary_attribute_present(maximum_boundary_attribute, 0);
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    boundary_attribute_present[mesh.GetBdrAttribute(be) - 1] = 1;
  }
  Mpi::GlobalMax(static_cast<int>(boundary_attribute_present.size()),
                 boundary_attribute_present.data(), mesh.GetComm());
  for (int attribute = 1; attribute <= maximum_boundary_attribute; attribute++)
  {
    if (!boundary_attribute_present[attribute - 1] ||
        attribute_conditions.find(attribute) != attribute_conditions.end() ||
        interface_attributes.find(attribute) != interface_attributes.end())
    {
      continue;
    }
    auto tree = BuildSupportTree(mesh, std::vector<int>{attribute}, false, true);
    if (tree)
    {
      truncation_support.emplace(attribute, std::move(tree));
    }
  }

  using Point = std::array<double, 3>;
  using PointKey = std::array<long long int, 3>;
  const double coordinate_tolerance = std::sqrt(tolerance_squared);
  auto GetPointKey = [&](const Point &point)
  {
    PointKey key;
    for (int d = 0; d < 3; d++)
    {
      key[d] = std::llround((point[d] - bbmin[d]) / coordinate_tolerance);
    }
    return key;
  };
  std::map<Point, std::size_t> vertex_indices;
  auto GetVertex = [&](const Point &point)
  {
    auto [it, inserted] =
        vertex_indices.try_emplace(point, static_cast<std::size_t>(result.vertices.size()));
    if (inserted)
    {
      result.vertices.push_back({point, {}, MetalEdgeVertexType::REGULAR});
    }
    return it->second;
  };

  result.segments.reserve(perimeter_segments.size());
  for (const auto &perimeter : perimeter_segments)
  {
    MetalEdgeSegment segment;
    segment.vertices = {GetVertex(perimeter.p0), GetVertex(perimeter.p1)};

    for (const auto &[attribute, tree] : attribute_support)
    {
      if (tree && IsCoincident(perimeter, *tree, tolerance_squared))
      {
        segment.metal_attributes.push_back(attribute);
        const auto &conditions = attribute_conditions.at(attribute);
        segment.conditions.insert(segment.conditions.end(), conditions.begin(),
                                  conditions.end());
      }
    }
    SortAndUnique(segment.conditions);
    MFEM_VERIFY(!segment.conditions.empty(),
                "Unable to classify an automatically extracted metal edge segment!");

    for (const auto &support : interface_support)
    {
      if (!IsCoincident(perimeter, *support.tree, tolerance_squared))
      {
        continue;
      }
      switch (support.type)
      {
        case InterfaceDielectric::SA:
          segment.sa_interfaces.push_back(support.index);
          break;
        case InterfaceDielectric::MS:
          segment.ms_interfaces.push_back(support.index);
          break;
        case InterfaceDielectric::MA:
          segment.ma_interfaces.push_back(support.index);
          break;
        case InterfaceDielectric::DEFAULT:
          break;
      }
    }

    for (const auto &[attribute, tree] : truncation_support)
    {
      if (IsCoincident(perimeter, *tree, tolerance_squared))
      {
        segment.truncation_attributes.push_back(attribute);
      }
    }
    if (!segment.truncation_attributes.empty())
    {
      segment.type = MetalEdgeSegmentType::TRUNCATION;
    }

    const std::size_t index = result.segments.size();
    result.vertices[segment.vertices[0]].segments.push_back(index);
    result.vertices[segment.vertices[1]].segments.push_back(index);
    result.segments.push_back(std::move(segment));
  }

  if (surface.classify_components || surface.retain_faces)
  {
    std::vector<MetalSurfaceFace> local_faces;
    mesh::MeshEdgeSegmentCache edge_segment_cache(mesh);
    mfem::Array<int> vertices, edges, orientations;
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      const int attribute = mesh.GetBdrAttribute(be);
      if (attribute <= 0 || attribute > metal_marker.Size() || !metal_marker[attribute - 1])
      {
        continue;
      }
      mesh.GetBdrElementVertices(be, vertices);
      MFEM_VERIFY(vertices.Size() >= 3 && vertices.Size() <= 4,
                  "Automatic metal-surface extraction supports triangular and "
                  "quadrilateral boundary elements!");
      MetalSurfaceFace face;
      std::vector<Point> corner_vertices(vertices.Size());
      for (int i = 0; i < vertices.Size(); i++)
      {
        const double *point = mesh.GetVertex(vertices[i]);
        std::copy_n(point, 3, corner_vertices[i].begin());
      }

      if (surface.retain_faces)
      {
        mesh.GetBdrElementEdges(be, edges, orientations);
        MFEM_VERIFY(edges.Size() == vertices.Size() && orientations.Size() == edges.Size(),
                    "Unexpected metal boundary face topology!");
        for (int i = 0; i < edges.Size(); i++)
        {
          auto edge_segments = edge_segment_cache.Get(edges[i]);
          if (orientations[i] < 0)
          {
            std::reverse(edge_segments.begin(), edge_segments.end());
            for (auto &segment : edge_segments)
            {
              std::swap(segment.p0, segment.p1);
            }
          }
          for (const auto &segment : edge_segments)
          {
            if (face.vertices.empty())
            {
              face.vertices.push_back(segment.p0);
            }
            else
            {
              double distance_squared = 0.0;
              for (int d = 0; d < 3; d++)
              {
                const double delta = face.vertices.back()[d] - segment.p0[d];
                distance_squared += delta * delta;
              }
              MFEM_VERIFY(distance_squared <= tolerance_squared,
                          "Metal boundary face edges do not form an ordered loop!");
            }
            face.vertices.push_back(segment.p1);
          }
        }
        double closure_distance_squared = 0.0;
        for (int d = 0; d < 3; d++)
        {
          const double delta = face.vertices.front()[d] - face.vertices.back()[d];
          closure_distance_squared += delta * delta;
        }
        MFEM_VERIFY(closure_distance_squared <= tolerance_squared,
                    "Metal boundary face edges do not form a closed loop!");
        face.vertices.pop_back();

        // Preserve the original representation for a geometrically linear face.
        if (face.vertices.size() == corner_vertices.size())
        {
          face.vertices = std::move(corner_vertices);
        }
      }
      else
      {
        face.vertices = std::move(corner_vertices);
      }
      local_faces.push_back(std::move(face));
    }

    std::vector<std::size_t> local_parent(local_faces.size());
    std::iota(local_parent.begin(), local_parent.end(), 0);
    auto FindLocal = [&](std::size_t item)
    {
      std::size_t root = item;
      while (local_parent[root] != root)
      {
        root = local_parent[root];
      }
      while (local_parent[item] != item)
      {
        const std::size_t next = local_parent[item];
        local_parent[item] = root;
        item = next;
      }
      return root;
    };
    auto UnionLocal = [&](std::size_t first, std::size_t second)
    {
      first = FindLocal(first);
      second = FindLocal(second);
      if (first != second)
      {
        local_parent[std::max(first, second)] = std::min(first, second);
      }
    };
    std::map<PointKey, std::size_t> face_by_vertex;
    for (std::size_t face = 0; face < local_faces.size(); face++)
    {
      for (const auto &point : local_faces[face].vertices)
      {
        auto [vertex, inserted] = face_by_vertex.emplace(GetPointKey(point), face);
        if (!inserted)
        {
          UnionLocal(face, vertex->second);
        }
      }
    }

    std::map<std::size_t, int> local_component_by_root;
    std::vector<int> local_face_components(local_faces.size());
    for (std::size_t face = 0; face < local_faces.size(); face++)
    {
      auto [component, inserted] =
          local_component_by_root.emplace(FindLocal(face), local_component_by_root.size());
      (void)inserted;
      local_face_components[face] = component->second;
    }

    const int local_component_count = static_cast<int>(local_component_by_root.size());
    std::vector<int> component_counts(Mpi::Size(mesh.GetComm()));
    Mpi::Allgather(1, &local_component_count, component_counts.data(), mesh.GetComm());
    std::vector<long long int> component_offsets(component_counts.size() + 1);
    for (std::size_t rank = 0; rank < component_counts.size(); rank++)
    {
      component_offsets[rank + 1] = component_offsets[rank] + component_counts[rank];
    }
    const long long int component_offset = component_offsets[Mpi::Rank(mesh.GetComm())];

    std::map<PointKey, long long int> local_component_by_vertex;
    for (std::size_t face = 0; face < local_faces.size(); face++)
    {
      const long long int component = component_offset + local_face_components[face];
      for (const auto &point : local_faces[face].vertices)
      {
        auto [vertex, inserted] =
            local_component_by_vertex.emplace(GetPointKey(point), component);
        MFEM_VERIFY(inserted || vertex->second == component,
                    "A metal-surface vertex belongs to multiple local components!");
      }
    }

    const int local_vertex_count = static_cast<int>(local_component_by_vertex.size());
    std::vector<int> vertex_counts(component_counts.size());
    Mpi::Allgather(1, &local_vertex_count, vertex_counts.data(), mesh.GetComm());
    std::vector<int> vertex_offsets(vertex_counts.size()), key_counts(vertex_counts.size()),
        key_offsets(vertex_counts.size());
    int total_vertices = 0;
    for (std::size_t rank = 0; rank < vertex_counts.size(); rank++)
    {
      vertex_offsets[rank] = total_vertices;
      key_counts[rank] = 3 * vertex_counts[rank];
      key_offsets[rank] = 3 * total_vertices;
      total_vertices += vertex_counts[rank];
    }
    std::vector<long long int> local_keys;
    std::vector<long long int> local_components;
    local_keys.reserve(3 * local_vertex_count);
    local_components.reserve(local_vertex_count);
    for (const auto &[key, component] : local_component_by_vertex)
    {
      local_keys.insert(local_keys.end(), key.begin(), key.end());
      local_components.push_back(component);
    }
    std::vector<long long int> global_keys(3 * total_vertices);
    std::vector<long long int> global_components(total_vertices);
    Mpi::Allgatherv(static_cast<int>(local_keys.size()), local_keys.data(),
                    global_keys.data(), key_counts.data(), key_offsets.data(),
                    mesh.GetComm());
    Mpi::Allgatherv(local_vertex_count, local_components.data(), global_components.data(),
                    vertex_counts.data(), vertex_offsets.data(), mesh.GetComm());

    const long long int total_components = component_offsets.back();
    std::vector<long long int> parent(total_components);
    std::iota(parent.begin(), parent.end(), 0);
    auto Find = [&](long long int item)
    {
      long long int root = item;
      while (parent[root] != root)
      {
        root = parent[root];
      }
      while (parent[item] != item)
      {
        const long long int next = parent[item];
        parent[item] = root;
        item = next;
      }
      return root;
    };
    auto Union = [&](long long int first, long long int second)
    {
      first = Find(first);
      second = Find(second);
      if (first != second)
      {
        parent[std::max(first, second)] = std::min(first, second);
      }
    };

    std::map<PointKey, long long int> provisional_component_by_vertex;
    for (int vertex = 0; vertex < total_vertices; vertex++)
    {
      PointKey key;
      std::copy_n(global_keys.data() + 3 * vertex, 3, key.begin());
      auto [entry, inserted] =
          provisional_component_by_vertex.emplace(key, global_components[vertex]);
      if (!inserted)
      {
        Union(entry->second, global_components[vertex]);
      }
    }
    std::map<long long int, int> component_by_root;
    for (long long int component = 0; component < total_components; component++)
    {
      component_by_root.try_emplace(Find(component), component_by_root.size());
    }
    result.metal_components = static_cast<int>(component_by_root.size());

    std::map<PointKey, int> component_by_vertex;
    for (const auto &[key, provisional] : provisional_component_by_vertex)
    {
      component_by_vertex.emplace(key, component_by_root.at(Find(provisional)));
    }
    auto FindPointComponent = [&](const Point &point)
    {
      const auto component = component_by_vertex.find(GetPointKey(point));
      return component == component_by_vertex.end() ? -1 : component->second;
    };

    for (auto &segment : result.segments)
    {
      const int first = FindPointComponent(result.vertices[segment.vertices[0]].coordinate);
      const int second =
          FindPointComponent(result.vertices[segment.vertices[1]].coordinate);
      MFEM_VERIFY(first >= 0 && second >= 0 && first == second,
                  "Unable to associate a metal perimeter edge with one supporting "
                  "metal surface!");
      segment.metal_component = first;
    }

    if (surface.retain_faces)
    {
      for (auto &face : local_faces)
      {
        face.component = FindPointComponent(face.vertices.front());
        MFEM_VERIFY(
            face.component >= 0 &&
                std::all_of(face.vertices.begin(), face.vertices.end(),
                            [&](const auto &point)
                            { return FindPointComponent(point) == face.component; }),
            "A metal boundary face belongs to multiple connected metal surfaces!");
      }
      result.surface_faces = std::move(local_faces);
    }
  }

  auto LabelComponents = [&](bool physical)
  {
    std::vector<bool> visited(result.segments.size(), false);
    int components = 0;
    for (std::size_t seed = 0; seed < result.segments.size(); seed++)
    {
      if (visited[seed] ||
          (physical && result.segments[seed].type != MetalEdgeSegmentType::PHYSICAL))
      {
        continue;
      }
      std::queue<std::size_t> queue;
      queue.push(seed);
      visited[seed] = true;
      while (!queue.empty())
      {
        const std::size_t current = queue.front();
        queue.pop();
        if (physical)
        {
          result.segments[current].physical_component = components;
        }
        else
        {
          result.segments[current].component = components;
        }
        for (const std::size_t vertex : result.segments[current].vertices)
        {
          for (const std::size_t neighbor : result.vertices[vertex].segments)
          {
            if (!visited[neighbor] && (!physical || result.segments[neighbor].type ==
                                                        MetalEdgeSegmentType::PHYSICAL))
            {
              visited[neighbor] = true;
              queue.push(neighbor);
            }
          }
        }
      }
      components++;
    }
    return components;
  };
  result.components = LabelComponents(false);
  result.physical_components = LabelComponents(true);

  // Boundary meshes commonly represent smooth layout curves by short polygonal facets.
  // Treat modest local turns as part of the same smooth chain so that chain topology does
  // not depend on the curve tessellation. Sharper turns remain explicit corner vertices
  // for separate corner treatment by a surface-response model.
  constexpr double corner_angle_tolerance_degrees = 30.0;
  const double straight_dot_tolerance =
      -std::cos(corner_angle_tolerance_degrees * std::acos(-1.0) / 180.0);
  auto ClassifyVertex = [&](std::size_t vertex_index,
                            bool physical) -> std::optional<MetalEdgeVertexType>
  {
    auto &vertex = result.vertices[vertex_index];
    std::vector<std::size_t> segments;
    segments.reserve(vertex.segments.size());
    for (const std::size_t segment : vertex.segments)
    {
      if (!physical || result.segments[segment].type == MetalEdgeSegmentType::PHYSICAL)
      {
        segments.push_back(segment);
      }
    }
    if (segments.empty())
    {
      return std::nullopt;
    }
    if (segments.size() == 1)
    {
      return MetalEdgeVertexType::ENDPOINT;
    }
    if (segments.size() > 2)
    {
      return MetalEdgeVertexType::JUNCTION;
    }

    std::array<std::array<double, 3>, 2> directions{};
    for (int i = 0; i < 2; i++)
    {
      const auto &edge = result.segments[segments[i]];
      const std::size_t other =
          edge.vertices[0] == vertex_index ? edge.vertices[1] : edge.vertices[0];
      double norm_squared = 0.0;
      for (int d = 0; d < 3; d++)
      {
        directions[i][d] = result.vertices[other].coordinate[d] - vertex.coordinate[d];
        norm_squared += directions[i][d] * directions[i][d];
      }
      MFEM_VERIFY(norm_squared > 0.0, "Metal edge graph contains a zero-length segment!");
      const double inverse_norm = 1.0 / std::sqrt(norm_squared);
      for (double &value : directions[i])
      {
        value *= inverse_norm;
      }
    }
    double dot = 0.0;
    for (int d = 0; d < 3; d++)
    {
      dot += directions[0][d] * directions[1][d];
    }
    return dot <= straight_dot_tolerance ? MetalEdgeVertexType::REGULAR
                                         : MetalEdgeVertexType::CORNER;
  };
  for (std::size_t vertex_index = 0; vertex_index < result.vertices.size(); vertex_index++)
  {
    auto &vertex = result.vertices[vertex_index];
    vertex.type = *ClassifyVertex(vertex_index, false);
    vertex.physical_type = ClassifyVertex(vertex_index, true);
    vertex.on_truncation_boundary = std::any_of(
        vertex.segments.begin(), vertex.segments.end(), [&](std::size_t segment)
        { return result.segments[segment].type == MetalEdgeSegmentType::TRUNCATION; });
  }

  // A physical chain is a maximal path which can pass through regular (locally straight)
  // vertices but stops at corners, endpoints, and junctions. This grouping is independent
  // of the finite-element subdivision along a straight fabricated edge.
  std::vector<bool> chain_visited(result.segments.size(), false);
  for (std::size_t seed = 0; seed < result.segments.size(); seed++)
  {
    if (chain_visited[seed] || result.segments[seed].type != MetalEdgeSegmentType::PHYSICAL)
    {
      continue;
    }
    std::queue<std::size_t> queue;
    queue.push(seed);
    chain_visited[seed] = true;
    while (!queue.empty())
    {
      const std::size_t current = queue.front();
      queue.pop();
      result.segments[current].physical_chain = result.physical_chains;
      for (const std::size_t vertex_index : result.segments[current].vertices)
      {
        const auto &vertex = result.vertices[vertex_index];
        if (vertex.physical_type != MetalEdgeVertexType::REGULAR)
        {
          continue;
        }
        for (const std::size_t neighbor : vertex.segments)
        {
          if (!chain_visited[neighbor] &&
              result.segments[neighbor].type == MetalEdgeSegmentType::PHYSICAL)
          {
            chain_visited[neighbor] = true;
            queue.push(neighbor);
          }
        }
      }
    }
    result.physical_chains++;
  }
  return result;
}

std::vector<std::array<double, 3>>
BuildMetalEdgeProcessNormals(const mfem::ParMesh &mesh, const MetalEdgeGeometry &geometry,
                             const std::vector<std::size_t> &segment_indices,
                             const std::function<double(int)> &material_score,
                             const std::optional<std::array<double, 3>> &fallback)
{
  MFEM_VERIFY(mesh.Dimension() == 3 && mesh.SpaceDimension() == 3,
              "Automatic metal edge frames require a three-dimensional mesh!");
  MFEM_VERIFY(!segment_indices.empty(),
              "Cannot infer process normals for an empty metal edge selection!");

  using Point = std::array<double, 3>;
  using PointKey = std::array<std::int64_t, 3>;
  using SegmentKey = std::pair<PointKey, PointKey>;

  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
  double extent = 0.0;
  for (int d = 0; d < 3; d++)
  {
    extent = std::max(extent, bbmax[d] - bbmin[d]);
  }
  MFEM_VERIFY(extent > 0.0, "Degenerate geometry for automatic metal edge frames!");
  const double coordinate_tolerance = 1.0e-10 * extent;
  auto GetPointKey = [&](const Point &point)
  {
    PointKey key;
    for (int d = 0; d < 3; d++)
    {
      key[d] = std::llround((point[d] - bbmin[d]) / coordinate_tolerance);
    }
    return key;
  };
  auto GetSegmentKey = [&](Point p0, Point p1)
  {
    PointKey k0 = GetPointKey(p0);
    PointKey k1 = GetPointKey(p1);
    if (k1 < k0)
    {
      std::swap(k0, k1);
    }
    return SegmentKey{k0, k1};
  };

  std::map<SegmentKey, std::size_t> selected_segments;
  std::set<int> selected_attributes;
  for (std::size_t i = 0; i < segment_indices.size(); i++)
  {
    const std::size_t segment_index = segment_indices[i];
    MFEM_VERIFY(segment_index < geometry.segments.size(),
                "Invalid metal edge segment index for process-normal inference!");
    const auto &segment = geometry.segments[segment_index];
    MFEM_VERIFY(segment.type == MetalEdgeSegmentType::PHYSICAL,
                "Cannot infer a process normal for a truncation segment!");
    const auto key = GetSegmentKey(geometry.vertices[segment.vertices[0]].coordinate,
                                   geometry.vertices[segment.vertices[1]].coordinate);
    MFEM_VERIFY(selected_segments.emplace(key, i).second,
                "Duplicate physical metal edge geometry!");
    selected_attributes.insert(segment.metal_attributes.begin(),
                               segment.metal_attributes.end());
  }

  struct Candidate
  {
    std::size_t segment;
    double score;
    std::array<double, 3> normal;
  };
  std::vector<Candidate> candidates;
  mesh::MeshEdgeSegmentCache edge_segment_cache(mesh);
  auto &mutable_mesh = const_cast<mfem::ParMesh &>(mesh);
  mutable_mesh.ExchangeFaceNbrData();
  mfem::Array<int> edges, orientations;
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  mfem::Vector normal(3);
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int attribute = mesh.GetBdrAttribute(be);
    if (selected_attributes.find(attribute) == selected_attributes.end())
    {
      continue;
    }

    auto *T = mutable_mesh.GetBdrElementTransformation(be);
    const auto &ip = mfem::Geometries.GetCenter(T->GetGeometryType());
    T->SetIntPoint(&ip);
    const bool orientation =
        BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(be, mesh, FET, T1,
                                                                         T2, &ip);
    BdrGridFunctionCoefficient::GetNormal(*T, normal, orientation);

    std::array<double, 3> normal_into_element_1{};
    std::copy_n(normal.GetData(), 3, normal_into_element_1.begin());
    mesh.GetBdrElementEdges(be, edges, orientations);
    for (const int edge : edges)
    {
      for (const auto &edge_segment : edge_segment_cache.Get(edge))
      {
        const auto selected =
            selected_segments.find(GetSegmentKey(edge_segment.p0, edge_segment.p1));
        if (selected == selected_segments.end())
        {
          continue;
        }
        const auto &segment = geometry.segments[segment_indices[selected->second]];
        if (std::find(segment.metal_attributes.begin(), segment.metal_attributes.end(),
                      attribute) == segment.metal_attributes.end())
        {
          continue;
        }

        candidates.push_back({selected->second, material_score(FET.Elem1->Attribute),
                              normal_into_element_1});
        if (FET.Elem2)
        {
          auto normal_into_element_2 = normal_into_element_1;
          for (double &value : normal_into_element_2)
          {
            value *= -1.0;
          }
          candidates.push_back({selected->second, material_score(FET.Elem2->Attribute),
                                normal_into_element_2});
        }
      }
    }
  }

  std::vector<double> minimum_score(segment_indices.size(), mfem::infinity());
  std::vector<double> maximum_score(segment_indices.size(), -mfem::infinity());
  for (const auto &candidate : candidates)
  {
    minimum_score[candidate.segment] =
        std::min(minimum_score[candidate.segment], candidate.score);
    maximum_score[candidate.segment] =
        std::max(maximum_score[candidate.segment], candidate.score);
  }
  Mpi::GlobalMin(static_cast<int>(minimum_score.size()), minimum_score.data(),
                 mesh.GetComm());
  Mpi::GlobalMax(static_cast<int>(maximum_score.size()), maximum_score.data(),
                 mesh.GetComm());

  std::vector<double> normal_sum(3 * segment_indices.size(), 0.0);
  for (const auto &candidate : candidates)
  {
    const double minimum = minimum_score[candidate.segment];
    const double maximum = maximum_score[candidate.segment];
    const double tolerance =
        1.0e-10 * std::max({1.0, std::abs(minimum), std::abs(maximum)});
    double sign = 1.0;
    if (maximum - minimum > tolerance)
    {
      sign = candidate.score >= 0.5 * (minimum + maximum) ? 1.0 : -1.0;
    }
    else if (fallback)
    {
      double dot = 0.0;
      for (int d = 0; d < 3; d++)
      {
        dot += candidate.normal[d] * (*fallback)[d];
      }
      sign = dot >= 0.0 ? 1.0 : -1.0;
    }
    else
    {
      int dominant = 0;
      for (int d = 1; d < 3; d++)
      {
        if (std::abs(candidate.normal[d]) > std::abs(candidate.normal[dominant]))
        {
          dominant = d;
        }
      }
      sign = candidate.normal[dominant] >= 0.0 ? 1.0 : -1.0;
    }
    for (int d = 0; d < 3; d++)
    {
      normal_sum[3 * candidate.segment + d] += sign * candidate.normal[d];
    }
  }
  Mpi::GlobalSum(static_cast<int>(normal_sum.size()), normal_sum.data(), mesh.GetComm());

  std::vector<std::array<double, 3>> process_normals(segment_indices.size());
  for (std::size_t i = 0; i < segment_indices.size(); i++)
  {
    MFEM_VERIFY(std::isfinite(minimum_score[i]) && std::isfinite(maximum_score[i]),
                "No supporting metal face was found for an automatic edge segment!");
    auto &process_normal = process_normals[i];
    std::copy_n(normal_sum.data() + 3 * i, 3, process_normal.begin());

    const auto &segment = geometry.segments[segment_indices[i]];
    const auto &p0 = geometry.vertices[segment.vertices[0]].coordinate;
    const auto &p1 = geometry.vertices[segment.vertices[1]].coordinate;
    std::array<double, 3> tangent{};
    double tangent_norm_squared = 0.0;
    for (int d = 0; d < 3; d++)
    {
      tangent[d] = p1[d] - p0[d];
      tangent_norm_squared += tangent[d] * tangent[d];
    }
    MFEM_VERIFY(tangent_norm_squared > 0.0,
                "Cannot infer a process normal for a zero-length edge segment!");
    double normal_tangent = 0.0;
    for (int d = 0; d < 3; d++)
    {
      normal_tangent += process_normal[d] * tangent[d];
    }
    for (int d = 0; d < 3; d++)
    {
      process_normal[d] -= normal_tangent * tangent[d] / tangent_norm_squared;
    }

    double norm_squared = 0.0;
    for (double value : process_normal)
    {
      norm_squared += value * value;
    }
    if (norm_squared <= 1.0e-20 && fallback)
    {
      process_normal = *fallback;
      normal_tangent = 0.0;
      for (int d = 0; d < 3; d++)
      {
        normal_tangent += process_normal[d] * tangent[d];
      }
      for (int d = 0; d < 3; d++)
      {
        process_normal[d] -= normal_tangent * tangent[d] / tangent_norm_squared;
      }
      norm_squared = 0.0;
      for (double value : process_normal)
      {
        norm_squared += value * value;
      }
    }
    MFEM_VERIFY(norm_squared > 1.0e-20,
                "Unable to infer a process normal for an automatic edge segment!");
    const double inverse_norm = 1.0 / std::sqrt(norm_squared);
    for (double &value : process_normal)
    {
      value *= inverse_norm;
    }
  }
  return process_normals;
}

std::vector<std::array<double, 3>>
BuildMetalEdgeGapDirections(const mfem::ParMesh &mesh, const MetalEdgeGeometry &geometry,
                            const std::vector<std::size_t> &segment_indices,
                            const std::vector<std::array<double, 3>> &process_normals)
{
  MFEM_VERIFY(mesh.Dimension() == 3 && mesh.SpaceDimension() == 3,
              "Automatic metal edge frames require a three-dimensional mesh!");
  MFEM_VERIFY(!segment_indices.empty() && process_normals.size() == segment_indices.size(),
              "Automatic metal edge gap directions require matching nonempty segment "
              "and process-normal lists!");

  using Point = std::array<double, 3>;
  using PointKey = std::array<std::int64_t, 3>;
  using SegmentKey = std::pair<PointKey, PointKey>;

  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
  double extent = 0.0;
  for (int d = 0; d < 3; d++)
  {
    extent = std::max(extent, bbmax[d] - bbmin[d]);
  }
  MFEM_VERIFY(extent > 0.0, "Degenerate geometry for automatic metal edge frames!");
  const double coordinate_tolerance = 1.0e-10 * extent;
  auto GetPointKey = [&](const Point &point)
  {
    PointKey key;
    for (int d = 0; d < 3; d++)
    {
      key[d] = std::llround((point[d] - bbmin[d]) / coordinate_tolerance);
    }
    return key;
  };
  auto GetSegmentKey = [&](Point p0, Point p1)
  {
    PointKey k0 = GetPointKey(p0);
    PointKey k1 = GetPointKey(p1);
    if (k1 < k0)
    {
      std::swap(k0, k1);
    }
    return SegmentKey{k0, k1};
  };

  std::map<SegmentKey, std::size_t> selected_segments;
  std::set<int> selected_attributes;
  for (std::size_t i = 0; i < segment_indices.size(); i++)
  {
    const std::size_t segment_index = segment_indices[i];
    MFEM_VERIFY(segment_index < geometry.segments.size(),
                "Invalid metal edge segment index for gap-direction inference!");
    const auto &segment = geometry.segments[segment_index];
    MFEM_VERIFY(segment.type == MetalEdgeSegmentType::PHYSICAL,
                "Cannot infer a gap direction for a truncation segment!");
    const auto key = GetSegmentKey(geometry.vertices[segment.vertices[0]].coordinate,
                                   geometry.vertices[segment.vertices[1]].coordinate);
    MFEM_VERIFY(selected_segments.emplace(key, i).second,
                "Duplicate physical metal edge geometry!");
    selected_attributes.insert(segment.metal_attributes.begin(),
                               segment.metal_attributes.end());
  }

  std::vector<double> inward_sum(3 * segment_indices.size(), 0.0);
  mesh::MeshEdgeSegmentCache edge_segment_cache(mesh);
  auto &mutable_mesh = const_cast<mfem::ParMesh &>(mesh);
  mfem::Array<int> edges, orientations;
  mfem::Vector center(3);
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int attribute = mesh.GetBdrAttribute(be);
    if (selected_attributes.find(attribute) == selected_attributes.end())
    {
      continue;
    }
    auto *T = mutable_mesh.GetBdrElementTransformation(be);
    const auto &ip = mfem::Geometries.GetCenter(T->GetGeometryType());
    T->Transform(ip, center);
    mesh.GetBdrElementEdges(be, edges, orientations);
    for (const int edge : edges)
    {
      for (const auto &edge_segment : edge_segment_cache.Get(edge))
      {
        const Point &p0 = edge_segment.p0;
        const Point &p1 = edge_segment.p1;
        const auto selected = selected_segments.find(GetSegmentKey(p0, p1));
        if (selected == selected_segments.end())
        {
          continue;
        }
        const std::size_t local_index = selected->second;
        const auto &segment = geometry.segments[segment_indices[local_index]];
        if (std::find(segment.metal_attributes.begin(), segment.metal_attributes.end(),
                      attribute) == segment.metal_attributes.end())
        {
          continue;
        }

        const auto &normal = process_normals[local_index];
        Point tangent{}, inward{};
        double tangent_norm_squared = 0.0;
        for (int d = 0; d < 3; d++)
        {
          tangent[d] = p1[d] - p0[d];
          tangent_norm_squared += tangent[d] * tangent[d];
          inward[d] = center[d] - 0.5 * (p0[d] + p1[d]);
        }
        MFEM_VERIFY(tangent_norm_squared > 0.0,
                    "Cannot infer a gap direction for a zero-length edge segment!");
        double inward_tangent = 0.0, inward_normal = 0.0;
        for (int d = 0; d < 3; d++)
        {
          inward_tangent += inward[d] * tangent[d];
          inward_normal += inward[d] * normal[d];
        }
        double inward_norm_squared = 0.0;
        for (int d = 0; d < 3; d++)
        {
          inward[d] -= inward_tangent * tangent[d] / tangent_norm_squared +
                       inward_normal * normal[d];
          inward_norm_squared += inward[d] * inward[d];
        }
        if (inward_norm_squared <= coordinate_tolerance * coordinate_tolerance)
        {
          continue;
        }
        const double inverse_norm = 1.0 / std::sqrt(inward_norm_squared);
        for (int d = 0; d < 3; d++)
        {
          inward_sum[3 * local_index + d] += inward[d] * inverse_norm;
        }
      }
    }
  }
  Mpi::GlobalSum(static_cast<int>(inward_sum.size()), inward_sum.data(), mesh.GetComm());

  std::vector<std::array<double, 3>> gap_directions(segment_indices.size());
  for (std::size_t i = 0; i < segment_indices.size(); i++)
  {
    auto &gap = gap_directions[i];
    double norm_squared = 0.0;
    for (int d = 0; d < 3; d++)
    {
      gap[d] = -inward_sum[3 * i + d];
      norm_squared += gap[d] * gap[d];
    }
    MFEM_VERIFY(norm_squared > 1.0e-20,
                "Unable to infer the metal-to-gap direction for an automatic edge "
                "segment!");
    const double inverse_norm = 1.0 / std::sqrt(norm_squared);
    for (double &value : gap)
    {
      value *= inverse_norm;
    }
  }
  return gap_directions;
}

}  // namespace palace
