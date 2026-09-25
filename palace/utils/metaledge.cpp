// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "metaledge.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <set>
#include <string>
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

struct SegmentVectorContribution
{
  std::size_t segment;
  std::array<double, 3> vector;
};

// Sum per-segment vector contributions from every rank in a canonical order (segment
// index, then the contribution itself). The per-face contributions of a segment are a
// property of the global mesh while their distribution over ranks is not, so an
// MPI_Allreduce of rank-local partial sums would make the result depend on the partition
// through floating-point roundoff.
std::vector<double>
SumSegmentVectorsInCanonicalOrder(MPI_Comm comm, std::size_t segment_count,
                                  const std::vector<SegmentVectorContribution> &local)
{
  std::vector<double> local_records;
  local_records.reserve(4 * local.size());
  for (const auto &contribution : local)
  {
    local_records.push_back(static_cast<double>(contribution.segment));
    local_records.insert(local_records.end(), contribution.vector.begin(),
                         contribution.vector.end());
  }
  MFEM_VERIFY(local_records.size() <=
                  static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "Local metal edge frame data exceeds the MPI count limit!");
  const int local_count = static_cast<int>(local_records.size());
  std::vector<int> counts(Mpi::Size(comm));
  Mpi::Allgather(1, &local_count, counts.data(), comm);
  std::vector<int> offsets(counts.size());
  int total = 0;
  for (std::size_t rank = 0; rank < counts.size(); rank++)
  {
    offsets[rank] = total;
    MFEM_VERIFY(counts[rank] <= std::numeric_limits<int>::max() - total,
                "Global metal edge frame data exceeds the MPI count limit!");
    total += counts[rank];
  }
  std::vector<double> records(total);
  Mpi::Allgatherv(local_count, local_records.data(), records.data(), counts.data(),
                  offsets.data(), comm);

  std::vector<SegmentVectorContribution> contributions(total / 4);
  for (std::size_t i = 0; i < contributions.size(); i++)
  {
    contributions[i].segment = static_cast<std::size_t>(std::llround(records[4 * i]));
    std::copy_n(records.data() + 4 * i + 1, 3, contributions[i].vector.begin());
  }
  std::sort(contributions.begin(), contributions.end(), [](const auto &a, const auto &b)
            { return std::tie(a.segment, a.vector) < std::tie(b.segment, b.vector); });

  std::vector<double> sum(3 * segment_count, 0.0);
  for (const auto &contribution : contributions)
  {
    MFEM_VERIFY(contribution.segment < segment_count,
                "Invalid segment index in a gathered metal edge frame contribution!");
    for (int d = 0; d < 3; d++)
    {
      sum[3 * contribution.segment + d] += contribution.vector[d];
    }
  }
  return sum;
}

// Flat buffers for the broadcast of the perimeter built on the root and the scatter of
// every rank's retained facets: integers and reals, read back in the order written.
struct GeometryBuffers
{
  std::vector<long long int> ints;
  std::vector<double> reals;
};

void PackInts(GeometryBuffers &b, const std::vector<int> &values)
{
  b.ints.push_back(static_cast<long long int>(values.size()));
  b.ints.insert(b.ints.end(), values.begin(), values.end());
}

std::vector<int> UnpackInts(const GeometryBuffers &b, std::size_t &pos)
{
  const auto n = static_cast<std::size_t>(b.ints[pos++]);
  std::vector<int> values(n);
  for (std::size_t i = 0; i < n; i++)
  {
    values[i] = static_cast<int>(b.ints[pos++]);
  }
  return values;
}

void PackFace(const MetalSurfaceFace &face, GeometryBuffers &b)
{
  b.ints.push_back(face.component);
  b.ints.push_back(face.attribute);
  b.ints.push_back(face.on_bounding_box ? 1 : 0);
  b.ints.push_back(static_cast<long long int>(face.vertices.size()));
  b.reals.insert(b.reals.end(), face.normal.begin(), face.normal.end());
  for (const auto &p : face.vertices)
  {
    b.reals.insert(b.reals.end(), p.begin(), p.end());
  }
}

MetalSurfaceFace UnpackFace(const GeometryBuffers &b, std::size_t &int_pos,
                            std::size_t &real_pos)
{
  MetalSurfaceFace face;
  face.component = static_cast<int>(b.ints[int_pos++]);
  face.attribute = static_cast<int>(b.ints[int_pos++]);
  face.on_bounding_box = b.ints[int_pos++] != 0;
  const auto n = static_cast<std::size_t>(b.ints[int_pos++]);
  std::copy_n(b.reals.data() + real_pos, 3, face.normal.begin());
  real_pos += 3;
  face.vertices.resize(n);
  for (std::size_t i = 0; i < n; i++)
  {
    std::copy_n(b.reals.data() + real_pos, 3, face.vertices[i].begin());
    real_pos += 3;
  }
  return face;
}

// The perimeter without faces: scalars, vertices, segments.
void PackGeometry(const MetalEdgeGeometry &g, GeometryBuffers &b)
{
  b.reals.insert(b.reals.end(), g.layer_normal.begin(), g.layer_normal.end());
  b.ints.push_back(g.components);
  b.ints.push_back(g.physical_components);
  b.ints.push_back(g.physical_chains);
  b.ints.push_back(g.metal_components);
  b.ints.push_back(static_cast<long long int>(g.vertices.size()));
  for (const auto &v : g.vertices)
  {
    b.reals.insert(b.reals.end(), v.coordinate.begin(), v.coordinate.end());
    b.ints.push_back(static_cast<long long int>(v.segments.size()));
    for (const std::size_t s : v.segments)
    {
      b.ints.push_back(static_cast<long long int>(s));
    }
    b.ints.push_back(static_cast<long long int>(v.type));
    b.ints.push_back(v.physical_type ? static_cast<long long int>(*v.physical_type) : -1);
    b.ints.push_back(v.on_truncation_boundary ? 1 : 0);
    b.ints.push_back(v.on_port_boundary ? 1 : 0);
  }
  b.ints.push_back(static_cast<long long int>(g.segments.size()));
  for (const auto &seg : g.segments)
  {
    b.ints.push_back(static_cast<long long int>(seg.vertices[0]));
    b.ints.push_back(static_cast<long long int>(seg.vertices[1]));
    b.ints.push_back(seg.component);
    b.ints.push_back(seg.physical_component);
    b.ints.push_back(seg.physical_chain);
    b.ints.push_back(seg.metal_component);
    b.ints.push_back(static_cast<long long int>(seg.type));
    PackInts(b, seg.metal_attributes);
    b.ints.push_back(static_cast<long long int>(seg.conditions.size()));
    for (const auto &c : seg.conditions)
    {
      b.ints.push_back(static_cast<long long int>(c.type));
      b.ints.push_back(c.index);
    }
    PackInts(b, seg.sa_interfaces);
    PackInts(b, seg.ms_interfaces);
    PackInts(b, seg.ma_interfaces);
    PackInts(b, seg.truncation_attributes);
    PackInts(b, seg.port_attributes);
    b.ints.push_back(seg.face_count);
    b.ints.push_back(static_cast<long long int>(seg.face_normals.size()));
    for (const auto &n : seg.face_normals)
    {
      b.reals.insert(b.reals.end(), n.begin(), n.end());
    }
    PackInts(b, seg.side_attributes);
    b.ints.push_back(seg.on_bounding_box ? 1 : 0);
  }
}

void UnpackGeometry(const GeometryBuffers &b, MetalEdgeGeometry &g)
{
  std::size_t ip = 0, rp = 0;
  std::copy_n(b.reals.data() + rp, 3, g.layer_normal.begin());
  rp += 3;
  g.components = static_cast<int>(b.ints[ip++]);
  g.physical_components = static_cast<int>(b.ints[ip++]);
  g.physical_chains = static_cast<int>(b.ints[ip++]);
  g.metal_components = static_cast<int>(b.ints[ip++]);
  g.vertices.resize(static_cast<std::size_t>(b.ints[ip++]));
  for (auto &v : g.vertices)
  {
    std::copy_n(b.reals.data() + rp, 3, v.coordinate.begin());
    rp += 3;
    v.segments.resize(static_cast<std::size_t>(b.ints[ip++]));
    for (auto &s : v.segments)
    {
      s = static_cast<std::size_t>(b.ints[ip++]);
    }
    v.type = static_cast<MetalEdgeVertexType>(b.ints[ip++]);
    const long long int physical = b.ints[ip++];
    v.physical_type = physical < 0 ? std::nullopt
                                   : std::optional<MetalEdgeVertexType>(
                                         static_cast<MetalEdgeVertexType>(physical));
    v.on_truncation_boundary = b.ints[ip++] != 0;
    v.on_port_boundary = b.ints[ip++] != 0;
  }
  g.segments.resize(static_cast<std::size_t>(b.ints[ip++]));
  for (auto &seg : g.segments)
  {
    seg.vertices[0] = static_cast<std::size_t>(b.ints[ip++]);
    seg.vertices[1] = static_cast<std::size_t>(b.ints[ip++]);
    seg.component = static_cast<int>(b.ints[ip++]);
    seg.physical_component = static_cast<int>(b.ints[ip++]);
    seg.physical_chain = static_cast<int>(b.ints[ip++]);
    seg.metal_component = static_cast<int>(b.ints[ip++]);
    seg.type = static_cast<MetalEdgeSegmentType>(b.ints[ip++]);
    seg.metal_attributes = UnpackInts(b, ip);
    seg.conditions.resize(static_cast<std::size_t>(b.ints[ip++]));
    for (auto &c : seg.conditions)
    {
      c.type = static_cast<MetalBoundaryConditionType>(b.ints[ip++]);
      c.index = static_cast<int>(b.ints[ip++]);
    }
    seg.sa_interfaces = UnpackInts(b, ip);
    seg.ms_interfaces = UnpackInts(b, ip);
    seg.ma_interfaces = UnpackInts(b, ip);
    seg.truncation_attributes = UnpackInts(b, ip);
    seg.port_attributes = UnpackInts(b, ip);
    seg.face_count = static_cast<int>(b.ints[ip++]);
    seg.face_normals.resize(static_cast<std::size_t>(b.ints[ip++]));
    for (auto &n : seg.face_normals)
    {
      std::copy_n(b.reals.data() + rp, 3, n.begin());
      rp += 3;
    }
    seg.side_attributes = UnpackInts(b, ip);
    seg.on_bounding_box = b.ints[ip++] != 0;
  }
  MFEM_VERIFY(ip == b.ints.size() && rp == b.reals.size(),
              "Broadcast metal perimeter buffers were not consumed exactly!");
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
  // Stage lines (counts and wall time on the root): the extraction is replicated on every
  // rank and a chip-scale mesh takes minutes here before the identification starts.
  const auto extraction_started = std::chrono::steady_clock::now();
  auto StageLine = [&](const std::string &text)
  {
    Mpi::Print(
        "  Metal perimeter {}: ({:.2f} s)\n", text,
        std::chrono::duration<double>(std::chrono::steady_clock::now() - extraction_started)
            .count());
  };
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
  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
  double extent = 0.0;
  for (int d = 0; d < 3; d++)
  {
    extent = std::max(extent, bbmax[d] - bbmin[d]);
  }
  MFEM_VERIFY(extent > 0.0, "Degenerate mesh geometry for metal edge extraction!");
  const double tolerance_squared = 1.0e-18 * extent * extent;
  const double coordinate_tolerance = 1.0e-10 * extent;

  using Point = std::array<double, 3>;
  using PointKey = std::array<long long int, 3>;

  // (1) Every rank collects its metal boundary faces: attribute, the two adjacent domain
  // element attributes (-1 when the face is exterior or the neighbour is absent), and the
  // ordered polygon of the face with high-order edges sampled as in the perimeter.
  std::vector<double> local_face_data;
  std::size_t local_face_count = 0;
  {
    auto &mutable_mesh = const_cast<mfem::ParMesh &>(mesh);
    mutable_mesh.ExchangeFaceNbrData();
    mfem::L2_FECollection material_fec(0, mesh.Dimension());
    mfem::ParFiniteElementSpace material_fespace(&mutable_mesh, &material_fec);
    mfem::ParGridFunction material_attribute(&material_fespace);
    mfem::Array<int> dofs;
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      material_fespace.GetElementDofs(element, dofs);
      material_attribute[dofs[0]] = mesh.GetAttribute(element);
    }
    material_attribute.ExchangeFaceNbrData();

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
      const int face = mesh.GetBdrElementFaceIndex(be);
      const auto face_info = mesh.GetFaceInformation(face);
      int side_1 = mesh.GetAttribute(face_info.element[0].index);
      int side_2 = -1;
      if (face_info.element[1].location == mfem::Mesh::ElementLocation::Local)
      {
        side_2 = mesh.GetAttribute(face_info.element[1].index);
      }
      else if (face_info.element[1].location == mfem::Mesh::ElementLocation::FaceNbr)
      {
        material_fespace.GetFaceNbrElementVDofs(face_info.element[1].index, dofs);
        side_2 = static_cast<int>(std::llround(material_attribute.FaceNbrData()[dofs[0]]));
      }

      std::vector<Point> loop;
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
          if (loop.empty())
          {
            loop.push_back(segment.p0);
          }
          else
          {
            double distance_squared = 0.0;
            for (int d = 0; d < 3; d++)
            {
              const double delta = loop.back()[d] - segment.p0[d];
              distance_squared += delta * delta;
            }
            MFEM_VERIFY(distance_squared <= tolerance_squared,
                        "Metal boundary face edges do not form an ordered loop!");
          }
          loop.push_back(segment.p1);
        }
      }
      double closure_distance_squared = 0.0;
      for (int d = 0; d < 3; d++)
      {
        const double delta = loop.front()[d] - loop.back()[d];
        closure_distance_squared += delta * delta;
      }
      MFEM_VERIFY(closure_distance_squared <= tolerance_squared,
                  "Metal boundary face edges do not form a closed loop!");
      loop.pop_back();
      if (loop.size() == static_cast<std::size_t>(vertices.Size()))
      {
        // Geometrically linear face: keep the exact vertex representation.
        for (int i = 0; i < vertices.Size(); i++)
        {
          loop[i] = mesh::GetVertexCoordinates(mesh, vertices[i]);
        }
      }
      local_face_data.push_back(static_cast<double>(attribute));
      local_face_data.push_back(static_cast<double>(side_1));
      local_face_data.push_back(static_cast<double>(side_2));
      local_face_data.push_back(static_cast<double>(loop.size()));
      for (const auto &point : loop)
      {
        local_face_data.insert(local_face_data.end(), point.begin(), point.end());
      }
      local_face_count++;
    }
  }

  // (2) Gather the metal faces on the root only (decision 82 infrastructure): the global
  // perimeter is built once, on the root, and the compact result (segments, vertices,
  // each rank's own retained facets) is broadcast / scattered below. The gathered faces,
  // canonical points, face incidence and global faces exist on the root alone (a
  // chip-scale mesh replicated ~3 GB of them on every rank: 560 GiB at np192).
  MPI_Comm comm = mesh.GetComm();
  const bool root = Mpi::Root(comm);
  MFEM_VERIFY(local_face_data.size() <=
                  static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "Local metal face data exceeds the MPI count limit!");
  const int local_count = static_cast<int>(local_face_data.size());
  std::vector<int> counts(Mpi::Size(comm));
  Mpi::Allgather(1, &local_count, counts.data(), comm);
  std::vector<int> offsets(counts.size());
  int total = 0;
  for (std::size_t rank = 0; rank < counts.size(); rank++)
  {
    offsets[rank] = total;
    MFEM_VERIFY(counts[rank] <= std::numeric_limits<int>::max() - total,
                "Global metal face data exceeds the MPI count limit!");
    total += counts[rank];
  }
  std::vector<double> face_data(root ? total : 0);
  MPI_Gatherv(local_face_data.data(), local_count, mpi::DataType<double>(), face_data.data(),
              counts.data(), offsets.data(), mpi::DataType<double>(), 0, comm);
  local_face_data.clear();
  local_face_data.shrink_to_fit();
  if (total == 0)
  {
    return result;
  }

  struct GatheredFace
  {
    int attribute;
    std::array<int, 2> sides;
    std::vector<Point> loop;
  };
  std::vector<GatheredFace> gathered;
  std::vector<int> gathered_rank;  // source rank of every gathered face (root only)
  {
    std::size_t rank = 0;
    for (std::size_t offset = 0; offset < face_data.size();)
    {
      while (rank + 1 < counts.size() &&
             offset >= static_cast<std::size_t>(offsets[rank] + counts[rank]))
      {
        rank++;
      }
      GatheredFace face;
      face.attribute = static_cast<int>(std::llround(face_data[offset]));
      face.sides = {static_cast<int>(std::llround(face_data[offset + 1])),
                    static_cast<int>(std::llround(face_data[offset + 2]))};
      const auto n = static_cast<std::size_t>(std::llround(face_data[offset + 3]));
      offset += 4;
      face.loop.resize(n);
      for (std::size_t i = 0; i < n; i++)
      {
        std::copy_n(face_data.data() + offset + 3 * i, 3, face.loop[i].begin());
      }
      offset += 3 * n;
      gathered.push_back(std::move(face));
      gathered_rank.push_back(static_cast<int>(rank));
    }
  }
  face_data.clear();
  face_data.shrink_to_fit();
  StageLine(std::to_string(local_face_count) + " local metal faces, " +
            std::to_string(gathered.size()) + " gathered on the root");

  // Collective preliminaries of the perimeter classification (every rank takes part in
  // the gathers inside; the trees are used on the root only): the dielectric-interface
  // support trees and the simulation-cut (truncation) surfaces.
  // An interface defined on metal attributes (MS / MA on the metal itself) supports every
  // perimeter segment of those attributes; an interface on nonmetal attributes (the SA
  // sheet) supports the segments coincident with its own perimeter.
  struct InterfaceSupport
  {
    int index;
    InterfaceDielectric type;
    std::set<int> metal_attributes;
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
    InterfaceSupport support{index, dielectric.type, {}, nullptr};
    std::vector<int> nonmetal;
    for (const int attribute : dielectric.attributes)
    {
      if (attribute_conditions.find(attribute) != attribute_conditions.end())
      {
        support.metal_attributes.insert(attribute);
      }
      else
      {
        nonmetal.push_back(attribute);
      }
    }
    std::sort(nonmetal.begin(), nonmetal.end());
    nonmetal.erase(std::unique(nonmetal.begin(), nonmetal.end()), nonmetal.end());
    if (!nonmetal.empty())
    {
      auto [tree_it, inserted] = interface_support_trees.try_emplace(nonmetal);
      if (inserted)
      {
        tree_it->second = BuildSupportTree(mesh, nonmetal, true);
      }
      support.tree = tree_it->second;
    }
    if (support.tree || !support.metal_attributes.empty())
    {
      interface_support.push_back(std::move(support));
    }
  }

  // An exterior, nonmetal boundary surface which is not itself a dielectric interface is
  // a simulation cut surface. Internal lumped ports and sources are deliberately omitted.
  // Where an exterior face edge coincides with the metal perimeter, the metal has been
  // cut by the simulation domain (for example at a wave port or an outer box).
  // Ports are not metal (decision 82(5)): the metal perimeter bordering a LumpedPort /
  // WavePort boundary face (the port attributes of the configuration, whatever their names)
  // is a cut like a truncation, reported as the Port exclusion.
  std::map<int, std::shared_ptr<const EdgeDistanceTree>> port_support;
  {
    std::set<int> port_attributes;
    for (const auto &[index, port] : boundaries.lumpedport)
    {
      (void)index;
      for (const auto &element : port.elements)
      {
        port_attributes.insert(element.attributes.begin(), element.attributes.end());
      }
    }
    for (const auto &[index, port] : boundaries.waveport)
    {
      (void)index;
      port_attributes.insert(port.attributes.begin(), port.attributes.end());
    }
    for (const int attribute : port_attributes)
    {
      if (attribute_conditions.find(attribute) != attribute_conditions.end())
      {
        continue;  // a metal attribute is metal, not a port face
      }
      auto tree = BuildSupportTree(mesh, std::vector<int>{attribute}, false, false);
      if (tree)
      {
        port_support.emplace(attribute, std::move(tree));
      }
    }
  }
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


  // Global face of every gathered face (root only): the per-rank retained facets below.
  std::vector<std::size_t> gathered_face;
  // Per-rank retained facets (root only), packed for the scatter; this rank's share.
  std::vector<GeometryBuffers> rank_faces(counts.size());
  GeometryBuffers local_faces;
  if (root)
  {
  // (3) Canonical points: coordinates within coordinate_tolerance of an already registered
  // point (the 27 neighbouring grid cells are searched so a crack copy straddling a cell
  // boundary still merges) map to that point. The representative is the lexicographically
  // smallest copy; the numbering is fixed afterwards by sorting the distinct points.
  std::vector<Point> canonical_points;
  std::map<PointKey, std::vector<std::size_t>> point_cells;
  auto GetPointKey = [&](const Point &point)
  {
    PointKey key;
    for (int d = 0; d < 3; d++)
    {
      key[d] = std::llround((point[d] - bbmin[d]) / coordinate_tolerance);
    }
    return key;
  };
  auto CanonicalPoint = [&](const Point &point) -> std::size_t
  {
    const PointKey key = GetPointKey(point);
    for (long long int dx = -1; dx <= 1; dx++)
    {
      for (long long int dy = -1; dy <= 1; dy++)
      {
        for (long long int dz = -1; dz <= 1; dz++)
        {
          const auto cell = point_cells.find({key[0] + dx, key[1] + dy, key[2] + dz});
          if (cell == point_cells.end())
          {
            continue;
          }
          for (const std::size_t candidate : cell->second)
          {
            double distance_squared = 0.0;
            for (int d = 0; d < 3; d++)
            {
              const double delta = canonical_points[candidate][d] - point[d];
              distance_squared += delta * delta;
            }
            if (distance_squared <= coordinate_tolerance * coordinate_tolerance)
            {
              // Representative = the lexicographically smallest copy merged into the
              // point (a function of the set of copies, not of their gathered order; the
              // cell of the representative moves by at most one, within the search).
              if (point < canonical_points[candidate])
              {
                canonical_points[candidate] = point;
              }
              return candidate;
            }
          }
        }
      }
    }
    canonical_points.push_back(point);
    point_cells[key].push_back(canonical_points.size() - 1);
    return canonical_points.size() - 1;
  };

  // Deduplicate coincident faces (crack copies, duplicate boundary elements of one
  // attribute): a face is identified by the set of its canonical vertices.
  struct GlobalFace
  {
    int attribute = -1;
    std::set<int> sides;
    std::vector<std::size_t> loop;  // canonical point indices, ordered
    Point normal{};
    Point centroid{};
    double area = 0.0;
    bool on_bounding_box = false;
  };
  std::vector<GlobalFace> faces;
  std::map<std::vector<std::size_t>, std::size_t> face_by_vertices;
  {
    // Register every point (first pass, in the gathered order), then renumber the
    // canonical points by their sorted quantized coordinates: the point, segment and vertex
    // numbering of the perimeter (and hence of the manifest) is a function of the set of
    // distinct points only, not of the partition or of the order in which the crack copies
    // of a face were gathered (the first-encounter numbering permuted 2,647 segment-table
    // entries between np8 and np192 on DS-SCT-002). Ties on the quantized key (distinct
    // points closer than sqrt(3) tolerance in one cell) are broken by the raw coordinates.
    std::vector<std::vector<std::size_t>> gathered_loops(gathered.size());
    for (std::size_t g = 0; g < gathered.size(); g++)
    {
      auto &loop = gathered_loops[g];
      loop.resize(gathered[g].loop.size());
      for (std::size_t i = 0; i < loop.size(); i++)
      {
        loop[i] = CanonicalPoint(gathered[g].loop[i]);
      }
    }
    std::vector<std::size_t> order(canonical_points.size());
    std::iota(order.begin(), order.end(), 0);
    std::vector<PointKey> keys(canonical_points.size());
    for (std::size_t p = 0; p < canonical_points.size(); p++)
    {
      keys[p] = GetPointKey(canonical_points[p]);
    }
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b)
              {
                return std::tie(keys[a], canonical_points[a]) <
                       std::tie(keys[b], canonical_points[b]);
              });
    std::vector<std::size_t> new_index(canonical_points.size());
    std::vector<Point> sorted_points(canonical_points.size());
    for (std::size_t rank_of_point = 0; rank_of_point < order.size(); rank_of_point++)
    {
      new_index[order[rank_of_point]] = rank_of_point;
      sorted_points[rank_of_point] = canonical_points[order[rank_of_point]];
    }
    canonical_points = std::move(sorted_points);
    for (auto &[cell, members] : point_cells)
    {
      (void)cell;
      for (auto &member : members)
      {
        member = new_index[member];
      }
      std::sort(members.begin(), members.end());
    }
    // Faces in the order of their sorted canonical vertex sets (then loop, attribute,
    // sides): the face numbering, the metal components and the face incidence order of
    // every segment follow the point numbering.
    std::vector<std::vector<std::size_t>> face_keys(gathered.size());
    for (std::size_t g = 0; g < gathered.size(); g++)
    {
      for (auto &p : gathered_loops[g])
      {
        p = new_index[p];
      }
      face_keys[g] = gathered_loops[g];
      std::sort(face_keys[g].begin(), face_keys[g].end());
      face_keys[g].erase(std::unique(face_keys[g].begin(), face_keys[g].end()),
                         face_keys[g].end());
      MFEM_VERIFY(face_keys[g].size() >= 3, "A metal boundary face degenerates to fewer "
                                            "than three distinct vertices!");
    }
    std::vector<std::size_t> face_order(gathered.size());
    std::iota(face_order.begin(), face_order.end(), 0);
    std::sort(face_order.begin(), face_order.end(),
              [&](std::size_t a, std::size_t b)
              {
                return std::tie(face_keys[a], gathered_loops[a], gathered[a].attribute,
                                gathered[a].sides) <
                       std::tie(face_keys[b], gathered_loops[b], gathered[b].attribute,
                                gathered[b].sides);
              });
    gathered_face.resize(gathered.size());
    for (const std::size_t g : face_order)
    {
      const auto &source = gathered[g];
      auto [entry, inserted] = face_by_vertices.try_emplace(face_keys[g], faces.size());
      gathered_face[g] = entry->second;
      if (inserted)
      {
        GlobalFace face;
        face.attribute = source.attribute;
        face.loop = std::move(gathered_loops[g]);
        faces.push_back(std::move(face));
      }
      else
      {
        auto &face = faces[entry->second];
        if (face.attribute != source.attribute)
        {
          const auto &p = source.loop.front();
          MFEM_ABORT("Coincident metal boundary elements carry different metal attributes "
                     << face.attribute << " and " << source.attribute << " (face at ("
                     << p[0] << ", " << p[1] << ", " << p[2]
                     << ")); drop or merge the duplicate boundary elements!");
        }
      }
      auto &face = faces[entry->second];
      for (const int side : source.sides)
      {
        if (side > 0)
        {
          face.sides.insert(side);
        }
      }
    }
  }
  StageLine(std::to_string(faces.size()) + " distinct faces, " +
            std::to_string(canonical_points.size()) + " canonical points");

  // Face normals (Newell), centroids and areas on the canonical loops; the layer normal is
  // the area-weighted principal direction of the face normals. Metal faces lying on the
  // bounding box of the mesh (a PEC simulation box) are not process metal and do not vote.
  auto OnBoundingBox = [&](const GlobalFace &face)
  {
    for (int d = 0; d < 3; d++)
    {
      for (const double bound : {bbmin[d], bbmax[d]})
      {
        if (std::all_of(face.loop.begin(), face.loop.end(), [&](std::size_t p)
                        { return std::abs(canonical_points[p][d] - bound) <= coordinate_tolerance; }))
        {
          return true;
        }
      }
    }
    return false;
  };
  std::array<double, 6> normal_tensor{};  // xx, yy, zz, xy, xz, yz
  std::array<double, 6> box_normal_tensor{};
  for (auto &face : faces)
  {
    Point normal{};
    Point centroid{};
    for (std::size_t i = 0; i < face.loop.size(); i++)
    {
      const auto &a = canonical_points[face.loop[i]];
      const auto &b = canonical_points[face.loop[(i + 1) % face.loop.size()]];
      normal[0] += (a[1] - b[1]) * (a[2] + b[2]);
      normal[1] += (a[2] - b[2]) * (a[0] + b[0]);
      normal[2] += (a[0] - b[0]) * (a[1] + b[1]);
      for (int d = 0; d < 3; d++)
      {
        centroid[d] += a[d] / static_cast<double>(face.loop.size());
      }
    }
    double norm_squared = 0.0;
    for (const double value : normal)
    {
      norm_squared += value * value;
    }
    MFEM_VERIFY(norm_squared > 0.0, "Degenerate metal boundary face!");
    const double norm = std::sqrt(norm_squared);
    face.area = 0.5 * norm;
    for (double &value : normal)
    {
      value /= norm;
    }
    // Sign-canonical: the largest-magnitude component positive.
    int dominant = 0;
    for (int d = 1; d < 3; d++)
    {
      if (std::abs(normal[d]) > std::abs(normal[dominant]) + 1.0e-12)
      {
        dominant = d;
      }
    }
    if (normal[dominant] < 0.0)
    {
      for (double &value : normal)
      {
        value = -value;
      }
    }
    face.normal = normal;
    face.centroid = centroid;
    face.on_bounding_box = OnBoundingBox(face);
    auto &tensor = face.on_bounding_box ? box_normal_tensor : normal_tensor;
    tensor[0] += face.area * normal[0] * normal[0];
    tensor[1] += face.area * normal[1] * normal[1];
    tensor[2] += face.area * normal[2] * normal[2];
    tensor[3] += face.area * normal[0] * normal[1];
    tensor[4] += face.area * normal[0] * normal[2];
    tensor[5] += face.area * normal[1] * normal[2];
  }
  if (normal_tensor[0] + normal_tensor[1] + normal_tensor[2] <= 0.0)
  {
    normal_tensor = box_normal_tensor;  // only box metal: nothing better to vote
  }
  {
    // Power iteration for the dominant eigenvector of the symmetric 3 x 3 tensor.
    const double trace = normal_tensor[0] + normal_tensor[1] + normal_tensor[2];
    Point v{};
    int start = 0;
    for (int d = 1; d < 3; d++)
    {
      if (normal_tensor[d] > normal_tensor[start])
      {
        start = d;
      }
    }
    v[start] = 1.0;
    for (int iteration = 0; iteration < 200; iteration++)
    {
      Point w = {normal_tensor[0] * v[0] + normal_tensor[3] * v[1] + normal_tensor[4] * v[2],
                 normal_tensor[3] * v[0] + normal_tensor[1] * v[1] + normal_tensor[5] * v[2],
                 normal_tensor[4] * v[0] + normal_tensor[5] * v[1] + normal_tensor[2] * v[2]};
      // Shift keeps the iteration well conditioned for a rank-one tensor.
      for (int d = 0; d < 3; d++)
      {
        w[d] += 1.0e-3 * trace * v[d];
      }
      double norm = 0.0;
      for (const double value : w)
      {
        norm += value * value;
      }
      norm = std::sqrt(norm);
      MFEM_VERIFY(norm > 0.0, "Degenerate metal face normal tensor!");
      for (int d = 0; d < 3; d++)
      {
        v[d] = w[d] / norm;
      }
    }
    int dominant = 0;
    for (int d = 1; d < 3; d++)
    {
      if (std::abs(v[d]) > std::abs(v[dominant]) + 1.0e-12)
      {
        dominant = d;
      }
    }
    if (v[dominant] < 0.0)
    {
      for (double &value : v)
      {
        value = -value;
      }
    }
    result.layer_normal = v;
  }

  // (4) Atomic edge segments of the faces, split at canonical points lying on them when
  // the mesh is nonconforming (crack sides with different hanging-node subdivisions), and
  // the metal faces supporting each of them with their in-plane inward directions.
  using SegmentKey = std::pair<std::size_t, std::size_t>;
  struct FaceIncidence
  {
    std::size_t face;
    Point inward;
  };
  std::map<SegmentKey, std::vector<FaceIncidence>> incidence;
  {
    std::map<PointKey, std::vector<std::size_t>> cells;  // for the nonconforming split
    if (mesh.Nonconforming())
    {
      for (std::size_t p = 0; p < canonical_points.size(); p++)
      {
        cells[GetPointKey(canonical_points[p])].push_back(p);
      }
    }
    auto PointsOnSegment = [&](std::size_t a, std::size_t b)
    {
      std::vector<std::pair<double, std::size_t>> found;
      const auto &p0 = canonical_points[a];
      const auto &p1 = canonical_points[b];
      Point direction{};
      double length_squared = 0.0;
      PointKey lo, hi;
      for (int d = 0; d < 3; d++)
      {
        direction[d] = p1[d] - p0[d];
        length_squared += direction[d] * direction[d];
        lo[d] = std::llround((std::min(p0[d], p1[d]) - bbmin[d]) / coordinate_tolerance) - 1;
        hi[d] = std::llround((std::max(p0[d], p1[d]) - bbmin[d]) / coordinate_tolerance) + 1;
      }
      MFEM_VERIFY(length_squared > 0.0, "Degenerate metal face edge!");
      for (auto it = cells.lower_bound({lo[0], lo[1], lo[2]});
           it != cells.end() && it->first[0] <= hi[0]; ++it)
      {
        if (it->first[1] < lo[1] || it->first[1] > hi[1] || it->first[2] < lo[2] ||
            it->first[2] > hi[2])
        {
          continue;
        }
        for (const std::size_t candidate : it->second)
        {
          if (candidate == a || candidate == b)
          {
            continue;
          }
          const auto &point = canonical_points[candidate];
          double projection = 0.0;
          for (int d = 0; d < 3; d++)
          {
            projection += (point[d] - p0[d]) * direction[d];
          }
          const double t = projection / length_squared;
          if (t <= 0.0 || t >= 1.0)
          {
            continue;
          }
          double distance_squared = 0.0;
          for (int d = 0; d < 3; d++)
          {
            const double delta = point[d] - p0[d] - t * direction[d];
            distance_squared += delta * delta;
          }
          if (distance_squared <= coordinate_tolerance * coordinate_tolerance)
          {
            found.emplace_back(t, candidate);
          }
        }
      }
      std::sort(found.begin(), found.end());
      return found;
    };
    for (std::size_t f = 0; f < faces.size(); f++)
    {
      const auto &face = faces[f];
      for (std::size_t i = 0; i < face.loop.size(); i++)
      {
        const std::size_t a = face.loop[i];
        const std::size_t b = face.loop[(i + 1) % face.loop.size()];
        if (a == b)
        {
          continue;
        }
        std::vector<std::size_t> chain = {a};
        if (mesh.Nonconforming())
        {
          for (const auto &[t, p] : PointsOnSegment(a, b))
          {
            (void)t;
            chain.push_back(p);
          }
        }
        chain.push_back(b);
        for (std::size_t k = 0; k + 1 < chain.size(); k++)
        {
          const std::size_t u = chain[k];
          const std::size_t v = chain[k + 1];
          const auto &p0 = canonical_points[u];
          const auto &p1 = canonical_points[v];
          Point tangent{}, inward{};
          double tangent_norm_squared = 0.0;
          for (int d = 0; d < 3; d++)
          {
            tangent[d] = p1[d] - p0[d];
            tangent_norm_squared += tangent[d] * tangent[d];
            inward[d] = face.centroid[d] - 0.5 * (p0[d] + p1[d]);
          }
          double inward_tangent = 0.0;
          for (int d = 0; d < 3; d++)
          {
            inward_tangent += inward[d] * tangent[d];
          }
          double inward_norm_squared = 0.0;
          for (int d = 0; d < 3; d++)
          {
            inward[d] -= inward_tangent * tangent[d] / tangent_norm_squared;
            inward_norm_squared += inward[d] * inward[d];
          }
          MFEM_VERIFY(inward_norm_squared > 0.0,
                      "Degenerate metal boundary face (centroid on an edge)!");
          const double inverse_norm = 1.0 / std::sqrt(inward_norm_squared);
          for (double &value : inward)
          {
            value *= inverse_norm;
          }
          incidence[{std::min(u, v), std::max(u, v)}].push_back({f, inward});
        }
      }
    }
  }

  StageLine(std::to_string(incidence.size()) + " face edges");

  // (5) Metal components through shared face edges.
  std::vector<std::size_t> face_parent(faces.size());
  std::iota(face_parent.begin(), face_parent.end(), 0);
  auto FindFace = [&](std::size_t item)
  {
    std::size_t root = item;
    while (face_parent[root] != root)
    {
      root = face_parent[root];
    }
    while (face_parent[item] != item)
    {
      const std::size_t next = face_parent[item];
      face_parent[item] = root;
      item = next;
    }
    return root;
  };
  for (const auto &[key, supports] : incidence)
  {
    (void)key;
    for (std::size_t i = 1; i < supports.size(); i++)
    {
      const std::size_t first = FindFace(supports[0].face);
      const std::size_t second = FindFace(supports[i].face);
      if (first != second)
      {
        face_parent[std::max(first, second)] = std::min(first, second);
      }
    }
  }
  std::map<std::size_t, int> component_by_root;
  std::vector<int> face_component(faces.size());
  for (std::size_t f = 0; f < faces.size(); f++)
  {
    auto [entry, inserted] = component_by_root.try_emplace(FindFace(f), component_by_root.size());
    (void)inserted;
    face_component[f] = entry->second;
  }
  result.metal_components = static_cast<int>(component_by_root.size());

  // (6) Perimeter classification. Faces whose inward directions coincide (crack copies with
  // different subdivisions, coarse and fine sides) form one direction class; a segment with
  // an opposite pair of classes and nothing else is interior to the metal.
  // The direction grid is the classification's 1e-12 direction quantum.
  constexpr double direction_quantum = 1.0e-12;
  auto QuantizeCosine = [](double cosine) { return std::round(cosine / direction_quantum); };
  const double same_cosine = QuantizeCosine(1.0 - 1.0e-8);
  const double opposite_cosine = QuantizeCosine(-1.0 + 1.0e-8);
  std::map<std::size_t, std::size_t> vertex_by_point;
  auto GetVertex = [&](std::size_t point)
  {
    auto [it, inserted] = vertex_by_point.try_emplace(point, result.vertices.size());
    if (inserted)
    {
      result.vertices.push_back(
          {canonical_points[point], {}, MetalEdgeVertexType::REGULAR});
    }
    return it->second;
  };
  for (const auto &[key, supports] : incidence)
  {
    // Direction classes.
    std::vector<Point> classes;
    std::vector<std::size_t> class_faces;
    for (const auto &support : supports)
    {
      bool known = false;
      for (const auto &direction : classes)
      {
        double dot = 0.0;
        for (int d = 0; d < 3; d++)
        {
          dot += direction[d] * support.inward[d];
        }
        if (QuantizeCosine(dot) >= same_cosine)
        {
          known = true;
          break;
        }
      }
      if (!known)
      {
        classes.push_back(support.inward);
      }
    }
    bool opposite_pair = false;
    for (std::size_t i = 0; i < classes.size() && !opposite_pair; i++)
    {
      for (std::size_t j = i + 1; j < classes.size(); j++)
      {
        double dot = 0.0;
        for (int d = 0; d < 3; d++)
        {
          dot += classes[i][d] * classes[j][d];
        }
        if (QuantizeCosine(dot) <= opposite_cosine)
        {
          opposite_pair = true;
          break;
        }
      }
    }
    MetalEdgeSegmentType type;
    if (classes.size() == 1)
    {
      type = MetalEdgeSegmentType::PHYSICAL;
    }
    else if (classes.size() == 2)
    {
      if (opposite_pair)
      {
        continue;  // metal continues on both sides: interior edge
      }
      type = MetalEdgeSegmentType::FOLD;
    }
    else
    {
      type = MetalEdgeSegmentType::NONMANIFOLD;
    }

    MetalEdgeSegment segment;
    segment.type = type;
    segment.vertices = {GetVertex(key.first), GetVertex(key.second)};
    std::set<std::size_t> distinct_faces;
    std::set<int> side_attributes;
    for (const auto &support : supports)
    {
      distinct_faces.insert(support.face);
    }
    segment.face_count = static_cast<int>(distinct_faces.size());
    segment.on_bounding_box =
        std::all_of(distinct_faces.begin(), distinct_faces.end(),
                    [&](std::size_t f) { return faces[f].on_bounding_box; });
    std::set<int> attributes;
    for (const std::size_t f : distinct_faces)
    {
      attributes.insert(faces[f].attribute);
      side_attributes.insert(faces[f].sides.begin(), faces[f].sides.end());
      segment.face_normals.push_back(faces[f].normal);
    }
    std::sort(segment.face_normals.begin(), segment.face_normals.end());
    segment.face_normals.erase(
        std::unique(segment.face_normals.begin(), segment.face_normals.end()),
        segment.face_normals.end());
    segment.metal_attributes.assign(attributes.begin(), attributes.end());
    segment.side_attributes.assign(side_attributes.begin(), side_attributes.end());
    for (const int attribute : segment.metal_attributes)
    {
      const auto &conditions = attribute_conditions.at(attribute);
      segment.conditions.insert(segment.conditions.end(), conditions.begin(),
                                conditions.end());
    }
    SortAndUnique(segment.conditions);
    segment.metal_component = face_component[*distinct_faces.begin()];

    const mesh::BoundaryEdgeSegment perimeter{canonical_points[key.first],
                                              canonical_points[key.second]};
    for (const auto &support : interface_support)
    {
      const bool by_attribute = std::any_of(
          segment.metal_attributes.begin(), segment.metal_attributes.end(),
          [&](int attribute)
          { return support.metal_attributes.find(attribute) != support.metal_attributes.end(); });
      if (!by_attribute &&
          !(support.tree && IsCoincident(perimeter, *support.tree, tolerance_squared)))
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
    if (type == MetalEdgeSegmentType::PHYSICAL)
    {
      for (const auto &[attribute, tree] : port_support)
      {
        if (IsCoincident(perimeter, *tree, tolerance_squared))
        {
          segment.port_attributes.push_back(attribute);
        }
      }
      for (const auto &[attribute, tree] : truncation_support)
      {
        if (IsCoincident(perimeter, *tree, tolerance_squared))
        {
          segment.truncation_attributes.push_back(attribute);
        }
      }
      // A port face on the simulation boundary (a wave port) is a port cut.
      if (!segment.port_attributes.empty())
      {
        segment.type = MetalEdgeSegmentType::PORT;
      }
      else if (!segment.truncation_attributes.empty())
      {
        segment.type = MetalEdgeSegmentType::TRUNCATION;
      }
    }

    const std::size_t index = result.segments.size();
    result.vertices[segment.vertices[0]].segments.push_back(index);
    result.vertices[segment.vertices[1]].segments.push_back(index);
    result.segments.push_back(std::move(segment));
  }
  StageLine(std::to_string(result.segments.size()) + " perimeter segments, " +
            std::to_string(result.vertices.size()) + " vertices, " +
            std::to_string(result.metal_components) + " metal components");

  // (7) Retained faces: the rank-local facets (with their global component) and the
  // global deduplicated faces. A rank retains every distinct geometric face it owns once,
  // with the canonical (global) vertex coordinates, so that the crack copies and duplicate
  // boundary elements of one face — whose own vertex coordinates differ by roundoff — are
  // one facet with identical coordinates on every rank (the plan-view canonicalisation
  // deduplicates facets on a 1e-9 R grid).
  if (surface.retain_faces && !result.segments.empty())
  {
    // The gathered faces of every rank in their local order (the rank's boundary element
    // order), each distinct global face once per rank, with the loop's canonical points.
    std::vector<std::set<std::size_t>> retained(counts.size());
    for (std::size_t g = 0; g < gathered.size(); g++)
    {
      const std::size_t f = gathered_face[g];
      const auto rank = static_cast<std::size_t>(gathered_rank[g]);
      if (!retained[rank].insert(f).second)
      {
        continue;
      }
      MetalSurfaceFace face;
      face.component = face_component[f];
      face.attribute = faces[f].attribute;
      face.normal = faces[f].normal;
      face.on_bounding_box = faces[f].on_bounding_box;
      face.vertices.reserve(gathered[g].loop.size());
      for (const auto &p : gathered[g].loop)
      {
        face.vertices.push_back(canonical_points[CanonicalPoint(p)]);
      }
      PackFace(face, rank_faces[rank]);
    }
  }
  gathered.clear();
  gathered.shrink_to_fit();
  if (surface.retain_global_faces && !result.segments.empty())
  {
    result.global_faces.reserve(faces.size());
    for (std::size_t f = 0; f < faces.size(); f++)
    {
      MetalSurfaceFace face;
      face.component = face_component[f];
      face.attribute = faces[f].attribute;
      face.normal = faces[f].normal;
      face.on_bounding_box = faces[f].on_bounding_box;
      face.vertices.reserve(faces[f].loop.size());
      for (const std::size_t p : faces[f].loop)
      {
        face.vertices.push_back(canonical_points[p]);
      }
      result.global_faces.push_back(std::move(face));
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
  // The turn test compares direction cosines on a fixed 1e-12 grid so that a roundoff-level
  // perturbation of the vertex coordinates cannot flip a vertex between REGULAR and CORNER
  // (the same direction quantum as the classification's parallelism tests).
  auto QuantizeDirection = QuantizeCosine;
  const double quantized_straight_dot_tolerance = QuantizeDirection(straight_dot_tolerance);
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
    return QuantizeDirection(dot) <= quantized_straight_dot_tolerance
               ? MetalEdgeVertexType::REGULAR
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
    vertex.on_port_boundary = std::any_of(
        vertex.segments.begin(), vertex.segments.end(), [&](std::size_t segment)
        { return result.segments[segment].type == MetalEdgeSegmentType::PORT; });
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
  StageLine(std::to_string(result.physical_chains) + " physical chains, " +
            std::to_string(result.global_faces.size()) + " global faces retained");
  }  // root

  // Broadcast the compact perimeter (segments, vertices, scalars; no faces) and scatter
  // every rank's own retained facets. An empty perimeter is empty everywhere.
  {
    GeometryBuffers buffers;
    if (root)
    {
      PackGeometry(result, buffers);
    }
    std::array<long long int, 2> sizes = {static_cast<long long int>(buffers.ints.size()),
                                          static_cast<long long int>(buffers.reals.size())};
    Mpi::Broadcast(2, sizes.data(), 0, comm);
    if (!root)
    {
      buffers.ints.resize(static_cast<std::size_t>(sizes[0]));
      buffers.reals.resize(static_cast<std::size_t>(sizes[1]));
    }
    Mpi::BroadcastLarge(sizes[0], buffers.ints.data(), 0, comm);
    Mpi::BroadcastLarge(sizes[1], buffers.reals.data(), 0, comm);
    if (!root)
    {
      UnpackGeometry(buffers, result);
    }
  }
  if (surface.retain_faces)
  {
    for (const bool reals : {false, true})
    {
      std::vector<int> face_counts(counts.size(), 0), face_offsets(counts.size(), 0);
      std::vector<long long int> send_ints;
      std::vector<double> send_reals;
      if (root)
      {
        long long int offset = 0;
        for (std::size_t rank = 0; rank < counts.size(); rank++)
        {
          const std::size_t n =
              reals ? rank_faces[rank].reals.size() : rank_faces[rank].ints.size();
          MFEM_VERIFY(offset + static_cast<long long int>(n) <=
                          std::numeric_limits<int>::max(),
                      "Retained metal facet data exceeds the MPI count limit!");
          face_offsets[rank] = static_cast<int>(offset);
          face_counts[rank] = static_cast<int>(n);
          offset += static_cast<long long int>(n);
          if (reals)
          {
            send_reals.insert(send_reals.end(), rank_faces[rank].reals.begin(),
                              rank_faces[rank].reals.end());
            rank_faces[rank].reals.clear();
          }
          else
          {
            send_ints.insert(send_ints.end(), rank_faces[rank].ints.begin(),
                             rank_faces[rank].ints.end());
            rank_faces[rank].ints.clear();
          }
        }
      }
      int my_count = 0;
      MPI_Scatter(face_counts.data(), 1, mpi::DataType<int>(), &my_count, 1,
                  mpi::DataType<int>(), 0, comm);
      if (reals)
      {
        local_faces.reals.resize(static_cast<std::size_t>(my_count));
        MPI_Scatterv(send_reals.data(), face_counts.data(), face_offsets.data(),
                     mpi::DataType<double>(), local_faces.reals.data(), my_count,
                     mpi::DataType<double>(), 0, comm);
      }
      else
      {
        local_faces.ints.resize(static_cast<std::size_t>(my_count));
        MPI_Scatterv(send_ints.data(), face_counts.data(), face_offsets.data(),
                     mpi::DataType<long long int>(), local_faces.ints.data(), my_count,
                     mpi::DataType<long long int>(), 0, comm);
      }
    }
    std::size_t int_pos = 0, real_pos = 0;
    while (int_pos < local_faces.ints.size())
    {
      result.surface_faces.push_back(UnpackFace(local_faces, int_pos, real_pos));
    }
  }
  return result;
}

std::vector<std::array<double, 3>>
BuildMetalEdgeProcessNormals(const mfem::ParMesh &mesh, const MetalEdgeGeometry &geometry,
                             const std::vector<std::size_t> &segment_indices,
                             const std::function<double(int)> &material_score,
                             const std::optional<std::array<double, 3>> &fallback,
                             std::vector<bool> *ambiguous_side)
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

  // The material scores are global (GlobalMin / GlobalMax), so the ambiguity of a segment's
  // process side is the same on every rank.
  if (ambiguous_side)
  {
    ambiguous_side->assign(segment_indices.size(), false);
    for (std::size_t i = 0; i < segment_indices.size(); i++)
    {
      const double tolerance =
          1.0e-10 * std::max({1.0, std::abs(minimum_score[i]), std::abs(maximum_score[i])});
      (*ambiguous_side)[i] = !(maximum_score[i] - minimum_score[i] > tolerance);
    }
  }
  std::vector<SegmentVectorContribution> normal_contributions;
  normal_contributions.reserve(candidates.size());
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
    std::array<double, 3> signed_normal{};
    for (int d = 0; d < 3; d++)
    {
      signed_normal[d] = sign * candidate.normal[d];
    }
    normal_contributions.push_back({candidate.segment, signed_normal});
  }
  const auto normal_sum = SumSegmentVectorsInCanonicalOrder(
      mesh.GetComm(), segment_indices.size(), normal_contributions);

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

  std::vector<SegmentVectorContribution> inward_contributions;
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
        std::array<double, 3> inward_unit{};
        for (int d = 0; d < 3; d++)
        {
          inward_unit[d] = inward[d] * inverse_norm;
        }
        inward_contributions.push_back({local_index, inward_unit});
      }
    }
  }
  const auto inward_sum = SumSegmentVectorsInCanonicalOrder(
      mesh.GetComm(), segment_indices.size(), inward_contributions);

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
    if (norm_squared <= 1.0e-20)
    {
      const auto &segment = geometry.segments[segment_indices[i]];
      const auto &p0 = geometry.vertices[segment.vertices[0]].coordinate;
      const auto &p1 = geometry.vertices[segment.vertices[1]].coordinate;
      std::string attributes;
      for (const int attribute : segment.metal_attributes)
      {
        attributes += (attributes.empty() ? "" : ", ") + std::to_string(attribute);
      }
      MFEM_ABORT("Unable to infer the metal-to-gap direction for the automatic edge "
                 "segment ("
                 << p0[0] << ", " << p0[1] << ", " << p0[2] << ") - (" << p1[0] << ", "
                 << p1[1] << ", " << p1[2] << ") with process normal ("
                 << process_normals[i][0] << ", " << process_normals[i][1] << ", "
                 << process_normals[i][2] << ") on metal attributes {" << attributes
                 << "} (" << inward_contributions.size() << " face contributions)!");
    }
    const double inverse_norm = 1.0 / std::sqrt(norm_squared);
    for (double &value : gap)
    {
      value *= inverse_norm;
    }
  }
  return gap_directions;
}

}  // namespace palace
