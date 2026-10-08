// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponsemirror.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <mfem.hpp>
#include "utils/communication.hpp"
#include "utils/metaledge.hpp"

namespace palace
{

namespace
{

using Point3D = std::array<double, 3>;

Point3D Sub(const Point3D &a, const Point3D &b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

double Dot(const Point3D &a, const Point3D &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double Norm(const Point3D &a)
{
  return std::sqrt(Dot(a, a));
}

Point3D Normalize(const Point3D &a)
{
  const double n = Norm(a);
  return n > 0.0 ? Point3D{a[0] / n, a[1] / n, a[2] / n} : a;
}

Point3D Cross(const Point3D &a, const Point3D &b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// Reflect through a plane subset in canonical (ascending) order.
Point3D ReflectThrough(const Point3D &p, const std::vector<MirrorPlane> &planes,
                       const std::vector<int> &subset)
{
  Point3D q = p;
  for (const int k : subset)
  {
    q = planes[k].Reflect(q);
  }
  return q;
}

Point3D ReflectDirectionThrough(const Point3D &v, const std::vector<MirrorPlane> &planes,
                                const std::vector<int> &subset)
{
  Point3D w = v;
  for (const int k : subset)
  {
    w = planes[k].ReflectDirection(w);
  }
  return w;
}

}  // namespace

std::vector<MirrorPlane> FitTruncationPlanes(const mfem::ParMesh &mesh,
                                             const std::set<int> &truncation_attributes,
                                             const std::array<double, 3> &process_normal,
                                             double matching_radius)
{
  const auto comm = mesh.GetComm();
  std::vector<MirrorPlane> planes;
  if (truncation_attributes.empty())
  {
    return planes;
  }
  MFEM_VERIFY(mesh.Dimension() == 3 && mesh.SpaceDimension() == 3,
              "Truncation planes are fitted on three-dimensional meshes only!");
  // The coordinate scale of the quantisation grid: the mesh bounding box (global).
  mfem::Vector bb_min, bb_max;
  const_cast<mfem::ParMesh &>(mesh).GetBoundingBox(bb_min, bb_max);
  double scale = 0.0;
  for (int d = 0; d < 3; d++)
  {
    scale = std::max({scale, std::abs(bb_min(d)), std::abs(bb_max(d))});
  }
  Mpi::GlobalMax(1, &scale, comm);
  scale = std::max(scale, matching_radius);
  const double normal_quantum = 1.0e-9;
  const double offset_quantum = 1.0e-9 * scale;
  // Local faces: (attribute, quantised normal, quantised offset) -> count, deviation, the
  // normal and offset sums for the averaged plane.
  struct Group
  {
    int faces = 0;
    double max_deviation = 0.0;
    Point3D normal{};
    double offset = 0.0;
    Point3D box_min = {std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity()};
    Point3D box_max = {-std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity()};
  };
  std::map<std::array<long long int, 5>, Group> local;
  mfem::Vector normal(3), center(3);
  mfem::Array<int> vertices;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int attribute = mesh.GetBdrAttribute(be);
    if (!truncation_attributes.count(attribute))
    {
      continue;
    }
    int element, info;
    mesh.GetBdrElementAdjacentElement(be, element, info);
    auto &T = *const_cast<mfem::ParMesh &>(mesh).GetBdrElementTransformation(be);
    const mfem::IntegrationPoint &ip =
        mfem::Geometries.GetCenter(mesh.GetBdrElementGeometry(be));
    T.SetIntPoint(&ip);
    mfem::CalcOrtho(T.Jacobian(), normal);
    T.Transform(ip, center);
    const double norm = normal.Norml2();
    if (norm <= 0.0)
    {
      continue;
    }
    normal /= norm;
    // Outward: away from the adjacent element's centre.
    {
      mfem::Vector element_center(3);
      const_cast<mfem::ParMesh &>(mesh).GetElementCenter(element, element_center);
      double dot = 0.0;
      for (int d = 0; d < 3; d++)
      {
        dot += (center(d) - element_center(d)) * normal(d);
      }
      if (dot < 0.0)
      {
        normal *= -1.0;
      }
    }
    const double offset = normal * center;
    double deviation = 0.0;
    mesh.GetBdrElementVertices(be, vertices);
    Point3D face_min = {std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
    Point3D face_max = {-std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity()};
    for (const int v : vertices)
    {
      const double *x = mesh.GetVertex(v);
      deviation = std::max(deviation, std::abs(normal(0) * x[0] + normal(1) * x[1] +
                                               normal(2) * x[2] - offset));
      for (int d = 0; d < 3; d++)
      {
        face_min[d] = std::min(face_min[d], x[d]);
        face_max[d] = std::max(face_max[d], x[d]);
      }
    }
    const std::array<long long int, 5> key = {
        attribute, std::llround(normal(0) / normal_quantum),
        std::llround(normal(1) / normal_quantum), std::llround(normal(2) / normal_quantum),
        std::llround(offset / offset_quantum)};
    auto &group = local[key];
    group.faces++;
    group.max_deviation = std::max(group.max_deviation, deviation);
    for (int d = 0; d < 3; d++)
    {
      group.normal[d] += normal(d);
    }
    group.offset += offset;
    for (int d = 0; d < 3; d++)
    {
      group.box_min[d] = std::min(group.box_min[d], face_min[d]);
      group.box_max[d] = std::max(group.box_max[d], face_max[d]);
    }
  }
  // Gather every rank's groups (flat records) and merge by key.
  constexpr int record_size = 17;
  std::vector<double> send;
  for (const auto &[key, group] : local)
  {
    for (const auto k : key)
    {
      send.push_back(static_cast<double>(k));
    }
    send.push_back(group.faces);
    send.push_back(group.max_deviation);
    send.push_back(group.normal[0]);
    send.push_back(group.normal[1]);
    send.push_back(group.normal[2]);
    send.push_back(group.offset);
    for (int d = 0; d < 3; d++)
    {
      send.push_back(group.box_min[d]);
    }
    for (int d = 0; d < 3; d++)
    {
      send.push_back(group.box_max[d]);
    }
  }
  const int size = Mpi::Size(comm);
  std::vector<int> counts(size), displs(size);
  const int local_count = static_cast<int>(send.size());
  Mpi::Allgather(1, &local_count, counts.data(), comm);
  int total = 0;
  for (int r = 0; r < size; r++)
  {
    displs[r] = total;
    total += counts[r];
  }
  std::vector<double> all(total);
  Mpi::Allgatherv(local_count, send.data(), all.data(), counts.data(), displs.data(), comm);
  std::map<std::array<long long int, 5>, Group> merged;
  for (int i = 0; i + record_size <= total; i += record_size)
  {
    std::array<long long int, 5> key{};
    for (int k = 0; k < 5; k++)
    {
      key[k] = static_cast<long long int>(std::llround(all[i + k]));
    }
    auto &group = merged[key];
    group.faces += static_cast<int>(std::lround(all[i + 5]));
    group.max_deviation = std::max(group.max_deviation, all[i + 6]);
    group.normal[0] += all[i + 7];
    group.normal[1] += all[i + 8];
    group.normal[2] += all[i + 9];
    group.offset += all[i + 10];
    for (int d = 0; d < 3; d++)
    {
      group.box_min[d] = std::min(group.box_min[d], all[i + 11 + d]);
      group.box_max[d] = std::max(group.box_max[d], all[i + 14 + d]);
    }
  }
  // One plane per group (the map order = (attribute, normal, offset): identical on every
  // rank), then the per-attribute status.
  for (const auto &[key, group] : merged)
  {
    // The plane constants are the QUANTISED key's (the grid both the grouping and every
    // rank see), not the gathered sums' average: identical for every rank count and
    // partition (decision 473 MINOR-9; a sum in rank order differs in the last bits).
    MirrorPlane plane;
    plane.attribute = static_cast<int>(key[0]);
    plane.normal = Normalize({static_cast<double>(key[1]) * normal_quantum,
                              static_cast<double>(key[2]) * normal_quantum,
                              static_cast<double>(key[3]) * normal_quantum});
    plane.offset = static_cast<double>(key[4]) * offset_quantum;
    plane.faces = group.faces;
    plane.max_deviation = group.max_deviation;
    plane.box_min = group.box_min;
    plane.box_max = group.box_max;
    planes.push_back(plane);
  }
  const Point3D n_process = Normalize(process_normal);
  const double joint_cap = std::cos(kArcMaxJointTurnDegrees * std::acos(-1.0) / 180.0);
  std::set<int> non_planar;
  for (std::size_t i = 0; i < planes.size(); i++)
  {
    for (std::size_t j = i + 1; j < planes.size(); j++)
    {
      if (planes[i].attribute != planes[j].attribute)
      {
        continue;
      }
      const double cosine = Dot(planes[i].normal, planes[j].normal);
      // Two non-parallel plane groups of one attribute turning by less than the arc rule's
      // joint cap: chords of a curved cut surface.
      if (cosine > joint_cap && cosine < 1.0 - 1.0e-9)
      {
        non_planar.insert(planes[i].attribute);
      }
    }
  }
  for (auto &plane : planes)
  {
    if (non_planar.count(plane.attribute))
    {
      plane.status = "NonPlanar";
    }
    else if (std::abs(Dot(plane.normal, n_process)) > 1.0e-6)
    {
      plane.status = "Unsupported";
    }
    else
    {
      plane.status = "Natural";
    }
  }
  return planes;
}

std::optional<ReflectedPoint> ReflectIntoDomain(const Point3D &p,
                                                const std::vector<MirrorPlane> &planes,
                                                double band, double tolerance)
{
  ReflectedPoint result;
  result.point = p;
  for (std::size_t k = 0; k < planes.size(); k++)
  {
    const double inside = planes[k].Inside(result.point);
    // Any point strictly beyond the plane is reflected (a sample a roundoff beyond a cut
    // face is outside the mesh for the locator; its reflection lands the same roundoff
    // inside): the tolerance bounds only the face test.
    if (inside >= 0.0 || !planes[k].NearFaces(result.point, tolerance))
    {
      continue;  // inside, or beyond the infinite plane but off its faces
    }
    if (!planes[k].Mirrors() || -inside > band)
    {
      return std::nullopt;
    }
    result.point = planes[k].Reflect(result.point);
    result.planes.push_back(static_cast<int>(k));
  }
  // A non-convex domain: the composition may still leave the point beyond a plane.
  for (const auto &plane : planes)
  {
    if (plane.Inside(result.point) < -tolerance && plane.NearFaces(result.point, tolerance))
    {
      return std::nullopt;
    }
  }
  return result;
}

MirrorExtension
ExtendIdentificationInputAcrossMirrorPlanes(const IdentificationInput &input,
                                            const std::vector<MirrorPlane> &planes,
                                            double band_over_radius)
{
  MirrorExtension extension;
  extension.input = input;
  extension.real_segments = input.segments.size();
  extension.real_vertices = input.vertices.size();
  auto &out = extension.input;
  const double R = input.radius;
  const double band = band_over_radius * R;
  const double tolerance = kSignatureParameterToleranceOverRadius * R;
  std::vector<int> mirroring;
  for (std::size_t k = 0; k < planes.size(); k++)
  {
    if (planes[k].Mirrors())
    {
      mirroring.push_back(static_cast<int>(k));
    }
  }
  if (mirroring.empty())
  {
    return extension;
  }
  auto Physical = [&](const IdentificationSegment &segment)
  { return !segment.truncation && !segment.exclusion && segment.chain >= 0; };
  // The planes a vertex lies on (within the tolerance).
  auto PlanesOf = [&](const Point3D &p)
  {
    std::vector<int> on;
    for (const int k : mirroring)
    {
      if (std::abs(planes[k].Inside(p)) <= tolerance && planes[k].NearFaces(p, tolerance))
      {
        on.push_back(k);
      }
    }
    return on;
  };
  // The subsets of the near planes of a segment (non-empty, ascending members).
  auto Subsets = [](const std::vector<int> &near)
  {
    std::vector<std::vector<int>> subsets;
    const std::size_t n = near.size();
    for (std::size_t mask = 1; mask < (std::size_t(1) << n); mask++)
    {
      std::vector<int> subset;
      for (std::size_t i = 0; i < n; i++)
      {
        if (mask & (std::size_t(1) << i))
        {
          subset.push_back(near[i]);
        }
      }
      subsets.push_back(subset);
    }
    return subsets;
  };
  // Image vertex of (real vertex, subset): the subset minus the planes the vertex lies on;
  // an empty remainder is the real vertex itself (the joint).
  std::map<std::pair<std::size_t, std::vector<int>>, std::size_t> image_vertex;
  auto ImageVertex = [&](std::size_t v, const std::vector<int> &subset)
  {
    std::vector<int> reduced;
    const auto on = PlanesOf(input.vertices[v].coordinate);
    for (const int k : subset)
    {
      if (std::find(on.begin(), on.end(), k) == on.end())
      {
        reduced.push_back(k);
      }
    }
    if (reduced.empty())
    {
      return v;
    }
    auto [it, inserted] = image_vertex.try_emplace({v, reduced}, out.vertices.size());
    if (inserted)
    {
      IdentificationVertex vertex;
      vertex.coordinate = ReflectThrough(input.vertices[v].coordinate, planes, reduced);
      vertex.physical_type = input.vertices[v].physical_type;
      vertex.on_truncation_boundary = input.vertices[v].on_truncation_boundary;
      vertex.on_port_boundary = input.vertices[v].on_port_boundary;
      out.vertices.push_back(vertex);
      extension.image_vertices++;
    }
    return it->second;
  };
  // Image chain ids per (real chain, subset); a straight joint lends the real chain's id.
  int next_chain = 0;
  for (const auto &segment : input.segments)
  {
    next_chain = std::max(next_chain, segment.chain + 1);
  }
  std::map<std::pair<int, std::vector<int>>, int> image_chain;
  // The image segments, in real segment order then subset order (deterministic).
  struct ImageRecord
  {
    std::size_t real = 0;
    std::vector<int> subset;
    std::size_t image = 0;
  };
  std::vector<ImageRecord> images;
  for (std::size_t s = 0; s < input.segments.size(); s++)
  {
    const auto &segment = input.segments[s];
    if (!Physical(segment))
    {
      continue;
    }
    std::vector<int> near;
    for (const int k : mirroring)
    {
      const double d = std::min(planes[k].Inside(segment.p0), planes[k].Inside(segment.p1));
      if (d <= band && (planes[k].NearFaces(segment.p0, tolerance) ||
                        planes[k].NearFaces(segment.p1, tolerance)))
      {
        near.push_back(k);
      }
    }
    for (const auto &subset : Subsets(near))
    {
      IdentificationSegment image = segment;
      image.p0 = ReflectThrough(segment.p0, planes, subset);
      image.p1 = ReflectThrough(segment.p1, planes, subset);
      image.vertices = {ImageVertex(segment.vertices[0], subset),
                        ImageVertex(segment.vertices[1], subset)};
      image.gap_direction = ReflectDirectionThrough(segment.gap_direction, planes, subset);
      image.process_normal = segment.process_normal;  // vertical planes: unchanged
      image.image_of = static_cast<int>(s);
      image.mirror_planes = subset;
      auto [it, inserted] = image_chain.try_emplace({segment.chain, subset}, -1);
      if (inserted)
      {
        it->second = next_chain++;
      }
      image.chain = it->second;
      images.push_back({s, subset, out.segments.size()});
      out.segments.push_back(image);
      extension.image_segments++;
    }
  }
  if (images.empty())
  {
    return extension;
  }
  // Vertex incidence of the image segments.
  for (const auto &record : images)
  {
    for (const std::size_t v : out.segments[record.image].vertices)
    {
      out.vertices[v].segments.push_back(record.image);
    }
  }
  // Band cuts: an image vertex whose real counterpart has more physical segments than were
  // imaged with its subset ends the image chain where the band ends.
  for (const auto &[key, v] : image_vertex)
  {
    const auto &real_vertex = input.vertices[key.first];
    int real_physical = 0;
    for (const std::size_t s : real_vertex.segments)
    {
      real_physical += Physical(input.segments[s]) ? 1 : 0;
    }
    int imaged = 0;
    for (const std::size_t s : out.vertices[v].segments)
    {
      imaged += s >= extension.real_segments ? 1 : 0;
    }
    if (imaged < real_physical)
    {
      out.vertices[v].image_band_cut = true;
      out.vertices[v].on_truncation_boundary = true;
      out.vertices[v].physical_type = MetalEdgeVertexType::ENDPOINT;
      extension.band_cut_vertices.push_back(v);
    }
  }
  // Joints: a real truncation vertex that now carries image segments.
  const double noise_sagitta = kJointNoiseSagittaOverRadius * R;
  auto Direction = [&](const IdentificationSegment &segment, std::size_t from)
  {
    return Normalize(segment.vertices[0] == from ? Sub(segment.p1, segment.p0)
                                                 : Sub(segment.p0, segment.p1));
  };
  // The straight piece from a vertex along its segment s: consecutive collinear segments
  // (within the direction quantum of the identification) through two-segment vertices.
  auto Piece = [&](std::size_t s, std::size_t from)
  {
    double length = 0.0;
    std::size_t current = s, v = from;
    const Point3D direction = Direction(out.segments[s], from);
    std::set<std::size_t> seen;
    while (seen.insert(current).second)
    {
      const auto &segment = out.segments[current];
      length += Norm(Sub(segment.p1, segment.p0));
      const std::size_t other =
          segment.vertices[0] == v ? segment.vertices[1] : segment.vertices[0];
      std::vector<std::size_t> next;
      for (const std::size_t t : out.vertices[other].segments)
      {
        if (t != current && Physical(out.segments[t]))
        {
          next.push_back(t);
        }
      }
      if (next.size() != 1 ||
          std::abs(Dot(Direction(out.segments[next[0]], other), direction) - 1.0) > 1.0e-12)
      {
        break;
      }
      current = next[0];
      v = other;
    }
    return length;
  };
  for (std::size_t v = 0; v < extension.real_vertices; v++)
  {
    const auto &real_vertex = input.vertices[v];
    if (!real_vertex.on_truncation_boundary)
    {
      continue;
    }
    std::vector<std::size_t> real_physical, image_physical;
    for (const std::size_t s : out.vertices[v].segments)
    {
      if (!Physical(out.segments[s]))
      {
        continue;
      }
      (s < extension.real_segments ? real_physical : image_physical).push_back(s);
    }
    if (real_physical.empty() || image_physical.empty())
    {
      continue;
    }
    auto &vertex = out.vertices[v];
    const auto on = PlanesOf(real_vertex.coordinate);
    extension.joined_vertices.push_back(v);
    extension.joined_planes.push_back(on);
    bool straight = false;
    if (real_physical.size() == 1 && image_physical.size() == 1)
    {
      const Point3D a = Direction(out.segments[real_physical[0]], v);
      const Point3D b = Direction(out.segments[image_physical[0]], v);
      // The turn between the incoming real direction (-a) and the outgoing image direction.
      const double cosine = std::clamp(-Dot(a, b), -1.0, 1.0);
      const double turn = std::acos(cosine);
      const double piece =
          std::min(Piece(real_physical[0], v), Piece(image_physical[0], v));
      straight = JointIsNoise(turn, piece, noise_sagitta);
      vertex.physical_type =
          straight ? MetalEdgeVertexType::REGULAR : MetalEdgeVertexType::CORNER;
    }
    else
    {
      vertex.physical_type = MetalEdgeVertexType::JUNCTION;
    }
    extension.joined_straight.push_back(straight);
    vertex.mirror_joint = true;
    // The vertex leaves the truncation boundary when every plane it lies on mirrors.
    bool all_mirror = true;
    for (const int k : on)
    {
      all_mirror = all_mirror && planes[k].Mirrors();
    }
    // Planes the vertex lies on that are not in the mirroring list (NonPlanar /
    // Unsupported) keep it a cut.
    for (std::size_t k = 0; k < planes.size(); k++)
    {
      if (!planes[k].Mirrors() &&
          std::abs(planes[k].Inside(real_vertex.coordinate)) <= tolerance)
      {
        all_mirror = false;
      }
    }
    if (all_mirror)
    {
      vertex.on_truncation_boundary = false;
    }
    // A straight joint: the image chain continues the real chain (one chain id).
    if (straight)
    {
      const int real_chain = out.segments[real_physical[0]].chain;
      const int image_chain_id = out.segments[image_physical[0]].chain;
      if (real_chain != image_chain_id)
      {
        for (auto &segment : out.segments)
        {
          if (segment.image_of >= 0 && segment.chain == image_chain_id)
          {
            segment.chain = real_chain;
          }
        }
      }
    }
  }
  // Image faces (CrossLayer rule) within the band.
  const std::size_t real_faces = input.faces.size();
  for (std::size_t f = 0; f < real_faces; f++)
  {
    const auto &face = input.faces[f];
    std::vector<int> near;
    for (const int k : mirroring)
    {
      double d = std::numeric_limits<double>::infinity();
      bool near_faces = false;
      for (const auto &vertex : face.vertices)
      {
        d = std::min(d, planes[k].Inside(vertex));
        near_faces = near_faces || planes[k].NearFaces(vertex, tolerance);
      }
      if (d <= band && near_faces)
      {
        near.push_back(k);
      }
    }
    for (const auto &subset : Subsets(near))
    {
      IdentificationFace image;
      for (const auto &vertex : face.vertices)
      {
        image.vertices.push_back(ReflectThrough(vertex, planes, subset));
      }
      // The reflected polygon's orientation flips; the sign-canonical normal is what the
      // CrossLayer rule reads.
      image.normal = ReflectDirectionThrough(face.normal, planes, subset);
      std::reverse(image.vertices.begin(), image.vertices.end());
      out.faces.push_back(image);
      extension.image_faces++;
    }
  }
  return extension;
}

double HalfCornerArmStart(double angle_degrees, double matching_radius)
{
  const double theta = angle_degrees * std::acos(-1.0) / 180.0;
  const double exit =
      matching_radius / std::max(std::abs(std::cos(theta)), std::abs(std::sin(theta)));
  return 0.5 * (matching_radius + exit);
}

namespace
{

using Interval = std::pair<double, double>;

// Overlap length of two intervals.
double Overlap(const Interval &a, const Interval &b)
{
  return std::max(0.0, std::min(a.second, b.second) - std::max(a.first, b.first));
}

// The mirror-formed cluster REQUIREMENT CONTRACT (decision 557 (1) / 562; impl-B5
// CONTRACT.md section 2, version 1) of an unmerged mirror-formed SpatialEdgeCluster /
// CurvedEdge configuration of the extended run: {Version, Planes, Frame, RealPortions,
// RealLengthOverR, ImageLengthOverR, RealFeatures, ExtendedFeature, Rule}. RealPortions are
// the indices of the serialised Signature.Portions whose world portions (the feature's
// IdentifiedPortions on the extended segments) are ALL real (segment < real_segments),
// mapped EXACTLY by geometry: every world portion longer than the signature tolerance is
// projected into the signature's canonical frame and must lie on exactly one serialised
// portion: both ends within the arc-fit tolerance widened by the serialisation quanta (as
// VerifySpatialEdgesInSignatureFrame reads a model Edge) on a straight segment, and on a
// mesh chord of a fitted arc (IdentifiedSegment::arc) within that bound plus the arc's
// recorded largest chord sagitta (IdentifiedArc::max_sagitta_over_R: a world portion is a
// piece of a chord whose joints lie on the fitted circle within the fit tolerance, so any
// point of it lies within the sagitta of the circle - the geometry's own bound, no
// heuristic); every serialised portion must receive world portions of ONE class (the chain
// vertex on the plane is a vertex of the extended chain, never interior to a serialised
// portion). Anything else REFUSES the contract by name (`refused`): never a length
// heuristic. The Frame is the feature's canonical frame in the identification's units (the
// manifest writers scale it). A CurvedEdge key carries no Portions: RealPortions is empty
// and the Rule says so.
nlohmann::json BuildMirrorFormedContract(
    const IdentifiedFeature &feature, const IdentificationResult &extended,
    std::size_t real_segments, double R, const std::vector<int> &planes, double real_length,
    double image_length, const nlohmann::json &real_ids, std::string &refused)
{
  refused.clear();
  const auto &segments = extended.segments;
  const double tolerance = kSignatureParameterToleranceOverRadius * R;
  // Units of R: the arc-fit tolerance widened by two signature quanta and read inclusive
  // (the library load's reading of a model Edge against the Signature's portions).
  const double straight_tolerance =
      kSignatureParameterToleranceOverRadius + 2.5 * kSignatureLengthQuantumOverRadius;
  std::vector<int> real_portions;
  std::string rule =
      "decision 557 / round-3 DESIGN 4.4 (MF): the identification's mirror-formed cluster "
      "requirement; RealPortions index Signature.Portions (0-based) whose world portions "
      "lie on real segments (segment < real_segments); every other portion is an image; "
      "Frame = the signature's canonical frame in mesh units with its Chirality (Axes[2] = "
      "Chirality x (Axes[0] x Axes[1])); lengths over R";
  if (feature.type == "CurvedEdge")
  {
    rule += "; a CurvedEdge key carries no Portions (RealPortions empty): the real / image "
            "split is RealLengthOverR / ImageLengthOverR along the one arc (the 2D "
            "curvature family's consumer, decision 557)";
  }
  else
  {
    const auto portions = feature.signature.find("Portions");
    if (portions == feature.signature.end() || !portions->is_array() || portions->empty())
    {
      refused = "NoSerialisedPortions: the " + feature.type +
                " signature carries no Portions to map the real / image split onto";
      return nullptr;
    }
    auto Local = [&](const std::array<double, 3> &p)
    {
      const Point3D r = Sub(p, feature.origin);
      return std::array<double, 3>{Dot(r, feature.axes[0]) / R, Dot(r, feature.axes[1]) / R,
                                   Dot(r, feature.axes[2]) / R};
    };
    // Per serialised portion: 0 unclaimed, 1 real, 2 image, 3 mixed.
    std::vector<int> classes(portions->size(), 0);
    for (std::size_t p = 0; p < feature.portions.size(); p++)
    {
      const auto &portion = feature.portions[p];
      if (portion.s1 - portion.s0 <= tolerance)
      {
        continue;  // a sub-tolerance piece carries no class of its own
      }
      const auto &segment = segments[portion.segment];
      const double mapping_tolerance =
          straight_tolerance +
          (segment.arc >= 0 && static_cast<std::size_t>(segment.arc) < extended.arcs.size()
               ? extended.arcs[static_cast<std::size_t>(segment.arc)].max_sagitta_over_R
               : 0.0);
      std::array<std::array<double, 2>, 2> ends{};
      for (int k = 0; k < 2; k++)
      {
        const double s = k == 0 ? portion.s0 : portion.s1;
        std::array<double, 3> world{};
        for (int d = 0; d < 3; d++)
        {
          const double direction =
              segment.length > 0.0
                  ? (segment.key[1][d] - segment.key[0][d]) / segment.length
                  : 0.0;
          world[d] = segment.key[0][d] + direction * s;
        }
        const auto local = Local(world);
        if (std::abs(local[2]) > straight_tolerance)
        {
          refused = "PortionOffPlane: world portion " + std::to_string(p) + " lies " +
                    std::to_string(std::abs(local[2])) +
                    " R off the signature's canonical plane";
          return nullptr;
        }
        ends[k] = {local[0], local[1]};
      }
      std::vector<std::size_t> candidates;
      for (std::size_t i = 0; i < portions->size(); i++)
      {
        double distance = 0.0;
        for (const auto &q : ends)
        {
          distance = std::max(distance, DistanceToSerializedPortion((*portions)[i], q));
        }
        if (distance <= mapping_tolerance)
        {
          candidates.push_back(i);
        }
      }
      if (candidates.empty())
      {
        refused = "PortionUnmapped: world portion " + std::to_string(p) + " (segment " +
                  std::to_string(portion.segment) +
                  ") lies on no serialised Signature.Portions entry within " +
                  std::to_string(mapping_tolerance) + " R";
        return nullptr;
      }
      if (candidates.size() > 1)
      {
        refused = "PortionAmbiguous: world portion " + std::to_string(p) + " (segment " +
                  std::to_string(portion.segment) + ") lies on " +
                  std::to_string(candidates.size()) +
                  " serialised Signature.Portions entries";
        return nullptr;
      }
      const int cls = portion.segment >= real_segments ? 2 : 1;
      int &entry = classes[candidates.front()];
      entry = entry == 0 ? cls : (entry == cls ? cls : 3);
    }
    for (std::size_t i = 0; i < classes.size(); i++)
    {
      if (classes[i] == 0)
      {
        refused = "PortionUnclaimed: serialised Signature.Portions[" + std::to_string(i) +
                  "] receives no world portion of the configuration";
        return nullptr;
      }
      if (classes[i] == 3)
      {
        refused = "MixedPortion: serialised Signature.Portions[" + std::to_string(i) +
                  "] receives real AND image world portions (the plane vertex is interior "
                  "to it)";
        return nullptr;
      }
      if (classes[i] == 1)
      {
        real_portions.push_back(static_cast<int>(i));
      }
    }
    if (real_portions.empty())
    {
      refused = "NoRealPortion: every serialised portion is an image";
      return nullptr;
    }
    if (real_portions.size() == classes.size())
    {
      refused = "NoImagePortion: every serialised portion is real (not mirror-formed)";
      return nullptr;
    }
  }
  // The Frame is the identification's canonical frame AS IS (Features[].Frame: Axes[2] =
  // the process normal n for both handedness values) with its Chirality made explicit
  // (CONTRACT.md v3, decision 584 (2)): Axes[2] = Chirality x (Axes[0] x Axes[1]), so a
  // chirality -1 key's left-handed triple is a valid frame, read by the consumer with the
  // same rule; the placement uses this very frame (M = identity for a Signature-keyed
  // model).
  return nlohmann::json{{"Version", 1},
                        {"Planes", planes},
                        {"Frame",
                         {{"Origin", feature.origin},
                          {"Axes", {feature.axes[0], feature.axes[1], feature.axes[2]}},
                          {"Chirality", feature.chirality}}},
                        {"RealPortions", real_portions},
                        {"RealLengthOverR", real_length / R},
                        {"ImageLengthOverR", image_length / R},
                        {"RealFeatures", real_ids},
                        {"ExtendedFeature", feature.id},
                        {"Rule", rule}};
}

}  // namespace

MirrorMergeSummary MergeMirrorIdentification(const IdentificationResult &real,
                                             const IdentificationResult &extended,
                                             const MirrorExtension &extension,
                                             const std::vector<MirrorPlane> &planes,
                                             IdentificationResult &merged)
{
  MirrorMergeSummary summary;
  merged = real;
  merged.real_segments = extension.real_segments;
  const double R = real.radius;
  const double tolerance = kSignatureParameterToleranceOverRadius * R;
  const std::size_t n_real_segments = extension.real_segments;
  const std::size_t n_real_vertices = extension.real_vertices;
  MFEM_VERIFY(real.segments.size() == n_real_segments,
              "The unextended identification lists "
                  << real.segments.size() << " segments, the input " << n_real_segments
                  << "!");
  // The image segments of the extended run follow the real ones.
  for (std::size_t s = n_real_segments; s < extended.segments.size(); s++)
  {
    merged.segments.push_back(extended.segments[s]);
    merged.segments.back().portions.clear();
  }
  auto IsImage = [&](const IdentifiedPortion &portion)
  { return portion.segment >= n_real_segments; };
  auto PlanesOfFeature = [&](const IdentifiedFeature &feature)
  {
    std::set<int> used;
    for (const auto &portion : feature.portions)
    {
      if (IsImage(portion))
      {
        const auto &segment = extension.input.segments[portion.segment];
        used.insert(segment.mirror_planes.begin(), segment.mirror_planes.end());
      }
    }
    return std::vector<int>(used.begin(), used.end());
  };
  // Real features by (segment) for the overlap search.
  std::map<std::size_t, std::vector<std::size_t>> real_features_on_segment;
  for (std::size_t i = 0; i < merged.features.size(); i++)
  {
    for (const auto &portion : merged.features[i].portions)
    {
      real_features_on_segment[portion.segment].push_back(i);
    }
  }
  struct Formed
  {
    const IdentifiedFeature *feature = nullptr;
    std::vector<int> planes;
    double real_length = 0.0, image_length = 0.0;
  };
  std::vector<Formed> formed;
  // Unmerged mirror-formed configurations emitted as merged features with their contract
  // (decision 557): appended after the formed features, numbered after them.
  struct Configuration
  {
    IdentifiedFeature feature;
    std::size_t unmerged_index = 0;
  };
  std::vector<Configuration> configurations;
  std::map<std::size_t, nlohmann::json> continued;  // real feature index -> record
  // Whether a REAL portion of a feature of the extended run is identified IDENTICALLY by
  // the unextended run (decision 512 (b) / (c), DESIGN ERRATA-7), piece by piece: the
  // portion is cut at the boundaries of the real portions covering it; a piece is identical
  // when its covering real feature has the same type, signature key and side and the turn
  // agrees within the joint noise rule (the turn difference over the piece, read as a joint
  // on a piece of that length, implies a sagitta below kJointNoiseSagittaOverRadius x R;
  // both turns are read at their portions' mean density); a part covered by no real portion
  // differs ("Uncovered"). The split of the portions does not matter (a re-split at the
  // plane, a far end clipped by a virtual corner at another plane). Pieces shorter than the
  // signature tolerance are dropped.
  const double noise_sagitta = kJointNoiseSagittaOverRadius * R;
  struct IdentityPiece
  {
    double s0 = 0.0, s1 = 0.0;
    bool identical = false;
    std::string reason;  // "" when identical, else "Uncovered" / "Type" / "Side" / "Turn"
    int real_feature = -1;  // index into merged.features, -1 when uncovered
  };
  auto IdentityOf = [&](const IdentifiedFeature &feature, const IdentifiedPortion &portion)
  {
    std::vector<IdentityPiece> pieces;
    const double length = std::max(0.0, portion.s1 - portion.s0);
    const auto it = real_features_on_segment.find(portion.segment);
    if (it != real_features_on_segment.end())
    {
      for (const std::size_t i : it->second)
      {
        const auto &candidate = merged.features[i];
        for (const auto &other : candidate.portions)
        {
          if (other.segment != portion.segment)
          {
            continue;
          }
          IdentityPiece piece;
          piece.s0 = std::max(portion.s0, other.s0);
          piece.s1 = std::min(portion.s1, other.s1);
          if (piece.s1 - piece.s0 <= tolerance)
          {
            continue;
          }
          piece.real_feature = static_cast<int>(i);
          if (candidate.type != feature.type ||
              candidate.signature_key != feature.signature_key)
          {
            piece.reason = "Type";
          }
          else if (other.side != portion.side)
          {
            piece.reason = "Side";
          }
          else
          {
            const double other_length = std::max(0.0, other.s1 - other.s0);
            const double density = length > 0.0 ? portion.turn / length : 0.0;
            const double other_density =
                other_length > 0.0 ? other.turn / other_length : 0.0;
            const double overlap = piece.s1 - piece.s0;
            const double turn_difference = std::abs(density - other_density) * overlap;
            if (!JointIsNoise(turn_difference, overlap, noise_sagitta))
            {
              piece.reason = "Turn";
            }
          }
          piece.identical = piece.reason.empty();
          pieces.push_back(piece);
        }
      }
    }
    std::sort(pieces.begin(), pieces.end(),
              [](const IdentityPiece &a, const IdentityPiece &b) { return a.s0 < b.s0; });
    // The parts of the portion no real portion covers.
    std::vector<IdentityPiece> covered = std::move(pieces);
    pieces.clear();
    double cursor = portion.s0;
    for (const auto &piece : covered)
    {
      if (piece.s0 > cursor + tolerance)
      {
        pieces.push_back({cursor, piece.s0, false, "Uncovered", -1});
      }
      pieces.push_back(piece);
      cursor = std::max(cursor, piece.s1);
    }
    if (portion.s1 > cursor + tolerance)
    {
      pieces.push_back({cursor, portion.s1, false, "Uncovered", -1});
    }
    return pieces;
  };
  for (const auto &feature : extended.features)
  {
    bool touches_image = false;
    double real_length = 0.0, image_length = 0.0;
    for (const auto &portion : feature.portions)
    {
      const double length = std::max(0.0, portion.s1 - portion.s0);
      if (IsImage(portion))
      {
        touches_image = true;
        image_length += length;
      }
      else
      {
        real_length += length;
      }
    }
    // A vertex feature at a joint vertex (its window claims lie on real and image segments)
    // touches the image through its portions; a vertex feature whose vertex is an image
    // vertex is image-only even without portions.
    bool image_vertex = false;
    for (const std::size_t v : feature.vertices)
    {
      image_vertex = image_vertex || v >= n_real_vertices;
    }
    if (!touches_image && !image_vertex)
    {
      continue;
    }
    if (real_length <= tolerance && (feature.portions.empty() ? image_vertex : true))
    {
      summary.image_only_features++;
      continue;
    }
    // The real features whose portions it overlaps, and whether every real portion is
    // identified identically by one of them (piece by piece).
    std::set<std::size_t> overlapping;
    std::vector<std::vector<IdentityPiece>> identities;
    bool all_identical = true;
    for (const auto &portion : feature.portions)
    {
      if (IsImage(portion))
      {
        identities.emplace_back();
        continue;
      }
      identities.push_back(IdentityOf(feature, portion));
      for (const auto &piece : identities.back())
      {
        if (piece.real_feature >= 0)
        {
          overlapping.insert(static_cast<std::size_t>(piece.real_feature));
        }
        all_identical = all_identical && piece.identical;
      }
    }
    // Continued (decision 512 (b)): the extended run reads every real portion of this
    // feature exactly as the unextended run reads it, on ONE real feature - whatever the
    // split of the portions and whether that real feature extends further (a far end
    // clipped by a virtual corner at another plane is clipped by the formed-feature step
    // below; a far end taken by an unmerged configuration is that configuration's own
    // cells): the real feature keeps the unextended reading verbatim.
    const bool is_continued = overlapping.size() == 1 && all_identical;
    if (is_continued)
    {
      auto &record = continued[*overlapping.begin()];
      if (record.is_null())
      {
        record = {{"Planes", PlanesOfFeature(feature)},
                  {"RealLength", real_length},
                  {"ImageLength", image_length},
                  {"Status", "Continued"}};
      }
      else
      {
        // A chain continued through two planes (a straight joint at each end).
        std::set<int> used;
        for (const auto &k : record.at("Planes"))
        {
          used.insert(k.get<int>());
        }
        for (const int k : PlanesOfFeature(feature))
        {
          used.insert(k);
        }
        record["Planes"] = std::vector<int>(used.begin(), used.end());
        record["RealLength"] = record.at("RealLength").get<double>() + real_length;
        record["ImageLength"] = record.at("ImageLength").get<double>() + image_length;
      }
      summary.continued_features++;
      continue;
    }
    // Only the topologies the placement knows how to halve take part in the merge (a
    // virtual corner, a pair / strip with its image side); a mirror-formed stack, cluster
    // or curved pair has no mirror placement (DESIGN 2.2.5): the real features it would
    // replace stay as the unextended run read them, the configuration is read as a real
    // Missing feature (decision 481: its own cells DomainBoundary with their raw claims,
    // its neighbours untouched) and the formed key is recorded for the discovery.
    const bool mergeable =
        feature.type == "ConvexCorner" || feature.type == "ConcaveCorner" ||
        feature.type == "SameConductorGap" || feature.type == "DifferentConductorGap" ||
        feature.type == "SameConductorStrip";
    if (!mergeable)
    {
      // The real features the unmerged configuration touches (a portion whose
      // identification DIFFERS overlapped, or a real vertex of the configuration) are
      // recorded with Status "Unmerged" (a record, not a block); the configuration's
      // portions, real and image, are recorded in world coordinates (mesh units), every
      // real portion flagged Identical when the unextended run identifies it exactly as the
      // extended run does, so that the operator can read it as a Missing feature (decision
      // 481): its own cells - those on DIFFERING real portions - DomainBoundary, raw kept,
      // never applied with a single-sided model; the identically identified cells keep
      // their models (decision 512 (c): the half = full / 2 identity holds for them, and an
      // unmerged chain never sends its far cells raw).
      std::set<std::size_t> touched;
      double identical_length = 0.0, differing_length = 0.0;
      for (std::size_t p = 0; p < feature.portions.size(); p++)
      {
        for (const auto &piece : identities[p])
        {
          if (piece.identical)
          {
            identical_length += piece.s1 - piece.s0;
          }
          else
          {
            differing_length += piece.s1 - piece.s0;
            if (piece.real_feature >= 0)
            {
              touched.insert(static_cast<std::size_t>(piece.real_feature));
            }
          }
        }
      }
      for (std::size_t i = 0; i < merged.features.size(); i++)
      {
        for (const std::size_t v : merged.features[i].vertices)
        {
          if (std::find(feature.vertices.begin(), feature.vertices.end(), v) !=
              feature.vertices.end())
          {
            touched.insert(i);
          }
        }
      }
      nlohmann::json real_ids = nlohmann::json::array();
      for (const std::size_t i : touched)
      {
        auto &target = merged.features[i];
        real_ids.push_back(target.id);
        summary.touched_feature_ids.push_back(target.id);
        if (target.mirror.is_null())
        {
          target.mirror = {{"Planes", PlanesOfFeature(feature)},
                           {"RealLength", target.length},
                           {"ImageLength", 0.0},
                           {"Status", "Unmerged"},
                           {"UnmergedType", feature.type}};
        }
      }
      // The world portions: image portions whole; real portions piece by piece (cut at the
      // real portions' boundaries), each flagged Identical / Differs (and the real feature
      // identifying it), so that the operator's own-cell test reads the differing pieces
      // only.
      nlohmann::json world_portions = nlohmann::json::array();
      auto WorldPortion = [&](const IdentifiedPortion &portion, double s0, double s1)
      {
        const auto &segment = extended.segments[portion.segment];
        std::array<double, 3> p0{}, p1{};
        for (int d = 0; d < 3; d++)
        {
          const double direction =
              segment.length > 0.0
                  ? (segment.key[1][d] - segment.key[0][d]) / segment.length
                  : 0.0;
          p0[d] = segment.key[0][d] + direction * s0;
          p1[d] = segment.key[0][d] + direction * s1;
        }
        return nlohmann::json{{"Segment", portion.segment}, {"S0", s0}, {"S1", s1},
                              {"Image", IsImage(portion)},  {"P0", p0}, {"P1", p1},
                              {"Length", s1 - s0}};
      };
      for (std::size_t p = 0; p < feature.portions.size(); p++)
      {
        const auto &portion = feature.portions[p];
        if (IsImage(portion))
        {
          world_portions.push_back(WorldPortion(portion, portion.s0, portion.s1));
          continue;
        }
        for (const auto &piece : identities[p])
        {
          nlohmann::json entry = WorldPortion(portion, piece.s0, piece.s1);
          entry["Identical"] = piece.identical;
          if (!piece.identical)
          {
            entry["Differs"] = piece.reason;
          }
          if (piece.real_feature >= 0)
          {
            entry["RealFeature"] =
                merged.features[static_cast<std::size_t>(piece.real_feature)].id;
          }
          world_portions.push_back(std::move(entry));
        }
      }
      const std::size_t unmerged_index = summary.unmerged_features.size();
      summary.unmerged_features.push_back({{"Feature", feature.id},
                                           {"Type", feature.type},
                                           {"Key", feature.signature_key},
                                           {"RealLength", real_length},
                                           {"ImageLength", image_length},
                                           {"IdenticalRealLength", identical_length},
                                           {"DifferingRealLength", differing_length},
                                           {"Planes", PlanesOfFeature(feature)},
                                           {"RealFeatures", real_ids},
                                           {"Portions", world_portions},
                                           {"Status", "Unmerged"}});
      // The REQUIREMENT CONTRACT (decision 557 (1) / (4), 562; impl-B5 CONTRACT.md): the
      // configuration itself becomes a merged feature of its own type with ALL its
      // portions (real and image), its real length, its real vertices and its canonical
      // frame, carrying Mirror Status "Unmerged" and the contract (or the refusal, by
      // name), so that the operator's requirement records emit it as an ordinary Missing /
      // Exact requirement WITH the Frame and the real / image split and the placement can
      // apply a matched model on its real half. The real features it touches keep the
      // unextended reading (above); the configuration is never registered in the segments'
      // portion tables (the real reading's) and is numbered after every real and formed
      // feature.
      if (feature.type == "SpatialEdgeCluster" || feature.type == "CurvedEdge")
      {
        std::string refused;
        nlohmann::json contract = BuildMirrorFormedContract(
            feature, extended, n_real_segments, R, PlanesOfFeature(feature), real_length,
            image_length, real_ids, refused);
        IdentifiedFeature copy = feature;
        copy.length = real_length;
        copy.vertices.erase(std::remove_if(copy.vertices.begin(), copy.vertices.end(),
                                           [&](std::size_t v)
                                           { return v >= n_real_vertices; }),
                            copy.vertices.end());
        copy.mirror = {{"Planes", PlanesOfFeature(feature)},
                       {"RealLength", real_length},
                       {"ImageLength", image_length},
                       {"Status", "Unmerged"},
                       {"ExtendedFeature", feature.id},
                       {"UnmergedIndex", unmerged_index},
                       {"RealFeatures", real_ids}};
        if (refused.empty())
        {
          copy.mirror["Contract"] = std::move(contract);
        }
        else
        {
          copy.mirror["ContractRefused"] = refused;
          summary.unmerged_features.back()["ContractRefused"] = refused;
        }
        configurations.push_back({std::move(copy), unmerged_index});
      }
      continue;
    }
    formed.push_back({&feature, PlanesOfFeature(feature), real_length, image_length});
  }
  for (const auto &[index, record] : continued)
  {
    auto &target = merged.features[index];
    if (target.mirror.is_null())
    {
      target.mirror = record;
    }
    else
    {
      // Continued through one plane and touched by an unmerged configuration elsewhere
      // (its differing portions): the Unmerged record stands, the continued part noted.
      target.mirror["ContinuedRealLength"] = record.at("RealLength");
      target.mirror["ContinuedPlanes"] = record.at("Planes");
    }
  }
  std::sort(summary.touched_feature_ids.begin(), summary.touched_feature_ids.end());
  summary.touched_feature_ids.erase(
      std::unique(summary.touched_feature_ids.begin(), summary.touched_feature_ids.end()),
      summary.touched_feature_ids.end());
  // The joint vertices take the extended run's reading (the virtual corner's vertex entry,
  // MirrorJoint); the feature reference is remapped to the merged Id.
  auto UpdateJointVertices = [&](const std::map<int, int> &extended_to_merged)
  {
    for (const std::size_t v : extension.joined_vertices)
    {
      const auto extended_entry =
          std::find_if(extended.vertices.begin(), extended.vertices.end(),
                       [&](const IdentifiedVertex &entry) { return entry.vertex == v; });
      auto merged_entry =
          std::find_if(merged.vertices.begin(), merged.vertices.end(),
                       [&](const IdentifiedVertex &entry) { return entry.vertex == v; });
      if (extended_entry == extended.vertices.end())
      {
        if (merged_entry != merged.vertices.end())
        {
          merged.vertices.erase(merged_entry);
        }
        continue;
      }
      IdentifiedVertex entry = *extended_entry;
      if (entry.feature >= 0)
      {
        const auto remapped = extended_to_merged.find(entry.feature);
        entry.feature = remapped != extended_to_merged.end() ? remapped->second : -1;
      }
      if (merged_entry != merged.vertices.end())
      {
        *merged_entry = entry;
      }
      else
      {
        merged.vertices.push_back(entry);
      }
    }

    std::sort(merged.vertices.begin(), merged.vertices.end(),
              [](const IdentifiedVertex &a, const IdentifiedVertex &b)
              { return a.vertex < b.vertex; });
  };
  if (formed.empty() && configurations.empty())
  {
    UpdateJointVertices({});
    return summary;
  }
  // Mirror-formed features: clip the real features' portions they overlap, reuse the Id of
  // a wholly replaced real feature, else number after every real feature.
  int next_id = 0;
  for (const auto &feature : merged.features)
  {
    next_id = std::max(next_id, feature.id + 1);
  }
  std::vector<std::size_t> emptied;
  std::vector<IdentifiedFeature> added;
  std::map<int, int> extended_to_merged;
  for (const auto &item : formed)
  {
    const auto &feature = *item.feature;
    std::set<std::size_t> touched;
    for (const auto &portion : feature.portions)
    {
      if (IsImage(portion))
      {
        continue;
      }
      const Interval claim = {portion.s0, portion.s1};
      for (std::size_t i = 0; i < merged.features.size(); i++)
      {
        auto &target = merged.features[i];
        std::vector<IdentifiedPortion> kept;
        bool changed = false;
        for (const auto &other : target.portions)
        {
          if (other.segment != portion.segment ||
              Overlap(claim, {other.s0, other.s1}) <= tolerance)
          {
            kept.push_back(other);
            continue;
          }
          changed = true;
          touched.insert(i);
          // The parts of `other` outside the claim.
          if (other.s0 < claim.first - tolerance)
          {
            auto piece = other;
            piece.s1 = std::min(other.s1, claim.first);
            kept.push_back(piece);
          }
          if (other.s1 > claim.second + tolerance)
          {
            auto piece = other;
            piece.s0 = std::max(other.s0, claim.second);
            kept.push_back(piece);
          }
        }
        if (changed)
        {
          target.portions = std::move(kept);
          target.length = 0.0;
          for (const auto &kept_portion : target.portions)
          {
            target.length += kept_portion.s1 - kept_portion.s0;
          }
        }
      }
    }
    IdentifiedFeature copy = feature;
    // The Id: a wholly replaced real feature's (the smallest), else a new one.
    int id = -1;
    for (const std::size_t i : touched)
    {
      if (merged.features[i].portions.empty() &&
          std::find(emptied.begin(), emptied.end(), i) == emptied.end())
      {
        emptied.push_back(i);
        if (id < 0 || merged.features[i].id < id)
        {
          id = merged.features[i].id;
        }
      }
    }
    copy.id = id >= 0 ? id : next_id++;
    copy.length = item.real_length;
    // Real vertices only.
    copy.vertices.erase(std::remove_if(copy.vertices.begin(), copy.vertices.end(),
                                       [&](std::size_t v) { return v >= n_real_vertices; }),
                        copy.vertices.end());
    copy.mirror = {{"Planes", item.planes},
                   {"RealLength", item.real_length},
                   {"ImageLength", item.image_length},
                   {"Status", "MirrorFormed"}};
    summary.mirror_formed_ids.push_back(copy.id);
    summary.mirror_formed_features++;
    added.push_back(std::move(copy));
  }
  // Drop the emptied real features (their Ids were lent), append the formed ones, then the
  // unmerged configurations (decision 557) numbered after them.
  std::sort(emptied.begin(), emptied.end());
  for (auto it = emptied.rbegin(); it != emptied.rend(); ++it)
  {
    merged.features.erase(merged.features.begin() + static_cast<std::ptrdiff_t>(*it));
  }
  for (auto &feature : added)
  {
    merged.features.push_back(std::move(feature));
  }
  for (std::size_t i = 0; i < formed.size(); i++)
  {
    extended_to_merged[formed[i].feature->id] = summary.mirror_formed_ids[i];
  }
  for (auto &configuration : configurations)
  {
    auto &feature = configuration.feature;
    extended_to_merged[feature.id] = next_id;
    feature.id = next_id++;
    summary.unmerged_features[configuration.unmerged_index]["MergedFeature"] = feature.id;
    merged.features.push_back(std::move(feature));
  }
  // The segments' portion lists ({s0, s1, feature}) rebuilt from the merged features: the
  // REAL reading (an unmerged configuration's portions overlap the real features' and are
  // not registered).
  for (auto &segment : merged.segments)
  {
    segment.portions.clear();
  }
  for (const auto &feature : merged.features)
  {
    if (IsUnmergedMirrorConfiguration(feature))
    {
      continue;
    }
    for (const auto &portion : feature.portions)
    {
      merged.segments[portion.segment].portions.push_back(
          {portion.s0, portion.s1, static_cast<double>(feature.id)});
    }
  }
  for (auto &segment : merged.segments)
  {
    std::sort(segment.portions.begin(), segment.portions.end());
  }
  UpdateJointVertices(extended_to_merged);
  return summary;
}

nlohmann::json DescribeMirrorBand(const std::vector<MirrorPlane> &planes,
                                  const MirrorExtension &extension,
                                  const MirrorMergeSummary &summary,
                                  const IdentificationResult &merged,
                                  double band_over_radius, double coordinate_scale)
{
  nlohmann::json plane_list = nlohmann::json::array();
  for (const auto &plane : planes)
  {
    plane_list.push_back(
        {{"Attribute", plane.attribute},
         {"Normal", plane.normal},
         {"Offset", plane.offset * coordinate_scale},
         {"Faces", plane.faces},
         {"MaxDeviation", plane.max_deviation * coordinate_scale},
         {"Box",
          {{plane.box_min[0] * coordinate_scale, plane.box_min[1] * coordinate_scale,
            plane.box_min[2] * coordinate_scale},
           {plane.box_max[0] * coordinate_scale, plane.box_max[1] * coordinate_scale,
            plane.box_max[2] * coordinate_scale}}},
         {"Status", plane.status}});
  }
  nlohmann::json formed = nlohmann::json::array();
  for (const auto &feature : merged.features)
  {
    if (feature.mirror.is_null() || feature.mirror.value("Status", "") == "Continued")
    {
      continue;
    }
    formed.push_back({{"Feature", feature.id},
                      {"Type", feature.type},
                      {"Status", feature.mirror.value("Status", "")},
                      {"Key", feature.signature_key},
                      {"Planes", feature.mirror.value("Planes", nlohmann::json::array())}});
  }
  nlohmann::json joints = nlohmann::json::array();
  for (std::size_t i = 0; i < extension.joined_vertices.size(); i++)
  {
    joints.push_back({{"Vertex", extension.joined_vertices[i]},
                      {"Planes", extension.joined_planes[i]},
                      {"Straight", extension.joined_straight[i]}});
  }
  return {{"Planes", std::move(plane_list)},
          {"BandOverR", band_over_radius},
          {"ImageSegments", extension.image_segments},
          {"ImageVertices", extension.image_vertices},
          {"ImageFaces", extension.image_faces},
          {"JoinedVertices", std::move(joints)},
          {"BandCutVertices", static_cast<int>(extension.band_cut_vertices.size())},
          {"ContinuedFeatures", summary.continued_features},
          {"ImageOnlyFeatures", summary.image_only_features},
          {"MirrorFormedFeatures", std::move(formed)},
          {"UnmergedFeatures", summary.unmerged_features},
          {"TouchedRealFeatures", summary.touched_feature_ids},
          {"Rule",
           "boundary-cut DESIGN 2.2 (decisions 442 / 454 / 481 / 512): the metal perimeter "
           "within BandOverR x R of every planar NATURAL vertical truncation plane is "
           "reflected into the identification input (image segments joined to the real "
           "chain at the truncation vertex on the plane; a straight joint - collinear "
           "within the direction quantum or a sub-noise turn - continues the chain, an "
           "oblique one is a corner of 2 theta, a parallel edge at d < R a strip / gap of "
           "2 d); the result is merged onto the unextended run so that real features keep "
           "their Ids, portions and order wherever no new feature formed (Mirror Status "
           "Continued: every real portion of the extended chain identified identically - "
           "type, key, side, turn within the joint noise rule - whatever the split of the "
           "portions and whether the real chain extends further); a mirror-formed feature "
           "is placed on its real half only (a "
           "vertex coupon on the plane with weight 1 / 2 and the real arm's cells from "
           "s_half = (R + s) / 2; a pair's real side with its side factor) with the trace "
           "taken by even extension (mirror-point evaluation); Missing / out-of-range "
           "mirror-formed features and non-mirroring planes fall back to F-DB-a; a "
           "mirror-formed stack, cluster or curved pair (no mirror placement: Unmerged) is "
           "read as a real Missing feature (decision 481) - its own cells, those on its "
           "real portions whose identification DIFFERS between the unextended and the "
           "extended run (decision 512 (c)), DomainBoundary with their raw claims; the "
           "identically identified cells and its neighbours' cells Applied / Mirrored as "
           "read"}};
}

}  // namespace palace
