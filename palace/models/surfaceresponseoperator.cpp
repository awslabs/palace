// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponseoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/interpolator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/edgedistance.hpp"
#include "utils/enum_string.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/metaledge.hpp"
#include "utils/tablecsv.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

using MatrixEntry = std::tuple<int, int, double>;
using Point2D = std::array<double, 2>;
using Point3D = std::array<double, 3>;
using ResponseCorrectionData = config::ElectrostaticSolverData::ResponseCorrectionData;
using ResponseModelData = config::ElectrostaticSolverData::ResponseCorrectionModelData;
using ResponsePatchData = config::ElectrostaticSolverData::ResponseCorrectionPatchData;

enum class LibraryTopology : char
{
  ISOLATED_EDGE,
  SAME_CONDUCTOR_GAP,
  DIFFERENT_CONDUCTOR_GAP,
  SAME_CONDUCTOR_STRIP
};

struct LibraryInterface
{
  InterfaceDielectric type;
  int coupon;
};

struct LibraryModel
{
  std::string name;
  LibraryTopology topology;
  double separation = 0.0;
  double separation_tolerance = 0.0;
  double coupon_depth = 0.0;
  ResponseModelData response;
  std::array<double, 3> reference{};
  std::vector<LibraryInterface> interfaces;
};

struct ProcessLibrary
{
  std::string name;
  double matching_radius = 0.0;
  std::vector<LibraryModel> models;
};

struct EdgeGroup2D
{
  std::vector<int> edge_attributes;
  std::map<InterfaceDielectric, int> targets;
  std::array<double, 2> process_normal{};
  double matching_radius = 0.0;
};

struct EdgeSite2D
{
  Point2D point{};
  Point2D axis_u{};
  Point2D axis_v{};
  int conductor = std::numeric_limits<int>::min();
};

struct EdgeGroup3D
{
  std::vector<std::size_t> segments;
  std::map<InterfaceDielectric, int> targets;
  std::optional<Point3D> process_normal;
  double matching_radius = 0.0;
};

struct EdgeSegment3D
{
  std::size_t geometry_index = 0;
  Point3D p0{};
  Point3D p1{};
  Point3D tangent{};
  Point3D axis_u{};
  Point3D axis_v{};
  double length = 0.0;
  int conductor = std::numeric_limits<int>::min();
};

struct EdgePair3D
{
  std::size_t first = 0;
  std::size_t second = 0;
  double first_begin = 0.0;
  double first_end = 0.0;
  double second_begin = 0.0;
  double second_end = 0.0;
};

struct AttributedSegment2D
{
  Point2D p0{};
  Point2D p1{};
  int attribute = 0;
};

struct PendingPatch
{
  std::size_t library_model = 0;
  ResponsePatchData patch;
};

struct AutomaticResponseDiagnostics
{
  double matching_radius = 0.0;
  double minimum_wave_speed = mfem::infinity();
  double selected_length = 0.0;
  double matched_length = 0.0;
  double matched_corner_neighborhood_length = 0.0;
  double maximum_curvature_ratio = 0.0;
  double maximum_library_distance = 0.0;
};

std::string ResolveLibraryPath(const std::filesystem::path &directory,
                               const std::string &path)
{
  const std::filesystem::path input(path);
  return (input.is_absolute() ? input : directory / input).lexically_normal().string();
}

LibraryTopology ParseLibraryTopology(const std::string &topology)
{
  if (topology == "IsolatedEdge")
  {
    return LibraryTopology::ISOLATED_EDGE;
  }
  if (topology == "SameConductorGap")
  {
    return LibraryTopology::SAME_CONDUCTOR_GAP;
  }
  if (topology == "DifferentConductorGap")
  {
    return LibraryTopology::DIFFERENT_CONDUCTOR_GAP;
  }
  if (topology == "SameConductorStrip")
  {
    return LibraryTopology::SAME_CONDUCTOR_STRIP;
  }
  MFEM_ABORT("Unknown fabrication-process response topology \"" << topology << "\"!");
}

ProcessLibrary ReadProcessLibrary(const std::string &path, double coordinate_scale)
{
  std::ifstream input(path);
  MFEM_VERIFY(input,
              "Unable to open fabrication-process response library \"" << path << "\"!");
  nlohmann::json data;
  input >> data;
  MFEM_VERIFY(data.value("Version", 0) == 1, "Fabrication-process response library \""
                                                 << path << "\" has unsupported version!");

  ProcessLibrary library;
  library.name = data.value("Name", std::filesystem::path(path).stem().string());
  library.matching_radius = data.at("MatchingRadius").get<double>() / coordinate_scale;
  MFEM_VERIFY(std::isfinite(library.matching_radius) && library.matching_radius > 0.0,
              "Fabrication-process response-library matching radius must be positive!");

  const auto directory = std::filesystem::path(path).parent_path();
  double default_coupon_depth = 0.0;
  if (auto depth = data.find("CouponDepth"); depth != data.end())
  {
    default_coupon_depth = depth->get<double>() / coordinate_scale;
    MFEM_VERIFY(std::isfinite(default_coupon_depth) && default_coupon_depth > 0.0,
                "Fabrication-process response-library CouponDepth must be positive!");
  }
  const auto &models = data.at("Models");
  MFEM_VERIFY(models.is_array() && !models.empty(),
              "Fabrication-process response library must contain at least one model!");
  std::set<std::string> names;
  for (const auto &entry : models)
  {
    LibraryModel model;
    model.name = entry.at("Name").get<std::string>();
    MFEM_VERIFY(!model.name.empty() && names.insert(model.name).second,
                "Fabrication-process response model names must be nonempty and unique!");
    model.topology = ParseLibraryTopology(entry.at("Topology").get<std::string>());
    model.separation = entry.value("Separation", 0.0) / coordinate_scale;
    model.separation_tolerance = entry.value("SeparationTolerance", 0.0) / coordinate_scale;
    if (auto depth = entry.find("CouponDepth"); depth != entry.end())
    {
      model.coupon_depth = depth->get<double>() / coordinate_scale;
      MFEM_VERIFY(std::isfinite(model.coupon_depth) && model.coupon_depth > 0.0,
                  "Fabrication-process response-model CouponDepth must be positive!");
    }
    else
    {
      model.coupon_depth = default_coupon_depth;
    }
    if (model.topology == LibraryTopology::ISOLATED_EDGE)
    {
      MFEM_VERIFY(model.separation == 0.0,
                  "An isolated-edge response model cannot specify a separation!");
    }
    else
    {
      MFEM_VERIFY(std::isfinite(model.separation) && model.separation > 0.0 &&
                      std::isfinite(model.separation_tolerance) &&
                      model.separation_tolerance >= 0.0,
                  "Paired-edge response models require a positive separation and a "
                  "nonnegative separation tolerance!");
    }

    model.response.fabricated_matrix =
        ResolveLibraryPath(directory, entry.at("FabricatedMatrix").get<std::string>());
    model.response.thin_matrix =
        ResolveLibraryPath(directory, entry.at("ThinMatrix").get<std::string>());
    model.response.fabricated_surface_matrix =
        entry.contains("FabricatedSurfaceMatrix")
            ? ResolveLibraryPath(directory,
                                 entry.at("FabricatedSurfaceMatrix").get<std::string>())
            : std::string{};
    model.response.thin_surface_matrix =
        entry.contains("ThinSurfaceMatrix")
            ? ResolveLibraryPath(directory,
                                 entry.at("ThinSurfaceMatrix").get<std::string>())
            : std::string{};
    MFEM_VERIFY(model.response.fabricated_surface_matrix.empty() ==
                    model.response.thin_surface_matrix.empty(),
                "Fabricated and thin surface matrices must be specified together in a "
                "fabrication-process response model!");
    model.response.basis_points =
        ResolveLibraryPath(directory, entry.at("BasisPoints").get<std::string>());
    model.reference = entry.value("Reference", std::array<double, 3>{0.0, 0.0, 0.0});
    for (double &value : model.reference)
    {
      value /= coordinate_scale;
    }

    std::set<InterfaceDielectric> types;
    if (auto interfaces = entry.find("Interfaces"); interfaces != entry.end())
    {
      for (const auto &interface : *interfaces)
      {
        InterfaceDielectric type = InterfaceDielectric::DEFAULT;
        FromString(interface.at("Type").get<std::string>(), type);
        const int coupon = interface.at("Coupon");
        MFEM_VERIFY(type != InterfaceDielectric::DEFAULT && coupon > 0 &&
                        types.insert(type).second,
                    "Fabrication-process response-model interface mappings must have "
                    "unique explicit types and positive coupon indices!");
        model.interfaces.push_back({type, coupon});
      }
    }
    MFEM_VERIFY(model.response.fabricated_surface_matrix.empty() ||
                    !model.interfaces.empty(),
                "Surface response matrices require interface mappings in the "
                "fabrication-process response model!");
    library.models.push_back(std::move(model));
  }
  return library;
}

double Dot(const Point2D &a, const Point2D &b)
{
  return a[0] * b[0] + a[1] * b[1];
}

double Dot(const Point3D &a, const Point3D &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double Norm(const Point2D &a)
{
  return std::hypot(a[0], a[1]);
}

double Norm(const Point3D &a)
{
  return std::sqrt(Dot(a, a));
}

Point2D Normalize(const Point2D &a)
{
  const double norm = Norm(a);
  MFEM_VERIFY(norm > 0.0, "Cannot normalize a zero response-correction direction!");
  return {a[0] / norm, a[1] / norm};
}

Point3D Normalize(const Point3D &a)
{
  const double norm = Norm(a);
  MFEM_VERIFY(norm > 0.0, "Cannot normalize a zero response-correction direction!");
  return {a[0] / norm, a[1] / norm, a[2] / norm};
}

Point3D Subtract(const Point3D &a, const Point3D &b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Point3D Add(const Point3D &a, const Point3D &b)
{
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Point3D Scale(double scale, const Point3D &a)
{
  return {scale * a[0], scale * a[1], scale * a[2]};
}

Point3D Interpolate(const EdgeSegment3D &segment, double distance)
{
  return Add(segment.p0, Scale(distance, segment.tangent));
}

double Distance(const Point2D &a, const Point2D &b)
{
  return std::hypot(a[0] - b[0], a[1] - b[1]);
}

double Distance(const Point3D &a, const Point3D &b)
{
  return Norm(Subtract(a, b));
}

bool PointOnSegment(const Point2D &point, const Point2D &a, const Point2D &b, double tol)
{
  const Point2D direction = {b[0] - a[0], b[1] - a[1]};
  const double length_squared = Dot(direction, direction);
  if (length_squared == 0.0)
  {
    return Distance(point, a) <= tol;
  }
  const Point2D offset = {point[0] - a[0], point[1] - a[1]};
  const double t = Dot(offset, direction) / length_squared;
  if (t < -tol || t > 1.0 + tol)
  {
    return false;
  }
  const Point2D closest = {a[0] + t * direction[0], a[1] + t * direction[1]};
  return Distance(point, closest) <= tol;
}

std::optional<int> GetConductor(const config::BoundaryData &boundaries, int attribute)
{
  if (std::find(boundaries.pec.attributes.begin(), boundaries.pec.attributes.end(),
                attribute) != boundaries.pec.attributes.end())
  {
    return 0;
  }
  for (const auto &[index, terminal] : boundaries.terminal)
  {
    if (std::find(terminal.attributes.begin(), terminal.attributes.end(), attribute) !=
        terminal.attributes.end())
    {
      return index;
    }
  }
  for (const auto &[index, potential] : boundaries.prescribed_potential)
  {
    if (std::find(potential.attributes.begin(), potential.attributes.end(), attribute) !=
        potential.attributes.end())
    {
      return index;
    }
  }
  return std::nullopt;
}

std::vector<AttributedSegment2D> GetAttributedSegments(const mfem::ParMesh &mesh,
                                                       const std::vector<int> &attributes)
{
  std::vector<AttributedSegment2D> result;
  for (const int attribute : attributes)
  {
    auto marker = mesh::BdrAttrToMarker(mesh, std::vector<int>{attribute}, true);
    for (const auto &segment : mesh::GetBoundaryElementEdgeSegments(mesh, marker))
    {
      result.push_back(
          {{segment.p0[0], segment.p0[1]}, {segment.p1[0], segment.p1[1]}, attribute});
    }
  }
  return result;
}

std::vector<EdgeSite2D> ExtractEdgeSites(const mfem::ParMesh &mesh,
                                         const config::BoundaryData &boundaries,
                                         const EdgeGroup2D &group)
{
  auto marker = mesh::BdrAttrToMarker(mesh, group.edge_attributes, true);
  const auto endpoints = mesh::GetBoundaryEdgeSegments(mesh, marker);
  const auto segments = GetAttributedSegments(mesh, group.edge_attributes);

  std::set<int> metal_attributes(boundaries.pec.attributes.begin(),
                                 boundaries.pec.attributes.end());
  metal_attributes.insert(boundaries.auxpec.attributes.begin(),
                          boundaries.auxpec.attributes.end());
  for (const auto &[index, terminal] : boundaries.terminal)
  {
    (void)index;
    metal_attributes.insert(terminal.attributes.begin(), terminal.attributes.end());
  }
  for (const auto &[index, potential] : boundaries.prescribed_potential)
  {
    (void)index;
    metal_attributes.insert(potential.attributes.begin(), potential.attributes.end());
  }
  for (const auto &conductivity : boundaries.conductivity)
  {
    metal_attributes.insert(conductivity.attributes.begin(), conductivity.attributes.end());
  }
  for (const auto &impedance : boundaries.impedance)
  {
    metal_attributes.insert(impedance.attributes.begin(), impedance.attributes.end());
  }
  for (const auto &impedance : boundaries.rational_impedance)
  {
    metal_attributes.insert(impedance.attributes.begin(), impedance.attributes.end());
  }
  const std::vector<int> metal_attribute_list(metal_attributes.begin(),
                                              metal_attributes.end());
  const auto metal_marker = mesh::BdrAttrToMarker(mesh, metal_attribute_list, true);
  const auto physical_metal_endpoints = mesh::GetBoundaryEdgeSegments(mesh, metal_marker);

  std::set<int> nontruncation_attributes = metal_attributes;
  nontruncation_attributes.insert(boundaries.cracked_attributes.begin(),
                                  boundaries.cracked_attributes.end());
  mfem::Array<int> truncation_marker(mesh.bdr_attributes.Max());
  truncation_marker = 0;
  for (int i = 0; i < mesh.bdr_attributes.Size(); i++)
  {
    const int attribute = mesh.bdr_attributes[i];
    if (nontruncation_attributes.find(attribute) == nontruncation_attributes.end())
    {
      truncation_marker[attribute - 1] = 1;
    }
  }
  const auto exterior_segments =
      mesh::GetBoundaryElementEdgeSegments(mesh, truncation_marker, true);

  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
  double extent = 0.0;
  for (int d = 0; d < 2; d++)
  {
    extent = std::max(extent, bbmax[d] - bbmin[d]);
  }
  const double tolerance = 1.0e-9 * std::max(extent, group.matching_radius);

  std::vector<EdgeSite2D> sites;
  for (const auto &endpoint : endpoints)
  {
    const Point2D point = {endpoint.p0[0], endpoint.p0[1]};
    const bool physical_metal_edge = std::any_of(
        physical_metal_endpoints.begin(), physical_metal_endpoints.end(),
        [&](const auto &metal_endpoint)
        {
          return Distance(point, {metal_endpoint.p0[0], metal_endpoint.p0[1]}) <= tolerance;
        });
    if (!physical_metal_edge)
    {
      continue;
    }
    const bool truncated =
        std::any_of(exterior_segments.begin(), exterior_segments.end(),
                    [&](const auto &segment)
                    {
                      return PointOnSegment(point, {segment.p0[0], segment.p0[1]},
                                            {segment.p1[0], segment.p1[1]}, tolerance);
                    });
    if (truncated)
    {
      continue;
    }

    Point2D inward = {};
    std::set<int> conductors;
    for (const auto &segment : segments)
    {
      Point2D direction;
      if (Distance(point, segment.p0) <= tolerance)
      {
        direction = {segment.p1[0] - point[0], segment.p1[1] - point[1]};
      }
      else if (Distance(point, segment.p1) <= tolerance)
      {
        direction = {segment.p0[0] - point[0], segment.p0[1] - point[1]};
      }
      else
      {
        continue;
      }
      if (Norm(direction) > tolerance)
      {
        direction = Normalize(direction);
        inward[0] += direction[0];
        inward[1] += direction[1];
      }
      if (auto conductor = GetConductor(boundaries, segment.attribute))
      {
        conductors.insert(*conductor);
      }
    }
    MFEM_VERIFY(Norm(inward) > 0.0,
                "Unable to infer the in-plane direction of an automatically detected "
                "two-dimensional metal edge!");
    MFEM_VERIFY(conductors.size() == 1,
                "Unable to assign an automatically detected two-dimensional metal edge "
                "to exactly one electrostatic conductor!");

    EdgeSite2D site;
    site.point = point;
    inward = Normalize(inward);
    site.axis_u = {-inward[0], -inward[1]};
    Point2D normal = group.process_normal;
    normal[0] -= Dot(normal, site.axis_u) * site.axis_u[0];
    normal[1] -= Dot(normal, site.axis_u) * site.axis_u[1];
    MFEM_VERIFY(Norm(normal) > 1.0e-8,
                "Interface EdgeFrameNormal is parallel to a detected metal edge!");
    site.axis_v = Normalize(normal);
    site.conductor = *conductors.begin();
    sites.push_back(site);
  }
  return sites;
}

std::optional<std::size_t> FindLibraryModel(const ProcessLibrary &library,
                                            LibraryTopology topology, double separation)
{
  std::optional<std::size_t> best;
  double best_error = mfem::infinity();
  for (std::size_t i = 0; i < library.models.size(); i++)
  {
    const auto &model = library.models[i];
    if (model.topology != topology)
    {
      continue;
    }
    const double error = std::abs(model.separation - separation);
    const double tolerance =
        std::max(model.separation_tolerance,
                 1.0e-10 * std::max(library.matching_radius, separation));
    if (error <= tolerance && error < best_error)
    {
      best = i;
      best_error = error;
    }
  }
  return best;
}

std::string TopologyName(LibraryTopology topology)
{
  switch (topology)
  {
    case LibraryTopology::ISOLATED_EDGE:
      return "isolated edge";
    case LibraryTopology::SAME_CONDUCTOR_GAP:
      return "same-conductor gap";
    case LibraryTopology::DIFFERENT_CONDUCTOR_GAP:
      return "different-conductor gap";
    case LibraryTopology::SAME_CONDUCTOR_STRIP:
      return "same-conductor strip";
  }
  return "unknown";
}

ResponseCorrectionData BuildAutomaticResponseData2D(const IoData &iodata,
                                                    const LaplaceOperator &laplace_op,
                                                    const ResponseCorrectionData &request)
{
  const auto &mesh = laplace_op.GetH1Space().GetParMesh();
  MFEM_VERIFY(mesh.Dimension() == 2 && mesh.SpaceDimension() == 2,
              "Automatic fabrication-process response matching currently supports only "
              "two-dimensional electrostatic meshes!");
  const double coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
  const auto library = ReadProcessLibrary(request.library, coordinate_scale);

  std::set<int> target_filter(request.target_interfaces.begin(),
                              request.target_interfaces.end());
  std::map<std::vector<int>, EdgeGroup2D> groups_by_attributes;
  for (const auto &[index, dielectric] : iodata.boundaries.postpro.dielectric)
  {
    if ((!target_filter.empty() && target_filter.find(index) == target_filter.end()) ||
        dielectric.type == InterfaceDielectric::DEFAULT ||
        dielectric.edge_distances.empty())
    {
      continue;
    }
    MFEM_VERIFY(!dielectric.automatic_edges && !dielectric.edge_attributes.empty(),
                "Automatic two-dimensional response matching requires EdgeAttributes on "
                "every target dielectric interface!");
    const double radius = dielectric.edge_distances.back();
    MFEM_VERIFY(std::abs(radius - library.matching_radius) <=
                    1.0e-10 * std::max(radius, library.matching_radius),
                "The largest EdgeDistances value for target interface "
                    << index << " does not match the fabrication-process library radius!");

    auto &group = groups_by_attributes[dielectric.edge_attributes];
    if (group.edge_attributes.empty())
    {
      group.edge_attributes = dielectric.edge_attributes;
      group.matching_radius = radius;
    }
    MFEM_VERIFY(group.targets.emplace(dielectric.type, index).second,
                "Automatic response matching found multiple target interfaces of type "
                    << ToString(dielectric.type) << " with the same EdgeAttributes!");
    if (dielectric.edge_frame_normal)
    {
      const Point2D normal = {(*dielectric.edge_frame_normal)[0],
                              (*dielectric.edge_frame_normal)[1]};
      MFEM_VERIFY(Norm(normal) > 0.0,
                  "Two-dimensional response matching requires an in-plane "
                  "EdgeFrameNormal!");
      const auto normalized = Normalize(normal);
      if (Norm(group.process_normal) == 0.0)
      {
        group.process_normal = normalized;
      }
      else
      {
        MFEM_VERIFY(Dot(group.process_normal, normalized) > 1.0 - 1.0e-10,
                    "Target interfaces sharing EdgeAttributes must use the same "
                    "EdgeFrameNormal!");
      }
    }
  }
  MFEM_VERIFY(!groups_by_attributes.empty(),
              "Fabrication-process response matching found no target interfaces!");
  if (!target_filter.empty())
  {
    std::set<int> found;
    for (const auto &[attributes, group] : groups_by_attributes)
    {
      (void)attributes;
      for (const auto &[type, index] : group.targets)
      {
        (void)type;
        found.insert(index);
      }
    }
    MFEM_VERIFY(found == target_filter,
                "One or more response-correction TargetInterfaces is missing, untyped, "
                "or does not configure edge-distance postprocessing!");
  }

  ResponseCorrectionData result;
  result.unmatched_policy = request.unmatched_policy;
  int next_model_index = 1;
  int matched_clusters = 0;
  int matched_edges = 0;
  int unmatched_clusters = 0;
  for (auto &[attributes, group] : groups_by_attributes)
  {
    (void)attributes;
    if (Norm(group.process_normal) == 0.0)
    {
      const std::string reason =
          "target interfaces sharing EdgeAttributes require EdgeFrameNormal for "
          "automatic two-dimensional response matching";
      if (request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
      {
        MFEM_ABORT(reason);
      }
      Mpi::Warning("{}; correction is disabled for this interface group!\n", reason);
      unmatched_clusters++;
      continue;
    }

    const auto sites = ExtractEdgeSites(mesh, iodata.boundaries, group);
    MFEM_VERIFY(!sites.empty(),
                "Automatic response matching found no physical metal edges for target "
                "interface group!");

    std::vector<int> component(sites.size(), -1);
    int component_count = 0;
    for (std::size_t seed = 0; seed < sites.size(); seed++)
    {
      if (component[seed] >= 0)
      {
        continue;
      }
      std::vector<std::size_t> queue = {seed};
      component[seed] = component_count;
      for (std::size_t cursor = 0; cursor < queue.size(); cursor++)
      {
        const std::size_t current = queue[cursor];
        for (std::size_t neighbor = 0; neighbor < sites.size(); neighbor++)
        {
          if (component[neighbor] < 0 &&
              Distance(sites[current].point, sites[neighbor].point) <
                  2.0 * group.matching_radius *
                      (1.0 - 16.0 * std::numeric_limits<double>::epsilon()))
          {
            component[neighbor] = component_count;
            queue.push_back(neighbor);
          }
        }
      }
      component_count++;
    }

    std::vector<PendingPatch> pending;
    bool group_matched = true;
    for (int component_index = 0; component_index < component_count; component_index++)
    {
      std::vector<std::size_t> cluster;
      for (std::size_t i = 0; i < component.size(); i++)
      {
        if (component[i] == component_index)
        {
          cluster.push_back(i);
        }
      }

      LibraryTopology topology = LibraryTopology::ISOLATED_EDGE;
      double separation = 0.0;
      ResponsePatchData patch;
      if (cluster.size() == 1)
      {
        const auto &edge = sites[cluster.front()];
        patch.origin = {edge.point[0], edge.point[1], 0.0};
        patch.axis_u = {edge.axis_u[0], edge.axis_u[1], 0.0};
        patch.axis_v = {edge.axis_v[0], edge.axis_v[1], 0.0};
      }
      else if (cluster.size() == 2)
      {
        const auto &first = sites[cluster[0]];
        const auto &second = sites[cluster[1]];
        separation = Distance(first.point, second.point);
        Point2D direction = Normalize(
            Point2D{second.point[0] - first.point[0], second.point[1] - first.point[1]});
        const bool facing =
            Dot(first.axis_u, direction) > 0.95 && Dot(second.axis_u, direction) < -0.95;
        const bool outward =
            Dot(first.axis_u, direction) < -0.95 && Dot(second.axis_u, direction) > 0.95;
        const bool same_conductor = first.conductor == second.conductor;
        if (facing)
        {
          topology = same_conductor ? LibraryTopology::SAME_CONDUCTOR_GAP
                                    : LibraryTopology::DIFFERENT_CONDUCTOR_GAP;
        }
        else if (outward && same_conductor)
        {
          topology = LibraryTopology::SAME_CONDUCTOR_STRIP;
        }
        else
        {
          group_matched = false;
          Mpi::Warning(
              "No canonical paired-edge topology for two edges separated by {:.6e} mesh "
              "units; correction is disabled for this interface group!\n",
              separation * coordinate_scale);
          unmatched_clusters++;
          break;
        }
        MFEM_VERIFY(Dot(first.axis_v, second.axis_v) > 0.95,
                    "Nearby edges with opposing process normals require a dedicated "
                    "cross-layer response model!");
        patch.origin = {0.5 * (first.point[0] + second.point[0]),
                        0.5 * (first.point[1] + second.point[1]), 0.0};
        patch.axis_u = {direction[0], direction[1], 0.0};
        const Point2D normal = Normalize(Point2D{first.axis_v[0] + second.axis_v[0],
                                                 first.axis_v[1] + second.axis_v[1]});
        patch.axis_v = {normal[0], normal[1], 0.0};
      }
      else
      {
        group_matched = false;
        Mpi::Warning(
            "No fabrication-process response model for a cluster of {} nearby metal "
            "edges; correction is disabled for this interface group!\n",
            cluster.size());
        unmatched_clusters++;
        break;
      }

      const auto model_index = FindLibraryModel(library, topology, separation);
      if (!model_index)
      {
        group_matched = false;
        Mpi::Warning(
            "Fabrication-process response library \"{}\" has no {} model at separation "
            "{:.6e} mesh units; correction is disabled for this interface group!\n",
            library.name, TopologyName(topology), separation * coordinate_scale);
        unmatched_clusters++;
        break;
      }
      patch.reference = library.models[*model_index].reference;
      if (library.models[*model_index].coupon_depth > 0.0)
      {
        patch.weight = 1.0 / library.models[*model_index].coupon_depth;
      }
      pending.push_back({*model_index, patch});
    }

    if (!group_matched)
    {
      if (request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
      {
        MFEM_ABORT("Automatic fabrication-process response matching failed!");
      }
      continue;
    }

    std::map<std::size_t, int> runtime_models;
    for (auto &selection : pending)
    {
      auto [model_it, inserted] =
          runtime_models.emplace(selection.library_model, next_model_index);
      if (inserted)
      {
        const auto &source = library.models[selection.library_model];
        auto model = source.response;
        model.idx = next_model_index++;
        for (const auto &[type, target] : group.targets)
        {
          const auto interface_type = type;
          auto interface = std::find_if(source.interfaces.begin(), source.interfaces.end(),
                                        [interface_type](const auto &entry)
                                        { return entry.type == interface_type; });
          MFEM_VERIFY(interface != source.interfaces.end() ||
                          source.response.fabricated_surface_matrix.empty(),
                      "Fabrication-process response model \"" << source.name << "\" has no "
                                                              << ToString(type)
                                                              << " interface mapping!");
          if (interface != source.interfaces.end())
          {
            model.interfaces.push_back({target, interface->coupon});
          }
        }
        result.models.push_back(std::move(model));
      }
      selection.patch.model = model_it->second;
      result.patches.push_back(selection.patch);
      matched_clusters++;
    }
    matched_edges += static_cast<int>(sites.size());
  }

  MFEM_VERIFY(!result.models.empty() && !result.patches.empty(),
              "Fabrication-process response matching produced no usable correction "
              "patches!");
  Mpi::Print("\nAutomatic fabrication-process response matching:\n"
             " Library: {}\n"
             " Matched edge sites: {:d}\n"
             " Matched clusters: {:d}\n"
             " Unmatched interface groups: {:d}\n",
             library.name, matched_edges, matched_clusters, unmatched_clusters);
  return result;
}

bool SegmentsShareVertex(const MetalEdgeSegment &a, const MetalEdgeSegment &b)
{
  return a.vertices[0] == b.vertices[0] || a.vertices[0] == b.vertices[1] ||
         a.vertices[1] == b.vertices[0] || a.vertices[1] == b.vertices[1];
}

double SegmentDistanceSquared(const Point3D &p0, const Point3D &p1, const Point3D &q0,
                              const Point3D &q1)
{
  const Point3D u = Subtract(p1, p0);
  const Point3D v = Subtract(q1, q0);
  const Point3D w = Subtract(p0, q0);
  const double a = Dot(u, u);
  const double b = Dot(u, v);
  const double c = Dot(v, v);
  const double d = Dot(u, w);
  const double e = Dot(v, w);
  const double denominator = a * c - b * b;
  MFEM_VERIFY(a > 0.0 && c > 0.0,
              "Cannot measure distance to a zero-length metal edge segment!");

  double s_numerator, s_denominator = denominator;
  double t_numerator, t_denominator = denominator;
  if (denominator <= 1.0e-14 * a * c)
  {
    s_numerator = 0.0;
    s_denominator = 1.0;
    t_numerator = e;
    t_denominator = c;
  }
  else
  {
    s_numerator = b * e - c * d;
    t_numerator = a * e - b * d;
    if (s_numerator < 0.0)
    {
      s_numerator = 0.0;
      t_numerator = e;
      t_denominator = c;
    }
    else if (s_numerator > s_denominator)
    {
      s_numerator = s_denominator;
      t_numerator = e + b;
      t_denominator = c;
    }
  }
  if (t_numerator < 0.0)
  {
    t_numerator = 0.0;
    if (-d < 0.0)
    {
      s_numerator = 0.0;
    }
    else if (-d > a)
    {
      s_numerator = s_denominator;
    }
    else
    {
      s_numerator = -d;
      s_denominator = a;
    }
  }
  else if (t_numerator > t_denominator)
  {
    t_numerator = t_denominator;
    if (-d + b < 0.0)
    {
      s_numerator = 0.0;
    }
    else if (-d + b > a)
    {
      s_numerator = s_denominator;
    }
    else
    {
      s_numerator = -d + b;
      s_denominator = a;
    }
  }
  const double s = std::abs(s_numerator) <= 1.0e-30 ? 0.0 : s_numerator / s_denominator;
  const double t = std::abs(t_numerator) <= 1.0e-30 ? 0.0 : t_numerator / t_denominator;
  const Point3D delta = Add(w, Subtract(Scale(s, u), Scale(t, v)));
  return Dot(delta, delta);
}

ResponseCorrectionData BuildAutomaticResponseData3D(
    const IoData &iodata, const mfem::ParMesh &mesh, const MaterialOperator &mat_op,
    const ResponseCorrectionData &request, bool maxwell,
    AutomaticResponseDiagnostics *diagnostics = nullptr)
{
  MFEM_VERIFY(mesh.Dimension() == 3 && mesh.SpaceDimension() == 3,
              "Automatic three-dimensional fabrication-process response matching "
              "requires a three-dimensional mesh!");
  const double coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
  const auto library = ReadProcessLibrary(request.library, coordinate_scale);
  if (diagnostics)
  {
    diagnostics->matching_radius = library.matching_radius;
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      const int attribute = mesh.GetAttribute(element);
      diagnostics->minimum_wave_speed =
          std::min(diagnostics->minimum_wave_speed,
                   mat_op.GetLightSpeedMin(attribute));
    }
    Mpi::GlobalMin(1, &diagnostics->minimum_wave_speed, mesh.GetComm());
    MFEM_VERIFY(std::isfinite(diagnostics->minimum_wave_speed) &&
                    diagnostics->minimum_wave_speed > 0.0,
                "Unable to determine a positive wave speed for Maxwell surface-response "
                "confidence diagnostics!");
  }
  const auto geometry = ExtractMetalEdgeGeometry(mesh, iodata.boundaries);
  MFEM_VERIFY(!geometry.Empty(),
              "Fabrication-process response matching found no metal perimeter!");

  std::set<int> target_filter(request.target_interfaces.begin(),
                              request.target_interfaces.end());
  std::map<std::vector<std::size_t>, EdgeGroup3D> groups_by_segments;
  for (const auto &[index, dielectric] : iodata.boundaries.postpro.dielectric)
  {
    if ((!target_filter.empty() && target_filter.find(index) == target_filter.end()) ||
        dielectric.type == InterfaceDielectric::DEFAULT ||
        dielectric.edge_distances.empty())
    {
      continue;
    }
    MFEM_VERIFY(dielectric.automatic_edges,
                "Automatic three-dimensional response matching requires "
                "AutomaticEdges on every target dielectric interface!");
    const double radius = dielectric.edge_distances.back();
    MFEM_VERIFY(std::abs(radius - library.matching_radius) <=
                    1.0e-10 * std::max(radius, library.matching_radius),
                "The largest EdgeDistances value for target interface "
                    << index << " does not match the fabrication-process library radius!");
    auto segment_indices =
        GetInterfaceMetalEdgeSegmentIndices(geometry, index, dielectric.type);
    auto &group = groups_by_segments[segment_indices];
    if (group.segments.empty())
    {
      group.segments = std::move(segment_indices);
      group.matching_radius = radius;
    }
    MFEM_VERIFY(group.targets.emplace(dielectric.type, index).second,
                "Automatic response matching found multiple target interfaces of type "
                    << ToString(dielectric.type)
                    << " on the same three-dimensional metal perimeter!");
    if (dielectric.edge_frame_normal)
    {
      const Point3D normal = Normalize(*dielectric.edge_frame_normal);
      if (!group.process_normal)
      {
        group.process_normal = normal;
      }
      else
      {
        MFEM_VERIFY(Dot(*group.process_normal, normal) > 1.0 - 1.0e-10,
                    "Target interfaces on the same metal perimeter must use the same "
                    "EdgeFrameNormal!");
      }
    }
  }
  MFEM_VERIFY(!groups_by_segments.empty(),
              "Fabrication-process response matching found no target interfaces!");
  if (!target_filter.empty())
  {
    std::set<int> found;
    for (const auto &[segments, group] : groups_by_segments)
    {
      (void)segments;
      for (const auto &[type, index] : group.targets)
      {
        (void)type;
        found.insert(index);
      }
    }
    MFEM_VERIFY(found == target_filter,
                "One or more response-correction TargetInterfaces is missing, untyped, "
                "or does not configure edge-distance postprocessing!");
  }

  std::set<std::size_t> owned_segments;
  for (const auto &[segments, group] : groups_by_segments)
  {
    (void)group;
    for (const std::size_t segment : segments)
    {
      MFEM_VERIFY(owned_segments.insert(segment).second,
                  "Three-dimensional response-correction interface groups have "
                  "partially overlapping metal perimeters. Split or merge the target "
                  "interface definitions so each physical segment has one complete "
                  "SA/MS/MA mapping!");
    }
  }

  const auto &quadrature =
      mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * std::max(1, iodata.solver.order));
  ResponseCorrectionData result;
  result.unmatched_policy = request.unmatched_policy;
  int next_model_index = 1;
  int matched_intervals = 0;
  int matched_segments = 0;
  int unmatched_groups = 0;
  int nonregular_vertices = 0;
  std::optional<Point3D> layer_normal;
  double layer_offset = 0.0;
  mfem::Vector mesh_bbmin, mesh_bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, mesh_bbmin, mesh_bbmax);
  double mesh_extent = 0.0;
  for (int d = 0; d < 3; d++)
  {
    mesh_extent = std::max(mesh_extent, mesh_bbmax[d] - mesh_bbmin[d]);
  }
  for (const auto &segment_group : groups_by_segments)
  {
    const auto &segment_key = segment_group.first;
    const auto &group = segment_group.second;
    (void)segment_key;
    int group_matched_intervals = 0;
    double group_selected_length = 0.0;
    double group_corner_neighborhood_length = 0.0;
    double group_maximum_curvature_ratio = 0.0;
    double group_maximum_library_distance = 0.0;
    const auto process_normals = BuildMetalEdgeProcessNormals(
        mesh, geometry, group.segments,
        [&](int attribute) { return mat_op.GetLightSpeedMax(attribute); },
        group.process_normal);
    const auto gap_directions =
        BuildMetalEdgeGapDirections(mesh, geometry, group.segments, process_normals);

    std::vector<EdgeSegment3D> segments;
    segments.reserve(group.segments.size());
    std::map<std::size_t, std::size_t> local_indices;
    for (std::size_t i = 0; i < group.segments.size(); i++)
    {
      const std::size_t geometry_index = group.segments[i];
      const auto &source = geometry.segments[geometry_index];
      EdgeSegment3D segment;
      segment.geometry_index = geometry_index;
      segment.p0 = geometry.vertices[source.vertices[0]].coordinate;
      segment.p1 = geometry.vertices[source.vertices[1]].coordinate;
      segment.length = Distance(segment.p0, segment.p1);
      segment.tangent = Normalize(Subtract(segment.p1, segment.p0));
      segment.axis_u = gap_directions[i];
      segment.axis_v = process_normals[i];
      if (maxwell)
      {
        MFEM_VERIFY(
            !source.conditions.empty() &&
                std::all_of(source.conditions.begin(), source.conditions.end(),
                            [](const auto &condition)
                            {
                              return condition.type == MetalBoundaryConditionType::PEC;
                            }),
            "Maxwell surface-response correction currently supports PEC target metal "
            "boundaries only!");
        MFEM_VERIFY(source.component >= 0,
                    "Unable to determine connected PEC ownership for a Maxwell edge!");
        segment.conductor = source.component;
      }
      else
      {
        std::set<int> conductors;
        for (const int attribute : source.metal_attributes)
        {
          if (auto conductor = GetConductor(iodata.boundaries, attribute))
          {
            conductors.insert(*conductor);
          }
        }
        MFEM_VERIFY(conductors.size() == 1,
                    "Unable to assign an automatically detected three-dimensional metal "
                    "edge to exactly one electrostatic conductor!");
      segment.conductor = *conductors.begin();
      }
      local_indices.emplace(geometry_index, segments.size());
      segments.push_back(segment);
      group_selected_length += segment.length;
    }

    if (diagnostics)
    {
      diagnostics->selected_length += group_selected_length;
      const double layer_tolerance =
          1.0e-8 * std::max(mesh_extent, group.matching_radius);
      for (std::size_t i = 0; i < segments.size(); i++)
      {
        const auto &segment = segments[i];
        if (!layer_normal)
        {
          layer_normal = segment.axis_v;
          layer_offset = Dot(*layer_normal, segment.p0);
        }
        MFEM_VERIFY(Dot(*layer_normal, segment.axis_v) > 1.0 - 1.0e-6,
                    "Maxwell surface-response correction currently requires one planar "
                    "fabrication layer with a common process normal!");
        MFEM_VERIFY(
            std::abs(Dot(*layer_normal, segment.p0) - layer_offset) <= layer_tolerance &&
                std::abs(Dot(*layer_normal, segment.p1) - layer_offset) <= layer_tolerance,
            "Maxwell surface-response correction currently requires all target PEC edges "
            "to lie in one fabrication plane!");

        const auto &source = geometry.segments[segment.geometry_index];
        double corner_length = 0.0;
        for (const std::size_t vertex : source.vertices)
        {
          if (geometry.vertices[vertex].physical_type &&
              *geometry.vertices[vertex].physical_type != MetalEdgeVertexType::REGULAR)
          {
            corner_length += std::min(group.matching_radius, segment.length);
          }
        }
        group_corner_neighborhood_length += std::min(segment.length, corner_length);
      }

      std::set<std::size_t> group_vertices;
      for (const auto &segment : segments)
      {
        const auto &source = geometry.segments[segment.geometry_index];
        group_vertices.insert(source.vertices.begin(), source.vertices.end());
      }
      for (const std::size_t vertex : group_vertices)
      {
        if (geometry.vertices[vertex].physical_type != MetalEdgeVertexType::REGULAR)
        {
          continue;
        }
        std::vector<const EdgeSegment3D *> incident;
        for (const std::size_t geometry_segment : geometry.vertices[vertex].segments)
        {
          auto local = local_indices.find(geometry_segment);
          if (local != local_indices.end())
          {
            incident.push_back(&segments[local->second]);
          }
        }
        if (incident.size() != 2)
        {
          continue;
        }
        auto DirectionAway = [&](const EdgeSegment3D &segment)
        {
          const auto &source = geometry.segments[segment.geometry_index];
          return source.vertices[0] == vertex ? segment.tangent
                                              : Scale(-1.0, segment.tangent);
        };
        const Point3D first = DirectionAway(*incident[0]);
        const Point3D second = DirectionAway(*incident[1]);
        const double turning_angle =
            std::acos(std::clamp(-Dot(first, second), -1.0, 1.0));
        const double local_length =
            0.5 * (incident[0]->length + incident[1]->length);
        if (local_length > 0.0)
        {
          group_maximum_curvature_ratio =
              std::max(group_maximum_curvature_ratio,
                       group.matching_radius * turning_angle / local_length);
        }
      }
    }

    std::set<std::size_t> vertices;
    for (const auto &segment : segments)
    {
      const auto &source = geometry.segments[segment.geometry_index];
      vertices.insert(source.vertices.begin(), source.vertices.end());
    }
    for (const std::size_t vertex : vertices)
    {
      const auto type = geometry.vertices[vertex].physical_type;
      nonregular_vertices += type && *type != MetalEdgeVertexType::REGULAR;
    }

    bool group_matched = true;
    const double interaction_distance =
        2.0 * group.matching_radius * (1.0 - 16.0 * std::numeric_limits<double>::epsilon());
    const double interaction_distance_squared = interaction_distance * interaction_distance;

    // A nearby edge outside this interface group would require patch-specific interface
    // mappings. Fail closed until such cross-group coupled coupons are represented.
    for (const auto &segment : segments)
    {
      const auto &source = geometry.segments[segment.geometry_index];
      for (std::size_t other_index = 0; other_index < geometry.segments.size();
           other_index++)
      {
        const auto &other = geometry.segments[other_index];
        if (other.type != MetalEdgeSegmentType::PHYSICAL ||
            local_indices.find(other_index) != local_indices.end() ||
            other.physical_chain == source.physical_chain ||
            SegmentsShareVertex(source, other))
        {
          continue;
        }
        const auto &q0 = geometry.vertices[other.vertices[0]].coordinate;
        const auto &q1 = geometry.vertices[other.vertices[1]].coordinate;
        if (SegmentDistanceSquared(segment.p0, segment.p1, q0, q1) <
            interaction_distance_squared)
        {
          Mpi::Warning(
              "A three-dimensional target edge is within 2R of a physical metal edge "
              "outside its interface group; correction is disabled for this group!\n");
          group_matched = false;
          break;
        }
      }
      if (!group_matched)
      {
        break;
      }
    }

    std::vector<EdgePair3D> pairs;
    std::vector<std::vector<std::pair<double, double>>> paired_intervals(segments.size());
    for (std::size_t i = 0; group_matched && i < segments.size(); i++)
    {
      const auto &first_source = geometry.segments[segments[i].geometry_index];
      for (std::size_t j = i + 1; j < segments.size(); j++)
      {
        const auto &second_source = geometry.segments[segments[j].geometry_index];
        if (first_source.physical_chain == second_source.physical_chain ||
            SegmentsShareVertex(first_source, second_source) ||
            SegmentDistanceSquared(segments[i].p0, segments[i].p1, segments[j].p0,
                                   segments[j].p1) >= interaction_distance_squared)
        {
          continue;
        }

        const double tangent_dot = Dot(segments[i].tangent, segments[j].tangent);
        if (std::abs(tangent_dot) < 1.0 - 1.0e-8)
        {
          Mpi::Warning(
              "Nearby three-dimensional metal edges are not parallel; correction is "
              "disabled for this interface group!\n");
          group_matched = false;
          break;
        }
        const double second_s0 =
            Dot(Subtract(segments[j].p0, segments[i].p0), segments[i].tangent);
        const double second_s1 =
            Dot(Subtract(segments[j].p1, segments[i].p0), segments[i].tangent);
        const double first_begin = std::max(0.0, std::min(second_s0, second_s1));
        const double first_end =
            std::min(segments[i].length, std::max(second_s0, second_s1));
        const double tolerance = 1.0e-10 * std::max({segments[i].length, segments[j].length,
                                                     group.matching_radius});
        if (first_end - first_begin <= tolerance)
        {
          continue;
        }
        const Point3D first_mid = Interpolate(segments[i], 0.5 * (first_begin + first_end));
        double second_mid = Dot(Subtract(first_mid, segments[j].p0), segments[j].tangent);
        second_mid = std::clamp(second_mid, 0.0, segments[j].length);
        const double half_length = 0.5 * (first_end - first_begin);
        const double second_begin = second_mid - half_length;
        const double second_end = second_mid + half_length;
        MFEM_VERIFY(second_begin >= -tolerance &&
                        second_end <= segments[j].length + tolerance,
                    "Inconsistent overlap between nearby parallel metal edges!");
        pairs.push_back({i, j, first_begin, first_end, std::max(0.0, second_begin),
                         std::min(segments[j].length, second_end)});
        paired_intervals[i].emplace_back(first_begin, first_end);
        paired_intervals[j].emplace_back(std::max(0.0, second_begin),
                                         std::min(segments[j].length, second_end));
      }
    }

    for (auto &intervals : paired_intervals)
    {
      std::sort(intervals.begin(), intervals.end());
      for (std::size_t i = 1; i < intervals.size(); i++)
      {
        const double tolerance =
            1.0e-10 * std::max(group.matching_radius, intervals[i - 1].second);
        if (intervals[i].first < intervals[i - 1].second - tolerance)
        {
          Mpi::Warning("More than two three-dimensional metal edges interact within one 2R "
                       "neighborhood; correction is disabled for this interface group!\n");
          group_matched = false;
          break;
        }
      }
      if (!group_matched)
      {
        break;
      }
    }

    std::vector<PendingPatch> pending;
    auto AppendQuadrature = [&](std::size_t model_index, const EdgeSegment3D &first,
                                double begin, double end, const EdgeSegment3D *second)
    {
      const auto &source = library.models[model_index];
      MFEM_VERIFY(source.coupon_depth > 0.0,
                  "Three-dimensional response correction requires CouponDepth for every "
                  "selected fabrication-process response model!");
      for (int q = 0; q < quadrature.GetNPoints(); q++)
      {
        const auto &ip = quadrature.IntPoint(q);
        const double first_distance = begin + (end - begin) * ip.x;
        const Point3D first_point = Interpolate(first, first_distance);
        ResponsePatchData patch;
        if (second)
        {
          double second_distance = Dot(Subtract(first_point, second->p0), second->tangent);
          second_distance = std::clamp(second_distance, 0.0, second->length);
          const Point3D second_point = Interpolate(*second, second_distance);
          const Point3D direction = Normalize(Subtract(second_point, first_point));
          patch.origin = Scale(0.5, Add(first_point, second_point));
          patch.axis_u = direction;
          patch.axis_v = Normalize(Add(first.axis_v, second->axis_v));
        }
        else
        {
          patch.origin = first_point;
          patch.axis_u = first.axis_u;
          patch.axis_v = first.axis_v;
        }
        if (source.topology == LibraryTopology::SAME_CONDUCTOR_STRIP)
        {
          patch.maxwell_anchor = patch.origin;
        }
        else
        {
          patch.maxwell_anchor =
              Add(first_point, Scale(-group.matching_radius, first.axis_u));
        }
        patch.reference = source.reference;
        patch.weight = (end - begin) * ip.weight / source.coupon_depth;
        pending.push_back({model_index, patch});
      }
    };

    for (const auto &pair : pairs)
    {
      const auto &first = segments[pair.first];
      const auto &second = segments[pair.second];
      const Point3D first_mid =
          Interpolate(first, 0.5 * (pair.first_begin + pair.first_end));
      double second_distance = Dot(Subtract(first_mid, second.p0), second.tangent);
      second_distance = std::clamp(second_distance, 0.0, second.length);
      const Point3D second_mid = Interpolate(second, second_distance);
      const Point3D direction = Normalize(Subtract(second_mid, first_mid));
      if (Dot(first.axis_v, second.axis_v) <= 0.95)
      {
        Mpi::Warning(
            "Nearby three-dimensional edges have incompatible process normals; a "
            "dedicated cross-layer coupon is required and correction is disabled for "
            "this interface group!\n");
        group_matched = false;
        break;
      }
      const Point3D process_normal = Normalize(Add(first.axis_v, second.axis_v));
      if (std::abs(Dot(direction, process_normal)) > 1.0e-8)
      {
        Mpi::Warning(
            "Nearby three-dimensional edges are offset along the process normal; a "
            "dedicated cross-layer coupon is required and correction is disabled for "
            "this interface group!\n");
        group_matched = false;
        break;
      }
      const bool facing =
          Dot(first.axis_u, direction) > 0.95 && Dot(second.axis_u, direction) < -0.95;
      const bool outward =
          Dot(first.axis_u, direction) < -0.95 && Dot(second.axis_u, direction) > 0.95;
      const bool same_conductor = first.conductor == second.conductor;
      LibraryTopology topology;
      if (facing)
      {
        topology = same_conductor ? LibraryTopology::SAME_CONDUCTOR_GAP
                                  : LibraryTopology::DIFFERENT_CONDUCTOR_GAP;
      }
      else if (outward && same_conductor)
      {
        topology = LibraryTopology::SAME_CONDUCTOR_STRIP;
      }
      else
      {
        Mpi::Warning("No canonical paired-edge topology for nearby three-dimensional metal "
                     "edges; correction is disabled for this interface group!\n");
        group_matched = false;
        break;
      }
      const double separation = Distance(first_mid, second_mid);
      const auto model_index = FindLibraryModel(library, topology, separation);
      if (!model_index)
      {
        Mpi::Warning(
            "Fabrication-process response library \"{}\" has no {} model at separation "
            "{:.6e} mesh units; correction is disabled for this interface group!\n",
            library.name, TopologyName(topology), separation * coordinate_scale);
        group_matched = false;
        break;
      }
      const auto &matched_model = library.models[*model_index];
      const double model_tolerance =
          std::max(matched_model.separation_tolerance,
                   1.0e-10 * std::max(library.matching_radius, separation));
      group_maximum_library_distance =
          std::max(group_maximum_library_distance,
                   std::abs(matched_model.separation - separation) / model_tolerance);
      AppendQuadrature(*model_index, first, pair.first_begin, pair.first_end, &second);
      group_matched_intervals++;
    }

    const auto isolated_model =
        FindLibraryModel(library, LibraryTopology::ISOLATED_EDGE, 0.0);
    if (group_matched && !isolated_model)
    {
      Mpi::Warning(
          "Fabrication-process response library \"{}\" has no isolated-edge model; "
          "correction is disabled for this interface group!\n",
          library.name);
      group_matched = false;
    }
    for (std::size_t i = 0; group_matched && i < segments.size(); i++)
    {
      double begin = 0.0;
      for (const auto &[paired_begin, paired_end] : paired_intervals[i])
      {
        if (paired_begin > begin)
        {
          AppendQuadrature(*isolated_model, segments[i], begin, paired_begin, nullptr);
          group_matched_intervals++;
        }
        begin = std::max(begin, paired_end);
      }
      if (begin < segments[i].length)
      {
        AppendQuadrature(*isolated_model, segments[i], begin, segments[i].length, nullptr);
        group_matched_intervals++;
      }
    }

    if (!group_matched)
    {
      unmatched_groups++;
      if (request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
      {
        MFEM_ABORT("Automatic fabrication-process response matching failed!");
      }
      continue;
    }

    std::map<std::size_t, int> runtime_models;
    for (auto &selection : pending)
    {
      auto [model_it, inserted] =
          runtime_models.emplace(selection.library_model, next_model_index);
      if (inserted)
      {
        const auto &source = library.models[selection.library_model];
        auto model = source.response;
        model.idx = next_model_index++;
        for (const auto &[type, target] : group.targets)
        {
          const auto interface_type = type;
          auto interface = std::find_if(source.interfaces.begin(), source.interfaces.end(),
                                        [interface_type](const auto &entry)
                                        { return entry.type == interface_type; });
          MFEM_VERIFY(interface != source.interfaces.end() ||
                          source.response.fabricated_surface_matrix.empty(),
                      "Fabrication-process response model \"" << source.name << "\" has no "
                                                              << ToString(type)
                                                              << " interface mapping!");
          if (interface != source.interfaces.end())
          {
            model.interfaces.push_back({target, interface->coupon});
          }
        }
        result.models.push_back(std::move(model));
      }
      selection.patch.model = model_it->second;
      result.patches.push_back(selection.patch);
    }
    matched_intervals += group_matched_intervals;
    matched_segments += static_cast<int>(segments.size());
    if (diagnostics)
    {
      diagnostics->matched_length += group_selected_length;
      diagnostics->matched_corner_neighborhood_length +=
          group_corner_neighborhood_length;
      diagnostics->maximum_curvature_ratio =
          std::max(diagnostics->maximum_curvature_ratio,
                   group_maximum_curvature_ratio);
      diagnostics->maximum_library_distance =
          std::max(diagnostics->maximum_library_distance,
                   group_maximum_library_distance);
    }
  }

  MFEM_VERIFY(!result.models.empty() && !result.patches.empty(),
              "Fabrication-process response matching produced no usable correction "
              "patches!");
  Mpi::Print("\nAutomatic fabrication-process response matching:\n"
             " Library: {}\n"
             " Matched physical edge segments: {:d}\n"
             " Matched longitudinal intervals: {:d}\n"
             " Longitudinal quadrature patches: {:d}\n"
             " Unmatched interface groups: {:d}\n",
             library.name, matched_segments, matched_intervals,
             static_cast<int>(result.patches.size()), unmatched_groups);
  if (nonregular_vertices > 0)
  {
    Mpi::Warning(
        "The selected three-dimensional metal perimeter has {} corner, endpoint, or "
        "junction vertices. Straight-edge coupon response is integrated through these "
        "neighborhoods; validate or add a corner-response model when their contribution "
        "is significant.\n",
        nonregular_vertices);
  }
  return result;
}

ResponseCorrectionData BuildAutomaticResponseData(const IoData &iodata,
                                                  const LaplaceOperator &laplace_op,
                                                  const ResponseCorrectionData &request)
{
  const auto &mesh = laplace_op.GetH1Space().GetParMesh();
  if (mesh.Dimension() == 2 && mesh.SpaceDimension() == 2)
  {
    return BuildAutomaticResponseData2D(iodata, laplace_op, request);
  }
  if (mesh.Dimension() == 3 && mesh.SpaceDimension() == 3)
  {
    return BuildAutomaticResponseData3D(iodata, mesh, laplace_op.GetMaterialOp(), request,
                                        false);
  }
  MFEM_ABORT("Automatic fabrication-process response matching requires a 2D or 3D "
             "electrostatic mesh!");
}

std::vector<std::array<double, 3>> ReadBasisPoints(const std::string &path)
{
  std::ifstream input(path);
  MFEM_VERIFY(input,
              "Unable to open response-correction basis point file \"" << path << "\"!");

  std::vector<std::array<double, 3>> points;
  std::string line;
  while (std::getline(input, line))
  {
    auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#')
    {
      continue;
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    std::array<double, 3> point;
    if (!(row >> point[0] >> point[1] >> point[2]))
    {
      MFEM_VERIFY(points.empty(), "Could not parse response-correction basis point file \""
                                      << path << "\"!");
      continue;  // Optional header before the first data row.
    }
    MFEM_VERIFY(
        std::all_of(point.begin(), point.end(), [](double x) { return std::isfinite(x); }),
        "Non-finite coordinate in response-correction basis point file \"" << path
                                                                           << "\"!");
    points.push_back(point);
  }
  MFEM_VERIFY(!points.empty(),
              "Response-correction basis point file \"" << path << "\" is empty!");
  return points;
}

Table ReadTable(const std::string &path)
{
  TableWithCSVFile input(path, true);
  MFEM_VERIFY(!input.table.empty(),
              "Unable to read response matrix file \"" << path << "\"!");
  return std::move(input.table);
}

const Column &FindColumn(const Table &table, const std::string &header,
                         const std::string &path)
{
  for (auto it = table.cbegin(); it != table.cend(); ++it)
  {
    if (it->header_text == header)
    {
      return *it;
    }
  }
  MFEM_ABORT("Response matrix file \"" << path << "\" is missing column \"" << header
                                       << "\"!");
}

int ParseIndex(double value, const std::string &name, const std::string &path)
{
  const int idx = static_cast<int>(value);
  MFEM_VERIFY(idx > 0 && value == idx,
              "Invalid " << name << " index in response matrix file \"" << path << "\"!");
  return idx;
}

std::pair<int, std::vector<MatrixEntry>> ReadDomainResponseMatrix(const std::string &path)
{
  const Table table = ReadTable(path);
  const auto &basis_i = FindColumn(table, "basis_i", path);
  const auto &basis_j = FindColumn(table, "basis_j", path);
  const auto &q = FindColumn(table, "Q_ij (J)", path);
  MFEM_VERIFY(basis_i.n_rows() == basis_j.n_rows() && basis_i.n_rows() == q.n_rows(),
              "Response matrix columns have inconsistent lengths in \"" << path << "\"!");

  int size = 0;
  std::vector<MatrixEntry> entries;
  entries.reserve(q.n_rows());
  for (std::size_t row = 0; row < q.n_rows(); row++)
  {
    const int i = ParseIndex(basis_i.data[row], "basis_i", path);
    const int j = ParseIndex(basis_j.data[row], "basis_j", path);
    MFEM_VERIFY(j >= i && std::isfinite(q.data[row]),
                "Invalid response matrix entry in \"" << path << "\"!");
    entries.emplace_back(i - 1, j - 1, q.data[row]);
    size = std::max(size, j);
  }
  MFEM_VERIFY(size > 0, "Response matrix file \"" << path << "\" contains no entries!");
  return {size, std::move(entries)};
}

mfem::DenseMatrix BuildDenseMatrix(const std::vector<MatrixEntry> &entries, int size,
                                   const std::string &path)
{
  mfem::DenseMatrix matrix(size);
  matrix = 0.0;
  std::vector<bool> have(size * size, false);
  for (const auto &[i, j, value] : entries)
  {
    MFEM_VERIFY(i >= 0 && i < size && j >= i && j < size && !have[i * size + j],
                "Duplicate or out-of-range response matrix entry in \"" << path << "\"!");
    matrix(i, j) = matrix(j, i) = value;
    have[i * size + j] = have[j * size + i] = true;
  }
  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      MFEM_VERIFY(have[i * size + j], "Response matrix file \""
                                          << path
                                          << "\" is missing an upper-triangular entry!");
    }
  }
  return matrix;
}

struct DomainResponseMatrices
{
  mfem::DenseMatrix fabricated;
  mfem::DenseMatrix thin;
  mfem::DenseMatrix defect;
  mfem::DenseMatrix fixed_flux_transform;
};

DomainResponseMatrices BuildDomainResponseMatrices(const std::string &fabricated_path,
                                                   const std::string &thin_path,
                                                   int expected_size,
                                                   const Units &units)
{
  auto [fabricated_size, fabricated_entries] = ReadDomainResponseMatrix(fabricated_path);
  auto [thin_size, thin_entries] = ReadDomainResponseMatrix(thin_path);
  MFEM_VERIFY(fabricated_size == expected_size && thin_size == expected_size,
              "Response matrices and basis point file have inconsistent sizes!");

  auto fabricated = BuildDenseMatrix(fabricated_entries, expected_size, fabricated_path);
  auto thin = BuildDenseMatrix(thin_entries, expected_size, thin_path);

  // The CSV stores coupon energy Q in joules for basis traces measured in volts.
  // Internally, 1/2 xᵀ C x must equal the nondimensional energy defect, so
  // C = 2 V_scale² / E_scale Q.
  const double voltage_scale = units.GetScaleFactor<Units::ValueType::VOLTAGE>();
  const double energy_scale = units.GetScaleFactor<Units::ValueType::ENERGY>();
  const double scale = 2.0 * voltage_scale * voltage_scale / energy_scale;
  fabricated *= scale;
  thin *= scale;

  mfem::DenseMatrix defect(fabricated);
  defect.Add(-1.0, thin);
  mfem::DenseMatrix fixed_flux_transform;
  mfem::DenseMatrixInverse(fabricated, true).Mult(thin, fixed_flux_transform);
  return {std::move(fabricated), std::move(thin), std::move(defect),
          std::move(fixed_flux_transform)};
}

std::map<int, mfem::DenseMatrix> ReadSurfaceResponseMatrices(const std::string &path,
                                                             int expected_size)
{
  const Table table = ReadTable(path);
  const auto &interface_col = FindColumn(table, "interface", path);
  const auto &edge_col = FindColumn(table, "edge", path);
  const auto &basis_i = FindColumn(table, "basis_i", path);
  const auto &basis_j = FindColumn(table, "basis_j", path);
  const auto &q = FindColumn(table, "Q_total_ij (J)", path);
  const std::size_t rows = q.n_rows();
  MFEM_VERIFY(interface_col.n_rows() == rows && edge_col.n_rows() == rows &&
                  basis_i.n_rows() == rows && basis_j.n_rows() == rows,
              "Surface response matrix columns have inconsistent lengths in \"" << path
                                                                                << "\"!");

  // Q_total is repeated for every matching radius. Deduplicate those rows for each
  // interface/edge/basis pair, then sum all physical coupon edges belonging to an
  // interface. This makes one matrix represent either a one-edge or a coupled multi-edge
  // coupon.
  using Key = std::tuple<int, int, int, int>;
  std::map<Key, double> unique;
  for (std::size_t row = 0; row < rows; row++)
  {
    const int interface = ParseIndex(interface_col.data[row], "interface", path);
    const int edge = ParseIndex(edge_col.data[row], "edge", path);
    const int i = ParseIndex(basis_i.data[row], "basis_i", path);
    const int j = ParseIndex(basis_j.data[row], "basis_j", path);
    MFEM_VERIFY(i <= expected_size && j >= i && j <= expected_size &&
                    std::isfinite(q.data[row]),
                "Invalid surface response matrix entry in \"" << path << "\"!");
    auto [it, inserted] = unique.emplace(Key{interface, edge, i - 1, j - 1}, q.data[row]);
    if (!inserted)
    {
      const double scale =
          std::max({std::abs(it->second), std::abs(q.data[row]), 1.0e-300});
      MFEM_VERIFY(std::abs(it->second - q.data[row]) <= 1.0e-10 * scale,
                  "Inconsistent repeated Q_total entry in \"" << path << "\"!");
    }
  }

  std::map<int, std::vector<MatrixEntry>> entries;
  std::map<std::pair<int, int>, std::vector<bool>> have;
  for (const auto &[key, value] : unique)
  {
    const auto [interface, edge, i, j] = key;
    auto &edge_have = have[{interface, edge}];
    if (edge_have.empty())
    {
      edge_have.resize(expected_size * expected_size, false);
    }
    edge_have[i * expected_size + j] = true;
    edge_have[j * expected_size + i] = true;

    auto &interface_entries = entries[interface];
    const int basis_row = i;
    const int basis_col = j;
    auto it = std::find_if(
        interface_entries.begin(), interface_entries.end(),
        [basis_row, basis_col](const auto &entry)
        { return std::get<0>(entry) == basis_row && std::get<1>(entry) == basis_col; });
    if (it == interface_entries.end())
    {
      interface_entries.emplace_back(i, j, value);
    }
    else
    {
      std::get<2>(*it) += value;
    }
  }
  for (const auto &[key, edge_have] : have)
  {
    MFEM_VERIFY(
        std::all_of(edge_have.begin(), edge_have.end(), [](bool value) { return value; }),
        "Surface response matrix file \"" << path << "\" has an incomplete edge matrix!");
  }

  std::map<int, mfem::DenseMatrix> matrices;
  for (const auto &[interface, interface_entries] : entries)
  {
    matrices.emplace(interface, BuildDenseMatrix(interface_entries, expected_size, path));
  }
  return matrices;
}

struct SurfaceResponseMatrices
{
  std::map<int, mfem::DenseMatrix> fabricated;
  std::map<int, mfem::DenseMatrix> defects;
};

SurfaceResponseMatrices BuildSurfaceResponseMatrices(
    const config::ElectrostaticSolverData::ResponseCorrectionModelData &config,
    int expected_size, const Units &units)
{
  if (config.fabricated_surface_matrix.empty() && config.thin_surface_matrix.empty())
  {
    MFEM_VERIFY(
        config.interfaces.empty(),
        "Response-correction interface mappings require surface response matrices!");
    return {};
  }
  MFEM_VERIFY(!config.fabricated_surface_matrix.empty() &&
                  !config.thin_surface_matrix.empty() && !config.interfaces.empty(),
              "FabricatedSurfaceMatrix, ThinSurfaceMatrix, and Interfaces must be "
              "specified together for response-corrected surface participation!");

  auto fabricated =
      ReadSurfaceResponseMatrices(config.fabricated_surface_matrix, expected_size);
  auto thin = ReadSurfaceResponseMatrices(config.thin_surface_matrix, expected_size);
  const double voltage_scale = units.GetScaleFactor<Units::ValueType::VOLTAGE>();
  const double energy_scale = units.GetScaleFactor<Units::ValueType::ENERGY>();
  const double scale = voltage_scale * voltage_scale / energy_scale;

  SurfaceResponseMatrices result;
  for (const auto &mapping : config.interfaces)
  {
    const auto fabricated_it = fabricated.find(mapping.coupon);
    const auto thin_it = thin.find(mapping.coupon);
    MFEM_VERIFY(mapping.target > 0 && mapping.coupon > 0 &&
                    fabricated_it != fabricated.end() && thin_it != thin.end(),
                "Response-correction interface mapping refers to a missing coupon "
                "surface response!");
    MFEM_VERIFY(result.defects.find(mapping.target) == result.defects.end(),
                "Duplicate target interface in one response-correction model!");
    auto [fabricated_target, fabricated_inserted] =
        result.fabricated.emplace(mapping.target, fabricated_it->second);
    fabricated_target->second *= scale;
    auto [defect, defect_inserted] =
        result.defects.emplace(mapping.target, fabricated_it->second);
    defect->second.Add(-1.0, thin_it->second);
    defect->second *= scale;
  }
  return result;
}

double Orientation(const Point2D &a, const Point2D &b, const Point2D &c)
{
  return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

bool PointStrictlyInside(const Point2D &point, const std::vector<Point2D> &polygon,
                         double tol)
{
  bool inside = false;
  for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++)
  {
    const auto &a = polygon[j];
    const auto &b = polygon[i];
    const double cross = Orientation(a, b, point);
    if (std::abs(cross) <= tol && point[0] >= std::min(a[0], b[0]) - tol &&
        point[0] <= std::max(a[0], b[0]) + tol && point[1] >= std::min(a[1], b[1]) - tol &&
        point[1] <= std::max(a[1], b[1]) + tol)
    {
      return false;
    }
    if ((a[1] > point[1]) != (b[1] > point[1]))
    {
      const double x = a[0] + (point[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
      if (x > point[0] + tol)
      {
        inside = !inside;
      }
    }
  }
  return inside;
}

bool PolygonsOverlap(const std::vector<Point2D> &a, const std::vector<Point2D> &b)
{
  if (a.size() < 3 || b.size() < 3)
  {
    return false;
  }
  double coordinate_scale = 1.0;
  for (const auto &point : a)
  {
    coordinate_scale = std::max({coordinate_scale, std::abs(point[0]), std::abs(point[1])});
  }
  for (const auto &point : b)
  {
    coordinate_scale = std::max({coordinate_scale, std::abs(point[0]), std::abs(point[1])});
  }
  const double tol = 1.0e-12 * coordinate_scale;

  for (std::size_t i = 0; i < a.size(); i++)
  {
    const auto &a0 = a[i];
    const auto &a1 = a[(i + 1) % a.size()];
    for (std::size_t j = 0; j < b.size(); j++)
    {
      const auto &b0 = b[j];
      const auto &b1 = b[(j + 1) % b.size()];
      const double o0 = Orientation(a0, a1, b0);
      const double o1 = Orientation(a0, a1, b1);
      const double o2 = Orientation(b0, b1, a0);
      const double o3 = Orientation(b0, b1, a1);
      if (o0 * o1 < -tol * tol && o2 * o3 < -tol * tol)
      {
        return true;
      }
    }
  }
  auto HasInteriorSample =
      [tol](const std::vector<Point2D> &source, const std::vector<Point2D> &target)
  {
    for (std::size_t i = 0; i < source.size(); i++)
    {
      if (PointStrictlyInside(source[i], target, tol))
      {
        return true;
      }
      const auto &next = source[(i + 1) % source.size()];
      const Point2D midpoint = {0.5 * (source[i][0] + next[0]),
                                0.5 * (source[i][1] + next[1])};
      if (PointStrictlyInside(midpoint, target, tol))
      {
        return true;
      }
    }
    return false;
  };
  if (HasInteriorSample(a, b) || HasInteriorSample(b, a))
  {
    return true;
  }

  // Identical contours have no strict edge intersection or interior boundary sample.
  if (a.size() == b.size())
  {
    bool identical = true;
    for (std::size_t i = 0; i < a.size(); i++)
    {
      identical = identical && std::abs(a[i][0] - b[i][0]) <= tol &&
                  std::abs(a[i][1] - b[i][1]) <= tol;
    }
    if (identical)
    {
      return true;
    }
  }
  return false;
}

double QuadraticForm(const mfem::DenseMatrix &matrix, const Vector &x, Vector &workspace)
{
  workspace.SetSize(x.Size());
  matrix.Mult(x, workspace);
  return x * workspace;
}

}  // namespace

SurfaceResponseOperator::SurfaceResponseOperator(const IoData &iodata,
                                                 const LaplaceOperator &laplace_op)
  : Operator(laplace_op.GetH1Space().GetTrueVSize()), fespace(laplace_op.GetH1Space()),
    basis_size(0)
{
  const auto &request = iodata.solver.electrostatic.response_correction;
  MFEM_VERIFY(request, "Missing electrostatic surface response correction configuration!");
  const int dimension = fespace.Dimension();
  MFEM_VERIFY((dimension == 2 || dimension == 3) && fespace.SpaceDimension() == dimension,
              "Surface response correction requires a 2D or 3D electrostatic mesh!");
  MFEM_VERIFY(dimension == 2 || request->IsAutomatic(),
              "Three-dimensional surface response correction requires automatic "
              "fabrication-process library matching!");
  std::optional<ResponseCorrectionData> automatic_config;
  const ResponseCorrectionData *config = &*request;
  if (request->IsAutomatic())
  {
    automatic_config = BuildAutomaticResponseData(iodata, laplace_op, *request);
    config = &*automatic_config;
  }
  MFEM_VERIFY(!config->models.empty() && !config->patches.empty(),
              "Surface response correction requires at least one model and patch!");

#if defined(MFEM_USE_GSLIB)
  std::unordered_map<int, int> model_indices;
  std::vector<std::vector<std::array<double, 3>>> basis_points;
  models.reserve(config->models.size());
  basis_points.reserve(config->models.size());
  for (const auto &model_config : config->models)
  {
    MFEM_VERIFY(model_config.idx > 0 &&
                    model_indices.find(model_config.idx) == model_indices.end(),
                "Response-correction model indices must be positive and unique!");
    auto points = ReadBasisPoints(model_config.basis_points);
    ResponseModel model;
    model.idx = model_config.idx;
    model.basis_size = static_cast<int>(points.size());
    auto domain_response =
        BuildDomainResponseMatrices(model_config.fabricated_matrix,
                                    model_config.thin_matrix, model.basis_size,
                                    iodata.units);
    model.fabricated_domain = std::move(domain_response.fabricated);
    model.thin_domain = std::move(domain_response.thin);
    model.domain_defect = std::move(domain_response.defect);
    model.fixed_flux_transform = std::move(domain_response.fixed_flux_transform);
    auto surface_response =
        BuildSurfaceResponseMatrices(model_config, model.basis_size, iodata.units);
    model.fabricated_surfaces = std::move(surface_response.fabricated);
    model.surface_defects = std::move(surface_response.defects);
    model_indices.emplace(model.idx, static_cast<int>(models.size()));
    models.push_back(std::move(model));
    basis_points.push_back(std::move(points));
  }
  dbc_tdof_list = laplace_op.GetDbcTDofList();

  int point_count = 0;
  for (const auto &patch_config : config->patches)
  {
    const auto model_it = model_indices.find(patch_config.model);
    MFEM_VERIFY(model_it != model_indices.end(),
                "Response-correction patch refers to an unknown model index!");
    const auto &model = models[model_it->second];
    MFEM_VERIFY(std::isfinite(patch_config.weight) && patch_config.weight > 0.0,
                "Response-correction patch weights must be positive!");
    patches.push_back(
        Patch{model_it->second, point_count, basis_size, patch_config.weight});
    point_count += model.basis_size + 1;
    basis_size += model.basis_size;
  }

  mfem::Vector xyz(dimension * point_count);
  const double coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
  std::vector<std::vector<Point2D>> polygons;
  if (dimension == 2)
  {
    polygons.reserve(patches.size());
  }
  int point = 0;
  for (std::size_t patch_idx = 0; patch_idx < patches.size(); patch_idx++)
  {
    const auto &patch = patches[patch_idx];
    const auto &patch_config = config->patches[patch_idx];
    const auto &local_points = basis_points[patch.model];
    double norm_u = 0.0, norm_v = 0.0, dot = 0.0;
    for (int d = 0; d < dimension; d++)
    {
      norm_u += patch_config.axis_u[d] * patch_config.axis_u[d];
      norm_v += patch_config.axis_v[d] * patch_config.axis_v[d];
      dot += patch_config.axis_u[d] * patch_config.axis_v[d];
    }
    norm_u = std::sqrt(norm_u);
    norm_v = std::sqrt(norm_v);
    MFEM_VERIFY(std::abs(norm_u - 1.0) < 1.0e-10 && std::abs(norm_v - 1.0) < 1.0e-10 &&
                    std::abs(dot) < 1.0e-10,
                "Surface response correction AxisU and AxisV must be orthonormal in the "
                "coupon cross-section!");

    std::vector<Point2D> *polygon = nullptr;
    if (dimension == 2)
    {
      polygon = &polygons.emplace_back();
      polygon->reserve(local_points.size());
    }
    for (const auto &local : local_points)
    {
      MFEM_VERIFY(std::abs(local[2]) <= 1.0e-12,
                  "Response-correction basis points must lie in the local coupon "
                  "cross-section!");
      Point2D global{};
      for (int d = 0; d < dimension; d++)
      {
        const double coordinate =
            patch_config.origin[d] +
            (local[0] * patch_config.axis_u[d] + local[1] * patch_config.axis_v[d]) /
                coordinate_scale;
        xyz(d * point_count + point) = coordinate;
        if (d < 2)
        {
          global[d] = coordinate;
        }
      }
      if (polygon)
      {
        polygon->push_back(global);
      }
      point++;
    }
    for (int d = 0; d < dimension; d++)
    {
      xyz(d * point_count + point) = patch_config.origin[d] +
                                     patch_config.reference[0] * patch_config.axis_u[d] +
                                     patch_config.reference[1] * patch_config.axis_v[d];
    }
    point++;
  }
  for (std::size_t i = 0; i < polygons.size(); i++)
  {
    for (std::size_t j = i + 1; j < polygons.size(); j++)
    {
      MFEM_VERIFY(!PolygonsOverlap(polygons[i], polygons[j]),
                  "Response-correction patches " << i + 1 << " and " << j + 1
                                                 << " overlap. Replace nearby one-edge "
                                                    "patches with one coupled multi-edge "
                                                    "coupon model!");
    }
  }

  mfem::FindPointsGSLIB finder(fespace.GetComm());
  finder.Setup(const_cast<mfem::ParMesh &>(fespace.GetParMesh()), 0.01, 1.0e-12,
               point_count);
  finder.FindPoints(xyz, mfem::Ordering::byNODES);
  points.resize(point_count);
  const int rank = Mpi::Rank(fespace.GetComm());
  for (int i = 0; i < point_count; i++)
  {
    MFEM_VERIFY(finder.GetCode()[i] != 2,
                "Surface response contour point " << i << " could not be located!");
    auto &evaluation = points[i];
    evaluation.local = static_cast<int>(finder.GetProc()[i]) == rank;
    if (evaluation.local)
    {
      evaluation.element = static_cast<int>(finder.GetElem()[i]);
      if (dimension == 2)
      {
        evaluation.point.Set2(finder.GetReferencePosition()(2 * i),
                              finder.GetReferencePosition()(2 * i + 1));
      }
      else
      {
        evaluation.point.Set3(finder.GetReferencePosition()(3 * i),
                              finder.GetReferencePosition()(3 * i + 1),
                              finder.GetReferencePosition()(3 * i + 2));
      }
    }
  }

  Mpi::Print("\nConfigured surface response correction:\n"
             " Coupon models: {:d}\n"
             " Global patches: {:d}\n"
             " Total trace coefficients: {:d}\n",
             static_cast<int>(models.size()), static_cast<int>(patches.size()), basis_size);
#else
  MFEM_ABORT("Surface response correction requires MFEM_USE_GSLIB!");
#endif
}

SurfaceResponseOperator::SurfaceResponseOperator(const IoData &iodata,
                                                 const SpaceOperator &space_op)
  : Operator(space_op.GetNDSpace().GetTrueVSize()), fespace(space_op.GetNDSpace()),
    basis_size(0)
{
  const auto &request = iodata.solver.surface_response_correction;
  MFEM_VERIFY(request,
              "Missing Maxwell surface response correction configuration!");
  MFEM_VERIFY(request->IsAutomatic(),
              "Maxwell surface response correction requires automatic fabrication-"
              "process library matching!");
  MFEM_VERIFY(fespace.Dimension() == 3 && fespace.SpaceDimension() == 3,
              "Maxwell surface response correction currently requires a three-"
              "dimensional mesh!");

#if defined(MFEM_USE_GSLIB)
  AutomaticResponseDiagnostics diagnostics;
  const auto config =
      BuildAutomaticResponseData3D(iodata, fespace.GetParMesh(),
                                   space_op.GetMaterialOp(), *request, true, &diagnostics);
  matching_radius = diagnostics.matching_radius;
  minimum_wave_speed = diagnostics.minimum_wave_speed;
  matched_length_fraction =
      diagnostics.selected_length > 0.0
          ? diagnostics.matched_length / diagnostics.selected_length
          : 0.0;
  corner_neighborhood_fraction =
      diagnostics.matched_length > 0.0
          ? diagnostics.matched_corner_neighborhood_length /
                diagnostics.matched_length
          : 0.0;
  maximum_curvature_ratio = diagnostics.maximum_curvature_ratio;
  maximum_library_distance = diagnostics.maximum_library_distance;

  std::unordered_map<int, int> model_indices;
  std::vector<std::vector<std::array<double, 3>>> basis_points;
  models.reserve(config.models.size());
  basis_points.reserve(config.models.size());
  for (const auto &model_config : config.models)
  {
    MFEM_VERIFY(model_config.idx > 0 &&
                    model_indices.find(model_config.idx) == model_indices.end(),
                "Response-correction model indices must be positive and unique!");
    auto points = ReadBasisPoints(model_config.basis_points);
    MFEM_VERIFY(points.size() >= 3,
                "Maxwell response-correction contours require at least three points!");
    ResponseModel model;
    model.idx = model_config.idx;
    model.basis_size = static_cast<int>(points.size());
    auto domain_response =
        BuildDomainResponseMatrices(model_config.fabricated_matrix,
                                    model_config.thin_matrix, model.basis_size,
                                    iodata.units);
    model.fabricated_domain = std::move(domain_response.fabricated);
    model.thin_domain = std::move(domain_response.thin);
    model.domain_defect = std::move(domain_response.defect);
    model.fixed_flux_transform = std::move(domain_response.fixed_flux_transform);
    auto surface_response =
        BuildSurfaceResponseMatrices(model_config, model.basis_size, iodata.units);
    model.fabricated_surfaces = std::move(surface_response.fabricated);
    model.surface_defects = std::move(surface_response.defects);
    model_indices.emplace(model.idx, static_cast<int>(models.size()));
    models.push_back(std::move(model));
    basis_points.push_back(std::move(points));
  }

  maxwell_contours.reserve(config.patches.size());
  maxwell_anchors.reserve(config.patches.size());
  const double coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
  for (const auto &patch_config : config.patches)
  {
    const auto model_it = model_indices.find(patch_config.model);
    MFEM_VERIFY(model_it != model_indices.end(),
                "Response-correction patch refers to an unknown model index!");
    const auto &model = models[model_it->second];
    MFEM_VERIFY(std::isfinite(patch_config.weight) && patch_config.weight > 0.0,
                "Response-correction patch weights must be positive!");

    double norm_u = 0.0, norm_v = 0.0, dot = 0.0;
    for (int d = 0; d < 3; d++)
    {
      norm_u += patch_config.axis_u[d] * patch_config.axis_u[d];
      norm_v += patch_config.axis_v[d] * patch_config.axis_v[d];
      dot += patch_config.axis_u[d] * patch_config.axis_v[d];
    }
    MFEM_VERIFY(std::abs(std::sqrt(norm_u) - 1.0) < 1.0e-10 &&
                    std::abs(std::sqrt(norm_v) - 1.0) < 1.0e-10 &&
                    std::abs(dot) < 1.0e-10,
                "Surface response correction AxisU and AxisV must be orthonormal in the "
                "coupon cross-section!");

    auto &contour = maxwell_contours.emplace_back();
    contour.reserve(basis_points[model_it->second].size());
    for (const auto &local : basis_points[model_it->second])
    {
      MFEM_VERIFY(std::abs(local[2]) <= 1.0e-12,
                  "Response-correction basis points must lie in the local coupon "
                  "cross-section!");
      mfem::Vector point(3);
      for (int d = 0; d < 3; d++)
      {
        point[d] =
            patch_config.origin[d] +
            (local[0] * patch_config.axis_u[d] +
             local[1] * patch_config.axis_v[d]) /
                coordinate_scale;
      }
      contour.push_back(std::move(point));
    }
    MFEM_VERIFY(patch_config.maxwell_anchor,
                "Automatic Maxwell response patch is missing its PEC voltage anchor!");
    mfem::Vector anchor(3);
    std::copy(patch_config.maxwell_anchor->begin(),
              patch_config.maxwell_anchor->end(), anchor.GetData());
    maxwell_anchors.push_back(std::move(anchor));

    patches.push_back(
        Patch{model_it->second, 0, basis_size, patch_config.weight});
    basis_size += model.basis_size;
  }

  maxwell_quadrature_order = 2 * std::max(1, iodata.solver.order) + 2;
  maxwell_finder = std::make_unique<mfem::FindPointsGSLIB>(fespace.GetComm());
  fem::SetupInterpolator(*maxwell_finder,
                         const_cast<mfem::ParMesh &>(fespace.GetParMesh()));

  Mpi::Print("\nConfigured PEC Maxwell surface response postprocessing:\n"
             " Coupon models: {:d}\n"
             " Longitudinal quadrature patches: {:d}\n"
             " Total contour coefficients: {:d}\n"
             " Matched edge-length fraction: {:.6f}\n"
             " Corner-neighborhood fraction: {:.6f}\n"
             " Maximum R/rho: {:.3e}\n"
             " Maximum normalized library distance: {:.3e}\n",
             static_cast<int>(models.size()), static_cast<int>(patches.size()),
             basis_size, matched_length_fraction, corner_neighborhood_fraction,
             maximum_curvature_ratio, maximum_library_distance);
#else
  MFEM_ABORT("Maxwell surface response correction requires MFEM_USE_GSLIB!");
#endif
}

void SurfaceResponseOperator::EvaluatePoints(const Vector &x, Vector &values) const
{
  local_x.SetSize(fespace.GetVSize());
  fespace.GetProlongationMatrix()->Mult(x, local_x);
  local_x.HostRead();
  values.SetSize(static_cast<int>(points.size()));
  values = 0.0;
  auto *value_data = values.HostWrite();
  for (std::size_t i = 0; i < points.size(); i++)
  {
    const auto &evaluation = points[i];
    if (!evaluation.local)
    {
      continue;
    }
    const auto &fe = *fespace.Get().GetFE(evaluation.element);
    shape.SetSize(fe.GetDof());
    fe.CalcShape(evaluation.point, shape);
    fespace.Get().GetElementDofs(evaluation.element, element_dofs);
    element_values.SetSize(element_dofs.Size());
    local_x.GetSubVector(element_dofs, element_values);
    value_data[i] = shape * element_values;
  }
  Mpi::GlobalSum(values.Size(), value_data, fespace.GetComm());
}

void SurfaceResponseOperator::AddPointTranspose(int point, double value, Vector &y) const
{
  const auto &evaluation = points[point];
  if (!evaluation.local || value == 0.0)
  {
    return;
  }
  const auto &fe = *fespace.Get().GetFE(evaluation.element);
  shape.SetSize(fe.GetDof());
  fe.CalcShape(evaluation.point, shape);
  shape *= value;
  fespace.Get().GetElementDofs(evaluation.element, element_dofs);
  y.AddElementVector(element_dofs, shape);
}

void SurfaceResponseOperator::ApplyTrace(const Vector &x, Vector &values) const
{
  EvaluatePoints(x, correction);
  values.SetSize(basis_size);
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    const double reference = correction(patch.point_offset + model.basis_size);
    for (int i = 0; i < model.basis_size; i++)
    {
      values(patch.trace_offset + i) = correction(patch.point_offset + i) - reference;
    }
  }
}

void SurfaceResponseOperator::ApplyTraceTranspose(const Vector &values, Vector &y) const
{
  local_y.SetSize(fespace.GetVSize());
  local_y = 0.0;
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    double reference = 0.0;
    for (int i = 0; i < model.basis_size; i++)
    {
      const double value = values(patch.trace_offset + i);
      AddPointTranspose(patch.point_offset + i, value, local_y);
      reference -= value;
    }
    AddPointTranspose(patch.point_offset + model.basis_size, reference, local_y);
  }
  y.SetSize(fespace.GetTrueVSize());
  fespace.GetProlongationMatrix()->MultTranspose(local_y, y);
}

void SurfaceResponseOperator::ApplyUneliminated(const Vector &x, Vector &y) const
{
  ApplyTrace(x, trace);
  response.SetSize(trace.Size());
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    Vector patch_trace(trace.GetData() + patch.trace_offset, model.basis_size);
    Vector patch_response(response.GetData() + patch.trace_offset, model.basis_size);
    model.domain_defect.Mult(patch_trace, patch_response);
    patch_response *= patch.weight;
  }
  ApplyTraceTranspose(response, y);
}

void SurfaceResponseOperator::Mult(const Vector &x, Vector &y) const
{
  x_free.SetSize(x.Size());
  x_free = x;
  x_free.SetSubVector(dbc_tdof_list, 0.0);
  ApplyUneliminated(x_free, y);
  y.SetSubVector(dbc_tdof_list, 0.0);
}

void SurfaceResponseOperator::EliminateRHS(const Vector &x, Vector &rhs) const
{
  ApplyUneliminated(x, correction);
  correction.SetSubVector(dbc_tdof_list, 0.0);
  rhs.Add(-1.0, correction);
}

SurfaceResponseOperator::EnergyCorrection
SurfaceResponseOperator::GetEnergyCorrection(const Vector &x) const
{
  ApplyTrace(x, trace);
  EnergyCorrection energy;
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    Vector patch_trace(trace.GetData() + patch.trace_offset, model.basis_size);
    energy.domain +=
        0.5 * patch.weight * QuadraticForm(model.domain_defect, patch_trace, response);
    for (const auto &[interface, defect] : model.surface_defects)
    {
      energy.interfaces[interface] +=
          patch.weight * QuadraticForm(defect, patch_trace, response);
    }
  }
  return energy;
}

std::map<int, double>
SurfaceResponseOperator::GetFabricatedSurfaceEnergy(const Vector &x) const
{
  ApplyTrace(x, trace);
  std::map<int, double> energy;
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    Vector patch_trace(trace.GetData() + patch.trace_offset, model.basis_size);
    for (const auto &[interface, response] : model.fabricated_surfaces)
    {
      energy[interface] +=
          patch.weight * QuadraticForm(response, patch_trace, this->response);
    }
  }
  return energy;
}

SurfaceResponseOperator::MaxwellResponse
SurfaceResponseOperator::GetMaxwellResponse(const GridFunction &E,
                                            std::complex<double> omega) const
{
  MFEM_VERIFY(maxwell_finder && maxwell_contours.size() == patches.size(),
              "Maxwell surface response was not configured!");
  MFEM_ASSERT(maxwell_anchors.size() == patches.size(),
              "Inconsistent Maxwell response anchor count!");
  MFEM_VERIFY(E.HasImag(),
              "Maxwell surface response requires a complex electric field!");

  MaxwellResponse result;
  result.kR = std::abs(omega) * matching_radius / minimum_wave_speed;
  result.matched_length_fraction = matched_length_fraction;
  result.corner_neighborhood_fraction = corner_neighborhood_fraction;
  result.maximum_curvature_ratio = maximum_curvature_ratio;
  result.maximum_library_distance = maximum_library_distance;

  Vector trace_real, trace_imag, fixed_flux_real, fixed_flux_imag, workspace;
  for (std::size_t patch_index = 0; patch_index < patches.size(); patch_index++)
  {
    const auto &patch = patches[patch_index];
    const auto &model = models[patch.model];
    const auto &contour = maxwell_contours[patch_index];
    MFEM_ASSERT(static_cast<int>(contour.size()) == model.basis_size,
                "Inconsistent Maxwell response contour size!");

    trace_real.SetSize(model.basis_size);
    trace_imag.SetSize(model.basis_size);
    fixed_flux_real.SetSize(model.basis_size);
    fixed_flux_imag.SetSize(model.basis_size);
    trace_real = 0.0;
    trace_imag = 0.0;

    const auto &anchor = maxwell_anchors[patch_index];
    int start = 0;
    double start_distance = mfem::infinity();
    for (int i = 0; i < model.basis_size; i++)
    {
      mfem::Vector delta(contour[i]);
      delta -= anchor;
      if (delta.Norml2() < start_distance)
      {
        start = i;
        start_distance = delta.Norml2();
      }
    }
    if (start_distance > 1.0e-14 * std::max(1.0, matching_radius))
    {
      const std::complex<double> anchor_integral = {
          fem::ComputeLineIntegral(*maxwell_finder, anchor, contour[start], E.Real(),
                                   maxwell_quadrature_order),
          fem::ComputeLineIntegral(*maxwell_finder, anchor, contour[start], E.Imag(),
                                   maxwell_quadrature_order)};
      trace_real[start] = -anchor_integral.real();
      trace_imag[start] = -anchor_integral.imag();
    }

    std::complex<double> loop_integral = 0.0;
    double loop_scale = 0.0;
    for (int offset = 0; offset < model.basis_size; offset++)
    {
      const int i = (start + offset) % model.basis_size;
      const int next = (i + 1) % model.basis_size;
      const std::complex<double> integral = {
          fem::ComputeLineIntegral(*maxwell_finder, contour[i], contour[next], E.Real(),
                                   maxwell_quadrature_order),
          fem::ComputeLineIntegral(*maxwell_finder, contour[i], contour[next], E.Imag(),
                                   maxwell_quadrature_order)};
      loop_integral += integral;
      loop_scale += std::abs(integral);
      if (next != start)
      {
        trace_real[next] = trace_real[i] - integral.real();
        trace_imag[next] = trace_imag[i] - integral.imag();
      }
    }
    if (loop_scale > 0.0)
    {
      result.loop_residual =
          std::max(result.loop_residual, std::abs(loop_integral) / loop_scale);
    }

    model.fixed_flux_transform.Mult(trace_real, fixed_flux_real);
    model.fixed_flux_transform.Mult(trace_imag, fixed_flux_imag);
    auto HermitianForm = [&](const mfem::DenseMatrix &matrix, const Vector &real,
                             const Vector &imag)
    {
      return QuadraticForm(matrix, real, workspace) +
             QuadraticForm(matrix, imag, workspace);
    };
    result.domain_correction +=
        0.5 * patch.weight *
        HermitianForm(model.domain_defect, trace_real, trace_imag);
    result.domain_correction_fixed_flux +=
        0.5 * patch.weight *
        (HermitianForm(model.fabricated_domain, fixed_flux_real, fixed_flux_imag) -
         HermitianForm(model.thin_domain, trace_real, trace_imag));
    for (const auto &[interface, matrix] : model.fabricated_surfaces)
    {
      const double fixed_trace_energy =
          patch.weight * HermitianForm(matrix, trace_real, trace_imag);
      const double fixed_flux_energy =
          patch.weight * HermitianForm(matrix, fixed_flux_real, fixed_flux_imag);
      result.fabricated_surface_energy[interface] += fixed_trace_energy;
      result.fabricated_surface_energy_fixed_flux[interface] += fixed_flux_energy;
    }
  }

  constexpr double maximum_kR = 0.1;
  constexpr double maximum_loop_residual = 0.05;
  constexpr double maximum_corner_fraction = 0.1;
  constexpr double maximum_curvature = 0.25;
  constexpr double minimum_coverage = 1.0 - 1.0e-10;
  constexpr double maximum_library_match_distance = 0.8;
  // Fixed-trace and fixed-flux closures are both admissible in postprocessing-only
  // correction. A material difference between them means that the unresolved
  // fabricated field is not determined accurately enough by the thin-model trace.
  constexpr double maximum_trace_closure_spread = 0.05;
  for (const auto &[interface, fixed_trace] : result.fabricated_surface_energy)
  {
    const auto fixed_flux = result.fabricated_surface_energy_fixed_flux.find(interface);
    MFEM_ASSERT(fixed_flux != result.fabricated_surface_energy_fixed_flux.end(),
                "Missing fixed-flux fabricated surface response!");
    const double scale = std::max(std::abs(fixed_trace), std::abs(fixed_flux->second));
    if (scale > 0.0)
    {
      result.maximum_trace_closure_spread =
          std::max(result.maximum_trace_closure_spread,
                   std::abs(fixed_trace - fixed_flux->second) / scale);
    }
  }
  result.confident =
      result.kR <= maximum_kR && result.loop_residual <= maximum_loop_residual &&
      result.matched_length_fraction >= minimum_coverage &&
      result.corner_neighborhood_fraction <= maximum_corner_fraction &&
      result.maximum_curvature_ratio <= maximum_curvature &&
      result.maximum_library_distance <= maximum_library_match_distance &&
      result.maximum_trace_closure_spread <= maximum_trace_closure_spread;
  return result;
}

bool SurfaceResponseOperator::HasSurfaceResponse() const
{
  return std::any_of(models.begin(), models.end(),
                     [](const auto &model) { return !model.fabricated_surfaces.empty(); });
}

std::set<int> SurfaceResponseOperator::GetTargetInterfaces() const
{
  std::set<int> interfaces;
  for (const auto &model : models)
  {
    for (const auto &[interface, matrix] : model.fabricated_surfaces)
    {
      (void)matrix;
      interfaces.insert(interface);
    }
  }
  return interfaces;
}

double SurfaceResponseOperator::GetPatchWeight() const
{
  double weight = 0.0;
  for (const auto &patch : patches)
  {
    weight += patch.weight;
  }
  return weight;
}

}  // namespace palace
