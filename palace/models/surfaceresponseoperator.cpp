// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponseoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
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
#include "models/boundarymodeoperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/edgedistance.hpp"
#include "utils/enum_string.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/metaledge.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"
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

constexpr double maximum_trace_closure_spread = 0.05;
constexpr double maximum_trace_closure_response_failure_fraction = 0.01;

struct ElementBox
{
  std::array<double, 3> min;
  std::array<double, 3> max;

  ElementBox()
  {
    min.fill(mfem::infinity());
    max.fill(-mfem::infinity());
  }

  void Add(const ElementBox &box)
  {
    for (int d = 0; d < 3; d++)
    {
      min[d] = std::min(min[d], box.min[d]);
      max[d] = std::max(max[d], box.max[d]);
    }
  }

  bool Contains(const std::array<double, 3> &point, int dimension, double tolerance) const
  {
    for (int d = 0; d < dimension; d++)
    {
      if (point[d] < min[d] - tolerance || point[d] > max[d] + tolerance)
      {
        return false;
      }
    }
    return true;
  }

  bool IntersectsSegment(const std::array<double, 3> &p0, const std::array<double, 3> &p1,
                         int dimension, double tolerance) const
  {
    double begin = 0.0, end = 1.0;
    for (int d = 0; d < dimension; d++)
    {
      const double delta = p1[d] - p0[d];
      if (delta == 0.0)
      {
        if (p0[d] < min[d] - tolerance || p0[d] > max[d] + tolerance)
        {
          return false;
        }
        continue;
      }
      double first = (min[d] - tolerance - p0[d]) / delta;
      double second = (max[d] + tolerance - p0[d]) / delta;
      if (first > second)
      {
        std::swap(first, second);
      }
      begin = std::max(begin, first);
      end = std::min(end, second);
      if (end < begin)
      {
        return false;
      }
    }
    return true;
  }
};

// Lightweight local point locator for response-contour evaluation. Unlike
// FindPointsGSLIB, this stores no replicated global volume-search structure.
class ElementPointLocator
{
private:
  static constexpr std::size_t leaf_size = 8;

  struct Node
  {
    ElementBox box;
    std::size_t begin = 0;
    std::size_t end = 0;
    int left = -1;
    int right = -1;

    bool IsLeaf() const { return left < 0; }
  };

  mfem::ParMesh &mesh;
  int dimension;
  bool linear_mesh;
  std::vector<ElementBox> element_boxes;
  std::vector<int> indices;
  std::vector<Node> nodes;

  ElementBox GetElementBox(int element)
  {
    ElementBox box;
    mfem::DenseMatrix points;
    if (mesh.GetNodes())
    {
      auto &transformation = *mesh.GetElementTransformation(element);
      const int order = std::max(1, transformation.Order());
      const auto *refined =
          mfem::GlobGeometryRefiner.Refine(transformation.GetGeometryType(), order);
      transformation.Transform(refined->RefPts, points);
    }
    else
    {
      const auto &mesh_element = *mesh.GetElement(element);
      const int *vertices = mesh_element.GetVertices();
      points.SetSize(dimension, mesh_element.GetNVertices());
      for (int j = 0; j < points.Width(); j++)
      {
        const double *coordinate = mesh.GetVertex(vertices[j]);
        for (int d = 0; d < dimension; d++)
        {
          points(d, j) = coordinate[d];
        }
      }
    }
    for (int j = 0; j < points.Width(); j++)
    {
      for (int d = 0; d < dimension; d++)
      {
        box.min[d] = std::min(box.min[d], points(d, j));
        box.max[d] = std::max(box.max[d], points(d, j));
      }
    }
    for (int d = dimension; d < 3; d++)
    {
      box.min[d] = box.max[d] = 0.0;
    }
    return box;
  }

  int Build(std::size_t begin, std::size_t end)
  {
    Node node;
    node.begin = begin;
    node.end = end;
    for (std::size_t i = begin; i < end; i++)
    {
      node.box.Add(element_boxes[indices[i]]);
    }

    const int node_index = static_cast<int>(nodes.size());
    nodes.push_back(node);
    if (end - begin <= leaf_size)
    {
      return node_index;
    }

    int axis = 0;
    for (int d = 1; d < dimension; d++)
    {
      if (node.box.max[d] - node.box.min[d] > node.box.max[axis] - node.box.min[axis])
      {
        axis = d;
      }
    }
    const std::size_t mid = begin + (end - begin) / 2;
    std::nth_element(indices.begin() + begin, indices.begin() + mid, indices.begin() + end,
                     [this, axis](int a, int b)
                     {
                       const double center_a =
                           element_boxes[a].min[axis] + element_boxes[a].max[axis];
                       const double center_b =
                           element_boxes[b].min[axis] + element_boxes[b].max[axis];
                       return center_a < center_b || (center_a == center_b && a < b);
                     });
    const int left = Build(begin, mid);
    const int right = Build(mid, end);
    nodes[node_index].left = left;
    nodes[node_index].right = right;
    return node_index;
  }

  void FindCandidates(int node_index, const std::array<double, 3> &point, double tolerance,
                      std::vector<int> &candidates) const
  {
    const auto &node = nodes[node_index];
    if (!node.box.Contains(point, dimension, tolerance))
    {
      return;
    }
    if (node.IsLeaf())
    {
      for (std::size_t i = node.begin; i < node.end; i++)
      {
        const int element = indices[i];
        if (element_boxes[element].Contains(point, dimension, tolerance))
        {
          candidates.push_back(element);
        }
      }
      return;
    }
    FindCandidates(node.left, point, tolerance, candidates);
    FindCandidates(node.right, point, tolerance, candidates);
  }

  void FindSegmentCandidates(int node_index, const std::array<double, 3> &p0,
                             const std::array<double, 3> &p1, double tolerance,
                             std::vector<int> &candidates) const
  {
    const auto &node = nodes[node_index];
    if (!node.box.IntersectsSegment(p0, p1, dimension, tolerance))
    {
      return;
    }
    if (node.IsLeaf())
    {
      for (std::size_t i = node.begin; i < node.end; i++)
      {
        const int element = indices[i];
        if (element_boxes[element].IntersectsSegment(p0, p1, dimension, tolerance))
        {
          candidates.push_back(element);
        }
      }
      return;
    }
    FindSegmentCandidates(node.left, p0, p1, tolerance, candidates);
    FindSegmentCandidates(node.right, p0, p1, tolerance, candidates);
  }

  bool GetLinearSimplexReference(int element, const std::array<double, 3> &point,
                                 mfem::IntegrationPoint &reference) const
  {
    const auto &mesh_element = *mesh.GetElement(element);
    const auto geometry = mesh_element.GetGeometryType();
    if ((dimension == 2 && geometry != mfem::Geometry::TRIANGLE) ||
        (dimension == 3 && geometry != mfem::Geometry::TETRAHEDRON))
    {
      return false;
    }
    const int *vertices = mesh_element.GetVertices();
    const double *v0 = mesh.GetVertex(vertices[0]);
    const double *v1 = mesh.GetVertex(vertices[1]);
    const double *v2 = mesh.GetVertex(vertices[2]);
    const std::array<double, 3> a = {v1[0] - v0[0], v1[1] - v0[1],
                                     dimension == 3 ? v1[2] - v0[2] : 0.0};
    const std::array<double, 3> b = {v2[0] - v0[0], v2[1] - v0[1],
                                     dimension == 3 ? v2[2] - v0[2] : 0.0};
    const std::array<double, 3> q = {point[0] - v0[0], point[1] - v0[1],
                                     dimension == 3 ? point[2] - v0[2] : 0.0};
    if (dimension == 2)
    {
      const double determinant = a[0] * b[1] - a[1] * b[0];
      if (determinant == 0.0)
      {
        return false;
      }
      reference.Set2((q[0] * b[1] - q[1] * b[0]) / determinant,
                     (a[0] * q[1] - a[1] * q[0]) / determinant);
    }
    else
    {
      const double *v3 = mesh.GetVertex(vertices[3]);
      const std::array<double, 3> c = {v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2]};
      auto Determinant = [](const std::array<double, 3> &x, const std::array<double, 3> &y,
                            const std::array<double, 3> &z)
      {
        return x[0] * (y[1] * z[2] - y[2] * z[1]) - x[1] * (y[0] * z[2] - y[2] * z[0]) +
               x[2] * (y[0] * z[1] - y[1] * z[0]);
      };
      const double determinant = Determinant(a, b, c);
      if (determinant == 0.0)
      {
        return false;
      }
      reference.Set3(Determinant(q, b, c) / determinant, Determinant(a, q, c) / determinant,
                     Determinant(a, b, q) / determinant);
    }
    return true;
  }

  bool FindInLinearSimplex(int element, const std::array<double, 3> &point,
                           mfem::IntegrationPoint &reference) const
  {
    return GetLinearSimplexReference(element, point, reference) &&
           mfem::Geometry::CheckPoint(mesh.GetElement(element)->GetGeometryType(),
                                      reference, 1.0e-9);
  }

public:
  ElementPointLocator(mfem::ParMesh &mesh_, int dimension_)
    : mesh(mesh_), dimension(dimension_),
      linear_mesh(!mesh.GetNodes() ||
                  mesh.GetNodes()->FESpace()->GetMaxElementOrder() == 1),
      element_boxes(mesh.GetNE()), indices(mesh.GetNE())
  {
    MFEM_VERIFY(mesh.SpaceDimension() == dimension,
                "Surface-response point coordinates do not match the mesh dimension!");
    MFEM_VERIFY(mesh.GetNE() > 0, "Cannot locate points in an empty local mesh!");
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      element_boxes[element] = GetElementBox(element);
    }
    std::iota(indices.begin(), indices.end(), 0);
    nodes.reserve(2 * (indices.size() / leaf_size + 1));
    Build(0, indices.size());
  }

  const ElementBox &GetBounds() const { return nodes.front().box; }

  bool SupportsExactSegmentIntersections() const
  {
    if (!linear_mesh)
    {
      return false;
    }
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      const auto geometry = mesh.GetElement(element)->GetGeometryType();
      if ((dimension == 2 && geometry != mfem::Geometry::TRIANGLE) ||
          (dimension == 3 && geometry != mfem::Geometry::TETRAHEDRON))
      {
        return false;
      }
    }
    return true;
  }

  std::vector<ElementBox> GetRoutingBoxes(std::size_t count) const
  {
    std::vector<int> frontier = {0};
    while (frontier.size() < count)
    {
      auto split = std::max_element(frontier.begin(), frontier.end(),
                                    [this](int a, int b)
                                    {
                                      const auto &left = nodes[a];
                                      const auto &right = nodes[b];
                                      const std::size_t left_size =
                                          left.IsLeaf() ? 0 : left.end - left.begin;
                                      const std::size_t right_size =
                                          right.IsLeaf() ? 0 : right.end - right.begin;
                                      return left_size < right_size;
                                    });
      if (split == frontier.end() || nodes[*split].IsLeaf())
      {
        break;
      }
      const int node = *split;
      *split = nodes[node].left;
      frontier.push_back(nodes[node].right);
    }
    std::vector<ElementBox> boxes;
    boxes.reserve(frontier.size());
    for (const int node : frontier)
    {
      boxes.push_back(nodes[node].box);
    }
    return boxes;
  }

  bool Find(const std::array<double, 3> &point, double tolerance, int &element,
            mfem::IntegrationPoint &reference, std::vector<int> &candidates)
  {
    candidates.clear();
    FindCandidates(0, point, tolerance, candidates);
    std::sort(candidates.begin(), candidates.end());

    mfem::Vector physical(dimension);
    std::copy_n(point.data(), dimension, physical.HostWrite());
    mfem::InverseElementTransformation inverse;
    inverse.SetReferenceTol(1.0e-12);
    inverse.SetPhysicalRelTol(1.0e-12);
    for (const int candidate : candidates)
    {
      const auto geometry = mesh.GetElement(candidate)->GetGeometryType();
      const bool linear_simplex =
          linear_mesh && ((dimension == 2 && geometry == mfem::Geometry::TRIANGLE) ||
                          (dimension == 3 && geometry == mfem::Geometry::TETRAHEDRON));
      if (linear_simplex)
      {
        if (FindInLinearSimplex(candidate, point, reference))
        {
          element = candidate;
          return true;
        }
        continue;
      }
      auto &transformation = *mesh.GetElementTransformation(candidate);
      inverse.SetTransformation(transformation);
      inverse.SetInitialGuessType(transformation.Order() > 1
                                      ? mfem::InverseElementTransformation::ClosestPhysNode
                                      : mfem::InverseElementTransformation::Center);
      if (inverse.Transform(physical, reference) ==
          mfem::InverseElementTransformation::Inside)
      {
        element = candidate;
        return true;
      }
      if (transformation.Order() > 1)
      {
        inverse.SetInitialGuessType(mfem::InverseElementTransformation::ClosestRefNode);
        if (inverse.Transform(physical, reference) ==
            mfem::InverseElementTransformation::Inside)
        {
          element = candidate;
          return true;
        }
      }
    }
    return false;
  }

  void FindSegmentIntersections(const std::array<double, 3> &p0,
                                const std::array<double, 3> &p1, double tolerance,
                                std::vector<std::pair<double, double>> &intervals,
                                std::vector<int> &candidates) const
  {
    MFEM_ASSERT(linear_mesh, "Exact segment intersections require a linear mesh!");
    candidates.clear();
    FindSegmentCandidates(0, p0, p1, tolerance, candidates);
    std::sort(candidates.begin(), candidates.end());
    intervals.clear();
    constexpr double reference_tolerance = 64.0 * std::numeric_limits<double>::epsilon();
    for (const int element : candidates)
    {
      mfem::IntegrationPoint first, second;
      if (!GetLinearSimplexReference(element, p0, first) ||
          !GetLinearSimplexReference(element, p1, second))
      {
        continue;
      }
      std::array<double, 4> lambda_first = {1.0 - first.x - first.y -
                                                (dimension == 3 ? first.z : 0.0),
                                            first.x, first.y, first.z};
      std::array<double, 4> lambda_second = {1.0 - second.x - second.y -
                                                 (dimension == 3 ? second.z : 0.0),
                                             second.x, second.y, second.z};
      double begin = 0.0, end = 1.0;
      for (int i = 0; i < dimension + 1; i++)
      {
        const double value = lambda_first[i];
        const double slope = lambda_second[i] - value;
        if (slope == 0.0)
        {
          if (value < -reference_tolerance)
          {
            begin = 1.0;
            end = 0.0;
            break;
          }
          continue;
        }
        const double crossing = (-reference_tolerance - value) / slope;
        if (slope > 0.0)
        {
          begin = std::max(begin, crossing);
        }
        else
        {
          end = std::min(end, crossing);
        }
      }
      begin = std::clamp(begin, 0.0, 1.0);
      end = std::clamp(end, 0.0, 1.0);
      if (end > begin)
      {
        intervals.emplace_back(begin, end);
      }
    }
  }
};

enum class LibraryTopology : char
{
  ISOLATED_EDGE,
  SAME_CONDUCTOR_GAP,
  DIFFERENT_CONDUCTOR_GAP,
  SAME_CONDUCTOR_STRIP,
  PARALLEL_EDGE_CLUSTER,
  SPATIAL_EDGE_CLUSTER,
  CONVEX_CORNER,
  CONCAVE_CORNER,
  ENDPOINT,
  JUNCTION
};

struct LibraryInterface
{
  int slot = 0;
  InterfaceDielectric type;
  int coupon;
};

struct LibraryInterfaceLayer
{
  double thickness = 0.0;
  double permittivity = 0.0;
};

struct LibraryClusterEdge
{
  double offset = 0.0;
  int gap_direction = 0;
  int conductor = 0;
};

struct MetalBoundaryLaw
{
  MetalBoundaryConditionType type = MetalBoundaryConditionType::PEC;
  bool parameters_verified = true;
  std::vector<double> parameters;
  std::vector<double> numerator;
  std::vector<double> denominator;
};

struct LibrarySpatialEdge
{
  Point3D point{};
  Point3D gap_direction{};
  Point3D process_normal{};
  std::array<double, 2> interval{};
  int conductor = 0;
  int interface_slot = 0;
  MetalBoundaryLaw boundary_condition;
};

struct LibraryModel
{
  std::string name;
  LibraryTopology topology;
  double separation = 0.0;
  double separation_tolerance = 0.0;
  double angle = 0.0;
  double angle_tolerance = 0.0;
  double corner_radius = 0.0;
  double corner_radius_tolerance = 0.0;
  std::vector<double> arm_angles;
  double arm_angle_tolerance = 0.0;
  MetalBoundaryLaw boundary_condition;
  double coupon_depth = 0.0;
  ResponseModelData response;
  std::vector<std::array<double, 3>> conductor_references;
  std::vector<LibraryClusterEdge> cluster_edges;
  double cluster_offset_tolerance = 0.0;
  std::vector<LibrarySpatialEdge> spatial_edges;
  double spatial_position_tolerance = 0.0;
  double spatial_angle_tolerance = 0.0;
  bool boundary_law_physics_qualified = true;
  std::optional<std::string> plan_view_boundary;
  std::optional<std::string> mask_regularization;
  std::vector<LibraryInterface> interfaces;
};

struct ProcessLibrary
{
  int version = 0;
  std::string name;
  double matching_radius = 0.0;
  std::map<InterfaceDielectric, LibraryInterfaceLayer> interface_layers;
  std::vector<LibraryModel> models;
  std::set<std::pair<std::size_t, std::size_t>> corner_radius_interpolation;
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
  MetalBoundaryLaw boundary_condition;
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
  int metal_component = -1;
  std::map<InterfaceDielectric, int> targets;
  MetalBoundaryLaw boundary_condition;
};

// Deterministic broad-phase index for nearby segment queries. Exact distance and topology
// predicates remain at the call sites; this index only removes pairs whose expanded AABBs
// cannot interact.
class SegmentBoxIndex
{
private:
  static constexpr std::size_t leaf_size = 16;
  struct Box
  {
    Point3D min{};
    Point3D max{};
  };
  struct Node
  {
    Box box;
    std::size_t begin = 0;
    std::size_t end = 0;
    int left = -1;
    int right = -1;
    bool IsLeaf() const { return left < 0; }
  };

  std::vector<Box> boxes;
  std::vector<std::size_t> indices;
  std::vector<Node> nodes;

  static bool Intersects(const Box &a, const Box &b)
  {
    for (int d = 0; d < 3; d++)
    {
      if (a.max[d] < b.min[d] || b.max[d] < a.min[d])
      {
        return false;
      }
    }
    return true;
  }

  int Build(std::size_t begin, std::size_t end)
  {
    Node node;
    node.begin = begin;
    node.end = end;
    node.box.min.fill(mfem::infinity());
    node.box.max.fill(-mfem::infinity());
    for (std::size_t i = begin; i < end; i++)
    {
      const auto &box = boxes[indices[i]];
      for (int d = 0; d < 3; d++)
      {
        node.box.min[d] = std::min(node.box.min[d], box.min[d]);
        node.box.max[d] = std::max(node.box.max[d], box.max[d]);
      }
    }
    const int node_index = static_cast<int>(nodes.size());
    nodes.push_back(node);
    if (end - begin <= leaf_size)
    {
      return node_index;
    }
    int axis = 0;
    for (int d = 1; d < 3; d++)
    {
      if (node.box.max[d] - node.box.min[d] > node.box.max[axis] - node.box.min[axis])
      {
        axis = d;
      }
    }
    const std::size_t middle = begin + (end - begin) / 2;
    std::nth_element(indices.begin() + begin, indices.begin() + middle,
                     indices.begin() + end,
                     [&](std::size_t first, std::size_t second)
                     {
                       return boxes[first].min[axis] + boxes[first].max[axis] <
                                  boxes[second].min[axis] + boxes[second].max[axis] ||
                              (boxes[first].min[axis] + boxes[first].max[axis] ==
                                   boxes[second].min[axis] + boxes[second].max[axis] &&
                               first < second);
                     });
    nodes[node_index].left = Build(begin, middle);
    nodes[node_index].right = Build(middle, end);
    return node_index;
  }

  void Query(int node_index, const Box &query, std::vector<std::size_t> &matches) const
  {
    const auto &node = nodes[node_index];
    if (!Intersects(node.box, query))
    {
      return;
    }
    if (node.IsLeaf())
    {
      for (std::size_t i = node.begin; i < node.end; i++)
      {
        const std::size_t index = indices[i];
        if (Intersects(boxes[index], query))
        {
          matches.push_back(index);
        }
      }
      return;
    }
    Query(node.left, query, matches);
    Query(node.right, query, matches);
  }

public:
  explicit SegmentBoxIndex(const std::vector<std::pair<Point3D, Point3D>> &segments)
    : boxes(segments.size()), indices(segments.size())
  {
    MFEM_VERIFY(!segments.empty(), "Cannot build an empty segment proximity index!");
    for (std::size_t i = 0; i < segments.size(); i++)
    {
      for (int d = 0; d < 3; d++)
      {
        boxes[i].min[d] = std::min(segments[i].first[d], segments[i].second[d]);
        boxes[i].max[d] = std::max(segments[i].first[d], segments[i].second[d]);
      }
    }
    std::iota(indices.begin(), indices.end(), 0);
    nodes.reserve(2 * (segments.size() / leaf_size + 1));
    Build(0, segments.size());
  }

  std::vector<std::size_t> Query(const Point3D &p0, const Point3D &p1,
                                 double distance) const
  {
    Box query;
    for (int d = 0; d < 3; d++)
    {
      query.min[d] = std::min(p0[d], p1[d]) - distance;
      query.max[d] = std::max(p0[d], p1[d]) + distance;
    }
    std::vector<std::size_t> matches;
    Query(0, query, matches);
    std::sort(matches.begin(), matches.end());
    return matches;
  }

  std::vector<std::pair<std::size_t, std::size_t>> CandidatePairs(double distance) const
  {
    std::vector<std::pair<std::size_t, std::size_t>> result;
    for (std::size_t first = 0; first < boxes.size(); first++)
    {
      Point3D p0 = boxes[first].min;
      Point3D p1 = boxes[first].max;
      for (const std::size_t second : Query(p0, p1, distance))
      {
        if (second > first)
        {
          result.emplace_back(first, second);
        }
      }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
  }
};

std::vector<std::pair<Point3D, Point3D>>
SegmentGeometry(const std::vector<EdgeSegment3D> &segments)
{
  std::vector<std::pair<Point3D, Point3D>> result;
  result.reserve(segments.size());
  for (const auto &segment : segments)
  {
    result.emplace_back(segment.p0, segment.p1);
  }
  return result;
}

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

struct LibrarySelection
{
  struct WeightedModel
  {
    std::size_t index = 0;
    double weight = 1.0;
  };

  std::vector<WeightedModel> models;
  std::vector<std::array<double, 3>> conductor_references;
  double normalized_distance = 0.0;

  bool IsInterpolated() const { return models.size() > 1; }
};

struct ParallelClusterSelection
{
  LibrarySelection response;
  std::vector<std::size_t> ordered_edges;
  std::vector<std::size_t> reference_edges;
  Point2D axis_u{};
  Point2D axis_v{};
};

struct ParallelClusterSelection3D
{
  LibrarySelection response;
  std::vector<std::size_t> ordered_edges;
  std::vector<std::size_t> reference_edges;
  Point3D axis_u{};
  Point3D axis_v{};
};

struct ParallelClusterSpan3D
{
  ParallelClusterSelection3D selection;
  Point3D tangent{};
  double begin = 0.0;
  double end = 0.0;
};

struct UnmatchedParallelClusterSpan3D
{
  std::vector<std::size_t> edges;
  Point3D tangent{};
  double begin = 0.0;
  double end = 0.0;
};

struct ParallelClusterSpans3D
{
  std::vector<ParallelClusterSpan3D> matched;
  std::vector<UnmatchedParallelClusterSpan3D> unmatched;
};

struct SpatialEdgeSite3D
{
  int physical_chain = -1;
  std::size_t geometry_index = 0;
  std::size_t segment = 0;
  double distance = 0.0;
  std::array<double, 2> interval{};
  Point3D point{};
  Point3D gap_direction{};
  Point3D process_normal{};
  int conductor = std::numeric_limits<int>::min();
  int metal_component = -1;
  std::map<InterfaceDielectric, int> targets;
  MetalBoundaryLaw boundary_condition;
};

struct SpatialClusterSelection3D
{
  struct InteractionNeighborhood
  {
    int first_chain = -1;
    int second_chain = -1;
    Point3D center{};
  };

  LibrarySelection response;
  std::vector<SpatialEdgeSite3D> sites;
  std::vector<std::size_t> model_to_site;
  std::map<int, std::map<InterfaceDielectric, int>> targets_by_slot;
  std::vector<InteractionNeighborhood> interactions;
  Point3D origin{};
  std::array<Point3D, 3> axes{};
};

struct AutomaticResponseDiagnostics
{
  double matching_radius = 0.0;
  double minimum_wave_speed = mfem::infinity();
  double selected_length = 0.0;
  double matched_length = 0.0;
  double matched_corner_neighborhood_length = 0.0;
  std::map<int, double> selected_length_by_interface;
  std::map<int, double> matched_length_by_interface;
  std::map<int, double> matched_corner_neighborhood_length_by_interface;
  double maximum_curvature_ratio = 0.0;
  double maximum_library_distance = 0.0;
  bool boundary_law_verified = true;
};

struct AutomaticResponseStatistics
{
  long long int target_groups = 0;
  long long int edge_sites_2d = 0;
  long long int metal_vertices = 0;
  long long int metal_segments = 0;
  long long int metal_components = 0;
  long long int physical_components = 0;
  long long int physical_chains = 0;
  long long int surface_faces_local = 0;

  long long int pair_checks_global_spatial = 0;
  long long int pair_checks_external_conflict = 0;
  long long int pair_checks_group_spatial = 0;
  long long int pair_checks_safety = 0;
  long long int pair_checks_patch_construction = 0;
  long long int spatial_events = 0;

  long long int mask_gather_calls = 0;
  long long int mask_faces_scanned_local = 0;
  long long int mask_facets_packed_local = 0;
  long long int mask_payload_scalars_local = 0;
  long long int mask_gathered_scalars = 0;
};

long long int PairCount(std::size_t count)
{
  return count > 1 ? static_cast<long long int>(count) * (count - 1) / 2 : 0;
}

nlohmann::json BuildAutomaticStatistics(MPI_Comm comm,
                                        const AutomaticResponseStatistics &statistics)
{
  auto Replicated = [comm](long long int value, const char *name)
  {
    long long int minimum = value;
    long long int maximum = value;
    Mpi::GlobalMin(1, &minimum, comm);
    Mpi::GlobalMax(1, &maximum, comm);
    MFEM_VERIFY(minimum == maximum,
                "Rank-inconsistent surface-response statistic \"" << name << "\"!");
    return minimum;
  };
  auto Distribution = [comm](long long int value)
  {
    long long int total = value;
    long long int minimum = value;
    long long int maximum = value;
    int nonzero = value > 0 ? 1 : 0;
    Mpi::GlobalSum(1, &total, comm);
    Mpi::GlobalMin(1, &minimum, comm);
    Mpi::GlobalMax(1, &maximum, comm);
    Mpi::GlobalSum(1, &nonzero, comm);
    return nlohmann::json{{"Total", total},
                          {"Minimum", minimum},
                          {"Maximum", maximum},
                          {"NonzeroRanks", nonzero}};
  };

  const long long int global_pair_checks =
      statistics.pair_checks_global_spatial + statistics.pair_checks_external_conflict +
      statistics.pair_checks_group_spatial + statistics.pair_checks_safety +
      statistics.pair_checks_patch_construction;
  return {{"Version", 1},
          {"Geometry",
           {{"TargetGroups", Replicated(statistics.target_groups, "TargetGroups")},
            {"EdgeSites2D", Replicated(statistics.edge_sites_2d, "EdgeSites2D")},
            {"MetalVertices", Replicated(statistics.metal_vertices, "MetalVertices")},
            {"MetalSegments", Replicated(statistics.metal_segments, "MetalSegments")},
            {"MetalComponents", Replicated(statistics.metal_components, "MetalComponents")},
            {"PhysicalComponents",
             Replicated(statistics.physical_components, "PhysicalComponents")},
            {"PhysicalChains", Replicated(statistics.physical_chains, "PhysicalChains")},
            {"SurfaceFaces", Distribution(statistics.surface_faces_local)}}},
          {"Matching",
           {{"PairChecks",
             {{"GlobalSpatial", Replicated(statistics.pair_checks_global_spatial,
                                           "GlobalSpatialPairChecks")},
              {"ExternalConflict", Replicated(statistics.pair_checks_external_conflict,
                                              "ExternalConflictPairChecks")},
              {"GroupSpatial",
               Replicated(statistics.pair_checks_group_spatial, "GroupSpatialPairChecks")},
              {"SafetyClassification",
               Replicated(statistics.pair_checks_safety, "SafetyPairChecks")},
              {"PatchConstruction", Replicated(statistics.pair_checks_patch_construction,
                                               "PatchConstructionPairChecks")},
              {"Total", Replicated(global_pair_checks, "TotalPairChecks")}}},
            {"SpatialEvents", Replicated(statistics.spatial_events, "SpatialEvents")}}},
          {"Masks",
           {{"GatherCalls", Replicated(statistics.mask_gather_calls, "MaskGatherCalls")},
            {"FacesScanned", Distribution(statistics.mask_faces_scanned_local)},
            {"FacetsPacked", Distribution(statistics.mask_facets_packed_local)},
            {"PayloadScalars", Distribution(statistics.mask_payload_scalars_local)},
            {"GatheredScalars",
             Replicated(statistics.mask_gathered_scalars, "MaskGatheredScalars")}}}};
}

double Dot(const Point3D &a, const Point3D &b);
double Norm(const Point3D &a);

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
  if (topology == "ParallelEdgeCluster")
  {
    return LibraryTopology::PARALLEL_EDGE_CLUSTER;
  }
  if (topology == "SpatialEdgeCluster")
  {
    return LibraryTopology::SPATIAL_EDGE_CLUSTER;
  }
  if (topology == "ConvexCorner")
  {
    return LibraryTopology::CONVEX_CORNER;
  }
  if (topology == "ConcaveCorner")
  {
    return LibraryTopology::CONCAVE_CORNER;
  }
  if (topology == "Endpoint")
  {
    return LibraryTopology::ENDPOINT;
  }
  if (topology == "Junction")
  {
    return LibraryTopology::JUNCTION;
  }
  MFEM_ABORT("Unknown fabrication-process response topology \"" << topology << "\"!");
}

MetalBoundaryConditionType ParseLibraryBoundaryConditionType(const std::string &condition)
{
  if (condition == "PEC")
  {
    return MetalBoundaryConditionType::PEC;
  }
  if (condition == "Conductivity")
  {
    return MetalBoundaryConditionType::CONDUCTIVITY;
  }
  if (condition == "Impedance")
  {
    return MetalBoundaryConditionType::IMPEDANCE;
  }
  if (condition == "RationalImpedance")
  {
    return MetalBoundaryConditionType::RATIONAL_IMPEDANCE;
  }
  MFEM_ABORT("Unknown fabrication-process response boundary condition \""
             << condition
             << "\"; expected PEC, Conductivity, Impedance, or RationalImpedance!");
}

void NormalizeRationalLaw(MetalBoundaryLaw &law)
{
  auto TrimLeadingZeros = [](std::vector<double> &coefficients)
  {
    auto first = std::find_if(coefficients.begin(), coefficients.end(),
                              [](double value) { return value != 0.0; });
    coefficients.erase(coefficients.begin(), first);
  };
  TrimLeadingZeros(law.numerator);
  TrimLeadingZeros(law.denominator);
  MFEM_VERIFY(!law.numerator.empty() && !law.denominator.empty(),
              "Rational-impedance boundary-law polynomials must be nonzero!");
  const double scale = law.denominator.front();
  for (double &coefficient : law.numerator)
  {
    coefficient /= scale;
  }
  for (double &coefficient : law.denominator)
  {
    coefficient /= scale;
  }
}

MetalBoundaryLaw ParseLibraryBoundaryCondition(const nlohmann::json &condition,
                                               const Units &units, bool nondimensionalize)
{
  MetalBoundaryLaw law;
  if (condition.is_string())
  {
    law.type = ParseLibraryBoundaryConditionType(condition.get<std::string>());
    law.parameters_verified = law.type == MetalBoundaryConditionType::PEC;
    return law;
  }

  MFEM_VERIFY(condition.is_object() && condition.contains("Type"),
              "Fabrication-process response BoundaryCondition must be a string or an "
              "object containing Type!");
  law.type = ParseLibraryBoundaryConditionType(condition.at("Type").get<std::string>());
  const auto VerifyKeys = [&condition](std::initializer_list<const char *> allowed)
  {
    for (const auto &item : condition.items())
    {
      const auto &key = item.key();
      const bool valid =
          std::any_of(allowed.begin(), allowed.end(),
                      [&key](const char *candidate) { return key == candidate; });
      MFEM_VERIFY(valid, "Unknown fabrication-process response BoundaryCondition key \""
                             << key << "\"!");
    }
  };
  nlohmann::json parameters = condition;
  parameters.erase("Type");
  parameters["Attributes"] = {1};
  switch (law.type)
  {
    case MetalBoundaryConditionType::PEC:
      VerifyKeys({"Type"});
      MFEM_VERIFY(parameters.size() == 1,
                  "A PEC response BoundaryCondition cannot specify parameters!");
      break;
    case MetalBoundaryConditionType::CONDUCTIVITY:
      {
        VerifyKeys({"Type", "Conductivity", "Permeability", "Thickness", "External"});
        config::ConductivityData data(parameters);
        if (nondimensionalize)
        {
          config::Nondimensionalize(units, data);
        }
        MFEM_VERIFY(std::isfinite(data.sigma) && data.sigma > 0.0 &&
                        std::isfinite(data.mu_r) && data.mu_r > 0.0 &&
                        std::isfinite(data.h) && data.h >= 0.0,
                    "Conductivity response BoundaryCondition parameters are invalid!");
        law.parameters = {data.sigma, data.mu_r, data.external ? 2.0 * data.h : data.h};
        break;
      }
    case MetalBoundaryConditionType::IMPEDANCE:
      {
        VerifyKeys({"Type", "Rs", "Ls", "Cs"});
        config::ImpedanceData data(parameters);
        if (nondimensionalize)
        {
          config::Nondimensionalize(units, data);
        }
        MFEM_VERIFY(std::isfinite(data.Rs) && std::isfinite(data.Ls) &&
                        std::isfinite(data.Cs) &&
                        std::abs(data.Rs) + std::abs(data.Ls) + std::abs(data.Cs) > 0.0,
                    "Impedance response BoundaryCondition parameters are invalid!");
        law.parameters = {data.Rs, data.Ls, data.Cs};
        break;
      }
    case MetalBoundaryConditionType::RATIONAL_IMPEDANCE:
      {
        VerifyKeys({"Type", "Numerator", "Denominator"});
        config::RationalImpedanceData data(parameters);
        if (nondimensionalize)
        {
          config::Nondimensionalize(units, data);
        }
        law.numerator = std::move(data.num);
        law.denominator = std::move(data.den);
        NormalizeRationalLaw(law);
        break;
      }
  }
  return law;
}

bool CompatibleBoundaryLaw(const MetalBoundaryLaw &library, const MetalBoundaryLaw &actual)
{
  if (library.type != actual.type)
  {
    return false;
  }
  if (!library.parameters_verified)
  {
    return true;
  }
  const auto Compatible =
      [](const std::vector<double> &first, const std::vector<double> &second)
  {
    if (first.size() != second.size())
    {
      return false;
    }
    for (std::size_t i = 0; i < first.size(); i++)
    {
      const double scale = std::max({std::abs(first[i]), std::abs(second[i]), 1.0e-300});
      if (std::abs(first[i] - second[i]) > 1.0e-10 * scale)
      {
        return false;
      }
    }
    return true;
  };
  return Compatible(library.parameters, actual.parameters) &&
         Compatible(library.numerator, actual.numerator) &&
         Compatible(library.denominator, actual.denominator);
}

bool SameBoundaryLaw(const MetalBoundaryLaw &first, const MetalBoundaryLaw &second)
{
  return first.parameters_verified && second.parameters_verified &&
         CompatibleBoundaryLaw(first, second) && CompatibleBoundaryLaw(second, first);
}

bool IsBoundaryLawVerified(const LibraryModel &model)
{
  if (!model.boundary_law_physics_qualified)
  {
    return false;
  }
  if (model.topology == LibraryTopology::SPATIAL_EDGE_CLUSTER)
  {
    return std::all_of(model.spatial_edges.begin(), model.spatial_edges.end(),
                       [](const auto &edge)
                       { return edge.boundary_condition.parameters_verified; });
  }
  return model.boundary_condition.parameters_verified;
}

ProcessLibrary ReadProcessLibrary(const std::string &path, const Units &units,
                                  bool nondimensionalize, bool allow_empty_models = false)
{
  std::ifstream input(path);
  MFEM_VERIFY(input,
              "Unable to open fabrication-process response library \"" << path << "\"!");
  nlohmann::json data;
  input >> data;
  const int version = data.value("Version", 0);
  MFEM_VERIFY(version == 1 || version == 2 || version == 3,
              "Fabrication-process response library \"" << path
                                                        << "\" has unsupported version!");

  ProcessLibrary library;
  library.version = version;
  library.name = data.value("Name", std::filesystem::path(path).stem().string());
  const double coordinate_scale = units.GetMeshLengthRelativeScale();
  library.matching_radius = data.at("MatchingRadius").get<double>() / coordinate_scale;
  MFEM_VERIFY(std::isfinite(library.matching_radius) && library.matching_radius > 0.0,
              "Fabrication-process response-library matching radius must be positive!");
  const auto fabrication = data.find("Fabrication");
  const nlohmann::json *interface_layers = nullptr;
  if (fabrication != data.end())
  {
    MFEM_VERIFY(fabrication->is_object(),
                "Fabrication-process response-library Fabrication metadata must be an "
                "object!");
    const auto entries = fabrication->find("InterfaceLayers");
    if (entries != fabrication->end())
    {
      interface_layers = &*entries;
    }
  }
  MFEM_VERIFY(version < 3 || (interface_layers && interface_layers->is_object()),
              "Version-3 fabrication-process response libraries require "
              "Fabrication.InterfaceLayers metadata!");
  if (interface_layers)
  {
    MFEM_VERIFY(interface_layers->is_object(),
                "Fabrication-process response-library InterfaceLayers must be an object!");
    for (auto entry = interface_layers->begin(); entry != interface_layers->end(); ++entry)
    {
      InterfaceDielectric type = InterfaceDielectric::DEFAULT;
      FromString(entry.key(), type);
      MFEM_VERIFY(type != InterfaceDielectric::DEFAULT && entry.value().is_object(),
                  "Fabrication-process response-library InterfaceLayers entries must use "
                  "the explicit types MA, MS, or SA and contain layer properties!");
      LibraryInterfaceLayer layer;
      layer.thickness = entry.value().at("Thickness").get<double>() / coordinate_scale;
      layer.permittivity = entry.value().at("Permittivity").get<double>();
      MFEM_VERIFY(std::isfinite(layer.thickness) && layer.thickness > 0.0 &&
                      std::isfinite(layer.permittivity) && layer.permittivity > 0.0,
                  "Fabrication-process response-library InterfaceLayers thicknesses and "
                  "permittivities must be finite and positive!");
      MFEM_VERIFY(library.interface_layers.emplace(type, layer).second,
                  "Fabrication-process response-library InterfaceLayers types must be "
                  "unique!");
    }
  }

  const auto directory = std::filesystem::path(path).parent_path();
  double default_coupon_depth = 0.0;
  if (auto depth = data.find("CouponDepth"); depth != data.end())
  {
    default_coupon_depth = depth->get<double>() / coordinate_scale;
    MFEM_VERIFY(std::isfinite(default_coupon_depth) && default_coupon_depth > 0.0,
                "Fabrication-process response-library CouponDepth must be positive!");
  }
  const auto &models = data.at("Models");
  MFEM_VERIFY(models.is_array(),
              "Fabrication-process response-library Models must be an array!");
  MFEM_VERIFY(allow_empty_models || !models.empty(),
              "Fabrication-process response library must contain at least one model!");
  std::set<std::string> names;
  std::set<InterfaceDielectric> mapped_interface_types;
  for (const auto &entry : models)
  {
    LibraryModel model;
    model.name = entry.at("Name").get<std::string>();
    MFEM_VERIFY(!model.name.empty() && names.insert(model.name).second,
                "Fabrication-process response model names must be nonempty and unique!");
    model.topology = ParseLibraryTopology(entry.at("Topology").get<std::string>());
    model.separation = entry.value("Separation", 0.0) / coordinate_scale;
    model.separation_tolerance = entry.value("SeparationTolerance", 0.0) / coordinate_scale;
    model.angle = entry.value("Angle", 0.0) * std::acos(-1.0) / 180.0;
    model.angle_tolerance = entry.value("AngleTolerance", 0.0) * std::acos(-1.0) / 180.0;
    model.corner_radius = entry.value("CornerRadius", 0.0) / coordinate_scale;
    model.corner_radius_tolerance =
        entry.value("CornerRadiusTolerance", 0.0) / coordinate_scale;
    model.arm_angles = entry.value("ArmAngles", std::vector<double>{});
    model.arm_angle_tolerance =
        entry.value("ArmAngleTolerance", 0.0) * std::acos(-1.0) / 180.0;
    for (double &angle : model.arm_angles)
    {
      angle *= std::acos(-1.0) / 180.0;
    }
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
    const bool parallel_cluster = model.topology == LibraryTopology::PARALLEL_EDGE_CLUSTER;
    const bool spatial_cluster = model.topology == LibraryTopology::SPATIAL_EDGE_CLUSTER;
    if (auto edges = entry.find("Edges"); edges != entry.end())
    {
      MFEM_VERIFY(
          (parallel_cluster && edges->is_array() && edges->size() >= 3) ||
              (spatial_cluster && edges->is_array() && edges->size() >= 2),
          "Edges requires a ParallelEdgeCluster model with at least three entries or a "
          "SpatialEdgeCluster model with at least two entries!");
      if (parallel_cluster)
      {
        model.cluster_offset_tolerance =
            entry.value("EdgeOffsetTolerance", 0.0) / coordinate_scale;
        MFEM_VERIFY(std::isfinite(model.cluster_offset_tolerance) &&
                        model.cluster_offset_tolerance >= 0.0,
                    "ParallelEdgeCluster EdgeOffsetTolerance must be nonnegative!");
        int next_conductor = 1;
        std::set<int> conductors;
        for (const auto &edge : *edges)
        {
          LibraryClusterEdge cluster_edge;
          cluster_edge.offset = edge.at("Offset").get<double>() / coordinate_scale;
          cluster_edge.gap_direction = edge.at("GapDirection").get<int>();
          cluster_edge.conductor = edge.at("Conductor").get<int>();
          MFEM_VERIFY(
              std::isfinite(cluster_edge.offset) &&
                  (cluster_edge.gap_direction == -1 || cluster_edge.gap_direction == 1) &&
                  cluster_edge.conductor > 0,
              "ParallelEdgeCluster edges require a finite Offset, GapDirection "
              "equal to -1 or 1, and a positive Conductor!");
          if (conductors.insert(cluster_edge.conductor).second)
          {
            MFEM_VERIFY(cluster_edge.conductor == next_conductor++,
                        "ParallelEdgeCluster conductor labels must be canonical and "
                        "contiguous in order of first occurrence!");
          }
          model.cluster_edges.push_back(cluster_edge);
        }
        MFEM_VERIFY(std::abs(model.cluster_edges.front().offset) <=
                            1.0e-12 * library.matching_radius &&
                        std::adjacent_find(model.cluster_edges.begin(),
                                           model.cluster_edges.end(),
                                           [](const auto &first, const auto &second)
                                           { return first.offset >= second.offset; }) ==
                            model.cluster_edges.end(),
                    "ParallelEdgeCluster edge offsets must begin at zero and be strictly "
                    "increasing!");
      }
      else
      {
        model.spatial_position_tolerance =
            entry.value("EdgePositionTolerance", 0.0) / coordinate_scale;
        model.spatial_angle_tolerance =
            entry.value("EdgeAngleTolerance", 0.0) * std::acos(-1.0) / 180.0;
        MFEM_VERIFY(std::isfinite(model.spatial_position_tolerance) &&
                        model.spatial_position_tolerance >= 0.0 &&
                        std::isfinite(model.spatial_angle_tolerance) &&
                        model.spatial_angle_tolerance >= 0.0,
                    "SpatialEdgeCluster edge position and angle tolerances must be "
                    "nonnegative!");
        int next_conductor = 1;
        std::set<int> conductors;
        for (const auto &edge : *edges)
        {
          LibrarySpatialEdge spatial_edge;
          spatial_edge.point = edge.at("Point").get<Point3D>();
          spatial_edge.gap_direction = edge.at("GapDirection").get<Point3D>();
          spatial_edge.process_normal = edge.at("ProcessNormal").get<Point3D>();
          spatial_edge.interval = edge.at("Interval").get<std::array<double, 2>>();
          spatial_edge.conductor = edge.at("Conductor").get<int>();
          spatial_edge.interface_slot = edge.value("InterfaceSlot", 0);
          spatial_edge.boundary_condition = ParseLibraryBoundaryCondition(
              edge.value("BoundaryCondition", nlohmann::json("PEC")), units,
              nondimensionalize);
          for (double &value : spatial_edge.point)
          {
            value /= coordinate_scale;
          }
          for (double &value : spatial_edge.interval)
          {
            value /= coordinate_scale;
          }
          const double gap_norm = Norm(spatial_edge.gap_direction);
          const double normal_norm = Norm(spatial_edge.process_normal);
          MFEM_VERIFY(
              std::all_of(spatial_edge.point.begin(), spatial_edge.point.end(),
                          [](double value) { return std::isfinite(value); }) &&
                  std::isfinite(gap_norm) && std::isfinite(normal_norm) &&
                  std::abs(gap_norm - 1.0) <= 1.0e-10 &&
                  std::abs(normal_norm - 1.0) <= 1.0e-10 &&
                  std::abs(Dot(spatial_edge.gap_direction, spatial_edge.process_normal)) <=
                      1.0e-10 &&
                  std::isfinite(spatial_edge.interval[0]) &&
                  std::isfinite(spatial_edge.interval[1]) &&
                  spatial_edge.interval[0] <= 0.0 && spatial_edge.interval[1] >= 0.0 &&
                  spatial_edge.interval[1] > spatial_edge.interval[0] &&
                  spatial_edge.conductor > 0 && spatial_edge.interface_slot >= 0,
              "SpatialEdgeCluster edges require a finite Point, orthonormal unit "
              "GapDirection and ProcessNormal vectors, an Interval containing zero, "
              "a positive Conductor, and a nonnegative InterfaceSlot!");
          if (conductors.insert(spatial_edge.conductor).second)
          {
            MFEM_VERIFY(spatial_edge.conductor == next_conductor++,
                        "SpatialEdgeCluster conductor labels must be canonical and "
                        "contiguous in order of first occurrence!");
          }
          model.spatial_edges.push_back(spatial_edge);
        }
      }
    }
    const bool corner = model.topology == LibraryTopology::CONVEX_CORNER ||
                        model.topology == LibraryTopology::CONCAVE_CORNER;
    const bool endpoint = model.topology == LibraryTopology::ENDPOINT;
    const bool junction = model.topology == LibraryTopology::JUNCTION;
    const bool spatial_vertex = corner || endpoint || junction;
    const bool spatial_response = spatial_vertex || spatial_cluster;
    if (auto boundary = entry.find("PlanViewBoundary"); boundary != entry.end())
    {
      MFEM_VERIFY((spatial_cluster || endpoint || junction) && boundary->is_array() &&
                      !boundary->empty(),
                  "PlanViewBoundary requires a SpatialEdgeCluster, Endpoint, or Junction "
                  "model and must be a nonempty array!");
      for (const auto &component : *boundary)
      {
        MFEM_VERIFY(component.is_object() && component.contains("Conductor") &&
                        component["Conductor"].is_number_integer() &&
                        component["Conductor"].get<int>() > 0 &&
                        component.contains("Segments") &&
                        component["Segments"].is_array() && !component["Segments"].empty(),
                    "Invalid PlanViewBoundary component!");
        auto ValidSegment = [](const auto &segment)
        {
          return segment.is_array() && segment.size() == 2 &&
                 std::all_of(segment.begin(), segment.end(),
                             [](const auto &point)
                             {
                               return point.is_array() && point.size() == 3 &&
                                      std::all_of(
                                          point.begin(), point.end(),
                                          [](const auto &coordinate)
                                          { return coordinate.is_number_integer(); });
                             });
        };
        MFEM_VERIFY(std::all_of(component["Segments"].begin(), component["Segments"].end(),
                                ValidSegment),
                    "Invalid PlanViewBoundary segment!");
        if (auto continuation = component.find("ContinuationSegments");
            continuation != component.end())
        {
          MFEM_VERIFY(
              continuation->is_array() &&
                  std::all_of(continuation->begin(), continuation->end(), ValidSegment),
              "Invalid PlanViewBoundary continuation segment!");
          const std::set<nlohmann::json> segments(component["Segments"].begin(),
                                                  component["Segments"].end());
          MFEM_VERIFY(std::all_of(continuation->begin(), continuation->end(),
                                  [&](const auto &segment)
                                  { return segments.find(segment) != segments.end(); }),
                      "Every PlanViewBoundary continuation segment must also appear in "
                      "Segments!");
        }
      }
      model.plan_view_boundary = boundary->dump();
      if (auto regularization = entry.find("MaskRegularization");
          regularization != entry.end())
      {
        MFEM_VERIFY(
            regularization->is_object() && regularization->value("Version", 0) == 1 &&
                regularization->value("PhysicalBoundary", std::string{}) ==
                    "TaperAndRound" &&
                regularization->value("ContinuationBoundary", std::string{}) == "Vertical",
            "MaskRegularization must select the supported version-1 tapered physical "
            "boundary and vertical continuation policy!");
        model.mask_regularization = regularization->dump();
      }
    }
    MFEM_VERIFY(!entry.contains("MaskRegularization") || model.plan_view_boundary,
                "MaskRegularization requires PlanViewBoundary!");
    if (spatial_cluster)
    {
      MFEM_VERIFY(!entry.contains("BoundaryCondition"),
                  "SpatialEdgeCluster boundary conditions must be specified separately "
                  "for every Edges entry!");
    }
    else
    {
      model.boundary_condition = ParseLibraryBoundaryCondition(
          entry.value("BoundaryCondition", nlohmann::json("PEC")), units,
          nondimensionalize);
    }
    if (auto qualification = entry.find("BoundaryLawQualification");
        qualification != entry.end())
    {
      MFEM_VERIFY(
          qualification->is_object() && qualification->value("Version", 0) == 1 &&
              qualification->contains("Status") &&
              qualification->at("Status").is_string() &&
              (qualification->at("Status") == "Qualified" ||
               qualification->at("Status") == "Unqualified"),
          "BoundaryLawQualification must be a version-1 object with Status equal to "
          "Qualified or Unqualified!");
      model.boundary_law_physics_qualified = qualification->at("Status") == "Qualified";
    }
    if (model.topology == LibraryTopology::ISOLATED_EDGE || spatial_response ||
        parallel_cluster)
    {
      MFEM_VERIFY(model.separation == 0.0,
                  "An isolated-edge, parallel-cluster, or spatial response model cannot "
                  "specify a separation!");
    }
    else
    {
      MFEM_VERIFY(std::isfinite(model.separation) && model.separation > 0.0 &&
                      std::isfinite(model.separation_tolerance) &&
                      model.separation_tolerance >= 0.0,
                  "Paired-edge response models require a positive separation and a "
                  "nonnegative separation tolerance!");
    }
    if (corner)
    {
      MFEM_VERIFY(std::isfinite(model.angle) && model.angle > 0.0 &&
                      model.angle < std::acos(-1.0) &&
                      std::isfinite(model.angle_tolerance) && model.angle_tolerance >= 0.0,
                  "Corner response models require Angle strictly between zero and 180 "
                  "degrees and a nonnegative AngleTolerance!");
      MFEM_VERIFY(std::isfinite(model.corner_radius) && model.corner_radius >= 0.0 &&
                      model.corner_radius < library.matching_radius &&
                      std::isfinite(model.corner_radius_tolerance) &&
                      model.corner_radius_tolerance >= 0.0,
                  "Corner response models require CornerRadius in [0, MatchingRadius) "
                  "and a nonnegative CornerRadiusTolerance!");
      MFEM_VERIFY(model.arm_angles.empty(),
                  "Corner response models cannot specify ArmAngles!");
    }
    else if (junction)
    {
      const double full_angle = 2.0 * std::acos(-1.0);
      MFEM_VERIFY(
          model.angle == 0.0 && model.corner_radius == 0.0 &&
              model.arm_angles.size() >= 3 && std::isfinite(model.arm_angle_tolerance) &&
              model.arm_angle_tolerance >= 0.0 &&
              std::abs(model.arm_angles.front()) <= 1.0e-12 &&
              std::all_of(
                  model.arm_angles.begin(), model.arm_angles.end(),
                  [full_angle](double angle)
                  { return std::isfinite(angle) && angle >= 0.0 && angle < full_angle; }) &&
              std::adjacent_find(model.arm_angles.begin(), model.arm_angles.end(),
                                 [](double first, double second)
                                 { return first >= second; }) == model.arm_angles.end(),
          "Junction response models require strictly increasing ArmAngles in [0, 360) "
          "degrees, beginning with zero, and a nonnegative ArmAngleTolerance!");
    }
    else if (endpoint)
    {
      MFEM_VERIFY(model.angle == 0.0 && model.corner_radius == 0.0 &&
                      model.arm_angles.empty(),
                  "Endpoint response models cannot specify Angle, CornerRadius, or "
                  "ArmAngles!");
    }
    else
    {
      MFEM_VERIFY(model.angle == 0.0 && model.corner_radius == 0.0 &&
                      model.arm_angles.empty(),
                  "Straight-edge response models cannot specify Angle, CornerRadius, or "
                  "ArmAngles!");
    }
    MFEM_VERIFY(parallel_cluster == !model.cluster_edges.empty(),
                "ParallelEdgeCluster response models require Edges, and other "
                "topologies cannot specify them!");
    MFEM_VERIFY(spatial_cluster == !model.spatial_edges.empty(),
                "SpatialEdgeCluster response models require spatial Edges, and other "
                "topologies cannot specify them!");
    MFEM_VERIFY(spatial_cluster || endpoint || junction || !model.plan_view_boundary,
                "PlanViewBoundary is supported only by SpatialEdgeCluster, Endpoint, or "
                "Junction models!");
    if (spatial_response)
    {
      model.response.spatial_basis = true;
      model.response.contour_groups = entry.value("ContourGroups", std::vector<int>{});
      MFEM_VERIFY(std::all_of(model.response.contour_groups.begin(),
                              model.response.contour_groups.end(),
                              [](int size) { return size >= 3; }),
                  "Every spatial response ContourGroups entry must contain at "
                  "least three knots!");
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
    if (auto references = entry.find("ConductorReferences"); references != entry.end())
    {
      MFEM_VERIFY(version >= 2,
                  "ConductorReferences requires a version-2 fabrication-process "
                  "response library!");
      MFEM_VERIFY(!entry.contains("Reference") && references->is_array() &&
                      references->size() >= 2 &&
                      (model.topology == LibraryTopology::DIFFERENT_CONDUCTOR_GAP ||
                       parallel_cluster || spatial_cluster),
                  "ConductorReferences must contain at least two points and is supported "
                  "only by a DifferentConductorGap, ParallelEdgeCluster, or "
                  "SpatialEdgeCluster model without Reference!");
      model.conductor_references = references->get<std::vector<std::array<double, 3>>>();
      for (auto &reference : model.conductor_references)
      {
        for (double &value : reference)
        {
          value /= coordinate_scale;
          MFEM_VERIFY(std::isfinite(value),
                      "Fabrication-process response-model ConductorReferences must be "
                      "finite!");
        }
      }
      for (std::size_t i = 0; i < model.conductor_references.size(); i++)
      {
        for (std::size_t j = i + 1; j < model.conductor_references.size(); j++)
        {
          double reference_separation = 0.0;
          for (int d = 0; d < 3; d++)
          {
            const double delta =
                model.conductor_references[i][d] - model.conductor_references[j][d];
            reference_separation += delta * delta;
          }
          MFEM_VERIFY(reference_separation > 0.0,
                      "Fabrication-process response-model ConductorReferences must be "
                      "distinct!");
        }
      }
      model.response.conductor_state_count =
          static_cast<int>(model.conductor_references.size()) - 1;
    }
    else
    {
      auto reference = entry.value("Reference", std::array<double, 3>{0.0, 0.0, 0.0});
      for (double &value : reference)
      {
        value /= coordinate_scale;
      }
      model.conductor_references.push_back(reference);
    }
    if (parallel_cluster || spatial_cluster)
    {
      const int conductor_count =
          parallel_cluster
              ? std::max_element(model.cluster_edges.begin(), model.cluster_edges.end(),
                                 [](const auto &first, const auto &second)
                                 { return first.conductor < second.conductor; })
                    ->conductor
              : std::max_element(model.spatial_edges.begin(), model.spatial_edges.end(),
                                 [](const auto &first, const auto &second)
                                 { return first.conductor < second.conductor; })
                    ->conductor;
      MFEM_VERIFY(static_cast<int>(model.conductor_references.size()) == conductor_count,
                  "Edge-cluster models require one conductor reference for each "
                  "canonical conductor label!");
    }
    if (auto paths = entry.find("OpenContourPaths"); paths != entry.end())
    {
      MFEM_VERIFY(model.conductor_references.size() >= 2 &&
                      (model.topology == LibraryTopology::DIFFERENT_CONDUCTOR_GAP ||
                       parallel_cluster || spatial_cluster) &&
                      model.response.contour_groups.empty() && paths->is_array() &&
                      !paths->empty(),
                  "OpenContourPaths requires a version-2 DifferentConductorGap or "
                  "edge-cluster model with ConductorReferences and no ContourGroups!");
      std::set<int> point_indices;
      for (const auto &path : *paths)
      {
        ResponseModelData::OpenContourPathData path_data;
        path_data.indices = path.at("Indices").get<std::vector<int>>();
        path_data.start_conductor = path.at("StartConductor").get<int>() - 1;
        path_data.end_conductor = path.at("EndConductor").get<int>() - 1;
        MFEM_VERIFY(!path_data.indices.empty() && path_data.start_conductor >= 0 &&
                        path_data.start_conductor <
                            static_cast<int>(model.conductor_references.size()) &&
                        path_data.end_conductor >= 0 &&
                        path_data.end_conductor <
                            static_cast<int>(model.conductor_references.size()) &&
                        path_data.start_conductor != path_data.end_conductor,
                    "Every OpenContourPaths entry must contain at least one point and "
                    "connect two distinct valid conductor references!");
        for (int &index : path_data.indices)
        {
          MFEM_VERIFY(index > 0 && point_indices.insert(index).second,
                      "OpenContourPaths BasisPoints indices must be positive and unique!");
          index--;
        }
        model.response.open_contour_paths.push_back(std::move(path_data));
      }
      std::vector<bool> connected(model.conductor_references.size(), false);
      connected.front() = true;
      bool changed = true;
      while (changed)
      {
        changed = false;
        for (const auto &path : model.response.open_contour_paths)
        {
          if (connected[path.start_conductor] != connected[path.end_conductor])
          {
            connected[path.start_conductor] = true;
            connected[path.end_conductor] = true;
            changed = true;
          }
        }
      }
      MFEM_VERIFY(
          std::all_of(connected.begin(), connected.end(), [](bool value) { return value; }),
          "OpenContourPaths must connect every conductor reference!");
    }
    if (auto indices = entry.find("ZeroTraceIndices"); indices != entry.end())
    {
      MFEM_VERIFY(model.response.open_contour_paths.empty() && indices->is_array() &&
                      !indices->empty(),
                  "ZeroTraceIndices requires closed response-model contours!");
      std::set<int> point_indices;
      model.response.zero_trace_indices = indices->get<std::vector<int>>();
      for (int &index : model.response.zero_trace_indices)
      {
        MFEM_VERIFY(index > 0 && point_indices.insert(index).second,
                    "ZeroTraceIndices BasisPoints indices must be positive and unique!");
        index--;
      }
      std::sort(model.response.zero_trace_indices.begin(),
                model.response.zero_trace_indices.end());
    }
    MFEM_VERIFY(model.response.zero_trace_indices.empty() ||
                    (spatial_vertex &&
                     model.boundary_condition.type == MetalBoundaryConditionType::PEC) ||
                    (spatial_cluster &&
                     std::all_of(model.spatial_edges.begin(), model.spatial_edges.end(),
                                 [](const auto &edge)
                                 {
                                   return edge.boundary_condition.type ==
                                          MetalBoundaryConditionType::PEC;
                                 })),
                "Finite-impedance spatial response models cannot use ZeroTraceIndices!");

    std::set<std::pair<int, InterfaceDielectric>> interface_slots;
    if (auto interfaces = entry.find("Interfaces"); interfaces != entry.end())
    {
      for (const auto &interface : *interfaces)
      {
        InterfaceDielectric type = InterfaceDielectric::DEFAULT;
        FromString(interface.at("Type").get<std::string>(), type);
        const int slot = interface.value("Slot", 0);
        const int coupon = interface.at("Coupon");
        MFEM_VERIFY(type != InterfaceDielectric::DEFAULT && slot >= 0 && coupon > 0 &&
                        interface_slots.emplace(slot, type).second,
                    "Fabrication-process response-model interface mappings must have "
                    "unique Slot and Type pairs, nonnegative slots, and positive coupon "
                    "indices!");
        MFEM_VERIFY(spatial_cluster || slot == 0,
                    "Nonzero interface Slots are supported only by SpatialEdgeCluster "
                    "response models!");
        if (spatial_cluster)
        {
          MFEM_VERIFY(std::any_of(model.spatial_edges.begin(), model.spatial_edges.end(),
                                  [slot](const auto &edge)
                                  { return edge.interface_slot == slot; }),
                      "A SpatialEdgeCluster interface mapping refers to an unused "
                      "InterfaceSlot!");
        }
        model.interfaces.push_back({slot, type, coupon});
        mapped_interface_types.insert(type);
      }
    }
    MFEM_VERIFY(model.response.fabricated_surface_matrix.empty() ||
                    !model.interfaces.empty(),
                "Surface response matrices require interface mappings in the "
                "fabrication-process response model!");
    library.models.push_back(std::move(model));
  }
  if (version >= 3)
  {
    for (const auto type : mapped_interface_types)
    {
      MFEM_VERIFY(
          library.interface_layers.find(type) != library.interface_layers.end(),
          "Version-3 fabrication-process response library \""
              << library.name << "\" has a " << ToString(type)
              << " surface response but no matching Fabrication.InterfaceLayers entry!");
    }
  }
  if (auto spans = data.find("CornerRadiusInterpolation"); spans != data.end())
  {
    MFEM_VERIFY(version >= 3 && spans->is_array(),
                "CornerRadiusInterpolation requires a version-3 process library and "
                "must be an array!");
    std::map<std::string, std::size_t> model_indices;
    for (std::size_t i = 0; i < library.models.size(); i++)
    {
      model_indices.emplace(library.models[i].name, i);
    }
    for (const auto &span : *spans)
    {
      MFEM_VERIFY(span.is_object() && span.contains("LowerModel") &&
                      span.at("LowerModel").is_string() && span.contains("UpperModel") &&
                      span.at("UpperModel").is_string() && span.contains("Qualification") &&
                      span.at("Qualification").is_object(),
                  "Every CornerRadiusInterpolation entry requires string LowerModel and "
                  "UpperModel fields and a Qualification object!");
      const auto lower_name = span.at("LowerModel").get<std::string>();
      const auto upper_name = span.at("UpperModel").get<std::string>();
      const auto lower_index = model_indices.find(lower_name);
      const auto upper_index = model_indices.find(upper_name);
      MFEM_VERIFY(lower_index != model_indices.end() &&
                      upper_index != model_indices.end() &&
                      lower_index->second != upper_index->second,
                  "CornerRadiusInterpolation model names must refer to two distinct "
                  "models in the same process library!");
      const auto &lower = library.models[lower_index->second];
      const auto &upper = library.models[upper_index->second];
      const auto &qualification = span.at("Qualification");
      MFEM_VERIFY(qualification.value("Method", std::string{}) == "HeldOutCoupon" &&
                      qualification.value("Passed", false) &&
                      qualification.contains("HeldoutRadius") &&
                      qualification.at("HeldoutRadius").is_number(),
                  "CornerRadiusInterpolation Qualification requires Method=HeldOutCoupon, "
                  "Passed=true, and a numeric HeldoutRadius!");
      const double heldout_radius =
          qualification.at("HeldoutRadius").get<double>() / coordinate_scale;
      const bool corner = lower.topology == LibraryTopology::CONVEX_CORNER ||
                          lower.topology == LibraryTopology::CONCAVE_CORNER;
      MFEM_VERIFY(
          corner && upper.topology == lower.topology && lower.corner_radius > 0.0 &&
              lower.corner_radius < heldout_radius &&
              heldout_radius < upper.corner_radius &&
              std::abs(lower.angle - upper.angle) <=
                  1.0e-12 * std::max({lower.angle, upper.angle, 1.0}) &&
              CompatibleBoundaryLaw(lower.boundary_condition, upper.boundary_condition) &&
              CompatibleBoundaryLaw(upper.boundary_condition, lower.boundary_condition),
          "CornerRadiusInterpolation requires compatible same-angle corner models "
          "with positive radii bracketing the held-out radius!");
      MFEM_VERIFY(library.corner_radius_interpolation
                      .emplace(lower_index->second, upper_index->second)
                      .second,
                  "CornerRadiusInterpolation contains a duplicate model span!");
    }
  }
  return library;
}

void ValidateLibraryInterfaceLayers(
    const ProcessLibrary &library,
    const std::map<int, config::InterfaceDielectricData> &dielectrics,
    const std::set<int> &target_interfaces, double coordinate_scale)
{
  if (library.interface_layers.empty())
  {
    Mpi::Warning("Fabrication-process response library \"{}\" is version {} and has no "
                 "InterfaceLayers metadata; target dielectric thicknesses and "
                 "permittivities cannot be verified!\n",
                 library.name, library.version);
    return;
  }

  constexpr double relative_tolerance = 1.0e-10;
  const auto compatible = [](double actual, double expected)
  {
    return std::abs(actual - expected) <=
           relative_tolerance * std::max(std::abs(actual), std::abs(expected));
  };
  for (const int index : target_interfaces)
  {
    const auto dielectric = dielectrics.find(index);
    MFEM_VERIFY(dielectric != dielectrics.end(),
                "Response-correction target interface "
                    << index << " is not configured for dielectric postprocessing!");
    const auto layer = library.interface_layers.find(dielectric->second.type);
    MFEM_VERIFY(layer != library.interface_layers.end(),
                "Fabrication-process response library \""
                    << library.name << "\" has no InterfaceLayers metadata for target "
                    << "interface " << index << " (" << ToString(dielectric->second.type)
                    << ")!");
    MFEM_VERIFY(compatible(dielectric->second.t, layer->second.thickness),
                "Target dielectric interface "
                    << index << " (" << ToString(dielectric->second.type) << ") thickness "
                    << dielectric->second.t * coordinate_scale
                    << " does not match fabrication-process response library \""
                    << library.name << "\" thickness "
                    << layer->second.thickness * coordinate_scale
                    << " in mesh coordinate units!");
    MFEM_VERIFY(compatible(dielectric->second.epsilon_r, layer->second.permittivity),
                "Target dielectric interface "
                    << index << " (" << ToString(dielectric->second.type)
                    << ") permittivity " << dielectric->second.epsilon_r
                    << " does not match fabrication-process response library \""
                    << library.name << "\" permittivity " << layer->second.permittivity
                    << "!");
  }
}

void MapLibraryInterfaces(
    const LibraryModel &source,
    const std::map<int, std::map<InterfaceDielectric, int>> &targets_by_slot,
    ResponseModelData &model)
{
  for (const auto &[slot, targets] : targets_by_slot)
  {
    for (const auto &[type, target] : targets)
    {
      const int interface_slot = slot;
      const InterfaceDielectric interface_type = type;
      auto interface = std::find_if(
          source.interfaces.begin(), source.interfaces.end(),
          [interface_slot, interface_type](const auto &entry)
          { return entry.slot == interface_slot && entry.type == interface_type; });
      MFEM_VERIFY(interface != source.interfaces.end() ||
                      source.response.fabricated_surface_matrix.empty(),
                  "Fabrication-process response model \""
                      << source.name << "\" has no slot " << slot << " " << ToString(type)
                      << " interface mapping!");
      if (interface != source.interfaces.end())
      {
        model.interfaces.push_back({target, interface->coupon});
      }
    }
  }
}

bool MatchLibraryInterfaces(
    const LibraryModel &model,
    const std::map<int, std::map<InterfaceDielectric, int>> &targets_by_slot)
{
  if (model.response.fabricated_surface_matrix.empty())
  {
    return true;
  }

  std::map<int, std::set<InterfaceDielectric>> model_types_by_slot;
  for (const auto &interface : model.interfaces)
  {
    model_types_by_slot[interface.slot].insert(interface.type);
  }
  if (model_types_by_slot.size() != targets_by_slot.size())
  {
    return false;
  }
  for (const auto &[slot, targets] : targets_by_slot)
  {
    const auto model_types = model_types_by_slot.find(slot);
    if (model_types == model_types_by_slot.end())
    {
      return false;
    }
    for (const auto &[type, target] : targets)
    {
      (void)target;
      if (model_types->second.find(type) == model_types->second.end())
      {
        return false;
      }
    }
  }
  return true;
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

Point3D Cross(const Point3D &a, const Point3D &b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
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

struct PlanViewFacet
{
  int conductor = 0;
  std::vector<Point3D> points;
};

struct PlanViewGeometry
{
  std::vector<PlanViewFacet> facets;
  std::array<double, 2> lower{};
  std::array<double, 2> upper{};
  int process_axis = 1;
};

using QuantizedPoint3D = std::array<long long int, 3>;
using QuantizedSegment3D = std::pair<QuantizedPoint3D, QuantizedPoint3D>;

bool IntegerCrossIsZero(const QuantizedPoint3D &a, const QuantizedPoint3D &b)
{
  using WideInteger = __int128;
  return static_cast<WideInteger>(a[1]) * b[2] - static_cast<WideInteger>(a[2]) * b[1] ==
             0 &&
         static_cast<WideInteger>(a[2]) * b[0] - static_cast<WideInteger>(a[0]) * b[2] ==
             0 &&
         static_cast<WideInteger>(a[0]) * b[1] - static_cast<WideInteger>(a[1]) * b[0] == 0;
}

std::string CanonicalPlanViewBoundary(
    const std::vector<PlanViewFacet> &facets, double matching_radius, int process_axis = 1,
    const std::optional<std::pair<std::array<double, 2>, std::array<double, 2>>>
        &clip_bounds = std::nullopt)
{
  MFEM_VERIFY(matching_radius > 0.0,
              "Plan-view canonicalization requires a positive matching radius!");
  MFEM_VERIFY(process_axis >= 0 && process_axis < 3,
              "Plan-view canonicalization requires a valid process axis!");
  const std::array<int, 2> plan_axes = process_axis == 0   ? std::array<int, 2>{1, 2}
                                       : process_axis == 1 ? std::array<int, 2>{0, 2}
                                                           : std::array<int, 2>{0, 1};
  const double tolerance = 1.0e-9 * matching_radius;
  auto Quantize = [&](const Point3D &point)
  {
    QuantizedPoint3D key;
    for (int d = 0; d < 3; d++)
    {
      key[d] = std::llround(point[d] / tolerance);
    }
    return key;
  };
  auto SubtractInteger = [](const QuantizedPoint3D &a, const QuantizedPoint3D &b)
  { return QuantizedPoint3D{a[0] - b[0], a[1] - b[1], a[2] - b[2]}; };
  auto OrderedSegment = [](QuantizedPoint3D first, QuantizedPoint3D second)
  {
    if (second < first)
    {
      std::swap(first, second);
    }
    return QuantizedSegment3D{first, second};
  };

  using GroupKey = std::pair<int, long long int>;
  std::map<GroupKey, std::vector<std::vector<QuantizedPoint3D>>> polygons_by_group;
  for (const auto &facet : facets)
  {
    MFEM_VERIFY(facet.conductor > 0 && facet.points.size() >= 3,
                "Invalid plan-view facet!");
    std::vector<QuantizedPoint3D> ring;
    ring.reserve(facet.points.size());
    for (const auto &point : facet.points)
    {
      const auto key = Quantize(point);
      if (ring.empty() || ring.back() != key)
      {
        ring.push_back(key);
      }
    }
    if (ring.size() > 1 && ring.front() == ring.back())
    {
      ring.pop_back();
    }
    if (ring.size() < 3)
    {
      continue;
    }
    const long long int plane = ring.front()[process_axis];
    MFEM_VERIFY(std::all_of(ring.begin(), ring.end(), [=](const auto &point)
                            { return point[process_axis] == plane; }),
                "Plan-view facet is not on one process plane!");
    polygons_by_group[{facet.conductor, plane}].push_back(std::move(ring));
  }

  nlohmann::json result = nlohmann::json::array();
  for (auto &[group, polygons] : polygons_by_group)
  {
    std::set<std::vector<QuantizedPoint3D>> seen;
    std::vector<std::vector<QuantizedPoint3D>> unique_polygons;
    for (auto &polygon : polygons)
    {
      std::vector<QuantizedPoint3D> canonical;
      for (std::size_t start = 0; start < polygon.size(); start++)
      {
        for (const bool reverse : {false, true})
        {
          std::vector<QuantizedPoint3D> candidate;
          candidate.reserve(polygon.size());
          for (std::size_t step = 0; step < polygon.size(); step++)
          {
            const std::size_t index = reverse
                                          ? (start + polygon.size() - step) % polygon.size()
                                          : (start + step) % polygon.size();
            candidate.push_back(polygon[index]);
          }
          if (canonical.empty() || candidate < canonical)
          {
            canonical = std::move(candidate);
          }
        }
      }
      if (seen.insert(canonical).second)
      {
        unique_polygons.push_back(std::move(polygon));
      }
    }

    std::vector<QuantizedPoint3D> vertices;
    for (const auto &polygon : unique_polygons)
    {
      vertices.insert(vertices.end(), polygon.begin(), polygon.end());
    }
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());

    std::map<QuantizedSegment3D, int> counts;
    for (const auto &polygon : unique_polygons)
    {
      for (std::size_t i = 0; i < polygon.size(); i++)
      {
        const auto &begin = polygon[i];
        const auto &end = polygon[(i + 1) % polygon.size()];
        const auto direction = SubtractInteger(end, begin);
        if (begin == end)
        {
          continue;
        }
        std::vector<QuantizedPoint3D> split{begin, end};
        for (const auto &point : vertices)
        {
          if (point == begin || point == end)
          {
            continue;
          }
          const auto offset = SubtractInteger(point, begin);
          if (!IntegerCrossIsZero(direction, offset))
          {
            continue;
          }
          const long double coordinate =
              static_cast<long double>(offset[0]) * direction[0] +
              static_cast<long double>(offset[1]) * direction[1] +
              static_cast<long double>(offset[2]) * direction[2];
          const long double length_squared =
              static_cast<long double>(direction[0]) * direction[0] +
              static_cast<long double>(direction[1]) * direction[1] +
              static_cast<long double>(direction[2]) * direction[2];
          if (coordinate > 0.0L && coordinate < length_squared)
          {
            split.push_back(point);
          }
        }
        std::sort(split.begin(), split.end(),
                  [&](const auto &first, const auto &second)
                  {
                    const auto a = SubtractInteger(first, begin);
                    const auto b = SubtractInteger(second, begin);
                    const long double a_coordinate =
                        static_cast<long double>(a[0]) * direction[0] +
                        static_cast<long double>(a[1]) * direction[1] +
                        static_cast<long double>(a[2]) * direction[2];
                    const long double b_coordinate =
                        static_cast<long double>(b[0]) * direction[0] +
                        static_cast<long double>(b[1]) * direction[1] +
                        static_cast<long double>(b[2]) * direction[2];
                    return a_coordinate < b_coordinate;
                  });
        split.erase(std::unique(split.begin(), split.end()), split.end());
        for (std::size_t j = 1; j < split.size(); j++)
        {
          counts[OrderedSegment(split[j - 1], split[j])]++;
        }
      }
    }

    MFEM_VERIFY(std::none_of(counts.begin(), counts.end(),
                             [](const auto &entry) { return entry.second > 2; }),
                "Plan-view facets form a nonmanifold surface!");
    std::set<QuantizedSegment3D> boundary;
    for (const auto &[segment, count] : counts)
    {
      if (count % 2 == 1)
      {
        boundary.insert(segment);
      }
    }
    std::map<QuantizedPoint3D, int> degree;
    for (const auto &[first, second] : boundary)
    {
      degree[first]++;
      degree[second]++;
    }
    MFEM_VERIFY(std::none_of(degree.begin(), degree.end(),
                             [](const auto &entry) { return entry.second % 2 == 1; }),
                "Plan-view facet union has an open boundary!");

    while (true)
    {
      std::map<QuantizedPoint3D, std::set<QuantizedPoint3D>> adjacency;
      for (const auto &[first, second] : boundary)
      {
        adjacency[first].insert(second);
        adjacency[second].insert(first);
      }
      std::optional<std::tuple<QuantizedPoint3D, QuantizedPoint3D, QuantizedPoint3D>> merge;
      for (const auto &[vertex, neighbors] : adjacency)
      {
        if (neighbors.size() != 2)
        {
          continue;
        }
        auto neighbor = neighbors.begin();
        const auto first = *neighbor++;
        const auto second = *neighbor;
        if (IntegerCrossIsZero(SubtractInteger(first, vertex),
                               SubtractInteger(second, vertex)))
        {
          merge = std::make_tuple(vertex, first, second);
          break;
        }
      }
      if (!merge)
      {
        break;
      }
      const auto &[vertex, first, second] = *merge;
      boundary.erase(OrderedSegment(first, vertex));
      boundary.erase(OrderedSegment(vertex, second));
      boundary.insert(OrderedSegment(first, second));
    }
    MFEM_VERIFY(!boundary.empty(), "Plan-view facets have no union boundary!");

    nlohmann::json segments = nlohmann::json::array();
    nlohmann::json continuation_segments = nlohmann::json::array();
    std::array<QuantizedPoint3D, 2> quantized_bounds{};
    if (clip_bounds)
    {
      for (int side = 0; side < 2; side++)
      {
        Point3D point{};
        point[plan_axes[0]] = side == 0 ? clip_bounds->first[0] : clip_bounds->second[0];
        point[plan_axes[1]] = side == 0 ? clip_bounds->first[1] : clip_bounds->second[1];
        quantized_bounds[side] = Quantize(point);
      }
    }
    for (const auto &[first, second] : boundary)
    {
      segments.push_back({first, second});
      if (clip_bounds && ((first[plan_axes[0]] == second[plan_axes[0]] &&
                           (first[plan_axes[0]] == quantized_bounds[0][plan_axes[0]] ||
                            first[plan_axes[0]] == quantized_bounds[1][plan_axes[0]])) ||
                          (first[plan_axes[1]] == second[plan_axes[1]] &&
                           (first[plan_axes[1]] == quantized_bounds[0][plan_axes[1]] ||
                            first[plan_axes[1]] == quantized_bounds[1][plan_axes[1]]))))
      {
        continuation_segments.push_back({first, second});
      }
    }
    nlohmann::json component = {{"Conductor", group.first},
                                {"Segments", std::move(segments)}};
    if (clip_bounds)
    {
      component["ContinuationSegments"] = std::move(continuation_segments);
    }
    result.push_back(std::move(component));
  }
  std::sort(result.begin(), result.end(), [](const auto &first, const auto &second)
            { return first.dump() < second.dump(); });
  return result.dump();
}

bool HasClassifiedPlanViewBoundary(const std::string &boundary)
{
  const auto components = nlohmann::json::parse(boundary);
  return std::all_of(components.begin(), components.end(), [](const auto &component)
                     { return component.contains("ContinuationSegments"); });
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

std::optional<int> GetConductor(const config::BoundaryData &boundaries, int attribute,
                                bool pec_attribute_conductors = false)
{
  if (std::find(boundaries.pec.attributes.begin(), boundaries.pec.attributes.end(),
                attribute) != boundaries.pec.attributes.end())
  {
    return pec_attribute_conductors ? attribute : 0;
  }
  if (std::find(boundaries.auxpec.attributes.begin(), boundaries.auxpec.attributes.end(),
                attribute) != boundaries.auxpec.attributes.end())
  {
    return pec_attribute_conductors ? attribute : 0;
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
            potential.attributes.end() ||
        std::find(potential.terminal_attributes.begin(),
                  potential.terminal_attributes.end(),
                  attribute) != potential.terminal_attributes.end())
    {
      return index;
    }
  }
  if (pec_attribute_conductors)
  {
    const auto HasAttribute = [attribute](const auto &data)
    {
      return std::find(data.attributes.begin(), data.attributes.end(), attribute) !=
             data.attributes.end();
    };
    if (std::any_of(boundaries.conductivity.begin(), boundaries.conductivity.end(),
                    HasAttribute) ||
        std::any_of(boundaries.impedance.begin(), boundaries.impedance.end(),
                    HasAttribute) ||
        std::any_of(boundaries.rational_impedance.begin(),
                    boundaries.rational_impedance.end(), HasAttribute))
    {
      return attribute;
    }
  }
  return std::nullopt;
}

MetalBoundaryLaw GetBoundaryConditionLaw(const config::BoundaryData &boundaries,
                                         const MetalBoundaryCondition &condition)
{
  MetalBoundaryLaw law;
  law.type = condition.type;
  switch (condition.type)
  {
    case MetalBoundaryConditionType::PEC:
      break;
    case MetalBoundaryConditionType::CONDUCTIVITY:
      {
        MFEM_VERIFY(condition.index >= 0 &&
                        condition.index < static_cast<int>(boundaries.conductivity.size()),
                    "Invalid conductivity boundary-law index!");
        const auto &data = boundaries.conductivity[condition.index];
        law.parameters = {data.sigma, data.mu_r, data.external ? 2.0 * data.h : data.h};
        break;
      }
    case MetalBoundaryConditionType::IMPEDANCE:
      {
        MFEM_VERIFY(condition.index >= 0 &&
                        condition.index < static_cast<int>(boundaries.impedance.size()),
                    "Invalid impedance boundary-law index!");
        const auto &data = boundaries.impedance[condition.index];
        law.parameters = {data.Rs, data.Ls, data.Cs};
        break;
      }
    case MetalBoundaryConditionType::RATIONAL_IMPEDANCE:
      {
        MFEM_VERIFY(condition.index >= 0 &&
                        condition.index <
                            static_cast<int>(boundaries.rational_impedance.size()),
                    "Invalid rational-impedance boundary-law index!");
        const auto &data = boundaries.rational_impedance[condition.index];
        law.numerator = data.num;
        law.denominator = data.den;
        NormalizeRationalLaw(law);
        break;
      }
  }
  return law;
}

std::optional<MetalBoundaryLaw>
GetBoundaryConditionLaw(const config::BoundaryData &boundaries, int attribute)
{
  if (std::find(boundaries.pec.attributes.begin(), boundaries.pec.attributes.end(),
                attribute) != boundaries.pec.attributes.end() ||
      std::find(boundaries.auxpec.attributes.begin(), boundaries.auxpec.attributes.end(),
                attribute) != boundaries.auxpec.attributes.end())
  {
    return MetalBoundaryLaw{};
  }
  for (const auto &[index, terminal] : boundaries.terminal)
  {
    (void)index;
    if (std::find(terminal.attributes.begin(), terminal.attributes.end(), attribute) !=
        terminal.attributes.end())
    {
      return MetalBoundaryLaw{};
    }
  }
  for (const auto &[index, potential] : boundaries.prescribed_potential)
  {
    (void)index;
    if (std::find(potential.attributes.begin(), potential.attributes.end(), attribute) !=
            potential.attributes.end() ||
        std::find(potential.terminal_attributes.begin(),
                  potential.terminal_attributes.end(),
                  attribute) != potential.terminal_attributes.end())
    {
      return MetalBoundaryLaw{};
    }
  }
  const auto HasAttribute = [attribute](const auto &data)
  {
    return std::find(data.attributes.begin(), data.attributes.end(), attribute) !=
           data.attributes.end();
  };
  for (std::size_t i = 0; i < boundaries.conductivity.size(); i++)
  {
    if (HasAttribute(boundaries.conductivity[i]))
    {
      return GetBoundaryConditionLaw(
          boundaries, {MetalBoundaryConditionType::CONDUCTIVITY, static_cast<int>(i)});
    }
  }
  for (std::size_t i = 0; i < boundaries.impedance.size(); i++)
  {
    if (HasAttribute(boundaries.impedance[i]))
    {
      return GetBoundaryConditionLaw(
          boundaries, {MetalBoundaryConditionType::IMPEDANCE, static_cast<int>(i)});
    }
  }
  for (std::size_t i = 0; i < boundaries.rational_impedance.size(); i++)
  {
    if (HasAttribute(boundaries.rational_impedance[i]))
    {
      return GetBoundaryConditionLaw(
          boundaries,
          {MetalBoundaryConditionType::RATIONAL_IMPEDANCE, static_cast<int>(i)});
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
                                         const EdgeGroup2D &group,
                                         bool pec_attribute_conductors)
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
    metal_attributes.insert(potential.terminal_attributes.begin(),
                            potential.terminal_attributes.end());
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

  std::vector<int> metal_segment_conductors(segments.size(),
                                            std::numeric_limits<int>::max());
  if (pec_attribute_conductors)
  {
    int next_conductor = std::numeric_limits<int>::min();
    for (std::size_t seed = 0; seed < segments.size(); seed++)
    {
      if (metal_attributes.find(segments[seed].attribute) == metal_attributes.end() ||
          metal_segment_conductors[seed] != std::numeric_limits<int>::max())
      {
        continue;
      }
      std::vector<std::size_t> queue = {seed};
      metal_segment_conductors[seed] = next_conductor++;
      for (std::size_t cursor = 0; cursor < queue.size(); cursor++)
      {
        const auto &current = segments[queue[cursor]];
        for (std::size_t neighbor = 0; neighbor < segments.size(); neighbor++)
        {
          if (metal_attributes.find(segments[neighbor].attribute) ==
                  metal_attributes.end() ||
              metal_segment_conductors[neighbor] != std::numeric_limits<int>::max())
          {
            continue;
          }
          const auto &candidate = segments[neighbor];
          const bool connected = Distance(current.p0, candidate.p0) <= tolerance ||
                                 Distance(current.p0, candidate.p1) <= tolerance ||
                                 Distance(current.p1, candidate.p0) <= tolerance ||
                                 Distance(current.p1, candidate.p1) <= tolerance;
          if (connected)
          {
            metal_segment_conductors[neighbor] = metal_segment_conductors[seed];
            queue.push_back(neighbor);
          }
        }
      }
    }
  }

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
    std::optional<MetalBoundaryLaw> boundary_condition;
    for (std::size_t segment_index = 0; segment_index < segments.size(); segment_index++)
    {
      const auto &segment = segments[segment_index];
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
      std::optional<int> conductor;
      if (pec_attribute_conductors &&
          metal_attributes.find(segment.attribute) != metal_attributes.end())
      {
        MFEM_ASSERT(metal_segment_conductors[segment_index] !=
                        std::numeric_limits<int>::max(),
                    "Missing connected metal component!");
        conductor = metal_segment_conductors[segment_index];
      }
      else
      {
        conductor = GetConductor(boundaries, segment.attribute, pec_attribute_conductors);
      }
      if (conductor)
      {
        conductors.insert(*conductor);
        const auto segment_condition =
            GetBoundaryConditionLaw(boundaries, segment.attribute);
        MFEM_VERIFY(segment_condition,
                    "Unable to determine a two-dimensional metal boundary condition!");
        MFEM_VERIFY(
            !boundary_condition || SameBoundaryLaw(*boundary_condition, *segment_condition),
            "A two-dimensional metal edge cannot mix distinct metal boundary conditions!");
        boundary_condition = *segment_condition;
      }
    }
    MFEM_VERIFY(Norm(inward) > 0.0,
                "Unable to infer the in-plane direction of an automatically detected "
                "two-dimensional metal edge!");
    MFEM_VERIFY(conductors.size() == 1,
                "Unable to assign an automatically detected two-dimensional metal edge "
                "to exactly one conductor!");

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
    site.boundary_condition = boundary_condition.value_or(MetalBoundaryLaw{});
    sites.push_back(site);
  }
  return sites;
}

std::optional<LibrarySelection> FindLibraryModel(const ProcessLibrary &library,
                                                 LibraryTopology topology,
                                                 double separation,
                                                 const MetalBoundaryLaw &boundary_condition)
{
  std::optional<std::size_t> best;
  double best_distance = mfem::infinity();
  for (std::size_t i = 0; i < library.models.size(); i++)
  {
    const auto &model = library.models[i];
    if (model.topology != topology ||
        !CompatibleBoundaryLaw(model.boundary_condition, boundary_condition))
    {
      continue;
    }
    const double error = std::abs(model.separation - separation);
    const double tolerance =
        std::max(model.separation_tolerance,
                 1.0e-10 * std::max(library.matching_radius, separation));
    const double distance = error / tolerance;
    if (error > tolerance)
    {
      continue;
    }
    const bool prefer_conductor_state =
        best && distance == best_distance &&
        model.conductor_references.size() >
            library.models[*best].conductor_references.size();
    if (distance < best_distance || prefer_conductor_state)
    {
      best = i;
      best_distance = distance;
    }
  }
  if (best)
  {
    LibrarySelection selection;
    selection.models.push_back({*best, 1.0});
    selection.conductor_references = library.models[*best].conductor_references;
    selection.normalized_distance = best_distance;
    return selection;
  }

  if (topology == LibraryTopology::ISOLATED_EDGE)
  {
    return std::nullopt;
  }

  std::optional<std::pair<std::size_t, std::size_t>> bracket;
  double best_span = mfem::infinity();
  std::size_t best_reference_count = 0;
  for (std::size_t lower = 0; lower < library.models.size(); lower++)
  {
    const auto &lower_model = library.models[lower];
    if (lower_model.topology != topology ||
        !CompatibleBoundaryLaw(lower_model.boundary_condition, boundary_condition) ||
        lower_model.separation >= separation)
    {
      continue;
    }
    for (std::size_t upper = 0; upper < library.models.size(); upper++)
    {
      const auto &upper_model = library.models[upper];
      if (upper_model.topology != topology ||
          !CompatibleBoundaryLaw(upper_model.boundary_condition, boundary_condition) ||
          upper_model.separation <= separation ||
          lower_model.conductor_references.size() !=
              upper_model.conductor_references.size())
      {
        continue;
      }
      const double span = upper_model.separation - lower_model.separation;
      if (span < best_span ||
          (span == best_span &&
           lower_model.conductor_references.size() > best_reference_count))
      {
        bracket = std::pair{lower, upper};
        best_span = span;
        best_reference_count = lower_model.conductor_references.size();
      }
    }
  }
  if (!bracket)
  {
    return std::nullopt;
  }

  const auto &lower_model = library.models[bracket->first];
  const auto &upper_model = library.models[bracket->second];
  const double upper_weight = (separation - lower_model.separation) / best_span;
  const double lower_weight = 1.0 - upper_weight;
  LibrarySelection selection;
  selection.models = {{bracket->first, lower_weight}, {bracket->second, upper_weight}};
  selection.conductor_references.resize(lower_model.conductor_references.size());
  for (std::size_t i = 0; i < selection.conductor_references.size(); i++)
  {
    for (int d = 0; d < 3; d++)
    {
      selection.conductor_references[i][d] =
          lower_weight * lower_model.conductor_references[i][d] +
          upper_weight * upper_model.conductor_references[i][d];
    }
  }
  selection.normalized_distance = best_span / library.matching_radius;
  return selection;
}

std::optional<LibrarySelection>
FindCornerLibraryModel(const ProcessLibrary &library, LibraryTopology topology,
                       double angle, double radius,
                       const MetalBoundaryLaw &boundary_condition)
{
  std::optional<std::size_t> best;
  double best_error = mfem::infinity();
  for (std::size_t i = 0; i < library.models.size(); i++)
  {
    const auto &model = library.models[i];
    if (model.topology != topology ||
        !CompatibleBoundaryLaw(model.boundary_condition, boundary_condition))
    {
      continue;
    }
    const double error = std::abs(model.angle - angle);
    const double angle_tolerance =
        std::max(model.angle_tolerance, 1.0e-10 * std::max(model.angle, angle));
    if (error > angle_tolerance)
    {
      continue;
    }
    const double angle_distance = error / angle_tolerance;
    const double radius_error = std::abs(model.corner_radius - radius);
    const double radius_tolerance = std::max(
        model.corner_radius_tolerance, 1.0e-10 * std::max(library.matching_radius, radius));
    if (radius_error <= radius_tolerance)
    {
      const double normalized_error =
          std::max(angle_distance, radius_error / radius_tolerance);
      if (normalized_error < best_error)
      {
        best = i;
        best_error = normalized_error;
      }
    }
  }
  if (best)
  {
    LibrarySelection selection;
    selection.models.push_back({*best, 1.0});
    selection.conductor_references = library.models[*best].conductor_references;
    selection.normalized_distance = best_error;
    return selection;
  }

  // A sharp-corner model is not a radius-interpolation endpoint. Its singular local
  // geometry is qualitatively different from a resolved fillet, so rounded corners
  // require two positive-radius coupons and Palace never extrapolates beyond them.
  const double positive_radius_tolerance = 1.0e-10 * library.matching_radius;
  if (radius <= positive_radius_tolerance)
  {
    return std::nullopt;
  }
  std::optional<std::pair<std::size_t, std::size_t>> bracket;
  double bracket_distance = mfem::infinity();
  for (const auto &[lower_index, upper_index] : library.corner_radius_interpolation)
  {
    const auto &lower_model = library.models[lower_index];
    const auto &upper_model = library.models[upper_index];
    if (lower_model.topology != topology ||
        !CompatibleBoundaryLaw(lower_model.boundary_condition, boundary_condition) ||
        !CompatibleBoundaryLaw(upper_model.boundary_condition, boundary_condition) ||
        !(lower_model.corner_radius < radius && radius < upper_model.corner_radius))
    {
      continue;
    }
    const double lower_angle_tolerance =
        std::max(lower_model.angle_tolerance, 1.0e-10 * std::max(lower_model.angle, angle));
    const double upper_angle_tolerance =
        std::max(upper_model.angle_tolerance, 1.0e-10 * std::max(upper_model.angle, angle));
    const double lower_angle_error = std::abs(lower_model.angle - angle);
    const double upper_angle_error = std::abs(upper_model.angle - angle);
    if (lower_angle_error > lower_angle_tolerance ||
        upper_angle_error > upper_angle_tolerance)
    {
      continue;
    }
    const double span = upper_model.corner_radius - lower_model.corner_radius;
    const double normalized_distance = std::max({lower_angle_error / lower_angle_tolerance,
                                                 upper_angle_error / upper_angle_tolerance,
                                                 span / library.matching_radius});
    if (normalized_distance < bracket_distance)
    {
      bracket = std::make_pair(lower_index, upper_index);
      bracket_distance = normalized_distance;
    }
  }
  if (!bracket)
  {
    return std::nullopt;
  }

  const auto &lower_model = library.models[bracket->first];
  const auto &upper_model = library.models[bracket->second];
  const double span = upper_model.corner_radius - lower_model.corner_radius;
  MFEM_ASSERT(span > 0.0, "Invalid corner-radius interpolation bracket!");
  const double upper_weight = (radius - lower_model.corner_radius) / span;
  const double lower_weight = 1.0 - upper_weight;
  LibrarySelection selection;
  selection.models = {{bracket->first, lower_weight}, {bracket->second, upper_weight}};
  MFEM_VERIFY(lower_model.conductor_references.size() ==
                  upper_model.conductor_references.size(),
              "Corner-radius interpolation requires compatible conductor references!");
  selection.conductor_references.resize(lower_model.conductor_references.size());
  for (std::size_t i = 0; i < selection.conductor_references.size(); i++)
  {
    for (int d = 0; d < 3; d++)
    {
      selection.conductor_references[i][d] =
          lower_weight * lower_model.conductor_references[i][d] +
          upper_weight * upper_model.conductor_references[i][d];
    }
  }
  selection.normalized_distance = bracket_distance;
  return selection;
}

std::optional<ParallelClusterSelection>
FindParallelClusterLibraryModel(const ProcessLibrary &library,
                                const std::vector<EdgeSite2D> &sites,
                                const std::vector<std::size_t> &cluster)
{
  if (cluster.size() < 3)
  {
    return std::nullopt;
  }

  Point2D process_normal{};
  for (const std::size_t index : cluster)
  {
    process_normal[0] += sites[index].axis_v[0];
    process_normal[1] += sites[index].axis_v[1];
  }
  if (Norm(process_normal) == 0.0)
  {
    return std::nullopt;
  }
  process_normal = Normalize(process_normal);
  if (std::any_of(cluster.begin(), cluster.end(), [&](std::size_t index)
                  { return Dot(process_normal, sites[index].axis_v) <= 0.95; }))
  {
    return std::nullopt;
  }
  const auto boundary_condition = sites[cluster.front()].boundary_condition;
  if (std::any_of(cluster.begin(), cluster.end(),
                  [&](std::size_t index)
                  {
                    return !SameBoundaryLaw(sites[index].boundary_condition,
                                            boundary_condition);
                  }))
  {
    return std::nullopt;
  }

  std::optional<ParallelClusterSelection> best;
  double best_distance = mfem::infinity();
  for (const double orientation : {-1.0, 1.0})
  {
    const Point2D axis_u = {orientation * process_normal[1],
                            -orientation * process_normal[0]};
    std::vector<std::size_t> ordered(cluster);
    std::sort(
        ordered.begin(), ordered.end(), [&](std::size_t first, std::size_t second)
        { return Dot(sites[first].point, axis_u) < Dot(sites[second].point, axis_u); });
    const Point2D origin = sites[ordered.front()].point;
    std::vector<double> offsets;
    std::vector<int> gap_directions;
    std::vector<int> conductors;
    std::map<int, int> conductor_labels;
    std::vector<std::size_t> reference_edges;
    offsets.reserve(ordered.size());
    gap_directions.reserve(ordered.size());
    conductors.reserve(ordered.size());
    for (const std::size_t index : ordered)
    {
      const Point2D delta = {sites[index].point[0] - origin[0],
                             sites[index].point[1] - origin[1]};
      const double offset = Dot(delta, axis_u);
      if (std::abs(Dot(delta, process_normal)) > 1.0e-8 * library.matching_radius)
      {
        offsets.clear();
        break;
      }
      const double gap_dot = Dot(sites[index].axis_u, axis_u);
      if (std::abs(gap_dot) <= 0.95)
      {
        offsets.clear();
        break;
      }
      offsets.push_back(offset);
      gap_directions.push_back(gap_dot > 0.0 ? 1 : -1);
      auto [label, inserted] = conductor_labels.emplace(
          sites[index].conductor, static_cast<int>(conductor_labels.size()) + 1);
      if (inserted)
      {
        reference_edges.push_back(index);
      }
      conductors.push_back(label->second);
    }
    if (offsets.size() != ordered.size())
    {
      continue;
    }

    for (std::size_t model_index = 0; model_index < library.models.size(); model_index++)
    {
      const auto &model = library.models[model_index];
      if (model.topology != LibraryTopology::PARALLEL_EDGE_CLUSTER ||
          !CompatibleBoundaryLaw(model.boundary_condition, boundary_condition) ||
          model.cluster_edges.size() != ordered.size() ||
          model.conductor_references.size() != reference_edges.size())
      {
        continue;
      }
      const double tolerance =
          std::max(model.cluster_offset_tolerance, 1.0e-10 * library.matching_radius);
      double normalized_distance = 0.0;
      bool compatible = true;
      for (std::size_t i = 0; i < ordered.size(); i++)
      {
        const auto &edge = model.cluster_edges[i];
        const double error = std::abs(offsets[i] - edge.offset);
        compatible = compatible && error <= tolerance &&
                     gap_directions[i] == edge.gap_direction &&
                     conductors[i] == edge.conductor;
        normalized_distance = std::max(normalized_distance, error / tolerance);
      }
      if (!compatible || normalized_distance >= best_distance)
      {
        continue;
      }
      ParallelClusterSelection selection;
      selection.response.models.push_back({model_index, 1.0});
      selection.response.conductor_references = model.conductor_references;
      selection.response.normalized_distance = normalized_distance;
      selection.ordered_edges = ordered;
      selection.reference_edges = reference_edges;
      selection.axis_u = axis_u;
      selection.axis_v = process_normal;
      best = std::move(selection);
      best_distance = normalized_distance;
    }
  }
  return best;
}

Point3D InterpolateAtLongitudinalCoordinate(const EdgeSegment3D &segment,
                                            const Point3D &tangent, double coordinate)
{
  const double orientation = Dot(segment.tangent, tangent);
  MFEM_ASSERT(std::abs(orientation) > 1.0 - 1.0e-8,
              "A parallel-edge cluster contains incompatible tangents!");
  const double distance = (coordinate - Dot(segment.p0, tangent)) / orientation;
  return Interpolate(segment, std::clamp(distance, 0.0, segment.length));
}

std::optional<ParallelClusterSelection3D>
FindParallelClusterLibraryModel(const ProcessLibrary &library,
                                const std::vector<EdgeSegment3D> &segments,
                                const std::vector<std::size_t> &cluster,
                                const Point3D &tangent, double longitudinal_coordinate)
{
  if (cluster.size() < 3)
  {
    return std::nullopt;
  }

  Point3D process_normal{};
  for (const std::size_t index : cluster)
  {
    process_normal = Add(process_normal, segments[index].axis_v);
  }
  if (Norm(process_normal) == 0.0)
  {
    return std::nullopt;
  }
  process_normal = Normalize(process_normal);
  if (std::any_of(cluster.begin(), cluster.end(), [&](std::size_t index)
                  { return Dot(process_normal, segments[index].axis_v) <= 0.95; }))
  {
    return std::nullopt;
  }
  const auto boundary_condition = segments[cluster.front()].boundary_condition;
  if (std::any_of(cluster.begin(), cluster.end(),
                  [&](std::size_t index)
                  {
                    return !SameBoundaryLaw(segments[index].boundary_condition,
                                            boundary_condition);
                  }))
  {
    return std::nullopt;
  }

  const Point3D transverse = Normalize(Cross(process_normal, tangent));
  std::optional<ParallelClusterSelection3D> best;
  double best_distance = mfem::infinity();
  for (const double orientation : {-1.0, 1.0})
  {
    const Point3D axis_u = Scale(orientation, transverse);
    std::vector<std::size_t> ordered(cluster);
    std::sort(ordered.begin(), ordered.end(),
              [&](std::size_t first, std::size_t second)
              {
                const auto first_point = InterpolateAtLongitudinalCoordinate(
                    segments[first], tangent, longitudinal_coordinate);
                const auto second_point = InterpolateAtLongitudinalCoordinate(
                    segments[second], tangent, longitudinal_coordinate);
                return Dot(first_point, axis_u) < Dot(second_point, axis_u);
              });
    const Point3D origin = InterpolateAtLongitudinalCoordinate(
        segments[ordered.front()], tangent, longitudinal_coordinate);
    std::vector<double> offsets;
    std::vector<int> gap_directions;
    std::vector<int> conductors;
    std::map<int, int> conductor_labels;
    std::vector<std::size_t> reference_edges;
    offsets.reserve(ordered.size());
    gap_directions.reserve(ordered.size());
    conductors.reserve(ordered.size());
    for (const std::size_t index : ordered)
    {
      const Point3D point = InterpolateAtLongitudinalCoordinate(segments[index], tangent,
                                                                longitudinal_coordinate);
      const Point3D delta = Subtract(point, origin);
      if (std::abs(Dot(delta, process_normal)) > 1.0e-8 * library.matching_radius)
      {
        offsets.clear();
        break;
      }
      const double gap_dot = Dot(segments[index].axis_u, axis_u);
      if (std::abs(gap_dot) <= 0.95)
      {
        offsets.clear();
        break;
      }
      offsets.push_back(Dot(delta, axis_u));
      gap_directions.push_back(gap_dot > 0.0 ? 1 : -1);
      auto [label, inserted] = conductor_labels.emplace(
          segments[index].conductor, static_cast<int>(conductor_labels.size()) + 1);
      if (inserted)
      {
        reference_edges.push_back(index);
      }
      conductors.push_back(label->second);
    }
    if (offsets.size() != ordered.size())
    {
      continue;
    }

    for (std::size_t model_index = 0; model_index < library.models.size(); model_index++)
    {
      const auto &model = library.models[model_index];
      if (model.topology != LibraryTopology::PARALLEL_EDGE_CLUSTER ||
          !CompatibleBoundaryLaw(model.boundary_condition, boundary_condition) ||
          model.cluster_edges.size() != ordered.size() ||
          model.conductor_references.size() != reference_edges.size())
      {
        continue;
      }
      const double tolerance =
          std::max(model.cluster_offset_tolerance, 1.0e-10 * library.matching_radius);
      double normalized_distance = 0.0;
      bool compatible = true;
      for (std::size_t i = 0; i < ordered.size(); i++)
      {
        const auto &edge = model.cluster_edges[i];
        const double error = std::abs(offsets[i] - edge.offset);
        compatible = compatible && error <= tolerance &&
                     gap_directions[i] == edge.gap_direction &&
                     conductors[i] == edge.conductor;
        normalized_distance = std::max(normalized_distance, error / tolerance);
      }
      if (!compatible || normalized_distance >= best_distance)
      {
        continue;
      }
      ParallelClusterSelection3D selection;
      selection.response.models.push_back({model_index, 1.0});
      selection.response.conductor_references = model.conductor_references;
      selection.response.normalized_distance = normalized_distance;
      selection.ordered_edges = ordered;
      selection.reference_edges = reference_edges;
      selection.axis_u = axis_u;
      selection.axis_v = process_normal;
      best = std::move(selection);
      best_distance = normalized_distance;
    }
  }
  return best;
}

ParallelClusterSpans3D FindParallelClusterSpans(const ProcessLibrary &library,
                                                const std::vector<EdgeSegment3D> &segments,
                                                const std::vector<EdgePair3D> &pairs)
{
  ParallelClusterSpans3D result;
  std::vector<std::vector<std::size_t>> pairs_by_edge(segments.size());
  for (std::size_t pair_index = 0; pair_index < pairs.size(); pair_index++)
  {
    pairs_by_edge[pairs[pair_index].first].push_back(pair_index);
    pairs_by_edge[pairs[pair_index].second].push_back(pair_index);
  }

  std::vector<bool> visited_edge(segments.size(), false);
  for (std::size_t seed = 0; seed < segments.size(); seed++)
  {
    if (visited_edge[seed] || pairs_by_edge[seed].empty())
    {
      continue;
    }
    std::vector<std::size_t> component_edges = {seed};
    std::set<std::size_t> component_pairs;
    visited_edge[seed] = true;
    for (std::size_t cursor = 0; cursor < component_edges.size(); cursor++)
    {
      const std::size_t edge = component_edges[cursor];
      for (const std::size_t pair_index : pairs_by_edge[edge])
      {
        component_pairs.insert(pair_index);
        const auto &pair = pairs[pair_index];
        const std::size_t neighbor = pair.first == edge ? pair.second : pair.first;
        if (!visited_edge[neighbor])
        {
          visited_edge[neighbor] = true;
          component_edges.push_back(neighbor);
        }
      }
    }
    if (component_edges.size() < 3)
    {
      continue;
    }

    Point3D tangent = segments[component_edges.front()].tangent;
    for (double value : tangent)
    {
      if (std::abs(value) <= 1.0e-12)
      {
        continue;
      }
      if (value < 0.0)
      {
        tangent = Scale(-1.0, tangent);
      }
      break;
    }
    struct PairRange
    {
      std::size_t index;
      double begin;
      double end;
    };
    std::vector<PairRange> ranges;
    std::vector<double> events;
    ranges.reserve(component_pairs.size());
    events.reserve(2 * component_pairs.size());
    for (const std::size_t pair_index : component_pairs)
    {
      const auto &pair = pairs[pair_index];
      const auto &first = segments[pair.first];
      const auto &second = segments[pair.second];
      const double first_begin = std::min(Dot(first.p0, tangent), Dot(first.p1, tangent));
      const double first_end = std::max(Dot(first.p0, tangent), Dot(first.p1, tangent));
      const double second_begin =
          std::min(Dot(second.p0, tangent), Dot(second.p1, tangent));
      const double second_end = std::max(Dot(second.p0, tangent), Dot(second.p1, tangent));
      const double begin = std::max(first_begin, second_begin);
      const double end = std::min(first_end, second_end);
      if (end <= begin)
      {
        continue;
      }
      ranges.push_back({pair_index, begin, end});
      events.push_back(begin);
      events.push_back(end);
    }
    std::sort(events.begin(), events.end());
    events.erase(
        std::unique(
            events.begin(), events.end(), [&](double first, double second)
            { return std::abs(first - second) <= 1.0e-10 * library.matching_radius; }),
        events.end());
    for (std::size_t event = 1; event < events.size(); event++)
    {
      const double begin = events[event - 1];
      const double end = events[event];
      if (end <= begin)
      {
        continue;
      }
      const double midpoint = 0.5 * (begin + end);
      std::map<std::size_t, std::vector<std::size_t>> adjacency;
      for (const auto &range : ranges)
      {
        if (midpoint <= range.begin || midpoint >= range.end)
        {
          continue;
        }
        const auto &pair = pairs[range.index];
        adjacency[pair.first].push_back(pair.second);
        adjacency[pair.second].push_back(pair.first);
      }
      std::set<std::size_t> active_visited;
      for (const auto &[active_seed, neighbors] : adjacency)
      {
        (void)neighbors;
        if (!active_visited.insert(active_seed).second)
        {
          continue;
        }
        std::vector<std::size_t> active = {active_seed};
        for (std::size_t cursor = 0; cursor < active.size(); cursor++)
        {
          for (const std::size_t neighbor : adjacency[active[cursor]])
          {
            if (active_visited.insert(neighbor).second)
            {
              active.push_back(neighbor);
            }
          }
        }
        if (active.size() < 3)
        {
          continue;
        }
        const auto selection =
            FindParallelClusterLibraryModel(library, segments, active, tangent, midpoint);
        if (!selection)
        {
          result.unmatched.push_back({std::move(active), tangent, begin, end});
          continue;
        }
        result.matched.push_back({*selection, tangent, begin, end});
      }
    }
  }
  return result;
}

struct VertexLibrarySelection
{
  LibrarySelection response;
  std::size_t first_arm = 0;
};

std::optional<VertexLibrarySelection> FindVertexLibraryModel(
    const ProcessLibrary &library, LibraryTopology topology,
    const std::vector<Point3D> &directions, const Point3D &process_normal,
    const MetalBoundaryLaw &boundary_condition,
    const std::function<bool(const LibraryModel &, std::size_t)> &matches_plan_view = {})
{
  MFEM_ASSERT(topology == LibraryTopology::ENDPOINT ||
                  topology == LibraryTopology::JUNCTION,
              "Invalid spatial-vertex response topology!");
  if (topology == LibraryTopology::ENDPOINT && directions.size() != 1)
  {
    return std::nullopt;
  }

  std::optional<VertexLibrarySelection> best;
  double best_error = mfem::infinity();
  std::vector<std::size_t> model_indices(library.models.size());
  std::iota(model_indices.begin(), model_indices.end(), 0);
  std::stable_sort(model_indices.begin(), model_indices.end(),
                   [&](std::size_t first, std::size_t second)
                   {
                     return library.models[first].plan_view_boundary.has_value() >
                            library.models[second].plan_view_boundary.has_value();
                   });
  for (const std::size_t model_index : model_indices)
  {
    const auto &model = library.models[model_index];
    if (model.topology != topology ||
        !CompatibleBoundaryLaw(model.boundary_condition, boundary_condition))
    {
      continue;
    }
    if (topology == LibraryTopology::ENDPOINT)
    {
      if (matches_plan_view && !matches_plan_view(model, 0))
      {
        continue;
      }
      VertexLibrarySelection selection;
      selection.response.models.push_back({model_index, 1.0});
      selection.response.conductor_references = model.conductor_references;
      selection.response.normalized_distance = 0.0;
      return selection;
    }
    if (model.arm_angles.size() != directions.size())
    {
      continue;
    }

    const double tolerance = std::max(model.arm_angle_tolerance, 1.0e-10 * std::acos(-1.0));
    for (std::size_t first = 0; first < directions.size(); first++)
    {
      const Point3D axis_u = directions[first];
      const Point3D axis_v = Normalize(Cross(process_normal, axis_u));
      std::vector<double> angles;
      angles.reserve(directions.size());
      for (const auto &direction : directions)
      {
        double angle = std::atan2(Dot(direction, axis_v), Dot(direction, axis_u));
        if (angle < 0.0)
        {
          angle += 2.0 * std::acos(-1.0);
        }
        angles.push_back(angle);
      }
      std::sort(angles.begin(), angles.end());
      double normalized_error = 0.0;
      for (std::size_t i = 0; i < angles.size(); i++)
      {
        normalized_error = std::max(normalized_error,
                                    std::abs(angles[i] - model.arm_angles[i]) / tolerance);
      }
      if (normalized_error > 1.0 || normalized_error >= best_error)
      {
        continue;
      }
      if (matches_plan_view && !matches_plan_view(model, first))
      {
        continue;
      }

      VertexLibrarySelection selection;
      selection.response.models.push_back({model_index, 1.0});
      selection.response.conductor_references = model.conductor_references;
      selection.response.normalized_distance = normalized_error;
      selection.first_arm = first;
      best = std::move(selection);
      best_error = normalized_error;
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
    case LibraryTopology::PARALLEL_EDGE_CLUSTER:
      return "parallel-edge cluster";
    case LibraryTopology::SPATIAL_EDGE_CLUSTER:
      return "spatial edge cluster";
    case LibraryTopology::CONVEX_CORNER:
      return "convex corner";
    case LibraryTopology::CONCAVE_CORNER:
      return "concave corner";
    case LibraryTopology::ENDPOINT:
      return "endpoint";
    case LibraryTopology::JUNCTION:
      return "junction";
  }
  return "unknown";
}

std::string TopologyIdentifier(LibraryTopology topology)
{
  switch (topology)
  {
    case LibraryTopology::ISOLATED_EDGE:
      return "IsolatedEdge";
    case LibraryTopology::SAME_CONDUCTOR_GAP:
      return "SameConductorGap";
    case LibraryTopology::DIFFERENT_CONDUCTOR_GAP:
      return "DifferentConductorGap";
    case LibraryTopology::SAME_CONDUCTOR_STRIP:
      return "SameConductorStrip";
    case LibraryTopology::PARALLEL_EDGE_CLUSTER:
      return "ParallelEdgeCluster";
    case LibraryTopology::SPATIAL_EDGE_CLUSTER:
      return "SpatialEdgeCluster";
    case LibraryTopology::CONVEX_CORNER:
      return "ConvexCorner";
    case LibraryTopology::CONCAVE_CORNER:
      return "ConcaveCorner";
    case LibraryTopology::ENDPOINT:
      return "Endpoint";
    case LibraryTopology::JUNCTION:
      return "Junction";
  }
  return "Unknown";
}

std::string BoundaryConditionName(MetalBoundaryConditionType type)
{
  switch (type)
  {
    case MetalBoundaryConditionType::PEC:
      return "PEC";
    case MetalBoundaryConditionType::CONDUCTIVITY:
      return "Conductivity";
    case MetalBoundaryConditionType::IMPEDANCE:
      return "Impedance";
    case MetalBoundaryConditionType::RATIONAL_IMPEDANCE:
      return "RationalImpedance";
  }
  return "Unknown";
}

class AutomaticResponseRequirements
{
private:
  struct Aggregate
  {
    nlohmann::json requirement;
    int count = 0;
    double total_edge_length = 0.0;
  };

  std::map<std::string, Aggregate> requirements;
  std::string library_path;
  std::string library_name;
  double matching_radius = 0.0;
  double coordinate_scale = 1.0;
  nlohmann::json statistics;
  const Units &units;
  bool nondimensionalized = false;

  std::vector<double>
  DimensionalizeRationalCoefficients(const std::vector<double> &coefficients,
                                     bool numerator) const
  {
    if (!nondimensionalized)
    {
      return coefficients;
    }
    auto result = coefficients;
    const double impedance_scale =
        numerator ? units.GetScaleFactor<Units::ValueType::IMPEDANCE>() : 1.0;
    const double time_scale = 1.0e-9 * units.GetScaleFactor<Units::ValueType::TIME>();
    for (std::size_t i = 0; i < result.size(); i++)
    {
      const int degree = static_cast<int>(result.size() - 1 - i);
      result[i] *= impedance_scale * std::pow(time_scale, degree);
    }
    return result;
  }

  nlohmann::json BoundaryCondition(const MetalBoundaryLaw &law) const
  {
    MFEM_VERIFY(law.parameters_verified,
                "Surface-response preflight cannot export an unverified metal boundary "
                "law!");
    nlohmann::json result = {{"Type", BoundaryConditionName(law.type)}};
    switch (law.type)
    {
      case MetalBoundaryConditionType::PEC:
        break;
      case MetalBoundaryConditionType::CONDUCTIVITY:
        MFEM_VERIFY(law.parameters.size() == 3,
                    "Invalid conductivity boundary-law parameter count!");
        result["Conductivity"] =
            nondimensionalized
                ? units.Dimensionalize<Units::ValueType::CONDUCTIVITY>(law.parameters[0])
                : law.parameters[0];
        result["Permeability"] = law.parameters[1];
        result["Thickness"] = nondimensionalized
                                  ? law.parameters[2] * units.GetMeshLengthRelativeScale()
                                  : law.parameters[2];
        // The matcher stores only effective thickness, including the external-surface
        // factor. Export its canonical equivalent instead of inventing lost provenance.
        result["External"] = false;
        break;
      case MetalBoundaryConditionType::IMPEDANCE:
        MFEM_VERIFY(law.parameters.size() == 3,
                    "Invalid impedance boundary-law parameter count!");
        result["Rs"] =
            nondimensionalized
                ? units.Dimensionalize<Units::ValueType::IMPEDANCE>(law.parameters[0])
                : law.parameters[0];
        result["Ls"] =
            nondimensionalized
                ? units.Dimensionalize<Units::ValueType::INDUCTANCE>(law.parameters[1])
                : law.parameters[1];
        result["Cs"] =
            nondimensionalized
                ? units.Dimensionalize<Units::ValueType::CAPACITANCE>(law.parameters[2])
                : law.parameters[2];
        break;
      case MetalBoundaryConditionType::RATIONAL_IMPEDANCE:
        MFEM_VERIFY(!law.numerator.empty() && !law.denominator.empty(),
                    "Invalid rational-impedance boundary-law coefficients!");
        result["Numerator"] = DimensionalizeRationalCoefficients(law.numerator, true);
        result["Denominator"] = DimensionalizeRationalCoefficients(law.denominator, false);
        break;
    }
    return result;
  }

public:
  AutomaticResponseRequirements(const Units &units, bool nondimensionalized)
    : units(units), nondimensionalized(nondimensionalized)
  {
  }

  nlohmann::json DescribeBoundaryCondition(const MetalBoundaryLaw &law) const
  {
    return BoundaryCondition(law);
  }

  void SetLibrary(const std::string &path, const ProcessLibrary &library, double scale)
  {
    library_path = std::filesystem::absolute(path).lexically_normal().string();
    library_name = library.name;
    matching_radius = library.matching_radius * scale;
    coordinate_scale = scale;
  }

  double ScaleLength(double value) const
  {
    const double scaled = value * coordinate_scale;
    const double tolerance =
        std::max(1.0e-10 * matching_radius, 64.0 * std::numeric_limits<double>::epsilon());
    const double step = std::pow(10.0, std::floor(std::log10(tolerance)));
    const double snapped = std::round(scaled / step) * step;
    return snapped == 0.0 ? 0.0 : snapped;
  }

  double UnscaleLength(double value) const { return value / coordinate_scale; }

  double SnapDirection(double value) const
  {
    const double snapped = std::round(value * 1.0e12) * 1.0e-12;
    return snapped == 0.0 ? 0.0 : snapped;
  }

  double SnapAngleDegrees(double value) const
  {
    const double snapped = std::round(value * 1.0e8) * 1.0e-8;
    return snapped == 0.0 ? 0.0 : snapped;
  }

  void SetStatistics(nlohmann::json value) { statistics = std::move(value); }

  void Add(int dimension, LibraryTopology topology,
           const std::map<int, std::map<InterfaceDielectric, int>> &targets_by_slot,
           const MetalBoundaryLaw &boundary_condition, const nlohmann::json &geometry,
           const ProcessLibrary &library, const LibrarySelection *selection, double length,
           const std::string &reason = {})
  {
    nlohmann::json interfaces = nlohmann::json::array();
    for (const auto &[slot, targets] : targets_by_slot)
    {
      for (const auto &[type, target] : targets)
      {
        interfaces.push_back(
            {{"Slot", slot}, {"Type", ToString(type)}, {"Target", target}});
      }
    }

    nlohmann::json requirement = {
        {"Dimension", dimension},
        {"Topology", TopologyIdentifier(topology)},
        {"Status",
         selection ? (selection->IsInterpolated() ? "Interpolated" : "Exact") : "Missing"},
        {"Geometry", geometry},
        {"Interfaces", interfaces},
        {"BoundaryCondition", BoundaryCondition(boundary_condition)}};
    if (selection)
    {
      nlohmann::json models = nlohmann::json::array();
      for (const auto &weighted_model : selection->models)
      {
        const auto &model = library.models[weighted_model.index];
        models.push_back({{"Name", model.name},
                          {"Topology", TopologyIdentifier(model.topology)},
                          {"Weight", SnapDirection(weighted_model.weight)}});
      }
      requirement["SelectedModels"] = std::move(models);
      requirement["NormalizedLibraryDistance"] =
          SnapDirection(selection->normalized_distance);
    }
    if (!reason.empty())
    {
      requirement["Reason"] = reason;
    }

    // Lengths use a tolerance-scaled canonical representation, while angles and
    // directions come from the same production classifier used for model selection.
    const std::string key = requirement.dump();
    auto [it, inserted] = requirements.emplace(key, Aggregate{requirement, 0, 0.0});
    (void)inserted;
    it->second.count++;
    it->second.total_edge_length += ScaleLength(length);
  }

  void Add(int dimension, LibraryTopology topology,
           const std::map<InterfaceDielectric, int> &targets,
           const MetalBoundaryLaw &boundary_condition, const nlohmann::json &geometry,
           const ProcessLibrary &library, const LibrarySelection *selection, double length,
           const std::string &reason = {})
  {
    Add(dimension, topology, {{0, targets}}, boundary_condition, geometry, library,
        selection, length, reason);
  }

  nlohmann::json Build() const
  {
    nlohmann::json entries = nlohmann::json::array();
    std::map<std::string, int> counts = {{"Exact", 0}, {"Interpolated", 0}, {"Missing", 0}};
    std::map<std::string, double> lengths = {
        {"Exact", 0.0}, {"Interpolated", 0.0}, {"Missing", 0.0}};
    for (const auto &[key, aggregate] : requirements)
    {
      (void)key;
      auto entry = aggregate.requirement;
      entry["Count"] = aggregate.count;
      entry["TotalEdgeLength"] = aggregate.total_edge_length;
      counts[entry["Status"].get<std::string>()] += aggregate.count;
      lengths[entry["Status"].get<std::string>()] += aggregate.total_edge_length;
      entries.push_back(std::move(entry));
    }
    nlohmann::json result = {
        {"Version", 1},
        {"Complete", counts["Missing"] == 0},
        {"Library",
         {{"Path", library_path},
          {"Name", library_name},
          {"MatchingRadius", matching_radius}}},
        {"LengthUnit", "mesh"},
        {"Summary", {{"Counts", counts}, {"TotalEdgeLengths", lengths}}},
        {"Requirements", std::move(entries)}};
    if (!statistics.is_null())
    {
      result["Statistics"] = statistics;
    }
    return result;
  }
};

ResponseCorrectionData BuildAutomaticResponseData2D(
    const IoData &iodata, const mfem::ParMesh &mesh, const MaterialOperator &mat_op,
    const ResponseCorrectionData &request, bool pec_attribute_conductors = false,
    AutomaticResponseDiagnostics *diagnostics = nullptr,
    AutomaticResponseRequirements *requirements = nullptr,
    AutomaticResponseStatistics *statistics = nullptr)
{
  MFEM_VERIFY(mesh.Dimension() == 2 && mesh.SpaceDimension() == 2,
              "Automatic two-dimensional fabrication-process response matching requires "
              "a two-dimensional mesh!");
  const double coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
  const auto library =
      ReadProcessLibrary(request.library, iodata.units, iodata.InputsNondimensionalized(),
                         requirements != nullptr);
  if (requirements)
  {
    requirements->SetLibrary(request.library, library, coordinate_scale);
  }
  if (diagnostics)
  {
    diagnostics->matching_radius = library.matching_radius;
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      const int attribute = mesh.GetAttribute(element);
      diagnostics->minimum_wave_speed =
          std::min(diagnostics->minimum_wave_speed, mat_op.GetLightSpeedMin(attribute));
    }
    Mpi::GlobalMin(1, &diagnostics->minimum_wave_speed, mesh.GetComm());
    MFEM_VERIFY(std::isfinite(diagnostics->minimum_wave_speed) &&
                    diagnostics->minimum_wave_speed > 0.0,
                "Unable to determine a positive wave speed for Maxwell surface-response "
                "confidence diagnostics!");
  }

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
  if (statistics)
  {
    statistics->target_groups = groups_by_attributes.size();
  }
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
  if (!target_filter.empty())
  {
    MFEM_VERIFY(found == target_filter,
                "One or more response-correction TargetInterfaces is missing, untyped, "
                "or does not configure edge-distance postprocessing!");
  }
  ValidateLibraryInterfaceLayers(library, iodata.boundaries.postpro.dielectric, found,
                                 coordinate_scale);

  ResponseCorrectionData result;
  result.unmatched_policy = request.unmatched_policy;
  int next_model_index = 1;
  int matched_clusters = 0;
  int matched_edges = 0;
  int interpolated_paired_clusters = 0;
  int unmatched_clusters = 0;
  int next_interpolation_group = 1;
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

    const auto sites =
        ExtractEdgeSites(mesh, iodata.boundaries, group, pec_attribute_conductors);
    if (statistics)
    {
      statistics->edge_sites_2d += sites.size();
    }
    MFEM_VERIFY(!sites.empty(),
                "Automatic response matching found no physical metal edges for target "
                "interface group!");
    if (diagnostics)
    {
      diagnostics->selected_length += static_cast<double>(sites.size());
      for (const auto &[type, target] : group.targets)
      {
        (void)type;
        diagnostics->selected_length_by_interface[target] +=
            static_cast<double>(sites.size());
      }
    }

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
    int group_interpolated_paired_clusters = 0;
    double group_maximum_library_distance = 0.0;
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
      std::optional<LibrarySelection> model_selection;
      const auto boundary_condition = sites[cluster.front()].boundary_condition;
      auto ClusterGeometry = [&]()
      {
        nlohmann::json geometry = {{"EdgeCount", cluster.size()}};
        if (cluster.size() > 1)
        {
          geometry["Separation"] =
              requirements ? requirements->ScaleLength(
                                 Distance(sites[cluster[0]].point, sites[cluster[1]].point))
                           : 0.0;
        }
        if (cluster.size() > 2)
        {
          const auto &reference = sites[cluster.front()];
          std::map<int, int> conductor_ids;
          nlohmann::json edges = nlohmann::json::array();
          for (const std::size_t index : cluster)
          {
            const auto &site = sites[index];
            const Point2D offset = {site.point[0] - reference.point[0],
                                    site.point[1] - reference.point[1]};
            auto [conductor, inserted] =
                conductor_ids.emplace(site.conductor, conductor_ids.size() + 1);
            (void)inserted;
            edges.push_back(
                {{"Offset",
                  {requirements->ScaleLength(Dot(offset, reference.axis_u)),
                   requirements->ScaleLength(Dot(offset, reference.axis_v))}},
                 {"GapDirection",
                  {requirements->SnapDirection(Dot(site.axis_u, reference.axis_u)),
                   requirements->SnapDirection(Dot(site.axis_u, reference.axis_v))}},
                 {"Conductor", conductor->second}});
          }
          geometry["Edges"] = std::move(edges);
        }
        return geometry;
      };
      if (std::any_of(cluster.begin(), cluster.end(),
                      [&](std::size_t index)
                      {
                        return !SameBoundaryLaw(sites[index].boundary_condition,
                                                boundary_condition);
                      }))
      {
        group_matched = false;
        Mpi::Warning(
            "Nearby two-dimensional metal edges use different boundary conditions; "
            "correction is disabled for this interface group!\n");
        if (requirements)
        {
          requirements->Add(2, LibraryTopology::SPATIAL_EDGE_CLUSTER, group.targets,
                            boundary_condition, ClusterGeometry(), library, nullptr, 0.0,
                            "Nearby edges use different metal boundary conditions");
          unmatched_clusters++;
          continue;
        }
        unmatched_clusters++;
        break;
      }
      if (cluster.size() == 1)
      {
        const auto &edge = sites[cluster.front()];
        patch.origin = {edge.point[0], edge.point[1], 0.0};
        patch.axis_u = {edge.axis_u[0], edge.axis_u[1], 0.0};
        patch.axis_v = {edge.axis_v[0], edge.axis_v[1], 0.0};
        patch.maxwell_reference_is_pec =
            edge.boundary_condition.type == MetalBoundaryConditionType::PEC;
        if (diagnostics && !patch.maxwell_reference_is_pec)
        {
          patch.maxwell_conductor_anchors = {patch.origin};
        }
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
          if (requirements)
          {
            requirements->Add(2, LibraryTopology::SPATIAL_EDGE_CLUSTER, group.targets,
                              boundary_condition, ClusterGeometry(), library, nullptr, 0.0,
                              "No canonical paired-edge topology");
            unmatched_clusters++;
            continue;
          }
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
        patch.maxwell_reference_is_pec =
            boundary_condition.type == MetalBoundaryConditionType::PEC;
        if (diagnostics && !patch.maxwell_reference_is_pec)
        {
          patch.maxwell_conductor_anchors = {
              std::array<double, 3>{first.point[0], first.point[1], 0.0}};
        }
      }
      else
      {
        const auto cluster_selection =
            FindParallelClusterLibraryModel(library, sites, cluster);
        if (!cluster_selection)
        {
          group_matched = false;
          Mpi::Warning("Fabrication-process response library \"{}\" has no matching "
                       "ParallelEdgeCluster model for {} nearby metal edges; correction is "
                       "disabled for this interface group!\n",
                       library.name, cluster.size());
          if (requirements)
          {
            requirements->Add(2, LibraryTopology::PARALLEL_EDGE_CLUSTER, group.targets,
                              boundary_condition, ClusterGeometry(), library, nullptr, 0.0,
                              "No compatible parallel-edge cluster model");
            unmatched_clusters++;
            continue;
          }
          unmatched_clusters++;
          break;
        }
        model_selection = cluster_selection->response;
        const auto &first = sites[cluster_selection->ordered_edges.front()];
        patch.origin = {first.point[0], first.point[1], 0.0};
        patch.axis_u = {cluster_selection->axis_u[0], cluster_selection->axis_u[1], 0.0};
        patch.axis_v = {cluster_selection->axis_v[0], cluster_selection->axis_v[1], 0.0};
        patch.maxwell_reference_is_pec =
            std::all_of(cluster.begin(), cluster.end(),
                        [&](std::size_t index)
                        {
                          return sites[index].boundary_condition.type ==
                                 MetalBoundaryConditionType::PEC;
                        });
        patch.conductor_references = model_selection->conductor_references;
        if (diagnostics)
        {
          for (const std::size_t index : cluster_selection->reference_edges)
          {
            const auto &point = sites[index].point;
            patch.maxwell_conductor_anchors.push_back({point[0], point[1], 0.0});
          }
        }
      }

      if (!model_selection)
      {
        model_selection =
            FindLibraryModel(library, topology, separation, boundary_condition);
      }
      if (!model_selection)
      {
        group_matched = false;
        Mpi::Warning(
            "Fabrication-process response library \"{}\" has no {} model at separation "
            "{:.6e} mesh units; correction is disabled for this interface group!\n",
            library.name, TopologyName(topology), separation * coordinate_scale);
        if (requirements)
        {
          requirements->Add(2, topology, group.targets, boundary_condition,
                            ClusterGeometry(), library, nullptr, 0.0,
                            "No compatible process-library model");
          unmatched_clusters++;
          continue;
        }
        unmatched_clusters++;
        break;
      }
      if (requirements)
      {
        requirements->Add(
            2, cluster.size() > 2 ? LibraryTopology::PARALLEL_EDGE_CLUSTER : topology,
            group.targets, boundary_condition, ClusterGeometry(), library,
            &*model_selection, 0.0);
      }
      if (patch.conductor_references.empty())
      {
        patch.conductor_references = model_selection->conductor_references;
      }
      if (diagnostics && !patch.maxwell_reference_is_pec &&
          patch.conductor_references.size() == 2 && cluster.size() == 2)
      {
        const auto &second = sites[cluster[1]];
        patch.maxwell_conductor_anchors.push_back(
            std::array<double, 3>{second.point[0], second.point[1], 0.0});
      }
      group_maximum_library_distance =
          std::max(group_maximum_library_distance, model_selection->normalized_distance);
      if (model_selection->IsInterpolated())
      {
        patch.interpolation_group = next_interpolation_group++;
        group_interpolated_paired_clusters++;
      }
      for (const auto &weighted_model : model_selection->models)
      {
        const auto &source = library.models[weighted_model.index];
        auto weighted_patch = patch;
        weighted_patch.weight = weighted_model.weight;
        if (source.coupon_depth > 0.0)
        {
          weighted_patch.weight /= source.coupon_depth;
        }
        pending.push_back({weighted_model.index, std::move(weighted_patch)});
      }
    }

    if (!group_matched)
    {
      if (request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
      {
        MFEM_ABORT("Automatic fabrication-process response matching failed!");
      }
      continue;
    }

    if (diagnostics)
    {
      for (const auto &selection : pending)
      {
        diagnostics->boundary_law_verified &=
            IsBoundaryLawVerified(library.models[selection.library_model]);
      }
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
        model.name = source.name;
        model.topology = TopologyName(source.topology);
        MapLibraryInterfaces(source, {{0, group.targets}}, model);
        result.models.push_back(std::move(model));
      }
      selection.patch.model = model_it->second;
      result.patches.push_back(selection.patch);
    }
    matched_clusters += component_count;
    matched_edges += static_cast<int>(sites.size());
    interpolated_paired_clusters += group_interpolated_paired_clusters;
    if (diagnostics)
    {
      diagnostics->matched_length += static_cast<double>(sites.size());
      diagnostics->maximum_library_distance =
          std::max(diagnostics->maximum_library_distance, group_maximum_library_distance);
      for (const auto &[type, target] : group.targets)
      {
        (void)type;
        diagnostics->matched_length_by_interface[target] +=
            static_cast<double>(sites.size());
      }
    }
  }

  MFEM_VERIFY(requirements || (!result.models.empty() && !result.patches.empty()),
              "Fabrication-process response matching produced no usable correction "
              "patches!");
  Mpi::Print("\nAutomatic fabrication-process response matching:\n"
             " Library: {}\n"
             " Matched edge sites: {:d}\n"
             " Matched clusters: {:d}\n"
             " Interpolated paired clusters: {:d}\n"
             " Unmatched interface groups: {:d}\n",
             library.name, matched_edges, matched_clusters, interpolated_paired_clusters,
             unmatched_clusters);
  return result;
}

bool SegmentsShareVertex(const MetalEdgeSegment &a, const MetalEdgeSegment &b)
{
  return a.vertices[0] == b.vertices[0] || a.vertices[0] == b.vertices[1] ||
         a.vertices[1] == b.vertices[0] || a.vertices[1] == b.vertices[1];
}

struct SegmentClosestApproach
{
  double first = 0.0;
  double second = 0.0;
  double distance_squared = 0.0;
};

SegmentClosestApproach ClosestSegmentApproach(const Point3D &p0, const Point3D &p1,
                                              const Point3D &q0, const Point3D &q1)
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
  return {s, t, Dot(delta, delta)};
}

double SegmentDistanceSquared(const Point3D &p0, const Point3D &p1, const Point3D &q0,
                              const Point3D &q1)
{
  return ClosestSegmentApproach(p0, p1, q0, q1).distance_squared;
}

double PointSegmentDistanceSquared(const Point3D &point, const EdgeSegment3D &segment)
{
  const double distance =
      std::clamp(Dot(Subtract(point, segment.p0), segment.tangent), 0.0, segment.length);
  const Point3D delta = Subtract(point, Interpolate(segment, distance));
  return Dot(delta, delta);
}

Point3D TransformLocalVector(const std::array<Point3D, 3> &axes, const Point3D &local)
{
  Point3D result{};
  for (int d = 0; d < 3; d++)
  {
    result = Add(result, Scale(local[d], axes[d]));
  }
  return result;
}

Point3D TransformLocalPoint(const Point3D &origin, const std::array<Point3D, 3> &axes,
                            const Point3D &local)
{
  return Add(origin, TransformLocalVector(axes, local));
}

std::pair<Point3D, std::array<Point3D, 3>> AlignSpatialFrame(const Point3D &model_point,
                                                             const Point3D &model_gap,
                                                             const Point3D &model_normal,
                                                             const SpatialEdgeSite3D &site)
{
  const Point3D model_tangent = Normalize(Cross(model_gap, model_normal));
  const Point3D site_tangent = Normalize(Cross(site.gap_direction, site.process_normal));
  std::array<Point3D, 3> axes{};
  for (int d = 0; d < 3; d++)
  {
    Point3D local_axis{};
    local_axis[d] = 1.0;
    axes[d] = Add(Scale(Dot(model_gap, local_axis), site.gap_direction),
                  Add(Scale(Dot(model_normal, local_axis), site.process_normal),
                      Scale(Dot(model_tangent, local_axis), site_tangent)));
  }
  return {Subtract(site.point, TransformLocalVector(axes, model_point)), axes};
}

std::optional<SpatialClusterSelection3D> FindSpatialClusterLibraryModel(
    const ProcessLibrary &library, const std::vector<SpatialEdgeSite3D> &sites,
    const std::set<std::size_t> &excluded_models = {},
    const std::function<bool(const SpatialClusterSelection3D &, const LibraryModel &)>
        &matches_plan_view = {})
{
  if (sites.size() < 2)
  {
    return std::nullopt;
  }

  std::optional<SpatialClusterSelection3D> best;
  double best_distance = mfem::infinity();
  for (std::size_t model_index = 0; model_index < library.models.size(); model_index++)
  {
    const auto &model = library.models[model_index];
    if (excluded_models.find(model_index) != excluded_models.end() ||
        model.topology != LibraryTopology::SPATIAL_EDGE_CLUSTER ||
        model.spatial_edges.size() != sites.size())
    {
      continue;
    }
    const double position_tolerance =
        std::max(model.spatial_position_tolerance, 1.0e-10 * library.matching_radius);
    const double angle_tolerance =
        std::max(model.spatial_angle_tolerance, 1.0e-10 * std::acos(-1.0));
    const auto &model_anchor = model.spatial_edges.front();

    for (std::size_t site_anchor_index = 0; site_anchor_index < sites.size();
         site_anchor_index++)
    {
      const auto &site_anchor = sites[site_anchor_index];
      if (!CompatibleBoundaryLaw(model_anchor.boundary_condition,
                                 site_anchor.boundary_condition))
      {
        continue;
      }
      const auto aligned_frame =
          AlignSpatialFrame(model_anchor.point, model_anchor.gap_direction,
                            model_anchor.process_normal, site_anchor);
      const Point3D origin = aligned_frame.first;
      const std::array<Point3D, 3> axes = aligned_frame.second;

      struct Candidate
      {
        std::size_t site = 0;
        double distance = 0.0;
      };
      std::vector<std::vector<Candidate>> candidates(model.spatial_edges.size());
      for (std::size_t edge_index = 0; edge_index < model.spatial_edges.size();
           edge_index++)
      {
        const auto &edge = model.spatial_edges[edge_index];
        const Point3D point = TransformLocalPoint(origin, axes, edge.point);
        const Point3D gap = TransformLocalVector(axes, edge.gap_direction);
        const Point3D normal = TransformLocalVector(axes, edge.process_normal);
        for (std::size_t site_index = 0; site_index < sites.size(); site_index++)
        {
          if (edge_index == 0 && site_index != site_anchor_index)
          {
            continue;
          }
          const auto &site = sites[site_index];
          if (!CompatibleBoundaryLaw(edge.boundary_condition, site.boundary_condition))
          {
            continue;
          }
          const double position_error = Distance(point, site.point);
          const double gap_error =
              std::acos(std::clamp(Dot(gap, site.gap_direction), -1.0, 1.0));
          const double normal_error =
              std::acos(std::clamp(Dot(normal, site.process_normal), -1.0, 1.0));
          const double interval_error =
              std::max(std::abs(edge.interval[0] - site.interval[0]),
                       std::abs(edge.interval[1] - site.interval[1]));
          const double normalized_error = std::max(
              {position_error / position_tolerance, gap_error / angle_tolerance,
               normal_error / angle_tolerance, interval_error / position_tolerance});
          if (normalized_error <= 1.0)
          {
            candidates[edge_index].push_back({site_index, normalized_error});
          }
        }
        if (candidates[edge_index].empty())
        {
          break;
        }
      }
      if (std::any_of(candidates.begin(), candidates.end(),
                      [](const auto &candidate) { return candidate.empty(); }))
      {
        continue;
      }

      std::vector<std::size_t> order(model.spatial_edges.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](std::size_t first, std::size_t second)
                { return candidates[first].size() < candidates[second].size(); });
      std::vector<std::size_t> assignment(model.spatial_edges.size(),
                                          std::numeric_limits<std::size_t>::max());
      std::vector<bool> used_site(sites.size(), false);
      std::map<int, int> model_to_conductor;
      std::map<int, int> conductor_to_model;
      std::map<int, std::map<InterfaceDielectric, int>> targets_by_slot;
      std::map<std::vector<std::pair<InterfaceDielectric, int>>, int> slot_by_targets;
      std::function<void(std::size_t, double)> Match =
          [&](std::size_t depth, double normalized_distance)
      {
        if (normalized_distance >= best_distance)
        {
          return;
        }
        if (depth == order.size())
        {
          if (!MatchLibraryInterfaces(model, targets_by_slot))
          {
            return;
          }
          SpatialClusterSelection3D selection;
          selection.response.models.push_back({model_index, 1.0});
          selection.response.conductor_references = model.conductor_references;
          selection.response.normalized_distance = normalized_distance;
          selection.sites = sites;
          selection.model_to_site = assignment;
          selection.targets_by_slot = targets_by_slot;
          selection.origin = origin;
          selection.axes = axes;
          if (matches_plan_view && !matches_plan_view(selection, model))
          {
            return;
          }
          best = std::move(selection);
          best_distance = normalized_distance;
          return;
        }

        const std::size_t edge_index = order[depth];
        const int model_conductor = model.spatial_edges[edge_index].conductor;
        const int interface_slot = model.spatial_edges[edge_index].interface_slot;
        for (const auto &candidate : candidates[edge_index])
        {
          if (used_site[candidate.site])
          {
            continue;
          }
          const int site_conductor = sites[candidate.site].conductor;
          const auto &site_targets = sites[candidate.site].targets;
          const std::vector<std::pair<InterfaceDielectric, int>> target_signature(
              site_targets.begin(), site_targets.end());
          const auto mapped_site = model_to_conductor.find(model_conductor);
          const auto mapped_model = conductor_to_model.find(site_conductor);
          const auto mapped_targets = targets_by_slot.find(interface_slot);
          const auto mapped_slot = slot_by_targets.find(target_signature);
          if ((mapped_site != model_to_conductor.end() &&
               mapped_site->second != site_conductor) ||
              (mapped_model != conductor_to_model.end() &&
               mapped_model->second != model_conductor) ||
              (mapped_targets != targets_by_slot.end() &&
               mapped_targets->second != site_targets) ||
              (mapped_slot != slot_by_targets.end() &&
               mapped_slot->second != interface_slot))
          {
            continue;
          }
          const bool new_model = mapped_site == model_to_conductor.end();
          const bool new_site = mapped_model == conductor_to_model.end();
          const bool new_targets = mapped_targets == targets_by_slot.end();
          const bool new_slot = mapped_slot == slot_by_targets.end();
          if (new_model)
          {
            model_to_conductor.emplace(model_conductor, site_conductor);
          }
          if (new_site)
          {
            conductor_to_model.emplace(site_conductor, model_conductor);
          }
          if (new_targets)
          {
            targets_by_slot.emplace(interface_slot, site_targets);
          }
          if (new_slot)
          {
            slot_by_targets.emplace(target_signature, interface_slot);
          }
          assignment[edge_index] = candidate.site;
          used_site[candidate.site] = true;
          Match(depth + 1, std::max(normalized_distance, candidate.distance));
          used_site[candidate.site] = false;
          assignment[edge_index] = std::numeric_limits<std::size_t>::max();
          if (new_model)
          {
            model_to_conductor.erase(model_conductor);
          }
          if (new_site)
          {
            conductor_to_model.erase(site_conductor);
          }
          if (new_targets)
          {
            targets_by_slot.erase(interface_slot);
          }
          if (new_slot)
          {
            slot_by_targets.erase(target_signature);
          }
        }
      };
      Match(0, 0.0);
    }
  }
  return best;
}

ResponseCorrectionData
BuildAutomaticResponseData3D(const IoData &iodata, const mfem::ParMesh &mesh,
                             const MaterialOperator &mat_op,
                             const ResponseCorrectionData &request, bool maxwell,
                             AutomaticResponseDiagnostics *diagnostics = nullptr,
                             AutomaticResponseRequirements *requirements = nullptr,
                             AutomaticResponseStatistics *statistics = nullptr)
{
  MFEM_VERIFY(mesh.Dimension() == 3 && mesh.SpaceDimension() == 3,
              "Automatic three-dimensional fabrication-process response matching "
              "requires a three-dimensional mesh!");
  const double coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
  const auto library =
      ReadProcessLibrary(request.library, iodata.units, iodata.InputsNondimensionalized(),
                         requirements != nullptr);
  if (requirements)
  {
    requirements->SetLibrary(request.library, library, coordinate_scale);
  }
  if (diagnostics)
  {
    diagnostics->matching_radius = library.matching_radius;
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      const int attribute = mesh.GetAttribute(element);
      diagnostics->minimum_wave_speed =
          std::min(diagnostics->minimum_wave_speed, mat_op.GetLightSpeedMin(attribute));
    }
    Mpi::GlobalMin(1, &diagnostics->minimum_wave_speed, mesh.GetComm());
    MFEM_VERIFY(std::isfinite(diagnostics->minimum_wave_speed) &&
                    diagnostics->minimum_wave_speed > 0.0,
                "Unable to determine a positive wave speed for Maxwell surface-response "
                "confidence diagnostics!");
  }
  MetalSurfaceExtraction surface;
  surface.classify_components =
      requirements ||
      std::any_of(library.models.begin(), library.models.end(), [](const auto &model)
                  { return model.topology == LibraryTopology::SPATIAL_EDGE_CLUSTER; });
  surface.retain_faces =
      requirements ||
      std::any_of(library.models.begin(), library.models.end(),
                  [](const auto &model) { return model.plan_view_boundary.has_value(); });
  const auto geometry = ExtractMetalEdgeGeometry(mesh, iodata.boundaries, surface);
  MFEM_VERIFY(!geometry.Empty(),
              "Fabrication-process response matching found no metal perimeter!");
  if (statistics)
  {
    statistics->metal_vertices = geometry.vertices.size();
    statistics->metal_segments = geometry.segments.size();
    statistics->metal_components = geometry.metal_components;
    statistics->physical_components = geometry.physical_components;
    statistics->physical_chains = geometry.physical_chains;
    statistics->surface_faces_local = geometry.surface_faces.size();
  }

  std::set<int> target_filter(request.target_interfaces.begin(),
                              request.target_interfaces.end());
  struct TargetSelection
  {
    InterfaceDielectric type;
    int index;
    std::vector<std::size_t> segments;
    std::optional<Point3D> process_normal;
  };
  std::vector<TargetSelection> selections;
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
    ExcludeMetalEdgeSegmentIndices(mesh, geometry, dielectric.edge_exclude_attributes,
                                   segment_indices);
    selections.push_back({dielectric.type, index, std::move(segment_indices),
                          dielectric.edge_frame_normal ? std::optional<Point3D>(Normalize(
                                                             *dielectric.edge_frame_normal))
                                                       : std::nullopt});
  }
  MFEM_VERIFY(!selections.empty(),
              "Fabrication-process response matching found no target interfaces!");
  std::set<int> found;
  for (const auto &selection : selections)
  {
    found.insert(selection.index);
  }
  if (!target_filter.empty())
  {
    MFEM_VERIFY(found == target_filter,
                "One or more response-correction TargetInterfaces is missing, untyped, "
                "or does not configure edge-distance postprocessing!");
  }
  ValidateLibraryInterfaceLayers(library, iodata.boundaries.postpro.dielectric, found,
                                 coordinate_scale);

  // Partition the selected perimeter by the interface types available on each physical
  // segment. Existing SA attributes can legitimately omit a port footprint while the
  // corresponding MS and MA attributes continue along the PEC perimeter. A target
  // interface may therefore span several groups, but every physical segment belongs to
  // exactly one group and receives at most one target of each interface type.
  std::map<std::size_t, std::map<InterfaceDielectric, int>> targets_by_segment;
  std::map<std::pair<InterfaceDielectric, int>, std::optional<Point3D>> normals_by_target;
  for (const auto &selection : selections)
  {
    normals_by_target.emplace(std::make_pair(selection.type, selection.index),
                              selection.process_normal);
    for (const std::size_t segment : selection.segments)
    {
      const auto [target, inserted] =
          targets_by_segment[segment].emplace(selection.type, selection.index);
      MFEM_VERIFY(inserted || target->second == selection.index,
                  "Automatic response matching found multiple target interfaces of type "
                      << ToString(selection.type)
                      << " on the same three-dimensional metal-perimeter segment!");
    }
  }
  using TargetSignature = std::vector<std::pair<InterfaceDielectric, int>>;
  std::map<TargetSignature, EdgeGroup3D> groups_by_targets;
  for (const auto &[segment, targets] : targets_by_segment)
  {
    TargetSignature signature(targets.begin(), targets.end());
    auto &group = groups_by_targets[signature];
    group.segments.push_back(segment);
    group.targets = targets;
    group.matching_radius = library.matching_radius;
    for (const auto &[type, index] : targets)
    {
      const auto &normal = normals_by_target.at({type, index});
      if (!normal)
      {
        continue;
      }
      if (!group.process_normal)
      {
        group.process_normal = normal;
      }
      else
      {
        MFEM_VERIFY(Dot(*group.process_normal, *normal) > 1.0 - 1.0e-10,
                    "Target interfaces on the same metal-perimeter segment must use the "
                    "same EdgeFrameNormal!");
      }
    }
  }

  if (statistics)
  {
    statistics->target_groups = groups_by_targets.size();
  }
  const auto &quadrature =
      mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * std::max(1, iodata.solver.order));
  ResponseCorrectionData result;
  result.unmatched_policy = request.unmatched_policy;
  int next_model_index = 1;
  int matched_intervals = 0;
  int matched_segments = 0;
  int unmatched_groups = 0;
  int unmatched_rounded_corners = 0;
  int interpolated_paired_intervals = 0;
  int interpolated_rounded_corners = 0;
  int nonregular_vertices = 0;
  int matched_corner_patches = 0;
  int matched_endpoint_patches = 0;
  int matched_junction_patches = 0;
  int matched_spatial_cluster_patches = 0;
  int matched_nonregular_vertices = 0;

  std::map<std::size_t, EdgeSegment3D> segment_cache;
  auto BuildSegments =
      [&](const EdgeGroup3D &group, const std::vector<std::size_t> &geometry_indices)
  {
    std::vector<std::size_t> missing;
    for (const auto geometry_index : geometry_indices)
    {
      if (segment_cache.find(geometry_index) == segment_cache.end())
      {
        missing.push_back(geometry_index);
      }
    }
    if (!missing.empty())
    {
      const auto process_normals = BuildMetalEdgeProcessNormals(
          mesh, geometry, missing, [&](int attribute)
          { return mat_op.GetLightSpeedMax(attribute); }, group.process_normal);
      const auto gap_directions =
          BuildMetalEdgeGapDirections(mesh, geometry, missing, process_normals);
      for (std::size_t i = 0; i < missing.size(); i++)
      {
        const std::size_t geometry_index = missing[i];
        const auto &source = geometry.segments[geometry_index];
        EdgeSegment3D segment;
        segment.geometry_index = geometry_index;
        segment.p0 = geometry.vertices[source.vertices[0]].coordinate;
        segment.p1 = geometry.vertices[source.vertices[1]].coordinate;
        segment.length = Distance(segment.p0, segment.p1);
        segment.tangent = Normalize(Subtract(segment.p1, segment.p0));
        segment.axis_u = gap_directions[i];
        segment.axis_v = process_normals[i];
        segment.targets = group.targets;
        segment.metal_component = source.metal_component;
        if (maxwell)
        {
          MFEM_VERIFY(!source.conditions.empty(),
                      "Unable to determine the Maxwell metal boundary condition for an "
                      "automatically detected edge!");
          const auto boundary_condition =
              GetBoundaryConditionLaw(iodata.boundaries, source.conditions.front());
          MFEM_VERIFY(
              std::all_of(source.conditions.begin(), source.conditions.end(),
                          [&](const auto &condition)
                          {
                            return SameBoundaryLaw(
                                GetBoundaryConditionLaw(iodata.boundaries, condition),
                                boundary_condition);
                          }),
              "A Maxwell target edge cannot mix distinct metal boundary conditions!");
          MFEM_VERIFY(source.component >= 0,
                      "Unable to determine connected metal ownership for a Maxwell edge!");
          segment.conductor =
              source.metal_component >= 0 ? source.metal_component : source.component;
          segment.boundary_condition = boundary_condition;
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
          MFEM_VERIFY(
              conductors.size() == 1,
              "Unable to assign an automatically detected three-dimensional metal edge "
              "to exactly one electrostatic conductor!");
          segment.conductor = *conductors.begin();
          segment.boundary_condition = MetalBoundaryLaw{};
        }
        segment_cache.emplace(geometry_index, std::move(segment));
      }
    }

    std::vector<EdgeSegment3D> segments;
    segments.reserve(geometry_indices.size());
    for (const auto geometry_index : geometry_indices)
    {
      const auto &segment = segment_cache.at(geometry_index);
      MFEM_VERIFY(segment.targets == group.targets,
                  "Cached edge frame reused with incompatible target interfaces!");
      segments.push_back(segment);
    }
    return segments;
  };
  auto SpatialInterval = [&](const std::vector<EdgeSegment3D> &segments,
                             std::size_t segment_index, double segment_distance,
                             double radius)
  {
    const auto &segment = segments[segment_index];
    const int physical_chain = geometry.segments[segment.geometry_index].physical_chain;
    const Point3D point = Interpolate(segment, segment_distance);
    const Point3D tangent = Normalize(Cross(segment.axis_u, segment.axis_v));
    double begin = 0.0;
    double end = 0.0;
    for (const auto &candidate : segments)
    {
      if (geometry.segments[candidate.geometry_index].physical_chain != physical_chain)
      {
        continue;
      }
      begin = std::min({begin, Dot(Subtract(candidate.p0, point), tangent),
                        Dot(Subtract(candidate.p1, point), tangent)});
      end = std::max({end, Dot(Subtract(candidate.p0, point), tangent),
                      Dot(Subtract(candidate.p1, point), tangent)});
    }
    return std::array<double, 2>{std::max(-radius, begin), std::min(radius, end)};
  };

  auto SpatialGeometry = [&](const SpatialClusterSelection3D &selection)
  {
    nlohmann::json edges = nlohmann::json::array();
    std::map<int, int> conductor_ids;
    const auto &model = library.models[selection.response.models.empty()
                                           ? 0
                                           : selection.response.models.front().index];
    for (std::size_t model_edge = 0; model_edge < selection.model_to_site.size();
         model_edge++)
    {
      const auto &site = selection.sites[selection.model_to_site[model_edge]];
      auto [conductor, inserted] =
          conductor_ids.emplace(site.conductor, conductor_ids.size() + 1);
      (void)inserted;
      const auto relative = Subtract(site.point, selection.origin);
      nlohmann::json edge = {
          {"Point",
           {requirements->ScaleLength(Dot(relative, selection.axes[0])),
            requirements->ScaleLength(Dot(relative, selection.axes[1])),
            requirements->ScaleLength(Dot(relative, selection.axes[2]))}},
          {"GapDirection",
           {requirements->SnapDirection(Dot(site.gap_direction, selection.axes[0])),
            requirements->SnapDirection(Dot(site.gap_direction, selection.axes[1])),
            requirements->SnapDirection(Dot(site.gap_direction, selection.axes[2]))}},
          {"ProcessNormal",
           {requirements->SnapDirection(Dot(site.process_normal, selection.axes[0])),
            requirements->SnapDirection(Dot(site.process_normal, selection.axes[1])),
            requirements->SnapDirection(Dot(site.process_normal, selection.axes[2]))}},
          {"Conductor", conductor->second},
          {"BoundaryCondition",
           requirements->DescribeBoundaryCondition(site.boundary_condition)}};
      if (!selection.response.models.empty() && model_edge < model.spatial_edges.size())
      {
        edge["Interval"] = {
            requirements->ScaleLength(model.spatial_edges[model_edge].interval[0]),
            requirements->ScaleLength(model.spatial_edges[model_edge].interval[1])};
        edge["InterfaceSlot"] = model.spatial_edges[model_edge].interface_slot;
      }
      edges.push_back(std::move(edge));
    }
    return nlohmann::json{{"EdgeCount", edges.size()}, {"Edges", std::move(edges)}};
  };
  auto SpatialLength = [&](const SpatialClusterSelection3D &selection)
  {
    if (selection.response.models.empty())
    {
      return 0.0;
    }
    const auto &model = library.models[selection.response.models.front().index];
    return std::accumulate(model.spatial_edges.begin(), model.spatial_edges.end(), 0.0,
                           [](double length, const auto &edge)
                           { return length + edge.interval[1] - edge.interval[0]; });
  };
  auto ClipPlanViewPolygon =
      [](std::vector<Point3D> polygon, int axis, double bound, bool keep_greater)
  {
    if (polygon.empty())
    {
      return polygon;
    }
    std::vector<Point3D> clipped;
    clipped.reserve(polygon.size() + 2);
    auto IsInside = [&](const Point3D &point)
    { return keep_greater ? point[axis] >= bound : point[axis] <= bound; };
    Point3D previous = polygon.back();
    bool previous_inside = IsInside(previous);
    for (const auto &current : polygon)
    {
      const bool current_inside = IsInside(current);
      if (current_inside != previous_inside)
      {
        const double denominator = current[axis] - previous[axis];
        MFEM_ASSERT(std::abs(denominator) > 0.0,
                    "Invalid plan-view clipping intersection!");
        const double fraction = (bound - previous[axis]) / denominator;
        clipped.push_back(Add(previous, Scale(fraction, Subtract(current, previous))));
      }
      if (current_inside)
      {
        clipped.push_back(current);
      }
      previous = current;
      previous_inside = current_inside;
    }
    return clipped;
  };
  std::map<int, std::vector<const MetalSurfaceFace *>> surface_faces_by_component;
  for (const auto &face : geometry.surface_faces)
  {
    surface_faces_by_component[face.component].push_back(&face);
  }

  auto GatherPlanViewFacets = [&](const std::vector<SpatialEdgeSite3D> &sites,
                                  const Point3D &origin, const std::array<Point3D, 3> &axes,
                                  const std::map<int, int> &conductor_by_metal_component,
                                  int process_axis = 1)
  {
    if (statistics)
    {
      statistics->mask_gather_calls++;
    }
    MFEM_ASSERT(process_axis >= 0 && process_axis < 3,
                "Plan-view facet extraction requires a valid process axis!");
    const std::array<int, 2> plan_axes = process_axis == 0   ? std::array<int, 2>{1, 2}
                                         : process_axis == 1 ? std::array<int, 2>{0, 2}
                                                             : std::array<int, 2>{0, 1};
    std::array<double, 2> lower = {mfem::infinity(), mfem::infinity()};
    std::array<double, 2> upper = {-mfem::infinity(), -mfem::infinity()};
    std::map<int, std::vector<double>> planes_by_metal_component;
    for (const auto &site : sites)
    {
      Point3D point{}, gap{}, normal{};
      const Point3D relative = Subtract(site.point, origin);
      for (int d = 0; d < 3; d++)
      {
        point[d] = Dot(relative, axes[d]);
        gap[d] = Dot(site.gap_direction, axes[d]);
        normal[d] = Dot(site.process_normal, axes[d]);
      }
      const Point3D tangent = Normalize(Cross(gap, normal));
      double begin = site.interval[0];
      double end = site.interval[1];
      const double interval_tolerance = 1.0e-10 * library.matching_radius;
      if (begin <= -library.matching_radius + interval_tolerance)
      {
        begin -= 2.0 * library.matching_radius;
      }
      if (end >= library.matching_radius - interval_tolerance)
      {
        end += 2.0 * library.matching_radius;
      }
      for (const double coordinate : {begin, end})
      {
        const Point3D boundary = Add(point, Scale(coordinate, tangent));
        for (const double side : {-1.0, 1.0})
        {
          const Point3D sample = Add(boundary, Scale(side * library.matching_radius, gap));
          for (int d = 0; d < 2; d++)
          {
            lower[d] = std::min(lower[d], sample[plan_axes[d]]);
            upper[d] = std::max(upper[d], sample[plan_axes[d]]);
          }
        }
      }
      if (site.metal_component >= 0)
      {
        auto &planes = planes_by_metal_component[site.metal_component];
        if (std::none_of(planes.begin(), planes.end(),
                         [&](double plane)
                         {
                           return std::abs(plane - point[process_axis]) <=
                                  1.0e-8 * library.matching_radius;
                         }))
        {
          planes.push_back(point[process_axis]);
        }
      }
    }
    for (int d = 0; d < 2; d++)
    {
      lower[d] -= library.matching_radius;
      upper[d] += library.matching_radius;
    }

    std::vector<PlanViewFacet> local_facets;
    const double tolerance = 1.0e-9 * library.matching_radius;
    for (const auto &[component, conductor] : conductor_by_metal_component)
    {
      const auto component_faces = surface_faces_by_component.find(component);
      if (component_faces == surface_faces_by_component.end())
      {
        continue;
      }
      if (statistics)
      {
        statistics->mask_faces_scanned_local += component_faces->second.size();
      }
      for (const auto *face : component_faces->second)
      {
        std::vector<Point3D> polygon;
        polygon.reserve(face->vertices.size());
        for (const auto &vertex : face->vertices)
        {
          const Point3D relative = Subtract(vertex, origin);
          polygon.push_back(
              {Dot(relative, axes[0]), Dot(relative, axes[1]), Dot(relative, axes[2])});
        }
        const auto &planes = planes_by_metal_component.at(component);
        const auto plane = std::find_if(
            planes.begin(), planes.end(),
            [&](double candidate)
            {
              return std::all_of(
                  polygon.begin(), polygon.end(), [&](const auto &point)
                  { return std::abs(point[process_axis] - candidate) <= tolerance; });
            });
        if (plane == planes.end())
        {
          continue;
        }
        polygon = ClipPlanViewPolygon(std::move(polygon), plan_axes[0], lower[0], true);
        polygon = ClipPlanViewPolygon(std::move(polygon), plan_axes[0], upper[0], false);
        polygon = ClipPlanViewPolygon(std::move(polygon), plan_axes[1], lower[1], true);
        polygon = ClipPlanViewPolygon(std::move(polygon), plan_axes[1], upper[1], false);
        if (polygon.size() < 3)
        {
          continue;
        }
        std::vector<Point3D> unique;
        unique.reserve(polygon.size());
        for (const auto &point : polygon)
        {
          if (unique.empty() || Distance(point, unique.back()) > tolerance)
          {
            unique.push_back(point);
          }
        }
        if (unique.size() > 1 && Distance(unique.front(), unique.back()) <= tolerance)
        {
          unique.pop_back();
        }
        if (unique.size() >= 3)
        {
          local_facets.push_back({conductor, std::move(unique)});
        }
      }
    }

    std::vector<double> local_records;
    for (const auto &facet : local_facets)
    {
      local_records.push_back(static_cast<double>(facet.conductor));
      local_records.push_back(static_cast<double>(facet.points.size()));
      for (const auto &point : facet.points)
      {
        local_records.insert(local_records.end(), point.begin(), point.end());
      }
    }
    if (statistics)
    {
      statistics->mask_facets_packed_local += local_facets.size();
      statistics->mask_payload_scalars_local += local_records.size();
    }
    MFEM_VERIFY(local_records.size() <=
                    static_cast<std::size_t>(std::numeric_limits<int>::max()),
                "Local plan-view facet data exceeds the MPI count limit!");
    const int local_value_count = static_cast<int>(local_records.size());
    std::vector<int> value_counts(Mpi::Size(mesh.GetComm()));
    Mpi::Allgather(1, &local_value_count, value_counts.data(), mesh.GetComm());
    std::vector<int> value_offsets(value_counts.size());
    int total_values = 0;
    for (std::size_t rank = 0; rank < value_counts.size(); rank++)
    {
      value_offsets[rank] = total_values;
      MFEM_VERIFY(value_counts[rank] <= std::numeric_limits<int>::max() - total_values,
                  "Global plan-view facet data exceeds the MPI count limit!");
      total_values += value_counts[rank];
    }
    if (statistics)
    {
      statistics->mask_gathered_scalars += total_values;
    }
    std::vector<double> records(total_values);
    Mpi::Allgatherv(local_value_count, local_records.data(), records.data(),
                    value_counts.data(), value_offsets.data(), mesh.GetComm());

    PlanViewGeometry result;
    result.lower = lower;
    result.upper = upper;
    result.process_axis = process_axis;
    for (std::size_t rank = 0; rank < value_counts.size(); rank++)
    {
      std::size_t offset = value_offsets[rank];
      const std::size_t end = offset + value_counts[rank];
      while (offset < end)
      {
        MFEM_VERIFY(end - offset >= 2, "Truncated gathered plan-view facet record!");
        PlanViewFacet facet;
        facet.conductor = static_cast<int>(std::llround(records[offset++]));
        const auto point_count = static_cast<std::size_t>(std::llround(records[offset++]));
        MFEM_VERIFY(facet.conductor > 0 && point_count >= 3 &&
                        point_count <= (end - offset) / 3,
                    "Invalid gathered plan-view facet!");
        facet.points.resize(point_count);
        for (auto &point : facet.points)
        {
          std::copy_n(records.data() + offset, 3, point.begin());
          offset += 3;
        }
        result.facets.push_back(std::move(facet));
      }
      MFEM_VERIFY(offset == end, "Invalid gathered plan-view facet data!");
    }
    return result;
  };
  auto HasOverlappingConductorStrips = [&](const std::vector<SpatialEdgeSite3D> &sites)
  {
    const double tolerance = 1.0e-10 * library.matching_radius;
    auto Polygon = [&](const SpatialEdgeSite3D &site)
    {
      const Point3D tangent = Normalize(Cross(site.gap_direction, site.process_normal));
      double begin = site.interval[0];
      double end = site.interval[1];
      if (begin <= -library.matching_radius + tolerance)
      {
        begin -= 2.0 * library.matching_radius;
      }
      if (end >= library.matching_radius - tolerance)
      {
        end += 2.0 * library.matching_radius;
      }
      const Point3D p0 = Add(site.point, Scale(begin, tangent));
      const Point3D p1 = Add(site.point, Scale(end, tangent));
      return std::array<Point3D, 4>{
          p0, p1, Add(p1, Scale(-3.0 * library.matching_radius, site.gap_direction)),
          Add(p0, Scale(-3.0 * library.matching_radius, site.gap_direction))};
    };
    auto Overlap = [&](const std::array<Point3D, 4> &first,
                       const std::array<Point3D, 4> &second, const Point3D &axis_x,
                       const Point3D &axis_y)
    {
      std::array<Point2D, 4> first_2d{}, second_2d{};
      for (int i = 0; i < 4; i++)
      {
        first_2d[i] = {Dot(first[i], axis_x), Dot(first[i], axis_y)};
        second_2d[i] = {Dot(second[i], axis_x), Dot(second[i], axis_y)};
      }
      for (const auto *polygon : {&first_2d, &second_2d})
      {
        for (int i = 0; i < 4; i++)
        {
          const auto &start = (*polygon)[i];
          const auto &end = (*polygon)[(i + 1) % 4];
          Point2D axis = {start[1] - end[1], end[0] - start[0]};
          const double length = Norm(axis);
          if (length <= tolerance)
          {
            continue;
          }
          axis[0] /= length;
          axis[1] /= length;
          std::array<double, 4> first_projection{}, second_projection{};
          for (int point = 0; point < 4; point++)
          {
            first_projection[point] = Dot(first_2d[point], axis);
            second_projection[point] = Dot(second_2d[point], axis);
          }
          const auto [first_min, first_max] =
              std::minmax_element(first_projection.begin(), first_projection.end());
          const auto [second_min, second_max] =
              std::minmax_element(second_projection.begin(), second_projection.end());
          if (*first_max <= *second_min + tolerance ||
              *second_max <= *first_min + tolerance)
          {
            return false;
          }
        }
      }
      return true;
    };

    for (std::size_t i = 0; i < sites.size(); i++)
    {
      const auto &first = sites[i];
      const Point3D tangent = Normalize(Cross(first.gap_direction, first.process_normal));
      const auto first_polygon = Polygon(first);
      for (std::size_t j = i + 1; j < sites.size(); j++)
      {
        const auto &second = sites[j];
        if (first.conductor == second.conductor ||
            Dot(first.process_normal, second.process_normal) <= 1.0 - 1.0e-10 ||
            std::abs(Dot(Subtract(second.point, first.point), first.process_normal)) >
                tolerance)
        {
          continue;
        }
        if (Overlap(first_polygon, Polygon(second), tangent, first.gap_direction))
        {
          return true;
        }
      }
    }
    return false;
  };
  auto FindMatchingSpatialModel = [&](const std::vector<SpatialEdgeSite3D> &sites)
      -> std::optional<SpatialClusterSelection3D>
  {
    std::set<std::size_t> excluded;
    for (std::size_t i = 0; i < library.models.size(); i++)
    {
      if (!library.models[i].plan_view_boundary)
      {
        excluded.insert(i);
      }
    }
    auto MatchesPlanView =
        [&](const SpatialClusterSelection3D &selection, const LibraryModel &model)
    {
      std::map<int, int> conductor_by_metal_component;
      for (std::size_t model_edge = 0; model_edge < model.spatial_edges.size();
           model_edge++)
      {
        const auto &site = sites[selection.model_to_site[model_edge]];
        if (site.metal_component < 0)
        {
          continue;
        }
        auto [component, inserted] = conductor_by_metal_component.emplace(
            site.metal_component, model.spatial_edges[model_edge].conductor);
        MFEM_VERIFY(inserted ||
                        component->second == model.spatial_edges[model_edge].conductor,
                    "A connected metal surface maps to multiple plan-view conductors!");
      }
      std::set<int> expected_conductors;
      for (const auto &edge : model.spatial_edges)
      {
        expected_conductors.insert(edge.conductor);
      }
      const auto plan_view = GatherPlanViewFacets(sites, selection.origin, selection.axes,
                                                  conductor_by_metal_component);
      std::set<int> found_conductors;
      for (const auto &facet : plan_view.facets)
      {
        found_conductors.insert(facet.conductor);
      }
      const auto clip_bounds = std::make_pair(plan_view.lower, plan_view.upper);
      const std::optional<decltype(clip_bounds)> classified_bounds =
          HasClassifiedPlanViewBoundary(*model.plan_view_boundary)
              ? std::optional<decltype(clip_bounds)>(clip_bounds)
              : std::nullopt;
      return found_conductors == expected_conductors &&
             CanonicalPlanViewBoundary(plan_view.facets, library.matching_radius,
                                       plan_view.process_axis,
                                       classified_bounds) == *model.plan_view_boundary;
    };
    if (auto selection =
            FindSpatialClusterLibraryModel(library, sites, excluded, MatchesPlanView))
    {
      return selection;
    }

    if (HasOverlappingConductorStrips(sites))
    {
      return std::nullopt;
    }
    excluded.clear();
    for (std::size_t i = 0; i < library.models.size(); i++)
    {
      if (library.models[i].plan_view_boundary)
      {
        excluded.insert(i);
      }
    }
    return FindSpatialClusterLibraryModel(library, sites, excluded);
  };
  auto DescribeMissingSpatialGeometry = [&](const std::vector<SpatialEdgeSite3D> &sites)
      -> std::pair<nlohmann::json, std::map<int, std::map<InterfaceDielectric, int>>>
  {
    MFEM_ASSERT(requirements && sites.size() >= 2,
                "Missing spatial geometry requires at least two edge sites!");
    std::optional<
        std::pair<nlohmann::json, std::map<int, std::map<InterfaceDielectric, int>>>>
        best;
    std::string best_key;
    for (const auto &anchor : sites)
    {
      const Point3D tangent = Normalize(Cross(anchor.gap_direction, anchor.process_normal));
      const std::array<Point3D, 3> axes = {anchor.gap_direction, anchor.process_normal,
                                           tangent};
      struct Description
      {
        const SpatialEdgeSite3D *site = nullptr;
        Point3D point{};
        Point3D gap{};
        Point3D normal{};
        std::string key;
      };
      std::vector<Description> descriptions;
      descriptions.reserve(sites.size());
      for (const auto &site : sites)
      {
        const Point3D relative = Subtract(site.point, anchor.point);
        Description description;
        description.site = &site;
        for (int d = 0; d < 3; d++)
        {
          description.point[d] = requirements->ScaleLength(Dot(relative, axes[d]));
          description.gap[d] =
              requirements->SnapDirection(Dot(site.gap_direction, axes[d]));
          description.normal[d] =
              requirements->SnapDirection(Dot(site.process_normal, axes[d]));
        }
        description.key = nlohmann::json{{"Point", description.point},
                                         {"GapDirection", description.gap},
                                         {"ProcessNormal", description.normal}}
                              .dump();
        descriptions.push_back(std::move(description));
      }
      std::sort(descriptions.begin(), descriptions.end(),
                [](const auto &first, const auto &second)
                {
                  if (first.key != second.key)
                  {
                    return first.key < second.key;
                  }
                  return first.site->conductor < second.site->conductor;
                });

      std::map<int, int> conductor_ids;
      std::map<int, int> conductor_by_metal_component;
      std::map<TargetSignature, int> slots;
      std::map<int, std::map<InterfaceDielectric, int>> targets_by_slot;
      nlohmann::json edges = nlohmann::json::array();
      for (const auto &description : descriptions)
      {
        const auto &site = *description.site;
        auto [conductor, conductor_inserted] =
            conductor_ids.emplace(site.conductor, conductor_ids.size() + 1);
        (void)conductor_inserted;
        if (site.metal_component >= 0)
        {
          auto [metal_component, inserted] =
              conductor_by_metal_component.emplace(site.metal_component, conductor->second);
          MFEM_VERIFY(inserted || metal_component->second == conductor->second,
                      "A connected metal surface maps to multiple spatial-coupon "
                      "conductors!");
        }
        const TargetSignature signature(site.targets.begin(), site.targets.end());
        auto [slot, slot_inserted] = slots.emplace(signature, slots.size());
        if (slot_inserted)
        {
          targets_by_slot.emplace(slot->second, site.targets);
        }
        edges.push_back({{"Point", description.point},
                         {"GapDirection", description.gap},
                         {"ProcessNormal", description.normal},
                         {"Interval",
                          {requirements->ScaleLength(site.interval[0]),
                           requirements->ScaleLength(site.interval[1])}},
                         {"Conductor", conductor->second},
                         {"InterfaceSlot", slot->second},
                         {"BoundaryCondition", requirements->DescribeBoundaryCondition(
                                                   site.boundary_condition)}});
      }
      nlohmann::json spatial_geometry = {{"EdgeCount", edges.size()},
                                         {"Edges", std::move(edges)}};
      if (surface.retain_faces && !conductor_by_metal_component.empty())
      {
        const auto &model_anchor = descriptions.front();
        Point3D model_point{};
        for (int d = 0; d < 3; d++)
        {
          model_point[d] = requirements->UnscaleLength(model_anchor.point[d]);
        }
        const auto [model_origin, model_axes] = AlignSpatialFrame(
            model_point, model_anchor.gap, model_anchor.normal, *model_anchor.site);
        const auto plan_view = GatherPlanViewFacets(sites, model_origin, model_axes,
                                                    conductor_by_metal_component);
        std::vector<nlohmann::json> facets;
        std::set<int> found_conductors;
        for (const auto &facet : plan_view.facets)
        {
          std::vector<Point3D> scaled(facet.points.size());
          for (std::size_t i = 0; i < facet.points.size(); i++)
          {
            for (int d = 0; d < 3; d++)
            {
              scaled[i][d] = requirements->ScaleLength(facet.points[i][d]);
            }
          }
          auto Sequence = [&](std::size_t start, bool reverse)
          {
            nlohmann::json points = nlohmann::json::array();
            for (std::size_t step = 0; step < scaled.size(); step++)
            {
              const std::size_t index = reverse
                                            ? (start + scaled.size() - step) % scaled.size()
                                            : (start + step) % scaled.size();
              points.push_back(scaled[index]);
            }
            return points;
          };
          nlohmann::json canonical;
          std::string canonical_key;
          for (std::size_t start = 0; start < scaled.size(); start++)
          {
            for (const bool reverse : {false, true})
            {
              auto candidate = Sequence(start, reverse);
              const std::string key = candidate.dump();
              if (canonical.is_null() || key < canonical_key)
              {
                canonical = std::move(candidate);
                canonical_key = key;
              }
            }
          }
          facets.push_back(
              {{"Conductor", facet.conductor}, {"Points", std::move(canonical)}});
          found_conductors.insert(facet.conductor);
        }
        std::set<int> expected_conductors;
        for (const auto &[component, conductor] : conductor_by_metal_component)
        {
          (void)component;
          expected_conductors.insert(conductor);
        }
        if (found_conductors == expected_conductors)
        {
          std::sort(facets.begin(), facets.end(), [](const auto &first, const auto &second)
                    { return first.dump() < second.dump(); });
          facets.erase(std::unique(facets.begin(), facets.end(),
                                   [](const auto &first, const auto &second)
                                   { return first == second; }),
                       facets.end());
          spatial_geometry["PlanViewFacets"] = std::move(facets);
          spatial_geometry["PlanViewBoundary"] =
              nlohmann::json::parse(CanonicalPlanViewBoundary(
                  plan_view.facets, library.matching_radius, plan_view.process_axis,
                  std::make_pair(plan_view.lower, plan_view.upper)));
        }
      }
      const std::string key = spatial_geometry.dump();
      if (!best || key < best_key)
      {
        best = std::make_pair(std::move(spatial_geometry), std::move(targets_by_slot));
        best_key = key;
      }
    }
    MFEM_ASSERT(best, "Unable to describe a missing spatial edge cluster!");
    return std::move(*best);
  };

  std::vector<EdgeSegment3D> global_segments;
  const bool has_cross_interface_spatial_model =
      std::any_of(library.models.begin(), library.models.end(),
                  [](const auto &model)
                  {
                    if (model.topology != LibraryTopology::SPATIAL_EDGE_CLUSTER)
                    {
                      return false;
                    }
                    std::set<int> slots;
                    for (const auto &edge : model.spatial_edges)
                    {
                      slots.insert(edge.interface_slot);
                    }
                    return slots.size() > 1;
                  });
  if (has_cross_interface_spatial_model || requirements)
  {
    for (const auto &[signature, group] : groups_by_targets)
    {
      (void)signature;
      auto segments = BuildSegments(group, group.segments);
      global_segments.insert(global_segments.end(),
                             std::make_move_iterator(segments.begin()),
                             std::make_move_iterator(segments.end()));
    }
  }

  const double global_interaction_distance =
      2.0 * library.matching_radius * (1.0 - 16.0 * std::numeric_limits<double>::epsilon());
  const double global_interaction_distance_squared =
      global_interaction_distance * global_interaction_distance;
  struct GlobalSpatialInteractionEvent
  {
    std::size_t first = 0;
    std::size_t second = 0;
    double first_distance = 0.0;
    double second_distance = 0.0;
    double distance_squared = 0.0;
    Point3D center{};
  };
  std::map<std::pair<int, int>, GlobalSpatialInteractionEvent> global_chain_events;
  std::vector<std::pair<std::size_t, std::size_t>> global_candidate_pairs;
  if (!global_segments.empty())
  {
    const SegmentBoxIndex index(SegmentGeometry(global_segments));
    global_candidate_pairs = index.CandidatePairs(global_interaction_distance);
  }
  if (statistics)
  {
    statistics->pair_checks_global_spatial += global_candidate_pairs.size();
  }
  for (const auto &[i, j] : global_candidate_pairs)
  {
    const auto &first_source = geometry.segments[global_segments[i].geometry_index];
    const auto &second_source = geometry.segments[global_segments[j].geometry_index];
    if (global_segments[i].targets == global_segments[j].targets ||
        first_source.physical_chain == second_source.physical_chain ||
        SegmentsShareVertex(first_source, second_source))
    {
      continue;
    }
    const auto closest =
        ClosestSegmentApproach(global_segments[i].p0, global_segments[i].p1,
                               global_segments[j].p0, global_segments[j].p1);
    if (closest.distance_squared >= global_interaction_distance_squared)
    {
      continue;
    }
    int first_chain = first_source.physical_chain;
    int second_chain = second_source.physical_chain;
    if (first_chain > second_chain)
    {
      std::swap(first_chain, second_chain);
    }
    const auto key = std::make_pair(first_chain, second_chain);
    auto event = global_chain_events.find(key);
    if (event == global_chain_events.end() ||
        closest.distance_squared < event->second.distance_squared)
    {
      const double first_distance = closest.first * global_segments[i].length;
      const double second_distance = closest.second * global_segments[j].length;
      global_chain_events[key] = {
          i,
          j,
          first_distance,
          second_distance,
          closest.distance_squared,
          Scale(0.5, Add(Interpolate(global_segments[i], first_distance),
                         Interpolate(global_segments[j], second_distance)))};
    }
  }

  std::vector<GlobalSpatialInteractionEvent> global_spatial_events;
  global_spatial_events.reserve(global_chain_events.size());
  for (auto &[chains, event] : global_chain_events)
  {
    (void)chains;
    global_spatial_events.push_back(std::move(event));
  }
  if (statistics)
  {
    statistics->spatial_events += global_spatial_events.size();
  }
  std::vector<SpatialClusterSelection3D> cross_interface_selections;
  std::set<std::pair<int, int>> described_cross_interface_pairs;
  std::vector<bool> visited_global_event(global_spatial_events.size(), false);
  for (std::size_t seed = 0; seed < global_spatial_events.size(); seed++)
  {
    if (visited_global_event[seed])
    {
      continue;
    }
    std::vector<std::size_t> event_component = {seed};
    visited_global_event[seed] = true;
    for (std::size_t cursor = 0; cursor < event_component.size(); cursor++)
    {
      const std::size_t event_index = event_component[cursor];
      for (std::size_t candidate = 0; candidate < global_spatial_events.size(); candidate++)
      {
        if (!visited_global_event[candidate] &&
            Distance(global_spatial_events[event_index].center,
                     global_spatial_events[candidate].center) < global_interaction_distance)
        {
          visited_global_event[candidate] = true;
          event_component.push_back(candidate);
        }
      }
    }

    std::map<int, std::vector<Point3D>> points_by_chain;
    for (const std::size_t event_index : event_component)
    {
      const auto &event = global_spatial_events[event_index];
      const int first_chain =
          geometry.segments[global_segments[event.first].geometry_index].physical_chain;
      const int second_chain =
          geometry.segments[global_segments[event.second].geometry_index].physical_chain;
      points_by_chain[first_chain].push_back(
          Interpolate(global_segments[event.first], event.first_distance));
      points_by_chain[second_chain].push_back(
          Interpolate(global_segments[event.second], event.second_distance));
    }

    std::vector<SpatialEdgeSite3D> sites;
    sites.reserve(points_by_chain.size());
    for (const auto &[chain, points] : points_by_chain)
    {
      Point3D average{};
      for (const auto &point : points)
      {
        average = Add(average, point);
      }
      average = Scale(1.0 / points.size(), average);

      std::optional<std::size_t> best_segment;
      double best_segment_distance = 0.0;
      double best_distance_squared = mfem::infinity();
      for (std::size_t segment_index = 0; segment_index < global_segments.size();
           segment_index++)
      {
        const auto &source =
            geometry.segments[global_segments[segment_index].geometry_index];
        if (source.physical_chain != chain)
        {
          continue;
        }
        const double distance =
            std::clamp(Dot(Subtract(average, global_segments[segment_index].p0),
                           global_segments[segment_index].tangent),
                       0.0, global_segments[segment_index].length);
        const Point3D point = Interpolate(global_segments[segment_index], distance);
        const double distance_squared =
            Dot(Subtract(average, point), Subtract(average, point));
        if (distance_squared < best_distance_squared)
        {
          best_segment = segment_index;
          best_segment_distance = distance;
          best_distance_squared = distance_squared;
        }
      }
      MFEM_ASSERT(best_segment, "Missing segment for a cross-interface spatial chain!");
      const auto &segment = global_segments[*best_segment];
      const auto interval = SpatialInterval(global_segments, *best_segment,
                                            best_segment_distance, library.matching_radius);
      sites.push_back({chain, segment.geometry_index, *best_segment, best_segment_distance,
                       interval, Interpolate(segment, best_segment_distance),
                       segment.axis_u, segment.axis_v, segment.conductor,
                       segment.metal_component, segment.targets,
                       segment.boundary_condition});
    }

    auto selection = FindMatchingSpatialModel(sites);
    if (!selection || selection->targets_by_slot.size() < 2)
    {
      if (requirements)
      {
        auto [spatial_geometry, targets_by_slot] = DescribeMissingSpatialGeometry(sites);
        requirements->Add(3, LibraryTopology::SPATIAL_EDGE_CLUSTER, targets_by_slot,
                          sites.front().boundary_condition, spatial_geometry, library,
                          nullptr, 2.0 * library.matching_radius * sites.size(),
                          "No compatible cross-interface spatial-edge cluster model");
      }
      for (const std::size_t event_index : event_component)
      {
        const auto &event = global_spatial_events[event_index];
        int first_chain =
            geometry.segments[global_segments[event.first].geometry_index].physical_chain;
        int second_chain =
            geometry.segments[global_segments[event.second].geometry_index].physical_chain;
        if (first_chain > second_chain)
        {
          std::swap(first_chain, second_chain);
        }
        described_cross_interface_pairs.emplace(first_chain, second_chain);
      }
      continue;
    }
    for (const std::size_t event_index : event_component)
    {
      const auto &event = global_spatial_events[event_index];
      int first_chain =
          geometry.segments[global_segments[event.first].geometry_index].physical_chain;
      int second_chain =
          geometry.segments[global_segments[event.second].geometry_index].physical_chain;
      if (first_chain > second_chain)
      {
        std::swap(first_chain, second_chain);
      }
      described_cross_interface_pairs.emplace(first_chain, second_chain);
      selection->interactions.push_back({first_chain, second_chain, event.center});
    }
    cross_interface_selections.push_back(std::move(*selection));
  }

  for (const auto &selection : cross_interface_selections)
  {
    const auto &weighted_model = selection.response.models.front();
    const auto &source = library.models[weighted_model.index];
    if (requirements)
    {
      requirements->Add(3, LibraryTopology::SPATIAL_EDGE_CLUSTER, selection.targets_by_slot,
                        selection.sites.front().boundary_condition,
                        SpatialGeometry(selection), library, &selection.response,
                        SpatialLength(selection));
    }
    auto model = source.response;
    model.idx = next_model_index++;
    model.name = source.name;
    model.topology = TopologyName(source.topology);
    MapLibraryInterfaces(source, selection.targets_by_slot, model);
    result.models.push_back(std::move(model));

    ResponsePatchData patch;
    patch.model = result.models.back().idx;
    patch.origin = selection.origin;
    patch.axis_u = selection.axes[0];
    patch.axis_v = selection.axes[1];
    patch.axis_w = selection.axes[2];
    patch.conductor_references = selection.response.conductor_references;
    patch.weight = 1.0;
    patch.maxwell_reference_is_pec = std::all_of(
        selection.sites.begin(), selection.sites.end(), [](const auto &site)
        { return site.boundary_condition.type == MetalBoundaryConditionType::PEC; });
    for (const auto &reference : patch.conductor_references)
    {
      patch.maxwell_conductor_anchors.push_back(
          TransformLocalPoint(patch.origin, selection.axes, reference));
    }
    result.patches.push_back(std::move(patch));
  }
  matched_spatial_cluster_patches += static_cast<int>(cross_interface_selections.size());
  if (diagnostics)
  {
    for (const auto &selection : cross_interface_selections)
    {
      diagnostics->boundary_law_verified &=
          IsBoundaryLawVerified(library.models[selection.response.models.front().index]);
      diagnostics->maximum_library_distance = std::max(
          diagnostics->maximum_library_distance, selection.response.normalized_distance);
    }
  }

  auto IsCrossInterfaceSpatiallyMatched =
      [&](const MetalEdgeSegment &first, const Point3D &first_p0, const Point3D &first_p1,
          const MetalEdgeSegment &second, const Point3D &second_p0,
          const Point3D &second_p1)
  {
    int first_chain = first.physical_chain;
    int second_chain = second.physical_chain;
    if (first_chain > second_chain)
    {
      std::swap(first_chain, second_chain);
    }
    const auto closest = ClosestSegmentApproach(first_p0, first_p1, second_p0, second_p1);
    const Point3D center = Scale(
        0.5, Add(Add(first_p0, Scale(closest.first, Subtract(first_p1, first_p0))),
                 Add(second_p0, Scale(closest.second, Subtract(second_p1, second_p0)))));
    return std::any_of(cross_interface_selections.begin(), cross_interface_selections.end(),
                       [&](const auto &selection)
                       {
                         return std::any_of(
                             selection.interactions.begin(), selection.interactions.end(),
                             [&](const auto &interaction)
                             {
                               return interaction.first_chain == first_chain &&
                                      interaction.second_chain == second_chain &&
                                      Distance(interaction.center, center) <
                                          global_interaction_distance;
                             });
                       });
  };

  std::vector<std::pair<Point3D, Point3D>> geometry_segment_points;
  geometry_segment_points.reserve(geometry.segments.size());
  for (const auto &segment : geometry.segments)
  {
    geometry_segment_points.emplace_back(geometry.vertices[segment.vertices[0]].coordinate,
                                         geometry.vertices[segment.vertices[1]].coordinate);
  }
  const SegmentBoxIndex geometry_segment_index(geometry_segment_points);

  for (const auto &segment_group : groups_by_targets)
  {
    const auto &segment_key = segment_group.first;
    const auto &group = segment_group.second;
    (void)segment_key;
    int group_matched_intervals = 0;
    int group_interpolated_paired_intervals = 0;
    double group_selected_length = 0.0;
    double group_corner_neighborhood_length = 0.0;
    double group_modeled_corner_neighborhood_length = 0.0;
    double group_maximum_curvature_ratio = 0.0;
    double group_maximum_library_distance = 0.0;
    const double interaction_distance =
        2.0 * group.matching_radius * (1.0 - 16.0 * std::numeric_limits<double>::epsilon());
    const double interaction_distance_squared = interaction_distance * interaction_distance;
    const std::set<std::size_t> group_segment_indices(group.segments.begin(),
                                                      group.segments.end());
    std::vector<std::size_t> usable_segment_indices;
    usable_segment_indices.reserve(group.segments.size());
    double total_selected_length = 0.0;
    int externally_conflicted_segments = 0;
    for (const std::size_t geometry_index : group.segments)
    {
      const auto &source = geometry.segments[geometry_index];
      const auto &p0 = geometry.vertices[source.vertices[0]].coordinate;
      const auto &p1 = geometry.vertices[source.vertices[1]].coordinate;
      total_selected_length += Distance(p0, p1);
      std::map<int, std::pair<std::size_t, double>> conflicts;
      const auto nearby_geometry =
          geometry_segment_index.Query(p0, p1, interaction_distance);
      if (statistics)
      {
        statistics->pair_checks_external_conflict += nearby_geometry.size();
      }
      for (const std::size_t other_index : nearby_geometry)
      {
        const auto &other = geometry.segments[other_index];
        if (other.type != MetalEdgeSegmentType::PHYSICAL ||
            group_segment_indices.find(other_index) != group_segment_indices.end() ||
            other.physical_chain == source.physical_chain ||
            SegmentsShareVertex(source, other))
        {
          continue;
        }
        const auto &q0 = geometry.vertices[other.vertices[0]].coordinate;
        const auto &q1 = geometry.vertices[other.vertices[1]].coordinate;
        const double distance_squared = SegmentDistanceSquared(p0, p1, q0, q1);
        if (distance_squared < interaction_distance_squared &&
            !IsCrossInterfaceSpatiallyMatched(source, p0, p1, other, q0, q1))
        {
          auto conflict = conflicts.find(other.physical_chain);
          if (conflict == conflicts.end() || distance_squared < conflict->second.second)
          {
            conflicts[other.physical_chain] = {other_index, distance_squared};
          }
        }
      }
      if (!conflicts.empty())
      {
        externally_conflicted_segments++;
        if (requirements)
        {
          std::map<int, std::pair<std::size_t, double>> undescribed_conflicts;
          for (const auto &[physical_chain, conflict] : conflicts)
          {
            const auto ordered = std::minmax(source.physical_chain, physical_chain);
            const std::pair<int, int> chains = {ordered.first, ordered.second};
            if (described_cross_interface_pairs.find(chains) ==
                described_cross_interface_pairs.end())
            {
              undescribed_conflicts.emplace(physical_chain, conflict);
            }
          }
          if (undescribed_conflicts.empty())
          {
            continue;
          }
          std::map<int, std::map<InterfaceDielectric, int>> targets_by_slot = {
              {0, group.targets}};
          nlohmann::json interactions = nlohmann::json::array();
          int slot = 1;
          for (const auto &[physical_chain, conflict] : undescribed_conflicts)
          {
            (void)physical_chain;
            const auto other_index = conflict.first;
            const auto &other = geometry.segments[other_index];
            const auto &q0 = geometry.vertices[other.vertices[0]].coordinate;
            const auto &q1 = geometry.vertices[other.vertices[1]].coordinate;
            const auto closest = ClosestSegmentApproach(p0, p1, q0, q1);
            const auto source_tangent = Normalize(Subtract(p1, p0));
            const auto other_tangent = Normalize(Subtract(q1, q0));
            interactions.push_back(
                {{"Separation",
                  requirements->ScaleLength(std::sqrt(closest.distance_squared))},
                 {"AngleDegrees",
                  requirements->SnapAngleDegrees(
                      std::acos(std::clamp(std::abs(Dot(source_tangent, other_tangent)),
                                           0.0, 1.0)) *
                      180.0 / std::acos(-1.0))},
                 {"BoundaryCondition",
                  other.conditions.empty()
                      ? "Unknown"
                      : BoundaryConditionName(other.conditions.front().type)}});
            if (auto targets = targets_by_segment.find(other_index);
                targets != targets_by_segment.end())
            {
              targets_by_slot.emplace(slot, targets->second);
            }
            slot++;
          }
          const auto boundary_condition =
              maxwell && !source.conditions.empty()
                  ? GetBoundaryConditionLaw(iodata.boundaries, source.conditions.front())
                  : MetalBoundaryLaw{};
          requirements->Add(3, LibraryTopology::SPATIAL_EDGE_CLUSTER, targets_by_slot,
                            boundary_condition,
                            {{"EdgeCount", undescribed_conflicts.size() + 1},
                             {"Interactions", std::move(interactions)}},
                            library, nullptr, Distance(p0, p1),
                            "Nearby physical edges use a different interface mapping");
        }
      }
      else
      {
        usable_segment_indices.push_back(geometry_index);
      }
    }
    if (diagnostics)
    {
      diagnostics->selected_length += total_selected_length;
      for (const auto &[type, target] : group.targets)
      {
        (void)type;
        diagnostics->selected_length_by_interface[target] += total_selected_length;
      }
    }
    if (externally_conflicted_segments > 0)
    {
      if (request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
      {
        MFEM_ABORT("A three-dimensional target edge is within 2R of a physical metal "
                   "edge with a different interface mapping!");
      }
      Mpi::Warning(
          "Omitting {} of {} three-dimensional target edge segments which are within 2R "
          "of a physical metal edge with a different interface mapping.\n",
          externally_conflicted_segments, static_cast<int>(group.segments.size()));
    }
    if (usable_segment_indices.empty())
    {
      unmatched_groups++;
      continue;
    }
    auto segments = BuildSegments(group, usable_segment_indices);
    SegmentBoxIndex segment_index(SegmentGeometry(segments));
    auto nearby_segment_pairs = segment_index.CandidatePairs(interaction_distance);
    std::map<std::size_t, std::size_t> local_indices;
    for (std::size_t i = 0; i < segments.size(); i++)
    {
      local_indices.emplace(segments[i].geometry_index, i);
    }

    struct ConnectedVertexNeighborhood
    {
      Point3D point;
      std::set<int> physical_chains;
    };
    auto BuildConnectedVertexNeighborhoods = [&]()
    {
      std::vector<ConnectedVertexNeighborhood> neighborhoods;
      for (std::size_t vertex = 0; vertex < geometry.vertices.size(); vertex++)
      {
        if (geometry.vertices[vertex].physical_type == MetalEdgeVertexType::REGULAR)
        {
          continue;
        }
        ConnectedVertexNeighborhood neighborhood{geometry.vertices[vertex].coordinate, {}};
        for (const std::size_t geometry_segment : geometry.vertices[vertex].segments)
        {
          auto local = local_indices.find(geometry_segment);
          if (local != local_indices.end())
          {
            neighborhood.physical_chains.insert(
                geometry.segments[geometry_segment].physical_chain);
          }
        }
        if (neighborhood.physical_chains.size() > 1)
        {
          neighborhoods.push_back(std::move(neighborhood));
        }
      }
      return neighborhoods;
    };

    // Fail closed only where the local geometry cannot be represented by an available
    // one- or two-edge coupon. Screening before corner construction prevents a rejected
    // interaction from leaving behind corner patches on the same mesh segments.
    const auto candidate_connected_vertex_neighborhoods =
        BuildConnectedVertexNeighborhoods();
    auto ConnectedNearVertex = [&](std::size_t first, std::size_t second)
    {
      const auto &first_source = geometry.segments[segments[first].geometry_index];
      const auto &second_source = geometry.segments[segments[second].geometry_index];
      return std::any_of(
          candidate_connected_vertex_neighborhoods.begin(),
          candidate_connected_vertex_neighborhoods.end(),
          [&](const auto &neighborhood)
          {
            return neighborhood.physical_chains.find(first_source.physical_chain) !=
                       neighborhood.physical_chains.end() &&
                   neighborhood.physical_chains.find(second_source.physical_chain) !=
                       neighborhood.physical_chains.end() &&
                   PointSegmentDistanceSquared(neighborhood.point, segments[first]) <
                       interaction_distance_squared &&
                   PointSegmentDistanceSquared(neighborhood.point, segments[second]) <
                       interaction_distance_squared;
          });
    };

    struct SpatialInteractionEvent
    {
      std::size_t first = 0;
      std::size_t second = 0;
      double first_distance = 0.0;
      double second_distance = 0.0;
      double distance_squared = 0.0;
      Point3D center{};
    };
    std::map<std::pair<int, int>, SpatialInteractionEvent> closest_chain_events;
    if (statistics)
    {
      statistics->pair_checks_group_spatial += nearby_segment_pairs.size();
    }
    for (const auto &[i, j] : nearby_segment_pairs)
    {
      const auto &first_source = geometry.segments[segments[i].geometry_index];
      const auto &second_source = geometry.segments[segments[j].geometry_index];
      if (first_source.physical_chain == second_source.physical_chain ||
          SegmentsShareVertex(first_source, second_source) || ConnectedNearVertex(i, j))
      {
        continue;
      }
      const auto closest = ClosestSegmentApproach(segments[i].p0, segments[i].p1,
                                                  segments[j].p0, segments[j].p1);
      if (closest.distance_squared >= interaction_distance_squared)
      {
        continue;
      }
      const double first_distance = closest.first * segments[i].length;
      const double second_distance = closest.second * segments[j].length;
      int first_chain = first_source.physical_chain;
      int second_chain = second_source.physical_chain;
      if (first_chain > second_chain)
      {
        std::swap(first_chain, second_chain);
      }
      const auto key = std::make_pair(first_chain, second_chain);
      auto event = closest_chain_events.find(key);
      if (event == closest_chain_events.end() ||
          closest.distance_squared < event->second.distance_squared)
      {
        closest_chain_events[key] = {
            i,
            j,
            first_distance,
            second_distance,
            closest.distance_squared,
            Scale(0.5, Add(Interpolate(segments[i], first_distance),
                           Interpolate(segments[j], second_distance)))};
      }
    }
    std::vector<SpatialInteractionEvent> spatial_events;
    spatial_events.reserve(closest_chain_events.size());
    for (auto &[chains, event] : closest_chain_events)
    {
      (void)chains;
      spatial_events.push_back(std::move(event));
    }

    if (statistics)
    {
      statistics->spatial_events += spatial_events.size();
    }
    std::vector<SpatialClusterSelection3D> spatial_cluster_selections;
    std::set<std::pair<int, int>> described_spatial_pairs;
    std::vector<bool> visited_spatial_event(spatial_events.size(), false);
    for (std::size_t seed = 0; seed < spatial_events.size(); seed++)
    {
      if (visited_spatial_event[seed])
      {
        continue;
      }
      std::vector<std::size_t> event_component = {seed};
      visited_spatial_event[seed] = true;
      for (std::size_t cursor = 0; cursor < event_component.size(); cursor++)
      {
        const std::size_t event_index = event_component[cursor];
        for (std::size_t candidate = 0; candidate < spatial_events.size(); candidate++)
        {
          if (!visited_spatial_event[candidate] &&
              Distance(spatial_events[event_index].center,
                       spatial_events[candidate].center) < interaction_distance)
          {
            visited_spatial_event[candidate] = true;
            event_component.push_back(candidate);
          }
        }
      }

      std::map<int, std::vector<Point3D>> points_by_chain;
      for (const std::size_t event_index : event_component)
      {
        const auto &event = spatial_events[event_index];
        const int first_chain =
            geometry.segments[segments[event.first].geometry_index].physical_chain;
        const int second_chain =
            geometry.segments[segments[event.second].geometry_index].physical_chain;
        points_by_chain[first_chain].push_back(
            Interpolate(segments[event.first], event.first_distance));
        points_by_chain[second_chain].push_back(
            Interpolate(segments[event.second], event.second_distance));
      }
      if (points_by_chain.size() < 2)
      {
        continue;
      }
      std::vector<SpatialEdgeSite3D> sites;
      sites.reserve(points_by_chain.size());
      for (const auto &[chain, points] : points_by_chain)
      {
        Point3D average{};
        for (const auto &point : points)
        {
          average = Add(average, point);
        }
        average = Scale(1.0 / points.size(), average);

        std::optional<std::size_t> best_segment;
        double best_segment_distance = 0.0;
        double best_distance_squared = mfem::infinity();
        for (std::size_t segment_index = 0; segment_index < segments.size();
             segment_index++)
        {
          const auto &source = geometry.segments[segments[segment_index].geometry_index];
          if (source.physical_chain != chain)
          {
            continue;
          }
          const double distance =
              std::clamp(Dot(Subtract(average, segments[segment_index].p0),
                             segments[segment_index].tangent),
                         0.0, segments[segment_index].length);
          const double distance_squared =
              Dot(Subtract(average, Interpolate(segments[segment_index], distance)),
                  Subtract(average, Interpolate(segments[segment_index], distance)));
          if (distance_squared < best_distance_squared)
          {
            best_segment = segment_index;
            best_segment_distance = distance;
            best_distance_squared = distance_squared;
          }
        }
        MFEM_ASSERT(best_segment, "Missing segment for a spatial edge chain!");
        const auto &segment = segments[*best_segment];
        const auto interval = SpatialInterval(segments, *best_segment,
                                              best_segment_distance, group.matching_radius);
        sites.push_back({chain, segment.geometry_index, *best_segment,
                         best_segment_distance, interval,
                         Interpolate(segment, best_segment_distance), segment.axis_u,
                         segment.axis_v, segment.conductor, segment.metal_component,
                         segment.targets, segment.boundary_condition});
      }

      auto selection = FindMatchingSpatialModel(sites);
      if (!selection)
      {
        if (requirements)
        {
          auto [spatial_geometry, targets_by_slot] = DescribeMissingSpatialGeometry(sites);
          requirements->Add(3, LibraryTopology::SPATIAL_EDGE_CLUSTER, targets_by_slot,
                            sites.front().boundary_condition, spatial_geometry, library,
                            nullptr, 2.0 * group.matching_radius * sites.size(),
                            "No compatible spatial-edge cluster model");
        }
        for (const std::size_t event_index : event_component)
        {
          const auto &event = spatial_events[event_index];
          int first_chain =
              geometry.segments[segments[event.first].geometry_index].physical_chain;
          int second_chain =
              geometry.segments[segments[event.second].geometry_index].physical_chain;
          if (first_chain > second_chain)
          {
            std::swap(first_chain, second_chain);
          }
          described_spatial_pairs.emplace(first_chain, second_chain);
        }
        continue;
      }
      group_maximum_library_distance =
          std::max(group_maximum_library_distance, selection->response.normalized_distance);
      for (const std::size_t event_index : event_component)
      {
        const auto &event = spatial_events[event_index];
        int first_chain =
            geometry.segments[segments[event.first].geometry_index].physical_chain;
        int second_chain =
            geometry.segments[segments[event.second].geometry_index].physical_chain;
        if (first_chain > second_chain)
        {
          std::swap(first_chain, second_chain);
        }
        described_spatial_pairs.emplace(first_chain, second_chain);
        selection->interactions.push_back({first_chain, second_chain, event.center});
      }
      spatial_cluster_selections.push_back(std::move(*selection));
    }
    auto IsSpatiallyMatched =
        [&](std::size_t first, std::size_t second,
            const std::vector<SpatialClusterSelection3D::InteractionNeighborhood>
                &interactions)
    {
      int first_chain = geometry.segments[segments[first].geometry_index].physical_chain;
      int second_chain = geometry.segments[segments[second].geometry_index].physical_chain;
      if (first_chain > second_chain)
      {
        std::swap(first_chain, second_chain);
      }
      const auto closest = ClosestSegmentApproach(segments[first].p0, segments[first].p1,
                                                  segments[second].p0, segments[second].p1);
      const Point3D center = Scale(
          0.5,
          Add(Interpolate(segments[first], closest.first * segments[first].length),
              Interpolate(segments[second], closest.second * segments[second].length)));
      return std::any_of(interactions.begin(), interactions.end(),
                         [&](const auto &interaction)
                         {
                           return interaction.first_chain == first_chain &&
                                  interaction.second_chain == second_chain &&
                                  Distance(interaction.center, center) <
                                      interaction_distance;
                         });
    };
    std::vector<SpatialClusterSelection3D::InteractionNeighborhood>
        candidate_spatial_interactions;
    for (const auto &selection : spatial_cluster_selections)
    {
      candidate_spatial_interactions.insert(candidate_spatial_interactions.end(),
                                            selection.interactions.begin(),
                                            selection.interactions.end());
    }

    std::vector<EdgePair3D> candidate_pairs;
    std::vector<bool> candidate_pair_has_model;
    std::vector<LibraryTopology> candidate_pair_topologies;
    std::set<std::size_t> unsupported_segments;
    int nonparallel_interactions = 0;
    int incompatible_process_interactions = 0;
    int process_offset_interactions = 0;
    int unclassified_interactions = 0;
    int unclassified_same_conductor_interactions = 0;
    int unclassified_different_conductor_interactions = 0;
    int missing_library_interactions = 0;
    int multiedge_interactions = 0;
    auto RejectInteraction =
        [&](std::size_t first, std::size_t second, int &counter, const char *message,
            std::optional<LibraryTopology> requirement_topology = std::nullopt)
    {
      counter++;
      if (requirements)
      {
        const auto topology =
            requirement_topology.value_or(LibraryTopology::SPATIAL_EDGE_CLUSTER);
        int first_chain = geometry.segments[segments[first].geometry_index].physical_chain;
        int second_chain =
            geometry.segments[segments[second].geometry_index].physical_chain;
        if (first_chain > second_chain)
        {
          std::swap(first_chain, second_chain);
        }
        const bool already_described =
            topology == LibraryTopology::SPATIAL_EDGE_CLUSTER &&
            described_spatial_pairs.find({first_chain, second_chain}) !=
                described_spatial_pairs.end();
        if (!already_described)
        {
          const auto closest =
              ClosestSegmentApproach(segments[first].p0, segments[first].p1,
                                     segments[second].p0, segments[second].p1);
          const double angle = std::acos(std::clamp(
              std::abs(Dot(segments[first].tangent, segments[second].tangent)), 0.0, 1.0));
          nlohmann::json geometry = {{"EdgeCount", 2},
                                     {"Separation", requirements->ScaleLength(std::sqrt(
                                                        closest.distance_squared))}};
          if (topology == LibraryTopology::SPATIAL_EDGE_CLUSTER)
          {
            geometry["AngleDegrees"] =
                requirements->SnapAngleDegrees(angle * 180.0 / std::acos(-1.0));
          }
          requirements->Add(3, topology, group.targets, segments[first].boundary_condition,
                            geometry, library, nullptr,
                            std::min(segments[first].length, segments[second].length),
                            message);
        }
      }
      if (request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
      {
        MFEM_ABORT(message);
      }
      unsupported_segments.insert(first);
      unsupported_segments.insert(second);
    };
    if (statistics)
    {
      statistics->pair_checks_safety += nearby_segment_pairs.size();
    }
    for (const auto &[i, j] : nearby_segment_pairs)
    {
      const auto &first_source = geometry.segments[segments[i].geometry_index];
      const auto &second_source = geometry.segments[segments[j].geometry_index];
      if (first_source.physical_chain == second_source.physical_chain ||
          SegmentsShareVertex(first_source, second_source) ||
          SegmentDistanceSquared(segments[i].p0, segments[i].p1, segments[j].p0,
                                 segments[j].p1) >= interaction_distance_squared)
      {
        continue;
      }
      if (ConnectedNearVertex(i, j) ||
          IsSpatiallyMatched(i, j, candidate_spatial_interactions))
      {
        continue;
      }

      const double tangent_dot = Dot(segments[i].tangent, segments[j].tangent);
      if (std::abs(tangent_dot) < 1.0 - 1.0e-8)
      {
        RejectInteraction(i, j, nonparallel_interactions,
                          "Nearby three-dimensional metal edges are not parallel!");
        continue;
      }
      const double second_s0 =
          Dot(Subtract(segments[j].p0, segments[i].p0), segments[i].tangent);
      const double second_s1 =
          Dot(Subtract(segments[j].p1, segments[i].p0), segments[i].tangent);
      const double first_begin = std::max(0.0, std::min(second_s0, second_s1));
      const double first_end = std::min(segments[i].length, std::max(second_s0, second_s1));
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
      EdgePair3D pair{i,
                      j,
                      first_begin,
                      first_end,
                      std::max(0.0, second_begin),
                      std::min(segments[j].length, second_end)};

      const Point3D second_point = Interpolate(segments[j], second_mid);
      const Point3D direction = Normalize(Subtract(second_point, first_mid));
      if (Dot(segments[i].axis_v, segments[j].axis_v) <= 0.95)
      {
        RejectInteraction(
            i, j, incompatible_process_interactions,
            "Nearby three-dimensional edges have incompatible process normals!");
        continue;
      }
      const Point3D process_normal = Normalize(Add(segments[i].axis_v, segments[j].axis_v));
      if (std::abs(Dot(direction, process_normal)) > 1.0e-8)
      {
        RejectInteraction(
            i, j, process_offset_interactions,
            "Nearby three-dimensional edges are offset along the process normal!");
        continue;
      }
      const bool facing = Dot(segments[i].axis_u, direction) > 0.95 &&
                          Dot(segments[j].axis_u, direction) < -0.95;
      const bool outward = Dot(segments[i].axis_u, direction) < -0.95 &&
                           Dot(segments[j].axis_u, direction) > 0.95;
      const bool same_conductor = segments[i].conductor == segments[j].conductor;
      std::optional<LibraryTopology> topology;
      if (facing)
      {
        topology = same_conductor ? LibraryTopology::SAME_CONDUCTOR_GAP
                                  : LibraryTopology::DIFFERENT_CONDUCTOR_GAP;
      }
      else if (outward)
      {
        // The two gap directions point away from the interval between the edges, so
        // that interval is occupied by one physical metal strip. Perimeter loops on
        // opposite sides of the strip need not be connected in the edge graph.
        topology = LibraryTopology::SAME_CONDUCTOR_STRIP;
      }
      if (!topology)
      {
        if (same_conductor)
        {
          unclassified_same_conductor_interactions++;
        }
        else
        {
          unclassified_different_conductor_interactions++;
        }
        RejectInteraction(i, j, unclassified_interactions,
                          "No canonical paired-edge topology for nearby "
                          "three-dimensional metal edges!");
        continue;
      }
      if (!SameBoundaryLaw(segments[i].boundary_condition, segments[j].boundary_condition))
      {
        RejectInteraction(i, j, unclassified_interactions,
                          "Nearby three-dimensional edges use distinct metal boundary "
                          "conditions!");
        continue;
      }
      candidate_pairs.push_back(pair);
      candidate_pair_topologies.push_back(*topology);
      candidate_pair_has_model.push_back(FindLibraryModel(library, *topology,
                                                          Distance(first_mid, second_point),
                                                          segments[i].boundary_condition)
                                             .has_value());
    }

    const auto candidate_parallel_clusters =
        FindParallelClusterSpans(library, segments, candidate_pairs);
    const auto &candidate_parallel_cluster_spans = candidate_parallel_clusters.matched;
    auto PairCoveredByParallelCluster = [&](const EdgePair3D &pair)
    {
      std::vector<std::pair<double, double>> covered;
      const auto &first = segments[pair.first];
      for (const auto &span : candidate_parallel_cluster_spans)
      {
        const auto &cluster = span.selection.ordered_edges;
        if (std::find(cluster.begin(), cluster.end(), pair.first) == cluster.end() ||
            std::find(cluster.begin(), cluster.end(), pair.second) == cluster.end())
        {
          continue;
        }
        const double orientation = Dot(first.tangent, span.tangent);
        if (std::abs(orientation) <= 1.0 - 1.0e-8)
        {
          continue;
        }
        double begin = (span.begin - Dot(first.p0, span.tangent)) / orientation;
        double end = (span.end - Dot(first.p0, span.tangent)) / orientation;
        if (begin > end)
        {
          std::swap(begin, end);
        }
        covered.emplace_back(std::max(pair.first_begin, begin),
                             std::min(pair.first_end, end));
      }
      std::sort(covered.begin(), covered.end());
      const double tolerance = 1.0e-10 * std::max(group.matching_radius, pair.first_end);
      double end = pair.first_begin;
      for (const auto &[begin, interval_end] : covered)
      {
        if (begin > end + tolerance)
        {
          return false;
        }
        end = std::max(end, interval_end);
      }
      return end >= pair.first_end - tolerance;
    };
    for (std::size_t pair_index = 0; pair_index < candidate_pairs.size(); pair_index++)
    {
      const auto &pair = candidate_pairs[pair_index];
      if (!candidate_pair_has_model[pair_index] && !PairCoveredByParallelCluster(pair))
      {
        RejectInteraction(
            pair.first, pair.second, missing_library_interactions,
            "The fabrication-process response library has no model for a nearby "
            "three-dimensional edge pair outside an exact parallel-edge cluster!",
            candidate_pair_topologies[pair_index]);
      }
    }

    std::vector<std::vector<std::pair<std::pair<double, double>, std::size_t>>>
        candidate_intervals(segments.size());
    for (std::size_t pair_index = 0; pair_index < candidate_pairs.size(); pair_index++)
    {
      const auto &pair = candidate_pairs[pair_index];
      candidate_intervals[pair.first].push_back(
          {{pair.first_begin, pair.first_end}, pair_index});
      candidate_intervals[pair.second].push_back(
          {{pair.second_begin, pair.second_end}, pair_index});
    }
    for (auto &intervals : candidate_intervals)
    {
      std::sort(intervals.begin(), intervals.end());
      for (std::size_t i = 1; i < intervals.size(); i++)
      {
        const double tolerance =
            1.0e-10 * std::max(group.matching_radius, intervals[i - 1].first.second);
        if (intervals[i].first.first >= intervals[i - 1].first.second - tolerance)
        {
          continue;
        }
        multiedge_interactions++;
      }
    }
    if (!unsupported_segments.empty())
    {
      Mpi::Warning(
          "Omitting {} of {} three-dimensional target edge segments in unsupported local "
          "interaction neighborhoods (nonparallel: {}, incompatible process normal: {}, "
          "process-normal offset: {}, unclassified topology: {}, missing library model: "
          "{}, multi-edge: {}). Unclassified pairs by conductor ownership: same = {}, "
          "different = {}.\n",
          static_cast<int>(unsupported_segments.size()), static_cast<int>(segments.size()),
          nonparallel_interactions, incompatible_process_interactions,
          process_offset_interactions, unclassified_interactions,
          missing_library_interactions, multiedge_interactions,
          unclassified_same_conductor_interactions,
          unclassified_different_conductor_interactions);
      std::vector<EdgeSegment3D> supported_segments;
      supported_segments.reserve(segments.size() - unsupported_segments.size());
      for (std::size_t i = 0; i < segments.size(); i++)
      {
        if (unsupported_segments.find(i) == unsupported_segments.end())
        {
          supported_segments.push_back(segments[i]);
        }
      }
      segments = std::move(supported_segments);
      local_indices.clear();
      for (std::size_t i = 0; i < segments.size(); i++)
      {
        local_indices.emplace(segments[i].geometry_index, i);
      }
    }
    if (segments.empty())
    {
      unmatched_groups++;
      continue;
    }
    for (const auto &segment : segments)
    {
      group_selected_length += segment.length;
    }

    if (diagnostics)
    {
      for (std::size_t i = 0; i < segments.size(); i++)
      {
        const auto &segment = segments[i];
        const auto &source = geometry.segments[segment.geometry_index];
        double corner_length = 0.0;
        for (const std::size_t vertex : source.vertices)
        {
          if (geometry.vertices[vertex].physical_type &&
              *geometry.vertices[vertex].physical_type != MetalEdgeVertexType::REGULAR &&
              !geometry.vertices[vertex].on_truncation_boundary)
          {
            corner_length += std::min(group.matching_radius, segment.length);
          }
        }
        group_corner_neighborhood_length += std::min(segment.length, corner_length);
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
      nonregular_vertices += type && *type != MetalEdgeVertexType::REGULAR &&
                             !geometry.vertices[vertex].on_truncation_boundary;
    }
    std::vector<PendingPatch> pending;
    std::vector<std::vector<std::pair<double, double>>> vertex_excluded_intervals(
        segments.size());
    bool group_has_unmatched_rounded_corner = false;
    std::set<std::size_t> spatially_excluded_vertices;
    std::vector<SpatialClusterSelection3D::InteractionNeighborhood>
        active_spatial_interactions;
    auto ExcludeBeyondVertex =
        [&](std::size_t vertex, std::size_t previous_segment, double distance)
    {
      const double tolerance = 1.0e-10 * group.matching_radius;
      const int physical_chain =
          geometry.segments[segments[previous_segment].geometry_index].physical_chain;
      std::set<std::size_t> visited = {previous_segment};
      double remaining = distance;
      while (remaining > tolerance)
      {
        spatially_excluded_vertices.insert(vertex);
        std::optional<std::size_t> next_segment;
        for (const std::size_t geometry_segment : geometry.vertices[vertex].segments)
        {
          auto local = local_indices.find(geometry_segment);
          if (local == local_indices.end() || !visited.insert(local->second).second ||
              geometry.segments[geometry_segment].physical_chain != physical_chain)
          {
            continue;
          }
          MFEM_VERIFY(!next_segment,
                      "Metal physical chain branches within a spatial coupon interval!");
          next_segment = local->second;
        }
        if (!next_segment)
        {
          break;
        }

        const auto &segment = segments[*next_segment];
        const auto &source = geometry.segments[segment.geometry_index];
        MFEM_VERIFY(source.vertices[0] == vertex || source.vertices[1] == vertex,
                    "Inconsistent spatial coupon chain connectivity!");
        const double trim = std::min(remaining, segment.length);
        if (source.vertices[0] == vertex)
        {
          vertex_excluded_intervals[*next_segment].emplace_back(0.0, trim);
        }
        else
        {
          vertex_excluded_intervals[*next_segment].emplace_back(segment.length - trim,
                                                                segment.length);
        }
        remaining -= trim;
        if (trim < segment.length - tolerance)
        {
          break;
        }
        vertex = source.vertices[0] == vertex ? source.vertices[1] : source.vertices[0];
        previous_segment = *next_segment;
      }
    };
    auto ExcludeSpatialInterval =
        [&](const SpatialClusterSelection3D &selection, std::size_t model_edge_index)
    {
      const auto &model = library.models[selection.response.models.front().index];
      const auto &edge = model.spatial_edges[model_edge_index];
      const auto &site = selection.sites[selection.model_to_site[model_edge_index]];
      const auto local = local_indices.find(site.geometry_index);
      if (local == local_indices.end())
      {
        return;
      }
      const std::size_t segment_index = local->second;
      const auto &segment = segments[segment_index];
      const auto &source = geometry.segments[segment.geometry_index];
      const Point3D tangent = TransformLocalVector(
          selection.axes, Normalize(Cross(edge.gap_direction, edge.process_normal)));
      const double orientation = Dot(tangent, segment.tangent);
      MFEM_VERIFY(std::abs(orientation) > 1.0 - 1.0e-8,
                  "Matched spatial edge tangent is incompatible with the target edge!");
      const double center = std::clamp(
          Dot(Subtract(site.point, segment.p0), segment.tangent), 0.0, segment.length);
      double begin = center + orientation * edge.interval[0];
      double end = center + orientation * edge.interval[1];
      if (begin > end)
      {
        std::swap(begin, end);
      }
      vertex_excluded_intervals[segment_index].emplace_back(std::max(0.0, begin),
                                                            std::min(segment.length, end));
      const double tolerance = 1.0e-10 * group.matching_radius;
      if (begin <= tolerance)
      {
        spatially_excluded_vertices.insert(source.vertices[0]);
      }
      if (end >= segment.length - tolerance)
      {
        spatially_excluded_vertices.insert(source.vertices[1]);
      }
      if (begin < 0.0)
      {
        ExcludeBeyondVertex(source.vertices[0], segment_index, -begin);
      }
      if (end > segment.length)
      {
        ExcludeBeyondVertex(source.vertices[1], segment_index, end - segment.length);
      }
    };

    for (const auto &selection : cross_interface_selections)
    {
      const auto &model = library.models[selection.response.models.front().index];
      for (std::size_t edge = 0; edge < model.spatial_edges.size(); edge++)
      {
        ExcludeSpatialInterval(selection, edge);
      }
    }

    int matched_spatial_clusters = 0;
    for (auto &selection : spatial_cluster_selections)
    {
      const bool retained = std::all_of(
          selection.sites.begin(), selection.sites.end(), [&](const auto &site)
          { return local_indices.find(site.geometry_index) != local_indices.end(); });
      if (!retained)
      {
        continue;
      }
      for (auto &site : selection.sites)
      {
        site.segment = local_indices.at(site.geometry_index);
      }

      if (requirements)
      {
        requirements->Add(
            3, LibraryTopology::SPATIAL_EDGE_CLUSTER, selection.targets_by_slot,
            selection.sites.front().boundary_condition, SpatialGeometry(selection), library,
            &selection.response, SpatialLength(selection));
      }
      const auto &weighted_model = selection.response.models.front();
      ResponsePatchData patch;
      patch.origin = selection.origin;
      patch.axis_u = selection.axes[0];
      patch.axis_v = selection.axes[1];
      patch.axis_w = selection.axes[2];
      patch.conductor_references = selection.response.conductor_references;
      patch.weight = 1.0;
      patch.maxwell_reference_is_pec = std::all_of(
          selection.sites.begin(), selection.sites.end(), [](const auto &site)
          { return site.boundary_condition.type == MetalBoundaryConditionType::PEC; });
      for (const auto &reference : patch.conductor_references)
      {
        patch.maxwell_conductor_anchors.push_back(
            TransformLocalPoint(patch.origin, selection.axes, reference));
      }
      pending.push_back({weighted_model.index, std::move(patch)});
      for (std::size_t edge = 0;
           edge < library.models[weighted_model.index].spatial_edges.size(); edge++)
      {
        ExcludeSpatialInterval(selection, edge);
      }
      active_spatial_interactions.insert(active_spatial_interactions.end(),
                                         selection.interactions.begin(),
                                         selection.interactions.end());
      matched_spatial_clusters++;
    }
    auto ExcludeVertexArm =
        [&](std::size_t vertex, std::size_t segment_index, double distance)
    {
      const int physical_chain =
          geometry.segments[segments[segment_index].geometry_index].physical_chain;
      std::set<std::size_t> visited;
      double remaining = distance;
      const double tolerance = 1.0e-10 * group.matching_radius;
      while (remaining > tolerance)
      {
        MFEM_VERIFY(visited.insert(segment_index).second,
                    "Metal physical chain contains a cycle within a vertex neighborhood!");
        const auto &segment = segments[segment_index];
        const auto &source = geometry.segments[segment.geometry_index];
        MFEM_VERIFY(source.physical_chain == physical_chain &&
                        (source.vertices[0] == vertex || source.vertices[1] == vertex),
                    "Inconsistent metal physical-chain connectivity!");

        const double trim = std::min(remaining, segment.length);
        if (source.vertices[0] == vertex)
        {
          vertex_excluded_intervals[segment_index].emplace_back(0.0, trim);
        }
        else
        {
          vertex_excluded_intervals[segment_index].emplace_back(segment.length - trim,
                                                                segment.length);
        }
        remaining -= trim;
        if (trim < segment.length - tolerance)
        {
          break;
        }

        const std::size_t next_vertex =
            source.vertices[0] == vertex ? source.vertices[1] : source.vertices[0];
        std::optional<std::size_t> next_segment;
        for (const std::size_t geometry_segment : geometry.vertices[next_vertex].segments)
        {
          auto local = local_indices.find(geometry_segment);
          if (local == local_indices.end() || local->second == segment_index ||
              geometry.segments[geometry_segment].physical_chain != physical_chain)
          {
            continue;
          }
          MFEM_VERIFY(!next_segment,
                      "Metal physical chain branches within a vertex neighborhood!");
          next_segment = local->second;
        }
        if (!next_segment)
        {
          break;
        }
        vertex = next_vertex;
        segment_index = *next_segment;
      }
    };
    int matched_corners = 0;
    int matched_sharp_corners = 0;
    int matched_endpoints = 0;
    int matched_junctions = 0;
    std::set<std::size_t> modeled_curved_vertices;
    std::map<int, std::vector<std::size_t>> physical_chains;
    for (std::size_t i = 0; i < segments.size(); i++)
    {
      const int chain = geometry.segments[segments[i].geometry_index].physical_chain;
      if (chain >= 0)
      {
        physical_chains[chain].push_back(i);
      }
    }
    for (const auto &[chain, chain_segments] : physical_chains)
    {
      (void)chain;
      if (chain_segments.size() < 3)
      {
        continue;
      }
      std::map<std::size_t, std::vector<std::size_t>> adjacency;
      for (const std::size_t segment_index : chain_segments)
      {
        const auto &source = geometry.segments[segments[segment_index].geometry_index];
        for (const std::size_t vertex : source.vertices)
        {
          adjacency[vertex].push_back(segment_index);
        }
      }
      if (std::any_of(adjacency.begin(), adjacency.end(),
                      [](const auto &entry) { return entry.second.size() > 2; }))
      {
        continue;
      }

      auto endpoint = std::find_if(adjacency.begin(), adjacency.end(), [](const auto &entry)
                                   { return entry.second.size() == 1; });
      const std::size_t start_vertex =
          endpoint != adjacency.end() ? endpoint->first : adjacency.begin()->first;
      std::vector<std::size_t> ordered_segments;
      std::vector<std::size_t> ordered_vertices = {start_vertex};
      std::optional<std::size_t> previous_segment;
      std::size_t current_vertex = start_vertex;
      while (ordered_segments.size() < chain_segments.size())
      {
        const auto adjacent = adjacency.find(current_vertex);
        if (adjacent == adjacency.end())
        {
          break;
        }
        auto next = std::find_if(
            adjacent->second.begin(), adjacent->second.end(), [&](std::size_t segment_index)
            { return !previous_segment || segment_index != *previous_segment; });
        if (next == adjacent->second.end())
        {
          break;
        }
        const std::size_t segment_index = *next;
        ordered_segments.push_back(segment_index);
        const auto &source = geometry.segments[segments[segment_index].geometry_index];
        const std::size_t next_vertex =
            source.vertices[0] == current_vertex ? source.vertices[1] : source.vertices[0];
        ordered_vertices.push_back(next_vertex);
        previous_segment = segment_index;
        current_vertex = next_vertex;
        if (current_vertex == start_vertex)
        {
          break;
        }
      }
      if (ordered_segments.size() != chain_segments.size())
      {
        continue;
      }
      const bool cycle = ordered_vertices.back() == start_vertex;
      if (cycle)
      {
        ordered_vertices.pop_back();
      }
      const std::size_t vertex_count = ordered_vertices.size();
      if (vertex_count < 3)
      {
        continue;
      }

      std::vector<double> turning_angles(vertex_count);
      std::vector<bool> curved(vertex_count, false);
      for (std::size_t i = 0; i < vertex_count; i++)
      {
        if ((!cycle && (i == 0 || i + 1 == vertex_count)) ||
            geometry.vertices[ordered_vertices[i]].physical_type !=
                MetalEdgeVertexType::REGULAR)
        {
          continue;
        }
        const std::size_t previous = i == 0 ? vertex_count - 1 : i - 1;
        const std::size_t next = i + 1 == vertex_count ? 0 : i + 1;
        const Point3D incoming =
            Normalize(Subtract(geometry.vertices[ordered_vertices[i]].coordinate,
                               geometry.vertices[ordered_vertices[previous]].coordinate));
        const Point3D outgoing =
            Normalize(Subtract(geometry.vertices[ordered_vertices[next]].coordinate,
                               geometry.vertices[ordered_vertices[i]].coordinate));
        turning_angles[i] = std::acos(std::clamp(Dot(incoming, outgoing), -1.0, 1.0));
        curved[i] = turning_angles[i] > 1.0e-6;
      }

      auto MatchRoundedRun = [&](const std::vector<std::size_t> &run)
      {
        if (run.size() < 2)
        {
          return;
        }
        if (std::any_of(run.begin(), run.end(),
                        [&](std::size_t index)
                        {
                          const std::size_t vertex = ordered_vertices[index];
                          return geometry.vertices[vertex].on_truncation_boundary ||
                                 spatially_excluded_vertices.find(vertex) !=
                                     spatially_excluded_vertices.end();
                        }))
        {
          return;
        }
        const std::size_t first_vertex_index = run.front();
        const std::size_t last_vertex_index = run.back();
        if (!cycle && (first_vertex_index == 0 || last_vertex_index + 1 >= vertex_count))
        {
          return;
        }

        const std::size_t previous_vertex_index =
            first_vertex_index == 0 ? vertex_count - 1 : first_vertex_index - 1;
        const std::size_t next_vertex_index =
            last_vertex_index + 1 == vertex_count ? 0 : last_vertex_index + 1;
        const std::size_t incoming_segment_index =
            first_vertex_index == 0 ? ordered_segments.size() - 1 : first_vertex_index - 1;
        const std::size_t outgoing_segment_index =
            last_vertex_index == ordered_segments.size() ? 0 : last_vertex_index;
        if (incoming_segment_index >= ordered_segments.size() ||
            outgoing_segment_index >= ordered_segments.size())
        {
          return;
        }

        std::vector<std::size_t> arc_segments;
        std::size_t arc_index = first_vertex_index;
        while (arc_index != last_vertex_index)
        {
          if (arc_index >= ordered_segments.size())
          {
            arc_index = 0;
          }
          arc_segments.push_back(ordered_segments[arc_index]);
          arc_index = (arc_index + 1) % ordered_segments.size();
          if (arc_segments.size() > ordered_segments.size())
          {
            return;
          }
        }
        if (arc_segments.empty())
        {
          return;
        }

        const double total_angle = std::accumulate(
            run.begin(), run.end(), 0.0,
            [&](double value, std::size_t index) { return value + turning_angles[index]; });
        if (total_angle <= 1.0e-3)
        {
          return;
        }

        const std::size_t start_vertex_index = ordered_vertices[first_vertex_index];
        const std::size_t end_vertex_index = ordered_vertices[last_vertex_index];
        const Point3D start = geometry.vertices[start_vertex_index].coordinate;
        const Point3D end = geometry.vertices[end_vertex_index].coordinate;
        const Point3D incoming = Normalize(Subtract(
            start, geometry.vertices[ordered_vertices[previous_vertex_index]].coordinate));
        const Point3D outgoing = Normalize(Subtract(
            geometry.vertices[ordered_vertices[next_vertex_index]].coordinate, end));
        Point3D first_direction = Scale(-1.0, incoming);
        Point3D second_direction = outgoing;
        std::size_t first_segment_index = ordered_segments[incoming_segment_index];
        std::size_t second_segment_index = ordered_segments[outgoing_segment_index];
        const auto &first_segment = segments[first_segment_index];
        const auto &second_segment = segments[second_segment_index];
        if (maxwell && !SameBoundaryLaw(first_segment.boundary_condition,
                                        second_segment.boundary_condition))
        {
          return;
        }
        const auto boundary_condition =
            maxwell ? first_segment.boundary_condition : MetalBoundaryLaw{};
        const double corner_score = Dot(first_segment.axis_u, second_direction) +
                                    Dot(second_segment.axis_u, first_direction);
        if (std::abs(corner_score) <= 1.0e-8)
        {
          return;
        }
        const LibraryTopology topology = corner_score < 0.0
                                             ? LibraryTopology::CONVEX_CORNER
                                             : LibraryTopology::CONCAVE_CORNER;
        const double angle =
            std::acos(std::clamp(Dot(first_direction, second_direction), -1.0, 1.0));
        const double expected_turn = std::acos(-1.0) - angle;
        if (std::abs(total_angle - expected_turn) > 1.0e-3)
        {
          return;
        }

        MFEM_VERIFY(Dot(first_segment.axis_v, second_segment.axis_v) > 0.95,
                    "A rounded-corner response model requires compatible process "
                    "normals on both incident arms!");
        Point3D process_normal =
            Normalize(Add(first_segment.axis_v, second_segment.axis_v));
        const double denominator =
            Dot(Cross(first_direction, second_direction), process_normal);
        if (std::abs(denominator) <= 1.0e-8)
        {
          return;
        }
        const double first_offset =
            Dot(Cross(Subtract(end, start), second_direction), process_normal) /
            denominator;
        Point3D origin = Add(start, Scale(first_offset, first_direction));
        const double first_tangent_distance = Distance(origin, start);
        const double second_tangent_distance = Distance(origin, end);
        const double distance_tolerance = 1.0e-6 * group.matching_radius;
        if (first_tangent_distance >= group.matching_radius + distance_tolerance ||
            second_tangent_distance >= group.matching_radius + distance_tolerance)
        {
          return;
        }
        const double tangent_scale =
            std::max(first_tangent_distance, second_tangent_distance);
        if (std::abs(first_tangent_distance - second_tangent_distance) >
            std::max(distance_tolerance, 0.05 * tangent_scale))
        {
          return;
        }
        const double corner_radius = 0.5 *
                                     (first_tangent_distance + second_tangent_distance) *
                                     std::tan(0.5 * angle);
        if (!(corner_radius > 0.0 && corner_radius < group.matching_radius))
        {
          return;
        }
        const auto model_selection = FindCornerLibraryModel(
            library, topology, angle, corner_radius, boundary_condition);
        if (!model_selection)
        {
          if (requirements)
          {
            requirements->Add(
                3, topology, group.targets, boundary_condition,
                {{"AngleDegrees",
                  requirements->SnapAngleDegrees(angle * 180.0 / std::acos(-1.0))},
                 {"CornerRadius", requirements->ScaleLength(corner_radius)}},
                library, nullptr, 2.0 * group.matching_radius,
                "No compatible rounded-corner model or interpolation bracket");
          }
          unmatched_rounded_corners++;
          group_has_unmatched_rounded_corner = true;
          return;
        }
        if (requirements)
        {
          requirements->Add(3, topology, group.targets, boundary_condition,
                            {{"AngleDegrees", requirements->SnapAngleDegrees(
                                                  angle * 180.0 / std::acos(-1.0))},
                             {"CornerRadius", requirements->ScaleLength(corner_radius)}},
                            library, &*model_selection, 2.0 * group.matching_radius);
        }

        if (Dot(Cross(process_normal, first_direction), second_direction) < 0.0)
        {
          std::swap(first_segment_index, second_segment_index);
          std::swap(first_direction, second_direction);
        }
        ResponsePatchData patch;
        patch.origin = origin;
        patch.axis_u = first_direction;
        patch.axis_v = Normalize(Cross(process_normal, first_direction));
        patch.axis_w = process_normal;
        patch.conductor_references = model_selection->conductor_references;
        patch.weight = 1.0;
        patch.maxwell_reference_is_pec =
            boundary_condition.type == MetalBoundaryConditionType::PEC;
        for (const auto &reference : patch.conductor_references)
        {
          auto anchor = patch.origin;
          for (int d = 0; d < 3; d++)
          {
            anchor[d] += reference[0] * patch.axis_u[d] + reference[1] * patch.axis_v[d] +
                         reference[2] * patch.axis_w[d];
          }
          patch.maxwell_conductor_anchors.push_back(anchor);
        }
        for (const auto &weighted_model : model_selection->models)
        {
          auto weighted_patch = patch;
          weighted_patch.weight = weighted_model.weight;
          pending.push_back({weighted_model.index, std::move(weighted_patch)});
        }

        for (const std::size_t segment_index : arc_segments)
        {
          vertex_excluded_intervals[segment_index].emplace_back(
              0.0, segments[segment_index].length);
        }
        ExcludeVertexArm(start_vertex_index, ordered_segments[incoming_segment_index],
                         std::max(0.0, group.matching_radius - first_tangent_distance));
        ExcludeVertexArm(end_vertex_index, ordered_segments[outgoing_segment_index],
                         std::max(0.0, group.matching_radius - second_tangent_distance));
        for (const std::size_t vertex_index : run)
        {
          modeled_curved_vertices.insert(ordered_vertices[vertex_index]);
        }
        group_maximum_library_distance =
            std::max(group_maximum_library_distance, model_selection->normalized_distance);
        matched_corners++;
        interpolated_rounded_corners += model_selection->IsInterpolated();
      };

      std::vector<std::size_t> run;
      auto FlushRun = [&]()
      {
        MatchRoundedRun(run);
        run.clear();
      };
      if (cycle)
      {
        auto straight = std::find(curved.begin(), curved.end(), false);
        if (straight == curved.end())
        {
          continue;
        }
        const std::size_t start = std::distance(curved.begin(), straight);
        for (std::size_t step = 1; step <= vertex_count; step++)
        {
          const std::size_t index = (start + step) % vertex_count;
          if (curved[index])
          {
            run.push_back(index);
          }
          else
          {
            FlushRun();
          }
        }
      }
      else
      {
        for (std::size_t i = 0; i < vertex_count; i++)
        {
          if (curved[i])
          {
            run.push_back(i);
          }
          else
          {
            FlushRun();
          }
        }
        FlushRun();
      }
    }
    for (const std::size_t vertex : vertices)
    {
      if (geometry.vertices[vertex].physical_type != MetalEdgeVertexType::CORNER ||
          geometry.vertices[vertex].on_truncation_boundary ||
          spatially_excluded_vertices.find(vertex) != spatially_excluded_vertices.end())
      {
        continue;
      }
      std::vector<std::size_t> incident;
      for (const std::size_t geometry_segment : geometry.vertices[vertex].segments)
      {
        auto local = local_indices.find(geometry_segment);
        if (local != local_indices.end())
        {
          incident.push_back(local->second);
        }
      }
      if (incident.size() != 2)
      {
        continue;
      }
      auto DirectionAway = [&](std::size_t segment_index)
      {
        const auto &source = geometry.segments[segments[segment_index].geometry_index];
        return source.vertices[0] == vertex ? segments[segment_index].tangent
                                            : Scale(-1.0, segments[segment_index].tangent);
      };
      Point3D first_direction = DirectionAway(incident[0]);
      Point3D second_direction = DirectionAway(incident[1]);
      const auto &first_segment = segments[incident[0]];
      const auto &second_segment = segments[incident[1]];
      if (maxwell && !SameBoundaryLaw(first_segment.boundary_condition,
                                      second_segment.boundary_condition))
      {
        continue;
      }
      const auto boundary_condition =
          maxwell ? first_segment.boundary_condition : MetalBoundaryLaw{};
      const double corner_score = Dot(first_segment.axis_u, second_direction) +
                                  Dot(second_segment.axis_u, first_direction);
      if (std::abs(corner_score) <= 1.0e-8)
      {
        continue;
      }
      const LibraryTopology topology = corner_score < 0.0 ? LibraryTopology::CONVEX_CORNER
                                                          : LibraryTopology::CONCAVE_CORNER;
      const double angle =
          std::acos(std::clamp(Dot(first_direction, second_direction), -1.0, 1.0));
      const auto model_selection =
          FindCornerLibraryModel(library, topology, angle, 0.0, boundary_condition);
      if (!model_selection)
      {
        if (requirements)
        {
          requirements->Add(3, topology, group.targets, boundary_condition,
                            {{"AngleDegrees", requirements->SnapAngleDegrees(
                                                  angle * 180.0 / std::acos(-1.0))},
                             {"CornerRadius", 0.0}},
                            library, nullptr, 2.0 * group.matching_radius,
                            "No compatible sharp-corner model");
        }
        continue;
      }
      if (requirements)
      {
        requirements->Add(3, topology, group.targets, boundary_condition,
                          {{"AngleDegrees", requirements->SnapAngleDegrees(
                                                angle * 180.0 / std::acos(-1.0))},
                           {"CornerRadius", 0.0}},
                          library, &*model_selection, 2.0 * group.matching_radius);
      }
      MFEM_ASSERT(!model_selection->IsInterpolated(),
                  "A sharp corner cannot use radius interpolation!");

      MFEM_VERIFY(Dot(first_segment.axis_v, second_segment.axis_v) > 0.95,
                  "A corner-response model requires compatible process normals on both "
                  "incident arms!");
      Point3D process_normal = Normalize(Add(first_segment.axis_v, second_segment.axis_v));
      if (Dot(Cross(process_normal, first_direction), second_direction) < 0.0)
      {
        std::swap(incident[0], incident[1]);
        std::swap(first_direction, second_direction);
      }
      ResponsePatchData patch;
      patch.origin = geometry.vertices[vertex].coordinate;
      patch.axis_u = first_direction;
      patch.axis_v = Normalize(Cross(process_normal, first_direction));
      patch.axis_w = process_normal;
      const auto &weighted_model = model_selection->models.front();
      patch.conductor_references = model_selection->conductor_references;
      patch.weight = 1.0;
      patch.maxwell_reference_is_pec =
          boundary_condition.type == MetalBoundaryConditionType::PEC;
      for (const auto &reference : patch.conductor_references)
      {
        auto anchor = patch.origin;
        for (int d = 0; d < 3; d++)
        {
          anchor[d] += reference[0] * patch.axis_u[d] + reference[1] * patch.axis_v[d] +
                       reference[2] * patch.axis_w[d];
        }
        patch.maxwell_conductor_anchors.push_back(anchor);
      }
      pending.push_back({weighted_model.index, patch});

      for (const std::size_t segment_index : incident)
      {
        ExcludeVertexArm(vertex, segment_index, group.matching_radius);
      }
      group_maximum_library_distance =
          std::max(group_maximum_library_distance, model_selection->normalized_distance);
      matched_corners++;
      matched_sharp_corners++;
    }
    for (const std::size_t vertex : vertices)
    {
      const auto vertex_type = geometry.vertices[vertex].physical_type;
      if (!vertex_type ||
          (*vertex_type != MetalEdgeVertexType::ENDPOINT &&
           *vertex_type != MetalEdgeVertexType::JUNCTION) ||
          geometry.vertices[vertex].on_truncation_boundary ||
          spatially_excluded_vertices.find(vertex) != spatially_excluded_vertices.end())
      {
        continue;
      }

      std::vector<std::size_t> incident;
      for (const std::size_t geometry_segment : geometry.vertices[vertex].segments)
      {
        auto local = local_indices.find(geometry_segment);
        if (local != local_indices.end())
        {
          incident.push_back(local->second);
        }
      }
      const bool endpoint = *vertex_type == MetalEdgeVertexType::ENDPOINT;
      if ((endpoint && incident.size() != 1) || (!endpoint && incident.size() < 3))
      {
        continue;
      }

      auto DirectionAway = [&](std::size_t segment_index)
      {
        const auto &source = geometry.segments[segments[segment_index].geometry_index];
        return source.vertices[0] == vertex ? segments[segment_index].tangent
                                            : Scale(-1.0, segments[segment_index].tangent);
      };
      std::vector<Point3D> directions;
      directions.reserve(incident.size());
      Point3D process_normal{};
      bool compatible = true;
      const int conductor = segments[incident.front()].conductor;
      const auto boundary_condition = segments[incident.front()].boundary_condition;
      for (const std::size_t segment_index : incident)
      {
        const auto &segment = segments[segment_index];
        directions.push_back(DirectionAway(segment_index));
        compatible = compatible && segment.conductor == conductor &&
                     SameBoundaryLaw(segment.boundary_condition, boundary_condition);
        if (Norm(process_normal) == 0.0)
        {
          process_normal = segment.axis_v;
        }
        else
        {
          if (Dot(Normalize(process_normal), segment.axis_v) <= 0.95)
          {
            compatible = false;
          }
          process_normal = Add(process_normal, segment.axis_v);
        }
      }
      if (!compatible || Norm(process_normal) == 0.0)
      {
        continue;
      }
      process_normal = Normalize(process_normal);
      if (std::any_of(directions.begin(), directions.end(), [&](const auto &direction)
                      { return std::abs(Dot(direction, process_normal)) > 1.0e-8; }))
      {
        continue;
      }

      auto VertexAxes = [&](std::size_t first)
      {
        const Point3D axis_u = directions[first];
        Point3D axis_v = Normalize(Cross(process_normal, axis_u));
        if (endpoint && Dot(axis_v, segments[incident[first]].axis_u) < 0.0)
        {
          axis_v = Scale(-1.0, axis_v);
        }
        return std::array<Point3D, 3>{axis_u, axis_v, process_normal};
      };
      std::vector<std::optional<PlanViewGeometry>> vertex_plan_views(directions.size());
      if (surface.retain_faces)
      {
        std::vector<SpatialEdgeSite3D> vertex_sites;
        std::map<int, int> conductor_by_metal_component;
        for (std::size_t arm = 0; arm < incident.size(); arm++)
        {
          const auto &segment = segments[incident[arm]];
          const auto &source = geometry.segments[segment.geometry_index];
          const Point3D tangent = Normalize(Cross(segment.axis_u, segment.axis_v));
          const double orientation = Dot(directions[arm], tangent);
          MFEM_VERIFY(std::abs(std::abs(orientation) - 1.0) <= 1.0e-8,
                      "A vertex arm is inconsistent with its extracted edge frame!");
          vertex_sites.push_back(
              {source.physical_chain, segment.geometry_index, incident[arm], 0.0,
               orientation > 0.0 ? std::array<double, 2>{0.0, group.matching_radius}
                                 : std::array<double, 2>{-group.matching_radius, 0.0},
               geometry.vertices[vertex].coordinate, segment.axis_u, segment.axis_v,
               segment.conductor, segment.metal_component, segment.targets,
               segment.boundary_condition});
          if (segment.metal_component >= 0)
          {
            conductor_by_metal_component.emplace(segment.metal_component, 1);
          }
        }
        if (!conductor_by_metal_component.empty())
        {
          for (std::size_t first = 0; first < directions.size(); first++)
          {
            vertex_plan_views[first] =
                GatherPlanViewFacets(vertex_sites, geometry.vertices[vertex].coordinate,
                                     VertexAxes(first), conductor_by_metal_component, 2);
          }
        }
      }
      auto MatchesVertexPlanView = [&](const LibraryModel &model, std::size_t first)
      {
        if (!model.plan_view_boundary)
        {
          const bool exact_plan_view_available =
              first < vertex_plan_views.size() && vertex_plan_views[first] &&
              std::any_of(
                  vertex_plan_views[first]->facets.begin(),
                  vertex_plan_views[first]->facets.end(),
                  [](const auto &facet) { return facet.conductor == 1; });
          return !requirements || !exact_plan_view_available;
        }
        if (first >= vertex_plan_views.size() || !vertex_plan_views[first])
        {
          return false;
        }
        const auto &plan_view = *vertex_plan_views[first];
        const bool found_conductor =
            std::any_of(plan_view.facets.begin(), plan_view.facets.end(),
                        [](const auto &facet) { return facet.conductor == 1; });
        const auto clip_bounds = std::make_pair(plan_view.lower, plan_view.upper);
        const std::optional<decltype(clip_bounds)> classified_bounds =
            HasClassifiedPlanViewBoundary(*model.plan_view_boundary)
                ? std::optional<decltype(clip_bounds)>(clip_bounds)
                : std::nullopt;
        return found_conductor &&
               CanonicalPlanViewBoundary(plan_view.facets, library.matching_radius,
                                         plan_view.process_axis,
                                         classified_bounds) == *model.plan_view_boundary;
      };
      const LibraryTopology topology =
          endpoint ? LibraryTopology::ENDPOINT : LibraryTopology::JUNCTION;
      const auto model_selection =
          FindVertexLibraryModel(library, topology, directions, process_normal,
                                 boundary_condition, MatchesVertexPlanView);
      auto ArmAngles = [&](std::size_t first)
      {
        if (endpoint)
        {
          return std::vector<double>{};
        }
        const Point3D axis_u = directions[first];
        const Point3D axis_v = Normalize(Cross(process_normal, axis_u));
        std::vector<double> angles;
        for (const auto &direction : directions)
        {
          double angle = std::atan2(Dot(direction, axis_v), Dot(direction, axis_u));
          if (angle < 0.0)
          {
            angle += 2.0 * std::acos(-1.0);
          }
          angles.push_back(requirements->SnapAngleDegrees(angle * 180.0 / std::acos(-1.0)));
        }
        std::sort(angles.begin(), angles.end());
        return angles;
      };
      auto VertexGeometry = [&](std::size_t first)
      {
        const auto axes = VertexAxes(first);
        const auto &axis_u = axes[0];
        const auto &axis_v = axes[1];
        struct ArmDescription
        {
          double angle = 0.0;
          nlohmann::json data;
        };
        std::vector<ArmDescription> arm_descriptions;
        arm_descriptions.reserve(directions.size());
        for (std::size_t arm = 0; arm < directions.size(); arm++)
        {
          double angle =
              std::atan2(Dot(directions[arm], axis_v), Dot(directions[arm], axis_u));
          if (angle < 0.0)
          {
            angle += 2.0 * std::acos(-1.0);
          }
          const auto &segment = segments[incident[arm]];
          nlohmann::json arm_data = {
              {"Direction",
               {requirements->SnapDirection(Dot(directions[arm], axis_u)),
                requirements->SnapDirection(Dot(directions[arm], axis_v)), 0.0}},
              {"GapDirection",
               {requirements->SnapDirection(Dot(segment.axis_u, axis_u)),
                requirements->SnapDirection(Dot(segment.axis_u, axis_v)),
                requirements->SnapDirection(Dot(segment.axis_u, process_normal))}},
              {"ProcessNormal",
               {requirements->SnapDirection(Dot(segment.axis_v, axis_u)),
                requirements->SnapDirection(Dot(segment.axis_v, axis_v)),
                requirements->SnapDirection(Dot(segment.axis_v, process_normal))}},
              {"Interval", {0.0, requirements->ScaleLength(group.matching_radius)}},
              {"Conductor", 1},
              {"InterfaceSlot", 0},
              {"BoundaryCondition",
               requirements->DescribeBoundaryCondition(segment.boundary_condition)}};
          arm_descriptions.push_back({angle, std::move(arm_data)});
        }
        std::sort(arm_descriptions.begin(), arm_descriptions.end(),
                  [](const auto &first_arm, const auto &second_arm)
                  { return first_arm.angle < second_arm.angle; });
        nlohmann::json arms = nlohmann::json::array();
        for (auto &arm : arm_descriptions)
        {
          arms.push_back(std::move(arm.data));
        }
        nlohmann::json result = {{"SignatureVersion", 2},
                                 {"ArmCount", directions.size()},
                                 {"ArmAnglesDegrees", ArmAngles(first)},
                                 {"Arms", std::move(arms)}};
        if (first < vertex_plan_views.size() && vertex_plan_views[first])
        {
          const auto &plan_view = *vertex_plan_views[first];
          nlohmann::json facets = nlohmann::json::array();
          for (const auto &facet : plan_view.facets)
          {
            std::vector<Point3D> scaled(facet.points.size());
            for (std::size_t i = 0; i < facet.points.size(); i++)
            {
              for (int d = 0; d < 3; d++)
              {
                scaled[i][d] = requirements->ScaleLength(facet.points[i][d]);
              }
            }
            auto Sequence = [&](std::size_t start, bool reverse)
            {
              nlohmann::json points = nlohmann::json::array();
              for (std::size_t step = 0; step < scaled.size(); step++)
              {
                const std::size_t index =
                    reverse ? (start + scaled.size() - step) % scaled.size()
                            : (start + step) % scaled.size();
                points.push_back(scaled[index]);
              }
              return points;
            };
            nlohmann::json canonical;
            std::string canonical_key;
            for (std::size_t start = 0; start < scaled.size(); start++)
            {
              for (const bool reverse : {false, true})
              {
                auto candidate = Sequence(start, reverse);
                const std::string key = candidate.dump();
                if (canonical.is_null() || key < canonical_key)
                {
                  canonical = std::move(candidate);
                  canonical_key = key;
                }
              }
            }
            facets.push_back(
                {{"Conductor", facet.conductor}, {"Points", std::move(canonical)}});
          }
          std::sort(facets.begin(), facets.end(),
                    [](const auto &first_facet, const auto &second_facet)
                    { return first_facet.dump() < second_facet.dump(); });
          facets.erase(std::unique(facets.begin(), facets.end(),
                                   [](const auto &first_facet, const auto &second_facet)
                                   { return first_facet == second_facet; }),
                       facets.end());
          if (!facets.empty())
          {
            result["PlanViewFacets"] = std::move(facets);
            result["PlanViewBoundary"] = nlohmann::json::parse(CanonicalPlanViewBoundary(
                plan_view.facets, library.matching_radius, plan_view.process_axis,
                std::make_pair(plan_view.lower, plan_view.upper)));
          }
        }
        return result;
      };
      if (!model_selection)
      {
        if (requirements)
        {
          std::size_t canonical_first = 0;
          std::string canonical_key;
          for (std::size_t first = 0; first < directions.size(); first++)
          {
            const std::string key = VertexGeometry(first).dump();
            if (first == 0 || key < canonical_key)
            {
              canonical_first = first;
              canonical_key = key;
            }
          }
          requirements->Add(3, topology, group.targets, boundary_condition,
                            VertexGeometry(canonical_first), library, nullptr,
                            directions.size() * group.matching_radius,
                            endpoint ? "No compatible endpoint model"
                                     : "No compatible junction model");
        }
        continue;
      }
      if (requirements)
      {
        requirements->Add(3, topology, group.targets, boundary_condition,
                          VertexGeometry(model_selection->first_arm), library,
                          &model_selection->response,
                          directions.size() * group.matching_radius);
      }

      const std::size_t first_arm = model_selection->first_arm;
      ResponsePatchData patch;
      patch.origin = geometry.vertices[vertex].coordinate;
      patch.axis_u = directions[first_arm];
      patch.axis_v = Normalize(Cross(process_normal, patch.axis_u));
      if (endpoint && Dot(patch.axis_v, segments[incident[first_arm]].axis_u) < 0.0)
      {
        patch.axis_v = Scale(-1.0, patch.axis_v);
      }
      patch.axis_w = process_normal;
      patch.conductor_references = model_selection->response.conductor_references;
      patch.weight = 1.0;
      patch.maxwell_reference_is_pec =
          boundary_condition.type == MetalBoundaryConditionType::PEC;
      for (const auto &reference : patch.conductor_references)
      {
        auto anchor = patch.origin;
        if (patch.maxwell_reference_is_pec)
        {
          for (int d = 0; d < 3; d++)
          {
            anchor[d] += reference[0] * patch.axis_u[d] + reference[1] * patch.axis_v[d] +
                         reference[2] * patch.axis_w[d];
          }
        }
        patch.maxwell_conductor_anchors.push_back(anchor);
      }
      const auto &weighted_model = model_selection->response.models.front();
      pending.push_back({weighted_model.index, patch});

      for (const std::size_t segment_index : incident)
      {
        ExcludeVertexArm(vertex, segment_index, group.matching_radius);
      }
      group_maximum_library_distance = std::max(
          group_maximum_library_distance, model_selection->response.normalized_distance);
      if (endpoint)
      {
        matched_endpoints++;
      }
      else
      {
        matched_junctions++;
      }
    }
    if (diagnostics)
    {
      for (const std::size_t vertex : vertices)
      {
        if (geometry.vertices[vertex].physical_type != MetalEdgeVertexType::REGULAR ||
            modeled_curved_vertices.find(vertex) != modeled_curved_vertices.end())
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
        const double turning_angle = std::acos(std::clamp(-Dot(first, second), -1.0, 1.0));
        const double local_length = 0.5 * (incident[0]->length + incident[1]->length);
        if (local_length > 0.0)
        {
          group_maximum_curvature_ratio =
              std::max(group_maximum_curvature_ratio,
                       group.matching_radius * turning_angle / local_length);
        }
      }
    }
    for (auto intervals : vertex_excluded_intervals)
    {
      std::sort(intervals.begin(), intervals.end());
      double end = -mfem::infinity();
      for (const auto &[begin, interval_end] : intervals)
      {
        group_modeled_corner_neighborhood_length +=
            std::max(0.0, interval_end - std::max(begin, end));
        end = std::max(end, interval_end);
      }
    }
    bool group_matched = true;
    if (group_has_unmatched_rounded_corner && !requirements &&
        request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
    {
      group_matched = false;
    }

    std::vector<EdgePair3D> pairs;
    std::vector<std::vector<std::pair<double, double>>> paired_intervals(segments.size());
    std::vector<std::pair<std::size_t, std::size_t>> final_nearby_pairs;
    if (group_matched && !segments.empty())
    {
      const SegmentBoxIndex final_index(SegmentGeometry(segments));
      final_nearby_pairs = final_index.CandidatePairs(interaction_distance);
    }
    if (statistics)
    {
      statistics->pair_checks_patch_construction += final_nearby_pairs.size();
    }
    for (const auto &[i, j] : final_nearby_pairs)
    {
      const auto &first_source = geometry.segments[segments[i].geometry_index];
      const auto &second_source = geometry.segments[segments[j].geometry_index];
      if (first_source.physical_chain == second_source.physical_chain ||
          SegmentsShareVertex(first_source, second_source) ||
          SegmentDistanceSquared(segments[i].p0, segments[i].p1, segments[j].p0,
                                 segments[j].p1) >= interaction_distance_squared)
      {
        continue;
      }
      if (ConnectedNearVertex(i, j) ||
          IsSpatiallyMatched(i, j, active_spatial_interactions))
      {
        continue;
      }

      const double tangent_dot = Dot(segments[i].tangent, segments[j].tangent);
      if (std::abs(tangent_dot) < 1.0 - 1.0e-8)
      {
        Mpi::Warning("Nearby three-dimensional metal edges are not parallel; correction is "
                     "disabled for this interface group!\n");
        group_matched = false;
        break;
      }
      const double second_s0 =
          Dot(Subtract(segments[j].p0, segments[i].p0), segments[i].tangent);
      const double second_s1 =
          Dot(Subtract(segments[j].p1, segments[i].p0), segments[i].tangent);
      const double first_begin = std::max(0.0, std::min(second_s0, second_s1));
      const double first_end = std::min(segments[i].length, std::max(second_s0, second_s1));
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

    for (auto &intervals : paired_intervals)
    {
      std::sort(intervals.begin(), intervals.end());
    }
    auto parallel_clusters = FindParallelClusterSpans(library, segments, pairs);
    auto &parallel_cluster_spans = parallel_clusters.matched;
    const auto &unmatched_parallel_clusters = parallel_clusters.unmatched;
    auto DescribeParallelCluster = [&](const std::vector<std::size_t> &cluster,
                                       const Point3D &tangent, double coordinate,
                                       std::optional<Point3D> selected_axis = std::nullopt)
    {
      Point3D process_normal{};
      for (const std::size_t edge : cluster)
      {
        process_normal = Add(process_normal, segments[edge].axis_v);
      }
      process_normal = Normalize(process_normal);

      auto DescribeOrientation = [&](const Point3D &axis_u)
      {
        std::vector<std::size_t> ordered(cluster);
        std::sort(ordered.begin(), ordered.end(),
                  [&](std::size_t first, std::size_t second)
                  {
                    const auto first_point = InterpolateAtLongitudinalCoordinate(
                        segments[first], tangent, coordinate);
                    const auto second_point = InterpolateAtLongitudinalCoordinate(
                        segments[second], tangent, coordinate);
                    return Dot(first_point, axis_u) < Dot(second_point, axis_u);
                  });
        const Point3D origin = InterpolateAtLongitudinalCoordinate(
            segments[ordered.front()], tangent, coordinate);
        std::map<int, int> conductor_ids;
        nlohmann::json edges = nlohmann::json::array();
        for (const std::size_t edge_index : ordered)
        {
          const auto &edge = segments[edge_index];
          const Point3D point =
              InterpolateAtLongitudinalCoordinate(edge, tangent, coordinate);
          auto [conductor, inserted] =
              conductor_ids.emplace(edge.conductor, conductor_ids.size() + 1);
          (void)inserted;
          const Point3D offset = Subtract(point, origin);
          edges.push_back(
              {{"Offset",
                {requirements->ScaleLength(Dot(offset, axis_u)),
                 requirements->ScaleLength(Dot(offset, process_normal))}},
               {"GapDirection",
                {requirements->SnapDirection(Dot(edge.axis_u, axis_u)),
                 requirements->SnapDirection(Dot(edge.axis_u, process_normal))}},
               {"Conductor", conductor->second}});
        }
        return nlohmann::json{{"EdgeCount", edges.size()}, {"Edges", std::move(edges)}};
      };

      if (selected_axis)
      {
        return DescribeOrientation(*selected_axis);
      }
      const Point3D axis_u = Normalize(Cross(process_normal, tangent));
      auto forward = DescribeOrientation(axis_u);
      auto reverse = DescribeOrientation(Scale(-1.0, axis_u));
      return forward.dump() <= reverse.dump() ? forward : reverse;
    };
    if (!unmatched_parallel_clusters.empty())
    {
      Mpi::Warning(
          "Fabrication-process response library \"{}\" has no matching "
          "ParallelEdgeCluster model for {} three-dimensional longitudinal span(s); "
          "correction is disabled for this interface group!\n",
          library.name, unmatched_parallel_clusters.size());
      if (requirements)
      {
        for (const auto &span : unmatched_parallel_clusters)
        {
          requirements->Add(3, LibraryTopology::PARALLEL_EDGE_CLUSTER, group.targets,
                            segments[span.edges.front()].boundary_condition,
                            DescribeParallelCluster(span.edges, span.tangent,
                                                    0.5 * (span.begin + span.end)),
                            library, nullptr, span.end - span.begin,
                            "No compatible parallel-edge cluster model");
        }
      }
      group_matched = false;
    }
    for (const auto &span : parallel_cluster_spans)
    {
      group_maximum_library_distance = std::max(
          group_maximum_library_distance, span.selection.response.normalized_distance);
    }

    auto SubtractInterval = [](std::vector<std::pair<double, double>> &intervals,
                               double excluded_begin, double excluded_end)
    {
      std::vector<std::pair<double, double>> remainder;
      for (const auto &[begin, end] : intervals)
      {
        if (excluded_end <= begin || excluded_begin >= end)
        {
          remainder.emplace_back(begin, end);
          continue;
        }
        if (excluded_begin > begin)
        {
          remainder.emplace_back(begin, std::min(excluded_begin, end));
        }
        if (excluded_end < end)
        {
          remainder.emplace_back(std::max(excluded_end, begin), end);
        }
      }
      intervals = std::move(remainder);
    };

    auto AppendQuadrature = [&](const LibrarySelection &model_selection,
                                std::size_t first_index, double begin, double end,
                                std::optional<std::size_t> second_index)
    {
      MFEM_ASSERT(!model_selection.models.empty(), "Missing response-model selection!");
      const auto &representative = library.models[model_selection.models.front().index];
      const auto &first = segments[first_index];
      const EdgeSegment3D *second = second_index ? &segments[*second_index] : nullptr;
      std::vector<std::pair<double, double>> intervals = {{begin, end}};
      auto SubtractIntervals = [&](const std::vector<std::pair<double, double>> &excluded)
      {
        for (const auto &[excluded_begin, excluded_end] : excluded)
        {
          std::vector<std::pair<double, double>> remainder;
          for (const auto &[interval_begin, interval_end] : intervals)
          {
            if (excluded_end <= interval_begin || excluded_begin >= interval_end)
            {
              remainder.emplace_back(interval_begin, interval_end);
              continue;
            }
            if (excluded_begin > interval_begin)
            {
              remainder.emplace_back(interval_begin,
                                     std::min(excluded_begin, interval_end));
            }
            if (excluded_end < interval_end)
            {
              remainder.emplace_back(std::max(excluded_end, interval_begin), interval_end);
            }
          }
          intervals = std::move(remainder);
        }
      };
      SubtractIntervals(vertex_excluded_intervals[first_index]);
      if (second_index)
      {
        std::vector<std::pair<double, double>> mapped;
        for (const auto &[second_begin, second_end] :
             vertex_excluded_intervals[*second_index])
        {
          const double first_begin =
              Dot(Subtract(Interpolate(*second, second_begin), first.p0), first.tangent);
          const double first_end =
              Dot(Subtract(Interpolate(*second, second_end), first.p0), first.tangent);
          mapped.emplace_back(std::min(first_begin, first_end),
                              std::max(first_begin, first_end));
        }
        SubtractIntervals(mapped);
      }
      for (const auto &[interval_begin, interval_end] : intervals)
      {
        if (interval_end <= interval_begin)
        {
          continue;
        }
        if (requirements)
        {
          nlohmann::json geometry = {{"EdgeCount", second ? 2 : 1}};
          if (second)
          {
            const Point3D first_mid =
                Interpolate(first, 0.5 * (interval_begin + interval_end));
            double second_distance = Dot(Subtract(first_mid, second->p0), second->tangent);
            second_distance = std::clamp(second_distance, 0.0, second->length);
            geometry["Separation"] = requirements->ScaleLength(
                Distance(first_mid, Interpolate(*second, second_distance)));
          }
          requirements->Add(3, representative.topology, group.targets,
                            first.boundary_condition, geometry, library, &model_selection,
                            interval_end - interval_begin);
          continue;
        }
        for (int q = 0; q < quadrature.GetNPoints(); q++)
        {
          const auto &ip = quadrature.IntPoint(q);
          const double first_distance =
              interval_begin + (interval_end - interval_begin) * ip.x;
          const Point3D first_point = Interpolate(first, first_distance);
          std::optional<Point3D> paired_point;
          ResponsePatchData patch;
          if (second)
          {
            double second_distance =
                Dot(Subtract(first_point, second->p0), second->tangent);
            second_distance = std::clamp(second_distance, 0.0, second->length);
            paired_point = Interpolate(*second, second_distance);
            const Point3D direction = Normalize(Subtract(*paired_point, first_point));
            patch.origin = Scale(0.5, Add(first_point, *paired_point));
            patch.axis_u = direction;
            patch.axis_v = Normalize(Add(first.axis_v, second->axis_v));
          }
          else
          {
            patch.origin = first_point;
            patch.axis_u = first.axis_u;
            patch.axis_v = first.axis_v;
          }
          patch.maxwell_reference_is_pec =
              first.boundary_condition.type == MetalBoundaryConditionType::PEC &&
              (!second ||
               second->boundary_condition.type == MetalBoundaryConditionType::PEC);
          patch.conductor_references = model_selection.conductor_references;
          if (!patch.maxwell_reference_is_pec)
          {
            patch.maxwell_conductor_anchors = {first_point};
          }
          else if (representative.topology == LibraryTopology::SAME_CONDUCTOR_STRIP)
          {
            patch.maxwell_conductor_anchors = {patch.origin};
          }
          else
          {
            patch.maxwell_conductor_anchors = {
                Add(first_point, Scale(-group.matching_radius, first.axis_u))};
          }
          if (patch.conductor_references.size() > 1)
          {
            MFEM_VERIFY(paired_point && patch.conductor_references.size() == 2,
                        "A paired-edge response model requires exactly two conductor "
                        "references!");
            // The coupon references lie inside each metal, but a single line between
            // them crosses two metal-gap discontinuities. Use the physical edge points
            // as local conductor anchors so Maxwell quadrature spans only the dielectric
            // gap. For finite impedance this is a local quasi-electrostatic reference.
            patch.maxwell_conductor_anchors = {first_point, *paired_point};
          }
          const double quadrature_weight = (interval_end - interval_begin) * ip.weight;
          for (const auto &weighted_model : model_selection.models)
          {
            const auto &source = library.models[weighted_model.index];
            MFEM_VERIFY(
                source.coupon_depth > 0.0,
                "Three-dimensional response correction requires CouponDepth for every "
                "selected fabrication-process response model!");
            auto weighted_patch = patch;
            weighted_patch.weight =
                weighted_model.weight * quadrature_weight / source.coupon_depth;
            pending.push_back({weighted_model.index, std::move(weighted_patch)});
          }
        }
      }
    };

    auto AppendParallelClusterQuadrature = [&](const ParallelClusterSpan3D &span)
    {
      const auto &selection = span.selection;
      MFEM_ASSERT(!selection.response.models.empty() && !selection.ordered_edges.empty(),
                  "Missing parallel-cluster response-model selection!");
      const std::size_t first_index = selection.ordered_edges.front();
      const auto &first = segments[first_index];
      const double orientation = Dot(first.tangent, span.tangent);
      MFEM_ASSERT(std::abs(orientation) > 1.0 - 1.0e-8,
                  "A parallel-edge cluster contains incompatible tangents!");
      double begin = (span.begin - Dot(first.p0, span.tangent)) / orientation;
      double end = (span.end - Dot(first.p0, span.tangent)) / orientation;
      if (begin > end)
      {
        std::swap(begin, end);
      }
      std::vector<std::pair<double, double>> intervals = {
          {std::clamp(begin, 0.0, first.length), std::clamp(end, 0.0, first.length)}};
      for (const std::size_t edge_index : selection.ordered_edges)
      {
        const auto &edge = segments[edge_index];
        for (const auto &[excluded_begin, excluded_end] :
             vertex_excluded_intervals[edge_index])
        {
          const Point3D excluded_p0 = Interpolate(edge, excluded_begin);
          const Point3D excluded_p1 = Interpolate(edge, excluded_end);
          const double first_begin = Dot(Subtract(excluded_p0, first.p0), first.tangent);
          const double first_end = Dot(Subtract(excluded_p1, first.p0), first.tangent);
          SubtractInterval(intervals, std::min(first_begin, first_end),
                           std::max(first_begin, first_end));
        }
      }
      for (const auto &[interval_begin, interval_end] : intervals)
      {
        if (interval_end <= interval_begin)
        {
          continue;
        }
        if (requirements)
        {
          const double coordinate =
              0.5 * (interval_begin + interval_end) * Dot(first.tangent, span.tangent) +
              Dot(first.p0, span.tangent);
          requirements->Add(3, LibraryTopology::PARALLEL_EDGE_CLUSTER, group.targets,
                            first.boundary_condition,
                            DescribeParallelCluster(selection.ordered_edges, span.tangent,
                                                    coordinate, selection.axis_u),
                            library, &selection.response, interval_end - interval_begin);
          continue;
        }
        for (int q = 0; q < quadrature.GetNPoints(); q++)
        {
          const auto &ip = quadrature.IntPoint(q);
          const double first_distance =
              interval_begin + (interval_end - interval_begin) * ip.x;
          const Point3D first_point = Interpolate(first, first_distance);
          const double longitudinal_coordinate = Dot(first_point, span.tangent);

          ResponsePatchData patch;
          patch.origin = first_point;
          patch.axis_u = selection.axis_u;
          patch.axis_v = selection.axis_v;
          patch.conductor_references = selection.response.conductor_references;
          patch.maxwell_reference_is_pec =
              std::all_of(selection.ordered_edges.begin(), selection.ordered_edges.end(),
                          [&](std::size_t index)
                          {
                            return segments[index].boundary_condition.type ==
                                   MetalBoundaryConditionType::PEC;
                          });
          for (const std::size_t reference_edge : selection.reference_edges)
          {
            patch.maxwell_conductor_anchors.push_back(InterpolateAtLongitudinalCoordinate(
                segments[reference_edge], span.tangent, longitudinal_coordinate));
          }
          const double quadrature_weight = (interval_end - interval_begin) * ip.weight;
          for (const auto &weighted_model : selection.response.models)
          {
            const auto &source = library.models[weighted_model.index];
            MFEM_VERIFY(
                source.coupon_depth > 0.0,
                "Three-dimensional response correction requires CouponDepth for every "
                "selected fabrication-process response model!");
            auto weighted_patch = patch;
            weighted_patch.weight =
                weighted_model.weight * quadrature_weight / source.coupon_depth;
            pending.push_back({weighted_model.index, std::move(weighted_patch)});
          }
        }
      }
    };

    for (const auto &span : parallel_cluster_spans)
    {
      if (!group_matched)
      {
        break;
      }
      AppendParallelClusterQuadrature(span);
      group_matched_intervals++;
    }

    for (const auto &pair : pairs)
    {
      if (!group_matched)
      {
        break;
      }
      const auto &first = segments[pair.first];
      const auto &second = segments[pair.second];
      std::vector<std::pair<double, double>> pair_intervals = {
          {pair.first_begin, pair.first_end}};
      for (const auto &span : parallel_cluster_spans)
      {
        const auto &cluster = span.selection.ordered_edges;
        if (std::find(cluster.begin(), cluster.end(), pair.first) == cluster.end() ||
            std::find(cluster.begin(), cluster.end(), pair.second) == cluster.end())
        {
          continue;
        }
        const double orientation = Dot(first.tangent, span.tangent);
        MFEM_ASSERT(std::abs(orientation) > 1.0 - 1.0e-8,
                    "A parallel-edge cluster contains incompatible tangents!");
        double excluded_begin = (span.begin - Dot(first.p0, span.tangent)) / orientation;
        double excluded_end = (span.end - Dot(first.p0, span.tangent)) / orientation;
        if (excluded_begin > excluded_end)
        {
          std::swap(excluded_begin, excluded_end);
        }
        SubtractInterval(pair_intervals, excluded_begin, excluded_end);
      }
      for (const auto &[pair_begin, pair_end] : pair_intervals)
      {
        if (pair_end <= pair_begin)
        {
          continue;
        }
        const Point3D first_mid = Interpolate(first, 0.5 * (pair_begin + pair_end));
        double second_distance = Dot(Subtract(first_mid, second.p0), second.tangent);
        second_distance = std::clamp(second_distance, 0.0, second.length);
        const Point3D second_mid = Interpolate(second, second_distance);
        const Point3D direction = Normalize(Subtract(second_mid, first_mid));
        if (!SameBoundaryLaw(first.boundary_condition, second.boundary_condition))
        {
          Mpi::Warning(
              "Nearby three-dimensional edges use distinct metal boundary conditions; "
              "a dedicated spatial coupon is required and correction is disabled for "
              "this interface group!\n");
          group_matched = false;
          break;
        }
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
        else if (outward)
        {
          // Local outward-facing edges bound a physical strip even when its two perimeter
          // loops have different graph-component labels.
          topology = LibraryTopology::SAME_CONDUCTOR_STRIP;
        }
        else
        {
          Mpi::Warning(
              "No canonical paired-edge topology for nearby three-dimensional metal "
              "edges; correction is disabled for this interface group!\n");
          group_matched = false;
          break;
        }
        const double separation = Distance(first_mid, second_mid);
        const auto model_selection =
            FindLibraryModel(library, topology, separation, first.boundary_condition);
        if (!model_selection)
        {
          Mpi::Warning(
              "Fabrication-process response library \"{}\" has no {} model at separation "
              "{:.6e} mesh units; correction is disabled for this interface group!\n",
              library.name, TopologyName(topology), separation * coordinate_scale);
          if (requirements)
          {
            requirements->Add(
                3, topology, group.targets, first.boundary_condition,
                {{"EdgeCount", 2}, {"Separation", requirements->ScaleLength(separation)}},
                library, nullptr, pair_end - pair_begin,
                "No compatible paired-edge model or interpolation bracket");
          }
          group_matched = false;
          break;
        }
        group_maximum_library_distance =
            std::max(group_maximum_library_distance, model_selection->normalized_distance);
        AppendQuadrature(*model_selection, pair.first, pair_begin, pair_end, pair.second);
        group_matched_intervals++;
        group_interpolated_paired_intervals += model_selection->IsInterpolated();
      }
    }

    for (std::size_t i = 0; group_matched && i < segments.size(); i++)
    {
      std::vector<std::pair<double, double>> isolated_intervals;
      double begin = 0.0;
      for (const auto &[paired_begin, paired_end] : paired_intervals[i])
      {
        if (paired_begin > begin)
        {
          isolated_intervals.emplace_back(begin, paired_begin);
        }
        begin = std::max(begin, paired_end);
      }
      if (begin < segments[i].length)
      {
        isolated_intervals.emplace_back(begin, segments[i].length);
      }
      if (isolated_intervals.empty())
      {
        continue;
      }

      const auto isolated_model = FindLibraryModel(library, LibraryTopology::ISOLATED_EDGE,
                                                   0.0, segments[i].boundary_condition);
      if (!isolated_model)
      {
        Mpi::Warning(
            "Fabrication-process response library \"{}\" has no isolated-edge model for "
            "an unmatched longitudinal span using the selected metal boundary condition; "
            "correction is disabled for this interface group!\n",
            library.name);
        if (requirements)
        {
          for (const auto &[isolated_begin, isolated_end] : isolated_intervals)
          {
            requirements->Add(
                3, LibraryTopology::ISOLATED_EDGE, group.targets,
                segments[i].boundary_condition, {{"EdgeCount", 1}}, library, nullptr,
                isolated_end - isolated_begin,
                "No compatible isolated-edge model for this metal boundary condition");
          }
        }
        group_matched = false;
        break;
      }
      for (const auto &[isolated_begin, isolated_end] : isolated_intervals)
      {
        AppendQuadrature(*isolated_model, i, isolated_begin, isolated_end, std::nullopt);
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

    if (diagnostics)
    {
      for (const auto &selection : pending)
      {
        diagnostics->boundary_law_verified &=
            IsBoundaryLawVerified(library.models[selection.library_model]);
      }
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
        model.name = source.name;
        model.topology = TopologyName(source.topology);
        MapLibraryInterfaces(source, {{0, group.targets}}, model);
        result.models.push_back(std::move(model));
      }
      selection.patch.model = model_it->second;
      result.patches.push_back(selection.patch);
    }
    matched_intervals += group_matched_intervals;
    interpolated_paired_intervals += group_interpolated_paired_intervals;
    matched_segments += static_cast<int>(segments.size());
    matched_corner_patches += matched_corners;
    matched_endpoint_patches += matched_endpoints;
    matched_junction_patches += matched_junctions;
    matched_spatial_cluster_patches += matched_spatial_clusters;
    const int spatial_nonregular_vertices = static_cast<int>(std::count_if(
        spatially_excluded_vertices.begin(), spatially_excluded_vertices.end(),
        [&](std::size_t vertex)
        {
          return geometry.vertices[vertex].physical_type &&
                 *geometry.vertices[vertex].physical_type != MetalEdgeVertexType::REGULAR &&
                 !geometry.vertices[vertex].on_truncation_boundary;
        }));
    matched_nonregular_vertices += matched_sharp_corners + matched_endpoints +
                                   matched_junctions + spatial_nonregular_vertices;
    if (diagnostics)
    {
      const double unmatched_corner_neighborhood_length = std::max(
          0.0, group_corner_neighborhood_length - group_modeled_corner_neighborhood_length);
      diagnostics->matched_length += group_selected_length;
      diagnostics->matched_corner_neighborhood_length +=
          unmatched_corner_neighborhood_length;
      for (const auto &[type, target] : group.targets)
      {
        (void)type;
        diagnostics->matched_length_by_interface[target] += group_selected_length;
        diagnostics->matched_corner_neighborhood_length_by_interface[target] +=
            unmatched_corner_neighborhood_length;
      }
      diagnostics->maximum_curvature_ratio =
          std::max(diagnostics->maximum_curvature_ratio, group_maximum_curvature_ratio);
      diagnostics->maximum_library_distance =
          std::max(diagnostics->maximum_library_distance, group_maximum_library_distance);
    }
  }

  MFEM_VERIFY(requirements || (!result.models.empty() && !result.patches.empty()),
              "Fabrication-process response matching produced no usable correction "
              "patches!");
  Mpi::Print("\nAutomatic fabrication-process response matching:\n"
             " Library: {}\n"
             " Matched physical edge segments: {:d}\n"
             " Matched longitudinal intervals: {:d}\n"
             " Longitudinal quadrature patches: {:d}\n"
             " Matched corner patches: {:d}\n"
             " Matched endpoint patches: {:d}\n"
             " Matched junction patches: {:d}\n"
             " Matched spatial edge-cluster patches: {:d}\n"
             " Interpolated paired intervals: {:d}\n"
             " Interpolated rounded corners: {:d}\n"
             " Unmatched interface groups: {:d}\n",
             library.name, matched_segments, matched_intervals,
             static_cast<int>(result.patches.size()), matched_corner_patches,
             matched_endpoint_patches, matched_junction_patches,
             matched_spatial_cluster_patches, interpolated_paired_intervals,
             interpolated_rounded_corners, unmatched_groups);
  if (unmatched_rounded_corners > 0)
  {
    Mpi::Warning(
        "The selected three-dimensional metal perimeter has {} rounded corner "
        "neighborhoods with radius smaller than R but no compatible fabrication-process "
        "library model. Straight-edge coupon response is integrated through these "
        "neighborhoods; add a radius-aware corner model when their contribution is "
        "significant.\n",
        unmatched_rounded_corners);
  }
  const int unmatched_vertices = nonregular_vertices - matched_nonregular_vertices;
  if (unmatched_vertices > 0)
  {
    Mpi::Warning(
        "The selected three-dimensional metal perimeter has {} unmatched corner, endpoint, "
        "or junction vertices. Straight-edge coupon response is integrated through these "
        "neighborhoods; validate or add a matching spatial-vertex response model when "
        "their contribution is significant.\n",
        unmatched_vertices);
  }
  return result;
}

ResponseCorrectionData BuildAutomaticResponseData(const IoData &iodata,
                                                  const LaplaceOperator &laplace_op,
                                                  const ResponseCorrectionData &request,
                                                  AutomaticResponseStatistics *statistics)
{
  const auto &mesh = laplace_op.GetH1Space().GetParMesh();
  if (mesh.Dimension() == 2 && mesh.SpaceDimension() == 2)
  {
    return BuildAutomaticResponseData2D(iodata, mesh, laplace_op.GetMaterialOp(), request,
                                        false, nullptr, nullptr, statistics);
  }
  if (mesh.Dimension() == 3 && mesh.SpaceDimension() == 3)
  {
    return BuildAutomaticResponseData3D(iodata, mesh, laplace_op.GetMaterialOp(), request,
                                        false, nullptr, nullptr, statistics);
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

DomainResponseMatrices
BuildDomainResponseMatrices(const std::string &fabricated_path,
                            const std::string &thin_path, int expected_size,
                            const std::vector<int> &zero_trace_indices, const Units &units)
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
  mfem::DenseMatrix fixed_flux_transform(expected_size);
  fixed_flux_transform = 0.0;
  if (zero_trace_indices.empty())
  {
    mfem::DenseMatrixInverse(fabricated, true).Mult(thin, fixed_flux_transform);
  }
  else
  {
    std::vector<bool> constrained(expected_size, false);
    for (const int index : zero_trace_indices)
    {
      MFEM_VERIFY(index >= 0 && index < expected_size && !constrained[index],
                  "Invalid or duplicate zero-trace response basis index!");
      constrained[index] = true;
    }
    std::vector<int> free_indices;
    free_indices.reserve(expected_size - zero_trace_indices.size());
    for (int i = 0; i < expected_size; i++)
    {
      if (!constrained[i])
      {
        free_indices.push_back(i);
      }
    }
    MFEM_VERIFY(!free_indices.empty(),
                "Fixed-flux response requires at least one free trace basis function!");

    const int free_size = static_cast<int>(free_indices.size());
    mfem::DenseMatrix fabricated_free(free_size), thin_free(free_size);
    for (int i = 0; i < free_size; i++)
    {
      for (int j = 0; j < free_size; j++)
      {
        fabricated_free(i, j) = fabricated(free_indices[i], free_indices[j]);
        thin_free(i, j) = thin(free_indices[i], free_indices[j]);
      }
    }
    mfem::DenseMatrix fixed_flux_free;
    mfem::DenseMatrixInverse(fabricated_free, true).Mult(thin_free, fixed_flux_free);
    for (int i = 0; i < free_size; i++)
    {
      for (int j = 0; j < free_size; j++)
      {
        fixed_flux_transform(free_indices[i], free_indices[j]) = fixed_flux_free(i, j);
      }
    }
  }
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
    mfem::DenseMatrix fabricated_contribution(fabricated_it->second);
    fabricated_contribution *= scale;
    auto [fabricated_target, fabricated_inserted] =
        result.fabricated.emplace(mapping.target, fabricated_contribution);
    if (!fabricated_inserted)
    {
      fabricated_target->second += fabricated_contribution;
    }

    mfem::DenseMatrix defect_contribution(fabricated_it->second);
    defect_contribution.Add(-1.0, thin_it->second);
    defect_contribution *= scale;
    auto [defect, defect_inserted] =
        result.defects.emplace(mapping.target, defect_contribution);
    if (!defect_inserted)
    {
      defect->second += defect_contribution;
    }
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

void WriteSurfaceResponseRequirements(const IoData &iodata, const Mesh &mesh,
                                      const std::string &path)
{
  const auto &configured_request = iodata.problem.type == ProblemType::ELECTROSTATIC
                                       ? iodata.solver.electrostatic.response_correction
                                       : iodata.solver.surface_response_correction;
  MFEM_VERIFY(configured_request && configured_request->IsAutomatic(),
              "Surface-response preflight requires an automatic fabrication-process "
              "response library!");

  auto request = *configured_request;
  request.unmatched_policy = ResponseCorrectionData::UnmatchedPolicy::WARN;
  MaterialOperator mat_op(iodata, mesh);
  AutomaticResponseRequirements requirements(iodata.units,
                                             iodata.InputsNondimensionalized());
  AutomaticResponseStatistics statistics;
  const auto &parallel_mesh = mesh.Get();
  const bool maxwell = iodata.problem.type != ProblemType::ELECTROSTATIC;
  if (parallel_mesh.Dimension() == 2 && parallel_mesh.SpaceDimension() == 2)
  {
    BuildAutomaticResponseData2D(iodata, parallel_mesh, mat_op, request, maxwell, nullptr,
                                 &requirements, &statistics);
  }
  else if (parallel_mesh.Dimension() == 3 && parallel_mesh.SpaceDimension() == 3)
  {
    BuildAutomaticResponseData3D(iodata, parallel_mesh, mat_op, request, maxwell, nullptr,
                                 &requirements, &statistics);
  }
  else
  {
    MFEM_ABORT("Surface-response preflight requires a 2D or 3D solve mesh!");
  }

  requirements.SetStatistics(BuildAutomaticStatistics(parallel_mesh.GetComm(), statistics));
  if (Mpi::Root(parallel_mesh.GetComm()))
  {
    auto manifest = requirements.Build();
    manifest["MeshDimension"] = parallel_mesh.Dimension();
    manifest["Maxwell"] = maxwell;
    std::ofstream output(path);
    MFEM_VERIFY(output, "Unable to open surface-response requirements manifest \""
                            << path << "\"!");
    output << manifest.dump(2) << '\n';
  }
  Mpi::Barrier(parallel_mesh.GetComm());
  Mpi::Print(parallel_mesh.GetComm(),
             "\nSurface-response process-library preflight complete:\n Manifest: {}\n",
             path);
}

struct SurfaceResponseGeometry::Impl
{
  ResponseCorrectionData config;
  std::optional<AutomaticResponseDiagnostics> diagnostics;
  nlohmann::json statistics;
  bool maxwell = false;
  int dimension = 0;
};

SurfaceResponseOperator::SurfaceResponseOperator(
    const IoData &iodata, const LaplaceOperator &laplace_op,
    std::shared_ptr<const SurfaceResponseGeometry> *automatic_geometry)
  : Operator(laplace_op.GetH1Space().GetTrueVSize()), fespace(laplace_op.GetH1Space()),
    basis_size(0)
{
  BlockTimer setup_timer(Timer::CONSTRUCT_RESPONSE);
  const auto &request = iodata.solver.electrostatic.response_correction;
  MFEM_VERIFY(request, "Missing electrostatic surface response correction configuration!");
  const int dimension = fespace.Dimension();
  MFEM_VERIFY((dimension == 2 || dimension == 3) && fespace.SpaceDimension() == dimension,
              "Surface response correction requires a 2D or 3D electrostatic mesh!");
  MFEM_VERIFY(dimension == 2 || request->IsAutomatic(),
              "Three-dimensional surface response correction requires automatic "
              "fabrication-process library matching!");
  std::optional<ResponseCorrectionData> automatic_config;
  std::shared_ptr<const SurfaceResponseGeometry> cached_geometry;
  const ResponseCorrectionData *config = &*request;
  if (request->IsAutomatic())
  {
    if (automatic_geometry && *automatic_geometry)
    {
      cached_geometry = *automatic_geometry;
      MFEM_VERIFY(!cached_geometry->impl->maxwell &&
                      cached_geometry->impl->dimension == dimension,
                  "Cannot reuse Maxwell surface-response geometry for electrostatics!");
      config = &cached_geometry->impl->config;
      automatic_statistics = cached_geometry->impl->statistics;
    }
    else
    {
      BlockTimer geometry_timer(Timer::CONSTRUCT_RESPONSE_GEOMETRY);
      AutomaticResponseStatistics statistics;
      automatic_config =
          BuildAutomaticResponseData(iodata, laplace_op, *request, &statistics);
      automatic_statistics = BuildAutomaticStatistics(fespace.GetComm(), statistics);
      if (automatic_geometry)
      {
        auto impl = std::make_shared<SurfaceResponseGeometry::Impl>();
        impl->config = std::move(*automatic_config);
        impl->statistics = automatic_statistics;
        impl->dimension = dimension;
        cached_geometry = std::shared_ptr<const SurfaceResponseGeometry>(
            new SurfaceResponseGeometry(std::move(impl)));
        *automatic_geometry = cached_geometry;
        config = &cached_geometry->impl->config;
      }
      else
      {
        config = &*automatic_config;
      }
    }
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
    model.name =
        model_config.name.empty() ? fmt::format("model-{}", model.idx) : model_config.name;
    model.topology = model_config.topology.empty() ? "Explicit" : model_config.topology;
    model.contour_size = static_cast<int>(points.size());
    model.conductor_state_count = model_config.conductor_state_count;
    MFEM_VERIFY(model.conductor_state_count >= 0,
                "Response correction requires a nonnegative conductor-state count!");
    model.basis_size = model.contour_size + model.conductor_state_count;
    model.spatial_basis = model_config.spatial_basis;
    model.contour_groups = model_config.contour_groups;
    model.zero_trace_indices = model_config.zero_trace_indices;
    for (const auto &path : model_config.open_contour_paths)
    {
      model.open_contour_paths.push_back(
          {path.indices, path.start_conductor, path.end_conductor});
    }
    MFEM_VERIFY(model.contour_groups.empty() || model.open_contour_paths.empty(),
                "Response-correction models cannot combine closed ContourGroups with "
                "OpenContourPaths!");
    MFEM_VERIFY(model.zero_trace_indices.empty() || model.open_contour_paths.empty(),
                "Response-correction models cannot combine ZeroTraceIndices with "
                "OpenContourPaths!");
    if (model.contour_groups.empty() && model.open_contour_paths.empty())
    {
      model.contour_groups.push_back(model.contour_size);
    }
    if (!model.contour_groups.empty())
    {
      MFEM_VERIFY(std::accumulate(model.contour_groups.begin(), model.contour_groups.end(),
                                  0) == model.contour_size,
                  "Response-correction ContourGroups do not partition BasisPoints!");
      MFEM_VERIFY(std::all_of(model.zero_trace_indices.begin(),
                              model.zero_trace_indices.end(), [&](int index)
                              { return index >= 0 && index < model.contour_size; }),
                  "Response-correction ZeroTraceIndices contain an invalid BasisPoints "
                  "index!");
    }
    else
    {
      std::vector<bool> assigned(model.contour_size, false);
      for (const auto &path : model.open_contour_paths)
      {
        for (const int index : path.indices)
        {
          MFEM_VERIFY(index >= 0 && index < model.contour_size && !assigned[index],
                      "Response-correction OpenContourPaths contain an invalid or "
                      "duplicate BasisPoints index!");
          assigned[index] = true;
        }
      }
      MFEM_VERIFY(
          std::all_of(assigned.begin(), assigned.end(), [](bool value) { return value; }),
          "Response-correction OpenContourPaths do not partition BasisPoints!");
    }
    auto domain_response = BuildDomainResponseMatrices(
        model_config.fabricated_matrix, model_config.thin_matrix, model.basis_size,
        model.zero_trace_indices, iodata.units);
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

  const int rank = Mpi::Rank(fespace.GetComm());
  const int size = Mpi::Size(fespace.GetComm());
  int point_count = 0;
  std::vector<std::size_t> local_patch_indices;
  for (std::size_t patch_idx = 0; patch_idx < config->patches.size(); patch_idx++)
  {
    const auto &patch_config = config->patches[patch_idx];
    const auto model_it = model_indices.find(patch_config.model);
    MFEM_VERIFY(model_it != model_indices.end(),
                "Response-correction patch refers to an unknown model index!");
    const auto &model = models[model_it->second];
    MFEM_VERIFY(std::isfinite(patch_config.weight) && patch_config.weight > 0.0,
                "Response-correction patch weights must be positive!");
    MFEM_VERIFY(static_cast<int>(patch_config.conductor_references.size()) ==
                    model.conductor_state_count + 1,
                "Response-correction patch conductor references do not match its model!");
    if (static_cast<int>(patch_idx % size) != rank)
    {
      continue;
    }
    local_patch_indices.push_back(patch_idx);
    patches.push_back(
        Patch{model_it->second, point_count, basis_size, patch_config.weight});
    point_count += model.contour_size + 1 + model.conductor_state_count;
    basis_size += model.basis_size;
  }
  global_patch_count = static_cast<int>(patches.size());
  global_basis_size = basis_size;
  Mpi::GlobalSum(1, &global_patch_count, fespace.GetComm());
  Mpi::GlobalSum(1, &global_basis_size, fespace.GetComm());

  mfem::Vector xyz(dimension * point_count);
  const double coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
  std::vector<std::vector<Point2D>> polygons;
  if (dimension == 2)
  {
    polygons.reserve(config->patches.size());
    for (const auto &patch_config : config->patches)
    {
      const auto model_it = model_indices.find(patch_config.model);
      MFEM_ASSERT(model_it != model_indices.end(), "Unknown response model!");
      const auto &local_points = basis_points[model_it->second];
      auto &polygon = polygons.emplace_back();
      polygon.reserve(local_points.size());
      for (const auto &local : local_points)
      {
        polygon.push_back({patch_config.origin[0] + (local[0] * patch_config.axis_u[0] +
                                                     local[1] * patch_config.axis_v[0]) /
                                                        coordinate_scale,
                           patch_config.origin[1] + (local[0] * patch_config.axis_u[1] +
                                                     local[1] * patch_config.axis_v[1]) /
                                                        coordinate_scale});
      }
    }
  }
  int point = 0;
  for (std::size_t patch_idx = 0; patch_idx < patches.size(); patch_idx++)
  {
    const auto &patch = patches[patch_idx];
    const auto &patch_config = config->patches[local_patch_indices[patch_idx]];
    const auto &model = models[patch.model];
    const auto &local_points = basis_points[patch.model];
    double norm_u = 0.0, norm_v = 0.0, norm_w = 0.0;
    double dot_uv = 0.0, dot_uw = 0.0, dot_vw = 0.0;
    for (int d = 0; d < dimension; d++)
    {
      norm_u += patch_config.axis_u[d] * patch_config.axis_u[d];
      norm_v += patch_config.axis_v[d] * patch_config.axis_v[d];
      norm_w += patch_config.axis_w[d] * patch_config.axis_w[d];
      dot_uv += patch_config.axis_u[d] * patch_config.axis_v[d];
      dot_uw += patch_config.axis_u[d] * patch_config.axis_w[d];
      dot_vw += patch_config.axis_v[d] * patch_config.axis_w[d];
    }
    norm_u = std::sqrt(norm_u);
    norm_v = std::sqrt(norm_v);
    MFEM_VERIFY(std::abs(norm_u - 1.0) < 1.0e-10 && std::abs(norm_v - 1.0) < 1.0e-10 &&
                    std::abs(dot_uv) < 1.0e-10,
                "Surface response correction AxisU and AxisV must be orthonormal in the "
                "coupon cross-section!");
    if (model.spatial_basis)
    {
      MFEM_VERIFY(dimension == 3 && std::abs(std::sqrt(norm_w) - 1.0) < 1.0e-10 &&
                      std::abs(dot_uw) < 1.0e-10 && std::abs(dot_vw) < 1.0e-10,
                  "A spatial response-correction basis requires an orthonormal three-"
                  "dimensional coupon frame!");
    }

    for (const auto &local : local_points)
    {
      MFEM_VERIFY(model.spatial_basis || std::abs(local[2]) <= 1.0e-12,
                  "Response-correction basis points must lie in the local coupon "
                  "cross-section!");
      for (int d = 0; d < dimension; d++)
      {
        double coordinate = patch_config.origin[d] + (local[0] * patch_config.axis_u[d] +
                                                      local[1] * patch_config.axis_v[d]) /
                                                         coordinate_scale;
        if (model.spatial_basis)
        {
          coordinate += local[2] * patch_config.axis_w[d] / coordinate_scale;
        }
        xyz(d * point_count + point) = coordinate;
      }
      point++;
    }
    for (const auto &reference : patch_config.conductor_references)
    {
      for (int d = 0; d < dimension; d++)
      {
        xyz(d * point_count + point) =
            patch_config.origin[d] + reference[0] * patch_config.axis_u[d] +
            reference[1] * patch_config.axis_v[d] + reference[2] * patch_config.axis_w[d];
      }
      point++;
    }
  }
  for (std::size_t i = 0; i < polygons.size(); i++)
  {
    for (std::size_t j = i + 1; j < polygons.size(); j++)
    {
      const int interpolation_group = config->patches[i].interpolation_group;
      if (interpolation_group > 0 &&
          interpolation_group == config->patches[j].interpolation_group)
      {
        continue;
      }
      MFEM_VERIFY(!PolygonsOverlap(polygons[i], polygons[j]),
                  "Response-correction patches " << i + 1 << " and " << j + 1
                                                 << " overlap. Replace nearby one-edge "
                                                    "patches with one coupled multi-edge "
                                                    "coupon model!");
    }
  }

  {
    BlockTimer point_timer(Timer::CONSTRUCT_RESPONSE_POINTS);
    ConfigurePointCommunication(xyz, dimension);
  }

  Mpi::Print("\nConfigured surface response correction:\n"
             " Coupon models: {:d}\n"
             " Global patches: {:d}\n"
             " Total trace coefficients: {:d}\n",
             static_cast<int>(models.size()), global_patch_count, global_basis_size);
#else
  MFEM_ABORT("Surface response correction requires MFEM_USE_GSLIB!");
#endif
}

SurfaceResponseOperator::SurfaceResponseOperator(
    const IoData &iodata, const SpaceOperator &space_op,
    std::shared_ptr<const SurfaceResponseGeometry> *automatic_geometry)
  : Operator(space_op.GetNDSpace().GetTrueVSize()), fespace(space_op.GetNDSpace()),
    basis_size(0)
{
  ConfigureMaxwellResponse(iodata, space_op.GetMaterialOp(),
                           space_op.GetNDDbcTDofLists().back(), automatic_geometry);
}

SurfaceResponseOperator::SurfaceResponseOperator(
    const IoData &iodata, const BoundaryModeOperator &mode_op,
    std::shared_ptr<const SurfaceResponseGeometry> *automatic_geometry)
  : Operator(mode_op.GetNDSpace().GetTrueVSize()), fespace(mode_op.GetNDSpace()),
    basis_size(0)
{
  ConfigureMaxwellResponse(iodata, mode_op.GetMaterialOp(),
                           mode_op.GetNDDbcTDofLists().back(), automatic_geometry);
}

void SurfaceResponseOperator::ConfigureMaxwellResponse(
    const IoData &iodata, const MaterialOperator &mat_op,
    const mfem::Array<int> &essential_tdofs,
    std::shared_ptr<const SurfaceResponseGeometry> *automatic_geometry)
{
  BlockTimer setup_timer(Timer::CONSTRUCT_RESPONSE);
  const auto &request = iodata.solver.surface_response_correction;
  MFEM_VERIFY(request, "Missing Maxwell surface response correction configuration!");
  MFEM_VERIFY(request->IsAutomatic(),
              "Maxwell surface response correction requires automatic fabrication-"
              "process library matching!");
  const int dimension = fespace.Dimension();
  MFEM_VERIFY((dimension == 2 || dimension == 3) && fespace.SpaceDimension() == dimension,
              "Maxwell surface response correction requires a two- or three-dimensional "
              "mesh!");

#if defined(MFEM_USE_GSLIB)
  maxwell = true;
  AutomaticResponseDiagnostics local_diagnostics;
  AutomaticResponseStatistics local_statistics;
  std::optional<ResponseCorrectionData> local_config;
  std::shared_ptr<const SurfaceResponseGeometry> cached_geometry;
  const ResponseCorrectionData *config_ptr = nullptr;
  const AutomaticResponseDiagnostics *diagnostics_ptr = nullptr;
  if (automatic_geometry && *automatic_geometry)
  {
    cached_geometry = *automatic_geometry;
    MFEM_VERIFY(cached_geometry->impl->maxwell &&
                    cached_geometry->impl->diagnostics.has_value() &&
                    cached_geometry->impl->dimension == dimension,
                "Cannot reuse incompatible surface-response geometry for Maxwell!");
    config_ptr = &cached_geometry->impl->config;
    diagnostics_ptr = &*cached_geometry->impl->diagnostics;
    automatic_statistics = cached_geometry->impl->statistics;
  }
  else
  {
    BlockTimer geometry_timer(Timer::CONSTRUCT_RESPONSE_GEOMETRY);
    if (dimension == 2)
    {
      local_config =
          BuildAutomaticResponseData2D(iodata, fespace.GetParMesh(), mat_op, *request, true,
                                       &local_diagnostics, nullptr, &local_statistics);
    }
    else
    {
      local_config =
          BuildAutomaticResponseData3D(iodata, fespace.GetParMesh(), mat_op, *request, true,
                                       &local_diagnostics, nullptr, &local_statistics);
    }
    automatic_statistics = BuildAutomaticStatistics(fespace.GetComm(), local_statistics);
    if (automatic_geometry)
    {
      auto impl = std::make_shared<SurfaceResponseGeometry::Impl>();
      impl->config = std::move(*local_config);
      impl->diagnostics = local_diagnostics;
      impl->statistics = automatic_statistics;
      impl->maxwell = true;
      impl->dimension = dimension;
      cached_geometry = std::shared_ptr<const SurfaceResponseGeometry>(
          new SurfaceResponseGeometry(std::move(impl)));
      *automatic_geometry = cached_geometry;
      config_ptr = &cached_geometry->impl->config;
      diagnostics_ptr = &*cached_geometry->impl->diagnostics;
    }
    else
    {
      config_ptr = &*local_config;
      diagnostics_ptr = &local_diagnostics;
    }
  }
  const auto &config = *config_ptr;
  const auto &diagnostics = *diagnostics_ptr;
  matching_radius = diagnostics.matching_radius;
  minimum_wave_speed = diagnostics.minimum_wave_speed;
  MFEM_VERIFY(!diagnostics.selected_length_by_interface.empty(),
              "Missing interface-resolved Maxwell surface-response diagnostics!");
  matched_length_fraction = 1.0;
  corner_neighborhood_fraction = 0.0;
  for (const auto &[interface, selected_length] : diagnostics.selected_length_by_interface)
  {
    const double matched_length =
        diagnostics.matched_length_by_interface.count(interface)
            ? diagnostics.matched_length_by_interface.at(interface)
            : 0.0;
    const double interface_matched_fraction =
        selected_length > 0.0 ? matched_length / selected_length : 0.0;
    matched_length_fraction = std::min(matched_length_fraction, interface_matched_fraction);
    matched_length_fraction_by_interface.emplace(interface, interface_matched_fraction);
    const double corner_length =
        diagnostics.matched_corner_neighborhood_length_by_interface.count(interface)
            ? diagnostics.matched_corner_neighborhood_length_by_interface.at(interface)
            : 0.0;
    const double interface_corner_fraction =
        matched_length > 0.0 ? corner_length / matched_length : 0.0;
    corner_neighborhood_fraction =
        std::max(corner_neighborhood_fraction, interface_corner_fraction);
    corner_neighborhood_fraction_by_interface.emplace(interface, interface_corner_fraction);
  }
  maximum_curvature_ratio = diagnostics.maximum_curvature_ratio;
  maximum_library_distance = diagnostics.maximum_library_distance;
  boundary_law_verified = diagnostics.boundary_law_verified;

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
    model.name =
        model_config.name.empty() ? fmt::format("model-{}", model.idx) : model_config.name;
    model.topology = model_config.topology.empty() ? "Explicit" : model_config.topology;
    model.contour_size = static_cast<int>(points.size());
    model.conductor_state_count = model_config.conductor_state_count;
    MFEM_VERIFY(model.conductor_state_count >= 0,
                "Maxwell response correction requires a nonnegative conductor-state "
                "count!");
    model.basis_size = model.contour_size + model.conductor_state_count;
    model.spatial_basis = model_config.spatial_basis;
    MFEM_VERIFY(dimension == 3 || !model.spatial_basis,
                "Two-dimensional BoundaryMode response correction requires planar "
                "coupon models!");
    model.contour_groups = model_config.contour_groups;
    model.zero_trace_indices = model_config.zero_trace_indices;
    for (const auto &path : model_config.open_contour_paths)
    {
      MFEM_VERIFY(path.start_conductor >= 0 &&
                      path.start_conductor <= model.conductor_state_count &&
                      path.end_conductor >= 0 &&
                      path.end_conductor <= model.conductor_state_count &&
                      path.start_conductor != path.end_conductor,
                  "Maxwell OpenContourPaths refer to invalid conductor indices!");
      model.open_contour_paths.push_back(
          {path.indices, path.start_conductor, path.end_conductor});
    }
    MFEM_VERIFY(model.conductor_state_count == 0 || !model.open_contour_paths.empty(),
                "Maxwell response correction requires OpenContourPaths for every "
                "two-conductor coupon model!");
    MFEM_VERIFY(model.contour_groups.empty() || model.open_contour_paths.empty(),
                "Response-correction models cannot combine closed ContourGroups with "
                "OpenContourPaths!");
    MFEM_VERIFY(model.zero_trace_indices.empty() || model.open_contour_paths.empty(),
                "Response-correction models cannot combine ZeroTraceIndices with "
                "OpenContourPaths!");
    if (model.contour_groups.empty() && model.open_contour_paths.empty())
    {
      model.contour_groups.push_back(model.contour_size);
    }
    if (!model.contour_groups.empty())
    {
      MFEM_VERIFY(std::accumulate(model.contour_groups.begin(), model.contour_groups.end(),
                                  0) == model.contour_size,
                  "Response-correction ContourGroups do not partition BasisPoints!");
      MFEM_VERIFY(std::all_of(model.zero_trace_indices.begin(),
                              model.zero_trace_indices.end(), [&](int index)
                              { return index >= 0 && index < model.contour_size; }),
                  "Response-correction ZeroTraceIndices contain an invalid BasisPoints "
                  "index!");
    }
    else
    {
      std::vector<bool> assigned(model.contour_size, false);
      for (const auto &path : model.open_contour_paths)
      {
        for (const int index : path.indices)
        {
          MFEM_VERIFY(index >= 0 && index < model.contour_size && !assigned[index],
                      "Response-correction OpenContourPaths contain an invalid or "
                      "duplicate BasisPoints index!");
          assigned[index] = true;
        }
      }
      MFEM_VERIFY(
          std::all_of(assigned.begin(), assigned.end(), [](bool value) { return value; }),
          "Response-correction OpenContourPaths do not partition BasisPoints!");
    }
    auto domain_response = BuildDomainResponseMatrices(
        model_config.fabricated_matrix, model_config.thin_matrix, model.basis_size,
        model.zero_trace_indices, iodata.units);
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

  const int rank = Mpi::Rank(fespace.GetComm());
  const int size = Mpi::Size(fespace.GetComm());
  global_patch_count = static_cast<int>(config.patches.size());
  for (const auto &patch_config : config.patches)
  {
    const auto model_it = model_indices.find(patch_config.model);
    MFEM_VERIFY(model_it != model_indices.end(),
                "Response-correction patch refers to an unknown model index!");
    global_basis_size += models[model_it->second].basis_size;
  }
  const std::size_t local_patch_capacity =
      (config.patches.size() + static_cast<std::size_t>(size) - 1) / size;
  maxwell_contours.reserve(local_patch_capacity);
  maxwell_conductor_anchors.reserve(local_patch_capacity);
  maxwell_paths.reserve(local_patch_capacity);
  const double coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
  for (std::size_t patch_index = 0; patch_index < config.patches.size(); patch_index++)
  {
    const auto &patch_config = config.patches[patch_index];
    const auto model_it = model_indices.find(patch_config.model);
    MFEM_VERIFY(model_it != model_indices.end(),
                "Response-correction patch refers to an unknown model index!");
    const auto &model = models[model_it->second];
    MFEM_VERIFY(std::isfinite(patch_config.weight) && patch_config.weight > 0.0,
                "Response-correction patch weights must be positive!");
    const std::size_t anchor_count = patch_config.maxwell_conductor_anchors.empty()
                                         ? patch_config.conductor_references.size()
                                         : patch_config.maxwell_conductor_anchors.size();
    MFEM_VERIFY(static_cast<int>(anchor_count) == model.conductor_state_count + 1,
                "Maxwell response-correction patch conductor anchors do not match its "
                "model!");
    MFEM_VERIFY(patch_config.maxwell_reference_is_pec || model.zero_trace_indices.empty(),
                "Finite-impedance Maxwell response patches cannot use PEC-constrained "
                "ZeroTraceIndices!");
    if (static_cast<int>(patch_index % size) != rank)
    {
      continue;
    }

    double norm_u = 0.0, norm_v = 0.0, norm_w = 0.0;
    double dot_uv = 0.0, dot_uw = 0.0, dot_vw = 0.0;
    for (int d = 0; d < 3; d++)
    {
      norm_u += patch_config.axis_u[d] * patch_config.axis_u[d];
      norm_v += patch_config.axis_v[d] * patch_config.axis_v[d];
      norm_w += patch_config.axis_w[d] * patch_config.axis_w[d];
      dot_uv += patch_config.axis_u[d] * patch_config.axis_v[d];
      dot_uw += patch_config.axis_u[d] * patch_config.axis_w[d];
      dot_vw += patch_config.axis_v[d] * patch_config.axis_w[d];
    }
    MFEM_VERIFY(std::abs(std::sqrt(norm_u) - 1.0) < 1.0e-10 &&
                    std::abs(std::sqrt(norm_v) - 1.0) < 1.0e-10 &&
                    std::abs(dot_uv) < 1.0e-10,
                "Surface response correction AxisU and AxisV must be orthonormal in the "
                "coupon cross-section!");
    if (model.spatial_basis)
    {
      MFEM_VERIFY(std::abs(std::sqrt(norm_w) - 1.0) < 1.0e-10 &&
                      std::abs(dot_uw) < 1.0e-10 && std::abs(dot_vw) < 1.0e-10,
                  "A spatial response-correction basis requires an orthonormal three-"
                  "dimensional coupon frame!");
    }

    auto &contour = maxwell_contours.emplace_back();
    contour.reserve(basis_points[model_it->second].size());
    for (const auto &local : basis_points[model_it->second])
    {
      MFEM_VERIFY(model.spatial_basis || std::abs(local[2]) <= 1.0e-12,
                  "Response-correction basis points must lie in the local coupon "
                  "cross-section!");
      mfem::Vector point(3);
      for (int d = 0; d < 3; d++)
      {
        point[d] = patch_config.origin[d] +
                   (local[0] * patch_config.axis_u[d] + local[1] * patch_config.axis_v[d]) /
                       coordinate_scale;
        if (model.spatial_basis)
        {
          point[d] += local[2] * patch_config.axis_w[d] / coordinate_scale;
        }
      }
      contour.push_back(std::move(point));
    }
    auto &anchors = maxwell_conductor_anchors.emplace_back();
    anchors.reserve(anchor_count);
    for (std::size_t i = 0; i < anchor_count; i++)
    {
      mfem::Vector anchor(3);
      anchor = 0.0;
      if (!patch_config.maxwell_conductor_anchors.empty())
      {
        std::copy(patch_config.maxwell_conductor_anchors[i].begin(),
                  patch_config.maxwell_conductor_anchors[i].end(), anchor.GetData());
      }
      else
      {
        for (int d = 0; d < dimension; d++)
        {
          anchor[d] = patch_config.origin[d] +
                      patch_config.conductor_references[i][0] * patch_config.axis_u[d] +
                      patch_config.conductor_references[i][1] * patch_config.axis_v[d];
        }
      }
      anchors.push_back(std::move(anchor));
    }

    patches.push_back(Patch{model_it->second, 0, basis_size, patch_config.weight});
    basis_size += model.basis_size;
  }

  maxwell_quadrature_order = 2 * std::max(1, iodata.solver.order) + 2;
  std::vector<MaxwellLineGeometry> line_geometry;
  auto AppendLine = [&](const mfem::Vector &p0, const mfem::Vector &p1)
  {
    MaxwellLineGeometry geometry;
    std::array<double, 3> tangent;
    for (int d = 0; d < dimension; d++)
    {
      geometry.begin[d] = p0[d];
      geometry.end[d] = p1[d];
      tangent[d] = geometry.end[d] - geometry.begin[d];
    }
    MFEM_VERIFY(Norm(tangent) > 0.0,
                "Maxwell response contour contains a zero-length line segment!");
    line_geometry.push_back(geometry);
    maxwell_lines.emplace_back();
    return static_cast<int>(maxwell_lines.size()) - 1;
  };

  for (std::size_t patch_index = 0; patch_index < patches.size(); patch_index++)
  {
    const auto &patch = patches[patch_index];
    const auto &model = models[patch.model];
    const auto &contour = maxwell_contours[patch_index];
    const auto &anchors = maxwell_conductor_anchors[patch_index];
    const auto &anchor = anchors.front();
    MFEM_ASSERT(static_cast<int>(contour.size()) == model.contour_size,
                "Inconsistent Maxwell response contour size!");
    MFEM_ASSERT(static_cast<int>(anchors.size()) == model.conductor_state_count + 1,
                "Inconsistent Maxwell conductor-anchor count!");

    const int path_begin = static_cast<int>(maxwell_paths.size());

    int group_offset = 0;
    for (const int group_size : model.contour_groups)
    {
      const auto zero_begin = std::lower_bound(
          model.zero_trace_indices.begin(), model.zero_trace_indices.end(), group_offset);
      const auto zero_end =
          std::lower_bound(model.zero_trace_indices.begin(), model.zero_trace_indices.end(),
                           group_offset + group_size);
      if (zero_begin != zero_end)
      {
        // Each PEC knot has an exact zero trace. Reconstruct every free arc between
        // consecutive knots independently, omitting segments whose two endpoints are
        // both PEC constrained.
        for (auto zero = zero_begin; zero != zero_end; zero++)
        {
          const int start = *zero;
          const int end = std::next(zero) != zero_end ? *std::next(zero) : *zero_begin;
          MaxwellContourPath path;
          path.closed = false;
          std::vector<int> contour_indices = {start};
          int index = group_offset + (start - group_offset + 1) % group_size;
          while (index != end)
          {
            contour_indices.push_back(index);
            index = group_offset + (index - group_offset + 1) % group_size;
          }
          if (contour_indices.size() == 1)
          {
            continue;
          }
          path.trace_indices.reserve(contour_indices.size());
          for (const int contour_index : contour_indices)
          {
            path.trace_indices.push_back(patch.trace_offset + contour_index);
          }
          path.contour_line_offset = static_cast<int>(maxwell_lines.size());
          path.contour_line_count = static_cast<int>(contour_indices.size()) - 1;
          for (std::size_t i = 1; i < contour_indices.size(); i++)
          {
            AppendLine(contour[contour_indices[i - 1]], contour[contour_indices[i]]);
          }
          path.end_line = AppendLine(contour[contour_indices.back()], contour[end]);
          maxwell_paths.push_back(std::move(path));
        }
        group_offset += group_size;
        continue;
      }

      MaxwellContourPath path;
      double start_distance = mfem::infinity();
      int start = 0;
      for (int i = 0; i < group_size; i++)
      {
        mfem::Vector delta(contour[group_offset + i]);
        delta -= anchor;
        const double distance = delta.Norml2();
        if (distance < start_distance)
        {
          start = i;
          start_distance = distance;
        }
      }
      path.trace_indices.reserve(group_size);
      for (int offset = 0; offset < group_size; offset++)
      {
        path.trace_indices.push_back(patch.trace_offset + group_offset +
                                     (start + offset) % group_size);
      }
      if (start_distance > 1.0e-14 * std::max(1.0, matching_radius))
      {
        path.anchor_line = AppendLine(anchor, contour[group_offset + start]);
      }
      path.contour_line_offset = static_cast<int>(maxwell_lines.size());
      path.contour_line_count = group_size;
      for (int offset = 0; offset < group_size; offset++)
      {
        const int i = group_offset + (start + offset) % group_size;
        const int next = group_offset + (start + offset + 1) % group_size;
        AppendLine(contour[i], contour[next]);
      }
      maxwell_paths.push_back(std::move(path));
      group_offset += group_size;
    }
    MFEM_ASSERT(group_offset == model.contour_size,
                "Inconsistent Maxwell contour-group partition!");
    std::vector<int> open_path_indices;
    open_path_indices.reserve(model.open_contour_paths.size());
    for (const auto &path_config : model.open_contour_paths)
    {
      const auto &start_anchor = anchors[path_config.start_conductor];
      const auto &end_anchor = anchors[path_config.end_conductor];
      MaxwellContourPath path;
      path.closed = false;
      if (path_config.start_conductor > 0)
      {
        path.start_conductor_trace =
            patch.trace_offset + model.contour_size + path_config.start_conductor - 1;
      }
      if (path_config.end_conductor > 0)
      {
        path.end_conductor_trace =
            patch.trace_offset + model.contour_size + path_config.end_conductor - 1;
      }
      path.trace_indices.reserve(path_config.indices.size());
      for (const int index : path_config.indices)
      {
        path.trace_indices.push_back(patch.trace_offset + index);
      }
      const auto &first = contour[path_config.indices.front()];
      const auto &last = contour[path_config.indices.back()];
      mfem::Vector delta(first);
      delta -= start_anchor;
      if (delta.Norml2() > 1.0e-14 * std::max(1.0, matching_radius))
      {
        path.anchor_line = AppendLine(start_anchor, first);
      }
      path.contour_line_offset = static_cast<int>(maxwell_lines.size());
      path.contour_line_count = static_cast<int>(path_config.indices.size()) - 1;
      for (std::size_t i = 1; i < path_config.indices.size(); i++)
      {
        AppendLine(contour[path_config.indices[i - 1]], contour[path_config.indices[i]]);
      }
      delta = last;
      delta -= end_anchor;
      if (delta.Norml2() > 1.0e-14 * std::max(1.0, matching_radius))
      {
        path.end_line = AppendLine(last, end_anchor);
      }
      open_path_indices.push_back(static_cast<int>(maxwell_paths.size()));
      maxwell_paths.push_back(std::move(path));
    }
    std::vector<bool> connected(anchors.size(), false);
    connected.front() = true;
    bool changed = true;
    while (changed)
    {
      changed = false;
      for (std::size_t i = 0; i < model.open_contour_paths.size(); i++)
      {
        const auto &path = model.open_contour_paths[i];
        if (connected[path.start_conductor] == connected[path.end_conductor])
        {
          continue;
        }
        const int parent =
            connected[path.start_conductor] ? path.start_conductor : path.end_conductor;
        const int conductor =
            connected[path.start_conductor] ? path.end_conductor : path.start_conductor;
        const int parent_trace_offset =
            parent > 0 ? patch.trace_offset + model.contour_size + parent - 1 : -1;
        const int trace_offset = patch.trace_offset + model.contour_size + conductor - 1;
        const int integral_sign = parent == path.start_conductor ? -1 : 1;
        maxwell_conductor_paths.push_back(
            {open_path_indices[i], parent_trace_offset, trace_offset, integral_sign});
        connected[conductor] = true;
        changed = true;
      }
    }
    MFEM_VERIFY(
        std::all_of(connected.begin(), connected.end(), [](bool value) { return value; }),
        "Maxwell OpenContourPaths must connect every conductor reference!");
    maxwell_patch_paths.emplace_back(path_begin,
                                     static_cast<int>(maxwell_paths.size()) - path_begin);
  }

  {
    BlockTimer point_timer(Timer::CONSTRUCT_RESPONSE_POINTS);
    ConfigureMaxwellLines(line_geometry);
  }
  dbc_tdof_list = essential_tdofs;
  contour_line_count = maxwell_lines.size();
  int global_line_count = static_cast<int>(maxwell_lines.size());
  Mpi::GlobalSum(1, &global_line_count, fespace.GetComm());

  Mpi::Print("\nConfigured Maxwell surface response correction:\n"
             " Coupon models: {:d}\n"
             " Response patches: {:d}\n"
             " Total contour coefficients: {:d}\n"
             " Contour line functionals: {:d}\n"
             " Minimum interface matched edge-length fraction: {:.6f}\n"
             " Maximum interface unmodeled corner-neighborhood fraction: {:.6f}\n"
             " Maximum R/rho: {:.3e}\n"
             " Maximum normalized library distance: {:.3e}\n",
             static_cast<int>(models.size()), global_patch_count, global_basis_size,
             global_line_count, matched_length_fraction, corner_neighborhood_fraction,
             maximum_curvature_ratio, maximum_library_distance);
  Mpi::Print(" Interface-resolved geometric coverage:\n");
  for (const auto &[interface, matched_fraction] : matched_length_fraction_by_interface)
  {
    Mpi::Print("  {:d}: matched = {:.6f}, unmodeled corner = {:.6f}\n", interface,
               matched_fraction, corner_neighborhood_fraction_by_interface.at(interface));
  }
#else
  MFEM_ABORT("Maxwell surface response correction requires MFEM_USE_GSLIB!");
#endif
}

void SurfaceResponseOperator::ConfigureMaxwellLines(
    const std::vector<MaxwellLineGeometry> &line_geometry)
{
  MFEM_ASSERT(line_geometry.size() == maxwell_lines.size(),
              "Inconsistent Maxwell line geometry!");
  auto &mesh = const_cast<mfem::ParMesh &>(fespace.GetParMesh());
  const auto comm = fespace.GetComm();
  const int size = Mpi::Size(comm);
  const int dimension = fespace.Dimension();
  ElementPointLocator locator(mesh, dimension);

  auto SetOffsets = [](const std::vector<int> &counts, std::vector<int> &offsets)
  {
    offsets.resize(counts.size());
    int total = 0;
    for (std::size_t i = 0; i < counts.size(); i++)
    {
      offsets[i] = total;
      total += counts[i];
    }
    return total;
  };
  auto ScaleCommunicationPlan = [](const std::vector<int> &values, int scale)
  {
    std::vector<int> result(values);
    for (auto &value : result)
    {
      value *= scale;
    }
    return result;
  };

  double coordinate_scale = 0.0;
  const auto &bounds = locator.GetBounds();
  for (int d = 0; d < dimension; d++)
  {
    coordinate_scale = std::max({coordinate_scale, std::abs(bounds.min[d]),
                                 std::abs(bounds.max[d]), bounds.max[d] - bounds.min[d]});
  }
  Mpi::GlobalMax(1, &coordinate_scale, comm);
  const double box_tolerance =
      1.0e-11 * coordinate_scale +
      64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, coordinate_scale);
  constexpr double parameter_tolerance = 1.0e-10;

  int exact_intersections = locator.SupportsExactSegmentIntersections() ? 1 : 0;
  Mpi::GlobalMin(1, &exact_intersections, comm);
  std::vector<std::vector<std::pair<double, double>>> line_intervals(line_geometry.size());
  if (exact_intersections)
  {
    constexpr int routing_box_count = 8;
    constexpr int routing_box_values = 6;
    std::array<double, routing_box_count * routing_box_values> local_routing;
    for (int box = 0; box < routing_box_count; box++)
    {
      for (int d = 0; d < 3; d++)
      {
        local_routing[routing_box_values * box + d] = mfem::infinity();
        local_routing[routing_box_values * box + 3 + d] = -mfem::infinity();
      }
    }
    const auto routing_boxes = locator.GetRoutingBoxes(routing_box_count);
    for (std::size_t box = 0; box < routing_boxes.size(); box++)
    {
      for (int d = 0; d < 3; d++)
      {
        local_routing[routing_box_values * box + d] = routing_boxes[box].min[d];
        local_routing[routing_box_values * box + 3 + d] = routing_boxes[box].max[d];
      }
    }
    std::vector<double> global_routing(size * local_routing.size());
    Mpi::Allgather(static_cast<int>(local_routing.size()), local_routing.data(),
                   global_routing.data(), comm);
    auto RankIntersects = [&](int candidate_rank, const MaxwellLineGeometry &line)
    {
      const int rank_offset = candidate_rank * routing_box_count * routing_box_values;
      for (int box = 0; box < routing_box_count; box++)
      {
        ElementBox candidate;
        for (int d = 0; d < 3; d++)
        {
          candidate.min[d] = global_routing[rank_offset + routing_box_values * box + d];
          candidate.max[d] = global_routing[rank_offset + routing_box_values * box + 3 + d];
        }
        if (candidate.min[0] <= candidate.max[0] &&
            candidate.IntersectsSegment(line.begin, line.end, dimension, box_tolerance))
        {
          return true;
        }
      }
      return false;
    };

    std::vector<int> query_send_counts(size, 0);
    for (const auto &line : line_geometry)
    {
      for (int candidate_rank = 0; candidate_rank < size; candidate_rank++)
      {
        if (RankIntersects(candidate_rank, line))
        {
          query_send_counts[candidate_rank]++;
        }
      }
    }
    std::vector<int> query_receive_counts(size);
    Mpi::Alltoall(1, query_send_counts.data(), query_receive_counts.data(), comm);
    std::vector<int> query_send_offsets, query_receive_offsets;
    const int query_send_total = SetOffsets(query_send_counts, query_send_offsets);
    const int query_receive_total = SetOffsets(query_receive_counts, query_receive_offsets);
    std::vector<int> query_send_indices(query_send_total);
    std::vector<double> query_send_coordinates(6 * query_send_total);
    std::vector<int> query_cursor(query_send_offsets);
    for (std::size_t line_index = 0; line_index < line_geometry.size(); line_index++)
    {
      const auto &line = line_geometry[line_index];
      for (int candidate_rank = 0; candidate_rank < size; candidate_rank++)
      {
        if (!RankIntersects(candidate_rank, line))
        {
          continue;
        }
        const int packed = query_cursor[candidate_rank]++;
        query_send_indices[packed] = static_cast<int>(line_index);
        for (int d = 0; d < 3; d++)
        {
          query_send_coordinates[6 * packed + d] = line.begin[d];
          query_send_coordinates[6 * packed + 3 + d] = line.end[d];
        }
      }
    }
    std::vector<int> query_receive_indices(query_receive_total);
    Mpi::Alltoallv(query_send_indices.data(), query_send_counts.data(),
                   query_send_offsets.data(), query_receive_indices.data(),
                   query_receive_counts.data(), query_receive_offsets.data(), comm);
    const auto query_send_coordinate_counts = ScaleCommunicationPlan(query_send_counts, 6);
    const auto query_send_coordinate_offsets =
        ScaleCommunicationPlan(query_send_offsets, 6);
    const auto query_receive_coordinate_counts =
        ScaleCommunicationPlan(query_receive_counts, 6);
    const auto query_receive_coordinate_offsets =
        ScaleCommunicationPlan(query_receive_offsets, 6);
    std::vector<double> query_receive_coordinates(6 * query_receive_total);
    Mpi::Alltoallv(query_send_coordinates.data(), query_send_coordinate_counts.data(),
                   query_send_coordinate_offsets.data(), query_receive_coordinates.data(),
                   query_receive_coordinate_counts.data(),
                   query_receive_coordinate_offsets.data(), comm);

    std::vector<std::vector<int>> returned_indices(size);
    std::vector<std::vector<double>> returned_bounds(size);
    std::vector<std::pair<double, double>> intersections;
    std::vector<int> candidates;
    for (int source = 0; source < size; source++)
    {
      const int end = query_receive_offsets[source] + query_receive_counts[source];
      for (int packed = query_receive_offsets[source]; packed < end; packed++)
      {
        std::array<double, 3> begin, finish;
        for (int d = 0; d < 3; d++)
        {
          begin[d] = query_receive_coordinates[6 * packed + d];
          finish[d] = query_receive_coordinates[6 * packed + 3 + d];
        }
        locator.FindSegmentIntersections(begin, finish, box_tolerance, intersections,
                                         candidates);
        for (const auto &[interval_begin, interval_end] : intersections)
        {
          returned_indices[source].push_back(query_receive_indices[packed]);
          returned_bounds[source].push_back(interval_begin);
          returned_bounds[source].push_back(interval_end);
        }
      }
    }

    std::vector<int> interval_send_counts(size);
    for (int rank = 0; rank < size; rank++)
    {
      interval_send_counts[rank] = static_cast<int>(returned_indices[rank].size());
    }
    std::vector<int> interval_receive_counts(size);
    Mpi::Alltoall(1, interval_send_counts.data(), interval_receive_counts.data(), comm);
    std::vector<int> interval_send_offsets, interval_receive_offsets;
    const int interval_send_total = SetOffsets(interval_send_counts, interval_send_offsets);
    const int interval_receive_total =
        SetOffsets(interval_receive_counts, interval_receive_offsets);
    std::vector<int> interval_send_indices;
    std::vector<double> interval_send_bounds;
    interval_send_indices.reserve(interval_send_total);
    interval_send_bounds.reserve(2 * interval_send_total);
    for (int rank = 0; rank < size; rank++)
    {
      interval_send_indices.insert(interval_send_indices.end(),
                                   returned_indices[rank].begin(),
                                   returned_indices[rank].end());
      interval_send_bounds.insert(interval_send_bounds.end(), returned_bounds[rank].begin(),
                                  returned_bounds[rank].end());
    }
    std::vector<int> interval_receive_indices(interval_receive_total);
    Mpi::Alltoallv(interval_send_indices.data(), interval_send_counts.data(),
                   interval_send_offsets.data(), interval_receive_indices.data(),
                   interval_receive_counts.data(), interval_receive_offsets.data(), comm);
    const auto interval_send_bound_counts = ScaleCommunicationPlan(interval_send_counts, 2);
    const auto interval_send_bound_offsets =
        ScaleCommunicationPlan(interval_send_offsets, 2);
    const auto interval_receive_bound_counts =
        ScaleCommunicationPlan(interval_receive_counts, 2);
    const auto interval_receive_bound_offsets =
        ScaleCommunicationPlan(interval_receive_offsets, 2);
    std::vector<double> interval_receive_bounds(2 * interval_receive_total);
    Mpi::Alltoallv(interval_send_bounds.data(), interval_send_bound_counts.data(),
                   interval_send_bound_offsets.data(), interval_receive_bounds.data(),
                   interval_receive_bound_counts.data(),
                   interval_receive_bound_offsets.data(), comm);
    for (int i = 0; i < interval_receive_total; i++)
    {
      const int line = interval_receive_indices[i];
      MFEM_ASSERT(line >= 0 && line < static_cast<int>(line_intervals.size()),
                  "Invalid returned Maxwell line index!");
      line_intervals[line].emplace_back(interval_receive_bounds[2 * i],
                                        interval_receive_bounds[2 * i + 1]);
    }
  }
  else
  {
    Mpi::Warning(comm,
                 "Exact Maxwell response-contour integration requires a linear simplex "
                 "mesh; using composite line quadrature instead.\n");
    for (std::size_t line_index = 0; line_index < line_geometry.size(); line_index++)
    {
      const auto &line = line_geometry[line_index];
      std::array<double, 3> tangent{};
      for (int d = 0; d < dimension; d++)
      {
        tangent[d] = line.end[d] - line.begin[d];
      }
      const int segment_count =
          std::max(1, static_cast<int>(std::ceil(16.0 * Norm(tangent) / matching_radius)));
      for (int segment = 0; segment < segment_count; segment++)
      {
        line_intervals[line_index].emplace_back(
            static_cast<double>(segment) / segment_count,
            static_cast<double>(segment + 1) / segment_count);
      }
    }
  }

  struct PendingQuadraturePoint
  {
    std::array<double, 3> coordinate{};
    std::array<double, 3> weighted_tangent{};
  };
  std::vector<PendingQuadraturePoint> pending_points;
  const auto &line_rule =
      mfem::IntRules.Get(mfem::Geometry::SEGMENT, maxwell_quadrature_order);
  for (std::size_t line_index = 0; line_index < line_geometry.size(); line_index++)
  {
    const auto &line = line_geometry[line_index];
    auto &functional = maxwell_lines[line_index];
    functional.point_offset = static_cast<int>(pending_points.size());
    std::vector<std::pair<double, double>> integration_intervals;
    if (exact_intersections)
    {
      std::vector<double> breakpoints = {0.0, 1.0};
      for (auto &[begin, end] : line_intervals[line_index])
      {
        if (begin <= parameter_tolerance)
        {
          begin = 0.0;
        }
        if (end >= 1.0 - parameter_tolerance)
        {
          end = 1.0;
        }
        breakpoints.push_back(begin);
        breakpoints.push_back(end);
      }
      std::sort(breakpoints.begin(), breakpoints.end());
      breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end(),
                                    [](double a, double b)
                                    { return std::abs(a - b) <= parameter_tolerance; }),
                        breakpoints.end());
      for (std::size_t i = 0; i + 1 < breakpoints.size(); i++)
      {
        const double begin = breakpoints[i];
        const double end = breakpoints[i + 1];
        if (end - begin <= parameter_tolerance)
        {
          continue;
        }
        const double midpoint = 0.5 * (begin + end);
        const bool covered = std::any_of(
            line_intervals[line_index].begin(), line_intervals[line_index].end(),
            [&](const auto &interval)
            {
              return midpoint >= interval.first - parameter_tolerance &&
                     midpoint <= interval.second + parameter_tolerance;
            });
        MFEM_VERIFY(covered, "Maxwell response contour line "
                                 << line_index
                                 << " leaves the finite-element mesh over parameter "
                                 << begin << " to " << end << " (begin = " << line.begin[0]
                                 << ", " << line.begin[1] << ", " << line.begin[2]
                                 << "; end = " << line.end[0] << ", " << line.end[1] << ", "
                                 << line.end[2] << ")!");
        integration_intervals.emplace_back(begin, end);
      }
    }
    else
    {
      integration_intervals = std::move(line_intervals[line_index]);
    }

    std::array<double, 3> tangent{};
    for (int d = 0; d < dimension; d++)
    {
      tangent[d] = line.end[d] - line.begin[d];
    }
    for (const auto &[begin, end] : integration_intervals)
    {
      const double interval_length = end - begin;
      for (int q = 0; q < line_rule.GetNPoints(); q++)
      {
        const auto &ip = line_rule.IntPoint(q);
        const double parameter = begin + interval_length * ip.x;
        PendingQuadraturePoint point;
        for (int d = 0; d < dimension; d++)
        {
          point.coordinate[d] = line.begin[d] + parameter * tangent[d];
          point.weighted_tangent[d] = interval_length * ip.weight * tangent[d];
        }
        pending_points.push_back(point);
      }
    }
    functional.point_count =
        static_cast<int>(pending_points.size()) - functional.point_offset;
  }

  mfem::Vector xyz(dimension * pending_points.size());
  std::vector<std::array<double, 3>> weighted_tangents(pending_points.size());
  for (std::size_t i = 0; i < pending_points.size(); i++)
  {
    weighted_tangents[i] = pending_points[i].weighted_tangent;
    for (int d = 0; d < dimension; d++)
    {
      xyz(d * pending_points.size() + i) = pending_points[i].coordinate[d];
    }
  }
  ConfigurePointCommunication(xyz, dimension, &weighted_tangents);
}

void SurfaceResponseOperator::ConfigurePointCommunication(
    const mfem::Vector &xyz, int dimension,
    const std::vector<std::array<double, 3>> *weighted_tangents)
{
  MFEM_VERIFY(dimension == 2 || dimension == 3,
              "Surface response points require dimension two or three!");
  MFEM_VERIFY(xyz.Size() % dimension == 0,
              "Invalid surface-response point-coordinate array!");
  point_query_count = xyz.Size() / dimension;
  MFEM_VERIFY(!weighted_tangents ||
                  static_cast<int>(weighted_tangents->size()) == point_query_count,
              "Invalid surface-response point-tangent array!");

  auto &mesh = const_cast<mfem::ParMesh &>(fespace.GetParMesh());
  const auto comm = fespace.GetComm();
  const int size = Mpi::Size(fespace.GetComm());
  ElementPointLocator locator(mesh, dimension);

  constexpr int routing_box_count = 8;
  constexpr int routing_box_values = 6;
  std::array<double, routing_box_count * routing_box_values> local_routing;
  for (int box = 0; box < routing_box_count; box++)
  {
    for (int d = 0; d < 3; d++)
    {
      local_routing[routing_box_values * box + d] = mfem::infinity();
      local_routing[routing_box_values * box + 3 + d] = -mfem::infinity();
    }
  }
  const auto routing_boxes = locator.GetRoutingBoxes(routing_box_count);
  for (std::size_t box = 0; box < routing_boxes.size(); box++)
  {
    for (int d = 0; d < 3; d++)
    {
      local_routing[routing_box_values * box + d] = routing_boxes[box].min[d];
      local_routing[routing_box_values * box + 3 + d] = routing_boxes[box].max[d];
    }
  }
  std::vector<double> global_routing(size * local_routing.size());
  Mpi::Allgather(static_cast<int>(local_routing.size()), local_routing.data(),
                 global_routing.data(), comm);

  double coordinate_scale = 0.0;
  const auto &bounds = locator.GetBounds();
  for (int d = 0; d < dimension; d++)
  {
    coordinate_scale = std::max({coordinate_scale, std::abs(bounds.min[d]),
                                 std::abs(bounds.max[d]), bounds.max[d] - bounds.min[d]});
  }
  Mpi::GlobalMax(1, &coordinate_scale, comm);
  const double box_tolerance =
      1.0e-11 * coordinate_scale +
      64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, coordinate_scale);

  auto GetPoint = [&](int point)
  {
    std::array<double, 3> coordinate{};
    for (int d = 0; d < dimension; d++)
    {
      coordinate[d] = xyz(d * point_query_count + point);
    }
    return coordinate;
  };
  auto RankContains = [&](int candidate_rank, const std::array<double, 3> &point)
  {
    const int rank_offset = candidate_rank * routing_box_count * routing_box_values;
    for (int box = 0; box < routing_box_count; box++)
    {
      ElementBox bounds;
      for (int d = 0; d < 3; d++)
      {
        bounds.min[d] = global_routing[rank_offset + routing_box_values * box + d];
        bounds.max[d] = global_routing[rank_offset + routing_box_values * box + 3 + d];
      }
      if (bounds.Contains(point, dimension, box_tolerance))
      {
        return true;
      }
    }
    return false;
  };
  auto SetOffsets = [](const std::vector<int> &counts, std::vector<int> &offsets)
  {
    offsets.resize(counts.size());
    int total = 0;
    for (std::size_t i = 0; i < counts.size(); i++)
    {
      offsets[i] = total;
      total += counts[i];
    }
    return total;
  };
  auto ScaleCommunicationPlan = [](const std::vector<int> &values, int scale)
  {
    std::vector<int> result(values);
    for (auto &value : result)
    {
      value *= scale;
    }
    return result;
  };

  std::vector<int> candidate_send_counts(size, 0);
  for (int point = 0; point < point_query_count; point++)
  {
    const auto coordinate = GetPoint(point);
    for (int candidate_rank = 0; candidate_rank < size; candidate_rank++)
    {
      if (RankContains(candidate_rank, coordinate))
      {
        candidate_send_counts[candidate_rank]++;
      }
    }
  }
  std::vector<int> candidate_receive_counts(size);
  Mpi::Alltoall(1, candidate_send_counts.data(), candidate_receive_counts.data(), comm);
  std::vector<int> candidate_send_offsets, candidate_receive_offsets;
  const int candidate_send_total =
      SetOffsets(candidate_send_counts, candidate_send_offsets);
  candidate_query_count += candidate_send_total;
  const int candidate_receive_total =
      SetOffsets(candidate_receive_counts, candidate_receive_offsets);

  std::vector<int> candidate_query_indices(candidate_send_total);
  std::vector<double> candidate_send_coordinates(dimension * candidate_send_total);
  std::vector<int> candidate_cursor(candidate_send_offsets);
  for (int point = 0; point < point_query_count; point++)
  {
    const auto coordinate = GetPoint(point);
    for (int candidate_rank = 0; candidate_rank < size; candidate_rank++)
    {
      if (!RankContains(candidate_rank, coordinate))
      {
        continue;
      }
      const int packed = candidate_cursor[candidate_rank]++;
      candidate_query_indices[packed] = point;
      for (int d = 0; d < dimension; d++)
      {
        candidate_send_coordinates[dimension * packed + d] = coordinate[d];
      }
    }
  }
  const auto candidate_send_coordinate_counts =
      ScaleCommunicationPlan(candidate_send_counts, dimension);
  const auto candidate_send_coordinate_offsets =
      ScaleCommunicationPlan(candidate_send_offsets, dimension);
  const auto candidate_receive_coordinate_counts =
      ScaleCommunicationPlan(candidate_receive_counts, dimension);
  const auto candidate_receive_coordinate_offsets =
      ScaleCommunicationPlan(candidate_receive_offsets, dimension);
  std::vector<double> candidate_receive_coordinates(dimension * candidate_receive_total);
  Mpi::Alltoallv(candidate_send_coordinates.data(), candidate_send_coordinate_counts.data(),
                 candidate_send_coordinate_offsets.data(),
                 candidate_receive_coordinates.data(),
                 candidate_receive_coordinate_counts.data(),
                 candidate_receive_coordinate_offsets.data(), comm);

  std::vector<int> candidate_receive_elements(candidate_receive_total, -1);
  std::vector<double> candidate_receive_references(dimension * candidate_receive_total);
  std::vector<int> candidates;
  for (int i = 0; i < candidate_receive_total; i++)
  {
    std::array<double, 3> coordinate{};
    for (int d = 0; d < dimension; d++)
    {
      coordinate[d] = candidate_receive_coordinates[dimension * i + d];
    }
    mfem::IntegrationPoint reference;
    if (locator.Find(coordinate, box_tolerance, candidate_receive_elements[i], reference,
                     candidates))
    {
      reference.Get(candidate_receive_references.data() + dimension * i, dimension);
    }
  }

  std::vector<int> candidate_result_elements(candidate_send_total);
  Mpi::Alltoallv(candidate_receive_elements.data(), candidate_receive_counts.data(),
                 candidate_receive_offsets.data(), candidate_result_elements.data(),
                 candidate_send_counts.data(), candidate_send_offsets.data(), comm);
  std::vector<double> candidate_result_references(dimension * candidate_send_total);
  Mpi::Alltoallv(
      candidate_receive_references.data(), candidate_receive_coordinate_counts.data(),
      candidate_receive_coordinate_offsets.data(), candidate_result_references.data(),
      candidate_send_coordinate_counts.data(), candidate_send_coordinate_offsets.data(),
      comm);

  std::vector<int> point_owners(point_query_count, size);
  std::vector<int> point_elements(point_query_count, -1);
  std::vector<double> point_references(dimension * point_query_count);
  for (int candidate_rank = 0; candidate_rank < size; candidate_rank++)
  {
    const int begin = candidate_send_offsets[candidate_rank];
    const int end = begin + candidate_send_counts[candidate_rank];
    for (int packed = begin; packed < end; packed++)
    {
      if (candidate_result_elements[packed] < 0)
      {
        continue;
      }
      const int point = candidate_query_indices[packed];
      if (candidate_rank > point_owners[point] ||
          (candidate_rank == point_owners[point] &&
           candidate_result_elements[packed] >= point_elements[point]))
      {
        continue;
      }
      point_owners[point] = candidate_rank;
      point_elements[point] = candidate_result_elements[packed];
      for (int d = 0; d < dimension; d++)
      {
        point_references[dimension * point + d] =
            candidate_result_references[dimension * packed + d];
      }
    }
  }

  std::vector<int> fallback_indices;
  for (int point = 0; point < point_query_count; point++)
  {
    if (point_owners[point] == size)
    {
      fallback_indices.push_back(point);
    }
  }
  fallback_query_count += fallback_indices.size();
  int fallback_count = static_cast<int>(fallback_indices.size());
  Mpi::GlobalSum(1, &fallback_count, comm);
  if (fallback_count > 0)
  {
    Mpi::Warning(
        comm,
        "Distributed surface-response point location could not resolve {:d} contour "
        "points; falling back to FindPointsGSLIB!\n",
        fallback_count);
    mfem::Vector fallback_xyz(dimension * fallback_indices.size());
    for (std::size_t i = 0; i < fallback_indices.size(); i++)
    {
      for (int d = 0; d < dimension; d++)
      {
        fallback_xyz(d * fallback_indices.size() + i) =
            xyz(d * point_query_count + fallback_indices[i]);
      }
    }
    mfem::FindPointsGSLIB finder(comm);
    finder.Setup(mesh, 0.01, 1.0e-12, 256);
    finder.FindPoints(fallback_xyz, mfem::Ordering::byNODES);
    const auto &reference = finder.GetReferencePosition();
    for (std::size_t i = 0; i < fallback_indices.size(); i++)
    {
      const int point = fallback_indices[i];
      MFEM_VERIFY(finder.GetCode()[i] != 2,
                  "Surface-response contour point " << point << " could not be located!");
      point_owners[point] = static_cast<int>(finder.GetProc()[i]);
      point_elements[point] = static_cast<int>(finder.GetElem()[i]);
      for (int d = 0; d < dimension; d++)
      {
        point_references[dimension * point + d] = reference(dimension * i + d);
      }
    }
  }

  point_send_counts.assign(size, 0);
  for (int i = 0; i < point_query_count; i++)
  {
    const int owner = point_owners[i];
    MFEM_VERIFY(owner >= 0 && owner < size,
                "Surface-response contour point " << i << " has no owning rank!");
    point_send_counts[owner]++;
  }
  point_receive_counts.resize(size);
  Mpi::Alltoall(1, point_send_counts.data(), point_receive_counts.data(),
                fespace.GetComm());

  const int send_total = SetOffsets(point_send_counts, point_send_offsets);
  const int receive_total = SetOffsets(point_receive_counts, point_receive_offsets);
  point_send_peer_count = std::count_if(point_send_counts.begin(), point_send_counts.end(),
                                        [](int count) { return count > 0; });
  point_receive_peer_count =
      std::count_if(point_receive_counts.begin(), point_receive_counts.end(),
                    [](int count) { return count > 0; });
  point_send_item_count = send_total;
  point_receive_item_count = receive_total;
  MFEM_ASSERT(send_total == point_query_count, "Invalid point-query communication plan!");

  std::vector<int> send_elements(send_total);
  std::vector<double> send_references(dimension * send_total);
  std::vector<double> send_tangents(weighted_tangents ? 3 * send_total : 0);
  point_send_indices.resize(send_total);
  std::vector<int> cursor(point_send_offsets);
  for (int i = 0; i < point_query_count; i++)
  {
    const int owner = point_owners[i];
    const int packed = cursor[owner]++;
    point_send_indices[packed] = i;
    send_elements[packed] = point_elements[i];
    for (int d = 0; d < dimension; d++)
    {
      send_references[dimension * packed + d] = point_references[dimension * i + d];
    }
    if (weighted_tangents)
    {
      for (int d = 0; d < 3; d++)
      {
        send_tangents[3 * packed + d] =
            (*weighted_tangents)[i][static_cast<std::size_t>(d)];
      }
    }
  }

  std::vector<int> receive_elements(receive_total);
  std::vector<double> receive_references(dimension * receive_total);
  Mpi::Alltoallv(send_elements.data(), point_send_counts.data(), point_send_offsets.data(),
                 receive_elements.data(), point_receive_counts.data(),
                 point_receive_offsets.data(), fespace.GetComm());
  const auto send_reference_counts = ScaleCommunicationPlan(point_send_counts, dimension);
  const auto send_reference_offsets = ScaleCommunicationPlan(point_send_offsets, dimension);
  const auto receive_reference_counts =
      ScaleCommunicationPlan(point_receive_counts, dimension);
  const auto receive_reference_offsets =
      ScaleCommunicationPlan(point_receive_offsets, dimension);
  Mpi::Alltoallv(send_references.data(), send_reference_counts.data(),
                 send_reference_offsets.data(), receive_references.data(),
                 receive_reference_counts.data(), receive_reference_offsets.data(),
                 fespace.GetComm());

  std::vector<double> receive_tangents(weighted_tangents ? 3 * receive_total : 0);
  if (weighted_tangents)
  {
    const auto send_tangent_counts = ScaleCommunicationPlan(point_send_counts, 3);
    const auto send_tangent_offsets = ScaleCommunicationPlan(point_send_offsets, 3);
    const auto receive_tangent_counts = ScaleCommunicationPlan(point_receive_counts, 3);
    const auto receive_tangent_offsets = ScaleCommunicationPlan(point_receive_offsets, 3);
    Mpi::Alltoallv(send_tangents.data(), send_tangent_counts.data(),
                   send_tangent_offsets.data(), receive_tangents.data(),
                   receive_tangent_counts.data(), receive_tangent_offsets.data(),
                   fespace.GetComm());
  }

  point_dof_offsets.resize(receive_total + 1);
  point_dofs.clear();
  point_weights.clear();
  mfem::Array<int> element_dofs;
  mfem::DofTransformation dof_transform;
  mfem::DenseMatrix vector_shape;
  Vector shape;
  for (int i = 0; i < receive_total; i++)
  {
    mfem::IntegrationPoint point;
    if (dimension == 2)
    {
      point.Set2(receive_references[2 * i], receive_references[2 * i + 1]);
    }
    else
    {
      point.Set3(receive_references[3 * i], receive_references[3 * i + 1],
                 receive_references[3 * i + 2]);
    }
    const int element = receive_elements[i];
    const auto &fe = *fespace.Get().GetFE(element);
    if (weighted_tangents)
    {
      MFEM_ASSERT(fe.GetRangeType() == mfem::FiniteElement::VECTOR,
                  "Maxwell contour trace requires a vector finite element!");
      auto &transformation = *fespace.Get().GetElementTransformation(element);
      transformation.SetIntPoint(&point);
      vector_shape.SetSize(fe.GetDof(), fespace.SpaceDimension());
      fe.CalcVShape(transformation, vector_shape);
      shape.SetSize(fe.GetDof());
      shape = 0.0;
      for (int j = 0; j < fe.GetDof(); j++)
      {
        for (int d = 0; d < fespace.SpaceDimension(); d++)
        {
          shape[j] += vector_shape(j, d) * receive_tangents[3 * i + d];
        }
      }
      fespace.Get().GetElementVDofs(element, element_dofs, dof_transform);
      dof_transform.TransformDual(shape);
    }
    else
    {
      shape.SetSize(fe.GetDof());
      fe.CalcShape(point, shape);
      fespace.Get().GetElementDofs(element, element_dofs);
    }

    MFEM_ASSERT(element_dofs.Size() == shape.Size(),
                "Invalid surface-response point interpolation stencil!");
    point_dof_offsets[i] = static_cast<int>(point_dofs.size());
    for (int j = 0; j < element_dofs.Size(); j++)
    {
      const int signed_dof = element_dofs[j];
      point_dofs.push_back(signed_dof >= 0 ? signed_dof : -1 - signed_dof);
      point_weights.push_back((signed_dof >= 0 ? 1.0 : -1.0) * shape[j]);
    }
  }
  point_dof_offsets[receive_total] = static_cast<int>(point_dofs.size());
  stencil_nonzero_count = point_dofs.size();
  point_owned_values.resize(receive_total);
  point_packed_values.resize(point_query_count);
  point_owned_values_pair.resize(2 * receive_total);
  point_packed_values_pair.resize(2 * point_query_count);
  point_send_counts_pair = ScaleCommunicationPlan(point_send_counts, 2);
  point_send_offsets_pair = ScaleCommunicationPlan(point_send_offsets, 2);
  point_receive_counts_pair = ScaleCommunicationPlan(point_receive_counts, 2);
  point_receive_offsets_pair = ScaleCommunicationPlan(point_receive_offsets, 2);
}

void SurfaceResponseOperator::EvaluatePointValues(const Vector &x, Vector &values) const
{
  local_x.SetSize(fespace.GetVSize());
  fespace.GetProlongationMatrix()->Mult(x, local_x);
  const auto *local_data = local_x.HostRead();

  std::fill(point_owned_values.begin(), point_owned_values.end(), 0.0);
  for (std::size_t i = 0; i + 1 < point_dof_offsets.size(); i++)
  {
    for (int j = point_dof_offsets[i]; j < point_dof_offsets[i + 1]; j++)
    {
      point_owned_values[i] += point_weights[j] * local_data[point_dofs[j]];
    }
  }

  Mpi::Alltoallv(point_owned_values.data(), point_receive_counts.data(),
                 point_receive_offsets.data(), point_packed_values.data(),
                 point_send_counts.data(), point_send_offsets.data(), fespace.GetComm());
  values.SetSize(point_query_count);
  auto *value_data = values.HostWrite();
  for (int i = 0; i < point_query_count; i++)
  {
    value_data[point_send_indices[i]] = point_packed_values[i];
  }
}

void SurfaceResponseOperator::EvaluatePointValues(const Vector &xr, const Vector &xi,
                                                  Vector &vr, Vector &vi) const
{
  local_x.SetSize(fespace.GetVSize());
  local_x_imag.SetSize(fespace.GetVSize());
  fespace.GetProlongationMatrix()->Mult(xr, local_x);
  fespace.GetProlongationMatrix()->Mult(xi, local_x_imag);
  const auto *local_real = local_x.HostRead();
  const auto *local_imag = local_x_imag.HostRead();

  std::fill(point_owned_values_pair.begin(), point_owned_values_pair.end(), 0.0);
  for (std::size_t i = 0; i + 1 < point_dof_offsets.size(); i++)
  {
    for (int j = point_dof_offsets[i]; j < point_dof_offsets[i + 1]; j++)
    {
      const double weight = point_weights[j];
      const int dof = point_dofs[j];
      point_owned_values_pair[2 * i] += weight * local_real[dof];
      point_owned_values_pair[2 * i + 1] += weight * local_imag[dof];
    }
  }

  Mpi::Alltoallv(point_owned_values_pair.data(), point_receive_counts_pair.data(),
                 point_receive_offsets_pair.data(), point_packed_values_pair.data(),
                 point_send_counts_pair.data(), point_send_offsets_pair.data(),
                 fespace.GetComm());
  vr.SetSize(point_query_count);
  vi.SetSize(point_query_count);
  auto *real_data = vr.HostWrite();
  auto *imag_data = vi.HostWrite();
  for (int i = 0; i < point_query_count; i++)
  {
    const int point = point_send_indices[i];
    real_data[point] = point_packed_values_pair[2 * i];
    imag_data[point] = point_packed_values_pair[2 * i + 1];
  }
}

void SurfaceResponseOperator::AddPointValuesTranspose(const Vector &values, Vector &y) const
{
  MFEM_ASSERT(values.Size() == point_query_count,
              "Invalid surface-response point-functional vector size!");
  const auto *value_data = values.HostRead();
  for (int i = 0; i < point_query_count; i++)
  {
    point_packed_values[i] = value_data[point_send_indices[i]];
  }
  Mpi::Alltoallv(point_packed_values.data(), point_send_counts.data(),
                 point_send_offsets.data(), point_owned_values.data(),
                 point_receive_counts.data(), point_receive_offsets.data(),
                 fespace.GetComm());

  local_y.SetSize(fespace.GetVSize());
  local_y = 0.0;
  auto *local_data = local_y.HostWrite();
  for (std::size_t i = 0; i + 1 < point_dof_offsets.size(); i++)
  {
    const double value = point_owned_values[i];
    if (value == 0.0)
    {
      continue;
    }
    for (int j = point_dof_offsets[i]; j < point_dof_offsets[i + 1]; j++)
    {
      local_data[point_dofs[j]] += value * point_weights[j];
    }
  }
  y.SetSize(fespace.GetTrueVSize());
  fespace.GetProlongationMatrix()->MultTranspose(local_y, y);
}

void SurfaceResponseOperator::EvaluatePoints(const Vector &x, Vector &values) const
{
  EvaluatePointValues(x, values);
}

void SurfaceResponseOperator::EvaluateMaxwellLines(const Vector &x, Vector &values) const
{
  EvaluatePointValues(x, maxwell_point_values);
  const auto *point_data = maxwell_point_values.HostRead();
  values.SetSize(static_cast<int>(maxwell_lines.size()));
  values = 0.0;
  auto *value_data = values.HostWrite();
  for (std::size_t line_index = 0; line_index < maxwell_lines.size(); line_index++)
  {
    const auto &line = maxwell_lines[line_index];
    for (int q = 0; q < line.point_count; q++)
    {
      value_data[line_index] += point_data[line.point_offset + q];
    }
  }
}

void SurfaceResponseOperator::EvaluateMaxwellLines(const Vector &xr, const Vector &xi,
                                                   Vector &vr, Vector &vi) const
{
  Vector point_real, point_imag;
  EvaluatePointValues(xr, xi, point_real, point_imag);
  const auto *point_real_data = point_real.HostRead();
  const auto *point_imag_data = point_imag.HostRead();
  vr.SetSize(static_cast<int>(maxwell_lines.size()));
  vi.SetSize(static_cast<int>(maxwell_lines.size()));
  vr = 0.0;
  vi = 0.0;
  auto *real_data = vr.HostWrite();
  auto *imag_data = vi.HostWrite();
  for (std::size_t line_index = 0; line_index < maxwell_lines.size(); line_index++)
  {
    const auto &line = maxwell_lines[line_index];
    for (int q = 0; q < line.point_count; q++)
    {
      real_data[line_index] += point_real_data[line.point_offset + q];
      imag_data[line_index] += point_imag_data[line.point_offset + q];
    }
  }
}

void SurfaceResponseOperator::AddMaxwellLinesTranspose(const Vector &values,
                                                       Vector &y) const
{
  MFEM_ASSERT(values.Size() == static_cast<int>(maxwell_lines.size()),
              "Invalid Maxwell line-functional vector size!");
  maxwell_point_values.SetSize(point_query_count);
  maxwell_point_values = 0.0;
  auto *point_data = maxwell_point_values.HostWrite();
  const auto *value_data = values.HostRead();
  for (std::size_t line_index = 0; line_index < maxwell_lines.size(); line_index++)
  {
    const auto &line = maxwell_lines[line_index];
    for (int q = 0; q < line.point_count; q++)
    {
      point_data[line.point_offset + q] = value_data[line_index];
    }
  }
  AddPointValuesTranspose(maxwell_point_values, y);
}

void SurfaceResponseOperator::BuildMaxwellTrace(const Vector &line_values,
                                                Vector &values) const
{
  MFEM_ASSERT(line_values.Size() == static_cast<int>(maxwell_lines.size()),
              "Inconsistent Maxwell contour-path data!");
  values.SetSize(basis_size);
  values = 0.0;
  auto PathIntegral = [&](const MaxwellContourPath &path)
  {
    double value = 0.0;
    if (path.anchor_line >= 0)
    {
      value += line_values[path.anchor_line];
    }
    for (int i = 0; i < path.contour_line_count; i++)
    {
      value += line_values[path.contour_line_offset + i];
    }
    if (path.end_line >= 0)
    {
      value += line_values[path.end_line];
    }
    return value;
  };
  for (const auto &conductor_path : maxwell_conductor_paths)
  {
    const double parent = conductor_path.parent_trace_offset >= 0
                              ? values[conductor_path.parent_trace_offset]
                              : 0.0;
    values[conductor_path.trace_offset] =
        parent + conductor_path.integral_sign *
                     PathIntegral(maxwell_paths[conductor_path.contour_path]);
  }
  for (const auto &path : maxwell_paths)
  {
    MFEM_ASSERT(!path.trace_indices.empty(), "Empty Maxwell contour path!");
    double value = 0.0;
    if (!path.closed && path.start_conductor_trace >= 0)
    {
      value = values[path.start_conductor_trace];
    }
    if (path.anchor_line >= 0)
    {
      value -= line_values[path.anchor_line];
    }
    values[path.trace_indices.front()] = value;
    for (std::size_t i = 1; i < path.trace_indices.size(); i++)
    {
      value -= line_values[path.contour_line_offset + static_cast<int>(i) - 1];
      values[path.trace_indices[i]] = value;
    }
  }
}

void SurfaceResponseOperator::BuildMaxwellTraceTranspose(const Vector &values,
                                                         Vector &line_values) const
{
  MFEM_ASSERT(values.Size() == basis_size, "Inconsistent Maxwell contour-path data!");
  line_values.SetSize(static_cast<int>(maxwell_lines.size()));
  line_values = 0.0;
  maxwell_conductor_adjoint.SetSize(basis_size);
  maxwell_conductor_adjoint = 0.0;
  for (const auto &path : maxwell_conductor_paths)
  {
    maxwell_conductor_adjoint[path.trace_offset] = values[path.trace_offset];
  }
  for (const auto &path : maxwell_paths)
  {
    const int size = static_cast<int>(path.trace_indices.size());
    MFEM_ASSERT(size > 0, "Empty Maxwell contour path!");
    maxwell_path_adjoint.SetSize(size);
    for (int i = 0; i < size; i++)
    {
      maxwell_path_adjoint[i] = values[path.trace_indices[i]];
    }
    for (int i = size - 2; i >= 0; i--)
    {
      line_values[path.contour_line_offset + i] -= maxwell_path_adjoint[i + 1];
      maxwell_path_adjoint[i] += maxwell_path_adjoint[i + 1];
    }
    if (path.anchor_line >= 0)
    {
      line_values[path.anchor_line] -= maxwell_path_adjoint[0];
    }
    if (!path.closed && path.start_conductor_trace >= 0)
    {
      maxwell_conductor_adjoint[path.start_conductor_trace] += maxwell_path_adjoint[0];
    }
  }
  auto AddPathIntegralTranspose = [&](const MaxwellContourPath &path, double value)
  {
    if (path.anchor_line >= 0)
    {
      line_values[path.anchor_line] += value;
    }
    for (int i = 0; i < path.contour_line_count; i++)
    {
      line_values[path.contour_line_offset + i] += value;
    }
    if (path.end_line >= 0)
    {
      line_values[path.end_line] += value;
    }
  };
  for (auto path = maxwell_conductor_paths.rbegin(); path != maxwell_conductor_paths.rend();
       path++)
  {
    const double value = maxwell_conductor_adjoint[path->trace_offset];
    AddPathIntegralTranspose(maxwell_paths[path->contour_path],
                             path->integral_sign * value);
    if (path->parent_trace_offset >= 0)
    {
      maxwell_conductor_adjoint[path->parent_trace_offset] += value;
    }
  }
}

void SurfaceResponseOperator::ApplyTrace(const Vector &x, Vector &values) const
{
  trace_forward_count++;
  if (maxwell)
  {
    EvaluateMaxwellLines(x, correction);
    BuildMaxwellTrace(correction, values);
    return;
  }
  EvaluatePoints(x, correction);
  values.SetSize(basis_size);
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    const double reference = correction(patch.point_offset + model.contour_size);
    for (int i = 0; i < model.contour_size; i++)
    {
      values(patch.trace_offset + i) = correction(patch.point_offset + i) - reference;
    }
    for (int state = 0; state < model.conductor_state_count; state++)
    {
      values(patch.trace_offset + model.contour_size + state) =
          correction(patch.point_offset + model.contour_size + 1 + state) - reference;
    }
  }
}

void SurfaceResponseOperator::ApplyTraceTranspose(const Vector &values, Vector &y) const
{
  trace_transpose_count++;
  if (maxwell)
  {
    BuildMaxwellTraceTranspose(values, correction);
    AddMaxwellLinesTranspose(correction, y);
    return;
  }
  correction.SetSize(point_query_count);
  correction = 0.0;
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    double reference = 0.0;
    for (int i = 0; i < model.contour_size; i++)
    {
      const double value = values(patch.trace_offset + i);
      correction[patch.point_offset + i] = value;
      reference -= value;
    }
    for (int state = 0; state < model.conductor_state_count; state++)
    {
      const double value = values(patch.trace_offset + model.contour_size + state);
      correction[patch.point_offset + model.contour_size + 1 + state] = value;
      reference -= value;
    }
    correction[patch.point_offset + model.contour_size] = reference;
  }
  AddPointValuesTranspose(correction, y);
}

void SurfaceResponseOperator::ApplyUneliminated(const Vector &x, Vector &y) const
{
  BlockTimer timer(Timer::RESPONSE_APPLY);
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
  operator_mult_count++;
  x_free.SetSize(x.Size());
  x_free = x;
  x_free.SetSubVector(dbc_tdof_list, 0.0);
  ApplyUneliminated(x_free, y);
  y.SetSubVector(dbc_tdof_list, 0.0);
}

void SurfaceResponseOperator::EliminateRHS(const Vector &x, Vector &rhs) const
{
  eliminate_rhs_count++;
  ApplyUneliminated(x, correction);
  correction.SetSubVector(dbc_tdof_list, 0.0);
  rhs.Add(-1.0, correction);
}

SurfaceResponseOperator::EnergyCorrection
SurfaceResponseOperator::GetEnergyCorrection(const Vector &x) const
{
  ApplyTrace(x, trace);
  EnergyCorrection energy;
  for (const auto &model : models)
  {
    for (const auto &[interface, defect] : model.surface_defects)
    {
      (void)defect;
      energy.interfaces.try_emplace(interface, 0.0);
    }
  }
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
  std::vector<double> reduction;
  reduction.reserve(1 + energy.interfaces.size());
  reduction.push_back(energy.domain);
  for (const auto &[interface, value] : energy.interfaces)
  {
    (void)interface;
    reduction.push_back(value);
  }
  Mpi::GlobalSum(static_cast<int>(reduction.size()), reduction.data(), fespace.GetComm());
  energy.domain = reduction[0];
  std::size_t i = 1;
  for (auto &[interface, value] : energy.interfaces)
  {
    (void)interface;
    value = reduction[i++];
  }
  return energy;
}

std::map<int, double>
SurfaceResponseOperator::GetFabricatedSurfaceEnergy(const Vector &x) const
{
  ApplyTrace(x, trace);
  std::map<int, double> energy;
  for (const auto &model : models)
  {
    for (const auto &[interface, matrix] : model.fabricated_surfaces)
    {
      (void)matrix;
      energy.try_emplace(interface, 0.0);
    }
  }
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
  std::vector<double> reduction;
  reduction.reserve(energy.size());
  for (const auto &[interface, value] : energy)
  {
    (void)interface;
    reduction.push_back(value);
  }
  Mpi::GlobalSum(static_cast<int>(reduction.size()), reduction.data(), fespace.GetComm());
  std::size_t i = 0;
  for (auto &[interface, value] : energy)
  {
    (void)interface;
    value = reduction[i++];
  }
  return energy;
}

SurfaceResponseOperator::ElectrostaticResponse
SurfaceResponseOperator::GetElectrostaticResponse(const Vector &x,
                                                  bool include_fixed_flux) const
{
  ApplyTrace(x, trace);
  ElectrostaticResponse result;
  double weighted_trace_closure_spread_squared = 0.0;
  double trace_closure_response_weight = 0.0;
  double failed_trace_closure_response_weight = 0.0;
  Vector fixed_flux;
  for (const auto &model : models)
  {
    ModelContribution contribution;
    contribution.model = model.idx;
    for (const auto &[interface, matrix] : model.fabricated_surfaces)
    {
      (void)matrix;
      result.fabricated_surface_energy.try_emplace(interface, 0.0);
      contribution.fabricated_surface_energy.try_emplace(interface, 0.0);
      if (include_fixed_flux)
      {
        result.fabricated_surface_energy_fixed_flux.try_emplace(interface, 0.0);
        contribution.fabricated_surface_energy_fixed_flux.try_emplace(interface, 0.0);
      }
    }
    result.model_contributions.push_back(std::move(contribution));
  }
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    auto &contribution = result.model_contributions[patch.model];
    contribution.patch_count += 1.0;
    contribution.patch_weight += patch.weight;
    Vector patch_trace(trace.GetData() + patch.trace_offset, model.basis_size);
    const double domain_correction =
        0.5 * patch.weight * QuadraticForm(model.domain_defect, patch_trace, response);
    result.domain_correction += domain_correction;
    contribution.domain_correction += domain_correction;
    if (include_fixed_flux)
    {
      fixed_flux.SetSize(model.basis_size);
      model.fixed_flux_transform.Mult(patch_trace, fixed_flux);
      const double domain_correction_fixed_flux =
          0.5 * patch.weight *
          (QuadraticForm(model.fabricated_domain, fixed_flux, response) -
           QuadraticForm(model.thin_domain, patch_trace, response));
      result.domain_correction_fixed_flux += domain_correction_fixed_flux;
      contribution.domain_correction_fixed_flux += domain_correction_fixed_flux;
    }
    for (const auto &[interface, matrix] : model.fabricated_surfaces)
    {
      const double fixed_trace_energy =
          patch.weight * QuadraticForm(matrix, patch_trace, response);
      result.fabricated_surface_energy[interface] += fixed_trace_energy;
      contribution.fabricated_surface_energy[interface] += fixed_trace_energy;
      if (include_fixed_flux)
      {
        const double fixed_flux_energy =
            patch.weight * QuadraticForm(matrix, fixed_flux, response);
        result.fabricated_surface_energy_fixed_flux[interface] += fixed_flux_energy;
        contribution.fabricated_surface_energy_fixed_flux[interface] += fixed_flux_energy;
        const double weight =
            std::max(std::abs(fixed_trace_energy), std::abs(fixed_flux_energy));
        if (weight > 0.0)
        {
          const double spread = std::abs(fixed_trace_energy - fixed_flux_energy) / weight;
          weighted_trace_closure_spread_squared += weight * spread * spread;
          failed_trace_closure_response_weight +=
              spread > maximum_trace_closure_spread ? weight : 0.0;
          trace_closure_response_weight += weight;
        }
      }
    }
  }
  const std::size_t values_per_interface = include_fixed_flux ? 2 : 1;
  std::vector<double> reduction;
  reduction.reserve(
      1 + (include_fixed_flux ? 4 : 0) +
      values_per_interface * result.fabricated_surface_energy.size() +
      result.model_contributions.size() *
          (4 + values_per_interface * result.fabricated_surface_energy.size()));
  reduction.push_back(result.domain_correction);
  if (include_fixed_flux)
  {
    reduction.push_back(result.domain_correction_fixed_flux);
    reduction.push_back(weighted_trace_closure_spread_squared);
    reduction.push_back(trace_closure_response_weight);
    reduction.push_back(failed_trace_closure_response_weight);
  }
  for (const auto &[interface, fixed_trace] : result.fabricated_surface_energy)
  {
    reduction.push_back(fixed_trace);
    if (include_fixed_flux)
    {
      reduction.push_back(result.fabricated_surface_energy_fixed_flux.at(interface));
    }
  }
  for (const auto &contribution : result.model_contributions)
  {
    reduction.push_back(contribution.patch_count);
    reduction.push_back(contribution.patch_weight);
    reduction.push_back(contribution.domain_correction);
    if (include_fixed_flux)
    {
      reduction.push_back(contribution.domain_correction_fixed_flux);
    }
    for (const auto &[interface, fixed_trace] : contribution.fabricated_surface_energy)
    {
      (void)interface;
      reduction.push_back(fixed_trace);
      if (include_fixed_flux)
      {
        reduction.push_back(
            contribution.fabricated_surface_energy_fixed_flux.at(interface));
      }
    }
  }
  Mpi::GlobalSum(static_cast<int>(reduction.size()), reduction.data(), fespace.GetComm());
  std::size_t i = 0;
  result.domain_correction = reduction[i++];
  if (include_fixed_flux)
  {
    result.domain_correction_fixed_flux = reduction[i++];
    weighted_trace_closure_spread_squared = reduction[i++];
    trace_closure_response_weight = reduction[i++];
    failed_trace_closure_response_weight = reduction[i++];
    if (trace_closure_response_weight > 0.0)
    {
      result.response_weighted_trace_closure_spread =
          std::sqrt(weighted_trace_closure_spread_squared / trace_closure_response_weight);
      result.trace_closure_response_failure_fraction =
          failed_trace_closure_response_weight / trace_closure_response_weight;
    }
  }
  for (auto &[interface, fixed_trace] : result.fabricated_surface_energy)
  {
    fixed_trace = reduction[i++];
    if (include_fixed_flux)
    {
      result.fabricated_surface_energy_fixed_flux.at(interface) = reduction[i++];
    }
  }
  for (auto &contribution : result.model_contributions)
  {
    contribution.patch_count = reduction[i++];
    contribution.patch_weight = reduction[i++];
    contribution.domain_correction = reduction[i++];
    if (include_fixed_flux)
    {
      contribution.domain_correction_fixed_flux = reduction[i++];
    }
    for (auto &[interface, fixed_trace] : contribution.fabricated_surface_energy)
    {
      fixed_trace = reduction[i++];
      if (include_fixed_flux)
      {
        contribution.fabricated_surface_energy_fixed_flux.at(interface) = reduction[i++];
      }
    }
  }
  MFEM_ASSERT(i == reduction.size(), "Incorrect batched model-contribution reduction!");

  if (!include_fixed_flux)
  {
    return result;
  }
  for (const auto &[interface, fixed_trace] : result.fabricated_surface_energy)
  {
    const auto fixed_flux_energy =
        result.fabricated_surface_energy_fixed_flux.find(interface);
    MFEM_ASSERT(fixed_flux_energy != result.fabricated_surface_energy_fixed_flux.end(),
                "Missing fixed-flux fabricated surface response!");
    const double scale =
        std::max(std::abs(fixed_trace), std::abs(fixed_flux_energy->second));
    const double spread =
        scale > 0.0 ? std::abs(fixed_trace - fixed_flux_energy->second) / scale : 0.0;
    result.trace_closure_spread[interface] = spread;
    result.maximum_trace_closure_spread =
        std::max(result.maximum_trace_closure_spread, spread);
  }
  result.confident =
      result.maximum_trace_closure_spread <= maximum_trace_closure_spread &&
      result.response_weighted_trace_closure_spread <= maximum_trace_closure_spread &&
      result.trace_closure_response_failure_fraction <=
          maximum_trace_closure_response_failure_fraction;
  return result;
}

SurfaceResponseOperator::MaxwellResponse
SurfaceResponseOperator::GetMaxwellResponse(const GridFunction &E,
                                            std::complex<double> omega) const
{
  MFEM_VERIFY(maxwell && maxwell_contours.size() == patches.size(),
              "Maxwell surface response was not configured!");
  MFEM_ASSERT(maxwell_conductor_anchors.size() == patches.size(),
              "Inconsistent Maxwell response anchor count!");
  MFEM_VERIFY(E.HasImag(), "Maxwell surface response requires a complex electric field!");

  MaxwellResponse result;
  result.kR = std::abs(omega) * matching_radius / minimum_wave_speed;
  result.matched_length_fraction = matched_length_fraction;
  result.corner_neighborhood_fraction = corner_neighborhood_fraction;
  result.matched_length_fraction_by_interface = matched_length_fraction_by_interface;
  result.corner_neighborhood_fraction_by_interface =
      corner_neighborhood_fraction_by_interface;
  result.maximum_curvature_ratio = maximum_curvature_ratio;
  result.maximum_library_distance = maximum_library_distance;
  result.boundary_law_verified = boundary_law_verified;
  for (const auto &model : models)
  {
    for (const auto &[interface, matrix] : model.fabricated_surfaces)
    {
      (void)matrix;
      result.fabricated_surface_energy.try_emplace(interface, 0.0);
      result.fabricated_surface_energy_fixed_flux.try_emplace(interface, 0.0);
    }
  }

  Vector field_real, field_imag, line_real, line_imag, trace_real, trace_imag;
  E.Real().GetTrueDofs(field_real);
  E.Imag().GetTrueDofs(field_imag);
  EvaluateMaxwellLines(field_real, field_imag, line_real, line_imag);
  BuildMaxwellTrace(line_real, trace_real);
  BuildMaxwellTrace(line_imag, trace_imag);

  constexpr double maximum_path_loop_residual = 0.05;
  double weighted_loop_residual_squared = 0.0;
  double loop_response_weight = 0.0;
  double failed_loop_response_weight = 0.0;
  double weighted_trace_closure_spread_squared = 0.0;
  double trace_closure_response_weight = 0.0;
  double failed_trace_closure_response_weight = 0.0;
  Vector fixed_flux_real, fixed_flux_imag, workspace;
  for (std::size_t patch_index = 0; patch_index < patches.size(); patch_index++)
  {
    const auto &patch = patches[patch_index];
    const auto &model = models[patch.model];
    Vector patch_trace_real(trace_real.GetData() + patch.trace_offset, model.basis_size);
    Vector patch_trace_imag(trace_imag.GetData() + patch.trace_offset, model.basis_size);
    fixed_flux_real.SetSize(model.basis_size);
    fixed_flux_imag.SetSize(model.basis_size);

    std::vector<std::pair<double, double>> path_diagnostics;
    const auto [path_begin, path_count] = maxwell_patch_paths[patch_index];
    path_diagnostics.reserve(path_count);
    for (int path_index = path_begin; path_index < path_begin + path_count; path_index++)
    {
      const auto &path = maxwell_paths[path_index];
      std::complex<double> loop_integral = 0.0;
      double loop_scale = 0.0;
      auto AddLine = [&](int line, double scale = 1.0)
      {
        if (line < 0)
        {
          return;
        }
        const std::complex<double> integral = {line_real[line], line_imag[line]};
        loop_integral += scale * integral;
        loop_scale += std::abs(integral);
      };
      auto AddTrace = [&](int trace_index, double scale)
      {
        if (trace_index < 0)
        {
          return;
        }
        const std::complex<double> value = {trace_real[trace_index],
                                            trace_imag[trace_index]};
        loop_integral += scale * value;
        loop_scale += std::abs(value);
      };
      if (!path.closed)
      {
        AddLine(path.anchor_line);
      }
      for (int i = 0; i < path.contour_line_count; i++)
      {
        AddLine(path.contour_line_offset + i);
      }
      if (!path.closed)
      {
        AddLine(path.end_line);
        AddTrace(path.start_conductor_trace, -1.0);
        AddTrace(path.end_conductor_trace, 1.0);
      }
      if (loop_scale > 0.0)
      {
        const double residual = std::abs(loop_integral) / loop_scale;
        result.loop_residual = std::max(result.loop_residual, residual);
        path_diagnostics.emplace_back(residual, loop_scale * loop_scale);
      }
    }

    model.fixed_flux_transform.Mult(patch_trace_real, fixed_flux_real);
    model.fixed_flux_transform.Mult(patch_trace_imag, fixed_flux_imag);
    auto HermitianForm =
        [&](const mfem::DenseMatrix &matrix, const Vector &real, const Vector &imag)
    {
      return QuadraticForm(matrix, real, workspace) +
             QuadraticForm(matrix, imag, workspace);
    };
    result.domain_correction +=
        0.5 * patch.weight *
        HermitianForm(model.domain_defect, patch_trace_real, patch_trace_imag);
    result.domain_correction_fixed_flux +=
        0.5 * patch.weight *
        (HermitianForm(model.fabricated_domain, fixed_flux_real, fixed_flux_imag) -
         HermitianForm(model.thin_domain, patch_trace_real, patch_trace_imag));
    double patch_response_energy = 0.0;
    for (const auto &[interface, matrix] : model.fabricated_surfaces)
    {
      const double fixed_trace_energy =
          patch.weight * HermitianForm(matrix, patch_trace_real, patch_trace_imag);
      const double fixed_flux_energy =
          patch.weight * HermitianForm(matrix, fixed_flux_real, fixed_flux_imag);
      result.fabricated_surface_energy[interface] += fixed_trace_energy;
      result.fabricated_surface_energy_fixed_flux[interface] += fixed_flux_energy;
      const double weight =
          std::max(std::abs(fixed_trace_energy), std::abs(fixed_flux_energy));
      patch_response_energy += weight;
      if (weight > 0.0)
      {
        const double spread = std::abs(fixed_trace_energy - fixed_flux_energy) / weight;
        weighted_trace_closure_spread_squared += weight * spread * spread;
        trace_closure_response_weight += weight;
        if (spread > maximum_trace_closure_spread)
        {
          failed_trace_closure_response_weight += weight;
        }
      }
    }
    double path_scale_squared = 0.0;
    for (const auto &[residual, scale_squared] : path_diagnostics)
    {
      (void)residual;
      path_scale_squared += scale_squared;
    }
    if (patch_response_energy > 0.0 && path_scale_squared > 0.0)
    {
      for (const auto &[residual, scale_squared] : path_diagnostics)
      {
        const double weight = patch_response_energy * scale_squared / path_scale_squared;
        weighted_loop_residual_squared += weight * residual * residual;
        loop_response_weight += weight;
        if (residual > maximum_path_loop_residual)
        {
          failed_loop_response_weight += weight;
        }
      }
    }
  }
  std::vector<double> reduction;
  reduction.reserve(8 + 2 * result.fabricated_surface_energy.size());
  reduction.push_back(result.domain_correction);
  reduction.push_back(result.domain_correction_fixed_flux);
  reduction.push_back(weighted_loop_residual_squared);
  reduction.push_back(loop_response_weight);
  reduction.push_back(failed_loop_response_weight);
  reduction.push_back(weighted_trace_closure_spread_squared);
  reduction.push_back(trace_closure_response_weight);
  reduction.push_back(failed_trace_closure_response_weight);
  for (const auto &[interface, fixed_trace] : result.fabricated_surface_energy)
  {
    reduction.push_back(fixed_trace);
    reduction.push_back(result.fabricated_surface_energy_fixed_flux.at(interface));
  }
  Mpi::GlobalSum(static_cast<int>(reduction.size()), reduction.data(), fespace.GetComm());
  std::size_t i = 0;
  result.domain_correction = reduction[i++];
  result.domain_correction_fixed_flux = reduction[i++];
  weighted_loop_residual_squared = reduction[i++];
  loop_response_weight = reduction[i++];
  failed_loop_response_weight = reduction[i++];
  weighted_trace_closure_spread_squared = reduction[i++];
  trace_closure_response_weight = reduction[i++];
  failed_trace_closure_response_weight = reduction[i++];
  if (loop_response_weight > 0.0)
  {
    result.response_weighted_loop_residual =
        std::sqrt(weighted_loop_residual_squared / loop_response_weight);
    result.loop_response_failure_fraction =
        failed_loop_response_weight / loop_response_weight;
  }
  if (trace_closure_response_weight > 0.0)
  {
    result.response_weighted_trace_closure_spread =
        std::sqrt(weighted_trace_closure_spread_squared / trace_closure_response_weight);
    result.trace_closure_response_failure_fraction =
        failed_trace_closure_response_weight / trace_closure_response_weight;
  }
  Mpi::GlobalMax(1, &result.loop_residual, fespace.GetComm());
  for (auto &[interface, fixed_trace] : result.fabricated_surface_energy)
  {
    auto &fixed_flux = result.fabricated_surface_energy_fixed_flux.at(interface);
    fixed_trace = reduction[i++];
    fixed_flux = reduction[i++];
  }

  constexpr double maximum_kR = 0.1;
  constexpr double maximum_weighted_loop_residual = 0.05;
  constexpr double maximum_loop_response_failure_fraction = 0.01;
  constexpr double maximum_corner_fraction = 0.1;
  constexpr double maximum_curvature = 0.25;
  constexpr double minimum_coverage = 1.0 - 1.0e-10;
  constexpr double maximum_library_match_distance = 0.8;
  // Fixed-trace and fixed-flux closures are both admissible in postprocessing-only
  // correction. A material difference between them means that the unresolved
  // fabricated field is not determined accurately enough by the thin-model trace.
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
  result.closure_independent_confident =
      result.kR <= maximum_kR &&
      result.response_weighted_loop_residual <= maximum_weighted_loop_residual &&
      result.loop_response_failure_fraction <= maximum_loop_response_failure_fraction &&
      result.matched_length_fraction >= minimum_coverage &&
      result.corner_neighborhood_fraction <= maximum_corner_fraction &&
      result.maximum_curvature_ratio <= maximum_curvature &&
      result.maximum_library_distance <= maximum_library_match_distance &&
      result.boundary_law_verified;
  result.confident =
      result.closure_independent_confident &&
      result.maximum_trace_closure_spread <= maximum_trace_closure_spread &&
      result.response_weighted_trace_closure_spread <= maximum_trace_closure_spread &&
      result.trace_closure_response_failure_fraction <=
          maximum_trace_closure_response_failure_fraction;
  return result;
}

nlohmann::json SurfaceResponseOperator::GetStatistics() const
{
  const auto comm = fespace.GetComm();
  auto Replicated = [comm](long long int value, const char *name)
  {
    long long int minimum = value;
    long long int maximum = value;
    Mpi::GlobalMin(1, &minimum, comm);
    Mpi::GlobalMax(1, &maximum, comm);
    MFEM_VERIFY(minimum == maximum,
                "Rank-inconsistent surface-response statistic \"" << name << "\"!");
    return minimum;
  };
  auto Distribution = [comm](long long int value)
  {
    long long int total = value;
    long long int minimum = value;
    long long int maximum = value;
    int nonzero = value > 0 ? 1 : 0;
    Mpi::GlobalSum(1, &total, comm);
    Mpi::GlobalMin(1, &minimum, comm);
    Mpi::GlobalMax(1, &maximum, comm);
    Mpi::GlobalSum(1, &nonzero, comm);
    return nlohmann::json{{"Total", total},
                          {"Minimum", minimum},
                          {"Maximum", maximum},
                          {"NonzeroRanks", nonzero}};
  };

  nlohmann::json result = automatic_statistics.is_null() ? nlohmann::json{{"Version", 1}}
                                                         : automatic_statistics;
  std::vector<long long int> model_patch_counts(models.size(), 0);
  std::vector<double> model_patch_weights(models.size(), 0.0);
  for (const auto &patch : patches)
  {
    model_patch_counts[patch.model]++;
    model_patch_weights[patch.model] += patch.weight;
  }
  Mpi::GlobalSum(static_cast<int>(model_patch_counts.size()), model_patch_counts.data(),
                 comm);
  Mpi::GlobalSum(static_cast<int>(model_patch_weights.size()), model_patch_weights.data(),
                 comm);
  nlohmann::json model_catalog = nlohmann::json::array();
  for (std::size_t i = 0; i < models.size(); i++)
  {
    model_catalog.push_back({{"Index", models[i].idx},
                             {"Name", models[i].name},
                             {"Topology", models[i].topology},
                             {"BasisSize", models[i].basis_size},
                             {"PatchCount", model_patch_counts[i]},
                             {"PatchWeight", model_patch_weights[i]}});
  }
  result["ModelCatalog"] = std::move(model_catalog);
  result["Correction"] = {{"Models", Replicated(models.size(), "Models")},
                          {"Patches", global_patch_count},
                          {"TraceCoefficients", global_basis_size},
                          {"LocalPatches", Distribution(patches.size())},
                          {"LocalTraceCoefficients", Distribution(basis_size)},
                          {"ContourLines", Distribution(contour_line_count)}};
  const long long int stencil_rows =
      point_dof_offsets.empty() ? 0 : point_dof_offsets.size() - 1;
  result["Interpolation"] = {{"PointQueries", Distribution(point_query_count)},
                             {"StencilRows", Distribution(stencil_rows)},
                             {"StencilNonzeros", Distribution(stencil_nonzero_count)},
                             {"CandidateQueries", Distribution(candidate_query_count)},
                             {"FallbackQueries", Distribution(fallback_query_count)}};
  auto send_items = Distribution(point_send_item_count);
  auto receive_items = Distribution(point_receive_item_count);
  MFEM_VERIFY(send_items["Total"] == receive_items["Total"],
              "Surface-response point communication has mismatched global item counts!");
  result["Communication"] = {
      {"PointSendPeers", Distribution(point_send_peer_count)},
      {"PointReceivePeers", Distribution(point_receive_peer_count)},
      {"PointSendItems", std::move(send_items)},
      {"PointReceiveItems", std::move(receive_items)},
      {"ApplyPayloadScalarsPerDirection", Distribution(point_send_item_count)}};
  result["Runtime"] = {
      {"OperatorMultCalls", Replicated(operator_mult_count, "OperatorMultCalls")},
      {"EliminateRHSCalls", Replicated(eliminate_rhs_count, "EliminateRHSCalls")},
      {"TraceForwardCalls", Replicated(trace_forward_count, "TraceForwardCalls")},
      {"TraceTransposeCalls", Replicated(trace_transpose_count, "TraceTransposeCalls")}};
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
  Mpi::GlobalSum(1, &weight, fespace.GetComm());
  return weight;
}

}  // namespace palace
