// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularsurfacepostoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <fmt/format.h>

#include "fem/coefficient.hpp"
#include "fem/singularelements.hpp"
#include "fem/singularfield.hpp"
#include "fem/singulargeometry.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

struct TriangleBoundarySide
{
  int element;
  int attribute;
  mfem::IntegrationPoint point;
  std::array<int, 2> endpoint_nodes{-1, -1};
};

struct TetrahedronBoundarySide
{
  int element;
  int attribute;
  mfem::IntegrationPoint point;
};

struct TetrahedronFaceSingularity
{
  fem::singular::TetrahedronFaceSingularityType type;
  std::array<int, 2> nodes{-1, -1};
  double nu;
};

struct PhysicalFaceEdge
{
  std::array<int, 2> nodes{-1, -1};
  mfem::Vector start;
  mfem::Vector direction;
  double length_squared = 0.0;

  double Distance(const mfem::Vector &position) const
  {
    mfem::Vector displacement(position);
    displacement -= start;
    const double projection =
        std::clamp((displacement * direction) / length_squared, 0.0, 1.0);
    displacement.Add(-projection, direction);
    return displacement.Norml2();
  }
};

struct FaceCutoffRay
{
  double radius;
  int active_edge;
};

struct PhysicalSingularTip
{
  std::array<double, 2> position{};
  double nu = 1.0;
};

std::array<double, 3> GetTriangleBarycentricCoordinates(const mfem::IntegrationPoint &point)
{
  return {1.0 - point.x - point.y, point.x, point.y};
}

int FindTriangleVertex(const mfem::IntegrationPoint &point)
{
  const auto lambda = GetTriangleBarycentricCoordinates(point);
  const auto maximum = std::max_element(lambda.begin(), lambda.end());
  if (*maximum < 1.0 - 256.0 * std::numeric_limits<double>::epsilon())
  {
    return -1;
  }
  return static_cast<int>(std::distance(lambda.begin(), maximum));
}

std::array<int, 2> GetTriangleEdgeNodes(const std::array<mfem::IntegrationPoint, 2> &points)
{
  const auto first = GetTriangleBarycentricCoordinates(points[0]);
  const auto second = GetTriangleBarycentricCoordinates(points[1]);
  constexpr double tolerance = 1024.0 * std::numeric_limits<double>::epsilon();
  int opposite = -1;
  for (int node = 0; node < 3; node++)
  {
    if (std::abs(first[node]) <= tolerance && std::abs(second[node]) <= tolerance)
    {
      if (opposite >= 0)
      {
        throw std::runtime_error(
            "Singular surface quadrature mapped a degenerate boundary segment!");
      }
      opposite = node;
    }
  }
  if (opposite < 0)
  {
    throw std::runtime_error(
        "Singular surface quadrature could not identify the containing triangle edge!");
  }
  std::array<int, 2> nodes;
  int index = 0;
  for (int node = 0; node < 3; node++)
  {
    if (node != opposite)
    {
      nodes[index++] = node;
    }
  }
  const double direction = second[nodes[1]] - first[nodes[1]];
  if (std::abs(direction) <= tolerance)
  {
    throw std::runtime_error(
        "Singular surface quadrature mapped a zero-length boundary segment!");
  }
  if (direction < 0.0)
  {
    std::swap(nodes[0], nodes[1]);
  }
  return nodes;
}

std::array<double, 4>
GetTetrahedronBarycentricCoordinates(const mfem::IntegrationPoint &point)
{
  return {1.0 - point.x - point.y - point.z, point.x, point.y, point.z};
}

int FindTetrahedronVertex(const mfem::IntegrationPoint &point)
{
  const auto lambda = GetTetrahedronBarycentricCoordinates(point);
  const auto maximum = std::max_element(lambda.begin(), lambda.end());
  if (*maximum < 1.0 - 256.0 * std::numeric_limits<double>::epsilon())
  {
    return -1;
  }
  return static_cast<int>(std::distance(lambda.begin(), maximum));
}

std::array<int, 3>
GetContainingTetrahedronFaceNodes(const std::array<mfem::IntegrationPoint, 3> &points)
{
  constexpr double tolerance = 1024.0 * std::numeric_limits<double>::epsilon();
  std::array<std::array<double, 4>, 3> barycentric;
  for (int point = 0; point < 3; point++)
  {
    barycentric[point] = GetTetrahedronBarycentricCoordinates(points[point]);
    for (double coordinate : barycentric[point])
    {
      if (!std::isfinite(coordinate) || coordinate < -tolerance ||
          coordinate > 1.0 + tolerance)
      {
        throw std::runtime_error(
            "Singular surface quadrature mapped a boundary vertex outside its "
            "tetrahedron!");
      }
    }
  }

  int opposite = -1;
  for (int node = 0; node < 4; node++)
  {
    if (std::all_of(barycentric.begin(), barycentric.end(),
                    [node](const auto &lambda)
                    {
                      return std::abs(lambda[node]) <=
                             1024.0 * std::numeric_limits<double>::epsilon();
                    }))
    {
      if (opposite >= 0)
      {
        throw std::runtime_error(
            "Singular surface quadrature mapped a degenerate tetrahedron subface!");
      }
      opposite = node;
    }
  }
  if (opposite < 0)
  {
    throw std::runtime_error(
        "Singular surface quadrature could not identify the containing tetrahedron face!");
  }

  std::array<int, 3> nodes;
  int face_node = 0;
  for (int node = 0; node < 4; node++)
  {
    if (node != opposite)
    {
      nodes[face_node++] = node;
    }
  }
  return nodes;
}

mfem::DenseMatrix
GetBoundaryBernsteinControlPoints(mfem::ElementTransformation &transformation)
{
  if (transformation.GetGeometryType() != mfem::Geometry::SEGMENT ||
      transformation.GetDimension() != 1 || transformation.GetSpaceDim() != 2)
  {
    throw std::invalid_argument(
        "Singular surface cutoff clipping requires a physical planar segment!");
  }

  const int map_order = std::max(1, transformation.Order());
  const int coefficient_count = map_order + 1;
  mfem::DenseMatrix bernstein_values(coefficient_count);
  mfem::DenseMatrix physical_values(coefficient_count, 2);
  mfem::Vector bernstein_shape(coefficient_count);
  mfem::Vector physical_point(2);
  mfem::IntegrationPoint point;
  double coordinate_scale = 1.0;
  for (int q = 0; q < coefficient_count; q++)
  {
    point.x = static_cast<double>(q) / map_order;
    mfem::Poly_1D::CalcBernstein(map_order, point.x, bernstein_shape);
    bernstein_values.SetRow(q, bernstein_shape);
    transformation.Transform(point, physical_point);
    for (int d = 0; d < 2; d++)
    {
      if (!std::isfinite(physical_point[d]))
      {
        throw std::domain_error(
            "Singular surface cutoff clipping found a nonfinite boundary map!");
      }
      physical_values(q, d) = physical_point[d];
      coordinate_scale = std::max(coordinate_scale, std::abs(physical_point[d]));
    }
  }

  mfem::DenseMatrix control_points;
  mfem::DenseMatrixInverse(bernstein_values).Mult(physical_values, control_points);
  if (control_points.Height() != coefficient_count || control_points.Width() != 2)
  {
    throw std::runtime_error(
        "Singular surface cutoff clipping produced invalid Bernstein dimensions!");
  }

  const auto &rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * map_order + 2);
  mfem::Vector reconstructed(2);
  for (int q = -2; q < rule.GetNPoints(); q++)
  {
    point.x = q < 0 ? (q == -2 ? 0.0 : 1.0) : rule.IntPoint(q).x;
    transformation.Transform(point, physical_point);
    mfem::Poly_1D::CalcBernstein(map_order, point.x, bernstein_shape);
    reconstructed = 0.0;
    for (int i = 0; i < coefficient_count; i++)
    {
      for (int d = 0; d < 2; d++)
      {
        reconstructed[d] += bernstein_shape[i] * control_points(i, d);
      }
    }
    reconstructed -= physical_point;
    const double tolerance =
        1.0e-11 * std::max(1.0, physical_point.Norml2()) +
        4096.0 * std::numeric_limits<double>::epsilon() * coordinate_scale;
    if (!std::isfinite(reconstructed.Norml2()) || reconstructed.Norml2() > tolerance)
    {
      throw std::domain_error(
          "Singular surface cutoff clipping could not certify the high-order "
          "boundary map!");
    }
  }
  return control_points;
}

long double BinomialCoefficient(int n, int k)
{
  if (n < 0 || k < 0 || k > n)
  {
    throw std::invalid_argument(
        "Singular surface Bernstein product received invalid polynomial degrees!");
  }
  k = std::min(k, n - k);
  long double value = 1.0L;
  for (int i = 1; i <= k; i++)
  {
    value *= static_cast<long double>(n - k + i) / i;
  }
  return value;
}

std::vector<double> BuildSquaredDistanceBernstein(const mfem::DenseMatrix &control_points,
                                                  const PhysicalSingularTip &tip,
                                                  double cutoff)
{
  const int order = control_points.Height() - 1;
  if (order < 1 || control_points.Width() != 2 || !(cutoff > 0.0))
  {
    throw std::invalid_argument("Singular surface cutoff distance received invalid input!");
  }

  std::vector<double> coefficients(2 * order + 1, -cutoff * cutoff);
  for (int d = 0; d < 2; d++)
  {
    for (int k = 0; k <= 2 * order; k++)
    {
      long double product = 0.0L;
      const int begin = std::max(0, k - order);
      const int end = std::min(order, k);
      for (int i = begin; i <= end; i++)
      {
        const int j = k - i;
        product += BinomialCoefficient(order, i) * BinomialCoefficient(order, j) *
                   (control_points(i, d) - tip.position[d]) *
                   (control_points(j, d) - tip.position[d]);
      }
      coefficients[k] += static_cast<double>(product / BinomialCoefficient(2 * order, k));
    }
  }
  if (!std::all_of(coefficients.begin(), coefficients.end(),
                   [](double value) { return std::isfinite(value); }))
  {
    throw std::domain_error(
        "Singular surface cutoff distance produced nonfinite Bernstein coefficients!");
  }
  return coefficients;
}

void SplitBernsteinCoefficients(const std::vector<double> &coefficients,
                                std::vector<double> &left, std::vector<double> &right)
{
  if (coefficients.size() < 2)
  {
    throw std::invalid_argument(
        "Singular surface cutoff subdivision requires a nonconstant polynomial!");
  }
  const int order = static_cast<int>(coefficients.size()) - 1;
  left.resize(coefficients.size());
  right.resize(coefficients.size());
  auto work = coefficients;
  left[0] = work[0];
  right[order] = work[order];
  for (int level = 1; level <= order; level++)
  {
    for (int i = 0; i <= order - level; i++)
    {
      work[i] = 0.5 * (work[i] + work[i + 1]);
    }
    left[level] = work[0];
    right[order - level] = work[order - level];
  }
}

void AppendCutoffIntervals(const std::vector<double> &distance_squared, double lower,
                           double upper, double tolerance, int depth,
                           std::vector<std::pair<double, double>> &excluded)
{
  const auto [minimum, maximum] =
      std::minmax_element(distance_squared.begin(), distance_squared.end());
  if (*minimum > tolerance)
  {
    return;
  }
  if (*maximum < -tolerance)
  {
    excluded.emplace_back(lower, upper);
    return;
  }
  if (*minimum >= -tolerance && *maximum <= tolerance)
  {
    // The complete subinterval is inside the arithmetic uncertainty band
    // around the cutoff circle. Excluding it is conservative and prevents a
    // nearly tangent root from causing exponential subdivision.
    excluded.emplace_back(lower, upper);
    return;
  }

  constexpr int maximum_depth = 42;
  if (depth == maximum_depth)
  {
    // Conservatively remove the final roundoff-sized root enclosure. This
    // guarantees that no retained quadrature point lies inside the cutoff.
    excluded.emplace_back(lower, upper);
    return;
  }

  std::vector<double> left, right;
  SplitBernsteinCoefficients(distance_squared, left, right);
  const double midpoint = 0.5 * (lower + upper);
  AppendCutoffIntervals(left, lower, midpoint, tolerance, depth + 1, excluded);
  AppendCutoffIntervals(right, midpoint, upper, tolerance, depth + 1, excluded);
}

std::vector<std::pair<double, double>>
GetBoundaryRetainedIntervals(mfem::ElementTransformation &transformation,
                             const std::vector<PhysicalSingularTip> &tips, double cutoff)
{
  if (!(cutoff > 0.0))
  {
    throw std::invalid_argument(
        "Singular surface cutoff clipping requires a positive physical cutoff!");
  }
  if (tips.empty())
  {
    return {{0.0, 1.0}};
  }

  const auto control_points = GetBoundaryBernsteinControlPoints(transformation);
  std::vector<std::pair<double, double>> excluded;
  for (const auto &tip : tips)
  {
    const auto distance_squared =
        BuildSquaredDistanceBernstein(control_points, tip, cutoff);
    double scale = cutoff * cutoff;
    double coordinate_scale = cutoff;
    double distance_scale = cutoff;
    for (int i = 0; i < control_points.Height(); i++)
    {
      coordinate_scale = std::max({coordinate_scale, std::abs(control_points(i, 0)),
                                   std::abs(control_points(i, 1)),
                                   std::abs(tip.position[0]), std::abs(tip.position[1])});
      distance_scale =
          std::max(distance_scale, std::hypot(control_points(i, 0) - tip.position[0],
                                              control_points(i, 1) - tip.position[1]));
    }
    for (double coefficient : distance_squared)
    {
      scale = std::max(scale, std::abs(coefficient));
    }
    const double tolerance = 8192.0 * std::numeric_limits<double>::epsilon() *
                             std::max({scale, coordinate_scale * distance_scale,
                                       std::numeric_limits<double>::min()});
    if (std::all_of(distance_squared.begin(), distance_squared.end(),
                    [=](double value) { return std::abs(value) <= tolerance; }))
    {
      throw std::domain_error(
          "A curved boundary segment follows a singular EdgeCutoff circle; the "
          "retained interval cannot be certified!");
    }
    AppendCutoffIntervals(distance_squared, 0.0, 1.0, tolerance, 0, excluded);
  }
  if (excluded.empty())
  {
    return {{0.0, 1.0}};
  }

  std::sort(excluded.begin(), excluded.end());
  std::vector<std::pair<double, double>> merged;
  constexpr double coordinate_tolerance = 64.0 * std::numeric_limits<double>::epsilon();
  for (const auto &[lower, upper] : excluded)
  {
    if (merged.empty() || lower > merged.back().second + coordinate_tolerance)
    {
      merged.emplace_back(std::max(0.0, lower), std::min(1.0, upper));
    }
    else
    {
      merged.back().second = std::max(merged.back().second, std::min(1.0, upper));
    }
  }

  std::vector<std::pair<double, double>> retained;
  double coordinate = 0.0;
  for (const auto &[lower, upper] : merged)
  {
    if (lower > coordinate + coordinate_tolerance)
    {
      retained.emplace_back(coordinate, lower);
    }
    coordinate = std::max(coordinate, upper);
  }
  if (coordinate < 1.0 - coordinate_tolerance)
  {
    retained.emplace_back(coordinate, 1.0);
  }
  return retained;
}

double GetFaceEdgeCutoffCoordinate(mfem::ElementTransformation &transformation,
                                   const std::array<int, 2> &edge_nodes, double tangent,
                                   double cutoff)
{
  if (!(cutoff > 0.0) || !std::isfinite(tangent) || tangent < 0.0 || tangent > 1.0 ||
      edge_nodes[0] < 0 || edge_nodes[0] >= 3 || edge_nodes[1] < 0 || edge_nodes[1] >= 3 ||
      edge_nodes[0] == edge_nodes[1])
  {
    throw std::invalid_argument(
        "Singular surface face-cutoff mapping received invalid input!");
  }
  const int opposite_node = 3 - edge_nodes[0] - edge_nodes[1];
  const auto reference_point = [&](double radial)
  {
    std::array<double, 3> lambda{};
    lambda[opposite_node] = radial;
    lambda[edge_nodes[0]] = (1.0 - radial) * (1.0 - tangent);
    lambda[edge_nodes[1]] = (1.0 - radial) * tangent;
    mfem::IntegrationPoint point;
    point.Set2(lambda[1], lambda[2]);
    return point;
  };
  const auto vertex_position = [&](int node)
  {
    mfem::IntegrationPoint point;
    point.Set2(node == 1 ? 1.0 : 0.0, node == 2 ? 1.0 : 0.0);
    mfem::Vector position(transformation.GetSpaceDim());
    transformation.Transform(point, position);
    return position;
  };
  const auto edge_start = vertex_position(edge_nodes[0]);
  const auto edge_end = vertex_position(edge_nodes[1]);
  mfem::Vector edge_direction(edge_end);
  edge_direction -= edge_start;
  const double edge_length_squared = edge_direction * edge_direction;
  if (!std::isfinite(edge_length_squared) || !(edge_length_squared > 0.0))
  {
    throw std::domain_error(
        "Singular surface face-cutoff mapping found a degenerate feature edge!");
  }
  const double edge_length = std::sqrt(edge_length_squared);
  const auto distance = [&](double radial)
  {
    mfem::Vector position(transformation.GetSpaceDim());
    const auto point = reference_point(radial);
    transformation.Transform(point, position);
    position -= edge_start;
    const double projection =
        std::clamp((position * edge_direction) / edge_length_squared, 0.0, 1.0);
    position.Add(-projection, edge_direction);
    return position.Norml2();
  };

  const double edge_distance = distance(0.0);
  const double geometry_tolerance =
      4096.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, edge_length);
  if (!std::isfinite(edge_distance) || edge_distance > geometry_tolerance)
  {
    throw std::domain_error(
        "Physical EdgeCutoff on a curved face requires a geometrically straight "
        "singular edge!");
  }
  const double opposite_distance = distance(1.0);
  if (!std::isfinite(opposite_distance) || !(opposite_distance > cutoff))
  {
    throw std::domain_error(
        "Singular surface EdgeCutoff must be smaller than every incident face "
        "altitude!");
  }

  // Find the first intersection of this physical cross-section with the
  // cutoff tube, then verify that the retained ray does not re-enter it.
  constexpr int search_intervals = 64;
  double lower = 0.0;
  double upper = 1.0;
  bool bracketed = false;
  for (int interval = 1; interval <= search_intervals; interval++)
  {
    const double radial = static_cast<double>(interval) / search_intervals;
    const double radial_distance = distance(radial);
    if (!std::isfinite(radial_distance))
    {
      throw std::runtime_error(
          "Singular surface face-cutoff mapping produced a nonfinite distance!");
    }
    if (radial_distance >= cutoff)
    {
      lower = static_cast<double>(interval - 1) / search_intervals;
      upper = radial;
      bracketed = true;
      break;
    }
  }
  if (!bracketed)
  {
    throw std::domain_error(
        "Singular surface EdgeCutoff does not leave a retained curved face region!");
  }
  for (int iteration = 0; iteration < 64; iteration++)
  {
    const double midpoint = 0.5 * (lower + upper);
    if (distance(midpoint) < cutoff)
    {
      lower = midpoint;
    }
    else
    {
      upper = midpoint;
    }
  }
  const double radial_coordinate = 0.5 * (lower + upper);
  for (int interval = 1; interval <= search_intervals; interval++)
  {
    const double radial = radial_coordinate + (1.0 - radial_coordinate) *
                                                  static_cast<double>(interval) /
                                                  search_intervals;
    if (distance(radial) < cutoff * (1.0 - 1.0e-10))
    {
      throw std::domain_error(
          "Physical EdgeCutoff on a curved face requires a star-shaped retained "
          "cross-section!");
    }
  }
  return radial_coordinate;
}

PhysicalFaceEdge GetPhysicalFaceEdge(mfem::ElementTransformation &transformation,
                                     const std::array<int, 2> &edge_nodes)
{
  if (transformation.GetGeometryType() != mfem::Geometry::TRIANGLE ||
      transformation.GetDimension() != 2 || transformation.GetSpaceDim() != 3 ||
      edge_nodes[0] < 0 || edge_nodes[0] >= 3 || edge_nodes[1] < 0 || edge_nodes[1] >= 3 ||
      edge_nodes[0] == edge_nodes[1])
  {
    throw std::invalid_argument(
        "Singular surface physical-edge construction received invalid input!");
  }

  const auto reference_vertex = [](int node)
  {
    mfem::IntegrationPoint point;
    point.Set2(node == 1 ? 1.0 : 0.0, node == 2 ? 1.0 : 0.0);
    return point;
  };
  const int map_order = std::max(1, transformation.Order());
  mfem::H1_SegmentElement segment_element(map_order);
  mfem::IsoparametricTransformation edge_transformation;
  edge_transformation.SetFE(&segment_element);
  mfem::DenseMatrix physical_points(3, segment_element.GetDof());
  mfem::Vector physical_point(3);
  const auto start_point = reference_vertex(edge_nodes[0]);
  const auto end_point = reference_vertex(edge_nodes[1]);
  const auto &segment_nodes = segment_element.GetNodes();
  for (int i = 0; i < segment_nodes.GetNPoints(); i++)
  {
    const double tangent = segment_nodes.IntPoint(i).x;
    mfem::IntegrationPoint point;
    point.Set2((1.0 - tangent) * start_point.x + tangent * end_point.x,
               (1.0 - tangent) * start_point.y + tangent * end_point.y);
    transformation.Transform(point, physical_point);
    physical_points.SetCol(i, physical_point);
  }
  edge_transformation.SetPointMat(physical_points);
  if (!fem::singular::IsGeometricallyStraightSegmentTransformation(edge_transformation))
  {
    throw std::domain_error(
        "Physical EdgeCutoff on a curved face requires geometrically straight "
        "singular edges!");
  }

  PhysicalFaceEdge edge;
  edge.nodes = edge_nodes;
  edge.start.SetSize(3);
  transformation.Transform(start_point, edge.start);
  mfem::Vector end(3);
  transformation.Transform(end_point, end);
  edge.direction.SetSize(3);
  edge.direction = end;
  edge.direction -= edge.start;
  edge.length_squared = edge.direction * edge.direction;
  if (!std::isfinite(edge.length_squared) || !(edge.length_squared > 0.0))
  {
    throw std::domain_error(
        "Singular surface physical-edge construction found a degenerate edge!");
  }
  return edge;
}

FaceCutoffRay GetFaceCutoffRay(mfem::ElementTransformation &transformation,
                               const std::vector<PhysicalFaceEdge> &edges,
                               int boundary_edge, double tangent, double cutoff)
{
  if (edges.size() < 2 || boundary_edge < 0 || boundary_edge >= 3 ||
      !std::isfinite(tangent) || tangent < 0.0 || tangent > 1.0 || !std::isfinite(cutoff) ||
      !(cutoff > 0.0))
  {
    throw std::invalid_argument(
        "Singular surface multi-edge cutoff chart received invalid input!");
  }

  constexpr std::array<double, 3> anchor{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
  std::array<double, 3> boundary{};
  const int first = (boundary_edge + 1) % 3;
  const int second = (boundary_edge + 2) % 3;
  boundary[first] = 1.0 - tangent;
  boundary[second] = tangent;
  const auto minimum_distance = [&](double radius)
  {
    std::array<double, 3> lambda{};
    for (int node = 0; node < 3; node++)
    {
      lambda[node] = (1.0 - radius) * anchor[node] + radius * boundary[node];
    }
    mfem::IntegrationPoint point;
    point.Set2(lambda[1], lambda[2]);
    mfem::Vector position(3);
    transformation.Transform(point, position);
    double distance = std::numeric_limits<double>::infinity();
    int active_edge = -1;
    for (std::size_t edge = 0; edge < edges.size(); edge++)
    {
      const double edge_distance = edges[edge].Distance(position);
      if (edge_distance < distance)
      {
        distance = edge_distance;
        active_edge = static_cast<int>(edge);
      }
    }
    if (!std::isfinite(distance) || active_edge < 0)
    {
      throw std::runtime_error(
          "Singular surface multi-edge cutoff chart produced a nonfinite distance!");
    }
    return std::pair<double, int>{distance, active_edge};
  };

  const double anchor_distance = minimum_distance(0.0).first;
  if (!(anchor_distance > cutoff))
  {
    throw std::domain_error(
        "Physical EdgeCutoff leaves no certified retained curved-face neighborhood "
        "around the barycenter!");
  }

  constexpr int search_intervals = 128;
  double lower = 0.0;
  double upper = 1.0;
  bool bracketed = false;
  for (int interval = 1; interval <= search_intervals; interval++)
  {
    const double radius = static_cast<double>(interval) / search_intervals;
    if (minimum_distance(radius).first <= cutoff)
    {
      lower = static_cast<double>(interval - 1) / search_intervals;
      upper = radius;
      bracketed = true;
      break;
    }
  }
  if (!bracketed)
  {
    return {1.0, -1};
  }
  for (int iteration = 0; iteration < 64; iteration++)
  {
    const double midpoint = 0.5 * (lower + upper);
    if (minimum_distance(midpoint).first > cutoff)
    {
      lower = midpoint;
    }
    else
    {
      upper = midpoint;
    }
  }
  const double radius = 0.5 * (lower + upper);
  const int active_edge = minimum_distance(upper).second;

  const double distance_tolerance =
      1.0e-10 * cutoff + 4096.0 * std::numeric_limits<double>::epsilon() *
                             std::max(cutoff, std::sqrt(edges.front().length_squared));
  for (int interval = 1; interval <= search_intervals; interval++)
  {
    const double outer_radius =
        radius + (1.0 - radius) * static_cast<double>(interval) / search_intervals;
    if (minimum_distance(outer_radius).first > cutoff + distance_tolerance)
    {
      throw std::domain_error(
          "Physical EdgeCutoff on a curved face requires a star-shaped retained "
          "multi-edge region!");
    }
  }
  return {radius, active_edge};
}

bool IsVacuumSide(const MaterialOperator &material, int attribute)
{
  constexpr double threshold = 1.0 - 1.0e-6;
  return material.GetLightSpeedMax(attribute) >= threshold;
}

std::vector<TriangleBoundarySide> GetSelectedSides(const MaterialOperator &material,
                                                   mfem::FaceElementTransformations &face,
                                                   InterfaceDielectric type)
{
  std::vector<TriangleBoundarySide> sides;
  const auto add_side = [&](int element, mfem::ElementTransformation *transformation)
  {
    if (!transformation)
    {
      return;
    }
    const int attribute = transformation->Attribute;
    const bool selected =
        type == InterfaceDielectric::DEFAULT ||
        ((type == InterfaceDielectric::MA || type == InterfaceDielectric::SA) &&
         IsVacuumSide(material, attribute)) ||
        (type == InterfaceDielectric::MS && !IsVacuumSide(material, attribute));
    if (!selected)
    {
      return;
    }
    if (element < 0)
    {
      throw std::runtime_error(
          "Singular surface postprocessing received an invalid adjacent element!");
    }
    sides.push_back({element, attribute, transformation->GetIntPoint()});
  };
  if (face.Elem1)
  {
    add_side(face.Elem1No, face.Elem1);
  }
  if (face.Elem2)
  {
    add_side(face.Elem2No, face.Elem2);
  }
  return sides;
}

std::vector<TetrahedronBoundarySide>
GetSelectedTetrahedronSides(const MaterialOperator &material,
                            mfem::FaceElementTransformations &face,
                            InterfaceDielectric type)
{
  std::vector<TetrahedronBoundarySide> sides;
  const auto add_side = [&](int element, mfem::ElementTransformation *transformation)
  {
    if (!transformation)
    {
      return;
    }
    const int attribute = transformation->Attribute;
    const bool selected =
        type == InterfaceDielectric::DEFAULT ||
        ((type == InterfaceDielectric::MA || type == InterfaceDielectric::SA) &&
         IsVacuumSide(material, attribute)) ||
        (type == InterfaceDielectric::MS && !IsVacuumSide(material, attribute));
    if (!selected)
    {
      return;
    }
    if (element < 0)
    {
      throw std::runtime_error(
          "Tetrahedral singular surface postprocessing received an invalid adjacent "
          "element!");
    }
    sides.push_back({element, attribute, transformation->GetIntPoint()});
  };
  if (face.Elem1)
  {
    add_side(face.Elem1No, face.Elem1);
  }
  if (face.Elem2)
  {
    add_side(face.Elem2No, face.Elem2);
  }
  return sides;
}

double GetTriangleAltitude(mfem::ElementTransformation &transformation, int opposite_node)
{
  if (opposite_node < 0 || opposite_node >= 3)
  {
    throw std::invalid_argument(
        "Singular surface altitude received an invalid triangle node!");
  }
  std::array<mfem::Vector, 3> vertices{mfem::Vector(3), mfem::Vector(3), mfem::Vector(3)};
  for (int node = 0; node < 3; node++)
  {
    mfem::IntegrationPoint point;
    point.Set2(node == 1 ? 1.0 : 0.0, node == 2 ? 1.0 : 0.0);
    transformation.Transform(point, vertices[node]);
  }
  const int edge_node_0 = (opposite_node + 1) % 3;
  const int edge_node_1 = (opposite_node + 2) % 3;
  mfem::Vector edge(vertices[edge_node_1]);
  edge -= vertices[edge_node_0];
  mfem::Vector opposite(vertices[opposite_node]);
  opposite -= vertices[edge_node_0];
  mfem::Vector cross(3);
  cross[0] = edge[1] * opposite[2] - edge[2] * opposite[1];
  cross[1] = edge[2] * opposite[0] - edge[0] * opposite[2];
  cross[2] = edge[0] * opposite[1] - edge[1] * opposite[0];
  const double edge_length = edge.Norml2();
  if (!std::isfinite(edge_length) || !(edge_length > 0.0))
  {
    throw std::runtime_error(
        "Singular surface cutoff mapping found a degenerate triangle edge!");
  }
  const double altitude = cross.Norml2() / edge_length;
  if (!std::isfinite(altitude) || !(altitude > 0.0))
  {
    throw std::runtime_error(
        "Singular surface cutoff mapping found a degenerate triangle!");
  }
  return altitude;
}

using TracePower = std::pair<double, double>;

struct ComplexElectricTraceTerm
{
  fem::singular::Vector2 transverse_real{};
  fem::singular::Vector2 transverse_imaginary{};
  double longitudinal_real = 0.0;
  double longitudinal_imaginary = 0.0;
};

double Dot(const fem::singular::Vector2 &left, const fem::singular::Vector2 &right)
{
  return left[0] * right[0] + left[1] * right[1];
}

double Dot(const fem::singular::Vector3 &left, const fem::singular::Vector3 &right)
{
  return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

void AddVectorTerms(std::map<TracePower, ComplexElectricTraceTerm> &result,
                    const std::vector<fem::singular::TriangleVectorTraceTerm> &terms,
                    bool imaginary, double scale)
{
  for (const auto &term : terms)
  {
    auto &value = result[{term.exponents.left, term.exponents.right}];
    auto &field = imaginary ? value.transverse_imaginary : value.transverse_real;
    field[0] += scale * term.coefficient[0];
    field[1] += scale * term.coefficient[1];
  }
}

void AddScalarTerms(std::map<TracePower, ComplexElectricTraceTerm> &result,
                    const std::vector<fem::singular::TriangleScalarTraceTerm> &terms,
                    bool imaginary, double scale)
{
  for (const auto &term : terms)
  {
    auto &value = result[{term.exponents.left, term.exponents.right}];
    (imaginary ? value.longitudinal_imaginary : value.longitudinal_real) +=
        scale * term.coefficient;
  }
}

double IntegratePowerExpansion(
    const std::function<std::map<TracePower, double>(double)> &coefficient_expansion,
    double lower, double upper, const fem::singular::AdaptiveAssemblyOptions &options)
{
  if (options.quadrature_order < 1 || !std::isfinite(options.absolute_tolerance) ||
      options.absolute_tolerance < 0.0 || !std::isfinite(options.relative_tolerance) ||
      options.relative_tolerance < 0.0 ||
      !(options.absolute_tolerance > 0.0 || options.relative_tolerance > 0.0) ||
      options.maximum_subdivisions < 1)
  {
    throw std::invalid_argument(
        "Singular surface postprocessing received invalid adaptive quadrature "
        "options!");
  }
  if (!std::isfinite(lower) || !std::isfinite(upper) || lower < 0.0 || upper > 1.0 ||
      !(lower < upper))
  {
    throw std::invalid_argument(
        "Singular surface postprocessing received an invalid integration interval!");
  }

  std::set<TracePower> powers;
  for (double fraction : {0.25, 0.5, 0.75})
  {
    const double coordinate = lower + fraction * (upper - lower);
    for (const auto &[power, coefficient] : coefficient_expansion(coordinate))
    {
      if (coefficient != 0.0)
      {
        powers.insert(power);
      }
    }
  }
  for (const auto &[left, right] : powers)
  {
    if ((lower == 0.0 && !(left > -1.0)) || (upper == 1.0 && !(right > -1.0)))
    {
      throw std::domain_error(
          "Singular surface participation is nonintegrable at an endpoint with "
          "nu <= 1/2; use an explicit physical cutoff or surface-response model!");
    }
  }

  const auto integrate = [&](int order)
  {
    long double value = 0.0L;
    for (const auto &power : powers)
    {
      const double left = power.first;
      const double right = power.second;
      const auto integrate_mapped =
          [&](double local_left, double local_right, const auto &mapping)
      {
        const auto rule =
            fem::singular::BuildWeightedSegmentQuadrature(order, local_left, local_right);
        for (const auto &quadrature : rule)
        {
          const auto [coordinate, jacobian] = mapping(quadrature.coordinate);
          const auto coefficients = coefficient_expansion(coordinate);
          const auto coefficient = coefficients.find({left, right});
          if (coefficient != coefficients.end())
          {
            const double global_weight =
                std::pow(coordinate, left) * std::pow(1.0 - coordinate, right);
            const double local_weight = std::pow(quadrature.coordinate, local_left) *
                                        std::pow(1.0 - quadrature.coordinate, local_right);
            const double contribution = quadrature.weight * jacobian * coefficient->second *
                                        global_weight / local_weight;
            if (!std::isfinite(contribution))
            {
              throw std::runtime_error(
                  "Singular surface power quadrature produced a nonfinite "
                  "contribution!");
            }
            value += contribution;
          }
        }
      };
      const auto integrate_affine = [&](double begin, double end)
      {
        const double length = end - begin;
        const double local_left = begin == 0.0 ? left : 0.0;
        const double local_right = end == 1.0 ? right : 0.0;
        integrate_mapped(local_left, local_right, [=](double coordinate)
                         { return std::pair{begin + length * coordinate, length}; });
      };
      const auto integrate_log_left = [&](double begin, double end)
      {
        const double logarithmic_length = std::log(end / begin);
        const double local_right = end == 1.0 ? right : 0.0;
        integrate_mapped(0.0, local_right,
                         [=](double coordinate)
                         {
                           const double physical_coordinate =
                               begin * std::exp(logarithmic_length * coordinate);
                           return std::pair{physical_coordinate,
                                            logarithmic_length * physical_coordinate};
                         });
      };
      const auto integrate_log_right = [&](double begin, double end)
      {
        const double begin_distance = 1.0 - begin;
        const double end_distance = 1.0 - end;
        const double logarithmic_ratio = std::log(end_distance / begin_distance);
        const double local_left = begin == 0.0 ? left : 0.0;
        integrate_mapped(local_left, 0.0,
                         [=](double coordinate)
                         {
                           const double distance =
                               begin_distance * std::exp(logarithmic_ratio * coordinate);
                           return std::pair{1.0 - distance, -logarithmic_ratio * distance};
                         });
      };

      const bool logarithmic_left = lower > 0.0 && left < 0.0;
      const bool logarithmic_right = upper < 1.0 && right < 0.0;
      if (logarithmic_left && logarithmic_right)
      {
        const double midpoint = 0.5 * (lower + upper);
        integrate_log_left(lower, midpoint);
        integrate_log_right(midpoint, upper);
      }
      else if (logarithmic_left)
      {
        integrate_log_left(lower, upper);
      }
      else if (logarithmic_right)
      {
        integrate_log_right(lower, upper);
      }
      else
      {
        integrate_affine(lower, upper);
      }
    }
    return static_cast<double>(value);
  };

  constexpr int order_increment = 8;
  int comparison_order = std::max(8, 2 * options.quadrature_order);
  double comparison = integrate(comparison_order);
  double value = comparison;
  double error = std::numeric_limits<double>::infinity();
  double tolerance = 0.0;
  int order = comparison_order + order_increment;
  for (int refinement = 0; refinement <= options.maximum_subdivisions; refinement++)
  {
    value = integrate(order);
    const double scale = std::max({1.0, std::abs(value), std::abs(comparison)});
    error = 8.0 * std::abs(value - comparison) +
            128.0 * std::numeric_limits<double>::epsilon() * scale;
    tolerance = options.absolute_tolerance + options.relative_tolerance * std::abs(value);
    if (std::isfinite(value) && std::isfinite(comparison) && std::isfinite(error) &&
        error <= tolerance)
    {
      return value;
    }
    if (refinement < options.maximum_subdivisions)
    {
      comparison = value;
      comparison_order = order;
      order += order_increment;
    }
  }
  throw std::runtime_error(fmt::format(
      "Singular surface power quadrature did not meet tolerance: value = {:.17g}, "
      "comparison = {:.17g}, estimated absolute error = {:.17g}, tolerance = {:.17g}, "
      "orders = {}/{}, order refinements = {}!",
      value, comparison, error, tolerance, order, comparison_order,
      options.maximum_subdivisions));
}

}  // namespace

TriangleSingularSurfacePostOperator::TriangleSingularSurfacePostOperator(
    const config::BoundaryPostData &postpro, const MaterialOperator &material,
    mfem::ParFiniteElementSpace &fespace)
  : material(material), fespace(fespace)
{
  auto *mesh = fespace.GetParMesh();
  if (!mesh || mesh->Dimension() != 2 || mesh->SpaceDimension() != 2 ||
      (fespace.FEColl()->GetMapType(2) != mfem::FiniteElement::H_CURL &&
       fespace.FEColl()->GetMapType(2) != mfem::FiniteElement::VALUE))
  {
    throw std::invalid_argument(
        "Triangular singular surface postprocessing requires a two-dimensional "
        "H(curl) or H1 space!");
  }
  const int maximum_attribute =
      mesh->bdr_attributes.Size() ? mesh->bdr_attributes.Max() : 0;
  for (const auto &[index, data] : postpro.dielectric)
  {
    InterfaceData interface{data.type,        data.t,
                            data.epsilon_r,   data.tandelta,
                            data.edge_cutoff, mfem::Array<int>(maximum_attribute)};
    interface.attribute_marker = 0;
    for (int attribute : data.attributes)
    {
      if (attribute <= 0 || attribute > maximum_attribute)
      {
        throw std::invalid_argument(
            "Singular surface postprocessing received an unknown boundary attribute!");
      }
      interface.attribute_marker[attribute - 1] = 1;
    }
    if (!std::isfinite(interface.thickness) || !(interface.thickness > 0.0) ||
        !std::isfinite(interface.permittivity) || !(interface.permittivity > 0.0) ||
        !std::isfinite(interface.loss_tangent) || interface.loss_tangent < 0.0 ||
        !std::isfinite(interface.edge_cutoff) || interface.edge_cutoff < 0.0)
    {
      throw std::invalid_argument(
          "Singular surface postprocessing received invalid interface material data!");
    }
    interfaces.emplace(index, std::move(interface));
  }
}

double TriangleSingularSurfacePostOperator::IntegrateInterface(
    const InterfaceData &interface,
    fem::singular::TriangleEnrichedNDFieldEvaluator *real_evaluator,
    fem::singular::TriangleEnrichedNDFieldEvaluator *imaginary_evaluator,
    fem::singular::TriangleEnrichedH1FieldEvaluator *real_gradient_evaluator,
    fem::singular::TriangleEnrichedH1FieldEvaluator *imaginary_gradient_evaluator,
    fem::singular::TriangleEnrichedH1FieldEvaluator *real_longitudinal_evaluator,
    fem::singular::TriangleEnrichedH1FieldEvaluator *imaginary_longitudinal_evaluator,
    const fem::singular::AdaptiveAssemblyOptions &options) const
{
  auto &mesh = *fespace.GetParMesh();
  mfem::FaceElementTransformations face;
  mfem::IsoparametricTransformation element1, element2;
  const auto node_exponent = [&](int element, int node)
  {
    return real_evaluator
               ? real_evaluator->GetElementNodeSingularExponent(element, node)
               : real_gradient_evaluator->GetElementNodeSingularExponent(element, node);
  };

  // A physical cutoff belongs to a singular point, not to one mesh boundary
  // segment. Gather the complete tip set so refinement and repartitioning can
  // place arbitrarily many local segments inside one cutoff neighborhood.
  std::vector<double> local_tip_data;
  if (interface.edge_cutoff > 0.0)
  {
    for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
    {
      const int attribute = mesh.GetBdrAttribute(boundary);
      if (attribute <= 0 || attribute > interface.attribute_marker.Size() ||
          !interface.attribute_marker[attribute - 1])
      {
        continue;
      }
      auto *boundary_transformation = mesh.GetBdrElementTransformation(boundary);
      if (!boundary_transformation ||
          boundary_transformation->GetGeometryType() != mfem::Geometry::SEGMENT)
      {
        throw std::runtime_error(
            "Triangular singular surface postprocessing requires segment boundaries!");
      }
      for (int endpoint = 0; endpoint < 2; endpoint++)
      {
        mfem::IntegrationPoint point;
        point.x = static_cast<double>(endpoint);
        BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
            boundary, mesh, face, element1, element2, &point);
        const auto sides = GetSelectedSides(material, face, interface.type);
        double exponent = 1.0;
        for (const auto &side : sides)
        {
          const int node = FindTriangleVertex(side.point);
          if (node >= 0)
          {
            exponent = std::min(exponent, node_exponent(side.element, node));
          }
        }
        if (exponent < 1.0)
        {
          mfem::Vector position(2);
          boundary_transformation->Transform(point, position);
          if (!std::isfinite(position[0]) || !std::isfinite(position[1]) ||
              !std::isfinite(exponent) || !(exponent > 0.0))
          {
            throw std::domain_error(
                "Singular surface cutoff found invalid physical tip data!");
          }
          local_tip_data.insert(local_tip_data.end(), {position[0], position[1], exponent});
        }
      }
    }
  }

  std::vector<int> tip_counts(Mpi::Size(mesh.GetComm()));
  const int local_tip_count = static_cast<int>(local_tip_data.size());
  Mpi::Allgather(1, &local_tip_count, tip_counts.data(), mesh.GetComm());
  std::vector<int> tip_offsets(tip_counts.size());
  std::partial_sum(tip_counts.begin(), tip_counts.end() - 1, tip_offsets.begin() + 1);
  const int global_tip_count = std::accumulate(tip_counts.begin(), tip_counts.end(), 0);
  if (local_tip_count % 3 != 0 || global_tip_count % 3 != 0)
  {
    throw std::runtime_error(
        "Singular surface cutoff MPI tip records have invalid dimensions!");
  }
  std::vector<double> global_tip_data(global_tip_count);
  Mpi::Allgatherv(local_tip_count, local_tip_data.data(), global_tip_data.data(),
                  tip_counts.data(), tip_offsets.data(), mesh.GetComm());

  std::vector<PhysicalSingularTip> singular_tips;
  double tip_coordinate_scale = std::max(1.0, interface.edge_cutoff);
  for (int record = 0; record < global_tip_count; record += 3)
  {
    PhysicalSingularTip tip{{global_tip_data[record], global_tip_data[record + 1]},
                            global_tip_data[record + 2]};
    if (!std::isfinite(tip.position[0]) || !std::isfinite(tip.position[1]) ||
        !std::isfinite(tip.nu) || !(tip.nu > 0.0) || !(tip.nu < 1.0))
    {
      throw std::domain_error(
          "Singular surface cutoff received invalid gathered tip data!");
    }
    tip_coordinate_scale = std::max(
        {tip_coordinate_scale, std::abs(tip.position[0]), std::abs(tip.position[1])});
    singular_tips.push_back(tip);
  }
  std::sort(singular_tips.begin(), singular_tips.end(),
            [](const auto &left, const auto &right)
            {
              return left.position < right.position ||
                     (left.position == right.position && left.nu < right.nu);
            });
  const double duplicate_tolerance =
      4096.0 * std::numeric_limits<double>::epsilon() * tip_coordinate_scale;
  std::vector<PhysicalSingularTip> unique_singular_tips;
  for (const auto &tip : singular_tips)
  {
    if (!unique_singular_tips.empty())
    {
      const double dx = tip.position[0] - unique_singular_tips.back().position[0];
      const double dy = tip.position[1] - unique_singular_tips.back().position[1];
      if (std::hypot(dx, dy) <= duplicate_tolerance)
      {
        unique_singular_tips.back().nu = std::min(unique_singular_tips.back().nu, tip.nu);
        continue;
      }
    }
    unique_singular_tips.push_back(tip);
  }

  long double local_energy = 0.0L;
  for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
  {
    const int attribute = mesh.GetBdrAttribute(boundary);
    if (attribute <= 0 || attribute > interface.attribute_marker.Size() ||
        !interface.attribute_marker[attribute - 1])
    {
      continue;
    }
    auto *boundary_transformation = mesh.GetBdrElementTransformation(boundary);
    if (!boundary_transformation ||
        boundary_transformation->GetGeometryType() != mfem::Geometry::SEGMENT)
    {
      throw std::runtime_error(
          "Triangular singular surface postprocessing requires segment boundaries!");
    }

    std::map<int, std::array<mfem::IntegrationPoint, 2>> endpoint_points;
    std::map<int, std::array<bool, 2>> endpoint_seen;
    std::array<double, 2> endpoint_exponents{1.0, 1.0};
    bool has_selected_side = false;
    for (int endpoint = 0; endpoint < 2; endpoint++)
    {
      mfem::IntegrationPoint point;
      point.x = static_cast<double>(endpoint);
      BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
          boundary, mesh, face, element1, element2, &point);
      const auto sides = GetSelectedSides(material, face, interface.type);
      has_selected_side = has_selected_side || !sides.empty();
      for (const auto &side : sides)
      {
        auto [points, inserted] =
            endpoint_points.emplace(side.element, std::array<mfem::IntegrationPoint, 2>{});
        auto [seen, seen_inserted] =
            endpoint_seen.emplace(side.element, std::array<bool, 2>{false, false});
        MFEM_VERIFY(inserted == seen_inserted,
                    "Singular surface endpoint bookkeeping is inconsistent!");
        if (seen->second[endpoint] && (std::abs(points->second[endpoint].x - side.point.x) >
                                           256.0 * std::numeric_limits<double>::epsilon() ||
                                       std::abs(points->second[endpoint].y - side.point.y) >
                                           256.0 * std::numeric_limits<double>::epsilon()))
        {
          throw std::runtime_error(
              "Singular surface postprocessing found inconsistent boundary endpoint "
              "orientation!");
        }
        points->second[endpoint] = side.point;
        seen->second[endpoint] = true;
        const int node = FindTriangleVertex(side.point);
        if (node >= 0)
        {
          endpoint_exponents[endpoint] =
              std::min(endpoint_exponents[endpoint], node_exponent(side.element, node));
        }
      }
    }
    if (!has_selected_side)
    {
      continue;
    }
    if (interface.edge_cutoff == 0.0 &&
        (endpoint_exponents[0] <= 0.5 || endpoint_exponents[1] <= 0.5))
    {
      throw std::domain_error(
          "Singular surface participation is nonintegrable at an ideal thin-sheet "
          "endpoint with nu <= 1/2; specify EdgeCutoff or use a surface-response "
          "model!");
    }

    std::map<int, std::array<int, 2>> endpoint_nodes;
    for (const auto &[element, points] : endpoint_points)
    {
      const auto seen = endpoint_seen.find(element);
      if (seen == endpoint_seen.end() || !seen->second[0] || !seen->second[1])
      {
        throw std::runtime_error(
            "Singular surface postprocessing could not map both boundary endpoints into "
            "an adjacent triangle!");
      }
      endpoint_nodes.emplace(element, GetTriangleEdgeNodes(points));
    }
    const auto retained_intervals =
        interface.edge_cutoff > 0.0
            ? GetBoundaryRetainedIntervals(*boundary_transformation, unique_singular_tips,
                                           interface.edge_cutoff)
            : std::vector<std::pair<double, double>>{{0.0, 1.0}};

    const auto coefficient_expansion = [&](double coordinate)
    {
      mfem::IntegrationPoint point;
      point.x = coordinate;
      const bool inverted =
          BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
              boundary, mesh, face, element1, element2, &point);
      auto sides = GetSelectedSides(material, face, interface.type);
      if (sides.empty())
      {
        return std::map<TracePower, double>{};
      }
      for (auto &side : sides)
      {
        const auto endpoint = endpoint_nodes.find(side.element);
        if (endpoint == endpoint_nodes.end())
        {
          throw std::runtime_error(
              "Singular surface postprocessing changed adjacent elements along one "
              "boundary segment!");
        }
        side.endpoint_nodes = endpoint->second;
      }

      boundary_transformation->SetIntPoint(&point);
      mfem::Vector normal_vector(2);
      BdrGridFunctionCoefficient::GetNormal(*boundary_transformation, normal_vector,
                                            inverted);
      const fem::singular::Vector2 normal{normal_vector[0], normal_vector[1]};
      std::map<TracePower, ComplexElectricTraceTerm> field_terms;
      const double average_scale = 1.0 / static_cast<double>(sides.size());
      for (const auto &side : sides)
      {
        if (real_evaluator)
        {
          AddVectorTerms(field_terms,
                         real_evaluator->EvaluateTraceExpansion(side.element, side.point,
                                                                side.endpoint_nodes),
                         false, average_scale);
          AddVectorTerms(field_terms,
                         imaginary_evaluator->EvaluateTraceExpansion(
                             side.element, side.point, side.endpoint_nodes),
                         true, average_scale);
        }
        else
        {
          AddVectorTerms(field_terms,
                         real_gradient_evaluator->EvaluateGradientTraceExpansion(
                             side.element, side.point, side.endpoint_nodes),
                         false, -average_scale);
          AddVectorTerms(field_terms,
                         imaginary_gradient_evaluator->EvaluateGradientTraceExpansion(
                             side.element, side.point, side.endpoint_nodes),
                         true, -average_scale);
        }
        if (real_longitudinal_evaluator)
        {
          AddScalarTerms(field_terms,
                         real_longitudinal_evaluator->EvaluateValueTraceExpansion(
                             side.element, side.point, side.endpoint_nodes),
                         false, average_scale);
          AddScalarTerms(field_terms,
                         imaginary_longitudinal_evaluator->EvaluateValueTraceExpansion(
                             side.element, side.point, side.endpoint_nodes),
                         true, average_scale);
        }
      }

      std::map<TracePower, double> density_terms;
      const auto displacement_normal = [&](const fem::singular::Vector2 &field)
      {
        fem::singular::Vector2 displacement{};
        mfem::Vector field_vector(const_cast<double *>(field.data()), 2);
        mfem::Vector displacement_vector(displacement.data(), 2);
        material.GetPermittivityReal(sides.front().attribute)
            .Mult(field_vector, displacement_vector);
        return Dot(displacement, normal);
      };
      for (const auto &[left_power, left] : field_terms)
      {
        for (const auto &[right_power, right] : field_terms)
        {
          const double real_normal =
              Dot(left.transverse_real, normal) * Dot(right.transverse_real, normal);
          const double imaginary_normal = Dot(left.transverse_imaginary, normal) *
                                          Dot(right.transverse_imaginary, normal);
          const double normal_product = real_normal + imaginary_normal;
          const double transverse_product =
              Dot(left.transverse_real, right.transverse_real) +
              Dot(left.transverse_imaginary, right.transverse_imaginary);
          const double longitudinal_product =
              left.longitudinal_real * right.longitudinal_real +
              left.longitudinal_imaginary * right.longitudinal_imaginary;
          double density = 0.0;
          if (interface.type == InterfaceDielectric::DEFAULT)
          {
            density = 0.5 * interface.thickness * interface.permittivity *
                      (transverse_product + longitudinal_product);
          }
          else if (interface.type == InterfaceDielectric::MA)
          {
            density = 0.5 * (interface.thickness / interface.permittivity) * normal_product;
          }
          else if (interface.type == InterfaceDielectric::MS)
          {
            const double displacement_product =
                displacement_normal(left.transverse_real) *
                    displacement_normal(right.transverse_real) +
                displacement_normal(left.transverse_imaginary) *
                    displacement_normal(right.transverse_imaginary);
            density =
                0.5 * (interface.thickness / interface.permittivity) * displacement_product;
          }
          else
          {
            const double tangential_product =
                transverse_product - normal_product + longitudinal_product;
            density = 0.5 * interface.thickness *
                      (interface.permittivity * tangential_product +
                       normal_product / interface.permittivity);
          }
          const TracePower power{left_power.first + right_power.first,
                                 left_power.second + right_power.second};
          density_terms[power] += boundary_transformation->Weight() * density;
        }
      }
      return density_terms;
    };
    for (const auto &[lower, upper] : retained_intervals)
    {
      local_energy += IntegratePowerExpansion(coefficient_expansion, lower, upper, options);
    }
  }
  double energy = static_cast<double>(local_energy);
  Mpi::GlobalSum(1, &energy, mesh.GetComm());
  return energy;
}

std::vector<TriangleSingularSurfacePostOperator::Measurement>
TriangleSingularSurfacePostOperator::Measure(
    fem::singular::TriangleEnrichedNDFieldEvaluator &real_evaluator,
    fem::singular::TriangleEnrichedNDFieldEvaluator &imaginary_evaluator,
    double total_electric_energy, const fem::singular::AdaptiveAssemblyOptions &options,
    fem::singular::TriangleEnrichedH1FieldEvaluator *real_longitudinal_evaluator,
    fem::singular::TriangleEnrichedH1FieldEvaluator *imaginary_longitudinal_evaluator) const
{
  if (!std::isfinite(total_electric_energy) || !(total_electric_energy > 0.0))
  {
    throw std::invalid_argument(
        "Singular surface participation requires positive total electric energy!");
  }
  if ((real_longitudinal_evaluator == nullptr) !=
      (imaginary_longitudinal_evaluator == nullptr))
  {
    throw std::invalid_argument(
        "Singular surface participation requires both real and imaginary longitudinal "
        "field evaluators!");
  }
  std::vector<Measurement> measurements;
  measurements.reserve(interfaces.size());
  for (const auto &[index, interface] : interfaces)
  {
    const double energy = IntegrateInterface(
        interface, &real_evaluator, &imaginary_evaluator, nullptr, nullptr,
        real_longitudinal_evaluator, imaginary_longitudinal_evaluator, options);
    const double participation = energy / total_electric_energy;
    const double quality_factor = participation == 0.0 || interface.loss_tangent == 0.0
                                      ? mfem::infinity()
                                      : 1.0 / (participation * interface.loss_tangent);
    measurements.push_back(
        {index, energy, interface.loss_tangent, participation, quality_factor});
  }
  return measurements;
}

std::vector<TriangleSingularSurfacePostOperator::Measurement>
TriangleSingularSurfacePostOperator::MeasureElectrostatic(
    fem::singular::TriangleEnrichedH1FieldEvaluator &real_evaluator,
    fem::singular::TriangleEnrichedH1FieldEvaluator &imaginary_evaluator,
    double total_electric_energy,
    const fem::singular::AdaptiveAssemblyOptions &options) const
{
  if (!std::isfinite(total_electric_energy) || !(total_electric_energy > 0.0))
  {
    throw std::invalid_argument(
        "Singular surface participation requires positive total electric energy!");
  }
  std::vector<Measurement> measurements;
  measurements.reserve(interfaces.size());
  for (const auto &[index, interface] : interfaces)
  {
    const double energy =
        IntegrateInterface(interface, nullptr, nullptr, &real_evaluator,
                           &imaginary_evaluator, nullptr, nullptr, options);
    const double participation = energy / total_electric_energy;
    const double quality_factor = participation == 0.0 || interface.loss_tangent == 0.0
                                      ? mfem::infinity()
                                      : 1.0 / (participation * interface.loss_tangent);
    measurements.push_back(
        {index, energy, interface.loss_tangent, participation, quality_factor});
  }
  return measurements;
}

TetrahedronSingularSurfacePostOperator::TetrahedronSingularSurfacePostOperator(
    const config::BoundaryPostData &postpro, const MaterialOperator &material,
    mfem::ParFiniteElementSpace &fespace)
  : material(material), fespace(fespace)
{
  auto *mesh = fespace.GetParMesh();
  if (!mesh || mesh->Dimension() != 3 || mesh->SpaceDimension() != 3 ||
      (fespace.FEColl()->GetMapType(3) != mfem::FiniteElement::VALUE &&
       fespace.FEColl()->GetMapType(3) != mfem::FiniteElement::H_CURL))
  {
    throw std::invalid_argument(
        "Tetrahedral singular surface postprocessing requires a three-dimensional H1 "
        "or H(curl) space!");
  }
  const int maximum_attribute =
      mesh->bdr_attributes.Size() ? mesh->bdr_attributes.Max() : 0;
  for (const auto &[index, data] : postpro.dielectric)
  {
    InterfaceData interface{data.type,        data.t,
                            data.epsilon_r,   data.tandelta,
                            data.edge_cutoff, mfem::Array<int>(maximum_attribute)};
    interface.attribute_marker = 0;
    for (int attribute : data.attributes)
    {
      if (attribute <= 0 || attribute > maximum_attribute)
      {
        throw std::invalid_argument(
            "Tetrahedral singular surface postprocessing received an unknown boundary "
            "attribute!");
      }
      interface.attribute_marker[attribute - 1] = 1;
    }
    if (!std::isfinite(interface.thickness) || !(interface.thickness > 0.0) ||
        !std::isfinite(interface.permittivity) || !(interface.permittivity > 0.0) ||
        !std::isfinite(interface.loss_tangent) || interface.loss_tangent < 0.0 ||
        !std::isfinite(interface.edge_cutoff) || interface.edge_cutoff < 0.0)
    {
      throw std::invalid_argument(
          "Tetrahedral singular surface postprocessing received invalid interface "
          "material data!");
    }
    interfaces.emplace(index, std::move(interface));
  }
}

std::vector<std::vector<TetrahedronSingularSurfacePostOperator::IndexedInterface>>
TetrahedronSingularSurfacePostOperator::GetInterfaceGroups() const
{
  const auto same_geometry = [](const InterfaceData &left, const InterfaceData &right)
  {
    if (left.edge_cutoff != right.edge_cutoff ||
        left.attribute_marker.Size() != right.attribute_marker.Size())
    {
      return false;
    }
    for (int attribute = 0; attribute < left.attribute_marker.Size(); attribute++)
    {
      if (left.attribute_marker[attribute] != right.attribute_marker[attribute])
      {
        return false;
      }
    }
    return true;
  };

  std::vector<std::vector<IndexedInterface>> groups;
  for (const auto &[index, interface] : interfaces)
  {
    const auto *current = &interface;
    const auto group = std::find_if(
        groups.begin(), groups.end(), [current, &same_geometry](const auto &candidate)
        { return same_geometry(*current, *candidate.front().second); });
    if (group == groups.end())
    {
      groups.push_back({{index, &interface}});
    }
    else
    {
      group->push_back({index, &interface});
    }
  }
  return groups;
}

std::vector<double> TetrahedronSingularSurfacePostOperator::IntegrateInterfaces(
    const std::vector<const InterfaceData *> &interfaces,
    const std::vector<NDFieldEvaluatorPair> &field_evaluators,
    fem::singular::EnrichedH1FieldEvaluator *real_gradient_evaluator,
    fem::singular::EnrichedH1FieldEvaluator *imaginary_gradient_evaluator,
    const fem::singular::AdaptiveAssemblyOptions &options) const
{
  if (interfaces.empty() ||
      std::any_of(interfaces.begin(), interfaces.end(),
                  [](const auto *interface) { return interface == nullptr; }))
  {
    throw std::invalid_argument(
        "Tetrahedral singular surface postprocessing requires at least one interface!");
  }
  const bool full_wave = !field_evaluators.empty();
  const bool electrostatic =
      real_gradient_evaluator != nullptr && imaginary_gradient_evaluator != nullptr;
  if ((real_gradient_evaluator == nullptr) != (imaginary_gradient_evaluator == nullptr) ||
      full_wave == electrostatic)
  {
    throw std::invalid_argument(
        "Tetrahedral singular surface postprocessing requires exactly one complete "
        "full-wave batch or electrostatic evaluator pair!");
  }
  if (full_wave &&
      std::any_of(field_evaluators.begin(), field_evaluators.end(),
                  [](const auto &evaluators)
                  { return evaluators.first == nullptr || evaluators.second == nullptr; }))
  {
    throw std::invalid_argument(
        "Tetrahedral singular surface postprocessing received a null full-wave "
        "evaluator!");
  }
  const bool valid_fixed = options.fixed_subdivision && options.quadrature_order >= 1 &&
                           options.subdivisions >= 0 && options.subdivisions <= 8;
  const bool valid_adaptive =
      !options.fixed_subdivision && options.quadrature_order >= 1 &&
      std::isfinite(options.absolute_tolerance) && options.absolute_tolerance >= 0.0 &&
      std::isfinite(options.relative_tolerance) && options.relative_tolerance >= 0.0 &&
      (options.absolute_tolerance > 0.0 || options.relative_tolerance > 0.0) &&
      options.maximum_subdivisions >= 1;
  if (!valid_fixed && !valid_adaptive)
  {
    throw std::invalid_argument(
        "Tetrahedral singular surface postprocessing received invalid quadrature "
        "options!");
  }

  const auto &geometry = *interfaces.front();
  if (options.fixed_subdivision && geometry.edge_cutoff > 0.0)
  {
    throw std::invalid_argument(
        "FixedSubdivision singular surface postprocessing currently requires "
        "EdgeCutoff = 0!");
  }
  bool requires_vacuum = false;
  bool requires_substrate = false;
  for (const auto *interface : interfaces)
  {
    if (interface->edge_cutoff != geometry.edge_cutoff ||
        interface->attribute_marker.Size() != geometry.attribute_marker.Size())
    {
      throw std::invalid_argument(
          "Tetrahedral singular surface integration group has inconsistent geometry!");
    }
    for (int attribute = 0; attribute < geometry.attribute_marker.Size(); attribute++)
    {
      if (interface->attribute_marker[attribute] != geometry.attribute_marker[attribute])
      {
        throw std::invalid_argument(
            "Tetrahedral singular surface integration group has inconsistent "
            "attributes!");
      }
    }
    requires_vacuum = requires_vacuum || interface->type == InterfaceDielectric::DEFAULT ||
                      interface->type == InterfaceDielectric::MA ||
                      interface->type == InterfaceDielectric::SA;
    requires_substrate = requires_substrate ||
                         interface->type == InterfaceDielectric::DEFAULT ||
                         interface->type == InterfaceDielectric::MS;
  }

  auto &mesh = *fespace.GetParMesh();
  mfem::FaceElementTransformations face;
  mfem::IsoparametricTransformation element1, element2;
  const std::size_t mode_count = full_wave ? field_evaluators.size() : 1;
  const std::size_t component_count = mode_count * interfaces.size();
  std::vector<long double> local_energy(component_count, 0.0L);
  for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
  {
    const int attribute = mesh.GetBdrAttribute(boundary);
    if (attribute <= 0 || attribute > geometry.attribute_marker.Size() ||
        !geometry.attribute_marker[attribute - 1])
    {
      continue;
    }
    auto *boundary_transformation = mesh.GetBdrElementTransformation(boundary);
    if (!boundary_transformation ||
        boundary_transformation->GetGeometryType() != mfem::Geometry::TRIANGLE)
    {
      throw std::runtime_error(
          "Tetrahedral singular surface postprocessing requires triangle boundaries!");
    }

    std::map<int, std::array<mfem::IntegrationPoint, 3>> element_face_points;
    std::map<int, std::array<bool, 3>> element_face_points_seen;
    std::map<int, int> element_attributes;
    bool has_selected_side = false;
    for (int face_node = 0; face_node < 3; face_node++)
    {
      mfem::IntegrationPoint point;
      point.Set2(face_node == 1 ? 1.0 : 0.0, face_node == 2 ? 1.0 : 0.0);
      BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
          boundary, mesh, face, element1, element2, &point);
      const auto sides =
          GetSelectedTetrahedronSides(material, face, InterfaceDielectric::DEFAULT);
      for (const auto &side : sides)
      {
        const bool vacuum = IsVacuumSide(material, side.attribute);
        if ((vacuum && !requires_vacuum) || (!vacuum && !requires_substrate))
        {
          continue;
        }
        has_selected_side = true;
        auto [entry, inserted] = element_face_points.emplace(
            side.element, std::array<mfem::IntegrationPoint, 3>{});
        auto [seen, seen_inserted] = element_face_points_seen.emplace(
            side.element, std::array<bool, 3>{false, false, false});
        if (inserted != seen_inserted)
        {
          throw std::logic_error(
              "Tetrahedral singular surface face-point bookkeeping is inconsistent!");
        }
        auto [attribute_entry, attribute_inserted] =
            element_attributes.emplace(side.element, side.attribute);
        if (!attribute_inserted && attribute_entry->second != side.attribute)
        {
          throw std::runtime_error(
              "Tetrahedral singular surface postprocessing found inconsistent adjacent "
              "element attributes!");
        }
        constexpr double tolerance = 256.0 * std::numeric_limits<double>::epsilon();
        const auto &previous = entry->second[face_node];
        if (seen->second[face_node] && (std::abs(previous.x - side.point.x) > tolerance ||
                                        std::abs(previous.y - side.point.y) > tolerance ||
                                        std::abs(previous.z - side.point.z) > tolerance))
        {
          throw std::runtime_error(
              "Tetrahedral singular surface postprocessing found inconsistent face "
              "orientation!");
        }
        entry->second[face_node] = side.point;
        seen->second[face_node] = true;
      }
    }
    if (!has_selected_side)
    {
      continue;
    }
    if (element_attributes.size() != element_face_points.size() ||
        element_face_points_seen.size() != element_face_points.size())
    {
      throw std::runtime_error(
          "Tetrahedral singular surface postprocessing found incomplete adjacent "
          "element attributes!");
    }
    struct ElementFaceData
    {
      int element;
      std::array<mfem::IntegrationPoint, 3> points;
      std::array<int, 3> nodes;
      bool complete_face;
      int attribute;
    };
    std::vector<ElementFaceData> element_faces;
    element_faces.reserve(element_face_points.size());
    for (const auto &[element, points] : element_face_points)
    {
      const auto seen = element_face_points_seen.find(element);
      if (seen == element_face_points_seen.end() ||
          !std::all_of(seen->second.begin(), seen->second.end(),
                       [](bool value) { return value; }))
      {
        throw std::runtime_error(
            "Tetrahedral singular surface postprocessing could not map all boundary "
            "vertices into an adjacent element!");
      }
      const auto containing_face = GetContainingTetrahedronFaceNodes(points);
      auto mapped_points = points;
      std::array<int, 3> nodes;
      bool complete_face = true;
      for (int face_node = 0; face_node < 3; face_node++)
      {
        nodes[face_node] = FindTetrahedronVertex(points[face_node]);
        complete_face = complete_face && nodes[face_node] >= 0;
      }
      if (complete_face)
      {
        auto sorted = nodes;
        std::sort(sorted.begin(), sorted.end());
        if (sorted != containing_face)
        {
          throw std::runtime_error(
              "Tetrahedral singular surface postprocessing mapped an invalid complete "
              "element face!");
        }
        // Preserve the exact complete-face map used by conforming integration. Near-vertex
        // transformation roundoff is significant when evaluating a singular basis.
        for (int face_node = 0; face_node < 3; face_node++)
        {
          mapped_points[face_node].Set3(nodes[face_node] == 1 ? 1.0 : 0.0,
                                        nodes[face_node] == 2 ? 1.0 : 0.0,
                                        nodes[face_node] == 3 ? 1.0 : 0.0);
        }
      }
      else
      {
        nodes = containing_face;
      }
      element_faces.push_back(
          {element, mapped_points, nodes, complete_face, element_attributes.at(element)});
    }

    std::vector<TetrahedronFaceSingularity> singularities;
    const auto add_singularity =
        [&singularities](fem::singular::TetrahedronFaceSingularityType type,
                         std::array<int, 2> nodes, double nu)
    {
      if (type == fem::singular::TetrahedronFaceSingularityType::EDGE &&
          nodes[1] < nodes[0])
      {
        std::swap(nodes[0], nodes[1]);
      }
      const auto existing = std::find_if(
          singularities.begin(), singularities.end(), [type, &nodes](const auto &feature)
          { return feature.type == type && feature.nodes == nodes; });
      if (existing == singularities.end())
      {
        singularities.push_back({type, nodes, nu});
      }
      else if (existing->nu != nu)
      {
        throw std::runtime_error(
            "Tetrahedral singular surface face has inconsistent feature exponents!");
      }
    };
    for (const auto &element_face : element_faces)
    {
      const int element = element_face.element;
      const auto &face_nodes = element_face.nodes;
      const auto local_singularities =
          full_wave
              ? field_evaluators.front().first->GetElementFaceSingularities(element,
                                                                            face_nodes)
              : real_gradient_evaluator->GetElementFaceSingularities(element, face_nodes);
      if (!element_face.complete_face)
      {
        if (!local_singularities.empty())
        {
          throw std::runtime_error(
              "Tetrahedral singular surface postprocessing encountered a nonconforming "
              "subface on an enriched singular trace!");
        }
        continue;
      }
      const auto to_face_node = [&face_nodes](int local_node)
      {
        const auto location = std::find(face_nodes.begin(), face_nodes.end(), local_node);
        if (location == face_nodes.end())
        {
          throw std::runtime_error(
              "Tetrahedral singular surface feature is not contained in its face!");
        }
        return static_cast<int>(std::distance(face_nodes.begin(), location));
      };
      for (const auto &feature : local_singularities)
      {
        std::array<int, 2> nodes{to_face_node(feature.nodes[0]), -1};
        if (feature.type == fem::singular::TetrahedronFaceSingularityType::EDGE)
        {
          nodes[1] = to_face_node(feature.nodes[1]);
        }
        add_singularity(feature.type, nodes, feature.nu);
      }
    }

    std::array<double, 3> cutoffs{};
    const bool curved_face =
        !fem::singular::IsAffineElementTransformation(*boundary_transformation);
    bool curved_edge_cutoff = false;
    bool curved_multi_edge_cutoff = false;
    if (geometry.edge_cutoff > 0.0)
    {
      const int edge_count = static_cast<int>(std::count_if(
          singularities.begin(), singularities.end(), [](const auto &feature)
          { return feature.type == fem::singular::TetrahedronFaceSingularityType::EDGE; }));
      curved_edge_cutoff = curved_face && edge_count > 0;
      curved_multi_edge_cutoff = curved_edge_cutoff && edge_count > 1;
      for (const auto &feature : singularities)
      {
        if (feature.type != fem::singular::TetrahedronFaceSingularityType::EDGE ||
            curved_edge_cutoff)
        {
          continue;
        }
        const int opposite_node = 3 - feature.nodes[0] - feature.nodes[1];
        const double coordinate =
            geometry.edge_cutoff /
            GetTriangleAltitude(*boundary_transformation, opposite_node);
        cutoffs[opposite_node] = std::max(cutoffs[opposite_node], coordinate);
      }
    }
    const double cutoff_sum = cutoffs[0] + cutoffs[1] + cutoffs[2];
    if (!std::isfinite(cutoff_sum) || !(cutoff_sum < 1.0))
    {
      throw std::domain_error(
          "Tetrahedral singular surface EdgeCutoff neighborhoods overlap on a face!");
    }
    for (const auto &feature : singularities)
    {
      if (feature.type == fem::singular::TetrahedronFaceSingularityType::EDGE)
      {
        const int opposite_node = 3 - feature.nodes[0] - feature.nodes[1];
        if (geometry.edge_cutoff == 0.0 && cutoffs[opposite_node] == 0.0 &&
            !(feature.nu > 0.5))
        {
          throw std::domain_error(
              "Tetrahedral singular surface participation is nonintegrable at an ideal "
              "thin-sheet edge with nu <= 1/2; specify EdgeCutoff or use a "
              "surface-response model!");
        }
      }
    }
    std::vector<TetrahedronFaceSingularity> chart_singularities;
    for (const auto &feature : singularities)
    {
      const bool node_on_edge =
          feature.type == fem::singular::TetrahedronFaceSingularityType::NODE &&
          std::any_of(singularities.begin(), singularities.end(),
                      [&feature](const auto &candidate)
                      {
                        return candidate.type ==
                                   fem::singular::TetrahedronFaceSingularityType::EDGE &&
                               (candidate.nodes[0] == feature.nodes[0] ||
                                candidate.nodes[1] == feature.nodes[0]);
                      });
      if (!node_on_edge)
      {
        chart_singularities.push_back(feature);
      }
    }
    if (curved_edge_cutoff && !curved_multi_edge_cutoff &&
        (chart_singularities.size() != 1 ||
         chart_singularities.front().type !=
             fem::singular::TetrahedronFaceSingularityType::EDGE))
    {
      throw std::runtime_error(
          "Physical EdgeCutoff integration on a curved singular face currently "
          "requires one isolated singular-edge chart!");
    }
    std::vector<PhysicalFaceEdge> physical_cutoff_edges;
    std::array<std::vector<double>, 3> cutoff_tangent_breakpoints;
    if (curved_multi_edge_cutoff)
    {
      for (const auto &feature : singularities)
      {
        if (feature.type == fem::singular::TetrahedronFaceSingularityType::EDGE)
        {
          physical_cutoff_edges.push_back(
              GetPhysicalFaceEdge(*boundary_transformation, feature.nodes));
        }
      }
      constexpr int tangent_search_intervals = 256;
      for (int boundary_edge = 0; boundary_edge < 3; boundary_edge++)
      {
        auto &breakpoints = cutoff_tangent_breakpoints[boundary_edge];
        breakpoints.push_back(0.0);
        double previous_tangent = 0.0;
        int previous_active =
            GetFaceCutoffRay(*boundary_transformation, physical_cutoff_edges, boundary_edge,
                             previous_tangent, geometry.edge_cutoff)
                .active_edge;
        const auto locate_transitions = [&](auto &&self, double lower, int lower_active,
                                            double upper, int upper_active,
                                            int depth) -> void
        {
          if (lower_active == upper_active)
          {
            return;
          }
          if (depth >= 64 ||
              upper - lower <= 64.0 * std::numeric_limits<double>::epsilon() *
                                   std::max({1.0, std::abs(lower), std::abs(upper)}))
          {
            breakpoints.push_back(0.5 * (lower + upper));
            return;
          }
          const double midpoint = 0.5 * (lower + upper);
          const int midpoint_active =
              GetFaceCutoffRay(*boundary_transformation, physical_cutoff_edges,
                               boundary_edge, midpoint, geometry.edge_cutoff)
                  .active_edge;
          if (lower_active != midpoint_active)
          {
            self(self, lower, lower_active, midpoint, midpoint_active, depth + 1);
          }
          if (midpoint_active != upper_active)
          {
            self(self, midpoint, midpoint_active, upper, upper_active, depth + 1);
          }
        };
        for (int interval = 1; interval <= tangent_search_intervals; interval++)
        {
          const double tangent = static_cast<double>(interval) / tangent_search_intervals;
          const int active =
              GetFaceCutoffRay(*boundary_transformation, physical_cutoff_edges,
                               boundary_edge, tangent, geometry.edge_cutoff)
                  .active_edge;
          if (active != previous_active)
          {
            locate_transitions(locate_transitions, previous_tangent, previous_active,
                               tangent, active, 0);
          }
          previous_tangent = tangent;
          previous_active = active;
        }
        breakpoints.push_back(1.0);
        std::sort(breakpoints.begin(), breakpoints.end());
        breakpoints.erase(
            std::unique(breakpoints.begin(), breakpoints.end(),
                        [](double left, double right)
                        {
                          return std::abs(left - right) <=
                                 64.0 * std::numeric_limits<double>::epsilon() *
                                     std::max({1.0, std::abs(left), std::abs(right)});
                        }),
            breakpoints.end());
      }
    }

    struct AccumulatedField
    {
      fem::singular::Vector3 real{};
      fem::singular::Vector3 imaginary{};
      int count = 0;
      int attribute = -1;
    };
    std::vector<fem::singular::NDFieldValuePair> evaluated_fields(mode_count);
    std::vector<std::array<AccumulatedField, 3>> fields(mode_count);
    fem::singular::Vector3 affine_normal{};
    double affine_boundary_weight = 0.0;
    if (!curved_face)
    {
      mfem::IntegrationPoint point;
      point.Set2(1.0 / 3.0, 1.0 / 3.0);
      boundary_transformation->SetIntPoint(&point);
      double normal_data[3];
      mfem::Vector normal_vector(normal_data, 3);
      BdrGridFunctionCoefficient::GetNormal(*boundary_transformation, normal_vector);
      affine_normal = {normal_data[0], normal_data[1], normal_data[2]};
      affine_boundary_weight = boundary_transformation->Weight();
    }
    const auto accumulate_density = [&](const std::array<double, 3> &lambda,
                                        double quadrature_weight,
                                        std::vector<long double> &values)
    {
      if (values.size() != component_count)
      {
        throw std::logic_error(
            "Tetrahedral singular surface integration accumulator has invalid size!");
      }
      mfem::IntegrationPoint point;
      point.Set2(lambda[1], lambda[2]);
      std::fill(fields.begin(), fields.end(), std::array<AccumulatedField, 3>{});
      for (const auto &element_face : element_faces)
      {
        const int element = element_face.element;
        const auto &element_points = element_face.points;
        mfem::IntegrationPoint element_point;
        element_point.Set3(0.0, 0.0, 0.0);
        for (int face_node = 0; face_node < 3; face_node++)
        {
          element_point.x += lambda[face_node] * element_points[face_node].x;
          element_point.y += lambda[face_node] * element_points[face_node].y;
          element_point.z += lambda[face_node] * element_points[face_node].z;
        }
        if (full_wave)
        {
          field_evaluators.front().first->EvaluateValueClosureBatch(
              field_evaluators, element, element_point, evaluated_fields);
        }
        else
        {
          const auto real_value =
              real_gradient_evaluator->EvaluateClosure(element, element_point);
          const auto imaginary_value =
              imaginary_gradient_evaluator->EvaluateClosure(element, element_point);
          for (int d = 0; d < 3; d++)
          {
            evaluated_fields[0].first[d] = -real_value.gradient[d];
            evaluated_fields[0].second[d] = -imaginary_value.gradient[d];
          }
        }

        const int element_attribute = element_face.attribute;
        const int material_group = IsVacuumSide(material, element_attribute) ? 1 : 2;
        for (std::size_t mode = 0; mode < mode_count; mode++)
        {
          for (int group : {0, material_group})
          {
            auto &field = fields[mode][group];
            field.count++;
            if (field.attribute < 0)
            {
              field.attribute = element_attribute;
            }
            for (int d = 0; d < 3; d++)
            {
              field.real[d] += evaluated_fields[mode].first[d];
              field.imaginary[d] += evaluated_fields[mode].second[d];
            }
          }
        }
      }

      fem::singular::Vector3 normal = affine_normal;
      double boundary_weight = affine_boundary_weight;
      if (curved_face)
      {
        boundary_transformation->SetIntPoint(&point);
        double normal_data[3];
        mfem::Vector normal_vector(normal_data, 3);
        BdrGridFunctionCoefficient::GetNormal(*boundary_transformation, normal_vector);
        normal = {normal_data[0], normal_data[1], normal_data[2]};
        boundary_weight = boundary_transformation->Weight();
      }
      for (std::size_t mode = 0; mode < mode_count; mode++)
      {
        for (std::size_t interface_index = 0; interface_index < interfaces.size();
             interface_index++)
        {
          const auto &interface = *interfaces[interface_index];
          const int material_group =
              interface.type == InterfaceDielectric::DEFAULT
                  ? 0
                  : (interface.type == InterfaceDielectric::MS ? 2 : 1);
          const auto &field = fields[mode][material_group];
          if (field.count == 0)
          {
            continue;
          }
          const double average_scale = 1.0 / static_cast<double>(field.count);
          fem::singular::Vector3 electric_real = field.real;
          fem::singular::Vector3 electric_imaginary = field.imaginary;
          for (int d = 0; d < 3; d++)
          {
            electric_real[d] *= average_scale;
            electric_imaginary[d] *= average_scale;
          }
          const double real_normal = Dot(electric_real, normal);
          const double imaginary_normal = Dot(electric_imaginary, normal);
          const double normal_product =
              real_normal * real_normal + imaginary_normal * imaginary_normal;
          const double field_product = Dot(electric_real, electric_real) +
                                       Dot(electric_imaginary, electric_imaginary);

          double surface_density = 0.0;
          if (interface.type == InterfaceDielectric::DEFAULT)
          {
            surface_density =
                0.5 * interface.thickness * interface.permittivity * field_product;
          }
          else if (interface.type == InterfaceDielectric::MA)
          {
            surface_density =
                0.5 * (interface.thickness / interface.permittivity) * normal_product;
          }
          else if (interface.type == InterfaceDielectric::MS)
          {
            fem::singular::Vector3 displacement_real{};
            fem::singular::Vector3 displacement_imaginary{};
            mfem::Vector electric_real_vector(electric_real.data(), 3);
            mfem::Vector electric_imaginary_vector(electric_imaginary.data(), 3);
            mfem::Vector displacement_real_vector(displacement_real.data(), 3);
            mfem::Vector displacement_imaginary_vector(displacement_imaginary.data(), 3);
            const auto &permittivity = material.GetPermittivityReal(field.attribute);
            permittivity.Mult(electric_real_vector, displacement_real_vector);
            permittivity.Mult(electric_imaginary_vector, displacement_imaginary_vector);
            const double displacement_normal_real = Dot(displacement_real, normal);
            const double displacement_normal_imaginary =
                Dot(displacement_imaginary, normal);
            surface_density =
                0.5 * (interface.thickness / interface.permittivity) *
                (displacement_normal_real * displacement_normal_real +
                 displacement_normal_imaginary * displacement_normal_imaginary);
          }
          else
          {
            const double tangential_product = field_product - normal_product;
            surface_density = 0.5 * interface.thickness *
                              (interface.permittivity * tangential_product +
                               normal_product / interface.permittivity);
          }
          const double contribution = quadrature_weight * boundary_weight * surface_density;
          if (!std::isfinite(contribution))
          {
            throw std::runtime_error(
                "Tetrahedral singular surface quadrature produced a nonfinite "
                "contribution!");
          }
          values[mode * interfaces.size() + interface_index] += contribution;
        }
      }
    };

    const double simplex_scale = 1.0 - cutoff_sum;
    const auto feature_radius =
        [](const TetrahedronFaceSingularity &feature, const std::array<double, 3> &lambda)
    {
      if (feature.type == fem::singular::TetrahedronFaceSingularityType::NODE)
      {
        return 1.0 - lambda[feature.nodes[0]];
      }
      const int opposite_node = 3 - feature.nodes[0] - feature.nodes[1];
      return lambda[opposite_node];
    };
    const auto partition_weight =
        [&chart_singularities, &feature_radius](std::size_t selected,
                                                const std::array<double, 3> &lambda)
    {
      if (chart_singularities.size() == 1)
      {
        return 1.0;
      }
      std::vector<double> radii;
      radii.reserve(chart_singularities.size());
      double minimum = std::numeric_limits<double>::infinity();
      for (const auto &feature : chart_singularities)
      {
        const double radius = feature_radius(feature, lambda);
        if (!std::isfinite(radius) || !(radius > 0.0))
        {
          throw std::runtime_error(
              "Tetrahedral singular surface partition encountered a nonpositive "
              "feature radius!");
        }
        radii.push_back(radius);
        minimum = std::min(minimum, radius);
      }
      double denominator = 0.0;
      double numerator = 0.0;
      for (std::size_t feature = 0; feature < radii.size(); feature++)
      {
        const double score = std::pow(minimum / radii[feature],
                                      fem::singular::MultiFeatureDuffyPartitionPower);
        denominator += score;
        numerator = feature == selected ? score : numerator;
      }
      if (!std::isfinite(denominator) || !(denominator > 0.0) || !std::isfinite(numerator))
      {
        throw std::runtime_error(
            "Tetrahedral singular surface partition produced invalid weights!");
      }
      return numerator / denominator;
    };

    const auto integrate = [&](int order)
    {
      std::vector<long double> value(component_count, 0.0L);
      const auto standard_segment_rule =
          fem::singular::BuildWeightedSegmentQuadrature(order, 0.0, 0.0);
      if (curved_multi_edge_cutoff)
      {
        // The retained physical domain is the face minus the union of all edge
        // cutoff tubes. Partition the reference triangle into three sectors from
        // its barycenter and terminate each ray at its first tube intersection.
        // The sectors overlap only on measure-zero rays.
        constexpr double radial_power = 2.0;
        const auto integrate_interval = [&](int boundary_edge, int first, int second,
                                            double interval_lower, double interval_upper)
        {
          std::vector<long double> result(component_count, 0.0L);
          for (const auto &tangent_quadrature : standard_segment_rule)
          {
            const double tangent = interval_lower + (interval_upper - interval_lower) *
                                                        tangent_quadrature.coordinate;
            const double maximum_radius =
                GetFaceCutoffRay(*boundary_transformation, physical_cutoff_edges,
                                 boundary_edge, tangent, geometry.edge_cutoff)
                    .radius;
            for (const auto &radial_quadrature : standard_segment_rule)
            {
              const double outer_coordinate =
                  1.0 - std::pow(1.0 - radial_quadrature.coordinate, radial_power);
              const double radius = maximum_radius * outer_coordinate;
              const double dr_du =
                  maximum_radius * radial_power *
                  std::pow(1.0 - radial_quadrature.coordinate, radial_power - 1.0);
              std::array<double, 3> lambda{(1.0 - radius) / 3.0, (1.0 - radius) / 3.0,
                                           (1.0 - radius) / 3.0};
              lambda[first] += radius * (1.0 - tangent);
              lambda[second] += radius * tangent;
              const double quadrature_weight =
                  radial_quadrature.weight * tangent_quadrature.weight *
                  (interval_upper - interval_lower) * dr_du * radius / 3.0;
              accumulate_density(lambda, quadrature_weight, result);
            }
          }
          return result;
        };
        const auto integrate_tangent_interval =
            [&](auto &&self, int boundary_edge, int first, int second, double lower,
                double upper, const std::vector<long double> &coarse,
                int depth) -> std::vector<long double>
        {
          const double midpoint = 0.5 * (lower + upper);
          const auto left =
              integrate_interval(boundary_edge, first, second, lower, midpoint);
          const auto right =
              integrate_interval(boundary_edge, first, second, midpoint, upper);
          std::vector<long double> fine(component_count);
          for (std::size_t component = 0; component < fine.size(); component++)
          {
            fine[component] = left[component] + right[component];
          }
          const double interval_fraction = upper - lower;
          bool converged = true;
          std::size_t failed_component = 0;
          double failed_error = 0.0;
          double failed_tolerance = 0.0;
          for (std::size_t component = 0; component < fine.size(); component++)
          {
            const double fine_value = static_cast<double>(fine[component]);
            const double coarse_value = static_cast<double>(coarse[component]);
            const double scale = std::max(std::abs(fine_value), std::abs(coarse_value));
            const double error = 8.0 * std::abs(fine_value - coarse_value) +
                                 256.0 * std::numeric_limits<double>::epsilon() * scale;
            const double tolerance =
                0.125 * (options.absolute_tolerance * interval_fraction / 3.0 +
                         options.relative_tolerance * std::abs(fine_value));
            if (!std::isfinite(error) || error > tolerance)
            {
              converged = false;
              failed_component = component;
              failed_error = error;
              failed_tolerance = tolerance;
              break;
            }
          }
          if (converged)
          {
            return fine;
          }
          if (depth >= options.maximum_subdivisions)
          {
            throw std::runtime_error(fmt::format(
                "Tetrahedral curved-face multi-edge cutoff tangent quadrature did "
                "not meet tolerance: value = {:.17g}, estimated absolute error = "
                "{:.17g}, tolerance = {:.17g}, component = {}, interval = "
                "[{:.17g}, {:.17g}], order = {}, subdivisions = {}!",
                static_cast<double>(fine[failed_component]), failed_error, failed_tolerance,
                failed_component, lower, upper, order, options.maximum_subdivisions));
          }
          const auto refined_left =
              self(self, boundary_edge, first, second, lower, midpoint, left, depth + 1);
          const auto refined_right =
              self(self, boundary_edge, first, second, midpoint, upper, right, depth + 1);
          std::vector<long double> refined(component_count);
          for (std::size_t component = 0; component < refined.size(); component++)
          {
            refined[component] = refined_left[component] + refined_right[component];
          }
          return refined;
        };
        for (int boundary_edge = 0; boundary_edge < 3; boundary_edge++)
        {
          const int first = (boundary_edge + 1) % 3;
          const int second = (boundary_edge + 2) % 3;
          const auto &breakpoints = cutoff_tangent_breakpoints[boundary_edge];
          if (breakpoints.size() < 2)
          {
            throw std::runtime_error(
                "Tetrahedral curved-face multi-edge cutoff chart has no tangent "
                "partition!");
          }
          for (std::size_t interval = 1; interval < breakpoints.size(); interval++)
          {
            const double lower = breakpoints[interval - 1];
            const double upper = breakpoints[interval];
            const auto coarse =
                integrate_interval(boundary_edge, first, second, lower, upper);
            const auto refined =
                integrate_tangent_interval(integrate_tangent_interval, boundary_edge, first,
                                           second, lower, upper, coarse, 0);
            for (std::size_t component = 0; component < value.size(); component++)
            {
              value[component] += refined[component];
            }
          }
        }
        return value;
      }
      if (chart_singularities.empty())
      {
        const auto radial_rule =
            fem::singular::BuildWeightedSegmentQuadrature(order, 0.0, 0.0);
        for (const auto &radial : radial_rule)
        {
          for (const auto &tangent : standard_segment_rule)
          {
            const double r = radial.coordinate;
            const double t = tangent.coordinate;
            const std::array<double, 3> lambda{1.0 - r, r * (1.0 - t), r * t};
            accumulate_density(lambda, radial.weight * tangent.weight * r, value);
          }
        }
        return value;
      }

      for (std::size_t selected = 0; selected < chart_singularities.size(); selected++)
      {
        const auto &feature = chart_singularities[selected];
        const bool node =
            feature.type == fem::singular::TetrahedronFaceSingularityType::NODE;
        if (curved_edge_cutoff)
        {
          for (int endpoint = 0; endpoint < 2; endpoint++)
          {
            double endpoint_nu = 1.0;
            for (const auto &candidate : singularities)
            {
              if (candidate.type == fem::singular::TetrahedronFaceSingularityType::NODE &&
                  candidate.nodes[0] == feature.nodes[endpoint])
              {
                endpoint_nu = std::min(endpoint_nu, candidate.nu);
              }
            }
            const double tangent_power =
                std::abs(endpoint_nu - 2.0 / 3.0) <=
                        64.0 * std::numeric_limits<double>::epsilon()
                    ? 3.0
                    : (endpoint_nu < 1.0 ? fem::singular::H1DuffyRadialPower : 1.0);
            for (const auto &tangent : standard_segment_rule)
            {
              const double distance = 0.5 * std::pow(tangent.coordinate, tangent_power);
              const double dt_du =
                  0.5 * tangent_power * std::pow(tangent.coordinate, tangent_power - 1.0);
              const double t = endpoint == 0 ? distance : 1.0 - distance;
              const double floor = GetFaceEdgeCutoffCoordinate(
                  *boundary_transformation, feature.nodes, t, geometry.edge_cutoff);
              const double logarithmic_length = std::log(1.0 / floor);
              const auto radial_rule =
                  fem::singular::BuildWeightedSegmentQuadrature(order, 0.0, 0.0);
              for (const auto &radial : radial_rule)
              {
                const double r = floor * std::exp(logarithmic_length * radial.coordinate);
                const double dr_du = logarithmic_length * r;
                std::array<double, 3> lambda{};
                const int opposite_node = 3 - feature.nodes[0] - feature.nodes[1];
                lambda[opposite_node] = r;
                lambda[feature.nodes[0]] = (1.0 - r) * (1.0 - t);
                lambda[feature.nodes[1]] = (1.0 - r) * t;
                const double quadrature_weight =
                    radial.weight * tangent.weight * dt_du * (1.0 - r) * dr_du;
                accumulate_density(lambda, quadrature_weight, value);
              }
            }
          }
          continue;
        }
        int primary_node = feature.nodes[0];
        if (!node)
        {
          primary_node = 3 - feature.nodes[0] - feature.nodes[1];
        }
        const double floor =
            node ? cutoff_sum - cutoffs[primary_node] : cutoffs[primary_node];
        const double upper = floor + simplex_scale;
        if (!std::isfinite(floor) || floor < 0.0 || !std::isfinite(upper) ||
            !(upper > floor))
        {
          throw std::runtime_error(
              "Tetrahedral singular surface chart has an invalid radial interval!");
        }

        const auto accumulate = [&](double r, double radial_weight, double radial_jacobian,
                                    double radial_weight_function)
        {
          const auto accumulate_tangent = [&](double t, double tangent_weight)
          {
            std::array<double, 3> mu{};
            if (node)
            {
              mu[primary_node] = 1.0 - r;
              const int other_0 = (primary_node + 1) % 3;
              const int other_1 = (primary_node + 2) % 3;
              mu[other_0] = r * (1.0 - t);
              mu[other_1] = r * t;
            }
            else
            {
              mu[primary_node] = r;
              mu[feature.nodes[0]] = (1.0 - r) * (1.0 - t);
              mu[feature.nodes[1]] = (1.0 - r) * t;
            }
            std::array<double, 3> lambda{};
            for (int i = 0; i < 3; i++)
            {
              lambda[i] = cutoffs[i] + simplex_scale * mu[i];
            }
            const double chart_jacobian =
                simplex_scale * simplex_scale * (node ? r : 1.0 - r) * radial_jacobian;
            const double quadrature_weight =
                radial_weight * tangent_weight * chart_jacobian *
                partition_weight(selected, lambda) / radial_weight_function;
            accumulate_density(lambda, quadrature_weight, value);
          };
          if (node)
          {
            for (const auto &tangent : standard_segment_rule)
            {
              accumulate_tangent(tangent.coordinate, tangent.weight);
            }
            return;
          }

          for (int endpoint = 0; endpoint < 2; endpoint++)
          {
            double endpoint_nu = 1.0;
            for (const auto &candidate : singularities)
            {
              if (candidate.type == fem::singular::TetrahedronFaceSingularityType::NODE &&
                  candidate.nodes[0] == feature.nodes[endpoint])
              {
                endpoint_nu = std::min(endpoint_nu, candidate.nu);
              }
            }
            const double tangent_power =
                std::abs(endpoint_nu - 2.0 / 3.0) <=
                        64.0 * std::numeric_limits<double>::epsilon()
                    ? 3.0
                    : (endpoint_nu < 1.0 ? fem::singular::H1DuffyRadialPower : 1.0);
            for (const auto &tangent : standard_segment_rule)
            {
              const double distance = 0.5 * std::pow(tangent.coordinate, tangent_power);
              const double dt_du =
                  0.5 * tangent_power * std::pow(tangent.coordinate, tangent_power - 1.0);
              accumulate_tangent(endpoint == 0 ? distance : 1.0 - distance,
                                 tangent.weight * dt_du);
            }
          }
        };

        if (floor > 0.0)
        {
          const double logarithmic_length = std::log(upper / floor);
          const auto radial_rule =
              fem::singular::BuildWeightedSegmentQuadrature(order, 0.0, 0.0);
          for (const auto &radial : radial_rule)
          {
            const double radius = floor * std::exp(logarithmic_length * radial.coordinate);
            const double r = (radius - floor) / simplex_scale;
            const double dr_du = logarithmic_length * radius / simplex_scale;
            accumulate(r, radial.weight, dr_du, 1.0);
          }
        }
        else
        {
          // A cubic radial map polynomializes every one-third power in the
          // thick 270-degree edge expansion without pushing high-order
          // quadrature points unnecessarily close to the singular edge. Keep
          // the assembly convention's sixth-power map for other exponents.
          const double radial_power =
              !node && std::abs(feature.nu - 2.0 / 3.0) <=
                           64.0 * std::numeric_limits<double>::epsilon()
                  ? 3.0
                  : fem::singular::H1DuffyRadialPower;
          const double alpha = node ? 2.0 * radial_power * feature.nu - 1.0
                                    : radial_power * (2.0 * feature.nu - 1.0) - 1.0;
          const auto radial_rule =
              fem::singular::BuildWeightedSegmentQuadrature(order, alpha, 0.0);
          for (const auto &radial : radial_rule)
          {
            const double r = std::pow(radial.coordinate, radial_power);
            const double dr_du =
                radial_power * std::pow(radial.coordinate, radial_power - 1.0);
            const double radial_weight_function = std::pow(radial.coordinate, alpha);
            accumulate(r, radial.weight, dr_du, radial_weight_function);
          }
        }
      }
      return value;
    };

    if (options.fixed_subdivision)
    {
      const auto value = integrate(std::max(4, 2 * options.quadrature_order));
      for (std::size_t component = 0; component < local_energy.size(); component++)
      {
        local_energy[component] += value[component];
      }
      continue;
    }

    constexpr int order_increment = 4;
    int comparison_order = std::max(4, 2 * options.quadrature_order);
    auto comparison = integrate(comparison_order);
    auto value = comparison;
    std::size_t failed_component = 0;
    double failed_error = std::numeric_limits<double>::infinity();
    double failed_tolerance = 0.0;
    int order = comparison_order + order_increment;
    bool converged = false;
    for (int refinement = 0; refinement <= options.maximum_subdivisions; refinement++)
    {
      value = integrate(order);
      converged = true;
      for (std::size_t component = 0; component < value.size(); component++)
      {
        const double component_value = static_cast<double>(value[component]);
        const double component_comparison = static_cast<double>(comparison[component]);
        const double scale =
            std::max({1.0, std::abs(component_value), std::abs(component_comparison)});
        const double error = 8.0 * std::abs(component_value - component_comparison) +
                             256.0 * std::numeric_limits<double>::epsilon() * scale;
        const double tolerance = options.absolute_tolerance +
                                 options.relative_tolerance * std::abs(component_value);
        if (!std::isfinite(component_value) || !std::isfinite(component_comparison) ||
            !std::isfinite(error) || error > tolerance)
        {
          converged = false;
          failed_component = component;
          failed_error = error;
          failed_tolerance = tolerance;
          break;
        }
      }
      if (converged)
      {
        break;
      }
      if (refinement < options.maximum_subdivisions)
      {
        comparison = value;
        comparison_order = order;
        order += order_increment;
      }
    }
    if (!converged)
    {
      throw std::runtime_error(fmt::format(
          "Tetrahedral singular surface quadrature did not meet tolerance: value = "
          "{:.17g}, comparison = {:.17g}, estimated absolute error = {:.17g}, "
          "tolerance = {:.17g}, component = {}, orders = {}/{}, order refinements = "
          "{}!",
          static_cast<double>(value[failed_component]),
          static_cast<double>(comparison[failed_component]), failed_error, failed_tolerance,
          failed_component, order, comparison_order, options.maximum_subdivisions));
    }
    for (std::size_t component = 0; component < local_energy.size(); component++)
    {
      local_energy[component] += value[component];
    }
  }
  std::vector<double> energy(local_energy.size());
  std::transform(local_energy.begin(), local_energy.end(), energy.begin(),
                 [](long double value) { return static_cast<double>(value); });
  Mpi::GlobalSum(static_cast<int>(energy.size()), energy.data(), mesh.GetComm());
  return energy;
}

std::vector<TetrahedronSingularSurfacePostOperator::Measurement>
TetrahedronSingularSurfacePostOperator::Measure(
    fem::singular::EnrichedNDFieldEvaluator &real_evaluator,
    fem::singular::EnrichedNDFieldEvaluator &imaginary_evaluator,
    double total_electric_energy,
    const fem::singular::AdaptiveAssemblyOptions &options) const
{
  auto measurements =
      Measure({{&real_evaluator, &imaginary_evaluator}}, {total_electric_energy}, options);
  if (measurements.size() != 1)
  {
    throw std::logic_error(
        "Tetrahedral singular surface single-field measurement returned invalid "
        "dimensions!");
  }
  return std::move(measurements.front());
}

std::vector<std::vector<TetrahedronSingularSurfacePostOperator::Measurement>>
TetrahedronSingularSurfacePostOperator::Measure(
    const std::vector<NDFieldEvaluatorPair> &field_evaluators,
    const std::vector<double> &total_electric_energies,
    const fem::singular::AdaptiveAssemblyOptions &options) const
{
  if (field_evaluators.empty() ||
      total_electric_energies.size() != field_evaluators.size() ||
      std::any_of(total_electric_energies.begin(), total_electric_energies.end(),
                  [](double energy) { return !std::isfinite(energy) || !(energy > 0.0); }))
  {
    throw std::invalid_argument(
        "Tetrahedral singular surface participation requires matching fields and "
        "positive total electric energies!");
  }
  std::vector<std::map<int, double>> energies(field_evaluators.size());
  for (const auto &group : GetInterfaceGroups())
  {
    std::vector<const InterfaceData *> grouped_interfaces;
    grouped_interfaces.reserve(group.size());
    for (const auto &[index, interface] : group)
    {
      grouped_interfaces.push_back(interface);
    }
    const auto start = Timer::Now();
    const auto grouped_energies = IntegrateInterfaces(grouped_interfaces, field_evaluators,
                                                      nullptr, nullptr, options);
    double elapsed = Timer::Duration(Timer::Now() - start).count();
    Mpi::GlobalMax(1, &elapsed, fespace.GetComm());
    std::string indices = fmt::format("{}", group.front().first);
    for (std::size_t interface = 1; interface < group.size(); interface++)
    {
      indices += fmt::format(",{}", group[interface].first);
    }
    Mpi::Print(" Singular surface interface group {} for {} field{} integration (s): "
               "{:.3f}\n",
               indices, field_evaluators.size(), field_evaluators.size() == 1 ? "" : "s",
               elapsed);
    if (grouped_energies.size() != group.size() * field_evaluators.size())
    {
      throw std::logic_error(
          "Tetrahedral singular surface integration returned invalid dimensions!");
    }
    for (std::size_t field = 0; field < field_evaluators.size(); field++)
    {
      for (std::size_t interface = 0; interface < group.size(); interface++)
      {
        energies[field].emplace(group[interface].first,
                                grouped_energies[field * group.size() + interface]);
      }
    }
  }

  std::vector<std::vector<Measurement>> measurements(field_evaluators.size());
  for (std::size_t field = 0; field < field_evaluators.size(); field++)
  {
    measurements[field].reserve(interfaces.size());
    for (const auto &[index, interface] : interfaces)
    {
      const double energy = energies[field].at(index);
      const double participation = energy / total_electric_energies[field];
      const double quality_factor = participation == 0.0 || interface.loss_tangent == 0.0
                                        ? mfem::infinity()
                                        : 1.0 / (participation * interface.loss_tangent);
      measurements[field].push_back(
          {index, energy, interface.loss_tangent, participation, quality_factor});
    }
  }
  return measurements;
}

std::vector<TetrahedronSingularSurfacePostOperator::Measurement>
TetrahedronSingularSurfacePostOperator::MeasureElectrostatic(
    fem::singular::EnrichedH1FieldEvaluator &real_evaluator,
    fem::singular::EnrichedH1FieldEvaluator &imaginary_evaluator,
    double total_electric_energy,
    const fem::singular::AdaptiveAssemblyOptions &options) const
{
  if (!std::isfinite(total_electric_energy) || !(total_electric_energy > 0.0))
  {
    throw std::invalid_argument(
        "Tetrahedral singular surface participation requires positive total electric "
        "energy!");
  }
  std::map<int, double> energies;
  for (const auto &group : GetInterfaceGroups())
  {
    std::vector<const InterfaceData *> grouped_interfaces;
    grouped_interfaces.reserve(group.size());
    for (const auto &[index, interface] : group)
    {
      grouped_interfaces.push_back(interface);
    }
    const auto grouped_energies = IntegrateInterfaces(
        grouped_interfaces, {}, &real_evaluator, &imaginary_evaluator, options);
    if (grouped_energies.size() != group.size())
    {
      throw std::logic_error(
          "Tetrahedral singular surface integration returned invalid dimensions!");
    }
    for (std::size_t interface = 0; interface < group.size(); interface++)
    {
      energies.emplace(group[interface].first, grouped_energies[interface]);
    }
  }

  std::vector<Measurement> measurements;
  measurements.reserve(interfaces.size());
  for (const auto &[index, interface] : interfaces)
  {
    const double energy = energies.at(index);
    const double participation = energy / total_electric_energy;
    const double quality_factor = participation == 0.0 || interface.loss_tangent == 0.0
                                      ? mfem::infinity()
                                      : 1.0 / (participation * interface.loss_tangent);
    measurements.push_back(
        {index, energy, interface.loss_tangent, participation, quality_factor});
  }
  return measurements;
}

}  // namespace palace
