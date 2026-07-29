// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularsurfacepostoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <fmt/format.h>

#include "fem/coefficient.hpp"
#include "fem/singularelements.hpp"
#include "fem/singularfield.hpp"
#include "fem/singulargeometry.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

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

int GetTetrahedronVertex(const mfem::IntegrationPoint &point)
{
  const std::array<double, 4> lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                     point.z};
  const auto maximum = std::max_element(lambda.begin(), lambda.end());
  const int vertex = static_cast<int>(std::distance(lambda.begin(), maximum));
  if (*maximum < 1.0 - 256.0 * std::numeric_limits<double>::epsilon())
  {
    throw std::runtime_error(
        "Singular surface quadrature could not identify a tetrahedron face vertex!");
  }
  return vertex;
}

double GetBoundaryCutoffCoordinate(mfem::ElementTransformation &transformation,
                                   int endpoint, double cutoff)
{
  if (!(cutoff > 0.0) || (endpoint != 0 && endpoint != 1))
  {
    throw std::invalid_argument("Singular surface cutoff mapping received invalid input!");
  }
  mfem::IntegrationPoint endpoint_point;
  endpoint_point.x = static_cast<double>(endpoint);
  mfem::Vector endpoint_position(transformation.GetSpaceDim());
  transformation.Transform(endpoint_point, endpoint_position);
  const auto distance = [&](double coordinate)
  {
    mfem::IntegrationPoint point;
    point.x = coordinate;
    mfem::Vector position(transformation.GetSpaceDim());
    transformation.Transform(point, position);
    position -= endpoint_position;
    return position.Norml2();
  };
  const double segment_length = distance(static_cast<double>(1 - endpoint));
  if (!std::isfinite(segment_length) || !(segment_length > cutoff))
  {
    throw std::domain_error(
        "Singular surface EdgeCutoff must be smaller than every incident boundary "
        "segment!");
  }

  double lower = 0.0;
  double upper = 1.0;
  for (int iteration = 0; iteration < 64; iteration++)
  {
    const double midpoint = 0.5 * (lower + upper);
    const double coordinate = endpoint == 0 ? midpoint : 1.0 - midpoint;
    if (distance(coordinate) < cutoff)
    {
      lower = midpoint;
    }
    else
    {
      upper = midpoint;
    }
  }
  const double radial_coordinate = 0.5 * (lower + upper);
  return endpoint == 0 ? radial_coordinate : 1.0 - radial_coordinate;
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
          endpoint_exponents[endpoint] = std::min(
              endpoint_exponents[endpoint],
              real_evaluator
                  ? real_evaluator->GetElementNodeSingularExponent(side.element, node)
                  : real_gradient_evaluator->GetElementNodeSingularExponent(side.element,
                                                                            node));
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
    double lower = 0.0;
    double upper = 1.0;
    if (interface.edge_cutoff > 0.0 && endpoint_exponents[0] < 1.0)
    {
      lower =
          GetBoundaryCutoffCoordinate(*boundary_transformation, 0, interface.edge_cutoff);
    }
    if (interface.edge_cutoff > 0.0 && endpoint_exponents[1] < 1.0)
    {
      upper =
          GetBoundaryCutoffCoordinate(*boundary_transformation, 1, interface.edge_cutoff);
    }
    if (!(lower < upper))
    {
      throw std::domain_error(
          "Singular surface EdgeCutoff neighborhoods overlap on a boundary segment!");
    }

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
    local_energy += IntegratePowerExpansion(coefficient_expansion, lower, upper, options);
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

double TetrahedronSingularSurfacePostOperator::IntegrateInterface(
    const InterfaceData &interface, fem::singular::EnrichedNDFieldEvaluator *real_evaluator,
    fem::singular::EnrichedNDFieldEvaluator *imaginary_evaluator,
    fem::singular::EnrichedH1FieldEvaluator *real_gradient_evaluator,
    fem::singular::EnrichedH1FieldEvaluator *imaginary_gradient_evaluator,
    const fem::singular::AdaptiveAssemblyOptions &options) const
{
  if ((real_evaluator == nullptr) != (imaginary_evaluator == nullptr) ||
      (real_gradient_evaluator == nullptr) != (imaginary_gradient_evaluator == nullptr) ||
      (real_evaluator == nullptr) == (real_gradient_evaluator == nullptr))
  {
    throw std::invalid_argument(
        "Tetrahedral singular surface postprocessing requires exactly one complete "
        "full-wave or electrostatic evaluator pair!");
  }
  if (options.quadrature_order < 1 || !std::isfinite(options.absolute_tolerance) ||
      options.absolute_tolerance < 0.0 || !std::isfinite(options.relative_tolerance) ||
      options.relative_tolerance < 0.0 ||
      !(options.absolute_tolerance > 0.0 || options.relative_tolerance > 0.0) ||
      options.maximum_subdivisions < 1)
  {
    throw std::invalid_argument(
        "Tetrahedral singular surface postprocessing received invalid adaptive "
        "quadrature options!");
  }

  auto &mesh = *fespace.GetParMesh();
  mfem::FaceElementTransformations face;
  mfem::IsoparametricTransformation element1, element2;
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
        boundary_transformation->GetGeometryType() != mfem::Geometry::TRIANGLE)
    {
      throw std::runtime_error(
          "Tetrahedral singular surface postprocessing requires triangle boundaries!");
    }

    std::map<int, std::array<int, 3>> element_face_nodes;
    bool has_selected_side = false;
    for (int face_node = 0; face_node < 3; face_node++)
    {
      mfem::IntegrationPoint point;
      point.Set2(face_node == 1 ? 1.0 : 0.0, face_node == 2 ? 1.0 : 0.0);
      BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
          boundary, mesh, face, element1, element2, &point);
      const auto sides = GetSelectedTetrahedronSides(material, face, interface.type);
      has_selected_side = has_selected_side || !sides.empty();
      for (const auto &side : sides)
      {
        auto [entry, inserted] =
            element_face_nodes.emplace(side.element, std::array<int, 3>{-1, -1, -1});
        const int node = GetTetrahedronVertex(side.point);
        if (entry->second[face_node] >= 0 && entry->second[face_node] != node)
        {
          throw std::runtime_error(
              "Tetrahedral singular surface postprocessing found inconsistent face "
              "orientation!");
        }
        entry->second[face_node] = node;
      }
    }
    if (!has_selected_side)
    {
      continue;
    }
    for (const auto &[element, nodes] : element_face_nodes)
    {
      auto sorted = nodes;
      std::sort(sorted.begin(), sorted.end());
      if (sorted[0] < 0 || sorted[2] >= 4 ||
          std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
      {
        throw std::runtime_error(
            "Tetrahedral singular surface postprocessing could not map a complete "
            "adjacent element face!");
      }
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
    for (const auto &element_face : element_face_nodes)
    {
      const int element = element_face.first;
      const auto &face_nodes = element_face.second;
      const auto local_singularities =
          real_evaluator
              ? real_evaluator->GetElementFaceSingularities(element, face_nodes)
              : real_gradient_evaluator->GetElementFaceSingularities(element, face_nodes);
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
    if (interface.edge_cutoff > 0.0)
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
            interface.edge_cutoff /
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
        if (interface.edge_cutoff == 0.0 && cutoffs[opposite_node] == 0.0 &&
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
                             previous_tangent, interface.edge_cutoff)
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
                               boundary_edge, midpoint, interface.edge_cutoff)
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
                               boundary_edge, tangent, interface.edge_cutoff)
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

    const auto density = [&](const std::array<double, 3> &lambda)
    {
      mfem::IntegrationPoint point;
      point.Set2(lambda[1], lambda[2]);
      const bool inverted =
          BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
              boundary, mesh, face, element1, element2, &point);
      const auto sides = GetSelectedTetrahedronSides(material, face, interface.type);
      if (sides.empty())
      {
        return 0.0;
      }
      fem::singular::Vector3 electric_real{};
      fem::singular::Vector3 electric_imaginary{};
      const double average_scale = 1.0 / static_cast<double>(sides.size());
      for (const auto &side : sides)
      {
        if (real_evaluator)
        {
          const auto real_value = real_evaluator->EvaluateClosure(side.element, side.point);
          const auto imaginary_value =
              imaginary_evaluator->EvaluateClosure(side.element, side.point);
          for (int d = 0; d < 3; d++)
          {
            electric_real[d] += average_scale * real_value.value[d];
            electric_imaginary[d] += average_scale * imaginary_value.value[d];
          }
        }
        else
        {
          const auto real_value =
              real_gradient_evaluator->EvaluateClosure(side.element, side.point);
          const auto imaginary_value =
              imaginary_gradient_evaluator->EvaluateClosure(side.element, side.point);
          for (int d = 0; d < 3; d++)
          {
            electric_real[d] -= average_scale * real_value.gradient[d];
            electric_imaginary[d] -= average_scale * imaginary_value.gradient[d];
          }
        }
      }

      boundary_transformation->SetIntPoint(&point);
      mfem::Vector normal_vector(3);
      BdrGridFunctionCoefficient::GetNormal(*boundary_transformation, normal_vector,
                                            inverted);
      const fem::singular::Vector3 normal{normal_vector[0], normal_vector[1],
                                          normal_vector[2]};
      const double real_normal = Dot(electric_real, normal);
      const double imaginary_normal = Dot(electric_imaginary, normal);
      const double normal_product =
          real_normal * real_normal + imaginary_normal * imaginary_normal;
      const double field_product =
          Dot(electric_real, electric_real) + Dot(electric_imaginary, electric_imaginary);

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
        const auto &permittivity = material.GetPermittivityReal(sides.front().attribute);
        permittivity.Mult(electric_real_vector, displacement_real_vector);
        permittivity.Mult(electric_imaginary_vector, displacement_imaginary_vector);
        const double displacement_normal_real = Dot(displacement_real, normal);
        const double displacement_normal_imaginary = Dot(displacement_imaginary, normal);
        surface_density = 0.5 * (interface.thickness / interface.permittivity) *
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
      return boundary_transformation->Weight() * surface_density;
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
      const auto standard_segment_rule =
          fem::singular::BuildWeightedSegmentQuadrature(order, 0.0, 0.0);
      long double value = 0.0L;
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
          long double result = 0.0L;
          for (const auto &tangent_quadrature : standard_segment_rule)
          {
            const double tangent = interval_lower + (interval_upper - interval_lower) *
                                                        tangent_quadrature.coordinate;
            const double maximum_radius =
                GetFaceCutoffRay(*boundary_transformation, physical_cutoff_edges,
                                 boundary_edge, tangent, interface.edge_cutoff)
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
              const double contribution = radial_quadrature.weight *
                                          tangent_quadrature.weight *
                                          (interval_upper - interval_lower) * dr_du *
                                          radius / 3.0 * density(lambda);
              if (!std::isfinite(contribution))
              {
                throw std::runtime_error(
                    "Tetrahedral curved-face multi-edge cutoff quadrature produced "
                    "a nonfinite contribution!");
              }
              result += contribution;
            }
          }
          return result;
        };
        const auto integrate_tangent_interval =
            [&](auto &&self, int boundary_edge, int first, int second, double lower,
                double upper, long double coarse, int depth) -> long double
        {
          const double midpoint = 0.5 * (lower + upper);
          const long double left =
              integrate_interval(boundary_edge, first, second, lower, midpoint);
          const long double right =
              integrate_interval(boundary_edge, first, second, midpoint, upper);
          const long double fine = left + right;
          const double scale = std::max(std::abs(static_cast<double>(fine)),
                                        std::abs(static_cast<double>(coarse)));
          const double error = 8.0 * std::abs(static_cast<double>(fine - coarse)) +
                               256.0 * std::numeric_limits<double>::epsilon() * scale;
          const double interval_fraction = upper - lower;
          const double tolerance =
              0.125 * (options.absolute_tolerance * interval_fraction / 3.0 +
                       options.relative_tolerance * std::abs(static_cast<double>(fine)));
          if (std::isfinite(error) && error <= tolerance)
          {
            return fine;
          }
          if (depth >= options.maximum_subdivisions)
          {
            throw std::runtime_error(fmt::format(
                "Tetrahedral curved-face multi-edge cutoff tangent quadrature did "
                "not meet tolerance: value = {:.17g}, estimated absolute error = "
                "{:.17g}, tolerance = {:.17g}, interval = [{:.17g}, {:.17g}], "
                "order = {}, subdivisions = {}!",
                static_cast<double>(fine), error, tolerance, lower, upper, order,
                options.maximum_subdivisions));
          }
          return self(self, boundary_edge, first, second, lower, midpoint, left,
                      depth + 1) +
                 self(self, boundary_edge, first, second, midpoint, upper, right,
                      depth + 1);
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
            const long double coarse =
                integrate_interval(boundary_edge, first, second, lower, upper);
            value += integrate_tangent_interval(integrate_tangent_interval, boundary_edge,
                                                first, second, lower, upper, coarse, 0);
          }
        }
        return static_cast<double>(value);
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
            value += radial.weight * tangent.weight * r * density(lambda);
          }
        }
        return static_cast<double>(value);
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
                  *boundary_transformation, feature.nodes, t, interface.edge_cutoff);
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
                const double contribution = radial.weight * tangent.weight * dt_du *
                                            (1.0 - r) * dr_du * density(lambda);
                if (!std::isfinite(contribution))
                {
                  throw std::runtime_error(
                      "Tetrahedral curved-face cutoff quadrature produced a nonfinite "
                      "contribution!");
                }
                value += contribution;
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
            const double contribution = radial_weight * tangent_weight * chart_jacobian *
                                        partition_weight(selected, lambda) *
                                        density(lambda) / radial_weight_function;
            if (!std::isfinite(contribution))
            {
              throw std::runtime_error(
                  "Tetrahedral singular surface quadrature produced a nonfinite "
                  "contribution!");
            }
            value += contribution;
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
      return static_cast<double>(value);
    };

    constexpr int order_increment = 4;
    int comparison_order = std::max(4, 2 * options.quadrature_order);
    double comparison = integrate(comparison_order);
    double value = comparison;
    double error = std::numeric_limits<double>::infinity();
    double tolerance = 0.0;
    int order = comparison_order + order_increment;
    bool converged = false;
    for (int refinement = 0; refinement <= options.maximum_subdivisions; refinement++)
    {
      value = integrate(order);
      const double scale = std::max({1.0, std::abs(value), std::abs(comparison)});
      error = 8.0 * std::abs(value - comparison) +
              256.0 * std::numeric_limits<double>::epsilon() * scale;
      tolerance = options.absolute_tolerance + options.relative_tolerance * std::abs(value);
      if (std::isfinite(value) && std::isfinite(comparison) && std::isfinite(error) &&
          error <= tolerance)
      {
        converged = true;
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
          "tolerance = {:.17g}, orders = {}/{}, order refinements = {}!",
          value, comparison, error, tolerance, order, comparison_order,
          options.maximum_subdivisions));
    }
    local_energy += value;
  }
  double energy = static_cast<double>(local_energy);
  Mpi::GlobalSum(1, &energy, mesh.GetComm());
  return energy;
}

std::vector<TetrahedronSingularSurfacePostOperator::Measurement>
TetrahedronSingularSurfacePostOperator::Measure(
    fem::singular::EnrichedNDFieldEvaluator &real_evaluator,
    fem::singular::EnrichedNDFieldEvaluator &imaginary_evaluator,
    double total_electric_energy,
    const fem::singular::AdaptiveAssemblyOptions &options) const
{
  if (!std::isfinite(total_electric_energy) || !(total_electric_energy > 0.0))
  {
    throw std::invalid_argument(
        "Tetrahedral singular surface participation requires positive total electric "
        "energy!");
  }
  std::vector<Measurement> measurements;
  measurements.reserve(interfaces.size());
  for (const auto &[index, interface] : interfaces)
  {
    const double energy = IntegrateInterface(
        interface, &real_evaluator, &imaginary_evaluator, nullptr, nullptr, options);
    const double participation = energy / total_electric_energy;
    const double quality_factor = participation == 0.0 || interface.loss_tangent == 0.0
                                      ? mfem::infinity()
                                      : 1.0 / (participation * interface.loss_tangent);
    measurements.push_back(
        {index, energy, interface.loss_tangent, participation, quality_factor});
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
  std::vector<Measurement> measurements;
  measurements.reserve(interfaces.size());
  for (const auto &[index, interface] : interfaces)
  {
    const double energy = IntegrateInterface(interface, nullptr, nullptr, &real_evaluator,
                                             &imaginary_evaluator, options);
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
