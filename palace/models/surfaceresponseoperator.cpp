// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponseoperator.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
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
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/interpolator.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/cornertracebasis.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfaceresponseidentification.hpp"
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

  bool Contains(const ElementBox &other, int dimension) const
  {
    double scale = 1.0;
    for (int d = 0; d < dimension; d++)
    {
      scale = std::max({scale, std::abs(min[d]), std::abs(max[d]), std::abs(other.min[d]),
                        std::abs(other.max[d])});
    }
    const double tolerance = 1.0e-12 * scale;
    for (int d = 0; d < dimension; d++)
    {
      if (other.min[d] < min[d] - tolerance || other.max[d] > max[d] + tolerance)
      {
        return false;
      }
    }
    return true;
  }

  bool InteriorOverlaps(const ElementBox &other, int dimension) const
  {
    double scale = 1.0;
    for (int d = 0; d < dimension; d++)
    {
      scale = std::max({scale, std::abs(min[d]), std::abs(max[d]), std::abs(other.min[d]),
                        std::abs(other.max[d])});
    }
    const double tolerance = 1.0e-12 * scale;
    for (int d = 0; d < dimension; d++)
    {
      if (std::min(max[d], other.max[d]) - std::max(min[d], other.min[d]) <= tolerance)
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

  // Distance from a point to the nearest local element bounding box (0 inside one): a
  // lower bound of the distance to the local mesh, exact for an axis-aligned domain face.
  double BoxDistance(const std::array<double, 3> &point) const
  {
    auto Distance = [&](const ElementBox &box)
    {
      double distance2 = 0.0;
      for (int d = 0; d < dimension; d++)
      {
        const double excess = std::max({box.min[d] - point[d], point[d] - box.max[d], 0.0});
        distance2 += excess * excess;
      }
      return std::sqrt(distance2);
    };
    double best = mfem::infinity();
    std::vector<int> stack = {0};
    while (!stack.empty())
    {
      const auto &node = nodes[stack.back()];
      stack.pop_back();
      if (Distance(node.box) >= best)
      {
        continue;
      }
      if (node.IsLeaf())
      {
        for (std::size_t i = node.begin; i < node.end; i++)
        {
          best = std::min(best, Distance(element_boxes[indices[i]]));
        }
        continue;
      }
      // Nearer child last (visited first).
      const bool left_nearer =
          Distance(nodes[node.left].box) <= Distance(nodes[node.right].box);
      stack.push_back(left_nearer ? node.right : node.left);
      stack.push_back(left_nearer ? node.left : node.right);
    }
    return best;
  }

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

  // The largest extent of a local element bounding box along any axis.
  double MaxElementExtent() const
  {
    double extent = 0.0;
    for (const auto &box : element_boxes)
    {
      for (int d = 0; d < dimension; d++)
      {
        extent = std::max(extent, box.max[d] - box.min[d]);
      }
    }
    return extent;
  }
};

}  // namespace

// Point location in the distributed device mesh without a global search structure
// (decision 346 (b)): FindPointsGSLIB's Setup builds a uniform global hash of the mesh
// (hash_n = cbrt(32 N_tets)) whose crystal-router distribution overflows a 32-bit message
// at ~17-28 M tets of a graded AMR mesh. Here the ranks' routing boxes are gathered once
// per construction (one per AMR cycle, shared by the mortar-resolution probe and the
// response-point location); a query goes to the ranks whose boxes contain it, is searched
// in their rank-local meshes by the operator's ElementPointLocator, and is owned by the
// lowest rank / lowest element containing it. A point no candidate rank finds this way is
// searched by the candidate ranks with a rank-local FindPointsGSLIB (MPI_COMM_SELF: the
// hash of the local mesh only, no router across ranks) with gslib's border tolerance.
class DistributedPointLocator
{
public:
  struct Result
  {
    // Per query point: the owning rank (the communicator size when none was found), its
    // local element (-1 when none) and the reference coordinates (dimension per point).
    std::vector<int> owners;
    std::vector<int> elements;
    std::vector<double> references;
    // The owner's value of the optional per-element evaluation.
    std::vector<double> owner_values;
    long long int candidate_queries = 0;
    long long int fallback_queries = 0;
  };

private:
  static constexpr int routing_box_count = 8;
  static constexpr int routing_box_values = 6;

  mfem::ParMesh &mesh;
  int dimension;
  MPI_Comm comm;
  int size;
  ElementPointLocator local;
  std::vector<double> global_routing;
  double box_tolerance;
  // The margin by which a rank-local FindPointsGSLIB search may find a point outside the
  // rank's routing boxes: gslib expands every element bounding box by the relative slop
  // bb_t = 0.01 of its size and accepts a border point within sqrt(bdr_tol) = 1e-4 mesh
  // units (the maximum over the ranks).
  double fallback_margin;

  bool RankContains(int candidate_rank, const std::array<double, 3> &point,
                    double tolerance) const
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
      if (bounds.Contains(point, dimension, tolerance))
      {
        return true;
      }
    }
    return false;
  }

  static int SetOffsets(const std::vector<int> &counts, std::vector<int> &offsets)
  {
    offsets.resize(counts.size());
    int total = 0;
    for (std::size_t i = 0; i < counts.size(); i++)
    {
      offsets[i] = total;
      total += counts[i];
    }
    return total;
  }

  static std::vector<int> ScaleCommunicationPlan(const std::vector<int> &values, int scale)
  {
    std::vector<int> result(values);
    for (auto &value : result)
    {
      value *= scale;
    }
    return result;
  }

  // The candidate exchange: every query point of the local list `points` (3 coordinates
  // each) is sent to the ranks whose routing boxes contain it within `tolerance`; the
  // candidate ranks answer with `answer_size` values per query (filled by `Answer` from
  // the received coordinates), gathered here per (point, candidate rank).
  struct CandidateExchange
  {
    std::vector<int> send_counts, send_offsets;
    std::vector<int> query_indices;  // the local point of each packed query
    std::vector<double> answers;     // answer_size values per packed query
  };
  template <typename Answer>
  CandidateExchange Exchange(const std::vector<std::array<double, 3>> &points,
                             double tolerance, int answer_size, Answer &&answer) const
  {
    CandidateExchange exchange;
    exchange.send_counts.assign(size, 0);
    for (const auto &point : points)
    {
      for (int candidate_rank = 0; candidate_rank < size; candidate_rank++)
      {
        if (RankContains(candidate_rank, point, tolerance))
        {
          exchange.send_counts[candidate_rank]++;
        }
      }
    }
    std::vector<int> receive_counts(size);
    Mpi::Alltoall(1, exchange.send_counts.data(), receive_counts.data(), comm);
    std::vector<int> receive_offsets;
    const int send_total = SetOffsets(exchange.send_counts, exchange.send_offsets);
    const int receive_total = SetOffsets(receive_counts, receive_offsets);

    exchange.query_indices.resize(send_total);
    std::vector<double> send_coordinates(dimension * send_total);
    std::vector<int> cursor(exchange.send_offsets);
    for (std::size_t point = 0; point < points.size(); point++)
    {
      for (int candidate_rank = 0; candidate_rank < size; candidate_rank++)
      {
        if (!RankContains(candidate_rank, points[point], tolerance))
        {
          continue;
        }
        const int packed = cursor[candidate_rank]++;
        exchange.query_indices[packed] = static_cast<int>(point);
        for (int d = 0; d < dimension; d++)
        {
          send_coordinates[dimension * packed + d] = points[point][d];
        }
      }
    }
    const auto send_coordinate_counts =
        ScaleCommunicationPlan(exchange.send_counts, dimension);
    const auto send_coordinate_offsets =
        ScaleCommunicationPlan(exchange.send_offsets, dimension);
    const auto receive_coordinate_counts =
        ScaleCommunicationPlan(receive_counts, dimension);
    const auto receive_coordinate_offsets =
        ScaleCommunicationPlan(receive_offsets, dimension);
    std::vector<double> receive_coordinates(dimension * receive_total);
    Mpi::Alltoallv(send_coordinates.data(), send_coordinate_counts.data(),
                   send_coordinate_offsets.data(), receive_coordinates.data(),
                   receive_coordinate_counts.data(), receive_coordinate_offsets.data(),
                   comm);

    std::vector<double> receive_answers(answer_size * receive_total);
    answer(receive_coordinates, receive_total, receive_answers);

    const auto send_answer_counts =
        ScaleCommunicationPlan(exchange.send_counts, answer_size);
    const auto send_answer_offsets =
        ScaleCommunicationPlan(exchange.send_offsets, answer_size);
    const auto receive_answer_counts = ScaleCommunicationPlan(receive_counts, answer_size);
    const auto receive_answer_offsets =
        ScaleCommunicationPlan(receive_offsets, answer_size);
    exchange.answers.resize(answer_size * send_total);
    Mpi::Alltoallv(receive_answers.data(), receive_answer_counts.data(),
                   receive_answer_offsets.data(), exchange.answers.data(),
                   send_answer_counts.data(), send_answer_offsets.data(), comm);
    return exchange;
  }

public:
  DistributedPointLocator(mfem::ParMesh &mesh_, int dimension_)
    : mesh(mesh_), dimension(dimension_), comm(mesh.GetComm()), size(Mpi::Size(comm)),
      local(mesh, dimension)
  {
    std::array<double, routing_box_count * routing_box_values> local_routing;
    for (int box = 0; box < routing_box_count; box++)
    {
      for (int d = 0; d < 3; d++)
      {
        local_routing[routing_box_values * box + d] = mfem::infinity();
        local_routing[routing_box_values * box + 3 + d] = -mfem::infinity();
      }
    }
    const auto routing_boxes = local.GetRoutingBoxes(routing_box_count);
    for (std::size_t box = 0; box < routing_boxes.size(); box++)
    {
      for (int d = 0; d < 3; d++)
      {
        local_routing[routing_box_values * box + d] = routing_boxes[box].min[d];
        local_routing[routing_box_values * box + 3 + d] = routing_boxes[box].max[d];
      }
    }
    global_routing.resize(size * local_routing.size());
    Mpi::Allgather(static_cast<int>(local_routing.size()), local_routing.data(),
                   global_routing.data(), comm);

    double coordinate_scale = 0.0;
    const auto &bounds = local.GetBounds();
    for (int d = 0; d < dimension; d++)
    {
      coordinate_scale = std::max({coordinate_scale, std::abs(bounds.min[d]),
                                   std::abs(bounds.max[d]), bounds.max[d] - bounds.min[d]});
    }
    Mpi::GlobalMax(1, &coordinate_scale, comm);
    box_tolerance =
        1.0e-11 * coordinate_scale +
        64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, coordinate_scale);

    fallback_margin = 0.01 * local.MaxElementExtent() + 1.0e-4;
    Mpi::GlobalMax(1, &fallback_margin, comm);
  }

  // Locates the local list of points (byNODES: every x, then every y, then every z). With
  // `element_value`, the owner evaluates it on the owning element.
  Result Locate(const mfem::Vector &xyz,
                const std::function<double(int element)> *element_value = nullptr)
  {
    MFEM_VERIFY(xyz.Size() % dimension == 0, "Invalid point-coordinate array!");
    const int point_count = xyz.Size() / dimension;
    std::vector<std::array<double, 3>> points(point_count);
    for (int point = 0; point < point_count; point++)
    {
      points[point].fill(0.0);
      for (int d = 0; d < dimension; d++)
      {
        points[point][d] = xyz(d * point_count + point);
      }
    }

    Result result;
    result.owners.assign(point_count, size);
    result.elements.assign(point_count, -1);
    result.references.assign(dimension * point_count, 0.0);
    result.owner_values.assign(element_value ? point_count : 0, 0.0);
    const int value_size = element_value ? 1 : 0;

    // The answer of a candidate rank: the local element (-1 when not found), the reference
    // coordinates and the optional element value.
    const int answer_size = 1 + dimension + value_size;
    const auto exchange = Exchange(
        points, box_tolerance, answer_size,
        [&](const std::vector<double> &coordinates, int count, std::vector<double> &answers)
        {
          std::vector<int> candidates;
          for (int i = 0; i < count; i++)
          {
            std::array<double, 3> coordinate{};
            for (int d = 0; d < dimension; d++)
            {
              coordinate[d] = coordinates[dimension * i + d];
            }
            int element = -1;
            mfem::IntegrationPoint reference;
            double *answer = answers.data() + answer_size * i;
            if (local.Find(coordinate, box_tolerance, element, reference, candidates))
            {
              reference.Get(answer + 1, dimension);
              if (element_value)
              {
                answer[1 + dimension] = (*element_value)(element);
              }
            }
            answer[0] = element;
          }
        });
    result.candidate_queries = static_cast<long long int>(exchange.query_indices.size());
    for (int candidate_rank = 0; candidate_rank < size; candidate_rank++)
    {
      const int begin = exchange.send_offsets[candidate_rank];
      const int end = begin + exchange.send_counts[candidate_rank];
      for (int packed = begin; packed < end; packed++)
      {
        const double *answer = exchange.answers.data() + answer_size * packed;
        const int element = static_cast<int>(answer[0]);
        if (element < 0)
        {
          continue;
        }
        const int point = exchange.query_indices[packed];
        if (candidate_rank > result.owners[point] ||
            (candidate_rank == result.owners[point] && element >= result.elements[point]))
        {
          continue;
        }
        result.owners[point] = candidate_rank;
        result.elements[point] = element;
        std::copy_n(answer + 1, dimension, result.references.data() + dimension * point);
        if (element_value)
        {
          result.owner_values[point] = answer[1 + dimension];
        }
      }
    }

    std::vector<int> fallback_indices;
    for (int point = 0; point < point_count; point++)
    {
      if (result.owners[point] == size)
      {
        fallback_indices.push_back(point);
      }
    }
    result.fallback_queries = static_cast<long long int>(fallback_indices.size());
    int fallback_count = static_cast<int>(fallback_indices.size());
    Mpi::GlobalSum(1, &fallback_count, comm);
    if (fallback_count == 0)
    {
      return result;
    }
    Mpi::Warning(comm,
                 "Distributed point location could not resolve {:d} points in the local "
                 "meshes of their candidate ranks; searching them with a rank-local "
                 "FindPointsGSLIB!\n",
                 fallback_count);

    // Every candidate rank (its routing boxes within its own fallback margin) searches the
    // unresolved points in its local mesh with gslib; the answer carries gslib's code
    // (0 inside, 1 on a border within the tolerance, 2 not found), the local element, the
    // squared distance, the reference coordinates and the optional element value. The
    // owner is chosen as gslib's global search does: an inside point before a border one,
    // then the smaller distance, then the lower rank.
    std::vector<std::array<double, 3>> fallback_points(fallback_indices.size());
    for (std::size_t i = 0; i < fallback_indices.size(); i++)
    {
      fallback_points[i] = points[fallback_indices[i]];
    }
    const int fallback_answer_size = 3 + dimension + value_size;
    const auto fallback = Exchange(
        fallback_points, box_tolerance + fallback_margin, fallback_answer_size,
        [&](const std::vector<double> &coordinates, int count, std::vector<double> &answers)
        {
          for (int i = 0; i < count; i++)
          {
            double *answer = answers.data() + fallback_answer_size * i;
            answer[0] = 2.0;
            answer[1] = -1.0;
            answer[2] = mfem::infinity();
          }
          if (count == 0)
          {
            return;
          }
#if defined(MFEM_USE_GSLIB)
          mfem::Vector fallback_xyz(dimension * count);
          for (int i = 0; i < count; i++)
          {
            for (int d = 0; d < dimension; d++)
            {
              fallback_xyz(d * count + i) = coordinates[dimension * i + d];
            }
          }
          mfem::FindPointsGSLIB finder(MPI_COMM_SELF);
          finder.Setup(mesh, 0.01, 1.0e-12, 256);
          finder.FindPoints(fallback_xyz, mfem::Ordering::byNODES);
          const auto &reference = finder.GetReferencePosition();
          for (int i = 0; i < count; i++)
          {
            if (finder.GetCode()[i] == 2)
            {
              continue;
            }
            double *answer = answers.data() + fallback_answer_size * i;
            const int element = static_cast<int>(finder.GetElem()[i]);
            answer[0] = finder.GetCode()[i];
            answer[1] = element;
            answer[2] = finder.GetDist()(i);
            for (int d = 0; d < dimension; d++)
            {
              answer[3 + d] = reference(dimension * i + d);
            }
            if (element_value)
            {
              answer[3 + dimension] = (*element_value)(element);
            }
          }
#else
          MFEM_ABORT("Rank-local point location fallback requires MFEM_USE_GSLIB!");
#endif
        });
    std::vector<double> best_code(fallback_indices.size(), 2.0);
    std::vector<double> best_distance(fallback_indices.size(), mfem::infinity());
    for (int candidate_rank = 0; candidate_rank < size; candidate_rank++)
    {
      const int begin = fallback.send_offsets[candidate_rank];
      const int end = begin + fallback.send_counts[candidate_rank];
      for (int packed = begin; packed < end; packed++)
      {
        const double *answer = fallback.answers.data() + fallback_answer_size * packed;
        if (answer[0] == 2.0)
        {
          continue;
        }
        const int i = fallback.query_indices[packed];
        if (answer[0] > best_code[i] ||
            (answer[0] == best_code[i] && answer[2] >= best_distance[i]))
        {
          continue;
        }
        best_code[i] = answer[0];
        best_distance[i] = answer[2];
        const int point = fallback_indices[i];
        result.owners[point] = candidate_rank;
        result.elements[point] = static_cast<int>(answer[1]);
        std::copy_n(answer + 3, dimension, result.references.data() + dimension * point);
        if (element_value)
        {
          result.owner_values[point] = answer[3 + dimension];
        }
      }
    }
    return result;
  }
};

namespace
{

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
  JUNCTION,
  // Curved classes of the identification (decision 73(1)(c)): models keyed by their
  // version-2 Signature (RadiusOverR), patched like their straight analogues along the
  // curved portions.
  CURVED_EDGE,
  CURVED_SAME_CONDUCTOR_GAP,
  CURVED_DIFFERENT_CONDUCTOR_GAP,
  CURVED_SAME_CONDUCTOR_STRIP
};

bool IsCurvedTopology(LibraryTopology topology)
{
  return topology == LibraryTopology::CURVED_EDGE ||
         topology == LibraryTopology::CURVED_SAME_CONDUCTOR_GAP ||
         topology == LibraryTopology::CURVED_DIFFERENT_CONDUCTOR_GAP ||
         topology == LibraryTopology::CURVED_SAME_CONDUCTOR_STRIP;
}

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
  std::vector<Point3D> support_points;
  double spatial_position_tolerance = 0.0;
  double spatial_angle_tolerance = 0.0;
  bool boundary_law_physics_qualified = true;
  std::optional<std::string> plan_view_boundary;
  std::optional<std::string> mask_regularization;
  std::vector<LibraryInterface> interfaces;

  // Version-2 identification signature (the feature's canonical Signature object,
  // dimensionless parameters) which the key-based matching pass looks up directly.
  std::optional<nlohmann::json> identification_signature;
  // Legacy-contract aliases (USER decision 283): contract-3 keys this legacy model serves.
  std::vector<LegacyContractAlias> legacy_contract_aliases;

  // Curvature family node (curved coupons built on axisymmetric (r, z) meshes): the
  // coupon's kappa = R / rho and whether the metal lies inside the circle (convex, a disk
  // edge) or outside (concave, a hole edge). A straight coupon of the analogous topology is
  // the kappa = 0 anchor of both convexities.
  std::optional<double> kappa;
  std::optional<bool> convex;

  // Corner coupon built on the angle-independent trace basis rule (record TraceBasis of
  // generate_corner_response.py): required of every node of an angle-interpolated corner
  // family (MatchCornerFamily), which checks the coupon's files against the rule at its
  // angle and constructs the runtime basis at the device angle from it.
  std::optional<CornerTraceBasisRule> trace_basis;
  // The segment connectivity of a corner family coupon (TraceBasis
  // ConnectivityAngleDegrees, corner-qualification block 2026-09-29): the angle whose rule
  // layout orders the band triangulation next to the metal rings; coupons sharing it form
  // one interpolation segment. Absent on a legacy coupon (perimeter order at its own
  // angle): exact matches only.
  std::optional<double> corner_connectivity_angle_degrees;
};

struct ProcessLibrary
{
  int version = 0;
  int trace_lift_version = 0;
  bool exhaustive_spatial_closure = false;
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

  // The supporting faces have the same material on both sides and the target configures
  // no EdgeFrameNormal: the process side of axis_v is a guess the identification excludes.
  bool ambiguous_process_side = false;
};

// Two legacy segments count as parallel up to a cosine deficit of 1e-8, i.e. an angle of
// sqrt(2e-8) rad, so the overlap of the first projected back onto the second may overshoot
// the second by that angle times the projected distances (DS-SCT-001: the two edges of a
// CPW gap along a 250 um bend at 4 um chords differ by 3e-5 to 1e-4 rad); the callers clamp
// the interval to the second segment. Anything beyond that bound is a genuine
// inconsistency.
void VerifyParallelOverlap(const EdgeSegment3D &first, const EdgeSegment3D &second,
                           double tangent_dot, double second_begin, double second_end,
                           double tolerance, double interaction_distance)
{
  const double overlap_tolerance =
      tolerance + std::sqrt(2.0e-8) * (first.length + second.length + interaction_distance);
  MFEM_VERIFY(second_begin >= -overlap_tolerance &&
                  second_end <= second.length + overlap_tolerance,
              "Inconsistent overlap between nearby parallel metal edges (second segment "
              "interval ["
                  << second_begin << ", " << second_end << "] of length " << second.length
                  << ", tangent dot " << tangent_dot << ")!");
}

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

// A runtime model interpolated in kappa = R / rho between the coupons of a curvature family
// (decision 92): the anchor (straight coupon, kappa 0) supplies the basis, interfaces and
// conductor references; the matrices are the Lagrange combination of the nodes' matrices,
// each scaled to the anchor's coupon depth (per unit edge length).
struct PendingBlend
{
  std::string name;
  std::string topology;
  std::size_t anchor = 0;
  std::vector<LibrarySelection::WeightedModel> nodes;
};

struct PendingPatch
{
  std::size_t library_model = 0;
  ResponsePatchData patch;
  std::optional<PendingBlend> blend;
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
  if (topology == "CurvedEdge")
  {
    return LibraryTopology::CURVED_EDGE;
  }
  if (topology == "CurvedSameConductorGap")
  {
    return LibraryTopology::CURVED_SAME_CONDUCTOR_GAP;
  }
  if (topology == "CurvedDifferentConductorGap")
  {
    return LibraryTopology::CURVED_DIFFERENT_CONDUCTOR_GAP;
  }
  if (topology == "CurvedSameConductorStrip")
  {
    return LibraryTopology::CURVED_SAME_CONDUCTOR_STRIP;
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

std::string TopologyIdentifier(LibraryTopology topology);
void VerifySpatialEdgesInSignatureFrame(const LibraryModel &model, double radius);

std::vector<std::array<double, 3>> ReadBasisPoints(const std::string &path);
std::vector<std::string> ModelInterfaceNames(const LibraryModel &model);
// The load-time check of an AllRingsFollowMetal corner coupon's files against the rule at
// ITS angle (corner-basis refinement 2026-09-30; fail closed in ReadProcessLibrary): the
// outer ring levels {-R, -R/3, -OveretchDepth, 0, MetalThickness, MetalThickness + k
// OveretchDepth (the rule's k), R/3, R}, every basis point at the rule's position, every
// knot's zero flag equal to its ZeroTraceIndices membership and the
// trace mesh (vertices, slave parents, triangle set) equal to the rule's. (A MetalRingsOnly
// coupon with a segment connectivity is checked the same way at match time,
// MatchCornerFamily.) Returns the reason, empty when the coupon is the rule's.
std::string CheckCornerRuleCouponFiles(const LibraryModel &model,
                                       double position_tolerance);

ProcessLibrary ReadProcessLibrary(const std::string &path, const Units &units,
                                  bool nondimensionalize, bool allow_empty_models = false,
                                  bool geometry_only = false)
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
  library.trace_lift_version = data.value("TraceLiftVersion", 0);
  library.exhaustive_spatial_closure = data.value("ExhaustiveSpatialClosure", false);
  MFEM_VERIFY(library.trace_lift_version >= 0,
              "Fabrication-process response-library TraceLiftVersion must be "
              "nonnegative!");
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
  std::set<std::string> legacy_alias_keys;
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
    // Corner models record Angle (degrees; the corner family's records also carry
    // AngleDegrees, accepted as the same value).
    model.angle =
        entry.value("Angle", entry.value("AngleDegrees", 0.0)) * std::acos(-1.0) / 180.0;
    model.angle_tolerance = entry.value("AngleTolerance", 0.0) * std::acos(-1.0) / 180.0;
    model.corner_radius = entry.value("CornerRadius", 0.0) / coordinate_scale;
    model.corner_radius_tolerance =
        entry.value("CornerRadiusTolerance", 0.0) / coordinate_scale;
    if (auto signature = entry.find("Signature"); signature != entry.end())
    {
      MFEM_VERIFY(
          signature->is_object() && signature->contains("Type"),
          "Fabrication-process response model Signature must be the identification's "
          "canonical Signature object (with Type)!");
      MFEM_VERIFY(signature->at("Type").get<std::string>() ==
                      TopologyIdentifier(model.topology),
                  "Fabrication-process response model \""
                      << model.name << "\" Signature.Type does not match its Topology!");
      // An UnboxableFeature key (decision 282: no box satisfies the face rules) is a
      // Missing placeholder no builder makes: a library may not serve it.
      MFEM_VERIFY(!signature->value("Unboxable", false),
                  "Fabrication-process response model \""
                      << model.name
                      << "\" is keyed by an UnboxableFeature signature (\"Unboxable\": "
                         "true): no coupon exists for such a key!");
      model.identification_signature = *signature;
    }
    if (auto aliases = entry.find("LegacyContractAliases"); aliases != entry.end())
    {
      MFEM_VERIFY(aliases->is_array(), "Fabrication-process response model \""
                                           << model.name
                                           << "\" LegacyContractAliases must be an array!");
      MFEM_VERIFY(model.topology == LibraryTopology::SPATIAL_EDGE_CLUSTER &&
                      model.identification_signature,
                  "Fabrication-process response model \""
                      << model.name
                      << "\" lists LegacyContractAliases but is not a SpatialEdgeCluster "
                         "model with a Signature (USER decision 283)!");
      const std::string own_hash =
          SignatureKeyAndHash(*model.identification_signature, "SpatialEdgeCluster").second;
      auto Hex64 = [](const std::string &text)
      {
        return text.size() == 64 &&
               std::all_of(text.begin(), text.end(), [](char c)
                           { return std::isxdigit(static_cast<unsigned char>(c)); });
      };
      for (const auto &alias_entry : *aliases)
      {
        LegacyContractAlias alias;
        alias.key = alias_entry.value("Key", std::string{});
        alias.context_digest = alias_entry.value("ContextDigest", std::string{});
        alias.reason = alias_entry.value("Reason", std::string{});
        alias.context = alias_entry.value("Context", nlohmann::json(nullptr));
        MFEM_VERIFY(Hex64(alias.key) && Hex64(alias.context_digest) &&
                        !alias.reason.empty(),
                    "Fabrication-process response model \""
                        << model.name
                        << "\" LegacyContractAliases entries need a 64-hex Key, a 64-hex "
                           "ContextDigest and a Reason!");
        MFEM_VERIFY(alias.key != own_hash,
                    "Fabrication-process response model \""
                        << model.name
                        << "\" lists its own key as a legacy-contract alias!");
        MFEM_VERIFY(legacy_alias_keys.insert(alias.key).second,
                    "Fabrication-process response library lists the legacy-contract alias "
                    "key "
                        << alias.key << " twice (model \"" << model.name << "\")!");
        model.legacy_contract_aliases.push_back(std::move(alias));
      }
    }
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
    if (auto kappa = entry.find("Kappa"); kappa != entry.end())
    {
      MFEM_VERIFY(IsCurvedTopology(model.topology),
                  "Fabrication-process response model \""
                      << model.name << "\" carries Kappa but is not a curved topology!");
      model.kappa = kappa->get<double>();
      MFEM_VERIFY(std::isfinite(*model.kappa) && *model.kappa > 0.0 && *model.kappa < 1.0,
                  "Fabrication-process response model Kappa = R / rho must lie in (0, 1)!");
      const std::string convexity = entry.at("Convexity").get<std::string>();
      MFEM_VERIFY(
          convexity == "Convex" || convexity == "Concave",
          "Fabrication-process response model Convexity must be Convex or Concave!");
      model.convex = convexity == "Convex";
      MFEM_VERIFY(model.coupon_depth > 0.0,
                  "Fabrication-process response model \""
                      << model.name
                      << "\" requires CouponDepth (the curved edge length 2 pi rho)!");
    }
    const bool parallel_cluster = model.topology == LibraryTopology::PARALLEL_EDGE_CLUSTER;
    const bool spatial_cluster = model.topology == LibraryTopology::SPATIAL_EDGE_CLUSTER;
    const bool trace_mesh_topology = spatial_cluster ||
                                     model.topology == LibraryTopology::CONVEX_CORNER ||
                                     model.topology == LibraryTopology::CONCAVE_CORNER ||
                                     model.topology == LibraryTopology::ENDPOINT ||
                                     model.topology == LibraryTopology::JUNCTION;
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
        if (auto support = entry.find("SupportPoints"); support != entry.end())
        {
          MFEM_VERIFY(support->is_array() && support->size() == 8,
                      "SpatialEdgeCluster SupportPoints must contain the eight corners "
                      "of its matching volume!");
          model.support_points = support->get<std::vector<Point3D>>();
          for (auto &point : model.support_points)
          {
            for (double &value : point)
            {
              MFEM_VERIFY(std::isfinite(value),
                          "SpatialEdgeCluster SupportPoints must be finite!");
              value /= coordinate_scale;
            }
          }
        }
        MFEM_VERIFY(geometry_only || library.trace_lift_version < 2 ||
                        !model.support_points.empty(),
                    "New SpatialEdgeCluster models require explicit SupportPoints "
                    "matching-volume metadata!");
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
    // A curved model has no version-1 geometry parameters: it is keyed by its Signature,
    // or it is a curvature-family node (Kappa + Convexity) interpolated in kappa.
    MFEM_VERIFY(!IsCurvedTopology(model.topology) || model.identification_signature ||
                    model.kappa,
                "Curved fabrication-process response model \""
                    << model.name
                    << "\" requires its identification Signature or a Kappa / Convexity "
                       "curvature-family record!");
    if (model.topology == LibraryTopology::ISOLATED_EDGE || spatial_response ||
        parallel_cluster || model.topology == LibraryTopology::CURVED_EDGE)
    {
      MFEM_VERIFY(model.separation == 0.0,
                  "An isolated-edge, curved-edge, parallel-cluster, or spatial response "
                  "model cannot specify a separation!");
    }
    else if (model.identification_signature && model.separation == 0.0)
    {
      // A paired-edge model keyed by its Signature needs no separation of its own.
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
      // Angle in (0, 180]: 180 is the straight edge through the corner box, the anchor of
      // the angle-interpolated corner family (USER decision 121 (C)); no device corner has
      // that angle (a corner is a joint that is not noise), so the anchor is only ever
      // reached through MatchCornerFamily.
      MFEM_VERIFY(std::isfinite(model.angle) && model.angle > 0.0 &&
                      model.angle <= std::acos(-1.0) + 1.0e-12 &&
                      std::isfinite(model.angle_tolerance) && model.angle_tolerance >= 0.0,
                  "Corner response models require Angle in (0, 180] degrees and a "
                  "nonnegative AngleTolerance!");
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
    // A cluster model keyed by its identification Signature is built in the feature's
    // canonical frame and needs no Edges of its own.
    MFEM_VERIFY(parallel_cluster == !model.cluster_edges.empty() ||
                    (parallel_cluster && model.identification_signature),
                "ParallelEdgeCluster response models require Edges (or a Signature), and "
                "other topologies cannot specify them!");
    MFEM_VERIFY(
        spatial_cluster == !model.spatial_edges.empty() ||
            (spatial_cluster && model.identification_signature),
        "SpatialEdgeCluster response models require spatial Edges (or a Signature), "
        "and other topologies cannot specify them!");
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
    if (auto trace_mesh = entry.find("TraceMesh"); trace_mesh != entry.end())
    {
      MFEM_VERIFY(version >= 3 && trace_mesh->is_object() &&
                      trace_mesh->contains("Vertices") &&
                      trace_mesh->contains("Triangles") && trace_mesh_topology,
                  "TraceMesh requires Vertices and Triangles and is supported by a "
                  "version-3 spatial response model!");
      model.response.trace_vertices =
          ResolveLibraryPath(directory, trace_mesh->at("Vertices").get<std::string>());
      model.response.trace_triangles =
          ResolveLibraryPath(directory, trace_mesh->at("Triangles").get<std::string>());
    }
    if (auto trace_basis = entry.find("TraceBasis"); trace_basis != entry.end())
    {
      MFEM_VERIFY(corner && spatial_response && trace_basis->is_object(),
                  "TraceBasis (the corner family's trace basis rule) is supported only by "
                  "spatial corner response models!");
      CornerTraceBasisRule rule;
      rule.ring_size = trace_basis->at("RingSize").get<int>();
      rule.metal_interior_knots = trace_basis->at("MetalInteriorKnots").get<int>();
      rule.free_knots = trace_basis->at("FreeKnots").get<int>();
      rule.fractions = trace_basis->at("Fractions").get<std::string>();
      if (auto layout = trace_basis->find("RingLayout"); layout != trace_basis->end())
      {
        const std::string name = layout->get<std::string>();
        MFEM_VERIFY(name == "MetalRingsOnly" || name == "AllRingsFollowMetal",
                    "Fabrication-process response model \""
                        << model.name << "\" has an unknown TraceBasis RingLayout \""
                        << name << "\" (MetalRingsOnly | AllRingsFollowMetal)!");
        rule.ring_layout = name == "AllRingsFollowMetal"
                               ? CornerRingLayout::ALL_RINGS_FOLLOW_METAL
                               : CornerRingLayout::METAL_RINGS_ONLY;
      }
      if (auto grading = trace_basis->find("FreeKnotGrading");
          grading != trace_basis->end())
      {
        MFEM_VERIFY(grading->is_array(), "Fabrication-process response model \""
                                             << model.name
                                             << "\" TraceBasis FreeKnotGrading must be an "
                                                "array of distances over R!");
        rule.free_knot_grading = grading->get<std::vector<double>>();
      }
      if (auto reference = trace_basis->find("FreeKnotGradingReferenceFreeArcOverR");
          reference != trace_basis->end())
      {
        // Written by the generator only where the scaling of the graded free knots is
        // active (an acute concave node, decision 318); an absent key is the default
        // reference kFreeKnotGradingReferenceFreeArcOverR, so the records of the family's
        // other nodes are unchanged and the family stays on one rule.
        MFEM_VERIFY(reference->is_number() && reference->get<double>() > 0.0 &&
                        !rule.free_knot_grading.empty(),
                    "Fabrication-process response model \""
                        << model.name
                        << "\" TraceBasis FreeKnotGradingReferenceFreeArcOverR must be a "
                           "positive free arc over R of a FreeKnotGrading layout!");
        rule.free_knot_grading_reference_free_arc_over_r = reference->get<double>();
      }
      if (auto extra = trace_basis->find("ExtraLevelsAboveOverOveretch");
          extra != trace_basis->end())
      {
        MFEM_VERIFY(extra->is_array(),
                    "Fabrication-process response model \""
                        << model.name
                        << "\" TraceBasis ExtraLevelsAboveOverOveretch must "
                           "be an array of multiples of OveretchDepth!");
        rule.extra_levels_above_over_overetch = extra->get<std::vector<double>>();
      }
      {
        const std::string reason = CheckCornerTraceBasisRule(rule);
        MFEM_VERIFY(reason.empty(), "Fabrication-process response model \""
                                        << model.name
                                        << "\" has an invalid TraceBasis rule (" << reason
                                        << ")!");
      }
      model.trace_basis = rule;
      if (auto connectivity = trace_basis->find("ConnectivityAngleDegrees");
          connectivity != trace_basis->end())
      {
        MFEM_VERIFY(!rule.AllRings(),
                    "Fabrication-process response model \""
                        << model.name
                        << "\" carries a TraceBasis ConnectivityAngleDegrees on the "
                           "AllRingsFollowMetal layout, which has no events!");
        MFEM_VERIFY(connectivity->is_number() && connectivity->get<double>() > 0.0 &&
                        connectivity->get<double>() < 180.0,
                    "Fabrication-process response model \""
                        << model.name
                        << "\" has an invalid TraceBasis ConnectivityAngleDegrees (in "
                           "(0, 180))!");
        model.corner_connectivity_angle_degrees = connectivity->get<double>();
      }
    }
    if (auto interior = entry.find("InteriorTraceCount"); interior != entry.end())
    {
      MFEM_VERIFY(interior->is_number_integer() && interior->get<int>() >= 0 &&
                      spatial_response && !model.response.trace_vertices.empty(),
                  "InteriorTraceCount requires a nonnegative count on a spatial response "
                  "model with an explicit TraceMesh (its cap-interior hats are trace "
                  "coefficients only through the trace mesh)!");
      model.response.interior_trace_count = interior->get<int>();
    }
    if (auto references = entry.find("ConductorReferences"); references != entry.end())
    {
      MFEM_VERIFY(version >= 2,
                  "ConductorReferences requires a version-2 fabrication-process "
                  "response library!");
      MFEM_VERIFY(!entry.contains("Reference") && references->is_array() &&
                      references->size() >= 2 &&
                      (model.topology == LibraryTopology::DIFFERENT_CONDUCTOR_GAP ||
                       model.topology == LibraryTopology::CURVED_DIFFERENT_CONDUCTOR_GAP ||
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
      MFEM_VERIFY(geometry_only || library.trace_lift_version >= 2,
                  "Multiconductor fabrication-process response model \""
                      << model.name
                      << "\" was generated without the corrected trace/conductor "
                         "boundary lift. Regenerate the library with TraceLiftVersion "
                         "equal to 2 or newer!");
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
    if ((parallel_cluster && !model.cluster_edges.empty()) ||
        (spatial_cluster && !model.spatial_edges.empty()))
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
                       model.topology == LibraryTopology::CURVED_DIFFERENT_CONDUCTOR_GAP ||
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
    if (corner && spatial_response && !model.response.zero_trace_indices.empty() &&
        !model.response.contour_groups.empty())
    {
      // The basis gate of a corner coupon (corner-family review 2026-09-29; fail closed):
      // every crossing of a metal arm with a box ring that meets the metal is a PEC knot
      // and no free knot lies on the metal, so that no free trace hat has support on the
      // PEC part of the box contour (the review's root cause: such a hat imposes
      // conflicting Dirichlet data on the coupon and puts a near-singular field into the
      // MS / MA layers, 6-7 orders in the fabricated Q_MS of that knot).
      const auto points = ReadBasisPoints(model.response.basis_points);
      const std::string reason = CheckCornerBasisCrossings(
          points, model.response.contour_groups, model.response.zero_trace_indices,
          model.angle, model.topology == LibraryTopology::CONVEX_CORNER,
          1.0e-6 * library.matching_radius * coordinate_scale);
      MFEM_VERIFY(reason.empty(),
                  "Fabrication-process corner response model \""
                      << model.name << "\" fails the trace basis gate: " << reason << "!");
    }
    if (corner && spatial_response && model.trace_basis && model.trace_basis->AllRings() &&
        model.corner_radius == 0.0)
    {
      // The coupon's files against its trace basis rule at its own angle (fail closed at
      // load; corner-basis refinement 2026-09-30): the outer ring levels, the basis points
      // and the trace mesh (the rule fixes all three under AllRingsFollowMetal). The
      // MetalRingsOnly coupons keep their recorded checks: the segment structure at load
      // below, the files at match time (MatchCornerFamily).
      const std::string reason = CheckCornerRuleCouponFiles(
          model, 1.0e-9 * library.matching_radius * coordinate_scale);
      MFEM_VERIFY(reason.empty(),
                  "Fabrication-process corner response model \""
                      << model.name << "\" is not its trace basis rule's coupon: " << reason
                      << "!");
    }

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
        if (spatial_cluster && !model.spatial_edges.empty())
        {
          MFEM_VERIFY(std::any_of(model.spatial_edges.begin(), model.spatial_edges.end(),
                                  [slot](const auto &edge)
                                  { return edge.interface_slot == slot; }),
                      "A SpatialEdgeCluster interface mapping refers to an unused "
                      "InterfaceSlot!");
        }
        else if (spatial_cluster)
        {
          // A Signature-only model's slots are the distinct interface sets of its
          // Signature's portions (slot k = the k-th distinct set in sorted order).
          std::set<std::string> portion_interface_sets;
          for (const auto &portion :
               model.identification_signature->value("Portions", nlohmann::json::array()))
          {
            portion_interface_sets.insert(
                portion.value("Interfaces", nlohmann::json::array()).dump());
          }
          MFEM_VERIFY(slot < static_cast<int>(portion_interface_sets.size()),
                      "A SpatialEdgeCluster interface mapping refers to InterfaceSlot "
                          << slot << " but the model's Signature has only "
                          << portion_interface_sets.size() << " interface slot(s)!");
        }
        model.interfaces.push_back({slot, type, coupon});
        mapped_interface_types.insert(type);
      }
    }
    if (spatial_cluster && model.identification_signature && !model.spatial_edges.empty())
    {
      // The builder's contract: Edges in the canonical frame of the Signature (placed on a
      // feature with the identity map), each carrying its portion's Conductor and the
      // InterfaceSlot mapped to its portion's interface set (after the Interfaces above).
      VerifySpatialEdgesInSignatureFrame(model, library.matching_radius);
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
  // The segment structure of every corner family (decision 137 (1), fail closed at load):
  // the sharp trace-basis-rule corner coupons of one topology, interface set and boundary
  // law (MatchCornerFamily's family) form segments by connectivity angle, each in one
  // event-free interval with its connectivity angle, one coupon per angle in a segment,
  // overlapping at shared node angles only (CheckCornerFamilySegments).
  {
    std::vector<bool> grouped(library.models.size(), false);
    for (std::size_t i = 0; i < library.models.size(); i++)
    {
      const auto &first = library.models[i];
      const bool corner = first.topology == LibraryTopology::CONVEX_CORNER ||
                          first.topology == LibraryTopology::CONCAVE_CORNER;
      if (grouped[i] || !corner || first.corner_radius != 0.0 || !first.trace_basis)
      {
        continue;
      }
      std::vector<CornerFamilyNode> family;
      for (std::size_t j = i; j < library.models.size(); j++)
      {
        const auto &model = library.models[j];
        if (grouped[j] || model.topology != first.topology || model.corner_radius != 0.0 ||
            !model.trace_basis ||
            ModelInterfaceNames(model) != ModelInterfaceNames(first) ||
            !CompatibleBoundaryLaw(model.boundary_condition, first.boundary_condition) ||
            !CompatibleBoundaryLaw(first.boundary_condition, model.boundary_condition))
        {
          continue;
        }
        grouped[j] = true;
        CornerFamilyNode node;
        node.angle_degrees = model.angle * 180.0 / std::acos(-1.0);
        node.connectivity_angle_degrees = model.corner_connectivity_angle_degrees;
        node.index = j;
        family.push_back(node);
      }
      const std::string reason = CheckCornerFamilySegments(
          family, first.topology == LibraryTopology::CONVEX_CORNER, *first.trace_basis,
          kSignatureAngleToleranceDegrees);
      MFEM_VERIFY(reason.empty(), "Fabrication-process response library \""
                                      << library.name << "\" corner family (" << first.name
                                      << " ...): " << reason << "!");
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

// Roundoff-robust comparisons for the classification decisions. The decisions compare a
// geometric distance with the interaction distance 2R (or a multiple of it) and layouts
// routinely place metal edges at exactly 2R or R apart, where the roundoff of the computed
// distance (~eps x the coordinate magnitude) exceeds any relative margin applied to 2R, so
// a roundoff-level perturbation of the input coordinates would flip the decision. Both
// operands are rounded to the same grid of kDecisionLengthQuantumRelativeToMatchingRadius
// x R before they are compared: the rule is symmetric in the operands and deterministic,
// and a distance within half a quantum of the threshold is AT the threshold (not strictly
// within it). Direction cosines (parallelism tests) are compared on the fixed
// kDecisionDirectionQuantum grid for the same reason. Both constants are recorded in the
// requirements output under Library.DecisionQuantization.
constexpr double kDecisionLengthQuantumRelativeToMatchingRadius = 1.0e-8;
constexpr double kDecisionDirectionQuantum = 1.0e-12;
// Corner class of the legacy per-group classifier (PatchConstruction "Legacy", comparison
// only): a two-segment vertex turning more than this is a corner, a smaller turn continues
// the chain (decision 73). The identification uses the geometric joint noise rule
// kJointNoiseSagittaOverRadius (metaledge.hpp; USER decision 121 (B)) instead.
constexpr double kLegacyCornerTurnToleranceDegrees = 30.0;

// Axisymmetric (r, z) edge sites: the gap direction of an edge on a curved (revolved) edge
// is radial; a site whose in-plane gap direction makes a cosine below this with the r axis
// (more than ~18 deg off radial) is not a revolved edge of the curvature family and is
// refused. The same cosine bounds the "parallel" tests of the 2D site pairing below.
constexpr double kAxisymmetricRadialGapCosine = 0.95;

class DecisionQuantizer
{
private:
  double length_quantum;

public:
  explicit DecisionQuantizer(double matching_radius)
    : length_quantum(kDecisionLengthQuantumRelativeToMatchingRadius * matching_radius)
  {
    MFEM_VERIFY(std::isfinite(length_quantum) && length_quantum > 0.0,
                "Invalid matching radius for quantized classification decisions!");
  }

  double Length(double length) const { return std::round(length / length_quantum); }
  bool LengthLess(double first, double second) const
  {
    return Length(first) < Length(second);
  }
  bool LengthAtMost(double first, double second) const
  {
    return Length(first) <= Length(second);
  }
  bool LengthSquaredLess(double first_squared, double second) const
  {
    return LengthLess(std::sqrt(first_squared), second);
  }

  static double Direction(double cosine)
  {
    return std::round(cosine / kDecisionDirectionQuantum);
  }
  static bool DirectionLess(double first, double second)
  {
    return Direction(first) < Direction(second);
  }
};

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

    for (const auto &[segment, count] : counts)
    {
      MFEM_VERIFY(count <= 2, "Plan-view facets form a nonmanifold surface (segment ["
                                  << segment.first[0] * tolerance << ", "
                                  << segment.first[1] * tolerance << ", "
                                  << segment.first[2] * tolerance << "] - ["
                                  << segment.second[0] * tolerance << ", "
                                  << segment.second[1] * tolerance << ", "
                                  << segment.second[2] * tolerance << "] shared by "
                                  << count << " facets of " << unique_polygons.size()
                                  << " in conductor " << group.first << ")!");
    }
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

std::string TopologyName(LibraryTopology topology);

LibraryTopology CurvedTopologyOf(LibraryTopology straight)
{
  switch (straight)
  {
    case LibraryTopology::ISOLATED_EDGE:
      return LibraryTopology::CURVED_EDGE;
    case LibraryTopology::SAME_CONDUCTOR_GAP:
      return LibraryTopology::CURVED_SAME_CONDUCTOR_GAP;
    case LibraryTopology::DIFFERENT_CONDUCTOR_GAP:
      return LibraryTopology::CURVED_DIFFERENT_CONDUCTOR_GAP;
    case LibraryTopology::SAME_CONDUCTOR_STRIP:
      return LibraryTopology::CURVED_SAME_CONDUCTOR_STRIP;
    default:
      return straight;
  }
}

// Curvature interpolation (decision 92): a curved feature of curvature kappa = R / rho and
// given convexity is modelled by the Lagrange combination of the family's coupons in kappa:
// the straight anchor (kappa 0, the model of the analogous straight topology) and the
// curved coupons (Kappa, Convexity) of the library. Rule, deterministic and recorded: an
// exact node (|kappa - node| <= 1e-9) is that coupon; kappa at or below the first-order
// curvature 1 / StraightBendRadiusOverR (= 0.1: the identification's straight-like
// threshold, decision 75) is the linear combination of the anchor and the family's node AT
// that curvature (the first-order correction dR/dkappa of straight-like features) — a
// library without that node refuses such a kappa (reason), never a linear rule on another
// node; otherwise the cubic Lagrange interpolant on the four nodes nearest to kappa (all
// nodes when fewer). A kappa beyond the largest node is refused (reason), never treated as
// straight.
struct CurvedFamilySelection
{
  std::size_t anchor = 0;
  std::vector<LibrarySelection::WeightedModel> nodes;
  std::string rule;
  double kappa_max = 0.0;
  double first_order_kappa = 1.0 / kStraightBendRadiusOverRadius;
};

std::optional<CurvedFamilySelection>
FindCurvedLibraryModel(const ProcessLibrary &library, LibraryTopology straight_topology,
                       double separation, bool convex, double kappa,
                       const MetalBoundaryLaw &boundary_condition, std::string &reason)
{
  MFEM_VERIFY(std::isfinite(kappa) && kappa > 0.0,
              "Curvature interpolation requires a positive finite kappa = R / rho!");
  const auto anchor =
      FindLibraryModel(library, straight_topology, separation, boundary_condition);
  if (!anchor || anchor->IsInterpolated())
  {
    reason = "no straight " + TopologyName(straight_topology) +
             " anchor coupon for the curvature family";
    return std::nullopt;
  }
  const LibraryTopology curved_topology = CurvedTopologyOf(straight_topology);
  std::vector<std::pair<double, std::size_t>> nodes;
  for (std::size_t i = 0; i < library.models.size(); i++)
  {
    const auto &model = library.models[i];
    if (model.topology != curved_topology || !model.kappa || !model.convex ||
        *model.convex != convex ||
        !CompatibleBoundaryLaw(model.boundary_condition, boundary_condition))
    {
      continue;
    }
    const double separation_tolerance =
        std::max(model.separation_tolerance,
                 1.0e-10 * std::max(library.matching_radius, separation));
    if (std::abs(model.separation - separation) > separation_tolerance)
    {
      continue;
    }
    nodes.emplace_back(*model.kappa, i);
  }
  const std::string convexity = convex ? "convex" : "concave";
  if (nodes.empty())
  {
    reason = "library has no " + convexity + " " + TopologyName(curved_topology) +
             " coupons (Kappa records)";
    return std::nullopt;
  }
  std::sort(nodes.begin(), nodes.end());
  for (std::size_t i = 1; i < nodes.size(); i++)
  {
    MFEM_VERIFY(nodes[i].first - nodes[i - 1].first > 1.0e-9,
                "Curvature family has two " << convexity << " "
                                            << TopologyName(curved_topology)
                                            << " coupons at the same Kappa!");
  }
  CurvedFamilySelection selection;
  selection.anchor = anchor->models.front().index;
  selection.kappa_max = nodes.back().first;
  if (kappa > selection.kappa_max * (1.0 + 1.0e-9))
  {
    reason = "kappa = " + std::to_string(kappa) + " exceeds the largest " + convexity +
             " coupon kappa " + std::to_string(selection.kappa_max);
    return std::nullopt;
  }
  for (const auto &[node_kappa, index] : nodes)
  {
    if (std::abs(kappa - node_kappa) <= 1.0e-9)
    {
      selection.nodes = {{index, 1.0}};
      selection.rule = "exact";
      return selection;
    }
  }
  // Abscissae: the anchor at kappa 0 then the curved nodes.
  std::vector<std::pair<double, std::size_t>> abscissae = {{0.0, selection.anchor}};
  abscissae.insert(abscissae.end(), nodes.begin(), nodes.end());
  std::size_t begin = 0, end = abscissae.size();
  if (kappa <= selection.first_order_kappa * (1.0 + 1.0e-9))
  {
    const auto first_order_node = std::find_if(
        nodes.begin(), nodes.end(), [&](const auto &node)
        { return std::abs(node.first - selection.first_order_kappa) <= 1.0e-9; });
    if (first_order_node == nodes.end())
    {
      reason = "kappa = " + std::to_string(kappa) +
               " is in the first-order regime (<= 1 / StraightBendRadiusOverR = " +
               std::to_string(selection.first_order_kappa) + ") but the library has no " +
               convexity + " coupon at that kappa (smallest node " +
               std::to_string(nodes.front().first) + ")";
      return std::nullopt;
    }
    abscissae = {{0.0, selection.anchor}, *first_order_node};
    end = 2;
    selection.rule = "linear";
  }
  else
  {
    std::size_t interval = 0;
    while (interval + 1 < abscissae.size() && abscissae[interval + 1].first < kappa)
    {
      interval++;
    }
    if (abscissae.size() > 4)
    {
      begin = std::min(interval > 0 ? interval - 1 : 0, abscissae.size() - 4);
      end = begin + 4;
    }
    selection.rule =
        (end - begin == 4) ? "cubic" : (end - begin == 3 ? "quadratic" : "linear");
  }
  for (std::size_t i = begin; i < end; i++)
  {
    double weight = 1.0;
    for (std::size_t j = begin; j < end; j++)
    {
      if (j != i)
      {
        weight *= (kappa - abscissae[j].first) / (abscissae[i].first - abscissae[j].first);
      }
    }
    selection.nodes.push_back({abscissae[i].second, weight});
  }
  return selection;
}

// Curvature family on the three-dimensional FEATURES path (design (b)7 / (e)): a curved
// feature (CurvedEdge, curved pair) carries RadiusOverR and Convexity in its signature;
// when no Signature-keyed model matches it exactly, the family of its straight analogue
// (the anchor at the feature's separation) interpolated at kappa = 1 / RadiusOverR with
// FindCurvedLibraryModel's rule models it. A feature the family cannot model (no anchor,
// no coupons of its convexity, kappa beyond the largest node, "Mixed" convexity: an S-bend
// inside one curved section) is reported unmatched with the reason — never straight.
struct FeatureCurvatureMatch
{
  CurvedFamilySelection selection;
  std::string name;  // runtime model: <anchor>@<convexity>-kappa<kappa>-<rule>
  std::string topology;
  double kappa = 0.0;
  bool convex = true;
};

std::optional<LibraryTopology> StraightTopologyOfCurvedFeature(const std::string &type)
{
  if (type == "CurvedEdge")
  {
    return LibraryTopology::ISOLATED_EDGE;
  }
  if (type == "CurvedSameConductorGap")
  {
    return LibraryTopology::SAME_CONDUCTOR_GAP;
  }
  if (type == "CurvedDifferentConductorGap")
  {
    return LibraryTopology::DIFFERENT_CONDUCTOR_GAP;
  }
  if (type == "CurvedSameConductorStrip")
  {
    return LibraryTopology::SAME_CONDUCTOR_STRIP;
  }
  return std::nullopt;
}

std::vector<std::string> ModelInterfaceNames(const LibraryModel &model)
{
  std::vector<std::string> names;
  for (const auto &interface : model.interfaces)
  {
    names.push_back(ToString(interface.type));
  }
  std::sort(names.begin(), names.end());
  names.erase(std::unique(names.begin(), names.end()), names.end());
  return names;
}

std::string CurvedRuntimeModelName(const std::string &anchor, bool convex, double kappa,
                                   const std::string &rule)
{
  std::ostringstream name;
  name << anchor << "@" << (convex ? "convex" : "concave") << "-kappa"
       << std::setprecision(9) << kappa << "-" << rule;
  return name.str();
}

std::optional<FeatureCurvatureMatch>
MatchCurvatureFamily(const ProcessLibrary &library, const IdentifiedFeature &feature,
                     const MetalBoundaryLaw &boundary_condition, std::string &reason)
{
  const auto straight = StraightTopologyOfCurvedFeature(feature.type);
  if (!straight)
  {
    reason = "no curvature family for " + feature.type;
    return std::nullopt;
  }
  const auto &sig = feature.signature;
  if (!sig.contains("RadiusOverR") || !sig.contains("Convexity"))
  {
    reason = "signature without RadiusOverR / Convexity";
    return std::nullopt;
  }
  const std::string convexity = sig["Convexity"].get<std::string>();
  if (convexity == "Mixed")
  {
    reason = "mixed convexity (bends of both senses in the curved regime inside one "
             "section)";
    return std::nullopt;
  }
  const bool convex = convexity == "Convex";
  const double radius_over_R = sig["RadiusOverR"].get<double>();
  MFEM_VERIFY(radius_over_R > 0.0, "A curved feature requires a positive RadiusOverR!");
  const double kappa = 1.0 / radius_over_R;
  const double separation =
      sig.contains("SeparationOverR")
          ? sig["SeparationOverR"].get<double>() * library.matching_radius
          : 0.0;
  auto selection = FindCurvedLibraryModel(library, *straight, separation, convex, kappa,
                                          boundary_condition, reason);
  if (!selection)
  {
    return std::nullopt;
  }
  const auto &anchor = library.models[selection->anchor];
  if (sig.contains("Interfaces") &&
      sig["Interfaces"].get<std::vector<std::string>>() != ModelInterfaceNames(anchor))
  {
    reason = "the family anchor \"" + anchor.name +
             "\" does not map the feature's interfaces " + sig["Interfaces"].dump();
    return std::nullopt;
  }
  FeatureCurvatureMatch match;
  match.selection = *selection;
  match.name = CurvedRuntimeModelName(anchor.name, convex, kappa, selection->rule);
  match.topology = TopologyName(CurvedTopologyOf(*straight));
  match.kappa = kappa;
  match.convex = convex;
  return match;
}

// Angle-interpolated corner family on the FEATURES path (USER decision 121 (C)): a SHARP
// corner feature (CornerRadiusOverR 0) of one convexity whose angle no library coupon
// matches within the signature tolerance is modelled by the family of sharp corner coupons
// of that convexity, interfaces and law, interpolated in the TURN t = 180 - AngleDegrees:
// the nodes are the coupons' turns (e.g. 90 / 75 / 60 / 45 / 30 / 15 deg for 90 / 105 / 120
// / 135 / 150 / 165 deg corners) and the anchor is the straight edge through the corner box
// (Angle 180, t = 0: the family's 180-deg anchor, built on the same basis). Rule: an exact
// node within kSignatureAngleToleranceDegrees -> that coupon; t at or below the smallest
// node turn -> linear between the anchor and that node (first order in the turn: the corner
// excess of a small turn is O(t), as the curvature family's first-order regime); otherwise
// cubic Lagrange on the four nearest abscissae (anchor included); t above the largest node
// turn (an angle sharper than the sharpest coupon) -> unmatched with the reason — never
// silently straight, never the nearest node. The turn is the interpolation variable because
// the coupon's response is smooth in the arm direction and t = 0 is the straight anchor
// where the first-order regime is anchored (the angle itself is the same variable shifted).
// The runtime model is the BLEND of the nodes' matrices (Lagrange weights may be negative;
// the corner coupons share one box basis, checked at match time) on the basis of the
// nearest node (its ZeroTraceIndices: the knots inside the metal at the nearest coupon
// angle; the corner coupons have no coupon depth, the weights apply as they are). Rounded
// corners keep their per-radius coupons (CornerRadiusInterpolation): no angle family.
// A spatial model's explicit trace triangulation (TraceMesh files, or constructed by the
// corner family's trace basis rule at the device angle). A vertex is a basis knot (basis =
// its 1-based BasisPoints index), a conductor vertex (basis 0, conductor > 0) or a slave
// vertex (basis 0, conductor 0, parents parent_a / parent_b 1-based with weight_a on
// parent_a: a box corner between two knots whose trace is the linear interpolation).
struct TraceMeshData
{
  struct Vertex
  {
    Point3D point{};
    int basis = 0;
    int conductor = 0;
    int parent_a = 0;
    int parent_b = 0;
    double weight_a = 0.0;
  };
  std::vector<Vertex> vertices;
  std::vector<std::array<int, 3>> triangles;

  bool HasSlaveVertices() const
  {
    return std::any_of(vertices.begin(), vertices.end(),
                       [](const Vertex &vertex) { return vertex.parent_a > 0; });
  }
};

TraceMeshData ReadTraceMesh(const std::string &vertex_path,
                            const std::string &triangle_path);

struct FeatureCornerMatch
{
  std::size_t base = 0;  // the nearest node: basis, zero trace indices, references
  std::vector<LibrarySelection::WeightedModel> nodes;
  std::string rule;
  std::string name;  // runtime model: <base>@corner-angle<deg>-<rule>
  std::string topology;
  double angle_degrees = 0.0;
  double turn_degrees = 0.0;
  double max_turn_degrees = 0.0;
  double first_order_turn_degrees = 0.0;
  // The segment connectivity angle of the stencil (the coupons' TraceBasis
  // ConnectivityAngleDegrees; absent for an exact legacy node).
  std::optional<double> connectivity_angle_degrees;
  // The runtime basis of an interpolated corner: the family's trace basis rule at the
  // feature's angle with the segment's connectivity (empty for an exact node, which uses
  // the node's own files).
  std::optional<ConstructedCornerTraceBasis> constructed;
};

std::string CornerRuntimeModelName(const std::string &base, double angle_degrees,
                                   const std::string &rule)
{
  std::ostringstream name;
  name << base << "@corner-angle" << std::setprecision(9) << angle_degrees << "-" << rule;
  return name.str();
}

std::optional<FeatureCornerMatch>
MatchCornerFamily(const ProcessLibrary &library, const IdentifiedFeature &feature,
                  const MetalBoundaryLaw &boundary_condition, std::string &reason)
{
  const auto topology = ParseLibraryTopology(feature.type);
  if (topology != LibraryTopology::CONVEX_CORNER &&
      topology != LibraryTopology::CONCAVE_CORNER)
  {
    reason = "no corner family for " + feature.type;
    return std::nullopt;
  }
  const auto &sig = feature.signature;
  if (!sig.contains("AngleDegrees") || !sig.contains("CornerRadiusOverR"))
  {
    reason = "signature without AngleDegrees / CornerRadiusOverR";
    return std::nullopt;
  }
  if (sig["CornerRadiusOverR"].get<double>() > 0.0)
  {
    reason = "rounded corner (CornerRadiusOverR > 0): per-radius coupons only, no angle "
             "family";
    return std::nullopt;
  }
  const double angle = sig["AngleDegrees"].get<double>();
  MFEM_VERIFY(angle > 0.0 && angle < 180.0,
              "A sharp corner feature requires AngleDegrees strictly between 0 and 180!");
  const double turn = 180.0 - angle;
  const std::vector<std::string> interfaces =
      sig.contains("Interfaces") ? sig["Interfaces"].get<std::vector<std::string>>()
                                 : std::vector<std::string>{};
  // The family: sharp coupons of the feature's topology, interfaces and law, by turn.
  std::vector<std::pair<double, std::size_t>> nodes;  // (turn, index); the anchor at 0
  for (std::size_t i = 0; i < library.models.size(); i++)
  {
    const auto &model = library.models[i];
    if (model.topology != topology || model.corner_radius != 0.0 ||
        !CompatibleBoundaryLaw(model.boundary_condition, boundary_condition) ||
        (sig.contains("Interfaces") && ModelInterfaceNames(model) != interfaces))
    {
      continue;
    }
    nodes.emplace_back(180.0 - model.angle * 180.0 / std::acos(-1.0), i);
  }
  const std::string convexity =
      topology == LibraryTopology::CONVEX_CORNER ? "convex" : "concave";
  if (nodes.empty())
  {
    reason = "library has no sharp " + convexity + " corner coupons for the interfaces " +
             sig["Interfaces"].dump();
    return std::nullopt;
  }
  std::sort(nodes.begin(), nodes.end());
  // One knot SEMANTICS for the whole family (the blend combines matrices entry by entry):
  // every coupon is built on the trace basis rule with the same parameters, the same basis
  // size, ContourGroups and ZeroTraceIndices, its rings that do not meet the metal at the
  // same points, and its rings that meet the metal at the rule's positions for its own
  // angle with its own segment connectivity (checked against the coupon's files). Like-to-
  // like free knots then sit at the same indices with positions that vary smoothly with the
  // angle. Two coupons at one angle are allowed with different segment connectivities (one
  // per side of a knot-corner passage) or a legacy coupon beside them.
  const bool convex = topology == LibraryTopology::CONVEX_CORNER;
  const auto &first = library.models[nodes.front().second];
  for (const auto &[node_turn, index] : nodes)
  {
    (void)node_turn;
    if (!library.models[index].trace_basis)
    {
      // Coupons of the lane-2 angle-independent layout (no TraceBasis record) form no
      // family: their free hats cross the metal at any other angle (the corner-family
      // review's root cause) — unmatched with the reason, never interpolated.
      reason =
          "sharp " + convexity + " corner coupon \"" + library.models[index].name +
          "\" has no TraceBasis rule (lane-2 layout): the angle-interpolated corner "
          "family requires coupons built on the trace basis rule (corner-family review "
          "2026-09-29)";
      return std::nullopt;
    }
  }
  const auto first_points = ReadBasisPoints(first.response.basis_points);
  const CornerBoxRings first_rings = DescribeCornerBoxRings(
      first_points, first.response.contour_groups, first.response.zero_trace_indices);
  const double position_tolerance = 1.0e-9 * first_rings.radius;
  std::vector<CornerFamilyNode> family;
  for (const auto &[node_turn, index] : nodes)
  {
    const auto &model = library.models[index];
    MFEM_VERIFY(model.trace_basis && *model.trace_basis == *first.trace_basis,
                "Corner family coupons \"" << model.name << "\" and \"" << first.name
                                           << "\" are not built on one TraceBasis rule!");
    const auto points = ReadBasisPoints(model.response.basis_points);
    MFEM_VERIFY(points.size() == first_points.size() &&
                    model.response.contour_groups == first.response.contour_groups &&
                    model.response.zero_trace_indices == first.response.zero_trace_indices,
                "Corner family coupons \""
                    << model.name << "\" and \"" << first.name
                    << "\" have different basis sizes, ContourGroups "
                       "or ZeroTraceIndices!");
    // The coupon's files against the rule at its own angle (its fixed rings and its metal
    // rings; the connectivity of its segment).
    std::optional<double> connectivity_radians;
    if (model.corner_connectivity_angle_degrees)
    {
      connectivity_radians =
          *model.corner_connectivity_angle_degrees * std::acos(-1.0) / 180.0;
    }
    const auto rule_basis = BuildCornerTraceBasis(
        points, model.response.contour_groups, model.response.zero_trace_indices,
        model.angle, convex, *model.trace_basis, connectivity_radians);
    for (std::size_t k = 0; k < points.size(); k++)
    {
      MFEM_VERIFY(Distance(points[k], rule_basis.knots[k]) <= position_tolerance,
                  "Corner family coupon \"" << model.name << "\" basis point " << k + 1
                                            << " is not at the trace basis rule's position "
                                               "for its angle!");
    }
    for (const auto &ring : first_rings.rings)
    {
      if (ring.metal || first.trace_basis->AllRings())
      {
        // Under AllRingsFollowMetal every ring follows the angle (checked knot by knot
        // above); only the ring structure is shared.
        continue;
      }
      for (int i = 0; i < ring.size; i++)
      {
        MFEM_VERIFY(Distance(points[ring.offset + i], first_points[ring.offset + i]) <=
                        position_tolerance,
                    "Corner family coupons \"" << model.name << "\" and \"" << first.name
                                               << "\" differ on a ring that does not meet "
                                                  "the metal!");
      }
    }
    // The coupon's trace triangulation against the rule's (the segment's connectivity): a
    // segment node whose bands are triangulated otherwise would blend across a jump of the
    // hats (fail closed; the recorded perimeter-ordered coupons stamped with a connectivity
    // angle pass exactly when no event lies between their angle and it).
    if (model.corner_connectivity_angle_degrees || model.trace_basis->AllRings())
    {
      MFEM_VERIFY(
          !model.response.trace_vertices.empty() && !model.response.trace_triangles.empty(),
          "Corner family coupon \"" << model.name
                                    << "\" has a rule-fixed trace triangulation but "
                                       "no TraceMesh!");
      const auto mesh =
          ReadTraceMesh(model.response.trace_vertices, model.response.trace_triangles);
      MFEM_VERIFY(mesh.vertices.size() == rule_basis.vertices.size() &&
                      mesh.triangles.size() == rule_basis.triangles.size(),
                  "Corner family coupon \""
                      << model.name
                      << "\" trace mesh does not have the rule's vertex / triangle counts "
                         "for its segment connectivity!");
      for (std::size_t v = 0; v < mesh.vertices.size(); v++)
      {
        const auto &vertex = mesh.vertices[v];
        const auto &rule_vertex = rule_basis.vertices[v];
        const std::array<double, 3> point = {vertex.point[0], vertex.point[1],
                                             vertex.point[2]};
        MFEM_VERIFY(Distance(point, rule_vertex.point) <= position_tolerance &&
                        vertex.basis - 1 == rule_vertex.basis &&
                        vertex.parent_a - 1 == rule_vertex.parent_a &&
                        vertex.parent_b - 1 == rule_vertex.parent_b,
                    "Corner family coupon \"" << model.name << "\" trace vertex " << v + 1
                                              << " differs from the rule's for its segment "
                                                 "connectivity!");
      }
      std::set<std::array<int, 3>> coupon_triangles, rule_triangles;
      for (auto triangle : mesh.triangles)
      {
        std::sort(triangle.begin(), triangle.end());
        coupon_triangles.insert(triangle);
      }
      for (auto triangle : rule_basis.triangles)
      {
        std::sort(triangle.begin(), triangle.end());
        rule_triangles.insert(triangle);
      }
      MFEM_VERIFY(
          coupon_triangles == rule_triangles,
          "Corner family coupon \""
              << model.name << "\" trace triangulation is not the rule's"
              << (model.corner_connectivity_angle_degrees
                      ? " for its segment connectivity angle " +
                            std::to_string(*model.corner_connectivity_angle_degrees) +
                            " deg (the band triangulation differs: the coupon "
                            "belongs to another segment)"
                      : " at its angle")
              << "!");
    }
    CornerFamilyNode node;
    node.angle_degrees = 180.0 - node_turn;
    node.connectivity_angle_degrees = model.corner_connectivity_angle_degrees;
    node.index = index;
    family.push_back(node);
  }
  FeatureCornerMatch match;
  match.angle_degrees = angle;
  match.turn_degrees = turn;
  match.topology = feature.type;
  match.max_turn_degrees = nodes.back().first;
  const bool has_anchor = std::abs(nodes.front().first) <= kSignatureAngleToleranceDegrees;
  match.first_order_turn_degrees =
      nodes.size() > (has_anchor ? 1u : 0u) ? nodes[has_anchor ? 1 : 0].first : 0.0;
  // The stencil: exact node, or Lagrange on the nodes of the segment (coupons sharing a
  // connectivity angle) containing the angle, never across a knot-corner passage of the
  // trace basis (SelectCornerFamilyStencil; legacy coupons are exact matches only). The
  // segment structure itself was verified at library load (ReadProcessLibrary).
  const auto stencil = SelectCornerFamilyStencil(family, angle, convex, *first.trace_basis,
                                                 kSignatureAngleToleranceDegrees);
  if (!stencil.reason.empty())
  {
    reason = stencil.reason;
    return std::nullopt;
  }
  match.base = stencil.base;
  match.rule = stencil.rule;
  match.connectivity_angle_degrees = stencil.connectivity_angle_degrees;
  for (const auto &[index, weight] : stencil.nodes)
  {
    match.nodes.push_back({index, weight});
  }
  match.name = CornerRuntimeModelName(library.models[match.base].name, angle, match.rule);
  if (match.rule == "exact")
  {
    return match;
  }
  // The runtime basis at the feature's angle: the base's fixed rings, the metal rings by
  // the rule (positions) with the segment's connectivity, the same zero set and contour
  // groups as every node.
  {
    const auto &base = library.models[match.base];
    std::optional<double> connectivity_radians;
    if (stencil.connectivity_angle_degrees)
    {
      connectivity_radians = *stencil.connectivity_angle_degrees * std::acos(-1.0) / 180.0;
    }
    match.constructed = BuildCornerTraceBasis(
        ReadBasisPoints(base.response.basis_points), base.response.contour_groups,
        base.response.zero_trace_indices, angle * std::acos(-1.0) / 180.0, convex,
        *base.trace_basis, connectivity_radians);
  }
  return match;
}

// First-order curvature term of a straight-like feature (windowed bend radius >= 10 R,
// described by its straight model): the family node AT kappa = 1 / StraightBendRadiusOverR
// of each convexity for the feature's straight model (its anchor). The term is evaluated
// per portion on its signed windowed turn (IdentifiedPortion::turn); a missing node is
// recorded (the feature keeps its straight model: the recorded straight-like rule).
struct FirstOrderNodes
{
  std::optional<std::size_t> convex, concave;
};

FirstOrderNodes FindFirstOrderNodes(const ProcessLibrary &library, std::size_t anchor_index)
{
  FirstOrderNodes nodes;
  const auto &anchor = library.models[anchor_index];
  const LibraryTopology curved_topology = CurvedTopologyOf(anchor.topology);
  if (curved_topology == anchor.topology)
  {
    return nodes;
  }
  const double first_order_kappa = 1.0 / kStraightBendRadiusOverRadius;
  for (std::size_t i = 0; i < library.models.size(); i++)
  {
    const auto &model = library.models[i];
    if (model.topology != curved_topology || !model.kappa || !model.convex ||
        std::abs(*model.kappa - first_order_kappa) > 1.0e-9 ||
        !CompatibleBoundaryLaw(model.boundary_condition, anchor.boundary_condition) ||
        ModelInterfaceNames(model) != ModelInterfaceNames(anchor))
    {
      continue;
    }
    const double separation_tolerance =
        std::max(model.separation_tolerance,
                 1.0e-10 * std::max(library.matching_radius, anchor.separation));
    if (std::abs(model.separation - anchor.separation) > separation_tolerance)
    {
      continue;
    }
    auto &slot = *model.convex ? nodes.convex : nodes.concave;
    MFEM_VERIFY(!slot, "Curvature family has two " << (*model.convex ? "convex" : "concave")
                                                   << " " << TopologyName(curved_topology)
                                                   << " coupons at the first-order kappa!");
    slot = i;
  }
  return nodes;
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
    case LibraryTopology::CURVED_EDGE:
      return "curved edge";
    case LibraryTopology::CURVED_SAME_CONDUCTOR_GAP:
      return "curved same-conductor gap";
    case LibraryTopology::CURVED_DIFFERENT_CONDUCTOR_GAP:
      return "curved different-conductor gap";
    case LibraryTopology::CURVED_SAME_CONDUCTOR_STRIP:
      return "curved same-conductor strip";
  }
  return "unknown";
}

// Longitudinal families integrated along the edge (per unit coupon depth); the curved
// classes are integrated along their curved portions in the same way.
bool IsTranslationalTopology(std::string_view topology)
{
  return topology == "isolated edge" || topology == "same-conductor gap" ||
         topology == "different-conductor gap" || topology == "same-conductor strip" ||
         topology == "parallel-edge cluster" || topology == "curved edge" ||
         topology == "curved same-conductor gap" ||
         topology == "curved different-conductor gap" ||
         topology == "curved same-conductor strip";
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
    case LibraryTopology::CURVED_EDGE:
      return "CurvedEdge";
    case LibraryTopology::CURVED_SAME_CONDUCTOR_GAP:
      return "CurvedSameConductorGap";
    case LibraryTopology::CURVED_DIFFERENT_CONDUCTOR_GAP:
      return "CurvedDifferentConductorGap";
    case LibraryTopology::CURVED_SAME_CONDUCTOR_STRIP:
      return "CurvedSameConductorStrip";
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
  // Version-2 contract: the identification's feature records replace the legacy per-pass
  // records, which are kept for comparison only under LegacyRequirements.
  std::map<std::string, Aggregate> legacy_requirements;
  nlohmann::json identification;
  bool identification_active = false;
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
  double CoordinateScale() const { return coordinate_scale; }

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

  // Switch the legacy Add() calls to the comparison table; the Requirements array is then
  // derived from the identification features through AddFeatureRecord().
  void ActivateIdentification() { identification_active = true; }
  bool IdentificationActive() const { return identification_active; }
  void SetIdentification(nlohmann::json value) { identification = std::move(value); }

  void AddFeatureRecord(nlohmann::json requirement, int count, double length)
  {
    const std::string key = requirement.dump();
    auto [it, inserted] = requirements.emplace(key, Aggregate{requirement, 0, 0.0});
    (void)inserted;
    it->second.count += count;
    it->second.total_edge_length += ScaleLength(length);
  }

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
    auto &table = identification_active ? legacy_requirements : requirements;
    auto [it, inserted] = table.emplace(key, Aggregate{requirement, 0, 0.0});
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
    std::map<std::string, int> counts = {{"Exact", 0}, {"Interpolated", 0}, {"Missing", 0}};
    std::map<std::string, double> lengths = {
        {"Exact", 0.0}, {"Interpolated", 0.0}, {"Missing", 0.0}};
    auto Entries = [](const std::map<std::string, Aggregate> &table,
                      std::map<std::string, int> *count_table,
                      std::map<std::string, double> *length_table)
    {
      nlohmann::json entries = nlohmann::json::array();
      for (const auto &[key, aggregate] : table)
      {
        (void)key;
        auto entry = aggregate.requirement;
        entry["Count"] = aggregate.count;
        entry["TotalEdgeLength"] = aggregate.total_edge_length;
        if (count_table)
        {
          (*count_table)[entry["Status"].get<std::string>()] += aggregate.count;
          (*length_table)[entry["Status"].get<std::string>()] +=
              aggregate.total_edge_length;
        }
        entries.push_back(std::move(entry));
      }
      return entries;
    };
    nlohmann::json entries = Entries(requirements, &counts, &lengths);
    nlohmann::json result = {
        {"Version", identification_active ? 2 : 1},
        {"Complete", counts["Missing"] == 0},
        {"Library",
         {{"Path", library_path},
          {"Name", library_name},
          {"MatchingRadius", matching_radius},
          {"DecisionQuantization",
           {{"LengthRelativeToMatchingRadius",
             kDecisionLengthQuantumRelativeToMatchingRadius},
            {"Direction", kDecisionDirectionQuantum}}}}},
        {"LengthUnit", "mesh"},
        {"Summary", {{"Counts", counts}, {"TotalEdgeLengths", lengths}}},
        {"Requirements", std::move(entries)}};
    if (identification_active)
    {
      result["Identification"] = identification;
      result["LegacyRequirements"] = Entries(legacy_requirements, nullptr, nullptr);
    }
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
  // Axisymmetric (r, z) device: every edge site at r = rho is the circular edge of a disk
  // (gap toward +r: metal inside, convex) or a hole (concave) with kappa = R / rho and edge
  // length 2 pi rho; a pair of sites is a concentric annular gap / strip with kappa =
  // R / rho_inner; both are corrected by the curvature family interpolated at kappa.
  const bool axisymmetric = mat_op.GetMesh().IsAxisymmetric();
  const double coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
  const auto library =
      ReadProcessLibrary(request.library, iodata.units, iodata.InputsNondimensionalized(),
                         requirements != nullptr, requirements != nullptr);
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
  result.translational_domain_correction = request.translational_domain_correction;
  result.trace_coupling = request.trace_coupling;
  result.mortar_oversampling = request.mortar_oversampling;
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
      if (!requirements &&
          request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
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
    const DecisionQuantizer quantizer(group.matching_radius);
    const double interaction_distance = 2.0 * group.matching_radius;
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
              quantizer.LengthLess(Distance(sites[current].point, sites[neighbor].point),
                                   interaction_distance))
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
      std::optional<CurvedFamilySelection> curved_selection;
      double curved_kappa = 0.0, curved_edge_length = 0.0;
      bool curved_convex = true;
      const auto boundary_condition = sites[cluster.front()].boundary_condition;
      auto ClusterGeometry = [&]()
      {
        nlohmann::json geometry = {{"EdgeCount", cluster.size()}};
        if (axisymmetric)
        {
          geometry["AxisymmetricRadius"] =
              requirements ? requirements->ScaleLength(sites[cluster.front()].point[0])
                           : 0.0;
          if (curved_kappa > 0.0)
          {
            geometry["Kappa"] = curved_kappa;
            geometry["Convexity"] = curved_convex ? "Convex" : "Concave";
          }
          if (curved_selection)
          {
            geometry["InterpolationRule"] = curved_selection->rule;
            geometry["KappaMax"] = curved_selection->kappa_max;
            geometry["FirstOrderKappa"] = curved_selection->first_order_kappa;
          }
        }
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
        if (axisymmetric)
        {
          MFEM_VERIFY(edge.point[0] > 0.0 &&
                          std::abs(edge.axis_u[0]) > kAxisymmetricRadialGapCosine,
                      "An axisymmetric edge site must lie off the axis with an in-plane "
                      "radial gap direction (|cos| > "
                          << kAxisymmetricRadialGapCosine << ")!");
          curved_kappa = group.matching_radius / edge.point[0];
          curved_convex = edge.axis_u[0] > 0.0;
          curved_edge_length = 2.0 * M_PI * edge.point[0];
        }
      }
      else if (axisymmetric && cluster.size() > 2)
      {
        group_matched = false;
        Mpi::Warning("Axisymmetric response matching supports isolated edges and pairs "
                     "only ({} nearby edges at r = {:.6e} mesh units); correction is "
                     "disabled for this interface group!\n",
                     cluster.size(), sites[cluster.front()].point[0] * coordinate_scale);
        if (requirements)
        {
          requirements->Add(2, LibraryTopology::SPATIAL_EDGE_CLUSTER, group.targets,
                            boundary_condition, ClusterGeometry(), library, nullptr, 0.0,
                            "Curved cluster coupons are not available");
          unmatched_clusters++;
          continue;
        }
        unmatched_clusters++;
        break;
      }
      else if (cluster.size() == 2)
      {
        const auto &first = sites[cluster[0]];
        const auto &second = sites[cluster[1]];
        separation = Distance(first.point, second.point);
        if (axisymmetric)
        {
          // Concentric pair (an annular gap / strip): the family's Kappa = R / rho_inner,
          // Convexity that of the model's first edge e1 = `first` (metal inside the bend
          // when its gap points to +r), CouponDepth = the centreline 2 pi (rho + s / 2)
          // (the pair measure: the mean of the two sides).
          MFEM_VERIFY(std::min(first.point[0], second.point[0]) > 0.0 &&
                          std::abs(first.axis_u[0]) > kAxisymmetricRadialGapCosine &&
                          std::abs(second.axis_u[0]) > kAxisymmetricRadialGapCosine,
                      "An axisymmetric pair must lie off the axis with in-plane radial gap "
                      "directions (|cos| > "
                          << kAxisymmetricRadialGapCosine << ")!");
          curved_kappa = group.matching_radius / std::min(first.point[0], second.point[0]);
          curved_convex = first.axis_u[0] > 0.0;
          curved_edge_length = M_PI * (first.point[0] + second.point[0]);
        }
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

      if (!model_selection && curved_kappa > 0.0)
      {
        // Never silently straight: an axisymmetric edge is corrected by the curvature
        // family interpolated at its kappa, or reported unmatched with the reason.
        std::string reason;
        curved_selection =
            FindCurvedLibraryModel(library, topology, separation, curved_convex,
                                   curved_kappa, boundary_condition, reason);
        if (!curved_selection)
        {
          group_matched = false;
          Mpi::Warning("Fabrication-process response library \"{}\" cannot model the {} "
                       "curved edge at r = {:.6e} mesh units (kappa = {:.4f}): {}; "
                       "correction is disabled for this interface group!\n",
                       library.name, curved_convex ? "convex" : "concave",
                       sites[cluster.front()].point[0] * coordinate_scale, curved_kappa,
                       reason);
          if (requirements)
          {
            requirements->Add(2, CurvedTopologyOf(topology), group.targets,
                              boundary_condition, ClusterGeometry(), library, nullptr, 0.0,
                              "Curvature family: " + reason);
            unmatched_clusters++;
            continue;
          }
          unmatched_clusters++;
          break;
        }
        LibrarySelection selection;
        selection.models = curved_selection->nodes;
        selection.conductor_references =
            library.models[curved_selection->anchor].conductor_references;
        selection.normalized_distance = 0.0;
        model_selection = selection;
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
        requirements->Add(2,
                          cluster.size() > 2 ? LibraryTopology::PARALLEL_EDGE_CLUSTER
                          : curved_selection ? CurvedTopologyOf(topology)
                                             : topology,
                          group.targets, boundary_condition, ClusterGeometry(), library,
                          &*model_selection, 0.0);
      }
      // The patch carries the selected model's conductor references for every cluster size
      // (a ResponsePatchData starts with the configuration default {{0, 0, 0}}, so an
      // "if empty" guard never fired and the 1- / 2-site patches kept the default: a
      // same-conductor gap pair read its reference potential in the gap, a
      // different-conductor gap pair aborted on the reference count).
      patch.conductor_references = model_selection->conductor_references;
      if (diagnostics && !patch.maxwell_reference_is_pec &&
          patch.conductor_references.size() == 2 && cluster.size() == 2)
      {
        const auto &second = sites[cluster[1]];
        patch.maxwell_conductor_anchors.push_back(
            std::array<double, 3>{second.point[0], second.point[1], 0.0});
      }
      group_maximum_library_distance =
          std::max(group_maximum_library_distance, model_selection->normalized_distance);
      if (curved_selection)
      {
        // One patch on a runtime model interpolated in kappa (Lagrange weights may be
        // negative, so the combination is formed on the matrices, not on patch weights);
        // the patch weight is the curved edge length over the anchor's coupon depth.
        const auto &anchor = library.models[curved_selection->anchor];
        MFEM_VERIFY(anchor.coupon_depth > 0.0,
                    "Curvature interpolation requires CouponDepth on the straight anchor!");
        std::ostringstream name;
        name << anchor.name << "@" << (curved_convex ? "convex" : "concave") << "-kappa"
             << std::setprecision(9) << curved_kappa << "-" << curved_selection->rule;
        PendingBlend blend{name.str(), TopologyName(CurvedTopologyOf(topology)),
                           curved_selection->anchor, curved_selection->nodes};
        auto weighted_patch = patch;
        weighted_patch.weight = curved_edge_length / anchor.coupon_depth;
        pending.push_back({curved_selection->anchor, std::move(weighted_patch), blend});
        continue;
      }
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
      if (!requirements &&
          request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
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
    std::map<std::string, int> runtime_models;
    for (auto &selection : pending)
    {
      const auto &source = library.models[selection.library_model];
      const std::string key = selection.blend ? selection.blend->name : source.name;
      auto [model_it, inserted] = runtime_models.emplace(key, next_model_index);
      if (inserted)
      {
        auto model = source.response;
        model.idx = next_model_index++;
        model.name = key;
        model.topology = TopologyName(source.topology);
        if (selection.blend)
        {
          // The anchor's basis and interfaces; matrices = sum of the nodes' matrices, each
          // per unit edge length and rescaled to the anchor's coupon depth.
          model.topology = selection.blend->topology;
          std::ostringstream nodes;
          for (const auto &node : selection.blend->nodes)
          {
            const auto &coupon = library.models[node.index];
            MFEM_VERIFY(coupon.coupon_depth > 0.0,
                        "Curvature interpolation requires CouponDepth on every coupon!");
            model.blend.push_back({node.weight * source.coupon_depth / coupon.coupon_depth,
                                   coupon.response.fabricated_matrix,
                                   coupon.response.thin_matrix,
                                   coupon.response.fabricated_surface_matrix,
                                   coupon.response.thin_surface_matrix});
            nodes << (model.blend.size() == 1 ? "" : ", ") << coupon.name << " x "
                  << std::setprecision(6) << node.weight;
          }
          Mpi::Print(" Curvature interpolation {}: {}\n", key, nodes.str());
        }
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

ElementBox SpatialSupportBox(const SpatialClusterSelection3D &selection,
                             const LibraryModel &model)
{
  ElementBox box;
  for (const auto &local : model.support_points)
  {
    const auto point = TransformLocalPoint(selection.origin, selection.axes, local);
    ElementBox point_box;
    point_box.min = point;
    point_box.max = point;
    box.Add(point_box);
  }
  return box;
}

std::set<int> SpatialTargetAttributes(const SpatialClusterSelection3D &selection)
{
  std::set<int> targets;
  for (const auto &[slot, interfaces] : selection.targets_by_slot)
  {
    (void)slot;
    for (const auto &[type, target] : interfaces)
    {
      (void)type;
      targets.insert(target);
    }
  }
  return targets;
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

  auto OwnershipPriority = [](const LibraryModel &model)
  {
    ElementBox support;
    for (const auto &point : model.support_points)
    {
      for (int d = 0; d < 3; d++)
      {
        support.min[d] = std::min(support.min[d], point[d]);
        support.max[d] = std::max(support.max[d], point[d]);
      }
    }
    double volume = 0.0;
    if (!model.support_points.empty())
    {
      volume = 1.0;
      for (int d = 0; d < 3; d++)
      {
        volume *= std::max(0.0, support.max[d] - support.min[d]);
      }
    }
    double interval_span = 0.0;
    for (const auto &edge : model.spatial_edges)
    {
      interval_span += edge.interval[1] - edge.interval[0];
    }
    return std::make_pair(volume, interval_span);
  };

  std::optional<SpatialClusterSelection3D> best;
  double best_distance = mfem::infinity();
  std::pair<double, double> best_priority{-mfem::infinity(), -mfem::infinity()};
  for (std::size_t model_index = 0; model_index < library.models.size(); model_index++)
  {
    const auto &model = library.models[model_index];
    if (excluded_models.find(model_index) != excluded_models.end() ||
        model.topology != LibraryTopology::SPATIAL_EDGE_CLUSTER ||
        model.spatial_edges.size() < sites.size())
    {
      continue;
    }
    const double position_tolerance =
        std::max(model.spatial_position_tolerance, 1.0e-10 * library.matching_radius);
    const double angle_tolerance =
        std::max(model.spatial_angle_tolerance, 1.0e-10 * std::acos(-1.0));
    const auto &model_anchor = model.spatial_edges.front();
    const auto ownership_priority = OwnershipPriority(model);

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
          // atan2 of the cross and dot products: acos(1 - eps) of a snapped, nearly-unit
          // model direction would report sqrt(2 eps) (1e-6 rad for 12-digit components).
          auto AngleBetween = [](const Point3D &a, const Point3D &b)
          { return std::atan2(Norm(Cross(a, b)), Dot(a, b)); };
          const double gap_error = AngleBetween(gap, site.gap_direction);
          const double normal_error = AngleBetween(normal, site.process_normal);
          // A unified owner can deliberately extend farther along a physical edge than
          // the interaction site which selects it. Require complete containment rather
          // than identical clipped intervals so that one enlarged coupon can own every
          // overlapping subordinate support.
          const double interval_error = std::max({0.0, edge.interval[0] - site.interval[0],
                                                  site.interval[1] - edge.interval[1]});
          const double normalized_error = std::max(
              {position_error / position_tolerance, gap_error / angle_tolerance,
               normal_error / angle_tolerance, interval_error / position_tolerance});
          if (normalized_error <= 1.0)
          {
            candidates[edge_index].push_back({site_index, normalized_error});
          }
        }
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
      std::function<void(std::size_t, std::size_t, double)> Match =
          [&](std::size_t depth, std::size_t assigned_sites, double normalized_distance)
      {
        if (normalized_distance > best_distance)
        {
          return;
        }
        if (depth == order.size())
        {
          if (assigned_sites != sites.size() ||
              !MatchLibraryInterfaces(model, targets_by_slot))
          {
            return;
          }
          const double distance_tolerance =
              1.0e-12 * std::max({1.0, normalized_distance, best_distance});
          if (best && std::abs(normalized_distance - best_distance) <= distance_tolerance &&
              ownership_priority <= best_priority)
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
          best_priority = ownership_priority;
          return;
        }

        const std::size_t edge_index = order[depth];
        const std::size_t remaining_edges = order.size() - depth - 1;
        const std::size_t remaining_sites = sites.size() - assigned_sites;
        if (remaining_edges >= remaining_sites)
        {
          Match(depth + 1, assigned_sites, normalized_distance);
        }
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
          Match(depth + 1, assigned_sites + 1,
                std::max(normalized_distance, candidate.distance));
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
      Match(0, 0, 0.0);
    }
  }
  return best;
}

// Library-model signatures computed from the stored model geometry with the
// identification's canonicalisation, so that the matching pass is a lookup of feature
// signatures.
// A model whose law parameters are not verified cannot equal a device law; give it a key
// that no feature produces instead of exporting the unverified parameters.
std::string LibraryLawKey(const MetalBoundaryLaw &law,
                          const AutomaticResponseRequirements &describer)
{
  return law.parameters_verified ? describer.DescribeBoundaryCondition(law).dump()
                                 : nlohmann::json{{"Type", BoundaryConditionName(law.type)},
                                                  {"Unverified", true}}
                                       .dump();
}

// Canonical cluster signature of a SpatialEdgeCluster model's stored edges; its frame
// (origin, axes in the model's own coordinates) maps the model onto a device feature with
// the same signature (the feature's frame is canonical in mesh coordinates).
std::optional<CanonicalSignature>
ModelClusterSignature(const LibraryModel &model, double radius,
                      const AutomaticResponseRequirements &describer)
{
  std::vector<SignaturePortion> portions;
  // Every interface type of a slot, sorted by name: the identification serialises a
  // portion's interfaces as a sorted set, and the canonical frame is chosen by that
  // serialisation, so the model must serialise the same way to land in the same frame.
  std::map<int, std::set<std::string>> slot_types;
  for (const auto &interface : model.interfaces)
  {
    slot_types[interface.slot].insert(ToString(interface.type));
  }
  for (const auto &edge : model.spatial_edges)
  {
    // The stored Interval runs along gap x normal: the convention of the spatial coupon
    // generator, of AlignSpatialFrame and of the legacy site geometry below (normal x gap
    // read the edges mirrored through their Point and canonicalised the model to a frame
    // other than the one it was built in).
    const Point3D tangent = Normalize(Cross(edge.gap_direction, edge.process_normal));
    std::vector<std::string> edge_interfaces;
    if (auto it = slot_types.find(edge.interface_slot); it != slot_types.end())
    {
      edge_interfaces.assign(it->second.begin(), it->second.end());
    }
    portions.push_back({Add(edge.point, Scale(edge.interval[0], tangent)),
                        Add(edge.point, Scale(edge.interval[1], tangent)),
                        edge.gap_direction, edge.conductor, edge_interfaces,
                        LibraryLawKey(edge.boundary_condition, describer)});
  }
  if (portions.empty())
  {
    return std::nullopt;
  }
  return CanonicalClusterSignature(portions, {}, model.spatial_edges.front().process_normal,
                                   radius);
}

// A version-2 SpatialEdgeCluster model (one carrying the identification's Signature) stores
// its Edges in the canonical frame of that Signature: the library builder's contract (Point
// = P x R, Interval along gap x normal, process normal +z; an arc portion chorded on its
// circle). The patch construction then maps the model onto a feature with the identity
// (the feature's canonical frame is the model frame), which is exact for straight and arc
// portions alike — ModelClusterSignature cannot rebuild an arc portion from chords, so a
// re-canonicalisation of the chorded Edges would land in another frame. The contract is
// verified here (fail closed): every edge lies on one portion of the Signature (both
// endpoints within the signature parameter tolerance) whose relabelled Conductor is the
// edge's Conductor and whose interface set is the set of the edge's InterfaceSlot — a
// geometrically mirror-symmetric cluster with asymmetric labels lands on the portions but
// would attach its coupons to the wrong conductors / interfaces.
void VerifySpatialEdgesInSignatureFrame(const LibraryModel &model, double radius)
{
  MFEM_VERIFY(model.identification_signature && !model.spatial_edges.empty(),
              "VerifySpatialEdgesInSignatureFrame needs a Signature and Edges!");
  const auto &signature = *model.identification_signature;
  MFEM_VERIFY(signature.contains("Portions") && signature["Portions"].is_array() &&
                  !signature["Portions"].empty(),
              "SpatialEdgeCluster model \"" << model.name
                                            << "\" Signature carries no Portions!");
  // In units of R: the arc-fit tolerance widened by two signature quanta and read inclusive
  // (block (b) DESIGN A1 (5), decision 288 (2) half-quantum rule): the builder's chord rows
  // are the REBUILT circle's chords through the serialised / snapped ends, which deviate
  // from the model signature's arc by at most the fit tolerance plus the rounding of the
  // ends; nothing is compared at a picometre bound that noise can cross.
  const double tolerance = kSignatureParameterToleranceOverRadius +
                           2.0 * kSignatureLengthQuantumOverRadius +
                           0.5 * kSignatureLengthQuantumOverRadius;
  // The interface set of every InterfaceSlot (the types mapped to it; empty without a
  // mapping), serialised as the Signature serialises a portion's Interfaces (sorted).
  std::map<int, std::set<std::string>> slot_types;
  for (const auto &interface : model.interfaces)
  {
    slot_types[interface.slot].insert(ToString(interface.type));
  }
  auto SlotInterfaceSet = [&](int slot)
  {
    std::set<std::string> types;
    if (auto it = slot_types.find(slot); it != slot_types.end())
    {
      types = it->second;
    }
    return types;
  };
  auto PortionInterfaceSet = [](const nlohmann::json &portion)
  {
    std::set<std::string> types;
    for (const auto &type : portion.value("Interfaces", nlohmann::json::array()))
    {
      types.insert(type.get<std::string>());
    }
    return types;
  };
  auto JoinSet = [](const std::set<std::string> &types)
  {
    std::string joined = "[";
    for (const auto &type : types)
    {
      joined += (joined.size() > 1 ? ", " : "") + type;
    }
    return joined + "]";
  };
  // Distance (in units of R) from a point q to a portion: a straight portion is the
  // segment P; an arc portion the arc from P[0] through the midpoint to P[1] (a closed
  // circle when the ends coincide).
  auto PortionDistance = [&](const nlohmann::json &portion, const std::array<double, 2> &q)
  {
    const auto P = portion.at("P").get<std::array<double, 4>>();
    const std::array<double, 2> a = {P[0], P[1]}, b = {P[2], P[3]};
    if (!portion.contains("Arc"))
    {
      const double dx = b[0] - a[0], dy = b[1] - a[1];
      const double length2 = dx * dx + dy * dy;
      double s = length2 > 0.0 ? ((q[0] - a[0]) * dx + (q[1] - a[1]) * dy) / length2 : 0.0;
      s = std::clamp(s, 0.0, 1.0);
      return std::hypot(q[0] - (a[0] + s * dx), q[1] - (a[1] + s * dy));
    }
    const auto arc = portion.at("Arc").get<std::array<double, 4>>();
    const std::array<double, 2> c = {arc[0], arc[1]}, m = {arc[2], arc[3]};
    const double r = std::hypot(a[0] - c[0], a[1] - c[1]);
    const double radial = std::abs(std::hypot(q[0] - c[0], q[1] - c[1]) - r);
    const bool closed = std::hypot(a[0] - b[0], a[1] - b[1]) <= 1.0e-9 * std::max(r, 1.0);
    if (closed)
    {
      return radial;
    }
    const double two_pi = 2.0 * std::acos(-1.0);
    auto Angle = [&](const std::array<double, 2> &p)
    { return std::atan2(p[1] - c[1], p[0] - c[0]); };
    const double ta = Angle(a);
    // The arc runs from a to b through m: counterclockwise when m lies on the
    // counterclockwise sweep from a to b, clockwise otherwise.
    const double ccw = std::fmod(Angle(b) - ta + two_pi, two_pi);
    const bool counterclockwise =
        std::fmod(Angle(m) - ta + two_pi, two_pi) <= ccw + 1.0e-12;
    const double sweep = counterclockwise ? ccw : ccw - two_pi;
    double tq = std::fmod(Angle(q) - ta + two_pi, two_pi);
    if (!counterclockwise)
    {
      tq = tq > 0.0 ? tq - two_pi : tq;
    }
    const bool within = counterclockwise ? tq <= sweep + 1.0e-12 : tq >= sweep - 1.0e-12;
    if (within)
    {
      return radial;
    }
    return std::min(std::hypot(q[0] - a[0], q[1] - a[1]),
                    std::hypot(q[0] - b[0], q[1] - b[1]));
  };
  for (std::size_t e = 0; e < model.spatial_edges.size(); e++)
  {
    const auto &edge = model.spatial_edges[e];
    MFEM_VERIFY(std::abs(edge.process_normal[2] - 1.0) <= 1.0e-9 &&
                    std::abs(edge.process_normal[0]) <= 1.0e-9 &&
                    std::abs(edge.process_normal[1]) <= 1.0e-9,
                "SpatialEdgeCluster model \""
                    << model.name
                    << "\" carries a Signature: its Edges must lie in the "
                       "Signature's canonical frame (ProcessNormal [0, 0, 1]); edge "
                    << e << " has ProcessNormal [" << edge.process_normal[0] << ", "
                    << edge.process_normal[1] << ", " << edge.process_normal[2] << "]!");
    const Point3D tangent = Normalize(Cross(edge.gap_direction, edge.process_normal));
    std::array<Point3D, 2> endpoints;
    for (std::size_t k = 0; k < 2; k++)
    {
      endpoints[k] = Add(edge.point, Scale(edge.interval[k], tangent));
      MFEM_VERIFY(std::abs(endpoints[k][2]) / radius <= tolerance,
                  "SpatialEdgeCluster model \""
                      << model.name << "\" carries a Signature but edge " << e
                      << " endpoint [" << endpoints[k][0] << ", " << endpoints[k][1] << ", "
                      << endpoints[k][2] << "] lies off the process plane (tolerance "
                      << tolerance
                      << " R): a version-2 model's Edges must be the Signature's portions "
                         "in its canonical frame!");
    }
    // The portion holding the edge: both endpoints within the tolerance (a builder edge is
    // one straight portion or one chord of one arc portion). The nearest portion by the
    // farther endpoint reports a geometric miss.
    double nearest = std::numeric_limits<double>::infinity();
    std::size_t nearest_portion = 0;
    bool labels_agree = false;
    for (std::size_t p = 0; p < signature["Portions"].size(); p++)
    {
      const auto &portion = signature["Portions"][p];
      double distance = 0.0;
      for (const auto &endpoint : endpoints)
      {
        const std::array<double, 2> q = {endpoint[0] / radius, endpoint[1] / radius};
        distance = std::max(distance, PortionDistance(portion, q));
      }
      if (distance < nearest)
      {
        nearest = distance;
        nearest_portion = p;
      }
      labels_agree =
          labels_agree ||
          (distance <= tolerance && portion.value("Conductor", 0) == edge.conductor &&
           PortionInterfaceSet(portion) == SlotInterfaceSet(edge.interface_slot));
    }
    MFEM_VERIFY(nearest <= tolerance,
                "SpatialEdgeCluster model \""
                    << model.name << "\" carries a Signature but edge " << e
                    << " endpoints [" << endpoints[0][0] << ", " << endpoints[0][1]
                    << "] / [" << endpoints[1][0] << ", " << endpoints[1][1] << "] lie "
                    << nearest << " R from the nearest Signature portion (tolerance "
                    << tolerance
                    << " R): a version-2 model's Edges must be the Signature's portions "
                       "in its canonical frame (Point = P x R, Interval along gap x "
                       "normal; arc portions chorded on their circle)!");
    const auto &portion = signature["Portions"][nearest_portion];
    MFEM_VERIFY(labels_agree,
                "SpatialEdgeCluster model \""
                    << model.name << "\" carries a Signature but edge " << e
                    << " (Conductor " << edge.conductor << ", InterfaceSlot "
                    << edge.interface_slot << " = "
                    << JoinSet(SlotInterfaceSet(edge.interface_slot))
                    << ") lies on Signature portion " << nearest_portion << " (Conductor "
                    << portion.value("Conductor", 0) << ", Interfaces "
                    << JoinSet(PortionInterfaceSet(portion))
                    << "): a version-2 model's Edges carry their portion's relabelled "
                       "Conductor and the InterfaceSlot mapped to its interface set!");
  }
}

// Library models by signature topology (decision 85(1)): a feature matches the model of its
// type whose topology key equals its own and whose continuous parameters lie within the
// signature tolerance (kSignatureParameterToleranceOverRadius /
// kSignatureAngleToleranceDegrees, both orientations of a translational signature); among
// several the nearest (the smallest normalised deviation), then the smallest rank (corner
// family coupons at one angle, corner-qualification block 2026-09-29: the legacy coupon
// with the tie triangulation at its own angle 0, the coupon of the lower-angle segment 1,
// of the upper-angle segment 2), ties by model name — a deterministic choice independent of
// the feature or model order.
class LibrarySignatureIndex
{
public:
  void Add(nlohmann::json signature, const std::string &type, const std::string &name,
           int rank = 0)
  {
    signature["Type"] = type;
    const std::size_t index = entries.size();
    for (const nlohmann::json &orientation :
         {signature, MirrorTranslationalSignature(signature)})
    {
      by_topology[SplitSignatureParameters(orientation).topology_key].push_back(index);
    }
    entries.push_back({std::move(signature), name, rank});
  }

  // The matched model, its normalised deviation (<= 1) and its signature, or nullopt.
  struct Matched
  {
    std::string name;
    double deviation = 0.0;
    const nlohmann::json *signature = nullptr;
  };
  std::optional<Matched> Match(const nlohmann::json &feature_signature) const
  {
    const auto it =
        by_topology.find(SplitSignatureParameters(feature_signature).topology_key);
    if (it == by_topology.end())
    {
      return std::nullopt;
    }
    std::optional<std::tuple<double, int, std::string, std::size_t>> best;
    std::set<std::size_t> seen;
    for (const std::size_t index : it->second)
    {
      if (!seen.insert(index).second)
      {
        continue;
      }
      const auto deviation =
          SignatureDeviation(feature_signature, entries[index].signature);
      if (!deviation || *deviation > 1.0)
      {
        continue;
      }
      const std::tuple<double, int, std::string, std::size_t> candidate{
          *deviation, entries[index].rank, entries[index].name, index};
      if (!best || candidate < *best)
      {
        best = candidate;
      }
    }
    if (!best)
    {
      return std::nullopt;
    }
    return Matched{std::get<2>(*best), std::get<0>(*best),
                   &entries[std::get<3>(*best)].signature};
  }

  // Two SpatialEdgeCluster models within 2 x kClusterQuantumNearMatchMaxQuanta quanta of
  // each other (half-quantum inclusive) are one geometry at the grid (block (b) DESIGN
  // section 4, MINOR-5): refused at load, so that no feature can lie within the near-match
  // of two models.
  void RefuseNearDuplicateClusters() const
  {
    for (const auto &[topology, indices] : by_topology)
    {
      (void)topology;
      for (std::size_t i = 0; i < indices.size(); i++)
      {
        const auto &a = entries[indices[i]];
        if (a.signature.value("Type", std::string{}) != "SpatialEdgeCluster")
        {
          continue;
        }
        for (std::size_t j = i + 1; j < indices.size(); j++)
        {
          const auto &b = entries[indices[j]];
          if (indices[i] == indices[j])
          {
            continue;
          }
          const auto difference =
              ClusterSignatureQuantumDifference(a.signature, b.signature);
          MFEM_VERIFY(!difference || !ClusterQuantumDuplicate(difference->max_delta_quanta),
                      "Library models \""
                          << a.name << "\" and \"" << b.name
                          << "\" are SpatialEdgeCluster keys of one topology within "
                          << (difference ? difference->max_delta_quanta : 0.0)
                          << " signature quanta of each other (<= "
                          << 2 * kClusterQuantumNearMatchMaxQuanta << " + "
                          << kClusterQuantumInclusiveMargin
                          << "): two models, one geometry at the 1e-6 R grid (block (b) "
                             "DESIGN section 4: keep the lexicographically smallest key's "
                             "model, list the other key under its NearKeys)!");
        }
      }
    }
  }

  std::size_t Size() const { return entries.size(); }

  // Legacy-contract aliases (USER decision 283), keyed by the aliased v3 hash: the legacy
  // model's name and Signature with the alias. Looked up by the exact hash only.
  struct AliasEntry
  {
    std::string model;
    nlohmann::json model_signature;
    LegacyContractAlias alias;
  };
  void AddAlias(const std::string &model, const nlohmann::json &model_signature,
                const LegacyContractAlias &alias)
  {
    aliases.emplace(alias.key, AliasEntry{model, model_signature, alias});
  }
  const AliasEntry *FindAlias(const std::string &hash) const
  {
    const auto it = aliases.find(hash);
    return it == aliases.end() ? nullptr : &it->second;
  }

private:
  struct Entry
  {
    nlohmann::json signature;
    std::string name;
    int rank = 0;
  };
  std::vector<Entry> entries;
  std::map<std::string, std::vector<std::size_t>> by_topology;
  std::map<std::string, AliasEntry> aliases;
};

}  // namespace

nlohmann::json ResolveLegacyContractAlias(const std::string &model_name,
                                          const nlohmann::json &model_signature,
                                          const LegacyContractAlias &alias,
                                          const IdentifiedFeature &feature)
{
  MFEM_VERIFY(feature.hash == alias.key,
              "Legacy-contract alias of model \""
                  << model_name << "\" resolved for a feature whose key " << feature.hash
                  << " is not the alias key " << alias.key << "!");
  const std::string digest = SpatialSupportContextDigest(feature.signature);
  MFEM_VERIFY(!digest.empty() && digest == alias.context_digest,
              "Legacy-contract alias "
                  << alias.key << " of model \"" << model_name
                  << "\": the feature's context digest "
                  << (digest.empty() ? std::string("(no Box)") : digest)
                  << " does not match the alias's ContextDigest " << alias.context_digest
                  << " (USER decision 283: fail closed)!");
  // The feature's claims-only key (the claims canonicalised alone, recorded by the
  // identification: the Portions of a contract-3 signature are serialised in the frame
  // that minimises Box + Context, so they cannot be re-keyed by removing those members).
  const std::string claims_hash =
      feature.spatial_support.is_object()
          ? feature.spatial_support.value("ClaimsKey", std::string{})
          : std::string{};
  const auto [model_key, model_hash] = SignatureKeyAndHash(model_signature, feature.type);
  (void)model_key;
  MFEM_VERIFY(!claims_hash.empty() && claims_hash == model_hash,
              "Legacy-contract alias "
                  << alias.key << " of model \"" << model_name
                  << "\": the feature's claims-only key " << claims_hash
                  << " is not the legacy model's key " << model_hash
                  << " (the alias names another geometry; fail closed)!");
  return nlohmann::json{
      {"Model", model_name},
      {"Key", alias.key},
      {"ContextDigest", alias.context_digest},
      {"ClaimsKey", claims_hash},
      {"Reason", alias.reason},
      {"Context",
       {{"Box", feature.signature.value("Box", nlohmann::json(nullptr))},
        {"Context", feature.signature.value("Context", nlohmann::json::array())}}},
      {"Rule", "USER decision 283: an explicit library-side alias from a contract-3 "
               "key to a legacy model (decision-236 straight-continuation coupon); "
               "usable ONLY for the listed key, never a fallback; the feature's "
               "context digest and claims-only key are verified against the alias "
               "and the model (mismatch aborts)"}};
}

namespace
{

LibrarySignatureIndex
LibrarySignatureKeys(const ProcessLibrary &library,
                     const AutomaticResponseRequirements &requirements)
{
  LibrarySignatureIndex index;
  const double R = library.matching_radius;
  auto Law = [&](const MetalBoundaryLaw &law) { return LibraryLawKey(law, requirements); };
  for (const auto &model : library.models)
  {
    std::vector<std::string> interfaces;
    for (const auto &interface : model.interfaces)
    {
      interfaces.push_back(ToString(interface.type));
    }
    std::sort(interfaces.begin(), interfaces.end());
    interfaces.erase(std::unique(interfaces.begin(), interfaces.end()), interfaces.end());
    const std::string law = Law(model.boundary_condition);
    // Corner family coupons at one angle (one per side of a knot-corner passage beside a
    // legacy tie coupon): the exact match prefers the legacy coupon, then the lower-angle
    // segment's (SelectCornerFamilyStencil's rule; corner-qualification block 2026-09-29).
    int rank = 0;
    if (model.corner_connectivity_angle_degrees)
    {
      rank =
          *model.corner_connectivity_angle_degrees < model.angle * 180.0 / std::acos(-1.0)
              ? 1
              : 2;
    }
    std::optional<std::pair<nlohmann::json, std::string>> key;  // signature, type
    if (model.identification_signature)
    {
      index.Add(*model.identification_signature,
                model.identification_signature->at("Type").get<std::string>(), model.name,
                rank);
      for (const auto &alias : model.legacy_contract_aliases)
      {
        index.AddAlias(model.name, *model.identification_signature, alias);
      }
      continue;
    }
    switch (model.topology)
    {
      case LibraryTopology::ISOLATED_EDGE:
        key = std::make_pair(nlohmann::json{{"Interfaces", interfaces}, {"Law", law}},
                             std::string("IsolatedEdge"));
        break;
      case LibraryTopology::SAME_CONDUCTOR_GAP:
      case LibraryTopology::DIFFERENT_CONDUCTOR_GAP:
      case LibraryTopology::SAME_CONDUCTOR_STRIP:
        {
          const bool strip = model.topology == LibraryTopology::SAME_CONDUCTOR_STRIP;
          const bool different = model.topology == LibraryTopology::DIFFERENT_CONDUCTOR_GAP;
          std::vector<TranslationalEdge> edges = {
              {0.0, strip ? -1 : 1, 1, interfaces, law},
              {model.separation, strip ? 1 : -1, different ? 2 : 1, interfaces, law}};
          key = std::make_pair(CanonicalTranslationalSignature(edges, R).signature,
                               TopologyIdentifier(model.topology));
          break;
        }
      case LibraryTopology::PARALLEL_EDGE_CLUSTER:
        {
          std::vector<TranslationalEdge> edges;
          for (const auto &edge : model.cluster_edges)
          {
            edges.push_back({edge.offset, edge.gap_direction >= 0 ? 1 : -1, edge.conductor,
                             interfaces, law});
          }
          if (edges.size() >= 2)
          {
            key = std::make_pair(CanonicalTranslationalSignature(edges, R).signature,
                                 std::string("ParallelEdgeCluster"));
          }
          break;
        }
      case LibraryTopology::CONVEX_CORNER:
      case LibraryTopology::CONCAVE_CORNER:
        key = std::make_pair(CanonicalCornerSignature(interfaces, law,
                                                      model.angle * 180.0 / std::acos(-1.0),
                                                      model.corner_radius / R),
                             TopologyIdentifier(model.topology));
        break;
      case LibraryTopology::ENDPOINT:
        key = std::make_pair(nlohmann::json{{"Interfaces", interfaces}, {"Law", law}},
                             std::string("Endpoint"));
        break;
      case LibraryTopology::JUNCTION:
        {
          // The model stores sorted absolute arm angles (radians); the identification
          // signature uses the consecutive differences in degrees.
          std::vector<double> sorted = model.arm_angles, differences;
          std::sort(sorted.begin(), sorted.end());
          for (std::size_t i = 0; i < sorted.size(); i++)
          {
            const double next =
                i + 1 < sorted.size() ? sorted[i + 1] : sorted[0] + 2.0 * std::acos(-1.0);
            differences.push_back((next - sorted[i]) * 180.0 / std::acos(-1.0));
          }
          if (!differences.empty())
          {
            key = std::make_pair(CanonicalJunctionSignature(interfaces, law, differences),
                                 std::string("Junction"));
          }
          break;
        }
      case LibraryTopology::SPATIAL_EDGE_CLUSTER:
        {
          if (const auto canonical = ModelClusterSignature(model, R, requirements))
          {
            nlohmann::json signature = canonical->signature;
            signature["EdgeCount"] = model.spatial_edges.size();
            key = std::make_pair(signature, std::string("SpatialEdgeCluster"));
          }
          break;
        }
      case LibraryTopology::CURVED_EDGE:
      case LibraryTopology::CURVED_SAME_CONDUCTOR_GAP:
      case LibraryTopology::CURVED_DIFFERENT_CONDUCTOR_GAP:
      case LibraryTopology::CURVED_SAME_CONDUCTOR_STRIP:
        break;  // keyed by the Signature only (verified at load)
    }
    if (key)
    {
      index.Add(key->first, key->second, model.name, rank);
    }
  }
  return index;
}

// Geometry identification (design: SURFACE-RESPONSE-IDENTIFICATION.md): pure function of
// the perimeter, the segment frames and R; then the key-based matching pass and the derived
// version-1 requirement records.
// The identification runs in every path (preflight and solve): its per-segment exclusions
// (non-planar, non-manifold, undetermined process side, cross-layer zones) also remove
// those segments from the legacy classification. The manifest records are written only
// when a requirements sink is given (preflight).
IdentificationResult RunGeometryIdentification(
    MPI_Comm comm, const MetalEdgeGeometry &geometry,
    const std::vector<EdgeSegment3D> &framed_segments, const ProcessLibrary &library,
    const AutomaticResponseRequirements &describer,
    AutomaticResponseRequirements *requirements, bool frame_normal_configured,
    const std::vector<ResponseCorrectionData::SpanCapAllowanceData> &span_cap_allowances,
    std::map<int, FeatureCurvatureMatch> *curved_matches = nullptr,
    std::map<int, FeatureCornerMatch> *corner_matches = nullptr)
{
  const bool root = Mpi::Root(comm);
  // A one-sided edge whose face is not parallel to the process plane (the area-weighted
  // principal direction of the metal normals) is non-planar metal: a wall, a staple leg, a
  // via; the classification's kParallelCosineTolerance.
  auto NonPlanarFace = [&](const MetalEdgeSegment &source)
  {
    for (const auto &normal : source.face_normals)
    {
      double dot = 0.0;
      for (int d = 0; d < 3; d++)
      {
        dot += normal[d] * geometry.layer_normal[d];
      }
      if (std::abs(dot) < 1.0 - 1.0e-8)
      {
        return true;
      }
    }
    return false;
  };
  IdentificationInput input;
  input.radius = library.matching_radius;
  for (const auto &allowance : span_cap_allowances)
  {
    // Per-case span-cap allowances (block (b) DESIGN A4): validated by the identification.
    input.span_cap_allowances.push_back({nlohmann::json::parse(allowance.claims_signature),
                                         allowance.span_cap_over_R, allowance.label,
                                         allowance.reason, allowance.approval});
  }
  std::map<std::size_t, const EdgeSegment3D *> framed;
  for (const auto &segment : framed_segments)
  {
    framed.emplace(segment.geometry_index, &segment);
  }
  input.segments.resize(geometry.segments.size());
  for (std::size_t i = 0; i < geometry.segments.size(); i++)
  {
    const auto &source = geometry.segments[i];
    auto &segment = input.segments[i];
    segment.p0 = geometry.vertices[source.vertices[0]].coordinate;
    segment.p1 = geometry.vertices[source.vertices[1]].coordinate;
    segment.vertices = source.vertices;
    segment.chain = source.physical_chain;
    segment.truncation = source.type == MetalEdgeSegmentType::TRUNCATION;
    segment.conductor = source.metal_component;
    if (source.type == MetalEdgeSegmentType::PORT)
    {
      segment.exclusion = std::make_pair(
          "Port", "metal perimeter bordering a port boundary face (LumpedPort / WavePort "
                  "attribute of the configuration): the port is not metal, the metal edge "
                  "along it is a cut (decision 82(5))");
    }
    else if (source.on_bounding_box && !segment.truncation)
    {
      segment.exclusion = std::make_pair(
          "SimulationBoundary",
          "edge of a metal face on the bounding box of the mesh (PEC simulation box)");
    }
    else if (source.type == MetalEdgeSegmentType::FOLD)
    {
      segment.exclusion = std::make_pair(
          "NonPlanar", "edge where two non-coplanar metal faces meet (the metal folds: "
                       "staple, box edge)");
    }
    else if (source.type == MetalEdgeSegmentType::NONMANIFOLD)
    {
      segment.exclusion = std::make_pair(
          "NonManifold", "edge shared by three or more metal face directions (a wall "
                         "standing on a sheet)");
    }
    else if (!segment.truncation && NonPlanarFace(source))
    {
      segment.exclusion = std::make_pair(
          "NonPlanar", "edge of a metal face not parallel to the process plane (wall, "
                       "staple, via)");
    }
    else if (auto it = framed.find(i); it != framed.end())
    {
      segment.conductor = it->second->conductor;
      segment.targets = it->second->targets;
      segment.gap_direction = it->second->axis_u;
      segment.process_normal = it->second->axis_v;
      segment.boundary_law =
          describer.DescribeBoundaryCondition(it->second->boundary_condition).dump();
      if (it->second->ambiguous_process_side)
      {
        segment.exclusion = std::make_pair(
            "UndeterminedProcessSide",
            "metal sheet with the same material on both sides and no EdgeFrameNormal on "
            "the target interface: the process side cannot be determined");
      }
    }
    else if (!segment.truncation && source.side_attributes.size() <= 1 &&
             !frame_normal_configured)
    {
      segment.exclusion = std::make_pair(
          "UndeterminedProcessSide",
          "metal sheet with the same material on both sides and no EdgeFrameNormal on "
          "the target interface: the process side cannot be determined");
    }
    else if (!segment.truncation)
    {
      segment.exclusion = std::make_pair(
          "Untargeted", "physical metal edge without a target interface (or excluded by "
                        "EdgeExcludeAttributes)");
    }
  }
  input.vertices.resize(geometry.vertices.size());
  for (std::size_t v = 0; v < geometry.vertices.size(); v++)
  {
    input.vertices[v].coordinate = geometry.vertices[v].coordinate;
    input.vertices[v].segments = geometry.vertices[v].segments;
    input.vertices[v].physical_type = geometry.vertices[v].physical_type;
    input.vertices[v].on_truncation_boundary = geometry.vertices[v].on_truncation_boundary;
    input.vertices[v].on_port_boundary = geometry.vertices[v].on_port_boundary;
  }
  // The identification runs on the root only (the global metal faces exist there alone;
  // decision 82 infrastructure) and the result is broadcast: every rank builds the same
  // patches from the same feature list, no rank replicates the identification state.
  IdentificationResult result;
  {
    std::string buffer;
    if (root)
    {
      input.faces.reserve(geometry.global_faces.size());
      for (const auto &face : geometry.global_faces)
      {
        if (!face.on_bounding_box)  // a PEC simulation box is not fabricated metal
        {
          input.faces.push_back({face.vertices, face.normal});
        }
      }
      // Stage counts, wall times and progress of the identification.
      input.log = [](const std::string &line) { Mpi::Print("{}", line); };
      result = IdentifyMetalPerimeter(input);
      input.faces.clear();
      input.faces.shrink_to_fit();
      buffer = SerializeIdentificationResult(result);
    }
    std::int64_t size = static_cast<std::int64_t>(buffer.size());
    Mpi::Broadcast(1, &size, 0, comm);
    if (!root)
    {
      buffer.resize(static_cast<std::size_t>(size));
    }
    Mpi::BroadcastLarge(size, buffer.data(), 0, comm);
    if (!root)
    {
      result = DeserializeIdentificationResult(buffer);
    }
  }

  {
    std::string exclusions;
    for (const auto &exclusion : result.exclusions)
    {
      exclusions += fmt::format("{}{}: {:d} ({:.6e})", exclusions.empty() ? "" : ", ",
                                exclusion.cls, exclusion.count, exclusion.length);
    }
    Mpi::Print(
        "Geometry identification: {:d} features on {:d} segments; exclusions {{{}}}\n",
        static_cast<int>(result.features.size()), static_cast<int>(result.segments.size()),
        exclusions);
  }
  int point_contacts = 0;
  for (const auto &vertex : result.vertices)
  {
    if (vertex.point_contact)
    {
      point_contacts++;
      const auto &p = geometry.vertices[vertex.vertex].coordinate;
      Mpi::Warning("Metal of different edge-connected components meets at the point ({:.6e}"
                   ", {:.6e}, {:.6e}) (a point contact, degenerate geometry; vertex type "
                   "{})!\n",
                   p[0], p[1], p[2], vertex.type);
    }
  }
  // Matching pass (the solve path consumes it as well): the model of the feature's
  // topology within the signature parameter tolerance, the nearest one (decision 85(1)).
  const auto library_keys = LibrarySignatureKeys(library, describer);
  std::map<int, nlohmann::json> curvature_records;  // per feature, for the manifest
  std::map<int, nlohmann::json> corner_records;     // per feature, for the manifest
  library_keys.RefuseNearDuplicateClusters();
  for (auto &feature : result.features)
  {
    if (const auto match = library_keys.Match(feature.signature))
    {
      feature.matched_model = match->name;
      feature.match_deviation = match->deviation;
      if (feature.type == "SpatialEdgeCluster" && match->deviation > 0.0)
      {
        // The quantum near-match (block (b) DESIGN section 4, decision 303): the model's
        // key differs from the feature's by <= kClusterQuantumNearMatchMaxQuanta quanta;
        // the feature is placed in its own canonical frame exactly as an exact match.
        const auto difference =
            ClusterSignatureQuantumDifference(feature.signature, *match->signature);
        MFEM_VERIFY(difference.has_value(),
                    "A near-matched cluster without a quantum difference record!");
        const auto [model_key, model_hash] =
            SignatureKeyAndHash(*match->signature, feature.type);
        (void)model_key;
        feature.quantum_near_match =
            QuantumNearMatchRecord(model_hash, feature.hash, *difference);
        feature.match_note =
            "quantum near-match (block (b) DESIGN section 4): key " +
            feature.hash.substr(0, 12) + " resolved to model \"" + match->name +
            "\" (key " + model_hash.substr(0, 12) + ", " +
            std::to_string(difference->differing_paths.size()) +
            " numbers differ by <= " + fmt::format("{:g}", difference->max_delta_quanta) +
            " quanta)";
      }
    }
    else if (const auto *alias = library_keys.FindAlias(feature.hash))
    {
      // Legacy-contract alias (USER decision 283): the library lists this contract-3 key
      // explicitly for a legacy model; verified (digest, claims) or aborted, recorded.
      feature.legacy_contract = ResolveLegacyContractAlias(
          alias->model, alias->model_signature, alias->alias, feature);
      // The legacy model's Edges and basis points live in the claims-only canonical frame:
      // the feature is placed in that frame (the contract-3 frame minimises Box + Context
      // and may differ by a rotation / reflection), so that the patches are the legacy
      // library's patches exactly.
      feature.origin = feature.claims_origin;
      feature.axes = feature.claims_axes;
      feature.chirality = feature.claims_chirality;
      feature.legacy_contract["PlacementFrame"] =
          "the claims-only canonical frame (Features[].ClaimsFrame), the legacy model's";
      feature.matched_model = alias->model;
      feature.match_deviation = 0.0;
      feature.match_note =
          "legacy contract (USER decision 283): key " + feature.hash.substr(0, 12) +
          " resolved through the alias of \"" + alias->model + "\" (context digest " +
          alias->alias.context_digest.substr(0, 12) + " verified)";
    }
    else if (StraightTopologyOfCurvedFeature(feature.type) && !feature.portions.empty())
    {
      // Curvature family (never silently straight): matched at its kappa or reported with
      // the reason. The same deterministic selection on every rank (library + signature).
      const auto it = framed.find(feature.portions.front().segment);
      MFEM_VERIFY(it != framed.end(),
                  "A curved feature portion lies on a segment without an edge frame!");
      std::string reason;
      const auto match =
          MatchCurvatureFamily(library, feature, it->second->boundary_condition, reason);
      if (match)
      {
        feature.matched_model = match->name;
        feature.match_deviation = 0.0;
        std::ostringstream note;
        note << "curvature family: " << match->selection.rule << " at kappa "
             << std::setprecision(6) << match->kappa << " ("
             << (match->convex ? "convex" : "concave") << ")";
        feature.match_note = note.str();
        nlohmann::json nodes = nlohmann::json::array();
        for (const auto &node : match->selection.nodes)
        {
          nodes.push_back({{"Name", library.models[node.index].name},
                           {"Kappa", library.models[node.index].kappa.value_or(0.0)},
                           {"Weight", node.weight}});
        }
        curvature_records[feature.id] = {
            {"Anchor", library.models[match->selection.anchor].name},
            {"Kappa", match->kappa},
            {"Convexity", match->convex ? "Convex" : "Concave"},
            {"InterpolationRule", match->selection.rule},
            {"KappaMax", match->selection.kappa_max},
            {"FirstOrderKappa", match->selection.first_order_kappa},
            {"Nodes", nodes}};
        if (curved_matches)
        {
          curved_matches->emplace(feature.id, *match);
        }
      }
      else
      {
        feature.match_note = "curvature family: " + reason;
      }
    }
    else if ((feature.type == "ConvexCorner" || feature.type == "ConcaveCorner") &&
             !feature.portions.empty())
    {
      // Angle-interpolated corner family (USER decision 121 (C)): a sharp corner without
      // an exact coupon is interpolated in its turn or reported with the reason (never
      // silently straight, never the nearest coupon).
      const auto it = framed.find(feature.portions.front().segment);
      MFEM_VERIFY(it != framed.end(),
                  "A corner feature portion lies on a segment without an edge frame!");
      std::string reason;
      const auto match =
          MatchCornerFamily(library, feature, it->second->boundary_condition, reason);
      if (match)
      {
        feature.matched_model = match->name;
        feature.match_deviation = 0.0;
        std::ostringstream note;
        note << "corner family: " << match->rule << " at " << std::setprecision(6)
             << match->angle_degrees << " deg (turn " << match->turn_degrees << ")";
        feature.match_note = note.str();
        nlohmann::json nodes = nlohmann::json::array();
        for (const auto &node : match->nodes)
        {
          nodes.push_back(
              {{"Name", library.models[node.index].name},
               {"AngleDegrees", library.models[node.index].angle * 180.0 / std::acos(-1.0)},
               {"Weight", node.weight}});
        }
        corner_records[feature.id] = {
            {"Base", library.models[match->base].name},
            {"AngleDegrees", match->angle_degrees},
            {"TurnDegrees", match->turn_degrees},
            {"Convexity", feature.type == "ConvexCorner" ? "Convex" : "Concave"},
            {"InterpolationRule", match->rule},
            {"MaxTurnDegrees", match->max_turn_degrees},
            {"FirstOrderTurnDegrees", match->first_order_turn_degrees},
            {"ConnectivityAngleDegrees",
             match->connectivity_angle_degrees
                 ? nlohmann::json(*match->connectivity_angle_degrees)
                 : nlohmann::json()},
            {"Nodes", nodes}};
        if (corner_matches)
        {
          corner_matches->emplace(feature.id, *match);
        }
      }
      else
      {
        feature.match_note = "corner family: " + reason;
      }
    }
  }
  if (!requirements || !root)
  {
    return result;  // the manifest is built and written on the root
  }

  // Version-1 records derived from the features. Features of one topology whose signature
  // parameters agree within the signature tolerance are ONE record (the library's coupon,
  // decision 85(1)): single-linkage grouping over the distinct signatures (order
  // independent), the record's Geometry / Hash from the group's representative signature
  // (RepresentativeSignature), Count and TotalEdgeLength summed, the group's parameter
  // spread recorded.
  requirements->ActivateIdentification();
  const double R = library.matching_radius;
  auto GeometryOf = [&](const std::string &type, const nlohmann::json &sig)
  {
    nlohmann::json geometry_json;
    if (type == "IsolatedEdge")
    {
      // The version-1 isolated-edge record carried {"EdgeCount": 1}; consumers of the
      // derived Requirements (prepare_surface_response_coupons.plan_from_manifest) read
      // Geometry.
      geometry_json["EdgeCount"] = 1;
    }
    else if (type == "CurvedEdge")
    {
      geometry_json["EdgeCount"] = 1;
      geometry_json["BendRadius"] =
          requirements->ScaleLength(sig["RadiusOverR"].get<double>() * R);
      geometry_json["Kappa"] = 1.0 / sig["RadiusOverR"].get<double>();
      if (sig.contains("Convexity"))
      {
        geometry_json["Convexity"] = sig["Convexity"];
      }
    }
    else if (type == "SameConductorGap" || type == "DifferentConductorGap" ||
             type == "SameConductorStrip" || type == "UnclassifiedParallelPair" ||
             type == "CurvedSameConductorGap" || type == "CurvedDifferentConductorGap" ||
             type == "CurvedSameConductorStrip" || type == "CurvedUnclassifiedParallelPair")
    {
      geometry_json["EdgeCount"] = 2;
      geometry_json["Separation"] =
          requirements->ScaleLength(sig["SeparationOverR"].get<double>() * R);
      if (sig.contains("RadiusOverR"))
      {
        geometry_json["BendRadius"] =
            requirements->ScaleLength(sig["RadiusOverR"].get<double>() * R);
        geometry_json["Kappa"] = 1.0 / sig["RadiusOverR"].get<double>();
      }
      if (sig.contains("Convexity"))
      {
        geometry_json["Convexity"] = sig["Convexity"];
      }
    }
    else if (type == "ParallelEdgeCluster" || type == "CurvedParallelEdgeCluster")
    {
      nlohmann::json edges = nlohmann::json::array();
      for (const auto &edge : sig["Edges"])
      {
        edges.push_back(
            {{"Offset",
              {requirements->ScaleLength(edge["OffsetOverR"].get<double>() * R), 0.0}},
             {"GapDirection", {static_cast<double>(edge["GapSide"].get<int>()), 0.0}},
             {"Conductor", edge["Conductor"]}});
      }
      geometry_json["Edges"] = edges;
      geometry_json["EdgeCount"] = sig["Edges"].size();
      if (sig.contains("RadiusOverR"))
      {
        geometry_json["BendRadius"] =
            requirements->ScaleLength(sig["RadiusOverR"].get<double>() * R);
      }
    }
    else if (type == "ConvexCorner" || type == "ConcaveCorner")
    {
      geometry_json["AngleDegrees"] = sig["AngleDegrees"];
      geometry_json["CornerRadius"] =
          requirements->ScaleLength(sig["CornerRadiusOverR"].get<double>() * R);
    }
    else if (type == "Junction")
    {
      geometry_json["ArmAnglesDegrees"] = sig["ArmAnglesDegrees"];
    }
    else if (type == "SpatialEdgeCluster")
    {
      geometry_json["EdgeCount"] = sig["EdgeCount"];
      geometry_json["Signature"] = sig;
    }
    return geometry_json;
  };
  struct Instance
  {
    nlohmann::json signature;
    int count = 0;     // version-1 Count: mesh segments (longitudinal) or features (vertex)
    int features = 0;  // feature instances (the library's coupon Instances)
    double length = 0.0;
    std::set<std::string> models;
    bool exact = true;
    nlohmann::json curvature_family;    // the family selection of a curved feature
    nlohmann::json corner_family;       // the family selection of an interpolated corner
    std::set<std::string> notes;        // matching notes (family refusals)
    nlohmann::json legacy_contract;     // the alias record (USER decision 283) + Features
    nlohmann::json span_cap_allowance;  // SpatialSupport.SpanCapAllowance (block (b) A4)
  };
  struct GroupBase
  {
    std::string type;
    nlohmann::json interfaces, law;
    std::map<std::string, Instance> instances;  // by signature key
    std::vector<std::string> order;
  };
  std::map<std::string, GroupBase> bases;
  for (const auto &feature : result.features)
  {
    std::set<std::map<InterfaceDielectric, int>> target_maps;
    for (const auto &portion : feature.portions)
    {
      target_maps.insert(input.segments[portion.segment].targets);
    }
    nlohmann::json interfaces = nlohmann::json::array();
    int slot = 0;
    for (const auto &targets : target_maps)
    {
      for (const auto &[type, target] : targets)
      {
        interfaces.push_back(
            {{"Slot", slot}, {"Type", ToString(type)}, {"Target", target}});
      }
      slot++;
    }
    nlohmann::json law =
        feature.signature.contains("Law")
            ? nlohmann::json::parse(feature.signature["Law"].get<std::string>())
            : nlohmann::json{{"Type", "PEC"}};
    if (!feature.portions.empty() && !feature.signature.contains("Law"))
    {
      law = nlohmann::json::parse(
          input.segments[feature.portions.front().segment].boundary_law);
    }
    const bool per_segment = feature.type != "ConvexCorner" &&
                             feature.type != "ConcaveCorner" &&
                             feature.type != "Junction" && feature.type != "Endpoint" &&
                             feature.type != "SpatialEdgeCluster";
    const int count = per_segment ? static_cast<int>(feature.portions.size()) : 1;
    const std::string base_key = feature.type + "|" + interfaces.dump() + "|" + law.dump() +
                                 "|" +
                                 SplitSignatureParameters(feature.signature).topology_key;
    // Both orientations of a translational signature belong to one base (the mirror's
    // topology key may be the smaller one).
    const std::string mirror_key =
        feature.type + "|" + interfaces.dump() + "|" + law.dump() + "|" +
        SplitSignatureParameters(MirrorTranslationalSignature(feature.signature))
            .topology_key;
    auto [it, inserted] = bases.emplace(std::min(base_key, mirror_key), GroupBase{});
    if (inserted)
    {
      it->second.type = feature.type;
      it->second.interfaces = interfaces;
      it->second.law = law;
    }
    auto [instance, new_instance] =
        it->second.instances.emplace(feature.signature_key, Instance{});
    if (new_instance)
    {
      instance->second.signature = feature.signature;
      it->second.order.push_back(feature.signature_key);
    }
    instance->second.count += count;
    instance->second.features++;
    instance->second.length += feature.length;
    instance->second.exact = instance->second.exact && feature.exact_parameters;
    if (feature.matched_model)
    {
      instance->second.models.insert(*feature.matched_model);
    }
    if (!feature.legacy_contract.is_null())
    {
      if (instance->second.legacy_contract.is_null())
      {
        instance->second.legacy_contract = feature.legacy_contract;
        instance->second.legacy_contract["Features"] = nlohmann::json::array();
      }
      instance->second.legacy_contract["Features"].push_back(feature.id);
    }
    if (feature.spatial_support.is_object() &&
        feature.spatial_support.contains("SpanCapAllowance") &&
        instance->second.span_cap_allowance.is_null())
    {
      instance->second.span_cap_allowance = feature.spatial_support["SpanCapAllowance"];
    }
    if (const auto record = curvature_records.find(feature.id);
        record != curvature_records.end())
    {
      instance->second.curvature_family = record->second;
    }
    else if (const auto corner = corner_records.find(feature.id);
             corner != corner_records.end())
    {
      instance->second.corner_family = corner->second;
    }
    else if (feature.match_note && !feature.matched_model)
    {
      instance->second.notes.insert(*feature.match_note);
    }
  }
  for (auto &[base_key, base] : bases)
  {
    (void)base_key;
    // Single linkage over the distinct signatures at the tolerance (union-find; the
    // result does not depend on the order).
    const std::size_t n = base.order.size();
    std::vector<std::size_t> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    auto Find = [&](std::size_t i)
    {
      while (parent[i] != i)
      {
        parent[i] = parent[parent[i]];
        i = parent[i];
      }
      return i;
    };
    for (std::size_t i = 0; i < n; i++)
    {
      for (std::size_t j = i + 1; j < n; j++)
      {
        const auto deviation =
            SignatureDeviation(base.instances.at(base.order[i]).signature,
                               base.instances.at(base.order[j]).signature);
        if (deviation && *deviation <= 1.0)
        {
          parent[Find(i)] = Find(j);
        }
      }
    }
    std::map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < n; i++)
    {
      groups[Find(i)].push_back(i);
    }
    for (const auto &[root, members] : groups)
    {
      (void)root;
      std::vector<nlohmann::json> signatures;
      int count = 0, feature_instances = 0;
      double length = 0.0;
      std::set<std::string> models, notes;
      nlohmann::json curvature_family, corner_family, span_cap_allowance;
      nlohmann::json legacy_contract = nlohmann::json::array();
      bool exact = true;
      for (const std::size_t i : members)
      {
        const Instance &instance = base.instances.at(base.order[i]);
        signatures.push_back(instance.signature);
        if (span_cap_allowance.is_null() && !instance.span_cap_allowance.is_null())
        {
          span_cap_allowance = instance.span_cap_allowance;
        }
        count += instance.count;
        feature_instances += instance.features;
        length += instance.length;
        models.insert(instance.models.begin(), instance.models.end());
        notes.insert(instance.notes.begin(), instance.notes.end());
        if (!instance.legacy_contract.is_null())
        {
          legacy_contract.push_back(instance.legacy_contract);
        }
        if (curvature_family.is_null() && !instance.curvature_family.is_null())
        {
          curvature_family = instance.curvature_family;
        }
        if (corner_family.is_null() && !instance.corner_family.is_null())
        {
          corner_family = instance.corner_family;
        }
        exact = exact && instance.exact;
      }
      const nlohmann::json representative = RepresentativeSignature(signatures);
      double spread = 0.0;
      for (const auto &signature : signatures)
      {
        spread =
            std::max(spread, SignatureDeviation(representative, signature).value_or(0.0));
      }
      const auto [key, hash] = SignatureKeyAndHash(representative, base.type);
      (void)key;
      // The keys of a near-matching cluster group's other members (block (b) DESIGN
      // section 4): the record's Hash is the lexicographically smallest member's.
      nlohmann::json near_keys = nlohmann::json::array();
      if (base.type == "SpatialEdgeCluster")
      {
        std::set<std::string> others;
        for (const auto &signature : signatures)
        {
          const std::string member = SignatureKeyAndHash(signature, base.type).second;
          if (member != hash)
          {
            others.insert(member);
          }
        }
        near_keys = others;
      }
      nlohmann::json record = {{"Dimension", 3},
                               {"Topology", base.type},
                               {"Status", models.empty() ? "Missing" : "Exact"},
                               {"Geometry", GeometryOf(base.type, representative)},
                               {"Interfaces", base.interfaces},
                               {"BoundaryCondition", base.law},
                               {"Hash", hash},
                               {"Signature", representative},
                               {"Instances", feature_instances},
                               {"DistinctSignatures", signatures.size()},
                               {"ParameterSpread", spread},
                               {"ExactParameters", exact}};
      if (!near_keys.empty())
      {
        record["NearKeys"] = near_keys;
      }
      if (!span_cap_allowance.is_null())
      {
        // The per-case span-cap allowance the cluster resolved (block (b) DESIGN section 3
        // (a) / A4): the planner passes its SpanCapOverR to the coupon generator.
        record["SpanCapAllowance"] = span_cap_allowance;
      }
      if (!models.empty())
      {
        nlohmann::json selected = nlohmann::json::array();
        for (const auto &model : models)
        {
          selected.push_back({{"Name", model}, {"Topology", base.type}, {"Weight", 1.0}});
        }
        record["SelectedModels"] = selected;
      }
      if (!curvature_family.is_null())
      {
        // The family selection (anchor, nodes, weights, rule) is the recorded
        // interpolation at the group's kappa (a group spans one signature tolerance);
        // Status "Interpolated" as on the two-dimensional path unless kappa is a node.
        record["CurvatureFamily"] = curvature_family;
        if (curvature_family["InterpolationRule"] != "exact")
        {
          record["Status"] = "Interpolated";
        }
      }
      if (!corner_family.is_null())
      {
        // The corner family selection (base, nodes, weights, rule) at the group's angle
        // (a group spans one signature tolerance).
        record["CornerFamily"] = corner_family;
        if (corner_family["InterpolationRule"] != "exact")
        {
          record["Status"] = "Interpolated";
        }
      }
      if (!notes.empty())
      {
        record["Notes"] = notes;
      }
      if (!legacy_contract.empty())
      {
        // Matched through the library's explicit legacy-contract alias (USER decision
        // 283): Status Exact (the legacy coupon is applied), flagged with the record.
        record["LegacyContract"] = legacy_contract;
      }
      requirements->AddFeatureRecord(record, count, length);
    }
  }
  requirements->SetIdentification(result.ToJson(requirements->CoordinateScale()));
  return result;
}

// Features-driven patch construction (SURFACE-RESPONSE-IDENTIFICATION.md (e)): every
// feature of the identification matched by signature becomes its patches, built from the
// feature's own portions, vertices and canonical frame; an unmatched feature is omitted
// alone (never an interface group, chain or neighbour); excluded segments carry no feature
// and are never corrected. Every patch records its feature and the segment portion it
// integrates so that the patch dry run can be audited against the manifest.
struct FeaturePatchSummary
{
  int matched_features = 0;
  int unmatched_features = 0;
  double matched_length = 0.0;
  double unmatched_length = 0.0;
  std::map<std::string, std::pair<int, double>> unmatched_by_type;  // count, length
  std::map<std::string, int> patches_by_type;
  // Matched pairs whose side geometry disagrees with their separation (a sample's foot on
  // the partner side farther than the tolerance from the signature's separation): omitted
  // with a warning, never patched silently.
  std::vector<int> inconsistent_features;
  double inconsistent_length = 0.0;
  // Curvature (design (b)7): features patched by a curvature family blend, and the
  // straight-like features with a bend whose first-order term was evaluated (both family
  // nodes at 1 / StraightBendRadiusOverR found) or could not be (a node missing: the
  // feature keeps its straight model; count and bent length recorded, never silent).
  int curved_family_features = 0;
  int corner_family_features = 0;  // corners patched by an angle-interpolated blend
  int first_order_features = 0;
  int first_order_missing_features = 0;
  double first_order_missing_turn = 0.0;  // |turn| (radians) left uncorrected
  std::set<std::string> first_order_missing_nodes;
  // A10 extended to the context (decision 282 section 3): context piece ends of the placed
  // contract-3 models verified to lie on device edges.
  std::size_t context_points_checked = 0;
};

}  // namespace

// Distance from a point to a straight segment a-b.
double SegmentDistance(const std::array<double, 3> &q, const std::array<double, 3> &a,
                       const std::array<double, 3> &b)
{
  const Point3D ab = Subtract(b, a), aq = Subtract(q, a);
  const double length2 = Dot(ab, ab);
  const double t = length2 > 0.0 ? std::clamp(Dot(aq, ab) / length2, 0.0, 1.0) : 0.0;
  return Norm(Subtract(q, Add(a, Scale(t, ab))));
}

// Distance from a point to the ARC of the fitted circle (centre C, radius rho) that the
// chord a-b subtends (block (b) DESIGN section 2 (a), A10-extended arc-aware): the point's
// projection into the circle's plane is tested against the chord's angular interval (the
// short way from a to b about C); inside it the distance is to the circle (radial residual
// and out-of-plane offset), outside it the distance to the nearer chord end. Exact and
// chord-independent: an arc context entry's end cut at a box face lies on the circle, up to
// the chord sagitta (1.6e-3..5.4e-2 R on the stage-2 census) from the device polyline. A
// degenerate chord (collinear with the centre) falls back to the straight distance.
double ArcChordDistance(const std::array<double, 3> &q, const std::array<double, 3> &a,
                        const std::array<double, 3> &b, const std::array<double, 3> &center,
                        double rho)
{
  const Point3D ra = Subtract(a, center), rb = Subtract(b, center);
  Point3D n = Cross(ra, rb);
  const double n_norm = Norm(n);
  if (n_norm <= 1.0e-12 * std::max(Norm(ra) * Norm(rb), 1.0e-300))
  {
    return SegmentDistance(q, a, b);
  }
  n = Scale(1.0 / n_norm, n);
  const Point3D rq = Subtract(q, center);
  const double out_of_plane = Dot(rq, n);
  const Point3D in_plane = Subtract(rq, Scale(out_of_plane, n));
  const bool within =
      Dot(Cross(ra, in_plane), n) >= 0.0 && Dot(Cross(in_plane, rb), n) >= 0.0;
  if (!within)
  {
    return std::min(Norm(Subtract(q, a)), Norm(Subtract(q, b)));
  }
  return std::hypot(std::abs(Norm(in_plane) - rho), out_of_plane);
}

// Distance from a point to the device perimeter of the identification: its segments' keys
// as straight chords, except that a segment lying on a fitted arc (Segments[].Arc) is read
// on the arc of its circle (ArcChordDistance). The A10 check extended to the context reads
// it for the placed context piece ends.
double DevicePerimeterDistance(const IdentificationResult &identification,
                               const std::array<double, 3> &q)
{
  double best = std::numeric_limits<double>::infinity();
  for (const auto &segment : identification.segments)
  {
    const auto &a = segment.key[0], &b = segment.key[1];
    const bool on_arc = segment.arc >= 0 &&
                        static_cast<std::size_t>(segment.arc) < identification.arcs.size();
    // Bounding-box rejection against the best distance so far (an arc bulges past its
    // chord's box by at most the recorded sagitta).
    const double bulge =
        on_arc ? identification.arcs[segment.arc].max_sagitta_over_R * identification.radius
               : 0.0;
    bool outside = false;
    for (int d = 0; d < 3 && !outside; d++)
    {
      outside = q[d] < std::min(a[d], b[d]) - best - bulge ||
                q[d] > std::max(a[d], b[d]) + best + bulge;
    }
    if (outside)
    {
      continue;
    }
    if (on_arc)
    {
      const auto &arc = identification.arcs[segment.arc];
      best = std::min(best, ArcChordDistance(q, a, b, arc.center, arc.radius));
    }
    else
    {
      best = std::min(best, SegmentDistance(q, a, b));
    }
  }
  return best;
}

namespace
{

// A pair's sample-to-partner distances must agree with its separation within twice the
// pair tolerance (the mean separation of a 5 % taper is within 5 % of every sample; the
// chord / inscribed readings add < 1 %).
constexpr double kPairPatchSeparationTolerance = 0.10;

// Longitudinal cells of a one-dimensional quadrature rule on the unit interval: the cell of
// point q is [sum of the weights of the points before it, that sum + w_q] in the order of
// increasing position, so the cells tile [0, 1] exactly and every cell contains its point
// (Gauss-Legendre). The translational surface mortar projects the device trace over the
// cell of each quadrature patch, mapped onto its portion.
std::vector<std::array<double, 2>>
LongitudinalQuadratureCells(const mfem::IntegrationRule &quadrature)
{
  std::vector<int> order(quadrature.GetNPoints());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int a, int b)
            { return quadrature.IntPoint(a).x < quadrature.IntPoint(b).x; });
  std::vector<std::array<double, 2>> cells(quadrature.GetNPoints());
  double cumulative = 0.0;
  for (const int q : order)
  {
    const double next = cumulative + quadrature.IntPoint(q).weight;
    cells[q] = {cumulative, next};
    cumulative = next;
  }
  MFEM_VERIFY(std::abs(cumulative - 1.0) < 1.0e-12,
              "Longitudinal quadrature weights do not sum to the portion length!");
  cells[order.back()][1] = 1.0;
  for (int q = 0; q < quadrature.GetNPoints(); q++)
  {
    const double x = quadrature.IntPoint(q).x;
    MFEM_VERIFY(cells[q][0] <= x && x <= cells[q][1],
                "A longitudinal quadrature cell does not contain its point!");
  }
  return cells;
}

// The smallest |cos| between a translational patch's AxisW and the tangent of the segment
// its sample lies on. AxisW = AxisU x AxisV is parallel to the sample's own segment for an
// isolated edge, but for a pair or stack AxisU points at the sample's closest foot on the
// PARTNER side, so AxisW follows the partner's chord (or the sample-to-vertex direction
// where the foot is clamped at a joint): the identification classifies slow tapers and
// sub-noise polyline bends as pairs (1 - |cos| of 1e-6..1e-3 there: a 6 um / 100 um taper
// 4.5e-6, a 1.6 deg chord joint 3.9e-4). The tolerance is the pair regime's own facing
// threshold (the paired-edge topology's Dot(axis_u, direction) > 0.95 in the 2D and legacy
// 3D groupings, ~18 deg); anything below it is not a translational frame and fails closed.
constexpr double kLongitudinalAxisCosineTolerance = 0.95;

// The longitudinal cell of one quadrature patch on a straight portion [a, b] of a segment
// parametrised from p0 along `tangent`: offsets from the patch origin (at parameter t_q)
// along the patch's AxisW (= AxisU x AxisV), mesh units, ordered begin <= end. The arc
// cell is projected onto AxisW by Dot(tangent, axis_w) (its sign orients the cell, its
// magnitude shortens it): the cross-section perpendicular to AxisW at the projected offset
// contains the segment point at that arc offset (the sample's foot lies in the plane
// perpendicular to AxisW through the origin), so the slices sweep exactly the sample's
// cell of its own segment whether or not AxisW is parallel to it.
std::array<double, 2> LongitudinalCellOffsets(const std::array<double, 2> &unit_cell,
                                              double a, double b, double t_q,
                                              const Point3D &tangent, const Point3D &axis_w)
{
  const double projection = Dot(tangent, axis_w);
  MFEM_VERIFY(std::abs(projection) >= kLongitudinalAxisCosineTolerance,
              "A translational patch's AxisW must be within |cos| >= "
                  << kLongitudinalAxisCosineTolerance << " of its segment (found "
                  << std::abs(projection) << ")!");
  const double begin = projection * (a + (b - a) * unit_cell[0] - t_q);
  const double end = projection * (a + (b - a) * unit_cell[1] - t_q);
  return {std::min(begin, end), std::max(begin, end)};
}

FeaturePatchSummary BuildFeaturePatches(
    const ProcessLibrary &library, const IdentificationResult &identification,
    const std::vector<EdgeSegment3D> &framed_segments,
    const mfem::IntegrationRule &quadrature, const AutomaticResponseRequirements &describer,
    AutomaticResponseDiagnostics *diagnostics, ResponseCorrectionData &result,
    const std::map<int, FeatureCurvatureMatch> &curved_matches = {},
    const std::map<int, FeatureCornerMatch> &corner_matches = {})
{
  FeaturePatchSummary summary;
  const double R = library.matching_radius;
  std::map<std::size_t, const EdgeSegment3D *> framed;
  for (const auto &segment : framed_segments)
  {
    framed.emplace(segment.geometry_index, &segment);
  }
  std::map<std::string, std::size_t> model_by_name;
  for (std::size_t i = 0; i < library.models.size(); i++)
  {
    model_by_name.emplace(library.models[i].name, i);
  }

  // A feature portion in the parametrisation of its framed segment (from p0); the manifest
  // portion [s0, s1) runs from the segment's canonical key origin.
  struct FramedPortion
  {
    const EdgeSegment3D *segment = nullptr;
    double a = 0.0, b = 0.0;
    std::size_t geometry_index = 0;
    double s0 = 0.0, s1 = 0.0;
    int side = 0;
    double turn = 0.0;  // signed windowed turn toward the metal (radians)
    int stretch = -1;   // the chain stretch of the portion (IdentifiedPortion::stretch)
  };
  auto Frame = [&](const IdentifiedPortion &portion)
  {
    const auto it = framed.find(portion.segment);
    MFEM_VERIFY(it != framed.end(),
                "A feature portion lies on a segment without an edge frame!");
    const EdgeSegment3D &segment = *it->second;
    const bool forward = identification.segments[portion.segment].key[0] == segment.p0;
    FramedPortion fp;
    fp.segment = &segment;
    fp.a =
        std::clamp(forward ? portion.s0 : segment.length - portion.s1, 0.0, segment.length);
    fp.b =
        std::clamp(forward ? portion.s1 : segment.length - portion.s0, 0.0, segment.length);
    fp.geometry_index = portion.segment;
    fp.s0 = portion.s0;
    fp.s1 = portion.s1;
    fp.side = portion.side;
    fp.turn = portion.turn;
    fp.stretch = portion.stretch;
    return fp;
  };
  auto IsPec = [](const EdgeSegment3D &segment)
  { return segment.boundary_condition.type == MetalBoundaryConditionType::PEC; };

  // One runtime model per (library model, target interfaces by slot); a family blend
  // (curvature family, corner family) is one runtime model per (blend name, slots): the
  // base model's basis and interfaces with matrices = the weighted sum of the nodes'
  // matrices (Lagrange weights may be negative, so the combination is formed on the
  // matrices, as on the two-dimensional path). A curvature blend rescales every node to the
  // anchor's coupon depth (per unit edge length); a corner blend has no depth (one patch of
  // weight one per corner) and combines the coupons' matrices as they are.
  struct FamilyBlend
  {
    std::string name;
    std::string topology;
    std::size_t base = 0;  // the runtime model's basis / interfaces / references
    std::vector<LibrarySelection::WeightedModel> nodes;
    bool per_unit_length = true;
    std::string label;
    // A corner blend's basis constructed at the feature's angle (the trace basis rule).
    const ConstructedCornerTraceBasis *constructed = nullptr;
  };
  std::map<std::pair<std::string, std::string>, int> runtime_models;
  int next_model_index = 1;
  auto RuntimeModel =
      [&](std::size_t model_index,
          const std::map<int, std::map<InterfaceDielectric, int>> &targets_by_slot,
          const FamilyBlend *blend = nullptr)
  {
    std::string key;
    for (const auto &[slot, targets] : targets_by_slot)
    {
      key += fmt::format("{}:", slot);
      for (const auto &[type, target] : targets)
      {
        key += fmt::format("{}={},", ToString(type), target);
      }
      key += ";";
    }
    const auto &source = library.models[model_index];
    MFEM_VERIFY(!blend || blend->base == model_index,
                "A family blend must be built on its base model!");
    auto [it, inserted] =
        runtime_models.emplace(std::make_pair(blend ? blend->name : source.name, key), 0);
    if (inserted)
    {
      auto model = source.response;
      model.idx = next_model_index++;
      model.name = blend ? blend->name : source.name;
      model.topology = blend ? blend->topology : TopologyName(source.topology);
      if (blend)
      {
        MFEM_VERIFY(!blend->per_unit_length || source.coupon_depth > 0.0,
                    "Curvature interpolation requires CouponDepth on the straight anchor!");
        std::ostringstream nodes;
        for (const auto &node : blend->nodes)
        {
          const auto &coupon = library.models[node.index];
          double weight = node.weight;
          if (blend->per_unit_length)
          {
            MFEM_VERIFY(coupon.coupon_depth > 0.0,
                        "Curvature interpolation requires CouponDepth on every coupon!");
            weight *= source.coupon_depth / coupon.coupon_depth;
          }
          model.blend.push_back({weight, coupon.response.fabricated_matrix,
                                 coupon.response.thin_matrix,
                                 coupon.response.fabricated_surface_matrix,
                                 coupon.response.thin_surface_matrix});
          nodes << (model.blend.size() == 1 ? "" : ", ") << coupon.name << " x "
                << std::setprecision(6) << node.weight;
        }
        Mpi::Print(" {} {}: {}\n", blend->label, blend->name, nodes.str());
        if (blend->constructed)
        {
          const auto &constructed = *blend->constructed;
          model.constructed_basis_points = constructed.knots;
          model.constructed_trace_vertices.clear();
          for (const auto &vertex : constructed.vertices)
          {
            ResponseModelData::ConstructedTraceVertex data;
            data.point = vertex.point;
            data.basis = vertex.basis >= 0 ? vertex.basis + 1 : 0;
            data.conductor = vertex.basis >= 0 && constructed.zero[vertex.basis] ? 1 : 0;
            data.parent_a = vertex.basis < 0 ? vertex.parent_a + 1 : 0;
            data.parent_b = vertex.basis < 0 ? vertex.parent_b + 1 : 0;
            data.weight_a = vertex.basis < 0 ? vertex.weight_a : 0.0;
            model.constructed_trace_vertices.push_back(data);
          }
          model.constructed_trace_triangles = constructed.triangles;
        }
      }
      MapLibraryInterfaces(source, targets_by_slot, model);
      result.models.push_back(std::move(model));
      it->second = result.models.back().idx;
    }
    return it->second;
  };
  auto CurvatureBlend = [](const FeatureCurvatureMatch &match)
  {
    return FamilyBlend{
        match.name, match.topology,           match.selection.anchor, match.selection.nodes,
        true,       "Curvature interpolation"};
  };
  auto CornerBlend = [](const FeatureCornerMatch &match)
  {
    return FamilyBlend{match.name,
                       match.topology,
                       match.base,
                       match.nodes,
                       false,
                       "Corner interpolation",
                       match.constructed ? &*match.constructed : nullptr};
  };

  auto Emit = [&](ResponsePatchData patch, std::size_t model_index, int runtime,
                  const IdentifiedFeature &feature, double model_weight = 1.0)
  {
    patch.model = runtime;
    patch.provenance.feature = feature.id;
    patch.provenance.model_weight = model_weight;
    if (diagnostics)
    {
      diagnostics->boundary_law_verified &=
          IsBoundaryLawVerified(library.models[model_index]);
    }
    result.patches.push_back(std::move(patch));
    summary.patches_by_type[feature.type]++;
  };

  // First-order curvature term of a straight-like portion (design (b)7, the linear rule
  // keyed to StraightBendRadiusOverR evaluated where the curvature is): the portion's
  // local kappa = R |turn| / length (dimensionless, <= 1 / StraightBendRadiusOverR on a
  // straight-like chain) splits every quadrature point between the anchor (model weight
  // 1 - a) and the family node at the first-order kappa of the turn's convexity (weight
  // a = kappa / FirstOrderKappa): two co-located positive patches whose assembled matrices
  // are the linear interpolant. For a pair the kappa is that of the INNER edge (the
  // family's Kappa = R / rho_inner): an outer-side sample reads rho_inner = rho - s.
  struct FirstOrderSplit
  {
    std::size_t node = 0;  // library index of the first-order node
    int runtime = 0;
    double a = 0.0;  // weight of the node (0: anchor only)
  };
  // Longitudinal quadrature over one framed portion: one patch per quadrature point,
  // weight = portion length x quadrature weight x side factor / coupon depth (times the
  // model weight of a first-order split); the patch's longitudinal cell is its quadrature
  // cell of the portion (the cells tile the portion; a split's co-located patches and the
  // sides of a pair carry the full cell whatever their weight factors).
  const auto quadrature_cells = LongitudinalQuadratureCells(quadrature);
  auto Quadrature = [&](const FramedPortion &fp, double side_factor,
                        std::size_t model_index, int runtime,
                        const IdentifiedFeature &feature, const auto &Place,
                        const std::optional<FirstOrderSplit> &split = std::nullopt)
  {
    const auto &model = library.models[model_index];
    MFEM_VERIFY(model.coupon_depth > 0.0,
                "Three-dimensional response correction requires CouponDepth for every "
                "selected fabrication-process response model (\""
                    << model.name << "\")!");
    std::vector<std::tuple<std::size_t, int, double>> terms = {{model_index, runtime, 1.0}};
    if (split && split->a > 0.0)
    {
      terms = {{model_index, runtime, 1.0 - split->a},
               {split->node, split->runtime, split->a}};
      MFEM_VERIFY(library.models[split->node].coupon_depth > 0.0,
                  "Curvature interpolation requires CouponDepth on every coupon!");
    }
    for (int q = 0; q < quadrature.GetNPoints(); q++)
    {
      const auto &ip = quadrature.IntPoint(q);
      const double t = fp.a + (fp.b - fp.a) * ip.x;
      const Point3D point = Interpolate(*fp.segment, t);
      for (const auto &[term_model, term_runtime, model_weight] : terms)
      {
        if (model_weight <= 0.0)
        {
          continue;
        }
        const auto &term = library.models[term_model];
        ResponsePatchData patch;
        Place(point, patch);
        patch.conductor_references = term.conductor_references;
        patch.weight =
            model_weight * side_factor * (fp.b - fp.a) * ip.weight / term.coupon_depth;
        patch.longitudinal_cell = LongitudinalCellOffsets(
            quadrature_cells[q], fp.a, fp.b, t, fp.segment->tangent, patch.axis_w);
        patch.provenance.segment = static_cast<int>(fp.geometry_index);
        patch.provenance.s0 = fp.s0;
        patch.provenance.s1 = fp.s1;
        patch.provenance.stretch = fp.stretch;
        // The side's own edge point (the sample on its segment) relative to the placed
        // origin (a single edge: 0; a pair: the midline; a stack: the first side).
        patch.provenance.edge_offset = Dot(Subtract(point, patch.origin), patch.axis_u);
        patch.provenance.quadrature_weight = ip.weight;
        patch.provenance.side_factor = side_factor;
        patch.provenance.coupon_depth = term.coupon_depth;
        Emit(std::move(patch), term_model, term_runtime, feature, model_weight);
      }
    }
  };
  // The first-order split of a portion of a straight-like feature: the node of the turn's
  // convexity (for a pair: of the edge the model's first edge e1 lands on — the far side
  // turns in the opposite sense) at a = kappa_inner / FirstOrderKappa, clamped to [0, 1]
  // (roundoff at the straight-like boundary); nullopt (anchor only, recorded) when the node
  // is missing or the portion carries no turn.
  const double first_order_kappa = 1.0 / kStraightBendRadiusOverRadius;
  auto SplitOf = [&](const FramedPortion &fp, const FirstOrderNodes &nodes,
                     const std::map<int, std::map<InterfaceDielectric, int>> &slots,
                     bool e1_side, bool strip, double separation,
                     std::size_t anchor_index) -> std::optional<FirstOrderSplit>
  {
    const double length = fp.b - fp.a;
    if (std::abs(fp.turn) <= 0.0 || length <= 0.0)
    {
      return std::nullopt;
    }
    const bool convex_e1 = e1_side ? fp.turn > 0.0 : fp.turn < 0.0;
    // Inner edge of a pair: the convex edge of a gap, the concave edge of a strip.
    const bool inner = separation <= 0.0 || ((fp.turn > 0.0) != strip);
    const double rho_side = length / std::abs(fp.turn);
    const double rho_inner = inner ? rho_side : rho_side - separation;
    if (rho_inner <= 0.0)
    {
      return std::nullopt;  // an outer edge tighter than the separation: not a pair bend
    }
    const double kappa = R / rho_inner;
    const auto node = convex_e1 ? nodes.convex : nodes.concave;
    if (!node)
    {
      summary.first_order_missing_turn += std::abs(fp.turn);
      summary.first_order_missing_nodes.insert(
          TopologyName(CurvedTopologyOf(library.models[anchor_index].topology)) +
          std::string(convex_e1 ? " convex" : " concave"));
      return std::nullopt;
    }
    FirstOrderSplit split;
    split.node = *node;
    split.runtime = RuntimeModel(*node, slots);
    split.a = std::clamp(kappa / first_order_kappa, 0.0, 1.0);
    return split;
  };

  // Closest point of p on the sub-segments of the other side of a pair.
  struct Foot
  {
    Point3D point{};
    const EdgeSegment3D *segment = nullptr;
  };
  auto ClosestFoot = [&](const Point3D &p, const std::vector<FramedPortion> &side)
  {
    Foot best;
    double best_distance = mfem::infinity();
    for (const auto &fp : side)
    {
      const double t =
          std::clamp(Dot(Subtract(p, fp.segment->p0), fp.segment->tangent), fp.a, fp.b);
      const Point3D q = Interpolate(*fp.segment, t);
      const double distance = Distance(p, q);
      if (distance < best_distance)
      {
        best_distance = distance;
        best = {q, fp.segment};
      }
    }
    MFEM_VERIFY(best.segment, "A paired feature has no partner side!");
    return best;
  };

  for (const auto &feature : identification.features)
  {
    double feature_length = 0.0;
    for (const auto &portion : feature.portions)
    {
      feature_length += portion.s1 - portion.s0;
    }
    const bool vertex_feature = feature.type == "ConvexCorner" ||
                                feature.type == "ConcaveCorner" ||
                                feature.type == "Endpoint" || feature.type == "Junction";
    if (diagnostics && feature.bend_radius_over_R && *feature.bend_radius_over_R > 0.0)
    {
      diagnostics->maximum_curvature_ratio =
          std::max(diagnostics->maximum_curvature_ratio, 1.0 / *feature.bend_radius_over_R);
    }
    if (!feature.matched_model)
    {
      summary.unmatched_features++;
      summary.unmatched_length += feature_length;
      auto &entry = summary.unmatched_by_type[feature.type];
      entry.first++;
      entry.second += feature_length;
      if (diagnostics && vertex_feature)
      {
        for (const auto &portion : feature.portions)
        {
          const double length = portion.s1 - portion.s0;
          diagnostics->matched_corner_neighborhood_length += length;
          for (const auto &[type, target] : framed.at(portion.segment)->targets)
          {
            (void)type;
            diagnostics->matched_corner_neighborhood_length_by_interface[target] += length;
          }
        }
      }
      continue;
    }
    // The feature's model: a library model by name, a curvature-family blend on its
    // anchor (the straight analogue) or a corner-family blend on its nearest node (the
    // runtime model carries the blended matrices).
    const auto curved_match = curved_matches.find(feature.id);
    const auto corner_match = corner_matches.find(feature.id);
    std::optional<FamilyBlend> family;
    const FeatureCurvatureMatch *curved_blend =
        curved_match != curved_matches.end() ? &curved_match->second : nullptr;
    std::size_t model_index = 0;
    if (curved_blend)
    {
      MFEM_VERIFY(curved_blend->name == *feature.matched_model,
                  "Curved feature " << feature.id << " matched \"" << *feature.matched_model
                                    << "\" but its curvature family is \""
                                    << curved_blend->name << "\"!");
      family = CurvatureBlend(*curved_blend);
      model_index = curved_blend->selection.anchor;
      summary.curved_family_features++;
    }
    else if (corner_match != corner_matches.end())
    {
      MFEM_VERIFY(corner_match->second.name == *feature.matched_model,
                  "Corner feature " << feature.id << " matched \"" << *feature.matched_model
                                    << "\" but its corner family is \""
                                    << corner_match->second.name << "\"!");
      family = CornerBlend(corner_match->second);
      model_index = corner_match->second.base;
      summary.corner_family_features++;
    }
    else
    {
      const auto model_it = model_by_name.find(*feature.matched_model);
      MFEM_VERIFY(model_it != model_by_name.end(),
                  "Matched library model \"" << *feature.matched_model << "\" not found!");
      model_index = model_it->second;
    }
    const auto &model = library.models[model_index];
    summary.matched_features++;
    summary.matched_length += feature_length;

    std::vector<FramedPortion> portions;
    std::set<std::map<InterfaceDielectric, int>> target_maps;
    bool all_pec = true;
    for (const auto &portion : feature.portions)
    {
      portions.push_back(Frame(portion));
      target_maps.insert(portions.back().segment->targets);
      all_pec = all_pec && IsPec(*portions.back().segment);
      if (diagnostics)
      {
        const double length = portion.s1 - portion.s0;
        diagnostics->matched_length += length;
        for (const auto &[type, target] : portions.back().segment->targets)
        {
          (void)type;
          diagnostics->matched_length_by_interface[target] += length;
        }
      }
    }
    MFEM_VERIFY(!portions.empty(), "A matched feature claims no perimeter portion!");
    // Interface slots in the sorted order of the distinct target maps (slot 0 = the first).
    std::map<int, std::map<InterfaceDielectric, int>> targets_by_slot;
    for (const auto &targets : target_maps)
    {
      targets_by_slot.emplace(static_cast<int>(targets_by_slot.size()), targets);
    }
    const int runtime =
        RuntimeModel(model_index, targets_by_slot, family ? &*family : nullptr);
    const Point3D n = feature.axes[2];
    // First-order curvature term of a straight-like feature with a bend (a straight
    // model on portions with a nonzero turn): the family nodes at the first-order kappa.
    const FeatureCurvatureMatch *blend = curved_blend;
    const bool straight_like_bend =
        !blend && feature.bend_radius_over_R &&
        (feature.type == "IsolatedEdge" || feature.type == "SameConductorGap" ||
         feature.type == "DifferentConductorGap" || feature.type == "SameConductorStrip") &&
        std::any_of(portions.begin(), portions.end(),
                    [](const FramedPortion &fp) { return fp.turn != 0.0; });
    FirstOrderNodes first_order_nodes;
    if (straight_like_bend)
    {
      first_order_nodes = FindFirstOrderNodes(library, model_index);
      summary.first_order_features++;
    }
    const double missing_turn_before = summary.first_order_missing_turn;

    if (feature.type == "IsolatedEdge" || feature.type == "CurvedEdge")
    {
      for (const auto &fp : portions)
      {
        const EdgeSegment3D &segment = *fp.segment;
        const auto split = straight_like_bend
                               ? SplitOf(fp, first_order_nodes, targets_by_slot, true,
                                         false, 0.0, model_index)
                               : std::nullopt;
        Quadrature(
            fp, 1.0, model_index, runtime, feature,
            [&](const Point3D &point, ResponsePatchData &patch)
            {
              patch.origin = point;
              patch.axis_u = segment.axis_u;
              patch.axis_v = segment.axis_v;
              patch.axis_w = Normalize(Cross(segment.axis_u, segment.axis_v));
              patch.maxwell_reference_is_pec = IsPec(segment);
              patch.maxwell_conductor_anchors = {patch.maxwell_reference_is_pec
                                                     ? Add(point, Scale(-R, segment.axis_u))
                                                     : point};
            },
            split);
      }
    }
    else if (feature.type == "SameConductorGap" ||
             feature.type == "DifferentConductorGap" ||
             feature.type == "SameConductorStrip" ||
             feature.type == "CurvedSameConductorGap" ||
             feature.type == "CurvedDifferentConductorGap" ||
             feature.type == "CurvedSameConductorStrip")
    {
      // Two sides (the identification's side labels in increasing lateral offset; both
      // sides may lie on one perimeter chain, e.g. the slot of a shorted CPW); the model's
      // first edge is side 0 (side 1 for chirality -1: the canonical orientation is the
      // mirror). Each side carries half of the longitudinal measure (the mean of the two
      // sides: exact for a straight pair, the centreline for a concentric one).
      std::map<int, std::vector<FramedPortion>> sides;
      for (const auto &fp : portions)
      {
        sides[fp.side].push_back(fp);
      }
      MFEM_VERIFY(sides.size() == 2 && sides.count(0) && sides.count(1),
                  "A paired feature must claim exactly two sides!");
      std::vector<std::pair<double, const std::vector<FramedPortion> *>> ordered = {
          {0.0, &sides.at(0)}, {1.0, &sides.at(1)}};
      if (feature.chirality < 0)
      {
        std::swap(ordered[0], ordered[1]);
      }
      const bool strip =
          model.topology == LibraryTopology::SAME_CONDUCTOR_STRIP ||
          model.topology == LibraryTopology::CURVED_SAME_CONDUCTOR_STRIP ||
          (blend &&
           blend->topology == TopologyName(LibraryTopology::CURVED_SAME_CONDUCTOR_STRIP));
      MFEM_VERIFY(model.conductor_references.size() <= 2,
                  "A paired-edge response model requires at most two conductor "
                  "references!");
      // Consistency of the two sides with the signature's separation before any patch.
      const double separation = feature.signature.at("SeparationOverR").get<double>() * R;
      bool consistent = true;
      for (int k = 0; consistent && k < 2; k++)
      {
        for (const auto &fp : *ordered[k].second)
        {
          for (const double t : {fp.a, 0.5 * (fp.a + fp.b), fp.b})
          {
            const Foot foot =
                ClosestFoot(Interpolate(*fp.segment, t), *ordered[1 - k].second);
            if (std::abs(Distance(Interpolate(*fp.segment, t), foot.point) - separation) >
                kPairPatchSeparationTolerance * separation)
            {
              consistent = false;
              break;
            }
          }
          if (!consistent)
          {
            break;
          }
        }
      }
      if (!consistent)
      {
        summary.inconsistent_features.push_back(feature.id);
        summary.inconsistent_length += feature_length;
        Mpi::Warning(
            "Feature {:d} ({}, {:.6e} length units): the two sides do not face each "
            "other at its separation {:.6e}; the pair is omitted (not patched).\n",
            feature.id, feature.type, feature_length, separation);
        continue;
      }
      for (int k = 0; k < 2; k++)
      {
        const auto &side = *ordered[k].second;
        const auto &other = *ordered[1 - k].second;
        for (const auto &fp : side)
        {
          const auto split = straight_like_bend
                                 ? SplitOf(fp, first_order_nodes, targets_by_slot, k == 0,
                                           strip, separation, model_index)
                                 : std::nullopt;
          Quadrature(
              fp, 0.5, model_index, runtime, feature,
              [&](const Point3D &point, ResponsePatchData &patch)
              {
                const Foot foot = ClosestFoot(point, other);
                const Point3D e1 = k == 0 ? point : foot.point;
                const Point3D e2 = k == 0 ? foot.point : point;
                const EdgeSegment3D &first = k == 0 ? *fp.segment : *foot.segment;
                const EdgeSegment3D &second = k == 0 ? *foot.segment : *fp.segment;
                patch.origin = Scale(0.5, Add(e1, e2));
                patch.axis_u = Normalize(Subtract(e2, e1));
                patch.axis_v = Normalize(Add(first.axis_v, second.axis_v));
                patch.axis_w = Normalize(Cross(patch.axis_u, patch.axis_v));
                patch.maxwell_reference_is_pec = IsPec(first) && IsPec(second);
                if (!patch.maxwell_reference_is_pec)
                {
                  patch.maxwell_conductor_anchors = {e1};
                }
                else if (model.conductor_references.size() > 1)
                {
                  // The physical edge points as local conductor anchors so that the
                  // Maxwell quadrature spans only the dielectric gap.
                  patch.maxwell_conductor_anchors = {e1, e2};
                }
                else if (strip)
                {
                  patch.maxwell_conductor_anchors = {patch.origin};
                }
                else
                {
                  patch.maxwell_conductor_anchors = {Add(e1, Scale(-R, first.axis_u))};
                }
              },
              split);
        }
      }
    }
    else if (feature.type == "ParallelEdgeCluster")
    {
      // Sides in increasing lateral offset (the identification's labels; reversed for
      // chirality -1) = the model's edges by offset; conductor labels by first appearance
      // in that order.
      std::map<int, std::vector<FramedPortion>> sides;
      for (const auto &fp : portions)
      {
        sides[fp.side].push_back(fp);
      }
      MFEM_VERIFY(sides.size() >= 3 && sides.begin()->first == 0 &&
                      sides.rbegin()->first + 1 == static_cast<int>(sides.size()),
                  "A parallel-edge cluster feature must claim at least three sides!");
      std::vector<std::pair<double, const std::vector<FramedPortion> *>> ordered;
      for (const auto &[side, side_portions] : sides)
      {
        ordered.emplace_back(static_cast<double>(side), &side_portions);
      }
      if (feature.chirality < 0)
      {
        std::reverse(ordered.begin(), ordered.end());
      }
      MFEM_VERIFY(model.cluster_edges.empty() ||
                      model.cluster_edges.size() == ordered.size(),
                  "A parallel-edge cluster model's edge count differs from the feature!");
      std::vector<std::size_t> reference_sides;  // first side of every conductor label
      {
        std::set<int> conductors;
        for (std::size_t k = 0; k < ordered.size(); k++)
        {
          if (conductors.insert(ordered[k].second->front().segment->conductor).second)
          {
            reference_sides.push_back(k);
          }
        }
      }
      MFEM_VERIFY(model.conductor_references.size() == reference_sides.size(),
                  "Parallel-edge cluster model \""
                      << model.name
                      << "\" requires one conductor reference per canonical conductor!");
      Point3D axis_v{};
      for (const auto &fp : portions)
      {
        axis_v = Add(axis_v, Scale(fp.b - fp.a, fp.segment->axis_v));
      }
      axis_v = Normalize(axis_v);
      // Local frame at every sample, as for the pairs: the origin is the sample's foot on
      // the canonical first side, the lateral axis points from there to its foot on the
      // last side, and the conductor anchors are the feet on the first side of every
      // conductor. (A feature-wide frame — the first run's tangent and the canonical
      // lateral — placed every patch of a stack that turns, e.g. a flux line meandering
      // through straight-like bends of radius 98 R over 2 mm, with its lateral axis up to
      // 49 deg off the local perpendicular; the placement audit A10 on DS-SCT-001 found it.
      // On a straight stack the local frame is the feature frame.)
      const std::vector<FramedPortion> &first_side = *ordered.front().second;
      const std::vector<FramedPortion> &last_side = *ordered.back().second;
      const double side_factor = 1.0 / static_cast<double>(ordered.size());
      for (std::size_t k = 0; k < ordered.size(); k++)
      {
        for (const auto &fp : *ordered[k].second)
        {
          Quadrature(
              fp, side_factor, model_index, runtime, feature,
              [&](const Point3D &point, ResponsePatchData &patch)
              {
                patch.origin = k == 0 ? point : ClosestFoot(point, first_side).point;
                const Point3D far = ClosestFoot(patch.origin, last_side).point;
                patch.axis_u = Normalize(Subtract(far, patch.origin));
                patch.axis_v = axis_v;
                patch.axis_w = Normalize(Cross(patch.axis_u, patch.axis_v));
                patch.maxwell_reference_is_pec = all_pec;
                for (const std::size_t reference : reference_sides)
                {
                  patch.maxwell_conductor_anchors.push_back(
                      reference == 0
                          ? patch.origin
                          : ClosestFoot(patch.origin, *ordered[reference].second).point);
                }
              });
        }
      }
    }
    else if (vertex_feature)
    {
      // One patch in the feature's canonical frame (x = the first arm, design (b) 4). A
      // legacy junction model stores its arms as absolute angles (arm 0 along its x axis,
      // counterclockwise): map its canonical first arm onto the feature's and its canonical
      // orientation onto the feature's (sigma = -1 mirrors the frame).
      ResponsePatchData patch;
      patch.origin = feature.origin;
      patch.axis_u = feature.axes[0];
      patch.axis_v = feature.axes[1];
      patch.axis_w = n;
      MFEM_VERIFY(Norm(patch.axis_u) > 0.0 && Norm(patch.axis_v) > 0.0,
                  "A vertex feature without a frame cannot be patched!");
      if (feature.type == "Junction" && !model.arm_angles.empty())
      {
        std::vector<double> sorted = model.arm_angles, differences;
        std::sort(sorted.begin(), sorted.end());
        for (std::size_t i = 0; i < sorted.size(); i++)
        {
          const double next =
              i + 1 < sorted.size() ? sorted[i + 1] : sorted[0] + 2.0 * std::acos(-1.0);
          differences.push_back((next - sorted[i]) * 180.0 / std::acos(-1.0));
        }
        JunctionCanonicalOrder model_order;
        CanonicalJunctionSignature({}, "", differences, {}, &model_order);
        const double theta = sorted[model_order.first_arm];
        const Point3D D = feature.axes[0];
        const Point3D W = Normalize(Cross(n, D));
        const bool device_ccw = Dot(feature.axes[1], W) > 0.0;
        const bool model_ccw = !model_order.reversed;
        const double sigma = device_ccw == model_ccw ? 1.0 : -1.0;
        patch.axis_u = Normalize(
            Subtract(Scale(std::cos(theta), D), Scale(sigma * std::sin(theta), W)));
        patch.axis_v = Scale(sigma, Normalize(Cross(n, patch.axis_u)));
      }
      patch.conductor_references = model.conductor_references;
      patch.weight = 1.0;
      patch.maxwell_reference_is_pec = all_pec;
      const std::array<Point3D, 3> axes = {patch.axis_u, patch.axis_v, patch.axis_w};
      for (const auto &reference : patch.conductor_references)
      {
        patch.maxwell_conductor_anchors.push_back(
            patch.maxwell_reference_is_pec || feature.type == "ConvexCorner" ||
                    feature.type == "ConcaveCorner"
                ? TransformLocalPoint(patch.origin, axes, reference)
                : patch.origin);
      }
      Emit(std::move(patch), model_index, runtime, feature);
    }
    else if (feature.type == "SpatialEdgeCluster")
    {
      // A model keyed by its Signature is built in the canonical frame of that Signature
      // (M = identity; its Edges, when stored, were verified against the Signature's
      // portions at load — VerifySpatialEdgesInSignatureFrame). A legacy model's edges are
      // expressed in its own frame; their canonical frame M and the feature's canonical
      // frame F have the same serialisation, so a model-frame point m maps to F.origin +
      // F.axes^T M.axes (m - M.origin).
      Point3D model_origin{};
      std::array<Point3D, 3> model_axes = {Point3D{1.0, 0.0, 0.0}, Point3D{0.0, 1.0, 0.0},
                                           Point3D{0.0, 0.0, 1.0}};
      if (!model.spatial_edges.empty() && !model.identification_signature)
      {
        const auto canonical = ModelClusterSignature(model, R, describer);
        MFEM_VERIFY(canonical, "Unable to canonicalise a spatial edge-cluster model!");
        model_origin = canonical->origin;
        model_axes = canonical->axes;
      }
      ResponsePatchData patch;
      std::array<Point3D, 3> axes{};
      for (int j = 0; j < 3; j++)
      {
        for (int k = 0; k < 3; k++)
        {
          axes[j] = Add(axes[j], Scale(model_axes[k][j], feature.axes[k]));
        }
      }
      patch.origin = feature.origin;
      for (int k = 0; k < 3; k++)
      {
        patch.origin = Subtract(patch.origin,
                                Scale(Dot(model_axes[k], model_origin), feature.axes[k]));
      }
      patch.axis_u = axes[0];
      patch.axis_v = axes[1];
      patch.axis_w = axes[2];
      patch.conductor_references = model.conductor_references;
      patch.weight = 1.0;
      patch.maxwell_reference_is_pec = all_pec;
      for (const auto &reference : patch.conductor_references)
      {
        patch.maxwell_conductor_anchors.push_back(
            TransformLocalPoint(patch.origin, axes, reference));
      }
      // The cluster's claimed portions (mesh segment and global ends): the ownership
      // record's continuation class (a translational stretch inside the box that continues
      // one of them through the claim cut).
      for (const auto &portion : feature.portions)
      {
        const auto fp = Frame(portion);
        patch.provenance.claims.push_back({static_cast<int>(fp.geometry_index),
                                           Interpolate(*fp.segment, fp.a),
                                           Interpolate(*fp.segment, fp.b)});
      }
      // A contract-3 (device-plan) model: its support box and continuation chain in the
      // patch frame (rule B4, the placement's vertex ownership), and the A10 check extended
      // to the context (decision 282 section 3): every context piece of the model's
      // Signature, placed by the patch frame, lies on a device run within the signature
      // tolerance — a mis-keyed library or a frame defect fails closed here.
      if (model.identification_signature &&
          model.identification_signature->contains("Box") &&
          feature.legacy_contract.is_null())
      {
        const auto &signature = *model.identification_signature;
        patch.provenance.support_box = signature.at("Box").get<std::array<double, 4>>();
        patch.provenance.has_support_box = true;
        MFEM_VERIFY(feature.signature.contains("Context") &&
                        signature.value("Context", nlohmann::json::array()) ==
                            feature.signature["Context"],
                    "The Context of model \"" << model.name
                                              << "\" differs from the matched feature's!");
        const double tolerance = kSignatureParameterToleranceOverRadius * R;
        auto Global = [&](double x, double y)
        {
          return Add(patch.origin,
                     Add(Scale(x * R, patch.axis_u), Scale(y * R, patch.axis_v)));
        };
        std::size_t checked = 0;
        for (const auto &piece : ContextPieceChords(signature))
        {
          if (piece.chain)
          {
            patch.provenance.chain.push_back(piece.P);
          }
        }
        // The ENDS of every context entry are device vertices or face crossings of device
        // edges (an arc entry's chords lie on the fitted circle, the device polyline up to
        // its chord sagitta away: its ends are the test).
        for (const auto &entry : signature.at("Context"))
        {
          const auto P = entry.at("P").get<std::array<double, 4>>();
          for (const auto &q : {Global(P[0], P[1]), Global(P[2], P[3])})
          {
            const double distance = DevicePerimeterDistance(identification, q);
            MFEM_VERIFY(distance <= tolerance,
                        "Context piece of model \""
                            << model.name << "\" placed for feature " << feature.id
                            << " lies " << distance / R
                            << " R from every device edge (A10 extended to the context: a "
                               "mis-keyed library or a placement frame defect)!");
            checked++;
          }
        }
        summary.context_points_checked += checked;
      }
      Emit(std::move(patch), model_index, runtime, feature);
    }
    else
    {
      MFEM_ABORT("No patch construction for the matched feature type \"" << feature.type
                                                                         << "\"!");
    }
    if (straight_like_bend && summary.first_order_missing_turn > missing_turn_before)
    {
      summary.first_order_missing_features++;
    }
  }
  return summary;
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
                         requirements != nullptr, requirements != nullptr);
  if (requirements)
  {
    requirements->SetLibrary(request.library, library, coordinate_scale);
  }
  const bool exhaustive_spatial_closure =
      library.exhaustive_spatial_closure ||
      (requirements &&
       request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::WARN);
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
  surface.classify_components = true;
  // Geometric joint noise rule at the library's matching radius (USER decision 121 (B)).
  surface.joint_noise_sagitta = kJointNoiseSagittaOverRadius * library.matching_radius;
  surface.retain_faces =
      requirements ||
      std::any_of(library.models.begin(), library.models.end(),
                  [](const auto &model) { return model.plan_view_boundary.has_value(); });
  // The identification needs every metal face for the decision-73(3) cross-layer zones.
  surface.retain_global_faces = true;
  // Wall time of the geometry steps before the identification (a chip-scale mesh spends
  // minutes here; the identification prints its own stage lines).
  const auto geometry_started = std::chrono::steady_clock::now();
  auto GeometryStageLine = [&](const std::string &text)
  {
    Mpi::Print(
        "  Metal perimeter {} ({:.2f} s)\n", text,
        std::chrono::duration<double>(std::chrono::steady_clock::now() - geometry_started)
            .count());
  };
  auto geometry = ExtractMetalEdgeGeometry(mesh, iodata.boundaries, surface);
  MFEM_VERIFY(!geometry.Empty(),
              "Fabrication-process response matching found no metal perimeter!");
  GeometryStageLine("extracted: " + std::to_string(geometry.segments.size()) +
                    " segments, " + std::to_string(geometry.global_faces.size()) +
                    " global faces");
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
    ExcludeCoincidentMetalEdgeSegmentIndices(geometry, dielectric.edge_exclude_segments,
                                             dielectric.edge_exclude_segment_tolerance,
                                             mesh.SpaceDimension(), segment_indices);
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
  const auto quadrature_cells = LongitudinalQuadratureCells(quadrature);
  ResponseCorrectionData result;
  result.unmatched_policy = request.unmatched_policy;
  result.translational_domain_correction = request.translational_domain_correction;
  result.trace_coupling = request.trace_coupling;
  result.mortar_oversampling = request.mortar_oversampling;
  result.matching_radius = library.matching_radius;
  int next_model_index = 1;
  int next_interpolation_group = 1;
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
      std::vector<bool> ambiguous_side;
      const auto process_normals = BuildMetalEdgeProcessNormals(
          mesh, geometry, missing,
          [&](int attribute) { return mat_op.GetLightSpeedMax(attribute); },
          group.process_normal, &ambiguous_side);
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
        // The process side of a sheet with one material on both sides comes from the
        // configured EdgeFrameNormal or is undetermined (decision 74 step 4).
        segment.ambiguous_process_side = ambiguous_side[i] && !group.process_normal;
        // Conductor identity is the edge-connected metal component (geometric, independent
        // of the attribute numbering and of the problem type).
        MFEM_VERIFY(source.metal_component >= 0,
                    "Unable to determine the connected metal component of an "
                    "automatically detected edge!");
        segment.conductor = source.metal_component;
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
          segment.boundary_condition = boundary_condition;
        }
        else
        {
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
      if (selection.model_to_site[model_edge] == std::numeric_limits<std::size_t>::max())
      {
        continue;
      }
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
                                  int process_axis = 1,
                                  const std::vector<Point3D> *support_points = nullptr)
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
    if (support_points)
    {
      for (const auto &point : *support_points)
      {
        for (int d = 0; d < 2; d++)
        {
          lower[d] = std::min(lower[d], point[plan_axes[d]]);
          upper[d] = std::max(upper[d], point[plan_axes[d]]);
        }
      }
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
    // The facets arrive in rank order; order them by geometry so that every consumer sees
    // the same sequence regardless of the partition.
    std::sort(result.facets.begin(), result.facets.end(),
              [](const PlanViewFacet &first, const PlanViewFacet &second)
              {
                return std::tie(first.conductor, first.points) <
                       std::tie(second.conductor, second.points);
              });
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
        if (selection.model_to_site[model_edge] == std::numeric_limits<std::size_t>::max())
        {
          continue;
        }
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
      const auto plan_view =
          GatherPlanViewFacets(sites, selection.origin, selection.axes,
                               conductor_by_metal_component, 1, &model.support_points);
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
  for (const auto &[signature, group] : groups_by_targets)
  {
    (void)signature;
    auto segments = BuildSegments(group, group.segments);
    global_segments.insert(global_segments.end(), std::make_move_iterator(segments.begin()),
                           std::make_move_iterator(segments.end()));
  }
  GeometryStageLine("framed: " + std::to_string(global_segments.size()) +
                    " targeted segments in " + std::to_string(groups_by_targets.size()) +
                    " interface groups");

  // Boundary-condition labels are diagnostics only: geometrically connected metal carrying
  // different Terminal / PrescribedPotential indices is reported, never split.
  {
    std::map<int, std::set<int>> labels_by_component;
    for (const auto &segment : global_segments)
    {
      for (const int attribute : geometry.segments[segment.geometry_index].metal_attributes)
      {
        if (auto label = GetConductor(iodata.boundaries, attribute))
        {
          labels_by_component[segment.metal_component].insert(*label);
        }
      }
    }
    for (const auto &[component, labels] : labels_by_component)
    {
      if (labels.size() > 1)
      {
        std::string text;
        for (const int label : labels)
        {
          text += (text.empty() ? "" : ", ") + std::to_string(label);
        }
        Mpi::Warning("Edge-connected metal component {:d} carries the distinct conductor "
                     "labels {{{}}} (Ground = 0, Terminal / PrescribedPotential index); "
                     "the identification treats it as one conductor.\n",
                     component, text);
      }
    }
  }

  const AutomaticResponseRequirements law_describer(iodata.units,
                                                    iodata.InputsNondimensionalized());
  const bool frame_normal_configured =
      std::any_of(iodata.boundaries.postpro.dielectric.begin(),
                  iodata.boundaries.postpro.dielectric.end(), [](const auto &entry)
                  { return entry.second.edge_frame_normal.has_value(); });
  std::map<int, FeatureCurvatureMatch> curved_matches;
  std::map<int, FeatureCornerMatch> corner_matches;
  const auto identification = RunGeometryIdentification(
      mesh.GetComm(), geometry, global_segments, library,
      requirements ? *requirements : law_describer, requirements, frame_normal_configured,
      request.span_cap_allowances, &curved_matches, &corner_matches);
  GeometryStageLine("identified and matched: " +
                    std::to_string(identification.features.size()) + " features");
  if (request.patch_construction == ResponseCorrectionData::PatchConstruction::FEATURES)
  {
    // Features-driven construction (default): the identification's feature list is the
    // contract; the legacy per-group classification below is kept behind
    // PatchConstruction = "Legacy" for comparison only. The preflight builds the same
    // patches (the patch dry run, written as surface-response-patches.csv).
    if (diagnostics)
    {
      for (const auto &segment : global_segments)
      {
        diagnostics->selected_length += segment.length;
        for (const auto &[type, target] : segment.targets)
        {
          (void)type;
          diagnostics->selected_length_by_interface[target] += segment.length;
        }
      }
    }
    const auto summary =
        BuildFeaturePatches(library, identification, global_segments, quadrature,
                            requirements ? *requirements : law_describer, diagnostics,
                            result, curved_matches, corner_matches);
    GeometryStageLine("patches built: " + std::to_string(result.patches.size()));
    // Legacy-contract aliases resolved by the matching pass (USER decision 283): one record
    // per alias with the features it served, carried into the operator record and the
    // geometry cache.
    result.legacy_contract.clear();
    for (const auto &feature : identification.features)
    {
      if (feature.legacy_contract.is_null())
      {
        continue;
      }
      const std::string key = feature.legacy_contract.at("Key").get<std::string>();
      auto it = std::find_if(result.legacy_contract.begin(), result.legacy_contract.end(),
                             [&](const auto &entry) { return entry.key == key; });
      if (it == result.legacy_contract.end())
      {
        result.legacy_contract.push_back(
            {feature.legacy_contract.at("Model").get<std::string>(),
             key,
             feature.legacy_contract.at("ContextDigest").get<std::string>(),
             feature.legacy_contract.at("Reason").get<std::string>(),
             {}});
        it = std::prev(result.legacy_contract.end());
      }
      it->features.push_back(feature.id);
    }
    // Quantum near-matches resolved by the matching pass (block (b) DESIGN section 4): one
    // record per feature key with the features it covered, carried into the operator
    // record and the geometry cache.
    result.quantum_near_match.clear();
    for (const auto &feature : identification.features)
    {
      if (feature.quantum_near_match.is_null())
      {
        continue;
      }
      const std::string feature_key =
          feature.quantum_near_match.at("FeatureKey").get<std::string>();
      auto it =
          std::find_if(result.quantum_near_match.begin(), result.quantum_near_match.end(),
                       [&](const auto &entry) { return entry.feature_key == feature_key; });
      if (it == result.quantum_near_match.end())
      {
        result.quantum_near_match.push_back(
            {feature.matched_model.value_or(std::string{}),
             feature.quantum_near_match.at("ModelKey").get<std::string>(),
             feature_key,
             feature.quantum_near_match.at("MaxDeltaQuanta").get<double>(),
             feature.quantum_near_match.at("DifferingNumbers").at("Count").get<int>(),
             {}});
        it = std::prev(result.quantum_near_match.end());
      }
      it->features.push_back(feature.id);
    }
    if (!result.quantum_near_match.empty())
    {
      std::string text;
      for (const auto &entry : result.quantum_near_match)
      {
        text += fmt::format("{}{} <- key {} ({:g} quanta, {:d} feature{})",
                            text.empty() ? "" : ", ", entry.model,
                            entry.feature_key.substr(0, 12), entry.max_delta_quanta,
                            entry.features.size(), entry.features.size() == 1 ? "" : "s");
      }
      Mpi::Print(mesh.GetComm(), " Quantum near-matches (block (b) DESIGN section 4): {}\n",
                 text);
    }
    if (!result.legacy_contract.empty())
    {
      std::string text;
      for (const auto &entry : result.legacy_contract)
      {
        text += fmt::format("{}{} <- key {} ({:d} feature{})", text.empty() ? "" : ", ",
                            entry.model, entry.key.substr(0, 12), entry.features.size(),
                            entry.features.size() == 1 ? "" : "s");
      }
      Mpi::Print(mesh.GetComm(), " Legacy-contract aliases (USER decision 283): {}\n",
                 text);
    }
    std::string unmatched;
    for (const auto &[type, entry] : summary.unmatched_by_type)
    {
      unmatched += fmt::format("{}{}: {:d} ({:.6e})", unmatched.empty() ? "" : ", ", type,
                               entry.first, entry.second * coordinate_scale);
    }
    std::string patches;
    for (const auto &[type, count] : summary.patches_by_type)
    {
      patches += fmt::format("{}{}: {:d}", patches.empty() ? "" : ", ", type, count);
    }
    Mpi::Print("\nAutomatic fabrication-process response matching (Features):\n"
               " Library: {}\n"
               " Matched features: {:d} ({:.6e})\n"
               " Unmatched features (omitted alone): {:d} ({:.6e}) {{{}}}\n"
               " Patches: {:d} {{{}}}\n"
               " Runtime models: {:d}\n"
               " Curvature: {:d} feature(s) on a curvature family, {:d} straight-like "
               "feature(s) with a first-order term\n"
               " Corners: {:d} feature(s) on an angle-interpolated corner family\n"
               " Context (A10 extended): {:d} placed context piece end(s) verified on "
               "device edges\n",
               library.name, summary.matched_features,
               summary.matched_length * coordinate_scale, summary.unmatched_features,
               summary.unmatched_length * coordinate_scale, unmatched,
               static_cast<int>(result.patches.size()), patches,
               static_cast<int>(result.models.size()), summary.curved_family_features,
               summary.first_order_features, summary.corner_family_features,
               static_cast<int>(summary.context_points_checked));
    if (summary.first_order_missing_features > 0)
    {
      std::string nodes;
      for (const auto &node : summary.first_order_missing_nodes)
      {
        nodes += (nodes.empty() ? "" : ", ") + node;
      }
      Mpi::Warning("Fabrication-process response library \"{}\" has no first-order "
                   "curvature node (kappa = 1 / StraightBendRadiusOverR) for {} "
                   "straight-like feature(s) with bends ({:.6e} rad of turn left at the "
                   "straight model; missing nodes: {}).\n",
                   library.name, summary.first_order_missing_features,
                   summary.first_order_missing_turn, nodes);
    }
    if (!summary.inconsistent_features.empty())
    {
      Mpi::Warning(
          "{} matched pair feature(s) ({:.6e} length units) whose sides do not face "
          "each other at their separation are omitted (identification defect, "
          "recorded; ids {} ...).\n",
          static_cast<int>(summary.inconsistent_features.size()),
          summary.inconsistent_length * coordinate_scale,
          summary.inconsistent_features.front());
    }
    if (summary.unmatched_features > 0 || !summary.inconsistent_features.empty())
    {
      if (summary.unmatched_features > 0)
      {
        Mpi::Warning("Fabrication-process response library \"{}\" has no model for {} "
                     "identified feature(s) ({:.6e} length units); correction is disabled "
                     "for these features only.\n",
                     library.name, summary.unmatched_features,
                     summary.unmatched_length * coordinate_scale);
      }
      if (!requirements &&
          request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
      {
        MFEM_ABORT("Automatic fabrication-process response matching failed: "
                   << summary.unmatched_features
                   << " identified feature(s) have no library model and "
                   << summary.inconsistent_features.size()
                   << " pair feature(s) are inconsistent (UnmatchedPolicy = Error)!");
      }
    }
    MFEM_VERIFY(requirements || (!result.models.empty() && !result.patches.empty()),
                "Fabrication-process response matching produced no usable correction "
                "patches!");
    return result;
  }
  // Legacy per-group classification (comparison only): it reads fillets as runs of REGULAR
  // sub-corner joints and pairs / neighbourhoods on physical chains broken at corners of
  // the former 30 deg class, so it re-reads the vertex classes and chains of the same
  // segments at that class (kLegacyCornerTurnToleranceDegrees; the segments, vertices,
  // faces and their numbering are those of the extraction above: only the vertex types and
  // chain ids differ; the face pointers of the plan-view mask index are rebuilt). The
  // identification and its manifest above use the joint noise threshold (USER decision
  // 117(4)).
  {
    MetalSurfaceExtraction legacy_surface = surface;
    legacy_surface.joint_noise_sagitta = 0.0;
    legacy_surface.corner_turn_tolerance_degrees = kLegacyCornerTurnToleranceDegrees;
    geometry = ExtractMetalEdgeGeometry(mesh, iodata.boundaries, legacy_surface);
    surface_faces_by_component.clear();
    for (const auto &face : geometry.surface_faces)
    {
      surface_faces_by_component[face.component].push_back(&face);
    }
  }
  std::set<std::size_t> identification_excluded_segments;
  for (std::size_t i = 0; i < identification.segments.size(); i++)
  {
    if (identification.segments[i].exclusion ||
        !identification.segments[i].excluded_portions.empty())
    {
      identification_excluded_segments.insert(i);
    }
  }

  const double global_interaction_distance = 2.0 * library.matching_radius;
  const DecisionQuantizer global_quantizer(library.matching_radius);
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
    if ((!exhaustive_spatial_closure &&
         global_segments[i].targets == global_segments[j].targets) ||
        first_source.physical_chain == second_source.physical_chain ||
        SegmentsShareVertex(first_source, second_source))
    {
      continue;
    }
    const auto closest =
        ClosestSegmentApproach(global_segments[i].p0, global_segments[i].p1,
                               global_segments[j].p0, global_segments[j].p1);
    if (!global_quantizer.LengthSquaredLess(closest.distance_squared,
                                            global_interaction_distance))
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
    for (std::size_t candidate = 0; candidate < global_spatial_events.size(); candidate++)
    {
      if (!visited_global_event[candidate] &&
          std::all_of(event_component.begin(), event_component.end(),
                      [&](std::size_t member)
                      {
                        return global_quantizer.LengthLess(
                            Distance(global_spatial_events[member].center,
                                     global_spatial_events[candidate].center),
                            (exhaustive_spatial_closure ? 4.0 : 1.0) *
                                global_interaction_distance);
                      }))
      {
        visited_global_event[candidate] = true;
        event_component.push_back(candidate);
      }
    }

    double event_diameter = 0.0;
    std::set<TargetSignature> component_targets;
    for (std::size_t i = 0; i < event_component.size(); i++)
    {
      const auto &event = global_spatial_events[event_component[i]];
      component_targets.emplace(global_segments[event.first].targets.begin(),
                                global_segments[event.first].targets.end());
      component_targets.emplace(global_segments[event.second].targets.begin(),
                                global_segments[event.second].targets.end());
      for (std::size_t j = i + 1; j < event_component.size(); j++)
      {
        event_diameter = std::max(
            event_diameter, Distance(global_spatial_events[event_component[i]].center,
                                     global_spatial_events[event_component[j]].center));
      }
    }
    // Build cross-interface components from every nearby interaction, including
    // same-interface pairs connected to the cross-interface neighborhood. Otherwise a
    // later per-interface model can overlap the cross-interface model and expose a new
    // unified coupon only after expensive generation. Purely local components remain the
    // responsibility of the per-interface pass below.
    if (component_targets.size() < 2)
    {
      continue;
    }
    const double maximum_global_event_diameter = 4.0 * global_interaction_distance;
    MFEM_VERIFY(
        global_quantizer.LengthAtMost(event_diameter, maximum_global_event_diameter),
        "Cross-interface spatial-event clustering produced an oversized "
        "component. Split the component or increase the matching radius!");

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

  struct ClaimedSpatialSupport
  {
    ElementBox box;
    ElementBox local_box;
    std::set<int> targets;
    std::string model;
    std::vector<SpatialEdgeSite3D> sites;
    Point3D origin{};
    std::array<Point3D, 3> axes{};
  };
  auto OwnershipSites =
      [](const SpatialClusterSelection3D &selection, const LibraryModel &model)
  {
    std::vector<SpatialEdgeSite3D> sites;
    sites.reserve(selection.model_to_site.size());
    for (std::size_t model_edge = 0; model_edge < selection.model_to_site.size();
         model_edge++)
    {
      if (selection.model_to_site[model_edge] == std::numeric_limits<std::size_t>::max())
      {
        continue;
      }
      auto site = selection.sites[selection.model_to_site[model_edge]];
      site.interval = model.spatial_edges[model_edge].interval;
      sites.push_back(std::move(site));
    }
    return sites;
  };
  auto MakeClaimedSupport =
      [&](const SpatialClusterSelection3D &selection, const LibraryModel &model)
  {
    ElementBox local_box;
    for (const auto &point : model.support_points)
    {
      for (int d = 0; d < 3; d++)
      {
        local_box.min[d] = std::min(local_box.min[d], point[d]);
        local_box.max[d] = std::max(local_box.max[d], point[d]);
      }
    }
    return ClaimedSpatialSupport{SpatialSupportBox(selection, model),
                                 local_box,
                                 SpatialTargetAttributes(selection),
                                 model.name,
                                 OwnershipSites(selection, model),
                                 selection.origin,
                                 selection.axes};
  };
  auto OwnsSpatialSupport =
      [&](const ClaimedSpatialSupport &claimed, const SpatialClusterSelection3D &candidate,
          const LibraryModel &candidate_model, const std::set<int> &candidate_targets)
  {
    if (!std::includes(claimed.targets.begin(), claimed.targets.end(),
                       candidate_targets.begin(), candidate_targets.end()))
    {
      return false;
    }
    const double tolerance = 1.0e-10 * library.matching_radius;
    const auto candidate_sites = OwnershipSites(candidate, candidate_model);
    for (const auto &site : candidate_sites)
    {
      const auto owner =
          std::find_if(claimed.sites.begin(), claimed.sites.end(), [&](const auto &entry)
                       { return entry.physical_chain == site.physical_chain; });
      if (owner == claimed.sites.end())
      {
        return false;
      }
      const Point3D owner_tangent =
          Normalize(Cross(owner->gap_direction, owner->process_normal));
      const Point3D site_tangent =
          Normalize(Cross(site.gap_direction, site.process_normal));
      for (const double coordinate : site.interval)
      {
        const Point3D endpoint = Add(site.point, Scale(coordinate, site_tangent));
        const double owner_coordinate =
            Dot(Subtract(endpoint, owner->point), owner_tangent);
        if (owner_coordinate < owner->interval[0] - tolerance ||
            owner_coordinate > owner->interval[1] + tolerance)
        {
          return false;
        }
      }
    }
    for (const auto &point : candidate_model.support_points)
    {
      const Point3D global = TransformLocalPoint(candidate.origin, candidate.axes, point);
      const Point3D relative = Subtract(global, claimed.origin);
      Point3D local{};
      for (int d = 0; d < 3; d++)
      {
        local[d] = Dot(relative, claimed.axes[d]);
      }
      if (!claimed.local_box.Contains(local, 3, tolerance))
      {
        return false;
      }
    }
    return true;
  };
  std::vector<ClaimedSpatialSupport> claimed_spatial_supports;
  for (const auto &selection : cross_interface_selections)
  {
    const auto &source = library.models[selection.response.models.front().index];
    if (!source.support_points.empty())
    {
      claimed_spatial_supports.push_back(MakeClaimedSupport(selection, source));
    }
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
                                      global_quantizer.LengthLess(
                                          Distance(interaction.center, center),
                                          global_interaction_distance);
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
    const double interaction_distance = 2.0 * group.matching_radius;
    const DecisionQuantizer quantizer(group.matching_radius);
    const std::set<std::size_t> group_segment_indices(group.segments.begin(),
                                                      group.segments.end());
    std::vector<std::size_t> usable_segment_indices;
    usable_segment_indices.reserve(group.segments.size());
    double total_selected_length = 0.0;
    int externally_conflicted_segments = 0;
    int identification_excluded = 0;
    for (const std::size_t geometry_index : group.segments)
    {
      const auto &source = geometry.segments[geometry_index];
      const auto &p0 = geometry.vertices[source.vertices[0]].coordinate;
      const auto &p1 = geometry.vertices[source.vertices[1]].coordinate;
      total_selected_length += Distance(p0, p1);
      if (identification_excluded_segments.find(geometry_index) !=
          identification_excluded_segments.end())
      {
        // Excluded by the identification (non-planar, undetermined process side, within a
        // cross-layer zone): recorded in the manifest, never corrected by the legacy path.
        identification_excluded++;
        continue;
      }
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
        if (quantizer.LengthSquaredLess(distance_squared, interaction_distance) &&
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
    if (identification_excluded > 0)
    {
      Mpi::Print("Omitting {:d} of {:d} three-dimensional target edge segments excluded by "
                 "the geometry identification (see the Exclusions of the manifest).\n",
                 identification_excluded, static_cast<int>(group.segments.size()));
    }
    if (externally_conflicted_segments > 0)
    {
      if (!requirements &&
          request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
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
                   quantizer.LengthSquaredLess(
                       PointSegmentDistanceSquared(neighborhood.point, segments[first]),
                       interaction_distance) &&
                   quantizer.LengthSquaredLess(
                       PointSegmentDistanceSquared(neighborhood.point, segments[second]),
                       interaction_distance);
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
      if (!quantizer.LengthSquaredLess(closest.distance_squared, interaction_distance))
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
    std::set<std::pair<int, int>> owned_spatial_pairs;
    std::vector<bool> visited_spatial_event(spatial_events.size(), false);
    for (std::size_t seed = 0; seed < spatial_events.size(); seed++)
    {
      if (visited_spatial_event[seed])
      {
        continue;
      }
      std::vector<std::size_t> event_component = {seed};
      visited_spatial_event[seed] = true;
      for (std::size_t candidate = 0; candidate < spatial_events.size(); candidate++)
      {
        if (!visited_spatial_event[candidate] &&
            std::all_of(event_component.begin(), event_component.end(),
                        [&](std::size_t member)
                        {
                          return quantizer.LengthLess(
                              Distance(spatial_events[member].center,
                                       spatial_events[candidate].center),
                              (exhaustive_spatial_closure ? 4.0 : 1.0) *
                                  interaction_distance);
                        }))
        {
          visited_spatial_event[candidate] = true;
          event_component.push_back(candidate);
        }
      }

      double event_diameter = 0.0;
      for (std::size_t i = 0; i < event_component.size(); i++)
      {
        for (std::size_t j = i + 1; j < event_component.size(); j++)
        {
          event_diameter =
              std::max(event_diameter, Distance(spatial_events[event_component[i]].center,
                                                spatial_events[event_component[j]].center));
        }
      }
      const double maximum_event_diameter = 4.0 * interaction_distance;
      MFEM_VERIFY(quantizer.LengthAtMost(event_diameter, maximum_event_diameter),
                  "Spatial-event clustering produced an oversized component. Split the "
                  "component or increase the matching radius!");

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
      bool support_owned = false;
      std::vector<SpatialEdgeSite3D> unified_sites;
      std::vector<std::string> overlapping_models;
      std::string candidate_model;
      if (selection)
      {
        const auto &source = library.models[selection->response.models.front().index];
        candidate_model = source.name;
        if (!source.support_points.empty())
        {
          const auto candidate_box = SpatialSupportBox(*selection, source);
          const auto candidate_targets = SpatialTargetAttributes(*selection);
          for (const auto &claimed : claimed_spatial_supports)
          {
            if (!candidate_box.InteriorOverlaps(claimed.box, 3))
            {
              continue;
            }
            if (OwnsSpatialSupport(claimed, *selection, source, candidate_targets))
            {
              support_owned = true;
              break;
            }
            if (!requirements)
            {
              MFEM_ABORT("Spatial response matching volumes for models \""
                         << claimed.model << "\" and \"" << source.name
                         << "\" overlap without one model owning the complete support "
                            "and interface mapping. Generate one unified spatial "
                            "cluster!");
            }
            overlapping_models.push_back(claimed.model);
            auto MergeSites = [&](const std::vector<SpatialEdgeSite3D> &additional)
            {
              for (const auto &site : additional)
              {
                auto existing = std::find_if(
                    unified_sites.begin(), unified_sites.end(), [&](const auto &candidate)
                    { return candidate.physical_chain == site.physical_chain; });
                if (existing == unified_sites.end())
                {
                  unified_sites.push_back(site);
                  continue;
                }
                const Point3D tangent =
                    Normalize(Cross(existing->gap_direction, existing->process_normal));
                const Point3D site_tangent =
                    Normalize(Cross(site.gap_direction, site.process_normal));
                const double offset = Dot(Subtract(site.point, existing->point), tangent);
                const double orientation = Dot(site_tangent, tangent);
                const double begin =
                    offset + (orientation >= 0.0 ? site.interval[0] : -site.interval[1]);
                const double end =
                    offset + (orientation >= 0.0 ? site.interval[1] : -site.interval[0]);
                existing->interval[0] = std::min(existing->interval[0], begin);
                existing->interval[1] = std::max(existing->interval[1], end);
              }
            };
            if (unified_sites.empty())
            {
              if (claimed.sites.size() >= selection->sites.size())
              {
                unified_sites = claimed.sites;
                MergeSites(OwnershipSites(*selection, source));
              }
              else
              {
                unified_sites = OwnershipSites(*selection, source);
                MergeSites(claimed.sites);
              }
            }
            else
            {
              MergeSites(claimed.sites);
            }
          }
        }
      }
      if (!support_owned && !unified_sites.empty())
      {
        auto [spatial_geometry, targets_by_slot] =
            DescribeMissingSpatialGeometry(unified_sites);
        std::string reason = "Overlapping spatial matching volumes for models ";
        for (std::size_t i = 0; i < overlapping_models.size(); i++)
        {
          reason += (i == 0 ? "\"" : ", \"") + overlapping_models[i] + "\"";
        }
        reason += " and \"" + candidate_model + "\" require one complete support owner";
        requirements->Add(3, LibraryTopology::SPATIAL_EDGE_CLUSTER, targets_by_slot,
                          unified_sites.front().boundary_condition, spatial_geometry,
                          library, nullptr,
                          2.0 * group.matching_radius * unified_sites.size(), reason);
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
      if (support_owned)
      {
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
          owned_spatial_pairs.emplace(first_chain, second_chain);
        }
        continue;
      }
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
      const auto &source = library.models[selection->response.models.front().index];
      if (!source.support_points.empty())
      {
        claimed_spatial_supports.push_back(MakeClaimedSupport(*selection, source));
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
      if (exhaustive_spatial_closure &&
          owned_spatial_pairs.find({first_chain, second_chain}) !=
              owned_spatial_pairs.end())
      {
        return true;
      }
      const auto closest = ClosestSegmentApproach(segments[first].p0, segments[first].p1,
                                                  segments[second].p0, segments[second].p1);
      const Point3D center = Scale(
          0.5,
          Add(Interpolate(segments[first], closest.first * segments[first].length),
              Interpolate(segments[second], closest.second * segments[second].length)));
      return std::any_of(
          interactions.begin(), interactions.end(),
          [&](const auto &interaction)
          {
            return interaction.first_chain == first_chain &&
                   interaction.second_chain == second_chain &&
                   (exhaustive_spatial_closure ||
                    quantizer.LengthLess(Distance(interaction.center, center),
                                         interaction_distance));
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
      if (!requirements &&
          request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
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
          !quantizer.LengthSquaredLess(
              SegmentDistanceSquared(segments[i].p0, segments[i].p1, segments[j].p0,
                                     segments[j].p1),
              interaction_distance))
      {
        continue;
      }
      if (ConnectedNearVertex(i, j) ||
          IsSpatiallyMatched(i, j, candidate_spatial_interactions))
      {
        continue;
      }

      const double tangent_dot = Dot(segments[i].tangent, segments[j].tangent);
      if (DecisionQuantizer::DirectionLess(std::abs(tangent_dot), 1.0 - 1.0e-8))
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
      VerifyParallelOverlap(segments[i], segments[j], tangent_dot, second_begin, second_end,
                            tolerance, interaction_distance);
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
      if (selection.model_to_site[model_edge_index] ==
          std::numeric_limits<std::size_t>::max())
      {
        return;
      }
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
        if (model_selection->IsInterpolated())
        {
          patch.interpolation_group = next_interpolation_group++;
        }
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
    // Segments in a non-parallel neighborhood that no matched cluster covers lose their own
    // correction only (decision 74(1): an unmatched feature disables itself, never the
    // interface group).
    std::set<std::size_t> nonparallel_omitted_segments;
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
          !quantizer.LengthSquaredLess(
              SegmentDistanceSquared(segments[i].p0, segments[i].p1, segments[j].p0,
                                     segments[j].p1),
              interaction_distance))
      {
        continue;
      }
      if (ConnectedNearVertex(i, j) ||
          IsSpatiallyMatched(i, j, active_spatial_interactions))
      {
        continue;
      }

      const double tangent_dot = Dot(segments[i].tangent, segments[j].tangent);
      if (DecisionQuantizer::DirectionLess(std::abs(tangent_dot), 1.0 - 1.0e-8))
      {
        nonparallel_omitted_segments.insert(i);
        nonparallel_omitted_segments.insert(j);
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
      VerifyParallelOverlap(segments[i], segments[j], tangent_dot, second_begin, second_end,
                            tolerance, interaction_distance);
      pairs.push_back({i, j, first_begin, first_end, std::max(0.0, second_begin),
                       std::min(segments[j].length, second_end)});
      paired_intervals[i].emplace_back(first_begin, first_end);
      paired_intervals[j].emplace_back(std::max(0.0, second_begin),
                                       std::min(segments[j].length, second_end));
    }

    if (!nonparallel_omitted_segments.empty())
    {
      // An omitted segment takes no part in any pair either: its whole correction is off.
      pairs.erase(std::remove_if(pairs.begin(), pairs.end(),
                                 [&](const EdgePair3D &pair)
                                 {
                                   return nonparallel_omitted_segments.count(pair.first) ||
                                          nonparallel_omitted_segments.count(pair.second);
                                 }),
                  pairs.end());
      for (auto &intervals : paired_intervals)
      {
        intervals.clear();
      }
      for (const auto &pair : pairs)
      {
        paired_intervals[pair.first].emplace_back(pair.first_begin, pair.first_end);
        paired_intervals[pair.second].emplace_back(pair.second_begin, pair.second_end);
      }
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
    // Spans, pair intervals and isolated intervals without a library model lose their own
    // correction only (decision 74(1): an unmatched feature disables itself, never the
    // interface group); counted for the summary warning.
    int omitted_parallel_cluster_spans = 0, omitted_pair_intervals = 0,
        omitted_isolated_intervals = 0;
    if (!unmatched_parallel_clusters.empty())
    {
      Mpi::Warning(
          "Fabrication-process response library \"{}\" has no matching "
          "ParallelEdgeCluster model for {} three-dimensional longitudinal span(s); "
          "correction is disabled for these spans only.\n",
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
      omitted_parallel_cluster_spans +=
          static_cast<int>(unmatched_parallel_clusters.size());
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
          patch.longitudinal_cell = LongitudinalCellOffsets(
              quadrature_cells[q], interval_begin, interval_end, first_distance,
              first.tangent, Normalize(Cross(patch.axis_u, patch.axis_v)));
          if (model_selection.IsInterpolated())
          {
            patch.interpolation_group = next_interpolation_group++;
          }
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
          patch.longitudinal_cell = LongitudinalCellOffsets(
              quadrature_cells[q], interval_begin, interval_end, first_distance,
              first.tangent, Normalize(Cross(patch.axis_u, patch.axis_v)));
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
      // A longitudinal span claimed by a parallel-edge cluster (matched or unmatched: an
      // unmatched cluster span is omitted, never re-described by its pairs) is not a pair.
      auto SubtractSpan = [&](const std::vector<std::size_t> &cluster,
                              const Point3D &tangent, double begin, double end)
      {
        if (std::find(cluster.begin(), cluster.end(), pair.first) == cluster.end() ||
            std::find(cluster.begin(), cluster.end(), pair.second) == cluster.end())
        {
          return;
        }
        const double orientation = Dot(first.tangent, tangent);
        MFEM_ASSERT(std::abs(orientation) > 1.0 - 1.0e-8,
                    "A parallel-edge cluster contains incompatible tangents!");
        double excluded_begin = (begin - Dot(first.p0, tangent)) / orientation;
        double excluded_end = (end - Dot(first.p0, tangent)) / orientation;
        if (excluded_begin > excluded_end)
        {
          std::swap(excluded_begin, excluded_end);
        }
        SubtractInterval(pair_intervals, excluded_begin, excluded_end);
      };
      for (const auto &span : parallel_cluster_spans)
      {
        SubtractSpan(span.selection.ordered_edges, span.tangent, span.begin, span.end);
      }
      for (const auto &span : unmatched_parallel_clusters)
      {
        SubtractSpan(span.edges, span.tangent, span.begin, span.end);
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
              "this pair interval only.\n");
          omitted_pair_intervals++;
          continue;
        }
        if (Dot(first.axis_v, second.axis_v) <= 0.95)
        {
          Mpi::Warning(
              "Nearby three-dimensional edges have incompatible process normals; a "
              "dedicated cross-layer coupon is required and correction is disabled for "
              "this pair interval only.\n");
          omitted_pair_intervals++;
          continue;
        }
        const Point3D process_normal = Normalize(Add(first.axis_v, second.axis_v));
        if (std::abs(Dot(direction, process_normal)) > 1.0e-8)
        {
          Mpi::Warning(
              "Nearby three-dimensional edges are offset along the process normal; a "
              "dedicated cross-layer coupon is required and correction is disabled for "
              "this pair interval only.\n");
          omitted_pair_intervals++;
          continue;
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
              "edges; correction is disabled for this pair interval only.\n");
          omitted_pair_intervals++;
          continue;
        }
        const double separation = Distance(first_mid, second_mid);
        const auto model_selection =
            FindLibraryModel(library, topology, separation, first.boundary_condition);
        if (!model_selection)
        {
          Mpi::Warning(
              "Fabrication-process response library \"{}\" has no {} model at separation "
              "{:.6e} mesh units; correction is disabled for this pair interval only.\n",
              library.name, TopologyName(topology), separation * coordinate_scale);
          if (requirements)
          {
            requirements->Add(
                3, topology, group.targets, first.boundary_condition,
                {{"EdgeCount", 2}, {"Separation", requirements->ScaleLength(separation)}},
                library, nullptr, pair_end - pair_begin,
                "No compatible paired-edge model or interpolation bracket");
          }
          omitted_pair_intervals++;
          continue;
        }
        group_maximum_library_distance =
            std::max(group_maximum_library_distance, model_selection->normalized_distance);
        AppendQuadrature(*model_selection, pair.first, pair_begin, pair_end, pair.second);
        group_matched_intervals++;
        group_interpolated_paired_intervals += model_selection->IsInterpolated();
      }
    }

    if (!nonparallel_omitted_segments.empty())
    {
      Mpi::Warning(
          "Omitting {} of {} three-dimensional target edge segments in non-parallel "
          "neighborhoods without a matched cluster model (correction disabled for "
          "these segments only).\n",
          static_cast<int>(nonparallel_omitted_segments.size()),
          static_cast<int>(segments.size()));
    }
    for (std::size_t i = 0; group_matched && i < segments.size(); i++)
    {
      if (nonparallel_omitted_segments.find(i) != nonparallel_omitted_segments.end())
      {
        continue;
      }
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
        omitted_isolated_intervals += static_cast<int>(isolated_intervals.size());
        continue;
      }
      for (const auto &[isolated_begin, isolated_end] : isolated_intervals)
      {
        AppendQuadrature(*isolated_model, i, isolated_begin, isolated_end, std::nullopt);
        group_matched_intervals++;
      }
    }
    if (omitted_isolated_intervals > 0)
    {
      Mpi::Warning(
          "Fabrication-process response library \"{}\" has no isolated-edge model for {} "
          "unmatched longitudinal span(s) using the selected metal boundary condition; "
          "correction is disabled for these spans only.\n",
          library.name, omitted_isolated_intervals);
    }
    if (omitted_parallel_cluster_spans + omitted_pair_intervals +
            omitted_isolated_intervals >
        0)
    {
      Mpi::Warning(
          "Omitting {} parallel-edge cluster span(s), {} paired-edge interval(s) and "
          "{} isolated interval(s) without a library model in this interface group "
          "(correction disabled for these features only).\n",
          omitted_parallel_cluster_spans, omitted_pair_intervals,
          omitted_isolated_intervals);
    }

    if (!group_matched)
    {
      unmatched_groups++;
      if (!requirements &&
          request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR)
      {
        MFEM_ABORT("Automatic fabrication-process response matching failed!");
      }
      continue;
    }
    if (!requirements &&
        request.unmatched_policy == ResponseCorrectionData::UnmatchedPolicy::ERROR &&
        (omitted_parallel_cluster_spans + omitted_pair_intervals +
             omitted_isolated_intervals +
             static_cast<int>(nonparallel_omitted_segments.size()) >
         0))
    {
      MFEM_ABORT("Automatic fabrication-process response matching failed: "
                 << omitted_parallel_cluster_spans << " parallel-edge cluster span(s), "
                 << omitted_pair_intervals << " paired-edge interval(s), "
                 << omitted_isolated_intervals << " isolated interval(s) and "
                 << nonparallel_omitted_segments.size()
                 << " non-parallel segment(s) have no library model (UnmatchedPolicy = "
                    "Error)!");
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

TraceMeshData ReadTraceMesh(const std::string &vertex_path,
                            const std::string &triangle_path)
{
  TraceMeshData mesh;
  std::ifstream vertices(vertex_path);
  MFEM_VERIFY(vertices,
              "Unable to open response trace vertex file \"" << vertex_path << "\"!");
  std::string line;
  while (std::getline(vertices, line))
  {
    auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#')
    {
      continue;
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    int index = 0;
    TraceMeshData::Vertex vertex;
    if (!(row >> index >> vertex.point[0] >> vertex.point[1] >> vertex.point[2] >>
          vertex.basis >> vertex.conductor))
    {
      MFEM_VERIFY(mesh.vertices.empty(),
                  "Could not parse response trace vertex file \"" << vertex_path << "\"!");
      continue;
    }
    // Optional slave-vertex columns (parent_a, parent_b, weight_a); a knot row carries
    // zeros there or nothing.
    if (!(row >> vertex.parent_a >> vertex.parent_b >> vertex.weight_a))
    {
      vertex.parent_a = vertex.parent_b = 0;
      vertex.weight_a = 0.0;
    }
    MFEM_VERIFY(index == static_cast<int>(mesh.vertices.size()) + 1 && vertex.basis >= 0 &&
                    vertex.conductor >= 0 &&
                    std::all_of(vertex.point.begin(), vertex.point.end(),
                                [](double value) { return std::isfinite(value); }),
                "Invalid response trace vertex in \"" << vertex_path << "\"!");
    MFEM_VERIFY((vertex.parent_a == 0 && vertex.parent_b == 0) ||
                    (vertex.parent_a > 0 && vertex.parent_b > 0 &&
                     vertex.parent_a != vertex.parent_b && vertex.basis == 0 &&
                     vertex.conductor == 0 && std::isfinite(vertex.weight_a) &&
                     vertex.weight_a >= 0.0 && vertex.weight_a <= 1.0),
                "Invalid slave response trace vertex (parents / weight) in \""
                    << vertex_path << "\"!");
    mesh.vertices.push_back(vertex);
  }
  MFEM_VERIFY(!mesh.vertices.empty(),
              "Response trace vertex file \"" << vertex_path << "\" is empty!");
  for (const auto &vertex : mesh.vertices)
  {
    for (const int parent : {vertex.parent_a, vertex.parent_b})
    {
      MFEM_VERIFY(parent == 0 || (parent <= static_cast<int>(mesh.vertices.size()) &&
                                  mesh.vertices[parent - 1].basis == parent),
                  "A slave response trace vertex must name two basis knots as parents (\""
                      << vertex_path << "\")!");
    }
  }

  std::ifstream triangles(triangle_path);
  MFEM_VERIFY(triangles,
              "Unable to open response trace triangle file \"" << triangle_path << "\"!");
  while (std::getline(triangles, line))
  {
    auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#')
    {
      continue;
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    int index = 0;
    std::array<int, 3> triangle{};
    if (!(row >> index >> triangle[0] >> triangle[1] >> triangle[2]))
    {
      MFEM_VERIFY(mesh.triangles.empty(), "Could not parse response trace triangle file \""
                                              << triangle_path << "\"!");
      continue;
    }
    MFEM_VERIFY(index == static_cast<int>(mesh.triangles.size()) + 1,
                "Invalid response trace triangle index in \"" << triangle_path << "\"!");
    for (int &vertex : triangle)
    {
      MFEM_VERIFY(vertex > 0 && vertex <= static_cast<int>(mesh.vertices.size()),
                  "Invalid response trace triangle vertex in \"" << triangle_path << "\"!");
      vertex--;
    }
    MFEM_VERIFY(triangle[0] != triangle[1] && triangle[1] != triangle[2] &&
                    triangle[2] != triangle[0],
                "Degenerate response trace triangle in \"" << triangle_path << "\"!");
    mesh.triangles.push_back(triangle);
  }
  MFEM_VERIFY(!mesh.triangles.empty(),
              "Response trace triangle file \"" << triangle_path << "\" is empty!");
  return mesh;
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

std::string CheckCornerRuleCouponFiles(const LibraryModel &model, double position_tolerance)
{
  const bool convex = model.topology == LibraryTopology::CONVEX_CORNER;
  const auto &rule = *model.trace_basis;
  const auto points = ReadBasisPoints(model.response.basis_points);
  const CornerBoxRings rings = DescribeCornerBoxRings(points, model.response.contour_groups,
                                                      model.response.zero_trace_indices);
  if (rule.AllRings())
  {
    // The outer ring levels: -R, -R/3, -d, 0, t, t + d, R/3, R with t the upper metal ring
    // and d the trench depth (the ring below z = 0); both metal rings at z = 0 and t.
    std::vector<double> metal_levels, levels;
    for (int r = 0; r < rings.outer_count; r++)
    {
      levels.push_back(rings.rings[r].z);
      if (rings.rings[r].metal)
      {
        metal_levels.push_back(rings.rings[r].z);
      }
    }
    const double radius = rings.radius, tolerance = 1.0e-9 * radius;
    if (metal_levels.size() != 2 || std::abs(metal_levels[0]) > tolerance ||
        !(metal_levels[1] > tolerance))
    {
      return "the metal rings are not at z = 0 and z = MetalThickness";
    }
    const double t = metal_levels[1];
    double d = 0.0;
    for (const double z : levels)
    {
      if (z < -tolerance && (d == 0.0 || z > -d))
      {
        d = -z;
      }
    }
    const std::vector<double> expected = CornerRuleLevels(rule, radius, t, d);
    bool equal = levels.size() == expected.size();
    for (std::size_t i = 0; equal && i < expected.size(); i++)
    {
      equal = std::abs(levels[i] - expected[i]) <= tolerance;
    }
    if (!equal)
    {
      std::ostringstream text;
      text << "the outer ring levels are not the AllRingsFollowMetal set {-R, -R/3, "
              "-OveretchDepth, 0, MetalThickness, MetalThickness + k OveretchDepth (k in "
              "ExtraLevelsAboveOverOveretch), R/3, R} (found";
      for (const double z : levels)
      {
        text << " " << z;
      }
      text << ")";
      return text.str();
    }
  }
  std::optional<double> connectivity_radians;
  if (model.corner_connectivity_angle_degrees)
  {
    connectivity_radians =
        *model.corner_connectivity_angle_degrees * std::acos(-1.0) / 180.0;
  }
  // Every knot's zero flag against ZeroTraceIndices: the rule's PEC knots are the zero
  // slots of the rings that meet the metal and nothing else (a spuriously zeroed free knot,
  // or a free metal-interior knot, is refused here rather than caught by the self-check
  // alone).
  {
    const std::set<int> zero_set(model.response.zero_trace_indices.begin(),
                                 model.response.zero_trace_indices.end());
    const auto zero_slots = CornerZeroSlots(convex, rule);
    for (const auto &ring : rings.rings)
    {
      for (int slot = 0; slot < ring.size; slot++)
      {
        const int k = ring.offset + slot;
        const bool rule_zero = ring.metal && std::find(zero_slots.begin(), zero_slots.end(),
                                                       slot) != zero_slots.end();
        if (rule_zero != (zero_set.count(k) > 0))
        {
          return "basis point " + std::to_string(k + 1) +
                 (rule_zero ? " is a PEC knot of the trace basis rule but is not in "
                              "ZeroTraceIndices"
                            : " is in ZeroTraceIndices but is a free knot of the trace "
                              "basis rule");
        }
      }
    }
  }
  const auto rule_basis = BuildCornerTraceBasis(
      points, model.response.contour_groups, model.response.zero_trace_indices, model.angle,
      convex, rule, connectivity_radians);
  for (std::size_t k = 0; k < points.size(); k++)
  {
    if (Distance(points[k], rule_basis.knots[k]) > position_tolerance)
    {
      return "basis point " + std::to_string(k + 1) +
             " is not at the trace basis rule's position for the coupon's angle";
    }
  }
  if (!(model.corner_connectivity_angle_degrees || rule.AllRings()))
  {
    return "";
  }
  if (model.response.trace_vertices.empty() || model.response.trace_triangles.empty())
  {
    return "the rule fixes the trace triangulation but the coupon has no TraceMesh";
  }
  const auto mesh =
      ReadTraceMesh(model.response.trace_vertices, model.response.trace_triangles);
  if (mesh.vertices.size() != rule_basis.vertices.size() ||
      mesh.triangles.size() != rule_basis.triangles.size())
  {
    return "the trace mesh does not have the rule's vertex / triangle counts";
  }
  for (std::size_t v = 0; v < mesh.vertices.size(); v++)
  {
    const auto &vertex = mesh.vertices[v];
    const auto &rule_vertex = rule_basis.vertices[v];
    const std::array<double, 3> point = {vertex.point[0], vertex.point[1], vertex.point[2]};
    if (Distance(point, rule_vertex.point) > position_tolerance ||
        vertex.basis - 1 != rule_vertex.basis ||
        vertex.parent_a - 1 != rule_vertex.parent_a ||
        vertex.parent_b - 1 != rule_vertex.parent_b ||
        (rule_vertex.basis < 0 &&
         std::abs(vertex.weight_a - rule_vertex.weight_a) > 1.0e-9))
    {
      return "trace vertex " + std::to_string(v + 1) + " differs from the rule's";
    }
  }
  std::set<std::array<int, 3>> coupon_triangles, rule_triangles;
  for (auto triangle : mesh.triangles)
  {
    std::sort(triangle.begin(), triangle.end());
    coupon_triangles.insert(triangle);
  }
  for (auto triangle : rule_basis.triangles)
  {
    std::sort(triangle.begin(), triangle.end());
    rule_triangles.insert(triangle);
  }
  if (coupon_triangles != rule_triangles)
  {
    return "the trace triangulation is not the rule's";
  }
  return "";
}

// A model's basis points: the constructed basis of an angle-interpolated corner model, else
// the BasisPoints file.
std::vector<std::array<double, 3>> ModelBasisPoints(
    const config::ElectrostaticSolverData::ResponseCorrectionModelData &model_config)
{
  if (!model_config.constructed_basis_points.empty())
  {
    return model_config.constructed_basis_points;
  }
  return ReadBasisPoints(model_config.basis_points);
}

bool HasExplicitTraceMesh(
    const config::ElectrostaticSolverData::ResponseCorrectionModelData &model_config)
{
  return !model_config.trace_vertices.empty() ||
         !model_config.constructed_trace_vertices.empty();
}

// A model's explicit trace mesh: the constructed mesh of an angle-interpolated corner
// model, else the TraceMesh files.
TraceMeshData ModelTraceMesh(
    const config::ElectrostaticSolverData::ResponseCorrectionModelData &model_config)
{
  if (!model_config.constructed_trace_vertices.empty())
  {
    TraceMeshData mesh;
    mesh.vertices.reserve(model_config.constructed_trace_vertices.size());
    for (const auto &vertex : model_config.constructed_trace_vertices)
    {
      mesh.vertices.push_back({vertex.point, vertex.basis, vertex.conductor,
                               vertex.parent_a, vertex.parent_b, vertex.weight_a});
    }
    mesh.triangles = model_config.constructed_trace_triangles;
    return mesh;
  }
  MFEM_VERIFY(!model_config.trace_triangles.empty(),
              "A response TraceMesh requires both vertex and triangle files!");
  return ReadTraceMesh(model_config.trace_vertices, model_config.trace_triangles);
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

mfem::DenseMatrix PositiveSemidefiniteInverseProduct(const mfem::DenseMatrix &matrix,
                                                     const mfem::DenseMatrix &rhs)
{
  MFEM_VERIFY(matrix.Height() == matrix.Width() && rhs.Height() == matrix.Height(),
              "Incompatible response matrices for quotient-space inverse!");
  const int size = matrix.Height();
  Eigen::MatrixXd A(size, size), B(size, rhs.Width());
  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      A(i, j) = matrix(i, j);
    }
    for (int j = 0; j < rhs.Width(); j++)
    {
      B(i, j) = rhs(i, j);
    }
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensystem(A);
  MFEM_VERIFY(eigensystem.info() == Eigen::Success,
              "Failed to diagonalize a response matrix on its active quotient space!");
  const auto eigenvalues = eigensystem.eigenvalues();
  const double scale = std::max(eigenvalues.cwiseAbs().maxCoeff(), 1.0e-300);
  const double negative_tolerance = 1.0e-9 * scale;
  const double active_tolerance = 1.0e-12 * scale;
  MFEM_VERIFY(eigenvalues.minCoeff() >= -negative_tolerance,
              "Fabricated response matrix has a negative-energy mode beyond roundoff!");
  Eigen::VectorXd inverse = Eigen::VectorXd::Zero(size);
  int numerical_rank = 0;
  for (int i = 0; i < size; i++)
  {
    if (eigenvalues(i) > active_tolerance)
    {
      inverse(i) = 1.0 / eigenvalues(i);
      numerical_rank++;
    }
  }
  MFEM_VERIFY(numerical_rank > 0, "Fabricated response matrix has no active energy modes!");
  const Eigen::MatrixXd result = eigensystem.eigenvectors() * inverse.asDiagonal() *
                                 eigensystem.eigenvectors().transpose() * B;
  mfem::DenseMatrix output(size, rhs.Width());
  for (int i = 0; i < output.Height(); i++)
  {
    for (int j = 0; j < output.Width(); j++)
    {
      output(i, j) = result(i, j);
    }
  }
  return output;
}

struct DomainResponseMatrices
{
  mfem::DenseMatrix fabricated;
  mfem::DenseMatrix thin;
  mfem::DenseMatrix defect;
  mfem::DenseMatrix fixed_flux_transform;
  mfem::DenseMatrix fixed_flux_defect;
};

// Read one domain response matrix file as a dense matrix of the expected size.
mfem::DenseMatrix ReadDenseDomainResponseMatrix(const std::string &path, int expected_size)
{
  auto [size, entries] = ReadDomainResponseMatrix(path);
  MFEM_VERIFY(size == expected_size,
              "Response matrix \"" << path
                                   << "\" and basis point file have inconsistent sizes!");
  return BuildDenseMatrix(entries, expected_size, path);
}

// A model's domain response matrix: the file itself, or for an interpolated model (config
// blend) the weighted sum of its sources' files (all on the model's basis).
mfem::DenseMatrix BlendedDomainResponseMatrix(
    const config::ElectrostaticSolverData::ResponseCorrectionModelData &config,
    bool fabricated, int expected_size)
{
  if (config.blend.empty())
  {
    return ReadDenseDomainResponseMatrix(
        fabricated ? config.fabricated_matrix : config.thin_matrix, expected_size);
  }
  mfem::DenseMatrix result(expected_size);
  result = 0.0;
  for (const auto &source : config.blend)
  {
    MFEM_VERIFY(std::isfinite(source.weight),
                "Interpolated response model weights must be finite!");
    result.Add(source.weight,
               ReadDenseDomainResponseMatrix(fabricated ? source.fabricated_matrix
                                                        : source.thin_matrix,
                                             expected_size));
  }
  return result;
}

DomainResponseMatrices BuildDomainResponseMatrices(
    const config::ElectrostaticSolverData::ResponseCorrectionModelData &config,
    int expected_size, const std::vector<int> &zero_trace_indices, const Units &units)
{
  auto fabricated = BlendedDomainResponseMatrix(config, true, expected_size);
  auto thin = BlendedDomainResponseMatrix(config, false, expected_size);

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
    fixed_flux_transform = PositiveSemidefiniteInverseProduct(fabricated, thin);
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
    mfem::DenseMatrix fixed_flux_free =
        PositiveSemidefiniteInverseProduct(fabricated_free, thin_free);
    for (int i = 0; i < free_size; i++)
    {
      for (int j = 0; j < free_size; j++)
      {
        fixed_flux_transform(free_indices[i], free_indices[j]) = fixed_flux_free(i, j);
      }
    }
  }
  mfem::DenseMatrix fabricated_times_transform(expected_size),
      fixed_flux_defect(expected_size);
  mfem::Mult(fabricated, fixed_flux_transform, fabricated_times_transform);
  mfem::MultAtB(fixed_flux_transform, fabricated_times_transform, fixed_flux_defect);
  fixed_flux_defect.Add(-1.0, thin);
  // At fixed coupon flux, v_f = F⁺ T v_t and the fabricated energy operator pulled
  // back to the thin trace is Aᵀ F A. This congruence is symmetric even when quotient-
  // space pseudoinverses are required.
  return {std::move(fabricated), std::move(thin), std::move(defect),
          std::move(fixed_flux_transform), std::move(fixed_flux_defect)};
}

// The coupon interface energy a model adds to the device. A translational (2D cross-
// section) coupon spans exactly the matching radius on either side of its edges, so its
// whole-box `Q_total_ij (J)` is the energy within R of its edges. A spatial (3D box)
// coupon's box extends beyond R of its edges (pad continuations, the far under-metal
// surface) and the device keeps its own raw energy there, so such a model must add the
// localized `Q_ij (J)` column at the matching radius (`within_radius`, in the file's SI
// units) instead: adding Q_total counted the box energy outside R twice (lane J, decision
// 112(a)).
std::map<int, mfem::DenseMatrix>
ReadSurfaceResponseMatrices(const std::string &path, int expected_size,
                            std::optional<double> within_radius = std::nullopt)
{
  const Table table = ReadTable(path);
  const auto &interface_col = FindColumn(table, "interface", path);
  const auto &edge_col = FindColumn(table, "edge", path);
  const auto &basis_i = FindColumn(table, "basis_i", path);
  const auto &basis_j = FindColumn(table, "basis_j", path);
  const Column *distance_col = nullptr;
  if (within_radius)
  {
    MFEM_VERIFY(std::isfinite(*within_radius) && *within_radius > 0.0,
                "Spatial response models require a positive matching radius to select "
                "their within-R surface response!");
    for (auto it = table.cbegin(); it != table.cend(); ++it)
    {
      if (it->header_text == "R (m)")
      {
        distance_col = &*it;
      }
    }
    MFEM_VERIFY(distance_col,
                "Surface response matrix file \""
                    << path
                    << "\" of a spatial (3D box) response model has no \"R (m)\" "
                       "column: such a model adds its within-R energy Q_ij, not the "
                       "whole-box Q_total; regenerate the compact matrix with the "
                       "localized column (finalize_corner_response.py) or publish the "
                       "coupon's surface-response-matrix.csv!");
  }
  const auto &q = FindColumn(table, within_radius ? "Q_ij (J)" : "Q_total_ij (J)", path);
  const std::size_t rows = q.n_rows();
  MFEM_VERIFY(interface_col.n_rows() == rows && edge_col.n_rows() == rows &&
                  basis_i.n_rows() == rows && basis_j.n_rows() == rows &&
                  (!distance_col || distance_col->n_rows() == rows),
              "Surface response matrix columns have inconsistent lengths in \"" << path
                                                                                << "\"!");

  // Q_total is repeated for every matching radius. Deduplicate those rows for each
  // interface/edge/basis pair, then sum all physical coupon edges belonging to an
  // interface. This makes one matrix represent either a one-edge or a coupled multi-edge
  // coupon. The within-R energy is read from the rows at the matching radius only.
  using Key = std::tuple<int, int, int, int>;
  std::map<Key, double> unique;
  std::size_t radius_rows = 0;
  for (std::size_t row = 0; row < rows; row++)
  {
    if (distance_col)
    {
      const double distance = distance_col->data[row];
      MFEM_VERIFY(std::isfinite(distance) && distance > 0.0,
                  "Invalid matching radius in surface response matrix file \"" << path
                                                                               << "\"!");
      if (std::abs(distance - *within_radius) > 1.0e-6 * *within_radius)
      {
        continue;
      }
      radius_rows++;
    }
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
  MFEM_VERIFY(!distance_col || radius_rows > 0,
              "Surface response matrix file \""
                  << path << "\" has no within-R rows at the matching radius "
                  << *within_radius << " m: a spatial (3D box) response model adds its "
                  << "energy within R of its edges (Q_ij), not the whole-box Q_total!");

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

// A model's per-coupon-interface surface response matrices: the file itself, or for an
// interpolated model (config blend) the weighted sum of its sources' files.
std::map<int, mfem::DenseMatrix> BlendedSurfaceResponseMatrices(
    const config::ElectrostaticSolverData::ResponseCorrectionModelData &config,
    bool fabricated, int expected_size, std::optional<double> within_radius)
{
  if (config.blend.empty())
  {
    return ReadSurfaceResponseMatrices(fabricated ? config.fabricated_surface_matrix
                                                  : config.thin_surface_matrix,
                                       expected_size, within_radius);
  }
  std::map<int, mfem::DenseMatrix> result;
  std::set<int> interfaces;
  for (const auto &source : config.blend)
  {
    auto matrices = ReadSurfaceResponseMatrices(
        fabricated ? source.fabricated_surface_matrix : source.thin_surface_matrix,
        expected_size, within_radius);
    std::set<int> source_interfaces;
    for (auto &[interface, matrix] : matrices)
    {
      source_interfaces.insert(interface);
      auto [it, inserted] = result.emplace(interface, mfem::DenseMatrix(expected_size));
      if (inserted)
      {
        it->second = 0.0;
      }
      it->second.Add(source.weight, matrix);
    }
    if (&source == &config.blend.front())
    {
      interfaces = source_interfaces;
    }
    MFEM_VERIFY(source_interfaces == interfaces,
                "Interpolated response model sources do not share the same coupon "
                "interfaces!");
  }
  return result;
}

struct SurfaceResponseMatrices
{
  std::map<int, mfem::DenseMatrix> fabricated;
  std::map<int, mfem::DenseMatrix> defects;
};

// `matching_radius` (mesh units) is the device's edge-localization radius R; a spatial
// (3D box) model reads its within-R surface response at that radius, a translational
// model its whole-box response (equal to within R by construction of the 2D coupon box).
// The device's edge-localization radius R (mesh units): the largest EdgeDistances value
// of the target dielectric interfaces, which the automatic library matching verifies
// against the process library's MatchingRadius. Zero when no target interface localizes
// its edge energy (explicit 2D configurations without spatial models).
double TargetInterfaceMatchingRadius(const IoData &iodata,
                                     const std::vector<int> &target_interfaces)
{
  const std::set<int> target_filter(target_interfaces.begin(), target_interfaces.end());
  double radius = 0.0;
  for (const auto &[index, dielectric] : iodata.boundaries.postpro.dielectric)
  {
    if ((!target_filter.empty() && target_filter.find(index) == target_filter.end()) ||
        dielectric.type == InterfaceDielectric::DEFAULT ||
        dielectric.edge_distances.empty())
    {
      continue;
    }
    const double candidate = dielectric.edge_distances.back();
    MFEM_VERIFY(radius == 0.0 ||
                    std::abs(candidate - radius) <= 1.0e-10 * std::max(candidate, radius),
                "Response-corrected target interfaces must share one largest EdgeDistances "
                "value (the library matching radius)!");
    radius = candidate;
  }
  return radius;
}

SurfaceResponseMatrices BuildSurfaceResponseMatrices(
    const config::ElectrostaticSolverData::ResponseCorrectionModelData &config,
    int expected_size, const Units &units, double matching_radius)
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

  std::optional<double> within_radius;
  if (config.spatial_basis)
  {
    MFEM_VERIFY(std::isfinite(matching_radius) && matching_radius > 0.0,
                "A spatial (3D box) response model requires the target interfaces' "
                "matching radius to select its within-R surface response!");
    within_radius = units.Dimensionalize<Units::ValueType::LENGTH>(matching_radius);
  }
  auto fabricated =
      BlendedSurfaceResponseMatrices(config, true, expected_size, within_radius);
  auto thin = BlendedSurfaceResponseMatrices(config, false, expected_size, within_radius);
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

std::optional<std::filesystem::path> ResponseGeometryCachePath()
{
  const char *value = std::getenv("PALACE_RESPONSE_GEOMETRY_CACHE");
  if (!value || std::string_view(value).empty())
  {
    return std::nullopt;
  }
  return std::filesystem::absolute(value).lexically_normal();
}

}  // namespace

void WriteResponseGeometryCache(const std::filesystem::path &path,
                                const ResponseCorrectionData &config)
{
  nlohmann::json models = nlohmann::json::array();
  for (const auto &model : config.models)
  {
    nlohmann::json interfaces = nlohmann::json::array();
    for (const auto &interface : model.interfaces)
    {
      interfaces.push_back({{"Target", interface.target}, {"Coupon", interface.coupon}});
    }
    nlohmann::json open_paths = nlohmann::json::array();
    for (const auto &entry : model.open_contour_paths)
    {
      open_paths.push_back({{"Indices", entry.indices},
                            {"StartConductor", entry.start_conductor},
                            {"EndConductor", entry.end_conductor}});
    }
    // A curvature-family runtime model is its blend (weights x the source coupons'
    // matrices); without it the reader would fall back to the anchor's straight matrices.
    nlohmann::json blend = nlohmann::json::array();
    for (const auto &source : model.blend)
    {
      blend.push_back({{"Weight", source.weight},
                       {"FabricatedMatrix", source.fabricated_matrix},
                       {"ThinMatrix", source.thin_matrix},
                       {"FabricatedSurfaceMatrix", source.fabricated_surface_matrix},
                       {"ThinSurfaceMatrix", source.thin_surface_matrix}});
    }
    // A corner-family runtime model interpolated between nodes carries its basis and trace
    // mesh constructed at the device angle (nothing on disk describes them).
    nlohmann::json constructed_vertices = nlohmann::json::array();
    for (const auto &vertex : model.constructed_trace_vertices)
    {
      constructed_vertices.push_back({{"Point", vertex.point},
                                      {"Basis", vertex.basis},
                                      {"Conductor", vertex.conductor},
                                      {"ParentA", vertex.parent_a},
                                      {"ParentB", vertex.parent_b},
                                      {"WeightA", vertex.weight_a}});
    }
    models.push_back({{"Index", model.idx},
                      {"Name", model.name},
                      {"Topology", model.topology},
                      {"FabricatedMatrix", model.fabricated_matrix},
                      {"ThinMatrix", model.thin_matrix},
                      {"FabricatedSurfaceMatrix", model.fabricated_surface_matrix},
                      {"ThinSurfaceMatrix", model.thin_surface_matrix},
                      {"ConstructedBasisPoints", model.constructed_basis_points},
                      {"ConstructedTraceVertices", constructed_vertices},
                      {"ConstructedTraceTriangles", model.constructed_trace_triangles},
                      {"BasisPoints", model.basis_points},
                      {"TraceVertices", model.trace_vertices},
                      {"TraceTriangles", model.trace_triangles},
                      {"SpatialBasis", model.spatial_basis},
                      {"ContourGroups", model.contour_groups},
                      {"ZeroTraceIndices", model.zero_trace_indices},
                      {"OpenContourPaths", std::move(open_paths)},
                      {"InteriorTraceCount", model.interior_trace_count},
                      {"ConductorStateCount", model.conductor_state_count},
                      {"Interfaces", std::move(interfaces)},
                      {"Blend", std::move(blend)}});
  }
  nlohmann::json patches = nlohmann::json::array();
  for (const auto &patch : config.patches)
  {
    nlohmann::json claims = nlohmann::json::array();
    for (const auto &claim : patch.provenance.claims)
    {
      claims.push_back({{"Segment", claim.segment}, {"P0", claim.p0}, {"P1", claim.p1}});
    }
    nlohmann::json support_box = nullptr;
    if (patch.provenance.has_support_box)
    {
      support_box = patch.provenance.support_box;
    }
    patches.push_back({{"Model", patch.model},
                       {"SupportBox", support_box},
                       {"Chain", patch.provenance.chain},
                       {"Origin", patch.origin},
                       {"AxisU", patch.axis_u},
                       {"AxisV", patch.axis_v},
                       {"AxisW", patch.axis_w},
                       {"ConductorReferences", patch.conductor_references},
                       {"Weight", patch.weight},
                       {"LongitudinalCell", patch.longitudinal_cell},
                       {"InterpolationGroup", patch.interpolation_group},
                       {"MaxwellConductorAnchors", patch.maxwell_conductor_anchors},
                       {"MaxwellReferenceIsPEC", patch.maxwell_reference_is_pec},
                       {"Feature", patch.provenance.feature},
                       {"Segment", patch.provenance.segment},
                       {"Stretch", patch.provenance.stretch},
                       {"EdgeOffset", patch.provenance.edge_offset},
                       {"Claims", claims}});
  }
  if (path.has_parent_path())
  {
    std::filesystem::create_directories(path.parent_path());
  }
  std::ofstream output(path);
  MFEM_VERIFY(output,
              "Unable to write response-geometry cache \"" << path.string() << "\"!");
  nlohmann::json cache = {{"Version", 8},
                          {"MatchingRadius", config.matching_radius},
                          {"Models", std::move(models)},
                          {"Patches", std::move(patches)}};
  if (!config.quantum_near_match.empty())
  {
    // Quantum near-matches resolved for this geometry (block (b) DESIGN section 4).
    nlohmann::json near_matches = nlohmann::json::array();
    for (const auto &entry : config.quantum_near_match)
    {
      near_matches.push_back({{"Model", entry.model},
                              {"ModelKey", entry.model_key},
                              {"FeatureKey", entry.feature_key},
                              {"MaxDeltaQuanta", entry.max_delta_quanta},
                              {"DifferingNumbers", entry.differing_numbers},
                              {"Features", entry.features}});
    }
    cache["QuantumNearMatch"] = std::move(near_matches);
  }
  if (!config.legacy_contract.empty())
  {
    // Legacy-contract aliases resolved for this geometry (USER decision 283).
    nlohmann::json aliases = nlohmann::json::array();
    for (const auto &entry : config.legacy_contract)
    {
      aliases.push_back({{"Model", entry.model},
                         {"Key", entry.key},
                         {"ContextDigest", entry.context_digest},
                         {"Reason", entry.reason},
                         {"Features", entry.features}});
    }
    cache["LegacyContract"] = std::move(aliases);
  }
  output << cache.dump(2) << '\n';
}

ResponseCorrectionData ReadResponseGeometryCache(const std::filesystem::path &path,
                                                 const ResponseCorrectionData &request)
{
  std::ifstream input(path);
  MFEM_VERIFY(input, "Unable to read response-geometry cache \"" << path.string() << "\"!");
  nlohmann::json data;
  input >> data;
  MFEM_VERIFY(
      data.value("Version", 0) == 8,
      "Unsupported response-geometry cache version "
          << data.value("Version", 0)
          << " (version 8 carries the feature, mesh segment, chain stretch and own-edge "
             "offset of every patch, the claims, support box and chain of every spatial "
             "cluster patch, the matching radius for the continuation and vertex "
             "ownership and the quantum near-match records of the matching pass; delete a "
             "stale cache)!");
  ResponseCorrectionData result = request;
  result.library.clear();
  result.models.clear();
  result.patches.clear();
  result.legacy_contract.clear();
  result.quantum_near_match.clear();
  result.matching_radius = data.at("MatchingRadius");
  for (const auto &entry : data.at("Models"))
  {
    ResponseModelData model;
    model.idx = entry.at("Index");
    model.name = entry.value("Name", std::string{});
    model.topology = entry.value("Topology", std::string{});
    model.fabricated_matrix = entry.at("FabricatedMatrix");
    model.thin_matrix = entry.at("ThinMatrix");
    model.fabricated_surface_matrix = entry.value("FabricatedSurfaceMatrix", std::string{});
    model.thin_surface_matrix = entry.value("ThinSurfaceMatrix", std::string{});
    model.basis_points = entry.at("BasisPoints");
    model.trace_vertices = entry.value("TraceVertices", std::string{});
    model.trace_triangles = entry.value("TraceTriangles", std::string{});
    model.constructed_basis_points =
        entry.value("ConstructedBasisPoints", std::vector<std::array<double, 3>>{});
    for (const auto &value :
         entry.value("ConstructedTraceVertices", nlohmann::json::array()))
    {
      ResponseModelData::ConstructedTraceVertex vertex;
      vertex.point = value.at("Point");
      vertex.basis = value.at("Basis");
      vertex.conductor = value.at("Conductor");
      vertex.parent_a = value.at("ParentA");
      vertex.parent_b = value.at("ParentB");
      vertex.weight_a = value.at("WeightA");
      model.constructed_trace_vertices.push_back(vertex);
    }
    model.constructed_trace_triangles =
        entry.value("ConstructedTraceTriangles", std::vector<std::array<int, 3>>{});
    model.spatial_basis = entry.value("SpatialBasis", false);
    model.contour_groups = entry.value("ContourGroups", std::vector<int>{});
    model.zero_trace_indices = entry.value("ZeroTraceIndices", std::vector<int>{});
    model.interior_trace_count = entry.value("InteriorTraceCount", 0);
    model.conductor_state_count = entry.value("ConductorStateCount", 0);
    for (const auto &value : entry.value("OpenContourPaths", nlohmann::json::array()))
    {
      model.open_contour_paths.push_back(
          {value.at("Indices"), value.at("StartConductor"), value.at("EndConductor")});
    }
    for (const auto &value : entry.value("Interfaces", nlohmann::json::array()))
    {
      model.interfaces.push_back({value.at("Target"), value.at("Coupon")});
    }
    // A blended (curvature-family) runtime model is only complete with its Blend (version
    // 2 of the cache; a version-1 cache, which could describe such a model only as its
    // anchor's straight matrices, is refused above).
    for (const auto &value : entry.at("Blend"))
    {
      ResponseModelData::BlendSourceData source;
      source.weight = value.at("Weight");
      source.fabricated_matrix = value.at("FabricatedMatrix");
      source.thin_matrix = value.at("ThinMatrix");
      source.fabricated_surface_matrix =
          value.value("FabricatedSurfaceMatrix", std::string{});
      source.thin_surface_matrix = value.value("ThinSurfaceMatrix", std::string{});
      model.blend.push_back(std::move(source));
    }
    result.models.push_back(std::move(model));
  }
  for (const auto &entry : data.at("Patches"))
  {
    ResponsePatchData patch;
    patch.model = entry.at("Model");
    patch.origin = entry.at("Origin");
    patch.axis_u = entry.at("AxisU");
    patch.axis_v = entry.at("AxisV");
    patch.axis_w = entry.at("AxisW");
    patch.conductor_references = entry.at("ConductorReferences");
    patch.weight = entry.at("Weight");
    patch.longitudinal_cell = entry.at("LongitudinalCell");
    patch.interpolation_group = entry.value("InterpolationGroup", 0);
    patch.provenance.feature = entry.at("Feature");
    patch.provenance.segment = entry.at("Segment");
    patch.provenance.stretch = entry.at("Stretch");
    patch.provenance.edge_offset = entry.at("EdgeOffset");
    for (const auto &claim : entry.at("Claims"))
    {
      patch.provenance.claims.push_back(
          {claim.at("Segment"), claim.at("P0"), claim.at("P1")});
    }
    if (const auto box = entry.find("SupportBox"); box != entry.end() && !box->is_null())
    {
      patch.provenance.support_box = box->get<std::array<double, 4>>();
      patch.provenance.has_support_box = true;
    }
    patch.provenance.chain = entry.value("Chain", std::vector<std::array<double, 4>>{});
    patch.maxwell_conductor_anchors =
        entry.value("MaxwellConductorAnchors", std::vector<std::array<double, 3>>{});
    patch.maxwell_reference_is_pec = entry.value("MaxwellReferenceIsPEC", true);
    result.patches.push_back(std::move(patch));
  }
  MFEM_VERIFY(!result.models.empty() && !result.patches.empty(),
              "Response-geometry cache contains no models or patches!");
  for (const auto &entry : data.value("LegacyContract", nlohmann::json::array()))
  {
    result.legacy_contract.push_back(
        {entry.at("Model").get<std::string>(), entry.at("Key").get<std::string>(),
         entry.at("ContextDigest").get<std::string>(), entry.value("Reason", std::string{}),
         entry.value("Features", std::vector<int>{})});
  }
  for (const auto &entry : data.value("QuantumNearMatch", nlohmann::json::array()))
  {
    result.quantum_near_match.push_back(
        {entry.at("Model").get<std::string>(), entry.at("ModelKey").get<std::string>(),
         entry.at("FeatureKey").get<std::string>(), entry.value("MaxDeltaQuanta", 0.0),
         entry.value("DifferingNumbers", 0), entry.value("Features", std::vector<int>{})});
  }
  return result;
}

namespace
{

double QuadraticForm(const mfem::DenseMatrix &matrix, const Vector &x, Vector &workspace)
{
  workspace.SetSize(x.Size());
  matrix.Mult(x, workspace);
  return x * workspace;
}

}  // namespace

// The patch dry run: every automatically constructed 3D patch with its feature, model,
// weights and the segment portion it integrates (lengths and coordinates in the manifest's
// units), so that the audit can gate the construction against the manifest without a
// field solve.
void WriteSurfaceResponsePatches(const ResponseCorrectionData &data,
                                 const AutomaticResponseRequirements &requirements,
                                 const std::string &path)
{
  std::ofstream output(path);
  MFEM_VERIFY(output, "Unable to open surface-response patch dry run \"" << path << "\"!");
  std::map<int, const ResponseModelData *> models;
  for (const auto &model : data.models)
  {
    models.emplace(model.idx, &model);
  }
  output << "Patch,Feature,Topology,Model,ModelIndex,Weight,ModelWeight,QuadratureWeight,"
            "SideFactor,CouponDepth,Segment,S0,S1,OriginX,OriginY,OriginZ,AxisUX,AxisUY,"
            "AxisUZ,AxisVX,AxisVY,AxisVZ,AxisWX,AxisWY,AxisWZ,StripBegin,StripEnd\n";
  output << std::setprecision(17);
  for (std::size_t i = 0; i < data.patches.size(); i++)
  {
    const auto &patch = data.patches[i];
    const auto model = models.find(patch.model);
    MFEM_VERIFY(model != models.end(), "A dry-run patch refers to an unknown model!");
    const auto &provenance = patch.provenance;
    // The weight of a longitudinal patch carries 1 / CouponDepth (mesh units); the exported
    // weight is dimensionless in either unit system when the depth is scaled with it.
    output << i << ',' << provenance.feature << ',' << model->second->topology << ','
           << model->second->name << ',' << patch.model << ',' << patch.weight << ','
           << provenance.model_weight << ',' << provenance.quadrature_weight << ','
           << provenance.side_factor << ','
           << requirements.ScaleLength(provenance.coupon_depth) << ',' << provenance.segment
           << ',' << requirements.ScaleLength(provenance.s0) << ','
           << requirements.ScaleLength(provenance.s1);
    for (const double value : patch.origin)
    {
      output << ',' << requirements.ScaleLength(value);
    }
    for (const auto *axis : {&patch.axis_u, &patch.axis_v, &patch.axis_w})
    {
      for (const double value : *axis)
      {
        output << ',' << value;
      }
    }
    // The longitudinal cell of the surface-mortar strip: offsets along AxisW from the
    // origin; the cells of one portion tile it.
    for (const double value : patch.longitudinal_cell)
    {
      output << ',' << requirements.ScaleLength(value);
    }
    output << '\n';
  }
  MFEM_VERIFY(output.good(), "Failed writing the surface-response patch dry run!");
}

void WriteSurfaceResponseRequirements(const IoData &iodata, const Mesh &mesh,
                                      const std::string &path)
{
  const auto &configured_request = iodata.problem.type == ProblemType::ELECTROSTATIC
                                       ? iodata.solver.electrostatic.response_correction
                                       : iodata.solver.surface_response_correction;
  MFEM_VERIFY(configured_request && configured_request->IsAutomatic(),
              "Surface-response preflight requires an automatic fabrication-process "
              "response library!");

  const auto request = *configured_request;
  MaterialOperator mat_op(iodata, mesh);
  AutomaticResponseRequirements requirements(iodata.units,
                                             iodata.InputsNondimensionalized());
  AutomaticResponseStatistics statistics;
  const auto &parallel_mesh = mesh.Get();
  const bool maxwell = iodata.problem.type != ProblemType::ELECTROSTATIC;
  ResponseCorrectionData patches;
  if (parallel_mesh.Dimension() == 2 && parallel_mesh.SpaceDimension() == 2)
  {
    BuildAutomaticResponseData2D(iodata, parallel_mesh, mat_op, request, maxwell, nullptr,
                                 &requirements, &statistics);
  }
  else if (parallel_mesh.Dimension() == 3 && parallel_mesh.SpaceDimension() == 3)
  {
    patches = BuildAutomaticResponseData3D(iodata, parallel_mesh, mat_op, request, maxwell,
                                           nullptr, &requirements, &statistics);
  }
  else
  {
    MFEM_ABORT("Surface-response preflight requires a 2D or 3D solve mesh!");
  }

  requirements.SetStatistics(BuildAutomaticStatistics(parallel_mesh.GetComm(), statistics));
  const auto patches_path =
      (std::filesystem::path(path).parent_path() / "surface-response-patches.csv").string();
  auto manifest = requirements.Build();
  manifest["MeshDimension"] = parallel_mesh.Dimension();
  manifest["Maxwell"] = maxwell;
  // The identification is recorded on the root only; its placement records below include a
  // collective containment test, so every rank follows the root's decision.
  int identification_recorded =
      parallel_mesh.Dimension() == 3 && manifest.contains("Identification") ? 1 : 0;
  Mpi::Broadcast(1, &identification_recorded, 0, parallel_mesh.GetComm());
  if (identification_recorded)
  {
    // The placement's records on the dry run, computed on every rank (the geometry-only
    // ownership is deterministic; the domain-boundary containment test is collective): the
    // ownership record (decision 236) on the spatial models whose basis points the library
    // provides (a signature placeholder of a Missing key has none), the margin overlaps
    // (decision 244) and the domain-boundary exclusions (decision 258).
    const double coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
    std::map<int, std::vector<std::array<double, 3>>> basis_points;
    std::map<int, bool> spatial_basis;
    std::map<int, std::string> model_names;
    for (const auto &model : patches.models)
    {
      spatial_basis.emplace(model.idx, model.spatial_basis);
      model_names.emplace(model.idx, model.name);
      if (!model.constructed_basis_points.empty() ||
          std::filesystem::is_regular_file(model.basis_points))
      {
        basis_points.emplace(model.idx, ModelBasisPoints(model));
      }
    }
    auto BasisPoints = [&](int model_idx) -> const std::vector<std::array<double, 3>> *
    {
      const auto it = basis_points.find(model_idx);
      return it == basis_points.end() ? nullptr : &it->second;
    };
    std::vector<std::string> skipped;
    const auto boxes = CollectSpatialSupports(
        patches,
        [&](int model_idx) -> const std::vector<std::array<double, 3>> *
        {
          const auto it = basis_points.find(model_idx);
          return it == basis_points.end() || !spatial_basis.at(model_idx) ? nullptr
                                                                          : &it->second;
        },
        coordinate_scale, 3, &skipped);
    const double continuation_tolerance =
        kSignatureParameterToleranceOverRadius * patches.matching_radius;
    const auto records = FindTranslationalStretchInsideSpatialSupport(
        patches.patches, boxes, 3, continuation_tolerance);
    // The placement's continuation ownership (decision 236 (2)) on the dry run: the
    // written patches are the placed ones (clipped cells, weight 0 inside the boxes).
    const auto ownership = ApplyContinuationOwnership(
        patches.patches, boxes, 3, continuation_tolerance, patches.matching_radius);
    auto &diagnostics = manifest["Identification"]["Diagnostics"];
    diagnostics["TranslationalStretchesInsideSpatialSupport"] =
        DescribeTranslationalOwnershipRecords(records, boxes, patches, coordinate_scale,
                                              skipped, &ownership);
    diagnostics["ContinuationOwnership"] =
        DescribeContinuationOwnership(ownership, boxes, patches, coordinate_scale);
    if (!records.empty())
    {
      Mpi::Warning("{}", DescribeTranslationalOwnershipWarning(
                             diagnostics["TranslationalStretchesInsideSpatialSupport"]));
    }
    if (!ownership.cells.empty())
    {
      Mpi::Print(
          "{}", DescribeContinuationOwnershipSummary(diagnostics["ContinuationOwnership"]));
    }
    // Coupon-vs-coupon margin overlaps (decision 244): recorded here; the placement
    // aborts on a claim inside the other cluster's claims.
    const auto margin_overlaps =
        FindSpatialSupportMarginOverlaps(boxes, 3, continuation_tolerance);
    if (!margin_overlaps.empty())
    {
      diagnostics["SpatialSupportMarginOverlaps"] =
          DescribeSpatialSupportMarginOverlaps(margin_overlaps, patches, coordinate_scale);
      Mpi::Warning("{}", DescribeSpatialSupportMarginOverlapWarning(
                             diagnostics["SpatialSupportMarginOverlaps"]));
    }
    // Domain-boundary exclusions (decision 258) on the placed dry run: the excluded
    // patches are written with weight 0; their uncorrected CELL length enters the
    // inventory next to Missing (not a library gap), the portion sum alongside.
    const auto exclusions = FindDomainBoundaryExclusions(
        const_cast<mfem::ParMesh &>(parallel_mesh), patches.patches, BasisPoints,
        [&](int model_idx) { return spatial_basis.at(model_idx); },
        [&](int model_idx) { return model_names.at(model_idx); }, coordinate_scale,
        patches.matching_radius, {});
    diagnostics["DomainBoundaryExclusions"] =
        DescribeDomainBoundaryExclusions(exclusions, patches, coordinate_scale);
    const auto &exclusion_diagnostics = diagnostics["DomainBoundaryExclusions"];
    manifest["Summary"]["Counts"]["DomainBoundary"] = exclusion_diagnostics["Count"];
    manifest["Summary"]["TotalEdgeLengths"]["DomainBoundary"] =
        exclusion_diagnostics["CellLength"];
    manifest["Summary"]["DomainBoundary"] = {
        {"Patches", exclusion_diagnostics["Count"]},
        {"Features", exclusion_diagnostics["Features"]},
        {"CellLength", exclusion_diagnostics["CellLength"]},
        {"PortionLength", exclusion_diagnostics["PortionLength"]}};
    // Legacy-contract aliases (USER decision 283) in the inventory: the features matched
    // through an explicit alias (Status Exact, the legacy coupon applied), their claimed
    // length and the aliases.
    {
      nlohmann::json aliases = nlohmann::json::array();
      int legacy_features = 0;
      double legacy_length = 0.0;
      for (const auto &feature : manifest["Identification"]["Features"])
      {
        if (!feature["Match"].contains("LegacyContract"))
        {
          continue;
        }
        legacy_features++;
        legacy_length += feature["Length"].get<double>();
        const auto &record = feature["Match"]["LegacyContract"];
        auto it = std::find_if(aliases.begin(), aliases.end(), [&](const nlohmann::json &a)
                               { return a["Key"] == record["Key"]; });
        if (it == aliases.end())
        {
          aliases.push_back({{"Model", record["Model"]},
                             {"Key", record["Key"]},
                             {"ContextDigest", record["ContextDigest"]},
                             {"Reason", record["Reason"]},
                             {"Features", nlohmann::json::array()}});
          it = std::prev(aliases.end());
        }
        (*it)["Features"].push_back(feature["Id"]);
      }
      manifest["Summary"]["Counts"]["LegacyContract"] = legacy_features;
      manifest["Summary"]["TotalEdgeLengths"]["LegacyContract"] = legacy_length;
      manifest["Summary"]["LegacyContract"] = {
          {"Features", legacy_features},
          {"Length", legacy_length},
          {"Aliases", aliases},
          {"Rule", "USER decision 283: features whose contract-3 key the library maps "
                   "explicitly (LegacyContractAliases: key + context digest, verified) to "
                   "a legacy model built under the decision-236 straight-continuation "
                   "contract; counted in Exact as well (the legacy coupon is applied); "
                   "never a fallback for any other key"}};
      if (legacy_features > 0)
      {
        Mpi::Warning("{:d} feature(s) matched through legacy-contract aliases (USER "
                     "decision 283); see Summary.LegacyContract!\n",
                     legacy_features);
      }
    }
    // Quantum near-matches (block (b) DESIGN section 4, decision 303) in the inventory:
    // the cluster features matched to a model within kClusterQuantumNearMatchMaxQuanta
    // quanta (Status Exact: the same geometry at the grid), their claimed length and the
    // (ModelKey, FeatureKey) pairs — the thin-run guard and the library tooling read the
    // FeatureKey -> ModelKey map from here.
    {
      nlohmann::json pairs = nlohmann::json::array();
      int near_features = 0;
      double near_length = 0.0;
      for (const auto &feature : manifest["Identification"]["Features"])
      {
        if (!feature["Match"].contains("QuantumNearMatch"))
        {
          continue;
        }
        near_features++;
        near_length += feature["Length"].get<double>();
        const auto &record = feature["Match"]["QuantumNearMatch"];
        auto it = std::find_if(pairs.begin(), pairs.end(), [&](const nlohmann::json &a)
                               { return a["FeatureKey"] == record["FeatureKey"]; });
        if (it == pairs.end())
        {
          pairs.push_back({{"Model", feature["Match"]["Model"]},
                           {"ModelKey", record["ModelKey"]},
                           {"FeatureKey", record["FeatureKey"]},
                           {"MaxDeltaQuanta", record["MaxDeltaQuanta"]},
                           {"DifferingNumbers", record["DifferingNumbers"]["Count"]},
                           {"Features", nlohmann::json::array()}});
          it = std::prev(pairs.end());
        }
        (*it)["Features"].push_back(feature["Id"]);
      }
      manifest["Summary"]["Counts"]["QuantumNearMatched"] = near_features;
      manifest["Summary"]["TotalEdgeLengths"]["QuantumNearMatched"] = near_length;
      manifest["Summary"]["QuantumNearMatch"] = {
          {"Features", near_features},
          {"Length", near_length},
          {"MaxQuanta", kClusterQuantumNearMatchMaxQuanta},
          {"Keys", pairs},
          {"Rule",
           "block (b) DESIGN section 4 (decision 303): SpatialEdgeCluster features "
           "matched to a library model of the same topology key whose numbers lie "
           "within MaxQuanta signature quanta (1e-6 R / 1e-6 deg) of their own — "
           "the same geometry at the grid, counted in Exact (the model's coupon is "
           "applied in the feature's own canonical frame); the FeatureKey -> "
           "ModelKey map is the thin-run guard's and the library tooling's record"}};
    }
    // Knife-edge keys (decision 287 (a)): a cluster whose face-rule or box-rule readings
    // sat within the band on ANY T2 pass has a key decided at a threshold; the band census
    // of Diagnostics.SpatialSupport sums them, the features name the readings.
    {
      std::vector<std::string> knife_edge_features;
      for (const auto &feature : manifest["Identification"]["Features"])
      {
        const auto support = feature.find("SpatialSupport");
        if (support == feature.end() || !support->is_object() ||
            !support->contains("FaceRules"))
        {
          continue;
        }
        const auto &rules = (*support)["FaceRules"];
        const int hits = rules.value("ThresholdBandHits", 0) +
                         rules.value("BoxRuleThresholdBandHits", 0);
        if (hits > 0)
        {
          knife_edge_features.push_back(
              fmt::format("feature {} ({} hit(s) over {} pass(es), growth steps {})",
                          feature["Id"].get<int>(), hits, rules.value("Passes", 1),
                          (*support)["Growth"]["Steps"].dump()));
        }
      }
      if (!knife_edge_features.empty())
      {
        Mpi::Warning("{:d} spatial cluster key(s) rest on face-rule / box-rule readings "
                     "within the knife-edge band (decision 287; see Features[]."
                     "SpatialSupport.FaceRules.ThresholdBandHits and Growth.Passes):\n  "
                     "{}\n",
                     static_cast<int>(knife_edge_features.size()),
                     fmt::join(knife_edge_features, "\n  "));
      }
    }
    // The conductor-consistency gate (decision 277) needs the device trace: it is evaluated
    // by the operator at solve time (every excitation, before any energy) and reported in
    // the operator's Diagnostics.ConductorConsistency (palace.json), never by the dry run.
    manifest["Summary"]["ConductorConsistency"] = {
        {"Evaluated", false},
        {"Tolerance", SurfaceResponseOperator::kConductorConsistencyTolerance},
        {"Note",
         "solve-time gate (decision 277): the device potential at every spatial "
         "coupon's metal cross-sections on its box faces (the trace mesh's conductor "
         "vertices on the process plane) against the conductor potential, relative "
         "to the patch's trace amplitude; above Tolerance the patch is excluded "
         "(weight 0) like a DomainBoundary cell and its claimed length is left "
         "uncorrected (the B1 consumer must add it to Missing and DomainBoundary). The "
         "dry run has no "
         "device trace: see the operator record SurfaceResponse.Diagnostics."
         "ConductorConsistency of the solve"}};
    Mpi::Print(
        parallel_mesh.GetComm(),
        " Domain-boundary containment test: {:d} patches / {:d} points in {:.3f} s\n",
        exclusions.tested_patches, exclusions.tested_points, exclusions.wall_time);
    if (!exclusions.patches.empty())
    {
      Mpi::Print("{}", DescribeDomainBoundaryExclusionSummary(exclusion_diagnostics));
    }
  }
  if (Mpi::Root(parallel_mesh.GetComm()))
  {
    std::ofstream output(path);
    MFEM_VERIFY(output, "Unable to open surface-response requirements manifest \""
                            << path << "\"!");
    output << manifest.dump(2) << '\n';
    if (parallel_mesh.Dimension() == 3)
    {
      WriteSurfaceResponsePatches(patches, requirements, patches_path);
    }
  }
  Mpi::Barrier(parallel_mesh.GetComm());
  Mpi::Print(parallel_mesh.GetComm(),
             "\nSurface-response process-library preflight complete:\n Manifest: {}\n",
             path);
  if (parallel_mesh.Dimension() == 3)
  {
    Mpi::Print(parallel_mesh.GetComm(), " Patch dry run: {} ({:d} patches)\n", patches_path,
               static_cast<int>(patches.patches.size()));
  }
}

struct SurfaceResponseGeometry::Impl
{
  ResponseCorrectionData config;
  std::optional<AutomaticResponseDiagnostics> diagnostics;
  nlohmann::json statistics;
  bool maxwell = false;
  int dimension = 0;
};

namespace
{

// One feature side's maximal contiguous stretch along its chain (patch provenance feature /
// stretch): its cells, their ends (patch units), the ends of the same cells on the side's
// OWN edge (the cell ends shifted by the provenance edge offset along AxisU: identical for
// a single edge, half the separation off the midline for a pair side, the side's offset
// off the first side for a stack side), the mesh segments they lie on, and the AxisW of
// the first cell.
struct TranslationalStretch
{
  std::vector<std::size_t> patches;
  std::vector<std::array<double, 3>> ends;
  std::vector<std::array<double, 3>> edge_ends;
  std::set<int> segments;
  std::array<double, 3> direction{};
  double length = 0.0;
};

std::map<std::pair<int, int>, TranslationalStretch>
CollectTranslationalStretches(const std::vector<ResponsePatchData> &patches, int dimension)
{
  std::map<std::pair<int, int>, TranslationalStretch> stretches;
  for (std::size_t patch_idx = 0; patch_idx < patches.size(); patch_idx++)
  {
    const auto &patch = patches[patch_idx];
    if (patch.longitudinal_cell[1] <= patch.longitudinal_cell[0] ||
        patch.provenance.feature < 0 || patch.provenance.stretch < 0)
    {
      continue;  // a point-in-z (2D, spatial, explicit) patch, or no provenance
    }
    auto &stretch =
        stretches[std::make_pair(patch.provenance.feature, patch.provenance.stretch)];
    if (stretch.patches.empty())
    {
      stretch.direction = patch.axis_w;
    }
    stretch.patches.push_back(patch_idx);
    stretch.segments.insert(patch.provenance.segment);
    stretch.length += patch.longitudinal_cell[1] - patch.longitudinal_cell[0];
    for (const double offset : patch.longitudinal_cell)
    {
      std::array<double, 3> end{}, edge_end{};
      for (int d = 0; d < dimension; d++)
      {
        end[d] = patch.origin[d] + offset * patch.axis_w[d];
        edge_end[d] = end[d] + patch.provenance.edge_offset * patch.axis_u[d];
      }
      stretch.ends.push_back(end);
      stretch.edge_ends.push_back(edge_end);
    }
  }
  return stretches;
}

// A stretch continues a claim when a cell lies on the claim's mesh segment, or when it
// runs parallel to the claim, one of the ends of its cells ON THE SIDE'S OWN EDGE (the
// cell ends shifted by the provenance edge offset: a pair's cells sit on its midline, a
// stack's on its first side) abuts a claim end along the chain within the tolerance and
// within the same tolerance transversely — the own edge is the claim's edge — AND it
// extends beyond that claim end: every cell end lies on the outward side of the abutting
// end (away from the claim's other end) within the tolerance. A side whose own edge is not
// the claimed edge (the unclaimed edge of a pair, another side of a stack) does not
// continue the claim whatever its cells' proximity to the claim end (decision 252); a
// parallel stretch lying alongside the claim over the claim's own range, one end aligned
// with a claim end, is foreign: it does not go through the claim cut. Returns the claim end
// the stretch continues from (the abutting end; on the segment branch the claim end nearest
// to the stretch), nullopt when it does not continue the claim.
std::optional<std::array<double, 3>>
ContinuedClaimEnd(const TranslationalStretch &stretch,
                  const ResponsePatchData::Provenance::Claim &claim, int dimension,
                  double continuation_tolerance)
{
  auto Distance2 = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
  {
    double distance2 = 0.0;
    for (int d = 0; d < dimension; d++)
    {
      distance2 += (a[d] - b[d]) * (a[d] - b[d]);
    }
    return distance2;
  };
  if (claim.segment >= 0 && stretch.segments.count(claim.segment))
  {
    double nearest0 = mfem::infinity(), nearest1 = mfem::infinity();
    for (const auto &end : stretch.ends)
    {
      nearest0 = std::min(nearest0, Distance2(end, claim.p0));
      nearest1 = std::min(nearest1, Distance2(end, claim.p1));
    }
    return nearest0 <= nearest1 ? claim.p0 : claim.p1;
  }
  const double parallel_cosine = std::cos(kSignatureAngleToleranceDegrees * M_PI / 180.0);
  std::array<double, 3> tangent{};
  double length2 = 0.0, cosine = 0.0;
  for (int d = 0; d < dimension; d++)
  {
    tangent[d] = claim.p1[d] - claim.p0[d];
    length2 += tangent[d] * tangent[d];
  }
  if (length2 <= 0.0)
  {
    return std::nullopt;
  }
  for (int d = 0; d < dimension; d++)
  {
    tangent[d] /= std::sqrt(length2);
    cosine += tangent[d] * stretch.direction[d];
  }
  if (std::abs(cosine) < parallel_cosine)
  {
    return std::nullopt;
  }
  for (const auto &end : stretch.edge_ends)
  {
    for (const auto *claim_end : {&claim.p0, &claim.p1})
    {
      double along = 0.0;
      for (int d = 0; d < dimension; d++)
      {
        along += (end[d] - (*claim_end)[d]) * tangent[d];
      }
      const double distance2 = Distance2(end, *claim_end);
      if (!(std::abs(along) <= continuation_tolerance &&
            distance2 - along * along <= continuation_tolerance * continuation_tolerance))
      {
        continue;
      }
      // Outward = from the claim's other end toward this abutting end.
      const double outward = (claim_end == &claim.p1) ? 1.0 : -1.0;
      const bool beyond = std::all_of(stretch.edge_ends.begin(), stretch.edge_ends.end(),
                                      [&](const std::array<double, 3> &other)
                                      {
                                        double along_outward = 0.0;
                                        for (int d = 0; d < dimension; d++)
                                        {
                                          along_outward += (other[d] - (*claim_end)[d]) *
                                                           tangent[d] * outward;
                                        }
                                        return along_outward >= -continuation_tolerance;
                                      });
      if (beyond)
      {
        return *claim_end;
      }
    }
  }
  return std::nullopt;
}

// The claim end of the support the stretch continues from (the first claim in the
// support's order that it continues), nullopt when it continues none.
std::optional<std::array<double, 3>>
ContinuedSupportClaimEnd(const TranslationalStretch &stretch,
                         const SpatialSupportBounds &support, int dimension,
                         double continuation_tolerance)
{
  for (const auto &claim : support.claims)
  {
    if (const auto end =
            ContinuedClaimEnd(stretch, claim, dimension, continuation_tolerance))
    {
      return end;
    }
  }
  return std::nullopt;
}

// The offsets along a cell's AxisW line (relative to the patch origin) strictly inside a
// box, intersected with the cell [c0, c1]; nullopt when the cell does not enter the box.
std::optional<std::array<double, 2>> CellInsideBox(const ResponsePatchData &patch,
                                                   const SpatialSupportBounds &box,
                                                   int dimension)
{
  double lo = patch.longitudinal_cell[0], hi = patch.longitudinal_cell[1];
  for (int d = 0; d < dimension; d++)
  {
    const double w = patch.axis_w[d];
    if (std::abs(w) <= 1.0e-14)
    {
      if (!(patch.origin[d] > box.min[d] && patch.origin[d] < box.max[d]))
      {
        return std::nullopt;
      }
      continue;
    }
    const double a = (box.min[d] - patch.origin[d]) / w;
    const double b = (box.max[d] - patch.origin[d]) / w;
    lo = std::max(lo, std::min(a, b));
    hi = std::min(hi, std::max(a, b));
  }
  const double scale =
      std::max(1.0, patch.longitudinal_cell[1] - patch.longitudinal_cell[0]);
  if (hi - lo <= 1.0e-12 * scale)
  {
    return std::nullopt;
  }
  return std::array<double, 2>{lo, hi};
}

}  // namespace

std::vector<TranslationalOwnershipRecord> FindTranslationalStretchInsideSpatialSupport(
    const std::vector<ResponsePatchData> &patches,
    const std::vector<SpatialSupportBounds> &supports, int dimension,
    double continuation_tolerance)
{
  const auto stretches = CollectTranslationalStretches(patches, dimension);
  std::vector<TranslationalOwnershipRecord> records;
  for (const auto &entry : stretches)
  {
    const auto &key = entry.first;
    const auto &stretch = entry.second;
    for (const auto &support : supports)
    {
      // Strictly inside: a stretch touching the box face is not inside.
      const bool inside =
          std::all_of(stretch.ends.begin(), stretch.ends.end(),
                      [&](const std::array<double, 3> &end)
                      {
                        for (int d = 0; d < dimension; d++)
                        {
                          if (!(end[d] > support.min[d] && end[d] < support.max[d]))
                          {
                            return false;
                          }
                        }
                        return true;
                      });
      if (!inside)
      {
        continue;
      }
      TranslationalOwnershipRecord record;
      record.feature = key.first;
      record.stretch = key.second;
      record.first_patch = stretch.patches.front();
      record.patch_count = stretch.patches.size();
      record.spatial_patch = support.patch;
      record.length = stretch.length;
      record.lo = record.hi = stretch.ends.front();
      for (const auto &end : stretch.ends)
      {
        for (int d = 0; d < 3; d++)
        {
          record.lo[d] = std::min(record.lo[d], end[d]);
          record.hi[d] = std::max(record.hi[d], end[d]);
        }
      }
      record.continuation =
          ContinuedSupportClaimEnd(stretch, support, dimension, continuation_tolerance)
              .has_value();
      records.push_back(record);
    }
  }
  return records;
}

std::vector<SpatialSupportBounds> CollectSpatialSupports(
    const ResponseCorrectionData &config,
    const std::function<const std::vector<std::array<double, 3>> *(int model_idx)>
        &basis_points,
    double coordinate_scale, int dimension, std::vector<std::string> *skipped)
{
  std::unordered_map<int, const ResponseModelData *> models;
  for (const auto &model : config.models)
  {
    models.emplace(model.idx, &model);
  }
  std::vector<SpatialSupportBounds> supports;
  for (std::size_t patch_idx = 0; patch_idx < config.patches.size(); patch_idx++)
  {
    const auto &patch = config.patches[patch_idx];
    const auto model_it = models.find(patch.model);
    MFEM_VERIFY(model_it != models.end(), "Unknown response model!");
    if (!model_it->second->spatial_basis)
    {
      continue;
    }
    const auto *points = basis_points(patch.model);
    SpatialSupportBounds support;
    support.patch = patch_idx;
    support.claims = patch.provenance.claims;
    support.has_support_box = patch.provenance.has_support_box;
    support.support_box = patch.provenance.support_box;
    support.chain = patch.provenance.chain;
    // A contract-3 placeholder without basis points (a signature-only library, the
    // preflight of a Missing key): the box is the Signature's support box placed by the
    // patch frame, R above and below the plane (the coupon's cap and substrate reach at
    // least R), so that the ownership records of the dry run are complete.
    std::vector<std::array<double, 3>> box_corners;
    if (!points && support.has_support_box)
    {
      const auto &b = support.support_box;
      for (const double x : {b[0], b[2]})
      {
        for (const double y : {b[1], b[3]})
        {
          for (const double z : {-1.0, 1.0})
          {
            box_corners.push_back({x * config.matching_radius, y * config.matching_radius,
                                   z * config.matching_radius});
          }
        }
      }
      points = &box_corners;
      support.from_signature_box = true;
    }
    if (!points)
    {
      if (skipped)
      {
        skipped->push_back(model_it->second->name);
      }
      continue;
    }
    bool first = true;
    for (const auto &local : *points)
    {
      for (int d = 0; d < dimension; d++)
      {
        const double coordinate =
            patch.origin[d] + (local[0] * patch.axis_u[d] + local[1] * patch.axis_v[d] +
                               local[2] * patch.axis_w[d]) /
                                  coordinate_scale;
        support.min[d] = first ? coordinate : std::min(support.min[d], coordinate);
        support.max[d] = first ? coordinate : std::max(support.max[d], coordinate);
      }
      first = false;
    }
    supports.push_back(std::move(support));
  }
  return supports;
}

nlohmann::json DescribeTranslationalOwnershipRecords(
    const std::vector<TranslationalOwnershipRecord> &records,
    const std::vector<SpatialSupportBounds> &supports, const ResponseCorrectionData &config,
    double coordinate_scale, const std::vector<std::string> &skipped,
    const ContinuationOwnership *ownership)
{
  std::unordered_map<int, const ResponseModelData *> models;
  for (const auto &model : config.models)
  {
    models.emplace(model.idx, &model);
  }
  auto Scaled = [&](const std::array<double, 3> &point)
  {
    return std::array<double, 3>{point[0] * coordinate_scale, point[1] * coordinate_scale,
                                 point[2] * coordinate_scale};
  };
  nlohmann::json entries = nlohmann::json::array();
  int continuation_count = 0, foreign_count = 0;
  double continuation_length = 0.0, foreign_length = 0.0, owned_length = 0.0;
  for (const auto &record : records)
  {
    const auto support = std::find_if(supports.begin(), supports.end(), [&](const auto &s)
                                      { return s.patch == record.spatial_patch; });
    MFEM_VERIFY(support != supports.end(), "An ownership record names no spatial support!");
    const auto &patch = config.patches[record.first_patch];
    const auto &spatial = config.patches[record.spatial_patch];
    const double length = record.length * coordinate_scale;
    (record.continuation ? continuation_count : foreign_count)++;
    (record.continuation ? continuation_length : foreign_length) += length;
    double owned = 0.0;
    if (ownership)
    {
      const auto it = ownership->owned_by_stretch.find(
          std::make_tuple(record.feature, record.stretch, record.spatial_patch));
      owned = it == ownership->owned_by_stretch.end() ? 0.0 : it->second * coordinate_scale;
    }
    owned_length += owned;
    entries.push_back(
        {{"Feature", record.feature},
         {"OwnedLength", owned},
         {"Stretch", record.stretch},
         {"Model", models.at(patch.model)->name},
         {"Topology", models.at(patch.model)->topology},
         {"FirstPatch", record.first_patch},
         {"Patches", record.patch_count},
         {"Length", length},
         {"Extent", {{"Min", Scaled(record.lo)}, {"Max", Scaled(record.hi)}}},
         {"SpatialPatch", record.spatial_patch},
         {"SpatialFeature", spatial.provenance.feature},
         {"SpatialModel", models.at(spatial.model)->name},
         {"Box", {{"Min", Scaled(support->min)}, {"Max", Scaled(support->max)}}},
         {"Class", record.continuation ? "Continuation" : "Foreign"}});
  }
  return {
      {"Count", static_cast<int>(records.size())},
      {"Length", continuation_length + foreign_length},
      {"Continuation", {{"Count", continuation_count}, {"Length", continuation_length}}},
      {"Foreign", {{"Count", foreign_count}, {"Length", foreign_length}}},
      // The continuation ownership's removed length: per record (the stretch's cells inside
      // that box, clipped exactly) and in total over every stretch, inside or straddling.
      {"OwnedLength", owned_length},
      {"OwnedLengthTotal", ownership ? ownership->owned_length * coordinate_scale : 0.0},
      {"SpatialSupports", static_cast<int>(supports.size())},
      {"SpatialSupportsWithoutBasisPoints", skipped},
      {"ContinuationToleranceOverR", kSignatureParameterToleranceOverRadius},
      {"Records", std::move(entries)},
      {"Rule",
       "decision 236 (2026-10-01): a translational stretch (one feature side's maximal "
       "contiguous stretch along its chain, the patch provenance Feature / Stretch) whose "
       "every longitudinal cell lies strictly inside the box of a spatial support (the "
       "bbox of the model's basis points placed by the patch frame: claims + 3R past every "
       "claim-cut end along its edge, the 2R continuation + R padding; 2R transversely) is "
       "recorded, never an abort. Continuation = the stretch continues a claimed portion "
       "of that cluster through the claim cut by the side's OWN edge (decision 252: a cell "
       "on the claim's mesh segment, or parallel to the claim within the signature angle "
       "tolerance with a cell end on the side's own edge (the cell end shifted by the "
       "provenance EdgeOffset along AxisU: a pair's cells sit on the midline, a stack's on "
       "the first side) abutting a claim end along the chain within "
       "ContinuationToleranceOverR x R and within the same tolerance transversely, AND "
       "extending beyond that claim end: every cell end on the outward side of the "
       "abutting end within the tolerance; the side of a pair or stack whose own edge is "
       "not claimed is Foreign whatever its cells' proximity; a parallel stretch alongside "
       "the claim over its own range is Foreign): the coupon continues "
       "the claim straight to the box face, so its defect and the stretch's own patches "
       "correct the same surface (the double count that the continuation ownership at "
       "placement removes). Foreign = any other stretch: absent "
       "from the coupon's twins, a second-order model mismatch of the coupon, not a double "
       "count. Judged per stretch, never per cell: the first cells of every stack portion "
       "adjacent to a cluster lie inside its box legitimately. Lengths and coordinates in "
       "mesh units"}};
}

std::string DescribeTranslationalOwnershipWarning(const nlohmann::json &diagnostics)
{
  std::string lines;
  for (const auto &entry : diagnostics["Records"])
  {
    lines += fmt::format(
        "  feature {} stretch {} ({}, {} patches from patch {}, {:.6e} mesh units) inside "
        "spatial patch {} ({}): {}\n",
        entry["Feature"].get<int>(), entry["Stretch"].get<int>(),
        entry["Model"].get<std::string>(), entry["Patches"].get<std::size_t>(),
        entry["FirstPatch"].get<std::size_t>() + 1, entry["Length"].get<double>(),
        entry["SpatialPatch"].get<std::size_t>() + 1,
        entry["SpatialModel"].get<std::string>(), entry["Class"].get<std::string>());
  }
  return fmt::format(
      "{:d} translational response-correction stretch(es) ({:.6e} mesh units; Continuation "
      "{:d} / Foreign {:d}) lie wholly inside a spatial patch's matching volume (decision "
      "236, recorded under Diagnostics.TranslationalStretchesInsideSpatialSupport):\n{}",
      diagnostics["Count"].get<int>(), diagnostics["Length"].get<double>(),
      diagnostics["Continuation"]["Count"].get<int>(),
      diagnostics["Foreign"]["Count"].get<int>(), lines);
}

ContinuationOwnership
ApplyContinuationOwnership(std::vector<ResponsePatchData> &patches,
                           const std::vector<SpatialSupportBounds> &supports, int dimension,
                           double continuation_tolerance, double matching_radius)
{
  ContinuationOwnership ownership;
  const auto stretches = CollectTranslationalStretches(patches, dimension);
  for (const auto &entry : stretches)
  {
    const auto &key = entry.first;
    const auto &stretch = entry.second;
    // The supports (by index) this stretch continues, with the claim end it continues from.
    std::vector<std::pair<std::size_t, std::array<double, 3>>> continued;
    for (std::size_t s = 0; s < supports.size(); s++)
    {
      if (const auto end = ContinuedSupportClaimEnd(stretch, supports[s], dimension,
                                                    continuation_tolerance))
      {
        continued.emplace_back(s, *end);
      }
    }
    if (continued.empty())
    {
      continue;
    }
    for (const std::size_t patch_idx : stretch.patches)
    {
      auto &patch = patches[patch_idx];
      const double c0 = patch.longitudinal_cell[0], c1 = patch.longitudinal_cell[1];
      const double cell_length = c1 - c0;
      // The inside interval per owning support, in ascending spatial patch order, with the
      // continued claim end's offset along the cell's AxisW line.
      struct Owner
      {
        std::size_t support = 0;
        std::array<double, 2> inside{};
        double claim_end = 0.0;
      };
      std::vector<Owner> inside;
      for (const auto &[s, claim_end] : continued)
      {
        if (const auto interval = CellInsideBox(patch, supports[s], dimension))
        {
          double offset = 0.0;
          for (int d = 0; d < dimension; d++)
          {
            offset += (claim_end[d] - patch.origin[d]) * patch.axis_w[d];
          }
          inside.push_back({s, *interval, offset});
        }
      }
      if (inside.empty())
      {
        continue;
      }
      std::sort(inside.begin(), inside.end(), [&](const Owner &a, const Owner &b)
                { return supports[a.support].patch < supports[b.support].patch; });
      // The union of the inside intervals must leave at most one kept piece: a box's inside
      // interval touches a cell end (the stretch starts at the claim cut inside the box, so
      // the box's inside interval along the stretch begins at the stretch's first end), and
      // two boxes continued by one straight stretch leave no gap between their inside
      // intervals (both contain that first end). Verified, not assumed.
      const double scale = std::max(1.0, cell_length);
      std::vector<std::array<double, 2>> sorted_inside;
      for (const auto &owner : inside)
      {
        sorted_inside.push_back(owner.inside);
      }
      std::sort(sorted_inside.begin(), sorted_inside.end());
      double removed_lo = sorted_inside.front()[0], removed_hi = sorted_inside.front()[1];
      for (const auto &interval : sorted_inside)
      {
        MFEM_VERIFY(interval[0] <= removed_hi + 1.0e-12 * scale,
                    "Continuation ownership: the boxes owning the longitudinal cell of "
                    "patch "
                        << patch_idx + 1 << " leave a gap inside the cell!");
        removed_hi = std::max(removed_hi, interval[1]);
      }
      const bool at_begin = removed_lo <= c0 + 1.0e-12 * scale;
      const bool at_end = removed_hi >= c1 - 1.0e-12 * scale;
      MFEM_VERIFY(at_begin || at_end,
                  "Continuation ownership: the box of spatial patch "
                      << supports[inside.front().support].patch + 1
                      << " lies strictly inside the longitudinal cell of patch "
                      << patch_idx + 1 << " (feature " << key.first << ", stretch "
                      << key.second
                      << "): a translational cell longer than a coupon box along its own "
                         "direction cannot be clipped into one piece!");
      const double owned_length = std::min(removed_hi - removed_lo, cell_length);
      const double kept_lo = at_begin ? removed_hi : c0;
      const double kept_hi = at_begin ? c1 : removed_lo;
      const double kept_length = std::max(0.0, kept_hi - kept_lo);
      const bool wholly = kept_length <= 1.0e-12 * scale;
      // Attribution per owner: the owned interval is split at the midpoints between the
      // consecutive continued claim ends along the cell (one owner takes the whole owned
      // length).
      ContinuationOwnership::Cell cell;
      cell.patch = patch_idx;
      cell.feature = key.first;
      cell.stretch = key.second;
      cell.cell_length = cell_length;
      cell.owned_length = owned_length;
      if (inside.size() == 1)
      {
        cell.owners.push_back(supports[inside.front().support].patch);
        cell.attributed.push_back(owned_length);
      }
      else
      {
        std::vector<std::pair<double, std::size_t>> claim_ends;
        for (std::size_t i = 0; i < inside.size(); i++)
        {
          claim_ends.emplace_back(inside[i].claim_end, i);
        }
        std::sort(claim_ends.begin(), claim_ends.end());
        std::vector<double> cuts = {removed_lo};
        for (std::size_t i = 1; i < claim_ends.size(); i++)
        {
          cuts.push_back(std::clamp(0.5 * (claim_ends[i - 1].first + claim_ends[i].first),
                                    removed_lo, removed_hi));
        }
        cuts.push_back(removed_hi);
        std::vector<double> share(inside.size(), 0.0);
        for (std::size_t i = 0; i < claim_ends.size(); i++)
        {
          share[claim_ends[i].second] = std::max(0.0, cuts[i + 1] - cuts[i]);
        }
        for (std::size_t i = 0; i < inside.size(); i++)
        {
          cell.owners.push_back(supports[inside[i].support].patch);
          cell.attributed.push_back(share[i]);
        }
        ownership.shared_cells++;
        ownership.shared_length += owned_length;
      }
      for (std::size_t i = 0; i < cell.owners.size(); i++)
      {
        ownership
            .owned_by_stretch[std::make_tuple(key.first, key.second, cell.owners[i])] +=
            cell.attributed[i];
        ownership.owned_by_support[cell.owners[i]] += cell.attributed[i];
      }
      ownership.owned_length += owned_length;
      // Clip: the kept interval becomes the cell of a patch at its midpoint (the quadrature
      // point of the kept piece), weight and quadrature weight scaled by kept / cell.
      const double fraction = wholly ? 0.0 : kept_length / cell_length;
      const double shift = wholly ? 0.0 : 0.5 * (kept_lo + kept_hi);
      for (int d = 0; d < 3; d++)
      {
        patch.origin[d] += shift * patch.axis_w[d];
        for (auto &anchor : patch.maxwell_conductor_anchors)
        {
          anchor[d] += shift * patch.axis_w[d];
        }
      }
      patch.longitudinal_cell =
          wholly ? std::array<double, 2>{0.0, 0.0}
                 : std::array<double, 2>{-0.5 * kept_length, 0.5 * kept_length};
      patch.weight *= fraction;
      patch.provenance.quadrature_weight *= fraction;
      (wholly ? ownership.wholly_owned_cells : ownership.clipped_cells)++;
      ownership.cells.push_back(std::move(cell));
    }
  }
  std::sort(ownership.cells.begin(), ownership.cells.end(),
            [](const auto &a, const auto &b) { return a.patch < b.patch; });

  // Vertex ownership (rule B4): the vertex patches against the chain piece ends of every
  // contract-3 support, in the support's local frame (units of R).
  const double R = matching_radius;
  const double snap = kSupportFaceSnapOverRadius;
  for (std::size_t patch_idx = 0; patch_idx < patches.size(); patch_idx++)
  {
    auto &patch = patches[patch_idx];
    const auto &provenance = patch.provenance;
    const bool vertex_patch = provenance.coupon_depth == 0.0 && provenance.claims.empty() &&
                              provenance.stretch < 0 && !provenance.has_support_box;
    if (!vertex_patch || patch.weight <= 0.0)
    {
      continue;
    }
    ContinuationOwnership::Vertex record;
    record.patch = patch_idx;
    record.feature = provenance.feature;
    for (const auto &support : supports)
    {
      if (!support.has_support_box || support.chain.empty())
      {
        continue;
      }
      const auto &owner = patches[support.patch];
      std::array<double, 3> local{};
      for (int d = 0; d < dimension; d++)
      {
        const double r = patch.origin[d] - owner.origin[d];
        local[0] += r * owner.axis_u[d] / R;
        local[1] += r * owner.axis_v[d] / R;
        local[2] += r * owner.axis_w[d] / R;
      }
      if (std::abs(local[2]) > kSignatureParameterToleranceOverRadius)
      {
        continue;  // another plane
      }
      const auto &box = support.support_box;
      if (local[0] < box[0] - snap || local[0] > box[2] + snap ||
          local[1] < box[1] - snap || local[1] > box[3] + snap)
      {
        continue;
      }
      double end_distance = std::numeric_limits<double>::infinity();
      for (const auto &piece : support.chain)
      {
        end_distance =
            std::min({end_distance, std::hypot(local[0] - piece[0], local[1] - piece[1]),
                      std::hypot(local[0] - piece[2], local[1] - piece[3])});
      }
      if (end_distance > kSignatureParameterToleranceOverRadius)
      {
        continue;
      }
      if (record.owners.empty())
      {
        record.face_distance_over_r = std::min(
            {local[0] - box[0], box[2] - local[0], local[1] - box[1], box[3] - local[1]});
        record.arm_outside_box = record.face_distance_over_r < 1.0;
        record.lost_arm_length_over_r = std::max(0.0, 1.0 - record.face_distance_over_r);
        record.chain_end_distance_over_r = end_distance;
      }
      record.owners.push_back(support.patch);
    }
    if (record.owners.empty())
    {
      continue;
    }
    std::sort(record.owners.begin(), record.owners.end());
    if (record.owners.size() > 1)
    {
      ownership.shared_vertices++;
    }
    patch.weight = 0.0;
    ownership.vertices.push_back(std::move(record));
  }
  return ownership;
}

std::vector<SpatialSupportMarginOverlap>
FindSpatialSupportMarginOverlaps(const std::vector<SpatialSupportBounds> &supports,
                                 int dimension, double continuation_tolerance)
{
  const double parallel_cosine = std::cos(kSignatureAngleToleranceDegrees * M_PI / 180.0);
  // The straight continuations of a support's claim CUT ends to its box face.
  struct Continuation
  {
    std::array<double, 3> p0{}, p1{};
  };
  auto Continuations = [&](const SpatialSupportBounds &support)
  {
    std::vector<Continuation> continuations;
    for (std::size_t i = 0; i < support.claims.size(); i++)
    {
      const auto &claim = support.claims[i];
      std::array<double, 3> tangent{};
      double length = 0.0;
      for (int d = 0; d < dimension; d++)
      {
        tangent[d] = claim.p1[d] - claim.p0[d];
        length += tangent[d] * tangent[d];
      }
      length = std::sqrt(length);
      if (length <= 0.0)
      {
        continue;
      }
      for (int d = 0; d < dimension; d++)
      {
        tangent[d] /= length;
      }
      for (const double sign : {-1.0, 1.0})
      {
        const auto &end = sign > 0.0 ? claim.p1 : claim.p0;
        bool shared = false;
        for (std::size_t j = 0; j < support.claims.size() && !shared; j++)
        {
          if (j == i)
          {
            continue;
          }
          for (const auto *other : {&support.claims[j].p0, &support.claims[j].p1})
          {
            double distance2 = 0.0;
            for (int d = 0; d < dimension; d++)
            {
              distance2 += (end[d] - (*other)[d]) * (end[d] - (*other)[d]);
            }
            shared |= distance2 <= continuation_tolerance * continuation_tolerance;
          }
        }
        if (shared)
        {
          continue;
        }
        // Distance along the outward tangent from the cut end to the box boundary.
        double reach = mfem::infinity();
        for (int d = 0; d < dimension; d++)
        {
          const double w = sign * tangent[d];
          if (std::abs(w) <= 1.0e-14)
          {
            continue;
          }
          for (const double face : {support.min[d], support.max[d]})
          {
            const double s = (face - end[d]) / w;
            if (s > 1.0e-12)
            {
              reach = std::min(reach, s);
            }
          }
        }
        if (!std::isfinite(reach))
        {
          continue;
        }
        Continuation continuation;
        continuation.p0 = end;
        for (int d = 0; d < 3; d++)
        {
          continuation.p1[d] = end[d] + reach * sign * tangent[d];
        }
        continuations.push_back(continuation);
      }
    }
    return continuations;
  };
  // Length of segment b lying on segment a: parallel within the signature angle tolerance,
  // within continuation_tolerance of a's line transversely, the overlap of the projections.
  auto Overlap = [&](const std::array<double, 3> &a0, const std::array<double, 3> &a1,
                     const std::array<double, 3> &b0, const std::array<double, 3> &b1)
  {
    std::array<double, 3> ta{}, tb{};
    double la = 0.0, lb = 0.0, cosine = 0.0;
    for (int d = 0; d < dimension; d++)
    {
      ta[d] = a1[d] - a0[d];
      tb[d] = b1[d] - b0[d];
      la += ta[d] * ta[d];
      lb += tb[d] * tb[d];
    }
    la = std::sqrt(la);
    lb = std::sqrt(lb);
    if (la <= 0.0 || lb <= 0.0)
    {
      return 0.0;
    }
    for (int d = 0; d < dimension; d++)
    {
      cosine += ta[d] * tb[d] / (la * lb);
    }
    if (std::abs(cosine) < parallel_cosine)
    {
      return 0.0;
    }
    double s0 = 0.0, s1 = 0.0;
    for (const auto *b : {&b0, &b1})
    {
      double along = 0.0, distance2 = 0.0;
      for (int d = 0; d < dimension; d++)
      {
        const double delta = (*b)[d] - a0[d];
        along += delta * ta[d] / la;
        distance2 += delta * delta;
      }
      if (distance2 - along * along > continuation_tolerance * continuation_tolerance)
      {
        return 0.0;
      }
      (b == &b0 ? s0 : s1) = along;
    }
    if (s0 > s1)
    {
      std::swap(s0, s1);
    }
    return std::max(0.0, std::min(s1, la) - std::max(s0, 0.0));
  };
  // A claim SEGMENT of a entering the bounding box of b's claims (the claims hull) by a
  // positive length, strictly inside it beyond continuation_tolerance (both ends outside
  // with the segment crossing the hull counts; an end strictly inside counts); in a
  // direction where b's claims are coplanar (the plan of a planar cluster) b's box.
  auto ClaimInsideHull = [&](const SpatialSupportBounds &a, const SpatialSupportBounds &b)
  {
    if (b.claims.empty())
    {
      return false;
    }
    std::array<double, 3> lo = b.claims.front().p0, hi = b.claims.front().p0;
    for (const auto &claim : b.claims)
    {
      for (const auto *p : {&claim.p0, &claim.p1})
      {
        for (int d = 0; d < 3; d++)
        {
          lo[d] = std::min(lo[d], (*p)[d]);
          hi[d] = std::max(hi[d], (*p)[d]);
        }
      }
    }
    for (int d = 0; d < dimension; d++)
    {
      if (hi[d] - lo[d] <= continuation_tolerance)
      {
        lo[d] = b.min[d];
        hi[d] = b.max[d];
      }
    }
    for (const auto &claim : a.claims)
    {
      // The parameter interval of the segment p0 + t (p1 - p0), t in [0, 1], inside the
      // open hull shrunk by the tolerance (slab clipping).
      double t_lo = 0.0, t_hi = 1.0;
      for (int d = 0; d < dimension && t_lo < t_hi; d++)
      {
        const double slab_lo = lo[d] + continuation_tolerance;
        const double slab_hi = hi[d] - continuation_tolerance;
        const double direction = claim.p1[d] - claim.p0[d];
        if (std::abs(direction) <= 1.0e-14 * std::max(1.0, std::abs(slab_hi - slab_lo)))
        {
          if (!(claim.p0[d] > slab_lo && claim.p0[d] < slab_hi))
          {
            t_lo = t_hi = 0.0;
          }
          continue;
        }
        const double ta = (slab_lo - claim.p0[d]) / direction;
        const double tb = (slab_hi - claim.p0[d]) / direction;
        t_lo = std::max(t_lo, std::min(ta, tb));
        t_hi = std::min(t_hi, std::max(ta, tb));
      }
      if (t_hi > t_lo)
      {
        return true;
      }
    }
    return false;
  };
  std::vector<SpatialSupportMarginOverlap> overlaps;
  for (std::size_t i = 0; i < supports.size(); i++)
  {
    for (std::size_t j = i + 1; j < supports.size(); j++)
    {
      const auto &a = supports[i];
      const auto &b = supports[j];
      if (a.claims.empty() || b.claims.empty())
      {
        continue;  // corner / vertex supports carry no claims: the cluster priority rule
      }
      double scale = 1.0;
      for (int d = 0; d < dimension; d++)
      {
        scale = std::max({scale, std::abs(a.min[d]), std::abs(a.max[d]), std::abs(b.min[d]),
                          std::abs(b.max[d])});
      }
      const double tolerance = 1.0e-12 * scale;
      SpatialSupportMarginOverlap overlap;
      bool interior = true;
      for (int d = 0; d < dimension; d++)
      {
        overlap.overlap_min[d] = std::max(a.min[d], b.min[d]);
        overlap.overlap_max[d] = std::min(a.max[d], b.max[d]);
        interior &= overlap.overlap_max[d] - overlap.overlap_min[d] > tolerance;
      }
      if (!interior)
      {
        continue;
      }
      overlap.first_patch = a.patch;
      overlap.second_patch = b.patch;
      overlap.claim_in_hull = ClaimInsideHull(a, b) || ClaimInsideHull(b, a);
      const auto a_continuations = Continuations(a);
      const auto b_continuations = Continuations(b);
      for (const auto &c : a_continuations)
      {
        for (const auto &claim : b.claims)
        {
          overlap.first_margin_over_second_claims +=
              Overlap(c.p0, c.p1, claim.p0, claim.p1);
        }
        for (const auto &other : b_continuations)
        {
          overlap.margin_over_margin += Overlap(c.p0, c.p1, other.p0, other.p1);
        }
      }
      for (const auto &c : b_continuations)
      {
        for (const auto &claim : a.claims)
        {
          overlap.second_margin_over_first_claims +=
              Overlap(c.p0, c.p1, claim.p0, claim.p1);
        }
      }
      overlaps.push_back(overlap);
    }
  }
  return overlaps;
}

nlohmann::json
DescribeContinuationOwnership(const ContinuationOwnership &ownership,
                              const std::vector<SpatialSupportBounds> &supports,
                              const ResponseCorrectionData &config, double coordinate_scale)
{
  std::unordered_map<int, const ResponseModelData *> models;
  for (const auto &model : config.models)
  {
    models.emplace(model.idx, &model);
  }
  nlohmann::json cells = nlohmann::json::array();
  for (const auto &cell : ownership.cells)
  {
    const auto &patch = config.patches[cell.patch];
    nlohmann::json owners = nlohmann::json::array();
    for (std::size_t i = 0; i < cell.owners.size(); i++)
    {
      const auto &spatial = config.patches[cell.owners[i]];
      owners.push_back({{"SpatialPatch", cell.owners[i]},
                        {"SpatialFeature", spatial.provenance.feature},
                        {"SpatialModel", models.at(spatial.model)->name},
                        {"Length", cell.attributed[i] * coordinate_scale}});
    }
    cells.push_back({{"Patch", cell.patch},
                     {"Feature", cell.feature},
                     {"Stretch", cell.stretch},
                     {"Segment", patch.provenance.segment},
                     {"S0", patch.provenance.s0 * coordinate_scale},
                     {"S1", patch.provenance.s1 * coordinate_scale},
                     {"Model", models.at(patch.model)->name},
                     {"CellLength", cell.cell_length * coordinate_scale},
                     {"OwnedLength", cell.owned_length * coordinate_scale},
                     {"Owners", std::move(owners)}});
  }
  nlohmann::json by_support = nlohmann::json::array();
  for (const auto &support : supports)
  {
    const auto it = ownership.owned_by_support.find(support.patch);
    const auto &spatial = config.patches[support.patch];
    by_support.push_back(
        {{"SpatialPatch", support.patch},
         {"SpatialFeature", spatial.provenance.feature},
         {"SpatialModel", models.at(spatial.model)->name},
         {"FromSignatureBox", support.from_signature_box},
         {"OwnedLength",
          (it == ownership.owned_by_support.end() ? 0.0 : it->second) * coordinate_scale}});
  }
  nlohmann::json vertices = nlohmann::json::array();
  for (const auto &vertex : ownership.vertices)
  {
    const auto &patch = config.patches[vertex.patch];
    nlohmann::json owners = nlohmann::json::array();
    for (const std::size_t owner : vertex.owners)
    {
      const auto &spatial = config.patches[owner];
      owners.push_back({{"SpatialPatch", owner},
                        {"SpatialFeature", spatial.provenance.feature},
                        {"SpatialModel", models.at(spatial.model)->name}});
    }
    vertices.push_back(
        {{"Kind", "Vertex"},
         {"Patch", vertex.patch},
         {"Feature", vertex.feature},
         {"Model", models.at(patch.model)->name},
         {"Origin", std::array<double, 3>{patch.origin[0] * coordinate_scale,
                                          patch.origin[1] * coordinate_scale,
                                          patch.origin[2] * coordinate_scale}},
         {"FaceDistanceOverR", vertex.face_distance_over_r},
         {"ChainEndDistanceOverR", vertex.chain_end_distance_over_r},
         {"ArmOutsideBox", vertex.arm_outside_box},
         {"LostArmLengthOverR", vertex.lost_arm_length_over_r},
         {"Owners", std::move(owners)}});
  }
  return {
      {"Cells", static_cast<int>(ownership.cells.size())},
      {"WhollyOwnedCells", ownership.wholly_owned_cells},
      {"ClippedCells", ownership.clipped_cells},
      {"OwnedLength", ownership.owned_length * coordinate_scale},
      {"Shared",
       {{"Cells", ownership.shared_cells},
        {"Length", ownership.shared_length * coordinate_scale}}},
      {"BySupport", std::move(by_support)},
      {"OwnedCells", std::move(cells)},
      {"Vertices",
       {{"Count", static_cast<int>(ownership.vertices.size())},
        {"Shared", ownership.shared_vertices},
        {"Records", std::move(vertices)},
        {"Rule",
         "decision 282 rule B4 / decision 285 (4): a vertex feature's patch (corner / "
         "junction / endpoint coupon) whose vertex lies on a chain piece END of a "
         "contract-3 coupon's continuation chain (the Signature's Chain context, placed by "
         "the patch frame) inside that coupon's support box is owned by the coupon (weight "
         "0, once; every owner listed); a vertex closer than R to a face has an arm partly "
         "outside the box (ArmOutsideBox; LostArmLengthOverR = R minus the face distance, "
         "in units of R: the part of the corner's R window beyond the first owner's box "
         "that no coupon corrects once the vertex patch has weight 0, R1 final review "
         "MINOR-2); a vertex on another cluster's claims is never on a chain. The vertex "
         "feature's claimed length is in the identification manifest "
         "(Features[].Length), not in the patch"}}},
      {"Rule",
       "decision 236 (2) / 244 (2026-10-02): a translational cell of a stretch that "
       "continues "
       "a claim of a spatial cluster (the Continuation criterion of "
       "TranslationalStretchesInsideSpatialSupport, judged on the whole stretch) is owned "
       "by "
       "that coupon inside its box: the kept part of the cell is the part outside EVERY "
       "box "
       "whose claims the stretch continues (symmetric for a cell on two coupons' "
       "continuations, exact, idempotent), the patch weight and quadrature weight scale by "
       "kept / cell (a cell wholly inside keeps weight 0 and is skipped), the origin moves "
       "to "
       "the kept interval's midpoint with the cell symmetric about it; the provenance "
       "portion [S0, S1) stays, so a portion's quadrature weights sum to 1 - OwnedLength / "
       "portion length. Foreign cells, cells outside the box and the stack-end cells of a "
       "stretch that continues no claim are untouched; a curved cell never continues a "
       "claim "
       "(residual double count on arc-continues-arc); a pair's or stack's SIDE is owned "
       "only where its OWN edge continues the claim (decision 252: the segment branch or "
       "the own-edge abutment of the Continuation criterion), the cells of an unclaimed "
       "side stay. A shared cell's owned length is attributed per coupon by the midpoint "
       "between the two continued claim ends. Lengths in mesh units"}};
}

nlohmann::json DescribeSpatialSupportMarginOverlaps(
    const std::vector<SpatialSupportMarginOverlap> &overlaps,
    const ResponseCorrectionData &config, double coordinate_scale)
{
  std::unordered_map<int, const ResponseModelData *> models;
  for (const auto &model : config.models)
  {
    models.emplace(model.idx, &model);
  }
  auto Scaled = [&](const std::array<double, 3> &point)
  {
    return std::array<double, 3>{point[0] * coordinate_scale, point[1] * coordinate_scale,
                                 point[2] * coordinate_scale};
  };
  nlohmann::json pairs = nlohmann::json::array();
  double margin_over_claims = 0.0, margin_over_margin = 0.0;
  int claim_in_hull = 0;
  for (const auto &overlap : overlaps)
  {
    const auto &first = config.patches[overlap.first_patch];
    const auto &second = config.patches[overlap.second_patch];
    margin_over_claims +=
        overlap.first_margin_over_second_claims + overlap.second_margin_over_first_claims;
    margin_over_margin += overlap.margin_over_margin;
    claim_in_hull += overlap.claim_in_hull ? 1 : 0;
    pairs.push_back(
        {{"FirstPatch", overlap.first_patch},
         {"FirstFeature", first.provenance.feature},
         {"FirstModel", models.at(first.model)->name},
         {"SecondPatch", overlap.second_patch},
         {"SecondFeature", second.provenance.feature},
         {"SecondModel", models.at(second.model)->name},
         {"Overlap",
          {{"Min", Scaled(overlap.overlap_min)}, {"Max", Scaled(overlap.overlap_max)}}},
         {"FirstMarginOverSecondClaims",
          overlap.first_margin_over_second_claims * coordinate_scale},
         {"SecondMarginOverFirstClaims",
          overlap.second_margin_over_first_claims * coordinate_scale},
         {"MarginOverMargin", overlap.margin_over_margin * coordinate_scale},
         {"ClaimInHull", overlap.claim_in_hull}});
  }
  return {
      {"Count", static_cast<int>(overlaps.size())},
      {"ClaimInHull", claim_in_hull},
      {"MarginOverClaimsLength", margin_over_claims * coordinate_scale},
      {"MarginOverMarginLength", margin_over_margin * coordinate_scale},
      {"DoubleCountedLength", (margin_over_claims + margin_over_margin) * coordinate_scale},
      {"Pairs", std::move(pairs)},
      {"Rule",
       "decision 244 (2026-10-02): two spatial cluster supports whose boxes overlap in "
       "their "
       "interiors are recorded when the overlap is margins only: no claim SEGMENT of "
       "either enters the bounding box of the other's claims beyond the tolerance, ends "
       "outside or not (decision 252); they abort otherwise. Each "
       "coupon's twins continue every claim CUT end (an end no other claim of the same "
       "cluster shares) straight to its own box face; the length of those continuations "
       "lying on the other coupon's claims (MarginOverClaims) or on the other coupon's "
       "continuations (MarginOverMargin) is corrected by both coupons: a double count the "
       "placement cannot remove (the dense coupon operator is not clippable; the "
       "translational cell on a shared continuation is removed once by the continuation "
       "ownership). Follow-up: shrinking the continuation from 2R to R (option D) reduces "
       "every margin double count, with library rebuilds. Lengths in mesh units"}};
}

std::string DescribeSpatialSupportMarginOverlapWarning(const nlohmann::json &diagnostics)
{
  std::string lines;
  for (const auto &entry : diagnostics["Pairs"])
  {
    lines += fmt::format(
        "  spatial patches {} ({}) and {} ({}): margin over the other's claims {:.6e} / "
        "{:.6e}, margin over margin {:.6e} mesh units{}\n",
        entry["FirstPatch"].get<std::size_t>() + 1, entry["FirstModel"].get<std::string>(),
        entry["SecondPatch"].get<std::size_t>() + 1,
        entry["SecondModel"].get<std::string>(),
        entry["FirstMarginOverSecondClaims"].get<double>(),
        entry["SecondMarginOverFirstClaims"].get<double>(),
        entry["MarginOverMargin"].get<double>(),
        entry["ClaimInHull"].get<bool>() ? " (a claim inside the other's claims)" : "");
  }
  return fmt::format(
      "{:d} pair(s) of spatial cluster supports overlap in their margins (decision 244, "
      "recorded under Diagnostics.SpatialSupportMarginOverlaps): {:.6e} mesh units of edge "
      "corrected by both coupons\n{}",
      diagnostics["Count"].get<int>(), diagnostics["DoubleCountedLength"].get<double>(),
      lines);
}

std::string DescribeContinuationOwnershipSummary(const nlohmann::json &diagnostics)
{
  std::string lines;
  for (const auto &entry : diagnostics["BySupport"])
  {
    if (entry["OwnedLength"].get<double>() <= 0.0)
    {
      continue;
    }
    lines += fmt::format("  spatial patch {} ({}): {:.6e} mesh units owned\n",
                         entry["SpatialPatch"].get<std::size_t>() + 1,
                         entry["SpatialModel"].get<std::string>(),
                         entry["OwnedLength"].get<double>());
  }
  return fmt::format(
      "Continuation ownership (decision 236 (2)): {:d} translational cell(s) owned by the "
      "spatial coupons whose claims they continue ({:d} wholly, {:d} clipped at a box "
      "face; {:d} shared by two coupons), {:.6e} mesh units removed\n{}",
      diagnostics["Cells"].get<int>(), diagnostics["WhollyOwnedCells"].get<int>(),
      diagnostics["ClippedCells"].get<int>(), diagnostics["Shared"]["Cells"].get<int>(),
      diagnostics["OwnedLength"].get<double>(), lines);
}

DomainBoundaryExclusions FindDomainBoundaryExclusions(
    mfem::ParMesh &mesh, std::vector<ResponsePatchData> &patches,
    const std::function<const std::vector<std::array<double, 3>> *(int model_idx)>
        &basis_points,
    const std::function<bool(int model_idx)> &spatial_basis,
    const std::function<std::string(int model_idx)> &model_name, double coordinate_scale,
    double matching_radius, const std::set<std::size_t> &skipped)
{
  const auto start = std::chrono::steady_clock::now();
  const int dimension = mesh.Dimension();
  MFEM_VERIFY(dimension == 3 && mesh.SpaceDimension() == 3,
              "The domain-boundary exclusion applies to three-dimensional devices!");
  const auto comm = mesh.GetComm();
  DomainBoundaryExclusions result;

  // The tested points of every candidate patch (the same list on every rank).
  const double inset = kSignatureParameterToleranceOverRadius * matching_radius;
  constexpr std::size_t no_reference = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> candidates;
  std::vector<std::size_t> point_offsets = {0};
  // Per candidate, the index of its first conductor reference at the origin section (the
  // metal-edge point, which lies on a mesh face for every correctly placed patch, a metal
  // edge on a chip-outline face included): its absence is a misplaced or mis-scaled
  // coupon, never a domain cut, and fails closed below (decision 260).
  std::vector<std::size_t> origin_references;
  std::vector<std::array<double, 3>> points;
  for (std::size_t patch_idx = 0; patch_idx < patches.size(); patch_idx++)
  {
    const auto &patch = patches[patch_idx];
    const auto *local_points = basis_points(patch.model);
    if (skipped.count(patch_idx) || patch.weight <= 0.0 || !local_points)
    {
      continue;
    }
    origin_references.push_back(patch.conductor_references.empty()
                                    ? no_reference
                                    : points.size() + local_points->size());
    const bool spatial = spatial_basis(patch.model);
    // The cross-sections: the origin and, for a longitudinal cell, both ends moved inward.
    std::vector<double> sections = {0.0};
    const auto &cell = patch.longitudinal_cell;
    if (!spatial && cell[1] > cell[0])
    {
      sections.push_back(std::min(0.0, cell[0] + inset));
      sections.push_back(std::max(0.0, cell[1] - inset));
    }
    for (const double section : sections)
    {
      for (const auto &local : *local_points)
      {
        std::array<double, 3> point{};
        for (int d = 0; d < 3; d++)
        {
          point[d] = patch.origin[d] +
                     (local[0] * patch.axis_u[d] + local[1] * patch.axis_v[d] +
                      (spatial ? local[2] * patch.axis_w[d] : 0.0)) /
                         coordinate_scale +
                     section * patch.axis_w[d];
        }
        points.push_back(point);
      }
      for (const auto &reference : patch.conductor_references)
      {
        std::array<double, 3> point{};
        for (int d = 0; d < 3; d++)
        {
          point[d] = patch.origin[d] + reference[0] * patch.axis_u[d] +
                     reference[1] * patch.axis_v[d] +
                     (reference[2] + section) * patch.axis_w[d];
        }
        points.push_back(point);
      }
    }
    candidates.push_back(patch_idx);
    point_offsets.push_back(points.size());
  }
  result.tested_patches = static_cast<long long int>(candidates.size());
  result.tested_points = static_cast<long long int>(points.size());
  MFEM_VERIFY(points.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "Too many domain-boundary test points for one reduction!");

  // Every rank locates every point in its local mesh with the operator's locator and its
  // tolerances; the found flags are OR-reduced, so the decision is the partition's union.
  ElementPointLocator locator(mesh, dimension);
  double locator_scale = 0.0;
  for (int d = 0; d < dimension; d++)
  {
    locator_scale = std::max({locator_scale, std::abs(locator.GetBounds().min[d]),
                              std::abs(locator.GetBounds().max[d]),
                              locator.GetBounds().max[d] - locator.GetBounds().min[d]});
  }
  Mpi::GlobalMax(1, &locator_scale, comm);
  const double box_tolerance =
      1.0e-11 * locator_scale +
      64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, locator_scale);
  std::vector<unsigned char> found(points.size(), 0);
  {
    std::vector<int> element_candidates;
    int element;
    mfem::IntegrationPoint reference;
    for (std::size_t i = 0; i < points.size(); i++)
    {
      if (locator.Find(points[i], box_tolerance, element, reference, element_candidates))
      {
        found[i] = 1;
      }
    }
  }
  if (!found.empty())
  {
    Mpi::GlobalMax(static_cast<int>(found.size()), found.data(), comm);
  }
  // Fail closed on a metal-edge reference outside the mesh (the found flags are identical
  // on every rank after the reduction, so every rank aborts together).
  auto DescribeMisplacedPatch = [&](std::size_t c, std::size_t point, const char *reason)
  {
    const auto &patch = patches[candidates[c]];
    return fmt::format(
        "Surface-response patch {:d} (0-based, as in the record and the dry run; model {}) "
        "is a misplaced or mis-scaled coupon: {} lies outside the device mesh (tested "
        "point {:d} of {:d} at ({:.9e}, {:.9e}, {:.9e}) mesh units)!",
        candidates[c], model_name(patch.model), reason, point - point_offsets[c],
        point_offsets[c + 1] - point_offsets[c], points[point][0] * coordinate_scale,
        points[point][1] * coordinate_scale, points[point][2] * coordinate_scale);
  };
  for (std::size_t c = 0; c < candidates.size(); c++)
  {
    const std::size_t reference = origin_references[c];
    if (reference != no_reference && !found[reference])
    {
      MFEM_ABORT(DescribeMisplacedPatch(c, reference, "its metal-edge reference"));
    }
  }
  // The distance of every outside point to the nearest element box, over the ranks.
  std::vector<std::size_t> outside;
  for (std::size_t i = 0; i < points.size(); i++)
  {
    if (!found[i])
    {
      outside.push_back(i);
    }
  }
  std::vector<double> distances(outside.size());
  for (std::size_t k = 0; k < outside.size(); k++)
  {
    distances[k] = locator.BoxDistance(points[outside[k]]);
  }
  if (!distances.empty())
  {
    Mpi::GlobalMin(static_cast<int>(distances.size()), distances.data(), comm);
  }

  std::size_t k = 0;
  for (std::size_t c = 0; c < candidates.size(); c++)
  {
    const std::size_t begin = point_offsets[c], end = point_offsets[c + 1];
    DomainBoundaryExclusion exclusion;
    exclusion.patch = candidates[c];
    exclusion.tested_points = static_cast<int>(end - begin);
    exclusion.nearest_distance = mfem::infinity();
    while (k < outside.size() && outside[k] < end)
    {
      exclusion.outside_points++;
      if (distances[k] < exclusion.nearest_distance)
      {
        exclusion.nearest_distance = distances[k];
        exclusion.nearest_outside_point = points[outside[k]];
      }
      k++;
    }
    if (exclusion.outside_points > 0)
    {
      // No tested point inside: not a cut through the coupon but a coupon off the mesh.
      if (exclusion.outside_points == exclusion.tested_points)
      {
        MFEM_ABORT(DescribeMisplacedPatch(c, begin, "every one of its tested points"));
      }
      patches[exclusion.patch].weight = 0.0;
      result.patches.push_back(exclusion);
    }
  }
  MFEM_VERIFY(candidates.empty() || result.patches.size() < candidates.size(),
              "The domain-boundary exclusion leaves no applied surface-response patch: all "
                  << candidates.size()
                  << " tested patches have placed coupon points outside the device mesh "
                     "(misplaced or mis-scaled coupons)!");
  result.wall_time =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  Mpi::GlobalMax(1, &result.wall_time, comm);
  return result;
}

nlohmann::json DescribeDomainBoundaryExclusions(const DomainBoundaryExclusions &exclusions,
                                                const ResponseCorrectionData &config,
                                                double coordinate_scale)
{
  std::unordered_map<int, const ResponseModelData *> models;
  for (const auto &model : config.models)
  {
    models.emplace(model.idx, &model);
  }
  nlohmann::json entries = nlohmann::json::array();
  std::map<int, std::pair<int, double>> by_feature;
  double cell_length = 0.0, portion_length = 0.0;
  for (const auto &exclusion : exclusions.patches)
  {
    const auto &patch = config.patches[exclusion.patch];
    const auto &provenance = patch.provenance;
    const double cell =
        (patch.longitudinal_cell[1] - patch.longitudinal_cell[0]) * coordinate_scale;
    const double portion =
        provenance.segment >= 0 ? (provenance.s1 - provenance.s0) * coordinate_scale : 0.0;
    cell_length += cell;
    portion_length += portion;
    auto &feature = by_feature[provenance.feature];
    feature.first++;
    feature.second += cell;
    entries.push_back(
        {{"Patch", exclusion.patch},
         {"Feature", provenance.feature},
         {"Model", models.at(patch.model)->name},
         {"Topology", models.at(patch.model)->topology},
         {"Origin",
          {patch.origin[0] * coordinate_scale, patch.origin[1] * coordinate_scale,
           patch.origin[2] * coordinate_scale}},
         {"Cell",
          {patch.longitudinal_cell[0] * coordinate_scale,
           patch.longitudinal_cell[1] * coordinate_scale}},
         {"CellLength", cell},
         {"Segment", provenance.segment},
         {"S0", provenance.s0 * coordinate_scale},
         {"S1", provenance.s1 * coordinate_scale},
         {"PortionLength", portion},
         {"TestedPoints", exclusion.tested_points},
         {"OutsidePoints", exclusion.outside_points},
         {"NearestOutsidePoint",
          {exclusion.nearest_outside_point[0] * coordinate_scale,
           exclusion.nearest_outside_point[1] * coordinate_scale,
           exclusion.nearest_outside_point[2] * coordinate_scale}},
         {"NearestDistance", exclusion.nearest_distance * coordinate_scale}});
  }
  nlohmann::json features = nlohmann::json::array();
  for (const auto &[feature, count_length] : by_feature)
  {
    features.push_back({{"Feature", feature},
                        {"Patches", count_length.first},
                        {"CellLength", count_length.second}});
  }
  return {
      {"Count", static_cast<int>(exclusions.patches.size())},
      {"Features", static_cast<int>(by_feature.size())},
      {"CellLength", cell_length},
      {"PortionLength", portion_length},
      {"TestedPatches", exclusions.tested_patches},
      {"TestedPoints", exclusions.tested_points},
      {"WallTime", exclusions.wall_time},
      {"CellEndInsetOverRadius", kSignatureParameterToleranceOverRadius},
      {"ByFeature", std::move(features)},
      {"Patches", std::move(entries)},
      {"Rule",
       "decision 258 (2026-10-02): a library-placed patch any of whose placed coupon "
       "points "
       "(its model's basis points and conductor references at the origin cross-section "
       "and, for a translational cell, at both cell ends moved CellEndInsetOverRadius x R "
       "inward) lies outside the device mesh is not applied (weight 0) and recorded here: "
       "the coupon's trace coupling is undefined beyond the device domain, which happens "
       "where a metal edge meets an artificial domain cut (a window cut, a chip outline) "
       "within ~R |sin theta| of it. One containment test (the operator's element point "
       "locator on every rank, OR-reduced) decides for the preflight and the solve, "
       "independent of the rank count and locator path; the applied patches' points still "
       "fail closed when not located. CellLength is the length left uncorrected (the "
       "inventory total next to Missing; not a library gap); PortionLength the sum of the "
       "excluded patches' portions (information). NearestDistance is to the nearest "
       "element bounding box. Lengths in mesh units"}};
}

std::string DescribeDomainBoundaryExclusionSummary(const nlohmann::json &diagnostics)
{
  std::string lines;
  for (const auto &entry : diagnostics["Patches"])
  {
    const auto &point = entry["NearestOutsidePoint"];
    lines += fmt::format(
        "  patch {} (0-based {}; feature {}, {}): cell {:.6e} mesh units of portion "
        "[{:.6e}, {:.6e}] on segment {}, {} / {} points outside, nearest ({:.6e}, {:.6e}, "
        "{:.6e}) at {:.3e}\n",
        entry["Patch"].get<std::size_t>() + 1, entry["Patch"].get<std::size_t>(),
        entry["Feature"].get<int>(), entry["Model"].get<std::string>(),
        entry["CellLength"].get<double>(), entry["S0"].get<double>(),
        entry["S1"].get<double>(), entry["Segment"].get<int>(),
        entry["OutsidePoints"].get<int>(), entry["TestedPoints"].get<int>(),
        point[0].get<double>(), point[1].get<double>(), point[2].get<double>(),
        entry["NearestDistance"].get<double>());
  }
  return fmt::format(
      "DomainBoundary exclusion (decision 258): {:d} patch(es) on {:d} feature(s) with "
      "placed coupon points outside the device mesh not applied, {:.6e} mesh units of "
      "cells left uncorrected (portions {:.6e}); {:d} patches / {:d} points tested in "
      "{:.3f} s (patch indices 1-based like the ownership summaries, 0-based in the "
      "record and the dry run)\n{}",
      diagnostics["Count"].get<int>(), diagnostics["Features"].get<int>(),
      diagnostics["CellLength"].get<double>(), diagnostics["PortionLength"].get<double>(),
      diagnostics["TestedPatches"].get<long long int>(),
      diagnostics["TestedPoints"].get<long long int>(),
      diagnostics["WallTime"].get<double>(), lines);
}

void SurfaceResponseOperator::ConfigureConductorConsistencyProbes(ResponseModel &model)
{
  model.plane_conductor_vertices.clear();
  model.off_plane_conductor_vertices.clear();
  model.adjacent_free_vertices.clear();
  const int vertex_count = static_cast<int>(model.mortar_vertices.size());
  // The process plane is w = 0 of the coupon frame (the metal bottom, where the device's
  // metal sheet lies); the tolerance is relative to the trace mesh's extent normal to it.
  double extent = 0.0;
  for (const auto &vertex : model.mortar_vertices)
  {
    extent = std::max(extent, std::abs(vertex.point[2]));
  }
  const double plane_tolerance = 1.0e-6 * extent;
  std::set<int> conductors, plane_conductors;
  for (int v = 0; v < vertex_count; v++)
  {
    const auto &vertex = model.mortar_vertices[v];
    if (vertex.conductor <= 0)
    {
      continue;
    }
    conductors.insert(vertex.conductor);
    if (std::abs(vertex.point[2]) <= plane_tolerance)
    {
      model.plane_conductor_vertices.push_back(v);
      plane_conductors.insert(vertex.conductor);
    }
    else
    {
      model.off_plane_conductor_vertices.push_back(v);
    }
  }
  for (const int conductor : conductors)
  {
    MFEM_VERIFY(plane_conductors.count(conductor),
                "Conductor-consistency gate (decision 277): conductor "
                    << conductor << " of the trace mesh of response model \"" << model.name
                    << "\" has no vertex on the process plane (w = 0), so its metal "
                       "cross-section on the box faces cannot be probed!");
  }
  if (model.plane_conductor_vertices.empty())
  {
    return;
  }
  // The adjacent free knots: sharing a trace-triangle edge with a plane conductor vertex,
  // in its column (the same in-plane position), a basis knot (no slave, no conductor) and
  // next to vertices of one conductor only.
  std::vector<std::set<int>> neighbors(vertex_count);
  for (const auto &triangle : model.mortar_triangles)
  {
    for (const int a : triangle.vertices)
    {
      for (const int b : triangle.vertices)
      {
        if (a != b)
        {
          neighbors[a].insert(b);
        }
      }
    }
  }
  std::map<int, std::set<int>> free_conductors;
  std::map<int, int> free_plane_vertex;
  for (const int p : model.plane_conductor_vertices)
  {
    const auto &plane = model.mortar_vertices[p];
    for (const int n : neighbors[p])
    {
      const auto &free = model.mortar_vertices[n];
      if (free.conductor > 0 || free.basis < 0 || free.second_basis >= 0 ||
          std::abs(free.point[0] - plane.point[0]) > plane_tolerance ||
          std::abs(free.point[1] - plane.point[1]) > plane_tolerance)
      {
        continue;
      }
      free_conductors[n].insert(plane.conductor);
      free_plane_vertex.try_emplace(n, p);
    }
  }
  for (const auto &[n, conductor_set] : free_conductors)
  {
    if (conductor_set.size() == 1)
    {
      model.adjacent_free_vertices.emplace_back(n, free_plane_vertex.at(n));
    }
  }
}

bool SurfaceResponseOperator::HasConductorConsistencyProbes() const
{
  int local = 0;
  for (const auto &patch : patches)
  {
    if (!models[patch.model].plane_conductor_vertices.empty())
    {
      local = 1;
      break;
    }
  }
  Mpi::GlobalMax(1, &local, fespace.GetComm());
  return local > 0;
}

nlohmann::json &SurfaceResponseOperator::ConductorConsistencyDiagnostics()
{
  auto &diagnostics = ownership_diagnostics["ConductorConsistency"];
  if (!diagnostics.contains("Rule"))
  {
    diagnostics["Tolerance"] = kConductorConsistencyTolerance;
    diagnostics["AmplitudeFloor"] = kConductorConsistencyAmplitudeFloor;
    diagnostics["Count"] = 0;
    diagnostics["ClaimLength"] = 0.0;
    diagnostics["CellLength"] = 0.0;
    diagnostics["TestedPatches"] = 0;
    diagnostics["ExcludedPatches"] = nlohmann::json::array();
    diagnostics["Records"] = nlohmann::json::array();
    diagnostics["UnprobedModels"] = nlohmann::json::array();
    diagnostics["Rule"] =
        "decision 277 (2026-10-03): solve-time gate on every applied spatial "
        "surface-mortar patch, per excitation. MaxRatio = max over the trace mesh's "
        "conductor vertices on the process plane (w = 0: the coupon's metal bottom on the "
        "device's metal sheet) of |V_device(knot) - V_device(conductor reference)| / "
        "Normalization, Normalization = max(Amplitude, AmplitudeFloor x "
        "ExcitationPotential), Amplitude = max |trace coefficient| of the patch incl. its "
        "conductor states (= the state for a two-conductor coupon without overshoots; "
        "State alongside), ExcitationPotential = max |V| of the excitation (its largest "
        "terminal potential); FloorApplied marks the patches whose amplitude lies under "
        "the floor (noise over noise: never excluded, at most AmplitudeFloor^2 of a "
        "unit-amplitude patch's energy). Real metal reads ~0 (the knot lies on a "
        "Dirichlet surface of the device), coupon metal where the device has gap reads "
        "the gap potential. MaxRatio > Tolerance excludes the patch (weight 0 from this "
        "excitation on, like a DomainBoundary cell; the patches CSV is rewritten); "
        "ClaimLength is the excluded spatial clusters' claimed portion length (left "
        "uncorrected: the B1 consumer must add it to Missing and DomainBoundary), "
        "CellLength the excluded longitudinal cells (0 for spatial patches). Information "
        "only: OffPlaneMaxRatio (the metal top rows, in the device gap for a thin device) "
        "and AdjacentMaxRatio (the free knots adjacent to a plane conductor vertex in its "
        "column: their mortar coefficient vs the conductor; near-edge field x 50 nm, up "
        "to 0.27 of the amplitude on real metal). UnprobedModels lists the spatial models "
        "the gate cannot test (applied collocated, or a trace mesh without any conductor "
        "vertex: a finite-impedance coupon or a ring-path model without ZeroTraceIndices), "
        "recorded untestable, not excluded. Not evaluated by the preflight (no device "
        "trace). Lengths and coordinates in mesh-file units (the device coordinates: the "
        "mesh coordinates x the mesh coordinate scale), like the DomainBoundary record";
  }
  return diagnostics;
}

std::vector<SurfaceResponseOperator::ConductorConsistencyRecord>
SurfaceResponseOperator::ApplyConductorConsistencyGate(const Vector &x, int source)
{
  MFEM_VERIFY(!maxwell, "The conductor-consistency gate applies to electrostatic response "
                        "correction only!");
  // The trace walk also samples the probe points (the device potential at the conductor
  // vertices) into `correction`.
  ApplyTrace(x, trace);
  // The excitation's largest potential (the largest terminal potential: the maximum
  // principle) sets the normalization floor of the ratio.
  double excitation_potential = x.Size() > 0 ? x.Normlinf() : 0.0;
  Mpi::GlobalMax(1, &excitation_potential, fespace.GetComm());
  const double amplitude_floor = kConductorConsistencyAmplitudeFloor * excitation_potential;
  constexpr int record_size = 22;
  std::vector<double> local_records;
  for (auto &patch : patches)
  {
    const auto &model = models[patch.model];
    if (patch.probe_point_count == 0 || patch.weight <= 0.0)
    {
      continue;
    }
    MFEM_ASSERT(patch.probe_point_count ==
                    static_cast<int>(model.plane_conductor_vertices.size() +
                                     model.off_plane_conductor_vertices.size()),
                "Inconsistent conductor-consistency probe count!");
    const int reference_offset = patch.point_offset + patch.point_count -
                                 patch.probe_point_count - model.conductor_state_count - 1;
    const int probe_offset =
        patch.point_offset + patch.point_count - patch.probe_point_count;
    auto ConductorValue = [&](int conductor)
    { return correction(reference_offset + conductor - 1); };
    ConductorConsistencyRecord record;
    record.patch = patch.global_index;
    record.model = model.idx;
    record.source = source;
    record.plane_knots = static_cast<int>(model.plane_conductor_vertices.size());
    record.off_plane_knots = static_cast<int>(model.off_plane_conductor_vertices.size());
    record.adjacent_knots = static_cast<int>(model.adjacent_free_vertices.size());
    record.claim_length = patch.claim_length;
    record.cell_length = patch.cell_length;
    for (int i = 0; i < model.basis_size; i++)
    {
      record.amplitude =
          std::max(record.amplitude, std::abs(trace(patch.trace_offset + i)));
    }
    for (int state = 0; state < model.conductor_state_count; state++)
    {
      record.state = std::max(
          record.state, std::abs(trace(patch.trace_offset + model.contour_size + state)));
    }
    record.excitation_potential = excitation_potential;
    record.normalization = std::max(record.amplitude, amplitude_floor);
    auto Ratio = [&](double deviation)
    {
      if (record.normalization > 0.0)
      {
        return std::abs(deviation) / record.normalization;
      }
      return deviation == 0.0 ? 0.0 : mfem::infinity();
    };
    int probe = 0;
    for (const int v : model.plane_conductor_vertices)
    {
      const auto &vertex = model.mortar_vertices[v];
      const double deviation =
          correction(probe_offset + probe) - ConductorValue(vertex.conductor);
      if (probe == 0 || std::abs(deviation) > record.max_deviation)
      {
        record.max_deviation = std::abs(deviation);
        record.max_ratio = Ratio(deviation);
        record.worst_vertex = v;
        record.worst_conductor = vertex.conductor;
        record.worst_point = patch.probe_points[probe];
      }
      probe++;
    }
    for (const int v : model.off_plane_conductor_vertices)
    {
      const auto &vertex = model.mortar_vertices[v];
      const double deviation =
          correction(probe_offset + probe) - ConductorValue(vertex.conductor);
      record.off_plane_max_ratio = std::max(record.off_plane_max_ratio, Ratio(deviation));
      probe++;
    }
    for (const auto &[free, plane] : model.adjacent_free_vertices)
    {
      const int conductor = model.mortar_vertices[plane].conductor;
      // Trace coefficients are relative to conductor 1; conductor c >= 2 is its state.
      const double conductor_trace =
          conductor == 1 ? 0.0
                         : trace(patch.trace_offset + model.contour_size + conductor - 2);
      const double deviation =
          trace(patch.trace_offset + model.mortar_vertices[free].basis) - conductor_trace;
      record.adjacent_max_ratio = std::max(record.adjacent_max_ratio, Ratio(deviation));
    }
    record.excluded = record.max_ratio > kConductorConsistencyTolerance;
    if (record.excluded)
    {
      patch.weight = 0.0;
    }
    local_records.insert(local_records.end(), {static_cast<double>(record.patch),
                                               static_cast<double>(record.model),
                                               static_cast<double>(record.source),
                                               static_cast<double>(record.plane_knots),
                                               static_cast<double>(record.off_plane_knots),
                                               static_cast<double>(record.adjacent_knots),
                                               record.amplitude,
                                               record.state,
                                               record.excitation_potential,
                                               record.normalization,
                                               record.max_deviation,
                                               record.max_ratio,
                                               static_cast<double>(record.worst_vertex),
                                               static_cast<double>(record.worst_conductor),
                                               record.worst_point[0],
                                               record.worst_point[1],
                                               record.worst_point[2],
                                               record.off_plane_max_ratio,
                                               record.adjacent_max_ratio,
                                               record.claim_length,
                                               record.cell_length,
                                               record.excluded ? 1.0 : 0.0});
  }
  MFEM_VERIFY(local_records.size() <=
                  static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "Local conductor-consistency data exceeds the MPI count limit!");
  const int local_value_count = static_cast<int>(local_records.size());
  std::vector<int> value_counts(Mpi::Size(fespace.GetComm()));
  Mpi::Allgather(1, &local_value_count, value_counts.data(), fespace.GetComm());
  std::vector<int> value_offsets(value_counts.size());
  int total_values = 0;
  for (std::size_t rank = 0; rank < value_counts.size(); rank++)
  {
    value_offsets[rank] = total_values;
    total_values += value_counts[rank];
  }
  std::vector<double> values(total_values);
  Mpi::Allgatherv(local_value_count, local_records.data(), values.data(),
                  value_counts.data(), value_offsets.data(), fespace.GetComm());
  MFEM_VERIFY(values.size() % record_size == 0, "Truncated conductor-consistency record!");
  std::vector<ConductorConsistencyRecord> records;
  for (std::size_t offset = 0; offset < values.size(); offset += record_size)
  {
    const double *r = values.data() + offset;
    ConductorConsistencyRecord record;
    record.patch = static_cast<int>(std::llround(r[0]));
    record.model = static_cast<int>(std::llround(r[1]));
    record.source = static_cast<int>(std::llround(r[2]));
    record.plane_knots = static_cast<int>(std::llround(r[3]));
    record.off_plane_knots = static_cast<int>(std::llround(r[4]));
    record.adjacent_knots = static_cast<int>(std::llround(r[5]));
    record.amplitude = r[6];
    record.state = r[7];
    record.excitation_potential = r[8];
    record.normalization = r[9];
    record.max_deviation = r[10];
    record.max_ratio = r[11];
    record.worst_vertex = static_cast<int>(std::llround(r[12]));
    record.worst_conductor = static_cast<int>(std::llround(r[13]));
    record.worst_point = {r[14], r[15], r[16]};
    record.off_plane_max_ratio = r[17];
    record.adjacent_max_ratio = r[18];
    record.claim_length = r[19];
    record.cell_length = r[20];
    record.excluded = r[21] > 0.5;
    records.push_back(record);
  }
  std::sort(records.begin(), records.end(),
            [](const auto &a, const auto &b) { return a.patch < b.patch; });

  // The record (every rank holds the same): one entry per tested patch and excitation,
  // the excluded patches' claim / cell lengths totalled next to the DomainBoundary record.
  std::map<int, std::string> model_names = GetModelNames();
  auto &diagnostics = ConductorConsistencyDiagnostics();
  int excluded_count = 0;
  double excluded_claims = 0.0, excluded_cells = 0.0, max_ratio = 0.0;
  int max_patch = -1;
  for (const auto &record : records)
  {
    nlohmann::json entry = {{"Patch", record.patch},
                            {"Model", model_names.at(record.model)},
                            {"ModelIndex", record.model},
                            {"Source", record.source},
                            {"PlaneKnots", record.plane_knots},
                            {"OffPlaneKnots", record.off_plane_knots},
                            {"AdjacentKnots", record.adjacent_knots},
                            {"Amplitude", record.amplitude},
                            {"State", record.state},
                            {"ExcitationPotential", record.excitation_potential},
                            {"Normalization", record.normalization},
                            {"FloorApplied", record.normalization > record.amplitude},
                            {"MaxDeviation", record.max_deviation},
                            {"MaxRatio", record.max_ratio},
                            {"MaxRatioOverState", record.state > 0.0
                                                      ? record.max_deviation / record.state
                                                      : 0.0},
                            {"WorstVertex", record.worst_vertex},
                            {"WorstConductor", record.worst_conductor},
                            {"WorstPoint",
                             {record.worst_point[0] * mesh_coordinate_scale,
                              record.worst_point[1] * mesh_coordinate_scale,
                              record.worst_point[2] * mesh_coordinate_scale}},
                            {"OffPlaneMaxRatio", record.off_plane_max_ratio},
                            {"AdjacentMaxRatio", record.adjacent_max_ratio},
                            {"ClaimLength", record.claim_length * mesh_coordinate_scale},
                            {"CellLength", record.cell_length * mesh_coordinate_scale},
                            {"Excluded", record.excluded}};
    diagnostics["Records"].push_back(entry);
    if (record.excluded)
    {
      diagnostics["ExcludedPatches"].push_back(entry);
      excluded_count++;
      excluded_claims += record.claim_length * mesh_coordinate_scale;
      excluded_cells += record.cell_length * mesh_coordinate_scale;
      // The root's patch table (the patches CSV) follows the exclusion.
      for (auto &assignment : patch_assignments)
      {
        if (assignment.global_index == record.patch)
        {
          assignment.weight = 0.0;
        }
      }
    }
    if (record.max_ratio > max_ratio || max_patch < 0)
    {
      max_ratio = record.max_ratio;
      max_patch = record.patch;
    }
  }
  diagnostics["Count"] = diagnostics["Count"].get<int>() + excluded_count;
  diagnostics["ClaimLength"] = diagnostics["ClaimLength"].get<double>() + excluded_claims;
  diagnostics["CellLength"] = diagnostics["CellLength"].get<double>() + excluded_cells;
  diagnostics["TestedPatches"] =
      diagnostics["TestedPatches"].get<int>() + static_cast<int>(records.size());
  if (!records.empty())
  {
    Mpi::Print(" Conductor-consistency gate (decision 277), excitation {:d}: {:d} spatial "
               "patch(es) probed, max ratio {:.3e} (patch {:d}, 0-based), {:d} excluded "
               "(tolerance {:.3e}; amplitude floor {:.3e} x the excitation potential "
               "{:.6e}; claims {:.6e}, cells {:.6e} mesh-file units left uncorrected)\n",
               source, static_cast<int>(records.size()), max_ratio, max_patch,
               excluded_count, kConductorConsistencyTolerance,
               kConductorConsistencyAmplitudeFloor, excitation_potential, excluded_claims,
               excluded_cells);
    for (const auto &record : records)
    {
      if (record.excluded)
      {
        Mpi::Print(
            "  patch {:d} (0-based; {}): |dV| {:.6e} = {:.3e} of the normalization "
            "{:.6e} (amplitude {:.6e}, state {:.6e}) at trace vertex {:d} of conductor "
            "{:d}, point ({:.6e}, {:.6e}, {:.6e}) mesh-file units; claims {:.6e} mesh-file "
            "units EXCLUDED\n",
            record.patch, model_names.at(record.model), record.max_deviation,
            record.max_ratio, record.normalization, record.amplitude, record.state,
            record.worst_vertex, record.worst_conductor,
            record.worst_point[0] * mesh_coordinate_scale,
            record.worst_point[1] * mesh_coordinate_scale,
            record.worst_point[2] * mesh_coordinate_scale,
            record.claim_length * mesh_coordinate_scale);
      }
    }
  }
  return records;
}

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
  const double coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
  mesh_coordinate_scale = coordinate_scale;
  auto &response_mesh = const_cast<mfem::ParMesh &>(fespace.GetParMesh());
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
      const auto geometry_cache = ResponseGeometryCachePath();
      const char *write_cache_environment =
          std::getenv("PALACE_RESPONSE_GEOMETRY_CACHE_WRITE");
      const bool write_geometry_cache = write_cache_environment &&
                                        std::string_view(write_cache_environment) != "0" &&
                                        !std::string_view(write_cache_environment).empty();
      if (geometry_cache && !write_geometry_cache)
      {
        MFEM_VERIFY(std::filesystem::is_regular_file(*geometry_cache),
                    "PALACE_RESPONSE_GEOMETRY_CACHE does not name an existing cache. "
                    "Set PALACE_RESPONSE_GEOMETRY_CACHE_WRITE=1 in a serial run to "
                    "create it first!");
        automatic_config = ReadResponseGeometryCache(*geometry_cache, *request);
        automatic_statistics = {{"GeometryCache", geometry_cache->string()}};
        Mpi::Print("Loaded response-geometry cache: {}\n", geometry_cache->string());
      }
      else
      {
        AutomaticResponseStatistics statistics;
        automatic_config =
            BuildAutomaticResponseData(iodata, laplace_op, *request, &statistics);
        automatic_statistics = BuildAutomaticStatistics(fespace.GetComm(), statistics);
        if (geometry_cache)
        {
          MFEM_VERIFY(write_geometry_cache, "Response-geometry cache generation requires "
                                            "PALACE_RESPONSE_GEOMETRY_CACHE_WRITE=1!");
          if (Mpi::Root(fespace.GetComm()))
          {
            WriteResponseGeometryCache(*geometry_cache, *automatic_config);
          }
          Mpi::Barrier(fespace.GetComm());
          Mpi::Print("Wrote response-geometry cache: {}\n", geometry_cache->string());
        }
      }
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
  auto GetTranslationalDomainCorrectionMode = [&]()
  {
    using ConfigMode = ResponseCorrectionData::TranslationalDomainCorrection;
    switch (config->translational_domain_correction)
    {
      case ConfigMode::DISABLED:
        return DomainCorrectionMode::DISABLED;
      case ConfigMode::FIXED_TRACE:
        return DomainCorrectionMode::FIXED_TRACE;
      case ConfigMode::FIXED_FLUX:
        return DomainCorrectionMode::FIXED_FLUX;
    }
    MFEM_ABORT("Unknown translational response domain-correction mode!");
    return DomainCorrectionMode::FIXED_TRACE;
  };
  const double target_matching_radius =
      TargetInterfaceMatchingRadius(iodata, config->target_interfaces);
  for (const auto &model_config : config->models)
  {
    MFEM_VERIFY(model_config.idx > 0 &&
                    model_indices.find(model_config.idx) == model_indices.end(),
                "Response-correction model indices must be positive and unique!");
    auto points = ModelBasisPoints(model_config);
    ResponseModel model;
    model.idx = model_config.idx;
    model.name =
        model_config.name.empty() ? fmt::format("model-{}", model.idx) : model_config.name;
    model.topology = model_config.topology.empty() ? "Explicit" : model_config.topology;
    if (IsTranslationalTopology(model.topology))
    {
      model.domain_correction_mode = GetTranslationalDomainCorrectionMode();
      model.surface_mortar =
          config->trace_coupling == ResponseCorrectionData::TraceCoupling::SURFACE_MORTAR;
    }
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
    // Cap-interior hats (trailing basis points on no contour) are trace coefficients only
    // through the explicit trace mesh: collocated or surface-mortar spatial models.
    model.interior_trace_count = model_config.interior_trace_count;
    MFEM_VERIFY(model.interior_trace_count >= 0 &&
                    model.interior_trace_count < model.contour_size &&
                    (model.interior_trace_count == 0 ||
                     (model.spatial_basis && HasExplicitTraceMesh(model_config))),
                "InteriorTraceCount requires a spatial response model with an explicit "
                "TraceMesh and fewer interior points than BasisPoints!");
    const int contour_point_count = model.contour_size - model.interior_trace_count;
    if (model.contour_groups.empty() && model.open_contour_paths.empty())
    {
      model.contour_groups.push_back(contour_point_count);
    }
    if (!model.contour_groups.empty())
    {
      MFEM_VERIFY(std::accumulate(model.contour_groups.begin(), model.contour_groups.end(),
                                  0) == contour_point_count,
                  "Response-correction ContourGroups do not partition the contour "
                  "BasisPoints (BasisPoints minus the trailing InteriorTraceCount)!");
      MFEM_VERIFY(std::all_of(model.zero_trace_indices.begin(),
                              model.zero_trace_indices.end(), [&](int index)
                              { return index >= 0 && index < contour_point_count; }),
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
          MFEM_VERIFY(index >= 0 && index < contour_point_count && !assigned[index],
                      "Response-correction OpenContourPaths contain an invalid or "
                      "duplicate BasisPoints index!");
          assigned[index] = true;
        }
      }
      MFEM_VERIFY(std::all_of(assigned.begin(), assigned.begin() + contour_point_count,
                              [](bool value) { return value; }),
                  "Response-correction OpenContourPaths do not partition the contour "
                  "BasisPoints (BasisPoints minus the trailing InteriorTraceCount)!");
    }
    MFEM_VERIFY(config->trace_coupling !=
                        ResponseCorrectionData::TraceCoupling::SURFACE_MORTAR ||
                    !model.spatial_basis || model.open_contour_paths.empty() ||
                    HasExplicitTraceMesh(model_config),
                "SurfaceMortar requires an explicit TraceMesh for a spatial model with "
                "OpenContourPaths; regenerate or augment the process library!");
    if (config->trace_coupling == ResponseCorrectionData::TraceCoupling::SURFACE_MORTAR &&
        model.spatial_basis &&
        (model.open_contour_paths.empty() || HasExplicitTraceMesh(model_config)))
    {
      model.surface_mortar = true;
      model.spatial_mortar = true;
    }
    if (model.surface_mortar)
    {
      mfem::DenseMatrix mass(model.contour_size);
      mass = 0.0;
      if (model.spatial_mortar)
      {
        model.mortar_constant_load.SetSize(model.contour_size);
        model.mortar_constant_load = 0.0;
        model.mortar_conductor_loads.resize(model.conductor_state_count);
        for (auto &load : model.mortar_conductor_loads)
        {
          load.SetSize(model.contour_size);
          load = 0.0;
        }
        std::vector<bool> represented_basis(model.contour_size, false);
        const bool explicit_trace_mesh = HasExplicitTraceMesh(model_config);
        TraceMeshData trace_mesh;
        if (explicit_trace_mesh)
        {
          trace_mesh = ModelTraceMesh(model_config);
          model.mortar_vertices.reserve(trace_mesh.vertices.size());
          for (const auto &vertex : trace_mesh.vertices)
          {
            MFEM_VERIFY(
                vertex.basis <= model.contour_size &&
                    vertex.conductor <= model.conductor_state_count + 1 &&
                    (vertex.basis > 0 || vertex.conductor > 0 || vertex.parent_a > 0),
                "A response trace vertex has invalid basis/conductor ownership!");
            if (vertex.parent_a > 0)
            {
              // A slave vertex (corner-family trace basis rule): the trace is the linear
              // interpolation between its two parent knots.
              MFEM_VERIFY(vertex.parent_a <= model.contour_size &&
                              vertex.parent_b <= model.contour_size,
                          "A slave response trace vertex names a parent outside the "
                          "contour basis!");
              model.mortar_vertices.push_back({vertex.point, vertex.parent_a - 1, 0,
                                               vertex.parent_b - 1, vertex.weight_a});
              continue;
            }
            const int basis = vertex.basis - 1;
            if (basis >= 0)
            {
              MFEM_VERIFY(
                  !represented_basis[basis],
                  "A response TraceMesh maps one basis coefficient more than once!");
              represented_basis[basis] = true;
            }
            model.mortar_vertices.push_back({vertex.point, basis, vertex.conductor});
          }
        }
        else
        {
          model.mortar_vertices.reserve(points.size());
          std::set<int> constrained(model.zero_trace_indices.begin(),
                                    model.zero_trace_indices.end());
          for (int i = 0; i < model.contour_size; i++)
          {
            model.mortar_vertices.push_back({points[i], i, constrained.count(i) ? 1 : 0});
            represented_basis[i] = true;
          }
        }
        MFEM_VERIFY(std::all_of(represented_basis.begin(), represented_basis.end(),
                                [](bool represented) { return represented; }),
                    "A response TraceMesh must represent every coupon basis coefficient!");

        auto AddTriangle = [&](int first, int second, int third)
        {
          const auto &a = model.mortar_vertices[first].point;
          const auto &b = model.mortar_vertices[second].point;
          const auto &c = model.mortar_vertices[third].point;
          Point3D ab{}, ac{};
          double maximum_edge_squared = 0.0;
          for (int d = 0; d < 3; d++)
          {
            ab[d] = b[d] - a[d];
            ac[d] = c[d] - a[d];
          }
          for (const auto edge :
               {std::make_pair(first, second), std::make_pair(second, third),
                std::make_pair(third, first)})
          {
            double length_squared = 0.0;
            for (int d = 0; d < 3; d++)
            {
              const double delta = model.mortar_vertices[edge.first].point[d] -
                                   model.mortar_vertices[edge.second].point[d];
              length_squared += delta * delta;
            }
            maximum_edge_squared = std::max(maximum_edge_squared, length_squared);
          }
          const double area = 0.5 * Norm(Cross(ab, ac));
          const double scale = std::max(1.0, maximum_edge_squared);
          MFEM_VERIFY(area > 1.0e-14 * scale || !explicit_trace_mesh,
                      "A response TraceMesh contains a degenerate triangle!");
          if (area <= 1.0e-14 * scale)
          {
            return;
          }
          model.mortar_triangles.push_back(
              {{first, second, third}, area, std::sqrt(maximum_edge_squared)});
          // P1 mass on the triangle: the hat of knot k takes the value w_k(p) at vertex p
          // (1 at its own knot, the interpolation weight at a slave vertex, 0 elsewhere),
          // so int phi_i phi_j = area sum_{p, q} w_i(p) w_j(q) (p == q ? 1/6 : 1/12).
          constexpr double diagonal = 1.0 / 6.0;
          constexpr double off_diagonal = 1.0 / 12.0;
          for (int local_i = 0; local_i < 3; local_i++)
          {
            const auto &vertex_i =
                model.mortar_vertices[model.mortar_triangles.back().vertices[local_i]];
            vertex_i.ForEachBasis(
                [&](int basis_i, double weight_i)
                {
                  model.mortar_constant_load[basis_i] += weight_i * area / 3.0;
                  for (int local_j = 0; local_j < 3; local_j++)
                  {
                    const auto &vertex_j =
                        model.mortar_vertices[model.mortar_triangles.back()
                                                  .vertices[local_j]];
                    const double entry =
                        weight_i * area * (local_i == local_j ? diagonal : off_diagonal);
                    if (vertex_j.basis >= 0)
                    {
                      vertex_j.ForEachBasis(
                          [&](int basis_j, double weight_j)
                          { mass(basis_i, basis_j) += weight_j * entry; });
                    }
                    else if (vertex_j.conductor > 1)
                    {
                      model.mortar_conductor_loads[vertex_j.conductor - 2][basis_i] +=
                          entry;
                    }
                  }
                });
          }
        };
        if (explicit_trace_mesh)
        {
          for (const auto &triangle : trace_mesh.triangles)
          {
            AddTriangle(triangle[0], triangle[1], triangle[2]);
          }
        }
        else
        {
          MFEM_VERIFY(model.contour_groups.size() >= 2,
                      "A spatial surface mortar requires at least two contour rings!");
          const int ring_size = model.contour_groups.front();
          MFEM_VERIFY(ring_size >= 3 &&
                          std::all_of(model.contour_groups.begin(),
                                      model.contour_groups.end(), [ring_size](int count)
                                      { return count == ring_size; }),
                      "A spatial surface mortar requires equal nontrivial contour rings!");
          for (std::size_t ring = 0; ring + 1 < model.contour_groups.size(); ring++)
          {
            const int first_offset = static_cast<int>(ring) * ring_size;
            const int second_offset = first_offset + ring_size;
            for (int i = 0; i < ring_size; i++)
            {
              const int next = (i + 1) % ring_size;
              AddTriangle(first_offset + i, first_offset + next, second_offset + next);
              AddTriangle(first_offset + i, second_offset + next, second_offset + i);
            }
          }
          const int last_offset =
              (static_cast<int>(model.contour_groups.size()) - 1) * ring_size;
          for (int i = 1; i + 1 < ring_size; i++)
          {
            AddTriangle(0, i + 1, i);
            AddTriangle(last_offset, last_offset + i, last_offset + i + 1);
          }
        }
        MFEM_VERIFY(!model.mortar_triangles.empty(),
                    "A spatial surface mortar has no nondegenerate triangles!");
        ConfigureConductorConsistencyProbes(model);
      }
      else
      {
        auto AddSegment = [&](int begin, int end)
        {
          MFEM_VERIFY(begin >= 0 && begin < model.contour_size && end >= 0 &&
                          end < model.contour_size && begin != end,
                      "Surface-mortar contour contains an invalid segment!");
          double length_squared = 0.0;
          for (int d = 0; d < 3; d++)
          {
            const double delta = points[end][d] - points[begin][d];
            length_squared += delta * delta;
          }
          const double length = std::sqrt(length_squared);
          MFEM_VERIFY(length > 0.0,
                      "Surface-mortar contour contains a zero-length segment!");
          model.mortar_segments.push_back({begin, end, length, 1});
        };
        if (!model.open_contour_paths.empty())
        {
          for (const auto &path : model.open_contour_paths)
          {
            MFEM_VERIFY(
                path.indices.size() >= 2,
                "A translational surface-mortar path requires at least two points!");
            for (std::size_t i = 1; i < path.indices.size(); i++)
            {
              AddSegment(path.indices[i - 1], path.indices[i]);
            }
          }
        }
        else
        {
          int offset = 0;
          for (const int count : model.contour_groups)
          {
            MFEM_VERIFY(count >= 3,
                        "A closed translational surface-mortar contour requires at least "
                        "three points!");
            for (int i = 0; i < count; i++)
            {
              AddSegment(offset + i, offset + (i + 1) % count);
            }
            offset += count;
          }
        }
        for (const auto &segment : model.mortar_segments)
        {
          mass(segment.begin, segment.begin) += segment.length / 3.0;
          mass(segment.end, segment.end) += segment.length / 3.0;
          mass(segment.begin, segment.end) += segment.length / 6.0;
          mass(segment.end, segment.begin) += segment.length / 6.0;
        }
      }
      if (!model.spatial_mortar)
      {
        model.mortar_constant_load.SetSize(model.contour_size);
        for (int i = 0; i < model.contour_size; i++)
        {
          model.mortar_constant_load[i] = 0.0;
          for (int j = 0; j < model.contour_size; j++)
          {
            model.mortar_constant_load[i] += mass(i, j);
          }
        }
      }
      for (const int index : model.zero_trace_indices)
      {
        for (int j = 0; j < model.contour_size; j++)
        {
          mass(index, j) = 0.0;
          mass(j, index) = 0.0;
        }
        mass(index, index) = 1.0;
      }
      model.mortar_mass_inverse.SetSize(model.contour_size);
      mfem::DenseMatrixInverse(mass, true).GetInverseMatrix(model.mortar_mass_inverse);
    }
    auto domain_response = BuildDomainResponseMatrices(
        model_config, model.basis_size, model.zero_trace_indices, iodata.units);
    model.fabricated_domain = std::move(domain_response.fabricated);
    model.thin_domain = std::move(domain_response.thin);
    model.domain_defect = std::move(domain_response.defect);
    model.fixed_flux_transform = std::move(domain_response.fixed_flux_transform);
    model.fixed_flux_domain_defect = std::move(domain_response.fixed_flux_defect);
    auto surface_response = BuildSurfaceResponseMatrices(
        model_config, model.basis_size, iodata.units, target_matching_radius);
    model.fabricated_surfaces = std::move(surface_response.fabricated);
    model.surface_defects = std::move(surface_response.defects);
    // The conductor-consistency gate (decision 277) probes the trace mesh's conductor
    // vertices through the surface mortar only: a spatial model applied collocated keeps
    // its metal cross-sections unprobed, which is said once here.
    // Recorded untestable (Diagnostics.ConductorConsistency.UnprobedModels, decision 279
    // MINOR-2): a collocated spatial model, and a spatial surface-mortar trace mesh without
    // any conductor vertex (a finite-impedance coupon or a ring-path model without
    // ZeroTraceIndices: no metal cross-section at a fixed potential to probe).
    if (model.spatial_basis && !model.spatial_mortar && dimension == 3 &&
        HasExplicitTraceMesh(model_config))
    {
      Mpi::Warning(fespace.GetComm(),
                   "Conductor-consistency gate (decision 277) not evaluated for response "
                   "model \"{}\": its metal cross-sections are probed through the surface "
                   "mortar only (TraceCoupling \"SurfaceMortar\"), not collocated\n",
                   model.name);
      ConductorConsistencyDiagnostics()["UnprobedModels"].push_back(
          {{"Model", model.name},
           {"ModelIndex", model.idx},
           {"Reason", "collocated: the probes live on the surface-mortar trace mesh"}});
    }
    else if (model.spatial_mortar && dimension == 3 &&
             std::none_of(model.mortar_vertices.begin(), model.mortar_vertices.end(),
                          [](const auto &vertex) { return vertex.conductor > 0; }))
    {
      Mpi::Warning(
          fespace.GetComm(),
          "Conductor-consistency gate (decision 277) not evaluated for response "
          "model \"{}\": its trace mesh has no conductor vertex, so it has no metal "
          "cross-section at a fixed potential to probe\n",
          model.name);
      ConductorConsistencyDiagnostics()["UnprobedModels"].push_back(
          {{"Model", model.name},
           {"ModelIndex", model.idx},
           {"Reason", "no conductor vertex on the trace mesh: no metal cross-section at a "
                      "fixed potential to probe"}});
    }
    model_indices.emplace(model.idx, static_cast<int>(models.size()));
    models.push_back(std::move(model));
    basis_points.push_back(std::move(points));
  }
  dbc_tdof_list = laplace_op.GetDbcTDofList();

  // A three-dimensional spatial response represents one complete coupon volume and has
  // unit weight. Such volumes must never overlap: adding both dense defects would count
  // the shared physical domain twice. Use transformed matching-contour bounds as a
  // conservative fail-closed ownership check until an explicit partition-of-unity
  // representation is available.
  struct SpatialSupport
  {
    std::size_t patch = 0;
    ElementBox box;
  };
  std::vector<SpatialSupport> spatial_supports;
  std::set<std::size_t> spatially_owned_patches;
  // The placed patches: the configured ones with the continuation ownership applied (3D).
  std::vector<ResponsePatchData> placed_patches = config->patches;
  if (dimension == 3)
  {
    std::vector<std::string> skipped;
    const auto boxes = CollectSpatialSupports(
        *config, [&](int model_idx) -> const std::vector<std::array<double, 3>> *
        { return &basis_points[model_indices.at(model_idx)]; }, coordinate_scale, dimension,
        &skipped);
    const auto margin_overlaps = FindSpatialSupportMarginOverlaps(
        boxes, dimension, kSignatureParameterToleranceOverRadius * config->matching_radius);
    for (std::size_t patch_idx = 0; patch_idx < config->patches.size(); patch_idx++)
    {
      const auto &patch_config = config->patches[patch_idx];
      const auto model_it = model_indices.find(patch_config.model);
      MFEM_ASSERT(model_it != model_indices.end(), "Unknown response model!");
      const auto &model = models[model_it->second];
      if (!model.spatial_basis)
      {
        continue;
      }
      auto &support = spatial_supports.emplace_back();
      support.patch = patch_idx;
      for (const auto &local : basis_points[model_it->second])
      {
        ElementBox point_box;
        for (int d = 0; d < dimension; d++)
        {
          const double coordinate =
              patch_config.origin[d] +
              (local[0] * patch_config.axis_u[d] + local[1] * patch_config.axis_v[d] +
               local[2] * patch_config.axis_w[d]) /
                  coordinate_scale;
          point_box.min[d] = point_box.max[d] = coordinate;
        }
        support.box.Add(point_box);
      }
    }
    for (std::size_t i = 0; i < spatial_supports.size(); i++)
    {
      for (std::size_t j = i + 1; j < spatial_supports.size(); j++)
      {
        const auto first = spatial_supports[i].patch;
        const auto second = spatial_supports[j].patch;
        if (spatially_owned_patches.count(first) || spatially_owned_patches.count(second))
        {
          continue;
        }
        const int interpolation_group = config->patches[first].interpolation_group;
        if (interpolation_group > 0 &&
            interpolation_group == config->patches[second].interpolation_group)
        {
          continue;
        }
        if (!spatial_supports[i].box.InteriorOverlaps(spatial_supports[j].box, 3))
        {
          continue;
        }
        const auto &first_model = models[model_indices.at(config->patches[first].model)];
        const auto &second_model = models[model_indices.at(config->patches[second].model)];
        const bool first_cluster = first_model.topology == "spatial edge cluster";
        const bool second_cluster = second_model.topology == "spatial edge cluster";
        if (first_cluster != second_cluster)
        {
          // Exact spatial clusters own their complete local interaction neighborhood.
          // Give them deterministic priority over any overlapping corner/vertex patch;
          // retaining both would double count the intersection, while clipping the dense
          // coupon operator is not yet supported.
          const auto subordinate = first_cluster ? j : i;
          spatially_owned_patches.insert(spatial_supports[subordinate].patch);
          continue;
        }
        // Two clusters (decision 244): a margins-only overlap (no claim of either inside
        // the other's claims hull) is recorded below with its double-counted length; a
        // claim inside the other's hull is a true overlap of two coupon domains.
        const auto overlap = std::find_if(
            margin_overlaps.begin(), margin_overlaps.end(), [&](const auto &entry)
            { return entry.first_patch == first && entry.second_patch == second; });
        if (first_cluster && second_cluster && overlap != margin_overlaps.end() &&
            !overlap->claim_in_hull)
        {
          continue;
        }
        MFEM_ABORT("Three-dimensional response-correction matching volumes for patches "
                   << first + 1 << " (" << first_model.name << ") and " << second + 1
                   << " (" << second_model.name
                   << ") overlap without complete spatial-cluster ownership"
                   << (first_cluster && second_cluster
                           ? " (a claim of one cluster lies inside the other's claims)"
                           : "")
                   << ". Replace them with one coupled spatial model or an explicit "
                      "nonoverlapping partition!");
      }
    }
    if (!config->quantum_near_match.empty())
    {
      // Quantum near-matches (block (b) DESIGN section 4) in the operator record: the
      // cluster models applied for feature keys within the near-match of their own.
      nlohmann::json keys = nlohmann::json::array();
      for (const auto &entry : config->quantum_near_match)
      {
        keys.push_back({{"Model", entry.model},
                        {"ModelKey", entry.model_key},
                        {"FeatureKey", entry.feature_key},
                        {"MaxDeltaQuanta", entry.max_delta_quanta},
                        {"DifferingNumbers", entry.differing_numbers},
                        {"Features", entry.features}});
      }
      ownership_diagnostics["QuantumNearMatch"] = {
          {"Count", keys.size()},
          {"MaxQuanta", kClusterQuantumNearMatchMaxQuanta},
          {"Keys", std::move(keys)},
          {"Rule", "block (b) DESIGN section 4 (decision 303): a SpatialEdgeCluster model "
                   "applied for a feature key of the same topology within MaxQuanta "
                   "signature quanta of the model's key (the same geometry at the 1e-6 R "
                   "grid; the feature placed in its own canonical frame)"}};
    }
    if (!config->legacy_contract.empty())
    {
      // Legacy-contract aliases (USER decision 283) in the operator record: the legacy
      // models applied for contract-3 keys the library lists explicitly.
      nlohmann::json aliases = nlohmann::json::array();
      for (const auto &entry : config->legacy_contract)
      {
        aliases.push_back({{"Model", entry.model},
                           {"Key", entry.key},
                           {"ContextDigest", entry.context_digest},
                           {"Reason", entry.reason},
                           {"Features", entry.features}});
      }
      ownership_diagnostics["LegacyContract"] = {
          {"Count", aliases.size()},
          {"Aliases", std::move(aliases)},
          {"Rule", "USER decision 283: a legacy model (decision-236 straight-continuation "
                   "coupon) applied for a contract-3 key through the library's explicit "
                   "LegacyContractAliases entry (key + verified context digest); never a "
                   "fallback for any other key"}};
      Mpi::Warning("Legacy-contract aliases applied (USER decision 283): {:d} model(s) "
                   "serve contract-3 keys through explicit library aliases!\n",
                   static_cast<int>(config->legacy_contract.size()));
    }
    if (!margin_overlaps.empty())
    {
      ownership_diagnostics["SpatialSupportMarginOverlaps"] =
          DescribeSpatialSupportMarginOverlaps(margin_overlaps, *config, coordinate_scale);
      Mpi::Warning(fespace.GetComm(), "{}",
                   DescribeSpatialSupportMarginOverlapWarning(
                       ownership_diagnostics["SpatialSupportMarginOverlaps"]));
    }
    // Translational patches against the spatial supports (decisions 224 / 236 / 244):
    // every translational stretch wholly inside one spatial support is recorded
    // (Diagnostics + warning), never an abort — a Foreign one is a model mismatch of the
    // coupon (the S1p 41-edge loop end: three leads split by 1.738 um pieces of a 3-edge
    // stack between the cluster's claims, now absorbed by the identification); then the
    // continuation ownership clips the placed cells of every stretch continuing a
    // cluster's claims at that cluster's box face (the accepted transmon library: 112.4 um
    // of isolated-edge and strip cells on the continuations of its five coupons).
    {
      const auto records = FindTranslationalStretchInsideSpatialSupport(
          placed_patches, boxes, dimension,
          kSignatureParameterToleranceOverRadius * config->matching_radius);
      const auto ownership = ApplyContinuationOwnership(
          placed_patches, boxes, dimension,
          kSignatureParameterToleranceOverRadius * config->matching_radius,
          config->matching_radius);
      for (const auto &cell : ownership.cells)
      {
        if (placed_patches[cell.patch].weight <= 0.0)
        {
          spatially_owned_patches.insert(cell.patch);
        }
      }
      for (const auto &vertex : ownership.vertices)
      {
        spatially_owned_patches.insert(vertex.patch);
      }
      ownership_diagnostics["TranslationalStretchesInsideSpatialSupport"] =
          DescribeTranslationalOwnershipRecords(records, boxes, *config, coordinate_scale,
                                                skipped, &ownership);
      ownership_diagnostics["ContinuationOwnership"] =
          DescribeContinuationOwnership(ownership, boxes, *config, coordinate_scale);
      if (!records.empty())
      {
        Mpi::Warning(
            fespace.GetComm(), "{}",
            DescribeTranslationalOwnershipWarning(
                ownership_diagnostics["TranslationalStretchesInsideSpatialSupport"]));
      }
      if (!ownership.cells.empty())
      {
        Mpi::Print(fespace.GetComm(), "{}",
                   DescribeContinuationOwnershipSummary(
                       ownership_diagnostics["ContinuationOwnership"]));
      }
    }
    // Domain-boundary exclusion (decision 258): the placed patches with coupon points
    // outside the device mesh are not applied (weight 0, recorded). The same collective
    // test as the preflight's, on the placed patches not already owned; the mortar
    // resolution lookup and the point location below see applied patches only.
    {
      auto exclusions = FindDomainBoundaryExclusions(
          response_mesh, placed_patches,
          [&](int model_idx) -> const std::vector<std::array<double, 3>> *
          { return &basis_points[model_indices.at(model_idx)]; },
          [&](int model_idx) { return models[model_indices.at(model_idx)].spatial_basis; },
          [&](int model_idx) { return models[model_indices.at(model_idx)].name; },
          coordinate_scale, config->matching_radius, spatially_owned_patches);
      for (const auto &exclusion : exclusions.patches)
      {
        spatially_owned_patches.insert(exclusion.patch);
      }
      ownership_diagnostics["DomainBoundaryExclusions"] =
          DescribeDomainBoundaryExclusions(exclusions, *config, coordinate_scale);
      Mpi::Print(
          fespace.GetComm(),
          " Domain-boundary containment test: {:d} patches / {:d} points in {:.3f} s\n",
          exclusions.tested_patches, exclusions.tested_points, exclusions.wall_time);
      if (!exclusions.patches.empty())
      {
        Mpi::Print(fespace.GetComm(), "{}",
                   DescribeDomainBoundaryExclusionSummary(
                       ownership_diagnostics["DomainBoundaryExclusions"]));
      }
    }
  }

  const int rank = Mpi::Rank(fespace.GetComm());
  const int size = Mpi::Size(fespace.GetComm());
  int point_count = 0;
  std::vector<std::size_t> local_patch_indices;
  for (std::size_t patch_idx = 0; patch_idx < placed_patches.size(); patch_idx++)
  {
    if (spatially_owned_patches.count(patch_idx))
    {
      continue;
    }
    const auto &patch_config = placed_patches[patch_idx];
    const auto model_it = model_indices.find(patch_config.model);
    MFEM_VERIFY(model_it != model_indices.end(),
                "Response-correction patch refers to an unknown model index!");
    const auto &model = models[model_it->second];
    if (rank == 0)
    {
      patch_assignments.push_back({model.idx, patch_config.origin, patch_config.axis_u,
                                   patch_config.axis_v, patch_config.axis_w,
                                   patch_config.weight, static_cast<int>(patch_idx)});
    }
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
    MFEM_VERIFY(std::isfinite(patch_config.longitudinal_cell[0]) &&
                    std::isfinite(patch_config.longitudinal_cell[1]) &&
                    patch_config.longitudinal_cell[0] <= patch_config.longitudinal_cell[1],
                "Response-correction patch longitudinal cells must be ordered intervals!");
    patches.push_back(Patch{static_cast<int>(patch_idx), model_it->second, 0, basis_size, 0,
                            patch_config.longitudinal_cell, 1, 0.0, patch_config.weight});
    // The provenance of the conductor-consistency record (decision 277): the feature, the
    // claimed portions of a spatial cluster and the longitudinal cell, mesh units.
    auto &placed = patches.back();
    placed.feature = patch_config.provenance.feature;
    placed.cell_length =
        patch_config.longitudinal_cell[1] - patch_config.longitudinal_cell[0];
    for (const auto &claim : patch_config.provenance.claims)
    {
      double length_squared = 0.0;
      for (int d = 0; d < 3; d++)
      {
        const double delta = claim.p1[d] - claim.p0[d];
        length_squared += delta * delta;
      }
      placed.claim_length += std::sqrt(length_squared);
    }
    basis_size += model.basis_size;
  }

  // Surface-mortar quadrature follows the mesh resolution local to each matching surface,
  // not the smallest element anywhere in the distributed device. This keeps AMR in an
  // unrelated hotspot from refining every translational response patch.
  std::vector<std::size_t> mortar_patch_indices;
  for (std::size_t patch_idx = 0; patch_idx < patches.size(); patch_idx++)
  {
    if (models[patches[patch_idx].model].surface_mortar)
    {
      mortar_patch_indices.push_back(patch_idx);
    }
  }
  int global_mortar_patch_count = static_cast<int>(mortar_patch_indices.size());
  Mpi::GlobalSum(1, &global_mortar_patch_count, fespace.GetComm());
  // One point locator per construction (= per AMR cycle) for the mortar-resolution probe
  // and the response points below: no global search structure on the device mesh
  // (decision 346 (b)).
  DistributedPointLocator point_locator(response_mesh, dimension);
  if (global_mortar_patch_count > 0)
  {
    mfem::Vector centers(dimension * mortar_patch_indices.size());
    for (std::size_t i = 0; i < mortar_patch_indices.size(); i++)
    {
      const std::size_t patch_idx = mortar_patch_indices[i];
      const auto &patch = patches[patch_idx];
      const auto &patch_config = placed_patches[local_patch_indices[patch_idx]];
      const auto &point = basis_points[patch.model].front();
      for (int d = 0; d < dimension; d++)
      {
        centers(d * mortar_patch_indices.size() + i) =
            patch_config.origin[d] +
            (point[0] * patch_config.axis_u[d] + point[1] * patch_config.axis_v[d] +
             (models[patch.model].spatial_mortar ? point[2] * patch_config.axis_w[d]
                                                 : 0.0)) /
                coordinate_scale;
      }
    }
    // The owning element's size (the smallest singular value of its Jacobian).
    const std::function<double(int)> element_size = [&](int element)
    { return response_mesh.GetElementSize(element, 1); };
    const auto located = point_locator.Locate(centers, &element_size);
    const auto &local_resolution = located.owner_values;
    for (std::size_t i = 0; i < mortar_patch_indices.size(); i++)
    {
      // The applied patches only (decision 258): a patch whose placed coupon section
      // leaves the mesh was excluded above, so a failure here is an error of the mesh or
      // of the placement, named.
      const auto &patch = patches[mortar_patch_indices[i]];
      MFEM_VERIFY(
          located.owners[i] < Mpi::Size(fespace.GetComm()) &&
              std::isfinite(local_resolution[i]) && local_resolution[i] > 0.0,
          "Unable to determine a local surface-mortar mesh resolution at the first "
          "basis point ("
              << centers(0 * mortar_patch_indices.size() + i) << ", "
              << centers(1 * mortar_patch_indices.size() + i) << ", "
              << (dimension == 3 ? centers(2 * mortar_patch_indices.size() + i) : 0.0)
              << ") of patch " << patch.global_index + 1 << " (model "
              << models[patch.model].name << ")!");
      patches[mortar_patch_indices[i]].mortar_resolution =
          local_resolution[i] / config->mortar_oversampling;
    }
  }

  point_count = 0;
  for (auto &patch : patches)
  {
    const auto &model = models[patch.model];
    patch.point_offset = point_count;
    patch.point_count = model.contour_size + 1 + model.conductor_state_count;
    if (model.surface_mortar)
    {
      MFEM_ASSERT(patch.mortar_resolution > 0.0,
                  "Missing local surface-mortar resolution!");
      if (model.spatial_mortar)
      {
        patch.mortar_longitudinal_subdivisions = 1;
        int surface_sample_count = 0;
        for (const auto &triangle : model.mortar_triangles)
        {
          const int subdivisions = std::max(
              1, static_cast<int>(std::ceil(triangle.maximum_edge_length /
                                            coordinate_scale / patch.mortar_resolution)));
          surface_sample_count += 4 * subdivisions * subdivisions;
        }
        patch.point_count = surface_sample_count + model.conductor_state_count + 1;
        // The conductor-consistency probes (decision 277) follow the references.
        patch.probe_point_count =
            static_cast<int>(model.plane_conductor_vertices.size() +
                             model.off_plane_conductor_vertices.size());
        patch.point_count += patch.probe_point_count;
      }
      else
      {
        // The strip (the patch's longitudinal cell, mesh units) is sampled at the same
        // resolution as the contour, so that the projection is a surface integral over
        // the matching strip rather than a line integral in one cross-section.
        const double strip_length =
            patch.mortar_longitudinal_strip[1] - patch.mortar_longitudinal_strip[0];
        patch.mortar_longitudinal_subdivisions =
            dimension == 3 ? std::max(1, static_cast<int>(std::ceil(
                                             strip_length / patch.mortar_resolution)))
                           : 1;
        int contour_sample_count = 0;
        for (const auto &segment : model.mortar_segments)
        {
          const int subdivisions =
              std::max(1, static_cast<int>(std::ceil(segment.length / coordinate_scale /
                                                     patch.mortar_resolution)));
          contour_sample_count += 2 * subdivisions;
        }
        patch.point_count = patch.mortar_longitudinal_subdivisions *
                            (contour_sample_count + model.conductor_state_count + 1);
      }
    }
    point_count += patch.point_count;
  }
  global_patch_count = static_cast<int>(patches.size());
  global_basis_size = basis_size;
  long long int global_point_count = point_count;
  Mpi::GlobalSum(1, &global_point_count, fespace.GetComm());
  long long int global_probe_count = 0;
  for (const auto &patch : patches)
  {
    global_probe_count += patch.probe_point_count;
  }
  Mpi::GlobalSum(1, &global_probe_count, fespace.GetComm());
  // Runaway guard, not a memory model: every point stores its element stencil on the
  // owning rank (one int and one double per element dof: ~0.7 kB at p5 on tetrahedra,
  // 56 dofs; ~2.6 kB on hexahedra, 216 dofs) plus ~50 B of query bookkeeping, so the
  // bound corresponds to ~2.8 TB of stencils at p5 on tetrahedra spread over the ranks. The
  // translational strips sample their longitudinal cells at the local mortar resolution
  // (the element size at the patch's first basis point): the transmon device (17,395
  // translational patches, 35 mm of matched perimeter, 4 mm characteristic length) needs
  // 4.3e6 points on its initial mesh and 1.0e7 after 11 AMR cycles at p4 (measured).
  constexpr long long int maximum_experimental_mortar_points = 4000000000LL;
  MFEM_VERIFY(config->trace_coupling !=
                      ResponseCorrectionData::TraceCoupling::SURFACE_MORTAR ||
                  global_point_count <= maximum_experimental_mortar_points,
              "Experimental SurfaceMortar trace quadrature requires "
                  << global_point_count
                  << " distributed points, exceeding its current safety limit of "
                  << maximum_experimental_mortar_points
                  << ". Assemble and compress the mortar projection before continuing "
                     "this AMR level!");
  std::array<int, 3> domain_mode_patch_count{};
  for (const auto &patch : patches)
  {
    domain_mode_patch_count[static_cast<std::size_t>(
        models[patch.model].domain_correction_mode)]++;
  }
  Mpi::GlobalSum(1, &global_patch_count, fespace.GetComm());
  Mpi::GlobalSum(1, &global_basis_size, fespace.GetComm());
  Mpi::GlobalSum(static_cast<int>(domain_mode_patch_count.size()),
                 domain_mode_patch_count.data(), fespace.GetComm());

  mfem::Vector xyz(dimension * point_count);
  std::vector<std::vector<Point2D>> polygons;
  if (dimension == 2)
  {
    polygons.reserve(placed_patches.size());
    for (const auto &patch_config : placed_patches)
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
    const auto &patch_config = placed_patches[local_patch_indices[patch_idx]];
    const auto &model = models[patch.model];
    const auto &local_points = basis_points[patch.model];
    auto axis_w = patch_config.axis_w;
    if (model.surface_mortar && dimension == 3 &&
        std::all_of(axis_w.begin(), axis_w.end(),
                    [](double value) { return std::abs(value) < 1.0e-14; }))
    {
      axis_w = Cross(patch_config.axis_u, patch_config.axis_v);
    }
    double norm_u = 0.0, norm_v = 0.0, norm_w = 0.0;
    double dot_uv = 0.0, dot_uw = 0.0, dot_vw = 0.0;
    for (int d = 0; d < dimension; d++)
    {
      norm_u += patch_config.axis_u[d] * patch_config.axis_u[d];
      norm_v += patch_config.axis_v[d] * patch_config.axis_v[d];
      norm_w += axis_w[d] * axis_w[d];
      dot_uv += patch_config.axis_u[d] * patch_config.axis_v[d];
      dot_uw += patch_config.axis_u[d] * axis_w[d];
      dot_vw += patch_config.axis_v[d] * axis_w[d];
    }
    norm_u = std::sqrt(norm_u);
    norm_v = std::sqrt(norm_v);
    MFEM_VERIFY(std::abs(norm_u - 1.0) < 1.0e-10 && std::abs(norm_v - 1.0) < 1.0e-10 &&
                    std::abs(dot_uv) < 1.0e-10,
                "Surface response correction AxisU and AxisV must be orthonormal in the "
                "coupon cross-section!");
    if (model.spatial_basis || (model.surface_mortar && dimension == 3))
    {
      MFEM_VERIFY(dimension == 3 && std::abs(std::sqrt(norm_w) - 1.0) < 1.0e-10 &&
                      std::abs(dot_uw) < 1.0e-10 && std::abs(dot_vw) < 1.0e-10,
                  "A spatial or surface-mortar response basis requires an orthonormal "
                  "three-dimensional coupon frame!");
    }

    if (model.surface_mortar)
    {
      constexpr double gauss_offset = 0.5 / 1.7320508075688772935;
      if (model.spatial_mortar)
      {
        for (const auto &triangle : model.mortar_triangles)
        {
          const int subdivisions = std::max(
              1, static_cast<int>(std::ceil(triangle.maximum_edge_length /
                                            coordinate_scale / patch.mortar_resolution)));
          for (int first = 0; first < subdivisions; first++)
          {
            for (const double first_offset : {-gauss_offset, gauss_offset})
            {
              const double u =
                  (static_cast<double>(first) + 0.5 + first_offset) / subdivisions;
              for (int second = 0; second < subdivisions; second++)
              {
                for (const double second_offset : {-gauss_offset, gauss_offset})
                {
                  const double v =
                      (static_cast<double>(second) + 0.5 + second_offset) / subdivisions;
                  const std::array<double, 3> barycentric = {u, (1.0 - u) * v,
                                                             (1.0 - u) * (1.0 - v)};
                  Point3D local{};
                  for (int q = 0; q < 3; q++)
                  {
                    for (int d = 0; d < 3; d++)
                    {
                      local[d] += barycentric[q] *
                                  model.mortar_vertices[triangle.vertices[q]].point[d];
                    }
                  }
                  for (int d = 0; d < dimension; d++)
                  {
                    xyz(d * point_count + point) =
                        patch_config.origin[d] +
                        (local[0] * patch_config.axis_u[d] +
                         local[1] * patch_config.axis_v[d] + local[2] * axis_w[d]) /
                            coordinate_scale;
                  }
                  point++;
                }
              }
            }
          }
        }
        for (const auto &reference : patch_config.conductor_references)
        {
          for (int d = 0; d < dimension; d++)
          {
            xyz(d * point_count + point) =
                patch_config.origin[d] + reference[0] * patch_config.axis_u[d] +
                reference[1] * patch_config.axis_v[d] + reference[2] * axis_w[d];
          }
          point++;
        }
        // The conductor-consistency probes: the conductor vertices of the trace mesh, the
        // plane ones first (the gate), then those off the plane (recorded).
        auto &probes = patches[patch_idx].probe_points;
        probes.clear();
        for (const auto *vertices :
             {&model.plane_conductor_vertices, &model.off_plane_conductor_vertices})
        {
          for (const int v : *vertices)
          {
            const auto &local = model.mortar_vertices[v].point;
            std::array<double, 3> probe{};
            for (int d = 0; d < dimension; d++)
            {
              probe[d] = patch_config.origin[d] +
                         (local[0] * patch_config.axis_u[d] +
                          local[1] * patch_config.axis_v[d] + local[2] * axis_w[d]) /
                             coordinate_scale;
              xyz(d * point_count + point) = probe[d];
            }
            probes.push_back(probe);
            point++;
          }
        }
        MFEM_ASSERT(static_cast<int>(probes.size()) == patch.probe_point_count,
                    "Incorrect conductor-consistency probe count!");
      }
      else
      {
        for (int longitudinal = 0; longitudinal < patch.mortar_longitudinal_subdivisions;
             longitudinal++)
        {
          // Midpoint of the longitudinal slice within the strip (offset along axis_w).
          const double longitudinal_coordinate =
              dimension == 3 ? patch.mortar_longitudinal_strip[0] +
                                   (patch.mortar_longitudinal_strip[1] -
                                    patch.mortar_longitudinal_strip[0]) *
                                       (static_cast<double>(longitudinal) + 0.5) /
                                       patch.mortar_longitudinal_subdivisions
                             : 0.0;
          for (const auto &segment : model.mortar_segments)
          {
            const auto &begin = local_points[segment.begin];
            const auto &end = local_points[segment.end];
            const int transverse_subdivisions =
                std::max(1, static_cast<int>(std::ceil(segment.length / coordinate_scale /
                                                       patch.mortar_resolution)));
            for (int subdivision = 0; subdivision < transverse_subdivisions; subdivision++)
            {
              for (const double offset : {-gauss_offset, gauss_offset})
              {
                const double t = (static_cast<double>(subdivision) + 0.5 + offset) /
                                 transverse_subdivisions;
                for (int d = 0; d < dimension; d++)
                {
                  const double local_u = (1.0 - t) * begin[0] + t * end[0];
                  const double local_v = (1.0 - t) * begin[1] + t * end[1];
                  xyz(d * point_count + point) = patch_config.origin[d] +
                                                 (local_u * patch_config.axis_u[d] +
                                                  local_v * patch_config.axis_v[d]) /
                                                     coordinate_scale +
                                                 longitudinal_coordinate * axis_w[d];
                }
                point++;
              }
            }
          }
          for (const auto &reference : patch_config.conductor_references)
          {
            for (int d = 0; d < dimension; d++)
            {
              xyz(d * point_count + point) =
                  patch_config.origin[d] + reference[0] * patch_config.axis_u[d] +
                  reference[1] * patch_config.axis_v[d] +
                  (reference[2] + longitudinal_coordinate) * axis_w[d];
            }
            point++;
          }
        }
      }
    }
    else
    {
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
            coordinate += local[2] * axis_w[d] / coordinate_scale;
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
              reference[1] * patch_config.axis_v[d] + reference[2] * axis_w[d];
        }
        point++;
      }
    }
    MFEM_ASSERT(point == patch.point_offset + patch.point_count,
                "Incorrect surface-response patch point count!");
  }
  for (std::size_t i = 0; i < polygons.size(); i++)
  {
    for (std::size_t j = i + 1; j < polygons.size(); j++)
    {
      const int interpolation_group = placed_patches[i].interpolation_group;
      if (interpolation_group > 0 &&
          interpolation_group == placed_patches[j].interpolation_group)
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
    ConfigurePointCommunication(point_locator, xyz, dimension);
  }

  Mpi::Print(
      "\nConfigured surface response correction:\n"
      " Coupon models: {:d}\n"
      " Global patches: {:d}\n"
      " Domain-coupled patches (disabled/fixed-trace/fixed-flux): {:d}/{:d}/{:d}\n"
      " Trace quadrature points: {:d}\n"
      " Total trace coefficients: {:d}\n",
      static_cast<int>(models.size()), global_patch_count,
      domain_mode_patch_count[static_cast<std::size_t>(DomainCorrectionMode::DISABLED)],
      domain_mode_patch_count[static_cast<std::size_t>(DomainCorrectionMode::FIXED_TRACE)],
      domain_mode_patch_count[static_cast<std::size_t>(DomainCorrectionMode::FIXED_FLUX)],
      global_point_count - global_probe_count, global_basis_size);
  if (global_probe_count > 0)
  {
    Mpi::Print(" Conductor-consistency probe points (decision 277): {:d}\n",
               global_probe_count);
  }
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
  MFEM_VERIFY(request->trace_coupling !=
                  ResponseCorrectionData::TraceCoupling::SURFACE_MORTAR,
              "Maxwell SurfaceMortar response trace coupling is not implemented yet; "
              "use Collocated while the H(curl)-compatible mortar projection is under "
              "development!");
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
  const double target_matching_radius =
      TargetInterfaceMatchingRadius(iodata, config.target_interfaces);
  for (const auto &model_config : config.models)
  {
    MFEM_VERIFY(model_config.idx > 0 &&
                    model_indices.find(model_config.idx) == model_indices.end(),
                "Response-correction model indices must be positive and unique!");
    auto points = ModelBasisPoints(model_config);
    MFEM_VERIFY(points.size() >= 3,
                "Maxwell response-correction contours require at least three points!");
    // The corner family's trace basis rule (slave trace vertices at the box corners, or a
    // basis constructed at the device angle) is electrostatic only: the Maxwell contour
    // lines run between basis knots and would cut the box corners.
    MFEM_VERIFY(model_config.constructed_basis_points.empty() &&
                    (model_config.trace_vertices.empty() ||
                     !ModelTraceMesh(model_config).HasSlaveVertices()),
                "Maxwell surface response correction does not support trace meshes with "
                "slave vertices (the corner family's trace basis rule is electrostatic "
                "only)!");
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
    MFEM_VERIFY(model_config.interior_trace_count == 0,
                "Maxwell response correction does not support cap-interior trace "
                "coefficients (InteriorTraceCount): every coefficient must lie on a "
                "Maxwell contour!");
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
        model_config, model.basis_size, model.zero_trace_indices, iodata.units);
    model.fabricated_domain = std::move(domain_response.fabricated);
    model.thin_domain = std::move(domain_response.thin);
    model.domain_defect = std::move(domain_response.defect);
    model.fixed_flux_transform = std::move(domain_response.fixed_flux_transform);
    model.fixed_flux_domain_defect = std::move(domain_response.fixed_flux_defect);
    auto surface_response = BuildSurfaceResponseMatrices(
        model_config, model.basis_size, iodata.units, target_matching_radius);
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
    if (rank == 0)
    {
      patch_assignments.push_back({model.idx, patch_config.origin, patch_config.axis_u,
                                   patch_config.axis_v, patch_config.axis_w,
                                   patch_config.weight, static_cast<int>(patch_index)});
    }
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

    patches.push_back(Patch{static_cast<int>(patches.size()), model_it->second, 0,
                            basis_size, 0, patch_config.longitudinal_cell, 1, 0.0,
                            patch_config.weight});
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
  DistributedPointLocator point_locator(const_cast<mfem::ParMesh &>(fespace.GetParMesh()),
                                        dimension);
  ConfigurePointCommunication(point_locator, xyz, dimension, &weighted_tangents);
}

void SurfaceResponseOperator::ConfigurePointCommunication(
    DistributedPointLocator &locator, const mfem::Vector &xyz, int dimension,
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

  const int size = Mpi::Size(fespace.GetComm());
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

  // Fail closed on an unlocated point, naming its patch (decision 258: the domain-boundary
  // exclusion removes the patches whose placed coupon sections leave the mesh before this
  // point location; anything else not located is an error, never silently dropped).
  auto DescribeUnlocatedPoint = [&](int point)
  {
    std::array<double, 3> coordinate{};
    for (int d = 0; d < dimension; d++)
    {
      coordinate[d] = xyz(d * point_query_count + point);
    }
    std::string description =
        fmt::format("Surface-response contour point {:d} at ({:.9e}, {:.9e}, {:.9e}) (mesh "
                    "units) could not be located in the device mesh",
                    point, coordinate[0], coordinate[1], coordinate[2]);
    for (const auto &patch : patches)
    {
      if (patch.point_count > 0 && point >= patch.point_offset &&
          point < patch.point_offset + patch.point_count)
      {
        description += fmt::format(" (patch {:d}, model {})", patch.global_index + 1,
                                   models[patch.model].name);
        break;
      }
    }
    return description + "!";
  };
  auto located = locator.Locate(xyz);
  candidate_query_count += located.candidate_queries;
  fallback_query_count += located.fallback_queries;
  const auto &point_owners = located.owners;
  const auto &point_elements = located.elements;
  const auto &point_references = located.references;
  for (int point = 0; point < point_query_count; point++)
  {
    MFEM_VERIFY(point_owners[point] < size, DescribeUnlocatedPoint(point));
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
  constexpr double gauss_offset = 0.5 / 1.7320508075688772935;
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    if (model.surface_mortar)
    {
      mortar_load.SetSize(model.contour_size);
      mortar_load = 0.0;
      std::vector<double> references(model.conductor_state_count + 1, 0.0);
      int point = patch.point_offset;
      if (model.spatial_mortar)
      {
        for (const auto &triangle : model.mortar_triangles)
        {
          const int subdivisions = std::max(
              1,
              static_cast<int>(std::ceil(triangle.maximum_edge_length /
                                         mesh_coordinate_scale / patch.mortar_resolution)));
          for (int first = 0; first < subdivisions; first++)
          {
            for (const double first_offset : {-gauss_offset, gauss_offset})
            {
              const double u =
                  (static_cast<double>(first) + 0.5 + first_offset) / subdivisions;
              for (int second = 0; second < subdivisions; second++)
              {
                for (const double second_offset : {-gauss_offset, gauss_offset})
                {
                  const double v =
                      (static_cast<double>(second) + 0.5 + second_offset) / subdivisions;
                  const std::array<double, 3> barycentric = {u, (1.0 - u) * v,
                                                             (1.0 - u) * (1.0 - v)};
                  const double quadrature_weight =
                      triangle.area * (1.0 - u) / (2.0 * subdivisions * subdivisions);
                  const double value = correction(point++);
                  for (int q = 0; q < 3; q++)
                  {
                    model.mortar_vertices[triangle.vertices[q]].ForEachBasis(
                        [&](int basis, double weight)
                        {
                          mortar_load[basis] +=
                              weight * quadrature_weight * barycentric[q] * value;
                        });
                  }
                }
              }
            }
          }
        }
        for (double &reference : references)
        {
          reference = correction(point++);
        }
        for (int i = 0; i < model.contour_size; i++)
        {
          mortar_load[i] -= references[0] * model.mortar_constant_load[i];
          for (int state = 0; state < model.conductor_state_count; state++)
          {
            mortar_load[i] -= (references[state + 1] - references[0]) *
                              model.mortar_conductor_loads[state][i];
          }
        }
        for (const int index : model.zero_trace_indices)
        {
          mortar_load[index] = 0.0;
        }
      }
      else
      {
        const double longitudinal_weight = 1.0 / patch.mortar_longitudinal_subdivisions;
        for (int longitudinal = 0; longitudinal < patch.mortar_longitudinal_subdivisions;
             longitudinal++)
        {
          for (const auto &segment : model.mortar_segments)
          {
            const int transverse_subdivisions = std::max(
                1, static_cast<int>(std::ceil(segment.length / mesh_coordinate_scale /
                                              patch.mortar_resolution)));
            const double quadrature_weight =
                longitudinal_weight * segment.length / (2.0 * transverse_subdivisions);
            for (int subdivision = 0; subdivision < transverse_subdivisions; subdivision++)
            {
              for (const double offset : {-gauss_offset, gauss_offset})
              {
                const double t = (static_cast<double>(subdivision) + 0.5 + offset) /
                                 transverse_subdivisions;
                const double value = correction(point++);
                mortar_load[segment.begin] += quadrature_weight * (1.0 - t) * value;
                mortar_load[segment.end] += quadrature_weight * t * value;
              }
            }
          }
          for (double &reference : references)
          {
            reference += longitudinal_weight * correction(point++);
          }
        }
      }
      MFEM_ASSERT(point == patch.point_offset + patch.point_count - patch.probe_point_count,
                  "Incorrect surface-mortar point count!");
      mortar_coefficients.SetSize(model.contour_size);
      model.mortar_mass_inverse.Mult(mortar_load, mortar_coefficients);
      for (int i = 0; i < model.contour_size; i++)
      {
        values(patch.trace_offset + i) =
            mortar_coefficients[i] - (model.spatial_mortar ? 0.0 : references[0]);
      }
      for (int state = 0; state < model.conductor_state_count; state++)
      {
        values(patch.trace_offset + model.contour_size + state) =
            references[state + 1] - references[0];
      }
      continue;
    }
    const double reference = correction(patch.point_offset + model.contour_size);
    for (int i = 0; i < model.contour_size; i++)
    {
      values(patch.trace_offset + i) = correction(patch.point_offset + i) - reference;
    }
    // A PEC-constrained knot (ZeroTraceIndices) carries a zero trace by the coupon's
    // semantics whatever the device potential at its point (a knot on the fabricated slab
    // top lies in the air of the thin device): the collocated lift enforces it as the
    // surface mortar does (corner-family review 2026-09-29, m5).
    for (const int index : model.zero_trace_indices)
    {
      values(patch.trace_offset + index) = 0.0;
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
  constexpr double gauss_offset = 0.5 / 1.7320508075688772935;
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    if (model.surface_mortar)
    {
      Vector patch_values(const_cast<double *>(values.GetData()) + patch.trace_offset,
                          model.basis_size);
      Vector contour_values(patch_values.GetData(), model.contour_size);
      mortar_load.SetSize(model.contour_size);
      model.mortar_mass_inverse.MultTranspose(contour_values, mortar_load);
      std::vector<double> references(model.conductor_state_count + 1, 0.0);
      if (model.spatial_mortar)
      {
        for (const int index : model.zero_trace_indices)
        {
          mortar_load[index] = 0.0;
        }
        references[0] = -mfem::InnerProduct(model.mortar_constant_load, mortar_load);
        for (int state = 0; state < model.conductor_state_count; state++)
        {
          const double coupling =
              mfem::InnerProduct(model.mortar_conductor_loads[state], mortar_load);
          const double direct = patch_values[model.contour_size + state];
          references[0] += coupling - direct;
          references[state + 1] = direct - coupling;
        }
      }
      else
      {
        for (int i = 0; i < model.contour_size; i++)
        {
          references[0] -= patch_values[i];
        }
        for (int state = 0; state < model.conductor_state_count; state++)
        {
          references[state + 1] += patch_values[model.contour_size + state];
          references[0] -= patch_values[model.contour_size + state];
        }
      }
      int point = patch.point_offset;
      if (model.spatial_mortar)
      {
        for (const auto &triangle : model.mortar_triangles)
        {
          const int subdivisions = std::max(
              1,
              static_cast<int>(std::ceil(triangle.maximum_edge_length /
                                         mesh_coordinate_scale / patch.mortar_resolution)));
          for (int first = 0; first < subdivisions; first++)
          {
            for (const double first_offset : {-gauss_offset, gauss_offset})
            {
              const double u =
                  (static_cast<double>(first) + 0.5 + first_offset) / subdivisions;
              for (int second = 0; second < subdivisions; second++)
              {
                for (const double second_offset : {-gauss_offset, gauss_offset})
                {
                  const double v =
                      (static_cast<double>(second) + 0.5 + second_offset) / subdivisions;
                  const std::array<double, 3> barycentric = {u, (1.0 - u) * v,
                                                             (1.0 - u) * (1.0 - v)};
                  const double quadrature_weight =
                      triangle.area * (1.0 - u) / (2.0 * subdivisions * subdivisions);
                  double value = 0.0;
                  for (int q = 0; q < 3; q++)
                  {
                    model.mortar_vertices[triangle.vertices[q]].ForEachBasis(
                        [&](int basis, double weight)
                        { value += weight * barycentric[q] * mortar_load[basis]; });
                  }
                  correction[point++] = quadrature_weight * value;
                }
              }
            }
          }
        }
        for (const double reference : references)
        {
          correction[point++] = reference;
        }
      }
      else
      {
        const double longitudinal_weight = 1.0 / patch.mortar_longitudinal_subdivisions;
        for (int longitudinal = 0; longitudinal < patch.mortar_longitudinal_subdivisions;
             longitudinal++)
        {
          for (const auto &segment : model.mortar_segments)
          {
            const int transverse_subdivisions = std::max(
                1, static_cast<int>(std::ceil(segment.length / mesh_coordinate_scale /
                                              patch.mortar_resolution)));
            const double quadrature_weight =
                longitudinal_weight * segment.length / (2.0 * transverse_subdivisions);
            for (int subdivision = 0; subdivision < transverse_subdivisions; subdivision++)
            {
              for (const double offset : {-gauss_offset, gauss_offset})
              {
                const double t = (static_cast<double>(subdivision) + 0.5 + offset) /
                                 transverse_subdivisions;
                correction[point++] =
                    quadrature_weight *
                    ((1.0 - t) * mortar_load[segment.begin] + t * mortar_load[segment.end]);
              }
            }
          }
          for (const double reference : references)
          {
            correction[point++] = longitudinal_weight * reference;
          }
        }
      }
      MFEM_ASSERT(point == patch.point_offset + patch.point_count - patch.probe_point_count,
                  "Incorrect surface-mortar point count!");
      continue;
    }
    double reference = 0.0;
    std::vector<bool> constrained(model.contour_size, false);
    for (const int index : model.zero_trace_indices)
    {
      constrained[index] = true;
    }
    for (int i = 0; i < model.contour_size; i++)
    {
      // The transpose of the collocated lift with its PEC knots zeroed.
      const double value = constrained[i] ? 0.0 : values(patch.trace_offset + i);
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

void SurfaceResponseOperator::ApplyDomainDefect(const Vector &x, Vector &y,
                                                bool fixed_trace) const
{
  ApplyTrace(x, trace);
  response.SetSize(trace.Size());
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    Vector patch_trace(trace.GetData() + patch.trace_offset, model.basis_size);
    Vector patch_response(response.GetData() + patch.trace_offset, model.basis_size);
    const auto mode =
        fixed_trace ? DomainCorrectionMode::FIXED_TRACE : model.domain_correction_mode;
    switch (mode)
    {
      case DomainCorrectionMode::DISABLED:
        patch_response = 0.0;
        break;
      case DomainCorrectionMode::FIXED_TRACE:
        model.domain_defect.Mult(patch_trace, patch_response);
        break;
      case DomainCorrectionMode::FIXED_FLUX:
        model.fixed_flux_domain_defect.Mult(patch_trace, patch_response);
        break;
    }
    patch_response *= patch.weight;
  }
  ApplyTraceTranspose(response, y);
}

void SurfaceResponseOperator::ApplyUneliminated(const Vector &x, Vector &y) const
{
  BlockTimer timer(Timer::RESPONSE_APPLY);
  ApplyDomainDefect(x, y, false);
}

void SurfaceResponseOperator::FixedTraceDomainDefectMult(const Vector &x, Vector &y) const
{
  ApplyDomainDefect(x, y, true);
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
    if (model.domain_correction_mode == DomainCorrectionMode::FIXED_TRACE)
    {
      energy.domain +=
          0.5 * patch.weight * QuadraticForm(model.domain_defect, patch_trace, response);
    }
    else if (model.domain_correction_mode == DomainCorrectionMode::FIXED_FLUX)
    {
      energy.domain += 0.5 * patch.weight *
                       QuadraticForm(model.fixed_flux_domain_defect, patch_trace, response);
    }
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
    const double domain_correction_fixed_trace =
        0.5 * patch.weight * QuadraticForm(model.domain_defect, patch_trace, response);
    const bool evaluate_fixed_flux =
        include_fixed_flux ||
        model.domain_correction_mode == DomainCorrectionMode::FIXED_FLUX;
    double domain_correction_fixed_flux = 0.0;
    if (evaluate_fixed_flux)
    {
      fixed_flux.SetSize(model.basis_size);
      model.fixed_flux_transform.Mult(patch_trace, fixed_flux);
      domain_correction_fixed_flux =
          0.5 * patch.weight *
          (QuadraticForm(model.fabricated_domain, fixed_flux, response) -
           QuadraticForm(model.thin_domain, patch_trace, response));
    }
    if (include_fixed_flux)
    {
      result.domain_correction += domain_correction_fixed_trace;
      contribution.domain_correction += domain_correction_fixed_trace;
      result.domain_correction_fixed_flux += domain_correction_fixed_flux;
      contribution.domain_correction_fixed_flux += domain_correction_fixed_flux;
    }
    else if (model.domain_correction_mode == DomainCorrectionMode::FIXED_TRACE)
    {
      result.domain_correction += domain_correction_fixed_trace;
      contribution.domain_correction += domain_correction_fixed_trace;
    }
    else if (model.domain_correction_mode == DomainCorrectionMode::FIXED_FLUX)
    {
      result.domain_correction += domain_correction_fixed_flux;
      contribution.domain_correction += domain_correction_fixed_flux;
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

std::vector<SurfaceResponseOperator::PatchTrace>
SurfaceResponseOperator::GetSpatialPatchTraces(const Vector &x) const
{
  MFEM_VERIFY(!maxwell,
              "Spatial patch-trace export currently supports electrostatic response "
              "correction only!");
  ApplyTrace(x, trace);
  std::vector<double> local_records;
  for (const auto &patch : patches)
  {
    const auto &model = models[patch.model];
    if (!model.spatial_basis)
    {
      continue;
    }
    local_records.push_back(static_cast<double>(patch.global_index));
    local_records.push_back(static_cast<double>(model.idx));
    local_records.push_back(static_cast<double>(model.contour_size));
    local_records.push_back(static_cast<double>(model.basis_size));
    for (int i = 0; i < model.basis_size; i++)
    {
      local_records.push_back(trace(patch.trace_offset + i));
    }
  }

  MFEM_VERIFY(local_records.size() <=
                  static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "Local spatial patch-trace data exceeds the MPI count limit!");
  const int local_value_count = static_cast<int>(local_records.size());
  std::vector<int> value_counts(Mpi::Size(fespace.GetComm()));
  Mpi::Allgather(1, &local_value_count, value_counts.data(), fespace.GetComm());
  std::vector<int> value_offsets(value_counts.size());
  int total_values = 0;
  for (std::size_t rank = 0; rank < value_counts.size(); rank++)
  {
    value_offsets[rank] = total_values;
    MFEM_VERIFY(value_counts[rank] <= std::numeric_limits<int>::max() - total_values,
                "Global spatial patch-trace data exceeds the MPI count limit!");
    total_values += value_counts[rank];
  }
  std::vector<double> records(total_values);
  Mpi::Allgatherv(local_value_count, local_records.data(), records.data(),
                  value_counts.data(), value_offsets.data(), fespace.GetComm());

  std::vector<PatchTrace> traces;
  std::size_t offset = 0;
  while (offset < records.size())
  {
    MFEM_VERIFY(records.size() - offset >= 4, "Truncated spatial patch-trace record!");
    PatchTrace entry;
    entry.patch = static_cast<int>(std::llround(records[offset++]));
    entry.model = static_cast<int>(std::llround(records[offset++]));
    entry.contour_size = static_cast<int>(std::llround(records[offset++]));
    const auto basis = static_cast<std::size_t>(std::llround(records[offset++]));
    MFEM_VERIFY(entry.patch >= 0 && entry.model > 0 && entry.contour_size >= 0 &&
                    static_cast<std::size_t>(entry.contour_size) <= basis &&
                    basis <= records.size() - offset,
                "Invalid gathered spatial patch-trace record!");
    entry.coefficients.assign(records.begin() + offset, records.begin() + offset + basis);
    offset += basis;
    traces.push_back(std::move(entry));
  }
  std::sort(traces.begin(), traces.end(), [](const auto &first, const auto &second)
            { return first.patch < second.patch; });
  return traces;
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
                             {"SurfaceMortar", models[i].surface_mortar},
                             {"SpatialMortar", models[i].spatial_mortar},
                             {"PatchCount", model_patch_counts[i]},
                             {"PatchWeight", model_patch_weights[i]}});
  }
  result["ModelCatalog"] = std::move(model_catalog);
  if (!ownership_diagnostics.is_null())
  {
    for (const auto &[key, value] : ownership_diagnostics.items())
    {
      result["Diagnostics"][key] = value;
    }
  }
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
