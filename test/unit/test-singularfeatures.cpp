// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "fem/singularfeatures.hpp"
#include "fem/singulargeometry.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

mfem::Mesh RectangleSheetMesh(bool reverse_elements = false, bool split_attributes = true,
                              bool duplicate_face = false)
{
  mfem::Mesh mesh(3, 6, 4, duplicate_face ? 3 : 2, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(1.0, 1.0, 0.0);
  mesh.AddVertex(0.5, 0.5, 1.0);
  mesh.AddVertex(0.5, 0.5, -1.0);

  const std::array<std::array<int, 4>, 4> elements{
      std::array<int, 4>{0, 1, 2, 4}, std::array<int, 4>{0, 2, 1, 5},
      std::array<int, 4>{1, 3, 2, 4}, std::array<int, 4>{1, 2, 3, 5}};
  if (reverse_elements)
  {
    for (auto element = elements.rbegin(); element != elements.rend(); ++element)
    {
      mesh.AddTet(element->data(), 1);
    }
  }
  else
  {
    for (const auto &element : elements)
    {
      mesh.AddTet(element.data(), 1);
    }
  }
  mesh.AddBdrTriangle(0, 1, 2, 7);
  mesh.AddBdrTriangle(1, 3, 2, split_attributes ? 8 : 7);
  if (duplicate_face)
  {
    mesh.AddBdrTriangle(0, 1, 2, 7);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh SubdividedTriangleSheetMesh()
{
  mfem::Mesh mesh(3, 6, 4, 2, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(0.5, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(0.4, 0.3, 1.0);
  mesh.AddVertex(0.4, 0.3, -1.0);
  mesh.AddTet(0, 1, 3, 4, 1);
  mesh.AddTet(0, 3, 1, 5, 1);
  mesh.AddTet(1, 2, 3, 4, 1);
  mesh.AddTet(1, 3, 2, 5, 1);
  mesh.AddBdrTriangle(0, 1, 3, 7);
  mesh.AddBdrTriangle(1, 2, 3, 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh NonmanifoldSheetEdgeMesh()
{
  mfem::Mesh mesh(3, 11, 6, 3, 3);
  mesh.AddVertex(0.0, 0.0, -1.0);
  mesh.AddVertex(0.0, 0.0, 1.0);
  constexpr double pi = 3.14159265358979323846;
  for (int branch = 0; branch < 3; branch++)
  {
    const double angle = 2.0 * pi * branch / 3.0;
    const double x = std::cos(angle);
    const double y = std::sin(angle);
    const int sheet_vertex = 2 + 3 * branch;
    const int side_0 = sheet_vertex + 1;
    const int side_1 = sheet_vertex + 2;
    mesh.AddVertex(x, y, 0.0);
    mesh.AddVertex(0.5 * x - y, 0.5 * y + x, 0.0);
    mesh.AddVertex(0.5 * x + y, 0.5 * y - x, 0.0);
    mesh.AddTet(0, 1, sheet_vertex, side_0, 1);
    mesh.AddTet(0, sheet_vertex, 1, side_1, 1);
    mesh.AddBdrTriangle(0, 1, sheet_vertex, 7);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh PerimeterJunctionMesh()
{
  mfem::Mesh mesh(3, 13, 6, 3, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  constexpr double pi = 3.14159265358979323846;
  for (int branch = 0; branch < 3; branch++)
  {
    const double angle_0 = 2.0 * pi * branch / 3.0;
    const double angle_1 = angle_0 + pi / 3.0;
    const double x_0 = std::cos(angle_0);
    const double y_0 = std::sin(angle_0);
    const double x_1 = std::cos(angle_1);
    const double y_1 = std::sin(angle_1);
    const int face_0 = 1 + 4 * branch;
    const int face_1 = face_0 + 1;
    const int side_0 = face_0 + 2;
    const int side_1 = face_0 + 3;
    mesh.AddVertex(x_0, y_0, 0.0);
    mesh.AddVertex(x_1, y_1, 0.0);
    mesh.AddVertex(0.25 * (x_0 + x_1), 0.25 * (y_0 + y_1), 1.0);
    mesh.AddVertex(0.25 * (x_0 + x_1), 0.25 * (y_0 + y_1), -1.0);
    mesh.AddTet(0, face_0, face_1, side_0, 1);
    mesh.AddTet(0, face_1, face_0, side_1, 1);
    mesh.AddBdrTriangle(0, face_0, face_1, 7);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh FiniteMetalWedgeMesh(bool heterogeneous = false)
{
  // Three right-angle tetrahedral sectors fill a 3*pi/2 dielectric fan around
  // edge (0, 1). The omitted fourth sector is the finite metal.
  mfem::Mesh mesh(3, 6, 3, 2, 3);
  mesh.AddVertex(0.0, 0.0, -1.0);
  mesh.AddVertex(0.0, 0.0, 1.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(-1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, -1.0, 0.0);

  mesh.AddTet(0, 1, 2, 3, 1);
  mesh.AddTet(0, 1, 3, 4, 1);
  mesh.AddTet(0, 1, 4, 5, heterogeneous ? 2 : 1);
  mesh.AddBdrTriangle(0, 1, 2, 7);
  mesh.AddBdrTriangle(0, 5, 1, 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh SplitFiniteMetalWedgeMesh(bool axial_exponent_transition = false)
{
  // Split the same 3*pi/2 fan into two segments along edge (0, 1, 2).
  // With axial_exponent_transition, only the upper final sector changes
  // material, so the two collinear edge segments have different exponents.
  mfem::Mesh mesh(3, 7, 6, 4, 3);
  mesh.AddVertex(0.0, 0.0, -1.0);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 0.0, 1.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(-1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, -1.0, 0.0);

  mesh.AddTet(0, 1, 3, 4, 1);
  mesh.AddTet(1, 2, 3, 4, 1);
  mesh.AddTet(0, 1, 4, 5, 1);
  mesh.AddTet(1, 2, 4, 5, 1);
  mesh.AddTet(0, 1, 5, 6, 1);
  mesh.AddTet(1, 2, 5, 6, axial_exponent_transition ? 2 : 1);
  mesh.AddBdrTriangle(0, 1, 3, 7);
  mesh.AddBdrTriangle(1, 2, 3, 7);
  mesh.AddBdrTriangle(0, 6, 1, 7);
  mesh.AddBdrTriangle(1, 6, 2, 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh InternalLineTipMesh(bool bend = false, bool selected_external_edge = false,
                               bool straight_chain = false)
{
  mfem::Mesh mesh(2, 9, 8, (bend || straight_chain) ? 2 : 1, 2);
  mesh.AddVertex(-1.0, -1.0);
  mesh.AddVertex(0.0, -1.0);
  mesh.AddVertex(1.0, -1.0);
  mesh.AddVertex(-1.0, 0.0);
  mesh.AddVertex(0.0, 0.0);
  mesh.AddVertex(1.0, 0.0);
  mesh.AddVertex(-1.0, 1.0);
  mesh.AddVertex(0.0, 1.0);
  mesh.AddVertex(1.0, 1.0);

  mesh.AddTriangle(0, 1, 4, 1);
  mesh.AddTriangle(0, 4, 3, 1);
  mesh.AddTriangle(3, 4, 7, 1);
  mesh.AddTriangle(3, 7, 6, 1);
  mesh.AddTriangle(1, 2, 5, 1);
  mesh.AddTriangle(1, 5, 4, 1);
  mesh.AddTriangle(4, 5, 8, 1);
  mesh.AddTriangle(4, 8, 7, 1);

  if (selected_external_edge)
  {
    mesh.AddBdrSegment(0, 1, 7);
  }
  else
  {
    mesh.AddBdrSegment(3, 4, 7);
    if (bend)
    {
      mesh.AddBdrSegment(4, 7, 7);
    }
    else if (straight_chain)
    {
      mesh.AddBdrSegment(4, 5, 7);
    }
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh FiniteMetalCornerMesh(bool heterogeneous = true)
{
  // The omitted upper-right quadrant is metal. The six triangles form the
  // 3*pi/2 dielectric fan between the two selected PEC rays.
  mfem::Mesh mesh(2, 8, 6, 8, 2);
  mesh.AddVertex(0.0, 0.0);
  mesh.AddVertex(1.0, 0.0);
  mesh.AddVertex(1.0, -1.0);
  mesh.AddVertex(0.0, -1.0);
  mesh.AddVertex(-1.0, -1.0);
  mesh.AddVertex(-1.0, 0.0);
  mesh.AddVertex(-1.0, 1.0);
  mesh.AddVertex(0.0, 1.0);

  mesh.AddTriangle(0, 2, 1, 1);
  mesh.AddTriangle(0, 3, 2, 1);
  mesh.AddTriangle(0, 4, 3, 1);
  mesh.AddTriangle(0, 5, 4, 1);
  mesh.AddTriangle(0, 6, 5, heterogeneous ? 2 : 1);
  mesh.AddTriangle(0, 7, 6, heterogeneous ? 2 : 1);

  mesh.AddBdrSegment(0, 1, 7);
  mesh.AddBdrSegment(7, 0, 7);
  mesh.AddBdrSegment(1, 2, 8);
  mesh.AddBdrSegment(2, 3, 8);
  mesh.AddBdrSegment(3, 4, 8);
  mesh.AddBdrSegment(4, 5, 8);
  mesh.AddBdrSegment(5, 6, 8);
  mesh.AddBdrSegment(6, 7, 8);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh StraightFiniteMetalBoundaryMesh()
{
  mfem::Mesh mesh(2, 6, 4, 6, 2);
  mesh.AddVertex(0.0, 0.0);
  mesh.AddVertex(1.0, 0.0);
  mesh.AddVertex(2.0, 0.0);
  mesh.AddVertex(0.0, 1.0);
  mesh.AddVertex(1.0, 1.0);
  mesh.AddVertex(2.0, 1.0);

  mesh.AddTriangle(0, 1, 4, 1);
  mesh.AddTriangle(0, 4, 3, 1);
  mesh.AddTriangle(1, 2, 5, 1);
  mesh.AddTriangle(1, 5, 4, 1);

  mesh.AddBdrSegment(0, 1, 7);
  mesh.AddBdrSegment(1, 2, 7);
  mesh.AddBdrSegment(2, 5, 8);
  mesh.AddBdrSegment(5, 4, 8);
  mesh.AddBdrSegment(4, 3, 8);
  mesh.AddBdrSegment(3, 0, 8);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

void SetQuadraticGeometry(mfem::Mesh &mesh, bool curve_selected_feature)
{
  const int dimension = mesh.SpaceDimension();
  mesh.SetCurvature(2, false, dimension, mfem::Ordering::byVDIM);
  mfem::VectorFunctionCoefficient geometry(
      dimension,
      [dimension, curve_selected_feature](const mfem::Vector &x, mfem::Vector &value)
      {
        value.SetSize(dimension);
        value = x;
        if (dimension == 2)
        {
          value[1] += curve_selected_feature ? 0.1 * x[0] * x[0] : 0.1 * x[0] * x[1];
        }
        else
        {
          value[2] +=
              curve_selected_feature ? 0.1 * x[0] * (1.0 - x[0]) : 0.1 * x[0] * x[2];
        }
      });
  mesh.GetNodes()->ProjectCoefficient(geometry);
}

void SetFiniteMetalWedgeGeometry(mfem::Mesh &mesh, bool varying_angle)
{
  const int order = varying_angle ? 4 : 2;
  mesh.SetCurvature(order, false, 3, mfem::Ordering::byVDIM);
  mfem::VectorFunctionCoefficient geometry(
      3,
      [varying_angle](const mfem::Vector &x, mfem::Vector &value)
      {
        value.SetSize(3);
        value = x;
        if (varying_angle)
        {
          // This quartic shear leaves the z-axis feature straight and vanishes
          // at t = 1/4, 1/2, and 3/4 along it. A three-sample angle check would
          // therefore miss the variation between those points.
          value[0] += 0.4 * x[1] * x[2] * (x[2] * x[2] - 0.25);
        }
        else
        {
          // Isotropic radial scaling curves the tetrahedra while preserving
          // every cross-sectional wedge angle along the straight z-axis.
          const double scale = 1.0 + 0.1 * x[2];
          value[0] *= scale;
          value[1] *= scale;
        }
      });
  mesh.GetNodes()->ProjectCoefficient(geometry);
}

void SetQuadraticStraightNonAffineFeatureGeometry(mfem::Mesh &mesh)
{
  const int dimension = mesh.SpaceDimension();
  mesh.SetCurvature(2, false, dimension, mfem::Ordering::byVDIM);
  mfem::VectorFunctionCoefficient geometry(
      dimension,
      [dimension](const mfem::Vector &x, mfem::Vector &value)
      {
        value.SetSize(dimension);
        value = x;
        if (dimension == 2)
        {
          // The selected edge is y = 0, -1 <= x <= 0. Its image stays on
          // y = 0, but the quadratic x parametrization is nonuniform.
          value[0] += 0.2 * x[0] * (x[0] + 1.0);
          value[1] += 0.05 * x[0] * x[1];
        }
        else
        {
          // The selected sheet perimeter remains the unit-square boundary in
          // z = 0, with nonuniform quadratic maps along its straight edges.
          value[0] += 0.15 * x[0] * (1.0 - x[0]);
          value[1] += 0.1 * x[1] * (1.0 - x[1]);
          value[2] += 0.05 * x[0] * x[2];
        }
      });
  mesh.GetNodes()->ProjectCoefficient(geometry);
}

std::map<std::array<int, 4>,
         std::pair<std::vector<std::pair<std::size_t, std::array<int, 4>>>,
                   std::vector<std::pair<std::size_t, std::array<int, 4>>>>>
GlobalElementIncidence(const mfem::Mesh &mesh,
                       const fem::singular::FeatureTopology &topology)
{
  using Incidence = std::pair<std::vector<std::pair<std::size_t, std::array<int, 4>>>,
                              std::vector<std::pair<std::size_t, std::array<int, 4>>>>;
  std::map<std::array<int, 4>, Incidence> result;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const int *vertices = mesh.GetElement(element)->GetVertices();
    std::array<int, 4> element_key{vertices[0], vertices[1], vertices[2], vertices[3]};
    std::sort(element_key.begin(), element_key.end());
    auto &[nodes, edges] = result[element_key];
    for (const auto &node : topology.elements[element].nodes)
    {
      std::array<int, 4> global;
      for (int i = 0; i < 4; i++)
      {
        global[i] = vertices[node.canonical_nodes[i]];
      }
      nodes.push_back({node.vertex, global});
    }
    for (const auto &edge : topology.elements[element].edges)
    {
      std::array<int, 4> global;
      for (int i = 0; i < 4; i++)
      {
        global[i] = vertices[edge.canonical_nodes[i]];
      }
      edges.push_back({edge.feature, global});
    }
    std::sort(nodes.begin(), nodes.end());
    std::sort(edges.begin(), edges.end());
  }
  return result;
}

void CheckCanonicalElementTuples(const mfem::Mesh &mesh,
                                 const fem::singular::FeatureTopology &topology)
{
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const int *vertices = mesh.GetElement(element)->GetVertices();
    for (const auto &node : topology.elements[element].nodes)
    {
      CAPTURE(element, node.vertex, node.mesh_vertex, node.canonical_nodes);
      CHECK(vertices[node.canonical_nodes[0]] == node.mesh_vertex);
      CHECK(vertices[node.canonical_nodes[1]] < vertices[node.canonical_nodes[2]]);
      CHECK(vertices[node.canonical_nodes[2]] < vertices[node.canonical_nodes[3]]);
    }
    for (const auto &edge : topology.elements[element].edges)
    {
      CAPTURE(element, edge.feature, edge.mesh_edge, edge.canonical_nodes);
      CHECK(vertices[edge.canonical_nodes[0]] < vertices[edge.canonical_nodes[1]]);
      CHECK(vertices[edge.canonical_nodes[2]] < vertices[edge.canonical_nodes[3]]);
      const auto &segment = topology.segments[edge.segment];
      CHECK(vertices[edge.canonical_nodes[0]] == segment.mesh_vertices[0]);
      CHECK(vertices[edge.canonical_nodes[1]] == segment.mesh_vertices[1]);
    }
  }
}

}  // namespace

TEST_CASE("Straight high-order singular segments require globally regular maps",
          "[singularfeatures][geometry][Serial]")
{
  mfem::H1_SegmentElement segment_element(3);
  mfem::IsoparametricTransformation transformation;
  transformation.SetFE(&segment_element);
  mfem::DenseMatrix points(2, segment_element.GetDof());

  const auto set_map = [&](const auto &longitudinal, const auto &transverse)
  {
    const auto &nodes = segment_element.GetNodes();
    for (int i = 0; i < nodes.GetNPoints(); i++)
    {
      const double t = nodes.IntPoint(i).x;
      points(0, i) = longitudinal(t);
      points(1, i) = transverse(t);
    }
    transformation.SetPointMat(points);
  };

  set_map([](double t) { return t + 0.15 * t * (1.0 - t); }, [](double) { return 0.0; });
  CHECK(fem::singular::IsGeometricallyStraightSegmentTransformation(transformation));

  // The derivative is strictly positive, but its degree-two Bernstein control
  // polygon has a negative middle coefficient. This requires subdivision to
  // certify and guards against rejecting a valid map from one loose bound.
  constexpr double positive_floor = 0.01;
  constexpr double positive_integral = 1.0 / 12.0 + positive_floor;
  set_map(
      [](double t)
      {
        return (((t - 0.5) * (t - 0.5) * (t - 0.5) + 0.125) / 3.0 + positive_floor * t) /
               positive_integral;
      },
      [](double) { return 0.0; });
  CHECK(fem::singular::IsGeometricallyStraightSegmentTransformation(transformation));

  // This cubic map has a negative derivative only on (0.10, 0.18), between
  // the order-eight Gauss points used by the former sampled validator.
  constexpr double center = 0.14;
  constexpr double half_width = 0.04;
  constexpr double integral =
      ((1.0 - center) * (1.0 - center) * (1.0 - center) + center * center * center) / 3.0 -
      half_width * half_width;
  set_map(
      [](double t)
      {
        const double primitive =
            ((t - center) * (t - center) * (t - center) + center * center * center) / 3.0 -
            half_width * half_width * t;
        return primitive / integral;
      },
      [](double) { return 0.0; });
  CHECK_FALSE(fem::singular::IsGeometricallyStraightSegmentTransformation(transformation));

  set_map([](double t) { return t; }, [](double t) { return 1.0e-3 * t * (1.0 - t); });
  CHECK_FALSE(fem::singular::IsGeometricallyStraightSegmentTransformation(transformation));

  set_map([](double t) { return t * t; }, [](double) { return 0.0; });
  CHECK_FALSE(fem::singular::IsGeometricallyStraightSegmentTransformation(transformation));
}

TEST_CASE("Singular line features extract internal thin-sheet tips",
          "[singularfeatures][triangle][Serial]")
{
  auto mesh = InternalLineTipMesh();
  const auto topology = fem::singular::ExtractSerialLineTipFeatures(mesh, {7}, 0.5);
  REQUIRE(topology.selected_segments.size() == 1);
  CHECK(topology.selected_segments[0].boundary_attribute == 7);
  CHECK(topology.selected_segments[0].mesh_vertices == std::array<int, 2>{3, 4});
  REQUIRE(topology.vertices.size() == 1);
  CHECK(topology.vertices[0].id == 0);
  CHECK(topology.vertices[0].mesh_vertex == 4);
  CHECK(topology.vertices[0].selected_segments == std::vector<std::size_t>{0});
  CHECK(topology.vertices[0].nu == 0.5);
  REQUIRE(topology.elements.size() == static_cast<std::size_t>(mesh.GetNE()));

  std::size_t incidence_count = 0;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const int *vertices = mesh.GetElement(element)->GetVertices();
    const bool contains_tip =
        std::find(vertices, vertices + 3, topology.vertices[0].mesh_vertex) != vertices + 3;
    CAPTURE(element);
    CHECK(topology.elements[element].nodes.size() == (contains_tip ? 1 : 0));
    for (const auto &node : topology.elements[element].nodes)
    {
      incidence_count++;
      CHECK(node.vertex == 0);
      CHECK(vertices[node.canonical_nodes[0]] == 4);
      CHECK(vertices[node.canonical_nodes[1]] < vertices[node.canonical_nodes[2]]);
    }
  }
  CHECK(incidence_count == 6);

  const auto empty = fem::singular::ExtractSerialLineTipFeatures(mesh, {});
  CHECK(empty.Empty());
  CHECK(empty.elements.size() == static_cast<std::size_t>(mesh.GetNE()));

  auto quadratic_affine = InternalLineTipMesh();
  quadratic_affine.SetCurvature(2);
  CHECK(
      fem::singular::ExtractSerialLineTipFeatures(quadratic_affine, {7}).vertices.size() ==
      1);

  auto curved_elements = InternalLineTipMesh();
  SetQuadraticGeometry(curved_elements, false);
  CHECK(fem::singular::ExtractSerialLineTipFeatures(curved_elements, {7}).vertices.size() ==
        1);

  auto nonuniform_straight_feature = InternalLineTipMesh();
  SetQuadraticStraightNonAffineFeatureGeometry(nonuniform_straight_feature);
  CHECK(fem::singular::ExtractSerialLineTipFeatures(nonuniform_straight_feature, {7})
            .vertices.size() == 1);

  auto curved_feature = InternalLineTipMesh();
  SetQuadraticGeometry(curved_feature, true);
  CHECK_THROWS_AS(fem::singular::ExtractSerialLineTipFeatures(curved_feature, {7}),
                  std::invalid_argument);

  auto endpoint_degenerate_feature = InternalLineTipMesh();
  endpoint_degenerate_feature.SetCurvature(2, false, 2, mfem::Ordering::byVDIM);
  mfem::VectorFunctionCoefficient endpoint_degenerate_geometry(
      2,
      [](const mfem::Vector &x, mfem::Vector &value)
      {
        value.SetSize(2);
        value = x;
        value[0] = (x[0] + 1.0) * (x[0] + 1.0) - 1.0;
      });
  endpoint_degenerate_feature.GetNodes()->ProjectCoefficient(endpoint_degenerate_geometry);
  endpoint_degenerate_feature.NodesUpdated();
  CHECK_THROWS_AS(
      fem::singular::ExtractSerialLineTipFeatures(endpoint_degenerate_feature, {7}),
      std::invalid_argument);

  auto physically_bent_chain = InternalLineTipMesh(false, false, true);
  physically_bent_chain.SetCurvature(2, false, 2, mfem::Ordering::byVDIM);
  mfem::VectorFunctionCoefficient bent_geometry(
      2,
      [](const mfem::Vector &x, mfem::Vector &value)
      {
        value.SetSize(2);
        value = x;
        value[1] += 0.2 * (1.0 - std::abs(x[0]));
      });
  physically_bent_chain.GetNodes()->ProjectCoefficient(bent_geometry);
  physically_bent_chain.NodesUpdated();
  CHECK_THROWS_AS(fem::singular::ExtractSerialLineTipFeatures(physically_bent_chain, {7}),
                  std::invalid_argument);
}

TEST_CASE("Singular line feature extraction rejects unsupported topology",
          "[singularfeatures][triangle][Serial]")
{
  {
    auto mesh = InternalLineTipMesh(true);
    CHECK_THROWS_AS(fem::singular::ExtractSerialLineTipFeatures(mesh, {7}),
                    std::invalid_argument);
  }
  {
    auto mesh = InternalLineTipMesh(false, true);
    CHECK_THROWS_AS(fem::singular::ExtractSerialLineTipFeatures(mesh, {7}),
                    std::invalid_argument);
  }
  {
    auto mesh = InternalLineTipMesh();
    CHECK_THROWS_AS(fem::singular::ExtractSerialLineTipFeatures(mesh, {8}),
                    std::invalid_argument);
    CHECK_THROWS_AS(fem::singular::ExtractSerialLineTipFeatures(mesh, {0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(fem::singular::ExtractSerialLineTipFeatures(mesh, {7}, 1.0),
                    std::invalid_argument);
  }
}

TEST_CASE("Dirichlet material-wedge exponents recover homogeneous and transmission roots",
          "[singularfeatures][triangle][Serial]")
{
  constexpr double pi = 3.14159265358979323846;
  CHECK(fem::singular::ComputeDirichletWedgeExponent({{1, 1.5 * pi, 1.0}}) ==
        Catch::Approx(2.0 / 3.0).epsilon(2.0e-13));

  const std::vector<fem::singular::TriangleWedgeSector> sectors{{1, pi, 11.47},
                                                                {2, 0.5 * pi, 1.0}};
  const double nu = fem::singular::ComputeDirichletWedgeExponent(sectors);
  CHECK(nu == Catch::Approx(0.5255535).epsilon(2.0e-7));
  CHECK(std::abs(11.47 / std::tan(pi * nu) + 1.0 / std::tan(0.5 * pi * nu)) < 1.0e-11);

  CHECK(fem::singular::ComputeDirichletWedgeExponent({{1, pi, 1.0}}) == 1.0);
  CHECK_THROWS_AS(fem::singular::ComputeDirichletWedgeExponent({}), std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ComputeDirichletWedgeExponent({{1, 0.5 * pi, -1.0}}),
                  std::invalid_argument);
}

TEST_CASE("Singular line features extract one-sided finite-metal corners",
          "[singularfeatures][triangle][Serial]")
{
  constexpr double pi = 3.14159265358979323846;
  const std::vector<fem::singular::TriangleMaterial> materials{{1, 11.47}, {2, 1.0}};
  auto mesh = FiniteMetalCornerMesh();
  const auto topology = fem::singular::ExtractSerialLineFeatures(mesh, {7}, materials);
  REQUIRE(topology.selected_segments.size() == 2);
  REQUIRE(topology.vertices.size() == 1);
  const auto &corner = topology.vertices[0];
  CHECK(corner.mesh_vertex == 0);
  CHECK(corner.type == fem::singular::FeatureVertexType::CORNER);
  REQUIRE(corner.selected_segments.size() == 2);
  REQUIRE(corner.sectors.size() == 2);
  CHECK(corner.sectors[0].angle + corner.sectors[1].angle ==
        Catch::Approx(1.5 * pi).epsilon(2.0e-14));
  CHECK(corner.nu == Catch::Approx(0.5255535).epsilon(2.0e-7));

  std::size_t incidence_count = 0;
  for (const auto &element : topology.elements)
  {
    incidence_count += element.nodes.size();
  }
  CHECK(incidence_count == 6);

  auto homogeneous_mesh = FiniteMetalCornerMesh(false);
  const auto homogeneous =
      fem::singular::ExtractSerialLineFeatures(homogeneous_mesh, {7}, {{1, 1.0}});
  REQUIRE(homogeneous.vertices.size() == 1);
  CHECK(homogeneous.vertices[0].nu == Catch::Approx(2.0 / 3.0).epsilon(2.0e-13));
  REQUIRE(homogeneous.vertices[0].sectors.size() == 1);
  CHECK(homogeneous.vertices[0].sectors[0].angle ==
        Catch::Approx(1.5 * pi).epsilon(2.0e-14));

  auto straight_mesh = StraightFiniteMetalBoundaryMesh();
  const auto straight =
      fem::singular::ExtractSerialLineFeatures(straight_mesh, {7}, {{1, 1.0}});
  CHECK(straight.selected_segments.size() == 2);
  CHECK(straight.vertices.empty());

  CHECK_THROWS_AS(fem::singular::ExtractSerialLineFeatures(mesh, {7}, {{1, 11.47}}),
                  std::invalid_argument);
}

TEST_CASE("Serial singular line-tip blueprints broadcast sparsely",
          "[singularfeatures][triangle][Parallel]")
{
  auto mesh = InternalLineTipMesh();
  const auto expected = fem::singular::ExtractSerialLineTipFeatures(mesh, {7});
  auto broadcast =
      Mpi::Root(Mpi::World()) ? expected : fem::singular::TriangleFeatureTopology{};
  fem::singular::BroadcastSerialLineTipFeatures(broadcast, Mpi::World());

  REQUIRE(broadcast.vertices.size() == expected.vertices.size());
  REQUIRE(broadcast.selected_segments.size() == expected.selected_segments.size());
  REQUIRE(broadcast.elements.size() == expected.elements.size());
  for (std::size_t i = 0; i < expected.vertices.size(); i++)
  {
    CHECK(broadcast.vertices[i].id == expected.vertices[i].id);
    CHECK(broadcast.vertices[i].mesh_vertex == expected.vertices[i].mesh_vertex);
    CHECK(broadcast.vertices[i].selected_segments ==
          expected.vertices[i].selected_segments);
    CHECK(broadcast.vertices[i].nu == expected.vertices[i].nu);
  }
  for (std::size_t i = 0; i < expected.selected_segments.size(); i++)
  {
    CHECK(broadcast.selected_segments[i].boundary_element ==
          expected.selected_segments[i].boundary_element);
    CHECK(broadcast.selected_segments[i].mesh_edge ==
          expected.selected_segments[i].mesh_edge);
    CHECK(broadcast.selected_segments[i].mesh_vertices ==
          expected.selected_segments[i].mesh_vertices);
    CHECK(broadcast.selected_segments[i].boundary_attribute ==
          expected.selected_segments[i].boundary_attribute);
  }
  for (std::size_t element = 0; element < expected.elements.size(); element++)
  {
    REQUIRE(broadcast.elements[element].nodes.size() ==
            expected.elements[element].nodes.size());
    for (std::size_t node = 0; node < expected.elements[element].nodes.size(); node++)
    {
      CHECK(broadcast.elements[element].nodes[node].vertex ==
            expected.elements[element].nodes[node].vertex);
      CHECK(broadcast.elements[element].nodes[node].mesh_vertex ==
            expected.elements[element].nodes[node].mesh_vertex);
      CHECK(broadcast.elements[element].nodes[node].canonical_nodes ==
            expected.elements[element].nodes[node].canonical_nodes);
    }
  }
}

TEST_CASE("Finite-metal corner blueprints preserve material sectors across broadcast",
          "[singularfeatures][triangle][Parallel]")
{
  auto mesh = FiniteMetalCornerMesh();
  const auto expected =
      fem::singular::ExtractSerialLineFeatures(mesh, {7}, {{1, 11.47}, {2, 1.0}});
  auto broadcast =
      Mpi::Root(Mpi::World()) ? expected : fem::singular::TriangleFeatureTopology{};
  fem::singular::BroadcastSerialLineTipFeatures(broadcast, Mpi::World());

  REQUIRE(broadcast.vertices.size() == 1);
  REQUIRE(expected.vertices.size() == 1);
  const auto &actual_corner = broadcast.vertices[0];
  const auto &expected_corner = expected.vertices[0];
  CHECK(actual_corner.type == fem::singular::FeatureVertexType::CORNER);
  CHECK(actual_corner.selected_segments == expected_corner.selected_segments);
  CHECK(actual_corner.nu == expected_corner.nu);
  REQUIRE(actual_corner.sectors.size() == expected_corner.sectors.size());
  for (std::size_t i = 0; i < expected_corner.sectors.size(); i++)
  {
    CHECK(actual_corner.sectors[i].domain_attribute ==
          expected_corner.sectors[i].domain_attribute);
    CHECK(actual_corner.sectors[i].angle == expected_corner.sectors[i].angle);
    CHECK(actual_corner.sectors[i].permittivity == expected_corner.sectors[i].permittivity);
  }
}

TEST_CASE("Singular sheet features extract topological perimeter and element incidence",
          "[singularfeatures][Serial]")
{
  auto mesh = RectangleSheetMesh();
  const auto topology = fem::singular::ExtractSerialSheetFeatures(mesh, {8, 7}, 0.5);
  REQUIRE(topology.segments.size() == 4);
  REQUIRE(topology.vertices.size() == 4);
  REQUIRE(topology.features.size() == 4);
  REQUIRE(topology.sheet_faces.size() == 2);
  REQUIRE(topology.elements.size() == 4);
  CHECK(topology.sheet_faces[0].mesh_vertices ==
        std::array<fem::singular::GlobalVertexId, 3>{0, 1, 2});
  CHECK(topology.sheet_faces[0].boundary_attribute == 7);
  CHECK(topology.sheet_faces[1].mesh_vertices ==
        std::array<fem::singular::GlobalVertexId, 3>{1, 2, 3});
  CHECK(topology.sheet_faces[1].boundary_attribute == 8);

  const std::array<std::array<int, 2>, 4> expected_edges{
      std::array<int, 2>{0, 1}, std::array<int, 2>{0, 2}, std::array<int, 2>{1, 3},
      std::array<int, 2>{2, 3}};
  for (std::size_t i = 0; i < expected_edges.size(); i++)
  {
    CAPTURE(i);
    CHECK(topology.segments[i].mesh_vertices == expected_edges[i]);
    CHECK(topology.segments[i].feature == i);
    CHECK(topology.features[i].id == i);
    CHECK(topology.features[i].segments == std::vector<std::size_t>{i});
    CHECK(topology.features[i].nu == 0.5);
    CHECK_FALSE(topology.features[i].closed);
  }
  for (const auto &vertex : topology.vertices)
  {
    CHECK(vertex.type == fem::singular::FeatureVertexType::CORNER);
    CHECK(vertex.segments.size() == 2);
    CHECK(vertex.features.size() == 2);
    CHECK(vertex.id < topology.vertices.size());
    CHECK(vertex.nu == 0.5);
  }

  std::size_t node_incidence = 0;
  std::size_t edge_incidence = 0;
  for (const auto &element : topology.elements)
  {
    node_incidence += element.nodes.size();
    edge_incidence += element.edges.size();
    std::set<int> singular_vertices;
    for (const auto &node : element.nodes)
    {
      CHECK(singular_vertices.insert(node.mesh_vertex).second);
      CHECK(topology.vertices[node.vertex].mesh_vertex == node.mesh_vertex);
    }
    for (const auto &edge : element.edges)
    {
      const auto &segment = topology.segments[edge.segment];
      CHECK(singular_vertices.count(segment.mesh_vertices[0]) == 1);
      CHECK(singular_vertices.count(segment.mesh_vertices[1]) == 1);
    }
  }
  CHECK(node_incidence == 12);
  CHECK(edge_incidence == 8);
  CheckCanonicalElementTuples(mesh, topology);

  const auto one_triangle = fem::singular::ExtractSerialSheetFeatures(mesh, {7});
  REQUIRE(one_triangle.segments.size() == 3);
  CHECK(std::find_if(one_triangle.segments.begin(), one_triangle.segments.end(),
                     [](const auto &segment)
                     { return segment.mesh_vertices == std::array<int, 2>{1, 2}; }) !=
        one_triangle.segments.end());
}

TEST_CASE("Singular sheet features group collinear mesh edges into straight features",
          "[singularfeatures][Serial]")
{
  auto mesh = SubdividedTriangleSheetMesh();
  const auto topology = fem::singular::ExtractSerialSheetFeatures(mesh, {7}, 2.0 / 3.0);
  REQUIRE(topology.segments.size() == 4);
  REQUIRE(topology.features.size() == 3);

  const auto vertex =
      std::find_if(topology.vertices.begin(), topology.vertices.end(),
                   [](const auto &entry) { return entry.mesh_vertex == 1; });
  REQUIRE(vertex != topology.vertices.end());
  CHECK(vertex->type == fem::singular::FeatureVertexType::REGULAR);
  REQUIRE(vertex->features.size() == 1);
  const auto &straight = topology.features[vertex->features[0]];
  CHECK(straight.segments.size() == 2);
  CHECK(straight.mesh_vertices == std::vector<int>{0, 1, 2});
  CHECK(straight.nu == 2.0 / 3.0);
  CHECK_FALSE(straight.closed);
}

TEST_CASE("Finite-metal wedge features infer homogeneous and transmission exponents",
          "[singularfeatures][Serial]")
{
  constexpr double pi = 3.14159265358979323846;
  auto homogeneous_mesh = FiniteMetalWedgeMesh();
  const auto homogeneous = fem::singular::ExtractSerialSheetFeatures(
      homogeneous_mesh, {7}, std::vector<fem::singular::TriangleMaterial>{{1, 1.0}});
  REQUIRE(homogeneous.sheet_faces.size() == 2);
  REQUIRE(homogeneous.segments.size() == 1);
  REQUIRE(homogeneous.features.size() == 1);
  REQUIRE(homogeneous.vertices.size() == 2);
  CHECK(homogeneous.segments[0].mesh_vertices == std::array<int, 2>{0, 1});
  CHECK(homogeneous.segments[0].type ==
        fem::singular::FeatureSegmentType::TRANSMISSION_WEDGE);
  CHECK(homogeneous.features[0].nu == Catch::Approx(2.0 / 3.0).margin(2.0e-13));
  CHECK(homogeneous.features[0].segments == std::vector<std::size_t>{0});
  CHECK_FALSE(homogeneous.features[0].closed);
  for (const auto &vertex : homogeneous.vertices)
  {
    CHECK(vertex.type == fem::singular::FeatureVertexType::ENDPOINT);
    CHECK(vertex.nu == Catch::Approx(2.0 / 3.0).margin(2.0e-13));
  }
  for (const auto &element : homogeneous.elements)
  {
    REQUIRE(element.edges.size() == 1);
    REQUIRE(element.nodes.size() == 2);
    CHECK(element.edges[0].segment == 0);
  }

  auto curved_constant_angle_mesh = FiniteMetalWedgeMesh();
  SetFiniteMetalWedgeGeometry(curved_constant_angle_mesh, false);
  const auto curved_constant_angle = fem::singular::ExtractSerialSheetFeatures(
      curved_constant_angle_mesh, {7},
      std::vector<fem::singular::TriangleMaterial>{{1, 1.0}});
  REQUIRE(curved_constant_angle.features.size() == 1);
  CHECK(curved_constant_angle.features[0].nu == Catch::Approx(2.0 / 3.0).margin(2.0e-12));

  auto curved_varying_angle_mesh = FiniteMetalWedgeMesh();
  SetFiniteMetalWedgeGeometry(curved_varying_angle_mesh, true);
  CHECK_THROWS_WITH(
      fem::singular::ExtractSerialSheetFeatures(
          curved_varying_angle_mesh, {7},
          std::vector<fem::singular::TriangleMaterial>{{1, 1.0}}),
      Catch::Matchers::ContainsSubstring("cross-sectional angle which varies"));

  auto heterogeneous_mesh = FiniteMetalWedgeMesh(true);
  const std::vector<fem::singular::TriangleMaterial> materials{{1, 11.47}, {2, 1.0}};
  const auto heterogeneous =
      fem::singular::ExtractSerialSheetFeatures(heterogeneous_mesh, {7}, materials);
  REQUIRE(heterogeneous.features.size() == 1);
  const double expected =
      fem::singular::ComputeDirichletWedgeExponent({{1, pi, 11.47}, {2, 0.5 * pi, 1.0}});
  CHECK(heterogeneous.features[0].nu == Catch::Approx(expected).margin(2.0e-13));
  CHECK(heterogeneous.features[0].nu < 2.0 / 3.0);

  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(
                      heterogeneous_mesh, {7},
                      std::vector<fem::singular::TriangleMaterial>{{1, 11.47}}),
                  std::invalid_argument);

  auto split_mesh = SplitFiniteMetalWedgeMesh();
  const auto split = fem::singular::ExtractSerialSheetFeatures(
      split_mesh, {7}, std::vector<fem::singular::TriangleMaterial>{{1, 1.0}});
  REQUIRE(split.segments.size() == 2);
  REQUIRE(split.features.size() == 1);
  CHECK(split.features[0].segments == std::vector<std::size_t>{0, 1});
  CHECK(split.features[0].nu == Catch::Approx(2.0 / 3.0).margin(2.0e-13));

  auto transition_mesh = SplitFiniteMetalWedgeMesh(true);
  CHECK_THROWS_WITH(fem::singular::ExtractSerialSheetFeatures(
                        transition_mesh, {7},
                        std::vector<fem::singular::TriangleMaterial>{{1, 11.47}, {2, 1.0}}),
                    Catch::Matchers::ContainsSubstring("point-transition enrichment"));
}

TEST_CASE("Singular sheet feature IDs and basis tuples ignore element ordering",
          "[singularfeatures][Serial]")
{
  auto forward_mesh = RectangleSheetMesh(false);
  auto reverse_mesh = RectangleSheetMesh(true);
  const auto forward = fem::singular::ExtractSerialSheetFeatures(forward_mesh, {7, 8});
  const auto reverse = fem::singular::ExtractSerialSheetFeatures(reverse_mesh, {7, 8});

  REQUIRE(reverse.segments.size() == forward.segments.size());
  REQUIRE(reverse.features.size() == forward.features.size());
  for (std::size_t i = 0; i < forward.segments.size(); i++)
  {
    CHECK(reverse.segments[i].mesh_vertices == forward.segments[i].mesh_vertices);
    CHECK(reverse.segments[i].feature == forward.segments[i].feature);
    CHECK(reverse.segments[i].boundary_attributes ==
          forward.segments[i].boundary_attributes);
  }
  for (std::size_t i = 0; i < forward.features.size(); i++)
  {
    CHECK(reverse.features[i].segments == forward.features[i].segments);
    CHECK(reverse.features[i].mesh_vertices == forward.features[i].mesh_vertices);
  }
  CHECK(GlobalElementIncidence(reverse_mesh, reverse) ==
        GlobalElementIncidence(forward_mesh, forward));
}

TEST_CASE("Singular sheet feature extraction rejects unsupported inputs",
          "[singularfeatures][Serial]")
{
  auto mesh = RectangleSheetMesh();
  const auto empty = fem::singular::ExtractSerialSheetFeatures(mesh, {});
  CHECK(empty.Empty());
  CHECK(empty.vertices.empty());
  CHECK(empty.segments.empty());
  CHECK(empty.sheet_faces.empty());
  CHECK(empty.elements.size() == static_cast<std::size_t>(mesh.GetNE()));
  for (const auto &element : empty.elements)
  {
    CHECK(element.nodes.empty());
    CHECK(element.edges.empty());
  }

  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(mesh, {0}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(mesh, {99}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(mesh, {7}, 0.0),
                  std::invalid_argument);

  auto external = SingleTetMesh();
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(external, {1}),
                  std::invalid_argument);

  auto quadratic_affine = RectangleSheetMesh();
  quadratic_affine.SetCurvature(2);
  const auto quadratic_topology =
      fem::singular::ExtractSerialSheetFeatures(quadratic_affine, {7, 8});
  CHECK(quadratic_topology.segments.size() == 4);
  CHECK(quadratic_topology.elements.size() ==
        static_cast<std::size_t>(quadratic_affine.GetNE()));

  auto curved_elements = RectangleSheetMesh();
  SetQuadraticGeometry(curved_elements, false);
  CHECK(
      fem::singular::ExtractSerialSheetFeatures(curved_elements, {7, 8}).segments.size() ==
      4);

  auto nonuniform_straight_feature = RectangleSheetMesh();
  SetQuadraticStraightNonAffineFeatureGeometry(nonuniform_straight_feature);
  CHECK(fem::singular::ExtractSerialSheetFeatures(nonuniform_straight_feature, {7, 8})
            .segments.size() == 4);

  auto curved_feature = RectangleSheetMesh();
  SetQuadraticGeometry(curved_feature, true);
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(curved_feature, {7, 8}),
                  std::invalid_argument);

  auto hexahedron =
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(hexahedron, {1}),
                  std::invalid_argument);

  auto duplicate = RectangleSheetMesh(false, true, true);
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(duplicate, {7}),
                  std::invalid_argument);

  auto nonmanifold = NonmanifoldSheetEdgeMesh();
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(nonmanifold, {7}),
                  std::invalid_argument);

  auto junction = PerimeterJunctionMesh();
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(junction, {7}),
                  std::invalid_argument);
}

TEST_CASE("Serial singular sheet extraction rejects a multi-rank ParMesh",
          "[singularfeatures][Parallel]")
{
  auto serial_mesh = RectangleSheetMesh();
  mfem::ParMesh parallel_mesh(Mpi::World(), serial_mesh);
  if (parallel_mesh.GetNRanks() == 1)
  {
    SUCCEED("Multi-rank rejection is exercised by the [Parallel] test run.");
    return;
  }
  CHECK_THROWS_AS(fem::singular::ExtractSerialSheetFeatures(parallel_mesh, {7, 8}),
                  std::invalid_argument);
}

TEST_CASE("Serial singular feature blueprints broadcast sparsely",
          "[singularfeatures][Parallel]")
{
  auto mesh = RectangleSheetMesh();
  const auto expected = fem::singular::ExtractSerialSheetFeatures(mesh, {7, 8});
  auto broadcast = Mpi::Root(Mpi::World()) ? expected : fem::singular::FeatureTopology{};
  fem::singular::BroadcastSerialSheetFeatures(broadcast, Mpi::World());

  CHECK(broadcast.vertices.size() == expected.vertices.size());
  CHECK(broadcast.segments.size() == expected.segments.size());
  CHECK(broadcast.features.size() == expected.features.size());
  CHECK(broadcast.sheet_faces.size() == expected.sheet_faces.size());
  CHECK(broadcast.elements.size() == expected.elements.size());
  CHECK(GlobalElementIncidence(mesh, broadcast) == GlobalElementIncidence(mesh, expected));
  for (std::size_t i = 0; i < expected.segments.size(); i++)
  {
    CHECK(broadcast.segments[i].mesh_edge == expected.segments[i].mesh_edge);
    CHECK(broadcast.segments[i].mesh_vertices == expected.segments[i].mesh_vertices);
    CHECK(broadcast.segments[i].feature == expected.segments[i].feature);
    CHECK(broadcast.segments[i].boundary_attributes ==
          expected.segments[i].boundary_attributes);
    CHECK(broadcast.segments[i].type == expected.segments[i].type);
    CHECK(broadcast.segments[i].type == fem::singular::FeatureSegmentType::SHEET_EDGE);
  }
  for (std::size_t i = 0; i < expected.vertices.size(); i++)
  {
    CHECK(broadcast.vertices[i].id == expected.vertices[i].id);
    CHECK(broadcast.vertices[i].mesh_vertex == expected.vertices[i].mesh_vertex);
    CHECK(broadcast.vertices[i].type == expected.vertices[i].type);
    CHECK(broadcast.vertices[i].segments == expected.vertices[i].segments);
    CHECK(broadcast.vertices[i].features == expected.vertices[i].features);
    CHECK(broadcast.vertices[i].nu == expected.vertices[i].nu);
  }
  for (std::size_t i = 0; i < expected.sheet_faces.size(); i++)
  {
    CHECK(broadcast.sheet_faces[i].mesh_face == expected.sheet_faces[i].mesh_face);
    CHECK(broadcast.sheet_faces[i].mesh_vertices == expected.sheet_faces[i].mesh_vertices);
    CHECK(broadcast.sheet_faces[i].boundary_attribute ==
          expected.sheet_faces[i].boundary_attribute);
  }
}

}  // namespace palace
