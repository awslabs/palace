// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/mesh.hpp"
#include "fem/singularfield.hpp"
#include "models/materialoperator.hpp"
#include "models/singularsurfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

using namespace Catch::Matchers;

namespace palace
{

namespace
{

std::unique_ptr<mfem::Mesh> MakeSurfaceTestMesh()
{
  auto mesh = std::make_unique<mfem::Mesh>(2, 3, 1, 3, 2);
  const std::array<fem::singular::Vector2, 3> vertices{fem::singular::Vector2{0.0, 0.0},
                                                       fem::singular::Vector2{1.0, 0.0},
                                                       fem::singular::Vector2{0.0, 1.0}};
  for (const auto &vertex : vertices)
  {
    mesh->AddVertex(vertex.data());
  }
  mesh->AddTriangle(0, 1, 2, 1);
  mesh->AddBdrSegment(0, 1, 1);
  mesh->AddBdrSegment(1, 2, 1);
  mesh->AddBdrSegment(2, 0, 1);
  mesh->FinalizeTopology();
  mesh->Finalize(false, false);
  return mesh;
}

std::unique_ptr<mfem::Mesh> MakePartitionedSurfaceTestMesh()
{
  auto mesh = std::make_unique<mfem::Mesh>(2, 4, 2, 5, 2);
  mesh->AddVertex(0.0, 0.0);
  mesh->AddVertex(1.0, 0.0);
  mesh->AddVertex(0.0, 1.0);
  mesh->AddVertex(0.0, -1.0);

  // The first triangle owns the internal boundary after partitioning.
  mesh->AddTriangle(0, 1, 2, 2);
  mesh->AddTriangle(0, 3, 1, 1);
  mesh->AddBdrSegment(0, 1, 1);
  mesh->AddBdrSegment(1, 2, 2);
  mesh->AddBdrSegment(2, 0, 2);
  mesh->AddBdrSegment(0, 3, 2);
  mesh->AddBdrSegment(3, 1, 2);
  mesh->FinalizeTopology();
  mesh->Finalize(false, false);
  return mesh;
}

std::unique_ptr<mfem::Mesh>
MakeSplitCutoffSurfaceTestMesh(const std::vector<double> &breakpoints)
{
  if (breakpoints.size() < 2 || breakpoints.front() != 0.0 || breakpoints.back() != 1.0 ||
      !std::is_sorted(breakpoints.begin(), breakpoints.end()) ||
      std::adjacent_find(breakpoints.begin(), breakpoints.end()) != breakpoints.end())
  {
    throw std::invalid_argument("Invalid split surface test mesh breakpoints!");
  }
  const int intervals = static_cast<int>(breakpoints.size()) - 1;
  auto mesh = std::make_unique<mfem::Mesh>(2, 2 * (intervals + 1), 2 * intervals,
                                           2 * intervals + 2, 2);
  for (double x : breakpoints)
  {
    mesh->AddVertex(x, 0.0);
  }
  for (double x : breakpoints)
  {
    mesh->AddVertex(x, 1.0);
  }
  const int top_offset = intervals + 1;
  for (int interval = 0; interval < intervals; interval++)
  {
    const int bottom_left = interval;
    const int bottom_right = interval + 1;
    const int top_left = top_offset + interval;
    const int top_right = top_offset + interval + 1;
    mesh->AddTriangle(bottom_left, bottom_right, top_right, 1);
    mesh->AddTriangle(bottom_left, top_right, top_left, 1);
    mesh->AddBdrSegment(bottom_left, bottom_right, 1);
    mesh->AddBdrSegment(top_right, top_left, 2);
  }
  mesh->AddBdrSegment(top_offset, 0, 2);
  mesh->AddBdrSegment(intervals, top_offset + intervals, 2);
  mesh->FinalizeTopology();
  mesh->Finalize(false, false);
  return mesh;
}

std::unique_ptr<mfem::Mesh> MakeTetrahedronSurfaceTestMesh()
{
  auto mesh = std::make_unique<mfem::Mesh>(3, 4, 1, 4, 3);
  mesh->AddVertex(0.0, 0.0, 0.0);
  mesh->AddVertex(1.0, 0.0, 0.0);
  mesh->AddVertex(0.0, 1.0, 0.0);
  mesh->AddVertex(0.0, 0.0, 1.0);
  mesh->AddTet(0, 1, 2, 3, 1);
  mesh->AddBdrTriangle(1, 2, 3, 2);
  mesh->AddBdrTriangle(0, 3, 2, 2);
  mesh->AddBdrTriangle(0, 1, 3, 2);
  mesh->AddBdrTriangle(0, 2, 1, 1);
  mesh->FinalizeTopology();
  mesh->Finalize(true, false);
  return mesh;
}

std::unique_ptr<mfem::Mesh> MakePartitionedTetrahedronSurfaceTestMesh()
{
  auto mesh = std::make_unique<mfem::Mesh>(3, 5, 2, 7, 3);
  mesh->AddVertex(0.0, 0.0, 0.0);
  mesh->AddVertex(1.0, 0.0, 0.0);
  mesh->AddVertex(0.0, 1.0, 0.0);
  mesh->AddVertex(0.0, 0.0, 1.0);
  mesh->AddVertex(0.0, 0.0, -1.0);

  // The internal boundary orientation assigns it to the first (vacuum)
  // tetrahedron. The MS-selected substrate trace is therefore off-rank for
  // the partition used by the parallel test below.
  mesh->AddTet(0, 1, 2, 3, 2);
  mesh->AddTet(0, 2, 1, 4, 1);
  mesh->AddBdrTriangle(0, 2, 1, 1);
  mesh->AddBdrTriangle(1, 2, 3, 2);
  mesh->AddBdrTriangle(0, 3, 2, 2);
  mesh->AddBdrTriangle(0, 1, 3, 2);
  mesh->AddBdrTriangle(2, 1, 4, 2);
  mesh->AddBdrTriangle(0, 4, 1, 2);
  mesh->AddBdrTriangle(0, 2, 4, 2);
  mesh->FinalizeTopology();
  mesh->Finalize(true, false);
  return mesh;
}

config::BoundaryPostData MakeMSPostprocessing(double edge_cutoff = 0.0)
{
  config::InterfaceDielectricData data;
  data.type = InterfaceDielectric::MS;
  data.t = 0.01;
  data.epsilon_r = 2.0;
  data.tandelta = 0.125;
  data.edge_cutoff = edge_cutoff;
  data.attributes = {1};
  config::BoundaryPostData postpro;
  postpro.dielectric.emplace(7, data);
  return postpro;
}

config::BoundaryPostData MakeDefaultPostprocessing(double edge_cutoff = 0.0)
{
  config::InterfaceDielectricData data;
  data.type = InterfaceDielectric::DEFAULT;
  data.t = 0.01;
  data.epsilon_r = 2.0;
  data.tandelta = 0.125;
  data.edge_cutoff = edge_cutoff;
  data.attributes = {1};
  config::BoundaryPostData postpro;
  postpro.dielectric.emplace(7, data);
  return postpro;
}

config::BoundaryPostData MakeMSMAPostprocessing(bool include_ms, bool include_ma,
                                                double edge_cutoff = 0.0)
{
  config::BoundaryPostData postpro;
  config::InterfaceDielectricData data;
  data.t = 0.01;
  data.epsilon_r = 2.0;
  data.tandelta = 0.125;
  data.edge_cutoff = edge_cutoff;
  data.attributes = {1};
  if (include_ms)
  {
    data.type = InterfaceDielectric::MS;
    postpro.dielectric.emplace(7, data);
  }
  if (include_ma)
  {
    data.type = InterfaceDielectric::MA;
    postpro.dielectric.emplace(8, data);
  }
  return postpro;
}

void SetCurvedTetrahedronGeometry(mfem::Mesh &mesh, double curvature)
{
  mesh.SetCurvature(2, false, 3, mfem::Ordering::byVDIM);
  mfem::VectorFunctionCoefficient geometry(
      3,
      [curvature](const mfem::Vector &point, mfem::Vector &value)
      {
        value.SetSize(3);
        value[0] = point[0];
        value[1] = point[1];
        value[2] = point[2] + curvature * point[0] * point[1];
      });
  mesh.GetNodes()->ProjectCoefficient(geometry);
  mesh.NodesUpdated();
}

void SetBubbleCurvedTetrahedronGeometry(mfem::Mesh &mesh, double curvature)
{
  mesh.SetCurvature(3, false, 3, mfem::Ordering::byVDIM);
  mfem::VectorFunctionCoefficient geometry(
      3,
      [curvature](const mfem::Vector &point, mfem::Vector &value)
      {
        value.SetSize(3);
        value[0] = point[0];
        value[1] = point[1];
        value[2] = point[2] + curvature * point[0] * point[1] * (1.0 - point[0] - point[1]);
      });
  mesh.GetNodes()->ProjectCoefficient(geometry);
  mesh.NodesUpdated();
}

double CurvedFaceAreaOutsideCutoff(double curvature, double cutoff)
{
  const auto retained_height = [=](double x)
  { return 1.0 - x - cutoff / std::sqrt(1.0 + curvature * curvature * x * x); };
  REQUIRE(retained_height(0.0) > 0.0);
  REQUIRE(retained_height(1.0) < 0.0);
  double lower = 0.0;
  double upper = 1.0;
  for (int iteration = 0; iteration < 80; iteration++)
  {
    const double midpoint = 0.5 * (lower + upper);
    if (retained_height(midpoint) > 0.0)
    {
      lower = midpoint;
    }
    else
    {
      upper = midpoint;
    }
  }
  const double x_upper = 0.5 * (lower + upper);
  const auto rule = fem::singular::BuildWeightedSegmentQuadrature(64, 0.0, 0.0);
  long double area = 0.0L;
  for (const auto &x_quadrature : rule)
  {
    const double x = x_upper * x_quadrature.coordinate;
    const double y_lower = cutoff / std::sqrt(1.0 + curvature * curvature * x * x);
    const double y_upper = 1.0 - x;
    for (const auto &y_quadrature : rule)
    {
      const double y = y_lower + (y_upper - y_lower) * y_quadrature.coordinate;
      const double surface_jacobian =
          std::sqrt(1.0 + curvature * curvature * (x * x + y * y));
      area += x_quadrature.weight * y_quadrature.weight * x_upper * (y_upper - y_lower) *
              surface_jacobian;
    }
  }
  return static_cast<double>(area);
}

double CurvedFaceAreaOutsideTwoCutoffs(double curvature, double cutoff)
{
  const auto retained_height = [=](double x)
  {
    if (!(x > 0.0))
    {
      return -std::numeric_limits<double>::infinity();
    }
    const double first_cutoff = cutoff / std::sqrt(1.0 + curvature * curvature * x * x);
    double second_cutoff = 0.0;
    if (x < cutoff)
    {
      if (curvature == 0.0)
      {
        return -std::numeric_limits<double>::infinity();
      }
      second_cutoff =
          std::sqrt((cutoff * cutoff / (x * x) - 1.0) / (curvature * curvature));
    }
    return 1.0 - x - std::max(first_cutoff, second_cutoff);
  };

  constexpr int search_intervals = 4096;
  double lower = 0.0;
  double upper = 0.0;
  bool inside = false;
  for (int interval = 1; interval <= search_intervals; interval++)
  {
    const double x = static_cast<double>(interval) / search_intervals;
    const bool retained = retained_height(x) > 0.0;
    if (retained && !inside)
    {
      double left = static_cast<double>(interval - 1) / search_intervals;
      double right = x;
      for (int iteration = 0; iteration < 80; iteration++)
      {
        const double midpoint = 0.5 * (left + right);
        if (retained_height(midpoint) > 0.0)
        {
          right = midpoint;
        }
        else
        {
          left = midpoint;
        }
      }
      lower = 0.5 * (left + right);
      inside = true;
    }
    if (!retained && inside)
    {
      double left = static_cast<double>(interval - 1) / search_intervals;
      double right = x;
      for (int iteration = 0; iteration < 80; iteration++)
      {
        const double midpoint = 0.5 * (left + right);
        if (retained_height(midpoint) > 0.0)
        {
          left = midpoint;
        }
        else
        {
          right = midpoint;
        }
      }
      upper = 0.5 * (left + right);
      break;
    }
  }
  REQUIRE(inside);
  REQUIRE(upper > lower);

  const auto cutoff_difference = [=](double x)
  {
    const double first_cutoff = cutoff / std::sqrt(1.0 + curvature * curvature * x * x);
    const double second_cutoff =
        std::sqrt((cutoff * cutoff / (x * x) - 1.0) / (curvature * curvature));
    return first_cutoff - second_cutoff;
  };
  double crossover_lower = lower;
  double crossover_upper = cutoff;
  REQUIRE(cutoff_difference(crossover_lower) < 0.0);
  REQUIRE(cutoff_difference(crossover_upper) > 0.0);
  for (int iteration = 0; iteration < 80; iteration++)
  {
    const double midpoint = 0.5 * (crossover_lower + crossover_upper);
    if (cutoff_difference(midpoint) < 0.0)
    {
      crossover_lower = midpoint;
    }
    else
    {
      crossover_upper = midpoint;
    }
  }
  const double crossover = 0.5 * (crossover_lower + crossover_upper);

  const auto rule = fem::singular::BuildWeightedSegmentQuadrature(96, 0.0, 0.0);
  long double area = 0.0L;
  const std::array<double, 3> breakpoints{lower, crossover, upper};
  for (std::size_t interval = 1; interval < breakpoints.size(); interval++)
  {
    const double x_lower = breakpoints[interval - 1];
    const double x_upper = breakpoints[interval];
    for (const auto &x_quadrature : rule)
    {
      const double x = x_lower + (x_upper - x_lower) * x_quadrature.coordinate;
      const double first_cutoff = cutoff / std::sqrt(1.0 + curvature * curvature * x * x);
      const double second_cutoff =
          x < cutoff
              ? std::sqrt((cutoff * cutoff / (x * x) - 1.0) / (curvature * curvature))
              : 0.0;
      const double y_lower = std::max(first_cutoff, second_cutoff);
      const double y_upper = 1.0 - x;
      REQUIRE(y_upper >= y_lower);
      for (const auto &y_quadrature : rule)
      {
        const double y = y_lower + (y_upper - y_lower) * y_quadrature.coordinate;
        const double surface_jacobian =
            std::sqrt(1.0 + curvature * curvature * (x * x + y * y));
        area += x_quadrature.weight * y_quadrature.weight * (x_upper - x_lower) *
                (y_upper - y_lower) * surface_jacobian;
      }
    }
  }
  return static_cast<double>(area);
}

double BubbleCurvedFaceAreaOutsideThreeCutoffs(double curvature, double cutoff,
                                               int quadrature_order)
{
  const auto height = [=](double x, double y) { return curvature * x * y * (1.0 - x - y); };
  const auto minimum_edge_distance = [&](double x, double y)
  {
    const double z = height(x, y);
    return std::sqrt(std::min(
        {y * y + z * z, x * x + z * z, 0.5 * (1.0 - x - y) * (1.0 - x - y) + z * z}));
  };
  REQUIRE(minimum_edge_distance(1.0 / 3.0, 1.0 / 3.0) > cutoff);

  // This oracle independently integrates direct barycenter-to-boundary rays
  // with a high-order tensor rule and a pointwise physical-distance solve.
  const auto rule =
      fem::singular::BuildWeightedSegmentQuadrature(quadrature_order, 0.0, 0.0);
  constexpr std::array<double, 3> anchor{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
  long double area = 0.0L;
  for (int boundary_edge = 0; boundary_edge < 3; boundary_edge++)
  {
    const int first = (boundary_edge + 1) % 3;
    const int second = (boundary_edge + 2) % 3;
    for (const auto &tangent_quadrature : rule)
    {
      std::array<double, 3> boundary{};
      boundary[first] = 1.0 - tangent_quadrature.coordinate;
      boundary[second] = tangent_quadrature.coordinate;
      const auto distance = [&](double radius)
      {
        std::array<double, 3> lambda{};
        for (int node = 0; node < 3; node++)
        {
          lambda[node] = (1.0 - radius) * anchor[node] + radius * boundary[node];
        }
        return minimum_edge_distance(lambda[1], lambda[2]);
      };
      REQUIRE(distance(0.0) > cutoff);
      REQUIRE(distance(1.0) < cutoff);
      double lower = 0.0;
      double upper = 1.0;
      for (int iteration = 0; iteration < 80; iteration++)
      {
        const double midpoint = 0.5 * (lower + upper);
        if (distance(midpoint) > cutoff)
        {
          lower = midpoint;
        }
        else
        {
          upper = midpoint;
        }
      }
      const double maximum_radius = 0.5 * (lower + upper);
      for (const auto &radial_quadrature : rule)
      {
        const double radius = maximum_radius * radial_quadrature.coordinate;
        std::array<double, 3> lambda{};
        for (int node = 0; node < 3; node++)
        {
          lambda[node] = (1.0 - radius) * anchor[node] + radius * boundary[node];
        }
        const double x = lambda[1];
        const double y = lambda[2];
        const double dz_dx = curvature * y * (1.0 - 2.0 * x - y);
        const double dz_dy = curvature * x * (1.0 - x - 2.0 * y);
        const double surface_jacobian = std::sqrt(1.0 + dz_dx * dz_dx + dz_dy * dz_dy);
        area += tangent_quadrature.weight * radial_quadrature.weight * maximum_radius *
                radius / 3.0 * surface_jacobian;
      }
    }
  }
  return static_cast<double>(area);
}

MaterialOperator MakeSubstrateMaterial(const Mesh &mesh)
{
  config::MaterialData substrate;
  substrate.attributes = {1};
  substrate.epsilon_r.s = {4.0, 4.0, 4.0};
  config::MaterialData vacuum;
  vacuum.attributes = {2};
  vacuum.epsilon_r.s = {1.0, 1.0, 1.0};
  return MaterialOperator({substrate, vacuum}, config::PeriodicBoundaryData(),
                          ProblemType::BOUNDARYMODE, mesh);
}

fem::singular::ParallelDofNumbering MakeNumbering(bool enriched)
{
  fem::singular::ParallelDofNumbering numbering;
  if (enriched)
  {
    numbering.nd.global_local_size = 1;
    numbering.nd.local_size = 1;
    numbering.nd.global_size = 1;
    numbering.nd.owned_size = 1;
    numbering.nd.owner = {0};
    numbering.nd.local_to_true = {0};
  }
  return numbering;
}

fem::singular::ParallelDofNumbering MakeH1Numbering(int enrichment_size)
{
  auto numbering = fem::singular::ParallelDofNumbering{};
  if (enrichment_size > 0)
  {
    numbering.h1.global_local_size = enrichment_size;
    numbering.h1.local_size = enrichment_size;
    numbering.h1.global_size = enrichment_size;
    numbering.h1.owned_size = enrichment_size;
    numbering.h1.owner.assign(enrichment_size, 0);
    numbering.h1.local_to_true.resize(enrichment_size);
    std::iota(numbering.h1.local_to_true.begin(), numbering.h1.local_to_true.end(), 0);
  }
  return numbering;
}

fem::singular::DofKey MakeTestDofKey(fem::singular::HigherOrderBasisFamily family,
                                     const fem::singular::EntityKey &singular_entity)
{
  fem::singular::DofKey key;
  key.family = family;
  key.order = 1;
  key.singular_entity = singular_entity;
  key.support_entity = {2, {0, 1, -1, -1}};
  key.component_entity = {2, {0, 1, -1, -1}};
  return key;
}

}  // namespace

TEST_CASE("Singular surface postprocessor reproduces constant-field MS energy",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = MakeSurfaceTestMesh();
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));
  auto material = MakeSubstrateMaterial(mesh);

  mfem::ND_FECollection collection(1, 2);
  mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
  mfem::VectorFunctionCoefficient coefficient(2,
                                              [](const mfem::Vector &, mfem::Vector &value)
                                              {
                                                value.SetSize(2);
                                                value[0] = 0.0;
                                                value[1] = 2.0;
                                              });
  mfem::ParGridFunction field(&fespace);
  field.ProjectCoefficient(coefficient);
  mfem::Vector true_dofs(fespace.GetTrueVSize());
  field.ParallelProject(true_dofs);
  mfem::Vector zero(true_dofs.Size());
  zero = 0.0;

  fem::singular::TriangleDofTopology topology;
  topology.elements.resize(1);
  auto numbering = MakeNumbering(false);
  fem::singular::TriangleEnrichedNDFieldEvaluator real(topology, numbering, fespace);
  fem::singular::TriangleEnrichedNDFieldEvaluator imaginary(topology, numbering, fespace);
  real.SetFromTrueDofs(true_dofs);
  imaginary.SetFromTrueDofs(zero);

  TriangleSingularSurfacePostOperator postoperator(MakeMSPostprocessing(), material,
                                                   fespace);
  const auto measurements =
      postoperator.Measure(real, imaginary, 1.0, {4, 1.0e-12, 1.0e-11, 8});
  REQUIRE(measurements.size() == 1);
  CHECK(measurements[0].index == 7);
  // Two unit-normal projections contribute: the x-normal edge has zero field,
  // while the y-normal edge has |epsilon E_n|^2 = 64 and unit length.
  const double exact = 0.16 + 0.08 * std::sqrt(2.0);
  CHECK_THAT(measurements[0].energy, WithinAbs(exact, 3.0e-13));
  CHECK_THAT(measurements[0].participation, WithinAbs(exact, 3.0e-13));
  CHECK_THAT(measurements[0].quality_factor, WithinAbs(1.0 / (0.125 * exact), 3.0e-11));
}

TEST_CASE("Singular surface postprocessor uses integrable corner weight",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = MakeSurfaceTestMesh();
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));
  auto material = MakeSubstrateMaterial(mesh);

  mfem::ND_FECollection collection(1, 2);
  mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
  mfem::Vector combined(fespace.GetTrueVSize() + 1);
  combined = 0.0;
  combined[fespace.GetTrueVSize()] = 1.0;
  mfem::Vector zero(combined.Size());
  zero = 0.0;

  const auto measure = [&](double nu, int quadrature_order, double edge_cutoff = 0.0)
  {
    const fem::singular::TriangleBasis basis{
        fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, {0, 1, 2}, 1, nu};
    fem::singular::TriangleDofTopology topology;
    topology.nd_dofs.resize(1);
    topology.nd_dofs[0].family = basis.family;
    topology.nd_dofs[0].order = basis.order;
    topology.elements.resize(1);
    topology.elements[0].nd = {{0, basis}};
    auto numbering = MakeNumbering(true);
    fem::singular::TriangleEnrichedNDFieldEvaluator real(topology, numbering, fespace);
    fem::singular::TriangleEnrichedNDFieldEvaluator imaginary(topology, numbering, fespace);
    real.SetFromTrueDofs(combined);
    imaginary.SetFromTrueDofs(zero);
    TriangleSingularSurfacePostOperator postoperator(MakeMSPostprocessing(edge_cutoff),
                                                     material, fespace);
    return postoperator
        .Measure(real, imaginary, 1.0, {quadrature_order, 1.0e-12, 2.0e-9, 12})
        .front()
        .energy;
  };

  constexpr double transmission_nu = 0.525553491856;
  const double order4 = measure(transmission_nu, 4);
  const double order8 = measure(transmission_nu, 8);
  const double order12 = measure(transmission_nu, 12);
  CHECK(std::isfinite(order12));
  CHECK(order12 > 0.0);
  CHECK_THAT(order4, WithinRel(order12, 5.0e-9));
  CHECK_THAT(order8, WithinRel(order12, 5.0e-9));
  CHECK_THROWS_AS(measure(0.5, 8), std::domain_error);
  const double cutoff_large_order4 = measure(0.5, 4, 1.0e-2);
  const double cutoff_large = measure(0.5, 8, 1.0e-2);
  const double cutoff_small_order4 = measure(0.5, 4, 1.0e-3);
  const double cutoff_small = measure(0.5, 8, 1.0e-3);
  CHECK(std::isfinite(cutoff_large));
  CHECK(cutoff_large > 0.0);
  CHECK(cutoff_small > cutoff_large);
  CHECK_THAT(cutoff_large_order4, WithinRel(cutoff_large, 5.0e-9));
  CHECK_THAT(cutoff_small_order4, WithinRel(cutoff_small, 5.0e-9));
}

TEST_CASE("Singular surface cutoff spans refined boundary segments",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr double cutoff = 0.05;
  const auto measure = [=](const std::vector<double> &breakpoints)
  {
    auto serial_mesh = MakeSplitCutoffSurfaceTestMesh(breakpoints);
    auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
    Mesh mesh(std::move(par_mesh));
    auto material = MakeSubstrateMaterial(mesh);

    mfem::H1_FECollection h1_collection(1, 2);
    mfem::ND_FECollection nd_collection(1, 2);
    mfem::ParFiniteElementSpace h1_fespace(&mesh.Get(), &h1_collection);
    mfem::ParFiniteElementSpace nd_fespace(&mesh.Get(), &nd_collection);
    mfem::FunctionCoefficient potential_coefficient([](const mfem::Vector &point)
                                                    { return 2.0 * point[1]; });
    mfem::VectorFunctionCoefficient field_coefficient(
        2,
        [](const mfem::Vector &, mfem::Vector &value)
        {
          value.SetSize(2);
          value = 0.0;
          value[1] = 2.0;
        });
    mfem::ParGridFunction potential(&h1_fespace);
    mfem::ParGridFunction field(&nd_fespace);
    potential.ProjectCoefficient(potential_coefficient);
    field.ProjectCoefficient(field_coefficient);
    mfem::Vector h1_standard(h1_fespace.GetTrueVSize());
    mfem::Vector nd_standard(nd_fespace.GetTrueVSize());
    potential.ParallelProject(h1_standard);
    field.ParallelProject(nd_standard);

    mfem::Array<int> element_vertices;
    mesh.Get().GetElementVertices(0, element_vertices);
    REQUIRE(element_vertices.Size() == 3);
    const int *tip_location =
        std::find(element_vertices.begin(), element_vertices.end(), 0);
    REQUIRE(tip_location != element_vertices.end());
    const int tip_node = static_cast<int>(tip_location - element_vertices.begin());
    std::array<int, 3> basis_nodes{tip_node, -1, -1};
    int basis_node = 1;
    for (int node = 0; node < 3; node++)
    {
      if (node != tip_node)
      {
        basis_nodes[basis_node++] = node;
      }
    }
    const fem::singular::TriangleBasis basis{
        fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, basis_nodes, 1, 0.5};
    fem::singular::TriangleDofTopology topology;
    topology.h1_dofs.resize(1);
    topology.h1_dofs[0].family = basis.family;
    topology.h1_dofs[0].order = basis.order;
    topology.nd_dofs.resize(1);
    topology.nd_dofs[0].family = basis.family;
    topology.nd_dofs[0].order = basis.order;
    topology.h1_to_nd = {0};
    topology.elements.resize(h1_fespace.GetNE());
    topology.elements[0].h1 = {{0, basis}};
    topology.elements[0].nd = {{0, basis}};

    const auto h1_numbering = MakeH1Numbering(1);
    const auto nd_numbering = MakeNumbering(true);
    mfem::Vector h1_combined(h1_standard.Size() + 1);
    mfem::Vector nd_combined(nd_standard.Size() + 1);
    h1_combined = 0.0;
    nd_combined = 0.0;
    for (int dof = 0; dof < h1_standard.Size(); dof++)
    {
      h1_combined[dof] = h1_standard[dof];
    }
    for (int dof = 0; dof < nd_standard.Size(); dof++)
    {
      nd_combined[dof] = nd_standard[dof];
    }
    mfem::Vector h1_zero(h1_combined.Size());
    mfem::Vector nd_zero(nd_combined.Size());
    h1_zero = 0.0;
    nd_zero = 0.0;
    fem::singular::TriangleEnrichedH1FieldEvaluator h1_real(topology, h1_numbering,
                                                            h1_fespace);
    fem::singular::TriangleEnrichedH1FieldEvaluator h1_imaginary(topology, h1_numbering,
                                                                 h1_fespace);
    fem::singular::TriangleEnrichedNDFieldEvaluator nd_real(topology, nd_numbering,
                                                            nd_fespace);
    fem::singular::TriangleEnrichedNDFieldEvaluator nd_imaginary(topology, nd_numbering,
                                                                 nd_fespace);
    h1_real.SetFromTrueDofs(h1_combined);
    h1_imaginary.SetFromTrueDofs(h1_zero);
    nd_real.SetFromTrueDofs(nd_combined);
    nd_imaginary.SetFromTrueDofs(nd_zero);

    TriangleSingularSurfacePostOperator h1_postoperator(MakeMSPostprocessing(cutoff),
                                                        material, h1_fespace);
    TriangleSingularSurfacePostOperator nd_postoperator(MakeMSPostprocessing(cutoff),
                                                        material, nd_fespace);
    const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-12, 2.0e-9, 12};
    const double h1_energy =
        h1_postoperator.MeasureElectrostatic(h1_real, h1_imaginary, 1.0, options)
            .front()
            .energy;
    const double nd_energy =
        nd_postoperator.Measure(nd_real, nd_imaginary, 1.0, options).front().energy;
    return std::pair{h1_energy, nd_energy};
  };

  const auto coarse = measure({0.0, 1.0});
  const auto refined = measure({0.0, 0.02, 0.04, 0.10, 1.0});
  constexpr double exact = 0.16 * (1.0 - cutoff);
  CHECK_THAT(coarse.first, WithinAbs(exact, 3.0e-11));
  CHECK_THAT(coarse.second, WithinAbs(exact, 3.0e-11));
  CHECK_THAT(refined.first, WithinAbs(exact, 3.0e-11));
  CHECK_THAT(refined.second, WithinAbs(exact, 3.0e-11));
}

TEST_CASE("Singular electrostatic surface postprocessor uses combined H1 gradient",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = MakeSurfaceTestMesh();
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));
  auto material = MakeSubstrateMaterial(mesh);

  mfem::H1_FECollection collection(1, 2);
  mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
  mfem::FunctionCoefficient coefficient([](const mfem::Vector &point)
                                        { return 2.0 * point[1]; });
  mfem::ParGridFunction potential(&fespace);
  potential.ProjectCoefficient(coefficient);
  mfem::Vector true_dofs(fespace.GetTrueVSize());
  potential.ParallelProject(true_dofs);
  mfem::Vector zero(true_dofs.Size());
  zero = 0.0;

  fem::singular::TriangleDofTopology topology;
  topology.elements.resize(1);
  auto numbering = MakeH1Numbering(false);
  fem::singular::TriangleEnrichedH1FieldEvaluator real(topology, numbering, fespace);
  fem::singular::TriangleEnrichedH1FieldEvaluator imaginary(topology, numbering, fespace);
  real.SetFromTrueDofs(true_dofs);
  imaginary.SetFromTrueDofs(zero);

  TriangleSingularSurfacePostOperator postoperator(MakeMSPostprocessing(), material,
                                                   fespace);
  const auto measurements =
      postoperator.MeasureElectrostatic(real, imaginary, 1.0, {4, 1.0e-12, 1.0e-11, 8});
  REQUIRE(measurements.size() == 1);
  const double exact = 0.16 + 0.08 * std::sqrt(2.0);
  CHECK_THAT(measurements[0].energy, WithinAbs(exact, 3.0e-13));
  CHECK_THAT(measurements[0].participation, WithinAbs(exact, 3.0e-13));
}

TEST_CASE("Singular surface postprocessor handles nonconforming interface segments",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const auto measure = [](bool refine)
  {
    auto serial_mesh = MakePartitionedSurfaceTestMesh();
    if (refine)
    {
      serial_mesh->EnsureNCMesh(true);
    }
    auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
    if (refine)
    {
      const int initial_boundary_elements = par_mesh->GetNBE();
      mfem::Array<int> marked(1);
      marked[0] = 0;
      par_mesh->GeneralRefinement(marked, -1, 1);
      REQUIRE(par_mesh->GetNBE() > initial_boundary_elements);
    }
    Mesh mesh(std::move(par_mesh));
    auto material = MakeSubstrateMaterial(mesh);

    mfem::H1_FECollection collection(1, 2);
    mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
    mfem::FunctionCoefficient coefficient([](const mfem::Vector &point)
                                          { return 2.0 * point[1]; });
    mfem::ParGridFunction potential(&fespace);
    potential.ProjectCoefficient(coefficient);
    mfem::Vector true_dofs(fespace.GetTrueVSize());
    potential.ParallelProject(true_dofs);
    mfem::Vector zero(true_dofs.Size());
    zero = 0.0;

    fem::singular::TriangleDofTopology topology;
    topology.elements.resize(fespace.GetNE());
    const auto numbering = MakeH1Numbering(0);
    fem::singular::TriangleEnrichedH1FieldEvaluator real(topology, numbering, fespace);
    fem::singular::TriangleEnrichedH1FieldEvaluator imaginary(topology, numbering, fespace);
    real.SetFromTrueDofs(true_dofs);
    imaginary.SetFromTrueDofs(zero);

    TriangleSingularSurfacePostOperator postoperator(MakeMSPostprocessing(), material,
                                                     fespace);
    return postoperator.MeasureElectrostatic(real, imaginary, 1.0, {4, 1.0e-12, 1.0e-11, 8})
        .front()
        .energy;
  };

  const double conforming = measure(false);
  const double nonconforming = measure(true);
  CHECK_THAT(nonconforming, WithinAbs(conforming, 3.0e-13));
}

TEST_CASE("Tetrahedral singular surface postprocessor handles nonconforming interface "
          "triangles",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const auto measure = [](bool refine)
  {
    auto serial_mesh = MakePartitionedTetrahedronSurfaceTestMesh();
    if (refine)
    {
      serial_mesh->EnsureNCMesh(true);
    }
    auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
    if (refine)
    {
      const int initial_boundary_elements = par_mesh->GetNBE();
      mfem::Array<int> marked(1);
      marked[0] = 0;
      par_mesh->GeneralRefinement(marked, -1, 1);
      REQUIRE(par_mesh->GetNBE() > initial_boundary_elements);
    }
    Mesh mesh(std::move(par_mesh));
    auto material = MakeSubstrateMaterial(mesh);

    mfem::H1_FECollection collection(1, 3);
    mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
    mfem::FunctionCoefficient coefficient([](const mfem::Vector &point)
                                          { return 2.0 * point[2]; });
    mfem::ParGridFunction potential(&fespace);
    potential.ProjectCoefficient(coefficient);
    mfem::Vector true_dofs(fespace.GetTrueVSize());
    potential.ParallelProject(true_dofs);
    mfem::Vector zero(true_dofs.Size());
    zero = 0.0;

    fem::singular::DofTopology topology;
    topology.elements.resize(fespace.GetNE());
    const auto numbering = MakeH1Numbering(0);
    fem::singular::EnrichedH1FieldEvaluator real(topology, numbering, fespace);
    fem::singular::EnrichedH1FieldEvaluator imaginary(topology, numbering, fespace);
    real.SetFromTrueDofs(true_dofs);
    imaginary.SetFromTrueDofs(zero);

    TetrahedronSingularSurfacePostOperator postoperator(MakeMSPostprocessing(), material,
                                                        fespace);
    return postoperator.MeasureElectrostatic(real, imaginary, 1.0, {4, 1.0e-12, 1.0e-11, 8})
        .front()
        .energy;
  };

  const double conforming = measure(false);
  const double nonconforming = measure(true);
  CHECK_THAT(nonconforming, WithinAbs(conforming, 3.0e-13));
}

TEST_CASE("Tetrahedral singular electrostatic surface postprocessor reproduces a "
          "constant field",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = MakeTetrahedronSurfaceTestMesh();
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));
  auto material = MakeSubstrateMaterial(mesh);

  mfem::H1_FECollection collection(1, 3);
  mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
  mfem::FunctionCoefficient coefficient([](const mfem::Vector &point)
                                        { return 2.0 * point[2]; });
  mfem::ParGridFunction potential(&fespace);
  potential.ProjectCoefficient(coefficient);
  mfem::Vector true_dofs(fespace.GetTrueVSize());
  potential.ParallelProject(true_dofs);
  mfem::Vector zero(true_dofs.Size());
  zero = 0.0;

  fem::singular::DofTopology topology;
  topology.elements.resize(1);
  const auto numbering = MakeH1Numbering(false);
  fem::singular::EnrichedH1FieldEvaluator real(topology, numbering, fespace);
  fem::singular::EnrichedH1FieldEvaluator imaginary(topology, numbering, fespace);
  real.SetFromTrueDofs(true_dofs);
  imaginary.SetFromTrueDofs(zero);

  TetrahedronSingularSurfacePostOperator postoperator(MakeMSPostprocessing(), material,
                                                      fespace);
  const auto measurements =
      postoperator.MeasureElectrostatic(real, imaginary, 1.0, {4, 1.0e-12, 1.0e-11, 8});
  REQUIRE(measurements.size() == 1);
  CHECK_THAT(measurements[0].energy, WithinAbs(0.08, 3.0e-13));
  CHECK_THAT(measurements[0].participation, WithinAbs(0.08, 3.0e-13));
  CHECK_THAT(measurements[0].quality_factor, WithinAbs(100.0, 3.0e-11));
}

TEST_CASE("Tetrahedral singular full-wave surface postprocessor reproduces a "
          "constant field",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = MakeTetrahedronSurfaceTestMesh();
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));
  auto material = MakeSubstrateMaterial(mesh);

  mfem::ND_FECollection collection(1, 3);
  mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
  mfem::VectorFunctionCoefficient coefficient(3,
                                              [](const mfem::Vector &, mfem::Vector &value)
                                              {
                                                value.SetSize(3);
                                                value = 0.0;
                                                value[2] = 2.0;
                                              });
  mfem::ParGridFunction field(&fespace);
  field.ProjectCoefficient(coefficient);
  mfem::Vector true_dofs(fespace.GetTrueVSize());
  field.ParallelProject(true_dofs);
  mfem::Vector zero(true_dofs.Size());
  zero = 0.0;

  fem::singular::DofTopology topology;
  topology.elements.resize(1);
  const auto numbering = MakeNumbering(false);
  fem::singular::EnrichedNDFieldEvaluator real(topology, numbering, fespace);
  fem::singular::EnrichedNDFieldEvaluator imaginary(topology, numbering, fespace);
  real.SetFromTrueDofs(true_dofs);
  imaginary.SetFromTrueDofs(zero);

  TetrahedronSingularSurfacePostOperator postoperator(MakeMSPostprocessing(), material,
                                                      fespace);
  const auto measurements =
      postoperator.Measure(real, imaginary, 1.0, {4, 1.0e-12, 1.0e-11, 8});
  REQUIRE(measurements.size() == 1);
  CHECK_THAT(measurements[0].energy, WithinAbs(0.08, 3.0e-13));
  CHECK_THAT(measurements[0].participation, WithinAbs(0.08, 3.0e-13));
  CHECK_THAT(measurements[0].quality_factor, WithinAbs(100.0, 3.0e-11));
}

TEST_CASE("Tetrahedral singular surface postprocessor fuses compatible interfaces",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto serial_mesh = MakePartitionedTetrahedronSurfaceTestMesh();
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));
  auto material = MakeSubstrateMaterial(mesh);
  const fem::singular::AdaptiveAssemblyOptions options{4, 1.0e-12, 1.0e-11, 8};

  const auto check_measurements = [](const auto &grouped, const auto &ms, const auto &ma)
  {
    REQUIRE(grouped.size() == 2);
    REQUIRE(ms.size() == 1);
    REQUIRE(ma.size() == 1);
    REQUIRE(grouped[0].index == 7);
    REQUIRE(grouped[1].index == 8);
    CHECK_THAT(grouped[0].energy, WithinAbs(ms[0].energy, 3.0e-13));
    CHECK_THAT(grouped[1].energy, WithinAbs(ma[0].energy, 3.0e-13));
  };

  SECTION("Electrostatic H1 gradient")
  {
    mfem::H1_FECollection collection(1, 3);
    mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
    mfem::FunctionCoefficient coefficient([](const mfem::Vector &point)
                                          { return 2.0 * point[2]; });
    mfem::ParGridFunction potential(&fespace);
    potential.ProjectCoefficient(coefficient);
    mfem::Vector true_dofs(fespace.GetTrueVSize());
    potential.ParallelProject(true_dofs);
    mfem::Vector zero(true_dofs.Size());
    zero = 0.0;

    fem::singular::DofTopology topology;
    topology.elements.resize(fespace.GetNE());
    const auto numbering = MakeH1Numbering(0);
    fem::singular::EnrichedH1FieldEvaluator real(topology, numbering, fespace);
    fem::singular::EnrichedH1FieldEvaluator imaginary(topology, numbering, fespace);
    real.SetFromTrueDofs(true_dofs);
    imaginary.SetFromTrueDofs(zero);

    TetrahedronSingularSurfacePostOperator grouped_postoperator(
        MakeMSMAPostprocessing(true, true), material, fespace);
    TetrahedronSingularSurfacePostOperator ms_postoperator(
        MakeMSMAPostprocessing(true, false), material, fespace);
    TetrahedronSingularSurfacePostOperator ma_postoperator(
        MakeMSMAPostprocessing(false, true), material, fespace);
    check_measurements(
        grouped_postoperator.MeasureElectrostatic(real, imaginary, 1.0, options),
        ms_postoperator.MeasureElectrostatic(real, imaginary, 1.0, options),
        ma_postoperator.MeasureElectrostatic(real, imaginary, 1.0, options));
  }

  SECTION("Full-wave Nedelec field")
  {
    mfem::ND_FECollection collection(1, 3);
    mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
    mfem::VectorFunctionCoefficient coefficient(
        3,
        [](const mfem::Vector &, mfem::Vector &value)
        {
          value.SetSize(3);
          value = 0.0;
          value[2] = 2.0;
        });
    mfem::ParGridFunction field(&fespace);
    field.ProjectCoefficient(coefficient);
    mfem::Vector true_dofs(fespace.GetTrueVSize());
    field.ParallelProject(true_dofs);
    mfem::Vector zero(true_dofs.Size());
    zero = 0.0;

    fem::singular::DofTopology topology;
    topology.elements.resize(fespace.GetNE());
    const auto numbering = MakeNumbering(false);
    fem::singular::EnrichedNDFieldEvaluator real(topology, numbering, fespace);
    fem::singular::EnrichedNDFieldEvaluator imaginary(topology, numbering, fespace);
    real.SetFromTrueDofs(true_dofs);
    imaginary.SetFromTrueDofs(zero);

    TetrahedronSingularSurfacePostOperator grouped_postoperator(
        MakeMSMAPostprocessing(true, true), material, fespace);
    TetrahedronSingularSurfacePostOperator ms_postoperator(
        MakeMSMAPostprocessing(true, false), material, fespace);
    TetrahedronSingularSurfacePostOperator ma_postoperator(
        MakeMSMAPostprocessing(false, true), material, fespace);
    const auto grouped = grouped_postoperator.Measure(real, imaginary, 1.0, options);
    check_measurements(grouped, ms_postoperator.Measure(real, imaginary, 1.0, options),
                       ma_postoperator.Measure(real, imaginary, 1.0, options));
    const auto batched = grouped_postoperator.Measure(
        {{&real, &imaginary}, {&imaginary, &real}}, {1.0, 1.0}, options);
    REQUIRE(batched.size() == 2);
    REQUIRE(batched[0].size() == grouped.size());
    REQUIRE(batched[1].size() == grouped.size());
    for (std::size_t interface = 0; interface < grouped.size(); interface++)
    {
      CHECK_THAT(batched[0][interface].energy,
                 WithinAbs(grouped[interface].energy, 3.0e-13));
      CHECK_THAT(batched[1][interface].energy,
                 WithinAbs(grouped[interface].energy, 3.0e-13));
    }
  }
}

TEST_CASE("Tetrahedral singular surface postprocessor evaluates off-rank traces",
          "[singularsurface][Parallel]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 2);
  auto serial_mesh = MakePartitionedTetrahedronSurfaceTestMesh();
  const std::array<int, 2> partition{0, 1};
  auto par_mesh =
      std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh, partition.data());
  par_mesh->ExchangeFaceNbrData();
  Mesh mesh(std::move(par_mesh));
  auto material = MakeSubstrateMaterial(mesh);

  int local_interface_boundaries = 0;
  for (int boundary = 0; boundary < mesh.Get().GetNBE(); boundary++)
  {
    local_interface_boundaries += mesh.Get().GetBdrAttribute(boundary) == 1 ? 1 : 0;
  }
  std::array<int, 2> interface_boundaries{};
  Mpi::Allgather(1, &local_interface_boundaries, interface_boundaries.data(), Mpi::World());
  REQUIRE(interface_boundaries == std::array<int, 2>{1, 0});
  REQUIRE(mesh.Get().GetNE() == 1);
  REQUIRE(mesh.Get().GetAttribute(0) == (Mpi::Rank(Mpi::World()) == 0 ? 2 : 1));

  SECTION("Electrostatic H1 gradient")
  {
    mfem::H1_FECollection collection(1, 3);
    mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
    mfem::FunctionCoefficient coefficient([](const mfem::Vector &point)
                                          { return 2.0 * point[2]; });
    mfem::ParGridFunction potential(&fespace);
    potential.ProjectCoefficient(coefficient);
    mfem::Vector true_dofs(fespace.GetTrueVSize());
    potential.ParallelProject(true_dofs);
    mfem::Vector zero(true_dofs.Size());
    zero = 0.0;

    fem::singular::DofTopology topology;
    topology.elements.resize(fespace.GetNE());
    const auto numbering = MakeH1Numbering(0);
    fem::singular::EnrichedH1FieldEvaluator real(topology, numbering, fespace);
    fem::singular::EnrichedH1FieldEvaluator imaginary(topology, numbering, fespace);
    real.SetFromTrueDofs(true_dofs);
    imaginary.SetFromTrueDofs(zero);

    TetrahedronSingularSurfacePostOperator postoperator(MakeMSPostprocessing(), material,
                                                        fespace);
    const auto measurements =
        postoperator.MeasureElectrostatic(real, imaginary, 1.0, {4, 1.0e-12, 1.0e-11, 8});
    REQUIRE(measurements.size() == 1);
    CHECK_THAT(measurements[0].energy, WithinAbs(0.08, 3.0e-13));
  }

  SECTION("Full-wave Nedelec field")
  {
    mfem::ND_FECollection collection(1, 3);
    mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
    mfem::VectorFunctionCoefficient coefficient(
        3,
        [](const mfem::Vector &, mfem::Vector &value)
        {
          value.SetSize(3);
          value = 0.0;
          value[2] = 2.0;
        });
    mfem::ParGridFunction field(&fespace);
    field.ProjectCoefficient(coefficient);
    mfem::Vector true_dofs(fespace.GetTrueVSize());
    field.ParallelProject(true_dofs);
    mfem::Vector zero(true_dofs.Size());
    zero = 0.0;

    fem::singular::DofTopology topology;
    topology.elements.resize(fespace.GetNE());
    const auto numbering = MakeNumbering(false);
    fem::singular::EnrichedNDFieldEvaluator real(topology, numbering, fespace);
    fem::singular::EnrichedNDFieldEvaluator imaginary(topology, numbering, fespace);
    real.SetFromTrueDofs(true_dofs);
    imaginary.SetFromTrueDofs(zero);

    TetrahedronSingularSurfacePostOperator postoperator(MakeMSPostprocessing(), material,
                                                        fespace);
    const auto measurements =
        postoperator.Measure(real, imaginary, 1.0, {4, 1.0e-12, 1.0e-11, 8});
    REQUIRE(measurements.size() == 1);
    CHECK_THAT(measurements[0].energy, WithinAbs(0.08, 3.0e-13));
  }

  SECTION("Singular H1 and Nedelec enrichment")
  {
    constexpr double nu = 2.0 / 3.0;
    const auto bases =
        fem::singular::EnumerateHigherOrderEdgeGradientBases({0, 1, 2, 3}, 1, nu);
    REQUIRE_FALSE(bases.empty());
    const auto basis = bases.front();
    const auto key =
        MakeTestDofKey(basis.family, fem::singular::EntityKey{2, {0, 1, -1, -1}});

    fem::singular::DofTopology topology;
    topology.elements.resize(1);
    if (Mpi::Rank(Mpi::World()) == 1)
    {
      topology.h1_dofs = {key};
      topology.nd_dofs = {key};
      topology.h1_to_nd = {0};
      topology.elements[0].h1 = {{0, basis}};
      topology.elements[0].nd = {{0, basis}};
    }
    const auto numbering = fem::singular::BuildParallelDofNumbering(Mpi::World(), topology);
    REQUIRE(numbering.h1.owned_size == (Mpi::Rank(Mpi::World()) == 1 ? 1 : 0));
    REQUIRE(numbering.nd.owned_size == (Mpi::Rank(Mpi::World()) == 1 ? 1 : 0));

    mfem::H1_FECollection h1_collection(1, 3);
    mfem::ND_FECollection nd_collection(1, 3);
    mfem::ParFiniteElementSpace h1_fespace(&mesh.Get(), &h1_collection);
    mfem::ParFiniteElementSpace nd_fespace(&mesh.Get(), &nd_collection);
    mfem::Vector h1_combined(h1_fespace.GetTrueVSize() +
                             static_cast<int>(numbering.h1.owned_size));
    mfem::Vector nd_combined(nd_fespace.GetTrueVSize() +
                             static_cast<int>(numbering.nd.owned_size));
    h1_combined = 0.0;
    nd_combined = 0.0;
    if (Mpi::Rank(Mpi::World()) == 1)
    {
      h1_combined[h1_fespace.GetTrueVSize()] = 1.0;
      nd_combined[nd_fespace.GetTrueVSize()] = 1.0;
    }
    mfem::Vector h1_zero(h1_combined.Size());
    mfem::Vector nd_zero(nd_combined.Size());
    h1_zero = 0.0;
    nd_zero = 0.0;

    fem::singular::EnrichedH1FieldEvaluator h1_real(topology, numbering, h1_fespace);
    fem::singular::EnrichedH1FieldEvaluator h1_imaginary(topology, numbering, h1_fespace);
    fem::singular::EnrichedNDFieldEvaluator nd_real(topology, numbering, nd_fespace);
    fem::singular::EnrichedNDFieldEvaluator nd_imaginary(topology, numbering, nd_fespace);
    h1_real.SetFromTrueDofs(h1_combined);
    h1_imaginary.SetFromTrueDofs(h1_zero);
    nd_real.SetFromTrueDofs(nd_combined);
    nd_imaginary.SetFromTrueDofs(nd_zero);

    TetrahedronSingularSurfacePostOperator h1_postoperator(MakeMSPostprocessing(), material,
                                                           h1_fespace);
    TetrahedronSingularSurfacePostOperator nd_postoperator(MakeMSPostprocessing(), material,
                                                           nd_fespace);
    const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-12, 2.0e-9, 12};
    const double parallel_h1 =
        h1_postoperator.MeasureElectrostatic(h1_real, h1_imaginary, 1.0, options)
            .front()
            .energy;
    const double parallel_nd =
        nd_postoperator.Measure(nd_real, nd_imaginary, 1.0, options).front().energy;

    auto reference_serial_mesh = MakePartitionedTetrahedronSurfaceTestMesh();
    auto reference_par_mesh =
        std::make_unique<mfem::ParMesh>(MPI_COMM_SELF, *reference_serial_mesh);
    Mesh reference_mesh(std::move(reference_par_mesh));
    auto reference_material = MakeSubstrateMaterial(reference_mesh);
    mfem::H1_FECollection reference_collection(1, 3);
    mfem::ParFiniteElementSpace reference_fespace(&reference_mesh.Get(),
                                                  &reference_collection);
    fem::singular::DofTopology reference_topology;
    reference_topology.h1_dofs = {key};
    reference_topology.elements.resize(2);
    reference_topology.elements[1].h1 = {{0, basis}};
    const auto reference_numbering = MakeH1Numbering(1);
    mfem::Vector reference_combined(reference_fespace.GetTrueVSize() + 1);
    reference_combined = 0.0;
    reference_combined[reference_fespace.GetTrueVSize()] = 1.0;
    mfem::Vector reference_zero(reference_combined.Size());
    reference_zero = 0.0;
    fem::singular::EnrichedH1FieldEvaluator reference_real(
        reference_topology, reference_numbering, reference_fespace);
    fem::singular::EnrichedH1FieldEvaluator reference_imaginary(
        reference_topology, reference_numbering, reference_fespace);
    reference_real.SetFromTrueDofs(reference_combined);
    reference_imaginary.SetFromTrueDofs(reference_zero);
    TetrahedronSingularSurfacePostOperator reference_postoperator(
        MakeMSPostprocessing(), reference_material, reference_fespace);
    const double reference =
        reference_postoperator
            .MeasureElectrostatic(reference_real, reference_imaginary, 1.0, options)
            .front()
            .energy;

    REQUIRE(std::isfinite(reference));
    REQUIRE(reference > 0.0);
    CHECK_THAT(parallel_h1, WithinRel(reference, 2.0e-12));
    CHECK_THAT(parallel_nd, WithinRel(reference, 2.0e-12));
  }
}

TEST_CASE("Triangular singular surface postprocessor evaluates off-rank traces",
          "[singularsurface][Parallel]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 2);
  auto serial_mesh = MakePartitionedSurfaceTestMesh();
  const std::array<int, 2> partition{0, 1};
  auto par_mesh =
      std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh, partition.data());
  par_mesh->ExchangeFaceNbrData();
  Mesh mesh(std::move(par_mesh));
  auto material = MakeSubstrateMaterial(mesh);

  int local_interface_boundaries = 0;
  for (int boundary = 0; boundary < mesh.Get().GetNBE(); boundary++)
  {
    local_interface_boundaries += mesh.Get().GetBdrAttribute(boundary) == 1 ? 1 : 0;
  }
  std::array<int, 2> interface_boundaries{};
  Mpi::Allgather(1, &local_interface_boundaries, interface_boundaries.data(), Mpi::World());
  REQUIRE(interface_boundaries == std::array<int, 2>{1, 0});
  REQUIRE(mesh.Get().GetNE() == 1);
  REQUIRE(mesh.Get().GetAttribute(0) == (Mpi::Rank(Mpi::World()) == 0 ? 2 : 1));

  SECTION("Electrostatic H1 gradient")
  {
    mfem::H1_FECollection collection(1, 2);
    mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
    mfem::FunctionCoefficient coefficient([](const mfem::Vector &point)
                                          { return 2.0 * point[1]; });
    mfem::ParGridFunction potential(&fespace);
    potential.ProjectCoefficient(coefficient);
    mfem::Vector true_dofs(fespace.GetTrueVSize());
    potential.ParallelProject(true_dofs);
    mfem::Vector zero(true_dofs.Size());
    zero = 0.0;

    fem::singular::TriangleDofTopology topology;
    topology.elements.resize(fespace.GetNE());
    const auto numbering = MakeH1Numbering(0);
    fem::singular::TriangleEnrichedH1FieldEvaluator real(topology, numbering, fespace);
    fem::singular::TriangleEnrichedH1FieldEvaluator imaginary(topology, numbering, fespace);
    real.SetFromTrueDofs(true_dofs);
    imaginary.SetFromTrueDofs(zero);

    TriangleSingularSurfacePostOperator postoperator(MakeMSPostprocessing(), material,
                                                     fespace);
    const auto measurements =
        postoperator.MeasureElectrostatic(real, imaginary, 1.0, {4, 1.0e-12, 1.0e-11, 8});
    REQUIRE(measurements.size() == 1);
    CHECK_THAT(measurements[0].energy, WithinAbs(0.16, 3.0e-13));
  }

  SECTION("Full-wave Nedelec field")
  {
    mfem::ND_FECollection collection(1, 2);
    mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
    mfem::VectorFunctionCoefficient coefficient(
        2,
        [](const mfem::Vector &, mfem::Vector &value)
        {
          value.SetSize(2);
          value = 0.0;
          value[1] = 2.0;
        });
    mfem::ParGridFunction field(&fespace);
    field.ProjectCoefficient(coefficient);
    mfem::Vector true_dofs(fespace.GetTrueVSize());
    field.ParallelProject(true_dofs);
    mfem::Vector zero(true_dofs.Size());
    zero = 0.0;

    fem::singular::TriangleDofTopology topology;
    topology.elements.resize(fespace.GetNE());
    const auto numbering = MakeNumbering(false);
    fem::singular::TriangleEnrichedNDFieldEvaluator real(topology, numbering, fespace);
    fem::singular::TriangleEnrichedNDFieldEvaluator imaginary(topology, numbering, fespace);
    real.SetFromTrueDofs(true_dofs);
    imaginary.SetFromTrueDofs(zero);

    TriangleSingularSurfacePostOperator postoperator(MakeMSPostprocessing(), material,
                                                     fespace);
    const auto measurements =
        postoperator.Measure(real, imaginary, 1.0, {4, 1.0e-12, 1.0e-11, 8});
    REQUIRE(measurements.size() == 1);
    CHECK_THAT(measurements[0].energy, WithinAbs(0.16, 3.0e-13));
  }

  SECTION("BoundaryMode longitudinal H1 field")
  {
    mfem::ND_FECollection nd_collection(1, 2);
    mfem::H1_FECollection h1_collection(1, 2);
    mfem::ParFiniteElementSpace nd_fespace(&mesh.Get(), &nd_collection);
    mfem::ParFiniteElementSpace h1_fespace(&mesh.Get(), &h1_collection);

    mfem::Vector nd_zero(nd_fespace.GetTrueVSize());
    nd_zero = 0.0;
    mfem::FunctionCoefficient coefficient([](const mfem::Vector &) { return 2.0; });
    mfem::ParGridFunction longitudinal_field(&h1_fespace);
    longitudinal_field.ProjectCoefficient(coefficient);
    mfem::Vector h1_true_dofs(h1_fespace.GetTrueVSize());
    longitudinal_field.ParallelProject(h1_true_dofs);
    mfem::Vector h1_zero(h1_true_dofs.Size());
    h1_zero = 0.0;

    fem::singular::TriangleDofTopology topology;
    topology.elements.resize(nd_fespace.GetNE());
    const auto numbering = MakeH1Numbering(0);
    fem::singular::TriangleEnrichedNDFieldEvaluator nd_real(topology, numbering,
                                                            nd_fespace);
    fem::singular::TriangleEnrichedNDFieldEvaluator nd_imaginary(topology, numbering,
                                                                 nd_fespace);
    fem::singular::TriangleEnrichedH1FieldEvaluator h1_real(topology, numbering,
                                                            h1_fespace);
    fem::singular::TriangleEnrichedH1FieldEvaluator h1_imaginary(topology, numbering,
                                                                 h1_fespace);
    nd_real.SetFromTrueDofs(nd_zero);
    nd_imaginary.SetFromTrueDofs(nd_zero);
    h1_real.SetFromTrueDofs(h1_true_dofs);
    h1_imaginary.SetFromTrueDofs(h1_zero);

    TriangleSingularSurfacePostOperator postoperator(MakeDefaultPostprocessing(), material,
                                                     nd_fespace);
    const auto measurements = postoperator.Measure(
        nd_real, nd_imaginary, 1.0, {4, 1.0e-12, 1.0e-11, 8}, &h1_real, &h1_imaginary);
    REQUIRE(measurements.size() == 1);
    CHECK_THAT(measurements[0].energy, WithinAbs(0.04, 3.0e-13));
  }

  SECTION("Singular H1 and Nedelec enrichment")
  {
    constexpr double nu = 2.0 / 3.0;
    const fem::singular::TriangleBasis basis{
        fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, {0, 2, 1}, 1, nu};
    const auto key =
        MakeTestDofKey(basis.family, fem::singular::EntityKey{1, {0, -1, -1, -1}});

    fem::singular::TriangleDofTopology topology;
    topology.elements.resize(1);
    if (Mpi::Rank(Mpi::World()) == 1)
    {
      topology.h1_dofs = {key};
      topology.nd_dofs = {key};
      topology.h1_to_nd = {0};
      topology.elements[0].h1 = {{0, basis}};
      topology.elements[0].nd = {{0, basis}};
    }
    const auto numbering = fem::singular::BuildParallelDofNumbering(Mpi::World(), topology);
    REQUIRE(numbering.h1.owned_size == (Mpi::Rank(Mpi::World()) == 1 ? 1 : 0));
    REQUIRE(numbering.nd.owned_size == (Mpi::Rank(Mpi::World()) == 1 ? 1 : 0));

    mfem::H1_FECollection h1_collection(1, 2);
    mfem::ND_FECollection nd_collection(1, 2);
    mfem::ParFiniteElementSpace h1_fespace(&mesh.Get(), &h1_collection);
    mfem::ParFiniteElementSpace nd_fespace(&mesh.Get(), &nd_collection);
    mfem::Vector h1_combined(h1_fespace.GetTrueVSize() +
                             static_cast<int>(numbering.h1.owned_size));
    mfem::Vector nd_combined(nd_fespace.GetTrueVSize() +
                             static_cast<int>(numbering.nd.owned_size));
    h1_combined = 0.0;
    nd_combined = 0.0;
    if (Mpi::Rank(Mpi::World()) == 1)
    {
      h1_combined[h1_fespace.GetTrueVSize()] = 1.0;
      nd_combined[nd_fespace.GetTrueVSize()] = 1.0;
    }
    mfem::Vector h1_zero(h1_combined.Size());
    mfem::Vector nd_zero(nd_combined.Size());
    h1_zero = 0.0;
    nd_zero = 0.0;

    fem::singular::TriangleEnrichedH1FieldEvaluator h1_real(topology, numbering,
                                                            h1_fespace);
    fem::singular::TriangleEnrichedH1FieldEvaluator h1_imaginary(topology, numbering,
                                                                 h1_fespace);
    fem::singular::TriangleEnrichedNDFieldEvaluator nd_real(topology, numbering,
                                                            nd_fespace);
    fem::singular::TriangleEnrichedNDFieldEvaluator nd_imaginary(topology, numbering,
                                                                 nd_fespace);
    fem::singular::TriangleEnrichedNDFieldEvaluator nd_zero_real(topology, numbering,
                                                                 nd_fespace);
    fem::singular::TriangleEnrichedNDFieldEvaluator nd_zero_imaginary(topology, numbering,
                                                                      nd_fespace);
    h1_real.SetFromTrueDofs(h1_combined);
    h1_imaginary.SetFromTrueDofs(h1_zero);
    nd_real.SetFromTrueDofs(nd_combined);
    nd_imaginary.SetFromTrueDofs(nd_zero);
    nd_zero_real.SetFromTrueDofs(nd_zero);
    nd_zero_imaginary.SetFromTrueDofs(nd_zero);

    TriangleSingularSurfacePostOperator h1_postoperator(MakeMSPostprocessing(), material,
                                                        h1_fespace);
    TriangleSingularSurfacePostOperator nd_postoperator(MakeMSPostprocessing(), material,
                                                        nd_fespace);
    TriangleSingularSurfacePostOperator longitudinal_postoperator(
        MakeDefaultPostprocessing(), material, nd_fespace);
    const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-12, 2.0e-9, 12};
    const double parallel_h1 =
        h1_postoperator.MeasureElectrostatic(h1_real, h1_imaginary, 1.0, options)
            .front()
            .energy;
    const double parallel_nd =
        nd_postoperator.Measure(nd_real, nd_imaginary, 1.0, options).front().energy;
    const double parallel_longitudinal =
        longitudinal_postoperator
            .Measure(nd_zero_real, nd_zero_imaginary, 1.0, options, &h1_real, &h1_imaginary)
            .front()
            .energy;

    auto reference_serial_mesh = MakePartitionedSurfaceTestMesh();
    auto reference_par_mesh =
        std::make_unique<mfem::ParMesh>(MPI_COMM_SELF, *reference_serial_mesh);
    Mesh reference_mesh(std::move(reference_par_mesh));
    auto reference_material = MakeSubstrateMaterial(reference_mesh);
    mfem::H1_FECollection reference_collection(1, 2);
    mfem::ND_FECollection reference_nd_collection(1, 2);
    mfem::ParFiniteElementSpace reference_fespace(&reference_mesh.Get(),
                                                  &reference_collection);
    mfem::ParFiniteElementSpace reference_nd_fespace(&reference_mesh.Get(),
                                                     &reference_nd_collection);
    fem::singular::TriangleDofTopology reference_topology;
    reference_topology.h1_dofs = {key};
    reference_topology.elements.resize(2);
    reference_topology.elements[1].h1 = {{0, basis}};
    const auto reference_numbering = MakeH1Numbering(1);
    mfem::Vector reference_combined(reference_fespace.GetTrueVSize() + 1);
    reference_combined = 0.0;
    reference_combined[reference_fespace.GetTrueVSize()] = 1.0;
    mfem::Vector reference_zero(reference_combined.Size());
    reference_zero = 0.0;
    fem::singular::TriangleEnrichedH1FieldEvaluator reference_real(
        reference_topology, reference_numbering, reference_fespace);
    fem::singular::TriangleEnrichedH1FieldEvaluator reference_imaginary(
        reference_topology, reference_numbering, reference_fespace);
    reference_real.SetFromTrueDofs(reference_combined);
    reference_imaginary.SetFromTrueDofs(reference_zero);
    mfem::Vector reference_nd_zero(reference_nd_fespace.GetTrueVSize());
    reference_nd_zero = 0.0;
    fem::singular::TriangleEnrichedNDFieldEvaluator reference_nd_real(
        reference_topology, reference_numbering, reference_nd_fespace);
    fem::singular::TriangleEnrichedNDFieldEvaluator reference_nd_imaginary(
        reference_topology, reference_numbering, reference_nd_fespace);
    reference_nd_real.SetFromTrueDofs(reference_nd_zero);
    reference_nd_imaginary.SetFromTrueDofs(reference_nd_zero);
    TriangleSingularSurfacePostOperator reference_postoperator(
        MakeMSPostprocessing(), reference_material, reference_fespace);
    TriangleSingularSurfacePostOperator reference_longitudinal_postoperator(
        MakeDefaultPostprocessing(), reference_material, reference_nd_fespace);
    const double reference =
        reference_postoperator
            .MeasureElectrostatic(reference_real, reference_imaginary, 1.0, options)
            .front()
            .energy;
    const double reference_longitudinal =
        reference_longitudinal_postoperator
            .Measure(reference_nd_real, reference_nd_imaginary, 1.0, options,
                     &reference_real, &reference_imaginary)
            .front()
            .energy;

    REQUIRE(std::isfinite(reference));
    REQUIRE(reference > 0.0);
    REQUIRE(std::isfinite(reference_longitudinal));
    REQUIRE(reference_longitudinal > 0.0);
    CHECK_THAT(parallel_h1, WithinRel(reference, 2.0e-12));
    CHECK_THAT(parallel_nd, WithinRel(reference, 2.0e-12));
    CHECK_THAT(parallel_longitudinal, WithinRel(reference_longitudinal, 2.0e-12));
  }
}

TEST_CASE("Tetrahedral singular surface postprocessor uses edge-aligned quadrature",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const auto measure = [](double nu, int quadrature_order, double edge_cutoff = 0.0,
                          bool complete_set = false, bool fixed_subdivision = false)
  {
    auto serial_mesh = MakeTetrahedronSurfaceTestMesh();
    auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
    Mesh mesh(std::move(par_mesh));
    auto material = MakeSubstrateMaterial(mesh);

    mfem::H1_FECollection collection(1, 3);
    mfem::ParFiniteElementSpace fespace(&mesh.Get(), &collection);
    mfem::Array<int> element_vertices;
    mesh.Get().GetElementVertices(0, element_vertices);
    REQUIRE(element_vertices.Size() == 4);
    const auto local_node = [&element_vertices](int mesh_vertex)
    {
      const int *location =
          std::find(element_vertices.begin(), element_vertices.end(), mesh_vertex);
      REQUIRE(location != element_vertices.end());
      return static_cast<int>(location - element_vertices.begin());
    };
    const std::array<int, 3> face_nodes{local_node(0), local_node(1), local_node(2)};
    const std::array<int, 4> canonical_nodes{face_nodes[0], face_nodes[1], face_nodes[2],
                                             local_node(3)};
    auto edge_bases =
        fem::singular::EnumerateHigherOrderEdgeGradientBases(canonical_nodes, 1, nu);
    REQUIRE_FALSE(edge_bases.empty());
    std::vector<fem::singular::HigherOrderBasis> bases{edge_bases.front()};
    if (complete_set)
    {
      const auto append = [&bases](const auto &additional)
      { bases.insert(bases.end(), additional.begin(), additional.end()); };
      append(fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, 1, nu));
      append(fem::singular::EnumerateHigherOrderNodeGradientBases(
          {canonical_nodes[1], canonical_nodes[0], canonical_nodes[2], canonical_nodes[3]},
          1, nu));
      bases.insert(bases.end(), edge_bases.begin() + 1, edge_bases.end());
    }
    fem::singular::DofTopology topology;
    topology.h1_dofs.resize(bases.size());
    topology.elements.resize(1);
    for (std::size_t dof = 0; dof < bases.size(); dof++)
    {
      topology.h1_dofs[dof].family = bases[dof].family;
      topology.h1_dofs[dof].order = bases[dof].order;
      topology.elements[0].h1.push_back({dof, bases[dof]});
    }

    auto numbering = MakeH1Numbering(static_cast<int>(bases.size()));
    mfem::Vector combined(fespace.GetTrueVSize() + static_cast<int>(bases.size()));
    combined = 0.0;
    for (std::size_t dof = 0; dof < bases.size(); dof++)
    {
      combined[fespace.GetTrueVSize() + static_cast<int>(dof)] =
          1.0 / static_cast<double>(dof + 1);
    }
    mfem::Vector zero(combined.Size());
    zero = 0.0;
    fem::singular::EnrichedH1FieldEvaluator real(topology, numbering, fespace);
    fem::singular::EnrichedH1FieldEvaluator imaginary(topology, numbering, fespace);
    real.SetFromTrueDofs(combined);
    imaginary.SetFromTrueDofs(zero);
    const auto face_singularities = real.GetElementFaceSingularities(0, face_nodes);
    REQUIRE(face_singularities.size() == (complete_set ? 3 : 1));
    const auto edge = std::find_if(
        face_singularities.begin(), face_singularities.end(), [](const auto &feature)
        { return feature.type == fem::singular::TetrahedronFaceSingularityType::EDGE; });
    REQUIRE(edge != face_singularities.end());
    REQUIRE(edge->nodes[0] == std::min(face_nodes[0], face_nodes[1]));
    REQUIRE(edge->nodes[1] == std::max(face_nodes[0], face_nodes[1]));

    TetrahedronSingularSurfacePostOperator postoperator(MakeMSPostprocessing(edge_cutoff),
                                                        material, fespace);
    return postoperator
        .MeasureElectrostatic(real, imaginary, 1.0,
                              {quadrature_order, 1.0e-12, 2.0e-9, 12, fixed_subdivision, 6})
        .front()
        .energy;
  };

  const double integrable_order4 = measure(2.0 / 3.0, 4);
  const double integrable_order8 = measure(2.0 / 3.0, 8);
  CHECK(std::isfinite(integrable_order8));
  CHECK(integrable_order8 > 0.0);
  CHECK_THAT(integrable_order4, WithinRel(integrable_order8, 5.0e-9));
  CHECK_THROWS_AS(measure(0.5, 8), std::domain_error);

  const double cutoff_large_order4 = measure(0.5, 4, 1.0e-2);
  const double cutoff_large = measure(0.5, 8, 1.0e-2);
  const double cutoff_small_order4 = measure(0.5, 4, 1.0e-3);
  const double cutoff_small = measure(0.5, 8, 1.0e-3);
  CHECK(std::isfinite(cutoff_large));
  CHECK(cutoff_large > 0.0);
  CHECK(cutoff_small > cutoff_large);
  CHECK_THAT(cutoff_large_order4, WithinRel(cutoff_large, 5.0e-9));
  CHECK_THAT(cutoff_small_order4, WithinRel(cutoff_small, 5.0e-9));

  const double complete_integrable_order4 = measure(2.0 / 3.0, 4, 0.0, true);
  const double complete_integrable_order8 = measure(2.0 / 3.0, 8, 0.0, true);
  CHECK(std::isfinite(complete_integrable_order8));
  CHECK(complete_integrable_order8 > 0.0);
  CHECK_THAT(complete_integrable_order4, WithinRel(complete_integrable_order8, 5.0e-8));
  const double complete_integrable_fixed = measure(2.0 / 3.0, 8, 0.0, true, true);
  CHECK_THAT(complete_integrable_fixed, WithinRel(complete_integrable_order8, 2.0e-7));
  CHECK_THROWS_AS(measure(0.5, 8, 0.0, true), std::domain_error);
  const double complete_cutoff_large = measure(0.5, 8, 1.0e-2, true);
  const double complete_cutoff_small = measure(0.5, 8, 1.0e-3, true);
  CHECK(std::isfinite(complete_cutoff_large));
  CHECK(complete_cutoff_small > complete_cutoff_large);
}

TEST_CASE("Tetrahedral singular surface postprocessor maps a physical cutoff on a "
          "curved face",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr double curvature = 0.6;
  constexpr double cutoff = 0.05;
  auto serial_mesh = MakeTetrahedronSurfaceTestMesh();
  SetCurvedTetrahedronGeometry(*serial_mesh, curvature);
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));
  auto material = MakeSubstrateMaterial(mesh);

  mfem::H1_FECollection h1_collection(2, 3);
  mfem::ND_FECollection nd_collection(2, 3);
  mfem::ParFiniteElementSpace h1_fespace(&mesh.Get(), &h1_collection);
  mfem::ParFiniteElementSpace nd_fespace(&mesh.Get(), &nd_collection);
  mfem::FunctionCoefficient coefficient([](const mfem::Vector &point)
                                        { return 2.0 * point[2]; });
  mfem::ParGridFunction potential(&h1_fespace);
  potential.ProjectCoefficient(coefficient);
  mfem::Vector h1_standard_true_dofs(h1_fespace.GetTrueVSize());
  potential.ParallelProject(h1_standard_true_dofs);
  mfem::VectorFunctionCoefficient vector_coefficient(
      3,
      [](const mfem::Vector &, mfem::Vector &value)
      {
        value.SetSize(3);
        value = 0.0;
        value[2] = -2.0;
      });
  mfem::ParGridFunction field(&nd_fespace);
  field.ProjectCoefficient(vector_coefficient);
  mfem::Vector nd_standard_true_dofs(nd_fespace.GetTrueVSize());
  field.ParallelProject(nd_standard_true_dofs);

  mfem::Array<int> element_vertices;
  mesh.Get().GetElementVertices(0, element_vertices);
  REQUIRE(element_vertices.Size() == 4);
  const auto local_node = [&element_vertices](int mesh_vertex)
  {
    const int *location =
        std::find(element_vertices.begin(), element_vertices.end(), mesh_vertex);
    REQUIRE(location != element_vertices.end());
    return static_cast<int>(location - element_vertices.begin());
  };
  const std::array<int, 4> canonical_nodes{local_node(0), local_node(1), local_node(2),
                                           local_node(3)};
  const auto bases =
      fem::singular::EnumerateHigherOrderEdgeGradientBases(canonical_nodes, 1, 0.5);
  REQUIRE_FALSE(bases.empty());
  const auto basis = bases.front();
  REQUIRE(std::array<int, 2>{basis.nodes[0], basis.nodes[1]} ==
          std::array<int, 2>{local_node(0), local_node(1)});

  fem::singular::DofTopology topology;
  topology.h1_dofs.resize(1);
  topology.h1_dofs[0].family = basis.family;
  topology.h1_dofs[0].order = basis.order;
  topology.nd_dofs.resize(1);
  topology.nd_dofs[0].family = basis.family;
  topology.nd_dofs[0].order = basis.order;
  topology.elements.resize(1);
  topology.elements[0].h1 = {{0, basis}};
  topology.elements[0].nd = {{0, basis}};
  const auto h1_numbering = MakeH1Numbering(1);
  const auto nd_numbering = MakeNumbering(true);

  mfem::Vector h1_combined(h1_fespace.GetTrueVSize() + 1);
  h1_combined = 0.0;
  for (int dof = 0; dof < h1_standard_true_dofs.Size(); dof++)
  {
    h1_combined[dof] = h1_standard_true_dofs[dof];
  }
  mfem::Vector nd_combined(nd_fespace.GetTrueVSize() + 1);
  nd_combined = 0.0;
  for (int dof = 0; dof < nd_standard_true_dofs.Size(); dof++)
  {
    nd_combined[dof] = nd_standard_true_dofs[dof];
  }
  mfem::Vector h1_zero(h1_combined.Size());
  mfem::Vector nd_zero(nd_combined.Size());
  h1_zero = 0.0;
  nd_zero = 0.0;
  fem::singular::EnrichedH1FieldEvaluator h1_real(topology, h1_numbering, h1_fespace);
  fem::singular::EnrichedH1FieldEvaluator h1_imaginary(topology, h1_numbering, h1_fespace);
  fem::singular::EnrichedNDFieldEvaluator nd_real(topology, nd_numbering, nd_fespace);
  fem::singular::EnrichedNDFieldEvaluator nd_imaginary(topology, nd_numbering, nd_fespace);
  h1_real.SetFromTrueDofs(h1_combined);
  h1_imaginary.SetFromTrueDofs(h1_zero);
  nd_real.SetFromTrueDofs(nd_combined);
  nd_imaginary.SetFromTrueDofs(nd_zero);

  TetrahedronSingularSurfacePostOperator h1_postoperator(MakeDefaultPostprocessing(cutoff),
                                                         material, h1_fespace);
  TetrahedronSingularSurfacePostOperator nd_postoperator(MakeDefaultPostprocessing(cutoff),
                                                         material, nd_fespace);
  const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-12, 2.0e-10, 12};
  const auto h1_measurements =
      h1_postoperator.MeasureElectrostatic(h1_real, h1_imaginary, 1.0, options);
  const auto nd_measurements = nd_postoperator.Measure(nd_real, nd_imaginary, 1.0, options);
  REQUIRE(h1_measurements.size() == 1);
  REQUIRE(nd_measurements.size() == 1);
  const double exact = 0.04 * CurvedFaceAreaOutsideCutoff(curvature, cutoff);
  CHECK_THAT(h1_measurements[0].energy, WithinRel(exact, 2.0e-10));
  CHECK_THAT(nd_measurements[0].energy, WithinRel(exact, 2.0e-10));

  h1_combined[h1_fespace.GetTrueVSize()] = 0.37;
  nd_combined[nd_fespace.GetTrueVSize()] = -0.37;
  h1_real.SetFromTrueDofs(h1_combined);
  nd_real.SetFromTrueDofs(nd_combined);
  const double h1_singular_energy =
      h1_postoperator.MeasureElectrostatic(h1_real, h1_imaginary, 1.0, options)
          .front()
          .energy;
  const double nd_singular_energy =
      nd_postoperator.Measure(nd_real, nd_imaginary, 1.0, options).front().energy;
  REQUIRE(std::isfinite(h1_singular_energy));
  REQUIRE(h1_singular_energy > 0.0);
  CHECK_THAT(nd_singular_energy, WithinRel(h1_singular_energy, 2.0e-10));
}

TEST_CASE("Tetrahedral singular surface postprocessor maps intersecting physical "
          "cutoffs on a curved face",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr double curvature = 0.6;
  constexpr double cutoff = 0.05;
  auto serial_mesh = MakeTetrahedronSurfaceTestMesh();
  SetCurvedTetrahedronGeometry(*serial_mesh, curvature);
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));
  auto material = MakeSubstrateMaterial(mesh);

  mfem::H1_FECollection h1_collection(2, 3);
  mfem::ND_FECollection nd_collection(2, 3);
  mfem::ParFiniteElementSpace h1_fespace(&mesh.Get(), &h1_collection);
  mfem::ParFiniteElementSpace nd_fespace(&mesh.Get(), &nd_collection);
  mfem::FunctionCoefficient coefficient([](const mfem::Vector &point)
                                        { return 2.0 * point[2]; });
  mfem::ParGridFunction potential(&h1_fespace);
  potential.ProjectCoefficient(coefficient);
  mfem::Vector h1_standard_true_dofs(h1_fespace.GetTrueVSize());
  potential.ParallelProject(h1_standard_true_dofs);
  mfem::VectorFunctionCoefficient vector_coefficient(
      3,
      [](const mfem::Vector &, mfem::Vector &value)
      {
        value.SetSize(3);
        value = 0.0;
        value[2] = -2.0;
      });
  mfem::ParGridFunction field(&nd_fespace);
  field.ProjectCoefficient(vector_coefficient);
  mfem::Vector nd_standard_true_dofs(nd_fespace.GetTrueVSize());
  field.ParallelProject(nd_standard_true_dofs);

  mfem::Array<int> element_vertices;
  mesh.Get().GetElementVertices(0, element_vertices);
  REQUIRE(element_vertices.Size() == 4);
  const auto local_node = [&element_vertices](int mesh_vertex)
  {
    const int *location =
        std::find(element_vertices.begin(), element_vertices.end(), mesh_vertex);
    REQUIRE(location != element_vertices.end());
    return static_cast<int>(location - element_vertices.begin());
  };
  const std::array<int, 4> first_nodes{local_node(0), local_node(1), local_node(2),
                                       local_node(3)};
  const std::array<int, 4> second_nodes{local_node(0), local_node(2), local_node(1),
                                        local_node(3)};
  const auto first_bases =
      fem::singular::EnumerateHigherOrderEdgeGradientBases(first_nodes, 1, 0.5);
  const auto second_bases =
      fem::singular::EnumerateHigherOrderEdgeGradientBases(second_nodes, 1, 0.5);
  REQUIRE_FALSE(first_bases.empty());
  REQUIRE_FALSE(second_bases.empty());
  const std::array<fem::singular::HigherOrderBasis, 2> bases{first_bases.front(),
                                                             second_bases.front()};

  fem::singular::DofTopology topology;
  topology.h1_dofs.resize(bases.size());
  topology.nd_dofs.resize(bases.size());
  topology.elements.resize(1);
  for (std::size_t dof = 0; dof < bases.size(); dof++)
  {
    topology.h1_dofs[dof].family = bases[dof].family;
    topology.h1_dofs[dof].order = bases[dof].order;
    topology.nd_dofs[dof].family = bases[dof].family;
    topology.nd_dofs[dof].order = bases[dof].order;
    topology.elements[0].h1.push_back({dof, bases[dof]});
    topology.elements[0].nd.push_back({dof, bases[dof]});
  }
  const auto h1_numbering = MakeH1Numbering(static_cast<int>(bases.size()));
  auto nd_numbering = MakeNumbering(false);
  nd_numbering.nd.global_local_size = static_cast<int>(bases.size());
  nd_numbering.nd.local_size = static_cast<int>(bases.size());
  nd_numbering.nd.global_size = static_cast<int>(bases.size());
  nd_numbering.nd.owned_size = static_cast<int>(bases.size());
  nd_numbering.nd.owner.assign(bases.size(), 0);
  nd_numbering.nd.local_to_true.resize(bases.size());
  std::iota(nd_numbering.nd.local_to_true.begin(), nd_numbering.nd.local_to_true.end(), 0);

  mfem::Vector h1_combined(h1_fespace.GetTrueVSize() + static_cast<int>(bases.size()));
  h1_combined = 0.0;
  for (int dof = 0; dof < h1_standard_true_dofs.Size(); dof++)
  {
    h1_combined[dof] = h1_standard_true_dofs[dof];
  }
  mfem::Vector nd_combined(nd_fespace.GetTrueVSize() + static_cast<int>(bases.size()));
  nd_combined = 0.0;
  for (int dof = 0; dof < nd_standard_true_dofs.Size(); dof++)
  {
    nd_combined[dof] = nd_standard_true_dofs[dof];
  }
  mfem::Vector h1_zero(h1_combined.Size());
  mfem::Vector nd_zero(nd_combined.Size());
  h1_zero = 0.0;
  nd_zero = 0.0;
  fem::singular::EnrichedH1FieldEvaluator h1_real(topology, h1_numbering, h1_fespace);
  fem::singular::EnrichedH1FieldEvaluator h1_imaginary(topology, h1_numbering, h1_fespace);
  fem::singular::EnrichedNDFieldEvaluator nd_real(topology, nd_numbering, nd_fespace);
  fem::singular::EnrichedNDFieldEvaluator nd_imaginary(topology, nd_numbering, nd_fespace);
  h1_real.SetFromTrueDofs(h1_combined);
  h1_imaginary.SetFromTrueDofs(h1_zero);
  nd_real.SetFromTrueDofs(nd_combined);
  nd_imaginary.SetFromTrueDofs(nd_zero);

  TetrahedronSingularSurfacePostOperator h1_postoperator(MakeDefaultPostprocessing(cutoff),
                                                         material, h1_fespace);
  TetrahedronSingularSurfacePostOperator nd_postoperator(MakeDefaultPostprocessing(cutoff),
                                                         material, nd_fespace);
  const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-12, 2.0e-9, 12};
  const double exact = 0.04 * CurvedFaceAreaOutsideTwoCutoffs(curvature, cutoff);
  const double h1_energy =
      h1_postoperator.MeasureElectrostatic(h1_real, h1_imaginary, 1.0, options)
          .front()
          .energy;
  const double nd_energy =
      nd_postoperator.Measure(nd_real, nd_imaginary, 1.0, options).front().energy;
  CHECK_THAT(h1_energy, WithinRel(exact, 2.0e-8));
  CHECK_THAT(nd_energy, WithinRel(exact, 2.0e-8));
  TetrahedronSingularSurfacePostOperator rejected_postoperator(
      MakeDefaultPostprocessing(0.5), material, h1_fespace);
  CHECK_THROWS_AS(
      rejected_postoperator.MeasureElectrostatic(h1_real, h1_imaginary, 1.0, options),
      std::domain_error);

  h1_combined[h1_fespace.GetTrueVSize()] = 0.37;
  h1_combined[h1_fespace.GetTrueVSize() + 1] = -0.23;
  nd_combined[nd_fespace.GetTrueVSize()] = -0.37;
  nd_combined[nd_fespace.GetTrueVSize() + 1] = 0.23;
  h1_real.SetFromTrueDofs(h1_combined);
  nd_real.SetFromTrueDofs(nd_combined);
  const double h1_singular_energy =
      h1_postoperator.MeasureElectrostatic(h1_real, h1_imaginary, 1.0, options)
          .front()
          .energy;
  const double nd_singular_energy =
      nd_postoperator.Measure(nd_real, nd_imaginary, 1.0, options).front().energy;
  REQUIRE(std::isfinite(h1_singular_energy));
  REQUIRE(h1_singular_energy > 0.0);
  CHECK_THAT(nd_singular_energy, WithinRel(h1_singular_energy, 2.0e-9));
}

TEST_CASE("Tetrahedral singular surface postprocessor maps three physical cutoffs on "
          "a bubble-curved face",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr double curvature = 0.6;
  constexpr double cutoff = 0.05;
  auto serial_mesh = MakeTetrahedronSurfaceTestMesh();
  SetBubbleCurvedTetrahedronGeometry(*serial_mesh, curvature);
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));
  auto material = MakeSubstrateMaterial(mesh);

  mfem::H1_FECollection h1_collection(3, 3);
  mfem::ND_FECollection nd_collection(3, 3);
  mfem::ParFiniteElementSpace h1_fespace(&mesh.Get(), &h1_collection);
  mfem::ParFiniteElementSpace nd_fespace(&mesh.Get(), &nd_collection);
  mfem::FunctionCoefficient coefficient([](const mfem::Vector &point)
                                        { return 2.0 * point[2]; });
  mfem::ParGridFunction potential(&h1_fespace);
  potential.ProjectCoefficient(coefficient);
  mfem::Vector h1_standard_true_dofs(h1_fespace.GetTrueVSize());
  potential.ParallelProject(h1_standard_true_dofs);
  mfem::VectorFunctionCoefficient vector_coefficient(
      3,
      [](const mfem::Vector &, mfem::Vector &value)
      {
        value.SetSize(3);
        value = 0.0;
        value[2] = -2.0;
      });
  mfem::ParGridFunction field(&nd_fespace);
  field.ProjectCoefficient(vector_coefficient);
  mfem::Vector nd_standard_true_dofs(nd_fespace.GetTrueVSize());
  field.ParallelProject(nd_standard_true_dofs);

  mfem::Array<int> element_vertices;
  mesh.Get().GetElementVertices(0, element_vertices);
  REQUIRE(element_vertices.Size() == 4);
  const auto local_node = [&element_vertices](int mesh_vertex)
  {
    const int *location =
        std::find(element_vertices.begin(), element_vertices.end(), mesh_vertex);
    REQUIRE(location != element_vertices.end());
    return static_cast<int>(location - element_vertices.begin());
  };
  const std::array<std::array<int, 4>, 3> canonical_nodes{
      std::array<int, 4>{local_node(0), local_node(1), local_node(2), local_node(3)},
      std::array<int, 4>{local_node(0), local_node(2), local_node(1), local_node(3)},
      std::array<int, 4>{local_node(1), local_node(2), local_node(0), local_node(3)}};
  std::array<fem::singular::HigherOrderBasis, 3> bases;
  for (std::size_t edge = 0; edge < canonical_nodes.size(); edge++)
  {
    const auto edge_bases =
        fem::singular::EnumerateHigherOrderEdgeGradientBases(canonical_nodes[edge], 1, 0.5);
    REQUIRE_FALSE(edge_bases.empty());
    bases[edge] = edge_bases.front();
  }

  fem::singular::DofTopology topology;
  topology.h1_dofs.resize(bases.size());
  topology.nd_dofs.resize(bases.size());
  topology.elements.resize(1);
  for (std::size_t dof = 0; dof < bases.size(); dof++)
  {
    topology.h1_dofs[dof].family = bases[dof].family;
    topology.h1_dofs[dof].order = bases[dof].order;
    topology.nd_dofs[dof].family = bases[dof].family;
    topology.nd_dofs[dof].order = bases[dof].order;
    topology.elements[0].h1.push_back({dof, bases[dof]});
    topology.elements[0].nd.push_back({dof, bases[dof]});
  }
  const auto h1_numbering = MakeH1Numbering(static_cast<int>(bases.size()));
  auto nd_numbering = MakeNumbering(false);
  nd_numbering.nd.global_local_size = static_cast<int>(bases.size());
  nd_numbering.nd.local_size = static_cast<int>(bases.size());
  nd_numbering.nd.global_size = static_cast<int>(bases.size());
  nd_numbering.nd.owned_size = static_cast<int>(bases.size());
  nd_numbering.nd.owner.assign(bases.size(), 0);
  nd_numbering.nd.local_to_true.resize(bases.size());
  std::iota(nd_numbering.nd.local_to_true.begin(), nd_numbering.nd.local_to_true.end(), 0);

  mfem::Vector h1_combined(h1_fespace.GetTrueVSize() + static_cast<int>(bases.size()));
  h1_combined = 0.0;
  for (int dof = 0; dof < h1_standard_true_dofs.Size(); dof++)
  {
    h1_combined[dof] = h1_standard_true_dofs[dof];
  }
  mfem::Vector nd_combined(nd_fespace.GetTrueVSize() + static_cast<int>(bases.size()));
  nd_combined = 0.0;
  for (int dof = 0; dof < nd_standard_true_dofs.Size(); dof++)
  {
    nd_combined[dof] = nd_standard_true_dofs[dof];
  }
  mfem::Vector h1_zero(h1_combined.Size());
  mfem::Vector nd_zero(nd_combined.Size());
  h1_zero = 0.0;
  nd_zero = 0.0;
  fem::singular::EnrichedH1FieldEvaluator h1_real(topology, h1_numbering, h1_fespace);
  fem::singular::EnrichedH1FieldEvaluator h1_imaginary(topology, h1_numbering, h1_fespace);
  fem::singular::EnrichedNDFieldEvaluator nd_real(topology, nd_numbering, nd_fespace);
  fem::singular::EnrichedNDFieldEvaluator nd_imaginary(topology, nd_numbering, nd_fespace);
  h1_real.SetFromTrueDofs(h1_combined);
  h1_imaginary.SetFromTrueDofs(h1_zero);
  nd_real.SetFromTrueDofs(nd_combined);
  nd_imaginary.SetFromTrueDofs(nd_zero);

  TetrahedronSingularSurfacePostOperator h1_postoperator(MakeDefaultPostprocessing(cutoff),
                                                         material, h1_fespace);
  TetrahedronSingularSurfacePostOperator nd_postoperator(MakeDefaultPostprocessing(cutoff),
                                                         material, nd_fespace);
  const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-12, 2.0e-9, 12};
  const double oracle_coarse =
      0.04 * BubbleCurvedFaceAreaOutsideThreeCutoffs(curvature, cutoff, 768);
  const double oracle_fine =
      0.04 * BubbleCurvedFaceAreaOutsideThreeCutoffs(curvature, cutoff, 1536);
  const double h1_energy =
      h1_postoperator.MeasureElectrostatic(h1_real, h1_imaginary, 1.0, options)
          .front()
          .energy;
  const double nd_energy =
      nd_postoperator.Measure(nd_real, nd_imaginary, 1.0, options).front().energy;
  const double oracle_difference = std::abs(oracle_fine - oracle_coarse);
  CHECK(oracle_difference < 1.0e-7 * oracle_fine);
  CHECK(std::abs(h1_energy - oracle_fine) <= oracle_difference);
  CHECK(std::abs(nd_energy - oracle_fine) <= oracle_difference);

  const std::array<double, 3> coefficients{0.37, -0.23, 0.19};
  for (std::size_t dof = 0; dof < coefficients.size(); dof++)
  {
    h1_combined[h1_fespace.GetTrueVSize() + static_cast<int>(dof)] = coefficients[dof];
    nd_combined[nd_fespace.GetTrueVSize() + static_cast<int>(dof)] = -coefficients[dof];
  }
  h1_real.SetFromTrueDofs(h1_combined);
  nd_real.SetFromTrueDofs(nd_combined);
  const double h1_singular_energy =
      h1_postoperator.MeasureElectrostatic(h1_real, h1_imaginary, 1.0, options)
          .front()
          .energy;
  const double nd_singular_energy =
      nd_postoperator.Measure(nd_real, nd_imaginary, 1.0, options).front().energy;
  REQUIRE(std::isfinite(h1_singular_energy));
  REQUIRE(h1_singular_energy > 0.0);
  CHECK_THAT(nd_singular_energy, WithinRel(h1_singular_energy, 2.0e-9));
}

TEST_CASE("Tetrahedral H1 and ND gradient enrichments give identical surface energy",
          "[singularsurface][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const auto measure = [](double nu, double edge_cutoff)
  {
    auto serial_mesh = MakeTetrahedronSurfaceTestMesh();
    auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
    Mesh mesh(std::move(par_mesh));
    auto material = MakeSubstrateMaterial(mesh);

    mfem::Array<int> element_vertices;
    mesh.Get().GetElementVertices(0, element_vertices);
    REQUIRE(element_vertices.Size() == 4);
    const auto local_node = [&element_vertices](int mesh_vertex)
    {
      const int *location =
          std::find(element_vertices.begin(), element_vertices.end(), mesh_vertex);
      REQUIRE(location != element_vertices.end());
      return static_cast<int>(location - element_vertices.begin());
    };
    const std::array<int, 4> canonical_nodes{local_node(0), local_node(1), local_node(2),
                                             local_node(3)};
    const auto bases =
        fem::singular::EnumerateHigherOrderEdgeGradientBases(canonical_nodes, 1, nu);
    REQUIRE_FALSE(bases.empty());
    const auto basis = bases.front();

    fem::singular::DofTopology topology;
    topology.h1_dofs.resize(1);
    topology.h1_dofs[0].family = basis.family;
    topology.h1_dofs[0].order = basis.order;
    topology.nd_dofs.resize(1);
    topology.nd_dofs[0].family = basis.family;
    topology.nd_dofs[0].order = basis.order;
    topology.elements.resize(1);
    topology.elements[0].h1 = {{0, basis}};
    topology.elements[0].nd = {{0, basis}};

    auto h1_numbering = MakeH1Numbering(1);
    auto nd_numbering = MakeNumbering(true);
    mfem::H1_FECollection h1_collection(1, 3);
    mfem::ND_FECollection nd_collection(1, 3);
    mfem::ParFiniteElementSpace h1_fespace(&mesh.Get(), &h1_collection);
    mfem::ParFiniteElementSpace nd_fespace(&mesh.Get(), &nd_collection);
    mfem::Vector h1_combined(h1_fespace.GetTrueVSize() + 1);
    mfem::Vector nd_combined(nd_fespace.GetTrueVSize() + 1);
    h1_combined = 0.0;
    nd_combined = 0.0;
    h1_combined[h1_fespace.GetTrueVSize()] = 0.73;
    nd_combined[nd_fespace.GetTrueVSize()] = 0.73;
    mfem::Vector h1_zero(h1_combined.Size());
    mfem::Vector nd_zero(nd_combined.Size());
    h1_zero = 0.0;
    nd_zero = 0.0;

    fem::singular::EnrichedH1FieldEvaluator h1_real(topology, h1_numbering, h1_fespace);
    fem::singular::EnrichedH1FieldEvaluator h1_imaginary(topology, h1_numbering,
                                                         h1_fespace);
    fem::singular::EnrichedNDFieldEvaluator nd_real(topology, nd_numbering, nd_fespace);
    fem::singular::EnrichedNDFieldEvaluator nd_imaginary(topology, nd_numbering,
                                                         nd_fespace);
    h1_real.SetFromTrueDofs(h1_combined);
    h1_imaginary.SetFromTrueDofs(h1_zero);
    nd_real.SetFromTrueDofs(nd_combined);
    nd_imaginary.SetFromTrueDofs(nd_zero);

    TetrahedronSingularSurfacePostOperator h1_postoperator(
        MakeMSPostprocessing(edge_cutoff), material, h1_fespace);
    TetrahedronSingularSurfacePostOperator nd_postoperator(
        MakeMSPostprocessing(edge_cutoff), material, nd_fespace);
    const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-12, 2.0e-9, 12};
    const double h1_energy =
        h1_postoperator.MeasureElectrostatic(h1_real, h1_imaginary, 1.0, options)
            .front()
            .energy;
    const double nd_energy =
        nd_postoperator.Measure(nd_real, nd_imaginary, 1.0, options).front().energy;
    return std::pair{h1_energy, nd_energy};
  };

  for (const auto &[nu, edge_cutoff] : {std::pair{2.0 / 3.0, 0.0}, std::pair{0.5, 1.0e-3}})
  {
    const auto [h1_energy, nd_energy] = measure(nu, edge_cutoff);
    CHECK(std::isfinite(h1_energy));
    CHECK(h1_energy > 0.0);
    CHECK_THAT(nd_energy, WithinRel(h1_energy, 2.0e-12));
  }
}

}  // namespace palace
