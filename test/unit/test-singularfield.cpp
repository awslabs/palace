// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/singularfield.hpp"
#include "utils/communication.hpp"

using namespace Catch::Matchers;

namespace palace
{

namespace
{

using fem::singular::BarycentricPoint;
using fem::singular::H1FieldValue;
using fem::singular::Vector3;

constexpr std::array<Vector3, 4> kVertices{Vector3{0.1, -0.2, 0.3}, Vector3{1.4, 0.0, 0.4},
                                           Vector3{0.2, 1.1, 0.1}, Vector3{-0.1, 0.2, 1.6}};
constexpr BarycentricPoint kPhysicalBarycentricPoint{0.37, 0.21, 0.18, 0.24};
constexpr double kSingularCoefficient = -0.43;

void TranslateQuadraticGeometry(mfem::Mesh &mesh, const std::array<double, 3> &translation)
{
  const int dimension = mesh.SpaceDimension();
  mesh.SetCurvature(2, false, dimension, mfem::Ordering::byVDIM);
  mfem::VectorFunctionCoefficient geometry(
      dimension,
      [dimension, translation](const mfem::Vector &x, mfem::Vector &value)
      {
        value.SetSize(dimension);
        for (int d = 0; d < dimension; d++)
        {
          value[d] = x[d] + translation[d];
        }
      });
  mesh.GetNodes()->ProjectCoefficient(geometry);
  mesh.NodesUpdated();
}

double StandardPotential(const mfem::Vector &x)
{
  return 1.2 + 2.0 * x[0] - 3.0 * x[1] + 0.5 * x[2];
}

H1FieldValue EvaluatePermutation(const std::array<int, 4> &local_to_physical)
{
  mfem::Mesh serial_mesh(3, 4, 1, 0, 3);
  for (const auto &vertex : kVertices)
  {
    serial_mesh.AddVertex(vertex.data());
  }
  serial_mesh.AddTet(local_to_physical.data(), 1);
  serial_mesh.FinalizeTopology();
  serial_mesh.Finalize(false, false);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  std::array<int, 4> physical_to_local;
  BarycentricPoint lambda;
  for (int local = 0; local < 4; local++)
  {
    physical_to_local[local_to_physical[local]] = local;
    lambda[local] = kPhysicalBarycentricPoint[local_to_physical[local]];
  }
  const std::array<int, 4> canonical_nodes{physical_to_local[0], physical_to_local[1],
                                           physical_to_local[2], physical_to_local[3]};
  const auto bases =
      fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, 2, 0.5);
  REQUIRE_FALSE(bases.empty());
  const auto &basis = bases.front();

  fem::singular::DofTopology topology;
  topology.h1_dofs.resize(1);
  topology.h1_dofs[0].family = basis.family;
  topology.h1_dofs[0].order = basis.order;
  topology.elements.resize(1);
  topology.elements[0].h1.push_back({0, basis});

  fem::singular::ParallelDofNumbering numbering;
  numbering.h1.global_local_size = 1;
  numbering.h1.local_offset = 0;
  numbering.h1.local_size = 1;
  numbering.h1.global_size = 1;
  numbering.h1.owned_offset = 0;
  numbering.h1.owned_size = 1;
  numbering.h1.owner = {0};
  numbering.h1.local_to_true = {0};

  mfem::H1_FECollection collection(2, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &collection);
  mfem::FunctionCoefficient standard_coefficient(StandardPotential);
  mfem::ParGridFunction standard_field(&fespace);
  standard_field.ProjectCoefficient(standard_coefficient);
  mfem::Vector standard_true_dofs(fespace.GetTrueVSize());
  standard_field.ParallelProject(standard_true_dofs);

  mfem::Vector combined_true_dofs(standard_true_dofs.Size() + 1);
  for (int i = 0; i < standard_true_dofs.Size(); i++)
  {
    combined_true_dofs[i] = standard_true_dofs[i];
  }
  combined_true_dofs[standard_true_dofs.Size()] = kSingularCoefficient;

  fem::singular::EnrichedH1FieldEvaluator evaluator(topology, numbering, fespace);
  CHECK_THROWS_AS(evaluator.GetOwnedCoefficientDiagnostics(), std::logic_error);
  const auto face_singularities = evaluator.GetElementFaceSingularities(0, {0, 1, 2});
  REQUIRE(face_singularities.size() == 1);
  CHECK(face_singularities[0].type == fem::singular::TetrahedronFaceSingularityType::NODE);
  CHECK(face_singularities[0].nodes[0] == canonical_nodes[0]);
  CHECK(face_singularities[0].nodes[1] == -1);
  CHECK(face_singularities[0].nu == basis.nu);
  std::array<int, 3> opposite_face{};
  int opposite_node = 0;
  for (int node = 0; node < 4; node++)
  {
    if (node != canonical_nodes[0])
    {
      opposite_face[opposite_node++] = node;
    }
  }
  CHECK(evaluator.GetElementFaceSingularities(0, opposite_face).empty());
  CHECK_THROWS_AS(evaluator.GetElementFaceSingularities(0, {0, 0, 1}),
                  std::invalid_argument);
  evaluator.SetFromTrueDofs(combined_true_dofs);
  const auto coefficient_diagnostics = evaluator.GetOwnedCoefficientDiagnostics();
  REQUIRE(coefficient_diagnostics.size() == 1);
  CHECK(coefficient_diagnostics[0].true_dof == 0);
  CHECK(coefficient_diagnostics[0].key.family == basis.family);
  CHECK(coefficient_diagnostics[0].key.order == basis.order);
  CHECK(coefficient_diagnostics[0].exponent == basis.nu);
  CHECK(coefficient_diagnostics[0].coefficient == kSingularCoefficient);
  mfem::L2_FECollection sample_collection(2, 3, mfem::BasisType::GaussLegendre,
                                          mfem::FiniteElement::VALUE);
  mfem::ParFiniteElementSpace sampled_potential_space(&mesh, &sample_collection);
  mfem::ParFiniteElementSpace sampled_electric_space(&mesh, &sample_collection, 3,
                                                     mfem::Ordering::byVDIM);
  mfem::ParGridFunction sampled_potential(&sampled_potential_space);
  mfem::ParGridFunction sampled_electric(&sampled_electric_space);
  evaluator.ProjectToDiscontinuousGridFunctions(sampled_potential, sampled_electric);
  const auto &sample_points = sampled_potential_space.GetFE(0)->GetNodes();
  mfem::Vector sampled_electric_value(3);
  for (int i = 0; i < sample_points.GetNPoints(); i++)
  {
    const auto &sample_point = sample_points.IntPoint(i);
    const auto exact_sample = evaluator.Evaluate(0, sample_point);
    CHECK_THAT(sampled_potential.GetValue(0, sample_point),
               WithinAbs(exact_sample.potential, 2.0e-13));
    sampled_electric.GetVectorValue(0, sample_point, sampled_electric_value);
    for (int d = 0; d < 3; d++)
    {
      CHECK_THAT(sampled_electric_value[d], WithinAbs(-exact_sample.gradient[d], 2.0e-12));
    }
  }

  mfem::IntegrationPoint point;
  point.Set3(lambda[1], lambda[2], lambda[3]);
  const auto evaluated = evaluator.Evaluate(0, point);

  auto *transformation = mesh.GetElementTransformation(0);
  REQUIRE(transformation);
  transformation->SetIntPoint(&point);
  mfem::Vector physical_point(3);
  transformation->Transform(point, physical_point);
  double jacobian_determinant;
  const auto grad_lambda =
      fem::singular::GetAffineBarycentricGradients(*transformation, jacobian_determinant);
  const auto singular_gradient = fem::singular::EvaluateHigherOrderNodeGradient(
      lambda, grad_lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2], basis.nodes[3],
      basis.interpolation_indices, basis.order, basis.nu);
  const double singular_potential = fem::singular::EvaluateHigherOrderNodeGradientPotential(
      lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2], basis.nodes[3],
      basis.interpolation_indices, basis.order, basis.nu);

  CHECK_THAT(evaluated.potential, WithinAbs(StandardPotential(physical_point) +
                                                kSingularCoefficient * singular_potential,
                                            2.0e-13));
  const Vector3 standard_gradient{2.0, -3.0, 0.5};
  for (int d = 0; d < 3; d++)
  {
    CHECK_THAT(
        evaluated.gradient[d],
        WithinAbs(standard_gradient[d] + kSingularCoefficient * singular_gradient.value[d],
                  2.0e-12));
  }

  BarycentricPoint face_lambda{0.31, 0.27, 0.42, 0.0};
  mfem::IntegrationPoint face_point;
  face_point.Set3(face_lambda[1], face_lambda[2], face_lambda[3]);
  const auto face_value = evaluator.EvaluateClosure(0, face_point);
  transformation->SetIntPoint(&face_point);
  transformation->Transform(face_point, physical_point);
  const auto face_grad_lambda = fem::singular::GetBarycentricGradients(
      *transformation, face_point, jacobian_determinant);
  const auto face_singular =
      fem::singular::EvaluateHigherOrderBasis(face_lambda, face_grad_lambda, basis);
  const double face_singular_potential =
      fem::singular::EvaluateHigherOrderGradientPotential(face_lambda, basis);
  CHECK_THAT(face_value.potential,
             WithinAbs(StandardPotential(physical_point) +
                           kSingularCoefficient * face_singular_potential,
                       2.0e-13));
  for (int d = 0; d < 3; d++)
  {
    CHECK_THAT(
        face_value.gradient[d],
        WithinAbs(standard_gradient[d] + kSingularCoefficient * face_singular.value[d],
                  2.0e-12));
  }
  mfem::IntegrationPoint outside_point;
  outside_point.Set3(0.7, 0.7, 0.0);
  CHECK_THROWS_AS(evaluator.EvaluateClosure(0, outside_point), std::invalid_argument);
  return evaluated;
}

}  // namespace

TEST_CASE("Combined singular H1 field evaluation is affine and permutation covariant",
          "[singularfield][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const auto canonical = EvaluatePermutation({0, 1, 2, 3});
  const auto permuted = EvaluatePermutation({1, 2, 0, 3});
  CHECK_THAT(permuted.potential, WithinAbs(canonical.potential, 2.0e-13));
  for (int d = 0; d < 3; d++)
  {
    CHECK_THAT(permuted.gradient[d], WithinAbs(canonical.gradient[d], 2.0e-12));
  }
}

TEST_CASE("Combined singular tetrahedral ND field evaluates value and curl",
          "[singularfield][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  mfem::Mesh serial_mesh(3, 4, 1, 0, 3);
  for (const auto &vertex : kVertices)
  {
    serial_mesh.AddVertex(vertex.data());
  }
  serial_mesh.AddTet(0, 1, 2, 3, 1);
  serial_mesh.FinalizeTopology();
  serial_mesh.Finalize(false, false);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  const auto gradient_bases =
      fem::singular::EnumerateHigherOrderNodeGradientBases({0, 1, 2, 3}, 2, 0.5);
  const auto rotational_bases =
      fem::singular::EnumerateHigherOrderNodeRotationalBases({0, 1, 2, 3}, 2, 0.5);
  REQUIRE_FALSE(gradient_bases.empty());
  REQUIRE_FALSE(rotational_bases.empty());
  const std::array<fem::singular::HigherOrderBasis, 2> bases{gradient_bases.front(),
                                                             rotational_bases.front()};

  fem::singular::DofTopology topology;
  topology.nd_dofs.resize(bases.size());
  topology.elements.resize(1);
  for (std::size_t dof = 0; dof < bases.size(); dof++)
  {
    topology.nd_dofs[dof].family = bases[dof].family;
    topology.nd_dofs[dof].order = bases[dof].order;
    topology.elements[0].nd.push_back({dof, bases[dof]});
  }

  fem::singular::ParallelDofNumbering numbering;
  numbering.nd.global_local_size = bases.size();
  numbering.nd.local_offset = 0;
  numbering.nd.local_size = bases.size();
  numbering.nd.global_size = bases.size();
  numbering.nd.owned_offset = 0;
  numbering.nd.owned_size = bases.size();
  numbering.nd.owner = {0, 0};
  numbering.nd.local_to_true = {0, 1};

  mfem::ND_FECollection collection(2, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &collection);
  mfem::VectorFunctionCoefficient standard_coefficient(
      3,
      [](const mfem::Vector &x, mfem::Vector &value)
      {
        value.SetSize(3);
        value[0] = 1.0 + 2.0 * x[1];
        value[1] = -0.5 + 3.0 * x[2];
        value[2] = 0.25 + 4.0 * x[0];
      });
  mfem::ParGridFunction standard_field(&fespace);
  standard_field.ProjectCoefficient(standard_coefficient);
  mfem::Vector standard_true_dofs(fespace.GetTrueVSize());
  standard_field.ParallelProject(standard_true_dofs);

  constexpr std::array<double, 2> coefficients{-0.43, 0.27};
  mfem::Vector combined(standard_true_dofs.Size() + coefficients.size());
  for (int i = 0; i < standard_true_dofs.Size(); i++)
  {
    combined[i] = standard_true_dofs[i];
  }
  for (std::size_t i = 0; i < coefficients.size(); i++)
  {
    combined[standard_true_dofs.Size() + static_cast<int>(i)] = coefficients[i];
  }

  fem::singular::EnrichedNDFieldEvaluator evaluator(topology, numbering, fespace);
  const auto face_singularities = evaluator.GetElementFaceSingularities(0, {0, 1, 2});
  REQUIRE(face_singularities.size() == 1);
  CHECK(face_singularities[0].type == fem::singular::TetrahedronFaceSingularityType::NODE);
  CHECK(face_singularities[0].nodes == std::array<int, 2>{0, -1});
  CHECK(face_singularities[0].nu == 0.5);
  evaluator.SetFromTrueDofs(combined);

  mfem::IntegrationPoint point;
  point.Set3(kPhysicalBarycentricPoint[1], kPhysicalBarycentricPoint[2],
             kPhysicalBarycentricPoint[3]);
  const auto value = evaluator.Evaluate(0, point);
  auto *transformation = mesh.GetElementTransformation(0);
  REQUIRE(transformation);
  transformation->SetIntPoint(&point);
  mfem::Vector physical_point(3);
  transformation->Transform(point, physical_point);
  double jacobian_determinant;
  const auto grad_lambda =
      fem::singular::GetBarycentricGradients(*transformation, point, jacobian_determinant);
  fem::singular::Vector3 expected_value{1.0 + 2.0 * physical_point[1],
                                        -0.5 + 3.0 * physical_point[2],
                                        0.25 + 4.0 * physical_point[0]};
  fem::singular::Vector3 expected_curl{-3.0, -4.0, -2.0};
  for (std::size_t i = 0; i < bases.size(); i++)
  {
    const auto basis = fem::singular::EvaluateHigherOrderBasis(kPhysicalBarycentricPoint,
                                                               grad_lambda, bases[i]);
    for (int d = 0; d < 3; d++)
    {
      expected_value[d] += coefficients[i] * basis.value[d];
      expected_curl[d] += coefficients[i] * basis.curl[d];
    }
  }
  for (int d = 0; d < 3; d++)
  {
    CHECK_THAT(value.value[d], WithinAbs(expected_value[d], 3.0e-12));
    CHECK_THAT(value.curl[d], WithinAbs(expected_curl[d], 3.0e-12));
  }

  mfem::IntegrationPoint face_point;
  face_point.Set3(0.23, 0.31, 0.0);
  const auto face_value = evaluator.EvaluateClosure(0, face_point);
  for (double component : face_value.value)
  {
    CHECK(std::isfinite(component));
  }
  for (double component : face_value.curl)
  {
    CHECK(std::isfinite(component));
  }
}

TEST_CASE("Combined singular H1 energy uses exact quadrature away from enrichment",
          "[singularfield][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  mfem::Mesh serial_mesh(3, 8, 2, 0, 3);
  constexpr std::array<Vector3, 8> vertices{Vector3{0.0, 0.0, 0.0}, Vector3{1.0, 0.0, 0.0},
                                            Vector3{0.0, 1.0, 0.0}, Vector3{0.0, 0.0, 1.0},
                                            Vector3{2.0, 0.0, 0.0}, Vector3{3.0, 0.0, 0.0},
                                            Vector3{2.0, 1.0, 0.0}, Vector3{2.0, 0.0, 1.0}};
  for (const auto &vertex : vertices)
  {
    serial_mesh.AddVertex(vertex.data());
  }
  constexpr std::array<int, 4> first{0, 1, 2, 3};
  constexpr std::array<int, 4> second{4, 5, 6, 7};
  serial_mesh.AddTet(first.data(), 1);
  serial_mesh.AddTet(second.data(), 1);
  serial_mesh.FinalizeTopology();
  serial_mesh.Finalize(false, false);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  const auto basis =
      fem::singular::EnumerateHigherOrderNodeGradientBases({0, 1, 2, 3}, 1, 0.5).front();
  fem::singular::DofTopology topology;
  topology.h1_dofs.resize(1);
  topology.h1_dofs[0].family = basis.family;
  topology.h1_dofs[0].order = basis.order;
  topology.elements.resize(2);
  topology.elements[0].h1.push_back({0, basis});

  fem::singular::ParallelDofNumbering numbering;
  numbering.h1.global_local_size = 1;
  numbering.h1.local_offset = 0;
  numbering.h1.local_size = 1;
  numbering.h1.global_size = 1;
  numbering.h1.owned_offset = 0;
  numbering.h1.owned_size = 1;
  numbering.h1.owner = {0};
  numbering.h1.local_to_true = {0};

  mfem::H1_FECollection collection(2, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &collection);
  mfem::FunctionCoefficient standard_coefficient(StandardPotential);
  mfem::ParGridFunction standard_field(&fespace);
  standard_field.ProjectCoefficient(standard_coefficient);
  mfem::Vector standard_true_dofs(fespace.GetTrueVSize());
  standard_field.ParallelProject(standard_true_dofs);
  mfem::Vector combined_true_dofs(standard_true_dofs.Size() + 1);
  for (int i = 0; i < standard_true_dofs.Size(); i++)
  {
    combined_true_dofs[i] = standard_true_dofs[i];
  }
  combined_true_dofs[standard_true_dofs.Size()] = kSingularCoefficient;

  fem::singular::EnrichedH1FieldEvaluator evaluator(topology, numbering, fespace);
  const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-8, 1.0e-8, 4};
  CHECK_THROWS_AS(evaluator.GetOwnedCoefficientDiagnostics(), std::logic_error);
  CHECK_THROWS_AS(evaluator.IntegrateElementGradientEnergy(1, 3.25, options),
                  std::logic_error);
  evaluator.SetFromTrueDofs(combined_true_dofs);
  const auto integral = evaluator.IntegrateElementGradientEnergy(1, 3.25, options);
  CHECK(integral.converged);
  CHECK(integral.estimated_absolute_error == 0.0);
  CHECK(integral.leaf_count == 1);
  CHECK(integral.maximum_subdivision_depth == 0);
  CHECK_THAT(integral.value, WithinAbs(3.25 * (4.0 + 9.0 + 0.25) / 6.0, 2.0e-13));
}

TEST_CASE("Combined singular H1 edge-ray fit recovers the Meixner exponent",
          "[singularfield][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  mfem::Mesh serial_mesh(3, 4, 1, 0, 3);
  constexpr std::array<Vector3, 4> vertices{Vector3{0.0, 0.0, 0.0}, Vector3{2.0, 0.0, 0.0},
                                            Vector3{0.1, 1.0, 0.0},
                                            Vector3{-0.2, 0.1, 1.3}};
  for (const auto &vertex : vertices)
  {
    serial_mesh.AddVertex(vertex.data());
  }
  constexpr std::array<int, 4> tetrahedron{0, 1, 2, 3};
  serial_mesh.AddTet(tetrahedron.data(), 1);
  serial_mesh.FinalizeTopology();
  serial_mesh.Finalize(false, false);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  TranslateQuadraticGeometry(mesh, {4.0, -3.0, 2.0});

  const auto basis =
      fem::singular::EnumerateHigherOrderEdgeGradientBases({0, 1, 2, 3}, 1, 0.5).front();
  fem::singular::DofTopology topology;
  topology.h1_dofs.resize(1);
  topology.h1_dofs[0].family = basis.family;
  topology.h1_dofs[0].order = basis.order;
  topology.elements.resize(1);
  topology.elements[0].h1.push_back({0, basis});

  fem::singular::ParallelDofNumbering numbering;
  numbering.h1.global_local_size = 1;
  numbering.h1.local_size = 1;
  numbering.h1.global_size = 1;
  numbering.h1.owned_size = 1;
  numbering.h1.owner = {0};
  numbering.h1.local_to_true = {0};

  mfem::H1_FECollection collection(1, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &collection);
  mfem::Vector combined_true_dofs(fespace.GetTrueVSize() + 1);
  combined_true_dofs = 0.0;

  fem::singular::FeatureTopology features;
  features.segments.push_back({0, {0, 1}, 0, {}});
  features.features.push_back({0, {0}, {0, 1}, 0.5, false});
  features.elements.resize(1);
  features.elements[0].edges.push_back({0, 0, 0, {0, 1, 2, 3}});
  const std::vector<fem::singular::GlobalVertexId> source_vertices{0, 1, 2, 3};
  const std::vector<fem::singular::GlobalVertexId> source_elements{7};

  fem::singular::EnrichedH1FieldEvaluator evaluator(topology, numbering, fespace);
  CHECK_THROWS_AS(evaluator.FitEdgeSlopes(features, source_vertices, source_elements),
                  std::logic_error);
  evaluator.SetFromTrueDofs(combined_true_dofs);
  const auto zero_slopes =
      evaluator.FitEdgeSlopes(features, source_vertices, source_elements);
  REQUIRE(zero_slopes.size() == 1);
  CHECK_FALSE(zero_slopes.front().valid);
  CHECK(zero_slopes.front().fitted_slope == 0.0);
  CHECK(zero_slopes.front().r_squared == 0.0);
  CHECK(zero_slopes.front().minimum_distance == 0.0);
  CHECK(zero_slopes.front().maximum_distance == 0.0);
  CHECK(zero_slopes.front().field_norm_at_minimum_distance == 0.0);
  CHECK(zero_slopes.front().field_norm_at_maximum_distance == 0.0);

  combined_true_dofs[fespace.GetTrueVSize()] = 1.0;
  evaluator.SetFromTrueDofs(combined_true_dofs);
  const auto slopes = evaluator.FitEdgeSlopes(features, source_vertices, source_elements);
  REQUIRE(slopes.size() == 1);
  const auto &slope = slopes.front();
  CHECK(slope.source_element == 7);
  CHECK(slope.feature == 0);
  CHECK(slope.segment == 0);
  CHECK(std::equal(slope.canonical_vertices.begin(), slope.canonical_vertices.end(),
                   source_vertices.begin()));
  CHECK(slope.sample_count == 9);
  CHECK(slope.valid);
  CHECK(slope.exponent == 0.5);
  CHECK(slope.expected_slope == -0.5);
  CHECK_THAT(slope.fitted_slope, WithinAbs(-0.5, 1.0e-2));
  CHECK(slope.r_squared > 0.999);
  CHECK(slope.minimum_distance > 0.0);
  CHECK(slope.maximum_distance > slope.minimum_distance);
  CHECK_THAT(slope.maximum_distance / slope.minimum_distance, WithinRel(1024.0, 1.0e-9));
  CHECK(slope.field_norm_at_minimum_distance > slope.field_norm_at_maximum_distance);
}

TEST_CASE("Combined triangular singular H1 field evaluation and sampling",
          "[singularfield][triangle][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  mfem::Mesh serial_mesh(2, 3, 1, 0, 2);
  constexpr std::array<fem::singular::Vector2, 3> vertices{
      fem::singular::Vector2{0.2, -0.1}, fem::singular::Vector2{1.7, 0.3},
      fem::singular::Vector2{-0.2, 1.4}};
  for (const auto &vertex : vertices)
  {
    serial_mesh.AddVertex(vertex.data());
  }
  serial_mesh.AddTriangle(0, 1, 2, 1);
  serial_mesh.FinalizeTopology();
  serial_mesh.Finalize(false, false);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  const fem::singular::TriangleBasis basis{
      fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, {0, 1, 2}, 1, 0.5};
  fem::singular::TriangleDofTopology topology;
  topology.h1_dofs.resize(1);
  topology.h1_dofs[0].family = basis.family;
  topology.h1_dofs[0].order = basis.order;
  topology.elements.resize(1);
  topology.elements[0].h1.push_back({0, basis});

  fem::singular::ParallelDofNumbering numbering;
  numbering.h1.global_local_size = 1;
  numbering.h1.local_size = 1;
  numbering.h1.global_size = 1;
  numbering.h1.owned_size = 1;
  numbering.h1.owner = {0};
  numbering.h1.local_to_true = {0};

  mfem::H1_FECollection collection(1, 2);
  mfem::ParFiniteElementSpace fespace(&mesh, &collection);
  mfem::FunctionCoefficient standard_coefficient([](const mfem::Vector &x)
                                                 { return 0.7 + 1.5 * x[0] - 0.4 * x[1]; });
  mfem::ParGridFunction standard_field(&fespace);
  standard_field.ProjectCoefficient(standard_coefficient);
  mfem::Vector standard_true_dofs(fespace.GetTrueVSize());
  standard_field.ParallelProject(standard_true_dofs);
  mfem::Vector combined_true_dofs(standard_true_dofs.Size() + 1);
  for (int i = 0; i < standard_true_dofs.Size(); i++)
  {
    combined_true_dofs[i] = standard_true_dofs[i];
  }
  combined_true_dofs[standard_true_dofs.Size()] = -0.37;

  fem::singular::TriangleEnrichedH1FieldEvaluator evaluator(topology, numbering, fespace);
  CHECK_THROWS_AS(evaluator.GetOwnedCoefficientDiagnostics(), std::logic_error);
  evaluator.SetFromTrueDofs(combined_true_dofs);
  const auto diagnostics = evaluator.GetOwnedCoefficientDiagnostics();
  REQUIRE(diagnostics.size() == 1);
  CHECK(diagnostics[0].true_dof == 0);
  CHECK(diagnostics[0].exponent == 0.5);
  CHECK(diagnostics[0].coefficient == -0.37);

  mfem::IntegrationPoint point;
  point.Set2(0.23, 0.31);
  const fem::singular::TriangleBarycentricPoint lambda{0.46, 0.23, 0.31};
  auto *transformation = mesh.GetElementTransformation(0);
  REQUIRE(transformation);
  double jacobian_determinant;
  const auto grad_lambda = fem::singular::GetAffineTriangleBarycentricGradients(
      *transformation, jacobian_determinant);
  const auto singular =
      fem::singular::EvaluateTriangleNodeGradient(lambda, grad_lambda, 0, 1, 0.5);
  mfem::Vector physical_point(2);
  transformation->Transform(point, physical_point);
  const auto value = evaluator.Evaluate(0, point);
  CHECK_THAT(value.potential,
             WithinAbs(0.7 + 1.5 * physical_point[0] - 0.4 * physical_point[1] -
                           0.37 * fem::singular::EvaluateTriangleNodeGradientPotential(
                                      lambda, 0, 1, 0.5),
                       2.0e-13));
  CHECK_THAT(value.gradient[0], WithinAbs(1.5 - 0.37 * singular.value[0], 2.0e-12));
  CHECK_THAT(value.gradient[1], WithinAbs(-0.4 - 0.37 * singular.value[1], 2.0e-12));

  mfem::IntegrationPoint edge_point;
  edge_point.Set2(0.35, 0.0);
  const auto edge_value = evaluator.EvaluateClosure(0, edge_point);
  double expanded_potential = 0.0;
  fem::singular::Vector2 expanded_gradient{};
  for (const auto &term : evaluator.EvaluateValueTraceExpansion(0, edge_point, {0, 1}))
  {
    expanded_potential += std::pow(edge_point.x, term.exponents.left) *
                          std::pow(1.0 - edge_point.x, term.exponents.right) *
                          term.coefficient;
  }
  for (const auto &term : evaluator.EvaluateGradientTraceExpansion(0, edge_point, {0, 1}))
  {
    const double weight = std::pow(edge_point.x, term.exponents.left) *
                          std::pow(1.0 - edge_point.x, term.exponents.right);
    expanded_gradient[0] += weight * term.coefficient[0];
    expanded_gradient[1] += weight * term.coefficient[1];
  }
  CHECK_THAT(expanded_potential, WithinAbs(edge_value.potential, 2.0e-13));
  CHECK_THAT(expanded_gradient[0], WithinAbs(edge_value.gradient[0], 2.0e-12));
  CHECK_THAT(expanded_gradient[1], WithinAbs(edge_value.gradient[1], 2.0e-12));

  mfem::L2_FECollection sample_collection(2, 2, mfem::BasisType::GaussLegendre,
                                          mfem::FiniteElement::VALUE);
  mfem::ParFiniteElementSpace sampled_potential_space(&mesh, &sample_collection);
  mfem::ParFiniteElementSpace sampled_electric_space(&mesh, &sample_collection, 2,
                                                     mfem::Ordering::byVDIM);
  mfem::ParGridFunction sampled_potential(&sampled_potential_space);
  mfem::ParGridFunction sampled_electric(&sampled_electric_space);
  evaluator.ProjectToDiscontinuousGridFunctions(sampled_potential, sampled_electric);
  const auto &sample_points = sampled_potential_space.GetFE(0)->GetNodes();
  mfem::Vector sampled_electric_value(2);
  for (int i = 0; i < sample_points.GetNPoints(); i++)
  {
    const auto &sample_point = sample_points.IntPoint(i);
    const auto exact = evaluator.Evaluate(0, sample_point);
    CHECK_THAT(sampled_potential.GetValue(0, sample_point),
               WithinAbs(exact.potential, 2.0e-13));
    sampled_electric.GetVectorValue(0, sample_point, sampled_electric_value);
    for (int d = 0; d < 2; d++)
    {
      CHECK_THAT(sampled_electric_value[d], WithinAbs(-exact.gradient[d], 2.0e-12));
    }
  }
}

TEST_CASE("Combined triangular singular ND field evaluation curl and energy",
          "[singularfield][triangle][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  mfem::Mesh serial_mesh(2, 3, 1, 0, 2);
  constexpr std::array<fem::singular::Vector2, 3> vertices{
      fem::singular::Vector2{0.2, -0.1}, fem::singular::Vector2{1.7, 0.3},
      fem::singular::Vector2{-0.2, 1.4}};
  for (const auto &vertex : vertices)
  {
    serial_mesh.AddVertex(vertex.data());
  }
  serial_mesh.AddTriangle(0, 1, 2, 1);
  serial_mesh.FinalizeTopology();
  serial_mesh.Finalize(false, false);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  const fem::singular::TriangleBasis gradient_basis{
      fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, {0, 1, 2}, 1, 0.5};
  const fem::singular::TriangleBasis rotational_basis{
      fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL, {0, 1, 2}, 1, 0.5};
  fem::singular::TriangleDofTopology topology;
  topology.nd_dofs.resize(2);
  topology.nd_dofs[0].family = gradient_basis.family;
  topology.nd_dofs[0].order = gradient_basis.order;
  topology.nd_dofs[1].family = rotational_basis.family;
  topology.nd_dofs[1].order = rotational_basis.order;
  topology.elements.resize(1);
  topology.elements[0].nd = {{0, gradient_basis}, {1, rotational_basis}};

  fem::singular::ParallelDofNumbering numbering;
  numbering.nd.global_local_size = 2;
  numbering.nd.local_size = 2;
  numbering.nd.global_size = 2;
  numbering.nd.owned_size = 2;
  numbering.nd.owner = {0, 0};
  numbering.nd.local_to_true = {0, 1};

  mfem::ND_FECollection collection(1, 2);
  mfem::ParFiniteElementSpace fespace(&mesh, &collection);
  mfem::VectorFunctionCoefficient standard_coefficient(
      2,
      [](const mfem::Vector &x, mfem::Vector &value)
      {
        value.SetSize(2);
        value[0] = 0.8 - 0.4 * x[1];
        value[1] = -0.2 + 0.4 * x[0];
      });
  mfem::ParGridFunction standard_field(&fespace);
  standard_field.ProjectCoefficient(standard_coefficient);
  mfem::Vector standard_true_dofs(fespace.GetTrueVSize());
  standard_field.ParallelProject(standard_true_dofs);
  mfem::Vector combined_true_dofs(standard_true_dofs.Size() + 2);
  for (int i = 0; i < standard_true_dofs.Size(); i++)
  {
    combined_true_dofs[i] = standard_true_dofs[i];
  }
  combined_true_dofs[standard_true_dofs.Size()] = 0.31;
  combined_true_dofs[standard_true_dofs.Size() + 1] = -0.27;

  fem::singular::TriangleEnrichedNDFieldEvaluator evaluator(topology, numbering, fespace);
  mfem::IntegrationPoint point;
  point.Set2(0.23, 0.31);
  CHECK_THROWS_AS(evaluator.Evaluate(0, point), std::logic_error);
  CHECK_THROWS_AS(evaluator.GetOwnedCoefficientDiagnostics(), std::logic_error);
  evaluator.SetFromTrueDofs(combined_true_dofs);
  const auto coefficient_diagnostics = evaluator.GetOwnedCoefficientDiagnostics();
  REQUIRE(coefficient_diagnostics.size() == 2);
  CHECK(coefficient_diagnostics[0].true_dof == 0);
  CHECK(coefficient_diagnostics[0].key.family ==
        fem::singular::HigherOrderBasisFamily::NODE_GRADIENT);
  CHECK(coefficient_diagnostics[0].exponent == 0.5);
  CHECK(coefficient_diagnostics[0].coefficient == 0.31);
  CHECK(coefficient_diagnostics[1].true_dof == 1);
  CHECK(coefficient_diagnostics[1].key.family ==
        fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL);
  CHECK(coefficient_diagnostics[1].exponent == 0.5);
  CHECK(coefficient_diagnostics[1].coefficient == -0.27);

  const fem::singular::TriangleBarycentricPoint lambda{0.46, 0.23, 0.31};
  auto *transformation = mesh.GetElementTransformation(0);
  REQUIRE(transformation);
  double jacobian_determinant;
  const auto grad_lambda = fem::singular::GetAffineTriangleBarycentricGradients(
      *transformation, jacobian_determinant);
  const auto gradient =
      fem::singular::EvaluateTriangleNodeGradient(lambda, grad_lambda, 0, 1, 0.5);
  const auto rotational =
      fem::singular::EvaluateTriangleNodeRotational(lambda, grad_lambda, 0, 1, 2, 0.5);
  mfem::Vector physical_point(2);
  transformation->Transform(point, physical_point);
  const auto value = evaluator.Evaluate(0, point);
  CHECK_THAT(value.value[0],
             WithinAbs(0.8 - 0.4 * physical_point[1] + 0.31 * gradient.value[0] -
                           0.27 * rotational.value[0],
                       3.0e-12));
  CHECK_THAT(value.value[1],
             WithinAbs(-0.2 + 0.4 * physical_point[0] + 0.31 * gradient.value[1] -
                           0.27 * rotational.value[1],
                       3.0e-12));
  CHECK_THAT(value.curl, WithinAbs(0.8 - 0.27 * rotational.curl, 3.0e-12));
  CHECK(evaluator.GetElementNodeSingularExponent(0, 0) == 0.5);
  CHECK(evaluator.GetElementNodeSingularExponent(0, 1) == 1.0);

  mfem::IntegrationPoint edge_point;
  edge_point.Set2(0.35, 0.0);
  const fem::singular::TriangleBarycentricPoint edge_lambda{0.65, 0.35, 0.0};
  const auto edge_gradient =
      fem::singular::EvaluateTriangleNodeGradient(edge_lambda, grad_lambda, 0, 1, 0.5);
  const auto edge_rotational =
      fem::singular::EvaluateTriangleNodeRotational(edge_lambda, grad_lambda, 0, 1, 2, 0.5);
  transformation->Transform(edge_point, physical_point);
  const auto edge_value = evaluator.EvaluateClosure(0, edge_point);
  CHECK_THAT(edge_value.value[0],
             WithinAbs(0.8 - 0.4 * physical_point[1] + 0.31 * edge_gradient.value[0] -
                           0.27 * edge_rotational.value[0],
                       3.0e-12));
  CHECK_THAT(edge_value.value[1],
             WithinAbs(-0.2 + 0.4 * physical_point[0] + 0.31 * edge_gradient.value[1] -
                           0.27 * edge_rotational.value[1],
                       3.0e-12));
  CHECK_THAT(edge_value.curl, WithinAbs(0.8 - 0.27 * edge_rotational.curl, 3.0e-12));
  fem::singular::Vector2 expanded_edge_value{};
  for (const auto &term : evaluator.EvaluateTraceExpansion(0, edge_point, {0, 1}))
  {
    const double weight = std::pow(edge_point.x, term.exponents.left) *
                          std::pow(1.0 - edge_point.x, term.exponents.right);
    expanded_edge_value[0] += weight * term.coefficient[0];
    expanded_edge_value[1] += weight * term.coefficient[1];
  }
  CHECK_THAT(expanded_edge_value[0], WithinAbs(edge_value.value[0], 3.0e-12));
  CHECK_THAT(expanded_edge_value[1], WithinAbs(edge_value.value[1], 3.0e-12));
  CHECK_THROWS_AS(evaluator.Evaluate(0, edge_point), std::invalid_argument);

  mfem::L2_FECollection sample_collection(2, 2, mfem::BasisType::GaussLegendre,
                                          mfem::FiniteElement::VALUE);
  mfem::ParFiniteElementSpace sampled_field_space(&mesh, &sample_collection, 2,
                                                  mfem::Ordering::byVDIM);
  mfem::ParFiniteElementSpace sampled_curl_space(&mesh, &sample_collection);
  mfem::ParGridFunction sampled_field(&sampled_field_space);
  mfem::ParGridFunction sampled_curl(&sampled_curl_space);
  evaluator.ProjectToDiscontinuousGridFunctions(sampled_field, sampled_curl);
  const auto &sample_points = sampled_curl_space.GetFE(0)->GetNodes();
  mfem::Vector sampled_value(2);
  for (int i = 0; i < sample_points.GetNPoints(); i++)
  {
    const auto &sample_point = sample_points.IntPoint(i);
    const auto exact = evaluator.Evaluate(0, sample_point);
    sampled_field.GetVectorValue(0, sample_point, sampled_value);
    for (int d = 0; d < 2; d++)
    {
      CHECK_THAT(sampled_value[d], WithinAbs(exact.value[d], 3.0e-12));
    }
    CHECK_THAT(sampled_curl.GetValue(0, sample_point), WithinAbs(exact.curl, 3.0e-12));
  }

  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 4};
  const auto integrated = evaluator.IntegrateElementFieldEnergy(0, 2.4, options);
  REQUIRE(integrated.converged);
  REQUIRE(integrated.estimated_absolute_error >= 0.0);
  const double independent = fem::singular::IntegrateReferenceTriangleNodeDuffy(
      43, 0, fem::singular::H1DuffyRadialPower,
      [&](const fem::singular::TriangleBarycentricPoint &sample_lambda)
      {
        mfem::IntegrationPoint sample_point;
        sample_point.Set2(sample_lambda[1], sample_lambda[2]);
        const auto sample = evaluator.Evaluate(0, sample_point);
        return 2.4 * jacobian_determinant *
               (sample.value[0] * sample.value[0] + sample.value[1] * sample.value[1]);
      });
  CHECK(std::abs(integrated.value - independent) <=
        integrated.estimated_absolute_error + 2.0e-11);
}

TEST_CASE("Triangular BoundaryMode magnetic reconstruction follows Maxwell convention",
          "[singularfield][triangle][Serial]")
{
  const fem::singular::Vector2 et_real{0.7, -0.3};
  const fem::singular::Vector2 et_imag{-0.2, 0.5};
  const fem::singular::Vector2 grad_en_real{0.11, -0.17};
  const fem::singular::Vector2 grad_en_imag{-0.07, 0.13};
  const double curl_real = 0.41;
  const double curl_imag = -0.19;
  const std::complex<double> kn(1.3, -0.08);
  const double omega = 2.7;

  const auto field = fem::singular::ReconstructTriangleBoundaryModeMagneticField(
      et_real, et_imag, curl_real, curl_imag, grad_en_real, grad_en_imag, kn, omega);
  const std::complex<double> ex(et_real[0], et_imag[0]);
  const std::complex<double> ey(et_real[1], et_imag[1]);
  const std::complex<double> gx(grad_en_real[0], grad_en_imag[0]);
  const std::complex<double> gy(grad_en_real[1], grad_en_imag[1]);
  const std::complex<double> curl(curl_real, curl_imag);
  const std::complex<double> imaginary_unit(0.0, 1.0);
  const std::complex<double> expected_bz = curl / (imaginary_unit * omega);
  const std::complex<double> expected_bx =
      -(kn / omega) * (-ey) + gy / (imaginary_unit * omega);
  const std::complex<double> expected_by =
      -(kn / omega) * ex - gx / (imaginary_unit * omega);
  CHECK_THAT(field.normal_real, WithinAbs(expected_bz.real(), 1.0e-15));
  CHECK_THAT(field.normal_imag, WithinAbs(expected_bz.imag(), 1.0e-15));
  CHECK_THAT(field.transverse_real[0], WithinAbs(expected_bx.real(), 1.0e-15));
  CHECK_THAT(field.transverse_imag[0], WithinAbs(expected_bx.imag(), 1.0e-15));
  CHECK_THAT(field.transverse_real[1], WithinAbs(expected_by.real(), 1.0e-15));
  CHECK_THAT(field.transverse_imag[1], WithinAbs(expected_by.imag(), 1.0e-15));

  CHECK_THROWS_AS(
      fem::singular::ReconstructTriangleBoundaryModeMagneticField(
          et_real, et_imag, curl_real, curl_imag, grad_en_real, grad_en_imag, kn, 0.0),
      std::invalid_argument);
}

TEST_CASE("Triangular combined-field energy handles two singular tips",
          "[singularfield][triangle][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  mfem::Mesh serial_mesh(2, 3, 1, 0, 2);
  serial_mesh.AddVertex(0.0, 0.0);
  serial_mesh.AddVertex(1.3, 0.2);
  serial_mesh.AddVertex(-0.1, 1.1);
  serial_mesh.AddTriangle(0, 1, 2, 1);
  serial_mesh.FinalizeTopology();
  serial_mesh.Finalize(false, false);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  const fem::singular::TriangleBasis first{
      fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, {2, 0, 1}, 1, 0.5};
  const fem::singular::TriangleBasis second{
      fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, {1, 2, 0}, 1, 2.0 / 3.0};
  fem::singular::TriangleDofTopology topology;
  topology.h1_dofs.resize(2);
  topology.nd_dofs.resize(2);
  topology.h1_to_nd = {0, 1};
  for (int i = 0; i < 2; i++)
  {
    topology.h1_dofs[i].family = fem::singular::HigherOrderBasisFamily::NODE_GRADIENT;
    topology.h1_dofs[i].order = 1;
    topology.nd_dofs[i] = topology.h1_dofs[i];
  }
  topology.elements.resize(1);
  topology.elements[0].h1 = {{0, first}, {1, second}};
  topology.elements[0].nd = {{0, first}, {1, second}};

  fem::singular::ParallelDofNumbering numbering;
  numbering.h1.global_local_size = 2;
  numbering.h1.local_size = 2;
  numbering.h1.global_size = 2;
  numbering.h1.owned_size = 2;
  numbering.h1.owner = {0, 0};
  numbering.h1.local_to_true = {0, 1};

  mfem::H1_FECollection collection(1, 2);
  mfem::ParFiniteElementSpace fespace(&mesh, &collection);
  mfem::Vector combined_true_dofs(fespace.GetTrueVSize() + 2);
  combined_true_dofs = 0.0;
  combined_true_dofs[fespace.GetTrueVSize()] = 0.7;
  combined_true_dofs[fespace.GetTrueVSize() + 1] = -0.3;

  fem::singular::TriangleEnrichedH1FieldEvaluator evaluator(topology, numbering, fespace);
  evaluator.SetFromTrueDofs(combined_true_dofs);
  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 4};
  const auto integrated = evaluator.IntegrateElementGradientEnergy(0, 2.4, options);
  INFO("value = " << integrated.value
                  << ", error = " << integrated.estimated_absolute_error);
  REQUIRE(integrated.converged);
  REQUIRE(integrated.estimated_absolute_error >= 0.0);

  auto *transformation = mesh.GetElementTransformation(0);
  REQUIRE(transformation);
  double jacobian_determinant;
  const auto grad_lambda = fem::singular::GetAffineTriangleBarycentricGradients(
      *transformation, jacobian_determinant);
  const auto matrices = fem::singular::AssembleTriangleElementEnrichmentMatrices(
      topology.elements[0], grad_lambda, jacobian_determinant, options);
  mfem::Vector coefficients(2);
  coefficients[0] = 0.7;
  coefficients[1] = -0.3;
  mfem::Vector product(2);
  matrices.h1_diffusion.Mult(coefficients, product);
  const double algebraic = 2.4 * (coefficients * product);
  const double assembly_bound =
      2.4 * (std::abs(coefficients[0]) *
                 (std::abs(coefficients[0]) *
                      matrices.h1_diffusion_estimated_absolute_error(0, 0) +
                  std::abs(coefficients[1]) *
                      matrices.h1_diffusion_estimated_absolute_error(0, 1)) +
             std::abs(coefficients[1]) *
                 (std::abs(coefficients[0]) *
                      matrices.h1_diffusion_estimated_absolute_error(1, 0) +
                  std::abs(coefficients[1]) *
                      matrices.h1_diffusion_estimated_absolute_error(1, 1)));
  CHECK(std::abs(integrated.value - algebraic) <=
        integrated.estimated_absolute_error + assembly_bound + 1.0e-12);
}

TEST_CASE("Triangular tip-ray fit recovers the Meixner exponent",
          "[singularfield][triangle][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  mfem::Mesh serial_mesh(2, 3, 1, 0, 2);
  serial_mesh.AddVertex(0.0, 0.0);
  serial_mesh.AddVertex(2.0, 0.0);
  serial_mesh.AddVertex(0.1, 1.3);
  serial_mesh.AddTriangle(0, 1, 2, 1);
  serial_mesh.FinalizeTopology();
  serial_mesh.Finalize(false, false);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  TranslateQuadraticGeometry(mesh, {3.0, -2.0, 0.0});

  const fem::singular::TriangleBasis basis{
      fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, {0, 1, 2}, 1, 0.5};
  fem::singular::TriangleDofTopology topology;
  topology.h1_dofs.resize(1);
  topology.h1_dofs[0].family = basis.family;
  topology.h1_dofs[0].order = basis.order;
  topology.nd_dofs = topology.h1_dofs;
  topology.elements.resize(1);
  topology.elements[0].h1.push_back({0, basis});
  topology.elements[0].nd.push_back({0, basis});

  fem::singular::ParallelDofNumbering numbering;
  numbering.h1.global_local_size = 1;
  numbering.h1.local_size = 1;
  numbering.h1.global_size = 1;
  numbering.h1.owned_size = 1;
  numbering.h1.owner = {0};
  numbering.h1.local_to_true = {0};
  numbering.nd = numbering.h1;

  mfem::H1_FECollection collection(1, 2);
  mfem::ParFiniteElementSpace fespace(&mesh, &collection);
  mfem::Vector combined_true_dofs(fespace.GetTrueVSize() + 1);
  combined_true_dofs = 0.0;

  fem::singular::TriangleFeatureTopology features;
  features.vertices.push_back({0, 0, {0}, 0.5});
  features.selected_segments.push_back({0, 0, {0, 1}, 7});
  features.elements.resize(1);
  features.elements[0].nodes.push_back({0, 0, {0, 1, 2}});
  const std::vector<fem::singular::GlobalVertexId> source_vertices{0, 1, 2};
  const std::vector<fem::singular::GlobalVertexId> source_elements{9};

  fem::singular::TriangleEnrichedH1FieldEvaluator evaluator(topology, numbering, fespace);
  CHECK_THROWS_AS(evaluator.FitTipSlopes(features, source_vertices, source_elements),
                  std::logic_error);
  evaluator.SetFromTrueDofs(combined_true_dofs);
  const auto zero_slopes =
      evaluator.FitTipSlopes(features, source_vertices, source_elements);
  REQUIRE(zero_slopes.size() == 1);
  CHECK_FALSE(zero_slopes.front().valid);

  combined_true_dofs[fespace.GetTrueVSize()] = 1.0;
  evaluator.SetFromTrueDofs(combined_true_dofs);
  const auto slopes = evaluator.FitTipSlopes(features, source_vertices, source_elements);
  REQUIRE(slopes.size() == 1);
  const auto &slope = slopes.front();
  CHECK(slope.source_element == 9);
  CHECK(slope.feature == 0);
  CHECK(slope.selected_segment == 0);
  CHECK(slope.canonical_vertices == std::array<fem::singular::GlobalVertexId, 3>{0, 1, 2});
  CHECK(slope.valid);
  CHECK(slope.exponent == 0.5);
  CHECK(slope.expected_slope == -0.5);
  CHECK_THAT(slope.fitted_slope, WithinAbs(-0.5, 1.0e-2));
  CHECK(slope.r_squared > 0.999);
  CHECK_THAT(slope.maximum_distance / slope.minimum_distance, WithinRel(1024.0, 1.0e-9));
  CHECK(slope.field_norm_at_minimum_distance > slope.field_norm_at_maximum_distance);

  mfem::ND_FECollection nd_collection(1, 2);
  mfem::ParFiniteElementSpace nd_fespace(&mesh, &nd_collection);
  mfem::Vector real_true_dofs(nd_fespace.GetTrueVSize() + 1);
  mfem::Vector imaginary_true_dofs(nd_fespace.GetTrueVSize() + 1);
  real_true_dofs = 0.0;
  imaginary_true_dofs = 0.0;
  real_true_dofs[nd_fespace.GetTrueVSize()] = 1.0;
  imaginary_true_dofs[nd_fespace.GetTrueVSize()] = -0.4;
  fem::singular::TriangleEnrichedNDFieldEvaluator real_evaluator(topology, numbering,
                                                                 nd_fespace);
  fem::singular::TriangleEnrichedNDFieldEvaluator imaginary_evaluator(topology, numbering,
                                                                      nd_fespace);
  CHECK_THROWS_AS(real_evaluator.FitComplexTipSlopes(imaginary_evaluator, features,
                                                     source_vertices, source_elements),
                  std::logic_error);
  real_evaluator.SetFromTrueDofs(real_true_dofs);
  imaginary_evaluator.SetFromTrueDofs(imaginary_true_dofs);
  const auto complex_slopes = real_evaluator.FitComplexTipSlopes(
      imaginary_evaluator, features, source_vertices, source_elements);
  REQUIRE(complex_slopes.size() == 1);
  CHECK(complex_slopes.front().valid);
  CHECK_THAT(complex_slopes.front().fitted_slope, WithinAbs(-0.5, 1.0e-2));
  CHECK(complex_slopes.front().r_squared > 0.999);
}

}  // namespace palace
