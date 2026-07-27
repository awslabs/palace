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
  CHECK_THAT(slope.maximum_distance / slope.minimum_distance, WithinRel(1024.0, 1.0e-12));
  CHECK(slope.field_norm_at_minimum_distance > slope.field_norm_at_maximum_distance);
}

}  // namespace palace
