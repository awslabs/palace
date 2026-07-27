// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/mesh.hpp"
#include "fem/singularfeatures.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

using namespace Catch::Matchers;

namespace palace
{

namespace
{

mfem::Mesh InternalSheetMesh()
{
  mfem::Mesh mesh(3, 6, 4, 10, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(1.0, 1.0, 0.0);
  mesh.AddVertex(0.5, 0.5, 1.0);
  mesh.AddVertex(0.5, 0.5, -1.0);

  mesh.AddTet(0, 1, 2, 4, 1);
  mesh.AddTet(0, 2, 1, 5, 1);
  mesh.AddTet(1, 3, 2, 4, 1);
  mesh.AddTet(1, 2, 3, 5, 1);

  // Conforming internal PEC sheet.
  mesh.AddBdrTriangle(0, 1, 2, 7);
  mesh.AddBdrTriangle(1, 3, 2, 8);

  // Exterior grounded boundary.
  mesh.AddBdrTriangle(0, 1, 4, 1);
  mesh.AddBdrTriangle(2, 0, 4, 1);
  mesh.AddBdrTriangle(1, 3, 4, 1);
  mesh.AddBdrTriangle(3, 2, 4, 1);
  mesh.AddBdrTriangle(1, 0, 5, 1);
  mesh.AddBdrTriangle(0, 2, 5, 1);
  mesh.AddBdrTriangle(3, 1, 5, 1);
  mesh.AddBdrTriangle(2, 3, 5, 1);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh InternalLineTipMesh()
{
  mfem::Mesh mesh(2, 9, 8, 9, 2);
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

  mesh.AddBdrSegment(0, 1, 1);
  mesh.AddBdrSegment(1, 2, 1);
  mesh.AddBdrSegment(2, 5, 1);
  mesh.AddBdrSegment(5, 8, 1);
  mesh.AddBdrSegment(8, 7, 1);
  mesh.AddBdrSegment(7, 6, 1);
  mesh.AddBdrSegment(6, 3, 1);
  mesh.AddBdrSegment(3, 0, 1);
  mesh.AddBdrSegment(3, 4, 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

void SetQuadraticGeometry(mfem::Mesh &mesh)
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
          value[1] += 0.1 * x[0] * x[1];
        }
        else
        {
          value[2] += 0.1 * x[0] * x[2];
        }
      });
  mesh.GetNodes()->ProjectCoefficient(geometry);
}

IoData SingularElectrostaticData()
{
  config::json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "unused.mesh"}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}, {"Permittivity", 3.25}}}}}},
      {"Boundaries",
       {{"PEC", {{"Attributes", {1}}}},
        {"Terminal", {{{"Index", 1}, {"Attributes", {7, 8}}}}}}},
      {"Solver",
       {{"Order", 1},
        {"SingularElements",
         {{"Attributes", {7, 8}},
          {"Order", 1},
          {"QuadratureOrder", 8},
          {"AbsTol", 1.0e-3},
          {"RelTol", 1.0e-3},
          {"MaxSubdivisions", 6}}},
        {"Linear", {{"MGMaxLevels", 1}}}}}};
  return IoData(config, false);
}

IoData TriangleSingularElectrostaticData()
{
  config::json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "unused.mesh"}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}, {"Permittivity", 3.25}}}}}},
      {"Boundaries",
       {{"PEC", {{"Attributes", {1}}}},
        {"Terminal", {{{"Index", 1}, {"Attributes", {7}}}}}}},
      {"Solver",
       {{"Order", 1},
        {"SingularElements",
         {{"Attributes", {7}},
          {"Order", 1},
          {"QuadratureOrder", 8},
          {"AbsTol", 2.0e-6},
          {"RelTol", 2.0e-6},
          {"MaxSubdivisions", 6}}},
        {"Linear", {{"MGMaxLevels", 1}}}}}};
  return IoData(config, false);
}

void FillVector(Vector &x, double phase)
{
  for (int i = 0; i < x.Size(); i++)
  {
    x[i] = std::sin(phase + 0.37 * (i + 1));
  }
}

}  // namespace

void CheckSingularLaplaceOperator(const std::array<int, 4> &partition, bool curved = false)
{
  auto serial_mesh = InternalSheetMesh();
  if (curved)
  {
    SetQuadraticGeometry(serial_mesh);
  }
  const auto serial_features =
      fem::singular::ExtractSerialSheetFeatures(serial_mesh, {7, 8});
  auto parallel_mesh =
      std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh, partition.data());
  const auto source_vertex_ids = fem::singular::MapPartitionedSerialVertexIds(
      serial_mesh, *parallel_mesh, partition.data());
  std::vector<fem::singular::GlobalVertexId> source_element_ids;
  for (int element = 0; element < serial_mesh.GetNE(); element++)
  {
    if (partition[element] == Mpi::Rank(Mpi::World()))
    {
      source_element_ids.push_back(element);
    }
  }
  const auto local_features = fem::singular::DistributeSerialSheetFeatures(
      serial_features, *parallel_mesh, source_vertex_ids, source_element_ids);

  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(parallel_mesh)));
  auto iodata = SingularElectrostaticData();
  LaplaceOperator laplace(iodata, meshes, &local_features, &source_vertex_ids);
  auto K = laplace.GetStiffnessMatrix();
  const auto &G = laplace.GetGradMatrix();
  const int standard_h1_size = laplace.GetH1Space().GetTrueVSize();
  const int standard_nd_size = laplace.GetNDSpace().GetTrueVSize();
  const auto *hypre_K = dynamic_cast<const mfem::HypreParMatrix *>(K.get());
  const auto *hypre_G = dynamic_cast<const mfem::HypreParMatrix *>(&G);

  REQUIRE(K);
  REQUIRE(hypre_K);
  REQUIRE(hypre_G);
  CHECK(K->Height() == K->Width());
  CHECK(G.Width() == K->Width());
  CHECK(hypre_K->GetGlobalNumRows() > laplace.GetH1Space().Get().GlobalTrueVSize());
  CHECK(hypre_G->GetGlobalNumRows() > laplace.GetNDSpace().Get().GlobalTrueVSize());
  CHECK(laplace.GlobalTrueVSize() == hypre_K->GetGlobalNumRows());
  CHECK(K->Height() >= standard_h1_size);
  CHECK(G.Height() >= standard_nd_size);
  const auto &diagnostics = laplace.GetSingularDiagnostics();
  CHECK(diagnostics.convention_version ==
        fem::singular::ReferenceIntegral::ConventionVersion);
  CHECK(diagnostics.standard_order == 1);
  CHECK(diagnostics.singular_order == 1);
  CHECK(diagnostics.quadrature_order == 8);
  if (curved)
  {
    CHECK(diagnostics.quadrature_leaf_count > 0);
    CHECK(diagnostics.duffy_reference_table_maximum_entries == 0);
    CHECK(diagnostics.duffy_reference_cache_hits == 0);
  }
  else
  {
    CHECK(diagnostics.quadrature_leaf_count == 0);
    CHECK(diagnostics.quadrature_maximum_depth == 0);
    CHECK(diagnostics.duffy_reference_table_maximum_entries > 0);
    CHECK(diagnostics.duffy_reference_cache_hits > 0);
  }
  CHECK(diagnostics.h1_enrichment_dofs ==
        hypre_K->GetGlobalNumRows() - laplace.GetH1Space().Get().GlobalTrueVSize());
  CHECK(diagnostics.nd_enrichment_dofs ==
        hypre_G->GetGlobalNumRows() - laplace.GetNDSpace().Get().GlobalTrueVSize());
  CHECK(std::isfinite(diagnostics.standard_diagonal_spread));
  CHECK(diagnostics.standard_diagonal_spread >= 1.0);
  CHECK(std::isfinite(diagnostics.enrichment_diagonal_spread));
  CHECK(diagnostics.enrichment_diagonal_spread >= 1.0);
  CHECK(std::isfinite(diagnostics.combined_diagonal_spread));
  CHECK(diagnostics.combined_diagonal_spread >= 1.0);

  Vector x(K->Width()), y(K->Width()), Kx(K->Height()), Ky(K->Height());
  FillVector(x, 0.2);
  FillVector(y, 0.9);
  K->Mult(x, Kx);
  K->Mult(y, Ky);
  const double x_K_y = linalg::Dot(Mpi::World(), x, Ky);
  const double y_K_x = linalg::Dot(Mpi::World(), y, Kx);
  CHECK_THAT(x_K_y, WithinAbs(y_K_x, 5.0e-12 * std::max(1.0, std::abs(x_K_y))));
  CHECK(linalg::Dot(Mpi::World(), x, Kx) > 0.0);

  Vector V, RHS;
  laplace.GetExcitationVector(1, *K, V, RHS);
  REQUIRE(V.Size() == K->Width());
  REQUIRE(RHS.Size() == K->Height());
  for (int i = standard_h1_size; i < V.Size(); i++)
  {
    CHECK(V[i] == 0.0);
  }

  mfem::CGSolver cg(Mpi::World());
  cg.SetOperator(*K);
  cg.SetRelTol(1.0e-12);
  cg.SetAbsTol(0.0);
  cg.SetMaxIter(500);
  cg.SetPrintLevel(0);
  cg.Mult(RHS, V);
  CHECK(cg.GetConverged());

  Vector residual(RHS);
  K->AddMult(V, residual, -1.0);
  CHECK(linalg::Norml2(Mpi::World(), residual) <=
        1.0e-10 * linalg::Norml2(Mpi::World(), RHS));

  const auto &K_unconstrained = laplace.GetUnconstrainedStiffnessMatrix();
  REQUIRE(K_unconstrained.Height() == V.Size());
  Vector KV(V.Size());
  K_unconstrained.Mult(V, KV);
  const double capacitance = linalg::Dot(Mpi::World(), V, KV);
  CHECK(std::isfinite(capacitance));
  CHECK(capacitance > 0.0);
  mfem::Vector invalid_energy_vector(V.Size() - 1);
  CHECK_THROWS_AS(laplace.GetSingularStiffnessEnergyErrorBound(invalid_energy_vector),
                  std::invalid_argument);
  const double assembly_error = laplace.GetSingularStiffnessEnergyErrorBound(V);
  CHECK(std::isfinite(assembly_error));
  CHECK(assembly_error >= 0.0);

  auto field = laplace.GetSingularFieldEvaluator();
  mfem::IntegrationPoint center;
  center.Set3(0.25, 0.25, 0.25);
  CHECK_THROWS_AS(field->Evaluate(0, center), std::logic_error);
  mfem::Vector invalid_combined_vector(V.Size() - 1);
  CHECK_THROWS_AS(field->SetFromTrueDofs(invalid_combined_vector), std::invalid_argument);
  field->SetFromTrueDofs(V);
  const auto coefficient_diagnostics = field->GetOwnedCoefficientDiagnostics();
  HYPRE_BigInt coefficient_count = coefficient_diagnostics.size();
  Mpi::GlobalSum(1, &coefficient_count, Mpi::World());
  CHECK(coefficient_count == diagnostics.h1_enrichment_dofs);
  for (const auto &coefficient : coefficient_diagnostics)
  {
    CHECK(coefficient.true_dof >= 0);
    CHECK((coefficient.key.family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT ||
           coefficient.key.family == fem::singular::HigherOrderBasisFamily::EDGE_GRADIENT));
    CHECK(coefficient.key.order == diagnostics.singular_order);
    CHECK(std::isfinite(coefficient.exponent));
    CHECK(coefficient.exponent > 0.0);
    CHECK(coefficient.exponent < 1.0);
    CHECK(std::isfinite(coefficient.coefficient));
  }
  const auto edge_slopes =
      field->FitEdgeSlopes(local_features, source_vertex_ids, source_element_ids);
  HYPRE_BigInt edge_slope_count = edge_slopes.size();
  Mpi::GlobalSum(1, &edge_slope_count, Mpi::World());
  CHECK(edge_slope_count > 0);
  for (const auto &slope : edge_slopes)
  {
    CHECK(slope.source_element >= 0);
    CHECK(slope.sample_count == 9);
    CHECK(std::isfinite(slope.exponent));
    CHECK(slope.expected_slope == slope.exponent - 1.0);
    CHECK(std::isfinite(slope.fitted_slope));
    CHECK(std::isfinite(slope.r_squared));
    CHECK(slope.valid);
    CHECK(slope.minimum_distance > 0.0);
    CHECK(slope.maximum_distance > slope.minimum_distance);
    CHECK(slope.field_norm_at_minimum_distance > 0.0);
    CHECK(slope.field_norm_at_maximum_distance > 0.0);
  }
  mfem::IntegrationPoint feature_point;
  feature_point.Set3(0.0, 0.0, 0.0);
  CHECK_THROWS_AS(field->Evaluate(0, feature_point), std::invalid_argument);
  mfem::IntegrationPoint near_feature_point;
  near_feature_point.Set3(1.0e-10, 2.0e-10, 3.0e-10);
  const auto near_feature_sample = field->Evaluate(0, near_feature_point);
  CHECK(std::isfinite(near_feature_sample.potential));
  for (double value : near_feature_sample.gradient)
  {
    CHECK(std::isfinite(value));
  }
  double integrated_energy = 0.0;
  double estimated_integration_error = 0.0;
  const fem::singular::AdaptiveAssemblyOptions integration_options{8, 1.0e-4, 1.0e-4, 7};
  for (int element = 0; element < laplace.GetMesh().GetNE(); element++)
  {
    const auto sample = field->Evaluate(element, center);
    CHECK(std::isfinite(sample.potential));
    for (double value : sample.gradient)
    {
      CHECK(std::isfinite(value));
    }

    const auto integral =
        field->IntegrateElementGradientEnergy(element, 3.25, integration_options);
    REQUIRE(integral.converged);
    integrated_energy += integral.value;
    estimated_integration_error += integral.estimated_absolute_error;
  }
  Mpi::GlobalSum(1, &integrated_energy, Mpi::World());
  Mpi::GlobalSum(1, &estimated_integration_error, Mpi::World());
  const double floating_point_error = 256.0 * std::numeric_limits<double>::epsilon() *
                                      (std::abs(integrated_energy) + std::abs(capacitance));
  CAPTURE(integrated_energy, capacitance, estimated_integration_error, assembly_error,
          floating_point_error);
  CHECK(std::abs(integrated_energy - capacitance) <=
        estimated_integration_error + assembly_error + floating_point_error);
}

void CheckTriangleSingularLaplaceOperator(const std::array<int, 8> &partition,
                                          bool curved = false)
{
  auto serial_mesh = InternalLineTipMesh();
  if (curved)
  {
    SetQuadraticGeometry(serial_mesh);
  }
  const auto serial_features =
      fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  auto parallel_mesh =
      std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh, partition.data());
  const auto source_vertex_ids = fem::singular::MapPartitionedSerialVertexIds(
      serial_mesh, *parallel_mesh, partition.data());
  std::vector<fem::singular::GlobalVertexId> source_element_ids;
  for (int element = 0; element < serial_mesh.GetNE(); element++)
  {
    if (partition[element] == Mpi::Rank(Mpi::World()))
    {
      source_element_ids.push_back(element);
    }
  }
  const auto local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, *parallel_mesh, source_vertex_ids, source_element_ids);

  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(parallel_mesh)));
  auto iodata = TriangleSingularElectrostaticData();
  LaplaceOperator laplace(iodata, meshes, nullptr, &source_vertex_ids, &local_features);
  auto K = laplace.GetStiffnessMatrix();
  const auto &G = laplace.GetGradMatrix();
  const int standard_h1_size = laplace.GetH1Space().GetTrueVSize();
  const int standard_nd_size = laplace.GetNDSpace().GetTrueVSize();
  const auto *hypre_K = dynamic_cast<const mfem::HypreParMatrix *>(K.get());
  const auto *hypre_G = dynamic_cast<const mfem::HypreParMatrix *>(&G);

  REQUIRE(K);
  REQUIRE(hypre_K);
  REQUIRE(hypre_G);
  CHECK(laplace.HasTriangleSingularEnrichment());
  CHECK(K->Height() == K->Width());
  CHECK(G.Width() == K->Width());
  CHECK(hypre_K->GetGlobalNumRows() > laplace.GetH1Space().Get().GlobalTrueVSize());
  CHECK(hypre_G->GetGlobalNumRows() > laplace.GetNDSpace().Get().GlobalTrueVSize());
  CHECK(K->Height() >= standard_h1_size);
  CHECK(G.Height() >= standard_nd_size);

  const auto &diagnostics = laplace.GetSingularDiagnostics();
  CHECK(diagnostics.standard_order == 1);
  CHECK(diagnostics.singular_order == 1);
  CHECK(diagnostics.quadrature_leaf_count > 0);
  CHECK(diagnostics.quadrature_maximum_depth == 0);
  CHECK(diagnostics.duffy_reference_table_maximum_entries == 0);
  CHECK(diagnostics.feature_patch_count == 1);
  CHECK(diagnostics.feature_patch_minimum_enrichment_multiplicity == 1);
  CHECK(diagnostics.feature_patch_maximum_enrichment_multiplicity == 1);

  Vector x(K->Width()), y(K->Width()), Kx(K->Height()), Ky(K->Height());
  FillVector(x, 0.2);
  FillVector(y, 0.9);
  K->Mult(x, Kx);
  K->Mult(y, Ky);
  const double x_K_y = linalg::Dot(Mpi::World(), x, Ky);
  const double y_K_x = linalg::Dot(Mpi::World(), y, Kx);
  CHECK_THAT(x_K_y, WithinAbs(y_K_x, 5.0e-12 * std::max(1.0, std::abs(x_K_y))));
  CHECK(linalg::Dot(Mpi::World(), x, Kx) > 0.0);

  Vector V, RHS;
  laplace.GetExcitationVector(1, *K, V, RHS);
  REQUIRE(V.Size() == K->Width());
  REQUIRE(RHS.Size() == K->Height());
  for (int i = standard_h1_size; i < V.Size(); i++)
  {
    CHECK(V[i] == 0.0);
  }

  mfem::CGSolver cg(Mpi::World());
  cg.SetOperator(*K);
  cg.SetRelTol(1.0e-12);
  cg.SetAbsTol(0.0);
  cg.SetMaxIter(500);
  cg.SetPrintLevel(0);
  cg.Mult(RHS, V);
  CHECK(cg.GetConverged());

  Vector unconstrained_V(V.Size());
  laplace.GetUnconstrainedStiffnessMatrix().Mult(V, unconstrained_V);
  const double capacitance = linalg::Dot(Mpi::World(), V, unconstrained_V);
  CHECK(std::isfinite(capacitance));
  CHECK(capacitance > 0.0);
}

TEST_CASE("LaplaceOperator assembles and constrains singular H1 enrichment",
          "[laplaceoperator][singular][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  CheckSingularLaplaceOperator({0, 0, 0, 0});
}

TEST_CASE("LaplaceOperator singular H1 enrichment is partition conforming",
          "[laplaceoperator][singular][Parallel]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 2);
  CheckSingularLaplaceOperator({0, 0, 1, 1});
}

TEST_CASE("LaplaceOperator assembles singular H1 enrichment on curved tetrahedra",
          "[laplaceoperator][singular][curved][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  CheckSingularLaplaceOperator({0, 0, 0, 0}, true);
}

TEST_CASE("Curved tetrahedral singular H1 enrichment is partition conforming",
          "[laplaceoperator][singular][curved][Parallel]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 2);
  CheckSingularLaplaceOperator({0, 0, 1, 1}, true);
}

TEST_CASE("LaplaceOperator assembles triangular singular H1 enrichment",
          "[laplaceoperator][singular][triangle][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  CheckTriangleSingularLaplaceOperator({0, 0, 0, 0, 0, 0, 0, 0});
}

TEST_CASE("LaplaceOperator triangular singular H1 enrichment is partition conforming",
          "[laplaceoperator][singular][triangle][Parallel]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 2);
  CheckTriangleSingularLaplaceOperator({0, 0, 0, 0, 1, 1, 1, 1});
}

TEST_CASE("LaplaceOperator assembles singular H1 enrichment on curved triangles",
          "[laplaceoperator][singular][triangle][curved][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  CheckTriangleSingularLaplaceOperator({0, 0, 0, 0, 0, 0, 0, 0}, true);
}

TEST_CASE("Curved triangular singular H1 enrichment is partition conforming",
          "[laplaceoperator][singular][triangle][curved][Parallel]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 2);
  CheckTriangleSingularLaplaceOperator({0, 0, 0, 0, 1, 1, 1, 1}, true);
}

}  // namespace palace
