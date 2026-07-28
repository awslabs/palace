// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include "fem/mesh.hpp"
#include "fem/singularfeatures.hpp"
#include "linalg/divfree.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

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
  mesh.AddBdrTriangle(0, 1, 2, 7);
  mesh.AddBdrTriangle(1, 3, 2, 7);
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

void SetQuadraticGeometry(mfem::Mesh &mesh, bool curved)
{
  const int dimension = mesh.SpaceDimension();
  mesh.SetCurvature(2, false, dimension, mfem::Ordering::byVDIM);
  if (!curved)
  {
    return;
  }
  mfem::VectorFunctionCoefficient geometry(
      dimension,
      [dimension](const mfem::Vector &x, mfem::Vector &value)
      {
        value.SetSize(dimension);
        value = x;
        if (dimension == 2)
        {
          // The selected segment has y = 0 and remains exactly straight.
          value[1] += 0.1 * x[0] * x[1];
        }
        else
        {
          // The selected sheet has z = 0 and its perimeter remains straight.
          value[2] += 0.1 * x[0] * x[2];
        }
      });
  mesh.GetNodes()->ProjectCoefficient(geometry);
}

IoData SingularSpaceData(int dimension)
{
  config::json config = {
      {"Problem",
       {{"Type", dimension == 2 ? "BoundaryMode" : "Electrostatic"},
        {"Output", "test_output"}}},
      {"Model", {{"Mesh", "unused.mesh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permittivity", 2.3}, {"Permeability", 1.7}}}}}},
      {"Boundaries", {{"PEC", {{"Attributes", {1, 7}}}}}},
      {"Solver",
       {{"Order", 1},
        {"SingularElements",
         {{"Attributes", {7}},
          {"Order", 1},
          {"QuadratureOrder", 8},
          {"AbsTol", 1.0e-3},
          {"RelTol", 1.0e-3},
          {"MaxSubdivisions", 6}}},
        {"Linear", {{"MGMaxLevels", 1}, {"Tol", 1.0e-12}, {"MaxIts", 20}}}}}};
  IoData iodata(config, false);
  iodata.problem.type = ProblemType::EIGENMODE;
  return iodata;
}

void FillVector(Vector &vector, double phase)
{
  for (int i = 0; i < vector.Size(); i++)
  {
    vector[i] = std::sin(phase + 0.37 * (i + 1));
  }
}

double RelativeNorm(const Vector &vector, const Vector &reference)
{
  return linalg::Norml2(Mpi::World(), vector) /
         std::max(1.0, linalg::Norml2(Mpi::World(), reference));
}

void CheckSymmetricPositive(const mfem::HypreParMatrix &matrix)
{
  Vector x(matrix.Width()), y(matrix.Width()), Ax(matrix.Height()), Ay(matrix.Height());
  FillVector(x, 0.2);
  FillVector(y, 0.9);
  matrix.Mult(x, Ax);
  matrix.Mult(y, Ay);
  const double x_A_y = linalg::Dot(Mpi::World(), x, Ay);
  const double y_A_x = linalg::Dot(Mpi::World(), y, Ax);
  CHECK(std::abs(x_A_y - y_A_x) <=
        2.0e-10 * std::max({1.0, std::abs(x_A_y), std::abs(y_A_x)}));
  CHECK(linalg::Dot(Mpi::World(), x, Ax) > 0.0);
}

void CheckSpaceOperator(mfem::Mesh serial_mesh, bool curved, double loss_tangent = 0.0)
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int dimension = serial_mesh.Dimension();
  SetQuadraticGeometry(serial_mesh, curved);
  const bool tetrahedral = dimension == 3;
  const auto sheet_features =
      tetrahedral ? fem::singular::ExtractSerialSheetFeatures(serial_mesh, {7})
                  : fem::singular::FeatureTopology{};
  const auto line_features =
      tetrahedral ? fem::singular::TriangleFeatureTopology{}
                  : fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});

  std::vector<fem::singular::GlobalVertexId> source_vertex_ids(serial_mesh.GetNV());
  std::vector<fem::singular::GlobalVertexId> source_element_ids(serial_mesh.GetNE());
  std::iota(source_vertex_ids.begin(), source_vertex_ids.end(), 0);
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);
  auto parallel_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh);
  const auto local_sheet_features =
      tetrahedral
          ? fem::singular::DistributeSerialSheetFeatures(
                sheet_features, *parallel_mesh, source_vertex_ids, source_element_ids)
          : fem::singular::FeatureTopology{};
  const auto local_line_features =
      tetrahedral
          ? fem::singular::TriangleFeatureTopology{}
          : fem::singular::DistributeSerialLineTipFeatures(
                line_features, *parallel_mesh, source_vertex_ids, source_element_ids);

  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(parallel_mesh)));
  auto iodata = SingularSpaceData(dimension);
  iodata.domains.materials[0].tandelta.s = {loss_tangent, loss_tangent, loss_tangent};
  SpaceOperator space_op(iodata, meshes, tetrahedral ? &local_sheet_features : nullptr,
                         tetrahedral ? nullptr : &local_line_features, &source_vertex_ids);
  REQUIRE(space_op.HasSingularEnrichment());

  auto K = space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ONE);
  auto K_zero = space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<Operator>(Operator::DIAG_ONE);
  auto M_zero = space_op.GetMassMatrix<Operator>(Operator::DIAG_ZERO);
  const auto *hypre_K = dynamic_cast<const mfem::HypreParMatrix *>(K.get());
  const auto *hypre_K_zero = dynamic_cast<const mfem::HypreParMatrix *>(K_zero.get());
  const auto *hypre_M = dynamic_cast<const mfem::HypreParMatrix *>(M.get());
  const auto *hypre_M_zero = dynamic_cast<const mfem::HypreParMatrix *>(M_zero.get());
  const auto *hypre_G =
      dynamic_cast<const mfem::HypreParMatrix *>(&space_op.GetGradMatrix());
  REQUIRE(hypre_K);
  REQUIRE(hypre_K_zero);
  REQUIRE(hypre_M);
  REQUIRE(hypre_M_zero);
  REQUIRE(hypre_G);

  CHECK(K->Height() == space_op.GetNDTrueVSize());
  CHECK(M->Height() == space_op.GetNDTrueVSize());
  CHECK(hypre_G->Height() == space_op.GetNDTrueVSize());
  CHECK(hypre_G->Width() == space_op.GetH1TrueVSize());
  CHECK(K->Height() > space_op.GetNDSpace().GetTrueVSize());
  CHECK(hypre_G->Width() > space_op.GetH1Space().GetTrueVSize());
  CHECK(space_op.GlobalTrueVSize() == hypre_K->GetGlobalNumRows());
  CheckSymmetricPositive(*hypre_K);
  CheckSymmetricPositive(*hypre_M);

  auto complex_mass =
      space_op.GetMassMatrix<ComplexOperator>(Operator::DiagonalPolicy::DIAG_ZERO);
  const auto *complex_mass_real =
      dynamic_cast<const mfem::HypreParMatrix *>(complex_mass->Real());
  const auto *complex_mass_imag =
      dynamic_cast<const mfem::HypreParMatrix *>(complex_mass->Imag());
  REQUIRE(complex_mass_real);
  CHECK((complex_mass_imag != nullptr) == (loss_tangent > 0.0));
  if (complex_mass_imag)
  {
    Vector probe(complex_mass->Width()), real_action(complex_mass->Height()),
        imaginary_action(complex_mass->Height());
    FillVector(probe, 0.61);
    complex_mass_real->Mult(probe, real_action);
    complex_mass_imag->Mult(probe, imaginary_action);
    imaginary_action.Add(loss_tangent, real_action);
    CHECK(RelativeNorm(imaginary_action, real_action) < 2.0e-11);

    constexpr std::complex<double> coefficient(1.3, -0.4);
    auto system = space_op.GetSystemMatrix(
        std::complex<double>(0.0), std::complex<double>(0.0), coefficient,
        static_cast<const ComplexOperator *>(nullptr),
        static_cast<const ComplexOperator *>(nullptr), complex_mass.get());
    const auto *system_real = dynamic_cast<const mfem::HypreParMatrix *>(system->Real());
    const auto *system_imag = dynamic_cast<const mfem::HypreParMatrix *>(system->Imag());
    REQUIRE(system_real);
    REQUIRE(system_imag);
    Vector system_real_action(system->Height()), system_imaginary_action(system->Height()),
        expected_real(real_action), expected_imag(real_action);
    system_real->Mult(probe, system_real_action);
    system_imag->Mult(probe, system_imaginary_action);
    expected_real *= coefficient.real() + coefficient.imag() * loss_tangent;
    expected_imag *= coefficient.imag() - coefficient.real() * loss_tangent;
    system_real_action -= expected_real;
    system_imaginary_action -= expected_imag;
    CHECK(RelativeNorm(system_real_action, expected_real) < 2.0e-11);
    CHECK(RelativeNorm(system_imaginary_action, expected_imag) < 2.0e-11);
  }

  Vector h1(hypre_G->Width()), gradient(hypre_G->Height()),
      curl_gradient(hypre_K_zero->Height());
  FillVector(h1, 0.4);
  linalg::SetSubVector(h1, space_op.GetCombinedH1DbcTDofList(), 0.0);
  hypre_G->Mult(h1, gradient);
  hypre_K_zero->Mult(gradient, curl_gradient);
  CHECK(RelativeNorm(curl_gradient, gradient) < 5.0e-8);

  Vector nd(K->Width()), Knd(K->Height());
  FillVector(nd, 1.1);
  K->Mult(nd, Knd);
  for (const int dof : space_op.GetCombinedNDDbcTDofList())
  {
    CHECK(Knd[dof] == nd[dof]);
  }

  // The projector must act on both standard and enrichment coordinates.
  DivFreeSolver<ComplexVector> divfree(iodata, Mpi::World(), *hypre_M_zero, *hypre_G,
                                       space_op.GetCombinedH1DbcTDofList());
  ComplexVector electric(hypre_M_zero->Width()), projected(hypre_M_zero->Width()),
      mass_electric(hypre_M_zero->Height());
  FillVector(electric.Real(), 0.25);
  FillVector(electric.Imag(), 0.75);
  linalg::SetSubVector(electric, space_op.GetCombinedNDDbcTDofList(), 0.0);
  projected = electric;
  divfree.Mult(projected);

  ComplexVector weak_divergence(hypre_G->Width());
  hypre_M_zero->Mult(projected.Real(), mass_electric.Real());
  hypre_M_zero->Mult(projected.Imag(), mass_electric.Imag());
  hypre_G->MultTranspose(mass_electric.Real(), weak_divergence.Real());
  hypre_G->MultTranspose(mass_electric.Imag(), weak_divergence.Imag());
  linalg::SetSubVector(weak_divergence, space_op.GetCombinedH1DbcTDofList(), 0.0);
  CHECK(linalg::Norml2(Mpi::World(), weak_divergence) <=
        2.0e-10 * std::max(1.0, linalg::Norml2(Mpi::World(), mass_electric)));

  ComplexVector projected_twice = projected;
  divfree.Mult(projected_twice);
  linalg::AXPY(-1.0, projected, projected_twice);
  CHECK(linalg::Norml2(Mpi::World(), projected_twice) <=
        2.0e-10 * std::max(1.0, linalg::Norml2(Mpi::World(), projected)));
}

}  // namespace

TEST_CASE("Full-wave singular SpaceOperator preserves Maxwell algebra on high-order maps",
          "[spaceoperator][singularelements][curved][Serial]")
{
  SECTION("2D quadratic affine")
  {
    CheckSpaceOperator(InternalLineTipMesh(), false);
  }
  SECTION("2D genuinely curved")
  {
    CheckSpaceOperator(InternalLineTipMesh(), true);
  }
  SECTION("3D quadratic affine")
  {
    CheckSpaceOperator(InternalSheetMesh(), false);
  }
  SECTION("3D genuinely curved")
  {
    CheckSpaceOperator(InternalSheetMesh(), true);
  }
  SECTION("2D isotropic dielectric loss")
  {
    CheckSpaceOperator(InternalLineTipMesh(), false, 0.017);
  }
  SECTION("3D isotropic dielectric loss")
  {
    CheckSpaceOperator(InternalSheetMesh(), false, 0.017);
  }
}

}  // namespace palace
