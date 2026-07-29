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

IoData SingularSpaceData(int dimension, int order = 1, int mg_max_levels = 1)
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
  iodata.solver.order = order;
  iodata.solver.linear.mg_max_levels = mg_max_levels;
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

void CheckSpaceOperator(mfem::Mesh serial_mesh, bool curved, double loss_tangent = 0.0,
                        int order = 1, int mg_max_levels = 1)
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
  auto iodata = SingularSpaceData(dimension, order, mg_max_levels);
  iodata.domains.materials[0].tandelta.s = {loss_tangent, loss_tangent, loss_tangent};
  SpaceOperator space_op(iodata, meshes, tetrahedral ? &local_sheet_features : nullptr,
                         tetrahedral ? nullptr : &local_line_features, &source_vertex_ids);
  REQUIRE(space_op.HasSingularEnrichment());

  auto K = space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ONE);
  auto K_zero = space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<Operator>(Operator::DIAG_ONE);
  auto M_zero = space_op.GetMassMatrix<Operator>(Operator::DIAG_ZERO);
  auto M_bulk = space_op.GetBulkMassMatrix(Operator::DIAG_ZERO);
  const auto *hypre_K = dynamic_cast<const mfem::HypreParMatrix *>(K.get());
  const auto *hypre_K_zero = dynamic_cast<const mfem::HypreParMatrix *>(K_zero.get());
  const auto *hypre_M = dynamic_cast<const mfem::HypreParMatrix *>(M.get());
  const auto *hypre_M_zero = dynamic_cast<const mfem::HypreParMatrix *>(M_zero.get());
  const auto *hypre_M_bulk = dynamic_cast<const mfem::HypreParMatrix *>(M_bulk.get());
  const auto *hypre_G =
      dynamic_cast<const mfem::HypreParMatrix *>(&space_op.GetGradMatrix());
  REQUIRE(hypre_K);
  REQUIRE(hypre_K_zero);
  REQUIRE(hypre_M);
  REQUIRE(hypre_M_zero);
  REQUIRE(hypre_M_bulk);
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

  const auto nd_prolongations = space_op.GetCombinedNDProlongationOperators();
  const auto h1_prolongations = space_op.GetCombinedH1ProlongationOperators();
  const auto gradients = space_op.GetCombinedGradientOperators();
  const auto &nd_spaces = space_op.GetNDSpaces();
  const auto &h1_spaces = space_op.GetH1Spaces();
  const auto &nd_essential = space_op.GetCombinedNDDbcTDofLists();
  const auto &h1_essential = space_op.GetCombinedH1DbcTDofLists();
  REQUIRE(nd_spaces.GetNumLevels() == static_cast<std::size_t>(mg_max_levels));
  REQUIRE(h1_spaces.GetNumLevels() == nd_spaces.GetNumLevels());
  REQUIRE(nd_prolongations.size() + 1 == nd_spaces.GetNumLevels());
  REQUIRE(h1_prolongations.size() == nd_prolongations.size());
  REQUIRE(gradients.size() == nd_spaces.GetNumLevels());
  REQUIRE(nd_essential.size() == nd_spaces.GetNumLevels());
  REQUIRE(h1_essential.size() == h1_spaces.GetNumLevels());

  for (std::size_t level = 0; level < gradients.size(); level++)
  {
    const int standard_nd_size = nd_spaces.GetFESpaceAtLevel(level).GetTrueVSize();
    const int standard_h1_size = h1_spaces.GetFESpaceAtLevel(level).GetTrueVSize();
    CHECK(gradients[level]->Height() - standard_nd_size ==
          hypre_G->Height() - space_op.GetNDSpace().GetTrueVSize());
    CHECK(gradients[level]->Width() - standard_h1_size ==
          hypre_G->Width() - space_op.GetH1Space().GetTrueVSize());
    CHECK(std::is_sorted(nd_essential[level].begin(), nd_essential[level].end()));
    CHECK(std::is_sorted(h1_essential[level].begin(), h1_essential[level].end()));
    for (const int dof : nd_essential[level])
    {
      CHECK(dof >= 0);
      CHECK(dof < gradients[level]->Height());
    }
    for (const int dof : h1_essential[level])
    {
      CHECK(dof >= 0);
      CHECK(dof < gradients[level]->Width());
    }
  }

  for (std::size_t level = 0; level < nd_prolongations.size(); level++)
  {
    const auto *nd_prolongation =
        dynamic_cast<const mfem::HypreParMatrix *>(nd_prolongations[level]);
    const auto *h1_prolongation =
        dynamic_cast<const mfem::HypreParMatrix *>(h1_prolongations[level]);
    REQUIRE(nd_prolongation);
    REQUIRE(h1_prolongation);
    CHECK(nd_prolongation->Width() == gradients[level]->Height());
    CHECK(nd_prolongation->Height() == gradients[level + 1]->Height());
    CHECK(h1_prolongation->Width() == gradients[level]->Width());
    CHECK(h1_prolongation->Height() == gradients[level + 1]->Width());

    const int coarse_standard_nd = nd_spaces.GetFESpaceAtLevel(level).GetTrueVSize();
    const int fine_standard_nd = nd_spaces.GetFESpaceAtLevel(level + 1).GetTrueVSize();
    Vector coarse_enrichment(nd_prolongation->Width()),
        fine_enrichment(nd_prolongation->Height());
    coarse_enrichment = 0.0;
    for (int i = coarse_standard_nd; i < coarse_enrichment.Size(); i++)
    {
      coarse_enrichment[i] = 0.13 * (i - coarse_standard_nd + 1);
    }
    nd_prolongation->Mult(coarse_enrichment, fine_enrichment);
    for (int i = 0; i < fine_standard_nd; i++)
    {
      CHECK(fine_enrichment[i] == 0.0);
    }
    for (int i = fine_standard_nd; i < fine_enrichment.Size(); i++)
    {
      CHECK(fine_enrichment[i] ==
            coarse_enrichment[coarse_standard_nd + i - fine_standard_nd]);
    }

    Vector coarse_h1(gradients[level]->Width()),
        prolonged_h1(gradients[level + 1]->Width()),
        gradient_after_prolongation(gradients[level + 1]->Height()),
        coarse_gradient(gradients[level]->Height()),
        prolonged_gradient(gradients[level + 1]->Height());
    FillVector(coarse_h1, 0.83);
    h1_prolongation->Mult(coarse_h1, prolonged_h1);
    gradients[level + 1]->Mult(prolonged_h1, gradient_after_prolongation);
    gradients[level]->Mult(coarse_h1, coarse_gradient);
    nd_prolongation->Mult(coarse_gradient, prolonged_gradient);
    gradient_after_prolongation -= prolonged_gradient;
    CHECK(RelativeNorm(gradient_after_prolongation, prolonged_gradient) < 2.0e-12);
  }

  if (mg_max_levels > 1)
  {
    auto preconditioner = space_op.GetPreconditionerMatrix<ComplexOperator>(
        std::complex<double>(1.0), std::complex<double>(0.0), std::complex<double>(-0.04),
        0.2);
    const auto *hierarchy =
        dynamic_cast<const BaseMultigridOperator<ComplexOperator> *>(preconditioner.get());
    REQUIRE(hierarchy);
    REQUIRE(hierarchy->GetNumLevels() == gradients.size());
    REQUIRE(hierarchy->GetNumAuxiliaryLevels() == gradients.size());
    auto exact_stiffness =
        space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto exact_mass = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto exact_system = space_op.GetSystemMatrix(
        std::complex<double>(1.0), std::complex<double>(0.0), std::complex<double>(-0.04),
        exact_stiffness.get(), static_cast<const ComplexOperator *>(nullptr),
        exact_mass.get());
    ComplexVector probe(exact_system->Width()), exact_action(exact_system->Height()),
        preconditioner_action(exact_system->Height());
    FillVector(probe.Real(), 0.37);
    FillVector(probe.Imag(), 0.79);
    linalg::SetSubVector(probe, nd_essential.back(), 0.0);
    exact_system->Mult(probe, exact_action);
    hierarchy->GetOperatorAtLevel(hierarchy->GetNumLevels() - 1)
        .Mult(probe, preconditioner_action);
    linalg::AXPY(-1.0, exact_action, preconditioner_action);
    CHECK(linalg::Norml2(Mpi::World(), preconditioner_action) <
          2.0e-11 * std::max(1.0, linalg::Norml2(Mpi::World(), exact_action)));
    for (std::size_t level = 0; level < gradients.size(); level++)
    {
      CHECK(hierarchy->GetOperatorAtLevel(level).Height() == gradients[level]->Height());
      CHECK(hierarchy->GetAuxiliaryOperatorAtLevel(level).Height() ==
            gradients[level]->Width());
      CHECK((hierarchy->GetOperatorAtLevel(level).Imag() != nullptr) ==
            (loss_tangent > 0.0));
      CHECK((hierarchy->GetAuxiliaryOperatorAtLevel(level).Imag() != nullptr) ==
            (loss_tangent > 0.0));
    }
    for (std::size_t level = 0; level + 1 < gradients.size(); level++)
    {
      const auto *coarse = dynamic_cast<const mfem::HypreParMatrix *>(
          hierarchy->GetOperatorAtLevel(level).Real());
      const auto *fine = dynamic_cast<const mfem::HypreParMatrix *>(
          hierarchy->GetOperatorAtLevel(level + 1).Real());
      const auto *coarse_imaginary = dynamic_cast<const mfem::HypreParMatrix *>(
          hierarchy->GetOperatorAtLevel(level).Imag());
      const auto *fine_imaginary = dynamic_cast<const mfem::HypreParMatrix *>(
          hierarchy->GetOperatorAtLevel(level + 1).Imag());
      const auto *prolongation =
          dynamic_cast<const mfem::HypreParMatrix *>(nd_prolongations[level]);
      REQUIRE(coarse);
      REQUIRE(fine);
      REQUIRE(prolongation);
      REQUIRE((coarse_imaginary != nullptr) == (loss_tangent > 0.0));
      REQUIRE((fine_imaginary != nullptr) == (loss_tangent > 0.0));

      const int standard_size = nd_spaces.GetFESpaceAtLevel(level).GetTrueVSize();
      const auto galerkin_error = [&](const mfem::HypreParMatrix &coarse_operator,
                                      const mfem::HypreParMatrix &fine_operator,
                                      bool include_standard, bool include_enrichment)
      {
        Vector coarse_x(coarse_operator.Width()), coarse_action(coarse_operator.Height()),
            fine_x(fine_operator.Width()), fine_action(fine_operator.Height()),
            restricted_action(coarse_operator.Height());
        FillVector(coarse_x, 0.41 + level);
        for (int i = 0; i < coarse_x.Size(); i++)
        {
          if ((i < standard_size && !include_standard) ||
              (i >= standard_size && !include_enrichment))
          {
            coarse_x[i] = 0.0;
          }
        }
        linalg::SetSubVector(coarse_x, nd_essential[level], 0.0);
        coarse_operator.Mult(coarse_x, coarse_action);
        prolongation->Mult(coarse_x, fine_x);
        fine_operator.Mult(fine_x, fine_action);
        prolongation->MultTranspose(fine_action, restricted_action);
        // The recursively projected level is re-eliminated after P^T A P. The
        // V-cycle likewise zeros coarse essential residuals before correction,
        // so compare the Galerkin identity only on the free subspace.
        linalg::SetSubVector(restricted_action, nd_essential[level], 0.0);
        restricted_action -= coarse_action;
        return RelativeNorm(restricted_action, coarse_action);
      };
      const double standard_error = galerkin_error(*coarse, *fine, true, false);
      const double enrichment_error = galerkin_error(*coarse, *fine, false, true);
      const double combined_error = galerkin_error(*coarse, *fine, true, true);
      const double imaginary_error =
          coarse_imaginary ? galerkin_error(*coarse_imaginary, *fine_imaginary, true, true)
                           : 0.0;
      CAPTURE(level, standard_error, enrichment_error, combined_error, imaginary_error);
      CHECK(combined_error < 2.0e-8);
      CHECK(imaginary_error < 2.0e-8);
    }
  }

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
  DivFreeSolver<ComplexVector> divfree(iodata, Mpi::World(), *hypre_M_bulk, *hypre_G,
                                       space_op.GetCombinedH1DbcTDofList());
  ComplexVector electric(hypre_M_zero->Width()), projected(hypre_M_zero->Width()),
      mass_electric(hypre_M_zero->Height());
  FillVector(electric.Real(), 0.25);
  FillVector(electric.Imag(), 0.75);
  linalg::SetSubVector(electric, space_op.GetCombinedNDDbcTDofList(), 0.0);
  projected = electric;
  divfree.Mult(projected);

  ComplexVector weak_divergence(hypre_G->Width());
  hypre_M_bulk->Mult(projected.Real(), mass_electric.Real());
  hypre_M_bulk->Mult(projected.Imag(), mass_electric.Imag());
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
  SECTION("2D p = 2 combined hierarchy")
  {
    CheckSpaceOperator(InternalLineTipMesh(), false, 0.017, 2, 2);
  }
  SECTION("3D p = 2 combined hierarchy")
  {
    CheckSpaceOperator(InternalSheetMesh(), false, 0.017, 2, 2);
  }
}

}  // namespace palace
