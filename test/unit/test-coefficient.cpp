// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <map>
#include <memory>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/coefficient.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

using namespace palace;

TEST_CASE("BdrFieldCoefficient interior average", "[coefficient][Serial]")
{
  // 2x1 quad mesh: two elements sharing an interior face at x=1.
  mfem::Mesh mesh =
      mfem::Mesh::MakeCartesian2D(2, 1, mfem::Element::QUADRILATERAL, false, 2.0, 1.0);

  // Add interior boundary at x=1 (the shared face between elements).
  for (int f = 0; f < mesh.GetNumFaces(); f++)
  {
    int e1, e2;
    mesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0)  // Interior face.
    {
      auto *face_elem = mesh.GetFace(f)->Duplicate(&mesh);
      mesh.AddBdrElement(face_elem);
      mesh.SetBdrAttribute(mesh.GetNBE() - 1, 99);
    }
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  auto pmesh = mfem::ParMesh(MPI_COMM_WORLD, mesh);

  // H1 space, order 1.
  mfem::H1_FECollection fec(1, 2);
  mfem::ParFiniteElementSpace fes(&pmesh, &fec);
  mfem::ParGridFunction u(&fes);

  // Set u = x. At interior face (x=1): elem0 centroid ~0.5, elem1 centroid ~1.5.
  // Proper average = 1.0. Bug (comma op) would give 0.5 * (value from elem1).
  mfem::FunctionCoefficient x_coord([](const mfem::Vector &x) { return x(0); });
  u.ProjectCoefficient(x_coord);

  BdrFieldCoefficient coeff(u);

  // Evaluate on interior boundary (attr 99).
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    if (pmesh.GetBdrAttribute(be) == 99)
    {
      auto *T = pmesh.GetBdrElementTransformation(be);
      mfem::IntegrationPoint ip;
      ip.Set1w(0.5, 1.0);  // Midpoint of the face.
      double val = coeff.Eval(*T, ip);
      // At x=1 face, both elements should give u=1.0 (since u=x and face is at x=1).
      // Average is 1.0. With comma bug: 0.5 * 1.0 = 0.5.
      CHECK_THAT(val, Catch::Matchers::WithinAbs(1.0, 1e-12));
      return;
    }
  }
  FAIL("Interior boundary not found");
}

TEST_CASE("BdrSurfaceFluxCoefficient direction orients geometric normal",
          "[coefficient][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));
  auto &pmesh = mesh.Get();

  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::MAGNETOSTATIC, mesh);

  mfem::H1_FECollection fec(1, 3);
  mfem::ParFiniteElementSpace fes(&pmesh, &fec, 3);
  mfem::ParGridFunction B(&fes);

  auto *T = pmesh.GetBdrElementTransformation(0);
  const auto &ip = mfem::Geometries.GetCenter(T->GetGeometryType());
  T->SetIntPoint(&ip);
  mfem::Vector normal(3);
  BdrGridFunctionCoefficient::GetNormal(*T, normal);

  mfem::Vector tangent(3);
  tangent = 0.0;
  tangent[std::abs(normal[0]) < 0.9 ? 0 : 1] = 1.0;
  tangent.Add(-(tangent * normal), normal);
  tangent /= tangent.Norml2();

  // Choose a field for which B ⋅ n differs from B projected onto the nonparallel reference
  // direction. The reference should only choose the normal sign.
  mfem::Vector field(normal);
  field *= 2.0;
  field.Add(3.0, tangent);
  mfem::VectorConstantCoefficient field_coeff(field);
  B.ProjectCoefficient(field_coeff);
  B.ExchangeFaceNbrData();

  mfem::Vector reference(normal);
  reference += tangent;
  reference /= reference.Norml2();
  using FluxCoefficient = BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC>;
  FluxCoefficient positive(nullptr, &B, mat_op, false, reference,
                           FluxCoefficient::OrientationMode::DIRECTION_BASED);
  const double positive_flux = positive.Eval(*T, ip);
  CHECK_THAT(std::abs(positive_flux), Catch::Matchers::WithinAbs(2.0, 1.0e-12));
  CHECK(std::abs(positive_flux) != Catch::Approx(std::abs(field * reference)));

  reference.Neg();
  FluxCoefficient negative(nullptr, &B, mat_op, false, reference,
                           FluxCoefficient::OrientationMode::DIRECTION_BASED);
  CHECK_THAT(negative.Eval(*T, ip), Catch::Matchers::WithinAbs(-positive_flux, 1.0e-12));
}

TEST_CASE("SurfaceCurrentOperator owns element apertures and current weights",
          "[surfacecurrentoperator][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
  mfem::ParMesh mesh(Mpi::World(), *serial_mesh);

  config::json current_json = {
      {"Index", 1},
      {"Elements",
       config::json::array({{{"Attributes", {1}},
                             {"Direction", "+X"},
                             {"Aperture", {{"Attributes", {3}}, {"Direction", "+Y"}}}},
                            {{"Attributes", {2}},
                             {"Direction", "+X"},
                             {"Aperture", {{"Attributes", {4}}, {"Direction", "-Y"}}}}})}};
  config::SurfaceCurrentData current(current_json);
  SurfaceCurrentOperator current_op({{1, current}}, mesh);
  const auto &source = current_op.GetSource(1);
  REQUIRE(source.elements.size() == 2);
  for (const auto &element : source.elements)
  {
    CHECK(element.current_fraction == Catch::Approx(0.5));
    REQUIRE(element.aperture);
  }
  CHECK(source.elements[0].aperture->attributes[0] == 3);
  CHECK(source.elements[0].aperture->direction[1] == Catch::Approx(1.0));
  CHECK(source.elements[1].aperture->attributes[0] == 4);
  CHECK(source.elements[1].aperture->direction[1] == Catch::Approx(-1.0));

  config::SurfaceCurrentData invalid(
      config::json{{"Index", 2},
                   {"Attributes", {1}},
                   {"Direction", "+X"},
                   {"Aperture", {{"Attributes", {99}}, {"Direction", "+Z"}}}});
  CHECK_THROWS(SurfaceCurrentOperator({{2, invalid}}, mesh));
}
