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
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
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

TEST_CASE("InterfaceDielectricCoefficient recovered normal flux", "[coefficient][Serial]")
{
  // Two materials separated by an internal boundary at x=1. Element 1 is substrate
  // (epsilon=4), and element 2 is vacuum.
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian2D(2, 1, mfem::Element::QUADRILATERAL, false, 2.0, 1.0);
  serial_mesh.SetAttribute(0, 1);
  serial_mesh.SetAttribute(1, 2);
  for (int f = 0; f < serial_mesh.GetNumFaces(); f++)
  {
    int e1, e2;
    serial_mesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0)
    {
      serial_mesh.AddBdrElement(serial_mesh.GetFace(f)->Duplicate(&serial_mesh));
      serial_mesh.SetBdrAttribute(serial_mesh.GetNBE() - 1, 99);
    }
  }
  serial_mesh.FinalizeTopology();
  serial_mesh.Finalize();

  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh);
  Mesh mesh(std::move(par_mesh));
  config::MaterialData substrate, vacuum;
  substrate.attributes = {1};
  substrate.epsilon_r.s = {4.0, 4.0, 4.0};
  vacuum.attributes = {2};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({substrate, vacuum}, periodic, ProblemType::ELECTROSTATIC, mesh);

  mfem::ND_FECollection nd_fec(1, mesh.Dimension());
  mfem::RT_FECollection rt_fec(0, mesh.Dimension());
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  FiniteElementSpace rt_fespace(mesh, &rt_fec);
  GridFunction E(nd_fespace, true), D(rt_fespace, true);

  mfem::VectorFunctionCoefficient Er(2,
                                     [](const mfem::Vector &, mfem::Vector &v)
                                     {
                                       v[0] = 3.0;
                                       v[1] = 4.0;
                                     });
  mfem::VectorFunctionCoefficient Ei(2,
                                     [](const mfem::Vector &, mfem::Vector &v)
                                     {
                                       v[0] = 5.0;
                                       v[1] = 6.0;
                                     });
  mfem::VectorFunctionCoefficient Dr(2,
                                     [](const mfem::Vector &, mfem::Vector &v)
                                     {
                                       v[0] = 7.0;
                                       v[1] = 8.0;
                                     });
  mfem::VectorFunctionCoefficient Di(2,
                                     [](const mfem::Vector &, mfem::Vector &v)
                                     {
                                       v[0] = 9.0;
                                       v[1] = 10.0;
                                     });
  E.Real().ProjectCoefficient(Er);
  E.Imag().ProjectCoefficient(Ei);
  D.Real().ProjectCoefficient(Dr);
  D.Imag().ProjectCoefficient(Di);
  E.Real().ExchangeFaceNbrData();
  E.Imag().ExchangeFaceNbrData();
  D.Real().ExchangeFaceNbrData();
  D.Imag().ExchangeFaceNbrData();

  mfem::ElementTransformation *T = nullptr;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (mesh.Get().GetBdrAttribute(be) == 99)
    {
      T = mesh.Get().GetBdrElementTransformation(be);
      break;
    }
  }
  REQUIRE(T);
  mfem::IntegrationPoint ip;
  ip.Set1w(0.5, 1.0);
  T->SetIntPoint(&ip);

  constexpr double thickness = 2.0;
  constexpr double interface_epsilon = 5.0;
  constexpr double raw_normal_squared = 3.0 * 3.0 + 5.0 * 5.0;
  constexpr double recovered_normal_squared = 7.0 * 7.0 + 9.0 * 9.0;
  constexpr double tangent_squared = 4.0 * 4.0 + 6.0 * 6.0;

  InterfaceDielectricCoefficient<InterfaceDielectric::MA> ma_raw(E, mat_op, thickness,
                                                                 interface_epsilon);
  InterfaceDielectricCoefficient<InterfaceDielectric::MA> ma_recovered(
      E, mat_op, thickness, interface_epsilon, nullptr, &D);
  InterfaceDielectricCoefficient<InterfaceDielectric::MA> ma_normal(
      E, mat_op, thickness, interface_epsilon, nullptr, &D,
      InterfaceDielectricComponent::NORMAL);
  InterfaceDielectricCoefficient<InterfaceDielectric::MA> ma_tangential(
      E, mat_op, thickness, interface_epsilon, nullptr, &D,
      InterfaceDielectricComponent::TANGENTIAL);
  CHECK_THAT(ma_raw.Eval(*T, ip),
             Catch::Matchers::WithinAbs(
                 0.5 * thickness / interface_epsilon * raw_normal_squared, 1.0e-12));
  CHECK_THAT(ma_recovered.Eval(*T, ip),
             Catch::Matchers::WithinAbs(
                 0.5 * thickness / interface_epsilon * recovered_normal_squared, 1.0e-12));
  CHECK_THAT(ma_normal.Eval(*T, ip),
             Catch::Matchers::WithinAbs(ma_recovered.Eval(*T, ip), 1.0e-12));
  CHECK_THAT(ma_tangential.Eval(*T, ip), Catch::Matchers::WithinAbs(0.0, 1.0e-12));

  InterfaceDielectricCoefficient<InterfaceDielectric::MS> ms_raw(E, mat_op, thickness,
                                                                 interface_epsilon);
  InterfaceDielectricCoefficient<InterfaceDielectric::MS> ms_recovered(
      E, mat_op, thickness, interface_epsilon, nullptr, &D);
  CHECK_THAT(ms_raw.Eval(*T, ip),
             Catch::Matchers::WithinAbs(0.5 * thickness / interface_epsilon * 4.0 * 4.0 *
                                            raw_normal_squared,
                                        1.0e-12));
  CHECK_THAT(ms_recovered.Eval(*T, ip),
             Catch::Matchers::WithinAbs(
                 0.5 * thickness / interface_epsilon * recovered_normal_squared, 1.0e-12));

  InterfaceDielectricCoefficient<InterfaceDielectric::SA> sa_raw(E, mat_op, thickness,
                                                                 interface_epsilon);
  InterfaceDielectricCoefficient<InterfaceDielectric::SA> sa_recovered(
      E, mat_op, thickness, interface_epsilon, nullptr, &D);
  InterfaceDielectricCoefficient<InterfaceDielectric::SA> sa_normal(
      E, mat_op, thickness, interface_epsilon, nullptr, &D,
      InterfaceDielectricComponent::NORMAL);
  InterfaceDielectricCoefficient<InterfaceDielectric::SA> sa_tangential(
      E, mat_op, thickness, interface_epsilon, nullptr, &D,
      InterfaceDielectricComponent::TANGENTIAL);
  CHECK_THAT(sa_raw.Eval(*T, ip),
             Catch::Matchers::WithinAbs(0.5 * thickness *
                                            (interface_epsilon * tangent_squared +
                                             raw_normal_squared / interface_epsilon),
                                        1.0e-12));
  CHECK_THAT(sa_recovered.Eval(*T, ip),
             Catch::Matchers::WithinAbs(0.5 * thickness *
                                            (interface_epsilon * tangent_squared +
                                             recovered_normal_squared / interface_epsilon),
                                        1.0e-12));
  CHECK_THAT(sa_normal.Eval(*T, ip) + sa_tangential.Eval(*T, ip),
             Catch::Matchers::WithinAbs(sa_recovered.Eval(*T, ip), 1.0e-12));

  // The real-valued amplitudes used by batched response-matrix assembly reproduce the
  // ordinary quadratic coefficients exactly.
  GridFunction E_batch(nd_fespace), D_batch(rt_fespace);
  E_batch.Real().ProjectCoefficient(Er);
  D_batch.Real().ProjectCoefficient(Dr);
  E_batch.Real().ExchangeFaceNbrData();
  D_batch.Real().ExchangeFaceNbrData();
  mfem::Vector field;
  auto CheckBatched = [&](auto &coefficient)
  {
    coefficient.EvalEnergyField(*T, ip, field);
    CHECK_THAT(0.5 * (field * field),
               Catch::Matchers::WithinAbs(coefficient.Eval(*T, ip), 1.0e-12));
  };
  InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT> default_batch(
      E_batch, mat_op, thickness, interface_epsilon);
  InterfaceDielectricCoefficient<InterfaceDielectric::MA> ma_batch(
      E_batch, mat_op, thickness, interface_epsilon, nullptr, &D_batch);
  InterfaceDielectricCoefficient<InterfaceDielectric::MS> ms_batch(
      E_batch, mat_op, thickness, interface_epsilon, nullptr, &D_batch);
  InterfaceDielectricCoefficient<InterfaceDielectric::SA> sa_batch(
      E_batch, mat_op, thickness, interface_epsilon, nullptr, &D_batch);
  CheckBatched(default_batch);
  CheckBatched(ma_batch);
  CheckBatched(ms_batch);
  CheckBatched(sa_batch);
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

namespace
{

// Two unit squares sharing the interior face x = 1 (boundary attribute sheet_attr, NOT
// cracked: the face is shared by both elements, as on an uncracked sheet of the new default
// crack policy), with the element attributes given (1 = substrate eps_r 4, 2 = vacuum, 3 =
// substrate eps_r 4), the electric field E = (3, 4) left and (5, 4) right (a jump of its
// normal component), and the continuous H(div) flux D = (7, 8) left and (7, 9) right (its
// normal component single-valued across the face, as for any D on the continuous RT space).
constexpr int sheet_attr = 99;

struct TwoSidedSheet
{
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<MaterialOperator> mat_op;
  mfem::ND_FECollection nd_fec;
  mfem::RT_FECollection rt_fec;
  std::unique_ptr<FiniteElementSpace> nd_fespace, rt_fespace;
  std::unique_ptr<GridFunction> E, D;
  mfem::Array<int> sum_sides_marker;
  mfem::ElementTransformation *T = nullptr;
  mfem::IntegrationPoint ip;

  TwoSidedSheet(int attr_left, int attr_right) : nd_fec(1, 2), rt_fec(0, 2)
  {
    mfem::Mesh serial_mesh =
        mfem::Mesh::MakeCartesian2D(2, 1, mfem::Element::QUADRILATERAL, false, 2.0, 1.0);
    mfem::Vector c(2);
    for (int e = 0; e < serial_mesh.GetNE(); e++)
    {
      serial_mesh.GetElementCenter(e, c);
      serial_mesh.SetAttribute(e, (c(0) < 1.0) ? attr_left : attr_right);
    }
    for (int f = 0; f < serial_mesh.GetNumFaces(); f++)
    {
      int e1, e2;
      serial_mesh.GetFaceElements(f, &e1, &e2);
      if (e1 >= 0 && e2 >= 0)
      {
        serial_mesh.AddBdrElement(serial_mesh.GetFace(f)->Duplicate(&serial_mesh));
        serial_mesh.SetBdrAttribute(serial_mesh.GetNBE() - 1, sheet_attr);
      }
    }
    serial_mesh.FinalizeTopology();
    serial_mesh.Finalize();
    serial_mesh.SetAttributes();
    mesh =
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh));
    std::vector<config::MaterialData> materials;
    for (const int attr : {attr_left, attr_right})
    {
      if (std::any_of(materials.begin(), materials.end(),
                      [&](const auto &m) { return m.attributes.front() == attr; }))
      {
        continue;
      }
      config::MaterialData material;
      material.attributes = {attr};
      if (attr != 2)
      {
        material.epsilon_r.s = {4.0, 4.0, 4.0};
      }
      materials.push_back(material);
    }
    mat_op = std::make_unique<MaterialOperator>(materials, config::PeriodicBoundaryData(),
                                                ProblemType::ELECTROSTATIC, *mesh);
    nd_fespace = std::make_unique<FiniteElementSpace>(*mesh, &nd_fec);
    rt_fespace = std::make_unique<FiniteElementSpace>(*mesh, &rt_fec);
    E = std::make_unique<GridFunction>(*nd_fespace);
    D = std::make_unique<GridFunction>(*rt_fespace);
    mfem::VectorFunctionCoefficient E_coeff(2,
                                            [](const mfem::Vector &x, mfem::Vector &v)
                                            {
                                              v(0) = (x(0) < 1.0) ? 3.0 : 5.0;
                                              v(1) = 4.0;
                                            });
    mfem::VectorFunctionCoefficient D_coeff(2,
                                            [](const mfem::Vector &x, mfem::Vector &v)
                                            {
                                              v(0) = 7.0;
                                              v(1) = (x(0) < 1.0) ? 8.0 : 9.0;
                                            });
    E->Real().ProjectCoefficient(E_coeff);
    D->Real().ProjectCoefficient(D_coeff);
    E->Real().ExchangeFaceNbrData();
    D->Real().ExchangeFaceNbrData();
    sum_sides_marker.SetSize(sheet_attr);
    sum_sides_marker = 0;
    sum_sides_marker[sheet_attr - 1] = 1;
    for (int be = 0; be < mesh->GetNBE(); be++)
    {
      if (mesh->Get().GetBdrAttribute(be) == sheet_attr)
      {
        T = mesh->Get().GetBdrElementTransformation(be);
        break;
      }
    }
    REQUIRE(T);
    ip.Set1w(0.5, 1.0);
    T->SetIntPoint(&ip);
  }

  template <InterfaceDielectric Type>
  auto
  Coefficient(bool summed, bool with_flux,
              InterfaceDielectricComponent component = InterfaceDielectricComponent::TOTAL)
  {
    return InterfaceDielectricCoefficient<Type>(*E, *mat_op, thickness, interface_epsilon,
                                                summed ? &sum_sides_marker : nullptr,
                                                with_flux ? D.get() : nullptr, component);
  }

  static constexpr double thickness = 2.0, interface_epsilon = 5.0;
};

constexpr const char *summed_flux_refusal =
    "A recovered electric flux on a two-sided interface sheet is not supported yet";
constexpr const char *summed_batched_refusal =
    "Batched interface response fields cannot represent the summed energies of both sides";

}  // namespace

// On an interior sheet face whose two sides both qualify (upstream's SumSides rule: a
// marked sheet with vacuum on both sides for MA, substrate on both sides for MS, any marked
// sheet for DEFAULT), the energies of the two sides are summed, which is the sum over the
// two copies of a cracked mesh (R4 C8). A recovered flux D on the continuous H(div) space
// has a single-valued normal trace there, so that the per-side reads would count one D_n
// twice: the coefficient refuses such a D (P3-13; our method rule keeps FluxRecovery off
// interior sheets, R4 C9: fail closed, never average). Where only one side qualifies, D is
// read on that side, and SA is never summed.
TEST_CASE("InterfaceDielectricCoefficient summed sheet with recovered flux",
          "[coefficient][Serial]")
{
  constexpr double t = TwoSidedSheet::thickness, eps = TwoSidedSheet::interface_epsilon;
  constexpr double tol = 1.0e-12;
  SECTION("Metal-air, vacuum on both sides")
  {
    TwoSidedSheet sheet(2, 2);
    auto &T = *sheet.T;
    auto summed = sheet.Coefficient<InterfaceDielectric::MA>(true, false);
    auto averaged = sheet.Coefficient<InterfaceDielectric::MA>(false, false);
    auto normal = sheet.Coefficient<InterfaceDielectric::MA>(
        true, false, InterfaceDielectricComponent::NORMAL);
    auto tangential = sheet.Coefficient<InterfaceDielectric::MA>(
        true, false, InterfaceDielectricComponent::TANGENTIAL);
    CHECK_THAT(summed.Eval(T, sheet.ip),
               Catch::Matchers::WithinAbs(0.5 * t / eps * (3.0 * 3.0 + 5.0 * 5.0), tol));
    CHECK_THAT(averaged.Eval(T, sheet.ip),
               Catch::Matchers::WithinAbs(0.5 * t / eps * 4.0 * 4.0, tol));
    CHECK_THAT(normal.Eval(T, sheet.ip),
               Catch::Matchers::WithinAbs(summed.Eval(T, sheet.ip), tol));
    CHECK_THAT(tangential.Eval(T, sheet.ip), Catch::Matchers::WithinAbs(0.0, tol));
    for (auto component :
         {InterfaceDielectricComponent::TOTAL, InterfaceDielectricComponent::NORMAL,
          InterfaceDielectricComponent::TANGENTIAL})
    {
      auto with_flux = sheet.Coefficient<InterfaceDielectric::MA>(true, true, component);
      CHECK_THROWS_WITH(with_flux.Eval(T, sheet.ip),
                        Catch::Matchers::ContainsSubstring(summed_flux_refusal));
    }
    // Without the marker (no sheet sum), D is averaged as for any two-sided interface.
    auto averaged_flux = sheet.Coefficient<InterfaceDielectric::MA>(false, true);
    CHECK_THAT(averaged_flux.Eval(T, sheet.ip),
               Catch::Matchers::WithinAbs(0.5 * t / eps * 7.0 * 7.0, tol));
  }
  SECTION("Metal-substrate, substrate on both sides")
  {
    TwoSidedSheet sheet(1, 3);
    auto &T = *sheet.T;
    auto summed = sheet.Coefficient<InterfaceDielectric::MS>(true, false);
    CHECK_THAT(
        summed.Eval(T, sheet.ip),
        Catch::Matchers::WithinAbs(0.5 * t / eps * (12.0 * 12.0 + 20.0 * 20.0), tol));
    auto with_flux = sheet.Coefficient<InterfaceDielectric::MS>(true, true);
    CHECK_THROWS_WITH(with_flux.Eval(T, sheet.ip),
                      Catch::Matchers::ContainsSubstring(summed_flux_refusal));
    // Metal-air: no vacuum side, no energy.
    auto ma = sheet.Coefficient<InterfaceDielectric::MA>(true, true);
    CHECK_THAT(ma.Eval(T, sheet.ip), Catch::Matchers::WithinAbs(0.0, tol));
  }
  SECTION("One qualifying side")
  {
    // Substrate left, vacuum right: MA reads the vacuum side, MS the substrate side, with
    // or without the marker, and D replaces the field's normal part there.
    TwoSidedSheet sheet(1, 2);
    auto &T = *sheet.T;
    for (bool summed : {false, true})
    {
      INFO("sum-sides marker " << summed);
      CHECK_THAT(
          sheet.Coefficient<InterfaceDielectric::MA>(summed, false).Eval(T, sheet.ip),
          Catch::Matchers::WithinAbs(0.5 * t / eps * 5.0 * 5.0, tol));
      CHECK_THAT(sheet.Coefficient<InterfaceDielectric::MA>(summed, true).Eval(T, sheet.ip),
                 Catch::Matchers::WithinAbs(0.5 * t / eps * 7.0 * 7.0, tol));
      CHECK_THAT(
          sheet.Coefficient<InterfaceDielectric::MS>(summed, false).Eval(T, sheet.ip),
          Catch::Matchers::WithinAbs(0.5 * t / eps * 12.0 * 12.0, tol));
      CHECK_THAT(sheet.Coefficient<InterfaceDielectric::MS>(summed, true).Eval(T, sheet.ip),
                 Catch::Matchers::WithinAbs(0.5 * t / eps * 7.0 * 7.0, tol));
    }
  }
  SECTION("Generic and substrate-air")
  {
    TwoSidedSheet sheet(2, 2);
    auto &T = *sheet.T;
    // DEFAULT sums |E|^2 of both sides on a marked sheet; SA is never summed.
    CHECK_THAT(
        sheet.Coefficient<InterfaceDielectric::DEFAULT>(true, false).Eval(T, sheet.ip),
        Catch::Matchers::WithinAbs(0.5 * t * eps * (25.0 + 41.0), tol));
    CHECK_THAT(
        sheet.Coefficient<InterfaceDielectric::DEFAULT>(false, false).Eval(T, sheet.ip),
        Catch::Matchers::WithinAbs(0.5 * t * eps * 32.0, tol));
    for (bool with_flux : {false, true})
    {
      CHECK_THAT(
          sheet.Coefficient<InterfaceDielectric::SA>(true, with_flux).Eval(T, sheet.ip),
          Catch::Matchers::WithinAbs(
              sheet.Coefficient<InterfaceDielectric::SA>(false, with_flux)
                  .Eval(T, sheet.ip),
              0.0));
    }
  }
}

// The batched amplitudes (EvalEnergyField: one normal and sdim tangential amplitudes per
// sample, for the response-matrix Gram and the streamed rows) cannot represent the summed
// energies of two sides: a sample on a summed sheet face is refused (fail closed, S-5), for
// DEFAULT on any marked sheet face and for MA / MS where both sides qualify; elsewhere the
// amplitudes reproduce the quadratic coefficient.
TEST_CASE("InterfaceDielectricCoefficient batched amplitudes refuse a summed sheet",
          "[coefficient][Serial]")
{
  mfem::Vector field;
  auto Refused = [&](auto &&coefficient, TwoSidedSheet &sheet)
  {
    // The pre-solve test of the runs that use the amplitudes (P3-16) is the same rule.
    CHECK(coefficient.IsSummedSample(*sheet.T, sheet.ip));
    CHECK_THROWS_WITH(coefficient.EvalEnergyField(*sheet.T, sheet.ip, field),
                      Catch::Matchers::ContainsSubstring(summed_batched_refusal));
  };
  auto Reproduces = [&](auto &&coefficient, TwoSidedSheet &sheet)
  {
    CHECK(!coefficient.IsSummedSample(*sheet.T, sheet.ip));
    coefficient.EvalEnergyField(*sheet.T, sheet.ip, field);
    CHECK_THAT(0.5 * (field * field),
               Catch::Matchers::WithinAbs(coefficient.Eval(*sheet.T, sheet.ip), 1.0e-12));
    CHECK(field * field > 0.0);
  };
  TwoSidedSheet vacuum(2, 2), substrate(1, 3), mixed(1, 2);
  Refused(vacuum.Coefficient<InterfaceDielectric::DEFAULT>(true, false), vacuum);
  Refused(mixed.Coefficient<InterfaceDielectric::DEFAULT>(true, false), mixed);
  Refused(vacuum.Coefficient<InterfaceDielectric::MA>(true, false), vacuum);
  Refused(vacuum.Coefficient<InterfaceDielectric::MA>(true, true), vacuum);
  Refused(substrate.Coefficient<InterfaceDielectric::MS>(true, false), substrate);
  Refused(substrate.Coefficient<InterfaceDielectric::MS>(true, true), substrate);
  Reproduces(vacuum.Coefficient<InterfaceDielectric::DEFAULT>(false, false), vacuum);
  Reproduces(mixed.Coefficient<InterfaceDielectric::MA>(true, true), mixed);
  Reproduces(mixed.Coefficient<InterfaceDielectric::MS>(true, true), mixed);
  Reproduces(vacuum.Coefficient<InterfaceDielectric::SA>(true, true), vacuum);
  Reproduces(vacuum.Coefficient<InterfaceDielectric::SA>(true, false), vacuum);
}
