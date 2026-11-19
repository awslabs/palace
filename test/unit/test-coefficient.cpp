// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <fstream>
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
#include "fixtures.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

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

TEST_CASE("SA interface entry on an interior metal sheet in vacuum is warned at set-up",
          "[coefficient][Serial][Parallel]")
{
  // The USER's warning (decisions 2197 (4) / 2204): the SA coefficient averages the two
  // sides of an interior face, so on a metal sheet with vacuum on both sides (uncracked) a
  // substrate-air entry reads near-zero raw energy. 8 x 4 squares of [0, 4] x [0, 1]:
  // vacuum for x < 2 (attribute 1), a substrate (epsilon 4) beyond (attribute 2); interior
  // boundary faces at x = 1 (vacuum on both sides; a PEC sheet, attribute 10), at x = 2
  // (vacuum / substrate; a dielectric face, attribute 11) and at x = 3 (substrate on both
  // sides; a PEC sheet, attribute 12).
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian2D(8, 4, mfem::Element::QUADRILATERAL, false, 4.0, 1.0);
  for (int e = 0; e < serial_mesh.GetNE(); e++)
  {
    mfem::Vector center(2);
    serial_mesh.GetElementCenter(e, center);
    serial_mesh.SetAttribute(e, center(0) < 2.0 ? 1 : 2);
  }
  for (int f = 0; f < serial_mesh.GetNumFaces(); f++)
  {
    int e1, e2;
    serial_mesh.GetFaceElements(f, &e1, &e2);
    if (e1 < 0 || e2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    serial_mesh.GetFaceVertices(f, vertices);
    const double x0 = serial_mesh.GetVertex(vertices[0])[0];
    const double x1 = serial_mesh.GetVertex(vertices[1])[0];
    for (const int k : {1, 2, 3})
    {
      if (std::abs(x0 - k) < 1.0e-12 && std::abs(x1 - k) < 1.0e-12)
      {
        serial_mesh.AddBdrElement(serial_mesh.GetFace(f)->Duplicate(&serial_mesh));
        serial_mesh.SetBdrAttribute(serial_mesh.GetNBE() - 1, 9 + k);
      }
    }
  }
  serial_mesh.FinalizeTopology();
  serial_mesh.Finalize();
  serial_mesh.SetAttributes();
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh);
  Mesh mesh(std::move(par_mesh));
  config::MaterialData vacuum, substrate;
  vacuum.attributes = {1};
  substrate.attributes = {2};
  substrate.epsilon_r.s = {4.0, 4.0, 4.0};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, substrate}, periodic, ProblemType::ELECTROSTATIC, mesh);

  // The boundaries that model a metal sheet: not lumped ports.
  config::BoundaryData boundaries;
  boundaries.pec.attributes = {12, 10};
  config::LumpedPortData port;
  port.elements.emplace_back();
  port.elements.back().attributes = {30};
  boundaries.lumpedport.emplace(1, port);
  config::TerminalData terminal;
  terminal.attributes = {20};
  boundaries.terminal.emplace(1, terminal);
  const auto metal = boundaries.GetMetalSheetAttributes();
  CHECK(metal == std::vector<int>{10, 12, 20});

  auto Warned = [&](InterfaceDielectric type, const std::vector<int> &attributes)
  {
    return GetVacuumMetalSheetSubstrateAirAttributes(mesh.Get(), mat_op, type, attributes,
                                                     metal);
  };
  // Fires for the SA entry on the sheet in vacuum, also within a longer list.
  CHECK(Warned(InterfaceDielectric::SA, {10}) == std::vector<int>{10});
  CHECK(Warned(InterfaceDielectric::SA, {12, 11, 10}) == std::vector<int>{10});
  // Not for an SA entry on a dielectric face, nor on a sheet inside the substrate.
  CHECK(Warned(InterfaceDielectric::SA, {11}).empty());
  CHECK(Warned(InterfaceDielectric::SA, {12}).empty());
  // Not for a metal-air (or any other) entry on the sheet in vacuum.
  CHECK(Warned(InterfaceDielectric::MA, {10}).empty());
  CHECK(Warned(InterfaceDielectric::MS, {10}).empty());
  CHECK(Warned(InterfaceDielectric::DEFAULT, {10}).empty());

  // The one line names the entry and the attributes, what is read there and what to do.
  // On a perfect conductor (never cracked) it names no "Crack" setting.
  const auto line = VacuumMetalSheetSubstrateAirWarning(3, {10}, boundaries);
  CHECK(line.find("postprocessing index 3") != std::string::npos);
  CHECK(line.find("interior (unseparated) faces") != std::string::npos);
  CHECK(line.find("boundary attributes 10") != std::string::npos);
  CHECK(line.find("near-zero raw energy;") != std::string::npos);
  CHECK(line.find("uncracked") == std::string::npos);
  CHECK(line.find("metal-air (\"MA\") entry!") != std::string::npos);
  CHECK(line.find("Crack") == std::string::npos);
  CHECK(line.find('\n') == std::string::npos);

  // An Impedance (Ls) sheet has the "Crack" setting: the line names it, with the listed
  // attributes of that boundary only.
  config::ImpedanceData impedance;
  impedance.Ls = 1.0e-12;
  impedance.attributes = {14};
  boundaries.impedance.push_back(impedance);
  const auto ls_line = VacuumMetalSheetSubstrateAirWarning(3, {10, 14}, boundaries);
  CHECK(ls_line.find("interior (unseparated) faces") != std::string::npos);
  CHECK(ls_line.find("boundary attributes 10, 14") != std::string::npos);
  CHECK(ls_line.find("metal-air (\"MA\") entry (or separate the faces by the \"Crack\" "
                     "setting of the Impedance, RationalImpedance or Conductivity "
                     "boundary of attributes 14)!") != std::string::npos);
  CHECK(ls_line.find('\n') == std::string::npos);
  // Only the listed attributes decide: the PEC entry alone names no "Crack".
  CHECK(VacuumMetalSheetSubstrateAirWarning(3, {10}, boundaries) == line);
}

TEST_CASE("SA warning on a cracked Ls sheet lists only its unseparated faces",
          "[coefficient][Serial][Parallel]")
{
  // Review r2198 m-2: a crack can leave faces of a cracked sheet interior. 8 x 4 squares
  // (triangles) of [0, 4] x [0, 1] in vacuum, with interior boundary edges for
  // y in [0.25, 0.75] at x = 1 (attribute 9, an Ls sheet with "Crack": both of its edges
  // separated), at x = 2 (attribute 12, an Ls sheet without "Crack") and at x = 3 (a PEC
  // edge, attribute 10, below an Ls edge with "Crack", attribute 11: the crack has no
  // interior vertex there, so the Ls edge between the junction and the crack tip stays
  // interior). The mesh is cracked by mesh::ReadMesh, as in a run.
  test::SharedTempDir temp;
  const auto mesh_path = temp.temp_dir / "sa-cracked-ls-sheets.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    mfem::Mesh serial_mesh =
        mfem::Mesh::MakeCartesian2D(8, 4, mfem::Element::TRIANGLE, false, 4.0, 1.0);
    for (int f = 0; f < serial_mesh.GetNumFaces(); f++)
    {
      int e1, e2;
      serial_mesh.GetFaceElements(f, &e1, &e2);
      if (e1 < 0 || e2 < 0)
      {
        continue;
      }
      mfem::Array<int> vertices;
      serial_mesh.GetFaceVertices(f, vertices);
      const double *v0 = serial_mesh.GetVertex(vertices[0]);
      const double *v1 = serial_mesh.GetVertex(vertices[1]);
      const double x = v0[0], y = 0.5 * (v0[1] + v1[1]);
      if (std::abs(v1[0] - x) > 1.0e-12 || y < 0.25 || y > 0.75)
      {
        continue;
      }
      int attr = 0;
      if (std::abs(x - 1.0) < 1.0e-12)
      {
        attr = 9;
      }
      else if (std::abs(x - 2.0) < 1.0e-12)
      {
        attr = 12;
      }
      else if (std::abs(x - 3.0) < 1.0e-12)
      {
        attr = (y < 0.5) ? 10 : 11;
      }
      if (attr > 0)
      {
        serial_mesh.AddBdrElement(serial_mesh.GetFace(f)->Duplicate(&serial_mesh));
        serial_mesh.SetBdrAttribute(serial_mesh.GetNBE() - 1, attr);
      }
    }
    serial_mesh.FinalizeTopology();
    serial_mesh.Finalize();
    serial_mesh.SetAttributes();
    std::ofstream output(mesh_path);
    serial_mesh.Print(output);
  }
  Mpi::Barrier(Mpi::World());
  const config::json config = {
      {"Problem",
       {{"Type", "Eigenmode"}, {"Verbose", 0}, {"Output", temp.temp_dir.string()}}},
      {"Model", {{"Mesh", mesh_path.string()}}},
      {"Domains", {{"Materials", config::json::array({{{"Attributes", {1}}}})}}},
      {"Boundaries",
       {{"PEC", {{"Attributes", {10}}}},
        {"Impedance",
         config::json::array({{{"Attributes", {9, 11}}, {"Ls", 1.0e-12}, {"Crack", true}},
                              {{"Attributes", {12}}, {"Ls", 1.0e-12}}})}}},
      {"Solver", {{"Eigenmode", {{"Target", 1.0}}}}}};
  IoData iodata(config, false);
  Mesh mesh(mesh::ReadMesh(iodata, Mpi::World()));
  REQUIRE(iodata.boundaries.cracked_attributes.contains(9));
  REQUIRE(iodata.boundaries.cracked_attributes.contains(11));
  REQUIRE(!iodata.boundaries.cracked_attributes.contains(10));
  REQUIRE(!iodata.boundaries.cracked_attributes.contains(12));
  config::MaterialData vacuum;
  vacuum.attributes = {1};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum}, periodic, ProblemType::ELECTROSTATIC, mesh);

  const auto &boundaries = iodata.boundaries;
  const auto metal = boundaries.GetMetalSheetAttributes();
  auto Warned = [&](const std::vector<int> &attributes)
  {
    return GetVacuumMetalSheetSubstrateAirAttributes(
        mesh.Get(), mat_op, InterfaceDielectric::SA, attributes, metal);
  };
  // The fully separated Ls sheet is not listed; the uncracked PEC and Ls sheets are, and
  // so is the cracked Ls sheet with an unseparated face.
  CHECK(Warned({9}).empty());
  CHECK(Warned({11}) == std::vector<int>{11});
  CHECK(Warned({9, 10, 11, 12}) == std::vector<int>{10, 11, 12});

  // The texts: the "Crack" setting is named for the Ls sheets only.
  const auto pec_line = VacuumMetalSheetSubstrateAirWarning(2, Warned({10}), boundaries);
  CHECK(pec_line.find("interior (unseparated) faces") != std::string::npos);
  CHECK(pec_line.find("boundary attributes 10)") != std::string::npos);
  CHECK(pec_line.find("Crack") == std::string::npos);
  const auto ls_line = VacuumMetalSheetSubstrateAirWarning(2, Warned({11}), boundaries);
  CHECK(ls_line.find("interior (unseparated) faces") != std::string::npos);
  CHECK(ls_line.find("boundary attributes 11)") != std::string::npos);
  CHECK(ls_line.find("\"Crack\" setting of the Impedance, RationalImpedance or "
                     "Conductivity boundary of attributes 11)!") != std::string::npos);
  const auto all_line =
      VacuumMetalSheetSubstrateAirWarning(2, Warned({9, 10, 11, 12}), boundaries);
  CHECK(all_line.find("boundary attributes 10, 11, 12)") != std::string::npos);
  CHECK(all_line.find("boundary of attributes 11, 12)!") != std::string::npos);
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
