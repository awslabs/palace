// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <memory>
#include <set>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/edgedistance.hpp"
#include "utils/iodata.hpp"
#include "utils/metaledge.hpp"

namespace palace
{

namespace fs = std::filesystem;

using json = nlohmann::json;
using namespace Catch::Matchers;

TEST_CASE("SurfaceResponseOperator", "[surfaceresponseoperator][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  test::SharedTempDir temp;
  const auto points_path = temp.temp_dir / "basis-points.csv";
  const auto fabricated_path = temp.temp_dir / "fabricated.csv";
  const auto thin_path = temp.temp_dir / "thin.csv";
  const auto fabricated_surface_path = temp.temp_dir / "fabricated-surface.csv";
  const auto thin_surface_path = temp.temp_dir / "thin-surface.csv";
  const auto library_path = temp.temp_dir / "fabrication-process.json";
  const auto legacy_library_path = temp.temp_dir / "fabrication-process-legacy.json";
  const auto missing_layer_library_path =
      temp.temp_dir / "fabrication-process-missing-layer.json";
  const auto invalid_library_path =
      temp.temp_dir / "fabrication-process-invalid-depth.json";
  const auto library_3d_path = temp.temp_dir / "fabrication-process-3d.json";
  const auto exact_pair_library_2d_path =
      temp.temp_dir / "fabrication-process-exact-pair-2d.json";
  const auto interpolated_pair_library_2d_path =
      temp.temp_dir / "fabrication-process-interpolated-pair-2d.json";
  const auto coupled_library_3d_path =
      temp.temp_dir / "fabrication-process-coupled-3d.json";
  const auto interpolated_coupled_library_3d_path =
      temp.temp_dir / "fabrication-process-coupled-interpolated-3d.json";
  const auto coupled_fabricated_path = temp.temp_dir / "coupled-fabricated.csv";
  const auto coupled_thin_path = temp.temp_dir / "coupled-thin.csv";
  const auto coupled_fabricated_surface_path =
      temp.temp_dir / "coupled-fabricated-surface.csv";
  const auto coupled_thin_surface_path = temp.temp_dir / "coupled-thin-surface.csv";
  const auto zero_trace_path = temp.temp_dir / "zero-trace.csv";
  const auto corner_points_path = temp.temp_dir / "corner-basis-points.csv";
  const auto corner_fabricated_path = temp.temp_dir / "corner-fabricated.csv";
  const auto corner_thin_path = temp.temp_dir / "corner-thin.csv";
  const auto corner_fabricated_surface_path =
      temp.temp_dir / "corner-fabricated-surface.csv";
  const auto corner_thin_surface_path = temp.temp_dir / "corner-thin-surface.csv";
  const auto convex_library_3d_path = temp.temp_dir / "fabrication-process-convex-3d.json";
  const auto concave_library_3d_path =
      temp.temp_dir / "fabrication-process-concave-3d.json";
  const auto strip_library_3d_path =
      temp.temp_dir / "fabrication-process-strip-3d.json";
  const auto rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-3d.json";
  const auto rounded_concave_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-concave-3d.json";
  const auto interpolated_rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-interpolated-3d.json";
  if (Mpi::Root(Mpi::World()))
  {
    {
      std::ofstream output(points_path);
      output << "x,y,z\n"
             << "-0.08,-0.06,0.0\n"
             << "0.08,-0.06,0.0\n"
             << "0.08,0.06,0.0\n"
             << "-0.08,0.06,0.0\n";
    }
    const std::array<std::array<double, 4>, 4> fabricated = {
        {{{3.0e-12, 0.5e-12, 0.2e-12, 0.1e-12}},
         {{0.5e-12, 2.0e-12, 0.3e-12, 0.2e-12}},
         {{0.2e-12, 0.3e-12, 2.5e-12, 0.4e-12}},
         {{0.1e-12, 0.2e-12, 0.4e-12, 1.5e-12}}}};
    const std::array<std::array<double, 4>, 4> thin = {
        {{{1.0e-12, 0.1e-12, 0.05e-12, 0.02e-12}},
         {{0.1e-12, 0.5e-12, 0.08e-12, 0.04e-12}},
         {{0.05e-12, 0.08e-12, 0.7e-12, 0.09e-12}},
         {{0.02e-12, 0.04e-12, 0.09e-12, 0.4e-12}}}};
    auto write_domain_matrix = [](const auto &path, const auto &matrix)
    {
      std::ofstream output(path);
      output << "basis_i,basis_j,Q_ij (J)\n";
      for (std::size_t i = 0; i < matrix.size(); i++)
      {
        for (std::size_t j = i; j < matrix.size(); j++)
        {
          output << i + 1 << "," << j + 1 << "," << matrix[i][j] << "\n";
        }
      }
    };
    auto write_surface_matrix = [](const auto &path, const auto &matrix)
    {
      std::ofstream output(path);
      output << "interface,edge,basis_i,basis_j,Q_total_ij (J)\n";
      for (int edge = 1; edge <= 2; edge++)
      {
        for (std::size_t i = 0; i < matrix.size(); i++)
        {
          for (std::size_t j = i; j < matrix.size(); j++)
          {
            output << "1," << edge << "," << i + 1 << "," << j + 1 << ","
                   << 0.5 * matrix[i][j] << "\n";
          }
        }
      }
    };
    write_domain_matrix(fabricated_path, fabricated);
    write_domain_matrix(thin_path, thin);
    write_surface_matrix(fabricated_surface_path, fabricated);
    write_surface_matrix(thin_surface_path, thin);
    const json library = {{"Version", 3},
                          {"Name", "unit-test-process"},
                          {"MatchingRadius", 0.1},
                          {"Fabrication",
                           {{"InterfaceLayers",
                             {{"SA", {{"Thickness", 0.002}, {"Permittivity", 4.0}}},
                              {"MS", {{"Thickness", 0.002}, {"Permittivity", 11.47}}},
                              {"MA", {{"Thickness", 0.002}, {"Permittivity", 10.0}}}}}}},
                          {"Models",
                           {{{"Name", "isolated"},
                             {"Topology", "IsolatedEdge"},
                             {"FabricatedMatrix", fabricated_path.string()},
                             {"ThinMatrix", thin_path.string()},
                             {"FabricatedSurfaceMatrix", fabricated_surface_path.string()},
                             {"ThinSurfaceMatrix", thin_surface_path.string()},
                             {"BasisPoints", points_path.string()},
                             {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}}}}}};
    std::ofstream output(library_path);
    output << library.dump(2) << "\n";
    auto legacy_library = library;
    legacy_library["Version"] = 2;
    legacy_library.erase("Fabrication");
    std::ofstream legacy_output(legacy_library_path);
    legacy_output << legacy_library.dump(2) << "\n";
    auto missing_layer_library = library;
    missing_layer_library["Fabrication"]["InterfaceLayers"].erase("SA");
    std::ofstream missing_layer_output(missing_layer_library_path);
    missing_layer_output << missing_layer_library.dump(2) << "\n";
    auto exact_pair_library_2d = library;
    exact_pair_library_2d["Name"] = "unit-test-exact-pair-2d";
    exact_pair_library_2d["MatchingRadius"] = 0.3;
    auto exact_pair_model_2d = exact_pair_library_2d["Models"][0];
    exact_pair_model_2d["Name"] = "strip-0.5";
    exact_pair_model_2d["Topology"] = "SameConductorStrip";
    exact_pair_model_2d["Separation"] = 0.5;
    exact_pair_model_2d["SeparationTolerance"] = 1.0e-8;
    exact_pair_library_2d["Models"] = {exact_pair_model_2d};
    std::ofstream exact_pair_output_2d(exact_pair_library_2d_path);
    exact_pair_output_2d << exact_pair_library_2d.dump(2) << "\n";
    auto interpolated_pair_library_2d = exact_pair_library_2d;
    interpolated_pair_library_2d["Name"] = "unit-test-interpolated-pair-2d";
    auto lower_pair_model_2d = exact_pair_model_2d;
    lower_pair_model_2d["Name"] = "strip-0.4";
    lower_pair_model_2d["Separation"] = 0.4;
    auto upper_pair_model_2d = exact_pair_model_2d;
    upper_pair_model_2d["Name"] = "strip-0.6";
    upper_pair_model_2d["Separation"] = 0.6;
    interpolated_pair_library_2d["Models"] = {lower_pair_model_2d, upper_pair_model_2d};
    std::ofstream interpolated_pair_output_2d(interpolated_pair_library_2d_path);
    interpolated_pair_output_2d << interpolated_pair_library_2d.dump(2) << "\n";
    auto invalid_library = library;
    invalid_library["Models"][0]["CouponDepth"] = 0.0;
    std::ofstream invalid_output(invalid_library_path);
    invalid_output << invalid_library.dump(2) << "\n";
    auto library_3d = library;
    library_3d["Name"] = "unit-test-process-3d";
    library_3d["MatchingRadius"] = 2.0;
    library_3d["CouponDepth"] = 2.0;
    library_3d["Models"][0]["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}},
                                             {{"Type", "MS"}, {"Coupon", 1}},
                                             {{"Type", "MA"}, {"Coupon", 1}}};
    std::ofstream output_3d(library_3d_path);
    output_3d << library_3d.dump(2) << "\n";
    auto coupled_library_3d = library_3d;
    coupled_library_3d["Version"] = 3;
    coupled_library_3d["Name"] = "unit-test-process-coupled-3d";
    coupled_library_3d["MatchingRadius"] = 7.0;
    std::array<std::array<double, 5>, 5> coupled_fabricated{};
    std::array<std::array<double, 5>, 5> coupled_thin{};
    for (std::size_t i = 0; i < coupled_fabricated.size(); i++)
    {
      coupled_fabricated[i][i] = (i == 4 ? 3.0 : 1.0) * 1.0e-12;
      coupled_thin[i][i] = 1.0e-12;
    }
    write_domain_matrix(coupled_fabricated_path, coupled_fabricated);
    write_domain_matrix(coupled_thin_path, coupled_thin);
    write_surface_matrix(coupled_fabricated_surface_path, coupled_fabricated);
    write_surface_matrix(coupled_thin_surface_path, coupled_thin);
    {
      std::ofstream zero_trace(zero_trace_path);
      zero_trace << "x,y,z,V\n"
                 << "0.0,0.0,0.0,0.0\n"
                 << "0.5,0.5,0.0,0.0\n"
                 << "1.0,1.0,0.0,0.0\n";
    }
    auto coupled_model = coupled_library_3d["Models"][0];
    coupled_model["Name"] = "terminal-ground-gap-12um";
    coupled_model["Topology"] = "DifferentConductorGap";
    coupled_model["Separation"] = 12.0;
    coupled_model["SeparationTolerance"] = 1.0e-6;
    coupled_model["ConductorReferences"] = {{-13.0, 0.0, 0.0}, {13.0, 0.0, 0.0}};
    coupled_model["OpenContourPaths"] = {
        {{"Indices", {1, 2}}, {"StartConductor", 1}, {"EndConductor", 2}},
        {{"Indices", {4, 3}}, {"StartConductor", 1}, {"EndConductor", 2}}};
    coupled_model["FabricatedMatrix"] = coupled_fabricated_path.string();
    coupled_model["ThinMatrix"] = coupled_thin_path.string();
    coupled_model["FabricatedSurfaceMatrix"] = coupled_fabricated_surface_path.string();
    coupled_model["ThinSurfaceMatrix"] = coupled_thin_surface_path.string();
    auto interpolated_coupled_library_3d = coupled_library_3d;
    auto lower_coupled_model = coupled_model;
    lower_coupled_model["Name"] = "terminal-ground-gap-10um";
    lower_coupled_model["Separation"] = 10.0;
    lower_coupled_model["ConductorReferences"] = {{-12.0, 0.0, 0.0}, {12.0, 0.0, 0.0}};
    auto upper_coupled_model = coupled_model;
    upper_coupled_model["Name"] = "terminal-ground-gap-14um";
    upper_coupled_model["Separation"] = 14.0;
    upper_coupled_model["ConductorReferences"] = {{-14.0, 0.0, 0.0}, {14.0, 0.0, 0.0}};
    interpolated_coupled_library_3d["Name"] = "unit-test-process-coupled-interpolated-3d";
    interpolated_coupled_library_3d["Models"].push_back(lower_coupled_model);
    interpolated_coupled_library_3d["Models"].push_back(upper_coupled_model);
    std::ofstream interpolated_coupled_output_3d(interpolated_coupled_library_3d_path);
    interpolated_coupled_output_3d << interpolated_coupled_library_3d.dump(2) << "\n";
    coupled_library_3d["Models"].push_back(std::move(coupled_model));
    auto legacy_coupled_model = coupled_library_3d["Models"][0];
    legacy_coupled_model["Name"] = "legacy-terminal-ground-gap-20um";
    legacy_coupled_model["Topology"] = "DifferentConductorGap";
    legacy_coupled_model["Separation"] = 20.0;
    legacy_coupled_model["SeparationTolerance"] = 1.0e-6;
    legacy_coupled_model["Reference"] = {0.0, 0.0, 0.0};
    coupled_library_3d["Models"].push_back(std::move(legacy_coupled_model));
    std::ofstream coupled_output_3d(coupled_library_3d_path);
    coupled_output_3d << coupled_library_3d.dump(2) << "\n";

    {
      std::ofstream output(corner_points_path);
      output << "x,y,z\n";
      for (const double z : {-0.02, 0.0, 0.02})
      {
        output << "0.02,0.02," << z << "\n"
               << "0.08,0.02," << z << "\n"
               << "0.08,0.08," << z << "\n"
               << "0.02,0.08," << z << "\n";
      }
    }
    std::array<std::array<double, 12>, 12> corner_fabricated{};
    std::array<std::array<double, 12>, 12> corner_thin{};
    for (std::size_t i = 0; i < corner_fabricated.size(); i++)
    {
      for (std::size_t j = 0; j < corner_fabricated.size(); j++)
      {
        const double coupling =
            1.0 / (1.0 + std::abs(static_cast<int>(i) - static_cast<int>(j)));
        corner_fabricated[i][j] = (i == j ? 3.0 : 0.05 * coupling) * 1.0e-12;
        corner_thin[i][j] = (i == j ? 1.0 : 0.01 * coupling) * 1.0e-12;
      }
    }
    write_domain_matrix(corner_fabricated_path, corner_fabricated);
    write_domain_matrix(corner_thin_path, corner_thin);
    write_surface_matrix(corner_fabricated_surface_path, corner_fabricated);
    write_surface_matrix(corner_thin_surface_path, corner_thin);

    auto convex_library_3d = library_3d;
    convex_library_3d["Name"] = "unit-test-process-convex-3d";
    convex_library_3d["MatchingRadius"] = 0.2;
    convex_library_3d["CouponDepth"] = 0.2;
    auto corner_model = convex_library_3d["Models"][0];
    corner_model["Name"] = "convex-corner-90";
    corner_model["Topology"] = "ConvexCorner";
    corner_model["Angle"] = 90.0;
    corner_model["AngleTolerance"] = 1.0e-6;
    corner_model["FabricatedMatrix"] = corner_fabricated_path.string();
    corner_model["ThinMatrix"] = corner_thin_path.string();
    corner_model["FabricatedSurfaceMatrix"] = corner_fabricated_surface_path.string();
    corner_model["ThinSurfaceMatrix"] = corner_thin_surface_path.string();
    corner_model["BasisPoints"] = corner_points_path.string();
    corner_model["ContourGroups"] = {4, 4, 4};
    convex_library_3d["Models"].push_back(corner_model);
    std::ofstream convex_output_3d(convex_library_3d_path);
    convex_output_3d << convex_library_3d.dump(2) << "\n";

    auto concave_library_3d = convex_library_3d;
    concave_library_3d["Name"] = "unit-test-process-concave-3d";
    concave_library_3d["Models"][1]["Name"] = "concave-corner-90";
    concave_library_3d["Models"][1]["Topology"] = "ConcaveCorner";
    std::ofstream concave_output_3d(concave_library_3d_path);
    concave_output_3d << concave_library_3d.dump(2) << "\n";

    auto strip_library_3d = concave_library_3d;
    strip_library_3d["Name"] = "unit-test-process-strip-3d";
    auto strip_model = strip_library_3d["Models"][0];
    strip_model["Name"] = "same-conductor-strip-0.25";
    strip_model["Topology"] = "SameConductorStrip";
    strip_model["Separation"] = 0.25;
    strip_model["SeparationTolerance"] = 1.0e-8;
    strip_model["Reference"] = {0.0, 0.0, 0.0};
    strip_library_3d["Models"].push_back(std::move(strip_model));
    std::ofstream strip_output_3d(strip_library_3d_path);
    strip_output_3d << strip_library_3d.dump(2) << "\n";

    auto rounded_library_3d = convex_library_3d;
    rounded_library_3d["Name"] = "unit-test-process-rounded-3d";
    rounded_library_3d["Models"][1]["Name"] = "convex-corner-90-r0.125";
    rounded_library_3d["Models"][1]["CornerRadius"] = 0.125;
    rounded_library_3d["Models"][1]["CornerRadiusTolerance"] = 0.01;
    rounded_library_3d["Models"][1]["Reference"] = {0.125, 0.125, 0.0};
    rounded_library_3d["Models"][1]["ZeroTraceIndices"] = {5, 6, 7, 8};
    std::ofstream rounded_output_3d(rounded_library_3d_path);
    rounded_output_3d << rounded_library_3d.dump(2) << "\n";

    auto rounded_concave_library_3d = rounded_library_3d;
    rounded_concave_library_3d["Name"] = "unit-test-process-rounded-concave-3d";
    rounded_concave_library_3d["Models"][1]["Name"] = "concave-corner-90-r0.125";
    rounded_concave_library_3d["Models"][1]["Topology"] = "ConcaveCorner";
    rounded_concave_library_3d["Models"][1]["Reference"] = {0.0, 0.0, 0.0};
    std::ofstream rounded_concave_output_3d(rounded_concave_library_3d_path);
    rounded_concave_output_3d << rounded_concave_library_3d.dump(2) << "\n";

    auto interpolated_rounded_library_3d = rounded_library_3d;
    interpolated_rounded_library_3d["Name"] = "unit-test-process-rounded-interpolated-3d";
    interpolated_rounded_library_3d["Models"][1]["Name"] = "convex-corner-90-r0.1";
    interpolated_rounded_library_3d["Models"][1]["CornerRadius"] = 0.1;
    interpolated_rounded_library_3d["Models"][1]["CornerRadiusTolerance"] = 1.0e-4;
    interpolated_rounded_library_3d["Models"][1]["Reference"] = {0.1, 0.1, 0.0};
    auto upper_corner_model = interpolated_rounded_library_3d["Models"][1];
    upper_corner_model["Name"] = "convex-corner-90-r0.15";
    upper_corner_model["CornerRadius"] = 0.15;
    upper_corner_model["Reference"] = {0.15, 0.15, 0.0};
    interpolated_rounded_library_3d["Models"].push_back(std::move(upper_corner_model));
    std::ofstream interpolated_rounded_output_3d(interpolated_rounded_library_3d_path);
    interpolated_rounded_output_3d << interpolated_rounded_library_3d.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());

  json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", temp.temp_dir.string()}}},
      {"Model", {{"Mesh", "unused.msh"}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
      {"Boundaries",
       {{"Ground", {{"Attributes", {1, 3, 4}}}},
        {"Terminal", {{{"Index", 1}, {"Attributes", {2}}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Electrostatic",
         {{"ResponseCorrection",
           {{"Models",
             {{{"Index", 1},
               {"FabricatedMatrix", fabricated_path.string()},
               {"ThinMatrix", thin_path.string()},
               {"FabricatedSurfaceMatrix", fabricated_surface_path.string()},
               {"ThinSurfaceMatrix", thin_surface_path.string()},
               {"BasisPoints", points_path.string()},
               {"Interfaces", {{{"Target", 4}, {"Coupon", 1}}}}}}},
            {"Patches",
             {{{"Model", 1},
               {"Origin", {0.87, 0.35, 0.0}},
               {"AxisU", {1.0, 0.0, 0.0}},
               {"AxisV", {0.0, 1.0, 0.0}},
               {"Reference", {0.0, 0.0, 0.0}}}}}}}}}}}};
  IoData iodata(config, false);

  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0);
  while (serial_mesh.GetNE() < Mpi::Size(Mpi::World()))
  {
    serial_mesh.UniformRefinement();
  }

  auto parallel_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(parallel_mesh)));

  LaplaceOperator laplace_op(iodata, meshes);
  auto prescribed_config = config;
  prescribed_config["Boundaries"] = {{"Ground", {{"Attributes", {1, 4}}}},
                                     {"PrescribedPotential",
                                      {{{"Index", 1},
                                        {"Attributes", {2}},
                                        {"TerminalAttributes", {3}},
                                        {"DataFile", zero_trace_path.string()}}}}};
  prescribed_config["Solver"]["Electrostatic"] = json::object();
  IoData prescribed_iodata(prescribed_config, false);
  LaplaceOperator prescribed_laplace(prescribed_iodata, meshes);
  auto prescribed_stiffness = prescribed_laplace.GetStiffnessMatrix();
  Vector prescribed_excitation, prescribed_rhs;
  prescribed_laplace.GetExcitationVector(1, *prescribed_stiffness, prescribed_excitation,
                                         prescribed_rhs);
  double prescribed_max = prescribed_excitation.Normlinf();
  Mpi::GlobalMax(1, &prescribed_max, Mpi::World());
  CHECK_THAT(
      prescribed_max,
      WithinRel(1.0 / prescribed_iodata.units.GetScaleFactor<Units::ValueType::VOLTAGE>(),
                1.0e-12));

  SurfaceResponseOperator response(iodata, laplace_op);
  REQUIRE(response.GetBasisSize() == 4);
  REQUIRE(response.GetPatchCount() == 1);
  REQUIRE(response.HasSurfaceResponse());

  const int size = laplace_op.GetH1Space().GetTrueVSize();
  Vector x(size), y(size), Cx, Cy, Ctx;
  for (int i = 0; i < size; i++)
  {
    x(i) = 0.17 + 0.03 * (i + 1) * (Mpi::Rank(Mpi::World()) + 1);
    y(i) = -0.11 + 0.02 * (i + 2) * (Mpi::Rank(Mpi::World()) + 2);
  }
  response.Mult(x, Cx);
  response.Mult(y, Cy);
  response.MultTranspose(x, Ctx);

  double lhs = x * Cy;
  double rhs = Cx * y;
  double norm = Cx * Cx;
  Mpi::GlobalSum(1, &lhs, Mpi::World());
  Mpi::GlobalSum(1, &rhs, Mpi::World());
  Mpi::GlobalSum(1, &norm, Mpi::World());
  CHECK_THAT(lhs, WithinRel(rhs, 1.0e-12));
  CHECK(norm > 0.0);

  Ctx -= Cx;
  double transpose_error = Ctx * Ctx;
  Mpi::GlobalSum(1, &transpose_error, Mpi::World());
  CHECK(transpose_error == 0.0);

  Vector essential_values;
  Cx.GetSubVector(laplace_op.GetDbcTDofList(), essential_values);
  CHECK(essential_values.Normlinf() == 0.0);

  Vector prescribed(size), eliminated_rhs(size);
  prescribed = 0.0;
  const auto &essential = laplace_op.GetDbcTDofList();
  for (int i = 0; i < essential.Size(); i++)
  {
    prescribed(essential[i]) = 0.2 + 0.01 * (i + 1);
  }
  eliminated_rhs = 0.0;
  response.EliminateRHS(prescribed, eliminated_rhs);
  eliminated_rhs.GetSubVector(essential, essential_values);
  CHECK(essential_values.Normlinf() == 0.0);
  double rhs_norm = eliminated_rhs * eliminated_rhs;
  Mpi::GlobalSum(1, &rhs_norm, Mpi::World());
  CHECK(rhs_norm > 0.0);

  const auto energy = response.GetEnergyCorrection(x);
  REQUIRE(energy.interfaces.size() == 1);
  CHECK_THAT(energy.interfaces.at(4), WithinRel(energy.domain, 1.0e-12));
  const auto fabricated_energy = response.GetFabricatedSurfaceEnergy(x);
  REQUIRE(fabricated_energy.size() == 1);
  CHECK(fabricated_energy.at(4) > energy.interfaces.at(4));
  const auto local_electrostatic_response = response.GetElectrostaticResponse(x);
  CHECK_THAT(local_electrostatic_response.domain_correction,
             WithinRel(energy.domain, 1.0e-12));
  CHECK_THAT(local_electrostatic_response.fabricated_surface_energy.at(4),
             WithinRel(fabricated_energy.at(4), 1.0e-12));
  CHECK(std::isfinite(local_electrostatic_response.domain_correction_fixed_flux));
  CHECK(local_electrostatic_response.fabricated_surface_energy_fixed_flux.at(4) > 0.0);
  CHECK(local_electrostatic_response.maximum_trace_closure_spread > 0.0);

  auto separated_config = config;
  separated_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Patches"].push_back(
      {{"Model", 1},
       {"Origin", {0.35, 0.70, 0.0}},
       {"AxisU", {1.0, 0.0, 0.0}},
       {"AxisV", {0.0, 1.0, 0.0}},
       {"Reference", {0.0, 0.0, 0.0}}});
  IoData separated_iodata(separated_config, false);
  SurfaceResponseOperator separated_response(separated_iodata, laplace_op);
  CHECK(separated_response.GetPatchCount() == 2);
  CHECK(separated_response.GetBasisSize() == 8);

  auto overlapping_config = config;
  overlapping_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Patches"].push_back(
      {{"Model", 1},
       {"Origin", {0.90, 0.35, 0.0}},
       {"AxisU", {1.0, 0.0, 0.0}},
       {"AxisV", {0.0, 1.0, 0.0}},
       {"Reference", {0.0, 0.0, 0.0}}});
  IoData overlapping_iodata(overlapping_config, false);
  CHECK_THROWS_WITH(SurfaceResponseOperator(overlapping_iodata, laplace_op),
                    Catch::Matchers::ContainsSubstring("coupled multi-edge coupon model"));

  json automatic_config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", temp.temp_dir.string()}}},
      {"Model", {{"Mesh", "unused.msh"}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
      {"Boundaries",
       {{"Ground", {{"Attributes", {1, 3, 4, 9, 10}}}},
        {"Terminal", {{{"Index", 1}, {"Attributes", {2}}}}},
        {"Postprocessing",
         {{"Dielectric",
           {{{"Index", 4},
             {"Attributes", {9}},
             {"Type", "SA"},
             {"Thickness", 0.002},
             {"Permittivity", 4.0},
             {"EdgeAttributes", {9}},
             {"EdgeDistances", {0.1}},
             {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Electrostatic",
         {{"ResponseCorrection",
           {{"Library", library_path.string()}, {"UnmatchedPolicy", "Error"}}}}}}}};
  IoData automatic_iodata(automatic_config, false);
  automatic_iodata.boundaries.cracked_attributes.insert(9);
  automatic_iodata.boundaries.cracked_attributes.insert(10);

  mfem::Mesh automatic_serial =
      mfem::Mesh::MakeCartesian2D(8, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0);
  for (int face = 0; face < automatic_serial.GetNumFaces(); face++)
  {
    int element1, element2;
    automatic_serial.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    automatic_serial.GetFaceVertices(face, vertices);
    if (vertices.Size() != 2)
    {
      continue;
    }
    const double *p0 = automatic_serial.GetVertex(vertices[0]);
    const double *p1 = automatic_serial.GetVertex(vertices[1]);
    const double xmin = std::min(p0[0], p1[0]);
    const double xmax = std::max(p0[0], p1[0]);
    if (std::abs(p0[1] - 0.5) < 1.0e-12 && std::abs(p1[1] - 0.5) < 1.0e-12 &&
        xmin >= 0.25 - 1.0e-12 && xmax <= 0.75 + 1.0e-12)
    {
      automatic_serial.AddBdrElement(
          automatic_serial.GetFace(face)->Duplicate(&automatic_serial));
      const bool continuation = xmin >= 0.5 - 1.0e-12 && xmax <= 0.625 + 1.0e-12;
      automatic_serial.SetBdrAttribute(automatic_serial.GetNBE() - 1,
                                       continuation ? 10 : 9);
    }
  }
  automatic_serial.FinalizeTopology();
  automatic_serial.Finalize();
  while (automatic_serial.GetNE() < Mpi::Size(Mpi::World()))
  {
    automatic_serial.UniformRefinement();
  }
  auto automatic_parallel = std::make_unique<mfem::ParMesh>(Mpi::World(), automatic_serial);
  std::vector<std::unique_ptr<Mesh>> automatic_meshes;
  automatic_meshes.push_back(std::make_unique<Mesh>(std::move(automatic_parallel)));
  LaplaceOperator automatic_laplace(automatic_iodata, automatic_meshes);
  std::shared_ptr<const SurfaceResponseGeometry> automatic_geometry;
  SurfaceResponseOperator automatic_response(automatic_iodata, automatic_laplace,
                                             &automatic_geometry);
  REQUIRE(automatic_geometry);
  CHECK(automatic_response.GetPatchCount() == 2);
  CHECK(automatic_response.GetBasisSize() == 8);
  CHECK(automatic_response.HasSurfaceResponse());
  SurfaceResponseOperator cached_automatic_response(automatic_iodata, automatic_laplace,
                                                    &automatic_geometry);
  CHECK(cached_automatic_response.GetPatchCount() == automatic_response.GetPatchCount());
  CHECK(cached_automatic_response.GetBasisSize() == automatic_response.GetBasisSize());
  CHECK(cached_automatic_response.HasSurfaceResponse());

  auto thickness_mismatch_config = automatic_config;
  thickness_mismatch_config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Thickness"] =
      0.003;
  IoData thickness_mismatch_iodata(thickness_mismatch_config, false);
  thickness_mismatch_iodata.boundaries.cracked_attributes.insert(9);
  thickness_mismatch_iodata.boundaries.cracked_attributes.insert(10);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(thickness_mismatch_iodata, automatic_laplace),
      Catch::Matchers::ContainsSubstring("does not match fabrication-process response "
                                         "library \"unit-test-process\" thickness"));

  auto permittivity_mismatch_config = automatic_config;
  permittivity_mismatch_config["Boundaries"]["Postprocessing"]["Dielectric"][0]
                              ["Permittivity"] = 4.1;
  IoData permittivity_mismatch_iodata(permittivity_mismatch_config, false);
  permittivity_mismatch_iodata.boundaries.cracked_attributes.insert(9);
  permittivity_mismatch_iodata.boundaries.cracked_attributes.insert(10);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(permittivity_mismatch_iodata, automatic_laplace),
      Catch::Matchers::ContainsSubstring("does not match fabrication-process response "
                                         "library \"unit-test-process\" permittivity"));

  auto missing_layer_config = automatic_config;
  missing_layer_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      missing_layer_library_path.string();
  IoData missing_layer_iodata(missing_layer_config, false);
  missing_layer_iodata.boundaries.cracked_attributes.insert(9);
  missing_layer_iodata.boundaries.cracked_attributes.insert(10);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(missing_layer_iodata, automatic_laplace),
      Catch::Matchers::ContainsSubstring(
          "has a SA surface response but no matching Fabrication.InterfaceLayers entry"));

  auto legacy_config = automatic_config;
  legacy_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      legacy_library_path.string();
  IoData legacy_iodata(legacy_config, false);
  legacy_iodata.boundaries.cracked_attributes.insert(9);
  legacy_iodata.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator legacy_response(legacy_iodata, automatic_laplace);
  CHECK(legacy_response.GetPatchCount() == automatic_response.GetPatchCount());

  auto exact_pair_config_2d = automatic_config;
  exact_pair_config_2d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      exact_pair_library_2d_path.string();
  exact_pair_config_2d["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"] = {
      0.3};
  IoData exact_pair_iodata_2d(exact_pair_config_2d, false);
  exact_pair_iodata_2d.boundaries.cracked_attributes.insert(9);
  exact_pair_iodata_2d.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator exact_pair_response_2d(exact_pair_iodata_2d, automatic_laplace);
  CHECK(exact_pair_response_2d.GetPatchCount() == 1);
  CHECK(exact_pair_response_2d.GetBasisSize() == 4);

  auto interpolated_pair_config_2d = exact_pair_config_2d;
  interpolated_pair_config_2d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      interpolated_pair_library_2d_path.string();
  IoData interpolated_pair_iodata_2d(interpolated_pair_config_2d, false);
  interpolated_pair_iodata_2d.boundaries.cracked_attributes.insert(9);
  interpolated_pair_iodata_2d.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator interpolated_pair_response_2d(interpolated_pair_iodata_2d,
                                                        automatic_laplace);
  CHECK(interpolated_pair_response_2d.GetPatchCount() == 2);
  CHECK(interpolated_pair_response_2d.GetBasisSize() == 8);
  CHECK_THAT(interpolated_pair_response_2d.GetPatchWeight(),
             WithinRel(exact_pair_response_2d.GetPatchWeight(), 1.0e-12));

  Vector pair_probe(automatic_laplace.GetH1Space().GetTrueVSize());
  for (int i = 0; i < pair_probe.Size(); i++)
  {
    pair_probe(i) = std::sin(0.11 * (i + 1 + Mpi::Rank(Mpi::World())));
  }
  Vector exact_pair_correction, interpolated_pair_correction;
  exact_pair_response_2d.Mult(pair_probe, exact_pair_correction);
  interpolated_pair_response_2d.Mult(pair_probe, interpolated_pair_correction);
  interpolated_pair_correction.Add(-1.0, exact_pair_correction);
  CHECK(linalg::Norml2(Mpi::World(), interpolated_pair_correction) <=
        1.0e-12 * std::max(linalg::Norml2(Mpi::World(), exact_pair_correction), 1.0e-300));

  auto invalid_depth_config = automatic_config;
  invalid_depth_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      invalid_library_path.string();
  IoData invalid_depth_iodata(invalid_depth_config, false);
  invalid_depth_iodata.boundaries.cracked_attributes.insert(9);
  invalid_depth_iodata.boundaries.cracked_attributes.insert(10);
  CHECK_THROWS_WITH(SurfaceResponseOperator(invalid_depth_iodata, automatic_laplace),
                    Catch::Matchers::ContainsSubstring(
                        "Fabrication-process response-model CouponDepth must be positive"));

  const auto config_3d_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                              "examples/cpw3d_surface/cpw3d_surface_validation_thin.json";
  std::ifstream config_3d_input(config_3d_path);
  json config_3d = json::parse(config_3d_input);
  config_3d["Problem"]["Output"] = temp.temp_dir.string();
  config_3d["Model"]["Mesh"] =
      (fs::path(PALACE_TEST_DATA_DIR) / "mesh/cpw3d-surface-nc.msh").string();
  config_3d["Solver"]["Electrostatic"]["ResponseCorrection"] = {
      {"Library", library_3d_path.string()},
      {"TargetInterfaces", {1, 2, 3}},
      {"UnmatchedPolicy", "Error"}};
  IoData iodata_3d(config_3d, false);
  auto mesh_3d = mesh::ReadMesh(iodata_3d, Mpi::World());
  const auto geometry_3d = ExtractMetalEdgeGeometry(*mesh_3d, iodata_3d.boundaries);
  const auto segment_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry_3d, 1, InterfaceDielectric::SA);
  double physical_edge_length = 0.0;
  for (const std::size_t segment_index : segment_indices)
  {
    const auto &segment = geometry_3d.segments[segment_index];
    const auto &p0 = geometry_3d.vertices[segment.vertices[0]].coordinate;
    const auto &p1 = geometry_3d.vertices[segment.vertices[1]].coordinate;
    double length_squared = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length_squared += (p1[d] - p0[d]) * (p1[d] - p0[d]);
    }
    physical_edge_length += std::sqrt(length_squared);
  }
  std::vector<std::unique_ptr<Mesh>> meshes_3d;
  meshes_3d.push_back(std::make_unique<Mesh>(std::move(mesh_3d)));
  LaplaceOperator laplace_3d(iodata_3d, meshes_3d);
  SurfaceResponseOperator response_3d(iodata_3d, laplace_3d);
  const auto &line_rule =
      mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * iodata_3d.solver.order);
  CHECK(response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size()) * line_rule.GetNPoints());
  CHECK(response_3d.GetBasisSize() == 4 * response_3d.GetPatchCount());
  CHECK_THAT(response_3d.GetPatchWeight(), WithinRel(0.5 * physical_edge_length, 1.0e-12));
  CHECK(response_3d.HasSurfaceResponse());

  auto maxwell_config_3d = config_3d;
  maxwell_config_3d["Problem"]["Type"] = "Eigenmode";
  maxwell_config_3d["Boundaries"]["Ground"]["Attributes"] = {1, 2};
  maxwell_config_3d["Boundaries"].erase("Terminal");
  maxwell_config_3d["Solver"] = {{"Order", 1},
                                 {"Eigenmode", {{"Target", 1.0}}},
                                 {"SurfaceResponseCorrection",
                                  {{"Library", library_3d_path.string()},
                                   {"TargetInterfaces", {1, 2, 3}},
                                   {"UnmatchedPolicy", "Error"}}}};
  IoData maxwell_iodata_3d(maxwell_config_3d, false);
  auto maxwell_mesh_3d = mesh::ReadMesh(maxwell_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> maxwell_meshes_3d;
  maxwell_meshes_3d.push_back(std::make_unique<Mesh>(std::move(maxwell_mesh_3d)));
  SpaceOperator maxwell_space_3d(maxwell_iodata_3d, maxwell_meshes_3d);
  SurfaceResponseOperator maxwell_response_3d(maxwell_iodata_3d, maxwell_space_3d);
  CHECK(maxwell_response_3d.GetPatchCount() > 0);
  CHECK(maxwell_response_3d.HasSurfaceResponse());
  CHECK(maxwell_response_3d.GetTargetInterfaces() == std::set<int>{1, 2, 3});

  GridFunction maxwell_field(maxwell_space_3d.GetNDSpace(), true);
  mfem::Vector constant_field(3);
  constant_field[0] = 0.7;
  constant_field[1] = -0.4;
  constant_field[2] = 0.2;
  mfem::VectorConstantCoefficient field_coefficient(constant_field);
  maxwell_field.Real().ProjectCoefficient(field_coefficient);
  maxwell_field.Imag() = 0.0;
  const auto real_response = maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  CHECK(std::abs(real_response.domain_correction) > 0.0);
  CHECK(std::abs(real_response.domain_correction_fixed_flux) > 0.0);
  REQUIRE(real_response.fabricated_surface_energy.size() == 3);
  REQUIRE(real_response.fabricated_surface_energy_fixed_flux.size() == 3);
  for (const auto &[interface, energy] : real_response.fabricated_surface_energy_fixed_flux)
  {
    (void)interface;
    CHECK(energy > 0.0);
  }
  CHECK(real_response.loop_residual < 1.0e-10);
  CHECK(real_response.response_weighted_loop_residual < 1.0e-10);
  CHECK(real_response.loop_response_failure_fraction == 0.0);
  CHECK(real_response.kR == 0.0);
  CHECK(real_response.maximum_trace_closure_spread > 0.0);
  CHECK_THAT(real_response.matched_length_fraction, WithinAbs(1.0, 1.0e-12));
  REQUIRE(real_response.matched_length_fraction_by_interface.size() == 3);
  for (const auto &[interface, fraction] :
       real_response.matched_length_fraction_by_interface)
  {
    (void)interface;
    CHECK_THAT(fraction, WithinAbs(1.0, 1.0e-12));
  }

  // The contour reconstruction is also the trace operator used by the self-consistent
  // Maxwell mass correction. Verify its energy identity and transpose symmetry.
  Vector field_true, correction_true;
  maxwell_field.Real().GetTrueDofs(field_true);
  field_true.SetSubVector(maxwell_space_3d.GetNDDbcTDofLists().back(), 0.0);
  maxwell_field.Real().SetFromTrueDofs(field_true);
  const auto operator_response = maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  maxwell_response_3d.Mult(field_true, correction_true);
  CHECK_THAT(0.5 * linalg::Dot(Mpi::World(), field_true, correction_true),
             WithinRel(operator_response.domain_correction, 1.0e-10));

  Vector probe(field_true.Size()), correction_probe;
  auto *probe_data = probe.HostWrite();
  for (int i = 0; i < probe.Size(); i++)
  {
    probe_data[i] = std::sin(0.37 * (i + 1 + 11 * Mpi::Rank(Mpi::World())));
  }
  probe.SetSubVector(maxwell_space_3d.GetNDDbcTDofLists().back(), 0.0);
  maxwell_response_3d.Mult(probe, correction_probe);
  CHECK_THAT(linalg::Dot(Mpi::World(), probe, correction_true),
             WithinRel(linalg::Dot(Mpi::World(), field_true, correction_probe), 1.0e-10));

  // A Maxwell trace reconstructed from E = -grad(V) must reproduce the H1 coupon trace
  // relative to the PEC. This catches a missing contour-voltage gauge even when the
  // contour-loop residual and complex-field scaling tests pass.
  mfem::ParGridFunction potential(&laplace_3d.GetH1Space().Get());
  mfem::FunctionCoefficient potential_coefficient([](const mfem::Vector &x)
                                                  { return x[1]; });
  potential.ProjectCoefficient(potential_coefficient);
  Vector potential_true;
  potential.GetTrueDofs(potential_true);
  mfem::Vector normal_field(3);
  normal_field = 0.0;
  normal_field[1] = -1.0;
  mfem::VectorConstantCoefficient normal_field_coefficient(normal_field);
  maxwell_field.Real().ProjectCoefficient(normal_field_coefficient);
  maxwell_field.Imag() = 0.0;

  const auto electrostatic_correction = response_3d.GetEnergyCorrection(potential_true);
  const auto electrostatic_surfaces =
      response_3d.GetFabricatedSurfaceEnergy(potential_true);
  const auto electrostatic_response = response_3d.GetElectrostaticResponse(potential_true);
  const auto gradient_response = maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  CHECK_THAT(gradient_response.domain_correction,
             WithinRel(electrostatic_correction.domain, 1.0e-10));
  REQUIRE(gradient_response.fabricated_surface_energy.size() ==
          electrostatic_surfaces.size());
  for (const auto &[interface, energy] : electrostatic_surfaces)
  {
    CHECK_THAT(gradient_response.fabricated_surface_energy.at(interface),
               WithinRel(energy, 1.0e-10));
    CHECK_THAT(
        gradient_response.fabricated_surface_energy_fixed_flux.at(interface),
        WithinRel(electrostatic_response.fabricated_surface_energy_fixed_flux.at(interface),
                  1.0e-10));
  }
  CHECK_THAT(gradient_response.domain_correction_fixed_flux,
             WithinRel(electrostatic_response.domain_correction_fixed_flux, 1.0e-10));
  CHECK(gradient_response.loop_residual < 1.0e-10);

  maxwell_field.Real().ProjectCoefficient(field_coefficient);
  maxwell_field.Imag() = maxwell_field.Real();
  const auto complex_response = maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  CHECK_THAT(complex_response.domain_correction,
             WithinRel(2.0 * real_response.domain_correction, 1.0e-10));
  CHECK_THAT(complex_response.domain_correction_fixed_flux,
             WithinRel(2.0 * real_response.domain_correction_fixed_flux, 1.0e-10));
  for (const auto &[interface, energy] : real_response.fabricated_surface_energy)
  {
    CHECK_THAT(complex_response.fabricated_surface_energy.at(interface),
               WithinRel(2.0 * energy, 1.0e-10));
  }
  for (const auto &[interface, energy] : real_response.fabricated_surface_energy_fixed_flux)
  {
    CHECK_THAT(complex_response.fabricated_surface_energy_fixed_flux.at(interface),
               WithinRel(2.0 * energy, 1.0e-10));
  }

  mfem::VectorFunctionCoefficient rotational_field_coefficient(
      3,
      [](const mfem::Vector &x, mfem::Vector &field)
      {
        field[0] = -x[1];
        field[1] = 0.0;
        field[2] = 0.0;
      });
  maxwell_field.Real().ProjectCoefficient(rotational_field_coefficient);
  maxwell_field.Imag() = 0.0;
  const auto rotational_response =
      maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  CHECK(rotational_response.loop_residual > 0.05);
  CHECK(rotational_response.response_weighted_loop_residual > 0.0);
  CHECK(rotational_response.response_weighted_loop_residual <=
        rotational_response.loop_residual);
  CHECK(rotational_response.loop_response_failure_fraction > 0.0);
  CHECK(rotational_response.loop_response_failure_fraction <= 1.0);

  auto coupled_config_3d = config_3d;
  for (auto &interface : coupled_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 7.0};
  }
  coupled_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      coupled_library_3d_path.string();
  IoData coupled_iodata_3d(coupled_config_3d, false);
  coupled_iodata_3d.boundaries.cracked_attributes.insert(1);
  coupled_iodata_3d.boundaries.cracked_attributes.insert(2);
  auto coupled_mesh_3d = mesh::ReadMesh(coupled_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> coupled_meshes_3d;
  coupled_meshes_3d.push_back(std::make_unique<Mesh>(std::move(coupled_mesh_3d)));
  LaplaceOperator coupled_laplace_3d(coupled_iodata_3d, coupled_meshes_3d);
  SurfaceResponseOperator coupled_response_3d(coupled_iodata_3d, coupled_laplace_3d);
  CHECK(coupled_response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size() / 2) * line_rule.GetNPoints());
  CHECK(coupled_response_3d.GetBasisSize() == 5 * coupled_response_3d.GetPatchCount());
  CHECK_THAT(coupled_response_3d.GetPatchWeight(),
             WithinRel(0.25 * physical_edge_length, 1.0e-12));
  CHECK(coupled_response_3d.HasSurfaceResponse());

  auto interpolated_coupled_config_3d = coupled_config_3d;
  interpolated_coupled_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]
                                ["Library"] = interpolated_coupled_library_3d_path.string();
  IoData interpolated_coupled_iodata_3d(interpolated_coupled_config_3d, false);
  interpolated_coupled_iodata_3d.boundaries.cracked_attributes.insert(1);
  interpolated_coupled_iodata_3d.boundaries.cracked_attributes.insert(2);
  auto interpolated_coupled_mesh_3d =
      mesh::ReadMesh(interpolated_coupled_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> interpolated_coupled_meshes_3d;
  interpolated_coupled_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(interpolated_coupled_mesh_3d)));
  LaplaceOperator interpolated_coupled_laplace_3d(interpolated_coupled_iodata_3d,
                                                  interpolated_coupled_meshes_3d);
  SurfaceResponseOperator interpolated_coupled_response_3d(interpolated_coupled_iodata_3d,
                                                           interpolated_coupled_laplace_3d);
  CHECK(interpolated_coupled_response_3d.GetPatchCount() ==
        2 * coupled_response_3d.GetPatchCount());
  CHECK(interpolated_coupled_response_3d.GetBasisSize() ==
        2 * coupled_response_3d.GetBasisSize());
  CHECK_THAT(interpolated_coupled_response_3d.GetPatchWeight(),
             WithinRel(coupled_response_3d.GetPatchWeight(), 1.0e-12));

  auto coupled_maxwell_config_3d = maxwell_config_3d;
  for (auto &interface :
       coupled_maxwell_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 7.0};
  }
  coupled_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      coupled_library_3d_path.string();
  IoData coupled_maxwell_iodata_3d(coupled_maxwell_config_3d, false);
  auto coupled_maxwell_mesh_3d = mesh::ReadMesh(coupled_maxwell_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> coupled_maxwell_meshes_3d;
  coupled_maxwell_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(coupled_maxwell_mesh_3d)));
  SpaceOperator coupled_maxwell_space_3d(coupled_maxwell_iodata_3d,
                                         coupled_maxwell_meshes_3d);
  SurfaceResponseOperator coupled_maxwell_response_3d(coupled_maxwell_iodata_3d,
                                                      coupled_maxwell_space_3d);
  const auto &maxwell_line_rule = mfem::IntRules.Get(
      mfem::Geometry::SEGMENT, 2 * coupled_maxwell_iodata_3d.solver.order);
  CHECK(coupled_maxwell_response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size() / 2) * maxwell_line_rule.GetNPoints());

  auto interpolated_coupled_maxwell_config_3d = coupled_maxwell_config_3d;
  interpolated_coupled_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      interpolated_coupled_library_3d_path.string();
  IoData interpolated_coupled_maxwell_iodata_3d(interpolated_coupled_maxwell_config_3d,
                                                false);
  auto interpolated_coupled_maxwell_mesh_3d =
      mesh::ReadMesh(interpolated_coupled_maxwell_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> interpolated_coupled_maxwell_meshes_3d;
  interpolated_coupled_maxwell_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(interpolated_coupled_maxwell_mesh_3d)));
  SpaceOperator interpolated_coupled_maxwell_space_3d(
      interpolated_coupled_maxwell_iodata_3d, interpolated_coupled_maxwell_meshes_3d);
  SurfaceResponseOperator interpolated_coupled_maxwell_response_3d(
      interpolated_coupled_maxwell_iodata_3d, interpolated_coupled_maxwell_space_3d);
  CHECK(interpolated_coupled_maxwell_response_3d.GetPatchCount() ==
        2 * coupled_maxwell_response_3d.GetPatchCount());
  CHECK(interpolated_coupled_maxwell_response_3d.GetBasisSize() ==
        2 * coupled_maxwell_response_3d.GetBasisSize());

  GridFunction coupled_maxwell_field(coupled_maxwell_space_3d.GetNDSpace(), true);
  auto cpw_potential = [](const mfem::Vector &x)
  {
    const double distance = std::abs(x[0] - 62.0);
    if (distance <= 10.0)
    {
      return 1.0;
    }
    if (distance >= 22.0)
    {
      return 0.0;
    }
    const double t = (distance - 10.0) / 12.0;
    return 1.0 - 3.0 * t * t + 2.0 * t * t * t;
  };
  mfem::FunctionCoefficient transverse_potential_coefficient(cpw_potential);
  mfem::ParGridFunction maxwell_potential(&coupled_maxwell_space_3d.GetH1Space().Get());
  maxwell_potential.ProjectCoefficient(transverse_potential_coefficient);
  mfem::ParDiscreteLinearOperator gradient(&coupled_maxwell_space_3d.GetH1Space().Get(),
                                           &coupled_maxwell_space_3d.GetNDSpace().Get());
  gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
  gradient.Assemble();
  gradient.Mult(maxwell_potential, coupled_maxwell_field.Real());
  coupled_maxwell_field.Real() *= -1.0;
  coupled_maxwell_field.Imag() = 0.0;
  const auto coupled_maxwell_result =
      coupled_maxwell_response_3d.GetMaxwellResponse(coupled_maxwell_field, 0.0);
  // These matrices differ only in the appended V_B - V_A coefficient, so a nonzero
  // domain correction proves that the independent conductor state is active.
  CHECK(std::abs(coupled_maxwell_result.domain_correction) > 0.0);
  CHECK(coupled_maxwell_result.fabricated_surface_energy.size() == 3);
  CHECK(coupled_maxwell_result.loop_residual < 1.0e-10);

  GridFunction interpolated_coupled_maxwell_field(
      interpolated_coupled_maxwell_space_3d.GetNDSpace(), true);
  mfem::ParGridFunction interpolated_maxwell_potential(
      &interpolated_coupled_maxwell_space_3d.GetH1Space().Get());
  interpolated_maxwell_potential.ProjectCoefficient(transverse_potential_coefficient);
  mfem::ParDiscreteLinearOperator interpolated_gradient(
      &interpolated_coupled_maxwell_space_3d.GetH1Space().Get(),
      &interpolated_coupled_maxwell_space_3d.GetNDSpace().Get());
  interpolated_gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
  interpolated_gradient.Assemble();
  interpolated_gradient.Mult(interpolated_maxwell_potential,
                             interpolated_coupled_maxwell_field.Real());
  interpolated_coupled_maxwell_field.Real() *= -1.0;
  interpolated_coupled_maxwell_field.Imag() = 0.0;
  const auto interpolated_coupled_maxwell_result =
      interpolated_coupled_maxwell_response_3d.GetMaxwellResponse(
          interpolated_coupled_maxwell_field, 0.0);
  CHECK_THAT(interpolated_coupled_maxwell_result.maximum_library_distance,
             WithinRel(4.0 / 7.0, 1.0e-12));
  CHECK(interpolated_coupled_maxwell_result.loop_residual < 1.0e-10);

  mfem::ParGridFunction transverse_potential(&coupled_laplace_3d.GetH1Space().Get());
  transverse_potential.ProjectCoefficient(transverse_potential_coefficient);
  Vector transverse_potential_true;
  transverse_potential.GetTrueDofs(transverse_potential_true);
  const auto coupled_electrostatic_result =
      coupled_response_3d.GetElectrostaticResponse(transverse_potential_true);
  mfem::ParGridFunction interpolated_transverse_potential(
      &interpolated_coupled_laplace_3d.GetH1Space().Get());
  interpolated_transverse_potential.ProjectCoefficient(transverse_potential_coefficient);
  Vector interpolated_transverse_potential_true;
  interpolated_transverse_potential.GetTrueDofs(interpolated_transverse_potential_true);
  const auto interpolated_coupled_electrostatic_result =
      interpolated_coupled_response_3d.GetElectrostaticResponse(
          interpolated_transverse_potential_true);
  CHECK_THAT(interpolated_coupled_electrostatic_result.domain_correction,
             WithinRel(coupled_electrostatic_result.domain_correction, 1.0e-12));
  CHECK_THAT(interpolated_coupled_electrostatic_result.domain_correction_fixed_flux,
             WithinRel(coupled_electrostatic_result.domain_correction_fixed_flux, 1.0e-12));
  CHECK_THAT(interpolated_coupled_maxwell_result.domain_correction,
             WithinRel(coupled_maxwell_result.domain_correction, 1.0e-12));
  CHECK_THAT(interpolated_coupled_maxwell_result.domain_correction_fixed_flux,
             WithinRel(coupled_maxwell_result.domain_correction_fixed_flux, 1.0e-12));
  for (const auto &[interface, energy] :
       coupled_electrostatic_result.fabricated_surface_energy)
  {
    CHECK_THAT(
        interpolated_coupled_electrostatic_result.fabricated_surface_energy.at(interface),
        WithinRel(energy, 1.0e-12));
    CHECK_THAT(
        interpolated_coupled_electrostatic_result.fabricated_surface_energy_fixed_flux.at(
            interface),
        WithinRel(
            coupled_electrostatic_result.fabricated_surface_energy_fixed_flux.at(interface),
            1.0e-12));
    CHECK_THAT(
        interpolated_coupled_maxwell_result.fabricated_surface_energy.at(interface),
        WithinRel(coupled_maxwell_result.fabricated_surface_energy.at(interface), 1.0e-12));
    CHECK_THAT(
        interpolated_coupled_maxwell_result.fabricated_surface_energy_fixed_flux.at(
            interface),
        WithinRel(coupled_maxwell_result.fabricated_surface_energy_fixed_flux.at(interface),
                  1.0e-12));
  }
  CHECK_THAT(coupled_maxwell_result.domain_correction,
             WithinRel(coupled_electrostatic_result.domain_correction, 1.0e-10));
  CHECK_THAT(coupled_maxwell_result.domain_correction_fixed_flux,
             WithinRel(coupled_electrostatic_result.domain_correction_fixed_flux, 1.0e-10));
  for (const auto &[interface, energy] :
       coupled_electrostatic_result.fabricated_surface_energy)
  {
    CHECK_THAT(coupled_maxwell_result.fabricated_surface_energy.at(interface),
               WithinRel(energy, 1.0e-5));
    CHECK_THAT(
        coupled_maxwell_result.fabricated_surface_energy_fixed_flux.at(interface),
        WithinRel(
            coupled_electrostatic_result.fabricated_surface_energy_fixed_flux.at(interface),
            1.0e-5));
  }

  // Exercise line reconstruction with a nontrivial discrete gradient. A globally
  // constant field does not detect integration errors when contour lines cross elements.
  mfem::FunctionCoefficient curved_potential_coefficient(
      [cpw_potential](const mfem::Vector &x)
      { return cpw_potential(x) + x[1] * (1.0 + 0.02 * x[0] + 0.03 * x[2]); });
  maxwell_potential.ProjectCoefficient(curved_potential_coefficient);
  gradient.Mult(maxwell_potential, coupled_maxwell_field.Real());
  coupled_maxwell_field.Real() *= -1.0;
  const auto curved_maxwell_result =
      coupled_maxwell_response_3d.GetMaxwellResponse(coupled_maxwell_field, 0.0);

  mfem::ParGridFunction electrostatic_potential(&coupled_laplace_3d.GetH1Space().Get());
  electrostatic_potential.ProjectCoefficient(curved_potential_coefficient);
  Vector electrostatic_potential_true;
  electrostatic_potential.GetTrueDofs(electrostatic_potential_true);
  const auto curved_electrostatic_result =
      coupled_response_3d.GetElectrostaticResponse(electrostatic_potential_true);
  CHECK_THAT(curved_maxwell_result.domain_correction,
             WithinRel(curved_electrostatic_result.domain_correction, 1.0e-10));
  CHECK_THAT(curved_maxwell_result.domain_correction_fixed_flux,
             WithinRel(curved_electrostatic_result.domain_correction_fixed_flux, 1.0e-10));
  for (const auto &[interface, energy] :
       curved_electrostatic_result.fabricated_surface_energy)
  {
    CHECK_THAT(curved_maxwell_result.fabricated_surface_energy.at(interface),
               WithinRel(energy, 1.0e-2));
    CHECK_THAT(
        curved_maxwell_result.fabricated_surface_energy_fixed_flux.at(interface),
        WithinRel(
            curved_electrostatic_result.fabricated_surface_energy_fixed_flux.at(interface),
            1.0e-2));
  }

  // A closed rectangular PEC island supplies true in-plane corners, unlike the CPW
  // extrusion above whose longitudinal physical edges end on truncation boundaries.
  json island_config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", temp.temp_dir.string()}}},
      {"Model", {{"Mesh", "unused.msh"}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
      {"Boundaries",
       {{"Ground", {{"Attributes", {1, 2, 3, 4, 5, 6}}}},
        {"Terminal", {{{"Index", 1}, {"Attributes", {9}}}}},
        {"Postprocessing",
         {{"Dielectric",
           {{{"Index", 4},
             {"Attributes", {9}},
             {"Type", "SA"},
             {"Thickness", 0.002},
             {"Permittivity", 4.0},
             {"AutomaticEdges", true},
             {"EdgeDistances", {0.2}},
             {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Electrostatic",
         {{"ResponseCorrection",
           {{"Library", concave_library_3d_path.string()},
            {"TargetInterfaces", {4}},
            {"UnmatchedPolicy", "Error"}}}}}}}};
  auto MakeIslandMesh =
      [](bool rounded = false, bool tetrahedral = false, bool aperture = false,
         bool neighboring_island = false)
  {
    const double in_plane_extent = aperture || neighboring_island ? 2.0 : 1.0;
    const double center = 0.5 * in_plane_extent;
    const double half_width = 0.25;
    const int in_plane_elements =
        (rounded ? 16 : 8) * (aperture || neighboring_island ? 2 : 1);
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(in_plane_elements, 4, in_plane_elements,
                                                    tetrahedral ? mfem::Element::TETRAHEDRON
                                                                : mfem::Element::HEXAHEDRON,
                                                    in_plane_extent, 1.0, in_plane_extent);
    for (int face = 0; face < serial.GetNumFaces(); face++)
    {
      int element1, element2;
      serial.GetFaceElements(face, &element1, &element2);
      if (element1 < 0 || element2 < 0)
      {
        continue;
      }
      mfem::Array<int> vertices;
      serial.GetFaceVertices(face, vertices);
      bool on_plane = true;
      double xmin = in_plane_extent, xmax = 0.0;
      double zmin = in_plane_extent, zmax = 0.0;
      for (const int vertex : vertices)
      {
        const double *point = serial.GetVertex(vertex);
        on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
        xmin = std::min(xmin, point[0]);
        xmax = std::max(xmax, point[0]);
        zmin = std::min(zmin, point[2]);
        zmax = std::max(zmax, point[2]);
      }
      const double island_center_x = neighboring_island ? 0.625 : center;
      const bool inside_island =
          xmin >= island_center_x - half_width - 1.0e-12 &&
          xmax <= island_center_x + half_width + 1.0e-12 &&
          zmin >= center - half_width - 1.0e-12 && zmax <= center + half_width + 1.0e-12;
      constexpr double neighbor_center_x = 1.375;
      const bool inside_neighbor =
          neighboring_island && xmin >= neighbor_center_x - half_width - 1.0e-12 &&
          xmax <= neighbor_center_x + half_width + 1.0e-12 &&
          zmin >= center - half_width - 1.0e-12 && zmax <= center + half_width + 1.0e-12;
      if (on_plane &&
          (aperture ? !inside_island : (inside_island || inside_neighbor)))
      {
        serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
        serial.SetBdrAttribute(serial.GetNBE() - 1, inside_neighbor ? 10 : 9);
      }
    }
    if (rounded)
    {
      constexpr double radius = 0.125;
      constexpr double tolerance = 1.0e-12;
      for (int vertex = 0; vertex < serial.GetNV(); vertex++)
      {
        double *point = serial.GetVertex(vertex);
        if (std::abs(point[1] - 0.5) > tolerance)
        {
          continue;
        }
        for (const double sign_x : {-1.0, 1.0})
        {
          for (const double sign_z : {-1.0, 1.0})
          {
            const double corner_x = center + sign_x * half_width;
            const double corner_z = center + sign_z * half_width;
            const double center_x = corner_x - sign_x * radius;
            const double center_z = corner_z - sign_z * radius;
            const double local_x = sign_x * (point[0] - center_x);
            const double local_z = sign_z * (point[2] - center_z);
            if (local_x < -tolerance || local_x > radius + tolerance ||
                local_z < -tolerance || local_z > radius + tolerance)
            {
              continue;
            }

            double angle = 0.0;
            if (std::abs(point[2] - corner_z) <= tolerance)
            {
              angle = 0.5 * std::acos(-1.0) - 0.25 * std::acos(-1.0) * local_x / radius;
            }
            else if (std::abs(point[0] - corner_x) <= tolerance)
            {
              angle = 0.25 * std::acos(-1.0) * local_z / radius;
            }
            else
            {
              continue;
            }
            point[0] = center_x + sign_x * radius * std::cos(angle);
            point[2] = center_z + sign_z * radius * std::sin(angle);
          }
        }
      }
    }
    serial.FinalizeTopology();
    serial.Finalize();
    return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
  };

  IoData concave_island_iodata(island_config, false);
  concave_island_iodata.boundaries.cracked_attributes.insert(9);
  auto island_geometry_mesh = MakeIslandMesh();
  const auto island_geometry =
      ExtractMetalEdgeGeometry(*island_geometry_mesh, concave_island_iodata.boundaries);
  const auto island_segments =
      GetInterfaceMetalEdgeSegmentIndices(island_geometry, 4, InterfaceDielectric::SA);
  std::set<std::size_t> island_vertices;
  double island_perimeter = 0.0;
  for (const std::size_t segment_index : island_segments)
  {
    const auto &segment = island_geometry.segments[segment_index];
    island_vertices.insert(segment.vertices.begin(), segment.vertices.end());
    const auto &p0 = island_geometry.vertices[segment.vertices[0]].coordinate;
    const auto &p1 = island_geometry.vertices[segment.vertices[1]].coordinate;
    double length_squared = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length_squared += (p1[d] - p0[d]) * (p1[d] - p0[d]);
    }
    island_perimeter += std::sqrt(length_squared);
  }
  const int island_corners = static_cast<int>(
      std::count_if(island_vertices.begin(), island_vertices.end(),
                    [&](std::size_t vertex)
                    {
                      return island_geometry.vertices[vertex].physical_type ==
                             MetalEdgeVertexType::CORNER;
                    }));
  REQUIRE(island_corners == 4);
  CHECK_THAT(island_perimeter, WithinAbs(2.0, 1.0e-12));

  std::vector<std::unique_ptr<Mesh>> concave_island_meshes;
  concave_island_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
  LaplaceOperator concave_island_laplace(concave_island_iodata, concave_island_meshes);
  SurfaceResponseOperator concave_island_response(concave_island_iodata,
                                                  concave_island_laplace);
  const auto &island_line_rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2);
  CHECK(concave_island_response.GetPatchCount() ==
        static_cast<int>(island_segments.size()) * island_line_rule.GetNPoints());
  CHECK(concave_island_response.GetBasisSize() ==
        4 * concave_island_response.GetPatchCount());
  CHECK_THAT(concave_island_response.GetPatchWeight(),
             WithinRel(island_perimeter / 0.2, 1.0e-12));

  auto convex_island_config = island_config;
  convex_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      convex_library_3d_path.string();
  IoData convex_island_iodata(convex_island_config, false);
  convex_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> convex_island_meshes;
  convex_island_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
  LaplaceOperator convex_island_laplace(convex_island_iodata, convex_island_meshes);
  SurfaceResponseOperator convex_island_response(convex_island_iodata,
                                                 convex_island_laplace);
  const int removed_straight_patches = 2 * island_corners * island_line_rule.GetNPoints();
  CHECK(convex_island_response.GetPatchCount() == concave_island_response.GetPatchCount() -
                                                      removed_straight_patches +
                                                      island_corners);
  CHECK(convex_island_response.GetBasisSize() ==
        4 * (convex_island_response.GetPatchCount() - island_corners) +
            12 * island_corners);
  const double expected_convex_weight =
      (island_perimeter - 2.0 * island_corners * 0.2) / 0.2 + island_corners;
  CHECK_THAT(convex_island_response.GetPatchWeight(),
             WithinRel(expected_convex_weight, 1.0e-12));

  auto rounded_island_config = island_config;
  rounded_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      rounded_library_3d_path.string();
  IoData rounded_island_iodata(rounded_island_config, false);
  rounded_island_iodata.boundaries.cracked_attributes.insert(9);
  auto rounded_geometry_mesh = MakeIslandMesh(true);
  const auto rounded_geometry =
      ExtractMetalEdgeGeometry(*rounded_geometry_mesh, rounded_island_iodata.boundaries);
  const auto rounded_segments =
      GetInterfaceMetalEdgeSegmentIndices(rounded_geometry, 4, InterfaceDielectric::SA);
  std::set<std::size_t> rounded_vertices;
  for (const std::size_t segment_index : rounded_segments)
  {
    const auto &segment = rounded_geometry.segments[segment_index];
    rounded_vertices.insert(segment.vertices.begin(), segment.vertices.end());
  }
  CHECK(std::none_of(rounded_vertices.begin(), rounded_vertices.end(),
                     [&](std::size_t vertex)
                     {
                       return rounded_geometry.vertices[vertex].physical_type !=
                              MetalEdgeVertexType::REGULAR;
                     }));

  std::vector<std::unique_ptr<Mesh>> rounded_island_meshes;
  rounded_island_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh(true)));
  LaplaceOperator rounded_island_laplace(rounded_island_iodata, rounded_island_meshes);
  SurfaceResponseOperator rounded_island_response(rounded_island_iodata,
                                                  rounded_island_laplace);
  constexpr int rounded_corner_count = 4;
  constexpr int remaining_straight_intervals = 8;
  CHECK(rounded_island_response.GetPatchCount() ==
        remaining_straight_intervals * island_line_rule.GetNPoints() +
            rounded_corner_count);
  CHECK(rounded_island_response.GetBasisSize() ==
        4 * remaining_straight_intervals * island_line_rule.GetNPoints() +
            12 * rounded_corner_count);
  CHECK_THAT(rounded_island_response.GetPatchWeight(), WithinRel(6.0, 1.0e-12));

  auto rounded_concave_island_config = island_config;
  rounded_concave_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]
                               ["Library"] = rounded_concave_library_3d_path.string();
  rounded_concave_island_config["Boundaries"]["Ground"]["Attributes"].push_back(9);
  rounded_concave_island_config["Boundaries"].erase("Terminal");
  rounded_concave_island_config["Boundaries"]["Postprocessing"]["Dielectric"][0]
                               ["EdgeExcludeAttributes"] = {1, 2, 3, 4, 5, 6};
  IoData rounded_concave_island_iodata(rounded_concave_island_config, false);
  rounded_concave_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> rounded_concave_island_meshes;
  rounded_concave_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, false, true)));
  LaplaceOperator rounded_concave_island_laplace(rounded_concave_island_iodata,
                                                 rounded_concave_island_meshes);
  SurfaceResponseOperator rounded_concave_island_response(rounded_concave_island_iodata,
                                                          rounded_concave_island_laplace);
  CHECK(rounded_concave_island_response.GetPatchCount() ==
        remaining_straight_intervals * island_line_rule.GetNPoints() +
            rounded_corner_count);
  CHECK(rounded_concave_island_response.GetBasisSize() ==
        4 * remaining_straight_intervals * island_line_rule.GetNPoints() +
            12 * rounded_corner_count);
  CHECK_THAT(rounded_concave_island_response.GetPatchWeight(), WithinRel(6.0, 1.0e-12));

  auto interpolated_rounded_island_config = rounded_island_config;
  interpolated_rounded_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]
                                    ["Library"] =
                                        interpolated_rounded_library_3d_path.string();
  IoData interpolated_rounded_island_iodata(interpolated_rounded_island_config, false);
  interpolated_rounded_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> interpolated_rounded_island_meshes;
  interpolated_rounded_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true)));
  LaplaceOperator interpolated_rounded_island_laplace(interpolated_rounded_island_iodata,
                                                      interpolated_rounded_island_meshes);
  SurfaceResponseOperator interpolated_rounded_island_response(
      interpolated_rounded_island_iodata, interpolated_rounded_island_laplace);
  CHECK(interpolated_rounded_island_response.GetPatchCount() ==
        remaining_straight_intervals * island_line_rule.GetNPoints() +
            2 * rounded_corner_count);
  CHECK(interpolated_rounded_island_response.GetBasisSize() ==
        4 * remaining_straight_intervals * island_line_rule.GetNPoints() +
            2 * 12 * rounded_corner_count);
  CHECK_THAT(interpolated_rounded_island_response.GetPatchWeight(),
             WithinRel(6.0, 1.0e-12));
  Vector interpolation_probe(rounded_island_response.Height());
  auto *interpolation_probe_data = interpolation_probe.HostWrite();
  for (int i = 0; i < interpolation_probe.Size(); i++)
  {
    interpolation_probe_data[i] = std::cos(0.17 * (i + 1 + 3 * Mpi::Rank(Mpi::World())));
  }
  Vector rounded_correction, interpolated_rounded_correction;
  rounded_island_response.Mult(interpolation_probe, rounded_correction);
  interpolated_rounded_island_response.Mult(interpolation_probe,
                                            interpolated_rounded_correction);
  interpolated_rounded_correction.Add(-1.0, rounded_correction);
  CHECK(linalg::Norml2(Mpi::World(), interpolated_rounded_correction) <=
        1.0e-12 * std::max(linalg::Norml2(Mpi::World(), rounded_correction), 1.0e-300));

  auto convex_maxwell_island_config = convex_island_config;
  convex_maxwell_island_config["Problem"]["Type"] = "Eigenmode";
  convex_maxwell_island_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4,
                                                                        5, 6, 9};
  convex_maxwell_island_config["Boundaries"].erase("Terminal");
  convex_maxwell_island_config["Solver"] = {{"Order", 1},
                                            {"Eigenmode", {{"Target", 1.0}}},
                                            {"SurfaceResponseCorrection",
                                             {{"Library", convex_library_3d_path.string()},
                                              {"TargetInterfaces", {4}},
                                              {"UnmatchedPolicy", "Error"}}}};
  IoData convex_maxwell_island_iodata(convex_maxwell_island_config, false);
  convex_maxwell_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> convex_maxwell_island_meshes;
  convex_maxwell_island_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
  SpaceOperator convex_maxwell_island_space(convex_maxwell_island_iodata,
                                            convex_maxwell_island_meshes);
  SurfaceResponseOperator convex_maxwell_island_response(convex_maxwell_island_iodata,
                                                         convex_maxwell_island_space);
  CHECK(convex_maxwell_island_response.GetPatchCount() ==
        static_cast<int>(island_segments.size()) * island_line_rule.GetNPoints() -
            removed_straight_patches + island_corners);

  GridFunction island_field(convex_maxwell_island_space.GetNDSpace(), true);
  island_field.Real().ProjectCoefficient(field_coefficient);
  island_field.Imag() = 0.0;
  const auto constant_island_response =
      convex_maxwell_island_response.GetMaxwellResponse(island_field, 0.0);
  CHECK(constant_island_response.loop_residual < 1.0e-10);
  CHECK(constant_island_response.corner_neighborhood_fraction == 0.0);

  Vector island_true, island_correction, island_probe, island_probe_correction;
  island_field.Real().GetTrueDofs(island_true);
  auto *island_data = island_true.HostWrite();
  for (int i = 0; i < island_true.Size(); i++)
  {
    island_data[i] = std::cos(0.23 * (i + 1 + 7 * Mpi::Rank(Mpi::World())));
  }
  island_true.SetSubVector(convex_maxwell_island_space.GetNDDbcTDofLists().back(), 0.0);
  island_field.Real().SetFromTrueDofs(island_true);
  const auto random_island_response =
      convex_maxwell_island_response.GetMaxwellResponse(island_field, 0.0);
  convex_maxwell_island_response.Mult(island_true, island_correction);
  CHECK_THAT(0.5 * linalg::Dot(Mpi::World(), island_true, island_correction),
             WithinRel(random_island_response.domain_correction, 1.0e-10));
  island_probe.SetSize(island_true.Size());
  auto *island_probe_data = island_probe.HostWrite();
  for (int i = 0; i < island_probe.Size(); i++)
  {
    island_probe_data[i] = std::sin(0.31 * (i + 1 + 5 * Mpi::Rank(Mpi::World())));
  }
  island_probe.SetSubVector(convex_maxwell_island_space.GetNDDbcTDofLists().back(), 0.0);
  convex_maxwell_island_response.Mult(island_probe, island_probe_correction);
  CHECK_THAT(
      linalg::Dot(Mpi::World(), island_probe, island_correction),
      WithinRel(linalg::Dot(Mpi::World(), island_true, island_probe_correction), 1.0e-10));

  // A rounded-corner reference lies inside the PEC footprint. On a tetrahedral mesh,
  // exact line integration rejects an anchor path through that internal boundary. The
  // process-plane contour instead starts at a library-declared zero-trace knot.
  auto rounded_maxwell_island_config = rounded_island_config;
  rounded_maxwell_island_config["Problem"]["Type"] = "Eigenmode";
  rounded_maxwell_island_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4,
                                                                         5, 6, 9};
  rounded_maxwell_island_config["Boundaries"].erase("Terminal");
  rounded_maxwell_island_config["Solver"] = {
      {"Order", 1},
      {"Eigenmode", {{"Target", 1.0}}},
      {"SurfaceResponseCorrection",
       {{"Library", rounded_library_3d_path.string()},
        {"TargetInterfaces", {4}},
        {"UnmatchedPolicy", "Error"}}}};
  IoData rounded_maxwell_island_iodata(rounded_maxwell_island_config, false);
  rounded_maxwell_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> rounded_maxwell_island_meshes;
  rounded_maxwell_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, true)));
  SpaceOperator rounded_maxwell_island_space(rounded_maxwell_island_iodata,
                                             rounded_maxwell_island_meshes);
  SurfaceResponseOperator rounded_maxwell_island_response(rounded_maxwell_island_iodata,
                                                          rounded_maxwell_island_space);
  CHECK(rounded_maxwell_island_response.GetPatchCount() > rounded_corner_count);

  GridFunction rounded_island_field(rounded_maxwell_island_space.GetNDSpace(), true);
  rounded_island_field.Real().ProjectCoefficient(field_coefficient);
  rounded_island_field.Imag() = 0.0;
  const auto rounded_island_result =
      rounded_maxwell_island_response.GetMaxwellResponse(rounded_island_field, 0.0);
  CHECK(rounded_island_result.loop_residual < 1.0e-10);

  auto rounded_concave_maxwell_config = rounded_maxwell_island_config;
  rounded_concave_maxwell_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      rounded_concave_library_3d_path.string();
  rounded_concave_maxwell_config["Boundaries"]["Postprocessing"]["Dielectric"][0]
                                ["EdgeExcludeAttributes"] = {1, 2, 3, 4, 5, 6};
  IoData rounded_concave_maxwell_iodata(rounded_concave_maxwell_config, false);
  rounded_concave_maxwell_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> rounded_concave_maxwell_meshes;
  rounded_concave_maxwell_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, true, true)));
  SpaceOperator rounded_concave_maxwell_space(rounded_concave_maxwell_iodata,
                                              rounded_concave_maxwell_meshes);
  SurfaceResponseOperator rounded_concave_maxwell_response(
      rounded_concave_maxwell_iodata, rounded_concave_maxwell_space);
  CHECK(rounded_concave_maxwell_response.GetPatchCount() ==
        rounded_maxwell_island_response.GetPatchCount());

  GridFunction rounded_concave_maxwell_field(rounded_concave_maxwell_space.GetNDSpace(),
                                             true);
  rounded_concave_maxwell_field.Real().ProjectCoefficient(field_coefficient);
  rounded_concave_maxwell_field.Imag() = 0.0;
  const auto rounded_concave_maxwell_result =
      rounded_concave_maxwell_response.GetMaxwellResponse(rounded_concave_maxwell_field,
                                                         0.0);
  CHECK(rounded_concave_maxwell_result.loop_residual < 1.0e-10);
  CHECK(rounded_concave_maxwell_result.corner_neighborhood_fraction == 0.0);

  auto interpolated_rounded_maxwell_island_config = rounded_maxwell_island_config;
  interpolated_rounded_maxwell_island_config["Solver"]["SurfaceResponseCorrection"]
                                            ["Library"] =
                                                interpolated_rounded_library_3d_path
                                                    .string();
  IoData interpolated_rounded_maxwell_island_iodata(
      interpolated_rounded_maxwell_island_config, false);
  interpolated_rounded_maxwell_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> interpolated_rounded_maxwell_island_meshes;
  interpolated_rounded_maxwell_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, true)));
  SpaceOperator interpolated_rounded_maxwell_island_space(
      interpolated_rounded_maxwell_island_iodata,
      interpolated_rounded_maxwell_island_meshes);
  SurfaceResponseOperator interpolated_rounded_maxwell_island_response(
      interpolated_rounded_maxwell_island_iodata,
      interpolated_rounded_maxwell_island_space);
  GridFunction interpolated_rounded_island_field(
      interpolated_rounded_maxwell_island_space.GetNDSpace(), true);
  interpolated_rounded_island_field.Real().ProjectCoefficient(field_coefficient);
  interpolated_rounded_island_field.Imag() = 0.0;
  const auto interpolated_rounded_island_result =
      interpolated_rounded_maxwell_island_response.GetMaxwellResponse(
          interpolated_rounded_island_field, 0.0);
  CHECK(interpolated_rounded_island_result.loop_residual < 1.0e-10);
  CHECK_THAT(interpolated_rounded_island_result.maximum_library_distance,
             WithinRel(0.25, 1.0e-12));
  CHECK_THAT(interpolated_rounded_island_result.domain_correction,
             WithinRel(rounded_island_result.domain_correction, 1.0e-12));
  for (const auto &[interface, energy] : rounded_island_result.fabricated_surface_energy)
  {
    CHECK_THAT(interpolated_rounded_island_result.fabricated_surface_energy.at(interface),
               WithinRel(energy, 1.0e-12));
    CHECK_THAT(
        interpolated_rounded_island_result.fabricated_surface_energy_fixed_flux.at(
            interface),
        WithinRel(rounded_island_result.fabricated_surface_energy_fixed_flux.at(interface),
                  1.0e-12));
  }

  // A local interface-signature conflict must not invalidate every segment carrying that
  // signature. The two islands are separated by less than 2R only along their facing
  // sides, so Warn omits those local segments while Error remains strict.
  auto mixed_signature_config = convex_maxwell_island_config;
  mixed_signature_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5,
                                                                   6, 9, 10};
  mixed_signature_config["Boundaries"]["Postprocessing"]["Dielectric"] = {
      {{"Index", 4},
       {"Attributes", {9}},
       {"Type", "SA"},
       {"Thickness", 0.002},
       {"Permittivity", 4.0},
       {"AutomaticEdges", true},
       {"EdgeDistances", {0.2}},
       {"EdgeFrameNormal", {0.0, 1.0, 0.0}}},
      {{"Index", 5},
       {"Attributes", {10}},
       {"Type", "MS"},
       {"Thickness", 0.002},
       {"Permittivity", 11.47},
       {"AutomaticEdges", true},
       {"EdgeDistances", {0.2}},
       {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}};
  mixed_signature_config["Solver"]["SurfaceResponseCorrection"]["TargetInterfaces"] = {4,
                                                                                       5};
  mixed_signature_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] =
      "Warn";
  IoData mixed_signature_iodata(mixed_signature_config, false);
  mixed_signature_iodata.boundaries.cracked_attributes.insert(9);
  mixed_signature_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> mixed_signature_meshes;
  mixed_signature_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(false, false, false, true)));
  SpaceOperator mixed_signature_space(mixed_signature_iodata, mixed_signature_meshes);
  SurfaceResponseOperator mixed_signature_response(mixed_signature_iodata,
                                                   mixed_signature_space);
  CHECK(mixed_signature_response.GetPatchCount() > 0);

  GridFunction mixed_signature_field(mixed_signature_space.GetNDSpace(), true);
  mixed_signature_field.Real().ProjectCoefficient(field_coefficient);
  mixed_signature_field.Imag() = 0.0;
  const auto mixed_signature_result =
      mixed_signature_response.GetMaxwellResponse(mixed_signature_field, 0.0);
  CHECK(mixed_signature_result.matched_length_fraction > 0.0);
  CHECK(mixed_signature_result.matched_length_fraction < 1.0);
  CHECK(mixed_signature_result.fabricated_surface_energy.count(4) == 1);
  CHECK(mixed_signature_result.fabricated_surface_energy.count(5) == 1);

  auto local_interaction_config = mixed_signature_config;
  local_interaction_config["Boundaries"]["Postprocessing"]["Dielectric"] = {
      {{"Index", 4},
       {"Attributes", {9, 10}},
       {"Type", "SA"},
       {"Thickness", 0.002},
       {"Permittivity", 4.0},
       {"AutomaticEdges", true},
       {"EdgeDistances", {0.2}},
       {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}};
  local_interaction_config["Solver"]["SurfaceResponseCorrection"]["TargetInterfaces"] = {
      4};
  IoData local_interaction_iodata(local_interaction_config, false);
  local_interaction_iodata.boundaries.cracked_attributes.insert(9);
  local_interaction_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> local_interaction_meshes;
  local_interaction_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(false, false, false, true)));
  SpaceOperator local_interaction_space(local_interaction_iodata,
                                        local_interaction_meshes);
  SurfaceResponseOperator local_interaction_response(local_interaction_iodata,
                                                     local_interaction_space);
  CHECK(local_interaction_response.GetPatchCount() > 0);

  GridFunction local_interaction_field(local_interaction_space.GetNDSpace(), true);
  local_interaction_field.Real().ProjectCoefficient(field_coefficient);
  local_interaction_field.Imag() = 0.0;
  const auto local_interaction_result =
      local_interaction_response.GetMaxwellResponse(local_interaction_field, 0.0);
  CHECK(local_interaction_result.matched_length_fraction > 0.0);
  CHECK(local_interaction_result.matched_length_fraction < 1.0);
  CHECK(local_interaction_result.fabricated_surface_energy.count(4) == 1);

  // Two aperture perimeters are disconnected edge-graph components even though the
  // metal bridge between them is one physical PEC strip. Its outward-facing edges must
  // select a strip coupon from local geometry instead of being rejected as different
  // conductors.
  auto MakePairedApertureMesh = []()
  {
    constexpr double extent = 3.0;
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(
        24, 4, 24, mfem::Element::HEXAHEDRON, extent, 1.0, extent);
    for (int face = 0; face < serial.GetNumFaces(); face++)
    {
      int element1, element2;
      serial.GetFaceElements(face, &element1, &element2);
      if (element1 < 0 || element2 < 0)
      {
        continue;
      }
      mfem::Array<int> vertices;
      serial.GetFaceVertices(face, vertices);
      bool on_plane = true;
      double xmin = extent, xmax = 0.0;
      double zmin = extent, zmax = 0.0;
      for (const int vertex : vertices)
      {
        const double *point = serial.GetVertex(vertex);
        on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
        xmin = std::min(xmin, point[0]);
        xmax = std::max(xmax, point[0]);
        zmin = std::min(zmin, point[2]);
        zmax = std::max(zmax, point[2]);
      }
      const bool left_aperture =
          xmin >= 0.75 - 1.0e-12 && xmax <= 1.25 + 1.0e-12 &&
          zmin >= 0.5 - 1.0e-12 && zmax <= 2.5 + 1.0e-12;
      const bool right_aperture =
          xmin >= 1.5 - 1.0e-12 && xmax <= 2.0 + 1.0e-12 &&
          zmin >= 0.5 - 1.0e-12 && zmax <= 2.5 + 1.0e-12;
      if (on_plane && !left_aperture && !right_aperture)
      {
        serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
        serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
      }
    }
    serial.FinalizeTopology();
    serial.Finalize();
    return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
  };

  auto strip_aperture_config = convex_maxwell_island_config;
  strip_aperture_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      strip_library_3d_path.string();
  strip_aperture_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] =
      "Warn";
  strip_aperture_config["Boundaries"]["Postprocessing"]["Dielectric"][0]
                       ["EdgeExcludeAttributes"] = {1, 2, 3, 4, 5, 6};
  IoData strip_aperture_iodata(strip_aperture_config, false);
  strip_aperture_iodata.boundaries.cracked_attributes.insert(9);
  auto strip_geometry_mesh = MakePairedApertureMesh();
  const auto strip_geometry =
      ExtractMetalEdgeGeometry(*strip_geometry_mesh, strip_aperture_iodata.boundaries);
  auto strip_segments =
      GetInterfaceMetalEdgeSegmentIndices(strip_geometry, 4, InterfaceDielectric::SA);
  ExcludeMetalEdgeSegmentIndices(*strip_geometry_mesh, strip_geometry, {1, 2, 3, 4, 5, 6},
                                 strip_segments);
  std::set<int> strip_perimeter_components;
  for (const std::size_t segment : strip_segments)
  {
    strip_perimeter_components.insert(strip_geometry.segments[segment].component);
  }
  CHECK(strip_perimeter_components.size() == 2);

  std::vector<std::unique_ptr<Mesh>> strip_aperture_meshes;
  strip_aperture_meshes.push_back(
      std::make_unique<Mesh>(MakePairedApertureMesh()));
  SpaceOperator strip_aperture_space(strip_aperture_iodata, strip_aperture_meshes);
  SurfaceResponseOperator strip_aperture_response(strip_aperture_iodata,
                                                  strip_aperture_space);
  CHECK(strip_aperture_response.GetPatchCount() > 0);

  auto no_strip_aperture_config = strip_aperture_config;
  no_strip_aperture_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      concave_library_3d_path.string();
  IoData no_strip_aperture_iodata(no_strip_aperture_config, false);
  no_strip_aperture_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> no_strip_aperture_meshes;
  no_strip_aperture_meshes.push_back(
      std::make_unique<Mesh>(MakePairedApertureMesh()));
  SpaceOperator no_strip_aperture_space(no_strip_aperture_iodata,
                                        no_strip_aperture_meshes);
  SurfaceResponseOperator no_strip_aperture_response(no_strip_aperture_iodata,
                                                     no_strip_aperture_space);

  GridFunction strip_aperture_field(strip_aperture_space.GetNDSpace(), true);
  strip_aperture_field.Real().ProjectCoefficient(field_coefficient);
  strip_aperture_field.Imag() = 0.0;
  const auto strip_aperture_result =
      strip_aperture_response.GetMaxwellResponse(strip_aperture_field, 0.0);
  GridFunction no_strip_aperture_field(no_strip_aperture_space.GetNDSpace(), true);
  no_strip_aperture_field.Real().ProjectCoefficient(field_coefficient);
  no_strip_aperture_field.Imag() = 0.0;
  const auto no_strip_aperture_result =
      no_strip_aperture_response.GetMaxwellResponse(no_strip_aperture_field, 0.0);
  CHECK(strip_aperture_result.matched_length_fraction >
        no_strip_aperture_result.matched_length_fraction);

  auto strict_mixed_signature_config = mixed_signature_config;
  strict_mixed_signature_config["Solver"]["SurfaceResponseCorrection"]
                               ["UnmatchedPolicy"] = "Error";
  IoData strict_mixed_signature_iodata(strict_mixed_signature_config, false);
  strict_mixed_signature_iodata.boundaries.cracked_attributes.insert(9);
  strict_mixed_signature_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> strict_mixed_signature_meshes;
  strict_mixed_signature_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(false, false, false, true)));
  SpaceOperator strict_mixed_signature_space(strict_mixed_signature_iodata,
                                             strict_mixed_signature_meshes);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(strict_mixed_signature_iodata, strict_mixed_signature_space),
      Catch::Matchers::ContainsSubstring("different interface mapping"));
#endif
}

}  // namespace palace
