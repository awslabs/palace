// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"

#include <array>
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
  const auto invalid_library_path =
      temp.temp_dir / "fabrication-process-invalid-depth.json";
  const auto library_3d_path = temp.temp_dir / "fabrication-process-3d.json";
  const auto coupled_library_3d_path =
      temp.temp_dir / "fabrication-process-coupled-3d.json";
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
    auto write_domain_matrix =
        [](const auto &path, const std::array<std::array<double, 4>, 4> &matrix)
    {
      std::ofstream output(path);
      output << "basis_i,basis_j,Q_ij (J)\n";
      for (int i = 0; i < 4; i++)
      {
        for (int j = i; j < 4; j++)
        {
          output << i + 1 << "," << j + 1 << "," << matrix[i][j] << "\n";
        }
      }
    };
    auto write_surface_matrix =
        [](const auto &path, const std::array<std::array<double, 4>, 4> &matrix)
    {
      std::ofstream output(path);
      output << "interface,edge,basis_i,basis_j,Q_total_ij (J)\n";
      for (int edge = 1; edge <= 2; edge++)
      {
        for (int i = 0; i < 4; i++)
        {
          for (int j = i; j < 4; j++)
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
    const json library = {{"Version", 1},
                          {"Name", "unit-test-process"},
                          {"MatchingRadius", 0.1},
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
    coupled_library_3d["Name"] = "unit-test-process-coupled-3d";
    coupled_library_3d["MatchingRadius"] = 7.0;
    auto coupled_model = coupled_library_3d["Models"][0];
    coupled_model["Name"] = "terminal-ground-gap-12um";
    coupled_model["Topology"] = "DifferentConductorGap";
    coupled_model["Separation"] = 12.0;
    coupled_model["SeparationTolerance"] = 1.0e-6;
    coupled_library_3d["Models"].push_back(std::move(coupled_model));
    std::ofstream coupled_output_3d(coupled_library_3d_path);
    coupled_output_3d << coupled_library_3d.dump(2) << "\n";
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
  SurfaceResponseOperator automatic_response(automatic_iodata, automatic_laplace);
  CHECK(automatic_response.GetPatchCount() == 2);
  CHECK(automatic_response.GetBasisSize() == 8);
  CHECK(automatic_response.HasSurfaceResponse());

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
  maxwell_config_3d["Solver"] = {
      {"Order", 1},
      {"Eigenmode", {{"Target", 1.0}}},
      {"SurfaceResponseCorrection",
       {{"Library", library_3d_path.string()},
        {"TargetInterfaces", {1, 2, 3}},
        {"UnmatchedPolicy", "Error"}}}};
  IoData maxwell_iodata_3d(maxwell_config_3d, false);
  auto maxwell_mesh_3d = mesh::ReadMesh(maxwell_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> maxwell_meshes_3d;
  maxwell_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(maxwell_mesh_3d)));
  SpaceOperator maxwell_space_3d(maxwell_iodata_3d, maxwell_meshes_3d);
  SurfaceResponseOperator maxwell_response_3d(maxwell_iodata_3d,
                                              maxwell_space_3d);
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
  const auto real_response =
      maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  CHECK(std::abs(real_response.domain_correction) > 0.0);
  REQUIRE(real_response.fabricated_surface_energy.size() == 3);
  CHECK(real_response.loop_residual < 1.0e-10);
  CHECK(real_response.kR == 0.0);
  CHECK_THAT(real_response.matched_length_fraction, WithinAbs(1.0, 1.0e-12));

  // A Maxwell trace reconstructed from E = -grad(V) must reproduce the H1 coupon trace
  // relative to the PEC. This catches a missing contour-voltage gauge even when the
  // contour-loop residual and complex-field scaling tests pass.
  mfem::ParGridFunction potential(&laplace_3d.GetH1Space().Get());
  mfem::FunctionCoefficient potential_coefficient(
      [](const mfem::Vector &x) { return x[1]; });
  potential.ProjectCoefficient(potential_coefficient);
  Vector potential_true;
  potential.GetTrueDofs(potential_true);
  mfem::Vector normal_field(3);
  normal_field = 0.0;
  normal_field[1] = -1.0;
  mfem::VectorConstantCoefficient normal_field_coefficient(normal_field);
  maxwell_field.Real().ProjectCoefficient(normal_field_coefficient);
  maxwell_field.Imag() = 0.0;

  const auto electrostatic_correction =
      response_3d.GetEnergyCorrection(potential_true);
  const auto electrostatic_surfaces =
      response_3d.GetFabricatedSurfaceEnergy(potential_true);
  const auto gradient_response =
      maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  CHECK_THAT(gradient_response.domain_correction,
             WithinRel(electrostatic_correction.domain, 1.0e-10));
  REQUIRE(gradient_response.fabricated_surface_energy.size() ==
          electrostatic_surfaces.size());
  for (const auto &[interface, energy] : electrostatic_surfaces)
  {
    CHECK_THAT(gradient_response.fabricated_surface_energy.at(interface),
               WithinRel(energy, 1.0e-10));
  }
  CHECK(gradient_response.loop_residual < 1.0e-10);

  maxwell_field.Real().ProjectCoefficient(field_coefficient);
  maxwell_field.Imag() = maxwell_field.Real();
  const auto complex_response =
      maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  CHECK_THAT(complex_response.domain_correction,
             WithinRel(2.0 * real_response.domain_correction, 1.0e-10));
  for (const auto &[interface, energy] :
       real_response.fabricated_surface_energy)
  {
    CHECK_THAT(complex_response.fabricated_surface_energy.at(interface),
               WithinRel(2.0 * energy, 1.0e-10));
  }

  auto coupled_config_3d = config_3d;
  for (auto &interface : coupled_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 7.0};
  }
  coupled_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      coupled_library_3d_path.string();
  IoData coupled_iodata_3d(coupled_config_3d, false);
  auto coupled_mesh_3d = mesh::ReadMesh(coupled_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> coupled_meshes_3d;
  coupled_meshes_3d.push_back(std::make_unique<Mesh>(std::move(coupled_mesh_3d)));
  LaplaceOperator coupled_laplace_3d(coupled_iodata_3d, coupled_meshes_3d);
  SurfaceResponseOperator coupled_response_3d(coupled_iodata_3d, coupled_laplace_3d);
  CHECK(coupled_response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size() / 2) * line_rule.GetNPoints());
  CHECK(coupled_response_3d.GetBasisSize() == 4 * coupled_response_3d.GetPatchCount());
  CHECK_THAT(coupled_response_3d.GetPatchWeight(),
             WithinRel(0.25 * physical_edge_length, 1.0e-12));
  CHECK(coupled_response_3d.HasSurfaceResponse());

  auto coupled_maxwell_config_3d = maxwell_config_3d;
  for (auto &interface :
       coupled_maxwell_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 7.0};
  }
  coupled_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      coupled_library_3d_path.string();
  IoData coupled_maxwell_iodata_3d(coupled_maxwell_config_3d, false);
  auto coupled_maxwell_mesh_3d =
      mesh::ReadMesh(coupled_maxwell_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> coupled_maxwell_meshes_3d;
  coupled_maxwell_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(coupled_maxwell_mesh_3d)));
  SpaceOperator coupled_maxwell_space_3d(coupled_maxwell_iodata_3d,
                                         coupled_maxwell_meshes_3d);
  SurfaceResponseOperator coupled_maxwell_response_3d(
      coupled_maxwell_iodata_3d, coupled_maxwell_space_3d);
  const auto &maxwell_line_rule = mfem::IntRules.Get(
      mfem::Geometry::SEGMENT, 2 * coupled_maxwell_iodata_3d.solver.order);
  CHECK(coupled_maxwell_response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size() / 2) *
            maxwell_line_rule.GetNPoints());
  GridFunction coupled_maxwell_field(coupled_maxwell_space_3d.GetNDSpace(), true);
  coupled_maxwell_field.Real().ProjectCoefficient(normal_field_coefficient);
  coupled_maxwell_field.Imag() = 0.0;
  const auto coupled_maxwell_result =
      coupled_maxwell_response_3d.GetMaxwellResponse(coupled_maxwell_field, 0.0);
  CHECK(std::abs(coupled_maxwell_result.domain_correction) > 0.0);
  CHECK(coupled_maxwell_result.fabricated_surface_energy.size() == 3);
  CHECK(coupled_maxwell_result.loop_residual < 1.0e-10);
#endif
}

}  // namespace palace
