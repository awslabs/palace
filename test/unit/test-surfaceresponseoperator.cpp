// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <memory>
#include <set>
#include <string_view>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodeoperator.hpp"
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
  const auto compact_fabricated_surface_path =
      temp.temp_dir / "fabricated-surface-compact.csv";
  const auto compact_thin_surface_path = temp.temp_dir / "thin-surface-compact.csv";
  const auto library_path = temp.temp_dir / "fabrication-process.json";
  const auto legacy_library_path = temp.temp_dir / "fabrication-process-legacy.json";
  const auto missing_layer_library_path =
      temp.temp_dir / "fabrication-process-missing-layer.json";
  const auto invalid_library_path =
      temp.temp_dir / "fabrication-process-invalid-depth.json";
  const auto impedance_library_path = temp.temp_dir / "fabrication-process-impedance.json";
  const auto legacy_impedance_library_path =
      temp.temp_dir / "fabrication-process-impedance-legacy.json";
  const auto conductivity_library_path =
      temp.temp_dir / "fabrication-process-conductivity.json";
  const auto rational_impedance_library_path =
      temp.temp_dir / "fabrication-process-rational-impedance.json";
  const auto invalid_boundary_law_library_path =
      temp.temp_dir / "fabrication-process-invalid-boundary-law.json";
  const auto library_3d_path = temp.temp_dir / "fabrication-process-3d.json";
  const auto impedance_library_3d_path =
      temp.temp_dir / "fabrication-process-impedance-3d.json";
  const auto exact_pair_library_2d_path =
      temp.temp_dir / "fabrication-process-exact-pair-2d.json";
  const auto different_pair_library_2d_path =
      temp.temp_dir / "fabrication-process-different-pair-2d.json";
  const auto interpolated_pair_library_2d_path =
      temp.temp_dir / "fabrication-process-interpolated-pair-2d.json";
  const auto parallel_cluster_library_2d_path =
      temp.temp_dir / "fabrication-process-parallel-cluster-2d.json";
  const auto impedance_parallel_cluster_library_2d_path =
      temp.temp_dir / "fabrication-process-parallel-cluster-impedance-2d.json";
  const auto parallel_cluster_points_2d_path =
      temp.temp_dir / "parallel-cluster-basis-points-2d.csv";
  const auto coupled_library_3d_path =
      temp.temp_dir / "fabrication-process-coupled-3d.json";
  const auto missing_pair_library_3d_path =
      temp.temp_dir / "fabrication-process-missing-pair-3d.json";
  const auto interpolated_coupled_library_3d_path =
      temp.temp_dir / "fabrication-process-coupled-interpolated-3d.json";
  const auto parallel_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-parallel-cluster-3d.json";
  const auto parallel_cluster_only_library_3d_path =
      temp.temp_dir / "fabrication-process-parallel-cluster-only-3d.json";
  const auto disconnected_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-disconnected-cluster-3d.json";
  const auto coupled_fabricated_path = temp.temp_dir / "coupled-fabricated.csv";
  const auto coupled_thin_path = temp.temp_dir / "coupled-thin.csv";
  const auto coupled_fabricated_surface_path =
      temp.temp_dir / "coupled-fabricated-surface.csv";
  const auto coupled_thin_surface_path = temp.temp_dir / "coupled-thin-surface.csv";
  const auto cluster_fabricated_path = temp.temp_dir / "cluster-fabricated.csv";
  const auto cluster_thin_path = temp.temp_dir / "cluster-thin.csv";
  const auto cluster_fabricated_surface_path =
      temp.temp_dir / "cluster-fabricated-surface.csv";
  const auto cluster_thin_surface_path = temp.temp_dir / "cluster-thin-surface.csv";
  const auto zero_trace_path = temp.temp_dir / "zero-trace.csv";
  const auto corner_points_path = temp.temp_dir / "corner-basis-points.csv";
  const auto corner_fabricated_path = temp.temp_dir / "corner-fabricated.csv";
  const auto corner_thin_path = temp.temp_dir / "corner-thin.csv";
  const auto corner_constrained_perturbed_fabricated_path =
      temp.temp_dir / "corner-constrained-perturbed-fabricated.csv";
  const auto corner_constrained_perturbed_thin_path =
      temp.temp_dir / "corner-constrained-perturbed-thin.csv";
  const auto corner_fabricated_surface_path =
      temp.temp_dir / "corner-fabricated-surface.csv";
  const auto corner_thin_surface_path = temp.temp_dir / "corner-thin-surface.csv";
  const auto convex_library_3d_path = temp.temp_dir / "fabrication-process-convex-3d.json";
  const auto finite_impedance_convex_library_3d_path =
      temp.temp_dir / "fabrication-process-convex-finite-impedance-3d.json";
  const auto concave_library_3d_path =
      temp.temp_dir / "fabrication-process-concave-3d.json";
  const auto strip_library_3d_path = temp.temp_dir / "fabrication-process-strip-3d.json";
  const auto rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-3d.json";
  const auto constrained_perturbed_rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-constrained-perturbed-3d.json";
  const auto finite_impedance_rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-finite-impedance-3d.json";
  const auto rounded_concave_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-concave-3d.json";
  const auto interpolated_rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-interpolated-3d.json";
  const auto unqualified_interpolated_rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-interpolated-unqualified-3d.json";
  const auto endpoint_library_3d_path =
      temp.temp_dir / "fabrication-process-endpoint-3d.json";
  const auto junction_library_3d_path =
      temp.temp_dir / "fabrication-process-junction-3d.json";
  const auto spatial_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-3d.json";
  const auto cross_layer_spatial_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-cross-layer-3d.json";
  const auto incomplete_cross_layer_spatial_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-cross-layer-incomplete-3d.json";
  const auto spatial_cluster_position_mismatch_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-position-mismatch-3d.json";
  const auto spatial_cluster_orientation_mismatch_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-orientation-mismatch-3d.json";
  const auto spatial_cluster_interval_mismatch_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-interval-mismatch-3d.json";
  const auto spatial_cluster_extra_edge_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-extra-edge-3d.json";
  const auto spatial_cluster_impedance_mismatch_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-impedance-mismatch-3d.json";
  const auto spatial_cluster_mixed_impedance_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-mixed-impedance-3d.json";
  const auto spatial_cluster_parameter_mismatch_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-parameter-mismatch-3d.json";
  const auto spatial_cluster_points_path =
      temp.temp_dir / "spatial-cluster-basis-points.csv";
  const auto cross_layer_fabricated_surface_path =
      temp.temp_dir / "cross-layer-fabricated-surface.csv";
  const auto cross_layer_thin_surface_path = temp.temp_dir / "cross-layer-thin-surface.csv";
  const json impedance_law = {{"Type", "Impedance"}, {"Ls", 1.0e-13}};
  const json second_impedance_law = {{"Type", "Impedance"}, {"Ls", 2.0e-13}};
  const json conductivity_law = {{"Type", "Conductivity"},
                                 {"Conductivity", 5.8e7},
                                 {"Permeability", 1.2},
                                 {"Thickness", 1.0e-7},
                                 {"External", true}};
  const json rational_impedance_law = {{"Type", "RationalImpedance"},
                                       {"Numerator", {5.0e-8, 0.0}},
                                       {"Denominator", {5.0e-20, 1.0e-9, 50.0}}};
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
    auto write_compact_surface_matrix = [](const auto &path, const auto &matrix)
    {
      std::ofstream output(path);
      output << "interface,edge,basis_i,basis_j,Q_total_ij (J)\n";
      for (std::size_t i = 0; i < matrix.size(); i++)
      {
        for (std::size_t j = i; j < matrix.size(); j++)
        {
          output << "1,1," << i + 1 << "," << j + 1 << "," << matrix[i][j] << "\n";
        }
      }
    };
    write_domain_matrix(fabricated_path, fabricated);
    write_domain_matrix(thin_path, thin);
    write_surface_matrix(fabricated_surface_path, fabricated);
    write_surface_matrix(thin_surface_path, thin);
    write_compact_surface_matrix(compact_fabricated_surface_path, fabricated);
    write_compact_surface_matrix(compact_thin_surface_path, thin);
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
    auto impedance_library = library;
    impedance_library["Name"] = "unit-test-process-impedance";
    for (auto &model : impedance_library["Models"])
    {
      model["BoundaryCondition"] = impedance_law;
    }
    std::ofstream impedance_output(impedance_library_path);
    impedance_output << impedance_library.dump(2) << "\n";
    auto legacy_impedance_library = impedance_library;
    legacy_impedance_library["Name"] = "unit-test-process-impedance-legacy";
    for (auto &model : legacy_impedance_library["Models"])
    {
      model["BoundaryCondition"] = "Impedance";
    }
    std::ofstream legacy_impedance_output(legacy_impedance_library_path);
    legacy_impedance_output << legacy_impedance_library.dump(2) << "\n";
    auto conductivity_library = library;
    conductivity_library["Name"] = "unit-test-process-conductivity";
    for (auto &model : conductivity_library["Models"])
    {
      model["BoundaryCondition"] = conductivity_law;
    }
    std::ofstream conductivity_output(conductivity_library_path);
    conductivity_output << conductivity_library.dump(2) << "\n";
    auto rational_impedance_library = library;
    rational_impedance_library["Name"] = "unit-test-process-rational-impedance";
    for (auto &model : rational_impedance_library["Models"])
    {
      model["BoundaryCondition"] = {{"Type", "RationalImpedance"},
                                    {"Numerator", {1.0e-7, 0.0}},
                                    {"Denominator", {1.0e-19, 2.0e-9, 100.0}}};
    }
    std::ofstream rational_impedance_output(rational_impedance_library_path);
    rational_impedance_output << rational_impedance_library.dump(2) << "\n";
    auto invalid_boundary_law_library = impedance_library;
    invalid_boundary_law_library["Name"] = "unit-test-process-invalid-boundary-law";
    invalid_boundary_law_library["Models"][0]["BoundaryCondition"]["Inductance"] = 1.0e-13;
    std::ofstream invalid_boundary_law_output(invalid_boundary_law_library_path);
    invalid_boundary_law_output << invalid_boundary_law_library.dump(2) << "\n";
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
    auto different_pair_library_2d = exact_pair_library_2d;
    different_pair_library_2d["Name"] = "unit-test-different-pair-2d";
    different_pair_library_2d["MatchingRadius"] = 0.25;
    different_pair_library_2d["Models"][0]["Name"] = "different-gap-0.4";
    different_pair_library_2d["Models"][0]["Topology"] = "DifferentConductorGap";
    different_pair_library_2d["Models"][0]["Separation"] = 0.4;
    std::ofstream different_pair_output_2d(different_pair_library_2d_path);
    different_pair_output_2d << different_pair_library_2d.dump(2) << "\n";
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
    auto impedance_library_3d = library_3d;
    impedance_library_3d["Name"] = "unit-test-process-impedance-3d";
    for (auto &model : impedance_library_3d["Models"])
    {
      model["BoundaryCondition"] = impedance_law;
    }
    std::ofstream impedance_output_3d(impedance_library_3d_path);
    impedance_output_3d << impedance_library_3d.dump(2) << "\n";
    auto coupled_library_3d = library_3d;
    coupled_library_3d["Version"] = 3;
    coupled_library_3d["Name"] = "unit-test-process-coupled-3d";
    coupled_library_3d["MatchingRadius"] = 7.0;
    auto missing_pair_library_3d = coupled_library_3d;
    missing_pair_library_3d["Name"] = "unit-test-process-missing-pair-3d";
    std::ofstream missing_pair_output_3d(missing_pair_library_3d_path);
    missing_pair_output_3d << missing_pair_library_3d.dump(2) << "\n";
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

    auto parallel_cluster_library_3d = coupled_library_3d;
    parallel_cluster_library_3d["Name"] = "unit-test-process-parallel-cluster-3d";
    parallel_cluster_library_3d["MatchingRadius"] = 11.0;
    auto strip_20_model = library_3d["Models"][0];
    strip_20_model["Name"] = "trace-strip-20um";
    strip_20_model["Topology"] = "SameConductorStrip";
    strip_20_model["Separation"] = 20.0;
    strip_20_model["SeparationTolerance"] = 1.0e-6;
    auto parallel_cluster_model = coupled_library_3d["Models"][1];
    parallel_cluster_model["Name"] = "cpw-four-edge-cluster";
    parallel_cluster_model["Topology"] = "ParallelEdgeCluster";
    parallel_cluster_model.erase("Separation");
    parallel_cluster_model.erase("SeparationTolerance");
    parallel_cluster_model["EdgeOffsetTolerance"] = 1.0e-6;
    parallel_cluster_model["Edges"] = {
        {{"Offset", 0.0}, {"GapDirection", 1}, {"Conductor", 1}},
        {{"Offset", 12.0}, {"GapDirection", -1}, {"Conductor", 2}},
        {{"Offset", 32.0}, {"GapDirection", 1}, {"Conductor", 2}},
        {{"Offset", 44.0}, {"GapDirection", -1}, {"Conductor", 1}}};
    std::array<std::array<double, 6>, 6> cluster_fabricated{};
    std::array<std::array<double, 6>, 6> cluster_thin{};
    for (std::size_t i = 0; i < cluster_fabricated.size(); i++)
    {
      cluster_fabricated[i][i] = (i >= 4 ? 3.0 : 1.0) * 1.0e-12;
      cluster_thin[i][i] = 1.0e-12;
    }
    write_domain_matrix(cluster_fabricated_path, cluster_fabricated);
    write_domain_matrix(cluster_thin_path, cluster_thin);
    write_surface_matrix(cluster_fabricated_surface_path, cluster_fabricated);
    write_surface_matrix(cluster_thin_surface_path, cluster_thin);
    auto disconnected_ground_cluster_model = parallel_cluster_model;
    disconnected_ground_cluster_model["Name"] = "cpw-three-conductor-four-edge-cluster";
    disconnected_ground_cluster_model["ConductorReferences"] = {
        {0.0, 0.0, 0.0}, {12.0, 0.0, 0.0}, {44.0, 0.0, 0.0}};
    disconnected_ground_cluster_model["Edges"][3]["Conductor"] = 3;
    disconnected_ground_cluster_model["OpenContourPaths"] = {
        {{"Indices", {1, 2}}, {"StartConductor", 1}, {"EndConductor", 2}},
        {{"Indices", {3, 4}}, {"StartConductor", 2}, {"EndConductor", 3}}};
    disconnected_ground_cluster_model["FabricatedMatrix"] =
        cluster_fabricated_path.string();
    disconnected_ground_cluster_model["ThinMatrix"] = cluster_thin_path.string();
    disconnected_ground_cluster_model["FabricatedSurfaceMatrix"] =
        cluster_fabricated_surface_path.string();
    disconnected_ground_cluster_model["ThinSurfaceMatrix"] =
        cluster_thin_surface_path.string();
    parallel_cluster_library_3d["Models"].push_back(std::move(strip_20_model));
    parallel_cluster_library_3d["Models"].push_back(std::move(parallel_cluster_model));
    parallel_cluster_library_3d["Models"].push_back(
        std::move(disconnected_ground_cluster_model));
    std::ofstream parallel_cluster_output_3d(parallel_cluster_library_3d_path);
    parallel_cluster_output_3d << parallel_cluster_library_3d.dump(2) << "\n";
    auto parallel_cluster_only_library_3d = parallel_cluster_library_3d;
    parallel_cluster_only_library_3d["Name"] = "unit-test-process-parallel-cluster-only-3d";
    parallel_cluster_only_library_3d["Models"] = json::array();
    for (const auto &model : parallel_cluster_library_3d["Models"])
    {
      if (model.value("Topology", "") == "ParallelEdgeCluster")
      {
        parallel_cluster_only_library_3d["Models"].push_back(model);
      }
    }
    std::ofstream parallel_cluster_only_output_3d(parallel_cluster_only_library_3d_path);
    parallel_cluster_only_output_3d << parallel_cluster_only_library_3d.dump(2) << "\n";
    auto disconnected_cluster_library_3d = parallel_cluster_library_3d;
    disconnected_cluster_library_3d["Name"] = "unit-test-process-disconnected-cluster-3d";
    disconnected_cluster_library_3d["Models"].back()["OpenContourPaths"].erase(1);
    std::ofstream disconnected_cluster_output_3d(disconnected_cluster_library_3d_path);
    disconnected_cluster_output_3d << disconnected_cluster_library_3d.dump(2) << "\n";

    {
      std::ofstream output(parallel_cluster_points_2d_path);
      output << "x,y,z\n"
             << "0.05,-0.05,0.0\n"
             << "0.15,-0.05,0.0\n"
             << "0.45,-0.05,0.0\n"
             << "0.55,-0.05,0.0\n";
    }
    const json cpw_cluster_edges = {
        {{"Offset", 0.0}, {"GapDirection", 1}, {"Conductor", 1}},
        {{"Offset", 0.2}, {"GapDirection", -1}, {"Conductor", 2}},
        {{"Offset", 0.4}, {"GapDirection", 1}, {"Conductor", 2}},
        {{"Offset", 0.6}, {"GapDirection", -1}, {"Conductor", 1}}};
    auto two_conductor_cluster_model = parallel_cluster_library_3d["Models"][1];
    two_conductor_cluster_model["Name"] = "cpw-two-conductor-four-edge-cluster";
    two_conductor_cluster_model["Topology"] = "ParallelEdgeCluster";
    two_conductor_cluster_model.erase("Separation");
    two_conductor_cluster_model.erase("SeparationTolerance");
    two_conductor_cluster_model.erase("CouponDepth");
    two_conductor_cluster_model["BasisPoints"] = parallel_cluster_points_2d_path.string();
    two_conductor_cluster_model["ConductorReferences"] = {{0.0, 0.0, 0.0}, {0.2, 0.0, 0.0}};
    two_conductor_cluster_model["Edges"] = cpw_cluster_edges;
    two_conductor_cluster_model["EdgeOffsetTolerance"] = 1.0e-8;
    two_conductor_cluster_model["OpenContourPaths"] = {
        {{"Indices", {1, 2}}, {"StartConductor", 1}, {"EndConductor", 2}},
        {{"Indices", {3, 4}}, {"StartConductor", 2}, {"EndConductor", 1}}};
    auto three_conductor_cluster_model = parallel_cluster_library_3d["Models"].back();
    three_conductor_cluster_model["Name"] = "cpw-three-conductor-four-edge-cluster-2d";
    three_conductor_cluster_model.erase("CouponDepth");
    three_conductor_cluster_model["BasisPoints"] = parallel_cluster_points_2d_path.string();
    three_conductor_cluster_model["ConductorReferences"] = {
        {0.0, 0.0, 0.0}, {0.2, 0.0, 0.0}, {0.6, 0.0, 0.0}};
    three_conductor_cluster_model["Edges"] = cpw_cluster_edges;
    three_conductor_cluster_model["Edges"][3]["Conductor"] = 3;
    three_conductor_cluster_model["EdgeOffsetTolerance"] = 1.0e-8;
    auto parallel_cluster_library_2d = library;
    parallel_cluster_library_2d["Name"] = "unit-test-process-parallel-cluster-2d";
    parallel_cluster_library_2d["MatchingRadius"] = 0.25;
    parallel_cluster_library_2d["Models"] = {two_conductor_cluster_model,
                                             three_conductor_cluster_model};
    std::ofstream parallel_cluster_output_2d(parallel_cluster_library_2d_path);
    parallel_cluster_output_2d << parallel_cluster_library_2d.dump(2) << "\n";
    auto impedance_parallel_cluster_library_2d = parallel_cluster_library_2d;
    impedance_parallel_cluster_library_2d["Name"] =
        "unit-test-process-parallel-cluster-impedance-2d";
    for (auto &model : impedance_parallel_cluster_library_2d["Models"])
    {
      model["BoundaryCondition"] = impedance_law;
    }
    std::ofstream impedance_parallel_cluster_output_2d(
        impedance_parallel_cluster_library_2d_path);
    impedance_parallel_cluster_output_2d << impedance_parallel_cluster_library_2d.dump(2)
                                         << "\n";

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
    auto corner_constrained_perturbed_fabricated = corner_fabricated;
    auto corner_constrained_perturbed_thin = corner_thin;
    for (const std::size_t i : {4, 5, 6, 7})
    {
      corner_constrained_perturbed_fabricated[i][i] = 8.0e-12;
      corner_constrained_perturbed_thin[i][i] = 0.4e-12;
      for (std::size_t j = 0; j < corner_fabricated.size(); j++)
      {
        if (j == i)
        {
          continue;
        }
        const double coupling =
            1.0 / (1.0 + std::abs(static_cast<int>(i) - static_cast<int>(j)));
        corner_constrained_perturbed_fabricated[i][j] =
            corner_constrained_perturbed_fabricated[j][i] = 0.2e-12 * coupling;
        corner_constrained_perturbed_thin[i][j] = corner_constrained_perturbed_thin[j][i] =
            0.15e-12 * coupling;
      }
    }
    write_domain_matrix(corner_constrained_perturbed_fabricated_path,
                        corner_constrained_perturbed_fabricated);
    write_domain_matrix(corner_constrained_perturbed_thin_path,
                        corner_constrained_perturbed_thin);

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

    auto finite_impedance_convex_library_3d = convex_library_3d;
    finite_impedance_convex_library_3d["Name"] =
        "unit-test-process-convex-finite-impedance-3d";
    for (auto &model : finite_impedance_convex_library_3d["Models"])
    {
      model["BoundaryCondition"] = impedance_law;
    }
    finite_impedance_convex_library_3d["Models"][1]["Name"] =
        "convex-corner-90-finite-impedance";
    finite_impedance_convex_library_3d["Models"][1]["BoundaryCondition"] = impedance_law;
    finite_impedance_convex_library_3d["Models"][1]["Reference"] = {0.0, 0.0, 0.0};
    std::ofstream finite_impedance_convex_output_3d(
        finite_impedance_convex_library_3d_path);
    finite_impedance_convex_output_3d << finite_impedance_convex_library_3d.dump(2) << "\n";

    {
      std::ofstream output(spatial_cluster_points_path);
      output << "x,y,z\n"
             << "0.025,-0.025,0.02\n"
             << "0.050,-0.050,0.02\n"
             << "0.075,-0.075,-0.02\n"
             << "0.100,-0.100,-0.02\n";
    }
    auto spatial_cluster_library_3d = convex_library_3d;
    spatial_cluster_library_3d["Name"] = "unit-test-process-spatial-cluster-3d";
    auto spatial_cluster_model = coupled_library_3d["Models"][1];
    spatial_cluster_model["Name"] = "offset-corner-pair";
    spatial_cluster_model["Topology"] = "SpatialEdgeCluster";
    spatial_cluster_model.erase("Separation");
    spatial_cluster_model.erase("SeparationTolerance");
    spatial_cluster_model.erase("CouponDepth");
    spatial_cluster_model["BasisPoints"] = spatial_cluster_points_path.string();
    spatial_cluster_model["EdgePositionTolerance"] = 1.0e-6;
    spatial_cluster_model["EdgeAngleTolerance"] = 1.0e-6;
    spatial_cluster_model["ConductorReferences"] = {{0.0, 0.0, 0.0}, {0.125, -0.125, 0.0}};
    spatial_cluster_model["OpenContourPaths"] = {
        {{"Indices", {1, 2}}, {"StartConductor", 1}, {"EndConductor", 2}},
        {{"Indices", {4, 3}}, {"StartConductor", 1}, {"EndConductor", 2}}};
    spatial_cluster_model["Edges"] = {{{"Point", {0.0, 0.0, 0.0}},
                                       {"GapDirection", {0.0, -1.0, 0.0}},
                                       {"ProcessNormal", {0.0, 0.0, 1.0}},
                                       {"Interval", {0.0, 0.2}},
                                       {"Conductor", 1},
                                       {"BoundaryCondition", "PEC"}},
                                      {{"Point", {0.0, 0.0, 0.0}},
                                       {"GapDirection", {1.0, 0.0, 0.0}},
                                       {"ProcessNormal", {0.0, 0.0, 1.0}},
                                       {"Interval", {-0.2, 0.0}},
                                       {"Conductor", 1},
                                       {"BoundaryCondition", "PEC"}},
                                      {{"Point", {0.125, -0.125, 0.0}},
                                       {"GapDirection", {-1.0, 0.0, 0.0}},
                                       {"ProcessNormal", {0.0, 0.0, 1.0}},
                                       {"Interval", {-0.2, 0.0}},
                                       {"Conductor", 2},
                                       {"BoundaryCondition", "PEC"}},
                                      {{"Point", {0.125, -0.125, 0.0}},
                                       {"GapDirection", {0.0, 1.0, 0.0}},
                                       {"ProcessNormal", {0.0, 0.0, 1.0}},
                                       {"Interval", {0.0, 0.2}},
                                       {"Conductor", 2},
                                       {"BoundaryCondition", "PEC"}}};
    spatial_cluster_library_3d["Models"].push_back(std::move(spatial_cluster_model));
    std::ofstream spatial_cluster_output_3d(spatial_cluster_library_3d_path);
    spatial_cluster_output_3d << spatial_cluster_library_3d.dump(2) << "\n";

    auto write_cross_layer_surface_matrix = [](const auto &path, const auto &matrix)
    {
      std::ofstream output(path);
      output << "interface,edge,basis_i,basis_j,Q_total_ij (J)\n";
      for (int interface = 1; interface <= 2; interface++)
      {
        for (std::size_t i = 0; i < matrix.size(); i++)
        {
          for (std::size_t j = i; j < matrix.size(); j++)
          {
            output << interface << ",1," << i + 1 << "," << j + 1 << ","
                   << 0.5 * matrix[i][j] << "\n";
          }
        }
      }
    };
    write_cross_layer_surface_matrix(cross_layer_fabricated_surface_path,
                                     coupled_fabricated);
    write_cross_layer_surface_matrix(cross_layer_thin_surface_path, coupled_thin);
    auto cross_layer_spatial_cluster_library_3d = spatial_cluster_library_3d;
    cross_layer_spatial_cluster_library_3d["Name"] =
        "unit-test-process-spatial-cluster-cross-layer-3d";
    auto &cross_layer_model = cross_layer_spatial_cluster_library_3d["Models"].back();
    cross_layer_model["Name"] = "offset-corner-pair-cross-layer";
    cross_layer_model["FabricatedSurfaceMatrix"] =
        cross_layer_fabricated_surface_path.string();
    cross_layer_model["ThinSurfaceMatrix"] = cross_layer_thin_surface_path.string();
    cross_layer_model["Interfaces"] = {{{"Slot", 1}, {"Type", "SA"}, {"Coupon", 1}},
                                       {{"Slot", 2}, {"Type", "SA"}, {"Coupon", 2}}};
    for (std::size_t edge = 0; edge < cross_layer_model["Edges"].size(); edge++)
    {
      cross_layer_model["Edges"][edge]["InterfaceSlot"] = edge < 2 ? 1 : 2;
    }
    std::ofstream cross_layer_spatial_cluster_output_3d(
        cross_layer_spatial_cluster_library_3d_path);
    cross_layer_spatial_cluster_output_3d << cross_layer_spatial_cluster_library_3d.dump(2)
                                          << "\n";

    auto incomplete_cross_layer_library = cross_layer_spatial_cluster_library_3d;
    incomplete_cross_layer_library["Name"] =
        "unit-test-process-spatial-cluster-cross-layer-incomplete-3d";
    incomplete_cross_layer_library["Models"].back()["Interfaces"].erase(1);
    std::ofstream incomplete_cross_layer_output_3d(
        incomplete_cross_layer_spatial_cluster_library_3d_path);
    incomplete_cross_layer_output_3d << incomplete_cross_layer_library.dump(2) << "\n";

    auto position_mismatch_library = spatial_cluster_library_3d;
    position_mismatch_library["Name"] =
        "unit-test-process-spatial-cluster-position-mismatch-3d";
    position_mismatch_library["Models"].back()["Edges"][3]["Point"] =
        std::array<double, 3>{0.135, -0.125, 0.0};
    std::ofstream position_mismatch_output_3d(
        spatial_cluster_position_mismatch_library_3d_path);
    position_mismatch_output_3d << position_mismatch_library.dump(2) << "\n";

    auto orientation_mismatch_library = spatial_cluster_library_3d;
    orientation_mismatch_library["Name"] =
        "unit-test-process-spatial-cluster-orientation-mismatch-3d";
    orientation_mismatch_library["Models"].back()["Edges"][3]["GapDirection"] =
        std::array<double, 3>{0.6, 0.8, 0.0};
    std::ofstream orientation_mismatch_output_3d(
        spatial_cluster_orientation_mismatch_library_3d_path);
    orientation_mismatch_output_3d << orientation_mismatch_library.dump(2) << "\n";

    auto interval_mismatch_library = spatial_cluster_library_3d;
    interval_mismatch_library["Name"] =
        "unit-test-process-spatial-cluster-interval-mismatch-3d";
    interval_mismatch_library["Models"].back()["Edges"][3]["Interval"] =
        std::array<double, 2>{0.0, 0.15};
    std::ofstream interval_mismatch_output_3d(
        spatial_cluster_interval_mismatch_library_3d_path);
    interval_mismatch_output_3d << interval_mismatch_library.dump(2) << "\n";

    auto extra_edge_library = spatial_cluster_library_3d;
    extra_edge_library["Name"] = "unit-test-process-spatial-cluster-extra-edge-3d";
    extra_edge_library["Models"].back()["Edges"].erase(3);
    std::ofstream extra_edge_output_3d(spatial_cluster_extra_edge_library_3d_path);
    extra_edge_output_3d << extra_edge_library.dump(2) << "\n";

    auto impedance_mismatch_library = spatial_cluster_library_3d;
    impedance_mismatch_library["Name"] =
        "unit-test-process-spatial-cluster-impedance-mismatch-3d";
    for (auto &edge : impedance_mismatch_library["Models"].back()["Edges"])
    {
      edge["BoundaryCondition"] = "Impedance";
    }
    std::ofstream impedance_mismatch_output_3d(
        spatial_cluster_impedance_mismatch_library_3d_path);
    impedance_mismatch_output_3d << impedance_mismatch_library.dump(2) << "\n";

    auto mixed_impedance_library = spatial_cluster_library_3d;
    mixed_impedance_library["Name"] =
        "unit-test-process-spatial-cluster-mixed-impedance-3d";
    for (std::size_t edge = 2; edge < mixed_impedance_library["Models"][2]["Edges"].size();
         edge++)
    {
      mixed_impedance_library["Models"][2]["Edges"][edge]["BoundaryCondition"] =
          second_impedance_law;
    }
    auto mixed_impedance_isolated_model = mixed_impedance_library["Models"][0];
    mixed_impedance_isolated_model["Name"] = "isolated-impedance-ls2";
    mixed_impedance_isolated_model["BoundaryCondition"] = second_impedance_law;
    mixed_impedance_library["Models"].push_back(std::move(mixed_impedance_isolated_model));
    auto mixed_impedance_corner_model = mixed_impedance_library["Models"][1];
    mixed_impedance_corner_model["Name"] = "convex-corner-90-impedance-ls2";
    mixed_impedance_corner_model["BoundaryCondition"] = second_impedance_law;
    mixed_impedance_corner_model["Reference"] = {0.0, 0.0, 0.0};
    mixed_impedance_library["Models"].push_back(std::move(mixed_impedance_corner_model));
    std::ofstream mixed_impedance_output_3d(
        spatial_cluster_mixed_impedance_library_3d_path);
    mixed_impedance_output_3d << mixed_impedance_library.dump(2) << "\n";

    auto parameter_mismatch_library = mixed_impedance_library;
    parameter_mismatch_library["Name"] =
        "unit-test-process-spatial-cluster-parameter-mismatch-3d";
    parameter_mismatch_library["Models"][2]["Edges"][3]["BoundaryCondition"] = {
        {"Type", "Impedance"}, {"Ls", 3.0e-13}};
    std::ofstream parameter_mismatch_output_3d(
        spatial_cluster_parameter_mismatch_library_3d_path);
    parameter_mismatch_output_3d << parameter_mismatch_library.dump(2) << "\n";

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

    auto constrained_perturbed_rounded_library_3d = rounded_library_3d;
    constrained_perturbed_rounded_library_3d["Name"] =
        "unit-test-process-rounded-constrained-perturbed-3d";
    constrained_perturbed_rounded_library_3d["Models"][1]["FabricatedMatrix"] =
        corner_constrained_perturbed_fabricated_path.string();
    constrained_perturbed_rounded_library_3d["Models"][1]["ThinMatrix"] =
        corner_constrained_perturbed_thin_path.string();
    std::ofstream constrained_perturbed_rounded_output_3d(
        constrained_perturbed_rounded_library_3d_path);
    constrained_perturbed_rounded_output_3d
        << constrained_perturbed_rounded_library_3d.dump(2) << "\n";

    auto finite_impedance_rounded_library_3d = rounded_library_3d;
    finite_impedance_rounded_library_3d["Name"] =
        "unit-test-process-rounded-finite-impedance-3d";
    for (auto &model : finite_impedance_rounded_library_3d["Models"])
    {
      model["BoundaryCondition"] = "Impedance";
    }
    finite_impedance_rounded_library_3d["Models"][1]["Name"] =
        "convex-corner-90-r0.125-finite-impedance";
    finite_impedance_rounded_library_3d["Models"][1]["BoundaryCondition"] = "Impedance";
    const double rounded_reference = 0.125 * (1.0 - std::sqrt(0.5));
    finite_impedance_rounded_library_3d["Models"][1]["Reference"] = {
        rounded_reference, rounded_reference, 0.0};
    finite_impedance_rounded_library_3d["Models"][1].erase("ZeroTraceIndices");
    std::ofstream finite_impedance_rounded_output_3d(
        finite_impedance_rounded_library_3d_path);
    finite_impedance_rounded_output_3d << finite_impedance_rounded_library_3d.dump(2)
                                       << "\n";

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
    std::ofstream unqualified_interpolated_rounded_output_3d(
        unqualified_interpolated_rounded_library_3d_path);
    unqualified_interpolated_rounded_output_3d << interpolated_rounded_library_3d.dump(2)
                                               << "\n";
    interpolated_rounded_library_3d["CornerRadiusInterpolation"] = {
        {{"LowerModel", "convex-corner-90-r0.1"},
         {"UpperModel", "convex-corner-90-r0.15"},
         {"Qualification",
          {{"Method", "HeldOutCoupon"}, {"Passed", true}, {"HeldoutRadius", 0.125}}}}};
    std::ofstream interpolated_rounded_output_3d(interpolated_rounded_library_3d_path);
    interpolated_rounded_output_3d << interpolated_rounded_library_3d.dump(2) << "\n";

    auto endpoint_library_3d = library_3d;
    endpoint_library_3d["Name"] = "unit-test-process-endpoint-3d";
    auto endpoint_model = corner_model;
    endpoint_model["Name"] = "endpoint";
    endpoint_model["Topology"] = "Endpoint";
    endpoint_model.erase("Angle");
    endpoint_model.erase("AngleTolerance");
    endpoint_library_3d["Models"].push_back(std::move(endpoint_model));
    std::ofstream endpoint_output_3d(endpoint_library_3d_path);
    endpoint_output_3d << endpoint_library_3d.dump(2) << "\n";

    auto junction_library_3d = convex_library_3d;
    junction_library_3d["Name"] = "unit-test-process-junction-3d";
    auto junction_model = corner_model;
    junction_model["Name"] = "junction-4x90";
    junction_model["Topology"] = "Junction";
    junction_model.erase("Angle");
    junction_model.erase("AngleTolerance");
    junction_model["ArmAngles"] = {0.0, 90.0, 180.0, 270.0};
    junction_model["ArmAngleTolerance"] = 1.0e-6;
    junction_library_3d["Models"].push_back(junction_model);
    junction_model["Name"] = "junction-4x90-finite-impedance";
    junction_model["BoundaryCondition"] = "Impedance";
    junction_library_3d["Models"].push_back(std::move(junction_model));
    auto impedance_isolated_model = junction_library_3d["Models"][0];
    impedance_isolated_model["Name"] = "isolated-impedance";
    impedance_isolated_model["BoundaryCondition"] = "Impedance";
    junction_library_3d["Models"].push_back(std::move(impedance_isolated_model));
    std::ofstream junction_output_3d(junction_library_3d_path);
    junction_output_3d << junction_library_3d.dump(2) << "\n";
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
  auto compact_config = config;
  auto &compact_model =
      compact_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Models"][0];
  compact_model["FabricatedSurfaceMatrix"] = compact_fabricated_surface_path.string();
  compact_model["ThinSurfaceMatrix"] = compact_thin_surface_path.string();
  IoData compact_iodata(compact_config, false);
  LaplaceOperator compact_laplace_op(compact_iodata, meshes);
  SurfaceResponseOperator compact_response(compact_iodata, compact_laplace_op);

  const int size = laplace_op.GetH1Space().GetTrueVSize();
  Vector x(size), y(size), Cx, Cy, Ctx, compact_Cx;
  for (int i = 0; i < size; i++)
  {
    x(i) = 0.17 + 0.03 * (i + 1) * (Mpi::Rank(Mpi::World()) + 1);
    y(i) = -0.11 + 0.02 * (i + 2) * (Mpi::Rank(Mpi::World()) + 2);
  }
  response.Mult(x, Cx);
  response.Mult(y, Cy);
  response.MultTranspose(x, Ctx);
  compact_response.Mult(x, compact_Cx);

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
  compact_Cx -= Cx;
  double compact_error = compact_Cx * compact_Cx;
  Mpi::GlobalSum(1, &compact_error, Mpi::World());
  CHECK(compact_error == 0.0);

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
  const auto compact_fabricated_energy = compact_response.GetFabricatedSurfaceEnergy(x);
  CHECK_THAT(compact_fabricated_energy.at(4), WithinRel(fabricated_energy.at(4), 1.0e-12));
  const auto local_electrostatic_response = response.GetElectrostaticResponse(x);
  CHECK_THAT(local_electrostatic_response.domain_correction,
             WithinRel(energy.domain, 1.0e-12));
  CHECK_THAT(local_electrostatic_response.fabricated_surface_energy.at(4),
             WithinRel(fabricated_energy.at(4), 1.0e-12));
  CHECK(std::isfinite(local_electrostatic_response.domain_correction_fixed_flux));
  CHECK(local_electrostatic_response.fabricated_surface_energy_fixed_flux.at(4) > 0.0);
  CHECK(local_electrostatic_response.maximum_trace_closure_spread > 0.0);
  // One patch and one mapped interface have no averaging: the response-weighted local
  // spread is exactly the interface aggregate, and all response weight either passes or
  // fails the local 5% limit.
  CHECK_THAT(local_electrostatic_response.response_weighted_trace_closure_spread,
             WithinRel(local_electrostatic_response.maximum_trace_closure_spread, 1.0e-12));
  CHECK_THAT(local_electrostatic_response.trace_closure_response_failure_fraction,
             WithinAbs(local_electrostatic_response.maximum_trace_closure_spread > 0.05
                           ? 1.0
                           : 0.0,
                       1.0e-12));
  CHECK(local_electrostatic_response.confident ==
        (local_electrostatic_response.maximum_trace_closure_spread <= 0.05));

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

  const auto requirements_path = temp.temp_dir / "surface-response-requirements.json";
  WriteSurfaceResponseRequirements(automatic_iodata, *automatic_meshes.back(),
                                   requirements_path.string());
  std::ifstream requirements_input(requirements_path);
  REQUIRE(requirements_input);
  const json requirements = json::parse(requirements_input);
  CHECK(requirements["Version"] == 1);
  CHECK(requirements["Complete"]);
  CHECK(requirements["MeshDimension"] == 2);
  CHECK_FALSE(requirements["Maxwell"]);
  CHECK(requirements["Summary"]["Counts"]["Exact"] == 2);
  CHECK(requirements["Summary"]["Counts"]["Missing"] == 0);
  REQUIRE(requirements["Requirements"].size() == 1);
  CHECK(requirements["Requirements"][0]["Topology"] == "IsolatedEdge");
  CHECK(requirements["Requirements"][0]["Status"] == "Exact");
  CHECK(requirements["Requirements"][0]["Count"] == 2);
  CHECK(requirements["Requirements"][0]["Interfaces"][0]["Target"] == 4);

  std::ifstream empty_library_input(library_path);
  REQUIRE(empty_library_input);
  auto empty_library = json::parse(empty_library_input);
  empty_library["Name"] = "unit-test-empty-process";
  empty_library["Models"] = json::array();
  const auto empty_library_path = temp.temp_dir / "surface-process-empty.json";
  std::ofstream empty_library_output(empty_library_path);
  empty_library_output << empty_library.dump(2) << "\n";
  empty_library_output.close();
  auto empty_library_config = automatic_config;
  empty_library_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      empty_library_path.string();
  IoData empty_library_iodata(empty_library_config, false);
  empty_library_iodata.boundaries.cracked_attributes.insert(9);
  empty_library_iodata.boundaries.cracked_attributes.insert(10);
  const auto empty_requirements_path =
      temp.temp_dir / "surface-response-requirements-empty.json";
  WriteSurfaceResponseRequirements(empty_library_iodata, *automatic_meshes.back(),
                                   empty_requirements_path.string());
  std::ifstream empty_requirements_input(empty_requirements_path);
  REQUIRE(empty_requirements_input);
  const json empty_requirements = json::parse(empty_requirements_input);
  CHECK_FALSE(empty_requirements["Complete"]);
  CHECK(empty_requirements["Summary"]["Counts"]["Exact"] == 0);
  CHECK(empty_requirements["Summary"]["Counts"]["Missing"] == 2);
  REQUIRE(empty_requirements["Requirements"].size() == 1);
  CHECK(empty_requirements["Requirements"][0]["Topology"] == "IsolatedEdge");
  CHECK(empty_requirements["Requirements"][0]["Status"] == "Missing");
  CHECK(empty_requirements["Requirements"][0]["Count"] == 2);
  CHECK_THROWS(SurfaceResponseOperator(empty_library_iodata, automatic_laplace));

  auto boundary_mode_config = automatic_config;
  boundary_mode_config["Problem"]["Type"] = "BoundaryMode";
  boundary_mode_config["Boundaries"].erase("Ground");
  boundary_mode_config["Boundaries"].erase("Terminal");
  boundary_mode_config["Boundaries"]["PEC"] = {{"Attributes", {9, 10}}};
  boundary_mode_config["Solver"] = {
      {"Order", 1},
      {"BoundaryMode", {{"Freq", 5.0}}},
      {"SurfaceResponseCorrection",
       {{"Library", library_path.string()}, {"UnmatchedPolicy", "Error"}}}};
  IoData boundary_mode_iodata(boundary_mode_config, false);
  boundary_mode_iodata.boundaries.cracked_attributes.insert(9);
  boundary_mode_iodata.boundaries.cracked_attributes.insert(10);
  MaterialOperator boundary_mode_material(boundary_mode_iodata, *automatic_meshes.back());
  BoundaryModeOperator boundary_mode_op(boundary_mode_iodata, automatic_meshes,
                                        boundary_mode_material);
  SurfaceResponseOperator boundary_mode_response(boundary_mode_iodata, boundary_mode_op);
  CHECK(boundary_mode_response.GetPatchCount() == automatic_response.GetPatchCount());
  CHECK(boundary_mode_response.GetBasisSize() == automatic_response.GetBasisSize());
  CHECK(boundary_mode_response.HasSurfaceResponse());
  CHECK(boundary_mode_response.GetTargetInterfaces() == std::set<int>{4});

  GridFunction boundary_mode_field(boundary_mode_op.GetNDSpace(), true);
  mfem::Vector boundary_mode_field_value(2);
  boundary_mode_field_value[0] = 0.7;
  boundary_mode_field_value[1] = -0.4;
  mfem::VectorConstantCoefficient boundary_mode_field_coefficient(
      boundary_mode_field_value);
  boundary_mode_field.Real().ProjectCoefficient(boundary_mode_field_coefficient);
  boundary_mode_field.Imag() = 0.0;
  const auto boundary_mode_result =
      boundary_mode_response.GetMaxwellResponse(boundary_mode_field, 0.0);
  CHECK(boundary_mode_result.fabricated_surface_energy.at(4) > 0.0);
  CHECK(boundary_mode_result.fabricated_surface_energy_fixed_flux.at(4) > 0.0);
  CHECK(boundary_mode_result.loop_residual < 1.0e-10);
  CHECK_THAT(boundary_mode_result.matched_length_fraction, WithinAbs(1.0, 1.0e-12));

  // Finite-impedance boundary-mode edges use the local metal surface as their voltage
  // reference. For a conservative field normal to the sheet this is gauge-equivalent to
  // the PEC reference displaced into the metal.
  auto impedance_boundary_mode_config = boundary_mode_config;
  impedance_boundary_mode_config["Boundaries"].erase("PEC");
  impedance_boundary_mode_config["Boundaries"]["Impedance"] = {
      {{"Attributes", {9, 10}}, {"Ls", 1.0e-13}}};
  impedance_boundary_mode_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      impedance_library_path.string();
  IoData impedance_boundary_mode_iodata(impedance_boundary_mode_config, false);
  impedance_boundary_mode_iodata.boundaries.cracked_attributes.insert(9);
  impedance_boundary_mode_iodata.boundaries.cracked_attributes.insert(10);
  MaterialOperator impedance_boundary_mode_material(impedance_boundary_mode_iodata,
                                                    *automatic_meshes.back());
  BoundaryModeOperator impedance_boundary_mode_op(
      impedance_boundary_mode_iodata, automatic_meshes, impedance_boundary_mode_material);
  SurfaceResponseOperator impedance_boundary_mode_response(impedance_boundary_mode_iodata,
                                                           impedance_boundary_mode_op);
  CHECK(impedance_boundary_mode_response.GetPatchCount() ==
        boundary_mode_response.GetPatchCount());
  GridFunction boundary_mode_normal_field(boundary_mode_op.GetNDSpace(), true);
  mfem::Vector boundary_mode_normal_value(2);
  boundary_mode_normal_value[0] = 0.0;
  boundary_mode_normal_value[1] = -1.0;
  mfem::VectorConstantCoefficient boundary_mode_normal_coefficient(
      boundary_mode_normal_value);
  boundary_mode_normal_field.Real().ProjectCoefficient(boundary_mode_normal_coefficient);
  boundary_mode_normal_field.Imag() = 0.0;
  GridFunction impedance_boundary_mode_normal_field(impedance_boundary_mode_op.GetNDSpace(),
                                                    true);
  impedance_boundary_mode_normal_field.Real().ProjectCoefficient(
      boundary_mode_normal_coefficient);
  impedance_boundary_mode_normal_field.Imag() = 0.0;
  const auto boundary_mode_normal_result =
      boundary_mode_response.GetMaxwellResponse(boundary_mode_normal_field, 0.0);
  const auto impedance_boundary_mode_result =
      impedance_boundary_mode_response.GetMaxwellResponse(
          impedance_boundary_mode_normal_field, 0.0);
  CHECK_THAT(impedance_boundary_mode_result.domain_correction,
             WithinRel(boundary_mode_normal_result.domain_correction, 1.0e-10));
  CHECK_THAT(
      impedance_boundary_mode_result.fabricated_surface_energy.at(4),
      WithinRel(boundary_mode_normal_result.fabricated_surface_energy.at(4), 1.0e-10));
  CHECK(impedance_boundary_mode_result.loop_residual < 1.0e-10);
  CHECK(impedance_boundary_mode_result.boundary_law_verified);

  auto nondimensionalized_impedance_config = impedance_boundary_mode_config;
  nondimensionalized_impedance_config["Model"]["L0"] = 1.0e-6;
  nondimensionalized_impedance_config["Model"]["Lc"] = 1.0;
  IoData nondimensionalized_impedance_iodata(nondimensionalized_impedance_config, false);
  nondimensionalized_impedance_iodata.boundaries.cracked_attributes.insert(9);
  nondimensionalized_impedance_iodata.boundaries.cracked_attributes.insert(10);
  auto nondimensionalized_serial = std::make_unique<mfem::Mesh>(automatic_serial);
  nondimensionalized_impedance_iodata.NondimensionalizeInputs(nondimensionalized_serial);
  auto nondimensionalized_parallel =
      std::make_unique<mfem::ParMesh>(Mpi::World(), *nondimensionalized_serial);
  std::vector<std::unique_ptr<Mesh>> nondimensionalized_meshes;
  nondimensionalized_meshes.push_back(
      std::make_unique<Mesh>(std::move(nondimensionalized_parallel)));
  MaterialOperator nondimensionalized_impedance_material(
      nondimensionalized_impedance_iodata, *nondimensionalized_meshes.back());
  BoundaryModeOperator nondimensionalized_impedance_op(
      nondimensionalized_impedance_iodata, nondimensionalized_meshes,
      nondimensionalized_impedance_material);
  SurfaceResponseOperator nondimensionalized_impedance_response(
      nondimensionalized_impedance_iodata, nondimensionalized_impedance_op);
  GridFunction nondimensionalized_impedance_field(
      nondimensionalized_impedance_op.GetNDSpace(), true);
  nondimensionalized_impedance_field.Real().ProjectCoefficient(
      boundary_mode_normal_coefficient);
  nondimensionalized_impedance_field.Imag() = 0.0;
  const auto nondimensionalized_impedance_result =
      nondimensionalized_impedance_response.GetMaxwellResponse(
          nondimensionalized_impedance_field, 0.0);
  CHECK(nondimensionalized_impedance_result.boundary_law_verified);
  CHECK(nondimensionalized_impedance_result.loop_residual < 1.0e-10);

  auto impedance_mismatch_config = impedance_boundary_mode_config;
  impedance_mismatch_config["Boundaries"]["Impedance"][0]["Ls"] = 2.0e-13;
  IoData impedance_mismatch_iodata(impedance_mismatch_config, false);
  impedance_mismatch_iodata.boundaries.cracked_attributes.insert(9);
  impedance_mismatch_iodata.boundaries.cracked_attributes.insert(10);
  MaterialOperator impedance_mismatch_material(impedance_mismatch_iodata,
                                               *automatic_meshes.back());
  BoundaryModeOperator impedance_mismatch_op(impedance_mismatch_iodata, automatic_meshes,
                                             impedance_mismatch_material);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(impedance_mismatch_iodata, impedance_mismatch_op),
      Catch::Matchers::ContainsSubstring(
          "Automatic fabrication-process response matching failed"));

  auto invalid_boundary_law_config = impedance_boundary_mode_config;
  invalid_boundary_law_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      invalid_boundary_law_library_path.string();
  IoData invalid_boundary_law_iodata(invalid_boundary_law_config, false);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(invalid_boundary_law_iodata, impedance_boundary_mode_op),
      Catch::Matchers::ContainsSubstring(
          "Unknown fabrication-process response BoundaryCondition key \"Inductance\""));

  auto legacy_impedance_config = impedance_boundary_mode_config;
  legacy_impedance_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      legacy_impedance_library_path.string();
  IoData legacy_impedance_iodata(legacy_impedance_config, false);
  legacy_impedance_iodata.boundaries.cracked_attributes.insert(9);
  legacy_impedance_iodata.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator legacy_impedance_response(legacy_impedance_iodata,
                                                    impedance_boundary_mode_op);
  const auto legacy_impedance_result = legacy_impedance_response.GetMaxwellResponse(
      impedance_boundary_mode_normal_field, 0.0);
  CHECK_FALSE(legacy_impedance_result.boundary_law_verified);
  CHECK_FALSE(legacy_impedance_result.closure_independent_confident);
  CHECK_FALSE(legacy_impedance_result.confident);

  auto conductivity_boundary_mode_config = boundary_mode_config;
  conductivity_boundary_mode_config["Boundaries"].erase("PEC");
  conductivity_boundary_mode_config["Boundaries"]["Conductivity"] = {
      {{"Attributes", {9, 10}},
       {"Conductivity", 5.8e7},
       {"Permeability", 1.2},
       {"Thickness", 1.0e-7},
       {"External", true}}};
  conductivity_boundary_mode_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      conductivity_library_path.string();
  IoData conductivity_boundary_mode_iodata(conductivity_boundary_mode_config, false);
  conductivity_boundary_mode_iodata.boundaries.cracked_attributes.insert(9);
  conductivity_boundary_mode_iodata.boundaries.cracked_attributes.insert(10);
  MaterialOperator conductivity_boundary_mode_material(conductivity_boundary_mode_iodata,
                                                       *automatic_meshes.back());
  BoundaryModeOperator conductivity_boundary_mode_op(conductivity_boundary_mode_iodata,
                                                     automatic_meshes,
                                                     conductivity_boundary_mode_material);
  SurfaceResponseOperator conductivity_boundary_mode_response(
      conductivity_boundary_mode_iodata, conductivity_boundary_mode_op);
  GridFunction conductivity_boundary_mode_field(conductivity_boundary_mode_op.GetNDSpace(),
                                                true);
  conductivity_boundary_mode_field.Real().ProjectCoefficient(
      boundary_mode_normal_coefficient);
  conductivity_boundary_mode_field.Imag() = 0.0;
  const auto conductivity_boundary_mode_result =
      conductivity_boundary_mode_response.GetMaxwellResponse(
          conductivity_boundary_mode_field, 0.0);
  CHECK(conductivity_boundary_mode_result.boundary_law_verified);
  CHECK(conductivity_boundary_mode_result.loop_residual < 1.0e-10);

  auto rational_boundary_mode_config = boundary_mode_config;
  rational_boundary_mode_config["Boundaries"].erase("PEC");
  rational_boundary_mode_config["Boundaries"]["RationalImpedance"] = {
      {{"Attributes", {9, 10}},
       {"Numerator", rational_impedance_law["Numerator"]},
       {"Denominator", rational_impedance_law["Denominator"]}}};
  rational_boundary_mode_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      rational_impedance_library_path.string();
  IoData rational_boundary_mode_iodata(rational_boundary_mode_config, false);
  rational_boundary_mode_iodata.boundaries.cracked_attributes.insert(9);
  rational_boundary_mode_iodata.boundaries.cracked_attributes.insert(10);
  MaterialOperator rational_boundary_mode_material(rational_boundary_mode_iodata,
                                                   *automatic_meshes.back());
  BoundaryModeOperator rational_boundary_mode_op(
      rational_boundary_mode_iodata, automatic_meshes, rational_boundary_mode_material);
  SurfaceResponseOperator rational_boundary_mode_response(rational_boundary_mode_iodata,
                                                          rational_boundary_mode_op);
  GridFunction rational_boundary_mode_field(rational_boundary_mode_op.GetNDSpace(), true);
  rational_boundary_mode_field.Real().ProjectCoefficient(boundary_mode_normal_coefficient);
  rational_boundary_mode_field.Imag() = 0.0;
  const auto rational_boundary_mode_result =
      rational_boundary_mode_response.GetMaxwellResponse(rational_boundary_mode_field, 0.0);
  CHECK(rational_boundary_mode_result.boundary_law_verified);
  CHECK(rational_boundary_mode_result.loop_residual < 1.0e-10);

  auto different_pair_config = boundary_mode_config;
  different_pair_config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeAttributes"] =
      {9, 10};
  different_pair_config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"] =
      {0.25};
  different_pair_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      different_pair_library_2d_path.string();
  IoData different_pair_iodata(different_pair_config, false);
  different_pair_iodata.boundaries.cracked_attributes.insert(9);
  different_pair_iodata.boundaries.cracked_attributes.insert(10);

  mfem::Mesh different_pair_serial =
      mfem::Mesh::MakeCartesian2D(10, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0);
  for (int face = 0; face < different_pair_serial.GetNumFaces(); face++)
  {
    int element1, element2;
    different_pair_serial.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    different_pair_serial.GetFaceVertices(face, vertices);
    if (vertices.Size() != 2)
    {
      continue;
    }
    const double *p0 = different_pair_serial.GetVertex(vertices[0]);
    const double *p1 = different_pair_serial.GetVertex(vertices[1]);
    const double xmin = std::min(p0[0], p1[0]);
    const double xmax = std::max(p0[0], p1[0]);
    if (std::abs(p0[1] - 0.5) < 1.0e-12 && std::abs(p1[1] - 0.5) < 1.0e-12 &&
        (xmax <= 0.3 + 1.0e-12 || xmin >= 0.7 - 1.0e-12))
    {
      different_pair_serial.AddBdrElement(
          different_pair_serial.GetFace(face)->Duplicate(&different_pair_serial));
      different_pair_serial.SetBdrAttribute(different_pair_serial.GetNBE() - 1,
                                            xmax <= 0.3 + 1.0e-12 ? 9 : 10);
    }
  }
  different_pair_serial.FinalizeTopology();
  different_pair_serial.Finalize();
  while (different_pair_serial.GetNE() < Mpi::Size(Mpi::World()))
  {
    different_pair_serial.UniformRefinement();
  }
  auto different_pair_parallel =
      std::make_unique<mfem::ParMesh>(Mpi::World(), different_pair_serial);
  std::vector<std::unique_ptr<Mesh>> different_pair_meshes;
  different_pair_meshes.push_back(
      std::make_unique<Mesh>(std::move(different_pair_parallel)));
  MaterialOperator different_pair_material(different_pair_iodata,
                                           *different_pair_meshes.back());
  BoundaryModeOperator different_pair_mode(different_pair_iodata, different_pair_meshes,
                                           different_pair_material);
  SurfaceResponseOperator different_pair_response(different_pair_iodata,
                                                  different_pair_mode);
  CHECK(different_pair_response.GetPatchCount() == 1);
  CHECK(different_pair_response.GetBasisSize() == 4);

  auto parallel_cluster_config_2d = automatic_config;
  parallel_cluster_config_2d["Boundaries"]["Ground"]["Attributes"] = {1, 3, 4, 9};
  parallel_cluster_config_2d["Boundaries"]["Terminal"][0]["Attributes"] = {2, 10};
  parallel_cluster_config_2d["Boundaries"]["Postprocessing"]["Dielectric"][0]
                            ["EdgeAttributes"] = {9, 10};
  parallel_cluster_config_2d["Boundaries"]["Postprocessing"]["Dielectric"][0]
                            ["EdgeDistances"] = {0.25};
  parallel_cluster_config_2d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      parallel_cluster_library_2d_path.string();
  IoData parallel_cluster_iodata_2d(parallel_cluster_config_2d, false);
  parallel_cluster_iodata_2d.boundaries.cracked_attributes.insert(9);
  parallel_cluster_iodata_2d.boundaries.cracked_attributes.insert(10);

  mfem::Mesh parallel_cluster_serial_2d =
      mfem::Mesh::MakeCartesian2D(10, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0);
  for (int face = 0; face < parallel_cluster_serial_2d.GetNumFaces(); face++)
  {
    int element1, element2;
    parallel_cluster_serial_2d.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    parallel_cluster_serial_2d.GetFaceVertices(face, vertices);
    if (vertices.Size() != 2)
    {
      continue;
    }
    const double *p0 = parallel_cluster_serial_2d.GetVertex(vertices[0]);
    const double *p1 = parallel_cluster_serial_2d.GetVertex(vertices[1]);
    const double xmin = std::min(p0[0], p1[0]);
    const double xmax = std::max(p0[0], p1[0]);
    const bool left_ground = xmax <= 0.2 + 1.0e-12;
    const bool trace = xmin >= 0.4 - 1.0e-12 && xmax <= 0.6 + 1.0e-12;
    const bool right_ground = xmin >= 0.8 - 1.0e-12;
    if (std::abs(p0[1] - 0.5) < 1.0e-12 && std::abs(p1[1] - 0.5) < 1.0e-12 &&
        (left_ground || trace || right_ground))
    {
      parallel_cluster_serial_2d.AddBdrElement(
          parallel_cluster_serial_2d.GetFace(face)->Duplicate(&parallel_cluster_serial_2d));
      parallel_cluster_serial_2d.SetBdrAttribute(parallel_cluster_serial_2d.GetNBE() - 1,
                                                 trace ? 10 : 9);
    }
  }
  parallel_cluster_serial_2d.FinalizeTopology();
  parallel_cluster_serial_2d.Finalize();
  while (parallel_cluster_serial_2d.GetNE() < Mpi::Size(Mpi::World()))
  {
    parallel_cluster_serial_2d.UniformRefinement();
  }
  auto parallel_cluster_parallel_2d =
      std::make_unique<mfem::ParMesh>(Mpi::World(), parallel_cluster_serial_2d);
  std::vector<std::unique_ptr<Mesh>> parallel_cluster_meshes_2d;
  parallel_cluster_meshes_2d.push_back(
      std::make_unique<Mesh>(std::move(parallel_cluster_parallel_2d)));
  LaplaceOperator parallel_cluster_laplace_2d(parallel_cluster_iodata_2d,
                                              parallel_cluster_meshes_2d);
  SurfaceResponseOperator parallel_cluster_response_2d(parallel_cluster_iodata_2d,
                                                       parallel_cluster_laplace_2d);
  CHECK(parallel_cluster_response_2d.GetPatchCount() == 1);
  CHECK(parallel_cluster_response_2d.GetBasisSize() == 5);

  auto parallel_cluster_boundary_config_2d = parallel_cluster_config_2d;
  parallel_cluster_boundary_config_2d["Problem"]["Type"] = "BoundaryMode";
  parallel_cluster_boundary_config_2d["Boundaries"].erase("Ground");
  parallel_cluster_boundary_config_2d["Boundaries"].erase("Terminal");
  parallel_cluster_boundary_config_2d["Boundaries"]["PEC"] = {{"Attributes", {9, 10}}};
  parallel_cluster_boundary_config_2d["Solver"] = {
      {"Order", 1},
      {"BoundaryMode", {{"Freq", 5.0}}},
      {"SurfaceResponseCorrection",
       {{"Library", parallel_cluster_library_2d_path.string()},
        {"UnmatchedPolicy", "Error"}}}};
  IoData parallel_cluster_boundary_iodata_2d(parallel_cluster_boundary_config_2d, false);
  parallel_cluster_boundary_iodata_2d.boundaries.cracked_attributes.insert(9);
  parallel_cluster_boundary_iodata_2d.boundaries.cracked_attributes.insert(10);
  MaterialOperator parallel_cluster_boundary_material_2d(
      parallel_cluster_boundary_iodata_2d, *parallel_cluster_meshes_2d.back());
  BoundaryModeOperator parallel_cluster_boundary_op_2d(
      parallel_cluster_boundary_iodata_2d, parallel_cluster_meshes_2d,
      parallel_cluster_boundary_material_2d);
  SurfaceResponseOperator parallel_cluster_boundary_response_2d(
      parallel_cluster_boundary_iodata_2d, parallel_cluster_boundary_op_2d);
  CHECK(parallel_cluster_boundary_response_2d.GetPatchCount() == 1);
  CHECK(parallel_cluster_boundary_response_2d.GetBasisSize() == 6);

  // Finite-impedance conductors use the same connected-component ownership as PEC.
  // In particular, the two disconnected ground strips share attribute 9 but remain
  // distinct conductors, selecting the three-conductor cluster response.
  auto impedance_cluster_boundary_config_2d = parallel_cluster_boundary_config_2d;
  impedance_cluster_boundary_config_2d["Boundaries"].erase("PEC");
  impedance_cluster_boundary_config_2d["Boundaries"]["Impedance"] = {
      {{"Attributes", {9, 10}}, {"Ls", 1.0e-13}}};
  impedance_cluster_boundary_config_2d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      impedance_parallel_cluster_library_2d_path.string();
  IoData impedance_cluster_boundary_iodata_2d(impedance_cluster_boundary_config_2d, false);
  impedance_cluster_boundary_iodata_2d.boundaries.cracked_attributes.insert(9);
  impedance_cluster_boundary_iodata_2d.boundaries.cracked_attributes.insert(10);
  MaterialOperator impedance_cluster_boundary_material_2d(
      impedance_cluster_boundary_iodata_2d, *parallel_cluster_meshes_2d.back());
  BoundaryModeOperator impedance_cluster_boundary_op_2d(
      impedance_cluster_boundary_iodata_2d, parallel_cluster_meshes_2d,
      impedance_cluster_boundary_material_2d);
  SurfaceResponseOperator impedance_cluster_boundary_response_2d(
      impedance_cluster_boundary_iodata_2d, impedance_cluster_boundary_op_2d);
  CHECK(impedance_cluster_boundary_response_2d.GetPatchCount() == 1);
  CHECK(impedance_cluster_boundary_response_2d.GetBasisSize() == 6);

  mfem::FunctionCoefficient cluster_potential_coefficient(
      [](const mfem::Vector &x)
      {
        const double distance = std::abs(x[0] - 0.5);
        if (distance <= 0.1)
        {
          return 1.0;
        }
        if (distance >= 0.3)
        {
          return 0.0;
        }
        const double t = (distance - 0.1) / 0.2;
        return 1.0 - 3.0 * t * t + 2.0 * t * t * t;
      });
  mfem::ParGridFunction parallel_cluster_potential_2d(
      &parallel_cluster_laplace_2d.GetH1Space().Get());
  parallel_cluster_potential_2d.ProjectCoefficient(cluster_potential_coefficient);
  Vector parallel_cluster_potential_true_2d;
  parallel_cluster_potential_2d.GetTrueDofs(parallel_cluster_potential_true_2d);
  const auto parallel_cluster_electrostatic_result_2d =
      parallel_cluster_response_2d.GetElectrostaticResponse(
          parallel_cluster_potential_true_2d);

  GridFunction parallel_cluster_boundary_field_2d(
      parallel_cluster_boundary_op_2d.GetNDSpace(), true);
  mfem::ParGridFunction parallel_cluster_boundary_potential_2d(
      &parallel_cluster_boundary_op_2d.GetH1Space().Get());
  parallel_cluster_boundary_potential_2d.ProjectCoefficient(cluster_potential_coefficient);
  mfem::ParDiscreteLinearOperator parallel_cluster_gradient_2d(
      &parallel_cluster_boundary_op_2d.GetH1Space().Get(),
      &parallel_cluster_boundary_op_2d.GetNDSpace().Get());
  parallel_cluster_gradient_2d.AddDomainInterpolator(new mfem::GradientInterpolator());
  parallel_cluster_gradient_2d.Assemble();
  parallel_cluster_gradient_2d.Mult(parallel_cluster_boundary_potential_2d,
                                    parallel_cluster_boundary_field_2d.Real());
  parallel_cluster_boundary_field_2d.Real() *= -1.0;
  parallel_cluster_boundary_field_2d.Imag() = 0.0;
  const auto parallel_cluster_boundary_result_2d =
      parallel_cluster_boundary_response_2d.GetMaxwellResponse(
          parallel_cluster_boundary_field_2d, 0.0);
  CHECK(parallel_cluster_boundary_result_2d.loop_residual < 1.0e-10);
  CHECK_THAT(
      parallel_cluster_boundary_result_2d.domain_correction,
      WithinRel(parallel_cluster_electrostatic_result_2d.domain_correction, 1.0e-10));
  CHECK_THAT(
      parallel_cluster_boundary_result_2d.domain_correction_fixed_flux,
      WithinRel(parallel_cluster_electrostatic_result_2d.domain_correction_fixed_flux,
                1.0e-10));
  CHECK_THAT(
      parallel_cluster_boundary_result_2d.fabricated_surface_energy.at(4),
      WithinRel(parallel_cluster_electrostatic_result_2d.fabricated_surface_energy.at(4),
                1.0e-10));

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

  auto disconnected_cluster_config = automatic_config;
  disconnected_cluster_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      disconnected_cluster_library_3d_path.string();
  IoData disconnected_cluster_iodata(disconnected_cluster_config, false);
  disconnected_cluster_iodata.boundaries.cracked_attributes.insert(9);
  disconnected_cluster_iodata.boundaries.cracked_attributes.insert(10);
  CHECK_THROWS_WITH(SurfaceResponseOperator(disconnected_cluster_iodata, automatic_laplace),
                    Catch::Matchers::ContainsSubstring(
                        "OpenContourPaths must connect every conductor reference"));

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
  std::set<std::size_t> physical_edge_vertices;
  for (const std::size_t segment_index : segment_indices)
  {
    const auto &segment = geometry_3d.segments[segment_index];
    physical_edge_vertices.insert(segment.vertices.begin(), segment.vertices.end());
    const auto &p0 = geometry_3d.vertices[segment.vertices[0]].coordinate;
    const auto &p1 = geometry_3d.vertices[segment.vertices[1]].coordinate;
    double length_squared = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length_squared += (p1[d] - p0[d]) * (p1[d] - p0[d]);
    }
    physical_edge_length += std::sqrt(length_squared);
  }
  const int physical_endpoint_count = static_cast<int>(std::count_if(
      physical_edge_vertices.begin(), physical_edge_vertices.end(),
      [&](std::size_t vertex)
      {
        return geometry_3d.vertices[vertex].physical_type == MetalEdgeVertexType::ENDPOINT;
      }));
  REQUIRE(physical_endpoint_count > 0);
  CHECK(static_cast<int>(
            std::count_if(physical_edge_vertices.begin(), physical_edge_vertices.end(),
                          [&](std::size_t vertex)
                          {
                            return geometry_3d.vertices[vertex].physical_type ==
                                       MetalEdgeVertexType::ENDPOINT &&
                                   geometry_3d.vertices[vertex].on_truncation_boundary;
                          })) == physical_endpoint_count);
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

  const auto endpoint_requirements_path =
      temp.temp_dir / "surface-response-requirements-endpoint.json";
  WriteSurfaceResponseRequirements(iodata_3d, *meshes_3d.back(),
                                   endpoint_requirements_path.string());
  std::ifstream endpoint_requirements_input(endpoint_requirements_path);
  REQUIRE(endpoint_requirements_input);
  const json endpoint_requirements = json::parse(endpoint_requirements_input);
  const auto endpoint_requirement = std::find_if(
      endpoint_requirements["Requirements"].begin(),
      endpoint_requirements["Requirements"].end(),
      [](const auto &requirement)
      {
        return requirement["Topology"] == "Endpoint" && requirement["Status"] == "Missing";
      });
  CHECK(endpoint_requirement == endpoint_requirements["Requirements"].end());

  auto endpoint_config_3d = config_3d;
  endpoint_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      endpoint_library_3d_path.string();
  IoData endpoint_iodata_3d(endpoint_config_3d, false);
  auto endpoint_mesh_3d = mesh::ReadMesh(endpoint_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> endpoint_meshes_3d;
  endpoint_meshes_3d.push_back(std::make_unique<Mesh>(std::move(endpoint_mesh_3d)));
  LaplaceOperator endpoint_laplace_3d(endpoint_iodata_3d, endpoint_meshes_3d);
  SurfaceResponseOperator endpoint_response_3d(endpoint_iodata_3d, endpoint_laplace_3d);
  CHECK(endpoint_response_3d.GetPatchCount() == response_3d.GetPatchCount());
  CHECK(endpoint_response_3d.GetBasisSize() == response_3d.GetBasisSize());
  CHECK_THAT(endpoint_response_3d.GetPatchWeight(),
             WithinRel(response_3d.GetPatchWeight(), 1.0e-12));

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

  auto endpoint_maxwell_config_3d = maxwell_config_3d;
  endpoint_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      endpoint_library_3d_path.string();
  IoData endpoint_maxwell_iodata_3d(endpoint_maxwell_config_3d, false);
  auto endpoint_maxwell_mesh_3d = mesh::ReadMesh(endpoint_maxwell_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> endpoint_maxwell_meshes_3d;
  endpoint_maxwell_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(endpoint_maxwell_mesh_3d)));
  SpaceOperator endpoint_maxwell_space_3d(endpoint_maxwell_iodata_3d,
                                          endpoint_maxwell_meshes_3d);
  SurfaceResponseOperator endpoint_maxwell_response_3d(endpoint_maxwell_iodata_3d,
                                                       endpoint_maxwell_space_3d);
  CHECK(endpoint_maxwell_response_3d.GetPatchCount() ==
        maxwell_response_3d.GetPatchCount());

  // A finite-impedance metal does not provide an interior equipotential anchor. The
  // automatic Maxwell trace instead references the local metal edge point and must still
  // reproduce a conservative quasi-electrostatic field.
  auto impedance_maxwell_config_3d = maxwell_config_3d;
  impedance_maxwell_config_3d["Boundaries"].erase("Ground");
  impedance_maxwell_config_3d["Boundaries"]["Impedance"] = {
      {{"Attributes", {1, 2}}, {"Ls", 1.0e-13}}};
  impedance_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      impedance_library_3d_path.string();
  IoData impedance_maxwell_iodata_3d(impedance_maxwell_config_3d, false);
  auto impedance_maxwell_mesh_3d =
      mesh::ReadMesh(impedance_maxwell_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> impedance_maxwell_meshes_3d;
  impedance_maxwell_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(impedance_maxwell_mesh_3d)));
  SpaceOperator impedance_maxwell_space_3d(impedance_maxwell_iodata_3d,
                                           impedance_maxwell_meshes_3d);
  SurfaceResponseOperator impedance_maxwell_response_3d(impedance_maxwell_iodata_3d,
                                                        impedance_maxwell_space_3d);
  CHECK(impedance_maxwell_response_3d.GetPatchCount() ==
        maxwell_response_3d.GetPatchCount());
  GridFunction impedance_maxwell_field(impedance_maxwell_space_3d.GetNDSpace(), true);
  mfem::Vector impedance_gradient(3);
  impedance_gradient = 0.0;
  impedance_gradient[1] = -1.0;
  mfem::VectorConstantCoefficient impedance_gradient_coefficient(impedance_gradient);
  impedance_maxwell_field.Real().ProjectCoefficient(impedance_gradient_coefficient);
  impedance_maxwell_field.Imag() = 0.0;
  const auto impedance_gradient_response =
      impedance_maxwell_response_3d.GetMaxwellResponse(impedance_maxwell_field, 0.0);
  CHECK(std::abs(impedance_gradient_response.domain_correction) > 0.0);
  CHECK(impedance_gradient_response.loop_residual < 1.0e-10);
  CHECK(impedance_gradient_response.boundary_law_verified);

  GridFunction maxwell_field(maxwell_space_3d.GetNDSpace(), true);
  mfem::Vector constant_field(3);
  constant_field[0] = 0.7;
  constant_field[1] = -0.4;
  constant_field[2] = 0.2;
  mfem::VectorConstantCoefficient field_coefficient(constant_field);
  maxwell_field.Real().ProjectCoefficient(field_coefficient);
  maxwell_field.Imag() = 0.0;
  GridFunction endpoint_maxwell_field(endpoint_maxwell_space_3d.GetNDSpace(), true);
  endpoint_maxwell_field.Real().ProjectCoefficient(field_coefficient);
  endpoint_maxwell_field.Imag() = 0.0;
  const auto endpoint_maxwell_result =
      endpoint_maxwell_response_3d.GetMaxwellResponse(endpoint_maxwell_field, 0.0);
  CHECK(endpoint_maxwell_result.loop_residual < 1.0e-10);
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
  CHECK(real_response.response_weighted_trace_closure_spread > 0.0);
  CHECK(real_response.trace_closure_response_failure_fraction >= 0.0);
  CHECK(real_response.trace_closure_response_failure_fraction <= 1.0);
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
  CHECK_THAT(
      gradient_response.response_weighted_trace_closure_spread,
      WithinRel(electrostatic_response.response_weighted_trace_closure_spread, 1.0e-10));
  CHECK_THAT(
      gradient_response.trace_closure_response_failure_fraction,
      WithinAbs(electrostatic_response.trace_closure_response_failure_fraction, 1.0e-12));
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
  CHECK_THAT(complex_response.response_weighted_trace_closure_spread,
             WithinRel(real_response.response_weighted_trace_closure_spread, 1.0e-10));
  CHECK_THAT(complex_response.trace_closure_response_failure_fraction,
             WithinAbs(real_response.trace_closure_response_failure_fraction, 1.0e-12));

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

  auto missing_pair_config_3d = coupled_config_3d;
  missing_pair_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      missing_pair_library_3d_path.string();
  missing_pair_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]
                        ["UnmatchedPolicy"] = "Warn";
  IoData missing_pair_iodata_3d(missing_pair_config_3d, false);
  missing_pair_iodata_3d.boundaries.cracked_attributes.insert(1);
  missing_pair_iodata_3d.boundaries.cracked_attributes.insert(2);
  auto missing_pair_mesh_3d = mesh::ReadMesh(missing_pair_iodata_3d, Mpi::World());
  Mesh missing_pair_mesh(std::move(missing_pair_mesh_3d));
  const auto missing_pair_requirements_path =
      temp.temp_dir / "surface-response-requirements-missing-pair.json";
  WriteSurfaceResponseRequirements(missing_pair_iodata_3d, missing_pair_mesh,
                                   missing_pair_requirements_path.string());
  std::ifstream missing_pair_requirements_input(missing_pair_requirements_path);
  REQUIRE(missing_pair_requirements_input);
  const json missing_pair_requirements = json::parse(missing_pair_requirements_input);
  CHECK_FALSE(missing_pair_requirements["Complete"]);
  CHECK(missing_pair_requirements["Summary"]["Counts"]["Missing"].get<int>() > 0);
  const auto missing_pair_requirement =
      std::find_if(missing_pair_requirements["Requirements"].begin(),
                   missing_pair_requirements["Requirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "DifferentConductorGap" &&
                            requirement["Status"] == "Missing";
                   });
  REQUIRE(missing_pair_requirement != missing_pair_requirements["Requirements"].end());
  CHECK((*missing_pair_requirement)["Geometry"]["EdgeCount"] == 2);
  CHECK_THAT((*missing_pair_requirement)["Geometry"]["Separation"].get<double>(),
             WithinAbs(12.0, 1.0e-9));

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

  auto parallel_cluster_config_3d = coupled_config_3d;
  for (auto &interface :
       parallel_cluster_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 11.0};
  }
  parallel_cluster_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      parallel_cluster_library_3d_path.string();
  IoData parallel_cluster_iodata_3d(parallel_cluster_config_3d, false);
  parallel_cluster_iodata_3d.boundaries.cracked_attributes.insert(1);
  parallel_cluster_iodata_3d.boundaries.cracked_attributes.insert(2);
  auto parallel_cluster_mesh_3d = mesh::ReadMesh(parallel_cluster_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> parallel_cluster_meshes_3d;
  parallel_cluster_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(parallel_cluster_mesh_3d)));
  LaplaceOperator parallel_cluster_laplace_3d(parallel_cluster_iodata_3d,
                                              parallel_cluster_meshes_3d);
  SurfaceResponseOperator parallel_cluster_response_3d(parallel_cluster_iodata_3d,
                                                       parallel_cluster_laplace_3d);
  CHECK(parallel_cluster_response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size() / 4) * line_rule.GetNPoints());
  CHECK(parallel_cluster_response_3d.GetBasisSize() ==
        5 * parallel_cluster_response_3d.GetPatchCount());
  CHECK_THAT(parallel_cluster_response_3d.GetPatchWeight(),
             WithinRel(0.125 * physical_edge_length, 1.0e-12));

  const auto parallel_cluster_requirements_path =
      temp.temp_dir / "surface-response-requirements-parallel-cluster.json";
  WriteSurfaceResponseRequirements(parallel_cluster_iodata_3d,
                                   *parallel_cluster_meshes_3d.back(),
                                   parallel_cluster_requirements_path.string());
  std::ifstream parallel_cluster_requirements_input(parallel_cluster_requirements_path);
  REQUIRE(parallel_cluster_requirements_input);
  const json parallel_cluster_requirements =
      json::parse(parallel_cluster_requirements_input);
  const auto parallel_cluster_requirement =
      std::find_if(parallel_cluster_requirements["Requirements"].begin(),
                   parallel_cluster_requirements["Requirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "ParallelEdgeCluster" &&
                            requirement["Status"] == "Exact";
                   });
  REQUIRE(parallel_cluster_requirement !=
          parallel_cluster_requirements["Requirements"].end());
  CHECK((*parallel_cluster_requirement)["Geometry"]["EdgeCount"] == 4);
  CHECK((*parallel_cluster_requirement)["Geometry"]["Edges"].size() == 4);
  std::set<int> parallel_cluster_conductors;
  for (const auto &edge : (*parallel_cluster_requirement)["Geometry"]["Edges"])
  {
    const int conductor = edge["Conductor"];
    CHECK(conductor > 0);
    parallel_cluster_conductors.insert(conductor);
  }
  CHECK(parallel_cluster_conductors == std::set<int>{1, 2});
  CHECK((*parallel_cluster_requirement)["TotalEdgeLength"].get<double>() > 0.0);

  // An exact multi-edge coupon is self-contained. It must not require redundant
  // two-edge models for every pair in the active cluster.
  auto parallel_cluster_only_config_3d = parallel_cluster_config_3d;
  parallel_cluster_only_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]
                                 ["Library"] =
                                     parallel_cluster_only_library_3d_path.string();
  IoData parallel_cluster_only_iodata_3d(parallel_cluster_only_config_3d, false);
  parallel_cluster_only_iodata_3d.boundaries.cracked_attributes.insert(1);
  parallel_cluster_only_iodata_3d.boundaries.cracked_attributes.insert(2);
  auto parallel_cluster_only_mesh_3d =
      mesh::ReadMesh(parallel_cluster_only_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> parallel_cluster_only_meshes_3d;
  parallel_cluster_only_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(parallel_cluster_only_mesh_3d)));
  LaplaceOperator parallel_cluster_only_laplace_3d(parallel_cluster_only_iodata_3d,
                                                   parallel_cluster_only_meshes_3d);
  SurfaceResponseOperator parallel_cluster_only_response_3d(
      parallel_cluster_only_iodata_3d, parallel_cluster_only_laplace_3d);
  CHECK(parallel_cluster_only_response_3d.GetPatchCount() ==
        parallel_cluster_response_3d.GetPatchCount());
  CHECK(parallel_cluster_only_response_3d.GetBasisSize() ==
        parallel_cluster_response_3d.GetBasisSize());

  auto parallel_cluster_maxwell_config_3d = coupled_maxwell_config_3d;
  for (auto &interface :
       parallel_cluster_maxwell_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 11.0};
  }
  parallel_cluster_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      parallel_cluster_library_3d_path.string();
  IoData parallel_cluster_maxwell_iodata_3d(parallel_cluster_maxwell_config_3d, false);
  auto parallel_cluster_maxwell_mesh_3d =
      mesh::ReadMesh(parallel_cluster_maxwell_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> parallel_cluster_maxwell_meshes_3d;
  parallel_cluster_maxwell_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(parallel_cluster_maxwell_mesh_3d)));
  SpaceOperator parallel_cluster_maxwell_space_3d(parallel_cluster_maxwell_iodata_3d,
                                                  parallel_cluster_maxwell_meshes_3d);
  SurfaceResponseOperator parallel_cluster_maxwell_response_3d(
      parallel_cluster_maxwell_iodata_3d, parallel_cluster_maxwell_space_3d);
  CHECK(parallel_cluster_maxwell_response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size() / 4) * maxwell_line_rule.GetNPoints());
  CHECK(parallel_cluster_maxwell_response_3d.GetBasisSize() ==
        6 * parallel_cluster_maxwell_response_3d.GetPatchCount());

  auto parallel_cluster_only_maxwell_config_3d = parallel_cluster_maxwell_config_3d;
  parallel_cluster_only_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]
                                         ["Library"] =
                                             parallel_cluster_only_library_3d_path.string();
  IoData parallel_cluster_only_maxwell_iodata_3d(parallel_cluster_only_maxwell_config_3d,
                                                 false);
  auto parallel_cluster_only_maxwell_mesh_3d =
      mesh::ReadMesh(parallel_cluster_only_maxwell_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> parallel_cluster_only_maxwell_meshes_3d;
  parallel_cluster_only_maxwell_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(parallel_cluster_only_maxwell_mesh_3d)));
  SpaceOperator parallel_cluster_only_maxwell_space_3d(
      parallel_cluster_only_maxwell_iodata_3d, parallel_cluster_only_maxwell_meshes_3d);
  SurfaceResponseOperator parallel_cluster_only_maxwell_response_3d(
      parallel_cluster_only_maxwell_iodata_3d, parallel_cluster_only_maxwell_space_3d);
  CHECK(parallel_cluster_only_maxwell_response_3d.GetPatchCount() ==
        parallel_cluster_maxwell_response_3d.GetPatchCount());
  CHECK(parallel_cluster_only_maxwell_response_3d.GetBasisSize() ==
        parallel_cluster_maxwell_response_3d.GetBasisSize());

  GridFunction parallel_cluster_maxwell_field(
      parallel_cluster_maxwell_space_3d.GetNDSpace(), true);
  mfem::ParGridFunction parallel_cluster_maxwell_potential(
      &parallel_cluster_maxwell_space_3d.GetH1Space().Get());
  parallel_cluster_maxwell_potential.ProjectCoefficient(transverse_potential_coefficient);
  mfem::ParDiscreteLinearOperator parallel_cluster_gradient(
      &parallel_cluster_maxwell_space_3d.GetH1Space().Get(),
      &parallel_cluster_maxwell_space_3d.GetNDSpace().Get());
  parallel_cluster_gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
  parallel_cluster_gradient.Assemble();
  parallel_cluster_gradient.Mult(parallel_cluster_maxwell_potential,
                                 parallel_cluster_maxwell_field.Real());
  parallel_cluster_maxwell_field.Real() *= -1.0;
  parallel_cluster_maxwell_field.Imag() = 0.0;
  const auto parallel_cluster_maxwell_result =
      parallel_cluster_maxwell_response_3d.GetMaxwellResponse(
          parallel_cluster_maxwell_field, 0.0);

  mfem::ParGridFunction parallel_cluster_electrostatic_potential(
      &parallel_cluster_laplace_3d.GetH1Space().Get());
  parallel_cluster_electrostatic_potential.ProjectCoefficient(
      transverse_potential_coefficient);
  Vector parallel_cluster_electrostatic_true;
  parallel_cluster_electrostatic_potential.GetTrueDofs(parallel_cluster_electrostatic_true);
  const auto parallel_cluster_electrostatic_result =
      parallel_cluster_response_3d.GetElectrostaticResponse(
          parallel_cluster_electrostatic_true);
  CHECK(std::abs(parallel_cluster_maxwell_result.domain_correction) > 0.0);
  CHECK(parallel_cluster_maxwell_result.loop_residual < 1.0e-10);
  CHECK_THAT(parallel_cluster_maxwell_result.domain_correction,
             WithinRel(parallel_cluster_electrostatic_result.domain_correction, 1.0e-10));
  CHECK_THAT(parallel_cluster_maxwell_result.domain_correction_fixed_flux,
             WithinRel(parallel_cluster_electrostatic_result.domain_correction_fixed_flux,
                       1.0e-10));
  for (const auto &[interface, energy] :
       parallel_cluster_electrostatic_result.fabricated_surface_energy)
  {
    CHECK_THAT(parallel_cluster_maxwell_result.fabricated_surface_energy.at(interface),
               WithinRel(energy, 1.0e-5));
    CHECK_THAT(
        parallel_cluster_maxwell_result.fabricated_surface_energy_fixed_flux.at(interface),
        WithinRel(
            parallel_cluster_electrostatic_result.fabricated_surface_energy_fixed_flux.at(
                interface),
            1.0e-4));
  }

  Vector cluster_state(parallel_cluster_maxwell_space_3d.GetNDSpace().GetTrueVSize());
  Vector cluster_probe(cluster_state.Size());
  auto *cluster_state_data = cluster_state.HostWrite();
  auto *cluster_probe_data = cluster_probe.HostWrite();
  for (int i = 0; i < cluster_state.Size(); i++)
  {
    cluster_state_data[i] = std::cos(0.19 * (i + 1 + 5 * Mpi::Rank(Mpi::World())));
    cluster_probe_data[i] = std::sin(0.27 * (i + 1 + 7 * Mpi::Rank(Mpi::World())));
  }
  cluster_state.SetSubVector(parallel_cluster_maxwell_space_3d.GetNDDbcTDofLists().back(),
                             0.0);
  cluster_probe.SetSubVector(parallel_cluster_maxwell_space_3d.GetNDDbcTDofLists().back(),
                             0.0);
  Vector cluster_correction, cluster_probe_correction;
  parallel_cluster_maxwell_response_3d.Mult(cluster_state, cluster_correction);
  parallel_cluster_maxwell_response_3d.Mult(cluster_probe, cluster_probe_correction);
  CHECK_THAT(linalg::Dot(Mpi::World(), cluster_probe, cluster_correction),
             WithinRel(linalg::Dot(Mpi::World(), cluster_state, cluster_probe_correction),
                       1.0e-10));

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
  auto MakeIslandMesh = [](bool rounded = false, bool tetrahedral = false,
                           bool aperture = false, bool neighboring_island = false,
                           bool second_layer = false, bool high_order_rounded = false)
  {
    const double in_plane_extent = aperture || neighboring_island ? 2.0 : 1.0;
    const double center = 0.5 * in_plane_extent;
    const double half_width = 0.25;
    const int in_plane_elements = (rounded && !high_order_rounded ? 16 : 8) *
                                  (aperture || neighboring_island ? 2 : 1);
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
      const double plane_y = serial.GetVertex(vertices[0])[1];
      double xmin = in_plane_extent, xmax = 0.0;
      double zmin = in_plane_extent, zmax = 0.0;
      for (const int vertex : vertices)
      {
        const double *point = serial.GetVertex(vertex);
        on_plane = on_plane && std::abs(point[1] - plane_y) < 1.0e-12;
        xmin = std::min(xmin, point[0]);
        xmax = std::max(xmax, point[0]);
        zmin = std::min(zmin, point[2]);
        zmax = std::max(zmax, point[2]);
      }
      const double island_center_x = neighboring_island ? 0.625 : center;
      const bool inside_island = xmin >= island_center_x - half_width - 1.0e-12 &&
                                 xmax <= island_center_x + half_width + 1.0e-12 &&
                                 zmin >= center - half_width - 1.0e-12 &&
                                 zmax <= center + half_width + 1.0e-12;
      constexpr double neighbor_center_x = 1.375;
      const bool inside_neighbor =
          neighboring_island && xmin >= neighbor_center_x - half_width - 1.0e-12 &&
          xmax <= neighbor_center_x + half_width + 1.0e-12 &&
          zmin >= center - half_width - 1.0e-12 && zmax <= center + half_width + 1.0e-12;
      const bool selected_plane = second_layer ? (std::abs(plane_y - 0.25) < 1.0e-12 ||
                                                  std::abs(plane_y - 0.75) < 1.0e-12)
                                               : std::abs(plane_y - 0.5) < 1.0e-12;
      if (on_plane && selected_plane &&
          (aperture ? !inside_island : (inside_island || inside_neighbor)))
      {
        serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
        serial.SetBdrAttribute(serial.GetNBE() - 1, second_layer && plane_y > 0.5
                                                        ? 10
                                                        : (inside_neighbor ? 10 : 9));
      }
    }
    auto RoundPoint = [&](const mfem::Vector &input, mfem::Vector &output)
    {
      output = input;
      constexpr double radius = 0.125;
      constexpr double tolerance = 1.0e-12;
      if (std::abs(input[1] - 0.5) > tolerance)
      {
        return;
      }
      const std::vector<double> island_centers = neighboring_island
                                                     ? std::vector<double>{0.625, 1.375}
                                                     : std::vector<double>{center};
      for (const double island_center_x : island_centers)
      {
        for (const double sign_x : {-1.0, 1.0})
        {
          for (const double sign_z : {-1.0, 1.0})
          {
            const double corner_x = island_center_x + sign_x * half_width;
            const double corner_z = center + sign_z * half_width;
            const double center_x = corner_x - sign_x * radius;
            const double center_z = corner_z - sign_z * radius;
            const double local_x = sign_x * (input[0] - center_x);
            const double local_z = sign_z * (input[2] - center_z);
            if (local_x < -tolerance || local_x > radius + tolerance ||
                local_z < -tolerance || local_z > radius + tolerance)
            {
              continue;
            }

            double angle;
            if (std::abs(input[2] - corner_z) <= tolerance)
            {
              angle = 0.5 * std::acos(-1.0) - 0.25 * std::acos(-1.0) * local_x / radius;
            }
            else if (std::abs(input[0] - corner_x) <= tolerance)
            {
              angle = 0.25 * std::acos(-1.0) * local_z / radius;
            }
            else
            {
              continue;
            }
            output[0] = center_x + sign_x * radius * std::cos(angle);
            output[2] = center_z + sign_z * radius * std::sin(angle);
            return;
          }
        }
      }
    };
    if (rounded && !high_order_rounded)
    {
      for (int vertex = 0; vertex < serial.GetNV(); vertex++)
      {
        mfem::Vector input(serial.GetVertex(vertex), 3);
        mfem::Vector output(3);
        RoundPoint(input, output);
        for (int d = 0; d < 3; d++)
        {
          serial.GetVertex(vertex)[d] = output[d];
        }
      }
    }
    serial.FinalizeTopology();
    serial.Finalize();
    if (rounded && high_order_rounded)
    {
      serial.SetCurvature(2);
      serial.Transform(RoundPoint);
      for (int d = 0; d < serial.SpaceDimension(); d++)
      {
        mfem::Vector values;
        serial.GetNodes()->GetNodalValues(values, d + 1);
        for (int vertex = 0; vertex < serial.GetNV(); vertex++)
        {
          serial.GetVertex(vertex)[d] = values[vertex];
        }
      }
    }
    return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
  };
  auto MakeTouchingIslandMesh = []()
  {
    constexpr double extent = 2.0;
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(16, 4, 16, mfem::Element::HEXAHEDRON,
                                                    extent, 1.0, extent);
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
      const bool lower_left = xmin >= 0.25 - 1.0e-12 && xmax <= 1.0 + 1.0e-12 &&
                              zmin >= 0.25 - 1.0e-12 && zmax <= 1.0 + 1.0e-12;
      const bool upper_right = xmin >= 1.0 - 1.0e-12 && xmax <= 1.75 + 1.0e-12 &&
                               zmin >= 1.0 - 1.0e-12 && zmax <= 1.75 + 1.0e-12;
      if (on_plane && (lower_left || upper_right))
      {
        serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
        serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
      }
    }
    serial.FinalizeTopology();
    serial.Finalize();
    return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
  };
  auto MakeOffsetCornerPairMesh = []()
  {
    constexpr double extent = 2.0;
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(16, 4, 16, mfem::Element::HEXAHEDRON,
                                                    extent, 1.0, extent);
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
      const bool first = xmin >= 0.25 - 1.0e-12 && xmax <= 1.0 + 1.0e-12 &&
                         zmin >= 0.25 - 1.0e-12 && zmax <= 0.75 + 1.0e-12;
      const bool second = xmin >= 1.125 - 1.0e-12 && xmax <= 1.625 + 1.0e-12 &&
                          zmin >= 0.875 - 1.0e-12 && zmax <= 1.75 + 1.0e-12;
      if (on_plane && (first || second))
      {
        serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
        serial.SetBdrAttribute(serial.GetNBE() - 1, first ? 9 : 10);
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

  auto touching_geometry_mesh = MakeTouchingIslandMesh();
  const auto touching_geometry =
      ExtractMetalEdgeGeometry(*touching_geometry_mesh, convex_island_iodata.boundaries);
  const auto touching_segments =
      GetInterfaceMetalEdgeSegmentIndices(touching_geometry, 4, InterfaceDielectric::SA);
  std::set<std::size_t> touching_vertices;
  for (const std::size_t segment_index : touching_segments)
  {
    const auto &segment = touching_geometry.segments[segment_index];
    touching_vertices.insert(segment.vertices.begin(), segment.vertices.end());
  }
  const int touching_junctions = static_cast<int>(
      std::count_if(touching_vertices.begin(), touching_vertices.end(),
                    [&](std::size_t vertex)
                    {
                      return touching_geometry.vertices[vertex].physical_type ==
                             MetalEdgeVertexType::JUNCTION;
                    }));
  REQUIRE(touching_junctions == 1);

  std::vector<std::unique_ptr<Mesh>> touching_island_meshes;
  touching_island_meshes.push_back(std::make_unique<Mesh>(MakeTouchingIslandMesh()));
  LaplaceOperator touching_island_laplace(convex_island_iodata, touching_island_meshes);
  SurfaceResponseOperator touching_island_response(convex_island_iodata,
                                                   touching_island_laplace);
  auto junction_island_config = convex_island_config;
  junction_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      junction_library_3d_path.string();
  IoData junction_island_iodata(junction_island_config, false);
  junction_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> junction_island_meshes;
  junction_island_meshes.push_back(std::make_unique<Mesh>(MakeTouchingIslandMesh()));
  LaplaceOperator junction_island_laplace(junction_island_iodata, junction_island_meshes);
  SurfaceResponseOperator junction_island_response(junction_island_iodata,
                                                   junction_island_laplace);
  const int removed_junction_straight_patches = 4 * island_line_rule.GetNPoints();
  CHECK(junction_island_response.GetPatchCount() ==
        touching_island_response.GetPatchCount() - removed_junction_straight_patches +
            touching_junctions);
  CHECK(junction_island_response.GetBasisSize() ==
        touching_island_response.GetBasisSize() - 4 * removed_junction_straight_patches +
            12 * touching_junctions);
  CHECK_THAT(junction_island_response.GetPatchWeight(),
             WithinRel(touching_island_response.GetPatchWeight() - 3.0, 1.0e-12));

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

  // The same fillet represented by coarse quadratic edges must select the same four
  // rounded-corner coupons and leave the same straight response intervals.
  std::vector<std::unique_ptr<Mesh>> high_order_rounded_island_meshes;
  high_order_rounded_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, false, false, false, false, true)));
  LaplaceOperator high_order_rounded_island_laplace(rounded_island_iodata,
                                                    high_order_rounded_island_meshes);
  SurfaceResponseOperator high_order_rounded_island_response(
      rounded_island_iodata, high_order_rounded_island_laplace);
  CHECK(high_order_rounded_island_response.GetPatchCount() ==
        rounded_island_response.GetPatchCount());
  CHECK(high_order_rounded_island_response.GetBasisSize() ==
        rounded_island_response.GetBasisSize());
  CHECK_THAT(high_order_rounded_island_response.GetPatchWeight(),
             WithinRel(rounded_island_response.GetPatchWeight(), 1.0e-12));

  // Nearby curved conductors require a spatial cluster rather than independent corner
  // responses. Preflight must retain sampled face loops without an eight-vertex cap, and
  // its exact mask must round-trip through the process library on the same geometry.
  auto high_order_spatial_config = island_config;
  high_order_spatial_config["Boundaries"]["Terminal"] = {
      {{"Index", 1}, {"Attributes", {9}}}, {{"Index", 2}, {"Attributes", {10}}}};
  high_order_spatial_config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] =
      {9, 10};
  high_order_spatial_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      spatial_cluster_library_3d_path.string();
  high_order_spatial_config["Solver"]["Electrostatic"]["ResponseCorrection"]
                           ["UnmatchedPolicy"] = "Warn";
  IoData high_order_spatial_iodata(high_order_spatial_config, false);
  high_order_spatial_iodata.boundaries.cracked_attributes.insert(9);
  high_order_spatial_iodata.boundaries.cracked_attributes.insert(10);
  auto high_order_spatial_mesh = MakeIslandMesh(true, false, false, true, false, true);
  const auto high_order_spatial_requirements_path =
      temp.temp_dir / "surface-response-requirements-high-order-spatial.json";
  WriteSurfaceResponseRequirements(high_order_spatial_iodata, *high_order_spatial_mesh,
                                   high_order_spatial_requirements_path.string());
  std::ifstream high_order_spatial_requirements_input(high_order_spatial_requirements_path);
  REQUIRE(high_order_spatial_requirements_input);
  const auto high_order_spatial_requirements =
      json::parse(high_order_spatial_requirements_input);
  const auto high_order_spatial_requirement =
      std::find_if(high_order_spatial_requirements["Requirements"].begin(),
                   high_order_spatial_requirements["Requirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "SpatialEdgeCluster" &&
                            requirement["Geometry"].contains("PlanViewFacets") &&
                            requirement["Geometry"].contains("PlanViewBoundary");
                   });
  REQUIRE(high_order_spatial_requirement !=
          high_order_spatial_requirements["Requirements"].end());
  const auto &high_order_spatial_facets =
      (*high_order_spatial_requirement)["Geometry"]["PlanViewFacets"];
  CHECK(std::any_of(high_order_spatial_facets.begin(), high_order_spatial_facets.end(),
                    [](const auto &facet) { return facet["Points"].size() > 8; }));

  std::ifstream exact_high_order_spatial_library_input(spatial_cluster_library_3d_path);
  REQUIRE(exact_high_order_spatial_library_input);
  auto exact_high_order_spatial_library =
      json::parse(exact_high_order_spatial_library_input);
  auto &exact_high_order_spatial_model = exact_high_order_spatial_library["Models"].back();
  exact_high_order_spatial_model["Name"] = "high-order-curved-spatial-exact-mask";
  exact_high_order_spatial_model["Edges"] =
      (*high_order_spatial_requirement)["Geometry"]["Edges"];
  for (auto &edge : exact_high_order_spatial_model["Edges"])
  {
    edge["BoundaryCondition"] = "PEC";
  }
  exact_high_order_spatial_model["PlanViewBoundary"] =
      (*high_order_spatial_requirement)["Geometry"]["PlanViewBoundary"];
  const auto exact_high_order_spatial_library_path =
      temp.temp_dir / "fabrication-process-high-order-spatial-exact-mask-3d.json";
  std::ofstream exact_high_order_spatial_library_output(
      exact_high_order_spatial_library_path);
  exact_high_order_spatial_library_output << exact_high_order_spatial_library.dump(2)
                                          << "\n";
  exact_high_order_spatial_library_output.close();

  high_order_spatial_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      exact_high_order_spatial_library_path.string();
  IoData exact_high_order_spatial_iodata(high_order_spatial_config, false);
  exact_high_order_spatial_iodata.boundaries.cracked_attributes.insert(9);
  exact_high_order_spatial_iodata.boundaries.cracked_attributes.insert(10);
  const auto exact_high_order_spatial_requirements_path =
      temp.temp_dir / "surface-response-requirements-high-order-spatial-exact.json";
  WriteSurfaceResponseRequirements(exact_high_order_spatial_iodata,
                                   *high_order_spatial_mesh,
                                   exact_high_order_spatial_requirements_path.string());
  std::ifstream exact_high_order_spatial_requirements_input(
      exact_high_order_spatial_requirements_path);
  REQUIRE(exact_high_order_spatial_requirements_input);
  const auto exact_high_order_spatial_requirements =
      json::parse(exact_high_order_spatial_requirements_input);
  const auto matched_high_order_spatial =
      std::find_if(exact_high_order_spatial_requirements["Requirements"].begin(),
                   exact_high_order_spatial_requirements["Requirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "SpatialEdgeCluster" &&
                            requirement["Status"] == "Exact" &&
                            requirement["SelectedModels"][0]["Name"] ==
                                "high-order-curved-spatial-exact-mask";
                   });
  CHECK(matched_high_order_spatial !=
        exact_high_order_spatial_requirements["Requirements"].end());

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
  auto unqualified_interpolated_rounded_island_config = rounded_island_config;
  unqualified_interpolated_rounded_island_config
      ["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
          unqualified_interpolated_rounded_library_3d_path.string();
  IoData unqualified_interpolated_rounded_island_iodata(
      unqualified_interpolated_rounded_island_config, false);
  unqualified_interpolated_rounded_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> unqualified_interpolated_rounded_island_meshes;
  unqualified_interpolated_rounded_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true)));
  LaplaceOperator unqualified_interpolated_rounded_island_laplace(
      unqualified_interpolated_rounded_island_iodata,
      unqualified_interpolated_rounded_island_meshes);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(unqualified_interpolated_rounded_island_iodata,
                              unqualified_interpolated_rounded_island_laplace),
      Catch::Matchers::ContainsSubstring(
          "Automatic fabrication-process response matching failed"));
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
  CHECK(constant_island_response.closure_independent_confident);
  CHECK(constant_island_response.maximum_trace_closure_spread > 0.05);
  CHECK(constant_island_response.response_weighted_trace_closure_spread > 0.0);
  CHECK(constant_island_response.trace_closure_response_failure_fraction > 0.0);
  CHECK_FALSE(constant_island_response.confident);

  auto impedance_convex_maxwell_config = convex_maxwell_island_config;
  impedance_convex_maxwell_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3,
                                                                           4, 5, 6};
  impedance_convex_maxwell_config["Boundaries"]["Impedance"] = {
      {{"Attributes", {9}}, {"Ls", 1.0e-13}}};
  impedance_convex_maxwell_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      finite_impedance_convex_library_3d_path.string();
  IoData impedance_convex_maxwell_iodata(impedance_convex_maxwell_config, false);
  impedance_convex_maxwell_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> impedance_convex_maxwell_meshes;
  impedance_convex_maxwell_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
  SpaceOperator impedance_convex_maxwell_space(impedance_convex_maxwell_iodata,
                                               impedance_convex_maxwell_meshes);
  SurfaceResponseOperator impedance_convex_maxwell_response(impedance_convex_maxwell_iodata,
                                                            impedance_convex_maxwell_space);
  CHECK(impedance_convex_maxwell_response.GetPatchCount() ==
        convex_maxwell_island_response.GetPatchCount());
  GridFunction impedance_convex_maxwell_field(impedance_convex_maxwell_space.GetNDSpace(),
                                              true);
  impedance_convex_maxwell_field.Real().ProjectCoefficient(field_coefficient);
  impedance_convex_maxwell_field.Imag() = 0.0;
  const auto impedance_convex_maxwell_result =
      impedance_convex_maxwell_response.GetMaxwellResponse(impedance_convex_maxwell_field,
                                                           0.0);
  CHECK(impedance_convex_maxwell_result.loop_residual < 1.0e-10);
  CHECK_THAT(impedance_convex_maxwell_result.matched_length_fraction,
             WithinAbs(1.0, 1.0e-12));
  CHECK(impedance_convex_maxwell_result.corner_neighborhood_fraction == 0.0);
  CHECK(impedance_convex_maxwell_result.boundary_law_verified);

  auto junction_maxwell_config = convex_maxwell_island_config;
  junction_maxwell_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      junction_library_3d_path.string();
  IoData junction_maxwell_iodata(junction_maxwell_config, false);
  junction_maxwell_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> junction_maxwell_meshes;
  junction_maxwell_meshes.push_back(std::make_unique<Mesh>(MakeTouchingIslandMesh()));
  SpaceOperator junction_maxwell_space(junction_maxwell_iodata, junction_maxwell_meshes);
  SurfaceResponseOperator junction_maxwell_response(junction_maxwell_iodata,
                                                    junction_maxwell_space);
  CHECK(junction_maxwell_response.GetPatchCount() ==
        junction_island_response.GetPatchCount());
  GridFunction junction_maxwell_field(junction_maxwell_space.GetNDSpace(), true);
  junction_maxwell_field.Real().ProjectCoefficient(field_coefficient);
  junction_maxwell_field.Imag() = 0.0;
  const auto junction_maxwell_result =
      junction_maxwell_response.GetMaxwellResponse(junction_maxwell_field, 0.0);
  CHECK(junction_maxwell_result.loop_residual < 1.0e-10);
  CHECK(junction_maxwell_result.corner_neighborhood_fraction == 0.0);

  const auto junction_requirements_path =
      temp.temp_dir / "surface-response-requirements-junction.json";
  WriteSurfaceResponseRequirements(junction_maxwell_iodata, *junction_maxwell_meshes.back(),
                                   junction_requirements_path.string());
  std::ifstream junction_requirements_input(junction_requirements_path);
  REQUIRE(junction_requirements_input);
  const json junction_requirements = json::parse(junction_requirements_input);
  const auto junction_requirement =
      std::find_if(junction_requirements["Requirements"].begin(),
                   junction_requirements["Requirements"].end(), [](const auto &requirement)
                   { return requirement["Topology"] == "Junction"; });
  REQUIRE(junction_requirement != junction_requirements["Requirements"].end());
  const auto &junction_geometry = (*junction_requirement)["Geometry"];
  CHECK((*junction_requirement)["Status"] == "Missing");
  CHECK(junction_geometry["SignatureVersion"] == 2);
  CHECK(junction_geometry["ArmCount"].get<int>() >= 3);
  CHECK(junction_geometry["Arms"].size() ==
        junction_geometry["ArmCount"].get<std::size_t>());
  for (const auto &arm : junction_geometry["Arms"])
  {
    CHECK(arm["Conductor"] == 1);
    CHECK(arm["InterfaceSlot"] == 0);
    CHECK(arm["BoundaryCondition"]["Type"] == "PEC");
    CHECK(arm["Direction"].size() == 3);
    CHECK(arm["GapDirection"].size() == 3);
    CHECK(arm["ProcessNormal"].size() == 3);
  }
  REQUIRE(junction_geometry.contains("PlanViewFacets"));
  REQUIRE(!junction_geometry["PlanViewFacets"].empty());
  REQUIRE(junction_geometry.contains("PlanViewBoundary"));
  for (const auto &component : junction_geometry["PlanViewBoundary"])
  {
    REQUIRE(component.contains("ContinuationSegments"));
  }
  std::ifstream exact_junction_library_input(junction_library_3d_path);
  REQUIRE(exact_junction_library_input);
  auto exact_junction_library = json::parse(exact_junction_library_input);
  auto exact_junction_model = *std::find_if(
      exact_junction_library["Models"].begin(), exact_junction_library["Models"].end(),
      [](const auto &model)
      {
        return model["Topology"] == "Junction" &&
               model.value("BoundaryCondition", json("PEC")) == "PEC";
      });
  exact_junction_model["Name"] = "junction-4x90-exact-mask";
  exact_junction_model["PlanViewBoundary"] = junction_geometry["PlanViewBoundary"];
  exact_junction_library["Models"].push_back(std::move(exact_junction_model));
  const auto exact_junction_library_path =
      temp.temp_dir / "surface-process-junction-exact-mask-3d.json";
  std::ofstream exact_junction_library_output(exact_junction_library_path);
  exact_junction_library_output << exact_junction_library.dump(2) << "\n";
  exact_junction_library_output.close();
  auto exact_junction_config = junction_maxwell_config;
  exact_junction_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      exact_junction_library_path.string();
  IoData exact_junction_iodata(exact_junction_config, false);
  exact_junction_iodata.boundaries.cracked_attributes.insert(9);
  const auto exact_junction_requirements_path =
      temp.temp_dir / "surface-response-requirements-junction-exact-mask.json";
  WriteSurfaceResponseRequirements(exact_junction_iodata, *junction_maxwell_meshes.back(),
                                   exact_junction_requirements_path.string());
  std::ifstream exact_junction_requirements_input(exact_junction_requirements_path);
  REQUIRE(exact_junction_requirements_input);
  const auto exact_junction_requirements = json::parse(exact_junction_requirements_input);
  const auto matched_exact_junction = std::find_if(
      exact_junction_requirements["Requirements"].begin(),
      exact_junction_requirements["Requirements"].end(),
      [](const auto &requirement)
      {
        return requirement["Topology"] == "Junction" &&
               requirement["SelectedModels"][0]["Name"] == "junction-4x90-exact-mask";
      });
  CHECK(matched_exact_junction != exact_junction_requirements["Requirements"].end());

  auto impedance_junction_config = junction_maxwell_config;
  impedance_junction_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6};
  impedance_junction_config["Boundaries"]["Impedance"] = {
      {{"Attributes", {9}}, {"Ls", 1.0e-13}}};
  IoData impedance_junction_iodata(impedance_junction_config, false);
  impedance_junction_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> impedance_junction_meshes;
  impedance_junction_meshes.push_back(std::make_unique<Mesh>(MakeTouchingIslandMesh()));
  SpaceOperator impedance_junction_space(impedance_junction_iodata,
                                         impedance_junction_meshes);
  SurfaceResponseOperator impedance_junction_response(impedance_junction_iodata,
                                                      impedance_junction_space);
  CHECK(impedance_junction_response.GetPatchCount() ==
        static_cast<int>(touching_segments.size()) * island_line_rule.GetNPoints() -
            removed_junction_straight_patches + touching_junctions);
  CHECK(impedance_junction_response.GetBasisSize() ==
        4 * (static_cast<int>(touching_segments.size()) * island_line_rule.GetNPoints() -
             removed_junction_straight_patches) +
            12 * touching_junctions);
  GridFunction impedance_junction_field(impedance_junction_space.GetNDSpace(), true);
  impedance_junction_field.Real().ProjectCoefficient(field_coefficient);
  impedance_junction_field.Imag() = 0.0;
  const auto impedance_junction_result =
      impedance_junction_response.GetMaxwellResponse(impedance_junction_field, 0.0);
  CHECK(impedance_junction_result.loop_residual < 1.0e-10);

  // Two disconnected island corners separated diagonally by less than 2R form one
  // localized four-edge neighborhood. It contains perpendicular and endpoint-adjacent
  // interactions, so no longitudinal pair or parallel-cluster coupon can represent it.
  auto spatial_cluster_config = convex_maxwell_island_config;
  spatial_cluster_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9, 10};
  spatial_cluster_config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {
      9, 10};
  spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      spatial_cluster_library_3d_path.string();
  IoData spatial_cluster_iodata(spatial_cluster_config, false);
  spatial_cluster_iodata.boundaries.cracked_attributes.insert(9);
  spatial_cluster_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> spatial_cluster_meshes;
  spatial_cluster_meshes.push_back(std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
  SpaceOperator spatial_cluster_space(spatial_cluster_iodata, spatial_cluster_meshes);
  SurfaceResponseOperator spatial_cluster_response(spatial_cluster_iodata,
                                                   spatial_cluster_space);
  CHECK(spatial_cluster_response.GetPatchCount() > 0);
  CHECK(spatial_cluster_response.GetTargetInterfaces() == std::set<int>{4});
  CHECK_THAT(spatial_cluster_response.GetPatchWeight(), WithinRel(17.25, 1.0e-12));

  GridFunction spatial_cluster_field(spatial_cluster_space.GetNDSpace(), true);
  spatial_cluster_field.Real().ProjectCoefficient(field_coefficient);
  spatial_cluster_field.Imag() = 0.0;
  const auto spatial_cluster_result =
      spatial_cluster_response.GetMaxwellResponse(spatial_cluster_field, 0.0);
  CHECK(std::abs(spatial_cluster_result.domain_correction) > 0.0);
  CHECK(spatial_cluster_result.loop_residual < 1.0e-10);
  CHECK_THAT(spatial_cluster_result.matched_length_fraction, WithinAbs(1.0, 1.0e-12));
  CHECK(spatial_cluster_result.corner_neighborhood_fraction == 0.0);

  auto missing_spatial_cluster_config = spatial_cluster_config;
  missing_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      convex_library_3d_path.string();
  missing_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] =
      "Warn";
  IoData missing_spatial_cluster_iodata(missing_spatial_cluster_config, false);
  missing_spatial_cluster_iodata.boundaries.cracked_attributes.insert(9);
  missing_spatial_cluster_iodata.boundaries.cracked_attributes.insert(10);
  const auto spatial_requirements_path =
      temp.temp_dir / "surface-response-requirements-spatial.json";
  WriteSurfaceResponseRequirements(missing_spatial_cluster_iodata,
                                   *spatial_cluster_meshes.back(),
                                   spatial_requirements_path.string());
  std::ifstream spatial_requirements_input(spatial_requirements_path);
  REQUIRE(spatial_requirements_input);
  const json spatial_requirements = json::parse(spatial_requirements_input);
  const auto spatial_requirement =
      std::find_if(spatial_requirements["Requirements"].begin(),
                   spatial_requirements["Requirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "SpatialEdgeCluster" &&
                            requirement["Status"] == "Missing" &&
                            requirement["Geometry"].contains("Edges");
                   });
  REQUIRE(spatial_requirement != spatial_requirements["Requirements"].end());
  const auto &spatial_edges = (*spatial_requirement)["Geometry"]["Edges"];
  REQUIRE(spatial_edges.size() >= 2);
  std::set<int> spatial_conductors;
  for (const auto &edge : spatial_edges)
  {
    const int conductor = edge["Conductor"];
    CHECK(conductor > 0);
    spatial_conductors.insert(conductor);
    CHECK(edge["Point"].size() == 3);
    CHECK(edge["GapDirection"].size() == 3);
    CHECK(edge["ProcessNormal"].size() == 3);
    REQUIRE(edge["Interval"].size() == 2);
    CHECK(std::isfinite(edge["Interval"][0].get<double>()));
    CHECK(std::isfinite(edge["Interval"][1].get<double>()));
    CHECK(edge["Interval"][0].get<double>() <= 0.0);
    CHECK(edge["Interval"][1].get<double>() >= 0.0);
    CHECK(edge["Interval"][1].get<double>() > edge["Interval"][0].get<double>());
    CHECK(edge["InterfaceSlot"].get<int>() >= 0);
    CHECK(edge["BoundaryCondition"]["Type"] == "PEC");
  }
  CHECK(*spatial_conductors.begin() == 1);
  REQUIRE((*spatial_requirement)["Geometry"].contains("PlanViewFacets"));
  const auto &spatial_facets = (*spatial_requirement)["Geometry"]["PlanViewFacets"];
  REQUIRE(!spatial_facets.empty());
  std::set<int> facet_conductors;
  for (const auto &facet : spatial_facets)
  {
    facet_conductors.insert(facet["Conductor"].get<int>());
    REQUIRE(facet["Points"].size() >= 3);
    for (const auto &point : facet["Points"])
    {
      REQUIRE(point.size() == 3);
      CHECK(std::isfinite(point[0].get<double>()));
      CHECK(std::isfinite(point[1].get<double>()));
      CHECK(std::isfinite(point[2].get<double>()));
    }
  }
  CHECK(facet_conductors == spatial_conductors);
  REQUIRE((*spatial_requirement)["Geometry"].contains("PlanViewBoundary"));
  const auto &plan_view_boundary = (*spatial_requirement)["Geometry"]["PlanViewBoundary"];
  REQUIRE(plan_view_boundary.is_array());
  REQUIRE(!plan_view_boundary.empty());

  std::ifstream empty_spatial_library_input(convex_library_3d_path);
  REQUIRE(empty_spatial_library_input);
  auto empty_spatial_library = json::parse(empty_spatial_library_input);
  empty_spatial_library["Name"] = "unit-test-empty-spatial-process";
  empty_spatial_library["Models"] = json::array();
  const auto empty_spatial_library_path =
      temp.temp_dir / "surface-process-empty-spatial-3d.json";
  std::ofstream empty_spatial_library_output(empty_spatial_library_path);
  empty_spatial_library_output << empty_spatial_library.dump(2) << "\n";
  empty_spatial_library_output.close();
  auto empty_spatial_config = missing_spatial_cluster_config;
  empty_spatial_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      empty_spatial_library_path.string();
  IoData empty_spatial_iodata(empty_spatial_config, false);
  empty_spatial_iodata.boundaries.cracked_attributes.insert(9);
  empty_spatial_iodata.boundaries.cracked_attributes.insert(10);
  const auto empty_spatial_requirements_path =
      temp.temp_dir / "surface-response-requirements-empty-spatial.json";
  WriteSurfaceResponseRequirements(empty_spatial_iodata, *spatial_cluster_meshes.back(),
                                   empty_spatial_requirements_path.string());
  std::ifstream empty_spatial_requirements_input(empty_spatial_requirements_path);
  REQUIRE(empty_spatial_requirements_input);
  const json empty_spatial_requirements = json::parse(empty_spatial_requirements_input);
  CHECK_FALSE(empty_spatial_requirements["Complete"]);
  CHECK(empty_spatial_requirements["Summary"]["Counts"]["Exact"] == 0);
  const auto empty_spatial_requirement =
      std::find_if(empty_spatial_requirements["Requirements"].begin(),
                   empty_spatial_requirements["Requirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "SpatialEdgeCluster" &&
                            requirement["Status"] == "Missing" &&
                            requirement["Geometry"].contains("PlanViewBoundary");
                   });
  REQUIRE(empty_spatial_requirement != empty_spatial_requirements["Requirements"].end());
  CHECK((*empty_spatial_requirement)["Geometry"].contains("PlanViewFacets"));

  std::ifstream spatial_cluster_library_input(spatial_cluster_library_3d_path);
  REQUIRE(spatial_cluster_library_input);
  auto exact_mask_library = json::parse(spatial_cluster_library_input);
  exact_mask_library["Name"] = "unit-test-process-spatial-cluster-exact-mask-3d";
  auto &exact_mask_model = exact_mask_library["Models"].back();
  exact_mask_model["Name"] = "offset-corner-pair-exact-mask";
  exact_mask_model["Edges"] = spatial_edges;
  for (auto &edge : exact_mask_model["Edges"])
  {
    edge["BoundaryCondition"] = "PEC";
  }
  exact_mask_model["PlanViewBoundary"] = plan_view_boundary;
  const auto exact_mask_library_path =
      temp.temp_dir / "surface-process-spatial-cluster-exact-mask-3d.json";
  std::ofstream exact_mask_output(exact_mask_library_path);
  exact_mask_output << exact_mask_library.dump(2) << "\n";
  exact_mask_output.close();

  auto exact_mask_config = spatial_cluster_config;
  exact_mask_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      exact_mask_library_path.string();
  IoData exact_mask_iodata(exact_mask_config, false);
  exact_mask_iodata.boundaries.cracked_attributes.insert(9);
  exact_mask_iodata.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator exact_mask_response(exact_mask_iodata, spatial_cluster_space);
  CHECK(exact_mask_response.GetPatchCount() == spatial_cluster_response.GetPatchCount());
  CHECK_THAT(exact_mask_response.GetPatchWeight(),
             WithinRel(spatial_cluster_response.GetPatchWeight(), 1.0e-12));

  auto CheckPlanViewMaskStatus = [&](const fs::path &library_path,
                                     std::string_view expected_status,
                                     std::string_view suffix)
  {
    auto config = spatial_cluster_config;
    config["Solver"]["SurfaceResponseCorrection"]["Library"] = library_path.string();
    config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] = "Warn";
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    iodata.boundaries.cracked_attributes.insert(10);
    const auto path = temp.temp_dir / ("surface-response-requirements-mask-" +
                                       std::string(suffix) + ".json");
    WriteSurfaceResponseRequirements(iodata, *spatial_cluster_meshes.back(), path.string());
    std::ifstream input(path);
    REQUIRE(input);
    const json manifest = json::parse(input);
    const auto requirement =
        std::find_if(manifest["Requirements"].begin(), manifest["Requirements"].end(),
                     [&](const auto &entry)
                     {
                       return entry["Topology"] == "SpatialEdgeCluster" &&
                              entry["Geometry"].contains("Edges") &&
                              entry["Geometry"]["Edges"].size() == spatial_edges.size();
                     });
    REQUIRE(requirement != manifest["Requirements"].end());
    CHECK((*requirement)["Status"] == expected_status);
  };
  CheckPlanViewMaskStatus(exact_mask_library_path, "Exact", "exact");

  auto mismatched_mask_library = exact_mask_library;
  mismatched_mask_library["Name"] = "unit-test-process-spatial-cluster-mismatched-mask-3d";
  auto &mismatched_coordinate =
      mismatched_mask_library["Models"].back()["PlanViewBoundary"][0]["Segments"][0][0][0];
  mismatched_coordinate = mismatched_coordinate.get<long long int>() + 1;
  const auto mismatched_mask_library_path =
      temp.temp_dir / "surface-process-spatial-cluster-mismatched-mask-3d.json";
  std::ofstream mismatched_mask_output(mismatched_mask_library_path);
  mismatched_mask_output << mismatched_mask_library.dump(2) << "\n";
  mismatched_mask_output.close();
  CheckPlanViewMaskStatus(mismatched_mask_library_path, "Missing", "mismatched");

  // Spatial clusters match each physical edge's complete boundary law. The second
  // island is finite impedance while the first remains PEC; its straight edges,
  // corners, and coupled four-edge neighborhood all have exact object-form models.
  auto mixed_impedance_spatial_cluster_config = spatial_cluster_config;
  mixed_impedance_spatial_cluster_config["Boundaries"]["Ground"]["Attributes"] = {
      1, 2, 3, 4, 5, 6, 9};
  mixed_impedance_spatial_cluster_config["Boundaries"]["Impedance"] = {
      {{"Attributes", {10}}, {"Ls", 2.0e-13}}};
  mixed_impedance_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      spatial_cluster_mixed_impedance_library_3d_path.string();
  IoData mixed_impedance_spatial_cluster_iodata(mixed_impedance_spatial_cluster_config,
                                                false);
  mixed_impedance_spatial_cluster_iodata.boundaries.cracked_attributes.insert(9);
  mixed_impedance_spatial_cluster_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> mixed_impedance_spatial_cluster_meshes;
  mixed_impedance_spatial_cluster_meshes.push_back(
      std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
  SpaceOperator mixed_impedance_spatial_cluster_space(
      mixed_impedance_spatial_cluster_iodata, mixed_impedance_spatial_cluster_meshes);
  SurfaceResponseOperator mixed_impedance_spatial_cluster_response(
      mixed_impedance_spatial_cluster_iodata, mixed_impedance_spatial_cluster_space);
  GridFunction mixed_impedance_spatial_cluster_field(
      mixed_impedance_spatial_cluster_space.GetNDSpace(), true);
  mixed_impedance_spatial_cluster_field.Real().ProjectCoefficient(field_coefficient);
  mixed_impedance_spatial_cluster_field.Imag() = 0.0;
  const auto mixed_impedance_spatial_cluster_result =
      mixed_impedance_spatial_cluster_response.GetMaxwellResponse(
          mixed_impedance_spatial_cluster_field, 0.0);
  CHECK_THAT(mixed_impedance_spatial_cluster_result.matched_length_fraction,
             WithinAbs(1.0, 1.0e-12));
  CHECK(mixed_impedance_spatial_cluster_result.corner_neighborhood_fraction == 0.0);
  CHECK(mixed_impedance_spatial_cluster_result.boundary_law_verified);

  auto parameter_mismatch_spatial_cluster_config = mixed_impedance_spatial_cluster_config;
  parameter_mismatch_spatial_cluster_config
      ["Solver"]["SurfaceResponseCorrection"]["Library"] =
          spatial_cluster_parameter_mismatch_library_3d_path.string();
  parameter_mismatch_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]
                                           ["UnmatchedPolicy"] = "Warn";
  IoData parameter_mismatch_spatial_cluster_iodata(
      parameter_mismatch_spatial_cluster_config, false);
  parameter_mismatch_spatial_cluster_iodata.boundaries.cracked_attributes.insert(9);
  parameter_mismatch_spatial_cluster_iodata.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator parameter_mismatch_spatial_cluster_response(
      parameter_mismatch_spatial_cluster_iodata, mixed_impedance_spatial_cluster_space);
  CHECK(parameter_mismatch_spatial_cluster_response.GetPatchCount() > 0);
  CHECK(parameter_mismatch_spatial_cluster_response.GetPatchCount() <
        mixed_impedance_spatial_cluster_response.GetPatchCount());
  const auto parameter_mismatch_spatial_cluster_result =
      parameter_mismatch_spatial_cluster_response.GetMaxwellResponse(
          mixed_impedance_spatial_cluster_field, 0.0);
  CHECK(parameter_mismatch_spatial_cluster_result.matched_length_fraction > 0.0);
  CHECK(parameter_mismatch_spatial_cluster_result.matched_length_fraction < 1.0);

  auto cross_layer_spatial_cluster_config = spatial_cluster_config;
  cross_layer_spatial_cluster_config["Boundaries"]["Postprocessing"]["Dielectric"] = {
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
       {"Type", "SA"},
       {"Thickness", 0.002},
       {"Permittivity", 4.0},
       {"AutomaticEdges", true},
       {"EdgeDistances", {0.2}},
       {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}};
  cross_layer_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      cross_layer_spatial_cluster_library_3d_path.string();
  cross_layer_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]
                                    ["TargetInterfaces"] = {4, 5};
  IoData cross_layer_spatial_cluster_iodata(cross_layer_spatial_cluster_config, false);
  cross_layer_spatial_cluster_iodata.boundaries.cracked_attributes.insert(9);
  cross_layer_spatial_cluster_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> cross_layer_spatial_cluster_meshes;
  cross_layer_spatial_cluster_meshes.push_back(
      std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
  SpaceOperator cross_layer_spatial_cluster_space(cross_layer_spatial_cluster_iodata,
                                                  cross_layer_spatial_cluster_meshes);
  SurfaceResponseOperator cross_layer_spatial_cluster_response(
      cross_layer_spatial_cluster_iodata, cross_layer_spatial_cluster_space);
  CHECK(cross_layer_spatial_cluster_response.GetPatchCount() > 0);
  CHECK(cross_layer_spatial_cluster_response.GetTargetInterfaces() == std::set<int>{4, 5});

  GridFunction cross_layer_spatial_cluster_field(
      cross_layer_spatial_cluster_space.GetNDSpace(), true);
  cross_layer_spatial_cluster_field.Real().ProjectCoefficient(field_coefficient);
  cross_layer_spatial_cluster_field.Imag() = 0.0;
  const auto cross_layer_spatial_cluster_result =
      cross_layer_spatial_cluster_response.GetMaxwellResponse(
          cross_layer_spatial_cluster_field, 0.0);
  CHECK(cross_layer_spatial_cluster_result.loop_residual < 1.0e-10);
  CHECK_THAT(cross_layer_spatial_cluster_result.matched_length_fraction,
             WithinAbs(1.0, 1.0e-12));
  CHECK(cross_layer_spatial_cluster_result.corner_neighborhood_fraction == 0.0);
  CHECK(cross_layer_spatial_cluster_result.fabricated_surface_energy.count(4) == 1);
  CHECK(cross_layer_spatial_cluster_result.fabricated_surface_energy.count(5) == 1);

  const auto cross_layer_requirements_path =
      temp.temp_dir / "surface-response-requirements-cross-layer.json";
  WriteSurfaceResponseRequirements(cross_layer_spatial_cluster_iodata,
                                   *cross_layer_spatial_cluster_meshes.back(),
                                   cross_layer_requirements_path.string());
  std::ifstream cross_layer_requirements_input(cross_layer_requirements_path);
  REQUIRE(cross_layer_requirements_input);
  const json cross_layer_requirements = json::parse(cross_layer_requirements_input);
  const auto cross_layer_requirement =
      std::find_if(cross_layer_requirements["Requirements"].begin(),
                   cross_layer_requirements["Requirements"].end(),
                   [](const auto &requirement)
                   {
                     if (requirement["Topology"] != "SpatialEdgeCluster" ||
                         !requirement["Geometry"].contains("Edges"))
                     {
                       return false;
                     }
                     std::set<int> slots;
                     for (const auto &edge : requirement["Geometry"]["Edges"])
                     {
                       slots.insert(edge["InterfaceSlot"].template get<int>());
                     }
                     return slots.size() == 2;
                   });
  REQUIRE(cross_layer_requirement != cross_layer_requirements["Requirements"].end());
  std::set<int> cross_layer_slots;
  for (const auto &edge : (*cross_layer_requirement)["Geometry"]["Edges"])
  {
    cross_layer_slots.insert(edge["InterfaceSlot"].get<int>());
    CHECK(edge["Conductor"].get<int>() > 0);
    CHECK(edge["BoundaryCondition"]["Type"] == "PEC");
  }
  CHECK(cross_layer_slots.size() == 2);

  // Cross-layer replacement is all-or-nothing. A library with no multi-slot model, or
  // with a multi-slot model that omits one physical target type, must leave the
  // interaction neighborhood unmatched instead of applying independent local coupons.
  auto CheckCrossLayerSpatialClusterMismatch = [&](const fs::path &library)
  {
    CAPTURE(library.string());
    auto mismatch_config = cross_layer_spatial_cluster_config;
    mismatch_config["Solver"]["SurfaceResponseCorrection"]["Library"] = library.string();
    mismatch_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] = "Warn";
    IoData mismatch_iodata(mismatch_config, false);
    mismatch_iodata.boundaries.cracked_attributes.insert(9);
    mismatch_iodata.boundaries.cracked_attributes.insert(10);
    SurfaceResponseOperator mismatch_response(mismatch_iodata,
                                              cross_layer_spatial_cluster_space);
    CHECK(mismatch_response.GetPatchCount() > 0);
    CHECK(mismatch_response.GetPatchCount() <
          cross_layer_spatial_cluster_response.GetPatchCount());
    const auto result =
        mismatch_response.GetMaxwellResponse(cross_layer_spatial_cluster_field, 0.0);
    CHECK(result.matched_length_fraction > 0.0);
    CHECK(result.matched_length_fraction < 1.0);
  };
  CheckCrossLayerSpatialClusterMismatch(spatial_cluster_library_3d_path);
  CheckCrossLayerSpatialClusterMismatch(
      incomplete_cross_layer_spatial_cluster_library_3d_path);

  // A spatial cluster is an exact local topology signature. Perturbing one edge's
  // position or orientation, changing its metal boundary condition, or presenting one
  // more physical edge than the coupon describes must omit the interaction neighborhood.
  auto CheckSpatialClusterMismatch = [&](const fs::path &library)
  {
    CAPTURE(library.string());
    auto mismatch_config = spatial_cluster_config;
    mismatch_config["Solver"]["SurfaceResponseCorrection"]["Library"] = library.string();
    mismatch_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] = "Warn";
    IoData mismatch_iodata(mismatch_config, false);
    mismatch_iodata.boundaries.cracked_attributes.insert(9);
    mismatch_iodata.boundaries.cracked_attributes.insert(10);
    SurfaceResponseOperator mismatch_response(mismatch_iodata, spatial_cluster_space);
    CHECK(mismatch_response.GetPatchCount() > 0);
    CHECK(mismatch_response.GetPatchCount() < spatial_cluster_response.GetPatchCount());
    const auto result = mismatch_response.GetMaxwellResponse(spatial_cluster_field, 0.0);
    CHECK(result.matched_length_fraction > 0.0);
    CHECK(result.matched_length_fraction < 1.0);
  };
  CheckSpatialClusterMismatch(spatial_cluster_position_mismatch_library_3d_path);
  CheckSpatialClusterMismatch(spatial_cluster_orientation_mismatch_library_3d_path);
  CheckSpatialClusterMismatch(spatial_cluster_interval_mismatch_library_3d_path);
  CheckSpatialClusterMismatch(spatial_cluster_extra_edge_library_3d_path);
  CheckSpatialClusterMismatch(spatial_cluster_impedance_mismatch_library_3d_path);

  // Separated fabrication planes are independent placements of the same local process.
  // They should both match unless their radius-R neighborhoods actually interact.
  auto multilayer_maxwell_config = convex_maxwell_island_config;
  multilayer_maxwell_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4,
                                                                     5, 6, 9, 10};
  auto second_interface =
      multilayer_maxwell_config["Boundaries"]["Postprocessing"]["Dielectric"][0];
  second_interface["Index"] = 5;
  second_interface["Attributes"] = {10};
  multilayer_maxwell_config["Boundaries"]["Postprocessing"]["Dielectric"].push_back(
      second_interface);
  multilayer_maxwell_config["Solver"]["SurfaceResponseCorrection"]["TargetInterfaces"] = {
      4, 5};
  IoData multilayer_maxwell_iodata(multilayer_maxwell_config, false);
  multilayer_maxwell_iodata.boundaries.cracked_attributes.insert(9);
  multilayer_maxwell_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> multilayer_maxwell_meshes;
  multilayer_maxwell_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(false, false, false, false, true)));
  SpaceOperator multilayer_maxwell_space(multilayer_maxwell_iodata,
                                         multilayer_maxwell_meshes);
  SurfaceResponseOperator multilayer_maxwell_response(multilayer_maxwell_iodata,
                                                      multilayer_maxwell_space);
  CHECK(multilayer_maxwell_response.GetPatchCount() ==
        2 * convex_maxwell_island_response.GetPatchCount());
  CHECK(multilayer_maxwell_response.GetTargetInterfaces() == std::set<int>{4, 5});
  GridFunction multilayer_field(multilayer_maxwell_space.GetNDSpace(), true);
  multilayer_field.Real().ProjectCoefficient(field_coefficient);
  multilayer_field.Imag() = 0.0;
  const auto multilayer_result =
      multilayer_maxwell_response.GetMaxwellResponse(multilayer_field, 0.0);
  CHECK_THAT(multilayer_result.matched_length_fraction, WithinAbs(1.0, 1.0e-12));
  CHECK(multilayer_result.fabricated_surface_energy.count(4) == 1);
  CHECK(multilayer_result.fabricated_surface_energy.count(5) == 1);

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

  // Fixed-flux closure acts only on the free trace subspace. Matrix rows associated
  // with exact PEC trace knots are calibration artifacts and must not change either
  // closure when their free-free blocks are unchanged.
  auto constrained_perturbed_rounded_config = rounded_maxwell_island_config;
  constrained_perturbed_rounded_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      constrained_perturbed_rounded_library_3d_path.string();
  IoData constrained_perturbed_rounded_iodata(constrained_perturbed_rounded_config, false);
  constrained_perturbed_rounded_iodata.boundaries.cracked_attributes.insert(9);
  SurfaceResponseOperator constrained_perturbed_rounded_response(
      constrained_perturbed_rounded_iodata, rounded_maxwell_island_space);
  const auto constrained_perturbed_rounded_result =
      constrained_perturbed_rounded_response.GetMaxwellResponse(rounded_island_field, 0.0);
  CHECK_THAT(constrained_perturbed_rounded_result.domain_correction,
             WithinRel(rounded_island_result.domain_correction, 1.0e-12));
  CHECK_THAT(constrained_perturbed_rounded_result.domain_correction_fixed_flux,
             WithinRel(rounded_island_result.domain_correction_fixed_flux, 1.0e-12));
  CHECK_THAT(constrained_perturbed_rounded_result.fabricated_surface_energy.at(4),
             WithinRel(rounded_island_result.fabricated_surface_energy.at(4), 1.0e-12));
  CHECK_THAT(
      constrained_perturbed_rounded_result.fabricated_surface_energy_fixed_flux.at(4),
      WithinRel(rounded_island_result.fabricated_surface_energy_fixed_flux.at(4), 1.0e-12));
  CHECK_THAT(
      constrained_perturbed_rounded_result.response_weighted_trace_closure_spread,
      WithinRel(rounded_island_result.response_weighted_trace_closure_spread, 1.0e-12));
  CHECK_THAT(
      constrained_perturbed_rounded_result.trace_closure_response_failure_fraction,
      WithinAbs(rounded_island_result.trace_closure_response_failure_fraction, 1.0e-12));

  auto impedance_rounded_maxwell_config = rounded_maxwell_island_config;
  impedance_rounded_maxwell_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3,
                                                                            4, 5, 6};
  impedance_rounded_maxwell_config["Boundaries"]["Impedance"] = {
      {{"Attributes", {9}}, {"Ls", 1.0e-13}}};
  impedance_rounded_maxwell_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      finite_impedance_rounded_library_3d_path.string();
  IoData impedance_rounded_maxwell_iodata(impedance_rounded_maxwell_config, false);
  impedance_rounded_maxwell_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> impedance_rounded_maxwell_meshes;
  impedance_rounded_maxwell_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, true)));
  SpaceOperator impedance_rounded_maxwell_space(impedance_rounded_maxwell_iodata,
                                                impedance_rounded_maxwell_meshes);
  SurfaceResponseOperator impedance_rounded_maxwell_response(
      impedance_rounded_maxwell_iodata, impedance_rounded_maxwell_space);
  CHECK(impedance_rounded_maxwell_response.GetPatchCount() ==
        rounded_maxwell_island_response.GetPatchCount());
  GridFunction impedance_rounded_maxwell_field(impedance_rounded_maxwell_space.GetNDSpace(),
                                               true);
  impedance_rounded_maxwell_field.Real().ProjectCoefficient(field_coefficient);
  impedance_rounded_maxwell_field.Imag() = 0.0;
  const auto impedance_rounded_maxwell_result =
      impedance_rounded_maxwell_response.GetMaxwellResponse(impedance_rounded_maxwell_field,
                                                            0.0);
  CHECK(impedance_rounded_maxwell_result.loop_residual < 1.0e-10);
  CHECK_THAT(impedance_rounded_maxwell_result.matched_length_fraction,
             WithinAbs(1.0, 1.0e-12));
  CHECK(impedance_rounded_maxwell_result.corner_neighborhood_fraction == 0.0);
  CHECK_FALSE(impedance_rounded_maxwell_result.boundary_law_verified);
  CHECK_FALSE(impedance_rounded_maxwell_result.closure_independent_confident);

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
  SurfaceResponseOperator rounded_concave_maxwell_response(rounded_concave_maxwell_iodata,
                                                           rounded_concave_maxwell_space);
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
  mixed_signature_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9, 10};
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
  mixed_signature_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] = "Warn";
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
  local_interaction_config["Solver"]["SurfaceResponseCorrection"]["TargetInterfaces"] = {4};
  IoData local_interaction_iodata(local_interaction_config, false);
  local_interaction_iodata.boundaries.cracked_attributes.insert(9);
  local_interaction_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> local_interaction_meshes;
  local_interaction_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(false, false, false, true)));
  SpaceOperator local_interaction_space(local_interaction_iodata, local_interaction_meshes);
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
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(24, 4, 24, mfem::Element::HEXAHEDRON,
                                                    extent, 1.0, extent);
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
      const bool left_aperture = xmin >= 0.75 - 1.0e-12 && xmax <= 1.25 + 1.0e-12 &&
                                 zmin >= 0.5 - 1.0e-12 && zmax <= 2.5 + 1.0e-12;
      const bool right_aperture = xmin >= 1.5 - 1.0e-12 && xmax <= 2.0 + 1.0e-12 &&
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
  strip_aperture_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] = "Warn";
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
  strip_aperture_meshes.push_back(std::make_unique<Mesh>(MakePairedApertureMesh()));
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
  no_strip_aperture_meshes.push_back(std::make_unique<Mesh>(MakePairedApertureMesh()));
  SpaceOperator no_strip_aperture_space(no_strip_aperture_iodata, no_strip_aperture_meshes);
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
  strict_mixed_signature_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] =
      "Error";
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
