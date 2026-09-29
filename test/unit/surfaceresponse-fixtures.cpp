// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/communication.hpp"
#include "utils/edgedistance.hpp"
#include "utils/iodata.hpp"
#include "utils/metaledge.hpp"

namespace palace::test
{

using namespace Catch::Matchers;

SurfaceResponseFiles::SurfaceResponseFiles()
{
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
    // Coupon surface response files carry the energy within R of the coupon edges
    // (`Q_ij (J)` at `R (m)`, the SI matching radius) and the whole-box `Q_total_ij (J)`.
    // A spatial (3D box) model adds the within-R energy; the fixtures make both equal
    // unless a test says otherwise (`box_scale`). Every 3D spatial library below uses
    // MatchingRadius 0.2, and these IoData are not nondimensionalized (Units(1, 1): mesh
    // units are metres), so the SI radius is 0.2 m.
    constexpr double spatial_radius_m = 0.2;
    auto write_surface_matrix = [](const auto &path, const auto &matrix,
                                   double box_scale = 1.0, double radius_m = 0.2)
    {
      std::ofstream output(path);
      output << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
      for (int edge = 1; edge <= 2; edge++)
      {
        for (std::size_t i = 0; i < matrix.size(); i++)
        {
          for (std::size_t j = i; j < matrix.size(); j++)
          {
            output << "1," << edge << "," << radius_m << "," << i + 1 << "," << j + 1 << ","
                   << 0.5 * matrix[i][j] << "," << 0.5 * box_scale * matrix[i][j] << "\n";
          }
        }
      }
    };
    auto write_compact_surface_matrix = [](const auto &path, const auto &matrix)
    {
      std::ofstream output(path);
      output << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
      for (std::size_t i = 0; i < matrix.size(); i++)
      {
        for (std::size_t j = i; j < matrix.size(); j++)
        {
          output << "1,1," << spatial_radius_m << "," << i + 1 << "," << j + 1 << ","
                 << matrix[i][j] << "," << matrix[i][j] << "\n";
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
                          {"TraceLiftVersion", 2},
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
    // Pair models with explicit, non-default conductor references (the 2D pair patches must
    // carry them: a ResponsePatchData starts with the configuration default {{0, 0, 0}}).
    auto reference_strip_library_2d = exact_pair_library_2d;
    reference_strip_library_2d["Name"] = "unit-test-reference-strip-2d";
    reference_strip_library_2d["Models"][0]["Reference"] = {0.05, -0.1, 0.0};
    std::ofstream reference_strip_output_2d(reference_strip_library_2d_path);
    reference_strip_output_2d << reference_strip_library_2d.dump(2) << "\n";

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
      std::ofstream shared_trace(shared_boundary_trace_path);
      shared_trace << "x,y,z,V\n"
                   << "1.0,0.0,0.0,0.0\n"
                   << "1.0,0.5,0.0,0.5\n"
                   << "1.0,1.0,0.0,1.0\n";
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
    // The two-reference gap model of the 2D conductor-reference regression (the coupled
    // model's matrices and contour paths, a 2D separation and references).
    auto reference_gap_library_2d = different_pair_library_2d;
    reference_gap_library_2d["Name"] = "unit-test-reference-gap-2d";
    reference_gap_library_2d["TraceLiftVersion"] = coupled_library_3d["TraceLiftVersion"];
    auto reference_gap_model_2d = different_pair_library_2d["Models"][0];
    reference_gap_model_2d["Name"] = "different-gap-0.4-referenced";
    reference_gap_model_2d["ConductorReferences"] = {{-0.2, 0.0, 0.0}, {0.2, 0.0, 0.0}};
    for (const char *key : {"OpenContourPaths", "FabricatedMatrix", "ThinMatrix",
                            "FabricatedSurfaceMatrix", "ThinSurfaceMatrix"})
    {
      reference_gap_model_2d[key] = coupled_model[key];
    }
    reference_gap_library_2d["Models"] = {reference_gap_model_2d};
    std::ofstream reference_gap_output_2d(reference_gap_library_2d_path);
    reference_gap_output_2d << reference_gap_library_2d.dump(2) << "\n";
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
    auto write_corner_variant_library =
        [&](const auto &library_path, const std::string &tag, double within_scale,
            double box_scale, double radius_m, bool legacy_compact)
    {
      const auto fabricated_surface =
          temp.temp_dir / ("corner-fabricated-surface-" + tag + ".csv");
      const auto thin_surface = temp.temp_dir / ("corner-thin-surface-" + tag + ".csv");
      for (const auto &[path, matrix] : {std::pair{fabricated_surface, corner_fabricated},
                                         std::pair{thin_surface, corner_thin}})
      {
        if (!legacy_compact)
        {
          // Same row structure (two half-edges, default stream precision) as the base
          // corner file so that equal within-R values print identically.
          auto scaled = matrix;
          for (auto &row : scaled)
          {
            for (auto &value : row)
            {
              value *= within_scale;
            }
          }
          write_surface_matrix(path, scaled, box_scale / within_scale, radius_m);
          continue;
        }
        std::ofstream output(path);
        output << "interface,edge,basis_i,basis_j,Q_total_ij (J)\n";
        for (std::size_t i = 0; i < matrix.size(); i++)
        {
          for (std::size_t j = i; j < matrix.size(); j++)
          {
            output << "1,1," << i + 1 << "," << j + 1 << "," << matrix[i][j] << "\n";
          }
        }
      }
      auto library = convex_library_3d;
      library["Name"] = "unit-test-process-convex-3d-" + tag;
      library["Models"][1]["FabricatedSurfaceMatrix"] = fabricated_surface.string();
      library["Models"][1]["ThinSurfaceMatrix"] = thin_surface.string();
      std::ofstream output(library_path);
      output << library.dump(2) << "\n";
    };
    write_corner_variant_library(inflated_box_convex_library_3d_path, "inflated-box", 1.0,
                                 3.0, spatial_radius_m, false);
    write_corner_variant_library(scaled_convex_library_3d_path, "scaled", 3.0, 3.0,
                                 spatial_radius_m, false);
    write_corner_variant_library(legacy_compact_convex_library_3d_path, "legacy-compact",
                                 1.0, 1.0, spatial_radius_m, true);
    write_corner_variant_library(other_radius_convex_library_3d_path, "other-radius", 1.0,
                                 1.0, 2.0 * spatial_radius_m, false);

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
    spatial_cluster_model["SupportPoints"] = {{-0.25, -0.25, -0.05}, {-0.25, -0.25, 0.05},
                                              {-0.25, 0.25, -0.05},  {-0.25, 0.25, 0.05},
                                              {0.25, -0.25, -0.05},  {0.25, -0.25, 0.05},
                                              {0.25, 0.25, -0.05},   {0.25, 0.25, 0.05}};
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

    // The same cluster with two cap-interior hats: basis points 5 and 6 lie on no contour
    // (the open paths partition the four ring knots only) and are declared by
    // InteriorTraceCount 2 through an explicit TraceMesh that represents every coefficient.
    // Matrices 7 x 7 (six points + the conductor state): the ring / conductor entries of
    // the coupled fixture, the hat rows zero (the response of the ring-only model) or a
    // diagonal hat energy equal in the fabricated and thin coupon (the surface energy
    // grows with the hats, the domain defect does not).
    {
      std::ofstream points(cap_hat_points_path);
      points << "x,y,z\n"
             << "0.025,-0.025,0.02\n"
             << "0.050,-0.050,0.02\n"
             << "0.075,-0.075,-0.02\n"
             << "0.100,-0.100,-0.02\n"
             << "0.060,-0.030,0.02\n"
             << "0.030,-0.060,-0.02\n";
      std::ofstream vertices(cap_hat_trace_vertices_path);
      vertices << "# index,x,y,z,basis,conductor\n"
               << "1,0.025,-0.025,0.02,1,0\n"
               << "2,0.050,-0.050,0.02,2,0\n"
               << "3,0.075,-0.075,-0.02,3,0\n"
               << "4,0.100,-0.100,-0.02,4,0\n"
               << "5,0.060,-0.030,0.02,5,0\n"
               << "6,0.030,-0.060,-0.02,6,0\n";
      std::ofstream triangles(cap_hat_trace_triangles_path);
      triangles << "# index,v1,v2,v3\n"
                << "1,1,2,5\n"
                << "2,2,3,5\n"
                << "3,3,4,6\n"
                << "4,2,3,6\n";
    }
    auto write_cap_hat_matrices = [&](double hat_energy, const auto &domain_fabricated_path,
                                      const auto &domain_thin_path,
                                      const auto &surface_fabricated_path,
                                      const auto &surface_thin_path)
    {
      std::array<std::array<double, 7>, 7> cap_hat_fabricated{};
      std::array<std::array<double, 7>, 7> cap_hat_thin{};
      for (std::size_t i = 0; i < 4; i++)
      {
        cap_hat_fabricated[i][i] = coupled_fabricated[i][i];
        cap_hat_thin[i][i] = coupled_thin[i][i];
      }
      cap_hat_fabricated[4][4] = cap_hat_fabricated[5][5] = hat_energy;
      cap_hat_thin[4][4] = cap_hat_thin[5][5] = hat_energy;
      cap_hat_fabricated[6][6] = coupled_fabricated[4][4];
      cap_hat_thin[6][6] = coupled_thin[4][4];
      write_domain_matrix(domain_fabricated_path, cap_hat_fabricated);
      write_domain_matrix(domain_thin_path, cap_hat_thin);
      write_surface_matrix(surface_fabricated_path, cap_hat_fabricated);
      write_surface_matrix(surface_thin_path, cap_hat_thin);
    };
    write_cap_hat_matrices(0.0, cap_hat_fabricated_path, cap_hat_thin_path,
                           cap_hat_fabricated_surface_path, cap_hat_thin_surface_path);
    write_cap_hat_matrices(1.0e-12, cap_hat_loaded_fabricated_path,
                           cap_hat_loaded_thin_path, cap_hat_loaded_fabricated_surface_path,
                           cap_hat_loaded_thin_surface_path);
    auto cap_hat_library_3d = spatial_cluster_library_3d;
    cap_hat_library_3d["Name"] = "unit-test-process-spatial-cluster-cap-hats-3d";
    auto &cap_hat_model = cap_hat_library_3d["Models"].back();
    cap_hat_model["Name"] = "offset-corner-pair-cap-hats";
    cap_hat_model["BasisPoints"] = cap_hat_points_path.string();
    cap_hat_model["TraceMesh"] = {{"Vertices", cap_hat_trace_vertices_path.string()},
                                  {"Triangles", cap_hat_trace_triangles_path.string()}};
    cap_hat_model["InteriorTraceCount"] = 2;
    cap_hat_model["FabricatedMatrix"] = cap_hat_fabricated_path.string();
    cap_hat_model["ThinMatrix"] = cap_hat_thin_path.string();
    cap_hat_model["FabricatedSurfaceMatrix"] = cap_hat_fabricated_surface_path.string();
    cap_hat_model["ThinSurfaceMatrix"] = cap_hat_thin_surface_path.string();
    std::ofstream cap_hat_output_3d(cap_hat_spatial_cluster_library_3d_path);
    cap_hat_output_3d << cap_hat_library_3d.dump(2) << "\n";

    auto cap_hat_loaded_library_3d = cap_hat_library_3d;
    cap_hat_loaded_library_3d["Name"] =
        "unit-test-process-spatial-cluster-cap-hats-loaded-3d";
    auto &cap_hat_loaded_model = cap_hat_loaded_library_3d["Models"].back();
    cap_hat_loaded_model["FabricatedMatrix"] = cap_hat_loaded_fabricated_path.string();
    cap_hat_loaded_model["ThinMatrix"] = cap_hat_loaded_thin_path.string();
    cap_hat_loaded_model["FabricatedSurfaceMatrix"] =
        cap_hat_loaded_fabricated_surface_path.string();
    cap_hat_loaded_model["ThinSurfaceMatrix"] = cap_hat_loaded_thin_surface_path.string();
    std::ofstream cap_hat_loaded_output_3d(cap_hat_loaded_spatial_cluster_library_3d_path);
    cap_hat_loaded_output_3d << cap_hat_loaded_library_3d.dump(2) << "\n";

    // Refused: open paths summing to BasisPoints (the hats on a contour), the key without
    // a TraceMesh, the key on a translational model.
    auto cap_hat_full_partition_library_3d = cap_hat_library_3d;
    cap_hat_full_partition_library_3d["Name"] =
        "unit-test-process-spatial-cluster-cap-hats-full-partition-3d";
    cap_hat_full_partition_library_3d["Models"].back()["OpenContourPaths"] = {
        {{"Indices", {1, 2}}, {"StartConductor", 1}, {"EndConductor", 2}},
        {{"Indices", {4, 3, 5, 6}}, {"StartConductor", 1}, {"EndConductor", 2}}};
    std::ofstream cap_hat_full_partition_output_3d(cap_hat_full_partition_library_3d_path);
    cap_hat_full_partition_output_3d << cap_hat_full_partition_library_3d.dump(2) << "\n";
    auto cap_hat_without_trace_mesh_library_3d = cap_hat_library_3d;
    cap_hat_without_trace_mesh_library_3d["Name"] =
        "unit-test-process-spatial-cluster-cap-hats-no-trace-mesh-3d";
    cap_hat_without_trace_mesh_library_3d["Models"].back().erase("TraceMesh");
    std::ofstream cap_hat_without_trace_mesh_output_3d(
        cap_hat_without_trace_mesh_library_3d_path);
    cap_hat_without_trace_mesh_output_3d << cap_hat_without_trace_mesh_library_3d.dump(2)
                                         << "\n";
    auto cap_hat_translational_library_3d = library_3d;
    cap_hat_translational_library_3d["Name"] =
        "unit-test-process-cap-hats-translational-3d";
    cap_hat_translational_library_3d["Models"][0]["InteriorTraceCount"] = 1;
    std::ofstream cap_hat_translational_output_3d(cap_hat_translational_library_3d_path);
    cap_hat_translational_output_3d << cap_hat_translational_library_3d.dump(2) << "\n";

    auto write_cross_layer_surface_matrix = [&](const auto &path, const auto &matrix)
    {
      std::ofstream output(path);
      output << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
      for (int interface = 1; interface <= 2; interface++)
      {
        for (std::size_t i = 0; i < matrix.size(); i++)
        {
          for (std::size_t j = i; j < matrix.size(); j++)
          {
            output << interface << ",1," << spatial_radius_m << "," << i + 1 << "," << j + 1
                   << "," << 0.5 * matrix[i][j] << "," << 0.5 * matrix[i][j] << "\n";
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
}

json SurfaceResponseFiles::AutomaticConfig2D() const
{
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
  return automatic_config;
}

json SurfaceResponseFiles::BoundaryModeConfig2D() const
{
  auto boundary_mode_config = AutomaticConfig2D();
  boundary_mode_config["Problem"]["Type"] = "BoundaryMode";
  boundary_mode_config["Boundaries"].erase("Ground");
  boundary_mode_config["Boundaries"].erase("Terminal");
  boundary_mode_config["Boundaries"]["PEC"] = {{"Attributes", {9, 10}}};
  boundary_mode_config["Solver"] = {
      {"Order", 1},
      {"BoundaryMode", {{"Freq", 5.0}}},
      {"SurfaceResponseCorrection",
       {{"Library", library_path.string()}, {"UnmatchedPolicy", "Error"}}}};
  return boundary_mode_config;
}

json SurfaceResponseFiles::ExactPairConfig2D() const
{
  auto exact_pair_config_2d = AutomaticConfig2D();
  exact_pair_config_2d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      exact_pair_library_2d_path.string();
  exact_pair_config_2d["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"] = {
      0.3};
  return exact_pair_config_2d;
}

json SurfaceResponseFiles::Cpw3dLegacyConfig() const
{
  const auto config_3d_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                              "examples/cpw3d_surface/cpw3d_surface_validation_thin.json";
  std::ifstream config_3d_input(config_3d_path);
  json config_3d = json::parse(config_3d_input);
  config_3d["Problem"]["Output"] = temp.temp_dir.string();
  config_3d["Model"]["Mesh"] =
      (fs::path(PALACE_TEST_DATA_DIR) / "mesh/cpw3d-surface-nc.msh").string();
  // The three-dimensional sections below exercise the legacy per-interface-group
  // classification (parametric matching with tolerances, interpolation brackets, plan-view
  // masks), kept behind PatchConstruction = "Legacy"; the Features-driven construction
  // (the default) is tested in its own section with a signature-keyed library.
  config_3d["Solver"]["Electrostatic"]["ResponseCorrection"] = {
      {"Library", library_3d_path.string()},
      {"TargetInterfaces", {1, 2, 3}},
      {"UnmatchedPolicy", "Error"},
      {"PatchConstruction", "Legacy"}};
  return config_3d;
}

json SurfaceResponseFiles::IslandConfig() const
{
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
            {"UnmatchedPolicy", "Error"},
            {"PatchConstruction", "Legacy"}}}}}}}};
  return island_config;
}

json SurfaceResponseFiles::ConvexIslandConfig() const
{
  auto convex_island_config = IslandConfig();
  convex_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      convex_library_3d_path.string();
  return convex_island_config;
}

json SurfaceResponseFiles::ConvexMaxwellIslandConfig() const
{
  auto convex_maxwell_island_config = ConvexIslandConfig();
  convex_maxwell_island_config["Problem"]["Type"] = "Eigenmode";
  convex_maxwell_island_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4,
                                                                        5, 6, 9};
  convex_maxwell_island_config["Boundaries"].erase("Terminal");
  convex_maxwell_island_config["Solver"] = {{"Order", 1},
                                            {"Eigenmode", {{"Target", 1.0}}},
                                            {"SurfaceResponseCorrection",
                                             {{"Library", convex_library_3d_path.string()},
                                              {"TargetInterfaces", {4}},
                                              {"UnmatchedPolicy", "Error"},
                                              {"PatchConstruction", "Legacy"}}}};
  return convex_maxwell_island_config;
}

json SurfaceResponseFiles::RoundedIslandConfig() const
{
  auto rounded_island_config = IslandConfig();
  rounded_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      rounded_library_3d_path.string();
  return rounded_island_config;
}

json SurfaceResponseFiles::HighOrderSpatialConfig() const
{
  auto high_order_spatial_config = IslandConfig();
  high_order_spatial_config["Boundaries"]["Terminal"] = {
      {{"Index", 1}, {"Attributes", {9}}}, {{"Index", 2}, {"Attributes", {10}}}};
  high_order_spatial_config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] =
      {9, 10};
  high_order_spatial_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      spatial_cluster_library_3d_path.string();
  high_order_spatial_config["Solver"]["Electrostatic"]["ResponseCorrection"]
                           ["UnmatchedPolicy"] = "Warn";
  return high_order_spatial_config;
}

void SurfaceResponseFiles::ShareCrackedAttributes(IoData &iodata, const IoData &reader)
{
  iodata.boundaries.cracked_attributes.insert(reader.boundaries.cracked_attributes.begin(),
                                              reader.boundaries.cracked_attributes.end());
}

mfem::Mesh SurfaceResponseFiles::MakeAutomatic2DMesh()
{
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
  return automatic_serial;
}

std::unique_ptr<mfem::ParMesh>
SurfaceResponseFiles::MakeIslandMesh(bool rounded, bool tetrahedral, bool aperture,
                                     bool neighboring_island, bool second_layer,
                                     bool high_order_rounded)
{
  const double in_plane_extent = aperture || neighboring_island ? 2.0 : 1.0;
  const double center = 0.5 * in_plane_extent;
  const double half_width = 0.25;
  const int in_plane_elements =
      (rounded && !high_order_rounded ? 16 : 8) * (aperture || neighboring_island ? 2 : 1);
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
}

std::unique_ptr<mfem::ParMesh> SurfaceResponseFiles::MakeTouchingIslandMesh()
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
}

std::unique_ptr<mfem::ParMesh> SurfaceResponseFiles::MakeOffsetCornerPairMesh()
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
}

InterfaceEdgeSummary SurfaceResponseFiles::SummarizeInterfaceEdges(
    const mfem::ParMesh &mesh, const config::BoundaryData &boundaries, int interface)
{
  const auto geometry =
      ExtractMetalEdgeGeometry(mesh, boundaries, JointNoiseExtractionFor(boundaries));
  InterfaceEdgeSummary summary;
  summary.segments =
      GetInterfaceMetalEdgeSegmentIndices(geometry, interface, InterfaceDielectric::SA);
  std::set<std::size_t> vertices;
  for (const std::size_t segment_index : summary.segments)
  {
    const auto &segment = geometry.segments[segment_index];
    vertices.insert(segment.vertices.begin(), segment.vertices.end());
    const auto &p0 = geometry.vertices[segment.vertices[0]].coordinate;
    const auto &p1 = geometry.vertices[segment.vertices[1]].coordinate;
    double length_squared = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length_squared += (p1[d] - p0[d]) * (p1[d] - p0[d]);
    }
    summary.length += std::sqrt(length_squared);
  }
  for (const std::size_t vertex : vertices)
  {
    const auto type = geometry.vertices[vertex].physical_type;
    summary.corners += type == MetalEdgeVertexType::CORNER ? 1 : 0;
    summary.junctions += type == MetalEdgeVertexType::JUNCTION ? 1 : 0;
  }
  return summary;
}

mfem::VectorConstantCoefficient SurfaceResponseFiles::ConstantFieldCoefficient()
{
  mfem::Vector constant_field(3);
  constant_field[0] = 0.7;
  constant_field[1] = -0.4;
  constant_field[2] = 0.2;
  return mfem::VectorConstantCoefficient(constant_field);
}

// The library builder's model Edges of a SpatialEdgeCluster Signature (canonical frame,
// library units): every straight portion one edge (Point = P0 x R, Interval along
// gap x normal), every arc portion chorded into n = max(ceil(sweep / 5 deg),
// ceil(arc length / 0.25 R), 1) chords with the radial gap at each chord's middle
// (signature_library.cluster_plan_view_edges). Conductor labels 1, 2, ... in order of
// first occurrence (the canonical serialisation relabels the same way); InterfaceSlot k =
// the k-th distinct portion interface set in sorted order, the slots SignatureInterfaces
// maps (every type of slot k to Coupon 1).
std::vector<std::string> SurfaceResponseFiles::SignatureInterfaceSets(const json &signature)
{
  std::set<std::string> sets;
  for (const auto &portion : signature["Portions"])
  {
    sets.insert(portion.value("Interfaces", json::array()).dump());
  }
  return std::vector<std::string>(sets.begin(), sets.end());
}

json SurfaceResponseFiles::SignatureInterfaces(const json &signature)
{
  json interfaces = json::array();
  const auto sets = SignatureInterfaceSets(signature);
  for (std::size_t slot = 0; slot < sets.size(); slot++)
  {
    for (const auto &type : json::parse(sets[slot]))
    {
      interfaces.push_back({{"Slot", slot}, {"Type", type}, {"Coupon", 1}});
    }
  }
  return interfaces;
}

json SurfaceResponseFiles::ChordedSignatureEdges(const json &signature, double radius)
{
  json edges = json::array();
  std::map<int, int> conductor_labels;
  const auto interface_sets = SignatureInterfaceSets(signature);
  int interface_slot = 0;
  auto AddEdge = [&](std::array<double, 2> a, std::array<double, 2> b,
                     std::array<double, 2> gap, int conductor)
  {
    const double gap_norm = std::hypot(gap[0], gap[1]);
    gap = {gap[0] / gap_norm, gap[1] / gap_norm};
    const double dx = (b[0] - a[0]) * radius, dy = (b[1] - a[1]) * radius;
    const double length = std::hypot(dx, dy);
    const bool forward = dx * gap[1] - dy * gap[0] > 0.0;  // along gap x (+z)
    edges.push_back({{"Point", {a[0] * radius, a[1] * radius, 0.0}},
                     {"GapDirection", {gap[0], gap[1], 0.0}},
                     {"ProcessNormal", {0.0, 0.0, 1.0}},
                     {"Interval", forward ? json{0.0, length} : json{-length, 0.0}},
                     {"Conductor", conductor},
                     {"InterfaceSlot", interface_slot},
                     {"BoundaryCondition", "PEC"}});
  };
  for (const auto &portion : signature["Portions"])
  {
    const int raw_conductor = portion["Conductor"].get<int>();
    if (!conductor_labels.count(raw_conductor))
    {
      conductor_labels[raw_conductor] = static_cast<int>(conductor_labels.size()) + 1;
    }
    const int conductor = conductor_labels.at(raw_conductor);
    interface_slot =
        static_cast<int>(std::find(interface_sets.begin(), interface_sets.end(),
                                   portion.value("Interfaces", json::array()).dump()) -
                         interface_sets.begin());
    const auto P = portion["P"].get<std::array<double, 4>>();
    const std::array<double, 2> a = {P[0], P[1]}, b = {P[2], P[3]};
    if (!portion.contains("Arc"))
    {
      AddEdge(a, b, portion["Gap"].get<std::array<double, 2>>(), conductor);
      continue;
    }
    const auto arc = portion["Arc"].get<std::array<double, 4>>();
    const std::array<double, 2> c = {arc[0], arc[1]}, m = {arc[2], arc[3]};
    const double r = std::hypot(a[0] - c[0], a[1] - c[1]);
    auto Angle = [&](const std::array<double, 2> &p)
    { return std::atan2(p[1] - c[1], p[0] - c[0]); };
    const double two_pi = 2.0 * std::acos(-1.0);
    const double ta = Angle(a);
    const double ccw = std::fmod(Angle(b) - ta + two_pi, two_pi);
    const double sweep =
        std::fmod(Angle(m) - ta + two_pi, two_pi) <= ccw + 1.0e-12 ? ccw : ccw - two_pi;
    const int chords = std::max(
        {static_cast<int>(std::ceil(std::abs(sweep) / (5.0 * two_pi / 360.0) - 1.0e-9)),
         static_cast<int>(std::ceil(r * std::abs(sweep) / 0.25 - 1.0e-9)), 1});
    const double gap_sign = portion["GapRadial"].get<int>();
    for (int k = 0; k < chords; k++)
    {
      const double t0 = ta + sweep * k / chords, t1 = ta + sweep * (k + 1) / chords;
      const double tm = 0.5 * (t0 + t1);
      AddEdge({c[0] + r * std::cos(t0), c[1] + r * std::sin(t0)},
              {c[0] + r * std::cos(t1), c[1] + r * std::sin(t1)},
              {gap_sign * std::cos(tm), gap_sign * std::sin(tm)}, conductor);
    }
  }
  return edges;
}

// Rows of the patch dry run (surface-response-patches.csv).
std::vector<std::vector<std::string>>
SurfaceResponseFiles::ReadPatchRows(const fs::path &path)
{
  std::vector<std::vector<std::string>> rows;
  std::ifstream input(path);
  REQUIRE(input);
  std::string line;
  std::getline(input, line);  // header
  while (std::getline(input, line))
  {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ','))
    {
      fields.push_back(field);
    }
    rows.push_back(std::move(fields));
  }
  return rows;
}

// Every endpoint of the model's Edges (library units = mesh units here), mapped through
// the feature's dry-run patch frame (Origin 13-15, AxisU / V / W 16-24), lies on the
// feature's claimed portions in the mesh: a straight portion endpoint, or a fitted arc
// (radially on the circle, within the angular range of the claimed chords).
std::size_t SurfaceResponseFiles::CheckPlacedModelEdges(
    const json &feature, const json &identification,
    const std::vector<std::vector<std::string>> &patch_rows, const json &model_edges,
    double tolerance)
{
  const auto &segments = identification["Segments"];
  const auto &arcs = identification["Arcs"];
  std::vector<std::array<double, 3>> straight_endpoints;
  struct ClaimedArc
  {
    std::array<double, 3> center;
    double radius;
    std::array<double, 3> u, v;  // in-plane basis: u toward the first claimed point
    double angle_min = 0.0, angle_max = 0.0;
  };
  std::map<int, ClaimedArc> claimed_arcs;
  const auto n = feature["Frame"]["Axes"][2].get<std::array<double, 3>>();
  auto ArcAngle = [&](const ClaimedArc &arc, const std::array<double, 3> &q)
  {
    double x = 0.0, y = 0.0;
    for (int d = 0; d < 3; d++)
    {
      x += (q[d] - arc.center[d]) * arc.u[d];
      y += (q[d] - arc.center[d]) * arc.v[d];
    }
    return std::atan2(y, x);
  };
  for (const auto &portion : feature["Portions"])
  {
    const auto &segment = segments[portion[0].get<int>()];
    const auto key = segment["Key"].get<std::array<std::array<double, 3>, 2>>();
    const double length = segment["Length"].get<double>();
    for (int end = 1; end <= 2; end++)
    {
      const double s = portion[end].get<double>() / length;
      const std::array<double, 3> point = {key[0][0] + s * (key[1][0] - key[0][0]),
                                           key[0][1] + s * (key[1][1] - key[0][1]),
                                           key[0][2] + s * (key[1][2] - key[0][2])};
      if (!segment.contains("Arc"))
      {
        straight_endpoints.push_back(point);
        continue;
      }
      const int arc_index = segment["Arc"].get<int>();
      auto it = claimed_arcs.find(arc_index);
      if (it == claimed_arcs.end())
      {
        ClaimedArc arc{arcs[arc_index]["Center"].get<std::array<double, 3>>(),
                       arcs[arc_index]["Radius"].get<double>()};
        std::array<double, 3> u{};
        double dot = 0.0;
        for (int d = 0; d < 3; d++)
        {
          dot += (point[d] - arc.center[d]) * n[d];
        }
        for (int d = 0; d < 3; d++)
        {
          u[d] = point[d] - arc.center[d] - dot * n[d];
        }
        const double norm = std::hypot(u[0], u[1], u[2]);
        arc.u = {u[0] / norm, u[1] / norm, u[2] / norm};
        arc.v = {n[1] * arc.u[2] - n[2] * arc.u[1], n[2] * arc.u[0] - n[0] * arc.u[2],
                 n[0] * arc.u[1] - n[1] * arc.u[0]};
        it = claimed_arcs.emplace(arc_index, arc).first;
      }
      const double angle = ArcAngle(it->second, point);
      it->second.angle_min = std::min(it->second.angle_min, angle);
      it->second.angle_max = std::max(it->second.angle_max, angle);
    }
  }
  const auto patch =
      std::find_if(patch_rows.begin(), patch_rows.end(), [&](const auto &row)
                   { return std::stoi(row[1]) == feature["Id"].get<int>(); });
  REQUIRE(patch != patch_rows.end());
  CHECK((*patch)[2] == "spatial edge cluster");  // the model topology label
  std::array<double, 3> origin{};
  std::array<std::array<double, 3>, 3> axes{};
  for (int d = 0; d < 3; d++)
  {
    origin[d] = std::stod((*patch)[13 + d]);
    for (int a = 0; a < 3; a++)
    {
      axes[a][d] = std::stod((*patch)[16 + 3 * a + d]);
    }
  }
  for (const auto &edge : model_edges)
  {
    const auto point = edge["Point"].get<std::array<double, 3>>();
    const auto gap = edge["GapDirection"].get<std::array<double, 3>>();
    const auto interval = edge["Interval"].get<std::array<double, 2>>();
    const std::array<double, 3> tangent = {gap[1], -gap[0], 0.0};  // gap x (+z)
    for (const double s : interval)
    {
      std::array<double, 3> mapped = origin;
      for (int d = 0; d < 3; d++)
      {
        const double local = point[d] + s * tangent[d];
        for (int k = 0; k < 3; k++)
        {
          mapped[k] += local * axes[d][k];
        }
      }
      double nearest = std::numeric_limits<double>::infinity();
      for (const auto &endpoint : straight_endpoints)
      {
        nearest =
            std::min(nearest, std::hypot(mapped[0] - endpoint[0], mapped[1] - endpoint[1],
                                         mapped[2] - endpoint[2]));
      }
      for (const auto &[index, arc] : claimed_arcs)
      {
        double out_of_plane = 0.0, in_plane2 = 0.0;
        for (int d = 0; d < 3; d++)
        {
          out_of_plane += (mapped[d] - arc.center[d]) * n[d];
        }
        for (int d = 0; d < 3; d++)
        {
          const double r = mapped[d] - arc.center[d] - out_of_plane * n[d];
          in_plane2 += r * r;
        }
        const double radial = std::abs(std::sqrt(in_plane2) - arc.radius);
        const double angle = ArcAngle(arc, mapped);
        const double angular_excess =
            std::max({0.0, arc.angle_min - angle, angle - arc.angle_max});
        nearest = std::min(nearest,
                           std::hypot(radial, angular_excess * arc.radius, out_of_plane));
      }
      INFO("feature " << feature["Id"] << " model edge point (" << point[0] << ", "
                      << point[1] << ") s " << s << " mapped to (" << mapped[0] << ", "
                      << mapped[1] << ", " << mapped[2] << ")");
      CHECK_THAT(nearest, WithinAbs(0.0, tolerance));
    }
  }
  return claimed_arcs.size();
}

GeometryCacheEnvGuard::GeometryCacheEnvGuard(const std::string &cache_path, bool write)
{
  setenv("PALACE_RESPONSE_GEOMETRY_CACHE", cache_path.c_str(), 1);
  if (write)
  {
    setenv("PALACE_RESPONSE_GEOMETRY_CACHE_WRITE", "1", 1);
  }
}

GeometryCacheEnvGuard::~GeometryCacheEnvGuard()
{
  unsetenv("PALACE_RESPONSE_GEOMETRY_CACHE_WRITE");
  unsetenv("PALACE_RESPONSE_GEOMETRY_CACHE");
}

void GeometryCacheEnvGuard::DisableWrite()
{
  unsetenv("PALACE_RESPONSE_GEOMETRY_CACHE_WRITE");
}

}  // namespace palace::test
