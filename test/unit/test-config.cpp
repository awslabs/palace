// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"
#include "utils/jsonschema.hpp"

using json = nlohmann::json;
using namespace palace;

TEST_CASE("Config Boundary Ports", "[config][Serial]")
{
  SECTION("Basic passing config with bool excitation")
  {
    json boundaries = {{"LumpedPort",
                        {{{"Attributes", {5}},
                          {"Index", 1},
                          {"R", 50.0},
                          {"Direction", "+Y"},
                          {"Excitation", true}},
                         {{"Attributes", {16}},
                          {"Index", 2},
                          {"R", 50.0},
                          {"Direction", "-Y"},
                          {"Excitation", false}},
                         {{"Active", false},
                          {"Index", 3},
                          {"Elements",
                           {{{"Attributes", {27}}, {"Direction", "+X"}},
                            {{"Attributes", {38}}, {"Direction", "-X"}}}},
                          {"R", 50.0}}}},
                       {"WavePort",
                        {{{"Attributes", {6}}, {"Index", 4}, {"Excitation", true}},
                         {{"Attributes", {17}}, {"Index", 5}, {"Excitation", false}},
                         {{"Attributes", {18}}, {"Index", 6}, {"Active", false}}}},
                       {"PEC", {{"Attributes", {3}}}},
                       {"Absorbing", {{"Attributes", {4}}}}};
    config::BoundaryData boundary_ex_bool(boundaries);

    CHECK(boundary_ex_bool.lumpedport.at(1).active);
    CHECK(boundary_ex_bool.lumpedport.at(3).active == false);
    CHECK(boundary_ex_bool.lumpedport.at(1).excitation != 0);
    CHECK(boundary_ex_bool.lumpedport.at(3).excitation == 0);
    CHECK(boundary_ex_bool.waveport.at(5).excitation == 0);
    CHECK(boundary_ex_bool.waveport.at(6).excitation == 0);

    // Equivalent config with int excitation.
    json boundaries_int = {{"LumpedPort",
                            {{{"Attributes", {5}},
                              {"Index", 1},
                              {"R", 50.0},
                              {"Direction", "+Y"},
                              {"Excitation", 1}},
                             {{"Attributes", {16}},
                              {"Index", 2},
                              {"R", 50.0},
                              {"Direction", "-Y"},
                              {"Excitation", 0}},
                             {{"Active", false},
                              {"Index", 3},
                              {"Elements",
                               {{{"Attributes", {27}}, {"Direction", "+X"}},
                                {{"Attributes", {38}}, {"Direction", "-X"}}}},
                              {"R", 50.0}}}},
                           {"WavePort",
                            {{{"Attributes", {6}}, {"Index", 4}, {"Excitation", 1}},
                             {{"Attributes", {16}}, {"Index", 5}, {"Excitation", 0}},
                             {{"Attributes", {18}}, {"Index", 6}, {"Active", false}}}},
                           {"PEC", {{"Attributes", {3}}}},
                           {"Absorbing", {{"Attributes", {4}}}}};
    config::BoundaryData boundary_ex_int(boundaries_int);

    REQUIRE(boundary_ex_bool.lumpedport.size() == boundary_ex_int.lumpedport.size());
    auto it_int = boundary_ex_int.lumpedport.begin();
    auto it_bool = boundary_ex_bool.lumpedport.begin();
    for (; it_int != boundary_ex_int.lumpedport.end(); it_int++, it_bool++)
    {
      CHECK(it_bool->first == it_int->first);
      CHECK(it_bool->second.excitation == it_int->second.excitation);
    }
  }

  SECTION("Negative excitation throws")
  {
    json boundaries1 = {
        {"LumpedPort",
         {{{"Attributes", {5}},
           {"Index", 1},
           {"R", 50},
           {"Direction", "+Y"},
           {"Excitation", -1}}}},
        {"WavePort", {{{"Attributes", {6}}, {"Index", 2}, {"Excitation", 1}}}}};
    CHECK_THROWS(config::BoundaryData(boundaries1));

    json boundaries2 = {
        {"LumpedPort",
         {{{"Attributes", {5}},
           {"Index", 1},
           {"R", 50},
           {"Direction", "+Y"},
           {"Excitation", 1}}}},
        {"WavePort", {{{"Attributes", {6}}, {"Index", 2}, {"Excitation", -1}}}}};
    CHECK_THROWS(config::BoundaryData(boundaries2));
  }

  SECTION("Repeated index within LumpedPort throws")
  {
    json boundaries = {
        {"LumpedPort",
         {{{"Attributes", {5}},
           {"Index", 1},
           {"R", 50},
           {"Direction", "+Y"},
           {"Excitation", true}},
          {{"Attributes", {16}},
           {"Index", 1},
           {"R", 50},
           {"Direction", "-Y"},
           {"Excitation", false}}}},
        {"WavePort", {{{"Attributes", {6}}, {"Index", 2}, {"Excitation", true}}}}};
    CHECK_THROWS(config::BoundaryData(boundaries));
  }

  SECTION("Repeated index within WavePort throws")
  {
    json boundaries = {{"LumpedPort",
                        {{{"Attributes", {5}},
                          {"Index", 1},
                          {"R", 50},
                          {"Direction", "+Y"},
                          {"Excitation", true}}}},
                       {"WavePort",
                        {{{"Attributes", {6}}, {"Index", 2}, {"Excitation", true}},
                         {{"Attributes", {7}}, {"Index", 2}, {"Excitation", true}}}}};
    CHECK_THROWS(config::BoundaryData(boundaries));
  }

  SECTION("Repeated index cross-array fails validation")
  {
    json boundaries = {
        {"LumpedPort",
         {{{"Attributes", {5}},
           {"Index", 1},
           {"R", 50},
           {"Direction", "+Y"},
           {"Excitation", true}},
          {{"Attributes", {16}},
           {"Index", 2},
           {"R", 50},
           {"Direction", "-Y"},
           {"Excitation", false}}}},
        {"WavePort", {{{"Attributes", {6}}, {"Index", 1}, {"Excitation", true}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(config::Validate(boundary_data).has_value());
  }

  SECTION("Mislabeled excitation index fails validation")
  {
    json boundaries1 = {
        {"LumpedPort",
         {{{"Attributes", {5}},
           {"Index", 1},
           {"R", 50},
           {"Direction", "+Y"},
           {"Excitation", 2}}}},
        {"WavePort", {{{"Attributes", {6}}, {"Index", 2}, {"Excitation", 1}}}}};
    config::BoundaryData bd1(boundaries1);
    CHECK(config::Validate(bd1).has_value());

    json boundaries2 = {
        {"LumpedPort",
         {{{"Attributes", {5}},
           {"Index", 1},
           {"R", 50},
           {"Direction", "+Y"},
           {"Excitation", 2}}}},
        {"WavePort", {{{"Attributes", {6}}, {"Index", 2}, {"Excitation", 0}}}}};
    config::BoundaryData bd2(boundaries2);
    CHECK(config::Validate(bd2).has_value());
  }

  SECTION("Upgrade excitation index 1")
  {
    json boundaries = {
        {"LumpedPort",
         {{{"Attributes", {5}}, {"Index", 1}, {"R", 50}, {"Direction", "+Y"}},
          {{"Attributes", {16}},
           {"Index", 2},
           {"R", 50},
           {"Direction", "-Y"},
           {"Excitation", true}}}},
        {"WavePort",
         {{{"Attributes", {6}}, {"Index", 4}}, {{"Attributes", {17}}, {"Index", 5}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(boundary_data.lumpedport.at(1).excitation == 0);
    CHECK(boundary_data.lumpedport.at(2).excitation == 2);
    CHECK(boundary_data.waveport.at(4).excitation == 0);
    CHECK(boundary_data.waveport.at(5).excitation == 0);
  }

  SECTION("Upgrade excitation index 2")
  {
    json boundaries = {
        {"LumpedPort",
         {{{"Attributes", {5}}, {"Index", 1}, {"R", 50}, {"Direction", "+Y"}},
          {{"Attributes", {16}},
           {"Index", 2},
           {"R", 50},
           {"Direction", "-Y"},
           {"Excitation", 1}}}},
        {"WavePort",
         {{{"Attributes", {6}}, {"Index", 4}}, {{"Attributes", {17}}, {"Index", 5}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(boundary_data.lumpedport.at(1).excitation == 0);
    CHECK(boundary_data.lumpedport.at(2).excitation == 2);
    CHECK(boundary_data.waveport.at(4).excitation == 0);
    CHECK(boundary_data.waveport.at(5).excitation == 0);
  }

  SECTION("Upgrade excitation index 3")
  {
    json boundaries = {
        {"LumpedPort",
         {{{"Attributes", {5}}, {"Index", 1}, {"R", 50}, {"Direction", "+Y"}},
          {{"Attributes", {16}}, {"Index", 2}, {"R", 50}, {"Direction", "-Y"}}}},
        {"WavePort",
         {{{"Attributes", {6}}, {"Index", 4}, {"Excitation", true}},
          {{"Attributes", {17}}, {"Index", 5}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(boundary_data.lumpedport.at(1).excitation == 0);
    CHECK(boundary_data.lumpedport.at(2).excitation == 0);
    CHECK(boundary_data.waveport.at(4).excitation == 4);
    CHECK(boundary_data.waveport.at(5).excitation == 0);
  }

  SECTION("Upgrade excitation index 4")
  {
    json boundaries = {
        {"LumpedPort",
         {{{"Attributes", {5}},
           {"Index", 1},
           {"R", 50},
           {"Direction", "+Y"},
           {"Excitation", true}},
          {{"Attributes", {16}}, {"Index", 2}, {"R", 50}, {"Direction", "-Y"}}}},
        {"WavePort",
         {{{"Attributes", {6}}, {"Index", 4}, {"Excitation", true}},
          {{"Attributes", {17}}, {"Index", 5}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(boundary_data.lumpedport.at(1).excitation == 1);
    CHECK(boundary_data.lumpedport.at(2).excitation == 0);
    CHECK(boundary_data.waveport.at(4).excitation == 1);
    CHECK(boundary_data.waveport.at(5).excitation == 0);
  }

  SECTION("Shared excitation valid")
  {
    json boundaries = {{"LumpedPort",
                        {{{"Attributes", {5}},
                          {"Index", 1},
                          {"R", 50},
                          {"Direction", "+Y"},
                          {"Excitation", 1}},
                         {{"Attributes", {6}},
                          {"Index", 2},
                          {"R", 50},
                          {"Direction", "-Y"},
                          {"Excitation", 1}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(!config::Validate(boundary_data).has_value());
    CHECK(boundary_data.lumpedport.at(1).excitation == 1);
    CHECK(boundary_data.lumpedport.at(2).excitation == 1);
  }

  SECTION("Multi-excitation valid")
  {
    json boundaries = {{"LumpedPort",
                        {{{"Attributes", {5}},
                          {"Index", 1},
                          {"R", 50},
                          {"Direction", "+Y"},
                          {"Excitation", 1}},
                         {{"Attributes", {6}},
                          {"Index", 2},
                          {"R", 50},
                          {"Direction", "-Y"},
                          {"Excitation", 2}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(!config::Validate(boundary_data).has_value());
    CHECK(boundary_data.lumpedport.at(1).excitation == 1);
    CHECK(boundary_data.lumpedport.at(2).excitation == 2);
  }

  SECTION("Terminal duplicate index fails validation")
  {
    json boundaries = {
        {"LumpedPort",
         {{{"Attributes", {5}}, {"Index", 1}, {"R", 50}, {"Direction", "+Y"}}}},
        {"Terminal", {{{"Attributes", {6}}, {"Index", 1}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(config::Validate(boundary_data).has_value());
  }

  SECTION("SurfaceCurrent duplicate index fails validation")
  {
    json boundaries = {
        {"LumpedPort",
         {{{"Attributes", {5}}, {"Index", 1}, {"R", 50}, {"Direction", "+Y"}}}},
        {"SurfaceCurrent", {{{"Attributes", {6}}, {"Index", 1}, {"Direction", "+X"}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(config::Validate(boundary_data).has_value());
  }
}

TEST_CASE("Config Driven Solver", "[config][Serial]")
{
  using namespace Catch::Matchers;
  constexpr double delta_eps = 1.0e-9;

  SECTION("Base uniform sample")
  {
    json driven = {{"MinFreq", 0.1},
                   {"MaxFreq", 1.1},
                   {"FreqStep", 0.1},
                   {"SaveStep", 2},
                   {"Restart", 3}};
    config::DrivenSolverData driven_solver(driven);

    auto sample_f = std::vector{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1};
    auto save_indices = std::vector<size_t>{0, 2, 4, 6, 8, 10};

    for (size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
    CHECK(driven_solver.save_indices == save_indices);
    CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
  }

  SECTION("Uniform freq step with Samples")
  {
    json driven = {
        {"MinFreq", 0.1},
        {"MaxFreq", 1.1},
        {"FreqStep", 0.1},
        {"SaveStep", 2},
        {"Samples",
         {{{"MinFreq", 0.1}, {"MaxFreq", 1.1}, {"FreqStep", 0.1}, {"SaveStep", 2}}}}};
    config::DrivenSolverData driven_solver(driven);

    auto sample_f = std::vector{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1};
    auto save_indices = std::vector<size_t>{0, 2, 4, 6, 8, 10};

    for (size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
    CHECK(driven_solver.save_indices == save_indices);
    CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
  }

  SECTION("Uniform NSample")
  {
    json driven = {
        {"Samples",
         {{{"MinFreq", 0.0}, {"MaxFreq", 1.0}, {"NSample", 9}, {"SaveStep", 2}}}}};
    config::DrivenSolverData driven_solver(driven);

    auto sample_f = std::vector{0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
    auto save_indices = std::vector<size_t>{0, 2, 4, 6, 8};

    for (size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
    CHECK(driven_solver.save_indices == save_indices);
    CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
  }

  SECTION("Paired uniform sample")
  {
    json driven = {
        {"Samples",
         {{{"MinFreq", 0.0}, {"MaxFreq", 1.0}, {"NSample", 5}, {"SaveStep", 2}},
          {{"MinFreq", 0.0}, {"MaxFreq", 10.0}, {"NSample", 5}, {"SaveStep", 1}}}}};
    config::DrivenSolverData driven_solver(driven);

    auto sample_f = std::vector{0.0, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 7.5, 10.0};
    auto save_indices = std::vector<size_t>{0, 2, 4, 5, 6, 7, 8};

    for (size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
    CHECK(driven_solver.save_indices == save_indices);
    CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
  }

  SECTION("Uniform with point samples")
  {
    json driven = {{"Samples",
                    {{{"MinFreq", 0.0}, {"MaxFreq", 1.0}, {"NSample", 9}, {"SaveStep", 2}},
                     {{"Freq", {0.15, 0.35, 0.55}}, {"SaveStep", 1}, {"AddToPROM", true}}}},
                   {"Save", {0.0, 0.35, 1.0}},
                   {"AdaptiveTol", 1e-3}};
    config::DrivenSolverData driven_solver(driven);

    auto sample_f = std::vector{0.0, 0.125, 0.15,  0.25, 0.35,  0.375,
                                0.5, 0.55,  0.625, 0.75, 0.875, 1.0};
    auto save_indices = std::vector<size_t>{0, 2, 3, 4, 6, 7, 9, 11};
    auto prom_indices = std::vector<size_t>{0, sample_f.size() - 1, 2, 4, 7};

    for (size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
    CHECK(driven_solver.save_indices == save_indices);
    CHECK(driven_solver.prom_indices == prom_indices);
  }

  SECTION("Log with point samples")
  {
    json driven = {{"Samples",
                    {{{"Type", "Log"},
                      {"MinFreq", 0.1},
                      {"MaxFreq", 1.0},
                      {"NSample", 5},
                      {"SaveStep", 2}},
                     {{"Freq", {0.15, 0.35, 0.55}}, {"SaveStep", 1}, {"AddToPROM", true}}}},
                   {"AdaptiveTol", 1e-3}};
    config::DrivenSolverData driven_solver(driven);

    auto sample_f = std::vector{0.1,  0.15, 0.1778279410038923, 0.31622776601683794,
                                0.35, 0.55, 0.5623413251903491, 1.0};
    auto save_indices = std::vector<size_t>{0, 1, 3, 4, 5, 7};
    auto prom_indices = std::vector<size_t>{0, sample_f.size() - 1, 1, 4, 5};

    for (size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
    CHECK(driven_solver.save_indices == save_indices);
    CHECK(driven_solver.prom_indices == prom_indices);
  }

  SECTION("Invalid configs throw")
  {
    // Empty
    CHECK_THROWS(config::DrivenSolverData(json::object()));

    // Mismatch type - Point with range
    json mismatch1 = {{"Samples",
                       {{{"Type", "Point"},
                         {"MinFreq", 0.0},
                         {"MaxFreq", 1.0},
                         {"NSample", 11},
                         {"SaveStep", 2}}}}};
    CHECK_THROWS(config::DrivenSolverData(mismatch1));

    // Mismatch type - Log with FreqStep
    json mismatch2 = {{"Samples",
                       {{{"Type", "Log"},
                         {"MinFreq", 0.1},
                         {"MaxFreq", 1.0},
                         {"FreqStep", 11},
                         {"SaveStep", 2}}}}};
    CHECK_THROWS(config::DrivenSolverData(mismatch2));

    // Mismatch type - Linear with Freq array
    json mismatch3 = {{"Samples", {{{"Type", "Linear"}, {"Freq", {0.11, 0.22, 0.34}}}}}};
    CHECK_THROWS(config::DrivenSolverData(mismatch3));

    // Invalid log range (MinFreq = 0)
    json invalid_log = {{"Samples",
                         {{{"Type", "Log"},
                           {"MinFreq", 0.0},
                           {"MaxFreq", 1.0},
                           {"NSample", 11},
                           {"SaveStep", 2}}}}};
    CHECK_THROWS(config::DrivenSolverData(invalid_log));

    // Invalid save frequency
    json invalid_save = {
        {"Samples",
         {{{"MinFreq", 0.0}, {"MaxFreq", 1.0}, {"NSample", 11}, {"SaveStep", 2}},
          {{"Freq", {0.15, 0.35, 0.55}}, {"SaveStep", 1}, {"AddToPROM", true}}}},
        {"Save", {0.05}},
        {"AdaptiveTol", 1e-3}};
    CHECK_THROWS(config::DrivenSolverData(invalid_save));
  }
}

TEST_CASE("Config Linear Solver MaxIts", "[config][Serial]")
{
  SECTION("Linear solver MaxIts = 0 should fail schema validation")
  {
    json linear = {{"MaxIts", 0}};
    std::string err = ValidateConfig(linear, "Linear");
    CHECK(!err.empty());
    CHECK(err.find("MaxIts") != std::string::npos);
  }

  SECTION("Linear solver with valid MaxIts should pass schema validation")
  {
    json linear = {{"MaxIts", 1}};
    std::string err = ValidateConfig(linear, "Linear");
    CHECK(err.empty());
  }
}

TEST_CASE("FarField", "[config][Serial]")
{
  constexpr double delta_eps = 1.0e-6;

  SECTION("Default constructor")
  {
    config::FarFieldPostData data;

    CHECK(data.attributes.empty());
    CHECK(data.thetaphis.empty());
  }

  SECTION("Basic setup with attributes only")
  {
    json farfield = {{"Attributes", {1, 3, 5}}};
    config::FarFieldPostData data(farfield);

    CHECK(data.attributes == std::vector<int>{1, 3, 5});
    CHECK(data.thetaphis.empty());
  }

  SECTION("ThetaPhis conversion to radians")
  {
    json farfield = {{"Attributes", {1}}, {"ThetaPhis", {{0.0, 0.0}, {90.0, 180.0}}}};
    config::FarFieldPostData data(farfield);

    CHECK(data.thetaphis.size() == 2);
    CHECK(data.thetaphis[0].first == Catch::Approx(0.0).margin(delta_eps));
    CHECK(data.thetaphis[0].second == Catch::Approx(0.0).margin(delta_eps));
    CHECK(data.thetaphis[1].first == Catch::Approx(M_PI / 2).margin(delta_eps));
    CHECK(data.thetaphis[1].second == Catch::Approx(M_PI).margin(delta_eps));
  }

  SECTION("Duplicate removal")
  {
    json farfield = {{"Attributes", {1}},
                     {"ThetaPhis", {{0.0, 0.0}, {90.0, 180.0}, {0.0, 0.0}, {90.0, 180.0}}}};
    config::FarFieldPostData data(farfield);

    CHECK(data.thetaphis.size() == 2);
  }

  SECTION("Spherical coordinate duplicate removal")
  {
    // Test pole singularity: (0°, any φ) should be treated as same point.
    {
      json farfield = {{"Attributes", {1}},
                       {"ThetaPhis", {{0.0, 0.0}, {0.0, 90.0}, {0.0, 180.0}}}};
      config::FarFieldPostData data(farfield);
      CHECK(data.thetaphis.size() == 1);
    }

    // Test phi periodicity: φ and φ+360° are same point.
    {
      json farfield = {{"Attributes", {1}}, {"ThetaPhis", {{45.0, 30.0}, {45.0, 390.0}}}};
      config::FarFieldPostData data(farfield);
      CHECK(data.thetaphis.size() == 1);
    }

    // Test theta reflection: (θ, φ) ≡ (180°-θ, φ+180°).
    {
      json farfield = {{"Attributes", {1}}, {"ThetaPhis", {{60.0, 45.0}, {120.0, 225.0}}}};
      config::FarFieldPostData data(farfield);
      CHECK(data.thetaphis.size() == 1);
    }
  }

  SECTION("Combined NSample and ThetaPhis")
  {
    json farfield = {{"Attributes", {1}}, {"NSample", 10}, {"ThetaPhis", {{33.0, 22.0}}}};
    config::FarFieldPostData data(farfield);

    CHECK(data.thetaphis.size() == 11);
  }

  SECTION("NSample sphere sampling")
  {
    for (int nsample : {2, 6, 10, 15, 20, 25, 64800})
    {
      json farfield = {{"Attributes", {1}}, {"NSample", nsample}};
      config::FarFieldPostData data(farfield);

      // Exact point count.
      CHECK(data.thetaphis.size() == nsample);

      // Poles always included.
      bool has_north_pole = false, has_south_pole = false;
      for (const auto &point : data.thetaphis)
      {
        if (std::abs(point.first) < delta_eps)
          has_north_pole = true;
        if (std::abs(point.first - M_PI) < delta_eps)
          has_south_pole = true;
      }
      CHECK(has_north_pole);
      CHECK(has_south_pole);

      if (nsample == 2)
        continue;

      int n_theta = std::max(1, static_cast<int>(std::sqrt(nsample - 2)));

      // XZ plane inclusion (excluding poles).
      int xz_plane_count = 0;
      for (const auto &point : data.thetaphis)
      {
        if ((std::abs(point.second) < delta_eps ||
             std::abs(point.second - M_PI) < delta_eps) &&
            std::abs(point.first) > delta_eps && std::abs(point.first - M_PI) > delta_eps)
        {
          xz_plane_count++;
        }
      }
      // Each ring contributes two points.
      CHECK(xz_plane_count >= 2 * n_theta);

      // Equator inclusion.
      int equator_count = 0;
      for (const auto &point : data.thetaphis)
      {
        if (std::abs(point.first - M_PI / 2.0) < delta_eps)
        {
          equator_count++;
        }
      }
      CHECK(equator_count >= std::max(1, (nsample - 2) / (2 * n_theta)));
    }
  }

  SECTION("NSample edge cases")
  {
    // NSample = 0 produces no points.
    {
      json farfield = {{"Attributes", {1}}, {"NSample", 0}};
      config::FarFieldPostData data(farfield);
      CHECK(data.thetaphis.empty());
    }

    // NSample = 1 produces two points (the poles).
    {
      json farfield = {{"Attributes", {1}}, {"NSample", 1}};
      config::FarFieldPostData data(farfield);
      CHECK(data.thetaphis.size() == 2);
    }
  }
}

TEST_CASE("ParseStringAsDirection", "[config][Serial]")
{
  SECTION("Cartesian")
  {
    auto s = GENERATE("x", "+x", "-x", "y", "+y", "-y", "z", "+z", "-z");
    REQUIRE_NOTHROW(config::ParseStringAsDirection(s, true));

    auto [dir, cs] = config::ParseStringAsDirection(s, true);
    CHECK(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] == 1.0);
    CHECK(cs == CoordinateSystem::CARTESIAN);
    auto val = (s[0] == '-' ? -1 : 1);
    CHECK((dir[0] == val || dir[1] == val || dir[2] == val));
  }

  SECTION("Cylindrical")
  {
    auto s = GENERATE("r", "+r", "-r");
    REQUIRE_NOTHROW(config::ParseStringAsDirection(s, true));

    auto [dir, cs] = config::ParseStringAsDirection(s, true);
    CHECK(dir[0] == (s[0] == '-' ? -1 : 1));
    CHECK(dir[1] == 0);
    CHECK(dir[2] == 0);
    CHECK(cs == CoordinateSystem::CYLINDRICAL);
  }

  CHECK_NOTHROW(config::ParseStringAsDirection("", false));
}
