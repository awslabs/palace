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
                             {{"Attributes", {17}}, {"Index", 5}, {"Excitation", 0}},
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

TEST_CASE("ConcretizeDefaults", "[config][Serial]")
{
  SECTION("Electrostatic resolves linear solver sentinels")
  {
    json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", json::object()}};

    IoData iodata(config, false);

    // CheckConfiguration should have resolved the sentinels.
    CHECK(iodata.solver.linear.type == LinearSolver::BOOMER_AMG);
    CHECK(iodata.solver.linear.krylov_solver == KrylovSolver::CG);
    CHECK(iodata.solver.linear.initial_guess == 1);
    CHECK(iodata.solver.linear.pc_mat_shifted == 0);
    CHECK(iodata.solver.linear.mg_smooth_aux == 0);
    CHECK(iodata.solver.linear.ams_singular_op == 0);
    CHECK(iodata.solver.linear.amg_agg_coarsen == 1);
    CHECK(iodata.solver.linear.ams_max_it == 1);
    CHECK(iodata.solver.linear.mg_cycle_it == 1);
    CHECK(iodata.solver.linear.pc_mat_sym == MatrixSymmetry::SPD);

    // ConcretizeDefaults should write resolved values back to JSON.
    config = IoData::ConcretizeDefaults(iodata, config);

    auto &j_linear = config["Solver"]["Linear"];
    CHECK(j_linear["Type"].get<std::string>() == "BoomerAMG");
    CHECK(j_linear["KSPType"].get<std::string>() == "CG");
    CHECK(j_linear["InitialGuess"].get<int>() == 1);
    CHECK(j_linear["PCMatShifted"].get<int>() == 0);
    CHECK(j_linear["MGAuxiliarySmoother"].get<int>() == 0);
    CHECK(j_linear["AMSSingularOperator"].get<int>() == 0);
    CHECK(j_linear["AMGAggressiveCoarsening"].get<int>() == 1);
    CHECK(j_linear["AMSMaxIts"].get<int>() == 1);
    CHECK(j_linear["MGCycleIts"].get<int>() == 1);
    CHECK(j_linear["PCMatSymmetry"].get<std::string>() == "SPD");
  }

  SECTION("Magnetostatic resolves to AMS with singular operator")
  {
    json config = {{"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", json::object()}};

    IoData iodata(config, false);
    CHECK(iodata.solver.linear.type == LinearSolver::AMS);
    CHECK(iodata.solver.linear.krylov_solver == KrylovSolver::CG);
    CHECK(iodata.solver.linear.ams_singular_op == 1);
    CHECK(iodata.solver.linear.ams_max_it == 1);
    CHECK(iodata.solver.linear.mg_cycle_it == 1);
    CHECK(iodata.solver.linear.pc_mat_sym == MatrixSymmetry::SPD);

    config = IoData::ConcretizeDefaults(iodata, config);
    CHECK(config["Solver"]["Linear"]["Type"].get<std::string>() == "AMS");
    CHECK(config["Solver"]["Linear"]["AMSSingularOperator"].get<int>() == 1);
    CHECK(config["Solver"]["Linear"]["AMSMaxIts"].get<int>() == 1);
    CHECK(config["Solver"]["Linear"]["MGCycleIts"].get<int>() == 1);
    CHECK(config["Solver"]["Linear"]["PCMatSymmetry"].get<std::string>() == "SPD");
  }

  SECTION("User-specified values survive concretization")
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries", json::object()},
        {"Solver", {{"Order", 3}, {"Linear", {{"Type", "BoomerAMG"}, {"KSPType", "CG"}}}}}};

    IoData iodata(config, false);
    config = IoData::ConcretizeDefaults(iodata, config);

    // User-specified values should be preserved.
    CHECK(config["Solver"]["Order"].get<int>() == 3);
    CHECK(config["Solver"]["Linear"]["Type"].get<std::string>() == "BoomerAMG");
    CHECK(config["Solver"]["Linear"]["KSPType"].get<std::string>() == "CG");
  }

  SECTION("EigenSolverBackend DEFAULT resolves to concrete backend")
  {
    json config = {{"Problem", {{"Type", "Eigenmode"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", {{"Eigenmode", {{"Target", 1.0}}}}}};

    IoData iodata(config, false);

    // DEFAULT is a compile-time alias — it should equal the concrete backend.
#if defined(PALACE_WITH_SLEPC)
    CHECK(iodata.solver.eigenmode.type == EigenSolverBackend::SLEPC);
#elif defined(PALACE_WITH_ARPACK)
    CHECK(iodata.solver.eigenmode.type == EigenSolverBackend::ARPACK);
#endif

    config = IoData::ConcretizeDefaults(iodata, config);

    // The resolved config should say the concrete backend name, not "Default".
    auto type_str = config["Solver"]["Eigenmode"]["Type"].get<std::string>();
#if defined(PALACE_WITH_SLEPC)
    CHECK(type_str == "SLEPc");
#elif defined(PALACE_WITH_ARPACK)
    CHECK(type_str == "ARPACK");
#endif
  }

  SECTION("WavePort eigen solver DEFAULT resolves to concrete backend")
  {
    json config = {
        {"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries",
         {{"WavePort", {{{"Index", 1}, {"Attributes", {2}}, {"Excitation", true}}}}}},
        {"Solver",
         {{"Driven",
           {{"Samples", {{{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.5}}}}}}}}};

    IoData iodata(config, false);
    config = IoData::ConcretizeDefaults(iodata, config);

    auto type_str = config["Boundaries"]["WavePort"][0]["SolverType"].get<std::string>();
#if defined(PALACE_WITH_SLEPC)
    CHECK(type_str == "SLEPc");
#elif defined(PALACE_WITH_ARPACK)
    CHECK(type_str == "ARPACK");
#endif
  }

  SECTION("mg_smooth_order resolves from solver order")
  {
    json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", {{"Order", 5}}}};

    IoData iodata(config, false);
    // max(2 * 5, 4) = 10
    CHECK(iodata.solver.linear.mg_smooth_order == 10);
    // ams_max_it defaults to solver.order = 5
    CHECK(iodata.solver.linear.ams_max_it == 5);

    config = IoData::ConcretizeDefaults(iodata, config);
    CHECK(config["Solver"]["Linear"]["MGSmoothOrder"].get<int>() == 10);
    CHECK(config["Solver"]["Linear"]["AMSMaxIts"].get<int>() == 5);
  }

  SECTION("max_size defaults to max_it")
  {
    json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", {{"Linear", {{"MaxIts", 200}}}}}};

    IoData iodata(config, false);
    CHECK(iodata.solver.linear.max_size == 200);

    config = IoData::ConcretizeDefaults(iodata, config);
    CHECK(config["Solver"]["Linear"]["MaxSize"].get<int>() == 200);
  }

  SECTION("Solver.Order default captured")
  {
    json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", json::object()}};

    IoData iodata(config, false);
    CHECK(iodata.solver.order == 1);

    config = IoData::ConcretizeDefaults(iodata, config);
    CHECK(config["Solver"]["Order"].get<int>() == 1);
  }

  SECTION("Linear default Tol, MaxIts, and orthogonalization captured")
  {
    json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", json::object()}};

    IoData iodata(config, false);
    config = IoData::ConcretizeDefaults(iodata, config);
    auto &j_linear = config["Solver"]["Linear"];

    CHECK(j_linear["Tol"].get<double>() == Catch::Approx(1.0e-6));
    CHECK(j_linear["MaxIts"].get<int>() == 100);
    CHECK(j_linear["GSOrthogonalization"].get<std::string>() == "MGS");
    // PCSide and ColumnOrdering default to the DEFAULT enum value (not a compile-time
    // alias); the concretized name records that the user accepted Palace's internal
    // default rather than naming a specific backend option.
    CHECK(j_linear["PCSide"].get<std::string>() == "Default");
    CHECK(j_linear["ColumnOrdering"].get<std::string>() == "Default");
  }

  SECTION("Transient resolves Type DEFAULT to GeneralizedAlpha")
  {
    json config = {
        {"Problem", {{"Type", "Transient"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries", json::object()},
        {"Solver",
         {{"Transient",
           {{"Excitation", "Sinusoidal"}, {"MaxTime", 1.0}, {"TimeStep", 0.01}}}}}};

    IoData iodata(config, false);
    // DEFAULT is a compile-time alias of GEN_ALPHA.
    CHECK(iodata.solver.transient.type == TimeSteppingScheme::GEN_ALPHA);

    config = IoData::ConcretizeDefaults(iodata, config);
    auto &j_transient = config["Solver"]["Transient"];
    // Must resolve to the concrete name, NOT "Default". This depends on the ordering
    // contract in PALACE_JSON_SERIALIZE_ENUM(TimeSteppingScheme, ...).
    CHECK(j_transient["Type"].get<std::string>() == "GeneralizedAlpha");
    CHECK(j_transient["Order"].get<int>() == 2);
  }

  SECTION("Eigenmode captures additional defaults")
  {
    json config = {{"Problem", {{"Type", "Eigenmode"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", {{"Eigenmode", {{"Target", 1.0}}}}}};

    IoData iodata(config, false);
    config = IoData::ConcretizeDefaults(iodata, config);
    auto &j_eigen = config["Solver"]["Eigenmode"];

    CHECK(j_eigen["Tol"].get<double>() == Catch::Approx(1.0e-6));
    CHECK(j_eigen["N"].get<int>() == 1);
    CHECK(j_eigen["Save"].get<int>() == 0);
    CHECK(j_eigen["PEPLinear"].get<bool>() == true);
    CHECK(j_eigen["Scaling"].get<bool>() == true);
    CHECK(j_eigen["StartVector"].get<bool>() == true);
    CHECK(j_eigen["StartVectorConstant"].get<bool>() == false);
    CHECK(j_eigen["MassOrthogonal"].get<bool>() == false);
  }

  SECTION("Round-trip: resolved config validates and re-parses identically")
  {
    // Start from a sparse Electrostatic config with one non-default (Order=3) and one
    // user-set Linear field (MaxIts=250). Everything else resolves from defaults.
    json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", {{"Order", 3}, {"Linear", {{"MaxIts", 250}}}}}};

    IoData iodata1(config, false);
    config = IoData::ConcretizeDefaults(iodata1, config);

    // The resolved config must pass schema validation; otherwise a user cannot
    // actually re-run Palace on the produced file.
    std::string err = ValidateConfig(config);
    INFO("schema validation error: " << err);
    CHECK(err.empty());

    IoData iodata2(config, false);

    CHECK(iodata2.solver.order == iodata1.solver.order);
    const auto &l1 = iodata1.solver.linear;
    const auto &l2 = iodata2.solver.linear;
    CHECK(l2.type == l1.type);
    CHECK(l2.krylov_solver == l1.krylov_solver);
    CHECK(l2.tol == l1.tol);
    CHECK(l2.max_it == l1.max_it);
    CHECK(l2.max_size == l1.max_size);
    CHECK(l2.initial_guess == l1.initial_guess);
    CHECK(l2.pc_mat_shifted == l1.pc_mat_shifted);
    CHECK(l2.mg_smooth_aux == l1.mg_smooth_aux);
    CHECK(l2.mg_smooth_order == l1.mg_smooth_order);
    CHECK(l2.mg_cycle_it == l1.mg_cycle_it);
    CHECK(l2.ams_singular_op == l1.ams_singular_op);
    CHECK(l2.ams_max_it == l1.ams_max_it);
    CHECK(l2.amg_agg_coarsen == l1.amg_agg_coarsen);
    CHECK(l2.pc_mat_sym == l1.pc_mat_sym);
    CHECK(l2.reorder_reuse == l1.reorder_reuse);
    CHECK(l2.pc_side == l1.pc_side);
    CHECK(l2.sym_factorization == l1.sym_factorization);
    CHECK(l2.gs_orthog == l1.gs_orthog);
  }

  SECTION("Round-trip: Eigenmode DEFAULT backend resolved concretely")
  {
    json config = {{"Problem", {{"Type", "Eigenmode"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", {{"Eigenmode", {{"Target", 1.0}}}}}};

    IoData iodata1(config, false);
    config = IoData::ConcretizeDefaults(iodata1, config);

    std::string err = ValidateConfig(config);
    INFO("schema validation error: " << err);
    CHECK(err.empty());

    IoData iodata2(config, false);
    CHECK(iodata2.solver.eigenmode.type == iodata1.solver.eigenmode.type);
    CHECK(iodata2.solver.eigenmode.type != EigenSolverBackend::DEFAULT);
    CHECK(iodata2.solver.eigenmode.target == iodata1.solver.eigenmode.target);
    CHECK(iodata2.solver.eigenmode.target_upper == iodata1.solver.eigenmode.target_upper);
    CHECK(iodata2.solver.eigenmode.max_it == iodata1.solver.eigenmode.max_it);
    CHECK(iodata2.solver.eigenmode.max_size == iodata1.solver.eigenmode.max_size);
  }

  SECTION("Round-trip: Transient DEFAULT scheme resolved concretely")
  {
    json config = {
        {"Problem", {{"Type", "Transient"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries", json::object()},
        {"Solver",
         {{"Transient",
           {{"Excitation", "Sinusoidal"}, {"MaxTime", 1.0}, {"TimeStep", 0.01}}}}}};

    IoData iodata1(config, false);
    config = IoData::ConcretizeDefaults(iodata1, config);

    std::string err = ValidateConfig(config);
    INFO("schema validation error: " << err);
    CHECK(err.empty());

    IoData iodata2(config, false);
    CHECK(iodata2.solver.transient.type == iodata1.solver.transient.type);
    CHECK(iodata2.solver.transient.type != TimeSteppingScheme::DEFAULT);
    CHECK(iodata2.solver.transient.excitation == iodata1.solver.transient.excitation);
    CHECK(iodata2.solver.transient.max_t == iodata1.solver.transient.max_t);
    CHECK(iodata2.solver.transient.delta_t == iodata1.solver.transient.delta_t);
  }
}
