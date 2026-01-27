// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"

using json = nlohmann::json;
using namespace palace;

TEST_CASE("Config Boundary Ports", "[config][Serial]")
{
  auto filename = fmt::format("{}/{}", PALACE_TEST_DIR, "config/boundary_configs.json");
  auto jsonstream = PreprocessFile(filename.c_str());  // Apply custom palace json
  auto config = json::parse(jsonstream);

  {
    // Basic passing config with bool excitation.
    config::BoundaryData boundary_ex_bool;
    REQUIRE_NOTHROW(boundary_ex_bool.SetUp(*config.find("boundaries_1_pass")));

    // Check simple parsing & defaults:
    CHECK(boundary_ex_bool.lumpedport.at(1).active);
    CHECK(boundary_ex_bool.lumpedport.at(3).active == false);
    CHECK(boundary_ex_bool.lumpedport.at(1).excitation != 0);
    CHECK(boundary_ex_bool.lumpedport.at(3).excitation == 0);
    CHECK(boundary_ex_bool.waveport.at(5).excitation == 0);
    CHECK(boundary_ex_bool.waveport.at(6).excitation == 0);

    // Equivalent config with int excitation.
    config::BoundaryData boundary_ex_int;
    REQUIRE_NOTHROW(
        boundary_ex_int.SetUp(*config.find("boundaries_1_pass_excitation_int")));

    // FUTURE: Default equality is C++20.
    // CHECK(boundary_ex_bool == boundary_ex_int);

    REQUIRE(boundary_ex_bool.lumpedport.size() == boundary_ex_bool.lumpedport.size());
    auto it_int = boundary_ex_bool.lumpedport.begin();
    auto it_bool = boundary_ex_bool.lumpedport.begin();
    for (; it_int != boundary_ex_bool.lumpedport.end(); it_int++, it_bool++)
    {
      CHECK(it_bool->first == it_int->first);  // Order is same as indices are
      CHECK(it_bool->second.excitation == it_int->second.excitation);
    }
  }
  // Excitation Specification.
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(boundary_data.SetUp(*config.find("boundaries_negative_excitation_1")));
  }
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(boundary_data.SetUp(*config.find("boundaries_negative_excitation_2")));
  }
  // Index Specification.
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(boundary_data.SetUp(*config.find("boundaries_repeated_index_lumped")));
  }
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(boundary_data.SetUp(*config.find("boundaries_repeated_index_wave")));
  }
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(boundary_data.SetUp(*config.find("boundaries_repeated_index_mixed")));
  }
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(boundary_data.SetUp(*config.find("boundaries_negative_index_1")));
  }
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(boundary_data.SetUp(*config.find("boundaries_negative_index_2")));
  }
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(boundary_data.SetUp(*config.find("boundaries_mislabeled_index_1")));
  }
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(boundary_data.SetUp(*config.find("boundaries_mislabeled_index_2")));
  }
  // Mark single excitation index.
  {
    config::BoundaryData boundary_data;
    CHECK_NOTHROW(boundary_data.SetUp(*config.find("boundaries_upgrade_excitation_idx_1")));
    CHECK(boundary_data.lumpedport.at(1).excitation == 0);
    CHECK(boundary_data.lumpedport.at(2).excitation == 2);
    CHECK(boundary_data.waveport.at(4).excitation == 0);
    CHECK(boundary_data.waveport.at(5).excitation == 0);
  }
  {
    config::BoundaryData boundary_data;
    CHECK_NOTHROW(boundary_data.SetUp(*config.find("boundaries_upgrade_excitation_idx_2")));
    CHECK(boundary_data.lumpedport.at(1).excitation == 0);
    CHECK(boundary_data.lumpedport.at(2).excitation == 2);
    CHECK(boundary_data.waveport.at(4).excitation == 0);
    CHECK(boundary_data.waveport.at(5).excitation == 0);
  }
  {
    config::BoundaryData boundary_data;
    CHECK_NOTHROW(boundary_data.SetUp(*config.find("boundaries_upgrade_excitation_idx_3")));
    CHECK(boundary_data.lumpedport.at(1).excitation == 0);
    CHECK(boundary_data.lumpedport.at(2).excitation == 0);
    CHECK(boundary_data.waveport.at(4).excitation == 4);
    CHECK(boundary_data.waveport.at(5).excitation == 0);
  }
  {
    config::BoundaryData boundary_data;
    CHECK_NOTHROW(boundary_data.SetUp(*config.find("boundaries_upgrade_excitation_idx_4")));
    CHECK(boundary_data.lumpedport.at(1).excitation == 1);
    CHECK(boundary_data.lumpedport.at(2).excitation == 0);
    CHECK(boundary_data.waveport.at(4).excitation == 1);
    CHECK(boundary_data.waveport.at(5).excitation == 0);
  }
}

TEST_CASE("Config Driven Solver", "[config][Serial]")
{
  auto filename = fmt::format("{}/{}", PALACE_TEST_DIR, "config/solver_configs.json");
  auto jsonstream = PreprocessFile(filename.c_str());  // Apply custom palace json
  auto config = json::parse(jsonstream);

  using namespace Catch::Matchers;

  constexpr double delta_eps = 1.0e-9;  // Precision in frequency comparisons (Hz)
  {
    auto sample_f = std::vector{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1};
    auto save_indices = std::vector<size_t>{0, 2, 4, 6, 8, 10};
    {
      // Top level configuration
      config::DrivenSolverData driven_solver;
      REQUIRE_NOTHROW(driven_solver.SetUp(*config.find("driven_base_uniform_sample")));

      for (size_t i = 0; i < sample_f.size(); ++i)
      {
        CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
      }
      CHECK(driven_solver.save_indices == save_indices);
      CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
    }
    {
      // Equivalent to top level from within Samples, deduplicates
      config::DrivenSolverData driven_solver;
      REQUIRE_NOTHROW(driven_solver.SetUp(*config.find("driven_uniform_freq_step")));

      for (size_t i = 0; i < sample_f.size(); ++i)
      {
        CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
      }
      CHECK(driven_solver.save_indices == save_indices);
      CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
    }
  }
  {
    // Specification through number of points rather than step size
    config::DrivenSolverData driven_solver;
    REQUIRE_NOTHROW(driven_solver.SetUp(*config.find("driven_uniform_nsample")));

    auto sample_f = std::vector{0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
    auto save_indices = std::vector<size_t>{0, 2, 4, 6, 8};

    for (size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
    CHECK(driven_solver.save_indices == save_indices);
    CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
  }
  {
    // Combining two different linear sample resolutions
    config::DrivenSolverData driven_solver;
    REQUIRE_NOTHROW(driven_solver.SetUp(*config.find("driven_paired_uniform_sample")));

    auto sample_f = std::vector{0.0, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 7.5, 10.0};
    auto save_indices =
        std::vector<size_t>{0, 2, 4, 5, 6, 7, 8};  // 0.0, 0.5, 1.0, 2.5, 5.0, 7.5, 10.0

    for (size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
    CHECK(driven_solver.save_indices == save_indices);
    CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
  }
  {
    // Combining two different linear sample resolutions
    config::DrivenSolverData driven_solver;
    REQUIRE_NOTHROW(driven_solver.SetUp(*config.find("driven_uniform_with_point")));

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
  {
    // Combining two different linear sample resolutions
    config::DrivenSolverData driven_solver;
    REQUIRE_NOTHROW(driven_solver.SetUp(*config.find("driven_log_with_point")));

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
  std::vector<std::string> invalid_configs = {"driven_empty",
                                              "driven_mismatch_type_1",
                                              "driven_mismatch_type_2",
                                              "driven_mismatch_type_3",
                                              "driven_invalid_log_range",
                                              "driven_uniform_with_point_invalid_save"};
  for (auto c : invalid_configs)
  {
    config::DrivenSolverData driven_solver;
    CHECK_THROWS(config::DrivenSolverData().SetUp(*config.find(c)));
  }
}

TEST_CASE("Config Linear Solver MaxIts", "[config][Serial]")
{
  auto filename = fmt::format("{}/{}", PALACE_TEST_DIR, "config/solver_configs.json");
  auto jsonstream = PreprocessFile(filename.c_str());
  auto config = json::parse(jsonstream);

  SECTION("Linear solver MaxIts = 0 should throw")
  {
    config::LinearSolverData linear_solver;
    CHECK_THROWS_WITH(linear_solver.SetUp(*config.find("linear_maxits_zero")),
                      Catch::Matchers::ContainsSubstring("MaxIts"));
  }

  SECTION("Linear solver with valid MaxIts should succeed")
  {
    config::LinearSolverData linear_solver;
    REQUIRE_NOTHROW(linear_solver.SetUp(*config.find("linear_maxits_one")));
  }
}

TEST_CASE("FarField", "[config][Serial]")
{
  constexpr double delta_eps = 1.0e-6;  // Precision in angle comparisons (rad)

  SECTION("Missing FarField section")
  {
    json postpro = json::object();
    config::FarFieldPostData data;

    data.SetUp(postpro);

    CHECK(data.attributes.empty());
    CHECK(data.thetaphis.empty());
  }

  SECTION("Basic setup with attributes only")
  {
    // This should produce a warning because there is no target point.
    json postpro = {{"FarField", {{"Attributes", {1, 3, 5}}}}};
    config::FarFieldPostData data;

    data.SetUp(postpro);

    CHECK(data.attributes == std::vector<int>{1, 3, 5});
    CHECK(data.thetaphis.empty());
  }

  SECTION("ThetaPhis conversion to radians")
  {
    json postpro = {
        {"FarField", {{"Attributes", {1}}, {"ThetaPhis", {{0.0, 0.0}, {90.0, 180.0}}}}}};
    config::FarFieldPostData data;

    data.SetUp(postpro);

    CHECK(data.thetaphis.size() == 2);
    CHECK(data.thetaphis[0].first == Catch::Approx(0.0).margin(delta_eps));
    CHECK(data.thetaphis[0].second == Catch::Approx(0.0).margin(delta_eps));
    CHECK(data.thetaphis[1].first == Catch::Approx(M_PI / 2).margin(delta_eps));
    CHECK(data.thetaphis[1].second == Catch::Approx(M_PI).margin(delta_eps));
  }

  SECTION("Duplicate removal")
  {
    json postpro = {
        {"FarField",
         {{"Attributes", {1}},
          {"ThetaPhis", {{0.0, 0.0}, {90.0, 180.0}, {0.0, 0.0}, {90.0, 180.0}}}}}};
    config::FarFieldPostData data;

    data.SetUp(postpro);

    CHECK(data.thetaphis.size() == 2);
  }

  SECTION("Spherical coordinate duplicate removal")
  {
    // Test pole singularity: (0°, any φ) should be treated as same point.
    {
      json postpro = {
          {"FarField",
           {{"Attributes", {1}}, {"ThetaPhis", {{0.0, 0.0}, {0.0, 90.0}, {0.0, 180.0}}}}}};
      config::FarFieldPostData data;
      data.SetUp(postpro);
      CHECK(data.thetaphis.size() == 1);  // All should collapse to one pole.
    }

    // Test phi periodicity: φ and φ+360° are same point.
    {
      json postpro = {
          {"FarField",
           {{"Attributes", {1}},
            {"ThetaPhis", {{45.0, 30.0}, {45.0, 390.0}}}}}};  // 390° = 30° + 360°
      config::FarFieldPostData data;
      data.SetUp(postpro);
      CHECK(data.thetaphis.size() == 1);  // Should be deduplicated.
    }

    // Test theta reflection: (θ, φ) ≡ (180°-θ, φ+180°).
    {
      json postpro = {
          {"FarField",
           {{"Attributes", {1}},
            {"ThetaPhis",
             {{60.0, 45.0}, {120.0, 225.0}}}}}};  // 120° = 180°-60°, 225° = 45°+180°
      config::FarFieldPostData data;
      data.SetUp(postpro);
      CHECK(data.thetaphis.size() == 1);  // Should be deduplicated.
    }
  }

  SECTION("Combined NSample and ThetaPhis")
  {
    json postpro = {
        {"FarField",
         {{"Attributes", {1}}, {"NSample", 10}, {"ThetaPhis", {{33.0, 22.0}}}}}};
    config::FarFieldPostData data;

    data.SetUp(postpro);

    CHECK(data.thetaphis.size() == 11);  // 10 from NSample + 1 from ThetaPhis.
  }

  SECTION("NSample sphere sampling")
  {
    for (int nsample : {2, 6, 10, 15, 20, 25, 64800})
    {
      json postpro = {{"FarField", {{"Attributes", {1}}, {"NSample", nsample}}}};
      config::FarFieldPostData data;
      data.SetUp(postpro);

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
      json postpro = {{"FarField", {{"Attributes", {1}}, {"NSample", 0}}}};
      config::FarFieldPostData data;
      data.SetUp(postpro);
      CHECK(data.thetaphis.empty());
    }

    // NSample = 1 produces two points (the poles).
    {
      json postpro = {{"FarField", {{"Attributes", {1}}, {"NSample", 1}}}};
      config::FarFieldPostData data;
      data.SetUp(postpro);
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

  SECTION("Invalid")
  {
    auto req = GENERATE(true, false);
    CHECK_THROWS(config::ParseStringAsDirection("a", req));
    CHECK_THROWS(config::ParseStringAsDirection("+a", req));
    CHECK_THROWS(config::ParseStringAsDirection("-a", req));
    CHECK_THROWS(config::ParseStringAsDirection("xx", req));
    CHECK_THROWS(config::ParseStringAsDirection("~x", req));
    CHECK_THROWS(config::ParseStringAsDirection("x+", req));
    CHECK_THROWS(config::ParseStringAsDirection("xy", req));
    CHECK_THROWS(config::ParseStringAsDirection("xyz", req));
    CHECK_THROWS(config::ParseStringAsDirection("abc", req));
  }

  CHECK_THROWS(config::ParseStringAsDirection("", true));
  CHECK_NOTHROW(config::ParseStringAsDirection("", false));
}
