// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>
#include <fmt/core.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"
#include "utils/jsonschema.hpp"

using json = nlohmann::json;
using namespace palace;

TEST_CASE("Config Boundary Ports", "[config][Serial]")
{
  auto filename = fmt::format("{}/{}", PALACE_TEST_DIR, "config/boundary_configs.json");
  auto jsonstream = PreprocessFile(filename.c_str());  // Apply custom palace json
  auto config = json::parse(jsonstream);

  {
    // Basic passing config with bool excitation.
    auto entry = config.find("boundaries_1_pass")->find("Boundaries");
    config::BoundaryData boundary_ex_bool(*entry);

    // Check simple parsing & defaults:
    CHECK(boundary_ex_bool.lumpedport.at(1).active);
    CHECK(boundary_ex_bool.lumpedport.at(3).active == false);
    CHECK(boundary_ex_bool.lumpedport.at(1).excitation != 0);
    CHECK(boundary_ex_bool.lumpedport.at(3).excitation == 0);
    CHECK(boundary_ex_bool.waveport.at(5).excitation == 0);
    CHECK(boundary_ex_bool.waveport.at(6).excitation == 0);

    // Equivalent config with int excitation.
    auto entry_int = config.find("boundaries_1_pass_excitation_int")->find("Boundaries");
    config::BoundaryData boundary_ex_int(*entry_int);

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
    auto entry = config.find("boundaries_negative_excitation_1")->find("Boundaries");
    CHECK_THROWS(config::BoundaryData(*entry));
  }
  {
    auto entry = config.find("boundaries_negative_excitation_2")->find("Boundaries");
    CHECK_THROWS(config::BoundaryData(*entry));
  }
  // Index Specification.
  {
    auto entry = config.find("boundaries_repeated_index_lumped")->find("Boundaries");
    CHECK_THROWS(config::BoundaryData(*entry));
  }
  {
    auto entry = config.find("boundaries_repeated_index_wave")->find("Boundaries");
    CHECK_THROWS(config::BoundaryData(*entry));
  }
  {
    auto entry = config.find("boundaries_repeated_index_mixed")->find("Boundaries");
    CHECK_THROWS(config::BoundaryData(*entry));
  }
  {
    auto entry = config.find("boundaries_mislabeled_index_1")->find("Boundaries");
    CHECK_THROWS(config::BoundaryData(*entry));
  }
  {
    auto entry = config.find("boundaries_mislabeled_index_2")->find("Boundaries");
    CHECK_THROWS(config::BoundaryData(*entry));
  }
  // Mark single excitation index.
  {
    auto entry = config.find("boundaries_upgrade_excitation_idx_1")->find("Boundaries");
    config::BoundaryData boundary_data(*entry);
    CHECK(boundary_data.lumpedport.at(1).excitation == 0);
    CHECK(boundary_data.lumpedport.at(2).excitation == 2);
    CHECK(boundary_data.waveport.at(4).excitation == 0);
    CHECK(boundary_data.waveport.at(5).excitation == 0);
  }
  {
    auto entry = config.find("boundaries_upgrade_excitation_idx_2")->find("Boundaries");
    config::BoundaryData boundary_data(*entry);
    CHECK(boundary_data.lumpedport.at(1).excitation == 0);
    CHECK(boundary_data.lumpedport.at(2).excitation == 2);
    CHECK(boundary_data.waveport.at(4).excitation == 0);
    CHECK(boundary_data.waveport.at(5).excitation == 0);
  }
  {
    auto entry = config.find("boundaries_upgrade_excitation_idx_3")->find("Boundaries");
    config::BoundaryData boundary_data(*entry);
    CHECK(boundary_data.lumpedport.at(1).excitation == 0);
    CHECK(boundary_data.lumpedport.at(2).excitation == 0);
    CHECK(boundary_data.waveport.at(4).excitation == 4);
    CHECK(boundary_data.waveport.at(5).excitation == 0);
  }
  {
    auto entry = config.find("boundaries_upgrade_excitation_idx_4")->find("Boundaries");
    config::BoundaryData boundary_data(*entry);
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
    auto save_indices = std::vector<std::size_t>{0, 2, 4, 6, 8, 10};
    {
      // Top level configuration
      auto driven = config.find("driven_base_uniform_sample")->find("Driven");
      config::DrivenSolverData driven_solver(*driven);

      for (std::size_t i = 0; i < sample_f.size(); ++i)
      {
        CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
      }
      CHECK(driven_solver.save_indices == save_indices);
      CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
    }
    {
      // Equivalent to top level from within Samples, deduplicates
      auto driven = config.find("driven_uniform_freq_step")->find("Driven");
      config::DrivenSolverData driven_solver(*driven);

      for (std::size_t i = 0; i < sample_f.size(); ++i)
      {
        CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
      }
      CHECK(driven_solver.save_indices == save_indices);
      CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
    }
  }
  {
    // Specification through number of points rather than step size
    auto driven = config.find("driven_uniform_nsample")->find("Driven");
    config::DrivenSolverData driven_solver(*driven);

    auto sample_f = std::vector{0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
    auto save_indices = std::vector<std::size_t>{0, 2, 4, 6, 8};

    for (std::size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
    CHECK(driven_solver.save_indices == save_indices);
    CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
  }
  {
    // Combining two different linear sample resolutions
    auto driven = config.find("driven_paired_uniform_sample")->find("Driven");
    config::DrivenSolverData driven_solver(*driven);

    auto sample_f = std::vector{0.0, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 7.5, 10.0};
    auto save_indices = std::vector<std::size_t>{
        0, 2, 4, 5, 6, 7, 8};  // 0.0, 0.5, 1.0, 2.5, 5.0, 7.5, 10.0

    for (std::size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
    CHECK(driven_solver.save_indices == save_indices);
    CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
  }
  {
    // Combining two different linear sample resolutions
    auto driven = config.find("driven_uniform_with_point")->find("Driven");
    config::DrivenSolverData driven_solver(*driven);

    auto sample_f = std::vector{0.0, 0.125, 0.15,  0.25, 0.35,  0.375,
                                0.5, 0.55,  0.625, 0.75, 0.875, 1.0};
    auto save_indices = std::vector<std::size_t>{0, 2, 3, 4, 6, 7, 9, 11};
    auto prom_indices = std::vector<std::size_t>{0, sample_f.size() - 1, 2, 4, 7};

    for (std::size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(driven_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
    CHECK(driven_solver.save_indices == save_indices);
    CHECK(driven_solver.prom_indices == prom_indices);
  }
  {
    // Combining two different linear sample resolutions
    auto driven = config.find("driven_log_with_point")->find("Driven");
    config::DrivenSolverData driven_solver(*driven);

    auto sample_f = std::vector{0.1,  0.15, 0.1778279410038923, 0.31622776601683794,
                                0.35, 0.55, 0.5623413251903491, 1.0};
    auto save_indices = std::vector<std::size_t>{0, 1, 3, 4, 5, 7};
    auto prom_indices = std::vector<std::size_t>{0, sample_f.size() - 1, 1, 4, 5};

    for (std::size_t i = 0; i < sample_f.size(); ++i)
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
    auto driven = config.find(c)->find("Driven");
    CHECK_THROWS(config::DrivenSolverData(*driven));
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
  constexpr double delta_eps = 1.0e-6;  // Precision in angle comparisons (rad)

  SECTION("Default constructor")
  {
    config::FarFieldPostData data;

    CHECK(data.attributes.empty());
    CHECK(data.thetaphis.empty());
  }

  SECTION("Basic setup with attributes only")
  {
    // This should produce a warning because there is no target point.
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
      CHECK(data.thetaphis.size() == 1);  // All should collapse to one pole.
    }

    // Test phi periodicity: φ and φ+360° are same point.
    {
      json farfield = {{"Attributes", {1}},
                       {"ThetaPhis", {{45.0, 30.0}, {45.0, 390.0}}}};  // 390° = 30° + 360°
      config::FarFieldPostData data(farfield);
      CHECK(data.thetaphis.size() == 1);  // Should be deduplicated.
    }

    // Test theta reflection: (θ, φ) ≡ (180°-θ, φ+180°).
    {
      json farfield = {
          {"Attributes", {1}},
          {"ThetaPhis",
           {{60.0, 45.0}, {120.0, 225.0}}}};  // 120° = 180°-60°, 225° = 45°+180°
      config::FarFieldPostData data(farfield);
      CHECK(data.thetaphis.size() == 1);  // Should be deduplicated.
    }
  }

  SECTION("Combined NSample and ThetaPhis")
  {
    json farfield = {{"Attributes", {1}}, {"NSample", 10}, {"ThetaPhis", {{33.0, 22.0}}}};
    config::FarFieldPostData data(farfield);

    CHECK(data.thetaphis.size() == 11);  // 10 from NSample + 1 from ThetaPhis.
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
