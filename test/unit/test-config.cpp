#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <fmt/format.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"

using json = nlohmann::json;
using namespace palace;

TEST_CASE("Config Boundary Ports", "[config]")
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
      CHECK(it_bool->first == it_int->first);  // Order is same as indicies are
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

TEST_CASE("Config Driven Solver", "[config]")
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

      CHECK_THAT(driven_solver.sample_f, Approx(sample_f).margin(delta_eps));
      CHECK(driven_solver.save_indices == save_indices);
      CHECK(driven_solver.prom_indices == std::vector{0, sample_f.size() - 1});
    }
    {
      // Equivalent to top level from within Samples, deduplicates
      config::DrivenSolverData driven_solver;
      REQUIRE_NOTHROW(driven_solver.SetUp(*config.find("driven_uniform_freq_step")));

      CHECK_THAT(driven_solver.sample_f, Approx(sample_f).margin(delta_eps));
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

    CHECK_THAT(driven_solver.sample_f, Approx(sample_f).margin(delta_eps));
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

    CHECK_THAT(driven_solver.sample_f, Approx(sample_f).margin(delta_eps));
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

    CHECK_THAT(driven_solver.sample_f, Approx(sample_f).margin(delta_eps));
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

    CHECK_THAT(driven_solver.sample_f, Approx(sample_f).margin(delta_eps));
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