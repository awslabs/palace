#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
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
    // Basic passing config with bool excitaiotn
    config::BoundaryData boundary_ex_bool;
    REQUIRE_NOTHROW(boundary_ex_bool.SetUp(*config.find("boundaries_1_pass")));

    // Check simple parsing & defaults:
    CHECK(boundary_ex_bool.lumpedport.at(1).active);
    CHECK(boundary_ex_bool.lumpedport.at(3).active == false);
    CHECK(boundary_ex_bool.lumpedport.at(1).excitation);
    CHECK(boundary_ex_bool.lumpedport.at(3).excitation == false);
    CHECK(boundary_ex_bool.waveport.at(5).excitation == false);
    CHECK(boundary_ex_bool.waveport.at(6).excitation == false);

    // Equivalent config with int excitaiotn
    config::BoundaryData boundary_ex_int;
    REQUIRE_NOTHROW(
        boundary_ex_int.SetUp(*config.find("boundaries_1_pass_excitation_int")));

    // FUTURE: Default equality is C++20
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
  // Excitation Specification
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(boundary_data.SetUp(*config.find("boundaries_negative_excitation_1")));
  }
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(boundary_data.SetUp(*config.find("boundaries_negative_excitation_2")));
  }
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(
        boundary_data.SetUp(*config.find("boundaries_excitation_no_mix_int_bool_1")));
  }
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(
        boundary_data.SetUp(*config.find("boundaries_excitation_no_mix_int_bool_2")));
  }
  {
    config::BoundaryData boundary_data;
    CHECK_THROWS(
        boundary_data.SetUp(*config.find("boundaries_excitation_no_mix_int_bool_3")));
  }
  // Index Specification
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
  // Mark single excitation index
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
    CHECK(boundary_data.lumpedport.at(2).excitation == 1);
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