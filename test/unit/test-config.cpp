// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <set>
#include <string>
#include <vector>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "embedded_schema.hpp"
#include "linalg/ksp.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"
#include "utils/jsonschema.hpp"

using json = nlohmann::json;
using namespace palace;

namespace
{

// Returns the names of properties declared in the given schema scope that the
// concretized config does NOT contain, excluding (a) properties marked `required`
// in the schema (the user must provide those at parse time, no default applies)
// and (b) any names in `skip` (intentional omissions — e.g. `Lc` is computed from
// the mesh at runtime, `MaterialAxes`/`Center` carry opt-in semantics where absence
// is meaningful).
//
// `schema_filename` is the key in `schema::GetSchemaMap()` (e.g. "config/model.json").
// `pointer` is a JSON Pointer into that schema document selecting the scope to walk
// (e.g. "" for the top of model.json, "/properties/Refinement" for the nested
// Refinement object, or "/properties/Conductivity/items" for an array element).
//
// Use this in round-trip tests to catch the case where someone adds a new optional
// schema property without wiring `ConcretizeDefaults` to emit it.
std::vector<std::string> SchemaCoverageGaps(const std::string &schema_filename,
                                            const std::string &pointer,
                                            const json &concretized_subtree,
                                            const std::set<std::string> &skip = {})
{
  const auto &schema_map = schema::GetSchemaMap();
  auto it = schema_map.find(schema_filename);
  REQUIRE(it != schema_map.end());
  const json schema = json::parse(it->second);
  const json scope = schema.at(json::json_pointer(pointer));
  REQUIRE(scope.contains("properties"));
  std::set<std::string> required_set;
  if (auto rit = scope.find("required"); rit != scope.end())
  {
    for (const auto &r : *rit)
    {
      required_set.insert(r.get<std::string>());
    }
  }
  std::vector<std::string> gaps;
  for (const auto &[name, _sub] : scope.at("properties").items())
  {
    if (required_set.count(name) || skip.count(name))
    {
      continue;
    }
    if (!concretized_subtree.is_object() || !concretized_subtree.contains(name))
    {
      gaps.push_back(name);
    }
  }
  return gaps;
}

}  // namespace

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
    CHECK(GetPreconditionerMatrixSymmetry(iodata) == MatrixSymmetry::SPD);

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
  }

  SECTION("Omitted Output resolves to default and concretizes (issue #745)")
  {
    // The schema marks Problem.Output optional, so a config without it must parse
    // without aborting and fall back to the documented "postpro" default.
    json config = {{"Problem", {{"Type", "Electrostatic"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", json::object()}};

    IoData iodata(config, false);
    CHECK(iodata.problem.output == "postpro");

    // ConcretizeDefaults must emit the resolved, non-empty value back to JSON.
    config = IoData::ConcretizeDefaults(iodata, config);
    CHECK(config["Problem"]["Output"].get<std::string>() == "postpro");
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
    CHECK(GetPreconditionerMatrixSymmetry(iodata) == MatrixSymmetry::SPD);

    config = IoData::ConcretizeDefaults(iodata, config);
    CHECK(config["Solver"]["Linear"]["Type"].get<std::string>() == "AMS");
    CHECK(config["Solver"]["Linear"]["AMSSingularOperator"].get<int>() == 1);
    CHECK(config["Solver"]["Linear"]["AMSMaxIts"].get<int>() == 1);
    CHECK(config["Solver"]["Linear"]["MGCycleIts"].get<int>() == 1);
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
    CHECK(iodata.solver.transient.type == TimeSteppingScheme::GEN_ALPHA);

    config = IoData::ConcretizeDefaults(iodata, config);
    auto &j_transient = config["Solver"]["Transient"];
    // Must resolve to the concrete name, NOT "Default". This depends on the ordering
    // contract in PALACE_JSON_SERIALIZE_ENUM(TimeSteppingScheme, ...).
    CHECK(j_transient["Type"].get<std::string>() == "GeneralizedAlpha");
    // Order is not consumed by GeneralizedAlpha/RungeKutta and triggers a warning at
    // parse time if emitted; ConcretizeTransient must omit it for these schemes.
    CHECK_FALSE(j_transient.contains("Order"));
    CHECK_FALSE(j_transient.contains("RelTol"));
    CHECK_FALSE(j_transient.contains("AbsTol"));
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
    CHECK(l2.reorder_reuse == l1.reorder_reuse);
    CHECK(l2.pc_side == l1.pc_side);
    CHECK(l2.sym_factorization == l1.sym_factorization);
    CHECK(l2.gs_orthog == l1.gs_orthog);
  }

  SECTION("Round-trip: Model, Refinement, and Boundaries are reproducible")
  {
    json config = {
        {"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
        {"Model",
         {{"Mesh", "test.msh"},
          {"L0", 1.0e-3},
          {"MakeSimplex", true},
          {"Partitioning", "parts.txt"},
          {"Refinement", {{"MaxIts", 4}, {"MaxSize", 1000}, {"UpdateFraction", 0.5}}}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries",
         {{"Absorbing", {{"Attributes", {2}}, {"Order", 2}}},
          {"Conductivity", {{{"Attributes", {3}}, {"Conductivity", 5.8e7}}}},
          {"Impedance",
           {{{"Attributes", {4}}, {"Rs", 50.0}, {"Ls", 1.0e-9}, {"Cs", 1.0e-12}}}},
          {"LumpedPort",
           {{{"Index", 1},
             {"Attributes", {5}},
             {"Direction", "+X"},
             {"R", 50.0},
             {"L", 1.0e-9},
             {"Excitation", true},
             {"Active", false}}}},
          {"Periodic",
           {{"FloquetWaveVector", {0.1, 0.2, 0.3}},
            {"BoundaryPairs",
             {{{"DonorAttributes", {10}}, {"ReceiverAttributes", {11}}}}}}},
          {"Postprocessing",
           {{"SurfaceFlux",
             {{{"Index", 1},
               {"Attributes", {6}},
               {"Type", "Electric"},
               {"TwoSided", true}}}},
            {"Dielectric",
             {{{"Index", 1},
               {"Attributes", {7}},
               {"Type", "MA"},
               {"Thickness", 1.0e-9},
               {"Permittivity", 4.0},
               {"LossTan", 0.002}}}}}}}},
        {"Solver", {{"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.1}}}}}};

    IoData iodata1(config, false);
    config = IoData::ConcretizeDefaults(iodata1, config);

    std::string err = ValidateConfig(config);
    INFO("schema validation error: " << err);
    CHECK(err.empty());

    IoData iodata2(config, false);

    const auto &m1 = iodata1.model;
    const auto &m2 = iodata2.model;
    CHECK(m2.L0 == m1.L0);
    CHECK(m2.Lc == m1.Lc);
    CHECK(m2.remove_curvature == m1.remove_curvature);
    CHECK(m2.make_simplex == m1.make_simplex);
    CHECK(m2.make_hex == m1.make_hex);
    CHECK(m2.reorder_elements == m1.reorder_elements);
    CHECK(m2.clean_unused_elements == m1.clean_unused_elements);
    CHECK(m2.crack_bdr_elements == m1.crack_bdr_elements);
    CHECK(m2.refine_crack_elements == m1.refine_crack_elements);
    CHECK(m2.crack_displ_factor == m1.crack_displ_factor);
    CHECK(m2.add_bdr_elements == m1.add_bdr_elements);
    CHECK(m2.export_prerefined_mesh == m1.export_prerefined_mesh);
    CHECK(m2.reorient_tet_mesh == m1.reorient_tet_mesh);
    CHECK(m2.partitioning == m1.partitioning);
    const auto &r1 = m1.refinement;
    const auto &r2 = m2.refinement;
    CHECK(r2.tol == r1.tol);
    CHECK(r2.max_it == r1.max_it);
    CHECK(r2.max_size == r1.max_size);
    CHECK(r2.nonconformal == r1.nonconformal);
    CHECK(r2.max_nc_levels == r1.max_nc_levels);
    CHECK(r2.update_fraction == r1.update_fraction);
    CHECK(r2.maximum_imbalance == r1.maximum_imbalance);
    CHECK(r2.save_adapt_iterations == r1.save_adapt_iterations);
    CHECK(r2.save_adapt_mesh == r1.save_adapt_mesh);
    CHECK(r2.uniform_ref_levels == r1.uniform_ref_levels);
    CHECK(r2.ser_uniform_ref_levels == r1.ser_uniform_ref_levels);
    CHECK(iodata2.boundaries.farfield.order == iodata1.boundaries.farfield.order);
    REQUIRE(iodata1.boundaries.conductivity.size() == 1);
    REQUIRE(iodata2.boundaries.conductivity.size() == 1);
    CHECK(iodata2.boundaries.conductivity[0].h == iodata1.boundaries.conductivity[0].h);
    CHECK(iodata2.boundaries.conductivity[0].external ==
          iodata1.boundaries.conductivity[0].external);
    REQUIRE(iodata1.boundaries.impedance.size() == 1);
    REQUIRE(iodata2.boundaries.impedance.size() == 1);
    CHECK(iodata2.boundaries.impedance[0].Rs == iodata1.boundaries.impedance[0].Rs);
    CHECK(iodata2.boundaries.impedance[0].Ls == iodata1.boundaries.impedance[0].Ls);
    CHECK(iodata2.boundaries.impedance[0].Cs == iodata1.boundaries.impedance[0].Cs);
    REQUIRE(iodata2.boundaries.lumpedport.count(1) == 1);
    const auto &lp1 = iodata1.boundaries.lumpedport.at(1);
    const auto &lp2 = iodata2.boundaries.lumpedport.at(1);
    CHECK(lp2.R == lp1.R);
    CHECK(lp2.L == lp1.L);
    CHECK(lp2.C == lp1.C);
    CHECK(lp2.Rs == lp1.Rs);
    CHECK(lp2.Ls == lp1.Ls);
    CHECK(lp2.Cs == lp1.Cs);
    CHECK(lp2.excitation == lp1.excitation);
    CHECK(lp2.active == lp1.active);
    CHECK(iodata2.boundaries.periodic.wave_vector ==
          iodata1.boundaries.periodic.wave_vector);
    REQUIRE(iodata2.boundaries.postpro.flux.count(1) == 1);
    CHECK(iodata2.boundaries.postpro.flux.at(1).two_sided ==
          iodata1.boundaries.postpro.flux.at(1).two_sided);
    REQUIRE(iodata2.boundaries.postpro.dielectric.count(1) == 1);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).type ==
          iodata1.boundaries.postpro.dielectric.at(1).type);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).tandelta ==
          iodata1.boundaries.postpro.dielectric.at(1).tandelta);

    // Coverage gates. Each schema scope this fixture exercises is checked here so a
    // future schema addition without matching Concretize emission fails this section.
    auto model_gaps = SchemaCoverageGaps("config/model.json", "", config["Model"],
                                         /*skip=*/{"Lc"});
    INFO("Model missing keys: " << json(model_gaps).dump());
    CHECK(model_gaps.empty());

    // Boxes/Spheres are opt-in arrays of explicit refinement regions; absence means
    // no per-region refinement and there is no scalar default to emit.
    auto ref_gaps = SchemaCoverageGaps("config/model.json", "/properties/Refinement",
                                       config["Model"]["Refinement"],
                                       /*skip=*/{"Boxes", "Spheres"});
    INFO("Model.Refinement missing keys: " << json(ref_gaps).dump());
    CHECK(ref_gaps.empty());

    auto abs_gaps = SchemaCoverageGaps("config/boundaries.json", "/properties/Absorbing",
                                       config["Boundaries"]["Absorbing"]);
    INFO("Boundaries.Absorbing missing keys: " << json(abs_gaps).dump());
    CHECK(abs_gaps.empty());

    auto cond_gaps =
        SchemaCoverageGaps("config/boundaries.json", "/properties/Conductivity/items",
                           config["Boundaries"]["Conductivity"][0]);
    INFO("Boundaries.Conductivity[] missing keys: " << json(cond_gaps).dump());
    CHECK(cond_gaps.empty());

    auto imp_gaps =
        SchemaCoverageGaps("config/boundaries.json", "/properties/Impedance/items",
                           config["Boundaries"]["Impedance"][0]);
    INFO("Boundaries.Impedance[] missing keys: " << json(imp_gaps).dump());
    CHECK(imp_gaps.empty());

    // LumpedPort: Direction/CoordinateSystem/Elements are opt-in alternatives to
    // Attributes for declaring the port geometry; concretize does not synthesize them.
    auto lp_gaps =
        SchemaCoverageGaps("config/boundaries.json", "/properties/LumpedPort/items",
                           config["Boundaries"]["LumpedPort"][0],
                           /*skip=*/{"CoordinateSystem", "Elements"});
    INFO("Boundaries.LumpedPort[] missing keys: " << json(lp_gaps).dump());
    CHECK(lp_gaps.empty());

    auto per_gaps = SchemaCoverageGaps("config/boundaries.json", "/properties/Periodic",
                                       config["Boundaries"]["Periodic"]);
    INFO("Boundaries.Periodic missing keys: " << json(per_gaps).dump());
    CHECK(per_gaps.empty());
  }

  SECTION("User-written \"Default\" is replaced with the resolved concrete value")
  {
    // If the user explicitly writes the sentinel string, we must still concretize —
    // otherwise the resolved config contains a default, defeating the whole feature.
    json config = {{"Problem", {{"Type", "Eigenmode"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver",
                    {{"Eigenmode", {{"Target", 1.0}, {"Type", "Default"}}},
                     {"Linear", {{"Type", "Default"}, {"KSPType", "Default"}}}}}};

    IoData iodata(config, false);
    config = IoData::ConcretizeDefaults(iodata, config);

    auto &j_eigen = config["Solver"]["Eigenmode"];
    auto &j_linear = config["Solver"]["Linear"];
    CHECK(j_eigen["Type"].get<std::string>() != "Default");
    CHECK(j_linear["Type"].get<std::string>() != "Default");
    CHECK(j_linear["KSPType"].get<std::string>() != "Default");
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
    // The emitted Type must be the concrete backend name, not the alias "Default".
    CHECK(config["Solver"]["Eigenmode"]["Type"].get<std::string>() != "Default");
    CHECK(iodata2.solver.eigenmode.target == iodata1.solver.eigenmode.target);
    CHECK(iodata2.solver.eigenmode.target_upper == iodata1.solver.eigenmode.target_upper);
    CHECK(iodata2.solver.eigenmode.max_it == iodata1.solver.eigenmode.max_it);
    CHECK(iodata2.solver.eigenmode.max_size == iodata1.solver.eigenmode.max_size);

    // Coverage gate: every optional schema property under Solver.Eigenmode must be
    // emitted by ConcretizeDefaults. A schema addition without a Concretize update
    // surfaces here.
    auto eigen_gaps = SchemaCoverageGaps("config/solver.json", "/properties/Eigenmode",
                                         config["Solver"]["Eigenmode"]);
    INFO("Solver.Eigenmode missing keys: " << json(eigen_gaps).dump());
    CHECK(eigen_gaps.empty());
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
    CHECK(config["Solver"]["Transient"]["Type"].get<std::string>() != "Default");
    CHECK(iodata2.solver.transient.excitation == iodata1.solver.transient.excitation);
    CHECK(iodata2.solver.transient.max_t == iodata1.solver.transient.max_t);
    CHECK(iodata2.solver.transient.delta_t == iodata1.solver.transient.delta_t);

    // Coverage gate. Order/RelTol/AbsTol are intentionally only emitted for the
    // adaptive CVODE/ARKODE schemes (GEN_ALPHA/RUNGE_KUTTA warn if present); this
    // test uses the default GEN_ALPHA, so we skip them here.
    auto trans_gaps = SchemaCoverageGaps("config/solver.json", "/properties/Transient",
                                         config["Solver"]["Transient"],
                                         /*skip=*/{"Order", "RelTol", "AbsTol"});
    INFO("Solver.Transient missing keys: " << json(trans_gaps).dump());
    CHECK(trans_gaps.empty());
  }

  SECTION("Round-trip: BoundaryMode defaulted fields are reproducible")
  {
    // Freq is required at parse; only defaulted fields need verification.
    json config = {{"Problem", {{"Type", "BoundaryMode"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", {{"BoundaryMode", {{"Freq", 10.0}}}}}};

    IoData iodata1(config, false);
    config = IoData::ConcretizeDefaults(iodata1, config);

    std::string err = ValidateConfig(config);
    INFO("schema validation error: " << err);
    CHECK(err.empty());

    IoData iodata2(config, false);
    const auto &b1 = iodata1.solver.boundary_mode;
    const auto &b2 = iodata2.solver.boundary_mode;
    CHECK(b2.n == b1.n);
    CHECK(b2.n_post == b1.n_post);
    CHECK(b2.target == b1.target);
    CHECK(b2.tol == b1.tol);
    CHECK(b2.max_size == b1.max_size);
    CHECK(b2.type == b1.type);
    CHECK(config["Solver"]["BoundaryMode"]["Type"].get<std::string>() != "Default");

    // Coverage gate. Attributes is opt-in (used to extract a 2D submesh from a 3D
    // mesh); leaving it absent means "use the parent mesh as-is", so we skip it.
    auto bm_gaps = SchemaCoverageGaps("config/solver.json", "/properties/BoundaryMode",
                                      config["Solver"]["BoundaryMode"],
                                      /*skip=*/{"Attributes"});
    INFO("Solver.BoundaryMode missing keys: " << json(bm_gaps).dump());
    CHECK(bm_gaps.empty());
  }

  SECTION("Round-trip: WavePort defaulted fields are reproducible")
  {
    // Index, Attributes are required at parse; only defaulted fields need verification.
    json config = {
        {"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries", {{"WavePort", {{{"Index", 1}, {"Attributes", {2}}}}}}},
        {"Solver", {{"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.1}}}}}};

    IoData iodata1(config, false);
    config = IoData::ConcretizeDefaults(iodata1, config);

    std::string err = ValidateConfig(config);
    INFO("schema validation error: " << err);
    CHECK(err.empty());

    IoData iodata2(config, false);
    REQUIRE(iodata2.boundaries.waveport.count(1) == 1);
    const auto &w1 = iodata1.boundaries.waveport.at(1);
    const auto &w2 = iodata2.boundaries.waveport.at(1);
    CHECK(w2.mode_idx == w1.mode_idx);
    CHECK(w2.d_offset == w1.d_offset);
    CHECK(w2.eigen_solver == w1.eigen_solver);
    CHECK(config["Boundaries"]["WavePort"][0]["SolverType"].get<std::string>() !=
          "Default");
    CHECK(w2.active == w1.active);
    CHECK(w2.ksp_max_its == w1.ksp_max_its);
    CHECK(w2.ksp_tol == w1.ksp_tol);
    CHECK(w2.eig_tol == w1.eig_tol);
    CHECK(w2.max_size == w1.max_size);
    CHECK(w2.verbose == w1.verbose);
    CHECK(w2.n_samples == w1.n_samples);

    // Coverage gate. VoltagePath is opt-in coordinate path for line integral
    // postprocessing on the port face; absence means no voltage line integral.
    auto wp_gaps =
        SchemaCoverageGaps("config/boundaries.json", "/properties/WavePort/items",
                           config["Boundaries"]["WavePort"][0],
                           /*skip=*/{"VoltagePath"});
    INFO("Boundaries.WavePort[] missing keys: " << json(wp_gaps).dump());
    CHECK(wp_gaps.empty());
  }

  SECTION("Round-trip: Material defaulted physical properties are reproducible")
  {
    // Attributes is required at parse; physical properties (mu_r, epsilon_r, tandelta,
    // sigma, lambda_L) all have scalar defaults that must survive concretize → reparse.
    json config = {
        {"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries", json::object()},
        {"Solver", {{"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.1}}}}}};

    IoData iodata1(config, false);
    config = IoData::ConcretizeDefaults(iodata1, config);

    std::string err = ValidateConfig(config);
    INFO("schema validation error: " << err);
    CHECK(err.empty());

    // The concretized JSON must mention every defaulted physical property.
    auto &j_mat = config["Domains"]["Materials"][0];
    CHECK(j_mat.contains("Permeability"));
    CHECK(j_mat.contains("Permittivity"));
    CHECK(j_mat.contains("LossTan"));
    CHECK(j_mat.contains("Conductivity"));
    CHECK(j_mat.contains("LondonDepth"));

    IoData iodata2(config, false);
    REQUIRE(iodata1.domains.materials.size() == 1);
    REQUIRE(iodata2.domains.materials.size() == 1);
    const auto &m1 = iodata1.domains.materials[0];
    const auto &m2 = iodata2.domains.materials[0];
    for (int k = 0; k < 3; ++k)
    {
      CHECK(m2.mu_r.s[k] == m1.mu_r.s[k]);
      CHECK(m2.epsilon_r.s[k] == m1.epsilon_r.s[k]);
      CHECK(m2.tandelta.s[k] == m1.tandelta.s[k]);
      CHECK(m2.sigma.s[k] == m1.sigma.s[k]);
    }
    CHECK(m2.lambda_L == m1.lambda_L);

    // Coverage gate. MaterialAxes is opt-in: omission means "diagonal in standard
    // basis". Concretize does not synthesize one.
    auto mat_gaps = SchemaCoverageGaps("config/domains.json", "/properties/Materials/items",
                                       config["Domains"]["Materials"][0],
                                       /*skip=*/{"MaterialAxes"});
    INFO("Domains.Materials[] missing keys: " << json(mat_gaps).dump());
    CHECK(mat_gaps.empty());
  }
}
