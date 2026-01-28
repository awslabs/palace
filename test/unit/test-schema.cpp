// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <filesystem>
#include <fstream>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "embedded_schema.hpp"
#include "fixtures.hpp"
#include "utils/iodata.hpp"
#include "utils/jsonschema.hpp"

using json = nlohmann::json;
using namespace palace;

namespace fs = std::filesystem;

TEST_CASE("Schema Validation - Embedded Schema Matches Source", "[schema][Serial]")
{
  // Verify embedded schemas match source files (catches stale builds). This
  // test only makes sense when PALACE_TEST_DATA_DIR is in the folder where
  // Palace is being developed. E.g., if Palace is being built with Cmake in a
  // palace_repo/build type of folder.
  std::string schema_dir = fmt::format("{}/../../scripts/schema", PALACE_TEST_DATA_DIR);

  if (std::filesystem::exists(schema_dir) && std::filesystem::is_directory(schema_dir))
  {
    for (const auto &[path, embedded_content] : schema::GetSchemaMap())
    {
      SECTION(path)
      {
        std::string full_path = schema_dir + "/" + path;
        std::ifstream f(full_path);
        if (!f.is_open())
        {
          SKIP("Schema source not found (installed build?): " << full_path);
        }
        // Parse both to compare (ignores whitespace differences).
        json embedded = json::parse(embedded_content);
        json source = json::parse(f);
        INFO("Schema file: " << full_path);
        CHECK(embedded == source);
      }
    }
  }
  else
  {
    SKIP("Schema source not found (installed build?): " << schema_dir);
  }
}


TEST_CASE("Schema Validation - Example Configs", "[schema][Serial]")
{
  // Schema directory is relative to test source directory.
  std::string examples_dir = fmt::format("{}/examples", PALACE_TEST_DATA_DIR);

  // Collect JSON config files directly in example subdirectories (not in postpro/output).
  std::vector<std::string> config_files;
  for (const auto &entry : fs::directory_iterator(examples_dir))
  {
    if (entry.is_directory())
    {
      for (const auto &file : fs::directory_iterator(entry.path()))
      {
        if (file.path().extension() == ".json")
        {
          config_files.push_back(file.path().string());
        }
      }
    }
  }

  REQUIRE(!config_files.empty());

  for (const auto &config_file : config_files)
  {
    SECTION(fs::path(config_file).filename().string())
    {
      // Preprocess and parse the config file.
      std::stringstream buffer = PreprocessFile(config_file.c_str());
      json config;
      REQUIRE_NOTHROW(config = json::parse(buffer));

      // Validate against schema.
      std::string err = ValidateConfig(config);
      INFO("Config: " << config_file);
      INFO("Error: " << err);
      CHECK(err.empty());
    }
  }
}

TEST_CASE_METHOD(palace::test::PerRankTempDir, "Schema Validation - Config with Comments",
                 "[schema][Serial]")
{
  // Test that preprocessing (comment stripping) works with schema validation.
  auto temp_path = temp_dir / "palace_test_comments.json";
  {
    std::ofstream f(temp_path);
    f << R"({
      // C++ style comment
      "Problem": { "Type": "Eigenmode" },
      /* C style comment */
      "Model": { "Mesh": "test.msh" },
      "Domains": { "Materials": [{ "Attributes": [1] }] },
      "Boundaries": {},
      "Solver": { "Eigenmode": { "Target": 1.0 } }
    })";
  }

  std::stringstream buffer = PreprocessFile(temp_path.c_str());
  json config;
  REQUIRE_NOTHROW(config = json::parse(buffer));

  std::string err = ValidateConfig(config);
  CHECK(err.empty());
}

TEST_CASE("Schema Validation - Invalid Config", "[schema][Serial]")
{

  SECTION("Missing required field")
  {
    json config = {{"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", json::object()}};
    // Missing "Problem" which is required.

    std::string err = ValidateConfig(config);
    CHECK(!err.empty());
  }

  SECTION("Invalid enum value")
  {
    json config = {{"Problem", {{"Type", "InvalidType"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", json::object()}};

    std::string err = ValidateConfig(config);
    CHECK(!err.empty());
  }

  SECTION("Additional property not allowed")
  {
    json config = {{"Problem", {{"Type", "Eigenmode"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", json::object()},
                   {"UnknownSection", {}}};

    std::string err = ValidateConfig(config);
    CHECK(!err.empty());
  }

  SECTION("Problem.Type requires matching Solver section")
  {
    // Driven type requires Solver.Driven section
    json driven_missing = {{"Problem", {{"Type", "Driven"}}},
                           {"Model", {{"Mesh", "test.msh"}}},
                           {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                           {"Boundaries", json::object()},
                           {"Solver", {{"Linear", {}}}}};
    CHECK(!ValidateConfig(driven_missing).empty());

    // Eigenmode type requires Solver.Eigenmode section
    json eigen_missing = {{"Problem", {{"Type", "Eigenmode"}}},
                          {"Model", {{"Mesh", "test.msh"}}},
                          {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                          {"Boundaries", json::object()},
                          {"Solver", {{"Linear", {}}}}};
    CHECK(!ValidateConfig(eigen_missing).empty());

    // Transient type requires Solver.Transient section
    json transient_missing = {{"Problem", {{"Type", "Transient"}}},
                              {"Model", {{"Mesh", "test.msh"}}},
                              {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                              {"Boundaries", json::object()},
                              {"Solver", {{"Linear", {}}}}};
    CHECK(!ValidateConfig(transient_missing).empty());

    // Electrostatic and Magnetostatic don't require matching sections (have defaults)
    json electro_ok = {{"Problem", {{"Type", "Electrostatic"}}},
                       {"Model", {{"Mesh", "test.msh"}}},
                       {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                       {"Boundaries", json::object()},
                       {"Solver", json::object()}};
    CHECK(ValidateConfig(electro_ok).empty());

    json magneto_ok = {{"Problem", {{"Type", "Magnetostatic"}}},
                       {"Model", {{"Mesh", "test.msh"}}},
                       {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                       {"Boundaries", json::object()},
                       {"Solver", json::object()}};
    CHECK(ValidateConfig(magneto_ok).empty());
  }
}

TEST_CASE("Schema Validation - Sub-schema by Key", "[schema][Serial]")
{

  SECTION("Valid LumpedPort")
  {
    json port = {{"Index", 1}, {"Attributes", {1, 2}}, {"R", 50.0}};
    std::string err = ValidateConfig(port, "LumpedPort");
    INFO("Error: " << err);
    CHECK(err.empty());
  }

  SECTION("Valid LumpedPort with optional fields")
  {
    json port = {{"Index", 1}, {"Attributes", {1}}, {"R", 50.0},
                 {"L", 1e-9},  {"C", 1e-12},        {"Excitation", true}};
    std::string err = ValidateConfig(port, "LumpedPort");
    INFO("Error: " << err);
    CHECK(err.empty());
  }

  SECTION("Invalid LumpedPort - negative Index")
  {
    json port = {{"Index", -1}, {"Attributes", {1}}};
    std::string err = ValidateConfig(port, "LumpedPort");
    CHECK(!err.empty());
  }

  SECTION("Invalid LumpedPort - negative Index")
  {
    json port = {{"Index", -1}, {"Attributes", {1}}, {"R", 50.0}, {"Direction", "+Y"}};
    std::string err = ValidateConfig(port, "LumpedPort");
    CHECK(!err.empty());
  }

  SECTION("Invalid WavePort - negative Index")
  {
    json port = {{"Index", -1}, {"Attributes", {1}}};
    std::string err = ValidateConfig(port, "WavePort");
    CHECK(!err.empty());
  }

  SECTION("Invalid Direction strings")
  {
    std::vector<std::string> invalid_dirs = {"a",  "+a", "-a",  "xx", "~x",
                                             "x+", "xy", "xyz", "abc"};
    for (const auto &dir : invalid_dirs)
    {
      json port = {{"Index", 1}, {"Attributes", {1}}, {"Direction", dir}};
      std::string err = ValidateConfig(port, "LumpedPort");
      INFO("Direction: " << dir);
      CHECK(!err.empty());
    }
  }

  SECTION("Valid Material")
  {
    json mat = {{"Attributes", {1}}, {"Permittivity", 2.0}};
    std::string err = ValidateConfig(mat, "Materials");
    INFO("Error: " << err);
    CHECK(err.empty());
  }

  SECTION("Invalid Material - bad Permittivity type")
  {
    json mat = {{"Attributes", {1}}, {"Permittivity", "not a number"}};
    std::string err = ValidateConfig(mat, "Materials");
    CHECK(!err.empty());
  }

  SECTION("Valid WavePort")
  {
    json port = {{"Index", 1}, {"Attributes", {1}}, {"Mode", 2}};
    std::string err = ValidateConfig(port, "WavePort");
    INFO("Error: " << err);
    CHECK(err.empty());
  }

  SECTION("Wrong key - LumpedPort validated as WavePort")
  {
    // LumpedPort has R/L/C which WavePort doesn't allow.
    json port = {{"Index", 1}, {"Attributes", {1}}, {"R", 50.0}};
    std::string err = ValidateConfig(port, "WavePort");
    CHECK(!err.empty());
  }

  SECTION("Wrong key - WavePort validated as LumpedPort")
  {
    // WavePort has Mode which LumpedPort doesn't have.
    json port = {{"Index", 1}, {"Attributes", {1}}, {"Mode", 2}};
    std::string err = ValidateConfig(port, "LumpedPort");
    CHECK(!err.empty());
  }

  SECTION("Nonexistent schema key")
  {
    json data = {{"Index", 1}};
    std::string err = ValidateConfig(data, "NonexistentKey");
    CHECK(!err.empty());
    CHECK(err.find("not found") != std::string::npos);
  }

  SECTION("Ambiguous schema key")
  {
    // "Postprocessing" exists in both boundaries.json and domains.json.
    json data = {{"Index", 1}};
    std::string err = ValidateConfig(data, "Postprocessing");
    CHECK(!err.empty());
    CHECK(err.find("not found") != std::string::npos);
  }
}

TEST_CASE("Schema Validation - Array Type Checks", "[schema][Serial]")
{

  SECTION("FloquetWaveVector must be array")
  {
    // Valid: array
    json periodic_valid = {{"FloquetWaveVector", {1.0, 0.0, 0.0}},
                           {"BoundaryPairs",
                            {{{"DonorAttributes", {1}},
                              {"ReceiverAttributes", {2}},
                              {"Translation", {1, 0, 0}}}}}};
    std::string err = ValidateConfig(periodic_valid, "Periodic");
    INFO("Error: " << err);
    CHECK(err.empty());

    // Invalid: not an array
    json periodic_invalid = {{"FloquetWaveVector", "not an array"},
                             {"BoundaryPairs",
                              {{{"DonorAttributes", {1}},
                                {"ReceiverAttributes", {2}},
                                {"Translation", {1, 0, 0}}}}}};
    err = ValidateConfig(periodic_invalid, "Periodic");
    CHECK(!err.empty());
  }

  SECTION("BoundaryPairs must be array")
  {
    // Invalid: not an array
    json periodic_invalid = {{"BoundaryPairs", "not an array"}};
    std::string err = ValidateConfig(periodic_invalid, "Periodic");
    CHECK(!err.empty());
  }

  SECTION("Translation must be array")
  {
    // Valid: array
    json periodic_valid = {{"BoundaryPairs",
                            {{{"DonorAttributes", {1}},
                              {"ReceiverAttributes", {2}},
                              {"Translation", {1, 0, 0}}}}}};
    std::string err = ValidateConfig(periodic_valid, "Periodic");
    INFO("Error: " << err);
    CHECK(err.empty());

    // Invalid: not an array
    json periodic_invalid = {{"BoundaryPairs",
                              {{{"DonorAttributes", {1}},
                                {"ReceiverAttributes", {2}},
                                {"Translation", "not array"}}}}};
    err = ValidateConfig(periodic_invalid, "Periodic");
    CHECK(!err.empty());
  }

  SECTION("AffineTransformation must be array")
  {
    // Valid: 16-element array (4x4 matrix)
    json periodic_valid = {
        {"BoundaryPairs",
         {{{"DonorAttributes", {1}},
           {"ReceiverAttributes", {2}},
           {"AffineTransformation", {1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}}}}}};
    std::string err = ValidateConfig(periodic_valid, "Periodic");
    INFO("Error: " << err);
    CHECK(err.empty());

    // Invalid: not an array
    json periodic_invalid = {{"BoundaryPairs",
                              {{{"DonorAttributes", {1}},
                                {"ReceiverAttributes", {2}},
                                {"AffineTransformation", "not array"}}}}};
    err = ValidateConfig(periodic_invalid, "Periodic");
    CHECK(!err.empty());
  }

  SECTION("ThetaPhis must be array of arrays")
  {
    // Valid: array of [theta, phi] pairs
    json farfield_valid = {{"Attributes", {1}}, {"ThetaPhis", {{0.0, 0.0}, {90.0, 45.0}}}};
    std::string err = ValidateConfig(farfield_valid, "FarField");
    INFO("Error: " << err);
    CHECK(err.empty());

    // Invalid: not an array
    json farfield_invalid = {{"Attributes", {1}}, {"ThetaPhis", "not array"}};
    err = ValidateConfig(farfield_invalid, "FarField");
    CHECK(!err.empty());

    // Invalid: array but inner elements not arrays
    json farfield_invalid2 = {{"Attributes", {1}}, {"ThetaPhis", {1.0, 2.0}}};
    err = ValidateConfig(farfield_invalid2, "FarField");
    CHECK(!err.empty());
  }
}

TEST_CASE_METHOD(palace::test::PerRankTempDir, "Schema Validation - Range Expansion",
                 "[schema][Serial]")
{
  // Test that integer range syntax (e.g., 1-5) is expanded before validation.
  auto temp_path = temp_dir / "palace_test_range.json";
  {
    std::ofstream f(temp_path);
    f << R"({
      "Problem": { "Type": "Eigenmode" },
      "Model": { "Mesh": "test.msh" },
      "Domains": { "Materials": [{ "Attributes": [1-3, 5, 7-9] }] },
      "Boundaries": {},
      "Solver": { "Eigenmode": { "Target": 1.0 } }
    })";
  }

  std::stringstream buffer = PreprocessFile(temp_path.c_str());
  json config = json::parse(buffer);

  // Verify range expansion worked.
  auto attrs = config["Domains"]["Materials"][0]["Attributes"];
  CHECK(attrs == json({1, 2, 3, 5, 7, 8, 9}));

  // Verify schema validation passes.
  std::string err = ValidateConfig(config);
  CHECK(err.empty());
}
