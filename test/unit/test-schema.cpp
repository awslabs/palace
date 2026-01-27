// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <filesystem>
#include <fstream>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "utils/iodata.hpp"

using json = nlohmann::json;
using namespace palace;

namespace fs = std::filesystem;

TEST_CASE("Schema Validation - Example Configs", "[schema][Serial]")
{
  // Schema directory is relative to test source directory.
  std::string schema_dir = fmt::format("{}/../../scripts/schema", PALACE_TEST_DIR);
  std::string examples_dir = fmt::format("{}/../../examples", PALACE_TEST_DIR);

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
      std::string err = ValidateConfig(config, schema_dir);
      INFO("Config: " << config_file);
      INFO("Error: " << err);
      CHECK(err.empty());
    }
  }
}

TEST_CASE("Schema Validation - Invalid Config", "[schema][Serial]")
{
  std::string schema_dir = fmt::format("{}/../../scripts/schema", PALACE_TEST_DIR);

  SECTION("Missing required field")
  {
    json config = {{"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", {}},
                   {"Solver", {}}};
    // Missing "Problem" which is required.

    std::string err = ValidateConfig(config, schema_dir);
    CHECK(!err.empty());
  }

  SECTION("Invalid enum value")
  {
    json config = {{"Problem", {{"Type", "InvalidType"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", {}},
                   {"Solver", {}}};

    std::string err = ValidateConfig(config, schema_dir);
    CHECK(!err.empty());
  }

  SECTION("Additional property not allowed")
  {
    json config = {{"Problem", {{"Type", "Eigenmode"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", {}},
                   {"Solver", {}},
                   {"UnknownSection", {}}};

    std::string err = ValidateConfig(config, schema_dir);
    CHECK(!err.empty());
  }
}
