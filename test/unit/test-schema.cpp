// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <filesystem>
#include <fstream>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "utils/iodata.hpp"
#include "utils/jsonschema.hpp"

using json = nlohmann::json;
using namespace palace;

namespace fs = std::filesystem;

TEST_CASE("Schema Validation - Example Configs", "[schema][Serial]")
{
  // Schema directory is relative to test source directory.
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
      std::string err = ValidateConfig(config);
      INFO("Config: " << config_file);
      INFO("Error: " << err);
      CHECK(err.empty());
    }
  }
}

TEST_CASE("Schema Validation - Invalid Config", "[schema][Serial]")
{

  SECTION("Missing required field")
  {
    json config = {{"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", {}},
                   {"Solver", {}}};
    // Missing "Problem" which is required.

    std::string err = ValidateConfig(config);
    CHECK(!err.empty());
  }

  SECTION("Invalid enum value")
  {
    json config = {{"Problem", {{"Type", "InvalidType"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", {}},
                   {"Solver", {}}};

    std::string err = ValidateConfig(config);
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

    std::string err = ValidateConfig(config);
    CHECK(!err.empty());
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

TEST_CASE("Schema Validation - Required Field Checks", "[schema][Serial]")
{

  SECTION("Periodic BoundaryPairs requires DonorAttributes and ReceiverAttributes")
  {
    // Valid: both present
    json periodic_valid = {{"BoundaryPairs",
                            {{{"DonorAttributes", {1}},
                              {"ReceiverAttributes", {2}},
                              {"Translation", {1, 0, 0}}}}}};
    std::string err = ValidateConfig(periodic_valid, "Periodic");
    INFO("Error: " << err);
    CHECK(err.empty());

    // Invalid: missing DonorAttributes
    json periodic_no_donor = {
        {"BoundaryPairs", {{{"ReceiverAttributes", {2}}, {"Translation", {1, 0, 0}}}}}};
    err = ValidateConfig(periodic_no_donor, "Periodic");
    CHECK(!err.empty());

    // Invalid: missing ReceiverAttributes
    json periodic_no_receiver = {
        {"BoundaryPairs", {{{"DonorAttributes", {1}}, {"Translation", {1, 0, 0}}}}}};
    err = ValidateConfig(periodic_no_receiver, "Periodic");
    CHECK(!err.empty());
  }

  SECTION("LumpedPort requires either Attributes or Elements")
  {
    // Valid: with Attributes
    json port_attrs = {{"Index", 1}, {"Attributes", {1}}, {"Direction", "+X"}};
    std::string err = ValidateConfig(port_attrs, "LumpedPort");
    INFO("Error: " << err);
    CHECK(err.empty());

    // Valid: with Elements
    json port_elems = {{"Index", 1},
                       {"Elements", {{{"Attributes", {1}}, {"Direction", "+X"}}}}};
    err = ValidateConfig(port_elems, "LumpedPort");
    INFO("Error: " << err);
    CHECK(err.empty());

    // Invalid: neither Attributes nor Elements
    json port_neither = {{"Index", 1}, {"R", 50.0}};
    err = ValidateConfig(port_neither, "LumpedPort");
    CHECK(!err.empty());
  }

  SECTION("SurfaceCurrent requires either Attributes or Elements")
  {
    // Valid: with Attributes
    json current_attrs = {{"Index", 1}, {"Attributes", {1}}, {"Direction", "+X"}};
    std::string err = ValidateConfig(current_attrs, "SurfaceCurrent");
    INFO("Error: " << err);
    CHECK(err.empty());

    // Valid: with Elements
    json current_elems = {{"Index", 1},
                          {"Elements", {{{"Attributes", {1}}, {"Direction", "+X"}}}}};
    err = ValidateConfig(current_elems, "SurfaceCurrent");
    INFO("Error: " << err);
    CHECK(err.empty());

    // Invalid: neither Attributes nor Elements
    json current_neither = {{"Index", 1}};
    err = ValidateConfig(current_neither, "SurfaceCurrent");
    CHECK(!err.empty());
  }
}

TEST_CASE("Schema Validation - Mutual Exclusion", "[schema][Serial]")
{

  SECTION("PEC and Ground are mutually exclusive")
  {
    // Valid: only PEC
    json boundaries_pec = {{"PEC", {{"Attributes", {1}}}}};
    std::string err = ValidateConfig(boundaries_pec, "Boundaries");
    CHECK(err.empty());

    // Valid: only Ground
    json boundaries_ground = {{"Ground", {{"Attributes", {1}}}}};
    err = ValidateConfig(boundaries_ground, "Boundaries");
    CHECK(err.empty());

    // Invalid: both PEC and Ground
    json boundaries_both = {{"PEC", {{"Attributes", {1}}}},
                            {"Ground", {{"Attributes", {2}}}}};
    err = ValidateConfig(boundaries_both, "Boundaries");
    CHECK(!err.empty());
  }

  SECTION("PMC and ZeroCharge are mutually exclusive")
  {
    // Valid: only PMC
    json boundaries_pmc = {{"PMC", {{"Attributes", {1}}}}};
    std::string err = ValidateConfig(boundaries_pmc, "Boundaries");
    CHECK(err.empty());

    // Valid: only ZeroCharge
    json boundaries_zeroq = {{"ZeroCharge", {{"Attributes", {1}}}}};
    err = ValidateConfig(boundaries_zeroq, "Boundaries");
    CHECK(err.empty());

    // Invalid: both PMC and ZeroCharge
    json boundaries_both = {{"PMC", {{"Attributes", {1}}}},
                            {"ZeroCharge", {{"Attributes", {2}}}}};
    err = ValidateConfig(boundaries_both, "Boundaries");
    CHECK(!err.empty());
  }
}

TEST_CASE("Schema Validation - Error Message Format", "[schema][Serial]")
{

  SECTION("Invalid enum value shows valid options")
  {
    json config = {{"Problem", {{"Type", "InvalidType"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", json::object()}};

    std::string err = ValidateConfig(config);
    CHECK(
        err ==
        "At [\"Problem\"][\"Type\"]: instance not found in required enum; valid values: "
        "\"Eigenmode\", \"Driven\", \"Transient\", \"Electrostatic\", \"Magnetostatic\"\n");
  }

  SECTION("Invalid enum in nested array")
  {
    json config = {
        {"Problem", {{"Type", "Driven"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries",
         {{"LumpedPort", {{{"Index", 1}, {"Attributes", {1}}, {"Direction", "BadDir"}}}}}},
        {"Solver", {{"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.1}}}}}};

    std::string err = ValidateConfig(config);
    // Direction uses anyOf (string enum or array), so error shows subschema failures.
    CHECK(err.find("[\"Boundaries\"][\"LumpedPort\"][0][\"Direction\"]") !=
          std::string::npos);
    CHECK(err.find("anyOf") != std::string::npos);
  }

  SECTION("Wrong type shows actual type")
  {
    json port = {{"Index", "not a number"}, {"Attributes", {1}}};
    std::string err = ValidateConfig(port, "LumpedPort");
    CHECK(err == "At [\"Index\"]: unexpected instance type (got string)\n");
  }

  SECTION("Value below minimum")
  {
    json port = {{"Index", -1}, {"Attributes", {1}}};
    std::string err = ValidateConfig(port, "LumpedPort");
    CHECK(err == "At [\"Index\"]: instance is below or equals minimum of 0\n");
  }

  SECTION("Missing required field shows oneOf options")
  {
    json config = {
        {"Problem", {{"Type", "Driven"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries", {{"LumpedPort", {{{"Index", 1}, {"R", 50.0}}}}}},
        {"Solver", {{"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.1}}}}}};

    std::string err = ValidateConfig(config);
    CHECK(err ==
          "At [\"Boundaries\"][\"LumpedPort\"][0]: no subschema has succeeded, but one of "
          "them is required to validate. Type: oneOf, number of failed subschemas: 2\n"
          "At [\"Boundaries\"][\"LumpedPort\"][0]: [combination: oneOf / case#0] required "
          "property 'Attributes' not found in object\n"
          "At [\"Boundaries\"][\"LumpedPort\"][0]: [combination: oneOf / case#1] required "
          "property 'Elements' not found in object\n");
  }

  SECTION("Additional property not allowed")
  {
    json config = {{"Problem", {{"Type", "Eigenmode"}, {"UnknownField", 123}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver", {{"Eigenmode", {{"Target", 1.0}}}}}};

    std::string err = ValidateConfig(config);
    CHECK(err.find("[\"Problem\"]") != std::string::npos);
    CHECK(err.find("UnknownField") != std::string::npos);
  }
}
