// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <functional>
#include <iterator>
#include <limits>
#include <numbers>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <fmt/format.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
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
// `pointer` is a JSON Pointer into the root schema document selecting the scope to walk
// (e.g. "/properties/Model" for the Model section,
// "/properties/Model/properties/Refinement" for a nested object, or
// "/properties/Boundaries/properties/Conductivity/items" for an array element).
//
// Use this in round-trip tests to catch the case where someone adds a new optional
// schema property without wiring `ConcretizeDefaults` to emit it.
std::vector<std::string> SchemaCoverageGaps(const std::string &pointer,
                                            const json &concretized_subtree,
                                            const std::set<std::string> &skip = {})
{
  const auto &schema_map = schema::GetSchemaMap();
  auto it = schema_map.find("config-schema.json");
  REQUIRE(it != schema_map.end());
  const json schema = json::parse(it->second);

  // Resolve an in-document `#/$defs/...` reference against the root schema. After the
  // hoist, named shapes live under `$defs` and are reached via `$ref`, so a raw JSON
  // Pointer lands on the ref node rather than the inline object. Mirrors the production
  // resolver in jsonschema.cpp.
  auto resolve_ref = [&schema](const json &node) -> const json &
  {
    auto rit = node.find("$ref");
    if (rit == node.end())
    {
      return node;
    }
    const std::string ref = rit->get<std::string>();
    const std::string defs_prefix = "#/$defs/";
    REQUIRE(ref.starts_with(defs_prefix));
    const std::string def_name = ref.substr(defs_prefix.size());
    REQUIRE(schema.contains("$defs"));
    REQUIRE(schema.at("$defs").contains(def_name));
    return schema.at("$defs").at(def_name);
  };

  // Gather the arms to inspect. A `oneOf`/`anyOf` of `$ref`s (e.g. LumpedPort, split by
  // the hoist into Attributes/Elements variants) is treated as the union of its arms'
  // properties, with `required` the *intersection* across arms — a field required in
  // only one arm is optional overall, matching the pre-hoist single-object schema.
  const json &resolved = resolve_ref(schema.at(json::json_pointer(pointer)));
  std::vector<json> arms;
  if (auto cit = resolved.find("oneOf"); cit != resolved.end())
  {
    for (const auto &arm : *cit)
    {
      arms.push_back(resolve_ref(arm));
    }
  }
  else if (auto ait = resolved.find("anyOf"); ait != resolved.end())
  {
    for (const auto &arm : *ait)
    {
      arms.push_back(resolve_ref(arm));
    }
  }
  else
  {
    arms.push_back(resolved);
  }
  REQUIRE(!arms.empty());

  std::set<std::string> property_names;
  std::set<std::string> required_set;
  for (std::size_t i = 0; i < arms.size(); ++i)
  {
    const json &arm = arms[i];
    REQUIRE(arm.contains("properties"));
    for (const auto &[name, _sub] : arm.at("properties").items())
    {
      property_names.insert(name);
    }
    std::set<std::string> arm_required;
    if (auto rit = arm.find("required"); rit != arm.end())
    {
      for (const auto &r : *rit)
      {
        arm_required.insert(r.get<std::string>());
      }
    }
    if (i == 0)
    {
      required_set = std::move(arm_required);
    }
    else
    {
      std::set<std::string> intersection;
      std::set_intersection(required_set.begin(), required_set.end(), arm_required.begin(),
                            arm_required.end(),
                            std::inserter(intersection, intersection.begin()));
      required_set = std::move(intersection);
    }
  }

  std::vector<std::string> gaps;
  for (const auto &name : property_names)
  {
    if (required_set.contains(name) || skip.contains(name))
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

TEST_CASE("Config Domain Postprocessing", "[config][Serial]")
{
  auto MakeDomains = [](const std::vector<int> &material_attributes,
                        const std::vector<int> &postprocessing_attributes)
  {
    return json{
        {"Materials", {{{"Attributes", material_attributes}}}},
        {"Postprocessing",
         {{"Energy", {{{"Index", 1}, {"Attributes", postprocessing_attributes}}}}}}};
  };

  SECTION("Postprocessing material domain is accepted")
  {
    CHECK_NOTHROW(config::DomainData(MakeDomains({1, 3}, {1, 3})));
  }

  SECTION("Postprocessing domain between material attributes is rejected")
  {
    // This specifically exercises a missing value whose lower bound is not end().
    CHECK_THROWS_WITH(config::DomainData(MakeDomains({1, 3}, {2})),
                      Catch::Matchers::ContainsSubstring(
                          "Domain postprocessing attribute 2 has no corresponding entry in "
                          "config[\"Domains\"][\"Materials\"]!"));
  }

  SECTION("Postprocessing domain above material attributes is rejected")
  {
    CHECK_THROWS_WITH(config::DomainData(MakeDomains({1, 3}, {4})),
                      Catch::Matchers::ContainsSubstring(
                          "Domain postprocessing attribute 4 has no corresponding entry in "
                          "config[\"Domains\"][\"Materials\"]!"));
  }
}

TEST_CASE("Config Interior Boundary Sheets", "[config][Serial]")
{
  // Sheets (for error estimation and interface dielectric postprocessing) and the
  // attributes along which the mesh is cracked, by default and with the deprecated global
  // option.
  json boundaries = {
      {"PEC", {{"Attributes", {1}}}},
      {"PMC", {{"Attributes", {2}}}},
      {"Absorbing", {{"Attributes", {3}}}},
      {"Impedance",
       {{{"Attributes", {4}}, {"Ls", 1.0e-12}},
        {{"Attributes", {5}}, {"Ls", 1.0e-12}, {"Crack", true}}}},
      {"Conductivity",
       {{{"Attributes", {6}}, {"Conductivity", 1.0e7}},
        {{"Attributes", {7}}, {"Conductivity", 1.0e7}, {"Crack", true}}}},
      {"RationalImpedance",
       {{{"Attributes", {8}}, {"Numerator", {1.0}}, {"Denominator", {1.0}}}}},
      {"LumpedPort",
       {{{"Index", 1}, {"Attributes", {9}}, {"R", 50.0}, {"Direction", "+X"}}}},
      {"SurfaceCurrent", {{{"Index", 2}, {"Attributes", {10}}, {"Direction", "+X"}}}},
      {"WavePortPEC", {{"Attributes", {11}}}}};
  config::BoundaryData data(boundaries);
  CHECK_FALSE(data.impedance[0].crack);
  CHECK(data.impedance[1].crack);
  CHECK_FALSE(data.conductivity[0].crack);
  CHECK(data.conductivity[1].crack);
  // Whether "Crack" is given is recorded, to tell an explicit false from the default.
  CHECK_FALSE(data.conductivity[0].crack_given);
  CHECK(data.conductivity[1].crack_given);
  const json explicit_false = {
      {"Attributes", {15}}, {"Conductivity", 1.0e7}, {"Crack", false}};
  const config::ConductivityData cond_false(explicit_false);
  CHECK(cond_false.crack_given);
  CHECK_FALSE(cond_false.crack);
  CHECK_FALSE(data.rational_impedance[0].crack);
  CHECK(data.GetSheetAttributes() == std::vector<int>{1, 4, 5, 6, 7, 8, 9, 10});
  // The metal sheets: not the lumped port (9) or the surface current source (10).
  CHECK(data.GetMetalSheetAttributes() == std::vector<int>{1, 4, 5, 6, 7, 8});
  // All impedance sheets are sheets, including purely capacitive ones.
  boundaries["Impedance"][0].erase("Ls");
  boundaries["Impedance"][0]["Cs"] = 1.0e-15;
  config::BoundaryData data_cs(boundaries);
  CHECK(data_cs.GetSheetAttributes() == std::vector<int>{1, 4, 5, 6, 7, 8, 9, 10});
  CHECK(data.GetMeshCrackAttributes(std::nullopt) == std::vector<int>{2, 3, 5, 7});
  CHECK(data.GetMeshCrackAttributes(false).empty());
  CHECK(data.GetMeshCrackAttributes(true) ==
        std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 10, 11});

  // Superconductor sheets and flux loop films carry surface currents, but are never
  // cracked (thin films), and flux loop holes are not sheets.
  json sc_boundaries = {
      {"PEC", {{"Attributes", {1}}}},
      {"Superconductor", {{{"Attributes", {12}}, {"KineticInductance", 1.0e-12}}}},
      {"FluxLoop",
       {{{"Index", 1},
         {"FilmAttributes", {13}},
         {"HoleAttributes", {14}},
         {"FluxAmounts", {1.0}}}}}};
  config::BoundaryData sc_data(sc_boundaries);
  CHECK(sc_data.GetSheetAttributes() == std::vector<int>{1, 12, 13});
  CHECK(sc_data.GetMetalSheetAttributes() == std::vector<int>{1, 12, 13});
  CHECK(sc_data.GetMeshCrackAttributes(std::nullopt).empty());
  CHECK(sc_data.GetMeshCrackAttributes(true) == std::vector<int>{1});

  // The option is not available for boundaries which are never cracked.
  json config = {{"Problem", {{"Type", "Eigenmode"}}},
                 {"Model", {{"Mesh", "mesh.msh"}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                 {"Boundaries", {{"PEC", {{"Attributes", {1}}}}}},
                 {"Solver", {{"Eigenmode", {{"Target", 1.0}}}}}};
  std::string err = ValidateConfig(config);
  INFO("schema validation error: " << err);
  CHECK(err.empty());
  config["Boundaries"]["PEC"]["Crack"] = true;
  CHECK_FALSE(ValidateConfig(config).empty());
  config["Boundaries"]["PEC"].erase("Crack");
  config["Problem"]["Type"] = "Magnetostatic";
  config["Boundaries"]["Superconductor"] = {
      {{"Attributes", {2}}, {"KineticInductance", 1.0e-12}}};
  CHECK(ValidateConfig(config).empty());
  config["Boundaries"]["Superconductor"][0]["Crack"] = true;
  CHECK_FALSE(ValidateConfig(config).empty());
}

TEST_CASE("Config metal sheet attributes are the one list of metal sheets",
          "[config][Serial]")
{
  // Review r2198 m-4: BoundaryData::GetMetalSheetAttributes is the one list of the
  // boundaries which model a metal sheet (the SA vacuum-sheet warning reads it), and
  // GetSheetAttributes is that list with the lumped ports and surface current sources. One
  // boundary of every kind (the parse accepts Terminal with PrescribedPotential; only the
  // validation refuses the pair).
  json boundaries = {
      {"PEC", {{"Attributes", {1}}}},
      {"PMC", {{"Attributes", {2}}}},
      {"Terminal", {{{"Index", 1}, {"Attributes", {3}}}}},
      {"PrescribedPotential",
       {{{"Index", 2},
         {"Attributes", {4}},
         {"DataFile", "trace.csv"},
         {"TerminalAttributes", {5}}}}},
      {"Impedance", {{{"Attributes", {6}}, {"Ls", 1.0e-12}, {"Crack", true}}}},
      {"RationalImpedance",
       {{{"Attributes", {7}}, {"Numerator", {1.0}}, {"Denominator", {1.0}}}}},
      {"Conductivity", {{{"Attributes", {8}}, {"Conductivity", 1.0e7}}}},
      {"Superconductor", {{{"Attributes", {9}}, {"KineticInductance", 1.0e-12}}}},
      {"FluxLoop",
       {{{"Index", 3},
         {"FilmAttributes", {10}},
         {"HoleAttributes", {11}},
         {"FluxAmounts", {1.0}}}}},
      {"LumpedPort",
       {{{"Index", 4}, {"Attributes", {12}}, {"R", 50.0}, {"Direction", "+X"}}}},
      {"SurfaceCurrent", {{{"Index", 5}, {"Attributes", {13}}, {"Direction", "+X"}}}},
      {"Absorbing", {{"Attributes", {14}}}}};
  config::BoundaryData data(boundaries);
  const auto metal = data.GetMetalSheetAttributes();
  CHECK(metal == std::vector<int>{1, 3, 4, 5, 6, 7, 8, 9, 10});
  auto sheets = metal;
  sheets.insert(sheets.end(), {12, 13});
  CHECK(data.GetSheetAttributes() == sheets);

  // Without lumped ports and surface current sources the two lists are equal.
  boundaries.erase("LumpedPort");
  boundaries.erase("SurfaceCurrent");
  config::BoundaryData metal_data(boundaries);
  CHECK(metal_data.GetMetalSheetAttributes() == metal);
  CHECK(metal_data.GetSheetAttributes() == metal);
}

TEST_CASE("Config Terminal and PrescribedPotential are exclusive",
          "[config][prescribedpotential][Serial]")
{
  // A configuration may not specify both Terminal and PrescribedPotential boundaries (the
  // electrostatic terminal excitations and the coupon's prescribed-potential basis): the
  // boundary validation refuses the pair, so IoData does; either one alone passes it.
  json boundaries = {{"Ground", {{"Attributes", {1}}}},
                     {"Terminal", {{{"Index", 1}, {"Attributes", {2}}}}},
                     {"PrescribedPotential",
                      {{{"Index", 2}, {"Attributes", {3}}, {"DataFile", "trace.csv"}}}}};
  const auto error = config::Validate(config::BoundaryData(boundaries));
  REQUIRE(error.has_value());
  CHECK_THAT(*error, Catch::Matchers::ContainsSubstring(
                         "\"Terminal\" and \"PrescribedPotential\" cannot both be "
                         "specified"));
  for (const auto *key : {"Terminal", "PrescribedPotential"})
  {
    INFO("without " << key);
    auto single = boundaries;
    single.erase(key);
    CHECK_FALSE(config::Validate(config::BoundaryData(single)).has_value());
  }
  json config = {{"Problem", {{"Type", "Electrostatic"}}},
                 {"Model", {{"Mesh", "mesh.msh"}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                 {"Boundaries", boundaries}};
  CHECK_THROWS_WITH(IoData(config, false),
                    Catch::Matchers::ContainsSubstring(
                        "\"Terminal\" and \"PrescribedPotential\" cannot both be "
                        "specified"));
}

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
    CHECK(boundary_ex_bool.waveport.at(4).target == 0.0);
    CHECK(config::WavePortData(json{{"Index", 7}, {"Attributes", {19}}, {"Target", 3.1}})
              .target == 3.1);

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

  SECTION("FloquetPort repeated index cross-array fails validation")
  {
    json boundaries = {
        {"LumpedPort",
         {{{"Attributes", {5}}, {"Index", 1}, {"R", 50}, {"Direction", "+Y"}}}},
        {"FloquetPort", {{{"Attributes", {6}}, {"Index", 1}}}}};
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

  SECTION("FloquetPort mislabeled excitation index fails validation")
  {
    json boundaries = {
        {"FloquetPort", {{{"Attributes", {6}}, {"Index", 3}, {"Excitation", 2}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(config::Validate(boundary_data).has_value());
  }

  SECTION("FloquetPort bool excitation upgrades to port index")
  {
    json boundaries = {
        {"FloquetPort", {{{"Attributes", {6}}, {"Index", 3}, {"Excitation", true}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(boundary_data.floquetport.at(3).excitation == 3);
    CHECK(!config::Validate(boundary_data).has_value());
  }

  SECTION("FloquetPort excitation is not upgraded for invalid excitation maps")
  {
    json boundaries = {
        {"LumpedPort",
         {{{"Attributes", {5}},
           {"Index", 2},
           {"R", 50},
           {"Direction", "+Y"},
           {"Excitation", 3}}}},
        {"FloquetPort", {{{"Attributes", {6}}, {"Index", 4}, {"Excitation", true}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(boundary_data.floquetport.at(4).excitation == 1);
    CHECK(config::Validate(boundary_data).has_value());
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

  SECTION("PrescribedPotential duplicate index and terminal mixing fail validation")
  {
    json boundaries = {
        {"LumpedPort",
         {{{"Attributes", {5}}, {"Index", 1}, {"R", 50}, {"Direction", "+Y"}}}},
        {"PrescribedPotential",
         {{{"Attributes", {6}}, {"Index", 1}, {"DataFile", "trace.csv"}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(config::Validate(boundary_data).has_value());

    boundaries.erase("LumpedPort");
    boundaries["Terminal"] = {{{"Attributes", {7}}, {"Index", 2}}};
    boundary_data = config::BoundaryData(boundaries);
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

  SECTION("SurfaceCurrent and FluxLoop combination parses")
  {
    json boundaries = {{"SurfaceCurrent",
                        {{{"Attributes", {6}},
                          {"Index", 1},
                          {"Direction", "+X"},
                          {"Aperture", {{"Attributes", {7}}, {"Direction", "+Z"}}}}}},
                       {"FluxLoop",
                        {{{"Index", 2},
                          {"FilmAttributes", {8}},
                          {"HoleAttributes", {9}},
                          {"FluxAmounts", {1.0}},
                          {"Direction", "+Z"}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(!config::Validate(boundary_data).has_value());
    CHECK(boundary_data.current.size() == 1);
    CHECK(boundary_data.fluxloop.size() == 1);
    REQUIRE(boundary_data.current.at(1).elements.at(0).aperture);
    CHECK(boundary_data.current.at(1).elements.at(0).aperture->attributes ==
          std::vector<int>{7});
  }

  SECTION("FluxLoop excitation defaults to the loop index")
  {
    json boundaries = {{"FluxLoop",
                        {{{"Index", 1},
                          {"FilmAttributes", {8}},
                          {"HoleAttributes", {9}},
                          {"FluxAmounts", {1.0}},
                          {"Direction", "+Z"}},
                         {{"Index", 3},
                          {"FilmAttributes", {8}},
                          {"HoleAttributes", {10}},
                          {"FluxAmounts", {1.0}},
                          {"Direction", "+Z"}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(!config::Validate(boundary_data).has_value());
    CHECK(boundary_data.fluxloop.at(1).excitation == 1);
    CHECK(boundary_data.fluxloop.at(3).excitation == 3);
  }

  SECTION("FluxLoop shared excitation groups the loops")
  {
    json boundaries = {{"FluxLoop",
                        {{{"Index", 1},
                          {"FilmAttributes", {8}},
                          {"HoleAttributes", {9}},
                          {"FluxAmounts", {1.0}},
                          {"Direction", "+Z"},
                          {"Excitation", 2}},
                         {{"Index", 3},
                          {"FilmAttributes", {8}},
                          {"HoleAttributes", {10}},
                          {"FluxAmounts", {-0.5}},
                          {"Direction", "+Z"},
                          {"Excitation", 2}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(!config::Validate(boundary_data).has_value());
    CHECK(boundary_data.fluxloop.at(1).excitation == 2);
    CHECK(boundary_data.fluxloop.at(3).excitation == 2);
  }

  SECTION("FluxLoop bool excitation upgrades to the loop index")
  {
    json boundaries = {{"FluxLoop",
                        {{{"Index", 1},
                          {"FilmAttributes", {8}},
                          {"HoleAttributes", {9}},
                          {"FluxAmounts", {1.0}},
                          {"Direction", "+Z"},
                          {"Excitation", false}},
                         {{"Index", 3},
                          {"FilmAttributes", {8}},
                          {"HoleAttributes", {10}},
                          {"FluxAmounts", {1.0}},
                          {"Direction", "+Z"},
                          {"Excitation", true}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(!config::Validate(boundary_data).has_value());
    CHECK(boundary_data.fluxloop.at(1).excitation == 0);
    CHECK(boundary_data.fluxloop.at(3).excitation == 3);
  }

  SECTION("FluxLoop mislabeled single excitation fails validation")
  {
    json boundaries = {{"FluxLoop",
                        {{{"Index", 1},
                          {"FilmAttributes", {8}},
                          {"HoleAttributes", {9}},
                          {"FluxAmounts", {1.0}},
                          {"Direction", "+Z"},
                          {"Excitation", 2}},
                         {{"Index", 3},
                          {"FilmAttributes", {8}},
                          {"HoleAttributes", {10}},
                          {"FluxAmounts", {1.0}},
                          {"Direction", "+Z"}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(config::Validate(boundary_data).has_value());
  }

  SECTION("Grouped FluxLoop excitation must not collide with a SurfaceCurrent index")
  {
    json boundaries = {{"SurfaceCurrent",
                        {{{"Attributes", {6}},
                          {"Index", 2},
                          {"Direction", "+X"},
                          {"Aperture", {{"Attributes", {7}}, {"Direction", "+Z"}}}}}},
                       {"FluxLoop",
                        {{{"Index", 1},
                          {"FilmAttributes", {8}},
                          {"HoleAttributes", {9}},
                          {"FluxAmounts", {1.0}},
                          {"Direction", "+Z"},
                          {"Excitation", 2}},
                         {{"Index", 3},
                          {"FilmAttributes", {8}},
                          {"HoleAttributes", {10}},
                          {"FluxAmounts", {1.0}},
                          {"Direction", "+Z"},
                          {"Excitation", 2}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(config::Validate(boundary_data).has_value());
    boundaries["SurfaceCurrent"][0]["Index"] = 4;
    CHECK(!config::Validate(config::BoundaryData(boundaries)).has_value());
    // A lone loop under a foreign excitation index next to the group collides too.
    boundaries["FluxLoop"].push_back({{"Index", 5},
                                      {"FilmAttributes", {8}},
                                      {"HoleAttributes", {11}},
                                      {"FluxAmounts", {1.0}},
                                      {"Direction", "+Z"},
                                      {"Excitation", 4}});
    auto error = config::Validate(config::BoundaryData(boundaries));
    REQUIRE(error.has_value());
    CHECK_THAT(*error,
               Catch::Matchers::ContainsSubstring("FluxLoop \"Excitation\" 4 collides"));
  }

  SECTION("Excited FluxLoop must have a nonzero FluxAmounts")
  {
    json loop = {{"Index", 1},
                 {"FilmAttributes", {8}},
                 {"HoleAttributes", {9}},
                 {"FluxAmounts", {0.0}},
                 {"Direction", "+Z"}};
    CHECK_THROWS(config::FluxLoopData(loop));
    loop["Excitation"] = 0;
    CHECK(config::FluxLoopData(loop).excitation == 0);
  }

  SECTION("SurfaceCurrent Aperture is element-owned, sorted, and normalized")
  {
    json boundaries = {
        {"SurfaceCurrent",
         {{{"Attributes", {6}}, {"Index", 1}, {"Direction", "+X"}},
          {{"Attributes", {7}},
           {"Index", 2},
           {"Direction", "+X"},
           {"Aperture", {{"Attributes", {7, 5, 6}}, {"Direction", {0.0, 0.0, 2.0}}}}},
          {{"Index", 3},
           {"Elements",
            {{{"Attributes", {8}},
              {"Direction", "+X"},
              {"Aperture", {{"Attributes", {10}}, {"Direction", "+Z"}}}},
             {{"Attributes", {9}},
              {"Direction", "-X"},
              {"Aperture", {{"Attributes", {12, 11}}, {"Direction", "-Z"}}}}}}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(!config::Validate(boundary_data).has_value());
    CHECK(!boundary_data.current.at(1).elements.at(0).aperture);

    const auto &single = boundary_data.current.at(2).elements.at(0).aperture;
    REQUIRE(single);
    CHECK(single->attributes == std::vector<int>{5, 6, 7});
    CHECK(single->direction[0] == 0.0);
    CHECK(single->direction[1] == 0.0);
    CHECK(single->direction[2] == 1.0);

    const auto &multi = boundary_data.current.at(3).elements;
    REQUIRE(multi.size() == 2);
    REQUIRE(multi[0].aperture);
    REQUIRE(multi[1].aperture);
    CHECK(multi[0].aperture->attributes == std::vector<int>{10});
    CHECK(multi[0].aperture->direction[2] == 1.0);
    CHECK(multi[1].aperture->attributes == std::vector<int>{11, 12});
    CHECK(multi[1].aperture->direction[2] == -1.0);
  }

  SECTION("SurfaceCurrent InactiveMode is optional and parses per port")
  {
    json boundaries = {{"SurfaceCurrent",
                        {{{"Attributes", {4}}, {"Index", 1}, {"Direction", "+X"}},
                         {{"Attributes", {5}},
                          {"Index", 2},
                          {"Direction", "+X"},
                          {"InactiveMode", "Open"}},
                         {{"Attributes", {6}},
                          {"Index", 3},
                          {"Direction", "+X"},
                          {"InactiveMode", "Short"}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(!config::Validate(boundary_data).has_value());
    // Unset when not specified so the global default applies for that port.
    CHECK(!boundary_data.current.at(1).inactive_port_mode.has_value());
    REQUIRE(boundary_data.current.at(2).inactive_port_mode.has_value());
    CHECK(boundary_data.current.at(2).inactive_port_mode.value() == InactivePortMode::OPEN);
    REQUIRE(boundary_data.current.at(3).inactive_port_mode.has_value());
    CHECK(boundary_data.current.at(3).inactive_port_mode.value() ==
          InactivePortMode::SHORT);
  }

  SECTION("WavePort ComplexCoarseSolve is optional and parses per port")
  {
    json boundaries = {
        {"WavePort",
         {{{"Attributes", {4}}, {"Index", 1}},
          {{"Attributes", {5}}, {"Index", 2}, {"ComplexCoarseSolve", true}},
          {{"Attributes", {6}}, {"Index", 3}, {"ComplexCoarseSolve", false}}}}};
    config::BoundaryData boundary_data(boundaries);
    CHECK(!config::Validate(boundary_data).has_value());
    // Unset when not specified so the global linear solver settings apply for that port.
    CHECK(!boundary_data.waveport.at(1).complex_coarse_solve.has_value());
    REQUIRE(boundary_data.waveport.at(2).complex_coarse_solve.has_value());
    CHECK(boundary_data.waveport.at(2).complex_coarse_solve.value());
    REQUIRE(boundary_data.waveport.at(3).complex_coarse_solve.has_value());
    CHECK_FALSE(boundary_data.waveport.at(3).complex_coarse_solve.value());
  }
}

TEST_CASE("Config LumpedPort Length override", "[config][Serial]")
{
  // Decision 929: the optional "Length" geometry override of a ragged planar element (mesh
  // units, along the element's Direction): parsed on the single-element and the Elements
  // forms, 0 (absent) by default, refused when not positive or on a coaxial element.
  json without = {
      {"LumpedPort",
       {{{"Attributes", {5}}, {"Index", 1}, {"R", 50.0}, {"Direction", "+Y"}}}}};
  config::BoundaryData data_without(without);
  CHECK(data_without.lumpedport.at(1).elements.size() == 1);
  CHECK(data_without.lumpedport.at(1).elements[0].length == 0.0);

  json with = {{"LumpedPort",
                {{{"Attributes", {5}},
                  {"Index", 1},
                  {"R", 50.0},
                  {"Direction", "+Y"},
                  {"Length", 0.745}},
                 {{"Index", 2},
                  {"L", 1.0e-9},
                  {"Elements",
                   {{{"Attributes", {27}}, {"Direction", "+X"}, {"Length", 0.5}},
                    {{"Attributes", {38}}, {"Direction", "-X"}}}}}}}};
  config::BoundaryData data_with(with);
  CHECK(data_with.lumpedport.at(1).elements[0].length == 0.745);
  CHECK(data_with.lumpedport.at(1).R == 50.0);
  REQUIRE(data_with.lumpedport.at(2).elements.size() == 2);
  CHECK(data_with.lumpedport.at(2).elements[0].length == 0.5);
  CHECK(data_with.lumpedport.at(2).elements[1].length == 0.0);
  CHECK(data_with.lumpedport.at(2).L == 1.0e-9);

  json zero = without;
  zero["LumpedPort"][0]["Length"] = 0.0;
  CHECK_THROWS_WITH(config::BoundaryData(zero),
                    Catch::Matchers::ContainsSubstring("\"Length\" must be positive"));
  json negative = without;
  negative["LumpedPort"][0]["Length"] = -1.0;
  CHECK_THROWS_WITH(config::BoundaryData(negative),
                    Catch::Matchers::ContainsSubstring("\"Length\" must be positive"));
  json coaxial = without;
  coaxial["LumpedPort"][0]["Direction"] = "+R";
  coaxial["LumpedPort"][0]["Length"] = 0.745;
  CHECK_THROWS_WITH(config::BoundaryData(coaxial),
                    Catch::Matchers::ContainsSubstring("Cartesian (uniform) element only"));
}

TEST_CASE("Config interface dielectric flux recovery", "[config][Serial]")
{
  json dielectric = {{"Attributes", {1}},
                     {"Thickness", 2.0e-3},
                     {"Permittivity", 4.0},
                     {"FluxRecovery", true}};
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));

  dielectric["Type"] = "SA";
  CHECK_NOTHROW(config::InterfaceDielectricData(dielectric));
}

TEST_CASE("Config electrostatic response matrix", "[config][Serial]")
{
  json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                 {"Model", {{"Mesh", "test.msh"}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                 {"Boundaries", json::object()},
                 {"Solver", {{"Electrostatic", {{"ResponseMatrix", true}}}}}};
  CHECK_THROWS(IoData(config, false));

  config["Boundaries"]["PrescribedPotential"] = {
      {{"Index", 1}, {"Attributes", {1}}, {"DataFile", "basis.csv"}}};
  CHECK_THROWS(IoData(config, false));

  config["Boundaries"]["Postprocessing"]["Dielectric"] = {
      {{"Index", 1},
       {"Attributes", {2}},
       {"Thickness", 2.0e-3},
       {"Permittivity", 4.0},
       {"EdgeAttributes", {1}},
       {"EdgeDistances", {0.2}},
       {"LocalizeEdgeEnergy", true},
       {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}};
  CHECK_NOTHROW(IoData(config, false));

  config["Solver"]["Electrostatic"]["ResponseMatrix"] = false;
  config["Solver"]["Electrostatic"]["AggregateResponseMatrix"] = true;
  CHECK_THROWS(IoData(config, false));
  config["Solver"]["Electrostatic"]["ResponseMatrix"] = true;

  config["Boundaries"]["PrescribedPotential"][0]["TerminalAttributes"] = {3};
  const IoData conductor_basis(config, false);
  CHECK(conductor_basis.boundaries.prescribed_potential.at(1).terminal_attributes ==
        std::vector<int>{3});

  config["Boundaries"]["PrescribedPotential"][0]["TerminalAttributes"] = {1};
  CHECK_THROWS(IoData(config, false));
}

// Decision 352 follow-up (2): the optional quadrature-level Region of an interface
// dielectric entry (a box and / or translational cells of segments), its defaults and its
// refused inputs; the resolved configuration records the default normal.
TEST_CASE("Config interface dielectric region", "[config][Serial]")
{
  json dielectric = {{"Attributes", {1}},
                     {"Type", "MS"},
                     {"Thickness", 2.0e-3},
                     {"Permittivity", 11.45},
                     {"Region",
                      {{"Segments", {{0.0, 0.0, 0.0, 10.0, 0.0, 0.0}}},
                       {"Distance", 5.0},
                       {"BoxMin", {-1.0, -1.0, -0.5}},
                       {"BoxMax", {11.0, 1.0, 0.5}}}}};
  CHECK_FALSE(config::InterfaceDielectricData(
                  json{{"Attributes", {1}}, {"Thickness", 2.0e-3}, {"Permittivity", 4.0}})
                  .region);
  const config::InterfaceDielectricData data(dielectric);
  REQUIRE(data.region);
  REQUIRE(data.region->box_min);
  CHECK(*data.region->box_min == std::array<double, 3>{-1.0, -1.0, -0.5});
  CHECK(*data.region->box_max == std::array<double, 3>{11.0, 1.0, 0.5});
  REQUIRE(data.region->segments.size() == 1);
  CHECK(data.region->segments[0] == std::array<double, 6>{0.0, 0.0, 0.0, 10.0, 0.0, 0.0});
  CHECK(data.region->distance == 5.0);
  CHECK(data.region->normal == std::array<double, 3>{0.0, 0.0, 1.0});

  // A box alone, segments alone, a custom normal.
  {
    auto box_only = dielectric;
    box_only["Region"].erase("Segments");
    box_only["Region"].erase("Distance");
    const config::InterfaceDielectricData parsed(box_only);
    REQUIRE(parsed.region);
    CHECK(parsed.region->segments.empty());
    CHECK(parsed.region->box_min);
    auto segments_only = dielectric;
    segments_only["Region"].erase("BoxMin");
    segments_only["Region"].erase("BoxMax");
    segments_only["Region"]["Normal"] = {1.0, 0.0, 0.0};
    const config::InterfaceDielectricData cells(segments_only);
    REQUIRE(cells.region);
    CHECK_FALSE(cells.region->box_min);
    CHECK(cells.region->normal == std::array<double, 3>{1.0, 0.0, 0.0});
  }
  // Refused: an empty region, one box bound, segments without a distance (and the
  // converse), a non-positive distance, unordered box bounds, a zero-length segment, a
  // zero normal.
  for (const auto &[reason, mutate] :
       std::vector<std::pair<const char *, std::function<void(json &)>>>{
           {"empty", [](json &r) { r = json::object(); }},
           {"one bound", [](json &r) { r.erase("BoxMax"); }},
           {"no distance", [](json &r) { r.erase("Distance"); }},
           {"distance alone",
            [](json &r)
            {
              r.erase("Segments");
              r.erase("BoxMin");
              r.erase("BoxMax");
            }},
           {"zero distance", [](json &r) { r["Distance"] = 0.0; }},
           {"unordered box", [](json &r) { r["BoxMin"] = {12.0, -1.0, -0.5}; }},
           {"zero segment",
            [](json &r) { r["Segments"] = {{1.0, 1.0, 0.0, 1.0, 1.0, 0.0}}; }},
           {"zero normal", [](json &r) { r["Normal"] = {0.0, 0.0, 0.0}; }}})
  {
    INFO(reason);
    auto refused = dielectric;
    mutate(refused["Region"]);
    CHECK_THROWS(config::InterfaceDielectricData(refused));
  }

  // The parsed region keeps the configured (mesh-unit) values (the scaling with the mesh
  // is config::Nondimensionalize's); the resolved configuration records the default normal.
  dielectric["Index"] = 1;
  json input = {{"Problem", {{"Type", "Electrostatic"}}},
                {"Model", {{"Mesh", "unused.msh"}}},
                {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                {"Boundaries",
                 {{"Ground", {{"Attributes", {2}}}},
                  {"Terminal", {{{"Index", 1}, {"Attributes", {1}}}}},
                  {"Postprocessing", {{"Dielectric", {dielectric}}}}}},
                {"Solver", {{"Electrostatic", {{"Save", 0}}}}}};
  IoData parsed(input, false);
  const auto &region = parsed.boundaries.postpro.dielectric.at(1).region;
  REQUIRE(region);
  CHECK(region->distance == 5.0);
  CHECK(region->segments[0][3] == 10.0);
  CHECK((*region->box_max)[0] == 11.0);
  CHECK(region->normal == std::array<double, 3>{0.0, 0.0, 1.0});
  const auto concrete = IoData::ConcretizeDefaults(parsed, input);
  const auto &entry = concrete["Boundaries"]["Postprocessing"]["Dielectric"][0];
  CHECK(entry.at("Region").at("Normal") == json({0.0, 0.0, 1.0}));
  CHECK(entry.at("Region").at("Distance") == 5.0);
}

TEST_CASE("Config electrostatic response correction", "[config][Serial]")
{
  const json correction = {{"FabricatedMatrix", "fabricated.csv"},
                           {"ThinMatrix", "thin.csv"},
                           {"BasisPoints", "points.csv"},
                           {"Edges",
                            {{{"Origin", {1.0, 2.0, 0.0}},
                              {"AxisU", {1.0, 0.0, 0.0}},
                              {"AxisV", {0.0, 1.0, 0.0}}}}}};
  const config::ElectrostaticSolverData electrostatic(
      json{{"ResponseCorrection", correction}});
  REQUIRE(electrostatic.response_correction);
  using CorrectionMode =
      config::ElectrostaticSolverData::ResponseCorrectionData::CorrectionMode;
  using DomainCorrection = config::ElectrostaticSolverData::ResponseCorrectionData::
      TranslationalDomainCorrection;
  // The v1 default path: post-processing only; the developer switch is off.
  CHECK(electrostatic.response_correction->correction_mode ==
        CorrectionMode::POSTPROCESS_ONLY);
  CHECK_FALSE(electrostatic.response_correction->developer);
  CHECK(electrostatic.response_correction->developer_options.empty());
  CHECK(electrostatic.response_correction->translational_domain_correction ==
        DomainCorrection::FIXED_TRACE);
  CHECK(electrostatic.response_correction->IncludesPostprocessing());
  CHECK_FALSE(electrostatic.response_correction->IncludesSelfConsistent());
  CHECK(electrostatic.response_correction->solve_tol == 1.0e-6);
  REQUIRE(electrostatic.response_correction->models.size() == 1);
  CHECK(electrostatic.response_correction->models[0].idx == 1);
  CHECK(electrostatic.response_correction->models[0].fabricated_matrix == "fabricated.csv");
  CHECK(electrostatic.response_correction->models[0].thin_matrix == "thin.csv");
  CHECK(electrostatic.response_correction->models[0].basis_points == "points.csv");
  REQUIRE(electrostatic.response_correction->patches.size() == 1);
  CHECK(electrostatic.response_correction->patches[0].model == 1);
  CHECK(electrostatic.response_correction->patches[0].origin ==
        std::array<double, 3>{1.0, 2.0, 0.0});

  json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                 {"Model", {{"Mesh", "test.msh"}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                 {"Boundaries", json::object()},
                 {"Solver", {{"Electrostatic", {{"ResponseCorrection", correction}}}}}};
  const IoData iodata(config, false);

  // Concretization must also work when the raw JSON does not already contain the nested
  // correction object, as happens for programmatically assembled configurations.
  json sparse = config;
  sparse["Solver"]["Electrostatic"].erase("ResponseCorrection");
  sparse = IoData::ConcretizeDefaults(iodata, sparse);
  auto concretized_correction = correction;
  concretized_correction["Developer"] = false;
  concretized_correction["CorrectionMode"] = "PostprocessOnly";
  concretized_correction["SolveTol"] = 1.0e-6;
  concretized_correction["TranslationalDomainCorrection"] = "FixedTrace";
  concretized_correction["TraceCoupling"] = "Collocated";
  concretized_correction["MortarOversampling"] = 2;
  concretized_correction["PatchEnergy"] = false;
  CHECK(sparse["Solver"]["Electrostatic"]["ResponseCorrection"] == concretized_correction);

  config["Problem"]["Type"] = "Magnetostatic";
  CHECK_THROWS(IoData(config, false));

  const json multi_model = {
      {"Models",
       {{{"Index", 7},
         {"FabricatedMatrix", "fabricated-domain.csv"},
         {"ThinMatrix", "thin-domain.csv"},
         {"FabricatedSurfaceMatrix", "fabricated-surface.csv"},
         {"ThinSurfaceMatrix", "thin-surface.csv"},
         {"BasisPoints", "pair-points.csv"},
         {"TraceMesh",
          {{"Vertices", "trace-vertices.csv"}, {"Triangles", "trace-triangles.csv"}}},
         {"Interfaces", {{{"Target", 4}, {"Coupon", 1}}}}}}},
      {"Patches",
       {{{"Model", 7},
         {"Origin", {3.0, 4.0, 0.0}},
         {"AxisU", {1.0, 0.0, 0.0}},
         {"AxisV", {0.0, 1.0, 0.0}},
         {"Reference", {-1.0, 0.0, 0.0}}}}}};
  const config::ElectrostaticSolverData modern(json{{"ResponseCorrection", multi_model}});
  REQUIRE(modern.response_correction);
  REQUIRE(modern.response_correction->models.size() == 1);
  CHECK(modern.response_correction->models[0].idx == 7);
  REQUIRE(modern.response_correction->models[0].interfaces.size() == 1);
  CHECK(modern.response_correction->models[0].interfaces[0].target == 4);
  CHECK(modern.response_correction->models[0].interfaces[0].coupon == 1);
  CHECK(modern.response_correction->models[0].trace_vertices == "trace-vertices.csv");
  CHECK(modern.response_correction->models[0].trace_triangles == "trace-triangles.csv");
  REQUIRE(modern.response_correction->patches.size() == 1);
  CHECK(modern.response_correction->patches[0].model == 7);
  CHECK(modern.response_correction->patches[0].conductor_references.front() ==
        std::array<double, 3>{-1.0, 0.0, 0.0});

  const json automatic_correction = {{"Library", "fabrication-process.json"},
                                     {"TargetInterfaces", {4, 5, 6}},
                                     {"UnmatchedPolicy", "Error"}};
  const config::ElectrostaticSolverData automatic(
      json{{"ResponseCorrection", automatic_correction}});
  REQUIRE(automatic.response_correction);
  CHECK(automatic.response_correction->IsAutomatic());
  CHECK(automatic.response_correction->library == "fabrication-process.json");
  CHECK(automatic.response_correction->target_interfaces == std::vector<int>{4, 5, 6});
  CHECK(automatic.response_correction->unmatched_policy ==
        config::ElectrostaticSolverData::ResponseCorrectionData::UnmatchedPolicy::ERROR);
  CHECK(automatic.response_correction->solve_tol == 1.0e-6);
  CHECK(automatic.response_correction->correction_mode == CorrectionMode::POSTPROCESS_ONLY);
  {
    // SolveTol configures the self-consistent solve alone: without it the key is not read
    // and is named among the ignored keys; with the (parked) solve it is read.
    auto solve_tol = automatic_correction;
    solve_tol["SolveTol"] = 1.0e-7;
    const config::ElectrostaticSolverData ignored(json{{"ResponseCorrection", solve_tol}});
    REQUIRE(ignored.response_correction);
    CHECK(ignored.response_correction->solve_tol == 1.0e-6);
    CHECK(ignored.response_correction->ignored_keys ==
          std::vector<std::string>{"\"SolveTol\" (no self-consistent solve)"});
    CHECK(ignored.response_correction->developer_options.empty());
    solve_tol["CorrectionMode"] = "Both";
    solve_tol["Developer"] = true;
    const config::ElectrostaticSolverData developer(
        json{{"ResponseCorrection", solve_tol}});
    REQUIRE(developer.response_correction);
    CHECK(developer.response_correction->solve_tol == 1.0e-7);
    CHECK(developer.response_correction->ignored_keys.empty());
    CHECK(developer.response_correction->developer_options ==
          std::vector<std::string>{"\"CorrectionMode\": \"Both\""});
    auto not_boolean = automatic_correction;
    not_boolean["Developer"] = "yes";
    CHECK_THROWS_WITH(
        config::ElectrostaticSolverData(json{{"ResponseCorrection", not_boolean}}),
        Catch::Matchers::ContainsSubstring("\"Developer\" must be a boolean"));
  }
  CHECK(automatic.response_correction->models.empty());
  CHECK(automatic.response_correction->patches.empty());
  // The Features-driven patch construction is the default; Legacy only by explicit request.
  CHECK(
      automatic.response_correction->patch_construction ==
      config::ElectrostaticSolverData::ResponseCorrectionData::PatchConstruction::FEATURES);
  auto legacy_construction = automatic_correction;
  legacy_construction["PatchConstruction"] = "Legacy";
  CHECK(config::ElectrostaticSolverData(json{{"ResponseCorrection", legacy_construction}})
            .response_correction->patch_construction ==
        config::ElectrostaticSolverData::ResponseCorrectionData::PatchConstruction::LEGACY);
  auto invalid_construction = automatic_correction;
  invalid_construction["PatchConstruction"] = "Groups";
  CHECK_THROWS(
      config::ElectrostaticSolverData(json{{"ResponseCorrection", invalid_construction}}));

  auto duplicate_targets = automatic_correction;
  duplicate_targets["TargetInterfaces"] = {4, 4};
  CHECK_THROWS(
      config::ElectrostaticSolverData(json{{"ResponseCorrection", duplicate_targets}}));

  auto invalid_policy = automatic_correction;
  invalid_policy["UnmatchedPolicy"] = "Ignore";
  CHECK_THROWS(
      config::ElectrostaticSolverData(json{{"ResponseCorrection", invalid_policy}}));

  // The knife-edge census flag (the identification fix round DESIGN section 9 (i); decision
  // 1081): absent = unset (the operator decides: on for the preflight, off for a solve),
  // an explicit boolean overrides, anything else is refused.
  CHECK_FALSE(automatic.response_correction->knife_edge_census.has_value());
  for (const bool census : {true, false})
  {
    auto census_config = automatic_correction;
    census_config["KnifeEdgeCensus"] = census;
    const config::ElectrostaticSolverData census_data(
        json{{"ResponseCorrection", census_config}});
    REQUIRE(census_data.response_correction);
    REQUIRE(census_data.response_correction->knife_edge_census.has_value());
    CHECK(*census_data.response_correction->knife_edge_census == census);
  }
  auto invalid_census = automatic_correction;
  invalid_census["KnifeEdgeCensus"] = "On";
  CHECK_THROWS(
      config::ElectrostaticSolverData(json{{"ResponseCorrection", invalid_census}}));

  for (const auto &[name, mode, postprocess, self_consistent] :
       std::vector<std::tuple<const char *, CorrectionMode, bool, bool>>{
           {"PostprocessOnly", CorrectionMode::POSTPROCESS_ONLY, true, false},
           {"SelfConsistent", CorrectionMode::SELF_CONSISTENT, false, true},
           {"Both", CorrectionMode::BOTH, true, true}})
  {
    // SelfConsistent and Both are parked: refused without the developer switch, accepted
    // and listed with it; PostprocessOnly needs no switch.
    auto mode_config = automatic_correction;
    mode_config["CorrectionMode"] = name;
    if (self_consistent)
    {
      CHECK_THROWS_WITH(
          config::ElectrostaticSolverData(json{{"ResponseCorrection", mode_config}}),
          Catch::Matchers::ContainsSubstring(fmt::format(
              "\"CorrectionMode\": \"{}\" is a parked developer option", name)) &&
              Catch::Matchers::ContainsSubstring("\"Developer\": true"));
      mode_config["Developer"] = true;
    }
    const config::ElectrostaticSolverData mode_data(
        json{{"ResponseCorrection", mode_config}});
    REQUIRE(mode_data.response_correction);
    CHECK(mode_data.response_correction->correction_mode == mode);
    CHECK(mode_data.response_correction->IncludesPostprocessing() == postprocess);
    CHECK(mode_data.response_correction->IncludesSelfConsistent() == self_consistent);
    CHECK(mode_data.response_correction->developer == self_consistent);
    CHECK(mode_data.response_correction->developer_options ==
          (self_consistent
               ? std::vector<std::string>{fmt::format("\"CorrectionMode\": \"{}\"", name)}
               : std::vector<std::string>{}));
  }
  auto invalid_mode = automatic_correction;
  invalid_mode["CorrectionMode"] = "postprocessonly";
  CHECK_THROWS(config::ElectrostaticSolverData(json{{"ResponseCorrection", invalid_mode}}));

  for (const auto &[name, mode] : std::vector<std::pair<const char *, DomainCorrection>>{
           {"Disabled", DomainCorrection::DISABLED},
           {"FixedTrace", DomainCorrection::FIXED_TRACE},
           {"FixedFlux", DomainCorrection::FIXED_FLUX}})
  {
    // A key of the self-consistent solve alone: read with that (parked) solve; without it
    // not read (the default FixedTrace stays) and named among the ignored keys.
    auto domain_config = automatic_correction;
    domain_config["TranslationalDomainCorrection"] = name;
    const config::ElectrostaticSolverData ignored(
        json{{"ResponseCorrection", domain_config}});
    REQUIRE(ignored.response_correction);
    CHECK(ignored.response_correction->translational_domain_correction ==
          DomainCorrection::FIXED_TRACE);
    CHECK(ignored.response_correction->ignored_keys ==
          std::vector<std::string>{
              "\"TranslationalDomainCorrection\" (no self-consistent solve)"});
    domain_config["CorrectionMode"] = "Both";
    domain_config["Developer"] = true;
    const config::ElectrostaticSolverData domain_data(
        json{{"ResponseCorrection", domain_config}});
    REQUIRE(domain_data.response_correction);
    CHECK(domain_data.response_correction->translational_domain_correction == mode);
  }
  auto invalid_domain_mode = automatic_correction;
  invalid_domain_mode["CorrectionMode"] = "Both";
  invalid_domain_mode["Developer"] = true;
  invalid_domain_mode["TranslationalDomainCorrection"] = "disabled";
  CHECK_THROWS(
      config::ElectrostaticSolverData(json{{"ResponseCorrection", invalid_domain_mode}}));

  for (const auto &[name, mode] :
       std::vector<std::pair<const char *, config::ElectrostaticSolverData::
                                               ResponseCorrectionData::TraceCoupling>>{
           {"Collocated", config::ElectrostaticSolverData::ResponseCorrectionData::
                              TraceCoupling::COLLOCATED},
           {"SurfaceMortar", config::ElectrostaticSolverData::ResponseCorrectionData::
                                 TraceCoupling::SURFACE_MORTAR}})
  {
    auto trace_config = automatic_correction;
    trace_config["TraceCoupling"] = name;
    const config::ElectrostaticSolverData trace_data(
        json{{"ResponseCorrection", trace_config}});
    REQUIRE(trace_data.response_correction);
    CHECK(trace_data.response_correction->trace_coupling == mode);
  }
  auto invalid_trace_coupling = automatic_correction;
  invalid_trace_coupling["TraceCoupling"] = "Mortar";
  CHECK_THROWS(config::ElectrostaticSolverData(
      json{{"ResponseCorrection", invalid_trace_coupling}}));
  auto oversampled_mortar = automatic_correction;
  oversampled_mortar["MortarOversampling"] = 2;
  const config::ElectrostaticSolverData oversampled_data(
      json{{"ResponseCorrection", oversampled_mortar}});
  REQUIRE(oversampled_data.response_correction);
  CHECK(oversampled_data.response_correction->mortar_oversampling == 2);
  oversampled_mortar["MortarOversampling"] = 0;
  CHECK_THROWS(
      config::ElectrostaticSolverData(json{{"ResponseCorrection", oversampled_mortar}}));
  // Decision 352 follow-up (1): the per-patch energy export, off by default; since decision
  // 922 B the Maxwell SurfaceResponseCorrection parses it too (the driven / eigenmode
  // per-patch export; the boundary-mode refusal is IoData's).
  CHECK_FALSE(oversampled_data.response_correction->patch_energy);
  auto patch_energy = automatic_correction;
  patch_energy["PatchEnergy"] = true;
  const config::ElectrostaticSolverData patch_energy_data(
      json{{"ResponseCorrection", patch_energy}});
  REQUIRE(patch_energy_data.response_correction);
  CHECK(patch_energy_data.response_correction->patch_energy);
  CHECK_NOTHROW(
      config::SolverData(json{{"SurfaceResponseCorrection", automatic_correction}}));
  const config::SolverData maxwell_patch_energy(
      json{{"SurfaceResponseCorrection", patch_energy}});
  REQUIRE(maxwell_patch_energy.surface_response_correction);
  CHECK(maxwell_patch_energy.surface_response_correction->patch_energy);
  // Diagnostics.UncoveredFeatureEnergy (decision 1567 (b)): off by default, parsed on the
  // Maxwell path, refused on the electrostatic path (no silent no-op).
  CHECK_FALSE(maxwell_patch_energy.surface_response_correction->uncovered_feature_energy);
  auto uncovered_feature = automatic_correction;
  uncovered_feature["Diagnostics"] = {{"UncoveredFeatureEnergy", true}};
  const config::SolverData maxwell_uncovered_feature(
      json{{"SurfaceResponseCorrection", uncovered_feature}});
  REQUIRE(maxwell_uncovered_feature.surface_response_correction);
  CHECK(maxwell_uncovered_feature.surface_response_correction->uncovered_feature_energy);
  CHECK_THROWS(config::SolverData(
      json{{"Electrostatic", {{"ResponseCorrection", uncovered_feature}}}}));
  uncovered_feature["Diagnostics"] = {{"UncoveredFeatureEnergies", true}};
  CHECK_THROWS(config::SolverData(json{{"SurfaceResponseCorrection", uncovered_feature}}));
  // Coverage (lane L1, decision 1634 (4)): the interim bars by default and FIXED (ruling
  // 1689 M-5: any other value is refused); Policy "Warn" or "MarkNotApplied" (lane L4,
  // decision 1829: the key absent = kDefaultCoveragePolicy = Warn; "Abort" and any other
  // value refused); no unknown key; the same block on both paths.
  CHECK(maxwell_patch_energy.surface_response_correction->coverage.max_uncovered_share ==
        0.05);
  CHECK(maxwell_patch_energy.surface_response_correction->coverage.build_target == 0.02);
  CHECK(
      maxwell_patch_energy.surface_response_correction->coverage.validated_regime.empty());
  auto coverage = automatic_correction;
  coverage["Coverage"] = {{"Policy", "Warn"},
                          {"MaxUncoveredShare", 0.05},
                          {"BuildTarget", 0.02},
                          {"ValidatedRegime", "validated-regime.json"}};
  const config::SolverData maxwell_coverage(json{{"SurfaceResponseCorrection", coverage}});
  REQUIRE(maxwell_coverage.surface_response_correction);
  CHECK(maxwell_coverage.surface_response_correction->coverage.max_uncovered_share == 0.05);
  CHECK(maxwell_coverage.surface_response_correction->coverage.build_target == 0.02);
  CHECK(maxwell_coverage.surface_response_correction->coverage.validated_regime ==
        "validated-regime.json");
  const config::ElectrostaticSolverData electrostatic_coverage(
      json{{"ResponseCorrection", coverage}});
  REQUIRE(electrostatic_coverage.response_correction);
  CHECK(electrostatic_coverage.response_correction->coverage.build_target == 0.02);
  using CoveragePolicy =
      config::ElectrostaticSolverData::ResponseCorrectionData::CoverageData::Policy;
  CHECK(config::ElectrostaticSolverData::ResponseCorrectionData::CoverageData::
            kDefaultCoveragePolicy == CoveragePolicy::WARN);
  CHECK(maxwell_patch_energy.surface_response_correction->coverage.policy ==
        CoveragePolicy::WARN);
  CHECK(maxwell_coverage.surface_response_correction->coverage.policy ==
        CoveragePolicy::WARN);
  {
    auto marked = coverage;
    marked["Coverage"]["Policy"] = "MarkNotApplied";
    const config::SolverData maxwell_marked(json{{"SurfaceResponseCorrection", marked}});
    REQUIRE(maxwell_marked.surface_response_correction);
    CHECK(maxwell_marked.surface_response_correction->coverage.policy ==
          CoveragePolicy::MARK_NOT_APPLIED);
    const config::ElectrostaticSolverData electrostatic_marked(
        json{{"ResponseCorrection", marked}});
    REQUIRE(electrostatic_marked.response_correction);
    CHECK(electrostatic_marked.response_correction->coverage.policy ==
          CoveragePolicy::MARK_NOT_APPLIED);
  }
  for (const auto &[key, value] :
       std::vector<std::pair<std::string, json>>{{"Policy", "Abort"},
                                                 {"Policy", "markNotApplied"},
                                                 {"Policy", 1},
                                                 {"BuildTarget", 0.2},
                                                 {"BuildTarget", 0.01},
                                                 {"MaxUncoveredShare", 1.5},
                                                 {"MaxUncoveredShare", 0.1},
                                                 {"MaxUncoveredShare", 1.0},
                                                 {"MaxMissingLengthShare", 0.01}})
  {
    auto invalid_coverage = coverage;
    invalid_coverage["Coverage"][key] = value;
    CHECK_THROWS(config::SolverData(json{{"SurfaceResponseCorrection", invalid_coverage}}));
  }

  auto invalid_solve_tol = automatic_correction;
  invalid_solve_tol["CorrectionMode"] = "Both";
  invalid_solve_tol["Developer"] = true;
  invalid_solve_tol["SolveTol"] = 0.0;
  CHECK_THROWS(
      config::ElectrostaticSolverData(json{{"ResponseCorrection", invalid_solve_tol}}));

  const config::SolverData maxwell_solver(
      json{{"SurfaceResponseCorrection", automatic_correction}});
  REQUIRE(maxwell_solver.surface_response_correction);
  CHECK(maxwell_solver.surface_response_correction->IsAutomatic());
  CHECK(maxwell_solver.surface_response_correction->library == "fabrication-process.json");

  json maxwell_config = {{"Problem", {{"Type", "Eigenmode"}, {"Output", "test_output"}}},
                         {"Model", {{"Mesh", "test.msh"}}},
                         {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                         {"Boundaries", json::object()},
                         {"Solver",
                          {{"Eigenmode", {{"Target", 1.0}}},
                           {"SurfaceResponseCorrection", automatic_correction}}}};
  CHECK_NOTHROW(IoData(maxwell_config, false));
  auto boundary_mode_config = maxwell_config;
  boundary_mode_config["Problem"]["Type"] = "BoundaryMode";
  boundary_mode_config["Solver"].erase("Eigenmode");
  boundary_mode_config["Solver"]["BoundaryMode"] = {{"Freq", 5.0}};
  CHECK_NOTHROW(IoData(boundary_mode_config, false));
  boundary_mode_config["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] =
      "PostprocessOnly";
  CHECK_NOTHROW(IoData(boundary_mode_config, false));
  boundary_mode_config["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] =
      "SelfConsistent";
  CHECK_THROWS_WITH(IoData(boundary_mode_config, false),
                    Catch::Matchers::ContainsSubstring("is a parked developer option"));
  boundary_mode_config["Solver"]["SurfaceResponseCorrection"]["Developer"] = true;
  CHECK_THROWS_WITH(IoData(boundary_mode_config, false),
                    Catch::Matchers::ContainsSubstring(
                        "Boundary mode does not support SurfaceResponseCorrection"));
  // OI-MAXWELL-TRANSLATIONAL-MODE-INERT (decisions 1237 / 1260): the Maxwell path never
  // reads TranslationalDomainCorrection, so only the default FixedTrace is accepted (an
  // explicit FixedTrace parses; Disabled / FixedFlux are refused by name).
  // TranslationalDomainCorrection is a key of the (parked) self-consistent solve alone:
  // without that solve it is not read, whatever its value.
  auto translational_maxwell = maxwell_config;
  for (const char *mode : {"FixedTrace", "Disabled", "FixedFlux"})
  {
    translational_maxwell["Solver"]["SurfaceResponseCorrection"]
                         ["TranslationalDomainCorrection"] = mode;
    const IoData ignored(translational_maxwell, false);
    CHECK(ignored.solver.surface_response_correction->ignored_keys ==
          std::vector<std::string>{
              "\"TranslationalDomainCorrection\" (no self-consistent solve)"});
  }
  translational_maxwell["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] = "Both";
  translational_maxwell["Solver"]["SurfaceResponseCorrection"]["Developer"] = true;
  translational_maxwell["Solver"]["SurfaceResponseCorrection"]
                       ["TranslationalDomainCorrection"] = "FixedTrace";
  CHECK_NOTHROW(IoData(translational_maxwell, false));
  for (const char *mode : {"Disabled", "FixedFlux"})
  {
    translational_maxwell["Solver"]["SurfaceResponseCorrection"]
                         ["TranslationalDomainCorrection"] = mode;
    CHECK_THROWS_WITH(
        IoData(translational_maxwell, false),
        Catch::Matchers::ContainsSubstring("OI-MAXWELL-TRANSLATIONAL-MODE-INERT"));
  }

  maxwell_config["Problem"]["Type"] = "Electrostatic";
  maxwell_config["Solver"].erase("Eigenmode");
  maxwell_config["Solver"]["Electrostatic"] = json::object();
  CHECK_THROWS(IoData(maxwell_config, false));

  auto explicit_maxwell = maxwell_config;
  explicit_maxwell["Problem"]["Type"] = "Eigenmode";
  explicit_maxwell["Solver"].erase("Electrostatic");
  explicit_maxwell["Solver"]["Eigenmode"] = {{"Target", 1.0}};
  explicit_maxwell["Solver"]["SurfaceResponseCorrection"] = multi_model;
  CHECK_THROWS(IoData(explicit_maxwell, false));
}

// Decision 1685 (lane diagdir): Diagnostics.Directory (one directory name inside
// Problem.Output for the per-patch tables) and Diagnostics.ConductorConsistencyRecords
// ("Metadata" | "Directory"), both absent by default = every output as before; the same
// keys on the Maxwell SurfaceResponseCorrection.
TEST_CASE("Config surface-response Diagnostics.Directory and ConductorConsistencyRecords",
          "[config][Serial]")
{
  using Records = config::ElectrostaticSolverData::ResponseCorrectionData::
      ConductorConsistencyRecordsOutput;
  const json automatic_correction = {{"Library", "fabrication-process.json"},
                                     {"TargetInterfaces", {4, 5, 6}},
                                     {"UnmatchedPolicy", "Error"}};
  const config::ElectrostaticSolverData defaults(
      json{{"ResponseCorrection", automatic_correction}});
  REQUIRE(defaults.response_correction);
  CHECK(defaults.response_correction->diagnostics_directory.empty());
  CHECK(defaults.response_correction->conductor_consistency_records == Records::METADATA);

  auto moved = automatic_correction;
  moved["Diagnostics"] = {{"Directory", "diagnostics"},
                          {"ConductorConsistencyRecords", "Directory"}};
  const config::ElectrostaticSolverData electrostatic(json{{"ResponseCorrection", moved}});
  REQUIRE(electrostatic.response_correction);
  CHECK(electrostatic.response_correction->diagnostics_directory == "diagnostics");
  CHECK(electrostatic.response_correction->conductor_consistency_records ==
        Records::DIRECTORY);
  const config::SolverData maxwell(json{{"SurfaceResponseCorrection", moved}});
  REQUIRE(maxwell.surface_response_correction);
  CHECK(maxwell.surface_response_correction->diagnostics_directory == "diagnostics");
  CHECK(maxwell.surface_response_correction->conductor_consistency_records ==
        Records::DIRECTORY);

  auto directory_only = automatic_correction;
  directory_only["Diagnostics"] = {{"Directory", "heavy.tables_1"}};
  const config::ElectrostaticSolverData only(json{{"ResponseCorrection", directory_only}});
  REQUIRE(only.response_correction);
  CHECK(only.response_correction->diagnostics_directory == "heavy.tables_1");
  CHECK(only.response_correction->conductor_consistency_records == Records::METADATA);

  for (const auto *name :
       {"", "a/b", "..", ".", ".hidden", "-dash", "iteration1", "iterations", "with space"})
  {
    CAPTURE(name);
    auto invalid = automatic_correction;
    invalid["Diagnostics"] = {{"Directory", name}};
    CHECK_THROWS_WITH(
        config::ElectrostaticSolverData(json{{"ResponseCorrection", invalid}}),
        Catch::Matchers::ContainsSubstring("Diagnostics.Directory"));
  }
  auto not_a_string = automatic_correction;
  not_a_string["Diagnostics"] = {{"Directory", 1}};
  CHECK_THROWS_WITH(
      config::ElectrostaticSolverData(json{{"ResponseCorrection", not_a_string}}),
      Catch::Matchers::ContainsSubstring("Diagnostics.Directory"));
  auto records_without_directory = automatic_correction;
  records_without_directory["Diagnostics"] = {{"ConductorConsistencyRecords", "Directory"}};
  CHECK_THROWS_WITH(config::ElectrostaticSolverData(
                        json{{"ResponseCorrection", records_without_directory}}),
                    Catch::Matchers::ContainsSubstring("requires"));
  auto invalid_records = automatic_correction;
  invalid_records["Diagnostics"] = {{"Directory", "diagnostics"},
                                    {"ConductorConsistencyRecords", "palace.json"}};
  CHECK_THROWS_WITH(
      config::ElectrostaticSolverData(json{{"ResponseCorrection", invalid_records}}),
      Catch::Matchers::ContainsSubstring("ConductorConsistencyRecords"));
  auto unknown = automatic_correction;
  unknown["Diagnostics"] = {{"Directories", "diagnostics"}};
  CHECK_THROWS_WITH(config::ElectrostaticSolverData(json{{"ResponseCorrection", unknown}}),
                    Catch::Matchers::ContainsSubstring("unknown key"));
}

TEST_CASE("Config axisymmetric model electrostatic only", "[config][Serial]")
{
  // The 2 pi r measure is carried by the libCEED geometry factors and the electrostatic
  // surface post-processing only (curved coupons on (r, z) meshes): every other problem
  // type is refused with Model.Axisymmetric.
  json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                 {"Model", {{"Mesh", "test.msh"}, {"Axisymmetric", true}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                 {"Boundaries", json::object()},
                 {"Solver", {{"Electrostatic", json::object()}}}};
  const IoData electrostatic(config, false);
  CHECK(electrostatic.model.axisymmetric);

  for (const auto &[type, solver] : std::vector<std::pair<std::string, json>>{
           {"Magnetostatic", json::object()},
           {"Eigenmode", {{"Target", 1.0}}},
           {"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.5}}},
           {"Transient",
            {{"Excitation", "Gaussian"}, {"MaxTime", 1.0}, {"TimeStep", 0.1}}}})
  {
    json refused = config;
    refused["Problem"]["Type"] = type;
    refused["Solver"] = {{type, solver}};
    CHECK_THROWS_WITH(IoData(refused, false),
                      Catch::Matchers::ContainsSubstring(
                          "\"Model.Axisymmetric\" is only supported for electrostatic"));
    refused["Model"]["Axisymmetric"] = false;
    CHECK_NOTHROW(IoData(refused, false));
  }
}

TEST_CASE("Config interface dielectric edge refinement", "[config][Serial]")
{
  json dielectric = {{"Attributes", {1}},
                     {"Thickness", 2.0e-3},
                     {"Permittivity", 4.0},
                     {"EdgeAttributes", {2}},
                     {"EdgeDistances", {0.2, 0.5}},
                     {"EdgeRefinement", {{"Radius", 0.2}, {"ElementsPerRadius", 3}}}};
  CHECK_NOTHROW(config::InterfaceDielectricData(dielectric));

  dielectric["EdgeRefinement"]["Radius"] = 0.3;
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
}

TEST_CASE("Config refinement MaxSizePredicted", "[config][Serial]")
{
  // Default = upstream's rule (USER decision 1636 (4)): MaxSize is checked only against the
  // last solved size.
  config::RefinementData unset({{"MaxIts", 4}, {"MaxSize", 100000000}});
  CHECK(unset.max_size == 100000000);
  CHECK_FALSE(unset.max_size_predicted);

  // The predictive stop (MaxSize checked against the predicted size of the refined mesh
  // before its solve) is selected explicitly, as the fabrication-response writers do.
  config::RefinementData predicted(
      {{"MaxIts", 4}, {"MaxSize", 100000000}, {"MaxSizePredicted", true}});
  CHECK(predicted.max_size_predicted);
  config::RefinementData upstream(
      {{"MaxIts", 4}, {"MaxSize", 100000000}, {"MaxSizePredicted", false}});
  CHECK_FALSE(upstream.max_size_predicted);
}

TEST_CASE("Config refinement MinElementSize", "[config][Serial]")
{
  // Default: no element-size floor (every element is eligible for marking).
  json unfloored = {{"MaxIts", 4}};
  CHECK(config::RefinementData(unfloored).min_element_size == 0.0);

  // A floor is a length in mesh units; parsed as written (nondimensionalized later).
  json floored = {{"MaxIts", 4}, {"MinElementSize", 0.02}};
  CHECK(config::RefinementData(floored).min_element_size == 0.02);

  // A negative floor is refused.
  json negative = {{"MaxIts", 4}, {"MinElementSize", -0.02}};
  CHECK_THROWS(config::RefinementData(negative));
}

TEST_CASE("Config automatic interface dielectric edges", "[config][Serial]")
{
  json dielectric = {{"Attributes", {1}},
                     {"Type", "SA"},
                     {"Thickness", 2.0e-3},
                     {"Permittivity", 4.0},
                     {"AutomaticEdges", true},
                     {"EdgeDistances", {0.2, 0.5}},
                     {"LocalizeEdgeEnergy", true},
                     {"EdgeRefinement", {{"Radius", 0.2}, {"ElementsPerRadius", 3}}}};
  const config::InterfaceDielectricData automatic(dielectric);
  CHECK(automatic.automatic_edges);
  CHECK(automatic.edge_attributes.empty());
  CHECK(automatic.save_local_edge_energy);

  dielectric["SaveLocalEdgeEnergy"] = false;
  const config::InterfaceDielectricData compact_output(dielectric);
  CHECK_FALSE(compact_output.save_local_edge_energy);
  dielectric.erase("SaveLocalEdgeEnergy");

  dielectric["EdgeAttributes"] = {2};
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
  dielectric.erase("EdgeAttributes");

  dielectric["EdgeExcludeAttributes"] = {3};
  const config::InterfaceDielectricData excluded(dielectric);
  CHECK(excluded.edge_exclude_attributes == std::vector<int>{3});
  dielectric.erase("EdgeExcludeAttributes");

  // EdgeExcludeSegments (decision 282 rule B5): explicit 3D segments with a default
  // tolerance of 1e-3 of the smallest matching radius.
  dielectric["EdgeExcludeSegments"] = {{0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                                       {0.0, 1.0, 0.0, 0.0, 1.0, 2.0}};
  const config::InterfaceDielectricData segments(dielectric);
  REQUIRE(segments.edge_exclude_segments.size() == 2);
  CHECK(segments.edge_exclude_segments[1] ==
        std::array<double, 6>{0.0, 1.0, 0.0, 0.0, 1.0, 2.0});
  CHECK_THAT(segments.edge_exclude_segment_tolerance,
             Catch::Matchers::WithinRel(1.0e-3 * 0.2, 1.0e-12));
  dielectric["EdgeExcludeSegmentTolerance"] = 0.01;
  CHECK(config::InterfaceDielectricData(dielectric).edge_exclude_segment_tolerance == 0.01);
  dielectric["EdgeExcludeSegmentTolerance"] = -0.01;
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
  dielectric.erase("EdgeExcludeSegmentTolerance");
  // EdgeExcludeSegmentSagittas (decision 657 (a)): one chord sagitta per listed segment,
  // finite and nonnegative; absent = every sagitta 0 (the plain rule).
  CHECK(segments.edge_exclude_segment_sagittas.empty());
  dielectric["EdgeExcludeSegmentSagittas"] = {0.0, 4.2e-3};
  CHECK(config::InterfaceDielectricData(dielectric).edge_exclude_segment_sagittas ==
        std::vector<double>{0.0, 4.2e-3});
  dielectric["EdgeExcludeSegmentSagittas"] = {4.2e-3};  // one per segment
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
  dielectric["EdgeExcludeSegmentSagittas"] = {0.0, -4.2e-3};
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
  dielectric.erase("EdgeExcludeSegmentSagittas");
  dielectric["EdgeExcludeSegments"] = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};  // zero length
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
  dielectric.erase("EdgeExcludeSegments");
  dielectric["EdgeExcludeSegmentTolerance"] = 0.01;  // without segments
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
  dielectric.erase("EdgeExcludeSegmentTolerance");
  CHECK(config::InterfaceDielectricData(dielectric).edge_exclude_segments.empty());

  dielectric["Type"] = "Default";
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
  dielectric["Type"] = "SA";

  dielectric.erase("EdgeDistances");
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
}

TEST_CASE("Config interface dielectric edge frame", "[config][Serial]")
{
  json dielectric = {{"Attributes", {1}},      {"Thickness", 2.0e-3},
                     {"Permittivity", 4.0},    {"EdgeAttributes", {2}},
                     {"EdgeDistances", {0.2}}, {"LocalizeEdgeEnergy", true}};
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));

  dielectric["EdgeFrameNormal"] = {0.0, 3.0, 4.0};
  const config::InterfaceDielectricData data(dielectric);
  REQUIRE(data.edge_frame_normal);
  CHECK_THAT((*data.edge_frame_normal)[0], Catch::Matchers::WithinAbs(0.0, 1.0e-12));
  CHECK_THAT((*data.edge_frame_normal)[1], Catch::Matchers::WithinAbs(0.6, 1.0e-12));
  CHECK_THAT((*data.edge_frame_normal)[2], Catch::Matchers::WithinAbs(0.8, 1.0e-12));

  dielectric["EdgeFrameNormal"] = {0.0, 0.0, 0.0};
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));

  dielectric.erase("EdgeAttributes");
  dielectric.erase("EdgeDistances");
  dielectric["LocalizeEdgeEnergy"] = false;
  dielectric["EdgeFrameNormal"] = {0.0, 1.0, 0.0};
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
}

TEST_CASE("Config interface dielectric quadrature ownership", "[config][Serial]")
{
  json dielectric = {{"Attributes", {1}},
                     {"Type", "SA"},
                     {"Thickness", 2.0e-3},
                     {"Permittivity", 4.0},
                     {"OwnershipDataFile", "ownership.csv"},
                     {"OwnershipGroup", 0},
                     {"OwnershipSlot", 2},
                     {"OwnershipQuadratureOrder", 20}};
  const config::InterfaceDielectricData data(dielectric);
  CHECK(data.ownership_data_file == "ownership.csv");
  CHECK(data.ownership_group == 0);
  CHECK(data.ownership_slot == 2);
  CHECK(data.ownership_quadrature_order == 20);

  for (const auto *key : {"OwnershipGroup", "OwnershipSlot", "OwnershipQuadratureOrder"})
  {
    auto missing = dielectric;
    missing.erase(key);
    CHECK_THROWS(config::InterfaceDielectricData(missing));
  }
  dielectric["OwnershipQuadratureOrder"] = 101;
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
  dielectric["OwnershipQuadratureOrder"] = 20;
  dielectric["Type"] = "Default";
  CHECK_THROWS(config::InterfaceDielectricData(dielectric));
  dielectric["Type"] = "SA";
  dielectric["Index"] = 1;
  json input = {{"Problem", {{"Type", "Electrostatic"}}},
                {"Model", {{"Mesh", "unused.msh"}}},
                {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                {"Boundaries",
                 {{"Ground", {{"Attributes", {2}}}},
                  {"Terminal", {{{"Index", 1}, {"Attributes", {1}}}}},
                  {"Postprocessing", {{"Dielectric", {dielectric}}}}}},
                {"Solver", {{"Electrostatic", {{"Save", 0}}}}}};
  IoData parsed(input, false);
  const auto concrete = IoData::ConcretizeDefaults(parsed, input);
  const auto &entry = concrete["Boundaries"]["Postprocessing"]["Dielectric"][0];
  for (const auto *key :
       {"OwnershipDataFile", "OwnershipGroup", "OwnershipSlot", "OwnershipQuadratureOrder"})
  {
    CHECK(entry.at(key) == dielectric.at(key));
  }
}

TEST_CASE("LumpedPort excitation weight", "[config][Serial]")
{
  using Catch::Matchers::ContainsSubstring;

  auto MakePort = [](const json &extra)
  {
    json port = {{"Index", 3},
                 {"Attributes", {5}},
                 {"R", 50.0},
                 {"Direction", "+Y"},
                 {"Excitation", true}};
    port.update(extra);
    return port;
  };

  SECTION("Defaults are unit amplitude and zero phase")
  {
    config::LumpedPortData excited(MakePort(json::object()));
    CHECK(excited.excitation_amplitude == 1.0);
    CHECK(excited.excitation_phase == 0.0);

    // A port without excitation also has the default weight.
    json passive = MakePort(json::object());
    passive.erase("Excitation");
    config::LumpedPortData unexcited(passive);
    CHECK(unexcited.excitation == 0);
    CHECK(unexcited.excitation_amplitude == 1.0);
    CHECK(unexcited.excitation_phase == 0.0);
  }

  SECTION("Parsed values")
  {
    config::LumpedPortData port(
        MakePort({{"ExcitationAmplitude", 2.5}, {"ExcitationPhase", -45.0}}));
    CHECK(port.excitation != 0);
    CHECK(port.excitation_amplitude == 2.5);
    CHECK(port.excitation_phase == -45.0);

    // Integer values and a multielement port with an integer excitation index.
    json elements = {{"Index", 4},
                     {"R", 75.0},
                     {"Excitation", 2},
                     {"ExcitationAmplitude", 3},
                     {"ExcitationPhase", 90},
                     {"Elements",
                      {{{"Attributes", {2}}, {"Direction", "+X"}},
                       {{"Attributes", {9}}, {"Direction", "+Y"}}}}};
    config::LumpedPortData multi(elements);
    CHECK(multi.excitation == 2);
    CHECK(multi.excitation_amplitude == 3.0);
    CHECK(multi.excitation_phase == 90.0);
  }

  SECTION("Amplitude must be positive and finite")
  {
    CHECK_THROWS_WITH(config::LumpedPortData(MakePort({{"ExcitationAmplitude", 0.0}})),
                      ContainsSubstring("ExcitationAmplitude"));
    CHECK_THROWS(config::LumpedPortData(MakePort({{"ExcitationAmplitude", -1.0}})));
    CHECK_THROWS(config::LumpedPortData(
        MakePort({{"ExcitationAmplitude", std::numeric_limits<double>::infinity()}})));
    CHECK_THROWS(config::LumpedPortData(
        MakePort({{"ExcitationAmplitude", std::numeric_limits<double>::quiet_NaN()}})));
  }

  SECTION("Phase must be finite")
  {
    CHECK_THROWS_WITH(config::LumpedPortData(MakePort(
                          {{"ExcitationPhase", std::numeric_limits<double>::infinity()}})),
                      ContainsSubstring("ExcitationPhase"));
  }

  SECTION("Weight keys require an excitation")
  {
    for (const auto &key : {"ExcitationAmplitude", "ExcitationPhase"})
    {
      CAPTURE(key);

      // No "Excitation" key.
      json missing = MakePort({{key, 1.0}});
      missing.erase("Excitation");
      CHECK_THROWS_WITH(config::LumpedPortData(missing),
                        ContainsSubstring("require a nonzero \"Excitation\""));

      // Explicitly not excited, as a boolean and as an index.
      CHECK_THROWS(config::LumpedPortData(MakePort({{key, 1.0}, {"Excitation", false}})));
      CHECK_THROWS(config::LumpedPortData(MakePort({{key, 1.0}, {"Excitation", 0}})));
    }
  }

  SECTION("Transient simulations reject a complex weight")
  {
    auto MakeConfig = [&](const json &extra)
    {
      return json{
          {"Problem", {{"Type", "Transient"}, {"Output", "test_output"}}},
          {"Model", {{"Mesh", "test.msh"}}},
          {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
          {"Boundaries", {{"LumpedPort", json::array({MakePort(extra)})}}},
          {"Solver",
           {{"Transient",
             {{"Excitation", "Sinusoidal"}, {"MaxTime", 1.0}, {"TimeStep", 0.01}}}}}};
    };

    CHECK_NOTHROW(IoData(MakeConfig(json::object()), false));
    CHECK_THROWS_WITH(IoData(MakeConfig({{"ExcitationPhase", 30.0}}), false),
                      ContainsSubstring("ExcitationPhase"));
    CHECK_THROWS_WITH(IoData(MakeConfig({{"ExcitationPhase", 30.0}}), false),
                      ContainsSubstring("multiple of 180"));
    CHECK_THROWS(IoData(MakeConfig({{"ExcitationPhase", -90.0}}), false));
    CHECK_THROWS(IoData(MakeConfig({{"ExcitationPhase", 270.0}}), false));

    // A multiple of 180 degrees gives a real weight, and the amplitude is unrestricted.
    CHECK_NOTHROW(IoData(MakeConfig({{"ExcitationPhase", 180.0}}), false));
    CHECK_NOTHROW(IoData(MakeConfig({{"ExcitationPhase", -180.0}}), false));
    CHECK_NOTHROW(IoData(MakeConfig({{"ExcitationPhase", 360.0}}), false));
    CHECK_NOTHROW(IoData(MakeConfig({{"ExcitationAmplitude", 2.0}}), false));

    // Driven simulations accept any phase.
    json driven = MakeConfig({{"ExcitationPhase", 30.0}});
    driven["Problem"]["Type"] = "Driven";
    driven["Solver"] = {
        {"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.1}}}};
    CHECK_NOTHROW(IoData(driven, false));
  }

  SECTION("Concretized configuration keeps the weight and re-parses")
  {
    json config = {
        {"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries",
         {{"LumpedPort",
           json::array({MakePort({{"ExcitationAmplitude", 2.0}, {"ExcitationPhase", 30.0}}),
                        json{{"Index", 4},
                             {"Attributes", {6}},
                             {"R", 50.0},
                             {"Direction", "-Y"}}})}}},
        {"Solver", {{"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.1}}}}}};

    IoData iodata1(config, false);
    config = IoData::ConcretizeDefaults(iodata1, config);

    // The weight is emitted for the excited port only, because the parser rejects it on a
    // port without an excitation.
    const auto &j_ports = config["Boundaries"]["LumpedPort"];
    CHECK(j_ports[0]["ExcitationAmplitude"].get<double>() == 2.0);
    CHECK(j_ports[0]["ExcitationPhase"].get<double>() == 30.0);
    CHECK_FALSE(j_ports[1].contains("ExcitationAmplitude"));
    CHECK_FALSE(j_ports[1].contains("ExcitationPhase"));

    std::string err = ValidateConfig(config);
    INFO("schema validation error: " << err);
    CHECK(err.empty());

    IoData iodata2(config, false);
    CHECK(iodata2.boundaries.lumpedport.at(3).excitation_amplitude == 2.0);
    CHECK(iodata2.boundaries.lumpedport.at(3).excitation_phase == 30.0);
    CHECK(iodata2.boundaries.lumpedport.at(4).excitation_amplitude == 1.0);
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

  SECTION("Save and PROM samples which differ from the sample frequencies by rounding")
  {
    // Same frequency with different rounding, as from two overlapping ranges.
    const double f_sweep = 7.85, f_prom = std::nextafter(7.85, 8.0);
    REQUIRE(f_prom != f_sweep);
    json driven = {{"Samples",
                    {{{"Freq", {6.0, 7.0, f_sweep, 8.0}}},
                     {{"Freq", {f_prom}}, {"SaveStep", 1}, {"AddToPROM", true}},
                     {{"Freq", {std::nextafter(8.0, 9.0)}}, {"AddToPROM", true}}}},
                   {"AdaptiveTol", 1e-3}};
    config::DrivenSolverData driven_solver(driven);

    CHECK(driven_solver.sample_f == std::vector{6.0, 7.0, f_sweep, 8.0});
    CHECK(driven_solver.save_indices == std::vector<size_t>{2});
    CHECK(driven_solver.prom_indices == std::vector<size_t>{0, 3, 2});
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

    // Nonpositive deprecated FreqStep
    json zero_step = {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.0}};
    CHECK_THROWS_WITH(config::DrivenSolverData(zero_step),
                      Catch::Matchers::ContainsSubstring("Invalid range"));
    json negative_step = {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", -0.5}};
    CHECK_THROWS(config::DrivenSolverData(negative_step));

    // "MaxFreq" below "MinFreq" with a step would sample only "MinFreq".
    json reversed_step = {{"MinFreq", 2.0}, {"MaxFreq", 1.0}, {"FreqStep", 0.5}};
    CHECK_THROWS_WITH(config::DrivenSolverData(reversed_step),
                      ContainsSubstring("Invalid range"));
    json reversed_linear = {
        {"Samples", {{{"MinFreq", 2.0}, {"MaxFreq", 1.0}, {"FreqStep", 0.5}}}}};
    CHECK_THROWS_WITH(config::DrivenSolverData(reversed_linear),
                      ContainsSubstring("Invalid range"));

    // A step too small for the number of steps to fit in an int.
    json tiny_step = {{"MinFreq", 1.0}, {"MaxFreq", 10.0}, {"FreqStep", 1.0e-9}};
    CHECK_THROWS_WITH(config::DrivenSolverData(tiny_step), ContainsSubstring("too small"));

    // Circuit synthesis needs a positive "AdaptiveTol" and two distinct samples, but only
    // for a Driven problem, so the parser accepts both.
    json synthesis_no_tol = {{"Samples", {{{"Freq", {1.0, 2.0}}}}},
                             {"AdaptiveCircuitSynthesis", true}};
    CHECK_NOTHROW(config::DrivenSolverData(synthesis_no_tol));
    json synthesis_single = {{"Samples", {{{"Freq", {3.0}}}}},
                             {"AdaptiveTol", 1e-3},
                             {"AdaptiveCircuitSynthesis", true}};
    CHECK_NOTHROW(config::DrivenSolverData(synthesis_single));
  }

  SECTION("Circuit synthesis requires two distinct samples")
  {
    auto Check = [](ProblemType type, const json &driven)
    {
      IoData iodata{Units(1.0, 1.0)};
      iodata.problem.type = type;
      iodata.solver.driven = config::DrivenSolverData(driven);
      iodata.CheckConfiguration();
    };
    json synthesis = {{"AdaptiveTol", 1e-3}, {"AdaptiveCircuitSynthesis", true}};
    json single = synthesis, near_duplicate = synthesis;
    single["Samples"] = {{{"Freq", {3.0}}}};
    // Within the sample deduplication tolerance, so a single distinct sample.
    near_duplicate["Samples"] = {{{"Freq", {3.0, 3.0 + 1.0e-10}}}};
    CHECK_THROWS_WITH(Check(ProblemType::DRIVEN, single),
                      ContainsSubstring("two distinct"));
    CHECK_THROWS_WITH(Check(ProblemType::DRIVEN, near_duplicate),
                      ContainsSubstring("two distinct"));
    // An unused "Driven" block of another problem type is not checked.
    CHECK_NOTHROW(Check(ProblemType::EIGENMODE, single));
  }

  SECTION("Single sample ranges")
  {
    // A single linear or log sample is the interval start (numpy.linspace convention).
    json linear_single = {
        {"Samples", {{{"MinFreq", 2.0}, {"MaxFreq", 5.0}, {"NSample", 1}}}}};
    config::DrivenSolverData linear_solver(linear_single);
    REQUIRE(linear_solver.sample_f.size() == 1);
    CHECK_THAT(linear_solver.sample_f[0], WithinAbs(2.0, delta_eps));

    json log_single = {
        {"Samples",
         {{{"Type", "Log"}, {"MinFreq", 0.1}, {"MaxFreq", 10.0}, {"NSample", 1}}}}};
    config::DrivenSolverData log_solver(log_single);
    REQUIRE(log_solver.sample_f.size() == 1);
    CHECK_THAT(log_solver.sample_f[0], WithinAbs(0.1, delta_eps));

    // A single-sample adaptive sweep without circuit synthesis parses; DrivenSolver
    // reverts it to a uniform sweep at runtime.
    json adaptive_single = {
        {"Samples", {{{"MinFreq", 2.0}, {"MaxFreq", 5.0}, {"NSample", 1}}}},
        {"AdaptiveTol", 1e-3}};
    CHECK_NOTHROW(config::DrivenSolverData(adaptive_single));

    // Deprecated top-level interface with a valid positive step.
    json legacy = {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.5}};
    config::DrivenSolverData legacy_solver(legacy);
    auto sample_f = std::vector{1.0, 1.5, 2.0};
    REQUIRE(legacy_solver.sample_f.size() == sample_f.size());
    for (size_t i = 0; i < sample_f.size(); ++i)
    {
      CHECK_THAT(legacy_solver.sample_f[i], WithinAbs(sample_f[i], delta_eps));
    }
  }
}

TEST_CASE("Config Eigenmode saved modes", "[config][Serial]")
{
  CHECK(config::EigenSolverData(json{{"Target", 1.0}, {"N", 3}, {"Save", 2}}).n_post == 2);
  CHECK(config::EigenSolverData(json{{"Target", 1.0}, {"N", 2}, {"Save", 3}}).n_post == 2);
  CHECK(config::EigenSolverData(json{{"Target", 1.0}, {"N", 2}, {"Save", -1}}).n_post == 0);
}

TEST_CASE("Config Eigenmode iteration limits", "[config][Serial]")
{
  // The quasi-Newton iterations of the hybrid nonlinear eigenvalue solver have their own
  // limit, independent of the (large by default) eigenvalue solver iteration limit.
  const config::EigenSolverData defaults(json{{"Target", 1.0}});
  CHECK(defaults.max_it == 1'000'000);
  CHECK(defaults.nonlinear_max_it == 100);

  const config::EigenSolverData linear(json{{"Target", 1.0}, {"MaxIts", 200}});
  CHECK(linear.max_it == 200);
  CHECK(linear.nonlinear_max_it == 100);

  const config::EigenSolverData both(
      json{{"Target", 1.0}, {"MaxIts", 200}, {"NonlinearMaxIts", 30}});
  CHECK(both.max_it == 200);
  CHECK(both.nonlinear_max_it == 30);

  const config::EigenSolverData nonlinear(json{{"Target", 1.0}, {"NonlinearMaxIts", 500}});
  CHECK(nonlinear.max_it == 1'000'000);
  CHECK(nonlinear.nonlinear_max_it == 500);
}

TEST_CASE("Config Magnetostatic InactivePorts", "[config][Serial]")
{
  auto MixedConfig = []()
  {
    return json{{"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
                {"Model", {{"Mesh", "test.msh"}}},
                {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                {"Boundaries",
                 {{"SurfaceCurrent",
                   {{{"Attributes", {4}},
                     {"Index", 1},
                     {"Direction", "+X"},
                     {"Aperture", {{"Attributes", {5}}, {"Direction", "+Z"}}}}}},
                  {"FluxLoop",
                   {{{"Index", 2},
                     {"FilmAttributes", {6}},
                     {"HoleAttributes", {7}},
                     {"FluxAmounts", {1.0}},
                     {"Direction", "+Z"}}}}}},
                {"Solver", json::object()}};
  };

  SECTION("Defaults to Open when unspecified")
  {
    config::MagnetostaticSolverData magnetostatic(json::object());
    CHECK(magnetostatic.inactive_port_mode == InactivePortMode::OPEN);
  }

  SECTION("Global InactivePorts parses Open and Short")
  {
    config::MagnetostaticSolverData open(json{{"InactivePorts", "Open"}});
    CHECK(open.inactive_port_mode == InactivePortMode::OPEN);
    config::MagnetostaticSolverData shorted(json{{"InactivePorts", "Short"}});
    CHECK(shorted.inactive_port_mode == InactivePortMode::SHORT);
  }

  SECTION("Per-port InactiveMode overrides the global default")
  {
    // Global default is Short; port 1 has no override (inherits Short), port 2 overrides
    // to Open. This mirrors the value_or(global) resolution done in the solver.
    config::MagnetostaticSolverData magnetostatic(json{{"InactivePorts", "Short"}});
    const InactivePortMode global_mode = magnetostatic.inactive_port_mode;

    json boundaries = {{"SurfaceCurrent",
                        {{{"Attributes", {4}}, {"Index", 1}, {"Direction", "+X"}},
                         {{"Attributes", {5}},
                          {"Index", 2},
                          {"Direction", "+X"},
                          {"InactiveMode", "Open"}}}}};
    config::BoundaryData boundary_data(boundaries);

    auto resolve = [&](int idx)
    { return boundary_data.current.at(idx).inactive_port_mode.value_or(global_mode); };
    CHECK(resolve(1) == InactivePortMode::SHORT);  // inherits global
    CHECK(resolve(2) == InactivePortMode::OPEN);   // per-port override wins
  }

  SECTION("Mixed extraction requires complete apertures")
  {
    auto config = MixedConfig();
    CHECK_NOTHROW(IoData(config, false));

    config["Boundaries"]["SurfaceCurrent"][0].erase("Aperture");
    CHECK_THROWS(IoData(config, false));
  }

  SECTION("Mixed extraction requires effectively Open current ports")
  {
    auto config = MixedConfig();
    config["Solver"]["Magnetostatic"] = {{"InactivePorts", "Short"}};
    CHECK_THROWS(IoData(config, false));

    config["Boundaries"]["SurfaceCurrent"][0]["InactiveMode"] = "Open";
    CHECK_NOTHROW(IoData(config, false));

    config["Boundaries"]["SurfaceCurrent"][0]["InactiveMode"] = "Short";
    CHECK_THROWS(IoData(config, false));
  }

  SECTION("Aperture is rejected outside mixed extraction")
  {
    auto config = MixedConfig();
    config["Boundaries"].erase("FluxLoop");
    CHECK_THROWS(IoData(config, false));
  }

  SECTION("Aperture direction must be finite, nonzero, and Cartesian")
  {
    CHECK_THROWS(config::SurfaceCurrentApertureData(
        json{{"Attributes", {1}}, {"Direction", {0.0, 0.0, 0.0}}}));
    CHECK_THROWS(
        config::SurfaceCurrentApertureData(json{{"Attributes", {1}}, {"Direction", "+R"}}));

    config::SurfaceCurrentApertureData aperture(
        json{{"Attributes", {1}}, {"Direction", {0.0, 3.0, 4.0}}});
    CHECK(aperture.direction[0] == 0.0);
    CHECK(aperture.direction[1] == Catch::Approx(0.6));
    CHECK(aperture.direction[2] == Catch::Approx(0.8));
  }
}

TEST_CASE("Config Magnetostatic auxiliary smoother", "[config][Serial]")
{
  using Catch::Matchers::ContainsSubstring;
  auto Check = [](int mg_smooth_aux, double london_pc_shift)
  {
    IoData iodata{Units(1.0, 1.0)};
    iodata.problem.type = ProblemType::MAGNETOSTATIC;
    iodata.solver.linear.mg_smooth_aux = mg_smooth_aux;
    iodata.solver.linear.london_pc_shift = london_pc_shift;
    iodata.CheckConfiguration();
  };
  CHECK_THROWS_WITH(Check(1, 0.0), ContainsSubstring("\"LondonPCShift\""));
  CHECK_NOTHROW(Check(1, 0.1));
  CHECK_NOTHROW(Check(-1, 0.0));
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
    CHECK(data.thetaphis[1].first == Catch::Approx(std::numbers::pi / 2).margin(delta_eps));
    CHECK(data.thetaphis[1].second == Catch::Approx(std::numbers::pi).margin(delta_eps));
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
        if (std::abs(point.first - std::numbers::pi) < delta_eps)
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
             std::abs(point.second - std::numbers::pi) < delta_eps) &&
            std::abs(point.first) > delta_eps &&
            std::abs(point.first - std::numbers::pi) > delta_eps)
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
        if (std::abs(point.first - std::numbers::pi / 2.0) < delta_eps)
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
    CHECK(config["Solver"]["Electrostatic"]["ResponseMatrix"].get<bool>() == false);
  }

  SECTION("Electrostatic response matrix option")
  {
    const config::ElectrostaticSolverData defaults(json::object());
    CHECK_FALSE(defaults.response_matrix);
    CHECK_FALSE(defaults.aggregate_response_matrix);

    const json electrostatic = {{"ResponseMatrix", true},
                                {"AggregateResponseMatrix", true}};
    const config::ElectrostaticSolverData enabled(electrostatic);
    CHECK(enabled.response_matrix);
    CHECK(enabled.aggregate_response_matrix);
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
    CHECK(config["Solver"]["Magnetostatic"]["InactivePorts"] == "Open");
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

  SECTION("ComplexCoarseSolve is disabled with AMS and BoomerAMG")
  {
    // AMS and BoomerAMG cannot be given the complex-valued coarse matrix as a real-valued
    // system of twice the size: an explicit true falls back to false, also when the type
    // is the resolved default, and the resolved configuration says so. Other types keep
    // it, and so does a boundary mode problem without multigrid, which replaces AMS and
    // BoomerAMG by a sparse direct solver.
    auto Resolve = [](const char *problem, json solver)
    {
      json config = {{"Problem", {{"Type", problem}, {"Output", "test_output"}}},
                     {"Model", {{"Mesh", "test.msh"}}},
                     {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                     {"Boundaries", json::object()},
                     {"Solver", std::move(solver)}};
      IoData iodata(config, false);
      config = IoData::ConcretizeDefaults(iodata, config);
      CHECK(config["Solver"]["Linear"]["ComplexCoarseSolve"].get<bool>() ==
            iodata.solver.linear.complex_coarse_solve);
      return iodata.solver.linear;
    };
    auto Linear = [](const char *type, int mg_max_levels = -1)
    {
      json linear = {{"Type", type}, {"ComplexCoarseSolve", true}};
      if (mg_max_levels > 0)
      {
        linear["MGMaxLevels"] = mg_max_levels;
      }
      return linear;
    };
    const json eigenmode = {{"Target", 1.0}}, boundary_mode = {{"Freq", 10.0}};
    for (const char *type : {"AMS", "BoomerAMG"})
    {
      INFO("type " << type);
      CHECK_FALSE(Resolve("Eigenmode", {{"Eigenmode", eigenmode}, {"Linear", Linear(type)}})
                      .complex_coarse_solve);
      CHECK_FALSE(
          Resolve("Eigenmode", {{"Eigenmode", eigenmode}, {"Linear", Linear(type, 1)}})
              .complex_coarse_solve);
      CHECK_FALSE(Resolve("BoundaryMode",
                          {{"BoundaryMode", boundary_mode}, {"Linear", Linear(type, 2)}})
                      .complex_coarse_solve);
      CHECK(Resolve("BoundaryMode",
                    {{"BoundaryMode", boundary_mode}, {"Linear", Linear(type)}})
                .complex_coarse_solve);
    }
    {
      // A real-valued problem does not use the option: it is left as it is, without a
      // warning.
      const auto linear = Resolve("Magnetostatic", {{"Linear", Linear("Default")}});
      CHECK(linear.type == LinearSolver::AMS);
      CHECK(linear.complex_coarse_solve);
    }
    for (const bool pc_mat_real : {false, true})
    {
      // The wave ports of a driven problem keep the inherited value for their mode solver,
      // which uses a sparse direct solver, unless PCMatReal overrides it. An explicit
      // per-port value is kept as it is.
      INFO("PCMatReal " << pc_mat_real);
      auto linear = Linear("AMS");
      linear["PCMatReal"] = pc_mat_real;
      json config = {
          {"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
          {"Model", {{"Mesh", "test.msh"}}},
          {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
          {"Boundaries",
           {{"WavePort",
             {{{"Index", 1}, {"Attributes", {2}}, {"Excitation", true}},
              {{"Index", 2}, {"Attributes", {3}}, {"ComplexCoarseSolve", false}}}}}},
          {"Solver",
           {{"Driven",
             {{"Samples", {{{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 0.5}}}}}},
            {"Linear", linear}}}};
      IoData iodata(config, false);
      CHECK_FALSE(iodata.solver.linear.complex_coarse_solve);
      const auto &port1 = iodata.boundaries.waveport.at(1).complex_coarse_solve;
      const auto &port2 = iodata.boundaries.waveport.at(2).complex_coarse_solve;
      CHECK(port1.has_value() == !pc_mat_real);
      CHECK(port1.value_or(true));
      CHECK(port2.has_value());
      CHECK_FALSE(port2.value_or(true));

      // The resolved configuration has the values which are used, so that it gives the
      // same simulation.
      config = IoData::ConcretizeDefaults(iodata, config);
      const auto &j_ports = config["Boundaries"]["WavePort"];
      CHECK_FALSE(config["Solver"]["Linear"]["ComplexCoarseSolve"].get<bool>());
      CHECK(j_ports[0].contains("ComplexCoarseSolve") == !pc_mat_real);
      CHECK(j_ports[0].value("ComplexCoarseSolve", true));
      CHECK_FALSE(j_ports[1].at("ComplexCoarseSolve").get<bool>());
      IoData iodata2(config, false);
      CHECK_FALSE(iodata2.solver.linear.complex_coarse_solve);
      CHECK(iodata2.boundaries.waveport.at(1).complex_coarse_solve == port1);
      CHECK(iodata2.boundaries.waveport.at(2).complex_coarse_solve == port2);
    }
    {
      // The default of a frequency domain problem is a sparse direct solver if there is
      // one in the build, AMS otherwise.
      const auto linear =
          Resolve("Eigenmode", {{"Eigenmode", eigenmode}, {"Linear", Linear("Default")}});
      CHECK(linear.complex_coarse_solve == (linear.type != LinearSolver::AMS));
    }
    CHECK(Resolve("Eigenmode", {{"Eigenmode", eigenmode}, {"Linear", Linear("Jacobi")}})
              .complex_coarse_solve);
    CHECK_FALSE(
        Resolve("Eigenmode", {{"Eigenmode", eigenmode}, {"Linear", {{"Type", "AMS"}}}})
            .complex_coarse_solve);
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
          {"Refinement",
           {{"MaxIts", 4},
            {"MaxSize", 1000},
            {"MaxSizePredicted", false},
            {"UpdateFraction", 0.5}}}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries",
         {{"Absorbing", {{"Attributes", {2}}, {"Order", 2}}},
          {"Conductivity", {{{"Attributes", {3}}, {"Conductivity", 5.8e7}}}},
          {"Impedance",
           {{{"Attributes", {4}}, {"Rs", 50.0}, {"Ls", 1.0e-9}, {"Cs", 1.0e-12}}}},
          {"RationalImpedance",
           {{{"Attributes", {8}},
             {"Numerator", {5.0e-8, 0.0}},
             {"Denominator", {5.0e-20, 1.0e-9, 50.0}}}}},
          {"LumpedPort",
           {{{"Index", 1},
             {"Attributes", {5}},
             {"Direction", "+X"},
             {"R", 50.0},
             {"L", 1.0e-9},
             {"Excitation", true},
             {"Active", false}}}},
          {"FloquetPort", {{{"Index", 2}, {"Attributes", {12}}}}},
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
               {"LossTan", 0.002},
               {"EdgeAttributes", {5}},
               {"EdgeExcludeAttributes", {6, 7}},
               {"EdgeExcludeSegments", {{0.0, 0.0, 0.0, 1.0e-6, 0.0, 0.0}}},
               {"EdgeExcludeSegmentTolerance", 1.0e-9},
               {"EdgeExcludeSegmentSagittas", {4.2e-9}},
               {"EdgeDistances", {1.0e-6, 2.0e-6}},
               {"EdgeDistanceSmoothing", 0.2},
               {"LocalizeEdgeEnergy", true},
               {"EdgeFrameNormal", {0.0, 0.0, 1.0}},
               {"EdgeRefinement",
                {{"Radius", 1.0e-6},
                 {"ElementsPerRadius", 3},
                 {"OuterRadiusFactor", 2.5},
                 {"CoreIndicatorWeight", 0.1}}},
               {"FluxRecovery", true}}}}}}}},
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
    CHECK(r2.max_size_predicted == r1.max_size_predicted);
    CHECK_FALSE(r2.max_size_predicted);
    CHECK(r2.nonconformal == r1.nonconformal);
    CHECK(r2.max_nc_levels == r1.max_nc_levels);
    CHECK(r2.update_fraction == r1.update_fraction);
    CHECK(r2.maximum_imbalance == r1.maximum_imbalance);
    CHECK(r2.save_adapt_iterations == r1.save_adapt_iterations);
    CHECK(r2.save_adapt_mesh == r1.save_adapt_mesh);
    CHECK(r2.uniform_ref_levels == r1.uniform_ref_levels);
    CHECK(r2.ser_uniform_ref_levels == r1.ser_uniform_ref_levels);
    CHECK(iodata2.boundaries.farfield.order == iodata1.boundaries.farfield.order);
    REQUIRE(iodata1.boundaries.postpro.dielectric.count(1) == 1);
    REQUIRE(iodata2.boundaries.postpro.dielectric.count(1) == 1);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).edge_exclude_attributes ==
          iodata1.boundaries.postpro.dielectric.at(1).edge_exclude_attributes);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).edge_exclude_segments ==
          iodata1.boundaries.postpro.dielectric.at(1).edge_exclude_segments);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).edge_exclude_segment_tolerance ==
          iodata1.boundaries.postpro.dielectric.at(1).edge_exclude_segment_tolerance);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).edge_exclude_segment_sagittas ==
          iodata1.boundaries.postpro.dielectric.at(1).edge_exclude_segment_sagittas);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).automatic_edges ==
          iodata1.boundaries.postpro.dielectric.at(1).automatic_edges);
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
    CHECK(iodata2.boundaries.impedance[0].crack == iodata1.boundaries.impedance[0].crack);
    CHECK(iodata2.boundaries.conductivity[0].crack ==
          iodata1.boundaries.conductivity[0].crack);
    REQUIRE(iodata1.boundaries.rational_impedance.size() == 1);
    REQUIRE(iodata2.boundaries.rational_impedance.size() == 1);
    CHECK(iodata2.boundaries.rational_impedance[0].attributes ==
          iodata1.boundaries.rational_impedance[0].attributes);
    CHECK(iodata2.boundaries.rational_impedance[0].num ==
          iodata1.boundaries.rational_impedance[0].num);
    CHECK(iodata2.boundaries.rational_impedance[0].den ==
          iodata1.boundaries.rational_impedance[0].den);
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
    REQUIRE(iodata2.boundaries.floquetport.count(2) == 1);
    const auto &fp1 = iodata1.boundaries.floquetport.at(2);
    const auto &fp2 = iodata2.boundaries.floquetport.at(2);
    CHECK(fp2.excitation == fp1.excitation);
    CHECK(fp2.inc_polarization == fp1.inc_polarization);
    CHECK(fp2.max_order == fp1.max_order);
    CHECK(iodata2.boundaries.periodic.wave_vector ==
          iodata1.boundaries.periodic.wave_vector);
    CHECK(iodata2.boundaries.periodic.floquet_reference_freq ==
          iodata1.boundaries.periodic.floquet_reference_freq);
    REQUIRE(iodata2.boundaries.postpro.flux.count(1) == 1);
    CHECK(iodata2.boundaries.postpro.flux.at(1).two_sided ==
          iodata1.boundaries.postpro.flux.at(1).two_sided);
    REQUIRE(iodata2.boundaries.postpro.dielectric.count(1) == 1);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).type ==
          iodata1.boundaries.postpro.dielectric.at(1).type);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).tandelta ==
          iodata1.boundaries.postpro.dielectric.at(1).tandelta);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).edge_attributes ==
          iodata1.boundaries.postpro.dielectric.at(1).edge_attributes);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).edge_distances ==
          iodata1.boundaries.postpro.dielectric.at(1).edge_distances);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).edge_distance_smoothing ==
          iodata1.boundaries.postpro.dielectric.at(1).edge_distance_smoothing);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).localize_edge_energy ==
          iodata1.boundaries.postpro.dielectric.at(1).localize_edge_energy);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).save_local_edge_energy ==
          iodata1.boundaries.postpro.dielectric.at(1).save_local_edge_energy);
    REQUIRE(iodata2.boundaries.postpro.dielectric.at(1).edge_refinement);
    REQUIRE(iodata1.boundaries.postpro.dielectric.at(1).edge_refinement);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).edge_refinement->radius ==
          iodata1.boundaries.postpro.dielectric.at(1).edge_refinement->radius);
    CHECK(
        iodata2.boundaries.postpro.dielectric.at(1).edge_refinement->elements_per_radius ==
        iodata1.boundaries.postpro.dielectric.at(1).edge_refinement->elements_per_radius);
    CHECK(
        iodata2.boundaries.postpro.dielectric.at(1).edge_refinement->outer_radius_factor ==
        iodata1.boundaries.postpro.dielectric.at(1).edge_refinement->outer_radius_factor);
    CHECK(
        iodata2.boundaries.postpro.dielectric.at(1)
            .edge_refinement->core_indicator_weight ==
        iodata1.boundaries.postpro.dielectric.at(1).edge_refinement->core_indicator_weight);
    CHECK(iodata2.boundaries.postpro.dielectric.at(1).flux_recovery ==
          iodata1.boundaries.postpro.dielectric.at(1).flux_recovery);

    // Coverage gates. Each schema scope this fixture exercises is checked here so a
    // future schema addition without matching Concretize emission fails this section.
    // CrackInternalBoundaryElements is deprecated and only emitted when specified.
    auto model_gaps = SchemaCoverageGaps("/properties/Model", config["Model"],
                                         /*skip=*/{"Lc", "CrackInternalBoundaryElements"});
    INFO("Model missing keys: " << json(model_gaps).dump());
    CHECK(model_gaps.empty());

    // Boxes/Spheres are opt-in arrays of explicit refinement regions; absence means
    // no per-region refinement and there is no scalar default to emit.
    auto ref_gaps = SchemaCoverageGaps("/properties/Model/properties/Refinement",
                                       config["Model"]["Refinement"],
                                       /*skip=*/{"Boxes", "Spheres"});
    INFO("Model.Refinement missing keys: " << json(ref_gaps).dump());
    CHECK(ref_gaps.empty());

    auto abs_gaps = SchemaCoverageGaps("/properties/Boundaries/properties/Absorbing",
                                       config["Boundaries"]["Absorbing"]);
    INFO("Boundaries.Absorbing missing keys: " << json(abs_gaps).dump());
    CHECK(abs_gaps.empty());

    auto cond_gaps =
        SchemaCoverageGaps("/properties/Boundaries/properties/Conductivity/items",
                           config["Boundaries"]["Conductivity"][0]);
    INFO("Boundaries.Conductivity[] missing keys: " << json(cond_gaps).dump());
    CHECK(cond_gaps.empty());

    auto imp_gaps = SchemaCoverageGaps("/properties/Boundaries/properties/Impedance/items",
                                       config["Boundaries"]["Impedance"][0]);
    INFO("Boundaries.Impedance[] missing keys: " << json(imp_gaps).dump());
    CHECK(imp_gaps.empty());

    auto rz_gaps =
        SchemaCoverageGaps("/properties/Boundaries/properties/RationalImpedance/items",
                           config["Boundaries"]["RationalImpedance"][0]);
    INFO("Boundaries.RationalImpedance[] missing keys: " << json(rz_gaps).dump());
    CHECK(rz_gaps.empty());

    // LumpedPort: Direction/CoordinateSystem/Elements are opt-in alternatives to
    // Attributes for declaring the port geometry, and Length (decision 929) an opt-in
    // geometry override without a default; concretize does not synthesize them.
    auto lp_gaps = SchemaCoverageGaps("/properties/Boundaries/properties/LumpedPort/items",
                                      config["Boundaries"]["LumpedPort"][0],
                                      /*skip=*/{"CoordinateSystem", "Elements", "Length"});
    INFO("Boundaries.LumpedPort[] missing keys: " << json(lp_gaps).dump());
    CHECK(lp_gaps.empty());

    auto fp_gaps = SchemaCoverageGaps("/properties/Boundaries/properties/FloquetPort/items",
                                      config["Boundaries"]["FloquetPort"][0],
                                      /*skip=*/{"MaxOrder"});
    INFO("Boundaries.FloquetPort[] missing keys: " << json(fp_gaps).dump());
    CHECK(fp_gaps.empty());

    auto per_gaps = SchemaCoverageGaps("/properties/Boundaries/properties/Periodic",
                                       config["Boundaries"]["Periodic"]);
    INFO("Boundaries.Periodic missing keys: " << json(per_gaps).dump());
    CHECK(per_gaps.empty());

    // Ownership is an opt-in, all-fields-together partition; inactive defaults
    // must not be synthesized into an otherwise unpartitioned interface. Region is the
    // opt-in quadrature-level filter (absent = unfiltered).
    auto dielectric_gaps = SchemaCoverageGaps(
        "/$defs/Dielectric", config["Boundaries"]["Postprocessing"]["Dielectric"][0],
        /*skip=*/
        {"OwnershipDataFile", "OwnershipGroup", "OwnershipSlot", "OwnershipQuadratureOrder",
         "Region"});
    INFO("Boundaries.Postprocessing.Dielectric[] missing keys: "
         << json(dielectric_gaps).dump());
    CHECK(dielectric_gaps.empty());

    auto edge_refinement_gaps = SchemaCoverageGaps(
        "/$defs/EdgeRefinement",
        config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeRefinement"]);
    INFO("Boundaries.Postprocessing.Dielectric[].EdgeRefinement missing keys: "
         << json(edge_refinement_gaps).dump());
    CHECK(edge_refinement_gaps.empty());
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
    CHECK(iodata2.solver.eigenmode.nonlinear_max_it ==
          iodata1.solver.eigenmode.nonlinear_max_it);

    // Coverage gate: every optional schema property under Solver.Eigenmode must be
    // emitted by ConcretizeDefaults. A schema addition without a Concretize update
    // surfaces here.
    auto eigen_gaps = SchemaCoverageGaps("/properties/Solver/properties/Eigenmode",
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
    auto trans_gaps = SchemaCoverageGaps("/properties/Solver/properties/Transient",
                                         config["Solver"]["Transient"],
                                         /*skip=*/{"Order", "RelTol", "AbsTol"});
    INFO("Solver.Transient missing keys: " << json(trans_gaps).dump());
    CHECK(trans_gaps.empty());
  }

  SECTION("Round-trip: Magnetostatic and FluxLoop defaulted fields are reproducible")
  {
    json config = {{"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries",
                    {{"FluxLoop",
                      {{{"Index", 3},
                        {"FilmAttributes", {8}},
                        {"HoleAttributes", {9}},
                        {"FluxAmounts", {1.0}},
                        {"Direction", "+Z"}}}}}},
                   {"Solver", json::object()}};

    IoData iodata1(config, false);
    config = IoData::ConcretizeDefaults(iodata1, config);

    std::string err = ValidateConfig(config);
    INFO("schema validation error: " << err);
    CHECK(err.empty());

    // The default excitation is the loop index, written out as such.
    CHECK(config["Boundaries"]["FluxLoop"][0]["Excitation"] == 3);
    IoData iodata2(config, false);
    CHECK(iodata2.boundaries.fluxloop.at(3).excitation ==
          iodata1.boundaries.fluxloop.at(3).excitation);
    CHECK(iodata2.boundaries.fluxloop.at(3).pec_lperp ==
          iodata1.boundaries.fluxloop.at(3).pec_lperp);

    auto ms_gaps = SchemaCoverageGaps("/properties/Solver/properties/Magnetostatic",
                                      config["Solver"]["Magnetostatic"], /*skip=*/{});
    INFO("Solver.Magnetostatic missing keys: " << json(ms_gaps).dump());
    CHECK(ms_gaps.empty());
    auto loop_gaps = SchemaCoverageGaps("/properties/Boundaries/properties/FluxLoop/items",
                                        config["Boundaries"]["FluxLoop"][0], /*skip=*/{});
    INFO("Boundaries.FluxLoop[] missing keys: " << json(loop_gaps).dump());
    CHECK(loop_gaps.empty());
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
    // mesh); leaving it absent means "use the parent mesh as-is". Target is resolved
    // from the mesh/materials at solve time, so omission remains the auto-selection path.
    auto bm_gaps = SchemaCoverageGaps("/properties/Solver/properties/BoundaryMode",
                                      config["Solver"]["BoundaryMode"],
                                      /*skip=*/{"Attributes", "Target"});
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
    // An omitted per-port ComplexCoarseSolve must stay omitted (see coverage gate below).
    CHECK_FALSE(w2.complex_coarse_solve.has_value());
    CHECK_FALSE(config["Boundaries"]["WavePort"][0].contains("ComplexCoarseSolve"));
    CHECK(w2.eig_tol == w1.eig_tol);
    CHECK(w2.max_size == w1.max_size);
    CHECK(w2.verbose == w1.verbose);

    // Coverage gate. VoltagePath is an opt-in coordinate path for line integral
    // postprocessing on the port face; absence means no voltage line integral.
    // PolarityAttributes is an opt-in [high, low] terminal pair for fixing the mode
    // polarity; absence means the internal polarity convention is used. Both are
    // opt-in features with no meaningful default value to emit (mutually exclusive),
    // so they are deliberately not concretized. ComplexCoarseSolve is an opt-in per-port
    // override: an explicit value also takes precedence over Solver.Linear.PCMatReal, so
    // absence (inherit both global settings) has no equivalent explicit value to emit.
    // Target is opt-in: absence means the eigenvalue target from the port materials.
    auto wp_gaps = SchemaCoverageGaps(
        "/properties/Boundaries/properties/WavePort/items",
        config["Boundaries"]["WavePort"][0],
        /*skip=*/{"VoltagePath", "PolarityAttributes", "ComplexCoarseSolve", "Target"});
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
    auto mat_gaps = SchemaCoverageGaps("/properties/Domains/properties/Materials/items",
                                       config["Domains"]["Materials"][0],
                                       /*skip=*/{"MaterialAxes"});
    INFO("Domains.Materials[] missing keys: " << json(mat_gaps).dump());
    CHECK(mat_gaps.empty());
  }
}

#if defined(MFEM_USE_CUDSS)
TEST_CASE("Linear solver Type cuDSS parses when built with cuDSS", "[config][Serial]")
{
  // With cuDSS support built in, Solver.Linear.Type = "cuDSS" must parse to the
  // CUDSS enum and not trip the "not built with cuDSS" configuration guard.
  json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                 {"Model", {{"Mesh", "test.msh"}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                 {"Boundaries", json::object()},
                 {"Solver", {{"Linear", {{"Type", "cuDSS"}}}}}};
  IoData iodata(config, false);
  CHECK(iodata.solver.linear.type == LinearSolver::CUDSS);
}
#endif
