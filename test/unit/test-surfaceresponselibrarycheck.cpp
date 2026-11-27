// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"
#include "test-helpers.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "driver.hpp"
#include "models/surfaceparticipationtable.hpp"
#include "models/surfaceresponselibrarycheck.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"
#include "utils/timer.hpp"

namespace palace
{

using json = nlohmann::json;
namespace fs = std::filesystem;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinRel;

namespace
{

// The library the rows are checked against: an island config's MS-free layers plus an
// MS layer, the substrate 11.45, R 0.2, "um" (the config's L0 1e-6).
LibraryCheckLibrary IslandLibrary()
{
  LibraryCheckLibrary library;
  library.name = "unit-library";
  library.path = "unit-library.json";
  library.matching_radius = 0.2;
  library.interface_layers[InterfaceDielectric::SA] = {0.002, 4.0};
  library.interface_layers[InterfaceDielectric::MS] = {0.002, 11.45};
  library.substrate_permittivity = 11.45;
  library.length_unit = "um";
  library.profile = {{"MetalThickness", 0.1}, {"OveretchDepth", 0.05}};
  return library;
}

// The island's SA target (Index 4) and an MS target (Index 5) on the same metal, with a
// substrate domain of 11.45 beside the vacuum.
json IslandTwoTargetConfig(const test::SurfaceResponseFiles &files)
{
  json config = files.ConvexIslandConfig(test::CrackPolicy::Default);
  auto &dielectric = config["Boundaries"]["Postprocessing"]["Dielectric"];
  json ms = dielectric[0];
  ms["Index"] = 5;
  ms["Type"] = "MS";
  ms["Permittivity"] = 11.45;
  dielectric.push_back(ms);
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"] = {4, 5};
  config["Domains"]["Materials"] = {{{"Attributes", {1}}},
                                    {{"Attributes", {2}}, {"Permittivity", 11.45}}};
  return config;
}

LibraryCheckReport Judge(const IoData &iodata, const LibraryCheckLibrary &library)
{
  const auto &request = *iodata.solver.electrostatic.response_correction;
  const auto targets = LibraryCheckTargets(iodata, request);
  return JudgeLibraryRows(FindLibraryLoadRows(iodata, request, targets, library), targets);
}

// The check at load as CheckResponseLibraryAtLoad runs it: the load rows, the substrate row
// and its line, the verdict and the message.
struct LoadCheck
{
  LibraryCheckReport report;
  LibrarySubstrateLine line;
  std::string message;
};

LoadCheck CheckAtLoad(const json &config, const LibraryCheckLibrary &library)
{
  IoData iodata(config, false);
  const auto &request = *iodata.solver.electrostatic.response_correction;
  const auto targets = LibraryCheckTargets(iodata, request);
  auto rows = FindLibraryLoadRows(iodata, request, targets, library);
  LoadCheck check;
  if (auto row = FindSubstrateLoadRow(iodata, library, targets, rows, &check.line))
  {
    rows.push_back(std::move(*row));
  }
  check.report = JudgeLibraryRows(std::move(rows), targets);
  check.message = LibraryCheckMessage(check.report, library, "at load", check.line);
  return check;
}

}  // namespace

TEST_CASE("LibraryMismatchFactor: a layer acts on its own Type, the substrate on all",
          "[surfaceresponselibrarycheck][Serial]")
{
  using Q = LibraryMismatchQuantity;
  using T = InterfaceDielectric;
  CHECK(LibraryMismatchFactor(Q::LAYER_PERMITTIVITY, T::MS, T::MS) == 1.0);
  CHECK(LibraryMismatchFactor(Q::LAYER_PERMITTIVITY, T::SA, T::SA) == 1.0);
  CHECK(LibraryMismatchFactor(Q::LAYER_PERMITTIVITY, T::MA, T::MA) == 1.0);
  CHECK(LibraryMismatchFactor(Q::LAYER_PERMITTIVITY, T::MS, T::SA) == 0.0);
  CHECK(LibraryMismatchFactor(Q::LAYER_THICKNESS, T::MS, T::MS) == 1.0);
  CHECK(LibraryMismatchFactor(Q::LAYER_THICKNESS, T::MA, T::MS) == 0.0);
  // The larger of two windows' readings (CTX-002, CTX-003), rounded up to one decimal.
  CHECK(LibraryMismatchFactor(Q::SUBSTRATE_PERMITTIVITY, T::DEFAULT, T::MS) == 1.9);
  CHECK(LibraryMismatchFactor(Q::SUBSTRATE_PERMITTIVITY, T::DEFAULT, T::MA) == 0.3);
  CHECK(LibraryMismatchFactor(Q::SUBSTRATE_PERMITTIVITY, T::DEFAULT, T::SA) == 0.2);
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "The check at load: one row per quantity inside and beyond its tolerance",
                 "[surfaceresponselibrarycheck][Serial]")
{
  const json base = IslandTwoTargetConfig(*this);
  const auto library = IslandLibrary();
  {
    IoData iodata(base, false);
    const auto report = Judge(iodata, library);
    CHECK(report.Empty());
    CHECK_FALSE(report.stop);
  }
  // The layers (permittivity, thickness): judged by the rule on the summed effect.
  struct Case
  {
    std::string key;
    int entry;
    double value, effect;
    bool stop;
  };
  for (const Case &c : {Case{"Permittivity", 1, 11.45 * 1.005, 0.005, false},
                        Case{"Permittivity", 1, 11.45 * 1.015, 0.015, true},
                        Case{"Permittivity", 0, 4.0 * 0.995, 0.005, false},
                        Case{"Permittivity", 0, 4.0 * 0.98, 0.02, true},
                        Case{"Thickness", 1, 0.002 * 1.005, 0.005, false},
                        Case{"Thickness", 1, 0.002 * 1.02, 0.02, true}})
  {
    CAPTURE(c.key, c.entry, c.value);
    json config = base;
    config["Boundaries"]["Postprocessing"]["Dielectric"][c.entry][c.key] = c.value;
    IoData iodata(config, false);
    const auto report = Judge(iodata, library);
    REQUIRE(report.rows.size() == 1);
    CHECK(report.rows[0].kind == LibraryCheckRow::Kind::EFFECT);
    const int index = c.entry == 0 ? 4 : 5;
    CHECK_THAT(report.target_effect.at(index), WithinRel(c.effect, 1.0e-9));
    CHECK(report.stop == c.stop);
    CHECK(report.beyond_rule == c.stop);
  }
  // Rows that stop at any difference, and the one that cannot verify.
  {
    json config = base;
    config["Boundaries"]["Postprocessing"]["Dielectric"][1]["EdgeDistances"] = {0.2001};
    IoData iodata(config, false);
    const auto report = Judge(iodata, library);
    REQUIRE(report.rows.size() == 1);
    CHECK(report.rows[0].quantity == "R (the largest EdgeDistances)");
    CHECK(report.stop);
    config["Boundaries"]["Postprocessing"]["Dielectric"][1]["EdgeDistances"] = {
        0.2 * (1.0 + 1.0e-12)};
    CHECK(Judge(IoData(config, false), library).Empty());
  }
  {
    IoData iodata(base, false);
    auto other = library;
    other.interface_layers.erase(InterfaceDielectric::MS);
    auto report = Judge(iodata, other);
    REQUIRE(report.rows.size() == 1);
    CHECK(report.rows[0].quantity == "MS InterfaceLayers entry");
    CHECK(report.stop);
    other.interface_layers.clear();
    report = Judge(iodata, other);
    REQUIRE(report.rows.size() == 1);
    CHECK(report.rows[0].kind == LibraryCheckRow::Kind::WARNING);
    CHECK_FALSE(report.stop);
  }
  {
    IoData iodata(base, false);
    auto other = library;
    other.length_unit = "nm";
    CHECK(Judge(iodata, other).stop);
    // The micro forms name 1e-6 m; an unknown word is named as such.
    for (const std::string unit : {"\u00b5m", "\u03bcm", "micron"})
    {
      CAPTURE(unit);
      other.length_unit = unit;
      CHECK(Judge(iodata, other).Empty());
    }
    other.length_unit = "microns-ish";
    const auto unknown = Judge(iodata, other);
    REQUIRE(unknown.rows.size() == 1);
    CHECK(unknown.stop);
    CHECK_THAT(unknown.rows[0].library_value, ContainsSubstring("unknown LengthUnit word"));
    other.length_unit.reset();
    CHECK(Judge(iodata, other).Empty());
    // Parked metal laws: every such model named without the switch, allowed with it.
    other.parked_models = {"cond-a", "cond-b"};
    iodata.solver.electrostatic.response_correction->developer = false;
    const auto report = Judge(iodata, other);
    CHECK(report.rows.size() == 2);
    CHECK(report.stop);
    // MINOR-4 (second review): the remedy then names the switch for the parked coupons.
    CHECK_THAT(LibraryCheckMessage(report, other, "at load", {}),
               ContainsSubstring("made for the configuration's values, or set "
                                 "\"Developer\": true for the parked coupons named above"));
    iodata.solver.electrostatic.response_correction->developer = true;
    CHECK(Judge(iodata, other).rows.empty());
    other.refusal = "Unable to open fabrication-process response library!";
    const auto refused = Judge(iodata, other);
    REQUIRE(refused.rows.size() == 1);
    CHECK(refused.rows[0].quantity == "library file");
    CHECK(refused.stop);
  }
  {
    json config = base;
    config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"] = {4, 5, 6};
    const auto report = Judge(IoData(config, false), library);
    REQUIRE(report.rows.size() == 1);
    CHECK(report.rows[0].entry == "TargetInterfaces entry 6");
    CHECK(report.stop);
  }
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "The substrate: at load a material must be able to pass; after the mesh "
                 "its effect adds to the layers'",
                 "[surfaceresponselibrarycheck][Serial]")
{
  const json base = IslandTwoTargetConfig(*this);
  const auto library = IslandLibrary();
  const auto Targets = [](const IoData &iodata)
  { return LibraryCheckTargets(iodata, *iodata.solver.electrostatic.response_correction); };
  {
    // At load: 11.47 (+0.175 %) could be the substrate; 11.6 (+1.3 %) cannot, and no other
    // material but vacuum is configured: the check after the mesh cannot pass.
    json config = base;
    config["Domains"]["Materials"][1]["Permittivity"] = 11.47;
    IoData near(config, false);
    LibrarySubstrateLine candidates;
    CHECK_FALSE(FindSubstrateLoadRow(near, library, Targets(near), {}, &candidates));
    CHECK_THAT(candidates.text,
               ContainsSubstring("could be the substrate: Domains.Materials[1]"));
    config["Domains"]["Materials"][1]["Permittivity"] = 11.6;
    IoData far(config, false);
    const auto row = FindSubstrateLoadRow(far, library, Targets(far), {}, &candidates);
    REQUIRE(row);
    CHECK_THAT(candidates.text, ContainsSubstring("could be the substrate: none"));
    CHECK(JudgeLibraryRows({*row}, Targets(far)).stop);
    // Vacuum alone with an MS target: no material can be the substrate (row (b) of the
    // table at load): a STOP row.
    config["Domains"]["Materials"] = {{{"Attributes", {1, 2}}}};
    IoData vacuum(config, false);
    const auto none =
        FindSubstrateLoadRow(vacuum, library, Targets(vacuum), {}, &candidates);
    REQUIRE(none);
    CHECK(none->kind == LibraryCheckRow::Kind::STOP);
  }
  {
    // After the mesh: the material next to the metal; MS 1.9 x, SA 0.2 x.
    json config = base;
    config["Domains"]["Materials"][1]["Permittivity"] = 11.45 * 1.003;
    IoData iodata(config, false);
    const auto targets = Targets(iodata);
    const auto rows = SubstrateMeshRows(iodata, {{1, {2}}}, 11.45, targets);
    REQUIRE(rows.size() == 1);
    auto report = JudgeLibraryRows(rows, targets);
    CHECK_THAT(report.target_effect.at(5), WithinRel(1.9 * 0.003, 1.0e-9));
    CHECK_THAT(report.target_effect.at(4), WithinRel(0.2 * 0.003, 1.0e-9));
    CHECK_FALSE(report.stop);
    // THE EFFECTS ADD (the ruling's example): an MS layer 0.8 % off with a substrate 0.4 %
    // off is beyond the rule (0.8 + 1.9 x 0.4 = 1.56 %), each alone is not.
    config["Domains"]["Materials"][1]["Permittivity"] = 11.45 * 1.004;
    config["Boundaries"]["Postprocessing"]["Dielectric"][1]["Permittivity"] = 11.45 * 1.008;
    IoData both(config, false);
    const auto &request = *both.solver.electrostatic.response_correction;
    auto layer_rows = FindLibraryLoadRows(both, request, Targets(both), library);
    REQUIRE(layer_rows.size() == 1);
    CHECK_FALSE(JudgeLibraryRows(layer_rows, Targets(both)).stop);
    auto substrate_rows = SubstrateMeshRows(both, {{1, {2}}}, 11.45, Targets(both));
    CHECK_FALSE(JudgeLibraryRows(substrate_rows, Targets(both)).stop);
    layer_rows.insert(layer_rows.end(), substrate_rows.begin(), substrate_rows.end());
    report = JudgeLibraryRows(layer_rows, Targets(both));
    CHECK_THAT(report.target_effect.at(5), WithinRel(0.008 + 1.9 * 0.004, 1.0e-9));
    CHECK(report.stop);
    const std::string message = LibraryCheckMessage(report, library, "after the mesh", {});
    CHECK_THAT(message, ContainsSubstring("MS 1.56 % (Index 5)"));
    CHECK_THAT(message, ContainsSubstring("STOP (beyond the rule): MS layer permittivity"));
    CHECK_THAT(message,
               ContainsSubstring("STOP (beyond the rule): substrate permittivity"));
    CHECK_THAT(message, ContainsSubstring("no override"));
    CHECK(message.find("parked coupons") == std::string::npos);
  }
  {
    // Per target the LARGEST substrate row counts (one patch lies on one substrate): two
    // materials next to the metal, +0.3 % and +0.4 %, read 1.9 x 0.4 % on MS, not the sum.
    json config = base;
    config["Domains"]["Materials"][1]["Permittivity"] = 11.45 * 1.003;
    config["Domains"]["Materials"].push_back(
        {{"Attributes", {3}}, {"Permittivity", 11.45 * 1.004}});
    IoData two(config, false);
    const auto rows = SubstrateMeshRows(two, {{1, {2}}, {2, {3}}}, 11.45, Targets(two));
    REQUIRE(rows.size() == 2);
    const auto report = JudgeLibraryRows(rows, Targets(two));
    CHECK_THAT(report.target_effect.at(5), WithinRel(1.9 * 0.004, 1.0e-9));
    CHECK_FALSE(report.stop);
  }
}

// MAJOR-1 of the code review (lane v1loadcheck): at load each candidate substrate is judged
// WITH the layer rows, so that a load message that tolerates is never followed by a stop of
// the rule after the mesh (one material but vacuum).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "The substrate at load is judged with the layer rows: a summed overrun "
                 "stops at load",
                 "[surfaceresponselibrarycheck][Serial]")
{
  const json base = IslandTwoTargetConfig(*this);
  const auto library = IslandLibrary();
  auto AtLoad = [&](double ms_layer, double substrate, LibrarySubstrateLine *candidates)
  {
    json config = base;
    config["Boundaries"]["Postprocessing"]["Dielectric"][1]["Permittivity"] = ms_layer;
    config["Domains"]["Materials"][1]["Permittivity"] = substrate;
    IoData iodata(config, false);
    const auto &request = *iodata.solver.electrostatic.response_correction;
    const auto targets = LibraryCheckTargets(iodata, request);
    auto rows = FindLibraryLoadRows(iodata, request, targets, library);
    if (auto row = FindSubstrateLoadRow(iodata, library, targets, rows, candidates))
    {
      rows.push_back(std::move(*row));
    }
    return JudgeLibraryRows(std::move(rows), targets);
  };
  LibrarySubstrateLine candidates;
  // MS 0.9 % off with the substrate 0.1 % off: each alone is inside the rule, together
  // 0.9 + 1.9 x 0.1 = 1.09 % on MS: STOP AT LOAD, the substrate row listed.
  auto report = AtLoad(11.45 * 1.009, 11.45 * 1.001, &candidates);
  CHECK(report.stop);
  REQUIRE(report.rows.size() == 2);
  CHECK(report.rows[1].substrate);
  CHECK_THAT(report.target_effect.at(5), WithinRel(0.009 + 1.9 * 0.001, 1.0e-9));
  CHECK_THAT(candidates.text, ContainsSubstring("could be the substrate: none"));
  CHECK_FALSE(AtLoad(11.45 * 1.009, 11.45, &candidates).stop);
  CHECK_FALSE(AtLoad(11.45, 11.45 * 1.001, &candidates).stop);
  // The reviewer's flip-chip combination: the designers' substrate 11.47 (+0.175 %) with MS
  // 11.55 (+0.87 %): 0.87 + 0.33 = 1.2 %: STOP AT LOAD. The designers' 11.47 on both:
  // 0.175 + 0.33 = 0.51 %: tolerated, the material named as the substrate.
  CHECK(AtLoad(11.55, 11.47, &candidates).stop);
  report = AtLoad(11.47, 11.47, &candidates);
  CHECK_FALSE(report.stop);
  CHECK(report.rows.size() == 1);
  CHECK_THAT(candidates.text,
             ContainsSubstring("could be the substrate: Domains.Materials[1]"));
  // What the load message tolerates passes after the mesh too (the same layer rows with
  // the material's own row).
  {
    json config = base;
    config["Boundaries"]["Postprocessing"]["Dielectric"][1]["Permittivity"] = 11.47;
    config["Domains"]["Materials"][1]["Permittivity"] = 11.47;
    IoData iodata(config, false);
    const auto &request = *iodata.solver.electrostatic.response_correction;
    const auto targets = LibraryCheckTargets(iodata, request);
    auto rows = FindLibraryLoadRows(iodata, request, targets, library);
    const auto substrate = SubstrateMeshRows(iodata, {{1, {2}}}, 11.45, targets);
    rows.insert(rows.end(), substrate.begin(), substrate.end());
    CHECK_FALSE(JudgeLibraryRows(rows, targets).stop);
  }
  // A substrate equal to the library's is no row: the layer rows alone stop the run.
  report = AtLoad(11.45 * 1.015, 11.45, &candidates);
  CHECK(report.stop);
  CHECK(report.rows.size() == 1);
  CHECK_FALSE(report.rows[0].substrate);
}

// MAJOR-1 and MINOR-2 (g) of the second review (lane v1loadcheck): every material that
// cannot pass the summed rule is named at load and the verdict is conditional (no stop at
// load: it may never touch the corrected metal); a thickness row adds at load too.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "At load a material that would stop the run after the mesh is named and "
                 "the verdict is conditional",
                 "[surfaceresponselibrarycheck][Serial]")
{
  json base = IslandTwoTargetConfig(*this);
  base["Domains"]["Materials"].push_back({{"Attributes", {3}}, {"Permittivity", 11.45}});
  const auto library = IslandLibrary();
  auto AtLoad = [&](const json &config, LibrarySubstrateLine *line)
  {
    IoData iodata(config, false);
    const auto &request = *iodata.solver.electrostatic.response_correction;
    const auto targets = LibraryCheckTargets(iodata, request);
    auto rows = FindLibraryLoadRows(iodata, request, targets, library);
    if (auto row = FindSubstrateLoadRow(iodata, library, targets, rows, line))
    {
      rows.push_back(std::move(*row));
    }
    return JudgeLibraryRows(std::move(rows), targets);
  };
  LibrarySubstrateLine line;
  // The reviewer's h32: MS +0.5 % with a second material at +0.3 % (0.5 + 1.9 x 0.3 =
  // 1.07 % on MS): Materials[1] can be the substrate, Materials[2] is named.
  json config = base;
  config["Boundaries"]["Postprocessing"]["Dielectric"][1]["Permittivity"] = 11.45 * 1.005;
  config["Domains"]["Materials"][2]["Permittivity"] = 11.45 * 1.003;
  auto report = AtLoad(config, &line);
  CHECK_FALSE(report.stop);
  CHECK(line.conditional == std::vector<std::string>{"Domains.Materials[2]"});
  CHECK_THAT(line.text, ContainsSubstring("could be the substrate: Domains.Materials[1]"));
  CHECK_THAT(line.text, ContainsSubstring("Domains.Materials[2] (attributes 3) 11.48435 "
                                          "(+0.3 %): beyond the rule with the layer rows: "
                                          "it stops the run after the mesh if it touches "
                                          "the corrected metal"));
  std::string message = LibraryCheckMessage(report, library, "at load", line);
  CHECK_THAT(message, ContainsSubstring("-> TOLERATED unless Domains.Materials[2] touches "
                                        "the corrected metal (then the run stops after "
                                        "the mesh): the run goes on"));
  // The reviewer's h14: no layer row, the second material at +1 % (1.9 % on MS).
  config = base;
  config["Domains"]["Materials"][2]["Permittivity"] = 11.45 * 1.01;
  report = AtLoad(config, &line);
  CHECK(report.Empty());
  message = LibraryCheckMessage(report, library, "at load", line);
  CHECK(message.find("no mismatch") == std::string::npos);
  CHECK_THAT(message, ContainsSubstring("no row at load\n  -> TOLERATED unless "
                                        "Domains.Materials[2] touches"));
  // MINOR-2 (g): an MS THICKNESS +0.9 % with the substrate +0.1 %: STOP at load.
  config = base;
  config["Boundaries"]["Postprocessing"]["Dielectric"][1]["Thickness"] = 0.002 * 1.009;
  config["Domains"]["Materials"][1]["Permittivity"] = 11.45 * 1.001;
  config["Domains"]["Materials"].erase(2);
  report = AtLoad(config, &line);
  CHECK(report.stop);
  REQUIRE(report.rows.size() == 2);
  CHECK(report.rows[1].substrate);
}

// MINOR-1 of the second review: an effect of exactly the limit, written in decimal, is
// tolerated at load and after the mesh alike (one test of the limit).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "A slip of exactly the limit is tolerated at load and after the mesh",
                 "[surfaceresponselibrarycheck][Serial]")
{
  json config = IslandTwoTargetConfig(*this);
  auto &dielectric = config["Boundaries"]["Postprocessing"]["Dielectric"];
  json ma = dielectric[0];
  ma["Index"] = 6;
  ma["Type"] = "MA";
  ma["Permittivity"] = 10.1;  // the library's 10.0: exactly 1 %, 1.0000000000000009 %
  dielectric.push_back(ma);
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"] = {4, 5, 6};
  auto library = IslandLibrary();
  library.interface_layers[InterfaceDielectric::MA] = {0.002, 10.0};
  auto report = Judge(IoData(config, false), library);
  REQUIRE(report.rows.size() == 1);
  CHECK_FALSE(report.stop);
  CHECK_THAT(LibraryCheckMessage(report, library, "at load", {}),
             ContainsSubstring("TOLERATED: MA layer permittivity"));
  config["Boundaries"]["Postprocessing"]["Dielectric"][2]["Permittivity"] = 10.1001;
  CHECK(Judge(IoData(config, false), library).stop);  // 1.001 %: beyond
  // The substrate exactly at its limit alone (1 / 1.9 % on MS): a candidate at load, and
  // tolerated after the mesh.
  config = IslandTwoTargetConfig(*this);
  const double at_limit = 11.45 * (1.0 + 0.01 / 1.9);
  config["Domains"]["Materials"][1]["Permittivity"] = at_limit;
  IoData iodata(config, false);
  const auto targets =
      LibraryCheckTargets(iodata, *iodata.solver.electrostatic.response_correction);
  LibrarySubstrateLine line;
  CHECK_FALSE(FindSubstrateLoadRow(iodata, IslandLibrary(), targets, {}, &line));
  CHECK_THAT(line.text, ContainsSubstring("could be the substrate: Domains.Materials[1]"));
  report = JudgeLibraryRows(SubstrateMeshRows(iodata, {{1, {2}}}, 11.45, targets), targets);
  CHECK_FALSE(report.stop);
}

// MINOR-3 of the second review (the supervisor's ruling): the vacuum side is Palace's own
// test; beside a true vacuum a material nearer to 1 than to the substrate (air 1.0006, 3.9)
// is neither: named at load, a STOP row after the mesh only when it touches the metal.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "The vacuum side: beside a vacuum, air 1.0006 and 3.9 are neither the "
                 "substrate nor vacuum",
                 "[surfaceresponselibrarycheck][Serial]")
{
  const auto library = IslandLibrary();
  for (const double permittivity : {1.0006, 3.9})
  {
    CAPTURE(permittivity);
    json config = IslandTwoTargetConfig(*this);
    config["Domains"]["Materials"].push_back(
        {{"Attributes", {3}}, {"Permittivity", permittivity}});
    IoData iodata(config, false);
    const auto targets =
        LibraryCheckTargets(iodata, *iodata.solver.electrostatic.response_correction);
    LibrarySubstrateLine line;
    CHECK_FALSE(FindSubstrateLoadRow(iodata, library, targets, {}, &line));
    const auto touching = SubstrateMeshRows(iodata, {{2, {3}}}, 11.45, targets);
    const auto report = JudgeLibraryRows(touching, targets);
    const std::string message = LibraryCheckMessage(report, library, "after the mesh", {});
    CHECK(message.find("-91") == std::string::npos);
    CHECK(line.conditional == std::vector<std::string>{"Domains.Materials[2]"});
    CHECK_THAT(line.text,
               ContainsSubstring(fmt::format("Domains.Materials[2] (attributes 3) "
                                             "{}: not the library's substrate "
                                             "and not vacuum (Palace takes",
                                             permittivity)));
    CHECK_THAT(line.text,
               ContainsSubstring("use a library made for it): it stops the run after the "
                                 "mesh if it touches the corrected metal"));
    REQUIRE(touching.size() == 1);
    CHECK(touching[0].kind == LibraryCheckRow::Kind::STOP);
    CHECK(report.stop);
    CHECK_THAT(message, ContainsSubstring("STOP: material next to the corrected metal"));
    CHECK(JudgeLibraryRows(SubstrateMeshRows(iodata, {}, 11.45, targets), targets).Empty());
  }
  // Air 1.0006 as the only air with an SA target: no vacuum side anywhere: STOP AT LOAD;
  // with MS targets only, no row; exact 1.0: no row.
  json config = IslandTwoTargetConfig(*this);
  config["Domains"]["Materials"][0]["Permittivity"] = 1.0006;
  auto report = Judge(IoData(config, false), library);
  REQUIRE(report.rows.size() == 1);
  CHECK(report.stop);
  CHECK(report.rows[0].quantity == "vacuum side (MA and SA targets)");
  CHECK(report.rows[0].config_value == "Domains.Materials[0] (attributes 1) 1.0006");
  CHECK_THAT(report.rows[0].tolerance,
             ContainsSubstring("set the air's permittivity to exactly 1.0"));
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"] = {5};
  CHECK(Judge(IoData(config, false), library).Empty());
  config = IslandTwoTargetConfig(*this);
  config["Domains"]["Materials"][0]["Permittivity"] = 1.0;
  CHECK(Judge(IoData(config, false), library).Empty());
}

// THE DECISION TABLE AT LOAD (reference.md), one case per row. Each configured material is
// vacuum (Palace's test), of the substrate class or neither. The island config: the vacuum
// Materials[0], the substrate 11.45 Materials[1], an SA (Index 4) and an MS (Index 5)
// target.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "The decision table at load row (a): an MA or SA target and no vacuum material "
    "stop the run at load",
    "[surfaceresponselibrarycheck][Serial]")
{
  const auto library = IslandLibrary();
  json config = IslandTwoTargetConfig(*this);
  config["Domains"]["Materials"][0]["Permittivity"] = 1.0006;
  auto check = CheckAtLoad(config, library);
  REQUIRE(check.report.rows.size() == 1);
  CHECK(check.report.stop);
  CHECK(check.report.rows[0].quantity == "vacuum side (MA and SA targets)");
  // MINOR-2 of the third review: no library can serve air that is not vacuum.
  CHECK_THAT(check.message,
             ContainsSubstring("  -> STOP: set the air's permittivity to exactly 1.0\n"));
  CHECK(check.message.find("library's value") == std::string::npos);
  // Another row beside it: its remedy follows.
  config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Permittivity"] = 4.0 * 1.05;
  check = CheckAtLoad(config, library);
  CHECK(check.report.rows.size() == 2);
  CHECK_THAT(check.message, ContainsSubstring("-> STOP: set the air's permittivity to "
                                              "exactly 1.0; set each other entry above "
                                              "to the library's value"));
  // MS targets only: no row (a).
  config = IslandTwoTargetConfig(*this);
  config["Domains"]["Materials"][0]["Permittivity"] = 1.0006;
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"] = {5};
  check = CheckAtLoad(config, library);
  CHECK(check.report.Empty());
  CHECK_THAT(check.message, ContainsSubstring("-> TOLERATED unless Domains.Materials[0]"));
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "The decision table at load row (b): an MS or SA target and no material of the "
    "substrate class stop the run at load",
    "[surfaceresponselibrarycheck][Serial]")
{
  const auto library = IslandLibrary();
  // MAJOR-1 of the third review: the substrate 3.9 (an 11.45 library) beside a vacuum and
  // an MS target: not a substrate slip, but no material can be the library's substrate.
  json config = IslandTwoTargetConfig(*this);
  config["Domains"]["Materials"][1]["Permittivity"] = 3.9;
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"] = {5};
  auto check = CheckAtLoad(config, library);
  REQUIRE(check.report.rows.size() == 1);
  CHECK(check.report.stop);
  CHECK(check.report.rows[0].kind == LibraryCheckRow::Kind::STOP);
  CHECK(check.report.rows[0].entry ==
        "Domains.Materials: no configured material can be the library's substrate");
  CHECK(check.report.rows[0].config_value == "Domains.Materials[1] (attributes 2) 3.9");
  CHECK_THAT(check.message, ContainsSubstring("STOP: substrate (any corrected target)"));
  CHECK(check.message.find("-65") == std::string::npos);
  CHECK(check.message.find("TOLERATED") == std::string::npos);
  // v1loadcheckfix3-review MINOR-1: the SA target alone stops at load the same way.
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"] = {4};
  check = CheckAtLoad(config, library);
  REQUIRE(check.report.rows.size() == 1);
  CHECK(check.report.stop);
  CHECK(check.report.rows[0].kind == LibraryCheckRow::Kind::STOP);
  CHECK(check.report.rows[0].entry ==
        "Domains.Materials: no configured material can be the library's substrate");
  // The SA target too; 6.0 (row (b)) and 7.0 (nearer to 11.45: row (c)) stop at load alike.
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"] = {4, 5};
  for (const double permittivity : {3.9, 6.0, 7.0})
  {
    CAPTURE(permittivity);
    config["Domains"]["Materials"][1]["Permittivity"] = permittivity;
    check = CheckAtLoad(config, library);
    REQUIRE(check.report.rows.size() == 1);
    CHECK(check.report.stop);
    CHECK(
        check.report.rows[0].kind ==
        (permittivity < 6.5 ? LibraryCheckRow::Kind::STOP : LibraryCheckRow::Kind::EFFECT));
  }
  // With no vacuum either (3.9 the only material): rows (a) and (b) both.
  config["Domains"]["Materials"] = {{{"Attributes", {1, 2}}, {"Permittivity", 3.9}}};
  CHECK(CheckAtLoad(config, library).report.rows.size() == 2);
  // The reviewer's sub1: the substrate set to 1.0 is vacuum; with an MS target no material
  // is of the substrate class: STOP at load (it ran to the end before).
  config = IslandTwoTargetConfig(*this);
  config["Domains"]["Materials"][1]["Permittivity"] = 1.0;
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"] = {5};
  check = CheckAtLoad(config, library);
  REQUIRE(check.report.rows.size() == 1);
  CHECK(check.report.stop);
  CHECK(check.report.rows[0].config_value == "none (vacuum only)");
}

// The supervisor's ruling on v1loadcheckfix3-review MINOR-2 (decision 2414): row (b) holds
// for every target Type; a metal edge with no substrate is not what the library was built
// for.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "The decision table at load row (b) for MA targets: no material of the substrate class "
    "stops the run at load",
    "[surfaceresponselibrarycheck][Serial]")
{
  json config = IslandTwoTargetConfig(*this);
  json ma = config["Boundaries"]["Postprocessing"]["Dielectric"][0];
  ma["Index"] = 6;
  ma["Type"] = "MA";
  ma["Permittivity"] = 10.0;
  config["Boundaries"]["Postprocessing"]["Dielectric"].push_back(ma);
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"] = {6};
  auto library = IslandLibrary();
  library.interface_layers[InterfaceDielectric::MA] = {0.002, 10.0};
  // Vacuum and a material of the substrate class: no row (b).
  CHECK(CheckAtLoad(config, library).report.Empty());
  // Every material vacuum (the substrate entered as 1.0), then vacuum beside 3.9 (neither):
  // no material can be the library's substrate: STOP at load (both ran before).
  for (const double permittivity : {1.0, 3.9})
  {
    CAPTURE(permittivity);
    config["Domains"]["Materials"][1]["Permittivity"] = permittivity;
    const auto check = CheckAtLoad(config, library);
    REQUIRE(check.report.rows.size() == 1);
    CHECK(check.report.stop);
    CHECK(check.report.rows[0].entry ==
          "Domains.Materials: no configured material can be the library's substrate");
    CHECK(check.report.rows[0].config_value ==
          (permittivity == 1.0 ? "none (vacuum only)"
                               : "Domains.Materials[1] (attributes 2) 3.9"));
  }
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "The decision table at load row (c): materials of the substrate class of which none "
    "passes the summed rule stop the run at load",
    "[surfaceresponselibrarycheck][Serial]")
{
  const auto library = IslandLibrary();
  // The MS layer +0.5 % (tolerated alone); the substrate +0.4 % (MS 0.76 %) and a second
  // one +1 %: no material passes with the layer row; the closest one's row stops the run.
  json config = IslandTwoTargetConfig(*this);
  config["Boundaries"]["Postprocessing"]["Dielectric"][1]["Permittivity"] = 11.45 * 1.005;
  config["Domains"]["Materials"][1]["Permittivity"] = 11.45 * 1.004;
  config["Domains"]["Materials"].push_back(
      {{"Attributes", {3}}, {"Permittivity", 11.45 * 1.01}});
  const auto check = CheckAtLoad(config, library);
  REQUIRE(check.report.rows.size() == 2);
  CHECK(check.report.stop);
  CHECK(check.report.beyond_rule);
  CHECK(check.report.rows[1].entry == "Domains.Materials[1] (attributes 2; the closest)");
  CHECK_THAT(check.line.text, ContainsSubstring("could be the substrate: none"));
  CHECK_THAT(check.message,
             ContainsSubstring("STOP (beyond the rule): substrate permittivity"));
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "The decision table at load row (d): a failing material beside a passing one is "
    "named and the verdict is conditional",
    "[surfaceresponselibrarycheck][Serial]")
{
  const auto library = IslandLibrary();
  json config = IslandTwoTargetConfig(*this);
  config["Domains"]["Materials"].push_back(
      {{"Attributes", {3}}, {"Permittivity", 11.45 * 1.01}});
  const auto check = CheckAtLoad(config, library);
  CHECK(check.report.Empty());
  CHECK(check.line.conditional == std::vector<std::string>{"Domains.Materials[2]"});
  CHECK_THAT(check.line.text,
             ContainsSubstring("Domains.Materials[2] (attributes 3) 11.5645 "
                               "(+1 %): beyond the rule with the layer "
                               "rows"));
  CHECK_THAT(check.message, ContainsSubstring("no row at load\n  -> TOLERATED unless "
                                              "Domains.Materials[2] touches"));
  // After the mesh, next to the metal: beyond the rule (MS 1.9 %).
  IoData iodata(config, false);
  const auto targets =
      LibraryCheckTargets(iodata, *iodata.solver.electrostatic.response_correction);
  const auto report =
      JudgeLibraryRows(SubstrateMeshRows(iodata, {{2, {3}}}, 11.45, targets), targets);
  CHECK(report.stop);
  CHECK(report.beyond_rule);
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "The decision table at load row (e): a material that is neither beside a vacuum and a "
    "passing substrate is named and the verdict is conditional",
    "[surfaceresponselibrarycheck][Serial]")
{
  const auto library = IslandLibrary();
  json config = IslandTwoTargetConfig(*this);
  config["Domains"]["Materials"].push_back({{"Attributes", {3}}, {"Permittivity", 3.9}});
  const auto check = CheckAtLoad(config, library);
  CHECK(check.report.Empty());
  CHECK(check.line.conditional == std::vector<std::string>{"Domains.Materials[2]"});
  // The ruled words (the third review's MAJOR-1): the air or a dielectric.
  CHECK_THAT(
      check.line.text,
      ContainsSubstring("3.9: not the library's substrate and not vacuum (Palace takes "
                        "the vacuum side of MA and SA interfaces only where the "
                        "material is vacuum: permittivity and permeability exactly 1, "
                        "no frequency dependence; if it is the air, set its "
                        "permittivity to exactly 1.0; if it is a dielectric, use a "
                        "library made for it)"));
  CHECK(check.line.text.find("set it to 1.0") == std::string::npos);
  CHECK_THAT(check.message, ContainsSubstring("-> TOLERATED unless Domains.Materials[2]"));
  // After the mesh: a plain STOP row only when it is next to the metal.
  IoData iodata(config, false);
  const auto targets =
      LibraryCheckTargets(iodata, *iodata.solver.electrostatic.response_correction);
  const auto touching = SubstrateMeshRows(iodata, {{2, {3}}}, 11.45, targets);
  REQUIRE(touching.size() == 1);
  CHECK(touching[0].kind == LibraryCheckRow::Kind::STOP);
  CHECK(touching[0].quantity == "material next to the corrected metal");
  // v1loadcheckfix3-review MINOR-4: its remedy names the air and a dielectric; the
  // library's value (11.45) is not offered for it.
  const auto printed = LibraryCheckMessage(JudgeLibraryRows(touching, targets), library,
                                           "after the mesh", {});
  CHECK_THAT(printed,
             ContainsSubstring("  -> STOP: if it is the air, set its permittivity to "
                               "exactly 1.0; if it is a dielectric, use a library "
                               "made for it\n"));
  CHECK(printed.find("set each") == std::string::npos);
  CHECK(JudgeLibraryRows(SubstrateMeshRows(iodata, {}, 11.45, targets), targets).Empty());
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "The check at load lists every mismatch in one message and stops before "
                 "the mesh in both forms",
                 "[surfaceresponselibrarycheck][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  json config = ConvexIslandConfig(test::CrackPolicy::Default);
  // No mesh file: a run the check let through would stop at the mesh load instead.
  config["Model"]["Mesh"] = (temp.temp_dir / "no-such-mesh.mesh").string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  auto &sa = config["Boundaries"]["Postprocessing"]["Dielectric"][0];
  sa["EdgeDistances"] = {0.25};     // R: the library's is 0.2
  sa["Permittivity"] = 4.0 * 1.05;  // +5 %: 4.35 % of an SA value
  sa["Thickness"] = 0.002 * 1.03;   // +3 %: 2.82 %
  for (const bool preflight : {false, true})
  {
    CAPTURE(preflight);
    json run = config;
    run["Problem"]["Output"] =
        (temp.temp_dir / (preflight ? "three-pre" : "three")).string();
    run["Solver"]["Electrostatic"]["ResponseCorrection"]["PreflightOnly"] = preflight;
    IoData iodata(run, false);
    iodata.solver.electrostatic.response_correction->developer = false;
    MakeOutputFolder(iodata, comm);
    std::string what;
    const std::string printed = test::CaptureStdout(
        [&]
        {
          try
          {
            if (preflight)
            {
              RunSurfaceResponsePreflight(iodata, comm, 1);
            }
            else
            {
              palace::Run(iodata, comm, 1);
            }
          }
          catch (const std::exception &e)
          {
            what = e.what();
          }
        });
    CHECK_THAT(what, ContainsSubstring("library check at load stops the run"));
    if (Mpi::Root(comm))
    {
      CHECK_THAT(printed, ContainsSubstring("3 row(s):"));
      CHECK_THAT(printed, ContainsSubstring("STOP: R (the largest EdgeDistances); "
                                            "Boundaries.Postprocessing.Dielectric Index 4 "
                                            "(SA; attributes 9)"));
      CHECK_THAT(printed, ContainsSubstring("STOP (beyond the rule): SA layer thickness"));
      CHECK_THAT(printed,
                 ContainsSubstring("STOP (beyond the rule): SA layer permittivity"));
      CHECK_THAT(printed, ContainsSubstring("SA 8 % (Index 4)"));
    }
  }
  // NO OVERRIDE (the supervisor's ruling of the fix round): an effect beyond the rule stops
  // the run at load with "Developer": true as without it.
  for (const bool developer : {false, true})
  {
    CAPTURE(developer);
    json run = config;
    run["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"] = {0.2};
    run["Boundaries"]["Postprocessing"]["Dielectric"][0]["Thickness"] = 0.002;
    IoData iodata(run, false);
    iodata.solver.electrostatic.response_correction->developer = developer;
    std::string what;
    const std::string printed = test::CaptureStdout(
        [&]
        {
          try
          {
            CheckResponseLibraryAtLoad(iodata, comm, false);
          }
          catch (const std::exception &e)
          {
            what = e.what();
          }
        });
    CHECK_THAT(what, ContainsSubstring("library check at load stops the run"));
    if (Mpi::Root(comm))
    {
      CHECK_THAT(printed, ContainsSubstring("1 row(s):"));
      CHECK_THAT(printed,
                 ContainsSubstring("STOP (beyond the rule): SA layer permittivity"));
      CHECK_THAT(printed, ContainsSubstring("no override"));
      CHECK(printed.find("TOLERATED") == std::string::npos);
    }
  }
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "A tolerated mismatch is named in the outputs; without one nothing is added",
    "[surfaceresponselibrarycheck][surfaceparticipation][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "library-check-island.mesh";
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(comm);
  json config = ConvexIslandConfig(test::CrackPolicy::Default);
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["CorrectionMode"] =
      "PostprocessOnly";
  const auto ReadText = [](const fs::path &path)
  {
    std::ifstream input(path);
    return std::string(std::istreambuf_iterator<char>(input), {});
  };
  for (const double permittivity : {4.0, 4.0 * 1.005})
  {
    for (const bool preflight : {false, true})
    {
      CAPTURE(permittivity, preflight);
      json run = config;
      const fs::path output =
          temp.temp_dir / fmt::format("mismatch-{}-{}", permittivity != 4.0, preflight);
      run["Problem"]["Output"] = output.string();
      run["Boundaries"]["Postprocessing"]["Dielectric"][0]["Permittivity"] = permittivity;
      run["Solver"]["Electrostatic"]["ResponseCorrection"]["PreflightOnly"] = preflight;
      IoData iodata(run, false);
      MakeOutputFolder(iodata, comm);
      BlockTimer::Reset();
      if (preflight)
      {
        CHECK(RunSurfaceResponsePreflight(iodata, comm, 1) == 0);
      }
      else
      {
        palace::Run(iodata, comm, 1);
      }
      Mpi::Barrier(comm);
      if (!Mpi::Root(comm))
      {
        continue;
      }
      const std::string table = ReadText(output / kSurfaceParticipationFile);
      const json meta = json::parse(ReadText(output / "palace.json"));
      REQUIRE(!table.empty());
      if (permittivity == 4.0)
      {
        CHECK(table.find(kParticipationLibraryMismatchFlag) == std::string::npos);
        CHECK_FALSE(meta.at("SurfaceResponse").contains("LibraryMismatch"));
        continue;
      }
      // Every SA row (the target and its Type total) names it; palace.json lists it.
      std::size_t rows = 0, flagged = 0;
      std::istringstream lines(table);
      std::string line;
      std::getline(lines, line);  // the header
      while (std::getline(lines, line))
      {
        rows += line.find(" SA,") != std::string::npos;
        flagged += line.find(" SA,") != std::string::npos &&
                   line.find(kParticipationLibraryMismatchFlag) != std::string::npos;
      }
      CHECK(rows > 0);
      CHECK(flagged == rows);
      const auto &block = meta.at("SurfaceResponse").at("LibraryMismatch");
      REQUIRE(block.at("Mismatches").size() == 1);
      CHECK(block.at("Mismatches")[0].at("Quantity") == "SA layer permittivity");
      CHECK_THAT(block.at("Mismatches")[0].at("Mismatch").get<double>(),
                 WithinRel(0.005, 1.0e-9));
      CHECK_THAT(block.at("EstimatedEffectPerType").at("SA").get<double>(),
                 WithinRel(0.005, 1.0e-9));
      CHECK(block.at("Interfaces") == json::array({4}));
      CHECK_FALSE(block.contains("BeyondRule"));  // a recorded mismatch is inside the rule
    }
  }
#endif
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "A substrate mismatch tolerated after the mesh is in the outputs of both forms",
    "[surfaceresponselibrarycheck][surfaceparticipation][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "library-check-substrate-island.mesh";
  const fs::path library_path = temp.temp_dir / "library-check-substrate-library.json";
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
    std::ifstream input(convex_library_3d_path);
    json library = json::parse(input);
    library["Fabrication"]["SubstratePermittivity"] = 11.45;
    std::ofstream library_output(library_path);
    library_output << library.dump(2) << "\n";
  }
  Mpi::Barrier(comm);
  json config = ConvexIslandConfig(test::CrackPolicy::Default);
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  config["Domains"]["Materials"][0]["Permittivity"] = 11.47;  // the island's one domain
  // A vacuum material (on an attribute the mesh lacks): an SA target needs one at load.
  config["Domains"]["Materials"].push_back({{"Attributes", {7}}});
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction["Library"] = library_path.string();
  correction["CorrectionMode"] = "PostprocessOnly";
  for (const bool preflight : {false, true})
  {
    CAPTURE(preflight);
    json run = config;
    const fs::path output = temp.temp_dir / fmt::format("substrate-{}", preflight);
    run["Problem"]["Output"] = output.string();
    run["Solver"]["Electrostatic"]["ResponseCorrection"]["PreflightOnly"] = preflight;
    IoData iodata(run, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    if (preflight)
    {
      CHECK(RunSurfaceResponsePreflight(iodata, comm, 1) == 0);
    }
    else
    {
      palace::Run(iodata, comm, 1);
    }
    Mpi::Barrier(comm);
    if (!Mpi::Root(comm))
    {
      continue;
    }
    std::ifstream input(output / "palace.json");
    const json meta = json::parse(input);
    const auto &block = meta.at("SurfaceResponse").at("LibraryMismatch");
    REQUIRE(block.at("Mismatches").size() == 1);
    CHECK(block.at("Mismatches")[0].at("Quantity") == "substrate permittivity");
    CHECK_THAT(block.at("EstimatedEffectPerType").at("SA").get<double>(),
               WithinRel(0.2 * (11.47 / 11.45 - 1.0), 1.0e-9));
  }
#endif
}

// MAJOR-1, MINOR-2 (h) and MINOR-3 of the second review through the run (both forms): the
// load message names the material that the check after the mesh then stops on, with the
// layer rows (h); a material neither the substrate nor vacuum stops the run only when it
// touches the corrected metal. The island's metal touches its one domain (attribute 1); a
// material on attribute 7 (not in the mesh) never touches it.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "A material named at load stops the run after the mesh only if it touches the metal",
    "[surfaceresponselibrarycheck][surfaceparticipation][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "library-check-named-island.mesh";
  const fs::path library_path = temp.temp_dir / "library-check-named-library.json";
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
    std::ifstream input(convex_library_3d_path);
    json library = json::parse(input);
    library["Fabrication"]["SubstratePermittivity"] = 11.45;
    std::ofstream library_output(library_path);
    library_output << library.dump(2) << "\n";
  }
  Mpi::Barrier(comm);
  json base = ConvexIslandConfig(test::CrackPolicy::Default);
  base["Model"]["Mesh"] = mesh_path.string();
  base["Model"]["L0"] = 1.0;
  base["Model"]["Lc"] = 1.0;
  auto &correction = base["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction["Library"] = library_path.string();
  correction["CorrectionMode"] = "PostprocessOnly";
  // The candidate substrate, equal to the library's, and a vacuum material (an SA target
  // needs one at load), on attributes the mesh lacks.
  base["Domains"]["Materials"].push_back({{"Attributes", {7}}, {"Permittivity", 11.45}});
  base["Domains"]["Materials"].push_back({{"Attributes", {8}}});
  int count = 0;
  const auto RunIt = [&](const json &config, bool preflight, std::string *what)
  {
    json run = config;
    run["Problem"]["Output"] = (temp.temp_dir / fmt::format("named-{}", count++)).string();
    run["Solver"]["Electrostatic"]["ResponseCorrection"]["PreflightOnly"] = preflight;
    IoData iodata(run, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    what->clear();
    return test::CaptureStdout(
        [&]
        {
          try
          {
            if (preflight)
            {
              RunSurfaceResponsePreflight(iodata, comm, 1);
            }
            else
            {
              palace::Run(iodata, comm, 1);
            }
          }
          catch (const std::exception &e)
          {
            *what = e.what();
          }
        });
  };
  for (const bool preflight : {false, true})
  {
    CAPTURE(preflight);
    std::string what;
    // (h): the domain at the substrate +1 % (SA 0.2 % alone) with the SA layer +0.9 %:
    // named at load, the verdict conditional; after the mesh 1.1 % on SA: STOP.
    json config = base;
    config["Domains"]["Materials"][0]["Permittivity"] = 11.45 * 1.01;
    config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Permittivity"] = 4.0 * 1.009;
    std::string printed = RunIt(config, preflight, &what);
    CHECK_THAT(what, ContainsSubstring("library check after the mesh stops the run"));
    if (Mpi::Root(comm))
    {
      CHECK_THAT(printed, ContainsSubstring("Domains.Materials[0] (attributes 1) 11.5645 "
                                            "(+1 %): beyond the rule with the layer rows"));
      CHECK_THAT(printed, ContainsSubstring("-> TOLERATED unless Domains.Materials[0] "
                                            "touches the corrected metal"));
      CHECK_THAT(printed,
                 ContainsSubstring("STOP (beyond the rule): substrate permittivity"));
    }
    // A material of 3.9 next to the metal: named at load, a plain STOP after the mesh.
    config = base;
    config["Domains"]["Materials"][0]["Permittivity"] = 3.9;
    printed = RunIt(config, preflight, &what);
    CHECK_THAT(what, ContainsSubstring("library check after the mesh stops the run"));
    if (Mpi::Root(comm))
    {
      CHECK_THAT(printed,
                 ContainsSubstring("3.9: not the library's substrate and not vacuum"));
      CHECK_THAT(printed, ContainsSubstring("STOP: material next to the corrected metal"));
      CHECK(printed.find("-91") == std::string::npos);
    }
  }
  // The same material where it does not touch the metal: named at load, no stop.
  json config = base;
  config["Domains"]["Materials"][0]["Permittivity"] = 11.45;
  config["Domains"]["Materials"][1]["Permittivity"] = 3.9;
  std::string what;
  const std::string printed = RunIt(config, true, &what);
  CHECK(what.empty());
  if (Mpi::Root(comm))
  {
    CHECK_THAT(printed,
               ContainsSubstring("-> TOLERATED unless Domains.Materials[1] touches"));
    CHECK(printed.find("library check after the mesh") == std::string::npos);
  }
#endif
}

}  // namespace palace
