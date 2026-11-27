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

LibraryCheckReport Judge(const IoData &iodata, const LibraryCheckLibrary &library,
                         bool developer = false)
{
  const auto &request = *iodata.solver.electrostatic.response_correction;
  const auto targets = LibraryCheckTargets(iodata, request);
  return JudgeLibraryRows(FindLibraryLoadRows(iodata, request, targets, library), targets,
                          developer);
}

}  // namespace

TEST_CASE("LibraryMismatchFactor: a layer acts on its own Type, the substrate on all",
          "[surfaceresponselibrarycheck][Serial]")
{
  using Q = LibraryMismatchQuantity;
  using T = InterfaceDielectric;
  CHECK(LibraryMismatchFactor(Q::LAYER_PERMITTIVITY, T::MS, T::MS) == 0.94);
  CHECK(LibraryMismatchFactor(Q::LAYER_PERMITTIVITY, T::SA, T::SA) == 0.87);
  CHECK(LibraryMismatchFactor(Q::LAYER_PERMITTIVITY, T::MA, T::MA) == 0.90);
  CHECK(LibraryMismatchFactor(Q::LAYER_PERMITTIVITY, T::MS, T::SA) == 0.0);
  CHECK(LibraryMismatchFactor(Q::LAYER_THICKNESS, T::MS, T::MS) == 0.95);
  CHECK(LibraryMismatchFactor(Q::LAYER_THICKNESS, T::MA, T::MS) == 0.0);
  CHECK(LibraryMismatchFactor(Q::SUBSTRATE_PERMITTIVITY, T::DEFAULT, T::MS) == 1.86);
  CHECK(LibraryMismatchFactor(Q::SUBSTRATE_PERMITTIVITY, T::DEFAULT, T::MA) == 0.16);
  CHECK(LibraryMismatchFactor(Q::SUBSTRATE_PERMITTIVITY, T::DEFAULT, T::SA) == 0.11);
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
  for (const Case &c : {Case{"Permittivity", 1, 11.45 * 1.005, 0.94 * 0.005, false},
                        Case{"Permittivity", 1, 11.45 * 1.015, 0.94 * 0.015, true},
                        Case{"Permittivity", 0, 4.0 * 0.99, 0.87 * 0.01, false},
                        Case{"Permittivity", 0, 4.0 * 0.98, 0.87 * 0.02, true},
                        Case{"Thickness", 1, 0.002 * 1.005, 0.95 * 0.005, false},
                        Case{"Thickness", 1, 0.002 * 1.02, 0.95 * 0.02, true}})
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
    // The developer switch tolerates an effect beyond the rule.
    CHECK_FALSE(Judge(iodata, library, true).stop);
  }
  // Rows that stop at any difference, and the one that cannot verify.
  {
    json config = base;
    config["Boundaries"]["Postprocessing"]["Dielectric"][1]["EdgeDistances"] = {0.2001};
    IoData iodata(config, false);
    const auto report = Judge(iodata, library, true);
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
    auto report = Judge(iodata, other, true);
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
    CHECK(Judge(iodata, other, true).stop);
    other.length_unit.reset();
    CHECK(Judge(iodata, other).Empty());
    // Parked metal laws: every such model named without the switch, allowed with it.
    other.parked_models = {"cond-a", "cond-b"};
    iodata.solver.electrostatic.response_correction->developer = false;
    const auto report = Judge(iodata, other);
    CHECK(report.rows.size() == 2);
    CHECK(report.stop);
    iodata.solver.electrostatic.response_correction->developer = true;
    CHECK(Judge(iodata, other, true).rows.empty());
    other.refusal = "Unable to open fabrication-process response library!";
    const auto refused = Judge(iodata, other, true);
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
    std::string candidates;
    CHECK_FALSE(FindSubstrateLoadRow(near, library, Targets(near), &candidates));
    CHECK_THAT(candidates,
               ContainsSubstring("could be the substrate: Domains.Materials[1]"));
    config["Domains"]["Materials"][1]["Permittivity"] = 11.6;
    IoData far(config, false);
    const auto row = FindSubstrateLoadRow(far, library, Targets(far), &candidates);
    REQUIRE(row);
    CHECK_THAT(candidates, ContainsSubstring("could be the substrate: none"));
    const auto report = JudgeLibraryRows({*row}, Targets(far), false);
    CHECK(report.stop);
    CHECK_FALSE(JudgeLibraryRows({*row}, Targets(far), true).stop);
    // Vacuum alone: nothing to judge at load (the check after the mesh passes).
    config["Domains"]["Materials"] = {{{"Attributes", {1, 2}}}};
    IoData vacuum(config, false);
    CHECK_FALSE(FindSubstrateLoadRow(vacuum, library, Targets(vacuum), &candidates));
  }
  {
    // After the mesh: the material next to the metal; MS 1.86 x, SA 0.11 x.
    json config = base;
    config["Domains"]["Materials"][1]["Permittivity"] = 11.45 * 1.003;
    IoData iodata(config, false);
    const auto targets = Targets(iodata);
    const auto rows = SubstrateMeshRows(iodata, {{1, {2}}}, 11.45, targets);
    REQUIRE(rows.size() == 1);
    auto report = JudgeLibraryRows(rows, targets, false);
    CHECK_THAT(report.target_effect.at(5), WithinRel(1.86 * 0.003, 1.0e-9));
    CHECK_THAT(report.target_effect.at(4), WithinRel(0.11 * 0.003, 1.0e-9));
    CHECK_FALSE(report.stop);
    // THE EFFECTS ADD (the ruling's example): an MS layer 0.8 % off with a substrate 0.4 %
    // off is beyond the rule (0.94 x 0.8 + 1.86 x 0.4 = 1.5 %), each alone is not.
    config["Domains"]["Materials"][1]["Permittivity"] = 11.45 * 1.004;
    config["Boundaries"]["Postprocessing"]["Dielectric"][1]["Permittivity"] = 11.45 * 1.008;
    IoData both(config, false);
    const auto &request = *both.solver.electrostatic.response_correction;
    auto layer_rows = FindLibraryLoadRows(both, request, Targets(both), library);
    REQUIRE(layer_rows.size() == 1);
    CHECK_FALSE(JudgeLibraryRows(layer_rows, Targets(both), false).stop);
    auto substrate_rows = SubstrateMeshRows(both, {{1, {2}}}, 11.45, Targets(both));
    CHECK_FALSE(JudgeLibraryRows(substrate_rows, Targets(both), false).stop);
    layer_rows.insert(layer_rows.end(), substrate_rows.begin(), substrate_rows.end());
    report = JudgeLibraryRows(layer_rows, Targets(both), false);
    CHECK_THAT(report.target_effect.at(5), WithinRel(0.94 * 0.008 + 1.86 * 0.004, 1.0e-9));
    CHECK(report.stop);
    const std::string message = LibraryCheckMessage(report, library, "after the mesh", "");
    CHECK_THAT(message, ContainsSubstring("MS 1.5 % (Index 5)"));
    CHECK_THAT(message, ContainsSubstring("STOP (beyond the rule): MS layer permittivity"));
    CHECK_THAT(message,
               ContainsSubstring("STOP (beyond the rule): substrate permittivity"));
  }
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
      CHECK_THAT(printed, ContainsSubstring("SA 7.17 % (Index 4)"));
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
                 WithinRel(0.87 * 0.005, 1.0e-9));
      CHECK(block.at("Interfaces") == json::array({4}));
      CHECK(block.at("BeyondRule") == false);
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
               WithinRel(0.11 * (11.47 / 11.45 - 1.0), 1.0e-9));
  }
#endif
}

}  // namespace palace
