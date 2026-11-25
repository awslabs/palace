// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "driver.hpp"
#include "fem/mesh.hpp"
#include "models/surfaceparticipationtable.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"
#include "utils/timer.hpp"

using namespace palace;
using namespace Catch::Matchers;
using json = nlohmann::json;
namespace fs = std::filesystem;

namespace
{

std::string ReadBytes(const fs::path &path)
{
  std::ifstream input(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

// The cells of a csv file, each without its padding, the header row first.
std::vector<std::vector<std::string>> CsvCells(const fs::path &path)
{
  std::vector<std::vector<std::string>> rows;
  std::ifstream input(path);
  std::string line;
  while (std::getline(input, line))
  {
    std::vector<std::string> cells;
    std::stringstream stream(line);
    std::string cell;
    while (std::getline(stream, cell, ','))
    {
      const auto first = cell.find_first_not_of(' ');
      cells.push_back(first == std::string::npos
                          ? std::string()
                          : cell.substr(first, cell.find_last_not_of(' ') - first + 1));
    }
    rows.push_back(std::move(cells));
  }
  return rows;
}

// The members of palace.json's SurfaceResponse.Resolution that a corrected run fixes before
// its first solve (the others are filled by its solves).
constexpr const char *kPreSolveResolutionMembers[] = {
    "Enabled", "Measure", "Order", "Radius", "Bars", "ShareRule", "Knots", "Preflight"};

}  // namespace

// The configuration key PreflightOnly (tracer fixes T1): the response-correction block of
// the configured problem type selects the preflight exactly as --surface-response-preflight
// does; a key in the block of another problem type, or absent, selects nothing.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "PreflightOnly selects the surface-response preflight from the block of "
                 "the configured problem type",
                 "[surfaceresponseoperator][preflightonly][Serial]")
{
  // P4 S3-0 (review r2084 F-1): both crack forms (the legacy key; upstream's default = the
  // key unset), as P3-23.
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  json electrostatic = ConvexIslandConfig(crack_policy);
  json maxwell = ConvexMaxwellIslandConfig(crack_policy);
  CHECK(SurfaceResponsePreflightKey(IoData(electrostatic, false)).empty());
  CHECK(SurfaceResponsePreflightKey(IoData(maxwell, false)).empty());

  json electrostatic_key = electrostatic;
  electrostatic_key["Solver"]["Electrostatic"]["ResponseCorrection"]["PreflightOnly"] =
      true;
  CHECK(SurfaceResponsePreflightKey(IoData(electrostatic_key, false)) ==
        "Solver.Electrostatic.ResponseCorrection.PreflightOnly");
  json maxwell_key = maxwell;
  maxwell_key["Solver"]["SurfaceResponseCorrection"]["PreflightOnly"] = true;
  CHECK(SurfaceResponsePreflightKey(IoData(maxwell_key, false)) ==
        "Solver.SurfaceResponseCorrection.PreflightOnly");
  json maxwell_false = maxwell;
  maxwell_false["Solver"]["SurfaceResponseCorrection"]["PreflightOnly"] = false;
  CHECK(SurfaceResponsePreflightKey(IoData(maxwell_false, false)).empty());

  // A block of the other problem type is refused with its keys by the existing
  // configuration check, so a key there is never silently ignored.
  json maxwell_other = maxwell;
  maxwell_other["Solver"]["Electrostatic"] = electrostatic_key["Solver"]["Electrostatic"];
  CHECK_THROWS(IoData(maxwell_other, false));

  json invalid = electrostatic;
  invalid["Solver"]["Electrostatic"]["ResponseCorrection"]["PreflightOnly"] = "Yes";
  CHECK_THROWS_WITH(IoData(invalid, false),
                    ContainsSubstring("\"PreflightOnly\" must be a boolean"));
}

// The key changes nothing inside the preflight: the same configuration with and without it
// writes byte-identical manifests and patch dry runs.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "PreflightOnly leaves the preflight's manifest and patch dry run "
                 "byte-identical",
                 "[surfaceresponseoperator][preflightonly][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // P4 S3-0 (review r2084 F-1): both crack forms (the legacy key; upstream's default = the
  // key unset), as P3-23.
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  MPI_Comm comm = Mpi::World();
  auto Preflight = [&](bool key, const fs::path &dir)
  {
    json config = ConvexIslandConfig(crack_policy);
    config["Problem"]["Output"] = dir.string();
    if (key)
    {
      config["Solver"]["Electrostatic"]["ResponseCorrection"]["PreflightOnly"] = true;
    }
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    MakeOutputFolder(iodata, comm);
    Mesh mesh(MakeIslandMesh());
    const auto outcome = WriteSurfaceResponseRequirements(
        iodata, mesh, (dir / "surface-response-requirements.json").string());
    Mpi::Barrier(comm);
    CHECK(outcome.fail_closed_summary.empty());
  };
  const fs::path without = temp.temp_dir / "preflight-without-key";
  const fs::path with = temp.temp_dir / "preflight-with-key";
  Preflight(false, without);
  Preflight(true, with);
  // The patch dry run byte for byte; the manifest equal but for its one run-metadata leaf
  // (the containment test's wall time).
  REQUIRE(fs::is_regular_file(without / "surface-response-patches.csv"));
  CHECK(ReadBytes(without / "surface-response-patches.csv") ==
        ReadBytes(with / "surface-response-patches.csv"));
  auto Manifest = [](const fs::path &dir)
  {
    json manifest = json::parse(ReadBytes(dir / "surface-response-requirements.json"));
    auto &diagnostics = manifest["Identification"]["Diagnostics"];
    if (diagnostics.contains("DomainBoundaryExclusions"))
    {
      diagnostics["DomainBoundaryExclusions"].erase("WallTime");
    }
    return manifest;
  };
  CHECK(Manifest(without) == Manifest(with));
#endif
}

// CUT-PLAN step 4b (lane v1step4b): PreflightOnly evaluates on the first mesh, without a
// solve, the resolution census a corrected run computes before its first solve: the
// Resolution block's pre-solve members and surface-response-resolution-patches.csv equal
// the solve form's on the same mesh (also after edge-tube refinement), no energy cell is
// measured, and surface-participation.csv holds the level-1 PREFLIGHT rows (no value).
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "PreflightOnly writes the resolution census a corrected run computes "
    "before its first solve",
    "[surfaceresponseoperator][preflightonly][surfaceparticipation][3d][Serial]"
    "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  std::string problem;
  bool edge_tubes = false;
  std::tie(problem, edge_tubes) =
      GENERATE(table<std::string, bool>({{"Electrostatic", false},
                                         {"Electrostatic", true},
                                         {"Eigenmode", false},
                                         {"Driven", false}}));
  INFO("problem " << problem << ", edge tubes " << edge_tubes);
  MPI_Comm comm = Mpi::World();
  const std::string stem = "preflight-census-" + problem + (edge_tubes ? "-tubes" : "");
  const fs::path mesh_path = temp.temp_dir / (stem + ".mesh");
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(comm);
  const bool electrostatic = problem == "Electrostatic";
  json config = electrostatic ? ConvexIslandConfig(test::CrackPolicy::Default)
                              : ConvexMaxwellIslandConfig(test::CrackPolicy::Default);
  config["Problem"]["Type"] = problem;
  config["Model"]["Mesh"] = mesh_path.string();
  // Mesh units are metres, unscaled, as the library's R (m).
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  const json::json_pointer correction_key(electrostatic
                                              ? "/Solver/Electrostatic/ResponseCorrection"
                                              : "/Solver/SurfaceResponseCorrection");
  config[correction_key]["CorrectionMode"] = "PostprocessOnly";
  config[correction_key]["PatchEnergy"] = true;
  if (edge_tubes)
  {
    config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeRefinement"] = {
        {"Radius", 0.2}, {"ElementsPerRadius", 2}};
  }
  if (problem == "Eigenmode")
  {
    config["Solver"]["Eigenmode"] = {{"Target", 0.05}, {"N", 2}};
  }
  if (problem == "Driven")
  {
    config["Solver"].erase("Eigenmode");
    config["Boundaries"]["Ground"]["Attributes"] = {2, 3, 4, 5, 6, 9};
    config["Boundaries"]["LumpedPort"] = {{{"Index", 1},
                                           {"Attributes", {1}},
                                           {"R", 50.0},
                                           {"Direction", "+X"},
                                           {"Excitation", 1}}};
    config["Solver"]["Driven"] = {{"Samples", {{{"Type", "Point"}, {"Freq", {0.04}}}}}};
  }
  auto Run = [&](bool preflight)
  {
    json run = config;
    const fs::path output = temp.temp_dir / (stem + (preflight ? "-preflight" : "-solve"));
    run["Problem"]["Output"] = output.string();
    run[correction_key]["PreflightOnly"] = preflight;
    IoData iodata(run, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    if (preflight)
    {
      CHECK(RunSurfaceResponsePreflight(iodata, comm, utils::ConfigureOmp()) == 0);
    }
    else
    {
      palace::Run(iodata, comm, utils::ConfigureOmp());
    }
    Mpi::Barrier(comm);
    return output;
  };
  const fs::path solve = Run(false);
  const fs::path preflight = Run(true);
  if (!Mpi::Root(comm))
  {
    return;
  }
  auto Block = [](const fs::path &dir, const char *key)
  { return json::parse(ReadBytes(dir / "palace.json")).at("SurfaceResponse").at(key); };
  const json solved = Block(solve, "Resolution"), census = Block(preflight, "Resolution");
  for (const char *member : kPreSolveResolutionMembers)
  {
    INFO("SurfaceResponse.Resolution." << member);
    CHECK(census.at(member) == solved.at(member));
  }
  REQUIRE(census.at("Knots").at("Located").get<long long int>() > 0);
  REQUIRE(!census.at("Preflight").empty());
  CHECK(solved.at("Cells").get<int>() > 0);
  CHECK(census.at("Cells") == 0);
  CHECK(census.at("NoCellsReason") == "PreflightOnly: no solve");
  CHECK(census.at("Complete") == true);
  CHECK(census.at("Final") == true);
  CHECK(census.at("Level") == 1);
  CHECK(census.at("Levels") == 1);
  const std::string patches = "surface-response-resolution-patches.csv";
  REQUIRE(fs::is_regular_file(solve / patches));
  CHECK(ReadBytes(preflight / patches) == ReadBytes(solve / patches));
  CHECK(fs::is_regular_file(solve / "surface-response-resolution.csv"));
  CHECK_FALSE(fs::exists(preflight / "surface-response-resolution.csv"));
  // The table: the Type total, then interface 4 (the island's SA target), of iteration 1
  // with the first solve's unknowns; the PREFLIGHT word of the Type; no value.
  const auto rows = CsvCells(preflight / kSurfaceParticipationFile);
  const auto solved_rows = CsvCells(solve / kSurfaceParticipationFile);
  REQUIRE(rows.size() == 3);
  REQUIRE(solved_rows.size() > 1);
  auto Index = [&](const std::string &name)
  {
    for (std::size_t i = 0; i < rows[0].size(); i++)
    {
      if (rows[0][i] == name || rows[0][i].rfind(name + " (", 0) == 0)
      {
        return i;
      }
    }
    FAIL("no column " << name);
    return std::size_t(0);
  };
  REQUIRE(solved_rows[0] == rows[0]);  // the same 27 columns
  std::string word;
  for (const auto &row : census.at("Preflight"))
  {
    if (row.at("Type") == "SA" && row.at("Plane").is_null())
    {
      word = "PREFLIGHT " + row.at("Status").get<std::string>();
    }
  }
  REQUIRE(!word.empty());
  for (std::size_t r = 1; r < rows.size(); r++)
  {
    INFO("row " << r);
    const auto &row = rows[r];
    REQUIRE(row.size() == rows[0].size());
    CHECK(row[Index("iteration")] == "1");
    CHECK(row[Index("unknowns")] == solved_rows[1][Index("unknowns")]);
    CHECK(row[Index("problem")] == problem);
    CHECK(row[Index("point")] == "0");
    CHECK(row[Index("exc")] == "0");
    CHECK(row[Index("interface")] == (r == 1 ? "0" : "4"));
    CHECK(row[Index("scope")] == (r == 1 ? "Type" : "interface"));
    CHECK(row[Index("resolution")] == word);
    for (const char *name :
         {"p raw", "p corrected", "corrected / raw", "Q raw", "Q corrected", "h/R m90",
          "h/R m50", "uncovered+pending share of W", "change corrected", "change raw",
          "change per doubling"})
    {
      INFO(name);
      CHECK(row[Index(name)] == "nan");
    }
    for (const char *name : {"coverage", "regime", "confidence", "convergence"})
    {
      INFO(name);
      CHECK(row[Index(name)] == kParticipationNotComputed);
    }
    CHECK(row[Index("flags")] == kParticipationPreflightFlags);
    CHECK(row[Index("reason")] == kParticipationPreflightReason);
  }
  CHECK(Block(preflight, "ParticipationTable").at("PreflightOnly") == true);
  CHECK_FALSE(Block(solve, "ParticipationTable").contains("PreflightOnly"));
#endif
}
