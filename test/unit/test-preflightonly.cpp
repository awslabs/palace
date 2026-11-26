// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"
#include "test-helpers.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <Eigen/Dense>
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
#include "utils/tablecsv.hpp"
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
                                         {"Eigenmode", true},
                                         {"Driven", false},
                                         {"Driven", true}}));
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
  REQUIRE(solved_rows[0] == rows[0]);  // the same 28 columns
  // ONE verdict (lane v1verdict1): the Type's pre-solve word is bare, the table's word of
  // the census weighted by length (LengthM90, h / R) with no convergence: never RESOLVED;
  // the retired share-rule words are gone.
  std::string word;
  for (const auto &row : census.at("Preflight"))
  {
    CHECK_FALSE(row.contains("EdgeStatus"));
    CHECK_FALSE(row.contains("SpatialStatus"));
    const auto &m90 = row.at("LengthM90");
    CHECK(row.at("Status") ==
          ParticipationPreflightWord(
              m90.is_null() ? std::numeric_limits<double>::quiet_NaN() : m90.get<double>(),
              census.at("Order").get<int>(), BuiltInParticipationVerdictBars(), problem));
    CHECK(row.at("Status") != "RESOLVED");
    if (row.at("Type") == "SA" && row.at("Plane").is_null())
    {
      word = row.at("Status").get<std::string>();
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
          "change per doubling", "S corner+cluster share"})
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
  // THE FLAG's guarantee on a real mesh (lane v1verdict1 fix round, BLOCKER-1): a Type's
  // all-planes PREFLIGHT row reads every cell its table census weighs, by the same h/R:
  // every solved row's h/R m90 is at most the row's LargestM. Order 1: no bar, no share.
  std::map<std::string, double> largest;
  for (const auto &row : census.at("Preflight"))
  {
    CHECK(row.at("UnderResolvedLengthShare").is_null());
    if (row.at("Plane").is_null())
    {
      REQUIRE(row.at("LargestM").is_number());
      largest[row.at("Type").get<std::string>()] = row.at("LargestM").get<double>();
    }
  }
  int m90_read = 0;
  for (std::size_t r = 1; r < solved_rows.size(); r++)
  {
    const double m90 = std::stod(solved_rows[r][Index("h/R m90")]);
    if (std::isfinite(m90))
    {
      INFO("solved row " << r);
      REQUIRE(largest.count(solved_rows[r][Index("Type")]) == 1);
      CHECK(m90 <= largest.at(solved_rows[r][Index("Type")]) * (1.0 + 1.0e-5));  // print
      m90_read++;
    }
  }
  CHECK(m90_read > 0);
  CHECK(Block(preflight, "ParticipationTable").at("PreflightOnly") == true);
  CHECK_FALSE(Block(solve, "ParticipationTable").contains("PreflightOnly"));
#endif
}

// CUT-PLAN step 4b (lane v1step4b; the supervisor's ruling (b) of 2026-11-25): a design
// whose preflight fails closed by default (GeometryFailurePolicy absent) still gets the
// census a solve-form run (RecordRaw, its default) computes before its first solve, and
// then ends with exit status 3; an explicit FailClosed builds none; an explicit RecordRaw
// completes with it. The geometry fallback: InterfaceMappingConflict (R2A), the two
// diagonal corners of the diagonal island whose arms carry the two SA targets 4 and 5
// (test-geometryfallback.cpp).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "PreflightOnly builds the census of a design whose preflight fails closed "
                 "by default",
                 "[surfaceresponseoperator][preflightonly][geometryfallback][3d][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "preflight-census-diagonal-island.mesh";
  const fs::path sa_library = temp.temp_dir / "preflight-census-convex-3d-sa.json";
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialDiagonalIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
    json library = json::parse(ReadBytes(convex_library_3d_path));
    library["Name"] = "unit-test-process-convex-3d-sa";
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
    }
    std::ofstream(sa_library) << library.dump(2) << "\n";
  }
  Mpi::Barrier(comm);
  json config = ConvexIslandConfig(test::CrackPolicy::Default);
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  config["Boundaries"]["Terminal"] = {{{"Index", 1}, {"Attributes", {9, 10}}}};
  json &dielectric = config["Boundaries"]["Postprocessing"]["Dielectric"];
  json second = dielectric[0];
  second["Index"] = 5;
  second["Attributes"] = {10};
  dielectric.push_back(second);
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction["TargetInterfaces"] = {4, 5};
  correction["Library"] = sa_library.string();
  correction["UnmatchedPolicy"] = "Warn";
  correction.erase("PatchConstruction");
  correction["CorrectionMode"] = "PostprocessOnly";
  correction["PatchEnergy"] = true;
  // policy: "" = the key absent.
  auto Run = [&](bool preflight, const std::string &policy)
  {
    json run = config;
    const fs::path output =
        temp.temp_dir /
        ("preflight-census-fallback-" + std::string(preflight ? "preflight" : "solve") +
         (policy.empty() ? "" : "-" + policy));
    run["Problem"]["Output"] = output.string();
    auto &request = run["Solver"]["Electrostatic"]["ResponseCorrection"];
    request["PreflightOnly"] = preflight;
    if (!policy.empty())
    {
      request["GeometryFailurePolicy"] = policy;
    }
    IoData iodata(run, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    int status = 0;
    std::string log;
    if (preflight)
    {
      log = test::CaptureStdout(
          [&]
          { status = RunSurfaceResponsePreflight(iodata, comm, utils::ConfigureOmp()); });
    }
    else
    {
      palace::Run(iodata, comm, utils::ConfigureOmp());
    }
    Mpi::Barrier(comm);
    return std::tuple{output, status, log};
  };
  const auto [solve, solve_status, solve_log] = Run(false, "");
  const auto [absent, absent_status, absent_log] = Run(true, "");
  const auto [closed, closed_status, closed_log] = Run(true, "FailClosed");
  const auto [raw, raw_status, raw_log] = Run(true, "RecordRaw");
  CHECK(absent_status == kSurfaceResponsePreflightFailClosedExitStatus);
  CHECK(closed_status == kSurfaceResponsePreflightFailClosedExitStatus);
  CHECK(raw_status == 0);
  if (!Mpi::Root(comm))
  {
    return;
  }
  // The order of the lines of the default fail-closed exit (R2): the preflight's summary,
  // then the census, then the line that names how the census was built.
  const auto summary =
      absent_log.find("Surface-response preflight ended: GeometryFailurePolicy FailClosed");
  const auto knots = absent_log.find("Resolution census: ");
  const auto built = absent_log.find("The resolution census above");
  REQUIRE(summary != std::string::npos);
  REQUIRE(knots != std::string::npos);
  REQUIRE(built != std::string::npos);
  CHECK(summary < knots);
  CHECK(knots < built);
  auto Response = [](const fs::path &dir)
  {
    return json::parse(ReadBytes(dir / "palace.json"))
        .value("SurfaceResponse", json::object());
  };
  const json solved = Response(solve).at("Resolution");
  REQUIRE(solved.at("Knots").at("Located").get<long long int>() > 0);
  const std::string patches = "surface-response-resolution-patches.csv";
  REQUIRE(fs::is_regular_file(solve / patches));
  for (const auto &dir : {absent, raw})
  {
    INFO("run " << dir.filename().string());
    const json census = Response(dir).at("Resolution");
    for (const char *member : kPreSolveResolutionMembers)
    {
      INFO("SurfaceResponse.Resolution." << member);
      CHECK(census.at(member) == solved.at(member));
    }
    CHECK(census.at("NoCellsReason") == "PreflightOnly: no solve");
    CHECK(ReadBytes(dir / patches) == ReadBytes(solve / patches));
    const auto rows = CsvCells(dir / kSurfaceParticipationFile);
    REQUIRE(rows.size() == 4);  // the SA Type total, interfaces 4 and 5
    CHECK(rows[1].back() == kParticipationPreflightReason);
  }
  // The explicit FailClosed: the manifest only, no census output.
  CHECK(fs::is_regular_file(closed / "surface-response-requirements.json"));
  CHECK_FALSE(Response(closed).contains("Resolution"));
  CHECK_FALSE(Response(closed).contains("ParticipationTable"));
  CHECK_FALSE(fs::exists(closed / patches));
  CHECK_FALSE(fs::exists(closed / kSurfaceParticipationFile));
#endif
}

// CUT-PLAN step 4b (lane v1step4b; the supervisor's ruling R1): where a solve-form run
// stops before its census, PreflightOnly builds none, prints one line and ends with the
// exit status the preflight had before the census existed (0 here). The convex island,
// PatchConstruction Features (the default): (a) a library without its corner model (the
// corners Missing) under "UnmatchedPolicy": "Error" (a solve form aborts; "Warn" and the
// key absent build the census); (b) a library with no model (the solve form's reader
// refuses it); (c) an explicit "GeometryFailurePolicy": "FailClosed" on a library that
// matches nothing (the solve form's NoMatchedFeature, recorded at the operator's
// construction only, ends it there; RecordRaw gives the raw-fallback line). (f)-(k): a
// library file the operator's model loop refuses (STEP 4c, lane v1verdict1: the Try forms
// of its readers) is predicted like (d) and (e), every refused model named.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "PreflightOnly builds no resolution census where a solve-form run stops "
                 "before it",
                 "[surfaceresponseoperator][preflightonly][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const std::string problem = GENERATE("Electrostatic", "Eigenmode");
  INFO("problem " << problem);
  MPI_Comm comm = Mpi::World();
  const bool electrostatic = problem == "Electrostatic";
  const std::string stem = "preflight-refusal-" + problem;
  const fs::path mesh_path = temp.temp_dir / (stem + ".mesh");
  const fs::path no_corner = temp.temp_dir / (stem + "-no-corner.json");
  const fs::path no_model = temp.temp_dir / (stem + "-no-model.json");
  const fs::path no_match = temp.temp_dir / (stem + "-no-match.json");
  const fs::path negative_matrix = temp.temp_dir / (stem + "-negative-fabricated.csv");
  const fs::path negative_mode = temp.temp_dir / (stem + "-negative-mode.json");
  const fs::path no_file = temp.temp_dir / (stem + "-no-file.json");
  const fs::path absent_matrix = temp.temp_dir / (stem + "-absent-fabricated.csv");
  const fs::path surface_absent = temp.temp_dir / (stem + "-surface-absent.json");
  const fs::path absent_surface = temp.temp_dir / (stem + "-absent-surface.csv");
  const fs::path points_absent = temp.temp_dir / (stem + "-points-absent.json");
  const fs::path absent_points = temp.temp_dir / (stem + "-absent-points.csv");
  const fs::path no_surface = temp.temp_dir / (stem + "-no-surface.json");
  const fs::path short_matrix = temp.temp_dir / (stem + "-short-fabricated.csv");
  const fs::path short_library = temp.temp_dir / (stem + "-short.json");
  const fs::path no_column_matrix = temp.temp_dir / (stem + "-no-column-fabricated.csv");
  const fs::path no_column_library = temp.temp_dir / (stem + "-no-column.json");
  const fs::path no_entry_matrix = temp.temp_dir / (stem + "-no-entry-fabricated.csv");
  const fs::path no_entry_library = temp.temp_dir / (stem + "-no-entry.json");
  // The fix round (MAJOR-1, MINOR-3, MINOR-4): a non-numeric domain cell, a binary domain
  // file, a header-only surface file, an absent thin surface file, two refused files of one
  // model.
  const fs::path text_cell_matrix = temp.temp_dir / (stem + "-text-cell-fabricated.csv");
  const fs::path text_cell_library = temp.temp_dir / (stem + "-text-cell.json");
  const fs::path binary_matrix = temp.temp_dir / (stem + "-binary-fabricated.csv");
  const fs::path binary_library = temp.temp_dir / (stem + "-binary.json");
  const fs::path header_surface = temp.temp_dir / (stem + "-header-only-surface.csv");
  const fs::path header_library = temp.temp_dir / (stem + "-header-only.json");
  const fs::path absent_thin_surface = temp.temp_dir / (stem + "-absent-thin-surface.csv");
  const fs::path thin_absent = temp.temp_dir / (stem + "-thin-surface-absent.json");
  const fs::path two_files = temp.temp_dir / (stem + "-two-files.json");
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
    const json library = json::parse(ReadBytes(convex_library_3d_path));
    CHECK(library.at("Models").size() == 2);
    json sa = library;
    for (auto &model : sa["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};  // the island's one target
    }
    json edge_only = sa, empty = sa, concave_only = sa;
    edge_only["Models"].erase(1);  // convex-corner-90
    empty["Models"] = json::array();
    concave_only["Models"].erase(0);
    concave_only["Models"][0]["Topology"] = "ConcaveCorner";
    std::ofstream(no_corner) << edge_only.dump(2) << "\n";
    std::ofstream(no_model) << empty.dump(2) << "\n";
    std::ofstream(no_match) << concave_only.dump(2) << "\n";
    // (d) / (e): the edge model's fabricated domain matrix with its lowest mode on the free
    // knots moved to -1e-7 x the largest |eigenvalue| there (beyond the operator's roundoff
    // 1e-9, which its matrix step refuses; a looser test on the preflight's side alone
    // would miss it), or a file that does not exist.
    {
      const json &edge = sa["Models"][0];
      std::ifstream input(edge["FabricatedMatrix"].get<std::string>());
      std::string header, line;
      std::getline(input, header);
      std::vector<std::tuple<int, int, double>> entries;
      int size = 0;
      while (std::getline(input, line))
      {
        std::stringstream fields(line);
        std::string i, j, q;
        std::getline(fields, i, ',');
        std::getline(fields, j, ',');
        std::getline(fields, q);
        entries.emplace_back(std::stoi(i) - 1, std::stoi(j) - 1, std::stod(q));
        size = std::max(size, std::stoi(j));
      }
      mfem::DenseMatrix matrix(size);
      matrix = 0.0;
      for (const auto &[i, j, q] : entries)
      {
        matrix(i, j) = q;
        matrix(j, i) = q;
      }
      const auto zero = edge.value("ZeroTraceIndices", std::vector<int>{});
      std::vector<int> free;
      for (int i = 0; i < size; i++)
      {
        if (std::find(zero.begin(), zero.end(), i) == zero.end())
        {
          free.push_back(i);
        }
      }
      const int n = static_cast<int>(free.size());
      Eigen::MatrixXd block(n, n);
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          block(i, j) = matrix(free[i], free[j]);
        }
      }
      const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensystem(block);
      REQUIRE(eigensystem.info() == Eigen::Success);
      const auto &values = eigensystem.eigenvalues();  // ascending
      const Eigen::VectorXd lowest = eigensystem.eigenvectors().col(0);
      const double shift = values(0) + 1.0e-7 * values.cwiseAbs().maxCoeff();
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          matrix(free[i], free[j]) -= shift * lowest(i) * lowest(j);
        }
      }
      std::ofstream output(negative_matrix);
      output << header << "\n" << std::setprecision(17);
      for (int i = 0; i < size; i++)
      {
        for (int j = i; j < size; j++)
        {
          output << i + 1 << "," << j + 1 << "," << matrix(i, j) << "\n";
        }
      }
      // STEP 4c (lane v1verdict1): the edge model's fabricated domain matrix that does not
      // fit its basis points (one basis function short), without its Q column, without
      // one upper-triangular entry.
      std::ofstream short_output(short_matrix), no_column_output(no_column_matrix),
          no_entry_output(no_entry_matrix);
      short_output << header << "\n" << std::setprecision(17);
      no_column_output << "basis_i,basis_j,Q (J)\n" << std::setprecision(17);
      no_entry_output << header << "\n" << std::setprecision(17);
      for (const auto &[i, j, q] : entries)
      {
        const std::string row = std::to_string(i + 1) + "," + std::to_string(j + 1) + ",";
        if (j + 1 < size)
        {
          short_output << row << q << "\n";
        }
        no_column_output << row << q << "\n";
        if (!(i == 0 && j == 1))
        {
          no_entry_output << row << q << "\n";
        }
      }
      // The fix round: the first Q cell not a number; 4 kB of bytes 0-255 (not text).
      std::ofstream text_cell_output(text_cell_matrix), binary_output(binary_matrix);
      text_cell_output << header << "\n" << std::setprecision(17);
      for (const auto &[i, j, q] : entries)
      {
        text_cell_output << i + 1 << "," << j + 1 << ",";
        if (&q == &std::get<2>(entries.front()))
        {
          text_cell_output << "abc\n";
        }
        else
        {
          text_cell_output << q << "\n";
        }
      }
      for (int k = 0; k < 4096; k++)
      {
        binary_output.put(static_cast<char>(k % 256));
      }
    }
    json negative = sa, unreadable = sa;
    negative["Models"][0]["FabricatedMatrix"] = negative_matrix.string();
    unreadable["Models"][0]["FabricatedMatrix"] = absent_matrix.string();
    std::ofstream(negative_mode) << negative.dump(2) << "\n";
    std::ofstream(no_file) << unreadable.dump(2) << "\n";
    for (const auto &[matrix, library] : {std::pair{short_matrix, short_library},
                                          std::pair{no_column_matrix, no_column_library},
                                          std::pair{no_entry_matrix, no_entry_library}})
    {
      json misfit = sa;
      misfit["Models"][0]["FabricatedMatrix"] = matrix.string();
      std::ofstream(library) << misfit.dump(2) << "\n";
    }
    // (f) / (g) / (h): every model's fabricated surface matrix or basis-point file absent;
    // every model without surface matrices, its interface mappings kept.
    json surface = sa, points = sa, mappings = sa;
    for (auto &model : surface["Models"])
    {
      model["FabricatedSurfaceMatrix"] = absent_surface.string();
    }
    for (auto &model : points["Models"])
    {
      model["BasisPoints"] = absent_points.string();
    }
    for (auto &model : mappings["Models"])
    {
      model.erase("FabricatedSurfaceMatrix");
      model.erase("ThinSurfaceMatrix");
    }
    std::ofstream(surface_absent) << surface.dump(2) << "\n";
    std::ofstream(points_absent) << points.dump(2) << "\n";
    std::ofstream(no_surface) << mappings.dump(2) << "\n";
    // The fix round: the edge model's fabricated domain matrix with a text cell or binary;
    // the edge model's fabricated surface matrix its header line only (it reads: no coupon
    // response); every model's thin surface matrix absent; the edge model's fabricated
    // domain matrix AND its thin surface matrix absent (one model, two refused files).
    json text_cell = sa, binary = sa, header_only = sa, thin = sa, two = sa;
    text_cell["Models"][0]["FabricatedMatrix"] = text_cell_matrix.string();
    binary["Models"][0]["FabricatedMatrix"] = binary_matrix.string();
    {
      std::ifstream input(sa["Models"][0]["FabricatedSurfaceMatrix"].get<std::string>());
      std::string header;
      std::getline(input, header);
      std::ofstream(header_surface) << header << "\n";
    }
    header_only["Models"][0]["FabricatedSurfaceMatrix"] = header_surface.string();
    for (auto &model : thin["Models"])
    {
      model["ThinSurfaceMatrix"] = absent_thin_surface.string();
    }
    two["Models"][0]["FabricatedMatrix"] = absent_matrix.string();
    two["Models"][0]["ThinSurfaceMatrix"] = absent_thin_surface.string();
    for (const auto &[library, data] :
         {std::pair{text_cell_library, &text_cell}, std::pair{binary_library, &binary},
          std::pair{header_library, &header_only}, std::pair{thin_absent, &thin},
          std::pair{two_files, &two}})
    {
      std::ofstream(library) << data->dump(2) << "\n";
    }
  }
  Mpi::Barrier(comm);
  json config = electrostatic ? ConvexIslandConfig(test::CrackPolicy::Default)
                              : ConvexMaxwellIslandConfig(test::CrackPolicy::Default);
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  if (!electrostatic)
  {
    config["Solver"]["Eigenmode"] = {{"Target", 0.05}, {"N", 2}};
  }
  const json::json_pointer correction_key(electrostatic
                                              ? "/Solver/Electrostatic/ResponseCorrection"
                                              : "/Solver/SurfaceResponseCorrection");
  config[correction_key].erase("PatchConstruction");
  config[correction_key].erase("Developer");
  config[correction_key]["CorrectionMode"] = "PostprocessOnly";
  config[correction_key]["PatchEnergy"] = true;
  config[correction_key]["PreflightOnly"] = true;
  struct Result
  {
    fs::path dir;
    int status = -1;
    std::string log;
  };
  // unmatched / geometry: "" = the key absent.
  auto Run = [&](const fs::path &library, const std::string &unmatched,
                 const std::string &geometry)
  {
    json run = config;
    Result result;
    result.dir = temp.temp_dir / (stem + "-" + library.stem().string() + "-u" + unmatched +
                                  "-g" + geometry);
    run["Problem"]["Output"] = result.dir.string();
    auto &request = run[correction_key];
    request["Library"] = library.string();
    request.erase("UnmatchedPolicy");
    if (!unmatched.empty())
    {
      request["UnmatchedPolicy"] = unmatched;
    }
    if (!geometry.empty())
    {
      request["GeometryFailurePolicy"] = geometry;
    }
    IoData iodata(run, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    result.log = test::CaptureStdout(
        [&]
        {
          result.status = RunSurfaceResponsePreflight(iodata, comm, utils::ConfigureOmp());
        });
    Mpi::Barrier(comm);
    return result;
  };
  auto Response = [](const fs::path &dir)
  {
    return json::parse(ReadBytes(dir / "palace.json"))
        .value("SurfaceResponse", json::object());
  };
  const std::string patches = "surface-response-resolution-patches.csv";
  auto NoCensus = [&](const Result &result, const std::string &line)
  {
    INFO("run " << result.dir.filename().string());
    CHECK(result.status == 0);
    if (!Mpi::Root(comm))
    {
      return;
    }
    CHECK(fs::is_regular_file(result.dir / "surface-response-requirements.json"));
    CHECK_THAT(result.log, ContainsSubstring(line));
    // The one field of every not-evaluated census (ruling (B)): the printed line.
    CHECK_THAT(
        "\n" + Response(result.dir).value("ResolutionCensusNotEvaluated", std::string()) +
            "\n",
        ContainsSubstring(line));
    CHECK_FALSE(Response(result.dir).contains("Resolution"));
    CHECK_FALSE(Response(result.dir).contains("ParticipationTable"));
    CHECK_FALSE(fs::exists(result.dir / patches));
    CHECK_FALSE(fs::exists(result.dir / kSurfaceParticipationFile));
  };
  // (a) UnmatchedPolicy Error with missing requirements: no census; Warn and absent: the
  // census.
  const auto error = Run(no_corner, "Error", "");
  NoCensus(error,
           "\nResolution census: not evaluated: UnmatchedPolicy \"Error\" with missing "
           "requirements (a solve-form run stops there); \"Warn\" gives it\n");
  if (Mpi::Root(comm))
  {
    const json manifest =
        json::parse(ReadBytes(error.dir / "surface-response-requirements.json"));
    CHECK(manifest.at("Complete") == false);
  }
  for (const std::string unmatched : {"Warn", ""})
  {
    const auto census = Run(no_corner, unmatched, "");
    INFO("UnmatchedPolicy \"" << unmatched << "\"");
    CHECK(census.status == 0);
    if (Mpi::Root(comm))
    {
      CHECK_THAT(census.log, !ContainsSubstring("Resolution census: not evaluated"));
      CHECK_FALSE(Response(census.dir).contains("ResolutionCensusNotEvaluated"));
      const json resolution = Response(census.dir).at("Resolution");
      CHECK(resolution.at("Knots").at("Located").get<long long int>() > 0);
      CHECK(resolution.at("NoCellsReason") == "PreflightOnly: no solve");
      CHECK(fs::is_regular_file(census.dir / patches));
    }
  }
  // (b) A library with no model: no census whatever the policy.
  NoCensus(Run(no_model, "", ""),
           "\nResolution census: not evaluated: the process library has no model (a "
           "solve-form run stops there)\n");
  // (c) An explicit FailClosed on a library that matches nothing; RecordRaw.
  NoCensus(Run(no_match, "", "FailClosed"),
           "geometry fallback(s) recorded by the response operator's construction (a "
           "solve-form run stops there); \"RecordRaw\" gives it\n");
  NoCensus(Run(no_match, "", "RecordRaw"),
           "\nResolution census: not evaluated: the correction falls back to the raw run "
           "(RawFallback: no patch applied)\n");
  // (d) A model whose fabricated domain matrix the operator's matrix step refuses:
  // predicted by the operator's own test (the supervisor's ruling (B)); no census, the line
  // naming the model and the operator's message, exit status 0; the solve form stops with
  // that message.
  const std::string negative_line =
      "\": Fabricated response matrix has a negative-energy mode beyond roundoff! (a "
      "solve-form run stops there)\n";
  for (const std::string unmatched : {"Warn", ""})
  {
    const auto refused = Run(negative_mode, unmatched, "");
    NoCensus(refused, negative_line);
    if (Mpi::Root(comm))
    {
      CHECK_THAT(refused.log,
                 ContainsSubstring("\nResolution census: not evaluated: response "
                                   "model \""));
    }
  }
  auto SolveFormStops = [&](const fs::path &library, const std::string &message)
  {
    json run = config;
    run["Problem"]["Output"] =
        (temp.temp_dir / (library.stem().string() + "-solve")).string();
    run[correction_key]["Library"] = library.string();
    run[correction_key]["PreflightOnly"] = false;
    IoData iodata(run, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    if (Mpi::Root(comm))
    {
      CHECK_THROWS_WITH(palace::Run(iodata, comm, utils::ConfigureOmp()),
                        ContainsSubstring(message));
    }
    else
    {
      CHECK_THROWS(palace::Run(iodata, comm, utils::ConfigureOmp()));
    }
    Mpi::Barrier(comm);
  };
  SolveFormStops(negative_mode,
                 "Fabricated response matrix has a negative-energy mode beyond roundoff!");
  // (e) A fabricated domain matrix file that does not exist (the solve form's matrix step
  // reads it): ReadTable's message.
  const std::string unreadable_message =
      "Unable to read response matrix file \"" + absent_matrix.string() + "\"!";
  NoCensus(Run(no_file, "", ""), unreadable_message + " (a solve-form run stops there)\n");
  SolveFormStops(no_file, unreadable_message);
  // (f)-(k) STEP 4c (lane v1verdict1; formerly the named exception): a library file the
  // operator's model loop refuses is predicted by the Try form of its reader: no census,
  // the line naming every refused model with the solve form's message, exit status 0; the
  // solve form stops with that message. (g) and (h) are the Maxwell configs of the second
  // step-4b review (on Electrostatic a missing basis-point file stops main's preflight
  // too).
  auto LibraryFileRefused = [&](const fs::path &library, const std::string &message,
                                int models, bool solve_form = true)
  {
    const auto refused = Run(library, "", "");
    NoCensus(refused, message + " (a solve-form run stops there)\n");
    if (Mpi::Root(comm))
    {
      // Every model the line names carries the message of its first refused file.
      const auto begin = refused.log.find("Resolution census: not evaluated: ");
      REQUIRE(begin != std::string::npos);
      const std::string line =
          refused.log.substr(begin, refused.log.find('\n', begin) - begin);
      auto Count = [&](const std::string &text)
      {
        int count = 0;
        for (auto at = line.find(text); at != std::string::npos;
             at = line.find(text, at + 1))
        {
          count++;
        }
        return count;
      };
      CHECK(Count("response model \"") == models);
      CHECK(Count("response model \"") == Count(message));
    }
    if (solve_form)
    {
      SolveFormStops(library, message);
    }
  };
  // (f) A surface-matrix file that cannot be read: ReadTable's message. (f), (g), (h) touch
  // both census models (the edge and the corner), (i)-(k) the edge model alone.
  LibraryFileRefused(
      surface_absent,
      "Unable to read response matrix file \"" + absent_surface.string() + "\"!", 2);
  // (i) / (j) / (k) A domain matrix that does not fit its basis points, without its Q
  // column, without an upper-triangular entry (the 30 sweep rows of the named exception
  // were (i)'s class).
  LibraryFileRefused(short_library,
                     "Response matrix \"" + short_matrix.string() +
                         "\" and basis point file have inconsistent sizes!",
                     1);
  LibraryFileRefused(no_column_library,
                     "Response matrix file \"" + no_column_matrix.string() +
                         "\" is missing column \"Q_ij (J)\"!",
                     1);
  LibraryFileRefused(no_entry_library,
                     "Response matrix file \"" + no_entry_matrix.string() +
                         "\" is missing an upper-triangular entry!",
                     1);
  // (l)-(p) The fix round. (l) / (m) A cell the Table parser refuses (a text cell, a binary
  // file): validated before parsing by the parser's own rules (a build without MFEM
  // exceptions cannot catch the parser's abort), the parser's message, printable.
  auto CellsMessage = [](const fs::path &path)
  {
    return "Response matrix file \"" + path.string() +
           "\": Could not parse CSV table from string! " +
           TableCellsError(ReadBytes(path)) + "!";
  };
  CHECK(TableCellsError(ReadBytes(text_cell_matrix)) ==
        "Could not parse CSV entry \"abc\" as double");
  const std::string binary_cells = TableCellsError(ReadBytes(binary_matrix));
  CHECK(binary_cells.rfind("Could not parse CSV entry \"\\x0d\\x0e", 0) == 0);
  CHECK(std::all_of(binary_cells.begin(), binary_cells.end(),
                    [](char c) { return c >= 0x20 && c < 0x7f; }));
  LibraryFileRefused(text_cell_library, CellsMessage(text_cell_matrix), 1);
  LibraryFileRefused(binary_library, CellsMessage(binary_matrix), 1);
  // (n) The edge model's fabricated surface matrix its header only: it reads, and the
  // coupon check after the surface matrices are read refuses it (a spatial model's reader
  // refuses such a file itself: no within-R rows).
  LibraryFileRefused(header_library,
                     "Response-correction interface mapping refers to a missing coupon "
                     "surface response!",
                     1);
  // (o) Every model's THIN surface matrix absent (MINOR-3: the prediction reads it too).
  LibraryFileRefused(
      thin_absent,
      "Unable to read response matrix file \"" + absent_thin_surface.string() + "\"!", 2);
  // (p) MINOR-4: one model, two refused files: both named, in the loop's order; the solve
  // form stops with the first.
  const std::string first = "Unable to read response matrix file \"" +
                            absent_matrix.string() + "\"!",
                    second = "Unable to read response matrix file \"" +
                             absent_thin_surface.string() + "\"!";
  LibraryFileRefused(two_files, first + " " + second, 1, false);
  SolveFormStops(two_files, first);
  // (q) MINOR-2: a model listed twice among the census models (as a Maxwell run lists a
  // runtime model's copies) is named once with its message.
  {
    json run = config;
    run[correction_key]["Library"] = no_file.string();
    IoData iodata(run, false);
    config::ElectrostaticSolverData::ResponseCorrectionModelData model;
    model.name = "listed-twice";
    model.basis_points = absent_points.string();
    SurfaceResponsePreflightOutcome outcome;
    outcome.census_models = {model, model};
    const std::string line = SurfaceResponsePreflightModelRefusal(outcome, iodata, comm);
    const std::string named = "response model \"listed-twice\": Unable to open ";
    CHECK(line.find(named) != std::string::npos);
    CHECK(line.find(named) == line.rfind(named));
  }
  if (!electrostatic)
  {
    // (g) A basis-point file that cannot be opened: the prediction on the root at any rank
    // count (the matrices' size is then unknown: not checked); the solve form serially only
    // (in parallel one rank's reader can stop while the others wait in a collective, which
    // a test program cannot observe and go on).
    const std::string points_message =
        "Unable to open response-correction basis point file \"" + absent_points.string() +
        "\"!";
    LibraryFileRefused(points_absent,
                       points_message +
                           " (its matrices are not checked without its basis points)",
                       2, false);
    if (Mpi::Size(comm) == 1)
    {
      SolveFormStops(points_absent, points_message);
    }
    // (h) Interface mappings without surface matrices (ComputeSurfaceResponseMatrices).
    LibraryFileRefused(no_surface,
                       "Response-correction interface mappings require surface "
                       "response matrices!",
                       2);
  }
#endif
}

// CUT-PLAN step 4b (lane v1step4b; the supervisor's ruling R4): on the three-dimensional
// Features path the census operator reads the manifest pass's identification instead of
// computing it again. The census is bit-equal with and without the sharing
// (PALACE_PREFLIGHT_SHARED_IDENTIFICATION=0) and equal to the solve form's level 1; after
// an edge-tube refinement of the first mesh the operator identifies that mesh itself.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "PreflightOnly's census on the manifest's identification equals the "
                 "census of its own identification",
                 "[surfaceresponseoperator][preflightonly][3d][Serial][Parallel]")
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
                                         {"Driven", false},
                                         {"Driven", true}}));
  INFO("problem " << problem << ", edge tubes " << edge_tubes);
  MPI_Comm comm = Mpi::World();
  const bool electrostatic = problem == "Electrostatic";
  const std::string stem = "preflight-shared-" + problem + (edge_tubes ? "-tubes" : "");
  const fs::path mesh_path = temp.temp_dir / (stem + ".mesh");
  const fs::path sa_library = temp.temp_dir / (stem + "-sa.json");
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
    json library = json::parse(ReadBytes(convex_library_3d_path));
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
    }
    std::ofstream(sa_library) << library.dump(2) << "\n";
  }
  Mpi::Barrier(comm);
  json config = electrostatic ? ConvexIslandConfig(test::CrackPolicy::Default)
                              : ConvexMaxwellIslandConfig(test::CrackPolicy::Default);
  config["Problem"]["Type"] = problem;
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  const json::json_pointer correction_key(electrostatic
                                              ? "/Solver/Electrostatic/ResponseCorrection"
                                              : "/Solver/SurfaceResponseCorrection");
  auto &correction = config[correction_key];
  correction.erase("PatchConstruction");
  correction.erase("Developer");
  correction["Library"] = sa_library.string();
  correction["CorrectionMode"] = "PostprocessOnly";
  correction["PatchEnergy"] = true;
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
  // mode: "solve", "shared" (the default), "own" (the sharing switched off).
  auto Run = [&](const std::string &mode)
  {
    json run = config;
    const fs::path output = temp.temp_dir / (stem + "-" + mode);
    run["Problem"]["Output"] = output.string();
    run[correction_key]["PreflightOnly"] = mode != "solve";
    IoData iodata(run, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    ScopedEnvironmentVariable share("PALACE_PREFLIGHT_SHARED_IDENTIFICATION",
                                    mode == "own" ? "0" : nullptr);
    const auto log = test::CaptureStdout(
        [&]
        {
          if (mode == "solve")
          {
            palace::Run(iodata, comm, utils::ConfigureOmp());
          }
          else
          {
            CHECK(RunSurfaceResponsePreflight(iodata, comm, utils::ConfigureOmp()) == 0);
          }
        });
    Mpi::Barrier(comm);
    return std::pair{output, log};
  };
  const auto [solve, solve_log] = Run("solve");
  const auto [shared, shared_log] = Run("shared");
  const auto [own, own_log] = Run("own");
  if (!Mpi::Root(comm))
  {
    return;
  }
  // The census lines of a log (the knots line and the PREFLIGHT line).
  auto CensusLines = [](const std::string &log)
  {
    std::vector<std::string> lines;
    std::istringstream stream(log);
    std::string line;
    while (std::getline(stream, line))
    {
      if (line.rfind(" Resolution census: ", 0) == 0 ||
          line.rfind("Resolution census PREFLIGHT", 0) == 0)
      {
        lines.push_back(line);
      }
    }
    return lines;
  };
  auto Count = [](const std::string &log, const std::string &text)
  {
    std::size_t count = 0;
    for (auto at = log.find(text); at != std::string::npos; at = log.find(text, at + 1))
    {
      count++;
    }
    return count;
  };
  // One identification with the sharing, two without it or after an edge-tube refinement.
  CHECK(Count(shared_log, "Geometry identification: ") == (edge_tubes ? 2 : 1));
  CHECK(Count(own_log, "Geometry identification: ") == 2);
  REQUIRE(CensusLines(shared_log).size() == 2);
  CHECK(CensusLines(shared_log) == CensusLines(own_log));
  CHECK(CensusLines(shared_log) == CensusLines(solve_log));
  const std::string patches = "surface-response-resolution-patches.csv";
  REQUIRE(fs::is_regular_file(shared / patches));
  CHECK(ReadBytes(shared / patches) == ReadBytes(own / patches));
  CHECK(ReadBytes(shared / patches) == ReadBytes(solve / patches));
  CHECK(ReadBytes(shared / kSurfaceParticipationFile) ==
        ReadBytes(own / kSurfaceParticipationFile));
  auto Block = [](const fs::path &dir)
  {
    return json::parse(ReadBytes(dir / "palace.json"))
        .at("SurfaceResponse")
        .at("Resolution");
  };
  CHECK(Block(shared) == Block(own));
  for (const char *member : kPreSolveResolutionMembers)
  {
    INFO("SurfaceResponse.Resolution." << member);
    CHECK(Block(shared).at(member) == Block(solve).at(member));
  }
#endif
}

// The fix round of lane v1verdict1 (MINOR-3 and BLOCKER-1 on a real mesh): a
// Resolution.BarsFile whose Verdict moves the Order-2 bar moves the PREFLIGHT word and the
// table's first-iteration word together (both from the one verdict data), and the flag of
// cells beyond twice the bar follows LargestM against the limit; where it is not raised,
// no table row reads UNDER-RESOLVED after the solve.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "PreflightOnly's word and flag follow a BarsFile's Verdict; without the flag "
    "the table reads no UNDER-RESOLVED",
    "[surfaceresponseoperator][preflightonly][surfaceparticipation][3d][Serial]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "preflight-bars.mesh";
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  json config = ConvexIslandConfig(test::CrackPolicy::Default);
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  config["Solver"]["Order"] = 2;
  auto &request = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  request["CorrectionMode"] = "PostprocessOnly";
  request["PatchEnergy"] = true;
  struct Words
  {
    json preflight;          // the SA all-planes Preflight entry (PreflightOnly)
    std::string flags;       // its PreflightOnly SA row's flags
    std::string table;       // the solved SA Type row's resolution
    double table_m90 = 0.0;  // and its h/R m90
    double radius = 0.0;     // the run's R (m)
  };
  auto Run = [&](const std::string &label, const fs::path &bars)
  {
    Words words;
    for (const bool preflight : {true, false})
    {
      json run = config;
      const fs::path output = temp.temp_dir / ("preflight-bars-" + label +
                                               (preflight ? "-preflight" : "-solve"));
      run["Problem"]["Output"] = output.string();
      run["Solver"]["Electrostatic"]["ResponseCorrection"]["PreflightOnly"] = preflight;
      if (!bars.empty())
      {
        run["Solver"]["Electrostatic"]["ResponseCorrection"]["Resolution"] = {
            {"BarsFile", bars.string()}};
      }
      IoData iodata(run, false);
      MakeOutputFolder(iodata, comm);
      BlockTimer::Reset();
      if (preflight)
      {
        REQUIRE(RunSurfaceResponsePreflight(iodata, comm, utils::ConfigureOmp()) == 0);
      }
      else
      {
        palace::Run(iodata, comm, utils::ConfigureOmp());
      }
      Mpi::Barrier(comm);
      if (!Mpi::Root(comm))
      {
        continue;
      }
      const auto rows = CsvCells(output / kSurfaceParticipationFile);
      REQUIRE(rows.size() > 1);
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
      for (std::size_t r = 1; r < rows.size(); r++)
      {
        if (rows[r][Index("Type")] == "SA" && rows[r][Index("scope")] == "Type")
        {
          (preflight ? words.flags : words.table) =
              rows[r][Index(preflight ? "flags" : "resolution")];
          words.table_m90 =
              preflight ? words.table_m90 : std::stod(rows[r][Index("h/R m90")]);
        }
      }
      const json census = json::parse(ReadBytes(output / "palace.json"))
                              .at("SurfaceResponse")
                              .at("Resolution");
      words.radius = census.at("Radius").get<double>();
      for (const auto &row : census.at("Preflight"))
      {
        if (preflight && row.at("Type") == "SA" && row.at("Plane").is_null())
        {
          words.preflight = row;
        }
      }
    }
    return words;
  };
  const Words built_in = Run("built-in", {});
  REQUIRE(built_in.preflight.at("LengthM90").is_number());
  const double length_m90 = built_in.preflight.at("LengthM90").get<double>();
  const double largest = built_in.preflight.at("LargestM").get<double>();
  REQUIRE(largest >= length_m90);
  // A bars file: the run's radius, census bars (their shares only), the built-in verdict
  // with the Order-2 bar B.
  auto Bars = [&](const std::string &label, double bar)
  {
    const fs::path path = temp.temp_dir / ("preflight-bars-" + label + ".json");
    json census_bars = json::object();
    for (const char *order : {"1", "2", "3", "4"})
    {
      census_bars[order] = {{"b5", nullptr}, {"b10", 1.0}};
    }
    // The built-in verdict's values (palace/models/surfaceparticipationverdict.json) with
    // the Order-2 bar B.
    const json verdict = {
        {"MinUnknownsRatio", 1.2},
        {"MaxCensusRatio", 0.95},
        {"MaxChangePerDoublingPercent", 1.0},
        {"CensusBarOverR", {{"1", nullptr}, {"2", bar}, {"3", 2.0}, {"4", 4.0}}},
        {"UnderResolvedFactor", 2.0},
        {"TestedProblems", {"Electrostatic"}},
        {"TestedScopes", {"Type"}},
        {"Interim", false},
        {"Provenance", "test-preflightonly: a moved Order-2 bar"}};
    std::ofstream(path) << json{{"RadiusMicrons", built_in.radius * 1.0e6},
                                {"Bars", census_bars},
                                {"Verdict", verdict}}
                               .dump(2)
                        << "\n";
    return path;
  };
  const std::string flag = kParticipationPreflightUnderResolvedCellsFlag;
  // (1) A bar so low that twice it lies below an eighth of the length m90: the PREFLIGHT
  // word and the table's first-iteration word read UNDER-RESOLVED together, flagged.
  const Words low = Run("low", Bars("low", length_m90 / 16.0));
  CHECK(low.preflight.at("Status") == "UNDER-RESOLVED");
  CHECK(low.table == "UNDER-RESOLVED");
  CHECK(low.preflight.at("UnderResolvedLengthShare").get<double>() > 0.0);
  CHECK(low.flags.find(flag) != std::string::npos);
  // (2) A bar at LargestM (twice it above every cell): both MARGINAL, NO flag (share 0),
  // and the table reads no UNDER-RESOLVED (its m90 within twice the bar).
  const Words high = Run("high", Bars("high", largest));
  CHECK(high.preflight.at("Status") == "MARGINAL");
  CHECK(high.table == "MARGINAL");
  CHECK(high.preflight.at("UnderResolvedLengthShare").get<double>() == 0.0);
  CHECK(high.flags.find(flag) == std::string::npos);
  CHECK(high.table_m90 <= 2.0 * largest);
  // (3) Twice the bar just below LargestM: flagged, LargestM above the limit; the table's
  // h/R m90 is at most LargestM (the PREFLIGHT reads every cell the table weighs).
  const Words edge = Run("edge", Bars("edge", 0.49 * largest));
  CHECK(edge.preflight.at("UnderResolvedLengthShare").get<double>() > 0.0);
  CHECK(edge.preflight.at("LargestM").get<double>() > 2.0 * 0.49 * largest);
  CHECK(edge.flags.find(flag) != std::string::npos);
  for (const auto *words : {&built_in, &low, &high, &edge})
  {
    CHECK(words->table_m90 <= largest * (1.0 + 1.0e-5));  // the printed m90's rounding
    CHECK(words->preflight.at("LargestM") == built_in.preflight.at("LargestM"));
  }
#endif
}
