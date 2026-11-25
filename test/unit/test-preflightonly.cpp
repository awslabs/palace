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
// construction only, ends it there; RecordRaw gives the raw-fallback line). (f)-(h) pin
// the class the supervisor ruled into the named exception (STEP 4c, a non-aborting reader
// of the model files): a library file the operator's model loop refuses stops
// PreflightOnly after the manifest is written, with the solve form's message.
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
    }
    json negative = sa, unreadable = sa;
    negative["Models"][0]["FabricatedMatrix"] = negative_matrix.string();
    unreadable["Models"][0]["FabricatedMatrix"] = absent_matrix.string();
    std::ofstream(negative_mode) << negative.dump(2) << "\n";
    std::ofstream(no_file) << unreadable.dump(2) << "\n";
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
  // (f)-(h) The class ruled into the named exception (STEP 4c): PreflightOnly stops after
  // the manifest is written, before any census, with the solve form's message (an abort
  // with a non-zero exit status in a production run; a preflight without the census ended
  // with 0). These rows pin that behaviour; (g) and (h) are the Maxwell configs of the
  // second review (on Electrostatic a missing basis-point file stops main's preflight too).
  auto PreflightStops = [&](const fs::path &library, const std::string &message)
  {
    json run = config;
    const fs::path dir = temp.temp_dir / (stem + "-" + library.stem().string() + "-stop");
    run["Problem"]["Output"] = dir.string();
    run[correction_key]["Library"] = library.string();
    run[correction_key].erase("UnmatchedPolicy");
    IoData iodata(run, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    bool thrown = false;
    std::string what;
    test::CaptureStdout(
        [&]
        {
          try
          {
            RunSurfaceResponsePreflight(iodata, comm, utils::ConfigureOmp());
          }
          catch (const std::exception &e)
          {
            thrown = true;
            what = e.what();
          }
        });
    Mpi::Barrier(comm);
    INFO("run " << dir.filename().string());
    CHECK(thrown);
    if (Mpi::Root(comm))
    {
      CHECK_THAT(what, ContainsSubstring(message));
      CHECK(fs::is_regular_file(dir / "surface-response-requirements.json"));
      CHECK_FALSE(fs::exists(dir / patches));
      CHECK_FALSE(fs::exists(dir / kSurfaceParticipationFile));
    }
    SolveFormStops(library, message);
  };
  // (f) A surface-matrix file that cannot be read: ReadTable's message.
  PreflightStops(surface_absent, "Unable to read response matrix file \"" +
                                     absent_surface.string() + "\"!");
  if (!electrostatic)
  {
    // (g) A basis-point file that cannot be opened. Serial only: in parallel one rank's
    // reader stops while the others wait in a collective (a production run ends there, the
    // stopped rank's abort ending the job), which a test program cannot observe and go on.
    if (Mpi::Size(comm) == 1)
    {
      PreflightStops(points_absent, "Unable to open response-correction basis point file");
    }
    // (h) Interface mappings without surface matrices (ComputeSurfaceResponseMatrices).
    PreflightStops(no_surface,
                   "Response-correction interface mappings require surface response "
                   "matrices!");
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
