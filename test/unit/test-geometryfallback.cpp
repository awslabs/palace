// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "models/laplaceoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"

using namespace palace;
using namespace Catch::Matchers;
using json = nlohmann::json;
namespace fs = std::filesystem;

// Geometry fallbacks (readiness hardening; USER decision 1000 (b), supervisor 1001; DESIGN
// readiness-20261029/hardening sections 1-5). Every case below ABORTS on main's production
// sources at the converted site (the failing-first record) and RECORDS on the branch.

namespace
{

json ReadJson(const fs::path &path)
{
  std::ifstream input(path);
  REQUIRE(input);
  return json::parse(input);
}

}  // namespace

// H1 (operator 13058 / 16386 "no usable correction patches" + the constructor's "at least
// one model and patch"): a valid library whose keys match no feature of the device. The
// island mesh (isolated edges + convex corners) against the fixture's cluster-keyed library
// (every feature Missing under UnmatchedPolicy Warn). RecordRaw (the solve default): the
// operator constructs flagged IsRawFallback with ONE NoMatchedFeature entry, every type's
// RawShareByEnergy 1 (no correction applied); the explicit FailClosed leaf aborts with the
// class named; the preflight (FailClosed default) records nothing: a Missing-only manifest
// is the normal discovery state and predicts the solve's raw run by Summary.Counts (Exact +
// Interpolated + LegacyContract + QuantumNearMatched == 0).
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "Geometry fallback H1 NoMatchedFeature: a library matching no feature runs "
    "raw with the record (decisions 1000 (b) / 1001)",
    "[surfaceresponseoperator][geometryfallback][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  json config = ConvexIslandConfig();
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction["Library"] = spatial_cluster_library_3d_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  correction.erase("PatchConstruction");
  correction["CorrectionMode"] = "PostprocessOnly";
  auto MakeIoData = [&](const json &c)
  {
    IoData iodata(c, false);
    iodata.boundaries.cracked_attributes.insert(9);
    return iodata;
  };

  SECTION("RecordRaw (the solve default): the operator is a recorded raw fallback")
  {
    IoData iodata = MakeIoData(config);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    CHECK(response.IsRawFallback());
    REQUIRE(response.GetGeometryFallbacks().size() == 1);
    const auto &entry = response.GetGeometryFallbacks().front();
    CHECK(entry.cls == "NoMatchedFeature");
    CHECK(entry.what == "run");
    CHECK(entry.site == "H1/BuildAutomaticResponseData3D:features");
    REQUIRE(entry.interfaces.size() == 1);
    CHECK(entry.interfaces.front() == 4);
    CHECK(entry.length > 0.0);
    const json record = response.GetGeometryFallbackRecord();
    REQUIRE(!record.is_null());
    CHECK(record.at("Version").get<int>() == 1);
    CHECK(record.at("Policy") == "RecordRaw");
    CHECK(record.at("Run") == "Solve");
    CHECK(record.at("Count").get<int>() == 1);
    CHECK(record.at("Classes").at("NoMatchedFeature").get<int>() == 1);
    REQUIRE(record.at("RawShareByLength").contains("SA"));
    CHECK_THAT(record.at("RawShareByLength").at("SA").at("Share").get<double>(),
               WithinAbs(1.0, 1.0e-12));
    REQUIRE(record.at("RawShareByEnergy").contains("SA"));
    CHECK(record.at("RawShareByEnergy").at("SA").at("Share").get<double>() == 1.0);
    // The statistics record the drivers write into palace.json carries the entry.
    const json statistics = response.GetStatistics();
    CHECK(statistics.at("Diagnostics").at("GeometryFallbacks").at("Count").get<int>() == 1);
    // The banner names the class and the 100 % raw share.
    const std::string banner = GeometryFallbackBanner(record, "palace.json");
    CHECK_THAT(banner, ContainsSubstring("LOW CONFIDENCE") &&
                           ContainsSubstring("NoMatchedFeature 1") &&
                           ContainsSubstring("SA 100 % (no correction applied)"));
  }

  SECTION("the explicit FailClosed leaf aborts at the end of the configuration, naming "
          "the class")
  {
    json fail_closed = config;
    fail_closed["Solver"]["Electrostatic"]["ResponseCorrection"]["GeometryFailurePolicy"] =
        "FailClosed";
    IoData iodata = MakeIoData(fail_closed);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    LaplaceOperator laplace(iodata, meshes);
    CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, laplace),
                      ContainsSubstring("GeometryFailurePolicy FailClosed") &&
                          ContainsSubstring("NoMatchedFeature 1") &&
                          ContainsSubstring("1 geometry fallback(s)"));
  }

  SECTION("the leaf is parsed, echoed only when set, and refuses other values")
  {
    IoData plain = MakeIoData(config);
    CHECK(!plain.solver.electrostatic.response_correction->geometry_failure_policy);
    json record_raw = config;
    record_raw["Solver"]["Electrostatic"]["ResponseCorrection"]["GeometryFailurePolicy"] =
        "RecordRaw";
    IoData explicit_raw = MakeIoData(record_raw);
    REQUIRE(explicit_raw.solver.electrostatic.response_correction->geometry_failure_policy);
    CHECK(*explicit_raw.solver.electrostatic.response_correction->geometry_failure_policy ==
          GeometryFailurePolicy::RECORD_RAW);
    CHECK(ResolveGeometryFailurePolicy(*plain.solver.electrostatic.response_correction,
                                       false) == GeometryFailurePolicy::RECORD_RAW);
    CHECK(ResolveGeometryFailurePolicy(*plain.solver.electrostatic.response_correction,
                                       true) == GeometryFailurePolicy::FAIL_CLOSED);
    CHECK(
        ResolveGeometryFailurePolicy(*explicit_raw.solver.electrostatic.response_correction,
                                     true) == GeometryFailurePolicy::RECORD_RAW);
    json invalid = config;
    invalid["Solver"]["Electrostatic"]["ResponseCorrection"]["GeometryFailurePolicy"] =
        "Ignore";
    CHECK_THROWS_WITH(IoData(invalid, false), ContainsSubstring("GeometryFailurePolicy") &&
                                                  ContainsSubstring("RecordRaw") &&
                                                  ContainsSubstring("FailClosed"));
  }

  SECTION("the preflight predicts the raw run by its Counts and records no fallback")
  {
    const fs::path dir = temp.temp_dir / "h1-preflight";
    json preflight_config = config;
    preflight_config["Problem"]["Output"] = dir.string();
    IoData iodata = MakeIoData(preflight_config);
    MakeOutputFolder(iodata, comm);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    const fs::path manifest_path = dir / "surface-response-requirements.json";
    WriteSurfaceResponseRequirements(iodata, *meshes.back(), manifest_path.string());
    Mpi::Barrier(comm);
    const json manifest = ReadJson(manifest_path);
    const auto &counts = manifest.at("Summary").at("Counts");
    CHECK(counts.at("Exact").get<int>() + counts.at("Interpolated").get<int>() +
              counts.value("LegacyContract", 0) + counts.value("QuantumNearMatched", 0) ==
          0);
    CHECK(counts.at("Missing").get<int>() > 0);
    CHECK(!manifest.at("Summary").contains("GeometryFallbacks"));
    CHECK(!manifest.at("Identification").at("Diagnostics").contains("GeometryFallbacks"));
  }
#endif
}

// H1 through the electrostatic driver: the run completes; every raw output (surface-Q,
// domain-E, terminal-C, surface-Q-edge) is BYTE-IDENTICAL to the same configuration without
// ResponseCorrection; palace.json carries SurfaceResponse.Diagnostics.GeometryFallbacks and
// no corrected CSV / patch table is written.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Geometry fallback H1 NoMatchedFeature: the electrostatic run is the thin "
                 "solve plus the record (decisions 1000 (b) / 1001)",
                 "[surfaceresponseoperator][geometryfallback][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "h1-island.mesh";
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(comm);
  json config = ConvexIslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Solver"]["Linear"] = {{"Tol", 1.0e-10}, {"MaxIts", 400}};
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction["Library"] = spatial_cluster_library_3d_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  correction.erase("PatchConstruction");
  correction["CorrectionMode"] = "Both";
  const fs::path fallback_dir = temp.temp_dir / "h1-run-fallback";
  const fs::path plain_dir = temp.temp_dir / "h1-run-plain";
  config["Problem"]["Output"] = fallback_dir.string();
  test::RunElectrostatic(config);
  json plain = config;
  plain["Problem"]["Output"] = plain_dir.string();
  plain["Solver"]["Electrostatic"].erase("ResponseCorrection");
  test::RunElectrostatic(plain);
  Mpi::Barrier(comm);
  if (Mpi::Root(comm))
  {
    for (const char *name : {"surface-Q.csv", "domain-E.csv", "terminal-C.csv",
                             "surface-Q-edge.csv", "terminal-V.csv"})
    {
      INFO(name);
      REQUIRE(fs::is_regular_file(fallback_dir / name));
      CHECK(test::ReadFile(fallback_dir / name) == test::ReadFile(plain_dir / name));
    }
    for (const char *name :
         {"surface-Q-corrected.csv", "surface-response-patches.csv",
          "surface-response-model-energy.csv", "surface-response-geometry-fallbacks.csv"})
    {
      INFO(name);
      CHECK(!fs::exists(fallback_dir / name));
    }
    const json metadata = ReadJson(fallback_dir / "palace.json");
    const auto &record =
        metadata.at("SurfaceResponse").at("Diagnostics").at("GeometryFallbacks");
    CHECK(record.at("Count").get<int>() == 1);
    CHECK(record.at("Classes").at("NoMatchedFeature").get<int>() == 1);
    CHECK(record.at("Policy") == "RecordRaw");
    CHECK(record.at("RawShareByEnergy").at("SA").at("Share").get<double>() == 1.0);
    CHECK(!ReadJson(plain_dir / "palace.json").contains("SurfaceResponse"));
  }
#endif
}

// H8 (operator FindDomainBoundaryExclusions "leaves no applied surface-response patch"): a
// mis-scaled coupon (basis points 2.0 beyond the unit island domain) alone. Without a
// ledger (the direct callers of record) the exclusion fails closed as before; with the
// ledger it records AllPatchesOutside (every tested patch, weight 0) and returns.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Geometry fallback H8 AllPatchesOutside: every coupon outside the mesh is "
                 "recorded instead of aborting (decisions 1000 (b) / 1001)",
                 "[surfaceresponseoperator][geometryfallback][domainboundary][3d][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  using ResponsePatchData = config::ElectrostaticSolverData::ResponseCorrectionPatchData;
  constexpr double R = 0.2;
  auto mesh = MakeIslandMesh();
  const std::vector<std::array<double, 3>> wide_basis_points = {
      {-2.0, -2.0, 0.0}, {2.0, -2.0, 0.0}, {2.0, 2.0, 0.0}, {-2.0, 2.0, 0.0}};
  auto BasisPoints = [&](int) -> const std::vector<std::array<double, 3>> *
  { return &wide_basis_points; };
  auto Spatial = [](int) { return false; };
  auto Name = [](int) { return std::string("wide"); };
  auto Patch = [&](const std::array<double, 3> &origin)
  {
    ResponsePatchData patch;
    patch.model = 1;
    patch.origin = origin;
    patch.axis_u = {-1.0, 0.0, 0.0};
    patch.axis_v = {0.0, 1.0, 0.0};
    patch.axis_w = {0.0, 0.0, 1.0};
    patch.longitudinal_cell = {-0.03, 0.03};
    return patch;
  };
  std::vector<ResponsePatchData> patches = {Patch({0.5, 0.5, 0.5}), Patch({0.5, 0.5, 0.3})};
  // The direct call of record (no ledger) keeps the abort.
  {
    auto copy = patches;
    CHECK_THROWS_WITH(
        FindDomainBoundaryExclusions(*mesh, copy, BasisPoints, Spatial, Name, 1.0, R, {}),
        ContainsSubstring("leaves no applied surface-response patch"));
  }
  // With the ledger: the record, every patch weight 0, no abort.
  {
    std::vector<GeometryFallbackData> entries;
    GeometryFallbackLedger ledger(GeometryFailurePolicy::RECORD_RAW, false, entries);
    auto copy = patches;
    const auto exclusions = FindDomainBoundaryExclusions(
        *mesh, copy, BasisPoints, Spatial, Name, 1.0, R, {}, {}, 0.0, {}, &ledger);
    CHECK(exclusions.tested_patches == 2);
    CHECK(exclusions.patches.size() == 2);
    CHECK(copy[0].weight == 0.0);
    CHECK(copy[1].weight == 0.0);
    REQUIRE(entries.size() == 1);
    CHECK(entries.front().cls == "AllPatchesOutside");
    CHECK(entries.front().what == "run");
    CHECK(entries.front().patches == std::vector<int>{0, 1});
    CHECK_THAT(entries.front().reason, ContainsSubstring("all 2 tested patches"));
    // The deferred FailClosed abort names the class; RecordRaw's Enforce is a no-op.
    ledger.Enforce("test");
    GeometryFallbackLedger fail_closed(GeometryFailurePolicy::FAIL_CLOSED, false, entries);
    CHECK_THROWS_WITH(fail_closed.Enforce("test"),
                      ContainsSubstring("FailClosed") &&
                          ContainsSubstring("AllPatchesOutside 1") &&
                          ContainsSubstring("(test)"));
  }
#endif
}
