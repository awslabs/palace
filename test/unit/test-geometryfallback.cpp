// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "models/laplaceoperator.hpp"
#include "models/surfaceresponseidentification.hpp"
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

// F-6 (readiness hardening-2, decision 1018 Q4; operator FindDomainBoundaryExclusions
// "misplaced or mis-scaled coupon": the census QUO-019 strip islands whose metal-edge
// reference lies 1 nm outside the tile wall). On the unit island mesh ([0, 1]^3) a cell
// whose reference sits on the x = 1 face: pushed 1e-9 outside (within the boundary
// tolerance kSignatureParameterToleranceOverRadius x R = 2e-4) it is read as located and
// the patch keeps the classification of the on-face one (DomainBoundary MirrorOff by its
// gap-side points; no record, with or without a ledger); pushed 1e-2 R outside, without a
// ledger the abort of record stands, with the ledger the patch is excluded with the reason
// MisplacedReference and one UnlocatedPatch record (what "patch", its model, feature,
// cell length and the reference point). A cell with no reference and every tested point
// outside: the abort of record, or AllPointsOutside + UnlocatedPatch under the ledger. The
// preflight and the solve call this one function with their ledgers (P == S).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Geometry fallback F-6 UnlocatedPatch: a metal-edge reference on the mesh "
                 "boundary to roundoff is located, one off the mesh is recorded (decision "
                 "1018 Q4)",
                 "[surfaceresponseoperator][geometryfallback][domainboundary][3d][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  using ResponsePatchData = config::ElectrostaticSolverData::ResponseCorrectionPatchData;
  constexpr double R = 0.2;
  const double tolerance = kSignatureParameterToleranceOverRadius * R;
  auto mesh = MakeIslandMesh();
  const std::vector<std::array<double, 3>> basis_points = {
      {-0.05, -0.05, 0.0}, {0.05, -0.05, 0.0}, {0.05, 0.05, 0.0}, {-0.05, 0.05, 0.0}};
  auto BasisPoints = [&](int) -> const std::vector<std::array<double, 3>> *
  { return &basis_points; };
  auto Spatial = [](int) { return false; };
  auto Name = [](int) { return std::string("edge"); };
  // A cell in the lead's frame (u from the metal into the gap = -x, w along the edge)
  // whose reference (the origin) sits at x = 1 + delta: the gap-side basis points at
  // x = 1.05 + delta are outside the unit mesh whatever delta.
  auto Patch = [&](const std::array<double, 3> &origin, int feature)
  {
    ResponsePatchData patch;
    patch.model = 0;
    patch.origin = origin;
    patch.axis_u = {-1.0, 0.0, 0.0};
    patch.axis_v = {0.0, 1.0, 0.0};
    patch.axis_w = {0.0, 0.0, 1.0};
    patch.longitudinal_cell = {-0.03, 0.03};
    patch.provenance.feature = feature;
    patch.provenance.coupon_depth = R;
    return patch;
  };
  const ResponsePatchData inside = Patch({0.25, 0.5, 0.3}, 3);
  auto Find = [&](std::vector<ResponsePatchData> &patches, GeometryFallbackLedger *ledger)
  {
    return FindDomainBoundaryExclusions(*mesh, patches, BasisPoints, Spatial, Name, 1.0, R,
                                        {}, {}, 0.0, {}, ledger);
  };
  // (1) The reference ON the face and 1e-9 outside it: the same classification
  // (DomainBoundary MirrorOff by the gap-side points), no record, no abort.
  for (const double delta : {0.0, 1.0e-9, 0.5 * tolerance})
  {
    for (const bool with_ledger : {false, true})
    {
      std::vector<GeometryFallbackData> entries;
      GeometryFallbackLedger ledger(GeometryFailurePolicy::RECORD_RAW, false, entries);
      std::vector<ResponsePatchData> patches = {inside, Patch({1.0 + delta, 0.5, 0.3}, 7)};
      const auto exclusions = Find(patches, with_ledger ? &ledger : nullptr);
      CAPTURE(delta, with_ledger);
      REQUIRE(exclusions.tested_patches == 2);
      REQUIRE(exclusions.patches.size() == 1);
      CHECK(exclusions.patches[0].patch == 1);
      CHECK(exclusions.patches[0].reason == "MirrorOff");
      CHECK(exclusions.patches[0].outside_points > 0);
      CHECK(exclusions.patches[0].outside_points < exclusions.patches[0].tested_points);
      CHECK(patches[0].weight == 1.0);
      CHECK(patches[1].weight == 0.0);
      CHECK(entries.empty());
    }
  }
  // (2) The reference 1e-2 R outside: beyond the tolerance.
  {
    const double delta = 1.0e-2 * R;
    std::vector<ResponsePatchData> patches = {inside, Patch({1.0 + delta, 0.5, 0.3}, 7)};
    // The direct call of record (no ledger) keeps the abort.
    {
      auto copy = patches;
      CHECK_THROWS_WITH(Find(copy, nullptr),
                        ContainsSubstring("misplaced or mis-scaled coupon") &&
                            ContainsSubstring("patch 1 (0-based") &&
                            ContainsSubstring("model edge") &&
                            ContainsSubstring("metal-edge reference"));
    }
    // With the ledger: excluded MisplacedReference, one UnlocatedPatch record.
    {
      std::vector<GeometryFallbackData> entries;
      GeometryFallbackLedger ledger(GeometryFailurePolicy::RECORD_RAW, false, entries);
      auto copy = patches;
      const auto exclusions = Find(copy, &ledger);
      REQUIRE(exclusions.tested_patches == 2);
      REQUIRE(exclusions.patches.size() == 1);
      CHECK(exclusions.patches[0].patch == 1);
      CHECK(exclusions.patches[0].reason == "MisplacedReference");
      CHECK(copy[0].weight == 1.0);
      CHECK(copy[1].weight == 0.0);
      REQUIRE(entries.size() == 1);
      const auto &entry = entries.front();
      CHECK(entry.cls == "UnlocatedPatch");
      CHECK(entry.site == "H9a/FindDomainBoundaryExclusions:reference");
      CHECK(entry.what == "patch");
      CHECK(entry.patches == std::vector<int>{1});
      CHECK(entry.models == std::vector<std::string>{"edge"});
      CHECK(entry.features == std::vector<int>{7});
      CHECK_THAT(entry.length, WithinAbs(0.06, 1.0e-12));
      CHECK_THAT(entry.example_origin[0], WithinAbs(1.0 + delta, 1.0e-12));
      CHECK_THAT(entry.example_origin[1], WithinAbs(0.5, 1.0e-12));
      CHECK_THAT(entry.example_origin[2], WithinAbs(0.3, 1.0e-12));
      CHECK_THAT(entry.reason, ContainsSubstring("misplaced or mis-scaled coupon") &&
                                   ContainsSubstring("metal-edge reference") &&
                                   ContainsSubstring("distance 2.000000e-03 mesh units") &&
                                   ContainsSubstring("boundary tolerance 2.000000e-04"));
      // The deferred FailClosed abort names the class; RecordRaw's Enforce is a no-op.
      ledger.Enforce("test");
      GeometryFallbackLedger fail_closed(GeometryFailurePolicy::FAIL_CLOSED, false,
                                         entries);
      CHECK_THROWS_WITH(fail_closed.Enforce("test"),
                        ContainsSubstring("FailClosed") &&
                            ContainsSubstring("UnlocatedPatch 1") &&
                            ContainsSubstring("(test)"));
    }
  }
  // (3) No reference and every tested point outside (a coupon off the mesh).
  {
    std::vector<ResponsePatchData> patches = {inside, Patch({2.0, 0.5, 0.3}, 9)};
    patches[1].conductor_references.clear();
    {
      auto copy = patches;
      CHECK_THROWS_WITH(Find(copy, nullptr),
                        ContainsSubstring("misplaced or mis-scaled coupon") &&
                            ContainsSubstring("patch 1 (0-based") &&
                            ContainsSubstring("every one of its tested points"));
    }
    {
      std::vector<GeometryFallbackData> entries;
      GeometryFallbackLedger ledger(GeometryFailurePolicy::RECORD_RAW, false, entries);
      auto copy = patches;
      const auto exclusions = Find(copy, &ledger);
      REQUIRE(exclusions.patches.size() == 1);
      CHECK(exclusions.patches[0].patch == 1);
      CHECK(exclusions.patches[0].reason == "AllPointsOutside");
      CHECK(exclusions.patches[0].outside_points == exclusions.patches[0].tested_points);
      CHECK(copy[0].weight == 1.0);
      CHECK(copy[1].weight == 0.0);
      REQUIRE(entries.size() == 1);
      CHECK(entries.front().cls == "UnlocatedPatch");
      CHECK(entries.front().site == "H9a/FindDomainBoundaryExclusions:all-points");
      CHECK(entries.front().patches == std::vector<int>{1});
      CHECK(entries.front().features == std::vector<int>{9});
      CHECK_THAT(entries.front().reason,
                 ContainsSubstring("every one of its tested points"));
    }
  }
  // Bitwise neutrality: the well-placed cell alone is applied with no record under either
  // ledger state, as before.
  for (const bool with_ledger : {false, true})
  {
    std::vector<GeometryFallbackData> entries;
    GeometryFallbackLedger ledger(GeometryFailurePolicy::RECORD_RAW, false, entries);
    std::vector<ResponsePatchData> patches = {inside};
    const auto exclusions = Find(patches, with_ledger ? &ledger : nullptr);
    CHECK(exclusions.tested_patches == 1);
    CHECK(exclusions.patches.empty());
    CHECK(patches[0].weight == 1.0);
    CHECK(entries.empty());
  }
#endif
}

// (iii) OI-MX-R2A-COMPARATOR-FORM (readiness hardening-2, decision 1018 Q6; operator
// MapLibraryInterfaces "has no slot 1 SA interface mapping"): a feature whose portions
// carry more distinct interface target maps than its model has slots. The diagonal island
// (one conductor, attributes 9 / 10 split by the island's diagonal) with one SA interface
// per attribute (Index 4 on 9, Index 5 on 10): the two diagonal corners (0.25, 0.75) and
// (0.75, 0.25) claim R = 0.2 along each arm, one arm on attribute 9 and the other crossing
// into attribute 10, so their target maps are {SA -> 4} and {SA -> 5} (slots 0 / 1) while
// the one-slot "convex-corner-90" coupon maps slot 0 alone; the other two corners and every
// straight cell carry one map. On the base the configuration aborts forming the runtime
// model ("has no slot 1 SA interface mapping"); here the two corners are NOT placed (never
// mapped through slot 0) and recorded InterfaceMappingConflict (what "feature", the corner
// id / model / both targets / its claimed length), their claims uncovered (raw), the other
// corners and the straight cells placed as on the same mesh with ONE interface over both
// attributes; the preflight predicts the solve (U-INV); the explicit FailClosed leaf aborts
// at the end of the configuration naming the class.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "Geometry fallback R2A InterfaceMappingConflict: a corner whose arms carry "
    "different interface targets is recorded raw instead of aborting "
    "(decision 1018 Q6)",
    "[surfaceresponseoperator][geometryfallback][invariant][3d][Serial]"
    "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  // The convex library keyed on SA only like the island's features (the domainboundary
  // fixture's form): the isolated-edge and the one-slot convex-corner-90 coupons.
  const fs::path sa_library = temp.temp_dir / "fabrication-process-convex-3d-sa.json";
  if (Mpi::Root(comm))
  {
    json library = ReadJson(convex_library_3d_path);
    library["Name"] = "unit-test-process-convex-3d-sa";
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
    }
    std::ofstream output(sa_library);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(comm);
  // The island config on the diagonal mesh: both attributes metal (one terminal), the SA
  // interface split by attribute (two targets) or over both attributes (one target).
  auto Config = [&](bool split_interfaces)
  {
    json config = ConvexIslandConfig();
    config["Boundaries"]["Terminal"] = {{{"Index", 1}, {"Attributes", {9, 10}}}};
    json &dielectric = config["Boundaries"]["Postprocessing"]["Dielectric"];
    auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
    if (split_interfaces)
    {
      json second = dielectric[0];
      second["Index"] = 5;
      second["Attributes"] = {10};
      dielectric.push_back(second);
      correction["TargetInterfaces"] = {4, 5};
    }
    else
    {
      dielectric[0]["Attributes"] = {9, 10};
    }
    correction["Library"] = sa_library.string();
    correction["UnmatchedPolicy"] = "Warn";
    correction.erase("PatchConstruction");
    correction["CorrectionMode"] = "PostprocessOnly";
    return config;
  };
  auto MakeIoData = [&](const json &config)
  {
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    iodata.boundaries.cracked_attributes.insert(10);
    return iodata;
  };
  auto Preflight = [&](json config, const fs::path &dir)
  {
    config["Problem"]["Output"] = dir.string();
    config["Solver"]["Electrostatic"]["ResponseCorrection"]["GeometryFailurePolicy"] =
        "RecordRaw";
    IoData iodata = MakeIoData(config);
    MakeOutputFolder(iodata, comm);
    Mesh mesh(MakeDiagonalIslandMesh());
    const fs::path manifest_path = dir / "surface-response-requirements.json";
    WriteSurfaceResponseRequirements(iodata, mesh, manifest_path.string());
    Mpi::Barrier(comm);
    return ReadJson(manifest_path);
  };
  // The reference: one interface over both attributes places all four corners.
  const json reference = Preflight(Config(false), temp.temp_dir / "r2a-reference");
  CHECK(!reference.at("Identification").at("Diagnostics").contains("GeometryFallbacks"));
  CHECK(reference.at("Summary").at("Counts").at("Missing").get<int>() == 0);
  auto CornerIds = [](const json &manifest)
  {
    std::vector<int> ids;
    for (const auto &feature : manifest.at("Identification").at("Features"))
    {
      if (feature.at("Type") == "ConvexCorner")
      {
        ids.push_back(feature.at("Id").get<int>());
      }
    }
    return ids;
  };
  REQUIRE(CornerIds(reference).size() == 4);
  CHECK(reference.at("Summary").at("Uncovered").at("Features").get<int>() == 0);

  // The split interfaces: two InterfaceMappingConflict entries on the preflight.
  const json manifest = Preflight(Config(true), temp.temp_dir / "r2a-split");
  const auto &record =
      manifest.at("Identification").at("Diagnostics").at("GeometryFallbacks");
  REQUIRE(record.at("Count").get<int>() == 2);
  CHECK(record.at("Classes").at("InterfaceMappingConflict").get<int>() == 2);
  CHECK(record.at("Policy") == "RecordRaw");
  std::set<int> rejected;
  for (const auto &entry : record.at("Entries"))
  {
    CHECK(entry.at("Class") == "InterfaceMappingConflict");
    CHECK(entry.at("Site") == "R2A/BuildFeaturePatches:interface-slots");
    CHECK(entry.at("What") == "feature");
    CHECK(entry.at("Topology") == "ConvexCorner");
    CHECK(entry.at("Models") == json::array({"convex-corner-90"}));
    CHECK(entry.at("Interfaces") == json::array({4, 5}));
    CHECK_THAT(entry.at("Reason").get<std::string>(),
               ContainsSubstring("cannot map the 2 distinct interface target map(s)") &&
                   ContainsSubstring("slot 0 {SA -> 4}; slot 1 {SA -> 5}"));
    // The corner's claims: R along each arm.
    CHECK_THAT(entry.at("Length").get<double>(), WithinAbs(0.4, 1.0e-9));
    REQUIRE(entry.at("Features").size() == 1);
    rejected.insert(entry.at("Features").front().get<int>());
  }
  CHECK(rejected.size() == 2);
  // The two diagonal corners (x + z = 1), by their recorded frames.
  const auto corner_ids = CornerIds(manifest);
  REQUIRE(corner_ids.size() == 4);
  for (const auto &feature : manifest.at("Identification").at("Features"))
  {
    if (feature.at("Type") != "ConvexCorner")
    {
      continue;
    }
    const auto origin = feature.at("Frame").at("Origin").get<std::array<double, 3>>();
    const bool diagonal = std::abs(origin[0] + origin[2] - 1.0) < 1.0e-9;
    CHECK(diagonal == (rejected.count(feature.at("Id").get<int>()) > 0));
  }
  // The rejected corners are read as Missing: their claims are the uncovered portions
  // (raw), the Summary mirror names the class and the raw share by length.
  CHECK(manifest.at("Summary").at("Uncovered").at("Features").get<int>() == 2);
  CHECK_THAT(manifest.at("Summary").at("Uncovered").at("Length").get<double>(),
             WithinAbs(0.8, 1.0e-9));
  CHECK(manifest.at("Summary").at("GeometryFallbacks").at("Count").get<int>() == 2);
  const auto &share = manifest.at("Summary").at("GeometryFallbacks").at("RawShareByLength");
  REQUIRE(share.contains("SA"));
  CHECK_THAT(share.at("SA").at("FallbackLength").get<double>(), WithinAbs(0.8, 1.0e-9));
  CHECK_THAT(share.at("SA").at("TypeLength").get<double>(), WithinAbs(2.0, 1.0e-9));
  // The operator (the solve default RecordRaw): the same two entries (U-INV), not raw as a
  // whole; the explicit FailClosed leaf aborts naming the class.
  {
    IoData iodata = MakeIoData(Config(true));
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeDiagonalIslandMesh()));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    CHECK(!response.IsRawFallback());
    REQUIRE(response.GetGeometryFallbacks().size() == 2);
    for (const auto &entry : response.GetGeometryFallbacks())
    {
      CHECK(entry.cls == "InterfaceMappingConflict");
      CHECK(entry.interfaces == std::vector<int>{4, 5});
    }
    const json solve_record = response.GetGeometryFallbackRecord();
    REQUIRE(!solve_record.is_null());
    CHECK(solve_record.at("Classes").at("InterfaceMappingConflict").get<int>() == 2);
    CHECK(CheckGeometryFallbackInvariant(manifest, solve_record) == 2);
    CHECK_THAT(GeometryFallbackBanner(solve_record, "palace.json"),
               ContainsSubstring("LOW CONFIDENCE") &&
                   ContainsSubstring("InterfaceMappingConflict 2"));
  }
  {
    json fail_closed = Config(true);
    fail_closed["Solver"]["Electrostatic"]["ResponseCorrection"]["GeometryFailurePolicy"] =
        "FailClosed";
    IoData iodata = MakeIoData(fail_closed);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeDiagonalIslandMesh()));
    LaplaceOperator laplace(iodata, meshes);
    CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, laplace),
                      ContainsSubstring("GeometryFailurePolicy FailClosed") &&
                          ContainsSubstring("InterfaceMappingConflict 2"));
  }
  // The reference configuration on the operator: no record (bitwise neutrality of the
  // predicate on a one-map feature set).
  {
    IoData iodata = MakeIoData(Config(false));
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeDiagonalIslandMesh()));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    CHECK(response.GetGeometryFallbacks().empty());
    CHECK(response.GetGeometryFallbackRecord().is_null());
  }
#endif
}

// H4a (operator ApplyPlacementOwnership "matching volumes overlap without complete
// spatial-cluster ownership"; the dry run's FindSpatialClusterPriorityExclusions carries
// the same rule): the S1p pair of the margin-overlap case (decision 244) with B's claim
// reaching into A's claims hull. Without the ledger the same-class pair is skipped (the
// record of decision 244, the operator aborts); with it the SUBORDINATE (the smaller
// claimed length: B, 13.675 vs 15.475) is owned with the OverlapExcluded record naming both
// patches / models; a margins-only pair records nothing. The operator's branch shares
// OverlapExcludedSubordinate.
TEST_CASE(
    "Geometry fallback H4a OverlapExcluded: an owned-less same-class overlap drops the "
    "subordinate support with the record (decisions 1000 (b) / 1001)",
    "[surfaceresponseoperator][geometryfallback][Serial]")
{
  using Patch = config::ElectrostaticSolverData::ResponseCorrectionPatchData;
  using Claim = Patch::Provenance::Claim;
  using Model = config::ElectrostaticSolverData::ResponseCorrectionModelData;
  const double R = 1.9;
  const double tolerance = kSignatureParameterToleranceOverRadius * R;
  const double y = -108.325;
  SpatialSupportBounds a, b;
  a.patch = 0;
  a.SetMeshAlignedBox({568.55, -128.2, 2.85}, {592.75, -102.625, 6.8});
  a.claims = {Claim{494, {583.25, y, 4.8}, {583.55, y, 4.8}},
              Claim{514, {583.55, y, 4.8}, {587.05, y, 4.8}},
              Claim{520, {583.25, y, 4.8}, {583.25, -120.0, 4.8}}};
  b.patch = 1;
  b.SetMeshAlignedBox({583.05, -128.2, 2.85}, {607.25, -102.625, 6.8});
  // B's second claim turns into A's hull (claim_in_hull: the true overlap).
  b.claims = {Claim{646, {586.0, y, 4.8}, {593.0, y, 4.8}},
              Claim{650, {586.0, y, 4.8}, {586.0, -115.0, 4.8}}};
  config::ElectrostaticSolverData::ResponseCorrectionData config;
  config.matching_radius = R;
  for (int k = 0; k < 2; k++)
  {
    Model model;
    model.idx = k + 1;
    model.name = k == 0 ? "cluster-a" : "cluster-b";
    model.topology = "spatial edge cluster";
    model.spatial_basis = true;
    config.models.push_back(model);
    Patch patch;
    patch.model = k + 1;
    patch.weight = 1.0;
    patch.provenance.feature = 10 + k;
    patch.provenance.claims = k == 0 ? a.claims : b.claims;
    patch.origin = k == 0 ? std::array<double, 3>{580.0, y, 4.8}
                          : std::array<double, 3>{595.0, y, 4.8};
    config.patches.push_back(patch);
  }
  const auto margin_overlaps = FindSpatialSupportMarginOverlaps({a, b}, 3, tolerance);
  REQUIRE(margin_overlaps.size() == 1);
  REQUIRE(margin_overlaps.front().claim_in_hull);
  // The rule of record without a ledger: the same-class pair is skipped (0 SCP drops).
  CHECK(FindSpatialClusterPriorityExclusions({a, b}, config, {}, 3).empty());
  CHECK(FindSpatialClusterPriorityExclusions({a, b}, config, {}, 3, &margin_overlaps)
            .empty());
  // The subordinate rule: B's claims 7.0 + 6.675 = 13.675 < A's 0.3 + 3.5 + 11.675
  // = 15.475.
  CHECK(OverlapExcludedSubordinate(config, 0, 1) == 1);
  // With the ledger: B is owned, one OverlapExcluded record.
  {
    std::vector<GeometryFallbackData> entries;
    GeometryFallbackLedger ledger(GeometryFailurePolicy::RECORD_RAW, false, entries);
    std::vector<std::pair<std::size_t, std::size_t>> overlap_excluded;
    const auto scp = FindSpatialClusterPriorityExclusions(
        {a, b}, config, {}, 3, &margin_overlaps, &ledger, &overlap_excluded);
    CHECK(scp.empty());
    REQUIRE(overlap_excluded.size() == 1);
    CHECK(overlap_excluded.front() == std::make_pair(std::size_t{1}, std::size_t{0}));
    REQUIRE(entries.size() == 1);
    const auto &entry = entries.front();
    CHECK(entry.cls == "OverlapExcluded");
    CHECK(entry.what == "patch");
    CHECK(entry.patches == std::vector<int>{1, 0});
    CHECK(entry.models == std::vector<std::string>{"cluster-b", "cluster-a"});
    CHECK(entry.features == std::vector<int>{11});
    CHECK_THAT(entry.length, WithinAbs(13.675, 1.0e-9));
    CHECK_THAT(entry.reason, ContainsSubstring("patches 1 and 2") &&
                                 ContainsSubstring("inside the other's claims"));
    // The ties rule: equal claimed lengths drop the larger index.
    auto tie = config;
    tie.patches[1].provenance.claims = tie.patches[0].provenance.claims;
    CHECK(OverlapExcludedSubordinate(tie, 0, 1) == 1);
    CHECK(OverlapExcludedSubordinate(tie, 1, 0) == 1);
  }
  // A margins-only pair (decision 244: no claim inside the other's hull) records nothing
  // even under the ledger.
  {
    SpatialSupportBounds margins = b;
    margins.claims = {Claim{646, {589.95, y, 4.8}, {593.0, y, 4.8}},
                      Claim{650, {593.0, y, 4.8}, {593.0, -120.0, 4.8}}};
    const auto margins_only = FindSpatialSupportMarginOverlaps({a, margins}, 3, tolerance);
    REQUIRE(margins_only.size() == 1);
    REQUIRE(!margins_only.front().claim_in_hull);
    std::vector<GeometryFallbackData> entries;
    GeometryFallbackLedger ledger(GeometryFailurePolicy::RECORD_RAW, false, entries);
    std::vector<std::pair<std::size_t, std::size_t>> overlap_excluded;
    CHECK(FindSpatialClusterPriorityExclusions({a, margins}, config, {}, 3, &margins_only,
                                               &ledger, &overlap_excluded)
              .empty());
    CHECK(overlap_excluded.empty());
    CHECK(entries.empty());
  }
}

// H6a / H6b (operator ApplyContinuationOwnership: the owning boxes leave a gap inside a
// cell / a box lies strictly inside a cell). By construction every continued support
// contains the stretch's first end, so every inside interval starts at the cell's begin
// and neither condition arises from a consistent placement (the code's "verified, not
// assumed"); the branch is reached here with a support whose claims lie OUTSIDE its own
// box (an inconsistent input the guard exists for). Without the ledger the abort of record
// stands; with it the cell is left untouched and ContinuationUnresolved is recorded.
TEST_CASE("Geometry fallback H6 ContinuationUnresolved: a gap inside a continued cell is "
          "recorded and the cell left untouched (decisions 1000 (b) / 1001)",
          "[surfaceresponseoperator][geometryfallback][Serial]")
{
  using Patch = config::ElectrostaticSolverData::ResponseCorrectionPatchData;
  using Claim = Patch::Provenance::Claim;
  auto Cell = [](int feature, int stretch, double begin, double end)
  {
    Patch patch;
    patch.origin = {0.5 * (begin + end), 0.0, 0.0};
    patch.axis_u = {0.0, 1.0, 0.0};
    patch.axis_v = {0.0, 0.0, 1.0};
    patch.axis_w = {1.0, 0.0, 0.0};
    patch.longitudinal_cell = {begin - patch.origin[0], end - patch.origin[0]};
    patch.provenance.feature = feature;
    patch.provenance.stretch = stretch;
    patch.provenance.segment = 7;
    patch.weight = 1.0;
    return patch;
  };
  const double R = 1.5;
  const double tolerance = kSignatureParameterToleranceOverRadius * R;
  // Support A: box x in [-3, 1.5], its claim on y = 0 ends at x = 0 (continued along +x).
  SpatialSupportBounds a;
  a.patch = 0;
  a.SetMeshAlignedBox({-3.0, -3.0, -2.0}, {1.5, 3.0, 2.0});
  a.claims = {Claim{3, {-2.5, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  // Support B: the same claim end at x = 0 (so the stretch continues it) but a box
  // x in [2.0, 4.5] that does not contain it: on the cell [0, 3] its inside interval
  // [2.0, 3.0] leaves the gap (1.5, 2.0) after A's [0, 1.5].
  SpatialSupportBounds b;
  b.patch = 1;
  b.SetMeshAlignedBox({2.0, -3.0, -2.0}, {4.5, 3.0, 2.0});
  b.claims = {Claim{3, {-2.5, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  std::vector<Patch> patches = {Cell(4, 0, 0.0, 3.0), Cell(4, 0, 3.0, 4.0)};
  {
    auto copy = patches;
    CHECK_THROWS_WITH(ApplyContinuationOwnership(copy, {a, b}, 3, tolerance, R),
                      ContainsSubstring("leave a gap inside the cell"));
  }
  {
    std::vector<GeometryFallbackData> entries;
    GeometryFallbackLedger ledger(GeometryFailurePolicy::RECORD_RAW, false, entries);
    auto copy = patches;
    const auto before = copy;
    const auto ownership =
        ApplyContinuationOwnership(copy, {a, b}, 3, tolerance, R, &ledger);
    REQUIRE(entries.size() == 1);
    CHECK(entries.front().cls == "ContinuationUnresolved");
    CHECK(entries.front().site == "H6a/ApplyContinuationOwnership:gap-inside-cell");
    CHECK(entries.front().what == "cell");
    CHECK(entries.front().features == std::vector<int>{4});
    CHECK(entries.front().patches == std::vector<int>{0, 0, 1});
    CHECK_THAT(entries.front().reason, ContainsSubstring("patch 1 leave a gap"));
    // The cell is untouched (weight, cell, quadrature weight), not owned.
    CHECK(copy[0].weight == before[0].weight);
    CHECK(copy[0].longitudinal_cell == before[0].longitudinal_cell);
    // The next cell [3, 4] lies wholly inside B's box: owned as usual (the fallback is
    // per cell).
    REQUIRE(ownership.cells.size() == 1);
    CHECK(ownership.cells.front().patch == 1);
  }
}

// U-INV (DESIGN section 4; decisions 1001 / 1008): the preflight predicts the solve's
// geometry fallbacks. For the same (mesh, config, library) the manifest written by
// WriteSurfaceResponseRequirements under the explicit RecordRaw leaf (P) and the operator's
// record (S) agree on (Class, Features, Patches, Models, Length) for every geometry class.
// H1 is the A10.1 form (P has no entry; Summary.Counts matched == 0 predicts the raw run).
// H8 fires end to end on the island with the convex library keyed on SA only (the
// domainboundary fixture's form) whose every model's BasisPoints are widened beyond the
// unit domain (the isolated edge at +-2.0, the corner's 12 points scaled 30x): every placed
// coupon leaves the mesh on the dry run and on the operator
// (the same FindDomainBoundaryExclusions on the same patch table), so both record ONE
// AllPatchesOutside entry over the same patches. H4a / H5 / H6 / H13-I1 and H7a: no
// end-to-end fixture here (H4a / H6 are exercised on the shared dry-run functions, A10.9 /
// A10.10; H5 is GUARDED by construction, A10.8; H13-I1 needs a metal-free perimeter; H7a's
// invariant is the role-swap alias case of test-surfaceresponseoperator.cpp).
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "Geometry fallback invariant U-INV: the preflight predicts the solve's "
    "geometry fallbacks (H1 by Counts, H8 by the entry; decisions 1001 / 1008)",
    "[surfaceresponseoperator][geometryfallback][invariant][domainboundary]"
    "[3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  // The preflight (the RecordRaw leaf: the complete prediction with no abort) and the
  // operator (the solve default) on the same configuration and mesh.
  auto Preflight = [&](json config, const fs::path &dir)
  {
    config["Problem"]["Output"] = dir.string();
    config["Solver"]["Electrostatic"]["ResponseCorrection"]["GeometryFailurePolicy"] =
        "RecordRaw";
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    MakeOutputFolder(iodata, comm);
    Mesh mesh(MakeIslandMesh());
    const fs::path manifest_path = dir / "surface-response-requirements.json";
    WriteSurfaceResponseRequirements(iodata, mesh, manifest_path.string());
    Mpi::Barrier(comm);
    return ReadJson(manifest_path);
  };
  auto Solve = [&](const json &config)
  {
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    return std::make_pair(response.IsRawFallback(), response.GetGeometryFallbackRecord());
  };

  SECTION("H1 NoMatchedFeature: no entry on the preflight, Counts predict the raw run")
  {
    json config = ConvexIslandConfig();
    auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
    correction["Library"] = spatial_cluster_library_3d_path.string();
    correction["UnmatchedPolicy"] = "Warn";
    correction.erase("PatchConstruction");
    correction["CorrectionMode"] = "PostprocessOnly";
    const json manifest = Preflight(config, temp.temp_dir / "uinv-h1");
    const auto [raw, record] = Solve(config);
    CHECK(raw);
    REQUIRE(!record.is_null());
    CHECK(record.at("Classes").at("NoMatchedFeature").get<int>() == 1);
    // MAJOR-A of the d1000 confirmation (decision 1011): the whole-run entry counts 100 %
    // of every Type by length on the solve record (H1 names its interfaces with the max
    // per-interface length; the share no longer depends on that).
    REQUIRE(record.at("RawShareByLength").contains("SA"));
    CHECK_THAT(record.at("RawShareByLength").at("SA").at("Share").get<double>(),
               WithinAbs(1.0, 1.0e-12));
    CHECK(record.at("RawShareByEnergy").at("SA").at("Share").get<double>() == 1.0);
    CHECK(CheckGeometryFallbackInvariant(manifest, record) == 0);
  }

  SECTION("H8 AllPatchesOutside: the same entry over the same patches on both")
  {
    const fs::path wide_points = temp.temp_dir / "uinv-wide-basis-points.csv";
    const fs::path wide_corner_points = temp.temp_dir / "uinv-wide-corner-points.csv";
    const fs::path wide_library = temp.temp_dir / "fabrication-process-convex-3d-wide.json";
    if (Mpi::Root(comm))
    {
      {
        std::ofstream output(wide_points);
        output << "x,y,z\n-2.0,-2.0,0.0\n2.0,-2.0,0.0\n2.0,2.0,0.0\n-2.0,2.0,0.0\n";
      }
      {
        // The fixture's corner points (0.02 / 0.08 in plane, z = -0.02 / 0 / 0.02) x 30:
        // the normal offset +-0.6 alone puts every point of the two outer rings outside
        // the unit box; the middle ring's 0.6 / 2.4 in-plane offsets leave it too.
        std::ofstream output(wide_corner_points);
        output << "x,y,z\n";
        for (const double z : {-0.6, 0.0, 0.6})
        {
          output << "0.6,0.6," << z << "\n2.4,0.6," << z << "\n2.4,2.4," << z
                 << "\n0.6,2.4," << z << "\n";
        }
      }
      json library = ReadJson(convex_library_3d_path);
      library["Name"] = "unit-test-process-convex-3d-wide";
      for (auto &model : library["Models"])
      {
        // Keyed on SA only like the island's features (the domainboundary fixture's form).
        model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
        model["BasisPoints"] = model["Topology"] == "IsolatedEdge"
                                   ? wide_points.string()
                                   : wide_corner_points.string();
      }
      std::ofstream output(wide_library);
      output << library.dump(2) << "\n";
    }
    Mpi::Barrier(comm);
    json config = ConvexIslandConfig();
    auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
    correction["Library"] = wide_library.string();
    correction["UnmatchedPolicy"] = "Warn";
    correction.erase("PatchConstruction");
    correction["TraceCoupling"] = "SurfaceMortar";
    const json manifest = Preflight(config, temp.temp_dir / "uinv-h8");
    const auto &preflight_record =
        manifest.at("Identification").at("Diagnostics").at("GeometryFallbacks");
    REQUIRE(preflight_record.at("Count").get<int>() == 1);
    CHECK(preflight_record.at("Classes").at("AllPatchesOutside").get<int>() == 1);
    CHECK(preflight_record.at("Policy") == "RecordRaw");
    const auto &entry = preflight_record.at("Entries").front();
    CHECK(entry.at("What") == "run");
    CHECK(entry.at("Patches").size() == manifest.at("Identification")
                                            .at("Diagnostics")
                                            .at("DomainBoundaryExclusions")
                                            .at("TestedPatches")
                                            .get<std::size_t>());
    CHECK(entry.at("Patches").size() > 1);
    // MAJOR-A of the d1000 confirmation (decision 1011): the whole-run fallback reads 100 %
    // of every selected Type BY LENGTH on the preflight record (the only record of an H8
    // run under the discovery's RecordRaw leaf; the census / stage D read it), its Summary
    // mirror and the solve record alike; the entry's Length is the total selected length
    // (the island selects SA only), so the U-INV Length comparison below is non-trivial.
    const auto &preflight_share = preflight_record.at("RawShareByLength");
    REQUIRE(preflight_share.size() == 1);
    REQUIRE(preflight_share.contains("SA"));
    const double type_length = preflight_share.at("SA").at("TypeLength").get<double>();
    CHECK(type_length > 0.0);
    CHECK_THAT(preflight_share.at("SA").at("FallbackLength").get<double>(),
               WithinAbs(type_length, 1.0e-12 * type_length));
    CHECK_THAT(preflight_share.at("SA").at("Share").get<double>(), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(entry.at("Length").get<double>(),
               WithinAbs(type_length, 1.0e-12 * type_length));
    CHECK_THAT(manifest.at("Summary")
                   .at("GeometryFallbacks")
                   .at("RawShareByLength")
                   .at("SA")
                   .at("Share")
                   .get<double>(),
               WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(GeometryFallbackBanner(preflight_record, "manifest"),
               ContainsSubstring("SA 100.00 % of the length") &&
                   ContainsSubstring("AllPatchesOutside 1"));
    const auto [raw, record] = Solve(config);
    CHECK(raw);
    REQUIRE(!record.is_null());
    CHECK(record.at("Classes").at("AllPatchesOutside").get<int>() == 1);
    CHECK_THAT(record.at("RawShareByLength").at("SA").at("Share").get<double>(),
               WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(record.at("Entries").front().at("Length").get<double>(),
               WithinAbs(type_length, 1.0e-12 * type_length));
    CHECK(record.at("RawShareByEnergy").at("SA").at("Share").get<double>() == 1.0);
    CHECK(CheckGeometryFallbackInvariant(manifest, record) == 1);
  }
#endif
}

// MAJOR-A of the d1000 confirmation (decision 1011) on the record function itself: a
// whole-run entry (What "run") counts the TOTAL selected length of EVERY type whatever the
// site recorded (H8 records length 0 and no interface; H1 records its interfaces with the
// MAX per-interface length), its Length is that total, and a feature-level entry keeps the
// per-type rule (its length on the types of the interfaces it names).
TEST_CASE(
    "Geometry fallback record: a whole-run entry counts 100 % of every Type by length "
    "(decision 1011)",
    "[surfaceresponseoperator][geometryfallback][Serial]")
{
  // Two SA interfaces (lengths 3 and 5) and one MS interface (length 4), the lengths in
  // mesh units; the record scaled by 2.
  const std::map<int, double> lengths = {{1, 3.0}, {2, 5.0}, {3, 4.0}};
  const std::map<int, InterfaceDielectric> types = {{1, InterfaceDielectric::SA},
                                                    {2, InterfaceDielectric::SA},
                                                    {3, InterfaceDielectric::MS}};
  const double scale = 2.0;
  auto Share = [](const json &record, const char *type)
  { return record.at("RawShareByLength").at(type); };

  SECTION("H8 form: What run, length 0, no interface -> 100 % of SA and MS")
  {
    GeometryFallbackData entry;
    entry.cls = "AllPatchesOutside";
    entry.site = "H8/FindDomainBoundaryExclusions:all-excluded";
    entry.what = "run";
    entry.patches = {0, 1};
    const json record = DescribeGeometryFallbacks(
        {entry}, GeometryFailurePolicy::RECORD_RAW, true, scale, lengths, types);
    CHECK(record.at("Run") == "Preflight");
    CHECK_THAT(Share(record, "SA").at("TypeLength").get<double>(),
               WithinAbs(8.0 * scale, 1.0e-12));
    CHECK_THAT(Share(record, "SA").at("FallbackLength").get<double>(),
               WithinAbs(8.0 * scale, 1.0e-12));
    CHECK_THAT(Share(record, "SA").at("Share").get<double>(), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(Share(record, "MS").at("FallbackLength").get<double>(),
               WithinAbs(4.0 * scale, 1.0e-12));
    CHECK_THAT(Share(record, "MS").at("Share").get<double>(), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(record.at("Entries").front().at("Length").get<double>(),
               WithinAbs(12.0 * scale, 1.0e-12));
    CHECK_THAT(GeometryFallbackBanner(record, "manifest"),
               ContainsSubstring("MS 100.00 % of the length") &&
                   ContainsSubstring("SA 100.00 % of the length"));
  }

  SECTION("H1 form: What run naming the interfaces with the max length -> 100 % of both")
  {
    GeometryFallbackData entry;
    entry.cls = "NoMatchedFeature";
    entry.site = "H1/BuildAutomaticResponseData3D:features";
    entry.what = "run";
    entry.interfaces = {1, 2, 3};
    entry.length = 5.0;
    const json record = DescribeGeometryFallbacks(
        {entry}, GeometryFailurePolicy::RECORD_RAW, false, scale, lengths, types);
    CHECK(record.at("Run") == "Solve");
    CHECK_THAT(Share(record, "SA").at("Share").get<double>(), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(Share(record, "MS").at("Share").get<double>(), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(record.at("Entries").front().at("Length").get<double>(),
               WithinAbs(12.0 * scale, 1.0e-12));
  }

  SECTION("a feature-level entry keeps the per-type rule beside a whole-run one")
  {
    GeometryFallbackData feature;
    feature.cls = "KeyMismatch";
    feature.site = "H7a";
    feature.what = "feature";
    feature.features = {7};
    feature.interfaces = {1};
    feature.length = 1.5;
    const json alone = DescribeGeometryFallbacks(
        {feature}, GeometryFailurePolicy::RECORD_RAW, false, scale, lengths, types);
    CHECK_THAT(Share(alone, "SA").at("FallbackLength").get<double>(),
               WithinAbs(1.5 * scale, 1.0e-12));
    CHECK_THAT(Share(alone, "SA").at("Share").get<double>(), WithinAbs(1.5 / 8.0, 1.0e-12));
    CHECK_THAT(Share(alone, "MS").at("Share").get<double>(), WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(alone.at("Entries").front().at("Length").get<double>(),
               WithinAbs(1.5 * scale, 1.0e-12));
    GeometryFallbackData run;
    run.cls = "NoPerimeter";
    run.site = "H13-I1";
    run.what = "run";
    const json both = DescribeGeometryFallbacks(
        {feature, run}, GeometryFailurePolicy::FAIL_CLOSED, true, scale, lengths, types);
    CHECK(both.at("Count").get<int>() == 2);
    CHECK_THAT(Share(both, "SA").at("Share").get<double>(), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(Share(both, "MS").at("Share").get<double>(), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(Share(both, "SA").at("FallbackLength").get<double>(),
               WithinAbs(8.0 * scale, 1.0e-12));
  }
}
