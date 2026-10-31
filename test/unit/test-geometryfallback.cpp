// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
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
#include "utils/metaledge.hpp"
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

// U-H13b (hardening-2 DESIGN row H13-I2, landed in the identification fix round's lane B):
// a PINCHED metal sheet on the plane y = 0.5 of a 3 x 1 x 3 box (hexahedra, h = 0.25): the
// sheet [0, 3]^2 minus the hole [1, 2]^2 and, when `pinched`, minus the corner cell [2,
// 3]^2, so that the hole's loop and the corner notch meet at the vertex (2, 0.5, 2) with
// FOUR perimeter segments: a chain vertex with more than two chain segments. The sheet
// reaches the box walls (its rim is a cut); the hole and the notch are physical perimeter.
// The sheet's defects (U-H10 / U-H11c alike): DUPLICATE_ATTRIBUTE adds the plane face at
// the hole's rim centred (0.875, 0.5, 1.375) a SECOND time under attribute 10 (coincident
// boundary elements of different metal attributes: H10), DEGENERATE_FACE moves two
// vertices of the plane face [0, 0.25] x [0.25, 0.5] onto its other two so that it keeps
// two distinct canonical vertices (H11c; the collapsed hexahedra are never integrated
// here).
enum class SheetDefect
{
  NONE,
  PINCH,
  DUPLICATE_ATTRIBUTE,
  DEGENERATE_FACE
};

mfem::Mesh MakeSerialSheetMesh(SheetDefect defect)
{
  mfem::Mesh serial =
      mfem::Mesh::MakeCartesian3D(12, 4, 12, mfem::Element::HEXAHEDRON, 3.0, 1.0, 3.0);
  if (defect == SheetDefect::DEGENERATE_FACE)
  {
    // The two x = 0.25 vertices of the plane face [0, 0.25] x [0.25, 0.5] onto its x = 0
    // vertices: the face keeps two distinct canonical points.
    for (int v = 0; v < serial.GetNV(); v++)
    {
      double *point = serial.GetVertex(v);
      if (std::abs(point[0] - 0.25) < 1.0e-12 && std::abs(point[1] - 0.5) < 1.0e-12 &&
          (std::abs(point[2] - 0.25) < 1.0e-12 || std::abs(point[2] - 0.5) < 1.0e-12))
      {
        point[0] = 0.0;
      }
    }
  }
  for (int face = 0; face < serial.GetNumFaces(); face++)
  {
    int element1, element2;
    serial.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    serial.GetFaceVertices(face, vertices);
    bool on_plane = true;
    double x = 0.0, z = 0.0;
    for (const int vertex : vertices)
    {
      const double *point = serial.GetVertex(vertex);
      on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
      x += 0.25 * point[0];
      z += 0.25 * point[2];
    }
    if (!on_plane)
    {
      continue;
    }
    const bool in_hole = x > 1.0 && x < 2.0 && z > 1.0 && z < 2.0;
    const bool in_corner = defect == SheetDefect::PINCH && x > 2.0 && z > 2.0;
    if (in_hole || in_corner)
    {
      continue;  // bare dielectric
    }
    serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
    serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
    if (defect == SheetDefect::DUPLICATE_ATTRIBUTE && std::abs(x - 0.875) < 1.0e-9 &&
        std::abs(z - 1.375) < 1.0e-9)
    {
      serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
      serial.SetBdrAttribute(serial.GetNBE() - 1, 10);
    }
  }
  serial.FinalizeTopology();
  serial.Finalize();
  return serial;
}

mfem::Mesh MakeSerialPinchedSheetMesh(bool pinched)
{
  return MakeSerialSheetMesh(pinched ? SheetDefect::PINCH : SheetDefect::NONE);
}

// The sheet tests' library: the concave library's models (concave corners, isolated edges)
// on the config's interface (Type SA, Coupon 1), written by the root.
fs::path WriteSheetLibrary(const fs::path &directory, const fs::path &template_library,
                           const std::string &suffix)
{
  const auto path = directory / ("fabrication-process-sheet-" + suffix + ".json");
  if (Mpi::Root(Mpi::World()))
  {
    json library = ReadJson(template_library);
    library["Name"] = "unit-test-process-sheet-" + suffix;
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
    }
    std::ofstream output(path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  return path;
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

// (iv) B2 (readiness hardening-2, decision 1018 Q1; motivated by decision 1007 /
// findings-from-exec-B F-B1, the S1b-rc concave pair, since MEASURED a true overlap of the
// coupon cubes recorded OverlapExcluded on both tests, decision 1030): the same-class
// overlap test decides on the ORIENTED boxes (the coupon domains) once the mesh-aligned
// boxes overlap. A 2 x 0.4 x 1 coupon box rotated 45 deg in the plane (A) and a
// mesh-aligned 0.8 x 0.8 x 1 box (B) whose mesh-aligned boxes overlap: B at (1, -1)
// lies 1.414 along A's short axis beyond A's 0.2 and B's 0.566 reach (disjoint: no overlap,
// no record); B at (0.5, -0.5) reaches A (overlap: OverlapExcluded as H4a). Mesh-aligned
// frames reproduce InteriorOverlaps (touching faces are no interior overlap).
TEST_CASE("Geometry fallback B2 oriented overlap: a same-class pair whose oriented boxes "
          "are disjoint is no overlap (decision 1018 Q1)",
          "[surfaceresponseoperator][geometryfallback][Serial]")
{
  using Patch = config::ElectrostaticSolverData::ResponseCorrectionPatchData;
  using Claim = Patch::Provenance::Claim;
  using Model = config::ElectrostaticSolverData::ResponseCorrectionModelData;
  const double c = std::sqrt(0.5);
  OrientedSupportBox a;
  a.axes = {{{c, c, 0.0}, {-c, c, 0.0}, {0.0, 0.0, 1.0}}};
  a.min = {-1.0, -0.2, -0.5};
  a.max = {1.0, 0.2, 0.5};
  auto MeshAligned =
      [](const std::array<double, 3> &centre, const std::array<double, 3> &half)
  {
    OrientedSupportBox box;
    for (int k = 0; k < 3; k++)
    {
      box.min[k] = centre[k] - half[k];
      box.max[k] = centre[k] + half[k];
    }
    return box;
  };
  const OrientedSupportBox far = MeshAligned({1.0, -1.0, 0.0}, {0.4, 0.4, 0.5});
  const OrientedSupportBox near = MeshAligned({0.5, -0.5, 0.0}, {0.4, 0.4, 0.5});
  // The mesh-aligned boxes: A's is +-(1 + 0.2) / sqrt(2) in x and y; both pairs overlap.
  struct Aabb
  {
    std::array<double, 3> min{}, max{};
  };
  const double a_reach = 1.2 * c;
  const Aabb a_aabb = {{-a_reach, -a_reach, -0.5}, {a_reach, a_reach, 0.5}};
  const Aabb far_aabb = {far.min, far.max};
  const Aabb near_aabb = {near.min, near.max};
  auto AabbInteriorOverlap = [](const Aabb &p, const Aabb &q)
  {
    for (int d = 0; d < 3; d++)
    {
      if (std::min(p.max[d], q.max[d]) - std::max(p.min[d], q.min[d]) <= 1.0e-12)
      {
        return false;
      }
    }
    return true;
  };
  REQUIRE(AabbInteriorOverlap(a_aabb, far_aabb));
  REQUIRE(AabbInteriorOverlap(a_aabb, near_aabb));
  CHECK(!OrientedBoxesInteriorOverlap(a, far, 3));
  CHECK(!OrientedBoxesInteriorOverlap(far, a, 3));
  CHECK(OrientedBoxesInteriorOverlap(a, near, 3));
  CHECK(OrientedBoxesInteriorOverlap(near, a, 3));
  // Mesh-aligned frames: the separating-axis test is the mesh-aligned test (overlapping,
  // touching along a face, disjoint), in 3D and in 2D.
  const OrientedSupportBox unit = MeshAligned({0.5, 0.5, 0.5}, {0.5, 0.5, 0.5});
  CHECK(
      OrientedBoxesInteriorOverlap(unit, MeshAligned({0.9, 0.5, 0.5}, {0.5, 0.5, 0.5}), 3));
  CHECK(!OrientedBoxesInteriorOverlap(unit, MeshAligned({1.5, 0.5, 0.5}, {0.5, 0.5, 0.5}),
                                      3));
  CHECK(!OrientedBoxesInteriorOverlap(unit, MeshAligned({2.5, 0.5, 0.5}, {0.5, 0.5, 0.5}),
                                      3));
  CHECK(
      OrientedBoxesInteriorOverlap(unit, MeshAligned({0.9, 0.5, 7.0}, {0.5, 0.5, 0.5}), 2));
  CHECK(!OrientedBoxesInteriorOverlap(unit, MeshAligned({0.5, 1.5, 0.5}, {0.5, 0.5, 0.5}),
                                      2));
  // A rotated box inside a mesh-aligned one overlaps it on every axis.
  CHECK(OrientedBoxesInteriorOverlap(a, MeshAligned({0.0, 0.0, 0.0}, {3.0, 3.0, 3.0}), 3));

  // The dry run's same-class branch (FindSpatialClusterPriorityExclusions) on two corner
  // supports with these boxes: the far pair records nothing under the ledger, the near
  // pair the OverlapExcluded record of H4a (B, the shorter claims, the subordinate).
  auto Support = [&](std::size_t patch, const OrientedSupportBox &box, const Aabb &aabb,
                     const std::vector<Claim> &claims)
  {
    SpatialSupportBounds support;
    support.patch = patch;
    support.min = aabb.min;
    support.max = aabb.max;
    support.axes = box.axes;
    support.frame_min = box.min;
    support.frame_max = box.max;
    support.claims = claims;
    return support;
  };
  const std::vector<Claim> a_claims = {Claim{1, {-c, -c, 0.0}, {c, c, 0.0}}};
  const std::vector<Claim> b_claims = {Claim{2, {0.5, -0.9, 0.0}, {0.5, -0.1, 0.0}}};
  config::ElectrostaticSolverData::ResponseCorrectionData config;
  config.matching_radius = 1.0;
  for (int k = 0; k < 2; k++)
  {
    Model model;
    model.idx = k + 1;
    model.name = k == 0 ? "corner-a" : "corner-b";
    model.topology = "convex corner";
    model.spatial_basis = true;
    config.models.push_back(model);
    Patch patch;
    patch.model = k + 1;
    patch.weight = 1.0;
    patch.provenance.feature = 20 + k;
    patch.provenance.claims = k == 0 ? a_claims : b_claims;
    config.patches.push_back(patch);
  }
  CHECK(OrientedBoxesInteriorOverlap(OrientedBoxOf(Support(0, a, a_aabb, a_claims)),
                                     OrientedBoxOf(Support(1, near, near_aabb, b_claims)),
                                     3));
  {
    std::vector<GeometryFallbackData> entries;
    GeometryFallbackLedger ledger(GeometryFailurePolicy::RECORD_RAW, false, entries);
    std::vector<std::pair<std::size_t, std::size_t>> overlap_excluded;
    const std::vector<SpatialSupportBounds> supports = {
        Support(0, a, a_aabb, a_claims), Support(1, far, far_aabb, b_claims)};
    CHECK(FindSpatialClusterPriorityExclusions(supports, config, {}, 3, nullptr, &ledger,
                                               &overlap_excluded)
              .empty());
    CHECK(overlap_excluded.empty());
    CHECK(entries.empty());
  }
  {
    std::vector<GeometryFallbackData> entries;
    GeometryFallbackLedger ledger(GeometryFailurePolicy::RECORD_RAW, false, entries);
    std::vector<std::pair<std::size_t, std::size_t>> overlap_excluded;
    const std::vector<SpatialSupportBounds> supports = {
        Support(0, a, a_aabb, a_claims), Support(1, near, near_aabb, b_claims)};
    CHECK(FindSpatialClusterPriorityExclusions(supports, config, {}, 3, nullptr, &ledger,
                                               &overlap_excluded)
              .empty());
    REQUIRE(overlap_excluded.size() == 1);
    CHECK(overlap_excluded.front() == std::make_pair(std::size_t{1}, std::size_t{0}));
    REQUIRE(entries.size() == 1);
    CHECK(entries.front().cls == "OverlapExcluded");
    CHECK(entries.front().patches == std::vector<int>{1, 0});
    CHECK_THAT(entries.front().length, WithinAbs(0.8, 1.0e-12));
  }
  // Without a ledger the same-class pairs are skipped as of record, oriented or not.
  {
    const std::vector<SpatialSupportBounds> supports = {
        Support(0, a, a_aabb, a_claims), Support(1, near, near_aabb, b_claims)};
    CHECK(FindSpatialClusterPriorityExclusions(supports, config, {}, 3).empty());
  }
}

// U-H4a-op (readiness hardening-2 (iv) B2 at the OPERATOR; the r1020 review MAJOR-2,
// decision 1035; the operator-level H4a fixture hardening-1 MINOR-2 named): the
// neighboring-island mesh rotated 45 deg about the plane normal (MakeRotatedIslandMesh)
// under the SA-keyed convex library at R = 0.1 (the facing edges 0.25 apart stay beyond
// 2 R: no cluster, 8 convex corners) whose convex-corner-90 BasisPoints span the cube
// [-h, h]^3 about the vertex (the coupon's matching volume, the U-INV H8 widened-points
// form). Every corner frame lies 45 deg off the mesh axes and all eight cubes are
// co-aligned, so a pair's ORIENTED cubes overlap iff every offset along the island axes is
// below 2 h, while the MESH-ALIGNED boxes (half-width h sqrt 2) overlap whenever the
// rotated offsets are below 2 h sqrt 2. The pair distances before the rotation: the facing
// pair 0.25 along x (2 pairs), adjacent corners 0.5 along an edge (8), the diagonal 0.5 /
// 0.5 (4), facing-to-outer 0.75 along x (4), cross 0.25 / 0.5 (2). With h = 0.1 the facing
// pairs' mesh-aligned boxes overlap (rotated offsets 0.177 < 0.283) and their oriented
// cubes do not (0.25 >= 0.2): NO record on the preflight or the operator, the 8 corners
// applied (the base records the 2 pairs OverlapExcluded on the mesh-aligned test: the
// failing-first of B2 at the operator). With h = 0.2 the facing pairs truly overlap
// (0.25 < 0.4) while 14 more pairs overlap on the mesh-aligned test alone (the adjacent
// 0.354 / 0.354 and the 0.75-along-x 0.53 / 0.53 rotated offsets against 0.566): exactly 2
// OverlapExcluded with the rotated frames on the preflight and the operator (U-INV), the
// subordinate by the larger patch index on equal claims.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Geometry fallback B2 at the operator U-H4a-op: rotated same-class corner "
                 "supports decide on the oriented boxes, preflight == operator (decisions "
                 "1018 Q1 / 1035)",
                 "[surfaceresponseoperator][geometryfallback][invariant][3d][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  constexpr double angle_degrees = 45.0;
  constexpr double R = 0.1;
  // The SA-keyed convex library at R with the corner coupon's basis points on the cube
  // [-h, h]^3 (three rings of four, the fixture's 12-point layout).
  auto WriteLibrary = [&](double h, const std::string &tag)
  {
    const fs::path points = temp.temp_dir / ("uh4aop-corner-cube-" + tag + ".csv");
    const fs::path library_path =
        temp.temp_dir / ("fabrication-process-convex-3d-cube-" + tag + ".json");
    if (Mpi::Root(comm))
    {
      {
        std::ofstream output(points);
        output << "x,y,z\n";
        for (const double z : {-h, 0.0, h})
        {
          output << -h << "," << -h << "," << z << "\n"
                 << h << "," << -h << "," << z << "\n"
                 << h << "," << h << "," << z << "\n"
                 << -h << "," << h << "," << z << "\n";
        }
      }
      // The fixture's surface matrices carry their within-R rows at R (m) = 0.2: the
      // copies re-key that column to this library's R.
      auto SurfaceMatrixAtR = [&](const std::string &source)
      {
        const fs::path copy =
            temp.temp_dir / ("uh4aop-" + fs::path(source).filename().string());
        std::ifstream input(source);
        std::ofstream output(copy);
        std::string line;
        std::getline(input, line);
        output << line << "\n";
        while (std::getline(input, line))
        {
          const auto first = line.find(',');
          const auto second = line.find(',', first + 1);
          const auto third = line.find(',', second + 1);
          REQUIRE(third != std::string::npos);
          output << line.substr(0, second + 1) << R << line.substr(third) << "\n";
        }
        return copy.string();
      };
      json library = ReadJson(convex_library_3d_path);
      library["Name"] = "unit-test-process-convex-3d-cube-" + tag;
      library["MatchingRadius"] = R;
      for (auto &model : library["Models"])
      {
        model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
        for (const char *key : {"FabricatedSurfaceMatrix", "ThinSurfaceMatrix"})
        {
          if (model.contains(key))
          {
            model[key] = SurfaceMatrixAtR(model[key].get<std::string>());
          }
        }
        if (model["Topology"] == "ConvexCorner")
        {
          model["BasisPoints"] = points.string();
        }
      }
      std::ofstream output(library_path);
      output << library.dump(2) << "\n";
    }
    Mpi::Barrier(comm);
    return library_path;
  };
  // The island config on the two-island mesh: both islands terminals, one SA interface
  // over both attributes, the solve's PostprocessOnly construction. The ground is the two
  // walls normal to y (the Cartesian box's attributes 2 / 4) alone: the rotated side walls
  // are off the mesh's bounding box, and as metal they would out-vote the islands in the
  // process-plane normal (the area-weighted principal direction of the metal normals off
  // the bounding box) and read every island edge NonPlanar.
  auto Config = [&](const fs::path &library_path)
  {
    json config = ConvexIslandConfig();
    config["Boundaries"]["Ground"]["Attributes"] = {2, 4};
    config["Boundaries"]["Terminal"] = {{{"Index", 1}, {"Attributes", {9}}},
                                        {{"Index", 2}, {"Attributes", {10}}}};
    config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {9, 10};
    config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"] = {R};
    auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
    correction["Library"] = library_path.string();
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
    Mesh mesh(MakeRotatedIslandMesh(angle_degrees));
    const fs::path manifest_path = dir / "surface-response-requirements.json";
    WriteSurfaceResponseRequirements(iodata, mesh, manifest_path.string());
    Mpi::Barrier(comm);
    return ReadJson(manifest_path);
  };
  auto Operator = [&](const json &config)
  {
    IoData iodata = MakeIoData(config);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeRotatedIslandMesh(angle_degrees)));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    return std::make_pair(response.IsRawFallback(), response.GetGeometryFallbackRecord());
  };
  // The placed convex corners of a manifest: (id, origin, frame axes), every frame off the
  // mesh axes by the rotation.
  struct Corner
  {
    int id;
    std::array<double, 3> origin;
    std::array<std::array<double, 3>, 3> axes;
  };
  auto PlacedCorners = [&](const json &manifest)
  {
    std::vector<Corner> corners;
    const double c = std::cos(angle_degrees * std::acos(-1.0) / 180.0);
    for (const auto &feature : manifest.at("Identification").at("Features"))
    {
      if (feature.at("Type") != "ConvexCorner")
      {
        continue;
      }
      REQUIRE(feature.at("Match").at("Status") == "Matched");
      CHECK(feature.at("Match").at("Model") == "convex-corner-90");
      Corner corner;
      corner.id = feature.at("Id").get<int>();
      corner.origin = feature.at("Frame").at("Origin").get<std::array<double, 3>>();
      corner.axes =
          feature.at("Frame").at("Axes").get<std::array<std::array<double, 3>, 3>>();
      for (int k = 0; k < 3; k++)
      {
        // An in-plane axis has both in-plane components at cos 45 deg; the normal is y.
        if (std::abs(corner.axes[k][1]) > 0.5)
        {
          CHECK_THAT(std::abs(corner.axes[k][1]), WithinAbs(1.0, 1.0e-9));
        }
        else
        {
          CHECK_THAT(std::abs(corner.axes[k][0]), WithinAbs(c, 1.0e-9));
          CHECK_THAT(std::abs(corner.axes[k][2]), WithinAbs(c, 1.0e-9));
        }
      }
      corners.push_back(corner);
    }
    return corners;
  };
  // The pairs whose MESH-ALIGNED boxes of the cube [-h, h]^3 in their frames
  // interior-overlap (the base's decision) and the pairs whose ORIENTED cubes do (the
  // offsets along the first corner's axes, the cubes co-aligned); the ElementBox rule,
  // tolerance 1e-12.
  auto OverlappingPairs = [](const std::vector<Corner> &corners, double h)
  {
    auto HalfWidth = [&](const Corner &corner, int d)
    {
      double width = 0.0;
      for (int k = 0; k < 3; k++)
      {
        width += h * std::abs(corner.axes[k][d]);
      }
      return width;
    };
    std::pair<std::size_t, std::size_t> pairs{0, 0};  // mesh-aligned, oriented
    for (std::size_t i = 0; i < corners.size(); i++)
    {
      for (std::size_t j = i + 1; j < corners.size(); j++)
      {
        bool aligned = true, oriented = true;
        for (int d = 0; d < 3; d++)
        {
          const double lo = std::max(corners[i].origin[d] - HalfWidth(corners[i], d),
                                     corners[j].origin[d] - HalfWidth(corners[j], d));
          const double hi = std::min(corners[i].origin[d] + HalfWidth(corners[i], d),
                                     corners[j].origin[d] + HalfWidth(corners[j], d));
          aligned = aligned && hi - lo > 1.0e-12;
          double offset = 0.0;
          for (int e = 0; e < 3; e++)
          {
            offset += corners[i].axes[d][e] * (corners[j].origin[e] - corners[i].origin[e]);
          }
          oriented = oriented && 2.0 * h - std::abs(offset) > 1.0e-12;
        }
        pairs.first += aligned ? 1 : 0;
        pairs.second += oriented ? 1 : 0;
      }
    }
    return pairs;
  };
  auto Distance = [](const Corner &a, const Corner &b)
  {
    double sum = 0.0;
    for (int d = 0; d < 3; d++)
    {
      sum += (a.origin[d] - b.origin[d]) * (a.origin[d] - b.origin[d]);
    }
    return std::sqrt(sum);
  };

  // h = 0.1: the facing pairs' mesh-aligned boxes overlap, no oriented pair: no record on
  // either side, the 8 corners applied.
  {
    const json config = Config(WriteLibrary(0.1, "h01"));
    const json manifest = Preflight(config, temp.temp_dir / "uh4aop-h01");
    const auto corners = PlacedCorners(manifest);
    REQUIRE(corners.size() == 8);
    CHECK(manifest.at("Summary").at("Counts").at("Missing").get<int>() == 0);
    const auto pairs = OverlappingPairs(corners, 0.1);
    CHECK(pairs.first == 2);
    CHECK(pairs.second == 0);
    CHECK(!manifest.at("Identification").at("Diagnostics").contains("GeometryFallbacks"));
    CHECK(!manifest.at("Summary").contains("GeometryFallbacks"));
    const auto [raw, record] = Operator(config);
    CHECK(!raw);
    CHECK(record.is_null());
  }
  // h = 0.2: the facing pairs' oriented cubes overlap (0.25 < 0.4), 14 more pairs on the
  // mesh-aligned test alone: exactly 2 OverlapExcluded on both, the facing corners.
  {
    const json config = Config(WriteLibrary(0.2, "h02"));
    const json manifest = Preflight(config, temp.temp_dir / "uh4aop-h02");
    const auto corners = PlacedCorners(manifest);
    REQUIRE(corners.size() == 8);
    const auto pairs = OverlappingPairs(corners, 0.2);
    CHECK(pairs.first == 16);
    CHECK(pairs.second == 2);
    const auto &preflight_record =
        manifest.at("Identification").at("Diagnostics").at("GeometryFallbacks");
    REQUIRE(preflight_record.at("Count").get<int>() == 2);
    CHECK(preflight_record.at("Classes").at("OverlapExcluded").get<int>() == 2);
    std::set<int> subordinates;
    for (const auto &entry : preflight_record.at("Entries"))
    {
      CHECK(entry.at("Class") == "OverlapExcluded");
      CHECK(entry.at("Site") == "H4a/ApplyPlacementOwnership:matching-volumes-overlap");
      CHECK(entry.at("What") == "patch");
      CHECK(entry.at("Models") == json::array({"convex-corner-90", "convex-corner-90"}));
      REQUIRE(entry.at("Patches").size() == 2);
      // Equal claimed lengths (R along each arm): the larger patch index is dropped.
      CHECK(entry.at("Patches")[0].get<int>() > entry.at("Patches")[1].get<int>());
      CHECK_THAT(entry.at("Length").get<double>(), WithinAbs(2.0 * R, 1.0e-9));
      REQUIRE(entry.at("Features").size() == 1);
      const int subordinate = entry.at("Features").front().get<int>();
      subordinates.insert(subordinate);
      // The subordinate is a facing corner: its nearest corner lies 0.25 away.
      const auto self = std::find_if(corners.begin(), corners.end(),
                                     [&](const Corner &c) { return c.id == subordinate; });
      REQUIRE(self != corners.end());
      double nearest = 1.0e300;
      for (const auto &other : corners)
      {
        if (other.id != subordinate)
        {
          nearest = std::min(nearest, Distance(*self, other));
        }
      }
      CHECK_THAT(nearest, WithinAbs(0.25, 1.0e-9));
    }
    CHECK(subordinates.size() == 2);
    CHECK(manifest.at("Summary").at("GeometryFallbacks").at("Count").get<int>() == 2);
    const auto [raw, record] = Operator(config);
    CHECK(!raw);
    REQUIRE(!record.is_null());
    CHECK(record.at("Classes").at("OverlapExcluded").get<int>() == 2);
    CHECK(CheckGeometryFallbackInvariant(manifest, record) == 2);
  }
#endif
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

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Geometry fallback H13-I2 BranchingPerimeter: a pinched perimeter is two "
                 "chains at the junction (GUARDED by the extraction's chain rule)",
                 "[surfaceresponseoperator][geometryfallback][h13][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // U-H13b (hardening-2 DESIGN row H13-I2): the identification's former abort "a perimeter
  // chain has a vertex with more than two chain segments" is UNREACHABLE by construction:
  // the extraction's physical chains stop at every non-regular vertex (corners, endpoints,
  // junctions), so the pinch vertex (2, 0.5, 2) where the hole's loop and the corner notch
  // meet with four perimeter segments ends FOUR chains of two segments there, and the
  // identification's per-chain incidence never exceeds two. The conversion in
  // RunGeometryIdentification (the per-chain pass: the incident segments excluded
  // BranchingPerimeter with a record, the chain cut) is the safety net for a chain the
  // extraction would leave branched; this fixture reads it NOT firing on the pinch (0
  // records, the construction continues under both policies), the control alike.
  json config = IslandConfig();
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["UnmatchedPolicy"] = "Warn";
  correction["Library"] =
      WriteSheetLibrary(temp.temp_dir, concave_library_3d_path, "h13").string();
  auto Build =
      [&](const json &run_config, bool pinched, std::vector<std::unique_ptr<Mesh>> &meshes)
  {
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    mfem::Mesh serial = MakeSerialPinchedSheetMesh(pinched);
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    return iodata;
  };
  for (const bool pinched : {false, true})
  {
    CAPTURE(pinched);
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(config, pinched, meshes);
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    CHECK(response.GetGeometryFallbacks().empty());
    CHECK(response.GetPatchCount() > 0);
    bool branching = false;
    for (const auto &portion : response.GetIdentificationExclusionPortions())
    {
      branching = branching || portion.topology == "BranchingPerimeter";
    }
    CHECK_FALSE(branching);
  }
  // FailClosed constructs too (nothing recorded).
  {
    json fail_closed = config;
    fail_closed["Solver"]["Electrostatic"]["ResponseCorrection"]["GeometryFailurePolicy"] =
        "FailClosed";
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(fail_closed, true, meshes);
    LaplaceOperator laplace(iodata, meshes);
    CHECK_NOTHROW(SurfaceResponseOperator(iodata, laplace));
  }
#endif
}

// U-H10 / U-H11c (hardening-2 DESIGN rows H10 / H11c, landed in the identification fix
// round's lane B): the extraction's two former aborts under record_fallbacks. Without the
// flag (the post-processing edge tools) they abort as today; with it the duplicate face of
// another attribute is kept once with the segments it supports flagged, the degenerate face
// is skipped, each recorded by class with its perimeter; the identification turns the
// flagged segments into DuplicateBoundary exclusions and the ledger carries both classes
// (RecordRaw continues, FailClosed aborts by name).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Geometry fallback H10 DuplicateBoundary / H11c DegenerateFace: the "
                 "extraction records instead of aborting under record_fallbacks",
                 "[surfaceresponseoperator][geometryfallback][h10][h11][3d][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json config = IslandConfig();
  config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 10};
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["UnmatchedPolicy"] = "Warn";
  correction["Library"] =
      WriteSheetLibrary(temp.temp_dir, concave_library_3d_path, "h10").string();
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  auto ParallelMesh = [&](SheetDefect defect)
  {
    mfem::Mesh serial = MakeSerialSheetMesh(defect);
    return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
  };
  MetalSurfaceExtraction surface = JointNoiseExtraction(0.2);
  surface.classify_components = true;
  surface.retain_global_faces = true;

  SECTION("H10: coincident boundary elements of different metal attributes")
  {
    auto mesh = ParallelMesh(SheetDefect::DUPLICATE_ATTRIBUTE);
    if (Mpi::Size(Mpi::World()) == 1)
    {
      // The former abort fires on the ROOT's gathered faces alone (the other ranks would
      // wait in the next collective): read in the serial run.
      CHECK_THROWS_WITH(ExtractMetalEdgeGeometry(*mesh, iodata.boundaries, surface),
                        ContainsSubstring("different metal attributes 9 and 10"));
    }
    MetalSurfaceExtraction recording = surface;
    recording.record_fallbacks = true;
    const auto geometry = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries, recording);
    REQUIRE(geometry.fallbacks.size() == 1);
    CHECK(geometry.fallbacks[0].cls == "DuplicateBoundary");
    CHECK(geometry.fallbacks[0].count == 1);
    CHECK_THAT(geometry.fallbacks[0].length, WithinRel(1.0, 1.0e-9));  // 4 x 0.25
    // The duplicate face sits on the hole's rim: its ONE perimeter edge (x = 1, z in
    // [1.25, 1.5], along the hole) is flagged; its interior edges form no segment.
    int flagged = 0;
    for (const auto &segment : geometry.segments)
    {
      flagged += segment.duplicate_boundary ? 1 : 0;
    }
    CHECK(flagged == 1);
    // The identification excludes it (DuplicateBoundary) and records the class; the
    // control sheet records nothing.
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(
        std::make_unique<Mesh>(ParallelMesh(SheetDefect::DUPLICATE_ATTRIBUTE)));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    REQUIRE(response.GetGeometryFallbacks().size() == 1);
    const auto &entry = response.GetGeometryFallbacks().front();
    CHECK(entry.cls == "DuplicateBoundary");
    CHECK(entry.site == "H10/metaledge:coincident-faces");
    CHECK(entry.what == "segment");
    bool listed = false;
    for (const auto &portion : response.GetIdentificationExclusionPortions())
    {
      listed = listed || portion.topology == "DuplicateBoundary";
    }
    CHECK(listed);
    json fail_closed = config;
    fail_closed["Solver"]["Electrostatic"]["ResponseCorrection"]["GeometryFailurePolicy"] =
        "FailClosed";
    IoData closed(fail_closed, false);
    closed.boundaries.cracked_attributes.insert(9);
    std::vector<std::unique_ptr<Mesh>> closed_meshes;
    closed_meshes.push_back(
        std::make_unique<Mesh>(ParallelMesh(SheetDefect::DUPLICATE_ATTRIBUTE)));
    LaplaceOperator closed_laplace(closed, closed_meshes);
    CHECK_THROWS_WITH(SurfaceResponseOperator(closed, closed_laplace),
                      ContainsSubstring("GeometryFailurePolicy FailClosed") &&
                          ContainsSubstring("DuplicateBoundary"));
  }

  SECTION("H11c: a face degenerating to fewer than three distinct vertices")
  {
    auto mesh = ParallelMesh(SheetDefect::DEGENERATE_FACE);
    if (Mpi::Size(Mpi::World()) == 1)
    {
      CHECK_THROWS_WITH(ExtractMetalEdgeGeometry(*mesh, iodata.boundaries, surface),
                        ContainsSubstring("fewer than three distinct vertices"));
    }
    MetalSurfaceExtraction recording = surface;
    recording.record_fallbacks = true;
    const auto geometry = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries, recording);
    REQUIRE(geometry.fallbacks.size() == 1);
    CHECK(geometry.fallbacks[0].cls == "DegenerateFace");
    CHECK(geometry.fallbacks[0].count >= 1);
    CHECK(geometry.fallbacks[0].length > 0.0);
    CHECK_FALSE(geometry.segments.empty());
  }

  SECTION("the control sheet: nothing recorded; H11a / H11b GUARDED")
  {
    MetalSurfaceExtraction recording = surface;
    recording.record_fallbacks = true;
    auto mesh = ParallelMesh(SheetDefect::NONE);
    const auto geometry = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries, recording);
    CHECK(geometry.fallbacks.empty());
    // H11a / H11b (UnorientedEdge): the frame builders' `unframed` form is bitwise the
    // abort form wherever a frame exists, and flags nothing on a regular sheet. No fixture
    // reaches the zero-norm sites: two coplanar faces of opposite inward direction form an
    // INTERIOR edge (no segment) in the extraction, a near-cancellation a fold; the sites
    // are GUARDED by the extraction's edge classification and recorded so (hardening-2
    // DESIGN row H11b's clause).
    std::vector<std::size_t> physical;
    for (std::size_t i = 0; i < geometry.segments.size(); i++)
    {
      if (geometry.segments[i].type == MetalEdgeSegmentType::PHYSICAL)
      {
        physical.push_back(i);
      }
    }
    REQUIRE_FALSE(physical.empty());
    auto Score = [](int) { return 1.0; };
    std::vector<unsigned char> unframed_normals, unframed_gaps;
    const auto normals = BuildMetalEdgeProcessNormals(*mesh, geometry, physical, Score,
                                                      std::nullopt, nullptr, nullptr);
    const auto flagged_normals = BuildMetalEdgeProcessNormals(
        *mesh, geometry, physical, Score, std::nullopt, nullptr, &unframed_normals);
    CHECK(flagged_normals == normals);
    CHECK(std::count(unframed_normals.begin(), unframed_normals.end(), 1) == 0);
    const auto gaps = BuildMetalEdgeGapDirections(*mesh, geometry, physical, normals);
    const auto flagged_gaps =
        BuildMetalEdgeGapDirections(*mesh, geometry, physical, normals, &unframed_gaps);
    CHECK(flagged_gaps == gaps);
    CHECK(std::count(unframed_gaps.begin(), unframed_gaps.end(), 1) == 0);
  }
#endif
}
