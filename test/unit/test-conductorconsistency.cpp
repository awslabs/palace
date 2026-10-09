// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "linalg/ksp.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace fs = std::filesystem;
using json = nlohmann::json;
using namespace Catch::Matchers;

namespace
{

// A two-conductor spatial coupon hand-placed across the gap between the two islands of
// MakeIslandMesh(neighboring_island): island 1 (attribute 9, terminal 1) occupies x in
// [0.375, 0.875], island 2 (attribute 10, terminal 2) x in [1.125, 1.625], both z in [0.75,
// 1.25] on the plane y = 0.5 of the 2 x 1 x 2 box (h = 0.125). The coupon frame: u = +z, v
// = +x, w = +y (the process normal), origin (1.0, 0.5, 1.0) = the gap centre on the plane;
// box u in [-0.2, 0.2], v in [-0.25, 0.25], w in [-0.2, 0.2]. Its trace mesh: a ten-point
// perimeter ring at w = -0.2 (free knots 1-10), at w = 0 (the metal cross-sections: eight
// conductor vertices — conductor 1 at v <= -0.125 on island 1's metal, conductor 2 at v >=
// 0.125 on island 2's — and the two gap-centre points (+-0.2, 0, 0)) and at w = +0.2 (free
// knots 11-20), two cap-centre hats (InteriorTraceCount 2). The CONSISTENT coupon keeps
// both gap-centre points free (knots 21, 22); the INCONSISTENT coupon declares the one at
// u = -0.2 a conductor-1 vertex (knot 21 is then the other): coupon metal at the island-1
// potential on a face where the device has gap at an intermediate potential, the S1p
// fictitious-block pattern (decisions 271-277). The OFF-PLANE variant moves every
// conductor-2 vertex to w = 0.1 (no plane vertex: the gate cannot probe it and fails
// closed). The NO-CONDUCTOR variant declares every plane ring point free (a coupon without
// metal cross-sections at a fixed potential, like a finite-impedance coupon or a ring-path
// model without ZeroTraceIndices: nothing to probe, recorded untestable under
// UnprobedModels and warned, not gated; decision 279 MINOR-2). The CHORD-OFFSET variant
// (decision 611 (2)) moves the two conductor-1 metal-edge knots on the faces u = -+0.2 by
// kChordOffset into the gap (v = -0.125 + 0.02: device x = 0.895, 0.02 off island 1's edge
// at x = 0.875): the 32b0 pattern of a coupon arc face end built on the fitted circle up to
// the chord sagitta outside the device's chord polyline, which reads the gap potential next
// to the edge (fail-before: excluded). The CHORD-OFFSET-SNAPPED variant is the same coupon
// whose placed patch carries the probe snaps of those two knots (displacement -0.02 along
// AxisV = -x, as BuildFeaturePatches records them for an arc chord): the probes are read
// on the device edge, the reading is the Dirichlet value, the snaps are recorded.
enum class CouponVariant
{
  CONSISTENT,
  INCONSISTENT,
  CONDUCTOR_OFF_PLANE,
  NO_CONDUCTOR,
  CHORD_OFFSET,
  CHORD_OFFSET_SNAPPED
};

struct CouponFiles
{
  fs::path basis_points, trace_vertices, trace_triangles, fabricated, thin, cache;
  int contour_size = 0;
  std::vector<std::array<double, 3>> perimeter;  // (u, v) of the ten ring points
  std::vector<int> offset_vertices;              // the chord-offset knots (0-based)
};

constexpr double kHalfU = 0.2, kHalfV = 0.25, kHalfW = 0.2, kMetalEdgeV = 0.125;
constexpr double kChordOffset = 0.02;

CouponFiles WriteCoupon(const fs::path &dir, CouponVariant variant)
{
  const std::string tag = variant == CouponVariant::CONSISTENT            ? "consistent"
                          : variant == CouponVariant::INCONSISTENT        ? "inconsistent"
                          : variant == CouponVariant::CONDUCTOR_OFF_PLANE ? "off-plane"
                          : variant == CouponVariant::NO_CONDUCTOR        ? "no-conductor"
                          : variant == CouponVariant::CHORD_OFFSET        ? "chord-offset"
                                                                          : "chord-snapped";
  const bool chord_offset = variant == CouponVariant::CHORD_OFFSET ||
                            variant == CouponVariant::CHORD_OFFSET_SNAPPED;
  CouponFiles files;
  files.basis_points = dir / ("conductor-consistency-" + tag + "-basis-points.csv");
  files.trace_vertices = dir / ("conductor-consistency-" + tag + "-trace-vertices.csv");
  files.trace_triangles = dir / ("conductor-consistency-" + tag + "-trace-triangles.csv");
  files.fabricated = dir / ("conductor-consistency-" + tag + "-fabricated.csv");
  files.thin = dir / ("conductor-consistency-" + tag + "-thin.csv");
  files.cache = dir / ("conductor-consistency-" + tag + "-geometry.json");
  // The perimeter ring (u, v), counter-clockwise seen from +w.
  files.perimeter = {{-kHalfU, -kHalfV},     {-kHalfU, -kMetalEdgeV},
                     {-kHalfU, 0.0},         {-kHalfU, kMetalEdgeV},
                     {-kHalfU, kHalfV},      {kHalfU, kHalfV},
                     {kHalfU, kMetalEdgeV},  {kHalfU, 0.0},
                     {kHalfU, -kMetalEdgeV}, {kHalfU, -kHalfV}};
  const int ring = static_cast<int>(files.perimeter.size());
  struct Vertex
  {
    std::array<double, 3> point;
    int basis;
    int conductor;
  };
  std::vector<Vertex> vertices;
  std::vector<std::array<double, 3>> basis_points;
  auto AddFree = [&](std::array<double, 3> point)
  {
    basis_points.push_back(point);
    vertices.push_back({point, static_cast<int>(basis_points.size()), 0});
    return static_cast<int>(vertices.size());
  };
  auto AddConductor = [&](std::array<double, 3> point, int conductor)
  {
    vertices.push_back({point, 0, conductor});
    return static_cast<int>(vertices.size());
  };
  std::vector<int> lower(ring), plane(ring), upper(ring);
  for (int i = 0; i < ring; i++)
  {
    lower[i] = AddFree({files.perimeter[i][0], files.perimeter[i][1], -kHalfW});
  }
  for (int i = 0; i < ring; i++)
  {
    upper[i] = AddFree({files.perimeter[i][0], files.perimeter[i][1], kHalfW});
  }
  std::vector<int> deferred_plane_free;
  for (int i = 0; i < ring; i++)
  {
    const double u = files.perimeter[i][0], v = files.perimeter[i][1];
    if (variant == CouponVariant::NO_CONDUCTOR)
    {
      plane[i] = AddFree({u, v, 0.0});
    }
    else if (v <= -kMetalEdgeV + 1.0e-12)
    {
      // The chord-offset knots: the metal-edge corners of the two faces, into the gap.
      const bool offset = chord_offset && std::abs(v + kMetalEdgeV) <= 1.0e-12;
      plane[i] = AddConductor({u, offset ? v + kChordOffset : v, 0.0}, 1);
      if (offset)
      {
        files.offset_vertices.push_back(plane[i] - 1);
      }
    }
    else if (v >= kMetalEdgeV - 1.0e-12)
    {
      plane[i] = AddConductor(
          {u, v, variant == CouponVariant::CONDUCTOR_OFF_PLANE ? 0.1 : 0.0}, 2);
    }
    else if (variant == CouponVariant::INCONSISTENT && u < 0.0)
    {
      plane[i] = AddConductor({u, v, 0.0}, 1);  // the fictitious metal
    }
    else
    {
      plane[i] = AddFree({u, v, 0.0});
    }
  }
  const int lower_cap = AddFree({0.0, 0.0, -kHalfW});
  const int upper_cap = AddFree({0.0, 0.0, kHalfW});
  std::vector<std::array<int, 3>> triangles;
  for (int i = 0; i < ring; i++)
  {
    const int next = (i + 1) % ring;
    triangles.push_back({lower[i], lower[next], plane[next]});
    triangles.push_back({lower[i], plane[next], plane[i]});
    triangles.push_back({plane[i], plane[next], upper[next]});
    triangles.push_back({plane[i], upper[next], upper[i]});
    triangles.push_back({lower_cap, lower[next], lower[i]});
    triangles.push_back({upper_cap, upper[i], upper[next]});
  }
  files.contour_size = static_cast<int>(basis_points.size());
  if (Mpi::Root(Mpi::World()))
  {
    {
      std::ofstream output(files.basis_points);
      output << "x,y,z\n";
      for (const auto &point : basis_points)
      {
        output << point[0] << "," << point[1] << "," << point[2] << "\n";
      }
    }
    {
      std::ofstream output(files.trace_vertices);
      output << "vertex,x,y,z,basis,conductor\n";
      for (std::size_t i = 0; i < vertices.size(); i++)
      {
        output << i + 1 << "," << vertices[i].point[0] << "," << vertices[i].point[1] << ","
               << vertices[i].point[2] << "," << vertices[i].basis << ","
               << vertices[i].conductor << "\n";
      }
    }
    {
      std::ofstream output(files.trace_triangles);
      output << "triangle,vertex_i,vertex_j,vertex_k\n";
      for (std::size_t i = 0; i < triangles.size(); i++)
      {
        output << i + 1 << "," << triangles[i][0] << "," << triangles[i][1] << ","
               << triangles[i][2] << "\n";
      }
    }
    // Diagonal response matrices (the gate reads the trace, not the energies): the
    // fabricated domain energy 2x the thin one.
    const int basis_size =
        files.contour_size + (variant == CouponVariant::NO_CONDUCTOR ? 0 : 1);
    for (const auto &[path, scale] :
         {std::make_pair(files.fabricated, 2.0e-12), std::make_pair(files.thin, 1.0e-12)})
    {
      std::ofstream output(path);
      output << "basis_i,basis_j,Q_ij (J)\n";
      for (int i = 1; i <= basis_size; i++)
      {
        for (int j = i; j <= basis_size; j++)
        {
          output << i << "," << j << "," << (i == j ? scale : 0.0) << "\n";
        }
      }
    }
    // The placed geometry (the response-geometry cache replaces the library matching).
    config::ElectrostaticSolverData::ResponseCorrectionData data;
    data.matching_radius = kHalfU;
    auto &model = data.models.emplace_back();
    model.idx = 1;
    model.name = "gap-junction-" + tag;
    model.topology = "SpatialEdgeCluster";
    model.fabricated_matrix = files.fabricated.string();
    model.thin_matrix = files.thin.string();
    model.basis_points = files.basis_points.string();
    model.trace_vertices = files.trace_vertices.string();
    model.trace_triangles = files.trace_triangles.string();
    model.spatial_basis = true;
    model.interior_trace_count = 2;
    model.conductor_state_count = variant == CouponVariant::NO_CONDUCTOR ? 0 : 1;
    auto &patch = data.patches.emplace_back();
    patch.model = 1;
    patch.origin = {1.0, 0.5, 1.0};
    patch.axis_u = {0.0, 0.0, 1.0};
    patch.axis_v = {1.0, 0.0, 0.0};
    patch.axis_w = {0.0, 1.0, 0.0};
    patch.conductor_references = {{0.0, -0.2, 0.0}, {0.0, 0.2, 0.0}};
    if (variant == CouponVariant::NO_CONDUCTOR)
    {
      patch.conductor_references.pop_back();
    }
    patch.weight = 1.0;
    patch.provenance.feature = 7;
    // Two claimed portions of 0.4 along the islands' facing edges (z in [0.8, 1.2]).
    patch.provenance.claims = {{3, {0.875, 0.5, 0.8}, {0.875, 0.5, 1.2}},
                               {5, {1.125, 0.5, 0.8}, {1.125, 0.5, 1.2}}};
    if (variant == CouponVariant::CHORD_OFFSET_SNAPPED)
    {
      // The placement's record for a knot kChordOffset off a device arc chord (segment 3
      // of arc 0): back onto the edge along -x (AxisV points +x).
      for (const int vertex : files.offset_vertices)
      {
        patch.provenance.probe_snaps.push_back(
            {vertex, 3, 0, {-kChordOffset, 0.0, 0.0}, kChordOffset});
      }
    }
    WriteResponseGeometryCache(files.cache, data);
  }
  Mpi::Barrier(Mpi::World());
  return files;
}

}  // namespace

// Decision 277 (A): the conductor-consistency gate compares, at solve time, the device
// potential at every conductor vertex of a spatial coupon's trace mesh on the process plane
// (the coupon's metal cross-sections on its box faces, on the device's metal sheet) with
// the potential at that conductor's reference, relative to the patch's trace amplitude.
// Real metal reads exactly the Dirichlet value (ratio ~1e-16); coupon metal where the
// device has gap reads the gap potential (here ~0.6 of the 1-V amplitude at the gap centre
// between the islands; S1p 0.20-0.43). Above kConductorConsistencyTolerance the patch is
// excluded like a DomainBoundary cell (weight 0: fixed-trace energies, self-consistent
// operator, patch table) and recorded with its claimed length; the decision is the same on
// 1 and 2 ranks
// ([Parallel]) and the record is replicated. The free knots adjacent to the plane
// conductor vertices are recorded only (their coefficient carries the near-edge field; see
// the calibration in conductor-consistency-20261003/gate/REPORT.md).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator conductor-consistency gate",
                 "[surfaceresponseoperator][conductorconsistency][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json config = HighOrderSpatialConfig();
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["TraceCoupling"] = "SurfaceMortar";
  correction["MortarOversampling"] = 2;
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh(false, false, false, true)));
  LaplaceOperator laplace(iodata, meshes);
  auto K = laplace.GetStiffnessMatrix();
  KspSolver ksp(iodata, laplace.GetH1Spaces());
  ksp.SetOperators(*K, *K);
  // The raw thin-metal fields of the two terminals (island 1 at 1 V, then island 2).
  std::vector<Vector> V(2);
  for (int source = 1; source <= 2; source++)
  {
    Vector RHS;
    laplace.GetExcitationVector(source, *K, V[source - 1], RHS);
    ksp.Mult(RHS, V[source - 1]);
    REQUIRE(ksp.GetConverged());
  }
  const auto tolerance = SurfaceResponseOperator::kConductorConsistencyTolerance;
  CHECK(tolerance == 0.02);
  const auto floor = SurfaceResponseOperator::kConductorConsistencyAmplitudeFloor;
  CHECK(floor == 1.0e-3);

  auto Replicated = [](double value)
  {
    double minimum = value, maximum = value;
    Mpi::GlobalMin(1, &minimum, Mpi::World());
    Mpi::GlobalMax(1, &maximum, Mpi::World());
    CHECK(minimum == maximum);
    return minimum;
  };

  SECTION("a consistent coupon is applied: its metal knots read the conductor potential")
  {
    const auto files = WriteCoupon(temp.temp_dir, CouponVariant::CONSISTENT);
    test::GeometryCacheEnvGuard cache_env(files.cache.string(), false);
    SurfaceResponseOperator response(iodata, laplace);
    REQUIRE(response.GetPatchCount() == 1);
    REQUIRE(response.GetBasisSize() == files.contour_size + 1);
    REQUIRE(response.HasConductorConsistencyProbes());
    const auto before = response.GetElectrostaticResponse(V[0]);
    CHECK(before.domain_correction != 0.0);
    for (int source = 1; source <= 2; source++)
    {
      const auto records = response.ApplyConductorConsistencyGate(V[source - 1], source);
      REQUIRE(records.size() == 1);
      const auto &record = records.front();
      CHECK(record.patch == 0);
      CHECK(record.model == 1);
      CHECK(record.source == source);
      CHECK(record.plane_knots == 8);  // the metal ring points, four per conductor
      CHECK(record.off_plane_knots == 0);
      // Every plane conductor vertex has its column neighbours at w = -+0.2 (sharing a
      // strip-triangle edge), each next to one conductor only: 16 adjacent knots.
      CHECK(record.adjacent_knots == 16);
      CHECK_THAT(record.amplitude, WithinRel(1.0, 1.0e-6));  // the 1-V terminal
      CHECK_THAT(record.state, WithinRel(1.0, 1.0e-6));
      CHECK_THAT(record.excitation_potential, WithinRel(1.0, 1.0e-12));
      CHECK(record.normalization == record.amplitude);  // far above the floor
      // The plane knots lie on the islands' Dirichlet sheets: the device potential there
      // is the terminal value to roundoff.
      CHECK(record.max_deviation < 1.0e-10);
      CHECK(record.max_ratio < 1.0e-10);
      CHECK_FALSE(record.excluded);
      CHECK(record.adjacent_max_ratio > 0.01);  // information: the far rows differ
      CHECK(record.adjacent_max_ratio <= 1.0 + 1.0e-12);
      CHECK_THAT(record.claim_length, WithinAbs(0.8, 1.0e-12));
      CHECK(record.cell_length == 0.0);
      CHECK_THAT(Replicated(record.max_ratio), WithinAbs(record.max_ratio, 0.0));
    }
    const auto after = response.GetElectrostaticResponse(V[0]);
    CHECK(after.domain_correction == before.domain_correction);
    const auto statistics = response.GetStatistics();
    const auto &record = statistics["Diagnostics"]["ConductorConsistency"];
    CHECK(record["Count"].get<int>() == 0);
    CHECK(record["TestedPatches"].get<int>() == 2);
    CHECK(record["Records"].size() == 2);
    CHECK(record["ExcludedPatches"].empty());
    CHECK(record["ClaimLength"].get<double>() == 0.0);
    CHECK_THAT(record["Tolerance"].get<double>(), WithinAbs(tolerance, 0.0));
    CHECK(statistics["ModelCatalog"][0]["PatchWeight"].get<double>() == 1.0);
    REQUIRE(response.GetPatchAssignments().size() == (Mpi::Root(Mpi::World()) ? 1 : 0));
    if (Mpi::Root(Mpi::World()))
    {
      CHECK(response.GetPatchAssignments()[0].weight == 1.0);
      CHECK(response.GetPatchAssignments()[0].global_index == 0);
    }
  }

  SECTION("coupon metal over a device gap at another potential is excluded and recorded")
  {
    const auto files = WriteCoupon(temp.temp_dir, CouponVariant::INCONSISTENT);
    test::GeometryCacheEnvGuard cache_env(files.cache.string(), false);
    SurfaceResponseOperator response(iodata, laplace);
    REQUIRE(response.GetPatchCount() == 1);
    REQUIRE(response.GetBasisSize() == files.contour_size + 1);
    const auto before = response.GetElectrostaticResponse(V[0]);
    CHECK(before.domain_correction != 0.0);
    const auto records = response.ApplyConductorConsistencyGate(V[0], 1);
    REQUIRE(records.size() == 1);
    const auto &record = records.front();
    CHECK(record.patch == 0);
    CHECK(record.plane_knots == 9);
    // The fictitious conductor-1 vertex at the gap centre of the u = -0.2 face (device
    // (1.0, 0.5, 0.8)) reads the gap potential between the 1-V island and the grounded
    // island, far from both the conductor's 1 V and the tolerance; its column neighbours
    // touch one conductor, so the adjacent set grows by two.
    CHECK(record.adjacent_knots == 18);
    CHECK(record.worst_conductor == 1);
    CHECK_THAT(record.worst_point[0], WithinAbs(1.0, 1.0e-9));
    CHECK_THAT(record.worst_point[1], WithinAbs(0.5, 1.0e-9));
    CHECK_THAT(record.worst_point[2], WithinAbs(0.8, 1.0e-9));
    CHECK(record.max_ratio > 0.2);
    CHECK(record.max_ratio < 1.0);
    CHECK(record.max_ratio > tolerance);
    CHECK(record.normalization == record.amplitude);
    CHECK_THAT(record.excitation_potential, WithinRel(1.0, 1.0e-12));
    CHECK(record.excluded);
    CHECK_THAT(Replicated(record.max_ratio), WithinAbs(record.max_ratio, 0.0));
    CHECK_THAT(record.claim_length, WithinAbs(0.8, 1.0e-12));
    // Excluded like a DomainBoundary cell: weight 0 everywhere the patch acts.
    const auto after = response.GetElectrostaticResponse(V[0]);
    CHECK(after.domain_correction == 0.0);
    const auto statistics = response.GetStatistics();
    const auto &diagnostics = statistics["Diagnostics"]["ConductorConsistency"];
    CHECK(diagnostics["Count"].get<int>() == 1);
    CHECK(diagnostics["TestedPatches"].get<int>() == 1);
    REQUIRE(diagnostics["ExcludedPatches"].size() == 1);
    const auto &entry = diagnostics["ExcludedPatches"][0];
    CHECK(entry["Patch"].get<int>() == 0);
    CHECK(entry["Model"] == "gap-junction-inconsistent");
    CHECK(entry["Source"].get<int>() == 1);
    CHECK(entry["Excluded"].get<bool>());
    CHECK_THAT(entry["MaxRatio"].get<double>(), WithinAbs(record.max_ratio, 0.0));
    CHECK_THAT(entry["MaxRatioOverState"].get<double>(),
               WithinRel(record.max_deviation / record.state, 1.0e-12));
    CHECK_THAT(entry["Normalization"].get<double>(), WithinAbs(record.amplitude, 0.0));
    CHECK_FALSE(entry["FloorApplied"].get<bool>());
    CHECK_THAT(diagnostics["AmplitudeFloor"].get<double>(), WithinAbs(floor, 0.0));
    CHECK_THAT(entry["ClaimLength"].get<double>(), WithinAbs(0.8, 1.0e-12));
    CHECK_THAT(diagnostics["ClaimLength"].get<double>(), WithinAbs(0.8, 1.0e-12));
    CHECK(diagnostics["CellLength"].get<double>() == 0.0);
    CHECK(statistics["ModelCatalog"][0]["PatchWeight"].get<double>() == 0.0);
    if (Mpi::Root(Mpi::World()))
    {
      REQUIRE(response.GetPatchAssignments().size() == 1);
      CHECK(response.GetPatchAssignments()[0].weight == 0.0);
    }
    // The self-consistent operator omits the excluded patch too.
    Vector y(V[0].Size());
    response.Mult(V[0], y);
    CHECK(linalg::Norml2(Mpi::World(), y) == 0.0);
    // The exclusion is sticky: the next excitation finds no applied probe patch.
    const auto later = response.ApplyConductorConsistencyGate(V[1], 2);
    CHECK(later.empty());
    CHECK(statistics["Diagnostics"]["ConductorConsistency"]["Records"].size() == 1);
  }

  SECTION("the ratio is scale invariant and floored at a fraction of the excitation")
  {
    const auto files = WriteCoupon(temp.temp_dir, CouponVariant::INCONSISTENT);
    test::GeometryCacheEnvGuard cache_env(files.cache.string(), false);
    // The unit excitation's reading (a fresh operator: the exclusion is sticky).
    double unit_ratio = 0.0, unit_deviation = 0.0;
    {
      SurfaceResponseOperator response(iodata, laplace);
      const auto records = response.ApplyConductorConsistencyGate(V[0], 1);
      REQUIRE(records.size() == 1);
      unit_ratio = records.front().max_ratio;
      unit_deviation = records.front().max_deviation;
      REQUIRE(unit_ratio > tolerance);
    }
    // A 1e-4-V terminal: the same ratio (the floor 1e-3 x 1e-4 lies under the amplitude),
    // the same decision.
    {
      Vector scaled(V[0]);
      scaled *= 1.0e-4;
      SurfaceResponseOperator response(iodata, laplace);
      const auto records = response.ApplyConductorConsistencyGate(scaled, 1);
      REQUIRE(records.size() == 1);
      const auto &record = records.front();
      CHECK_THAT(record.excitation_potential, WithinRel(1.0e-4, 1.0e-12));
      CHECK_THAT(record.amplitude, WithinRel(1.0e-4 * 1.018, 0.05));
      CHECK(record.normalization == record.amplitude);
      CHECK_THAT(record.max_deviation, WithinRel(1.0e-4 * unit_deviation, 1.0e-9));
      CHECK_THAT(record.max_ratio, WithinRel(unit_ratio, 1.0e-9));
      CHECK(record.excluded);
    }
    // The coupon under an excitation it barely sees: both islands within 1e-6 V of each
    // other at ~1 V (the device's largest potential). The trace amplitude is 1e-6 of the
    // excitation and the fictitious knot's deviation 0.6e-6 V: without the floor the
    // ratio would be the unit one (0.63) and the patch excluded on a 1e-6-V reading; with
    // the floor the normalization is 1e-3 V, the ratio 6e-4, recorded FloorApplied, not
    // excluded (energy ~ amplitude^2: at most 1e-12 of the unit patch), and the next
    // excitation that fields the coupon excludes it.
    {
      Vector shifted(V[0]);
      shifted *= 1.0e-6;
      shifted += 1.0;
      SurfaceResponseOperator response(iodata, laplace);
      const auto records = response.ApplyConductorConsistencyGate(shifted, 1);
      REQUIRE(records.size() == 1);
      const auto &record = records.front();
      CHECK_THAT(record.excitation_potential, WithinRel(1.0 + 1.0e-6, 1.0e-12));
      CHECK_THAT(record.amplitude, WithinRel(1.0e-6 * 1.018, 0.05));
      CHECK_THAT(record.normalization, WithinRel(floor * (1.0 + 1.0e-6), 1.0e-12));
      CHECK(record.normalization > record.amplitude);
      CHECK_THAT(record.max_deviation, WithinRel(1.0e-6 * unit_deviation, 1.0e-6));
      CHECK_THAT(record.max_ratio,
                 WithinRel(record.max_deviation / record.normalization, 1.0e-12));
      CHECK(record.max_ratio < tolerance);
      CHECK(record.max_ratio > 1.0e-4);
      CHECK_FALSE(record.excluded);
      CHECK_THAT(Replicated(record.max_ratio), WithinAbs(record.max_ratio, 0.0));
      const auto statistics = response.GetStatistics();
      const auto &entry = statistics["Diagnostics"]["ConductorConsistency"]["Records"][0];
      CHECK(entry["FloorApplied"].get<bool>());
      CHECK_THAT(entry["Normalization"].get<double>(),
                 WithinAbs(record.normalization, 0.0));
      CHECK_THAT(entry["ExcitationPotential"].get<double>(),
                 WithinAbs(record.excitation_potential, 0.0));
      CHECK(statistics["Diagnostics"]["ConductorConsistency"]["Count"].get<int>() == 0);
      CHECK(response.GetElectrostaticResponse(V[0]).domain_correction != 0.0);
      const auto later = response.ApplyConductorConsistencyGate(V[0], 2);
      REQUIRE(later.size() == 1);
      CHECK_THAT(later.front().max_ratio, WithinRel(unit_ratio, 1.0e-9));
      CHECK(later.front().excluded);
      CHECK(response.GetElectrostaticResponse(V[0]).domain_correction == 0.0);
    }
  }

  SECTION("a metal-edge knot off the device chord reads the gap; snapped onto the chord it "
          "reads the conductor (decision 611)")
  {
    // Fail-before (the gate as on main): the two conductor-1 face knots 0.02 into the gap
    // read the device potential 0.02 from island 1's edge - the near-edge gap potential,
    // far above the tolerance - so the coupon is excluded; no snap is recorded.
    double offset_ratio = 0.0;
    {
      const auto files = WriteCoupon(temp.temp_dir, CouponVariant::CHORD_OFFSET);
      REQUIRE(files.offset_vertices.size() == 2);
      test::GeometryCacheEnvGuard cache_env(files.cache.string(), false);
      SurfaceResponseOperator response(iodata, laplace);
      REQUIRE(response.GetPatchCount() == 1);
      const auto records = response.ApplyConductorConsistencyGate(V[0], 1);
      REQUIRE(records.size() == 1);
      const auto &record = records.front();
      CHECK(record.plane_knots == 8);
      CHECK(record.worst_conductor == 1);
      CHECK(std::find(files.offset_vertices.begin(), files.offset_vertices.end(),
                      record.worst_vertex) != files.offset_vertices.end());
      CHECK_THAT(record.worst_point[0], WithinAbs(0.875 + kChordOffset, 1.0e-9));
      CHECK_THAT(record.worst_point[1], WithinAbs(0.5, 1.0e-9));
      CHECK(record.max_ratio > tolerance);
      CHECK(record.max_ratio < 1.0);
      CHECK(record.excluded);
      CHECK(record.probe_snaps.empty());
      CHECK(record.max_probe_snap == 0.0);
      offset_ratio = Replicated(record.max_ratio);
      const auto statistics = response.GetStatistics();
      const auto &entry = statistics["Diagnostics"]["ConductorConsistency"]["Records"][0];
      CHECK(entry["ProbeSnaps"]["Count"].get<int>() == 0);
      CHECK(entry["ProbeSnaps"]["MaxDisplacement"].get<double>() == 0.0);
      CHECK(entry["ProbeSnaps"]["Knots"].empty());
      CHECK(response.GetElectrostaticResponse(V[0]).domain_correction == 0.0);
    }
    // The fix: the same coupon with the placement's probe snaps - the two knots are read
    // on the device edge (their mapped position plus the recorded displacement), the
    // reading is the Dirichlet value to roundoff, the patch stays applied, and the record
    // carries the snaps (count, largest displacement, the knots).
    {
      const auto files = WriteCoupon(temp.temp_dir, CouponVariant::CHORD_OFFSET_SNAPPED);
      REQUIRE(files.offset_vertices.size() == 2);
      test::GeometryCacheEnvGuard cache_env(files.cache.string(), false);
      SurfaceResponseOperator response(iodata, laplace);
      REQUIRE(response.GetPatchCount() == 1);
      for (int source = 1; source <= 2; source++)
      {
        const auto records = response.ApplyConductorConsistencyGate(V[source - 1], source);
        REQUIRE(records.size() == 1);
        const auto &record = records.front();
        CHECK(record.plane_knots == 8);
        CHECK(record.max_deviation < 1.0e-10);
        CHECK(record.max_ratio < 1.0e-10);
        CHECK(record.max_ratio < 1.0e-6 * offset_ratio);
        CHECK_FALSE(record.excluded);
        REQUIRE(record.probe_snaps.size() == 2);
        CHECK_THAT(record.max_probe_snap, WithinAbs(kChordOffset, 1.0e-15));
        for (const auto &snap : record.probe_snaps)
        {
          CHECK(std::find(files.offset_vertices.begin(), files.offset_vertices.end(),
                          snap.vertex) != files.offset_vertices.end());
          CHECK(snap.segment == 3);
          CHECK(snap.arc == 0);
          CHECK_THAT(snap.displacement, WithinAbs(kChordOffset, 1.0e-15));
        }
        CHECK_THAT(Replicated(record.max_ratio), WithinAbs(record.max_ratio, 0.0));
        CHECK_THAT(Replicated(record.max_probe_snap), WithinAbs(kChordOffset, 1.0e-15));
      }
      CHECK(response.GetElectrostaticResponse(V[0]).domain_correction != 0.0);
      const auto statistics = response.GetStatistics();
      const auto &diagnostics = statistics["Diagnostics"]["ConductorConsistency"];
      CHECK(diagnostics["Count"].get<int>() == 0);
      REQUIRE(diagnostics["Records"].size() == 2);
      CHECK_THAT(diagnostics["Rule"].get<std::string>(), ContainsSubstring("ProbeSnaps"));
      for (const auto &entry : diagnostics["Records"])
      {
        const auto &snaps = entry["ProbeSnaps"];
        CHECK(snaps["Count"].get<int>() == 2);
        CHECK_THAT(snaps["MaxDisplacement"].get<double>(),
                   WithinAbs(kChordOffset, 1.0e-15));
        REQUIRE(snaps["Knots"].size() == 2);
        std::vector<int> knots;
        for (const auto &knot : snaps["Knots"])
        {
          knots.push_back(knot["Knot"].get<int>());
          CHECK(knot["Segment"].get<int>() == 3);
          CHECK(knot["Arc"].get<int>() == 0);
          CHECK_THAT(knot["Displacement"].get<double>(), WithinAbs(kChordOffset, 1.0e-15));
        }
        std::sort(knots.begin(), knots.end());
        auto expected = files.offset_vertices;
        std::sort(expected.begin(), expected.end());
        CHECK(knots == expected);
      }
      CHECK(statistics["ModelCatalog"][0]["PatchWeight"].get<double>() == 1.0);
    }
    // A snap record naming a knot that is not a conductor vertex of the trace mesh (a
    // stale cache or another library) fails closed by the model's name.
    {
      const auto files = WriteCoupon(temp.temp_dir, CouponVariant::CHORD_OFFSET_SNAPPED);
      auto data = ReadResponseGeometryCache(
          files.cache, *iodata.solver.electrostatic.response_correction);
      REQUIRE(data.patches.size() == 1);
      data.patches[0].provenance.probe_snaps.push_back({0, 3, 0, {0.0, 0.0, 0.0}, 0.0});
      const auto stale = temp.temp_dir / "conductor-consistency-chord-stale.json";
      if (Mpi::Root(Mpi::World()))
      {
        WriteResponseGeometryCache(stale, data);
      }
      Mpi::Barrier(Mpi::World());
      test::GeometryCacheEnvGuard cache_env(stale.string(), false);
      CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, laplace),
                        ContainsSubstring("conductor probe snaps") &&
                            ContainsSubstring("gap-junction-chord-snapped"));
    }
  }

  SECTION("a patch placed without its probe snaps (ProbeSnapsComputed false) is refused "
          "at solve on every rank alike (decision 632)")
  {
    // The placement records ProbeSnapsComputed = false when the model's trace data was not
    // readable (a discovery placeholder: the geometry-only preflight places it and never
    // solves it, test "placement without readable trace data" below); a SOLVE on such a
    // cache fails closed by the model's name, collectively (the 617 form), never reading
    // the gate on un-snapped probes.
    const auto files = WriteCoupon(temp.temp_dir, CouponVariant::CHORD_OFFSET_SNAPPED);
    auto data = ReadResponseGeometryCache(files.cache,
                                          *iodata.solver.electrostatic.response_correction);
    REQUIRE(data.patches.size() == 1);
    REQUIRE(data.patches[0].provenance.probe_snaps_computed);
    data.patches[0].provenance.probe_snaps.clear();
    data.patches[0].provenance.probe_snaps_computed = false;
    const auto uncomputed = temp.temp_dir / "conductor-consistency-chord-uncomputed.json";
    if (Mpi::Root(Mpi::World()))
    {
      WriteResponseGeometryCache(uncomputed, data);
    }
    Mpi::Barrier(Mpi::World());
    test::GeometryCacheEnvGuard cache_env(uncomputed.string(), false);
    CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, laplace),
                      ContainsSubstring("carries no conductor probe snaps") &&
                          ContainsSubstring("gap-junction-chord-snapped") &&
                          ContainsSubstring("decision 632"));
  }

  SECTION("a spatial coupon without any conductor vertex is recorded untestable")
  {
    const auto files = WriteCoupon(temp.temp_dir, CouponVariant::NO_CONDUCTOR);
    test::GeometryCacheEnvGuard cache_env(files.cache.string(), false);
    SurfaceResponseOperator response(iodata, laplace);
    REQUIRE(response.GetPatchCount() == 1);
    CHECK_FALSE(response.HasConductorConsistencyProbes());
    const auto records = response.ApplyConductorConsistencyGate(V[0], 1);
    CHECK(records.empty());
    const auto statistics = response.GetStatistics();
    const auto &diagnostics = statistics["Diagnostics"]["ConductorConsistency"];
    CHECK(diagnostics["Count"].get<int>() == 0);
    CHECK(diagnostics["TestedPatches"].get<int>() == 0);
    CHECK_THAT(diagnostics["Tolerance"].get<double>(), WithinAbs(tolerance, 0.0));
    REQUIRE(diagnostics["UnprobedModels"].size() == 1);
    const auto &unprobed = diagnostics["UnprobedModels"][0];
    CHECK(unprobed["Model"] == "gap-junction-no-conductor");
    CHECK(unprobed["ModelIndex"].get<int>() == 1);
    CHECK_THAT(unprobed["Reason"].get<std::string>(),
               ContainsSubstring("no conductor vertex"));
    // Not gated: the coupon stays applied.
    CHECK(response.GetElectrostaticResponse(V[0]).domain_correction != 0.0);
    CHECK(statistics["ModelCatalog"][0]["PatchWeight"].get<double>() == 1.0);
  }

  SECTION("a conductor without a vertex on the process plane fails closed")
  {
    const auto files = WriteCoupon(temp.temp_dir, CouponVariant::CONDUCTOR_OFF_PLANE);
    test::GeometryCacheEnvGuard cache_env(files.cache.string(), false);
    CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, laplace),
                      ContainsSubstring("Conductor-consistency gate") &&
                          ContainsSubstring("conductor 2") &&
                          ContainsSubstring("gap-junction-off-plane") &&
                          ContainsSubstring("no vertex on the process plane"));
  }
#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator conductor-consistency placement without readable "
                 "trace data",
                 "[surfaceresponseoperator][conductorconsistency][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // Decision 632: a discovery PLACEHOLDER - a Signature-keyed spatial model whose
  // BasisPoints name a file that is never written (discover_surface_response_requirements
  // .placeholder_model) - is placed by the geometry-only preflight exactly as a real model
  // is; the probe snaps are simply not computed (ProbeSnapsComputed false), so the
  // discovery closure proceeds, and a SOLVE refuses such a patch (the gate case above). The
  // form of the high-order spatial round-trip case: the island's cluster feature keyed by
  // its canonical Signature; the same model with the real basis-points file is the control
  // - the two preflights write the SAME patches CSV (the placement reads no basis; the
  // snaps are no CSV column).
  json config = HighOrderSpatialConfig();
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  Mesh island_mesh(MakeIslandMesh(false, false, false, true));
  auto Preflight = [&](const fs::path &library, const std::string &tag)
  {
    auto preflight_config = config;
    preflight_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        library.string();
    IoData iodata(preflight_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    iodata.boundaries.cracked_attributes.insert(10);
    const auto dir = temp.temp_dir / ("preflight-" + tag);
    if (Mpi::Root(Mpi::World()))
    {
      fs::create_directories(dir);
    }
    Mpi::Barrier(Mpi::World());
    WriteSurfaceResponseRequirements(iodata, island_mesh,
                                     (dir / "surface-response-requirements.json").string());
    Mpi::Barrier(Mpi::World());
    std::ifstream manifest_input(dir / "surface-response-requirements.json");
    REQUIRE(manifest_input);
    std::ifstream csv(dir / "surface-response-patches.csv");
    REQUIRE(csv);
    return std::pair<json, std::string>(json::parse(manifest_input),
                                        std::string((std::istreambuf_iterator<char>(csv)),
                                                    std::istreambuf_iterator<char>()));
  };
  // 1. The identification: the island's cluster feature and its Signature (Missing against
  //    the fixture library, whose cluster model is Edges-keyed for another geometry).
  const auto [identified, identified_csv] =
      Preflight(spatial_cluster_library_3d_path, "id");
  const auto &features = identified["Identification"]["Features"];
  const auto cluster =
      std::find_if(features.begin(), features.end(), [](const auto &feature)
                   { return feature["Type"] == "SpatialEdgeCluster"; });
  REQUIRE(cluster != features.end());
  REQUIRE((*cluster)["Match"]["Status"] == "Missing");
  // 2. The Signature-keyed model: the placeholder (BasisPoints absent) and its control (the
  //    fixture's real basis-points file).
  const auto absent_basis = temp.temp_dir / "placeholder-absent-basis-points.csv";
  REQUIRE_FALSE(fs::exists(absent_basis));
  const auto placeholder_library = temp.temp_dir / "spatial-cluster-placeholder-3d.json";
  const auto control_library = temp.temp_dir / "spatial-cluster-signature-control-3d.json";
  if (Mpi::Root(Mpi::World()))
  {
    std::ifstream input(spatial_cluster_library_3d_path);
    REQUIRE(input);
    auto library = json::parse(input);
    auto &model = library["Models"].back();
    REQUIRE(model["Topology"] == "SpatialEdgeCluster");
    model["Name"] = "v2-signature-placeholder";
    model["Signature"] = (*cluster)["Signature"];
    model.erase("Edges");  // keyed and built by the Signature alone, as a placeholder is
    {
      std::ofstream output(control_library);
      output << library.dump(2) << "\n";
    }
    model["BasisPoints"] = absent_basis.string();
    std::ofstream output(placeholder_library);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  auto Placed = [](const json &manifest)
  {
    return std::any_of(manifest["Requirements"].begin(), manifest["Requirements"].end(),
                       [](const auto &requirement)
                       {
                         return requirement["Topology"] == "SpatialEdgeCluster" &&
                                requirement["Status"] == "Exact";
                       });
  };
  const auto [control, control_csv] = Preflight(control_library, "control");
  REQUIRE(Placed(control));
  REQUIRE(control_csv.find("v2-signature-placeholder") != std::string::npos);
  // 3. The placeholder preflight neither throws nor changes a byte of the placement.
  std::pair<json, std::string> placeholder;
  REQUIRE_NOTHROW(placeholder = Preflight(placeholder_library, "placeholder"));
  CHECK(Placed(placeholder.first));
  CHECK(placeholder.second == control_csv);
#endif
}

}  // namespace palace
