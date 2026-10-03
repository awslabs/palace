// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
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
// closed).
enum class CouponVariant
{
  CONSISTENT,
  INCONSISTENT,
  CONDUCTOR_OFF_PLANE
};

struct CouponFiles
{
  fs::path basis_points, trace_vertices, trace_triangles, fabricated, thin, cache;
  int contour_size = 0;
  std::vector<std::array<double, 3>> perimeter;  // (u, v) of the ten ring points
};

constexpr double kHalfU = 0.2, kHalfV = 0.25, kHalfW = 0.2, kMetalEdgeV = 0.125;

CouponFiles WriteCoupon(const fs::path &dir, CouponVariant variant)
{
  const std::string tag = variant == CouponVariant::CONSISTENT     ? "consistent"
                          : variant == CouponVariant::INCONSISTENT ? "inconsistent"
                                                                   : "off-plane";
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
    if (v <= -kMetalEdgeV + 1.0e-12)
    {
      plane[i] = AddConductor({u, v, 0.0}, 1);
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
    const int basis_size = files.contour_size + 1;
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
    model.conductor_state_count = 1;
    auto &patch = data.patches.emplace_back();
    patch.model = 1;
    patch.origin = {1.0, 0.5, 1.0};
    patch.axis_u = {0.0, 0.0, 1.0};
    patch.axis_v = {1.0, 0.0, 0.0};
    patch.axis_w = {0.0, 1.0, 0.0};
    patch.conductor_references = {{0.0, -0.2, 0.0}, {0.0, 0.2, 0.0}};
    patch.weight = 1.0;
    patch.provenance.feature = 7;
    // Two claimed portions of 0.4 along the islands' facing edges (z in [0.8, 1.2]).
    patch.provenance.claims = {{3, {0.875, 0.5, 0.8}, {0.875, 0.5, 1.2}},
                               {5, {1.125, 0.5, 0.8}, {1.125, 0.5, 1.2}}};
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

}  // namespace palace
