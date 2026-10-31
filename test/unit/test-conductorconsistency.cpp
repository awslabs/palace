// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <fmt/format.h>
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
  std::string tag;  // the variant's name in the file names and the model name
  fs::path basis_points, trace_vertices, trace_triangles, fabricated, thin, cache;
  fs::path fabricated_surface, thin_surface;  // the within-R pair of coupon interface 1
  int contour_size = 0;
  std::vector<std::array<double, 3>> perimeter;  // (u, v) of the ten ring points
  std::vector<int> offset_vertices;              // the chord-offset knots (0-based)
};

constexpr double kHalfU = 0.2, kHalfV = 0.25, kHalfW = 0.2, kMetalEdgeV = 0.125;
constexpr double kChordOffset = 0.02;

// The coupon model of a variant's files at library index `idx` (the response-geometry
// cache replaces the library matching).
config::ElectrostaticSolverData::ResponseCorrectionModelData
CouponModel(const CouponFiles &files, CouponVariant variant, int idx)
{
  config::ElectrostaticSolverData::ResponseCorrectionModelData model;
  model.idx = idx;
  model.name = "gap-junction-" + files.tag;
  model.topology = "SpatialEdgeCluster";
  model.fabricated_matrix = files.fabricated.string();
  model.thin_matrix = files.thin.string();
  model.fabricated_surface_matrix = files.fabricated_surface.string();
  model.thin_surface_matrix = files.thin_surface.string();
  model.interfaces = {{4, 1}};
  model.basis_points = files.basis_points.string();
  model.trace_vertices = files.trace_vertices.string();
  model.trace_triangles = files.trace_triangles.string();
  model.spatial_basis = true;
  model.interior_trace_count = 2;
  model.conductor_state_count = variant == CouponVariant::NO_CONDUCTOR ? 0 : 1;
  return model;
}

// The coupon of model `model` placed across the gap at height z (origin (1.0, 0.5, z); u =
// +z, v = +x, w = +y). The claimed portions: 0.4 along each island's facing edge (z in [z
// - 0.2, z + 0.2]) as the identification writes them, one portion per MESH segment (h =
// 0.125: the pieces [-0.2, -0.125], [-0.125, 0], [0, 0.125], [0.125, 0.2] about z;
// segments first_segment.. on island 1's edge x = 0.875, first_segment + 4.. on island 2's
// x = 1.125) - the postoperator serves a portion on the one perimeter mesh segment nearest
// to its midpoint.
config::ElectrostaticSolverData::ResponseCorrectionPatchData
CouponPatch(CouponVariant variant, int model, double z, int feature, int first_segment)
{
  config::ElectrostaticSolverData::ResponseCorrectionPatchData patch;
  patch.model = model;
  patch.origin = {1.0, 0.5, z};
  patch.axis_u = {0.0, 0.0, 1.0};
  patch.axis_v = {1.0, 0.0, 0.0};
  patch.axis_w = {0.0, 1.0, 0.0};
  patch.conductor_references = {{0.0, -0.2, 0.0}, {0.0, 0.2, 0.0}};
  if (variant == CouponVariant::NO_CONDUCTOR)
  {
    patch.conductor_references.pop_back();
  }
  patch.weight = 1.0;
  patch.provenance.feature = feature;
  for (int edge = 0; edge < 2; edge++)
  {
    const double x = edge == 0 ? 0.875 : 1.125;
    int segment = first_segment + 4 * edge;
    for (const auto &[u0, u1] : {std::make_pair(-0.2, -0.125), std::make_pair(-0.125, 0.0),
                                 std::make_pair(0.0, 0.125), std::make_pair(0.125, 0.2)})
    {
      patch.provenance.claims.push_back({segment++, {x, 0.5, z + u0}, {x, 0.5, z + u1}});
    }
  }
  return patch;
}

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
  files.tag = tag;
  files.basis_points = dir / ("conductor-consistency-" + tag + "-basis-points.csv");
  files.trace_vertices = dir / ("conductor-consistency-" + tag + "-trace-vertices.csv");
  files.trace_triangles = dir / ("conductor-consistency-" + tag + "-trace-triangles.csv");
  files.fabricated = dir / ("conductor-consistency-" + tag + "-fabricated.csv");
  files.thin = dir / ("conductor-consistency-" + tag + "-thin.csv");
  files.fabricated_surface =
      dir / ("conductor-consistency-" + tag + "-fabricated-surface.csv");
  files.thin_surface = dir / ("conductor-consistency-" + tag + "-thin-surface.csv");
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
    // fabricated domain energy 2x the thin one, and the within-R surface energy of coupon
    // interface 1 (R = kHalfU, the matching radius; the driver case maps it onto target
    // interface 4) likewise 2x, so an APPLIED coupon adds a positive model energy.
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
    for (const auto &[path, scale] : {std::make_pair(files.fabricated_surface, 1.0e-12),
                                      std::make_pair(files.thin_surface, 0.5e-12)})
    {
      std::ofstream output(path);
      output << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
      for (int i = 1; i <= basis_size; i++)
      {
        for (int j = i; j <= basis_size; j++)
        {
          output << "1,1," << kHalfU << "," << i << "," << j << ","
                 << (i == j ? scale : 0.0) << "," << (i == j ? scale : 0.0) << "\n";
        }
      }
    }
    // The placed geometry (the response-geometry cache replaces the library matching): the
    // coupon at the gap centre z = 1.0, feature 7, claims on segments 3-10.
    config::ElectrostaticSolverData::ResponseCorrectionData data;
    data.matching_radius = kHalfU;
    data.models.push_back(CouponModel(files, variant, 1));
    auto &patch = data.patches.emplace_back(CouponPatch(variant, 1, 1.0, 7, 3));
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

// The two-coupon cache of the swap-detecting driver run (decision 773, review d769
// MINOR-3): the INCONSISTENT coupon at z = 1.0 (patch 0, model 1, feature 7, segments 3-10:
// excluded by the gate, its claim served) and the CONSISTENT coupon at z = 2.0 (patch 1,
// model 2, feature 8, segments 11-18: applied, so the self-consistent solve differs from
// the raw one), 3R apart on the gap of the long islands (MakeSerialLongGapIslandMesh).
fs::path WriteTwoCouponCache(const fs::path &dir, const CouponFiles &inconsistent,
                             const CouponFiles &consistent)
{
  const fs::path cache = dir / "conductor-consistency-two-coupon-geometry.json";
  if (Mpi::Root(Mpi::World()))
  {
    config::ElectrostaticSolverData::ResponseCorrectionData data;
    data.matching_radius = kHalfU;
    data.models.push_back(CouponModel(inconsistent, CouponVariant::INCONSISTENT, 1));
    data.models.push_back(CouponModel(consistent, CouponVariant::CONSISTENT, 2));
    data.patches.push_back(CouponPatch(CouponVariant::INCONSISTENT, 1, 1.0, 7, 3));
    data.patches.push_back(CouponPatch(CouponVariant::CONSISTENT, 2, 2.0, 8, 11));
    WriteResponseGeometryCache(cache, data);
  }
  Mpi::Barrier(Mpi::World());
  return cache;
}

// The island pair of MakeSerialIslandMesh(neighboring_island) with both islands extended
// along the gap: a 2 x 1 x 3 box of 16 x 4 x 24 hexahedra (h = 0.125 in the plane), island
// 1 (attribute 9) x in [0.375, 0.875] and island 2 (attribute 10) x in [1.125, 1.625], both
// z in [0.75, 2.25] on the plane y = 0.5, so the gap carries the two coupons of
// WriteTwoCouponCache (z = 1.0 and 2.0) without overlap.
mfem::Mesh MakeSerialLongGapIslandMesh()
{
  constexpr double tolerance = 1.0e-12;
  mfem::Mesh serial =
      mfem::Mesh::MakeCartesian3D(16, 4, 24, mfem::Element::HEXAHEDRON, 2.0, 1.0, 3.0);
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
    double xmin = 2.0, xmax = 0.0, zmin = 3.0, zmax = 0.0;
    for (const int vertex : vertices)
    {
      const double *point = serial.GetVertex(vertex);
      on_plane = on_plane && std::abs(point[1] - 0.5) < tolerance;
      xmin = std::min(xmin, point[0]);
      xmax = std::max(xmax, point[0]);
      zmin = std::min(zmin, point[2]);
      zmax = std::max(zmax, point[2]);
    }
    if (!on_plane || zmin < 0.75 - tolerance || zmax > 2.25 + tolerance)
    {
      continue;
    }
    const bool island_1 = xmin >= 0.375 - tolerance && xmax <= 0.875 + tolerance;
    const bool island_2 = xmin >= 1.125 - tolerance && xmax <= 1.625 + tolerance;
    if (!island_1 && !island_2)
    {
      continue;
    }
    serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
    serial.SetBdrAttribute(serial.GetNBE() - 1, island_1 ? 9 : 10);
  }
  serial.FinalizeTopology();
  serial.Finalize();
  return serial;
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
    // S7diag R1 (decision 758): nothing excluded, nothing served - the record carries no
    // RawPortions entry (a window with 0 exclusions writes main's palace.json).
    CHECK_FALSE(record.contains("RawPortions"));
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
    // S7diag R1 (decision 758): before the gate nothing is served.
    CHECK_FALSE(response.GetStatistics()["Diagnostics"]["ConductorConsistency"].contains(
        "RawPortions"));
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
    // S7diag R1 (decision 758): the excluded patch's claim is served RAW - the record's
    // RawPortions (the DomainBoundary RawPortions shape) hold EXACTLY the coupon's eight
    // claimed portions (0.8 in all), replicated on every rank; the response operator stays
    // at weight 0.
    auto CheckServedClaim = [&](const json &raw)
    {
      REQUIRE(raw["Count"].get<int>() == 8);
      CHECK(raw["GeometricCells"].get<int>() == 0);
      CHECK(raw["DuplicatePatches"].get<int>() == 0);
      CHECK_THAT(raw["Length"].get<double>(), WithinAbs(0.8, 1.0e-12));
      CHECK(raw["ByType"]["SpatialEdgeCluster"]["Portions"].get<int>() == 8);
      CHECK_THAT(raw["ByType"]["SpatialEdgeCluster"]["Length"].get<double>(),
                 WithinAbs(0.8, 1.0e-12));
      REQUIRE(raw["ByFeature"].size() == 1);
      CHECK(raw["ByFeature"][0]["Feature"].get<int>() == 7);
      CHECK(raw["ByFeature"][0]["Portions"].get<int>() == 8);
      const double z_cuts[5] = {0.8, 0.875, 1.0, 1.125, 1.2};
      for (int k = 0; k < 8; k++)
      {
        const auto &portion = raw["Portions"][k];
        const double x = k < 4 ? 0.875 : 1.125;
        const std::array<double, 3> p0 = {x, 0.5, z_cuts[k % 4]};
        const std::array<double, 3> p1 = {x, 0.5, z_cuts[k % 4 + 1]};
        CHECK(portion["Feature"].get<int>() == 7);
        CHECK(portion["Type"] == "SpatialEdgeCluster");
        CHECK(portion["Segment"].get<int>() == 3 + k);
        CHECK(portion["Patches"] == json::array({0}));
        CHECK_THAT(portion["Length"].get<double>(), WithinAbs(p1[2] - p0[2], 1.0e-12));
        for (int d = 0; d < 3; d++)
        {
          CHECK_THAT(portion["P0"][d].get<double>(), WithinAbs(p0[d], 1.0e-12));
          CHECK_THAT(portion["P1"][d].get<double>(), WithinAbs(p1[d], 1.0e-12));
        }
      }
      CHECK_THAT(Replicated(raw["Count"].get<double>()), WithinAbs(8.0, 0.0));
      CHECK_THAT(Replicated(raw["Length"].get<double>()),
                 WithinAbs(raw["Length"].get<double>(), 0.0));
      CHECK_THAT(raw["Rule"].get<std::string>(),
                 ContainsSubstring("surface-response-conductor-consistency-energy.csv"));
    };
    REQUIRE(diagnostics.contains("RawPortions"));
    CheckServedClaim(diagnostics["RawPortions"]);
    // The exclusion is sticky: the next excitation finds no applied probe patch.
    const auto later = response.ApplyConductorConsistencyGate(V[1], 2);
    CHECK(later.empty());
    CHECK(statistics["Diagnostics"]["ConductorConsistency"]["Records"].size() == 1);
    // ... and serves the same two portions (no duplicate appended).
    const auto later_statistics = response.GetStatistics();
    CheckServedClaim(
        later_statistics["Diagnostics"]["ConductorConsistency"]["RawPortions"]);
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

  SECTION("a conductor without a vertex on the process plane is recorded untestable "
          "(H9b: a library property, no abort)")
  {
    // U-H9b (hardening-2 DESIGN row H9b, landed in the identification fix round's lane B):
    // the former abort at this site is a coupon property, not device geometry; the model
    // is listed UnprobedModels by name with the conductor, the gate skips it, the coupon
    // stays applied.
    const auto files = WriteCoupon(temp.temp_dir, CouponVariant::CONDUCTOR_OFF_PLANE);
    test::GeometryCacheEnvGuard cache_env(files.cache.string(), false);
    SurfaceResponseOperator response(iodata, laplace);
    REQUIRE(response.GetPatchCount() == 1);
    CHECK_FALSE(response.HasConductorConsistencyProbes());
    CHECK(response.ApplyConductorConsistencyGate(V[0], 1).empty());
    const auto statistics = response.GetStatistics();
    const auto &diagnostics = statistics["Diagnostics"]["ConductorConsistency"];
    REQUIRE(diagnostics["UnprobedModels"].size() == 1);
    const auto &unprobed = diagnostics["UnprobedModels"][0];
    CHECK(unprobed["Model"] == "gap-junction-off-plane");
    CHECK_THAT(unprobed["Reason"].get<std::string>(),
               ContainsSubstring("conductor 2") &&
                   ContainsSubstring("no vertex on the process plane"));
    CHECK(response.GetElectrostaticResponse(V[0]).domain_correction != 0.0);
    CHECK(statistics["ModelCatalog"][0]["PatchWeight"].get<double>() == 1.0);
  }
#endif
}

// S7diag R1 (decision 758): a patch the conductor-consistency gate excludes is never
// dropped from the corrected interface energies. Its raw claim (here the coupon's claimed
// portions, 0.4 on each island's facing edge) is served RAW from the excluding
// excitation on, through the uncovered-requirement machinery of the postoperator (one
// nearest perimeter segment per quadrature point), as its own ConductorConsistency share
// beside Missing (uncovered) and DomainBoundary: E_ft = E_out + models + uncovered +
// domainboundary + conductorconsistency, the same for E_ff, and E_sc on the corrected
// field. The response operator is unchanged (weight 0 stands). The INCONSISTENT coupon
// (excluded at source 1, sticky) against the CONSISTENT one (applied: 0 exclusions, the
// pre-existing identity E_ft = E_out + models, no ConductorConsistency table, no
// RawPortions record) on the same raw solve; the served term is bracketed by two
// independent decision-366 Region readings of the raw field (the nearest-point window of a
// portion that stops 0.05 short of its edge's corners lies between the within-R rectangle
// over the portion's central [0.95, 1.05], whose points are nearer to the portion than to
// any corner edge, and the within-R rectangle over the whole portion [0.8, 1.2], which
// reaches into the corner edges' nearest-point sectors).
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "Electrostatic conductor-consistency-excluded claims keep their raw energy "
    "(S7diag R1)",
    "[electrostaticsolver][surfaceresponseoperator][conductorconsistency][3d]"
    "[Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = kHalfU;
  const fs::path mesh_path = temp.temp_dir / "conductor-consistency-islands.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    mfem::Mesh serial = MakeSerialIslandMesh(false, false, false, true);
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(Mpi::World());
  const auto inconsistent = WriteCoupon(temp.temp_dir, CouponVariant::INCONSISTENT);
  const auto consistent = WriteCoupon(temp.temp_dir, CouponVariant::CONSISTENT);

  json config = HighOrderSpatialConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  // Mesh units are metres and stay unscaled (Lc = 1 x L0): the placed geometry of the
  // cache and the surface response files' R (m) = 0.2 are in mesh coordinates.
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  // The Region readings of the raw field (decision 366: half-open along, transverse <= R,
  // in the plane y = 0.5): 5 / 6 the central thirds of the two claimed portions, 7 / 8 the
  // whole claimed portions, on the islands' facing edges x = 0.875 / 1.125.
  auto &dielectric = config["Boundaries"]["Postprocessing"]["Dielectric"];
  const json target = dielectric[0];
  int index = 5;
  for (const std::array<double, 6> &segment :
       {std::array<double, 6>{0.875, 0.5, 0.95, 0.875, 0.5, 1.05},
        std::array<double, 6>{1.125, 0.5, 0.95, 1.125, 0.5, 1.05},
        std::array<double, 6>{0.875, 0.5, 0.8, 0.875, 0.5, 1.2},
        std::array<double, 6>{1.125, 0.5, 0.8, 1.125, 0.5, 1.2}})
  {
    json region = target;
    region["Index"] = index++;
    region.erase("AutomaticEdges");
    region.erase("EdgeDistances");
    region.erase("EdgeFrameNormal");
    region["Region"] = {
        {"Segments", {segment}}, {"Distance", R}, {"Normal", {0.0, 1.0, 0.0}}};
    dielectric.push_back(region);
  }
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["TraceCoupling"] = "SurfaceMortar";
  correction["MortarOversampling"] = 2;
  correction["CorrectionMode"] = "Both";
  config["Solver"]["Linear"] = {{"Tol", 1.0e-12}, {"MaxIts", 400}};
  const fs::path excluded_dir = temp.temp_dir / "cc-excluded";
  const fs::path applied_dir = temp.temp_dir / "cc-applied";
  {
    test::GeometryCacheEnvGuard cache_env(inconsistent.cache.string(), false);
    config["Problem"]["Output"] = excluded_dir.string();
    test::RunElectrostatic(config);
  }
  {
    test::GeometryCacheEnvGuard cache_env(consistent.cache.string(), false);
    config["Problem"]["Output"] = applied_dir.string();
    test::RunElectrostatic(config);
  }
  // The swap-detecting run (decision 773, review d769 MINOR-3): the INCONSISTENT coupon
  // beside an APPLIED consistent coupon 3R along the gap of the long islands, so the
  // self-consistent field differs from the raw one and the served term on the corrected
  // field from the served term on the raw one - each table column is read against its own
  // field (the Region bracket of the raw field holds on the fixed-trace row alone, the
  // bracket of the corrected field on the self-consistent row alone).
  const fs::path long_mesh_path = temp.temp_dir / "conductor-consistency-long-islands.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    mfem::Mesh serial = MakeSerialLongGapIslandMesh();
    std::ofstream output(long_mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(Mpi::World());
  const fs::path two_coupon_cache =
      WriteTwoCouponCache(temp.temp_dir, inconsistent, consistent);
  const fs::path two_coupon_dir = temp.temp_dir / "cc-two-coupon";
  {
    test::GeometryCacheEnvGuard cache_env(two_coupon_cache.string(), false);
    json long_config = config;
    long_config["Model"]["Mesh"] = long_mesh_path.string();
    long_config["Problem"]["Output"] = two_coupon_dir.string();
    test::RunElectrostatic(long_config);
  }
  Mpi::Barrier(Mpi::World());
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }

  // The row of source `source` (the terminal index column) of a per-source table.
  auto SourceRow = [&](const Table &table, const std::string &source_header, int source)
  {
    const Column &column = test::ColumnByHeader(table, source_header);
    for (std::size_t i = 0; i < column.data.size(); i++)
    {
      if (static_cast<int>(std::lround(column.data[i])) == source)
      {
        return i;
      }
    }
    FAIL("source " << source << " has no row");
    return std::size_t(0);
  };
  auto ModelEnergySum = [&](const fs::path &dir, int source, int evaluation)
  {
    const Table models = test::LoadCsv(dir / "surface-response-model-energy.csv");
    const Column &sources = test::ColumnByHeader(models, "source");
    const Column &eval = test::ColumnByHeader(models, "evaluation");
    const Column &energy = test::ColumnByHeader(models, "fabricated surface energy[4] (J)");
    double sum = 0.0;
    for (std::size_t i = 0; i < eval.data.size(); i++)
    {
      if (static_cast<int>(std::lround(sources.data[i])) == source &&
          static_cast<int>(std::lround(eval.data[i])) == evaluation)
      {
        sum += energy.data[i];
      }
    }
    return sum;
  };
  auto OutsideEnergy = [&](const fs::path &dir, int source)
  {
    const Table edge = test::LoadCsv(dir / "surface-Q-edge.csv");
    REQUIRE(edge.n_rows() == 2);
    return test::ColumnByHeader(edge, "E_out (J)").data[SourceRow(edge, "i", source)];
  };
  auto RawEnergy = [&](const fs::path &dir, int source, int interface_index)
  {
    const Table surface = test::LoadCsv(dir / "surface-Q.csv");
    const Table domain = test::LoadCsv(dir / "domain-E.csv");
    REQUIRE(surface.n_rows() == 2);
    REQUIRE(domain.n_rows() == 2);
    const Column &participation =
        test::ColumnByHeader(surface, fmt::format("p_surf[{}]", interface_index));
    const Column &energy = test::ColumnByHeader(domain, "E_elec (J)");
    return participation.data[SourceRow(surface, "i", source)] *
           energy.data[SourceRow(domain, "i", source)];
  };
  auto CorrectedEnergies = [&](const fs::path &dir, int source)
  {
    const Table corrected = test::LoadCsv(dir / "surface-Q-corrected.csv");
    REQUIRE(corrected.n_rows() == 2);
    const std::size_t row = SourceRow(corrected, "i", source);
    return std::array<double, 3>{
        test::ColumnByHeader(corrected, "E_surf postprocessed fixed-trace[4] (J)")
            .data[row],
        test::ColumnByHeader(corrected, "E_surf postprocessed fixed-flux[4] (J)").data[row],
        test::ColumnByHeader(corrected, "E_surf corrected[4] (J)").data[row]};
  };
  // The per-type energy table of a run keyed by (source, evaluation, type).
  struct TypeTable
  {
    std::map<std::tuple<int, int, std::string>, std::vector<std::string>> rows;
    std::size_t portions = 0, length = 0, energy = 0, share = 0, extra = 0;
  };
  auto ReadTypeTable =
      [&](const fs::path &path, const std::string &label, const std::string &extra)
  {
    const test::TextCsv csv = test::ReadTextCsv(path);
    TypeTable table;
    const std::size_t source = csv.Column("source"), type = csv.Column("type"),
                      evaluation = csv.Column("evaluation");
    table.portions = csv.Column("portions");
    table.length = csv.Column("length (m)");
    table.extra = csv.Column(extra);
    table.energy = csv.Column(label + " raw energy[4] (J)");
    table.share = csv.Column(label + " share[4]");
    for (const auto &row : csv.rows)
    {
      table.rows[{std::stoi(row[source]), std::stoi(row[evaluation]), row[type]}] = row;
    }
    return table;
  };
  auto Metadata = [&](const fs::path &dir)
  {
    std::ifstream input(dir / "palace.json");
    REQUIRE(input);
    return json::parse(input);
  };

  // The raw solve knows nothing of the coupon: the raw outputs of the two runs are
  // byte-identical.
  for (const char *file : {"terminal-C.csv", "terminal-V.csv", "domain-E.csv",
                           "surface-Q.csv", "surface-Q-edge.csv"})
  {
    INFO(file);
    CHECK(test::ReadFile(excluded_dir / file) == test::ReadFile(applied_dir / file));
  }

  // The CONSISTENT coupon: applied (0 exclusions), the pre-existing identity E_ft = E_out +
  // models, nothing served, no table, no RawPortions record.
  CHECK_FALSE(
      fs::exists(applied_dir / "surface-response-conductor-consistency-energy.csv"));
  CHECK_FALSE(fs::exists(applied_dir / "surface-response-uncovered-energy.csv"));
  CHECK_FALSE(fs::exists(applied_dir / "surface-response-domain-boundary-energy.csv"));
  {
    const auto metadata = Metadata(applied_dir);
    const auto &record =
        metadata.at("SurfaceResponse").at("Diagnostics").at("ConductorConsistency");
    CHECK(record.at("Count").get<int>() == 0);
    CHECK(record.at("TestedPatches").get<int>() == 2);
    CHECK_FALSE(record.contains("RawPortions"));
    for (int source = 1; source <= 2; source++)
    {
      INFO("source " << source);
      const auto energies = CorrectedEnergies(applied_dir, source);
      const double outside = OutsideEnergy(applied_dir, source);
      const double models_ft = ModelEnergySum(applied_dir, source, 0);
      CHECK(models_ft > 0.0);
      CHECK_THAT(energies[0], WithinRel(outside + models_ft, 1.0e-10));
      CHECK_THAT(energies[1],
                 WithinRel(outside + ModelEnergySum(applied_dir, source, 1), 1.0e-10));
    }
    // The 0-exclusion fingerprint of this run's corrected outputs (FNV-1a 64 of the
    // file bytes; palace.json's SurfaceResponse record without its WallTime entries and
    // the GeometryCache path of this run's temp directory), reported for the byte-identity
    // reading of a build of this test against main's sources (palace-unit-tests -s).
    auto Fingerprint = [](const std::string &bytes)
    {
      std::uint64_t hash = 14695981039346656037ull;
      for (const unsigned char c : bytes)
      {
        hash = (hash ^ c) * 1099511628211ull;
      }
      return hash;
    };
    std::function<void(json &)> EraseVolatile = [&](json &node)
    {
      if (node.is_object())
      {
        node.erase("WallTime");
        node.erase("GeometryCache");
        for (auto &[key, value] : node.items())
        {
          EraseVolatile(value);
        }
      }
      else if (node.is_array())
      {
        for (auto &value : node)
        {
          EraseVolatile(value);
        }
      }
    };
    json stripped = metadata.at("SurfaceResponse");
    EraseVolatile(stripped);
    SUCCEED(fmt::format(
        "S7diag R1 0-exclusion fingerprints: surface-Q-corrected.csv {:016x} "
        "surface-response-model-energy.csv {:016x} "
        "surface-response-patches.csv {:016x} palace.json SurfaceResponse(-WallTime,"
        "-GeometryCache) {:016x}",
        Fingerprint(test::ReadFile(applied_dir / "surface-Q-corrected.csv")),
        Fingerprint(test::ReadFile(applied_dir / "surface-response-model-energy.csv")),
        Fingerprint(test::ReadFile(applied_dir / "surface-response-patches.csv")),
        Fingerprint(stripped.dump())));
  }

  // The INCONSISTENT coupon: excluded at source 1 (sticky), its eight claimed portions
  // served as the ConductorConsistency share of both sources.
  REQUIRE(fs::is_regular_file(excluded_dir /
                              "surface-response-conductor-consistency-energy.csv"));
  CHECK_FALSE(fs::exists(excluded_dir / "surface-response-uncovered-energy.csv"));
  CHECK_FALSE(fs::exists(excluded_dir / "surface-response-domain-boundary-energy.csv"));
  {
    const auto metadata = Metadata(excluded_dir);
    const auto &record =
        metadata.at("SurfaceResponse").at("Diagnostics").at("ConductorConsistency");
    CHECK(record.at("Count").get<int>() == 1);
    CHECK(record.at("TestedPatches").get<int>() == 1);
    REQUIRE(record.contains("RawPortions"));
    const auto &raw = record.at("RawPortions");
    CHECK(raw.at("Count").get<int>() == 8);
    CHECK(raw.at("GeometricCells").get<int>() == 0);
    CHECK(raw.at("DuplicatePatches").get<int>() == 0);
    CHECK_THAT(raw.at("Length").get<double>(), WithinAbs(0.8, 1.0e-12));
    CHECK(raw.at("ByType").at("SpatialEdgeCluster").at("Portions").get<int>() == 8);
    const TypeTable cc =
        ReadTypeTable(excluded_dir / "surface-response-conductor-consistency-energy.csv",
                      "conductorconsistency", "geometric cells");
    // 2 sources x 3 evaluations x (SpatialEdgeCluster + Total).
    REQUIRE(cc.rows.size() == 12);
    for (int source = 1; source <= 2; source++)
    {
      INFO("source " << source);
      REQUIRE(cc.rows.count({source, 0, "SpatialEdgeCluster"}) == 1);
      const auto &row = cc.rows.at({source, 0, "SpatialEdgeCluster"});
      CHECK(std::stoi(row[cc.portions]) == 8);
      CHECK(std::stoi(row[cc.extra]) == 0);
      CHECK_THAT(std::stod(row[cc.length]), WithinAbs(0.8, 1.0e-9));
      const double cc_ft = std::stod(cc.rows.at({source, 0, "Total"})[cc.energy]);
      const double cc_ff = std::stod(cc.rows.at({source, 1, "Total"})[cc.energy]);
      const double cc_sc = std::stod(cc.rows.at({source, 2, "Total"})[cc.energy]);
      CHECK(cc_ft > 0.0);
      CHECK(std::stod(row[cc.energy]) == cc_ft);
      CHECK(cc_ff == cc_ft);  // the same raw field
      const auto energies = CorrectedEnergies(excluded_dir, source);
      const double outside = OutsideEnergy(excluded_dir, source);
      // The only patch is excluded before any energy is kept: models = 0.
      CHECK(ModelEnergySum(excluded_dir, source, 0) == 0.0);
      CHECK(ModelEnergySum(excluded_dir, source, 1) == 0.0);
      CHECK_THAT(energies[0], WithinRel(outside + cc_ft, 1.0e-10));
      CHECK_THAT(energies[1], WithinRel(outside + cc_ft, 1.0e-10));
      CHECK_THAT(std::stod(cc.rows.at({source, 0, "Total"})[cc.share]),
                 WithinRel(cc_ft / energies[0], 1.0e-9));
      // The self-consistent solve of a weight-0 operator is the raw solve: E_sc = E_ft and
      // the served term is the same, to the solver tolerance.
      REQUIRE(std::isfinite(energies[2]));
      CHECK_THAT(energies[2], WithinRel(energies[0], 1.0e-8));
      CHECK_THAT(cc_sc, WithinRel(cc_ft, 1.0e-8));
      CHECK_THAT(std::stod(cc.rows.at({source, 2, "Total"})[cc.share]),
                 WithinRel(cc_sc / energies[2], 1.0e-9));
      // The independent bracket of the served term by the Region readings of the raw
      // field: the central thirds' rectangles lie wholly inside the portions' nearest-point
      // windows, the whole portions' rectangles contain them.
      const double central =
          RawEnergy(excluded_dir, source, 5) + RawEnergy(excluded_dir, source, 6);
      const double whole =
          RawEnergy(excluded_dir, source, 7) + RawEnergy(excluded_dir, source, 8);
      INFO(fmt::format("E_out {:.6e} within-R {:.6e} served {:.6e} central {:.6e} whole "
                       "{:.6e} E_ft {:.6e}",
                       outside, RawEnergy(excluded_dir, source, 4) - outside, cc_ft,
                       central, whole, energies[0]));
      CHECK(central > 0.0);
      CHECK(cc_ft > central);
      CHECK(cc_ft < whole);
      // ... and by the within-R energy of the whole interface, which the served term
      // cannot exceed.
      CHECK(cc_ft < RawEnergy(excluded_dir, source, 4) - outside);
    }
  }

  // The two-coupon run (decision 773): the inconsistent coupon (patch 0) excluded at source
  // 1 and served, the consistent coupon (patch 1) applied at both sources.
  REQUIRE(fs::is_regular_file(two_coupon_dir /
                              "surface-response-conductor-consistency-energy.csv"));
  {
    const auto metadata = Metadata(two_coupon_dir);
    const auto &record =
        metadata.at("SurfaceResponse").at("Diagnostics").at("ConductorConsistency");
    CHECK(record.at("Count").get<int>() == 1);
    REQUIRE(record.at("ExcludedPatches").size() == 1);
    CHECK(record.at("ExcludedPatches")[0].at("Patch").get<int>() == 0);
    // Both coupons probed at source 1, the applied one again at source 2.
    CHECK(record.at("TestedPatches").get<int>() == 3);
    REQUIRE(record.contains("RawPortions"));
    const auto &raw = record.at("RawPortions");
    CHECK(raw.at("Count").get<int>() == 8);
    CHECK_THAT(raw.at("Length").get<double>(), WithinAbs(0.8, 1.0e-12));
    REQUIRE(raw.at("ByFeature").size() == 1);
    CHECK(raw.at("ByFeature")[0].at("Feature").get<int>() == 7);
    const TypeTable cc =
        ReadTypeTable(two_coupon_dir / "surface-response-conductor-consistency-energy.csv",
                      "conductorconsistency", "geometric cells");
    REQUIRE(cc.rows.size() == 12);
    // The Region energies of a field: the regions are not response targets, so their
    // columns of surface-Q-corrected.csv carry the plain interface energy of each field.
    const Table corrected = test::LoadCsv(two_coupon_dir / "surface-Q-corrected.csv");
    REQUIRE(corrected.n_rows() == 2);
    auto RegionEnergy = [&](const char *field, int source, int region)
    {
      return test::ColumnByHeader(corrected,
                                  fmt::format("E_surf {}[{}] (J)", field, region))
          .data[SourceRow(corrected, "i", source)];
    };
    for (int source = 1; source <= 2; source++)
    {
      INFO("source " << source);
      CHECK(std::stoi(cc.rows.at({source, 0, "SpatialEdgeCluster"})[cc.portions]) == 8);
      const double cc_ft = std::stod(cc.rows.at({source, 0, "Total"})[cc.energy]);
      const double cc_ff = std::stod(cc.rows.at({source, 1, "Total"})[cc.energy]);
      const double cc_sc = std::stod(cc.rows.at({source, 2, "Total"})[cc.energy]);
      const auto energies = CorrectedEnergies(two_coupon_dir, source);
      const double outside = OutsideEnergy(two_coupon_dir, source);
      const double models_ft = ModelEnergySum(two_coupon_dir, source, 0);
      const double models_ff = ModelEnergySum(two_coupon_dir, source, 1);
      const double models_sc = ModelEnergySum(two_coupon_dir, source, 2);
      REQUIRE(std::isfinite(energies[2]));
      INFO(fmt::format(
          "E_ft {:.6e} E_sc {:.6e} models ft {:.6e} sc {:.6e} served ft "
          "{:.6e} sc {:.6e} central ft {:.6e} sc {:.6e} whole ft {:.6e} sc "
          "{:.6e}",
          energies[0], energies[2], models_ft, models_sc, cc_ft, cc_sc,
          RegionEnergy("raw", source, 5) + RegionEnergy("raw", source, 6),
          RegionEnergy("corrected", source, 5) + RegionEnergy("corrected", source, 6),
          RegionEnergy("raw", source, 7) + RegionEnergy("raw", source, 8),
          RegionEnergy("corrected", source, 7) + RegionEnergy("corrected", source, 8)));
      // The consistent coupon is applied: the models' energy is kept in every evaluation.
      CHECK(models_ft > 0.0);
      CHECK(models_sc > 0.0);
      CHECK(cc_ft > 0.0);
      CHECK(cc_ff == cc_ft);
      // The fixed-trace / fixed-flux identities (the fixture's model energies dwarf its
      // field energies, so these are blind to the served term here: the swap detectors are
      // the own-field brackets below).
      CHECK_THAT(energies[0], WithinRel(outside + models_ft + cc_ft, 1.0e-10));
      CHECK_THAT(energies[1], WithinRel(outside + models_ff + cc_ff, 1.0e-10));
      // The applied coupon's self-consistent solve differs from the raw one, and the served
      // term on the corrected field from the served term on the raw one: the two maps are
      // distinguishable.
      CHECK(std::abs(energies[2] - energies[0]) > 1.0e-6 * energies[0]);
      CHECK(std::abs(cc_sc - cc_ft) > 1.0e-6 * std::max(cc_sc, cc_ft));
      // Each map against its own field: the Region bracket of the A1 form read on the
      // corrected field for the self-consistent row and on the raw field for the
      // fixed-trace row ...
      const double central_ft =
          RegionEnergy("raw", source, 5) + RegionEnergy("raw", source, 6);
      const double whole_ft =
          RegionEnergy("raw", source, 7) + RegionEnergy("raw", source, 8);
      const double central_sc =
          RegionEnergy("corrected", source, 5) + RegionEnergy("corrected", source, 6);
      const double whole_sc =
          RegionEnergy("corrected", source, 7) + RegionEnergy("corrected", source, 8);
      CHECK(central_ft > 0.0);
      CHECK(cc_ft > central_ft);
      CHECK(cc_ft < whole_ft);
      CHECK(central_sc > 0.0);
      CHECK(cc_sc > central_sc);
      CHECK(cc_sc < whole_sc);
      // ... and neither value lies in the other field's bracket: a swap of the two maps
      // fails the brackets above.
      CHECK_FALSE((cc_ft > central_sc && cc_ft < whole_sc));
      CHECK_FALSE((cc_sc > central_ft && cc_sc < whole_ft));
      // Each share against its own evaluation's interface energy.
      CHECK_THAT(std::stod(cc.rows.at({source, 0, "Total"})[cc.share]),
                 WithinRel(cc_ft / energies[0], 1.0e-9));
      CHECK_THAT(std::stod(cc.rows.at({source, 2, "Total"})[cc.share]),
                 WithinRel(cc_sc / energies[2], 1.0e-9));
    }
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
