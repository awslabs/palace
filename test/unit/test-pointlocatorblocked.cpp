// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"

using namespace palace;

namespace
{

// The unit cube in 6 x 6 x 6 x 6 tetrahedra; at two or more ranks rank 0 owns every fifth
// element, spread over the whole cube, so that its routing boxes cover (nearly) every
// point: the hot-rank form of the candidate exchange (decision 1501).
std::unique_ptr<mfem::ParMesh> MakeSpreadParMesh(MPI_Comm comm)
{
  mfem::Mesh serial =
      mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0);
  const int size = Mpi::Size(comm);
  std::vector<int> partitioning(serial.GetNE());
  for (int e = 0; e < serial.GetNE(); e++)
  {
    partitioning[e] = (size == 1 || e % 5 == 0) ? 0 : 1 + (e % (size - 1));
  }
  auto mesh = std::make_unique<mfem::ParMesh>(comm, serial, partitioning.data());
  mesh->EnsureNodes();  // the gslib fallback searches the nodal mesh (as Palace's meshes)
  return mesh;
}

// A ragged local list (byNODES): interior points, grid vertices and face points (shared by
// many elements and ranks: the lowest-rank / lowest-element tie-break), points within
// gslib's border tolerance outside the cube (the fallback exchange) and points far
// outside (unresolved).
mfem::Vector MakeQueryPoints(int count, int seed)
{
  std::vector<std::array<double, 3>> points;
  unsigned long long state = 0x9E3779B97F4A7C15ull * static_cast<unsigned>(seed + 1);
  auto Next = [&state]()
  {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    return static_cast<double>(state >> 11) / static_cast<double>(1ull << 53);
  };
  for (int i = 0; i < count; i++)
  {
    std::array<double, 3> p = {Next(), Next(), Next()};
    switch (i % 6)
    {
      case 1:  // a grid vertex
        for (auto &x : p)
        {
          x = static_cast<int>(x * 7.0) / 6.0;
        }
        break;
      case 2:  // on a grid plane
        p[1] = static_cast<int>(p[1] * 7.0) / 6.0;
        break;
      case 3:  // just outside the face x = 1 (within gslib's border tolerance)
        p[0] = 1.0 + 1.0e-7;
        break;
      case 4:
        if (i % 12 == 4)  // far outside
        {
          p[2] = 3.0;
        }
        break;
      default:
        break;
    }
    points.push_back(p);
  }
  mfem::Vector xyz(3 * count);
  for (int i = 0; i < count; i++)
  {
    for (int d = 0; d < 3; d++)
    {
      xyz(d * count + i) = points[i][d];
    }
  }
  return xyz;
}

// Byte equality of two equally long double arrays; std::memcmp is not called on the null
// data() of an empty vector (UBSan: 'null pointer passed as argument 1').
bool BytesEqual(const std::vector<double> &a, const std::vector<double> &b)
{
  return a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

void CheckBitwiseEqual(const DistributedPointLocation &a, const DistributedPointLocation &b)
{
  CHECK(a.owners == b.owners);
  CHECK(a.elements == b.elements);
  REQUIRE(a.references.size() == b.references.size());
  CHECK(BytesEqual(a.references, b.references));
  REQUIRE(a.owner_values.size() == b.owner_values.size());
  CHECK(BytesEqual(a.owner_values, b.owner_values));
  CHECK(a.candidate_queries == b.candidate_queries);
  CHECK(a.fallback_queries == b.fallback_queries);
}

}  // namespace

TEST_CASE("Blocked point-locator exchanges locate every point as the unblocked one",
          "[surfaceresponseoperator][pointlocator][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("The point locator's fallback requires MFEM_USE_GSLIB");
#else
  // Round 3 M1 (decision 1501): the candidate and fallback exchanges in rounds of at most
  // Q local points give every point the same owner, element, reference and owner value
  // as the single exchange, for ragged per-rank lists (incl. an empty one).
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), size = Mpi::Size(comm);
  auto mesh = MakeSpreadParMesh(comm);
  const std::function<double(int)> element_value = [rank](int element)
  { return 0.25 * element + 1000.0 * rank; };
  for (const bool empty_rank : {false, true})
  {
    const int count = (empty_rank && rank == size - 1 && size > 1) ? 0 : 37 + 23 * rank;
    const mfem::Vector xyz = MakeQueryPoints(count, rank);
    const auto unblocked = LocateDistributedPoints(
        *mesh, 3, xyz, std::numeric_limits<int>::max(), &element_value);
    // The cases exercise the fallback, the unresolved points and (at two or more ranks)
    // the multi-candidate routing.
    long long int totals[3] = {unblocked.fallback_queries, unblocked.candidate_queries,
                               static_cast<long long int>(count)};
    Mpi::GlobalSum(3, totals, comm);
    CHECK(totals[0] > 0);
    if (size > 1)
    {
      CHECK(totals[1] > totals[2]);
    }
    int unresolved = 0;
    for (const int owner : unblocked.owners)
    {
      unresolved += owner == size ? 1 : 0;
    }
    Mpi::GlobalSum(1, &unresolved, comm);
    CHECK(unresolved > 0);
    for (const int block : {1, 3, 7, 4096})
    {
      const auto blocked = LocateDistributedPoints(*mesh, 3, xyz, block, &element_value);
      CheckBitwiseEqual(blocked, unblocked);
    }
    const auto without_values =
        LocateDistributedPoints(*mesh, 3, xyz, std::numeric_limits<int>::max());
    const auto without_values_blocked = LocateDistributedPoints(*mesh, 3, xyz, 1);
    CheckBitwiseEqual(without_values_blocked, without_values);
    CHECK(without_values.owners == unblocked.owners);
  }
#endif
}

TEST_CASE("Unlocated trace points keep the clear error unless their patch is excluded late",
          "[surfaceresponseoperator][pointlocator][Serial][Parallel]")
{
  // Decisions 2156 / 2180 (the late domain-boundary exclusion's classification): planted
  // points outside the unit cube (spread over the ranks) with their exact distances to the
  // meshed domain, and one inside it. On this first-order simplex mesh a point of the
  // device metal (an edge point), a point beyond its coupon's contour radius and the point
  // inside the mesh (a search defect, decision 2174 m-2) keep the clear error; the others
  // are late exclusions. The same answer on every rank (collective).
  using Catch::Matchers::WithinAbs;
  using Decision = UnlocatedTracePointDecision;
  using Verdict = UnlocatedTracePointVerdict;
  MPI_Comm comm = Mpi::World();
  auto mesh = MakeSpreadParMesh(comm);
  const std::vector<std::array<double, 3>> points = {
      {1.05, 0.5, 0.5},     // 0.05 beyond the face x = 1
      {1.5, 0.5, 0.5},      // 0.5 beyond it: beyond the contour radius 0.2
      {1.05, 0.5, 0.5},     // 0.05 beyond it, a point of the device metal
      {1.03, 1.04, 0.5},    // 0.05 from the cube's edge x = y = 1
      {-0.06, -0.08, 1.0},  // 0.1 from the cube's corner (0, 0, 1)
      {0.5, 0.5, 0.5}};     // inside: distance 0
  const std::vector<bool> edge_points = {false, false, true, false, false, false};
  const std::vector<double> radii(points.size(), 0.2);
  const auto result = ClassifyUnlocatedTracePoints(*mesh, 3, points, radii);
  REQUIRE(result.distances.size() == points.size());
  REQUIRE(result.verdicts.size() == points.size());
  CHECK(result.first_order_simplex_mesh);
  CHECK(result.first_order_mesh);
  CHECK_THAT(result.distances[0], WithinAbs(0.05, 1.0e-14));
  CHECK(std::isinf(result.distances[1]));
  CHECK_THAT(result.distances[2], WithinAbs(0.05, 1.0e-14));
  CHECK_THAT(result.distances[3], WithinAbs(0.05, 1.0e-14));
  CHECK_THAT(result.distances[4], WithinAbs(0.1, 1.0e-14));
  CHECK(result.distances[5] == 0.0);
  const std::vector<Verdict> verdicts = {
      Verdict::WITHIN_CONTOUR_RADIUS, Verdict::BEYOND_CONTOUR_RADIUS,
      Verdict::WITHIN_CONTOUR_RADIUS, Verdict::WITHIN_CONTOUR_RADIUS,
      Verdict::WITHIN_CONTOUR_RADIUS, Verdict::INSIDE_MESH};
  const std::vector<Decision> decisions = {
      Decision::LATE_EXCLUSION, Decision::BEYOND_CONTOUR_RADIUS, Decision::OWN_EDGE_POINT,
      Decision::LATE_EXCLUSION, Decision::LATE_EXCLUSION,        Decision::INSIDE_MESH};
  for (std::size_t i = 0; i < points.size(); i++)
  {
    INFO("point " << i);
    CHECK(result.verdicts[i] == verdicts[i]);
    CHECK(DecideUnlocatedTracePoint(result.verdicts[i], edge_points[i], false,
                                    result.first_order_simplex_mesh,
                                    result.first_order_mesh) == decisions[i]);
  }
  // The root's distances on every rank (the ranks' minimum).
  std::vector<double> root = result.distances;
  for (auto &value : root)
  {
    value = std::isinf(value) ? -1.0 : value;
  }
  std::vector<double> mine = root;
  Mpi::Broadcast(static_cast<int>(root.size()), root.data(), 0, comm);
  CHECK(root == mine);

  // The same mesh curved (second-order nodes, decision 2180): the distances are element-box
  // lower bounds, so the inside point's 0 does not show that it is inside; it stays a late
  // exclusion, counted as undecided. The outside distances are the same here.
  mesh->SetCurvature(2);
  const auto curved = ClassifyUnlocatedTracePoints(*mesh, 3, points, radii);
  CHECK(!curved.first_order_simplex_mesh);
  CHECK(!curved.first_order_mesh);
  CHECK_THAT(curved.distances[0], WithinAbs(0.05, 1.0e-12));
  CHECK(std::isinf(curved.distances[1]));
  CHECK_THAT(curved.distances[3], WithinAbs(0.05, 1.0e-12));
  CHECK_THAT(curved.distances[4], WithinAbs(0.1, 1.0e-12));
  CHECK(curved.distances[5] == 0.0);
  CHECK(curved.verdicts[5] == Verdict::INSIDE_MESH);
  CHECK(DecideUnlocatedTracePoint(
            curved.verdicts[5], false, false, curved.first_order_simplex_mesh,
            curved.first_order_mesh) == Decision::LATE_EXCLUSION_UNDECIDED);
  // On the curved mesh a positive element-box distance is a lower bound: the point is
  // outside, a decided late exclusion (decision 2204).
  CHECK(curved.verdicts[0] == Verdict::WITHIN_CONTOUR_RADIUS);
  CHECK(DecideUnlocatedTracePoint(curved.verdicts[0], false, false,
                                  curved.first_order_simplex_mesh,
                                  curved.first_order_mesh) == Decision::LATE_EXCLUSION);

  // Decision 2204 (ERRATUM-2190 m-3): a first-order MIXED mesh, hexahedra in the unit cube
  // and tetrahedra in the cube shifted by 2 along x. A point inside a hexahedron has a
  // positive facet distance (a non-simplex element is taken to exclude the point), so on
  // such a mesh no distance decides inside or outside: every late point is undecided; the
  // point beyond its contour radius keeps the clear error.
  {
    mfem::Mesh hexahedra =
        mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
    mfem::Mesh tetrahedra =
        mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0);
    for (int v = 0; v < tetrahedra.GetNV(); v++)
    {
      tetrahedra.GetVertex(v)[0] += 2.0;
    }
    mfem::Mesh *pieces[2] = {&hexahedra, &tetrahedra};
    mfem::Mesh serial(pieces, 2);
    mfem::ParMesh mixed(comm, serial);
    const std::vector<std::array<double, 3>> mixed_points = {
        {0.6, 0.6, 0.6},   // inside a hexahedron: facet distance 0.1
        {1.05, 0.6, 0.6},  // 0.05 beyond the hexahedra's face x = 1
        {2.6, 0.6, 0.6},   // inside the tetrahedra: distance 0
        {1.5, 0.5, 0.5}};  // 0.5 from both blocks: beyond the contour radius 0.2
    const std::vector<double> mixed_radii(mixed_points.size(), 0.2);
    const auto first_order =
        ClassifyUnlocatedTracePoints(mixed, 3, mixed_points, mixed_radii);
    CHECK(!first_order.first_order_simplex_mesh);
    CHECK(first_order.first_order_mesh);
    CHECK_THAT(first_order.distances[0], WithinAbs(0.1, 1.0e-14));
    CHECK_THAT(first_order.distances[1], WithinAbs(0.05, 1.0e-14));
    CHECK(first_order.distances[2] == 0.0);
    CHECK(std::isinf(first_order.distances[3]));
    const std::vector<Decision> mixed_decisions = {
        Decision::LATE_EXCLUSION_UNDECIDED, Decision::LATE_EXCLUSION_UNDECIDED,
        Decision::LATE_EXCLUSION_UNDECIDED, Decision::BEYOND_CONTOUR_RADIUS};
    for (std::size_t i = 0; i < mixed_points.size(); i++)
    {
      INFO("mixed-mesh point " << i);
      CHECK(DecideUnlocatedTracePoint(first_order.verdicts[i], false, false,
                                      first_order.first_order_simplex_mesh,
                                      first_order.first_order_mesh) == mixed_decisions[i]);
    }
  }
}

TEST_CASE("Unlocated trace point decisions of the late domain-boundary exclusion",
          "[surfaceresponseoperator][pointlocator][Serial]")
{
  // Decision 2180: the decision of the late block for every combination of the point's
  // verdict, kind (a point of the device metal or a trace point), patch class (applied or
  // diagnostic-only) and mesh kind (first-order simplex, first-order with a non-simplex
  // element, or curved), and the reason the clear error gives. Rows: the gross mismatch
  // (BEYOND_CONTOUR_RADIUS), decision 2174 m-1 (a diagnostic-only patch keeps the clear
  // error) and m-2 (inside a first-order simplex mesh: the clear error; elsewhere
  // undecided, a late exclusion), and decision 2204 (ERRATUM-2190 m-3: on a first-order
  // mesh with a non-simplex element every late point is undecided).
  using Decision = UnlocatedTracePointDecision;
  using Verdict = UnlocatedTracePointVerdict;
  struct Row
  {
    Verdict verdict;
    bool edge_point, diagnostic_patch, first_order_simplex_mesh, first_order_mesh;
    Decision decision;
  };
  const Verdict W = Verdict::WITHIN_CONTOUR_RADIUS, B = Verdict::BEYOND_CONTOUR_RADIUS,
                I = Verdict::INSIDE_MESH;
  const Decision L = Decision::LATE_EXCLUSION, U = Decision::LATE_EXCLUSION_UNDECIDED,
                 E = Decision::OWN_EDGE_POINT, R = Decision::BEYOND_CONTOUR_RADIUS,
                 M = Decision::INSIDE_MESH, D = Decision::DIAGNOSTIC_PATCH;
  // Mesh kinds: (true, true) first-order simplex, (false, true) first-order with a
  // non-simplex element, (false, false) curved.
  const std::vector<Row> rows = {
      {W, false, false, true, true, L},   {W, false, false, false, true, U},
      {W, false, false, false, false, L}, {W, false, true, true, true, D},
      {W, false, true, false, true, D},   {W, false, true, false, false, D},
      {W, true, false, true, true, E},    {W, true, false, false, true, E},
      {W, true, false, false, false, E},  {W, true, true, true, true, E},
      {W, true, true, false, true, E},    {W, true, true, false, false, E},
      {B, false, false, true, true, R},   {B, false, false, false, true, R},
      {B, false, false, false, false, R}, {B, false, true, true, true, R},
      {B, false, true, false, true, R},   {B, false, true, false, false, R},
      {B, true, false, true, true, E},    {B, true, false, false, true, E},
      {B, true, false, false, false, E},  {B, true, true, true, true, E},
      {B, true, true, false, true, E},    {B, true, true, false, false, E},
      {I, false, false, true, true, M},   {I, false, false, false, true, U},
      {I, false, false, false, false, U}, {I, false, true, true, true, M},
      {I, false, true, false, true, D},   {I, false, true, false, false, D},
      {I, true, false, true, true, E},    {I, true, false, false, true, E},
      {I, true, false, false, false, E},  {I, true, true, true, true, E},
      {I, true, true, false, true, E},    {I, true, true, false, false, E}};
  REQUIRE(rows.size() == 36);
  for (std::size_t i = 0; i < rows.size(); i++)
  {
    const auto &row = rows[i];
    INFO("row " << i);
    CHECK(DecideUnlocatedTracePoint(row.verdict, row.edge_point, row.diagnostic_patch,
                                    row.first_order_simplex_mesh,
                                    row.first_order_mesh) == row.decision);
  }
  CHECK(UnlocatedTracePointErrorReason(Decision::LATE_EXCLUSION, 1.0).empty());
  CHECK(UnlocatedTracePointErrorReason(Decision::LATE_EXCLUSION_UNDECIDED, 1.0).empty());
  CHECK(UnlocatedTracePointErrorReason(Decision::OWN_EDGE_POINT, 1.0)
            .find("own metal edge") != std::string::npos);
  CHECK(UnlocatedTracePointErrorReason(Decision::BEYOND_CONTOUR_RADIUS, 2.687006)
            .find("contour radius 2.687006e+00 mesh units") != std::string::npos);
  CHECK(UnlocatedTracePointErrorReason(Decision::INSIDE_MESH, 1.0)
            .find("inside the meshed domain") != std::string::npos);
  CHECK(UnlocatedTracePointErrorReason(Decision::DIAGNOSTIC_PATCH, 1.0)
            .find("diagnostic-only patch") != std::string::npos);
}
