// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cstring>
#include <functional>
#include <limits>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
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

void CheckBitwiseEqual(const DistributedPointLocation &a, const DistributedPointLocation &b)
{
  CHECK(a.owners == b.owners);
  CHECK(a.elements == b.elements);
  REQUIRE(a.references.size() == b.references.size());
  CHECK(std::memcmp(a.references.data(), b.references.data(),
                    a.references.size() * sizeof(double)) == 0);
  REQUIRE(a.owner_values.size() == b.owner_values.size());
  CHECK(std::memcmp(a.owner_values.data(), b.owner_values.data(),
                    a.owner_values.size() * sizeof(double)) == 0);
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
