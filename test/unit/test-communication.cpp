// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstdint>
#include <numeric>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "utils/communication.hpp"

using namespace palace;

TEST_CASE("BroadcastLarge", "[communication][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int root = 0;

  SECTION("Chunked matches unchunked")
  {
    // Length chosen to exercise full chunks plus a partial remainder chunk.
    const std::int64_t len = 1000;
    std::vector<char> ref(len), buf(len);
    if (Mpi::Rank(comm) == root)
    {
      for (std::int64_t i = 0; i < len; i++)
      {
        ref[i] = static_cast<char>(i % 251);
      }
      buf = ref;
    }
    Mpi::Broadcast(static_cast<int>(len), ref.data(), root, comm);
    for (const std::int64_t chunk : {std::int64_t(1), std::int64_t(7), std::int64_t(333),
                                     std::int64_t(1000), std::int64_t(5000)})
    {
      std::vector<char> out = (Mpi::Rank(comm) == root) ? buf : std::vector<char>(len);
      Mpi::BroadcastLarge(len, out.data(), root, comm, chunk);
      REQUIRE(out == ref);
    }
  }

  SECTION("Zero length")
  {
    std::vector<char> buf;
    Mpi::BroadcastLarge(0, buf.data(), root, comm, 16);
    REQUIRE(buf.empty());
  }

  SECTION("Non-char type")
  {
    const std::int64_t len = 137;
    std::vector<double> buf(len, 0.0);
    if (Mpi::Rank(comm) == root)
    {
      std::iota(buf.begin(), buf.end(), 0.5);
    }
    Mpi::BroadcastLarge(len, buf.data(), root, comm, 16);
    std::vector<double> ref(len);
    std::iota(ref.begin(), ref.end(), 0.5);
    REQUIRE(buf == ref);
  }
}
