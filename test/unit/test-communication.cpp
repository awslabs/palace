// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <string_view>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "utils/communication.hpp"

using namespace palace;

TEST_CASE("BroadcastLarge", "[communication][Parallel]")
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

// OI-SCALE-GATHER (idfix DESIGN section 9 (ii)): the 64-bit-sized variable-length gathers /
// scatter / reduction in rounds are element-for-element the one-call forms (the P == S
// identity) at every chunk, including chunks smaller than any rank's piece, a chunk that
// splits the pieces unevenly, one above every piece, and ranks with empty pieces.
TEST_CASE("LargeMessageRounds", "[communication][scalegather][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), size = Mpi::Size(comm);
  // Rank r holds 13 r + 5 values (rank 0: 5; the last rank's piece the largest), tagged by
  // rank and position so that a misplaced element is visible.
  auto Local = [&](int r)
  {
    std::vector<double> values(13 * r + 5);
    for (std::size_t i = 0; i < values.size(); i++)
    {
      values[i] = 1000.0 * r + static_cast<double>(i) + 0.25;
    }
    return values;
  };
  const std::vector<double> local = Local(rank);
  std::vector<double> reference;
  std::vector<std::int64_t> reference_counts;
  for (int r = 0; r < size; r++)
  {
    const auto piece = Local(r);
    reference.insert(reference.end(), piece.begin(), piece.end());
    reference_counts.push_back(static_cast<std::int64_t>(piece.size()));
  }
  for (const std::int64_t chunk : {std::int64_t(1), std::int64_t(3), std::int64_t(7),
                                   std::int64_t(64), std::int64_t(1 << 20)})
  {
    Mpi::LargeGatherPlan plan;
    const auto all = Mpi::AllgatherLarge(local, comm, &plan, chunk);
    CHECK(all == reference);
    CHECK(plan.counts == reference_counts);
    CHECK(plan.total == static_cast<std::int64_t>(reference.size()));
    const auto gathered = Mpi::GatherLarge(local, 0, comm, nullptr, chunk);
    if (rank == 0)
    {
      CHECK(gathered == reference);
    }
    else
    {
      CHECK(gathered.empty());
    }
    // The scatter of the gathered picture returns every rank its own piece.
    const auto scattered = Mpi::ScattervLarge(
        gathered, rank == 0 ? reference_counts : std::vector<std::int64_t>{}, 0, comm,
        chunk);
    CHECK(scattered == local);
    // The chunked MAX reduction = the one-call reduction.
    std::vector<unsigned char> flags(257, 0);
    for (std::size_t i = rank; i < flags.size(); i += size + 1)
    {
      flags[i] = 1;
    }
    std::vector<unsigned char> expected = flags;
    Mpi::GlobalMax(static_cast<int>(expected.size()), expected.data(), comm);
    Mpi::GlobalMaxLarge(flags, comm, chunk);
    CHECK(flags == expected);
  }
  // An empty piece on every rank: empty results, a zero plan.
  Mpi::LargeGatherPlan empty_plan;
  CHECK(Mpi::AllgatherLarge(std::vector<double>{}, comm, &empty_plan, 4).empty());
  CHECK(empty_plan.total == 0);
  CHECK(Mpi::GatherLarge(std::vector<double>{}, 0, comm, nullptr, 4).empty());
  std::vector<unsigned char> none;
  Mpi::GlobalMaxLarge(none, comm, 4);
  CHECK(none.empty());
  // The default chunk is INT_MAX / 8 unless PALACE_LARGE_MESSAGE_CHUNK says otherwise.
  CHECK(Mpi::LargeChunk() > 0);
}

// Decision 346 (b): an exit() from a library in the middle of a run (gslib's die() on one
// rank) must abort every rank at once instead of blocking in MPI_Finalize from Mpi's
// destructor. The case only runs under PALACE_TEST_ABNORMAL_EXIT=1 (the ctest entry
// failfast-AbnormalExitAbortsAllRanks sets it and passes on the abort line within its
// timeout): in a sweep it is skipped, since passing means the launcher reports a failed
// job.
TEST_CASE("AbnormalExitAbortsAllRanks", "[communication][Parallel]")
{
  const char *enabled = std::getenv("PALACE_TEST_ABNORMAL_EXIT");
  if (!enabled || std::string_view(enabled) != "1")
  {
    SKIP("Set PALACE_TEST_ABNORMAL_EXIT=1 to exit one rank mid-run (the job must abort)");
  }
  MPI_Comm comm = Mpi::World();
  Mpi::Barrier(comm);
  if (Mpi::Rank(comm) == Mpi::Size(comm) - 1)
  {
    std::exit(1);
  }
  // Without the abort these ranks wait here forever (the R-B hang of PBS 56421).
  Mpi::Barrier(comm);
  FAIL("The abnormal exit of the last rank did not abort this rank!");
}
