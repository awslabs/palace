// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "utils/communication.hpp"
#include "utils/memoryreporting.hpp"

using namespace palace;
using namespace palace::memory_reporting;

TEST_CASE("Local Memory Queries", "[memoryreporting][Serial][Parallel]")
{
  auto current = GetCurrentMemory();
  auto peak = GetPeakMemory();

  CHECK(current > 0);
  CHECK(peak > 0);
  CHECK(peak >= current);
}

TEST_CASE("Peak Memory Is Non-Decreasing", "[memoryreporting][Serial][Parallel]")
{
  auto peak_before = GetPeakMemory();

  // Allocate ~1 MB and touch it to ensure it becomes resident.
  volatile std::vector<char> buf(1024 * 1024, 1);

  auto peak_after = GetPeakMemory();
  CHECK(peak_after >= peak_before);
}

TEST_CASE("Memory Stats Single Process", "[memoryreporting][Serial]")
{
  auto stats = GetCurrentMemoryStats(MPI_COMM_WORLD);

  CHECK(stats.min > 0);
  CHECK(stats.min == stats.max);
  CHECK(stats.sum == stats.min);
  CHECK(stats.avg == static_cast<double>(stats.min));
}

TEST_CASE("Memory Stats Multi Process", "[memoryreporting][Parallel]")
{
  auto comm = MPI_COMM_WORLD;
  auto np = Mpi::Size(comm);

  auto current_stats = GetCurrentMemoryStats(comm);
  CHECK(current_stats.min > 0);
  CHECK(current_stats.min <= current_stats.max);
  CHECK(current_stats.sum >= current_stats.min * np);
  CHECK(current_stats.sum <= current_stats.max * np);
  CHECK(current_stats.avg >= static_cast<double>(current_stats.min));
  CHECK(current_stats.avg <= static_cast<double>(current_stats.max));

  auto peak_stats = GetPeakMemoryStats(comm);
  CHECK(peak_stats.min > 0);
  CHECK(peak_stats.min <= peak_stats.max);
  CHECK(peak_stats.sum >= peak_stats.min * np);
  CHECK(peak_stats.sum <= peak_stats.max * np);
}
