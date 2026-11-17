// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cfloat>
#include <cmath>
#include <limits>
#include <catch2/catch_test_macros.hpp>
#include "utils/communication.hpp"
#include "utils/perfmetrics.hpp"

using namespace palace;

TEST_CASE("PerfMetrics ReduceStats", "[perfmetrics][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm);
  // Value 0: rank + 1 on every process. Value 1: set on even ranks only. Value 2: nowhere.
  std::vector<double> v = {rank + 1.0, 10.0 * rank, 5.0};
  std::vector<char> s = {1, static_cast<char>(rank % 2 == 0), 0};
  const auto r = perf::ReduceStats(v, s, comm);
  REQUIRE(r.size() == 3);
  CHECK(r[0].min == 1.0);
  CHECK(r[0].max == np);
  CHECK(r[0].max_rank == np - 1);
  CHECK(r[0].sum == np * (np + 1) / 2.0);
  CHECK(r[0].ranks == np);
  const int evens = (np + 1) / 2;
  CHECK(r[1].ranks == evens);
  CHECK(r[1].min == 0.0);
  CHECK(r[1].max == 10.0 * 2 * (evens - 1));
  CHECK(r[2].ranks == 0);
  CHECK(perf::ReduceStats({}, {}, comm).empty());
}
