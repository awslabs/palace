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

TEST_CASE("PerfMetrics Handles", "[perfmetrics][Serial]")
{
  using perf::detail::entries;
  perf::SumMetric a("Test/Handles/Sum", perf::Red::MAX);
  perf::SumMetric a2("Test/Handles/Sum", perf::Red::MAX);  // Same metric
  perf::MaxMetric m("Test/Handles/Max", perf::Red::MAX);
  perf::LastMetric l("Test/Handles/Last", perf::Red::MAX, perf::Run::ITERATION_ONLY);
  perf::SumMetric invalid;
  perf::detail::ResetValues();
  a.Add(2.0);
  a2.Add(3.0);
  m.Update(1.0);
  m.Update(std::nan(""));
  m.Update(std::numeric_limits<double>::infinity());
  l.Set(4.0);
  l.Set(7.0);
  invalid.Add(1.0);
  CHECK(!invalid.Valid());
  const auto &es = entries;
  auto Find = [&](const char *p) -> const perf::detail::Entry &
  {
    for (const auto &e : es)
    {
      if (e.path == p)
      {
        return e;
      }
    }
    FAIL("missing " << p);
    return es.front();
  };
  CHECK(Find("Test/Handles/Sum").iter == 5.0);
  CHECK(Find("Test/Handles/Max").iter == DBL_MAX);
  CHECK(Find("Test/Handles/Last").iter == 7.0);
  perf::detail::FoldIteration();
  a.Add(1.0);
  l.Set(1.0);
  CHECK(Find("Test/Handles/Sum").iter == 1.0);
  CHECK(Find("Test/Handles/Sum").run_value == 5.0);
  CHECK(!Find("Test/Handles/Last").run_set);  // ITERATION_ONLY
  perf::detail::ResetValues();
  a.Add(1.0);  // Handles survive a reset
  CHECK(Find("Test/Handles/Sum").iter == 1.0);
  CHECK(!Find("Test/Handles/Sum").run_set);
}

TEST_CASE("PerfMetrics Paths", "[perfmetrics][Serial]")
{
  CHECK(perf::detail::ValidPath("LinearSolver/Main/Its"));
  CHECK(perf::detail::ValidPath("Tallies/X"));
  CHECK(!perf::detail::ValidPath(""));
  CHECK(!perf::detail::ValidPath("/A"));
  CHECK(!perf::detail::ValidPath("A/"));
  CHECK(!perf::detail::ValidPath("A//B"));
  CHECK(!perf::detail::ValidPath("ElapsedTime/X"));
  CHECK(!perf::detail::ValidPath("Cumulative"));
  CHECK(!perf::detail::ValidPath("Problem/AMR/Marked"));
  CHECK(perf::detail::ValidPath("ProblemSize"));  // Component-wise, not string prefix
}

TEST_CASE("PerfMetrics Roles", "[perfmetrics][Serial][Parallel]")
{
  CHECK(perf::CurrentRole() == "Unassigned");
  {
    perf::ScopedRole r("Outer");
    CHECK(perf::CurrentRole() == "Outer");
    const auto id = perf::CurrentRoleId();
    {
      perf::ScopedRole r2("Inner");
      CHECK(perf::CurrentRole() == "Inner");
      perf::ScopedRole r3(id);
      CHECK(perf::CurrentRole() == "Outer");
    }
    CHECK(perf::CurrentRole() == "Outer");
  }
  CHECK(perf::CurrentRole() == "Unassigned");
}

TEST_CASE("PerfMetrics Validation", "[perfmetrics][Serial]")
{
  using perf::Acc, perf::Red, perf::Run;
  std::vector<perf::MetricInfo> global = {
      {"LinearSolver/Main/Its", Acc::SUM, Red::MAX, Run::ACCUMULATE}};
  auto Line = [](const char *p, char a = '0', char r = '2', char run = '0')
  { return std::string(p) + '\t' + a + r + run + '\n'; };
  std::string error;
  // New paths, a duplicate within the batch, and one already global: sorted, unique.
  auto out = perf::detail::ValidateNew(
      global, Line("B/X") + Line("A/Y") + Line("B/X") + Line("LinearSolver/Main/Its"),
      error);
  CHECK(error.empty());
  REQUIRE(out.size() == 2);
  CHECK(out[0].starts_with("A/Y\t"));
  CHECK(out[1].starts_with("B/X\t"));
  // Component-wise prefixes conflict; string prefixes do not.
  perf::detail::ValidateNew(global, Line("LinearSolver/Main/ItsHistogram/B1"), error);
  CHECK(error.empty());
  perf::detail::ValidateNew(global, Line("LinearSolver/Main/Its/B1"), error);
  CHECK(!error.empty());
  error.clear();
  perf::detail::ValidateNew(global, Line("LinearSolver/Main"), error);
  CHECK(!error.empty());
  error.clear();
  perf::detail::ValidateNew(global, Line("C") + Line("C/D"), error);
  CHECK(!error.empty());
  error.clear();
  // Different rules for one path.
  perf::detail::ValidateNew(global, Line("LinearSolver/Main/Its", '1'), error);
  CHECK(!error.empty());
  error.clear();
  perf::detail::ValidateNew(global, Line("E", '0') + Line("E", '0', '4'), error);
  CHECK(!error.empty());
}

TEST_CASE("PerfMetrics Sync", "[perfmetrics][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm);
  perf::detail::ResetValues();
  // Registered on every process but the root, in rank-dependent order.
  perf::SumMetric a, b;
  if (rank > 0)
  {
    if (rank % 2)
    {
      a = perf::SumMetric("Test/Sync/A", perf::Red::STATS);
      b = perf::SumMetric("Test/Sync/B", perf::Red::SUM);
    }
    else
    {
      b = perf::SumMetric("Test/Sync/B", perf::Red::SUM);
      a = perf::SumMetric("Test/Sync/A", perf::Red::STATS);
    }
    a.Add(rank);
    b.Add(1.0);
  }
  perf::Sync(comm);
  const auto &g = perf::GlobalMetrics();
  int ia = -1, ib = -1;
  for (std::size_t i = 0; i < g.size(); i++)
  {
    ia = (g[i].path == "Test/Sync/A") ? static_cast<int>(i) : ia;
    ib = (g[i].path == "Test/Sync/B") ? static_cast<int>(i) : ib;
  }
  if (np == 1)
  {
    CHECK(ia < 0);
    CHECK(ib < 0);
    return;
  }
  REQUIRE(ia >= 0);
  REQUIRE(ib >= 0);
  CHECK(g[ia].red == perf::Red::STATS);
  std::vector<double> v;
  std::vector<char> s;
  perf::PackMetrics(false, v, s);
  const auto r = perf::ReduceStats(v, s, comm);
  CHECK(r[ia].ranks == np - 1);
  CHECK(r[ia].min == 1.0);
  CHECK(r[ib].sum == np - 1);
}
