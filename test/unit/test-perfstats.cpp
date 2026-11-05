// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <thread>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/communication.hpp"
#include "utils/perfstats.hpp"
#include "utils/timer.hpp"

using namespace palace;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{

void Sleep(int ms)
{
  std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

int CounterIndex(const perf::Results &r, const std::string &path)
{
  for (std::size_t i = 0; i < r.counter_paths.size(); i++)
  {
    if (r.counter_paths[i] == path)
    {
      return static_cast<int>(i);
    }
  }
  return -1;
}

int TallyIndex(const perf::TimerResults &r, const std::string &name)
{
  for (std::size_t i = 0; i < r.tally_names.size(); i++)
  {
    if (r.tally_names[i] == name)
    {
      return static_cast<int>(i);
    }
  }
  return -1;
}

}  // namespace

TEST_CASE("PerfStats Reductions", "[perfstats][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm);
  BlockTimer::Reset();

  const int tally = perf::RegisterTally("Test/Tally", Timer::KSP_PRECONDITIONER);
  CHECK(perf::RegisterTally("Test/Tally", Timer::KSP_PRECONDITIONER) == tally);
  perf::AddTally(tally, 1.0 + rank, 2);
  const int stats = perf::RegisterCounter("Test/Rows", perf::Reduce::STATS);
  perf::Set(stats, 10.0 * (rank + 1));
  const int sum = perf::RegisterCounter("Test/Bytes", perf::Reduce::SUM);
  perf::Add(sum, 3.0);
  perf::Add(sum, 4.0);
  const int mx =
      perf::RegisterCounter("Test/MaxIts", perf::Reduce::MAX, perf::Accumulate::MAX);
  perf::SetMax(mx, rank);
  perf::SetMax(mx, 0.5);
  perf::RegisterCounter("Test/Unset");

  BlockTimer::Finalize(comm);
  const auto &r = perf::GetResults();
  CHECK(r.iteration == 1);
  const int t = TallyIndex(r.iter, "Test/Tally");
  REQUIRE(t >= 0);
  CHECK(r.iter.tally_parents[t] == Timer::KSP_PRECONDITIONER);
  CHECK_THAT(r.iter.tally_time[t].min, WithinAbs(1.0, 1e-12));
  CHECK_THAT(r.iter.tally_time[t].max, WithinAbs(np, 1e-12));
  CHECK_THAT(r.iter.tally_time[t].avg, WithinAbs(0.5 * (np + 1), 1e-12));
  CHECK(r.iter.tally_time[t].max_rank == np - 1);
  CHECK(r.iter.tally_count[t].max == 2);

  const int s = CounterIndex(r, "Test/Rows");
  REQUIRE(s >= 0);
  CHECK(r.counter[s].min == 10.0);
  CHECK(r.counter[s].max == 10.0 * np);
  CHECK(r.counter[s].max_rank == np - 1);
  CHECK(r.counter[s].n == np);
  CHECK(r.counter[CounterIndex(r, "Test/Bytes")].sum == 7.0 * np);
  CHECK(r.counter[CounterIndex(r, "Test/MaxIts")].max == std::max(np - 1.0, 0.5));
  CHECK(r.counter[CounterIndex(r, "Test/Unset")].n == 0);
  CHECK(r.peak_mem.n == np);
  CHECK(r.peak_node_mem.n == r.num_nodes);
  CHECK(r.num_nodes >= 1);
  BlockTimer::Reset();
}

TEST_CASE("PerfStats Inconsistent Registration", "[perfstats][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm);
  BlockTimer::Reset();
  if (rank == np - 1)
  {
    perf::AddTally(perf::RegisterTally(fmt::format("Test/OnlyRank{}", rank), Timer::IO),
                   2.0);
    perf::Set(perf::RegisterCounter("Test/OnlyLast", perf::Reduce::STATS), 5.0);
  }
  perf::Add(perf::RegisterCounter("Test/All", perf::Reduce::SUM), 1.0);
  BlockTimer::Finalize(comm);
  const auto &r = perf::GetResults();
  const int t = TallyIndex(r.iter, fmt::format("Test/OnlyRank{}", np - 1));
  REQUIRE(t >= 0);
  CHECK(r.iter.tally_time[t].max == 2.0);
  CHECK(r.iter.tally_time[t].min == ((np > 1) ? 0.0 : 2.0));
  CHECK(r.iter.tally_time[t].max_rank == np - 1);
  const int c = CounterIndex(r, "Test/OnlyLast");
  REQUIRE(c >= 0);
  CHECK(r.counter[c].n == 1);
  CHECK(r.counter[c].min == 5.0);
  CHECK(r.counter[CounterIndex(r, "Test/All")].sum == np);
  BlockTimer::Reset();
}

TEST_CASE("PerfStats Per-Iteration And Cumulative", "[perfstats][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int np = Mpi::Size(comm);
  BlockTimer::Reset();
  const int its =
      perf::RegisterCounter("Test/Its", perf::Reduce::MAX, perf::Accumulate::SUM);
  const int once = perf::RegisterCounter("Test/Once", perf::Reduce::SUM);
  const int size =
      perf::RegisterCounter("Test/Size", perf::Reduce::MAX, perf::Accumulate::NONE);
  const int tally = perf::RegisterTally("Test/Iter", Timer::IO);
  {
    BlockTimer bt(Timer::INIT);
    Sleep(5);
  }
  {
    BlockTimer bt(Timer::IO);
    perf::ScopedTally st(tally);
    Sleep(20);
  }
  perf::Add(its, 3);
  perf::Add(once, 1);
  perf::Set(size, 100);
  perf::SetInfo("Problem/Test/Info", "first");
  BlockTimer::Finalize(comm);
  nlohmann::json meta = {{"Problem", {{"MPISize", np}}}};
  perf::WriteJson(meta);
  const double io1 = perf::GetResults().iter.time[Timer::IO].avg;
  const double total1 = perf::GetResults().iter.time[Timer::TOTAL].avg;
  CHECK(io1 >= 0.019);
  CHECK(perf::GetResults().iter.count[Timer::STARTUP].max == 1);
  CHECK(meta["SchemaVersion"] == perf::SCHEMA_VERSION);
  CHECK(meta["Test"]["Once"] == np);
  CHECK(meta["Problem"]["Test"]["Info"] == "first");
  CHECK(meta["Problem"]["MPISize"] == np);
  CHECK(meta["ElapsedTime"]["Counts"]["IO"] == 1);
  CHECK(meta["ElapsedTime"]["Tallies"]["Test/Iter"]["Parent"] == "IO");

  BlockTimer::StartIteration();
  {
    BlockTimer bt(Timer::ADAPTATION);
    Sleep(5);
  }
  {
    BlockTimer bt(Timer::IO);
    perf::ScopedTally st(tally);
    Sleep(10);
  }
  perf::Add(its, 4);
  perf::Set(size, 200);
  BlockTimer::Finalize(comm);
  perf::WriteJson(meta);
  const auto &r = perf::GetResults();
  CHECK(r.iteration == 2);
  CHECK(r.iter.count[Timer::IO].max == 1);
  CHECK(r.cum.count[Timer::IO].max == 2);
  CHECK(r.iter.count[Timer::STARTUP].max == 0);
  CHECK(r.iter.count[Timer::INIT].max == 0);
  CHECK(r.cum.count[Timer::INIT].max == 1);
  CHECK_THAT(r.cum.time[Timer::IO].avg, WithinRel(io1 + r.iter.time[Timer::IO].avg, 1e-12));
  CHECK_THAT(r.cum.time[Timer::TOTAL].avg,
             WithinRel(total1 + r.iter.time[Timer::TOTAL].avg, 1e-12));
  CHECK(r.iter.time[Timer::TOTAL].avg < r.cum.time[Timer::TOTAL].avg);
  CHECK(r.cum_peak_mem.max >= r.peak_mem.max);
  CHECK(r.cum_peak_node_mem.max >= r.peak_node_mem.max);

  // Per-iteration values replace the previous iteration's, stale keys are removed, and
  // Cumulative holds the totals.
  CHECK(meta["Test"]["Its"] == 4);
  CHECK(meta["Test"]["Size"] == 200);
  CHECK(!meta["Test"].contains("Once"));
  CHECK(!meta["Problem"].contains("Test"));
  CHECK(meta["Cumulative"]["Test"]["Its"] == 7);
  CHECK(meta["Cumulative"]["Test"]["Once"] == np);
  CHECK(!meta["Cumulative"]["Test"].contains("Size"));
  CHECK(meta["ElapsedTime"]["Tallies"]["Test/Iter"]["Count"] == 1);
  CHECK(meta["Cumulative"]["ElapsedTime"]["Tallies"]["Test/Iter"]["Count"] == 2);
  CHECK(meta["ElapsedTime"]["Counts"]["ADAPTATION"] == 1);
  CHECK(meta["Cumulative"]["ElapsedTime"]["Counts"]["IO"] == 2);
  CHECK(meta["ElapsedTime"]["Durations"]["IO"].is_number());
  CHECK(meta["ElapsedTime"]["MaxRank"]["IO"].is_number_integer());
  CHECK(meta["PeakMemoryMegabytes"]["Max"].get<double>() <=
        meta["Cumulative"]["PeakMemoryMegabytes"]["Max"].get<double>());
  BlockTimer::Reset();
}

TEST_CASE("PerfStats Roles", "[perfstats][Serial]")
{
  CHECK(perf::CurrentRole() == "Main");
  {
    perf::ScopedRole role("DivFree");
    CHECK(perf::CurrentRole() == "DivFree");
    {
      perf::ScopedRole inner("Estimator");
      CHECK(perf::CurrentRole() == "Estimator");
    }
    CHECK(perf::CurrentRole() == "DivFree");
  }
  CHECK(perf::CurrentRole() == "Main");
  CHECK(perf::LinearSolverPath("Main", "Solves") == "LinearSolver/TotalSolves");
  CHECK(perf::LinearSolverPath("Main", "Its") == "LinearSolver/TotalIts");
  CHECK(perf::LinearSolverPath("Main", "MaxIts") == "LinearSolver/Main/MaxIts");
  CHECK(perf::LinearSolverPath("DivFree", "Solves") == "LinearSolver/DivFree/Solves");
}
