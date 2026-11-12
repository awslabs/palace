// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/communication.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/perfprobes.hpp"
#include "utils/perfstats.hpp"
#include "utils/timer.hpp"

using namespace palace;
using Catch::Matchers::WithinAbs;

namespace
{

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

}  // namespace

TEST_CASE("PerfProbes Memory Attribution", "[perfprobes][perfstats][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  BlockTimer::Reset();
  std::vector<char> keep;
  {
    BlockTimer bt(Timer::CONSTRUCT);
    perf::ScopedMemoryProbe probe("Test/Alloc");
    keep.assign(64 << 20, 1);            // Persistent
    std::vector<char> tmp(32 << 20, 1);  // Transient
  }
  std::vector<char> outside(16 << 20, 1);  // Outside any timer row
  BlockTimer::Finalize(comm);
  const auto &r = perf::GetResults();

  // Baseline + rows + unattributed = peak, per rank (sums over ranks).
  double rows = 0.0;
  for (int i = 0; i < Timer::NUM_TIMINGS; i++)
  {
    if (i != Timer::TOTAL)
    {
      rows += r.iter.rank_mem[i].sum;
    }
  }
  CHECK_THAT(r.mem_baseline.sum + rows + r.mem_unattributed.sum,
             WithinAbs(r.peak_mem.sum, 1.0));
  CHECK_THAT(r.node_mem_baseline.sum + r.node_mem_unattributed.sum + rows,
             WithinAbs(r.peak_node_mem.sum, 1.0));
#if defined(__linux__)
  const int persistent = CounterIndex(r, "Memory/Test/Alloc/PersistentMegabytes");
  const int events = CounterIndex(r, "Memory/Test/Alloc/Events");
  REQUIRE(persistent >= 0);
  REQUIRE(events >= 0);
  CHECK(r.counter[persistent].min > 50.0);
  CHECK(r.counter[events].max == 1.0);
#endif
  BlockTimer::Reset();
}

TEST_CASE("PerfProbes Peak Memory Of Later Iterations", "[perfprobes][perfstats][Serial]")
{
  // The peak RSS of the process is not reset at the start of an iteration (the reset would
  // also change the peak seen by getrusage and the accounting of job schedulers). The peak
  // of an iteration which does not raise it is flagged as a lower bound.
  MPI_Comm comm = Mpi::World();
  BlockTimer::Reset();
  static char *volatile escape = nullptr;  // Keeps the allocation from being elided
  {
    constexpr std::size_t n = 256 << 20;
    std::unique_ptr<char[]> big(new char[n]);
    std::memset(big.get(), 1, n);
    escape = big.get();
  }
  const long peak = memory_reporting::GetPeakMemory();
  if (peak < memory_reporting::GetCurrentMemory() + (64l << 20))
  {
    // For example on macOS, where the peak and the current memory are different measures.
    BlockTimer::Reset();
    SKIP("The peak RSS is not above the RSS after the allocation is freed");
  }
  BlockTimer::Finalize(comm);
  CHECK(perf::GetResults().peak_lower_bound == 0);
  BlockTimer::StartIteration();
  CHECK(memory_reporting::GetPeakMemory() >= peak);
  {
    BlockTimer bt(Timer::IO);
    std::vector<char> small(1 << 20, 1);
    escape = small.data();
  }
  BlockTimer::Finalize(comm);
  const auto &r = perf::GetResults();
  CHECK(r.peak_lower_bound == 1);
  CHECK(r.peak_mem.max < peak);
  CHECK(r.cum_peak_mem.max >= peak);
  nlohmann::json meta = nlohmann::json::object();
  perf::WriteJson(meta);
  CHECK(meta["PeakMemoryMegabytes"]["LowerBound"] == true);
  CHECK(meta["PeakNodeMemoryMegabytes"]["LowerBound"] == true);
  CHECK(!meta["Cumulative"]["PeakMemoryMegabytes"].contains("LowerBound"));
#if defined(__linux__)
  // The high-water mark of the kernel is kept as well.
  std::ifstream status("/proc/self/status");
  std::string line;
  long hwm = 0;
  while (std::getline(status, line))
  {
    if (line.starts_with("VmHWM:"))
    {
      hwm = std::stol(line.substr(6)) * 1024;
    }
  }
  CHECK(hwm >= peak);
#endif
  escape = nullptr;
  BlockTimer::Reset();
}

TEST_CASE("PerfProbes Peak Memory After Reset", "[perfprobes][perfstats][Serial]")
{
  // A peak of the process reached before a Reset (by a previous run in the same process)
  // is not the peak of the next run: it is flagged as a lower bound, bounded below by the
  // RSS sampled in the run, such as by a memory probe.
  MPI_Comm comm = Mpi::World();
  static char *volatile escape = nullptr;  // Keeps the allocation from being elided
  {
    constexpr std::size_t n = 256 << 20;
    std::unique_ptr<char[]> big(new char[n]);
    std::memset(big.get(), 1, n);
    escape = big.get();
  }
  const long peak = memory_reporting::GetPeakMemory();
  if (peak < memory_reporting::GetCurrentMemory() + (64l << 20))
  {
    SKIP("The peak RSS is not above the RSS after the allocation is freed");
  }
  BlockTimer::Reset();
  std::vector<char> keep;
  long rss = 0;
  {
    BlockTimer bt(Timer::CONSTRUCT);
    perf::ScopedMemoryProbe probe("Test/AfterReset");
    keep.assign(16 << 20, 1);
    rss = memory_reporting::GetCurrentMemory();
  }
  keep = std::vector<char>();
  BlockTimer::Finalize(comm);
  const auto &r = perf::GetResults();
  CHECK(r.peak_lower_bound == 1);
  CHECK(r.peak_mem.max < peak);
  CHECK(r.peak_mem.max >= 0.99 * rss);  // The RSS at the end of the probe
  escape = nullptr;
  BlockTimer::Reset();
}

TEST_CASE("PerfProbes Node Communicators", "[perfprobes][perfstats][Serial][Parallel]")
{
  // The node communicators are cached for the group of a communicator: one of another
  // group gets its own, also if it reuses the handle of a freed communicator.
  MPI_Comm world = Mpi::World();
  MPI_Comm dup, single;
  MPI_Comm_dup(world, &dup);
  const int node_size = Mpi::Size(memory_reporting::GetNodeComms(dup).node);
  MPI_Comm_free(&dup);
  MPI_Comm_split(world, Mpi::Rank(world), 0, &single);
  const auto &nc = memory_reporting::GetNodeComms(single);
  CHECK(Mpi::Size(nc.node) == 1);
  CHECK(nc.node_rank == 0);
  CHECK(nc.leaders != MPI_COMM_NULL);
  MPI_Comm_free(&single);
  CHECK(Mpi::Size(memory_reporting::GetNodeComms(world).node) == node_size);
}

TEST_CASE("PerfProbes Assembly Tally", "[perfprobes][Serial]")
{
  BlockTimer::Reset();
  int id;
  {
    BlockTimer bt(Timer::CONSTRUCT);
    id = perf::AssemblyTally(perf::AssemblyKind::PARTIAL);
    CHECK(perf::AssemblyTally(perf::AssemblyKind::PARTIAL) == id);
    CHECK(perf::AssemblyTally(perf::AssemblyKind::FULL) != id);
    {
      BlockTimer bt1(Timer::KSP_SETUP);
      CHECK(perf::AssemblyTally(perf::AssemblyKind::PARTIAL) != id);
    }
  }
  CHECK(id >= 0);
  BlockTimer::Reset();
}
