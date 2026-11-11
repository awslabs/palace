// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <thread>
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
    keep.assign(64 << 20, 1);  // Persistent
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

TEST_CASE("PerfProbes Heap Reuse", "[perfprobes][perfstats][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  constexpr double to_mb = 1.0 / (1024.0 * 1024.0);
  BlockTimer::Reset();
  BlockTimer::SetHeapBudget(1.0, 1.0);  // Many free chunks: each sample may take a while

  // Leave 32 MB of freed but resident heap: small chunks (below any mmap threshold), every
  // other one freed so that the freed chunks can neither coalesce nor shrink the heap. The
  // allocator may still return their pages to the operating system; the RSS check below
  // only applies if it did not.
  constexpr std::size_t chunk = 4 << 10, n = 16384;
  std::vector<std::unique_ptr<char[]>> chunks(n), reuse;
  for (auto &c : chunks)
  {
    c.reset(new char[chunk]);
    std::memset(c.get(), 1, chunk);
  }
  const long rss_before = memory_reporting::GetCurrentMemory();
  for (std::size_t i = 0; i < n; i += 2)
  {
    chunks[i].reset();
  }
  const bool resident =
      rss_before - memory_reporting::GetCurrentMemory() < static_cast<long>(8 << 20);
  reuse.reserve(n / 2);
  {
    BlockTimer bt(Timer::CONSTRUCT);
    perf::ScopedMemoryProbe probe("Test/Reuse");
    for (std::size_t i = 0; i < n / 2; i++)
    {
      reuse.emplace_back(new char[chunk]);
      std::memset(reuse.back().get(), 2, chunk);
    }
  }
  BlockTimer::Finalize(comm);
  const auto &r = perf::GetResults();

  // Baseline + rows + unattributed = heap in use at the end, per rank (sums over ranks).
  double rows = 0.0;
  for (int i = 0; i < Timer::NUM_TIMINGS; i++)
  {
    if (i != Timer::TOTAL)
    {
      rows += r.heap_growth[i].sum;
    }
  }
  CHECK_THAT(r.heap_baseline.sum + rows + r.heap_unattributed.sum,
             WithinAbs(r.heap_mem.sum, 1.0));
  CHECK_THAT(rows + r.heap_unattributed.sum,
             WithinAbs(r.heap_growth[Timer::TOTAL].sum, 1.0));
  CHECK(r.counter[CounterIndex(r, "Memory/Heap/Samples")].max > 0);
  if (!memory_reporting::HeapMemorySource().starts_with("rss"))
  {
    // The allocations of the phase show up in the heap, whatever memory serves them.
    CHECK(r.heap_growth[Timer::CONSTRUCT].min * to_mb > 30.0);
    const int heap = CounterIndex(r, "Memory/Test/Reuse/HeapPersistentMegabytes");
    REQUIRE(heap >= 0);
    CHECK(r.counter[heap].min > 30.0);
  }
#if defined(__GLIBC__)
  // They were served from the freed heap: the RSS did not grow.
  const int rss = CounterIndex(r, "Memory/Test/Reuse/PersistentMegabytes");
  REQUIRE(rss >= 0);
  if (resident)
  {
    CHECK(r.counter[rss].max < 8.0);
  }
  else
  {
    WARN("The freed heap was returned to the operating system: RSS check skipped");
  }
#else
  static_cast<void>(resident);
#endif
  nlohmann::json meta = nlohmann::json::object();
  perf::WriteJson(meta);
  CHECK(meta["HeapGrowthMegabytes"]["Min"]["CONSTRUCT"].is_number());
  CHECK(meta["HeapGrowthMegabytes"]["Sum"].contains("Baseline"));
  CHECK(meta["HeapGrowthMegabytes"]["Coverage"]["CONSTRUCT"] == 1.0);
  CHECK(meta["HeapMegabytes"]["Max"].get<double>() > 0.0);
  CHECK(meta["Memory"]["HeapSource"] == memory_reporting::HeapMemorySource());
  BlockTimer::Reset();
}

TEST_CASE("PerfProbes Heap Skipped Samples", "[perfprobes][Serial]")
{
  auto Sleep = [](int ms) { std::this_thread::sleep_for(std::chrono::milliseconds(ms)); };

  // A sample taken, without a gap: the interval of the phase is covered.
  Timer timer;
  Sleep(2);
  timer.MarkTime(Timer::IO);
  CHECK(timer.MarkHeap(Timer::IO, Timer::NUM_TIMINGS, true));
  CHECK(timer.HeapCoverage(Timer::IO) == 1.0);

  // With no budget, the first sample which can be skipped is taken, the next ones not.
  timer.SetHeapBudget(0.0, 0.0);
  timer.MarkTime(Timer::KSP);
  CHECK(timer.MarkHeap(Timer::KSP, Timer::NUM_TIMINGS));
  timer.MarkTime(Timer::KSP);
  CHECK(!timer.MarkHeap(Timer::KSP, Timer::NUM_TIMINGS));
  CHECK(timer.HeapSkipped() == 1);
  CHECK(timer.HeapQueryTime() > 0.0);

  // After skipped samples, the growth goes to the phase open throughout the gap, here
  // EPS, and the intervals of the other phases are not covered.
  Sleep(2);
  timer.MarkTime(Timer::IO);
  CHECK(!timer.MarkHeap(Timer::IO, Timer::NUM_TIMINGS));
  const long eps = timer.HeapData(Timer::EPS), io = timer.HeapData(Timer::IO);
  std::vector<char> keep(4 << 20, 1);
  Sleep(2);
  timer.MarkTime(Timer::EPS);
  CHECK(timer.MarkHeap(Timer::EPS, Timer::EPS, true));
  CHECK(timer.HeapData(Timer::IO) == io);
  if (!memory_reporting::HeapMemorySource().starts_with("rss"))
  {
    CHECK(timer.HeapData(Timer::EPS) - eps >= static_cast<long>(4 << 20));
  }
  CHECK(timer.HeapCoverage(Timer::IO) > 0.2);
  CHECK(timer.HeapCoverage(Timer::IO) < 0.8);
  CHECK(timer.HeapCoverage(Timer::EPS) == 1.0);
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

TEST_CASE("PerfProbes Hash and Verification", "[perfprobes][Serial]")
{
  const char a[] = "abc";
  CHECK(perf::HashBytes(a, 3) == 0xe71fa2190541574bull);  // FNV-1a 64-bit of "abc"
  CHECK(perf::HashBytes(a, 0) == 14695981039346656037ull);
  perf::SetVerification(1);
  const bool env = std::getenv("PALACE_PERF_VERIFY") != nullptr;
  CHECK(perf::Verification() == env);
  perf::SetVerification(2);
  CHECK(perf::Verification());
  perf::SetVerification(0);
}

TEST_CASE("PerfProbes Residual", "[perfprobes][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::QUADRILATERAL);
  mfem::ParMesh pmesh(comm, smesh);
  mfem::H1_FECollection fec(1, 2);
  mfem::ParFiniteElementSpace fes(&pmesh, &fec);
  mfem::ParBilinearForm m(&fes);
  m.AddDomainIntegrator(new mfem::MassIntegrator);
  m.Assemble();
  m.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> M(m.ParallelAssemble());
  mfem::Vector x(M->Height()), b(M->Height());
  x.Randomize(1);
  M->Mult(x, b);
  CHECK_THAT(perf::RelativeResidual(*M, b, x, comm), WithinAbs(0.0, 1.0e-14));
  x *= 2.0;
  CHECK_THAT(perf::RelativeResidual(*M, b, x, comm), WithinAbs(1.0, 1.0e-12));

  std::vector<long> tdofs = {fes.GetTrueVSize()};
  perf::RecordPartitionHash(pmesh, tdofs);  // Collective, must not hang
}
