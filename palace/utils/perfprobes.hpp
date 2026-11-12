// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_PERFPROBES_HPP
#define PALACE_UTILS_PERFPROBES_HPP

#include <string>
#include <string_view>

namespace mfem
{
class Solver;
}  // namespace mfem

namespace palace::perf
{

//
// Probes feeding the performance report: memory attribution of large consumers, resolved
// solver decisions and assembly tallies.
//

// Counter ids of a memory probe name (see MemoryProbe), registered once.
struct MemoryProbeCounters
{
  int persistent, growth, transient, heap, events;
  explicit MemoryProbeCounters(std::string_view name);
};

// RSS, peak RSS (high-water mark) and heap in use over one or more segments, recorded
// under Memory/<name>/: PersistentMegabytes (RSS growth), PeakGrowthMegabytes (high-water
// mark growth), TransientMegabytes (in-segment peak above the final RSS, when the segment
// set a new peak; otherwise zero) and HeapPersistentMegabytes (heap in use growth, which
// also counts allocations served from freed but still resident heap), each the largest
// over the events of the iteration, and Events. Each segment queries the RSS, the peak RSS
// and the heap twice: avoid probes in hot paths (see MemoryProbeSite).
class MemoryProbe
{
  long cur0 = 0, hwm0 = 0, heap0 = 0;
  double persistent = 0.0, growth = 0.0, transient = 0.0, heap = 0.0;
  bool active = false, used = false;

public:
  void Start();
  void Stop();
  void Record(std::string_view name) { Record(MemoryProbeCounters(name)); }
  void Record(const MemoryProbeCounters &counters);
  bool Used() const { return used; }
};

// A probe call site in a loop (for example the error estimate of each time step), measured
// in the first of its scopes of each reporting iteration only. Declare it static.
class MemoryProbeSite
{
  MemoryProbeCounters counters;
  long iteration = -1;  // perf::detail::iteration_serial of the last scope measured

public:
  explicit MemoryProbeSite(std::string_view name) : counters(name) {}
  // Whether to measure the current scope: true for the first one of the iteration.
  bool Sample();
  const MemoryProbeCounters &Counters() const { return counters; }
};

// Scoped single-segment MemoryProbe, of every scope or of the first scope per iteration
// of a call site.
class ScopedMemoryProbe
{
  MemoryProbe probe;
  std::string name;
  const MemoryProbeCounters *counters = nullptr;
  bool active = true;

public:
  explicit ScopedMemoryProbe(std::string name) : name(std::move(name)) { probe.Start(); }
  explicit ScopedMemoryProbe(MemoryProbeSite &site)
    : counters(&site.Counters()), active(site.Sample())
  {
    if (active)
    {
      probe.Start();
    }
  }
  ~ScopedMemoryProbe()
  {
    if (active)
    {
      probe.Stop();
      counters ? probe.Record(*counters) : probe.Record(name);
    }
  }
  ScopedMemoryProbe(const ScopedMemoryProbe &) = delete;
  ScopedMemoryProbe &operator=(const ScopedMemoryProbe &) = delete;
};

// Tally id for libCEED partial assembly, full (sparse) assembly or parallel (RAP) assembly
// within the current exclusive timer row: Assembly/<Kind>/<ROW>.
enum class AssemblyKind
{
  PARTIAL,
  FULL,
  PARALLEL
};
int AssemblyTally(AssemblyKind kind);

// Hierarchy summary of a hypre BoomerAMG or AMS solver after setup, written as counters
// under LinearSolver/<Role>/AMG or LinearSolver/<Role>/AMS/{G,Pi}.
void RecordHypreHierarchy(const mfem::Solver &solver, std::string_view role);

}  // namespace palace::perf

#endif  // PALACE_UTILS_PERFPROBES_HPP
