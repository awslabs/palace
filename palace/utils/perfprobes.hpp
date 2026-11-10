// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_PERFPROBES_HPP
#define PALACE_UTILS_PERFPROBES_HPP

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <mpi.h>

namespace mfem
{
class HypreParMatrix;
class ParMesh;
class Solver;
class Vector;
}  // namespace mfem

namespace palace::perf
{

//
// Probes feeding the performance report: memory attribution of large consumers, resolved
// solver decisions, assembly tallies and opt-in verification checks.
//

// Opt-in verification checks, enabled by Problem.Verbose >= 2 or PALACE_PERF_VERIFY=1.
// They may synchronize and cost a little, but do not change any numerics.
void SetVerification(int verbose);
bool Verification();

// RSS, peak RSS (high-water mark) and heap in use over one or more segments, recorded
// under Memory/<name>/: PersistentMegabytes (RSS growth), PeakGrowthMegabytes (high-water
// mark growth), TransientMegabytes (in-segment peak above the final RSS, when the segment
// set a new peak; otherwise zero) and HeapPersistentMegabytes (heap in use growth, which
// also counts allocations served from freed but still resident heap), each the largest
// over the events of the iteration, and Events.
class MemoryProbe
{
  long cur0 = 0, hwm0 = 0, heap0 = 0;
  double persistent = 0.0, growth = 0.0, transient = 0.0, heap = 0.0;
  bool active = false, used = false;

public:
  void Start();
  void Stop();
  void Record(std::string_view name);
  bool Used() const { return used; }
};

// Scoped single-segment MemoryProbe.
class ScopedMemoryProbe
{
  MemoryProbe probe;
  std::string name;

public:
  explicit ScopedMemoryProbe(std::string name) : name(std::move(name)) { probe.Start(); }
  ~ScopedMemoryProbe()
  {
    probe.Stop();
    probe.Record(name);
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

// Fill x with uniform random entries in [-1, 1], seeded by the rank and a stream index,
// without changing the state of the C library generator (unlike mfem::Vector::Randomize).
void SetRandom(mfem::Vector &x, MPI_Comm comm, int stream = 0);

// Relative residual ‖b - A x‖ / ‖b‖ (collective).
double RelativeResidual(const mfem::HypreParMatrix &A, const mfem::Vector &b,
                        const mfem::Vector &x);

// Rank-order-sensitive hashes of the mesh partition and of the true DOF counts per rank
// (collective), written as Info under Verification/.
void RecordPartitionHash(const mfem::ParMesh &mesh, const std::vector<long> &tdofs);

// FNV-1a over raw bytes, exposed for testing.
std::uint64_t HashBytes(const void *data, std::size_t n,
                        std::uint64_t h = 14695981039346656037ull);

}  // namespace palace::perf

#endif  // PALACE_UTILS_PERFPROBES_HPP
