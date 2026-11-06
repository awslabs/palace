// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstdlib>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <mfem.hpp>
#include "utils/communication.hpp"
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
  CHECK_THAT(perf::RelativeResidual(*M, b, x), WithinAbs(0.0, 1.0e-14));
  x *= 2.0;
  CHECK_THAT(perf::RelativeResidual(*M, b, x), WithinAbs(1.0, 1.0e-12));

  std::vector<long> tdofs = {fes.GetTrueVSize()};
  perf::RecordPartitionHash(pmesh, tdofs);  // Collective, must not hang
}
