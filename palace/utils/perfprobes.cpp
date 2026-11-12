// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "perfprobes.hpp"

#include <algorithm>
#include <array>
#include <_hypre_parcsr_ls.h>
#include <fmt/format.h>
#include <mfem.hpp>
#include "linalg/ams.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/perfstats.hpp"

namespace palace::perf
{

namespace
{

constexpr double to_mb = 1.0 / (1024.0 * 1024.0);

void RecordAmg(HYPRE_Solver solver, const std::string &base)
{
  auto *amg = reinterpret_cast<hypre_ParAMGData *>(solver);
  if (!amg || !hypre_ParAMGDataAArray(amg))
  {
    return;
  }
  const int num_levels = hypre_ParAMGDataNumLevels(amg);
  Set(base + "/Levels", num_levels, Reduce::MAX);
  for (int l = 0; l < num_levels; l++)
  {
    hypre_ParCSRMatrix *A = hypre_ParAMGDataAArray(amg)[l];
    if (!A)
    {
      continue;
    }
    const auto rows = static_cast<double>(hypre_ParCSRMatrixNumRows(A));
    Set(fmt::format("{}/RowsPerRank/L{}", base, l), rows, Reduce::STATS);
    if (l == num_levels - 1)
    {
      Set(base + "/CoarseRowsPerRank", rows, Reduce::STATS);
    }
    if (const auto *pkg = hypre_ParCSRMatrixCommPkg(A))
    {
      Set(fmt::format("{}/Neighbours/L{}", base, l),
          std::max(hypre_ParCSRCommPkgNumSends(pkg), hypre_ParCSRCommPkgNumRecvs(pkg)),
          Reduce::STATS);
    }
  }
}

}  // namespace

void MemoryProbe::Start()
{
  cur0 = memory_reporting::GetCurrentMemory();
  hwm0 = memory_reporting::GetPeakMemory();
  heap0 = memory_reporting::GetHeapMemory();
  BlockTimer::NoteRss(cur0);
  active = true;
}

void MemoryProbe::Stop()
{
  if (!active)
  {
    return;
  }
  const long heap1 = memory_reporting::GetHeapMemory();
  const long cur1 = memory_reporting::GetCurrentMemory();
  BlockTimer::NoteRss(cur1);
  const long hwm1 = memory_reporting::GetPeakMemory();
  persistent += static_cast<double>(cur1 - cur0);
  growth += static_cast<double>(hwm1 - hwm0);
  heap += static_cast<double>(heap1 - heap0);
  if (hwm1 > hwm0)
  {
    transient = std::max(transient, static_cast<double>(hwm1 - cur1));
  }
  active = false;
  used = true;
}

MemoryProbeCounters::MemoryProbeCounters(std::string_view name)
{
  const std::string base = "Memory/" + std::string(name);
  auto Max = [&base](const char *field)
  { return RegisterCounter(base + field, Reduce::STATS, Accumulate::MAX); };
  persistent = Max("/PersistentMegabytes");
  growth = Max("/PeakGrowthMegabytes");
  transient = Max("/TransientMegabytes");
  heap = Max("/HeapPersistentMegabytes");
  events = RegisterCounter(base + "/Events", Reduce::MAX, Accumulate::SUM);
}

void MemoryProbe::Record(const MemoryProbeCounters &c)
{
  if (!used)
  {
    return;
  }
  SetMax(c.persistent, persistent * to_mb);
  SetMax(c.growth, growth * to_mb);
  SetMax(c.transient, transient * to_mb);
  SetMax(c.heap, heap * to_mb);
  Add(c.events, 1.0);
  persistent = growth = transient = heap = 0.0;
  used = false;
}

bool MemoryProbeSite::Sample()
{
  if (iteration == detail::iteration_serial)
  {
    return false;
  }
  iteration = detail::iteration_serial;
  return true;
}

int AssemblyTally(AssemblyKind kind)
{
  static std::array<std::array<int, Timer::NUM_TIMINGS>, 3> ids = []()
  {
    std::array<std::array<int, Timer::NUM_TIMINGS>, 3> a;
    for (auto &row : a)
    {
      row.fill(-2);
    }
    return a;
  }();
  const auto idx = BlockTimer::Current();
  int &id = ids[static_cast<int>(kind)][idx];
  if (id == -2)
  {
    constexpr std::array<const char *, 3> kinds = {"Partial", "Full", "Parallel"};
    id = RegisterTally(
        fmt::format("Assembly/{}/{}", kinds[static_cast<int>(kind)], Timer::names[idx]));
  }
  return id;
}

void RecordHypreHierarchy(const mfem::Solver &solver, std::string_view role)
{
  if (const auto *amg = dynamic_cast<const mfem::HypreBoomerAMG *>(&solver))
  {
    RecordAmg(static_cast<HYPRE_Solver>(*amg), LinearSolverPath(role, "AMG"));
  }
  else if (const auto *ams = dynamic_cast<const HypreAmsSolver *>(&solver))
  {
    auto *data = reinterpret_cast<hypre_AMSData *>(static_cast<HYPRE_Solver>(*ams));
    if (!data)
    {
      return;
    }
    const std::string base = LinearSolverPath(role, "AMS");
    RecordAmg(hypre_AMSDataPoissonBetaAMG(data), base + "/G");
    RecordAmg(hypre_AMSDataPoissonAlphaAMG(data) ? hypre_AMSDataPoissonAlphaAMG(data)
                                                 : hypre_AMSDataPoissonAlphaXAMG(data),
              base + "/Pi");
  }
}

}  // namespace palace::perf
