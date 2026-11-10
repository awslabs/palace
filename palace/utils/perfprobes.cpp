// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "perfprobes.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <random>
#include <_hypre_parcsr_ls.h>
#include <fmt/format.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "linalg/ams.hpp"
#include "utils/communication.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/perfstats.hpp"

namespace palace::perf
{

namespace
{

bool verification = false;

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

void SetVerification(int verbose)
{
  const char *env = std::getenv("PALACE_PERF_VERIFY");
  verification = (verbose >= 2) || (env && std::string_view(env) != "0");
}

bool Verification()
{
  return verification;
}

void MemoryProbe::Start()
{
  cur0 = memory_reporting::GetCurrentMemory();
  hwm0 = memory_reporting::GetPeakMemory();
  heap0 = BlockTimer::QueryHeap();
  active = true;
}

void MemoryProbe::Stop()
{
  if (!active)
  {
    return;
  }
  const long heap1 = BlockTimer::QueryHeap();
  const long cur1 = memory_reporting::GetCurrentMemory();
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

void MemoryProbe::Record(std::string_view name)
{
  if (!used)
  {
    return;
  }
  const std::string base = "Memory/" + std::string(name);
  auto Max = [&](const char *field, double v)
  { SetMax(RegisterCounter(base + field, Reduce::STATS, Accumulate::MAX), v * to_mb); };
  Max("/PersistentMegabytes", persistent);
  Max("/PeakGrowthMegabytes", growth);
  Max("/TransientMegabytes", transient);
  Max("/HeapPersistentMegabytes", heap);
  Add(RegisterCounter(base + "/Events", Reduce::MAX, Accumulate::SUM), 1.0);
  persistent = growth = transient = heap = 0.0;
  used = false;
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
        fmt::format("Assembly/{}/{}", kinds[static_cast<int>(kind)], Timer::names[idx]),
        idx);
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

void SetRandom(mfem::Vector &x, MPI_Comm comm, int stream)
{
  std::seed_seq seq{Mpi::Rank(comm), stream};
  std::mt19937_64 gen(seq);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  double *X = x.HostWrite();
  for (int i = 0; i < x.Size(); i++)
  {
    X[i] = dist(gen);
  }
}

double RelativeResidual(const mfem::HypreParMatrix &A, const mfem::Vector &b,
                        const mfem::Vector &x)
{
  mfem::Vector r(b.Size());
  r.UseDevice(true);
  A.Mult(x, r);
  r -= b;
  MPI_Comm comm = A.GetComm();
  const double nr = std::sqrt(mfem::InnerProduct(comm, r, r));
  const double nb = std::sqrt(mfem::InnerProduct(comm, b, b));
  return (nb > 0.0) ? nr / nb : nr;
}

std::uint64_t HashBytes(const void *data, std::size_t n, std::uint64_t h)
{
  const auto *p = static_cast<const unsigned char *>(data);
  for (std::size_t i = 0; i < n; i++)
  {
    h = (h ^ p[i]) * 1099511628211ull;
  }
  return h;
}

void RecordPartitionHash(const mfem::ParMesh &mesh, const std::vector<long> &tdofs)
{
  MPI_Comm comm = mesh.GetComm();
  const int rank = Mpi::Rank(comm);
  // Per-rank hashes of the local elements (vertex coordinates, attributes) and of the
  // true DOF counts, combined across ranks in rank order.
  std::uint64_t h[2] = {HashBytes(&rank, sizeof(rank)), HashBytes(&rank, sizeof(rank))};
  mfem::Array<int> v;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const int attr = mesh.GetAttribute(e);
    h[0] = HashBytes(&attr, sizeof(attr), h[0]);
    mesh.GetElementVertices(e, v);
    for (int i : v)
    {
      h[0] = HashBytes(mesh.GetVertex(i), sizeof(double) * mesh.SpaceDimension(), h[0]);
    }
  }
  h[1] = HashBytes(tdofs.data(), sizeof(long) * tdofs.size(), h[1]);
  std::vector<std::uint64_t> all(2 * Mpi::Size(comm));
  MPI_Gather(h, 2, MPI_UINT64_T, all.data(), 2, MPI_UINT64_T, 0, comm);
  if (Mpi::Root(comm))
  {
    std::uint64_t g[2] = {14695981039346656037ull, 14695981039346656037ull};
    for (std::size_t i = 0; i < all.size(); i++)
    {
      g[i % 2] = HashBytes(&all[i], sizeof(std::uint64_t), g[i % 2]);
    }
    SetInfo("Verification/PartitionHash", fmt::format("{:016x}", g[0]));
    SetInfo("Verification/DofHash", fmt::format("{:016x}", g[1]));
  }
}

}  // namespace palace::perf
