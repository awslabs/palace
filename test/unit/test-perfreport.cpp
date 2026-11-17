// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <memory>
#include <thread>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/ksp.hpp"
#include "linalg/solver.hpp"
#include "utils/communication.hpp"
#include "utils/perfreport.hpp"
#include "utils/timer.hpp"

using namespace palace;
using Catch::Matchers::WithinAbs;

namespace
{
void Sleep(int ms)
{
  std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

class SpdTestOperator : public Operator
{
public:
  SpdTestOperator() : Operator(2) {}
  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(2);
    diag[0] = 4.0;
    diag[1] = 3.0;
  }
  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(2);
    y[0] = 4.0 * x[0] + x[1];
    y[1] = x[0] + 3.0 * x[1];
  }
};
}  // namespace

TEST_CASE("Timer Iterations", "[perfreport][Serial]")
{
  BlockTimer::Reset();
  Sleep(5);
  {
    BlockTimer bt(Timer::IO);
    Sleep(5);
  }
  const auto &t = BlockTimer::GlobalTimer();
  CHECK(t.Counts(Timer::STARTUP) == 1);
  CHECK(t.Data(Timer::STARTUP) >= 0.004);
  // A row open at the report is accounted up to the report, and not twice.
  {
    BlockTimer bt(Timer::POSTPRO);
    Sleep(5);
    BlockTimer::Snapshot();
    const double t1 = t.Data(Timer::POSTPRO), total1 = t.Data(Timer::TOTAL);
    CHECK(t1 >= 0.004);
    CHECK(t.Counts(Timer::POSTPRO) == 0);
    BlockTimer::StartIteration();
    Sleep(5);
  }
  BlockTimer::Snapshot();
  CHECK(t.Iterations() == 2);
  CHECK(t.Counts(Timer::POSTPRO) == 1);
  CHECK(t.Data(Timer::POSTPRO) >= 0.004);
  CHECK(t.Counts(Timer::STARTUP) == 0);
  CHECK(t.RunCounts(Timer::STARTUP) == 1);
  CHECK_THAT(t.RunData(Timer::POSTPRO) - t.Data(Timer::POSTPRO), WithinAbs(0.005, 0.004));
  CHECK(t.RunCounts(Timer::TOTAL) == 2);
  BlockTimer::Reset();
}

TEST_CASE("PerfReport Iterations And Cumulative", "[perfreport][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  BlockTimer::Reset();
  perf::SumMetric its("Test/Report/Its", perf::Red::COUNT);
  perf::LastMetric size("Test/Report/Size", perf::Red::STATS, perf::Run::ITERATION_ONLY);
  perf::SumMetric once("Test/Report/Once", perf::Red::MAX);
  perf::MaxMetric prune("LinearSolver/TestPrune/X", perf::Red::MAX);
  {
    BlockTimer bt(Timer::IO);
    Sleep(5);
  }
  its.Add(3);
  size.Set(100);
  once.Add(1);
  prune.Update(1);
  const auto r1 = perf::Report(comm);  // Copy
  nlohmann::json meta = {{"Problem", {{"Iteration", 1}}}, {"LinearSolver", {{"Other", 1}}}};
  perf::WriteJson(meta);
  CHECK(meta["MetadataVersion"] == perf::METADATA_VERSION);
  CHECK(meta["Problem"]["Iteration"] == 1);  // Other writers' keys are kept
  for (const char *k : {"Min", "Avg", "Max", "Sum", "MaxRank", "Ranks", "Count"})
  {
    CHECK(meta["ElapsedTime"]["IO"].contains(k));
    CHECK(meta["ElapsedTime"]["TOTAL"].contains(k));
  }
  CHECK(!meta["ElapsedTime"].contains("KSP"));  // Unused rows are omitted
  CHECK(meta["Test"]["Report"]["Once"] == 1);
  CHECK(meta["Cumulative"]["Test"]["Report"]["Its"] == meta["Test"]["Report"]["Its"]);
  CHECK(!meta["Cumulative"]["Test"]["Report"].contains("Size"));
  CHECK(meta["PeakMemoryGrowthMegabytes"].contains("Baseline"));
  CHECK(meta["LinearSolver"]["TestPrune"]["X"] == 1);
  CHECK(meta["PeakMemoryGrowthMegabytes"]["Baseline"] ==
        meta["Cumulative"]["PeakMemoryGrowthMegabytes"]["Baseline"]);  // Iteration 1
  CHECK(meta["Cumulative"]["ElapsedTime"]["IO"] == meta["ElapsedTime"]["IO"]);

  BlockTimer::StartIteration();
  {
    BlockTimer bt(Timer::IO);
    Sleep(5);
  }
  // Growth which stays alive until after the report.
  constexpr std::size_t alloc = 64 << 20;
  std::vector<char> buf;
  long grew;
  {
    BlockTimer bt(Timer::CONSTRUCT);
    const long before = memory_reporting::GetPeakMemory();
    buf.assign(alloc, 1);
    grew = memory_reporting::GetPeakMemory() - before;
  }
  Mpi::GlobalMin(1, &grew, comm);
  its.Add(4);
  size.Set(200);
  const auto &r2 = perf::Report(comm);
  perf::WriteJson(meta);
  CHECK(r2.iteration == 2);
  CHECK(meta["Test"]["Report"]["Its"] == 4);
  CHECK(meta["Cumulative"]["Test"]["Report"]["Its"] == 7);
  CHECK(!meta["Test"]["Report"].contains("Once"));  // Stale leaf removed
  CHECK(meta["Cumulative"]["Test"]["Report"]["Once"] == 1);
  CHECK(meta["Test"]["Report"]["Size"]["Max"] == 200);
  // Sums of times over the iteration reports add up to Cumulative.
  for (int i : {static_cast<int>(Timer::IO), static_cast<int>(Timer::TOTAL)})
  {
    CHECK_THAT(r2.cum.time[i].sum, WithinAbs(r1.iter.time[i].sum + r2.iter.time[i].sum,
                                             1e-9 * r2.cum.time[i].sum));
  }
  CHECK(!meta["LinearSolver"].contains("TestPrune"));  // Stale subtree removed
  CHECK(meta["LinearSolver"]["Other"] == 1);
  CHECK(buf[alloc / 2] == 1);
  if (grew >= static_cast<long>(alloc / 2))
  {
    CHECK(r2.iter.mem[Timer::CONSTRUCT].max > alloc / 2);
    CHECK(r2.iter.unattributed.min >= 0.0);
  }
  else
  {
    WARN("Peak memory did not grow with the allocation, skipping growth checks");
  }
  BlockTimer::Reset();
}

TEST_CASE("PerfReport Subcommunicator Role", "[perfreport][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm);
  if (np < 2)
  {
    SKIP("Needs at least 2 processes");
  }
  BlockTimer::Reset();
  // Two "ports" on subcommunicators {even ranks}, {odd ranks}: each does 5 solves.
  MPI_Comm sub;
  MPI_Comm_split(comm, rank % 2, rank, &sub);
  {
    perf::ScopedRole role("TestPort", sub);
    perf::SumMetric solves("Test/Port/Solves", perf::Red::COUNT);
    perf::MaxMetric max_its("Test/Port/MaxIts", perf::Red::MAX);
    solves.Add(5);
    max_its.Update(10 + rank % 2);
  }
  MPI_Comm_free(&sub);
  const auto &r = perf::Report(comm);
  CHECK(r.Metric("Test/Port/Solves", false)->sum == 10);  // Totals over ports
  CHECK(r.Metric("Test/Port/MaxIts", false)->max == 11);
  BlockTimer::Reset();
}

TEST_CASE("PerfReport Baseline After Reset", "[perfreport][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  {
    std::vector<char> buf(64 << 20, 1);
    CHECK(buf[buf.size() / 2] == 1);
  }
  BlockTimer::Reset();
  {
    BlockTimer bt(Timer::IO);
  }
  const auto &r = perf::Report(comm);
  nlohmann::json meta = nlohmann::json::object();
  perf::WriteJson(meta);
  CHECK(meta["PeakMemoryGrowthMegabytes"]["Baseline"]["Min"].get<double>() > 0.0);
  CHECK(meta["Cumulative"]["PeakMemoryGrowthMegabytes"]["Baseline"]["Min"].get<double>() >
        0.0);
  CHECK(meta["PeakNodeMemoryMegabytes"]["Nodes"].get<int>() >= 1);
  if (meta["PeakNodeMemoryMegabytes"]["Nodes"] == 1)
  {
    CHECK(r.iter.node_peak.sum == r.iter.peak.sum);
  }
  BlockTimer::Reset();
}

TEST_CASE("PerfReport Row On One Process", "[perfreport][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  if (Mpi::Size(comm) < 2)
  {
    SKIP("Needs at least 2 processes");
  }
  BlockTimer::Reset();
  if (Mpi::Root(comm))
  {
    BlockTimer bt(Timer::EPS);
    Sleep(2);
  }
  const auto &r = perf::Report(comm);
  CHECK(r.iter.time[Timer::EPS].ranks == 1);
  CHECK(r.iter.time[Timer::TOTAL].ranks == Mpi::Size(comm));
  BlockTimer::Reset();
}

TEST_CASE("PerfReport Linear Solver Statistics", "[perfreport][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  BlockTimer::Reset();
  SpdTestOperator A;
  auto make = [&](int max_it)
  {
    auto cg = std::make_unique<CgSolver<Operator>>(comm, 0);
    cg->SetRelTol(1.0e-12);
    cg->SetAbsTol(0.0);
    cg->SetMaxIter(max_it);
    return std::make_unique<KspSolver>(std::move(cg),
                                       std::make_unique<JacobiSmoother<Operator>>(comm));
  };
  std::unique_ptr<KspSolver> converging, truncated;
  {
    perf::ScopedRole role("TestMain");
    converging = make(10);
    truncated = make(1);  // One iteration cannot solve the 2x2 system
  }
  {
    // Updates go to the role captured at construction.
    perf::ScopedRole role("TestOther");
    Vector b(2), x(2);
    b = 0.0;
    b[0] = 1.0;
    for (auto *ksp : {converging.get(), truncated.get()})
    {
      x = 0.0;
      ksp->SetOperators(A, A);
      ksp->Mult(b, x);
    }
  }
  const auto &r = perf::Report(comm);
  auto max = [&](const char *field)
  {
    const auto *s = r.Metric(std::string("LinearSolver/TestMain/") + field, false);
    REQUIRE(s);
    return s->max;
  };
  CHECK(max("Solves") == 2);
  CHECK(max("Setups") == 2);
  CHECK(max("Its") >= 2);
  CHECK(max("MaxIts") >= 1);
  CHECK(max("NotConverged") == 1);
  CHECK(!r.Metric("LinearSolver/TestOther/Solves", false));
  BlockTimer::Reset();
}

TEST_CASE("PerfReport Lazy Coarse Setup", "[perfreport][Serial]")
{
  // 1D Laplacian on one process; BoomerAMG sets up at the first apply.
  constexpr int n = 10;
  mfem::SparseMatrix S(n, n);
  for (int i = 0; i < n; i++)
  {
    S.Add(i, i, 2.0);
    if (i > 0)
    {
      S.Add(i, i - 1, -1.0);
    }
    if (i < n - 1)
    {
      S.Add(i, i + 1, -1.0);
    }
  }
  S.Finalize();
  HYPRE_BigInt row_starts[2] = {0, n};
  mfem::HypreParMatrix A(Mpi::World(), n, row_starts, &S);

  BlockTimer::Reset();
  MfemWrapperSolver<Operator> pc(std::make_unique<mfem::HypreBoomerAMG>());
  pc.EnableSetupTimer();
  Vector b(n), x(n);
  b = 1.0;
  pc.SetOperator(A);
  CHECK(BlockTimer::GlobalTimer().Counts(Timer::KSP_COARSE_SETUP) == 0);
  pc.Mult(b, x);
  pc.Mult(b, x);
  CHECK(BlockTimer::GlobalTimer().Counts(Timer::KSP_COARSE_SETUP) == 1);
  pc.SetOperator(A);
  pc.Mult(b, x);
  CHECK(BlockTimer::GlobalTimer().Counts(Timer::KSP_COARSE_SETUP) == 2);
  BlockTimer::Reset();
}
