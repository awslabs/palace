// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <thread>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/ksp.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"
#include "utils/perfprobes.hpp"
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

// Exact inverse of a diagonal matrix: a direct coarse solver for the verification checks.
class DiagonalSolver : public mfem::Solver
{
  mfem::Vector d;

public:
  void SetOperator(const mfem::Operator &op) override
  {
    const auto &A = dynamic_cast<const mfem::HypreParMatrix &>(op);
    height = width = A.Height();
    A.GetDiag(d);
  }
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    const double *X = x.HostRead(), *D = d.HostRead();
    double *Y = y.HostWrite();
    for (int i = 0; i < x.Size(); i++)
    {
      Y[i] = X[i] / D[i];
    }
  }
};

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

TEST_CASE("PerfStats Non-Finite Values", "[perfstats][Serial][Parallel]")
{
  // Maxima and minima keep a value which is not a number, also over processes and in the
  // Cumulative block.
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm);
  BlockTimer::Reset();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const int mx =
      perf::RegisterCounter("Test/NaNMax", perf::Reduce::MAX, perf::Accumulate::MAX);
  const int mn =
      perf::RegisterCounter("Test/NaNMin", perf::Reduce::MAX, perf::Accumulate::MIN);
  const int rk =
      perf::RegisterCounter("Test/NaNRoot", perf::Reduce::STATS, perf::Accumulate::MAX);
  perf::SetMax(mx, 1.0);
  perf::SetMax(mx, nan);
  perf::SetMax(mx, 2.0);
  perf::SetMin(mn, 1.0);
  perf::SetMin(mn, nan);
  perf::SetMin(mn, 0.5);
  perf::SetMax(rk, (rank == 0) ? nan : 1.0);
  BlockTimer::Finalize(comm);
  {
    const auto &r = perf::GetResults();
    CHECK(std::isnan(r.counter[CounterIndex(r, "Test/NaNMax")].max));
    CHECK(std::isnan(r.counter[CounterIndex(r, "Test/NaNMin")].max));
    CHECK(std::isnan(r.counter[CounterIndex(r, "Test/NaNRoot")].max));
    CHECK(std::isnan(r.counter[CounterIndex(r, "Test/NaNRoot")].min));
  }
  BlockTimer::StartIteration();
  perf::SetMax(mx, 3.0);
  perf::SetMin(mn, 0.1);
  BlockTimer::Finalize(comm);
  {
    const auto &r = perf::GetResults();
    CHECK(r.counter[CounterIndex(r, "Test/NaNMax")].max == 3.0);
    CHECK(std::isnan(r.cum_counter[CounterIndex(r, "Test/NaNMax")].max));
    CHECK(std::isnan(r.cum_counter[CounterIndex(r, "Test/NaNMin")].max));
  }
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

TEST_CASE("PerfStats Stale Counters Of Other Ranks", "[perfstats][Parallel]")
{
  // A counter registered on another process only is erased by the root as well, once it
  // is no longer set.
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm);
  BlockTimer::Reset();
  if (rank == np - 1)
  {
    perf::Set(perf::RegisterCounter("Test/LastRankOnly", perf::Reduce::MAX), 3.0);
  }
  BlockTimer::Finalize(comm);
  nlohmann::json meta = nlohmann::json::object();
  perf::WriteJson(meta);
  CHECK(meta["Test"]["LastRankOnly"] == 3);
  BlockTimer::StartIteration();
  BlockTimer::Finalize(comm);
  perf::WriteJson(meta);
  CHECK(!meta.contains("Test"));
  CHECK(meta["Cumulative"]["Test"]["LastRankOnly"] == 3);
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

TEST_CASE("PerfStats Coarse Setup And Apply", "[perfstats][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm);
  BlockTimer::Reset();
  constexpr int n = 20;
  mfem::SparseMatrix diag(n);
  for (int i = 0; i < n; i++)
  {
    diag.Set(i, i, 2.0);
  }
  diag.Finalize();
  HYPRE_BigInt row_starts[2] = {static_cast<HYPRE_BigInt>(rank) * n,
                                static_cast<HYPRE_BigInt>(rank + 1) * n};
  mfem::HypreParMatrix A(comm, static_cast<HYPRE_BigInt>(np) * n, row_starts, &diag);
  auto amg = std::make_unique<mfem::HypreBoomerAMG>();
  amg->SetPrintLevel(0);
  MfemWrapperSolver<Operator> pc(std::move(amg));
  pc.EnableTimer(true, Timer::KSP_COARSE_SOLVE);
  Vector x(n), y(n);
  x = 1.0;
  {
    BlockTimer bt(Timer::KSP_SETUP);
    pc.SetOperator(A);
  }
  for (int k = 0; k < 3; k++)
  {
    BlockTimer bt(Timer::KSP_PRECONDITIONER);
    pc.Mult(x, y);
  }
  pc.SetOperator(A);
  pc.Mult(x, y);
  BlockTimer::Finalize(comm);
  const auto &r = perf::GetResults();
  // The first application after each SetOperator is the (lazy) setup.
  CHECK(r.iter.count[Timer::KSP_COARSE_SETUP].max == 2);
  CHECK(r.iter.count[Timer::KSP_COARSE_SOLVE].max == 2);
  CHECK(r.iter.count[Timer::KSP_PRECONDITIONER].max == 3);
  CHECK(r.counter[CounterIndex(r, "LinearSolver/Main/Coarse/Setups")].max == 2);
  CHECK(r.counter[CounterIndex(r, "LinearSolver/Main/Coarse/Applies")].max == 4);
  CHECK(r.counter[CounterIndex(r, "LinearSolver/Main/Coarse/Refactorizations")].max == 0);
  CHECK(r.counter[CounterIndex(r, "LinearSolver/Main/Coarse/Rows")].sum == np * n);
  CHECK(r.counter[CounterIndex(r, "LinearSolver/Main/Coarse/RowsPerRank")].min == n);
  BlockTimer::Reset();
}

TEST_CASE("PerfStats Verification Of Coarse Solve And Preconditioner",
          "[perfstats][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm);
  constexpr int n = 20;
  HYPRE_BigInt row_starts[2] = {static_cast<HYPRE_BigInt>(rank) * n,
                                static_cast<HYPRE_BigInt>(rank + 1) * n};
  std::vector<std::unique_ptr<mfem::SparseMatrix>> diags;
  auto Identity = [&](double s)
  {
    diags.push_back(std::make_unique<mfem::SparseMatrix>(n));
    for (int i = 0; i < n; i++)
    {
      diags.back()->Set(i, i, s);
    }
    diags.back()->Finalize();
    return std::make_unique<mfem::HypreParMatrix>(comm, static_cast<HYPRE_BigInt>(np) * n,
                                                  row_starts, diags.back().get());
  };
  auto A1 = Identity(1.0), A2 = Identity(2.0), A4 = Identity(4.0);
  auto Amg = []()
  {
    auto amg = std::make_unique<mfem::HypreBoomerAMG>();
    amg->SetPrintLevel(0);
    return amg;
  };
  auto Direct = []() { return std::make_unique<DiagonalSolver>(); };
  auto Value = [](const std::string &path)
  {
    const auto &r = perf::GetResults();
    const int i = CounterIndex(r, path);
    return (i >= 0 && r.counter[i].n > 0) ? r.counter[i].max : -1.0;
  };
  const std::string growth = "Verification/Preconditioner/Main/MaxErrorPropagatorGrowth",
                    checks = "Verification/Preconditioner/Main/Checks",
                    forward = "Verification/CoarseSolve/Main/MaxForwardError";

  SECTION("Off by default")
  {
    perf::SetVerification(0);
    BlockTimer::Reset();
    KspSolver ksp(std::make_unique<CgSolver<Operator>>(comm, 0),
                  std::make_unique<MfemWrapperSolver<Operator>>(Direct()));
    Vector b(n), x(n);
    b = 1.0;
    ksp.SetOperators(*A2, *A4);
    ksp.Mult(b, x);
    BlockTimer::Finalize(comm);
    CHECK(Value(growth) == -1.0);
    CHECK(Value(forward) == -1.0);
    const auto &r = perf::GetResults();
    CHECK(r.iter.count[Timer::VERIFICATION].max == 0);
    CHECK(r.iter.time[Timer::VERIFICATION].max == 0.0);
  }

  SECTION("Real-valued")
  {
    // A = 2 I and M = (4 I)^-1 exactly: I - M A = I / 2.
    perf::SetVerification(2);
    BlockTimer::Reset();
    auto pc = std::make_unique<MfemWrapperSolver<Operator>>(Direct());
    pc->EnableTimer(true, Timer::KSP_COARSE_SOLVE);
    KspSolver ksp(std::make_unique<CgSolver<Operator>>(comm, 0), std::move(pc));
    Vector b(n), x(n);
    b = 0.0;
    ksp.SetOperators(*A2, *A4);
    ksp.Mult(b, x);  // No iteration: not checked yet
    BlockTimer::Finalize(comm);
    CHECK(Value(checks) == -1.0);
    b = 1.0;
    ksp.Mult(b, x);
    ksp.SetOperators(*A2, *A4);
    ksp.Mult(b, x);  // Checked after the first solve with iterations only
    BlockTimer::Finalize(comm);
    CHECK_THAT(Value(growth), WithinRel(0.5, 1.0e-10));
    CHECK(Value(checks) == 1.0);
    CHECK_THAT(Value(forward), WithinAbs(0.0, 1.0e-12));

    // The 16 preconditioner applications of the check and the coarse solve of the forward
    // error check are not counted as the solver's, only in their own tallies.
    const auto &r = perf::GetResults();
    const double applies = Value("LinearSolver/Main/Coarse/Applies");
    CHECK(applies > 0.0);
    CHECK(applies < 16.0);
    CHECK(r.iter.count[Timer::KSP_COARSE_SOLVE].max +
              r.iter.count[Timer::KSP_COARSE_SETUP].max ==
          applies);
    CHECK(r.iter.tally_count[TallyIndex(r.iter, "Main/Verification/ErrorPropagator")].max ==
          1);
    CHECK(r.iter.tally_count[TallyIndex(r.iter, "Main/Verification/CoarseForwardError")]
              .max == 2);
    CHECK(r.iter.tally_parents[TallyIndex(r.iter, "Main/Verification/ErrorPropagator")] ==
          Timer::VERIFICATION);

    // One verification row scope per check: the probe, and the forward error and the
    // residual checks after each of the two setups.
    CHECK(r.iter.count[Timer::VERIFICATION].max == 5);
    CHECK(r.iter.time[Timer::VERIFICATION].min > 0.0);
  }

  SECTION("No forward error for hypre coarse solvers")
  {
    perf::SetVerification(2);
    BlockTimer::Reset();
    KspSolver ksp(std::make_unique<CgSolver<Operator>>(comm, 0),
                  std::make_unique<MfemWrapperSolver<Operator>>(Amg()));
    Vector b(n), x(n);
    b = 1.0;
    ksp.SetOperators(*A2, *A4);
    ksp.Mult(b, x);
    BlockTimer::Finalize(comm);
    CHECK(Value(forward) == -1.0);
    CHECK_THAT(Value(growth), WithinRel(0.5, 1.0e-6));
  }

  SECTION("Complex-valued with a real-valued coarse solve")
  {
    // A = (2 + i) I and M = ((2 + 1) I)^-1, the real-valued approximation: the error
    // propagators of the coarse solve and of the preconditioner are (1 - i) I / 3.
    perf::SetVerification(2);
    BlockTimer::Reset();
    double total = 0.0, sum = 0.0;
    {
      BlockTimer bt(Timer::KSP);  // All of the time is in the exclusive rows
      ComplexWrapperOperator A(A2.get(), A1.get());
      ComplexKspSolver ksp(std::make_unique<GmresSolver<ComplexOperator>>(comm, 0),
                           std::make_unique<MfemWrapperSolver<ComplexOperator>>(
                               Direct(), true, /*complex_matrix*/ false));
      ComplexVector b(n), x(n);
      b = 1.0;
      ksp.SetOperators(A, A);
      ksp.Mult(b, x);
      BlockTimer::Finalize(comm);
      const auto &r = perf::GetResults();
      for (int i = 0; i < Timer::NUM_TIMINGS; i++)
      {
        if (i == Timer::TOTAL)
        {
          total = r.iter.time[i].avg;
        }
        else
        {
          sum += r.iter.time[i].avg;
        }
      }
    }
    CHECK_THAT(Value(growth), WithinRel(std::sqrt(2.0) / 3.0, 1.0e-10));
    CHECK_THAT(Value(forward), WithinRel(std::sqrt(2.0) / 3.0, 1.0e-10));
    CHECK(perf::GetResults().iter.time[Timer::VERIFICATION].min > 0.0);
    CHECK_THAT(sum, WithinAbs(total, 1.0e-4));
  }
  perf::SetVerification(0);
  BlockTimer::Reset();
}

TEST_CASE("PerfStats Phase Tallies", "[perfstats][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  BlockTimer::Reset();
  const int id = perf::RegisterPhaseTally("Test/Phase");
  perf::AddTally(id, 1.0);
  {
    BlockTimer bt(Timer::IO);
    perf::AddTally(id, 2.0);
  }
  BlockTimer::Finalize(comm);
  const auto &r = perf::GetResults();
  const int all = TallyIndex(r.iter, "Test/Phase");
  const int io = TallyIndex(r.iter, "Test/Phase/IO");
  const int none = TallyIndex(r.iter, "Test/Phase/STARTUP");
  REQUIRE(all >= 0);
  REQUIRE(io >= 0);
  REQUIRE(none >= 0);
  CHECK(r.iter.tally_parents[all] == Timer::TOTAL);
  CHECK(r.iter.tally_parents[io] == Timer::IO);
  CHECK(r.iter.tally_time[all].max == 3.0);
  CHECK(r.iter.tally_time[io].max == 2.0);
  CHECK(r.iter.tally_time[none].max == 1.0);
  BlockTimer::Reset();
}

TEST_CASE("PerfStats Wait And Work", "[perfstats][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm);
  BlockTimer::Reset();
  const int imb = perf::RegisterTally("Test/Imbalance", Timer::KSP);
  perf::AddTally(imb, 1.0 + rank);
  const int root_only = perf::RegisterTally("Test/RootOnly", Timer::IO);
  if (rank == 0)
  {
    perf::AddTally(root_only, 1.0);
  }

  // A reduction inside a scope is its wait.
  const int red = perf::RegisterTally("Test/Reduce", Timer::KSP);
  {
    perf::ScopedTally t(red);
    double v = 1.0;
    Mpi::GlobalSum(1, &v, comm);
  }
  const auto &reg = perf::detail::registry;
  const auto &tr = reg.tallies[red];
  CHECK(tr.wait > 0.0);
  CHECK(tr.wait <= tr.time);

  // So is a barrier, behind a slow rank.
  const int bar = perf::RegisterTally("Test/Barrier", Timer::KSP);
  {
    perf::ScopedTally t(bar);
    if (rank == 0)
    {
      Sleep(20);
    }
    Mpi::Barrier(comm);
  }
  if (rank > 0)
  {
    CHECK(reg.tallies[bar].wait > 0.01);
  }
  CHECK(reg.tallies[reg.tally_ids.at("Barrier")].count >= 1);

  // Rank 0 works while the others block in MPI_Wait for its message: the wait lands in
  // their exchange, the work in rank 0's.
  const int ex = perf::RegisterTally("Test/Exchange", Timer::KSP);
  const auto it = reg.tally_ids.find("Wait/MPI");
  const long waits0 = (it != reg.tally_ids.end()) ? reg.tallies[it->second].count : 0;
  {
    BlockTimer bt(Timer::KSP);
    perf::ScopedTally t(ex);
    int msg = 0;
    if (rank == 0)
    {
      Sleep(50);
      std::vector<MPI_Request> req(np - 1);
      for (int p = 1; p < np; p++)
      {
        MPI_Isend(&msg, 1, MPI_INT, p, 0, comm, &req[p - 1]);
      }
      MPI_Waitall(np - 1, req.data(), MPI_STATUSES_IGNORE);
    }
    else
    {
      MPI_Request req;
      MPI_Irecv(&msg, 1, MPI_INT, 0, 0, comm, &req);
      MPI_Wait(&req, MPI_STATUS_IGNORE);
    }
  }
  if (perf::detail::MpiWaitTally() && rank > 0)
  {
    const auto &te = reg.tallies[ex];
    CHECK(te.wait > 0.04);
    CHECK(te.time - te.wait < 0.04);
    const auto &tw = reg.tallies[reg.tally_ids.at("Wait/MPI")];
    CHECK(tw.count >= waits0 + 1);
    CHECK(reg.tallies[reg.tally_ids.at("Wait/MPI/KSP")].time >= te.wait - 1.0e-12);
  }
#if defined(PALACE_WITH_MPI_WAIT_TALLY)
  CHECK(perf::detail::MpiWaitTally());  // The test executable has the wrappers
#endif
  BlockTimer::Finalize(comm);
  nlohmann::json meta = nlohmann::json::object();
  perf::WriteJson(meta);
  const auto &tallies = meta["ElapsedTime"]["Tallies"];
  CHECK_THAT(tallies["Test/Imbalance"]["Imbalance"].get<double>(),
             WithinRel(2.0 * np / (np + 1.0), 1e-12));
  CHECK(!tallies["Test/Imbalance"].contains("Wait"));
  CHECK(tallies["Test/RootOnly"]["Imbalance"] == 1.0);  // Over the ranks which ran it
  CHECK(tallies["Test/Reduce"]["Wait"]["Max"].get<double>() > 0.0);
  CHECK(tallies["Test/Reduce"]["Work"].contains("Imbalance"));
  CHECK(meta["Communication"]["MPIWaitTally"] == perf::detail::MpiWaitTally());
  CHECK(!tallies["Reduce"].contains("Work"));  // Pure wait
  if (perf::detail::MpiWaitTally() && np > 1)
  {
    const auto &e = tallies["Test/Exchange"];
    CHECK(e["Wait"]["Max"].get<double>() > 0.04);
    CHECK(e["Work"]["Max"].get<double>() >= 0.049);
    CHECK(e["Work"]["MaxRank"] == 0);
    CHECK(e["Work"]["Imbalance"].get<double>() > 1.0);
    CHECK(tallies["Wait/MPI/KSP"]["Parent"] == "KSP");
  }
  BlockTimer::Reset();
}

TEST_CASE("PerfStats Pause Wait And Nesting", "[perfstats][Serial]")
{
  // Blocking waits and nested scopes in a paused scope are not counted either.
  BlockTimer::Reset();
  perf::detail::AddWait(0.0);  // Registers the communication tallies
  const auto &reg = perf::detail::registry;
  const int id = perf::RegisterTally("Test/PausedNesting", Timer::KSP);
  const auto &wait = reg.tallies[reg.tally_ids.at("Wait/MPI")];
  const long waits = wait.count;
  {
    perf::ScopedTally outer(id);
    perf::ScopedPause pause;
    perf::ScopedTally inner(id);
    perf::detail::AddWait(0.5);
  }
  CHECK(wait.count == waits);
  CHECK(reg.tallies[id].nested == 0);
  CHECK(reg.tallies[id].count == 1);
  CHECK(reg.tallies[id].wait >= 0.5 - 1.0e-12);  // The clock runs for the enclosing scope
  BlockTimer::Reset();
}

TEST_CASE("PerfStats Numeric Path Tokens", "[perfstats][Serial]")
{
  BlockTimer::Reset();
  perf::Set("Test/Levels/0/Rows", 10.0);
  perf::Set("Test/Levels/1/Rows", 20.0);
  perf::SetInfo("Test/Info/2", "x");
  BlockTimer::Finalize(Mpi::World());
  nlohmann::json meta = nlohmann::json::object();
  perf::WriteJson(meta);
  CHECK(meta["Test"]["Levels"].is_object());
  CHECK(meta["Test"]["Levels"]["0"]["Rows"] == 10);
  CHECK(meta["Test"]["Levels"]["1"]["Rows"] == 20);
  CHECK(meta["Test"]["Info"]["2"] == "x");
  perf::WriteJson(meta);  // Erase and rewrite
  CHECK(meta["Test"]["Levels"]["1"]["Rows"] == 20);
  BlockTimer::Reset();
}
