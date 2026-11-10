// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <thread>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "linalg/chebyshev.hpp"
#include "linalg/iterative.hpp"
#include "linalg/operator.hpp"
#include "linalg/solverstats.hpp"
#include "utils/communication.hpp"
#include "utils/perfstats.hpp"

namespace palace
{

namespace
{

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

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

const perf::detail::TallyData &Tally(const std::string &name)
{
  const auto &reg = perf::detail::registry;
  return reg.tallies[reg.tally_ids.at(name)];
}

// A tallied estimate, and a second entry point with the same tally which calls it twice
// (as a complex-valued eigenvalue estimate built on the real-valued one).
void RealEstimate(int id)
{
  perf::ScopedTally t(id);
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
}

void ComplexEstimate(int id)
{
  perf::ScopedTally t(id);
  RealEstimate(id);
  RealEstimate(id);
}

}  // namespace

TEST_CASE("PerfTallies Reduce", "[perfstats][Serial][Parallel]")
{
  double buf[3] = {1.0, 2.0, 3.0};
  Mpi::GlobalSum(3, buf, Mpi::World());
  const auto &t = Tally("Reduce");
  const auto &reg = perf::detail::registry;
  const auto &c = reg.counters[reg.counter_ids.at("Communication/ReduceBytes")];
  const long count = t.count;
  const double bytes = c.value;
  Mpi::GlobalMax(3, buf, Mpi::World());
  CHECK(t.count == count + 1);
  CHECK(c.value == bytes + 3 * sizeof(double));
  CHECK(t.parent == Timer::TOTAL);
}

TEST_CASE("PerfTallies Krylov", "[perfstats][Serial][Parallel]")
{
  SpdTestOperator A;
  Vector b(2), x(2);
  b[0] = 1.0;
  b[1] = 2.0;

  CgSolver<Operator> cg(Mpi::World(), 0);
  cg.SetOperator(A);
  cg.SetRelTol(1.0e-12);
  cg.SetMaxIter(10);
  cg.SetPerfTallies("TestCG/Krylov", Timer::KSP);
  x = 0.0;
  cg.Mult(b, x);
  CHECK(cg.GetNumIterations() > 0);
  CHECK(Tally("TestCG/Krylov/ApplyBA").count >= cg.GetNumIterations());
  CHECK(Tally("TestCG/Krylov/ApplyBA").parent == Timer::KSP);

  GmresSolver<Operator> gmres(Mpi::World(), 0);
  gmres.SetOperator(A);
  gmres.SetRelTol(1.0e-12);
  gmres.SetMaxIter(10);
  gmres.SetPerfTallies("TestGMRES/Krylov", Timer::KSP);
  x = 0.0;
  gmres.Mult(b, x);
  const auto &apply = Tally("TestGMRES/Krylov/ApplyBA");
  const auto &orth = Tally("TestGMRES/Krylov/Orthogonalize");
  CHECK(apply.count == gmres.GetNumIterations());
  CHECK(orth.count == apply.count);
  CHECK(orth.time >= 0.0);
}

TEST_CASE("PerfTallies Nested Same Tally", "[perfstats][Serial]")
{
  BlockTimer::Reset();
  const int id = perf::LambdaMaxTally();
  const int outer = perf::RegisterTally("Test/Outer", Timer::KSP_SETUP);
  const auto &t = perf::detail::registry.tallies[id];
  const auto t0 = Timer::Now();
  RealEstimate(id);
  ComplexEstimate(id);
  const double elapsed = Timer::Duration(Timer::Now() - t0).count();
  CHECK(t.name == "Main/LambdaMax");
  CHECK(t.count == 2);
  CHECK(t.nested == 2);      // Visible, not counted
  CHECK(t.time <= elapsed);  // Counting the nested scopes would add >= 20 ms
  CHECK(t.time >= 0.03);

  // Scopes of other tallies nest as usual, and the tally counts again once closed.
  {
    perf::ScopedTally a(outer);
    RealEstimate(id);
  }
  CHECK(Tally("Test/Outer").count == 1);
  CHECK(t.count == 3);
  BlockTimer::Finalize(Mpi::World());
  nlohmann::json meta = nlohmann::json::object();
  perf::WriteJson(meta);
  CHECK(meta["ElapsedTime"]["Tallies"]["Main/LambdaMax"]["Nested"] == 2);
  CHECK(!meta["ElapsedTime"]["Tallies"]["Test/Outer"].contains("Nested"));
  BlockTimer::Reset();
}

TEST_CASE("PerfTallies Chebyshev", "[perfstats][Serial][Parallel]")
{
  SpdTestOperator A;
  constexpr int order = 3;
  ChebyshevSmoother<Operator> cheby(Mpi::World(), 1, order, 1.0);
  cheby.SetOperator(A);
  cheby.SetPerfTallies("TestCheby", Timer::KSP_PRECONDITIONER);
  Vector b(2), x(2);
  b = 1.0;
  cheby.SetInitialGuess(false);
  cheby.Mult(b, x);
  const auto &t = Tally("TestCheby/OpApply");
  CHECK(t.count == order - 1);
  cheby.SetInitialGuess(true);
  cheby.Mult(b, x);
  CHECK(t.count == 2 * order - 1);
  CHECK(t.parent == Timer::KSP_PRECONDITIONER);
}

}  // namespace palace
