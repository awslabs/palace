// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/solverstats.hpp"
#include "utils/communication.hpp"
#include "utils/perfstats.hpp"
#include "utils/timer.hpp"

using namespace palace;
using Catch::Matchers::WithinRel;

namespace
{

double Value(const perf::Results &r, const std::string &path, bool cum = false)
{
  for (std::size_t i = 0; i < r.counter_paths.size(); i++)
  {
    if (r.counter_paths[i] == path)
    {
      const auto &s = cum ? r.cum_counter[i] : r.counter[i];
      return (s.n > 0) ? s.max : -1.0;
    }
  }
  return -1.0;
}

}  // namespace

TEST_CASE("SolverStats Histogram Buckets", "[perfstats][Serial]")
{
  CHECK(perf::KspStats::Bucket(0) == 0);
  CHECK(perf::KspStats::Bucket(1) == 1);
  CHECK(perf::KspStats::Bucket(3) == 2);
  CHECK(perf::KspStats::Bucket(4) == 3);
  CHECK(perf::KspStats::Bucket(1023) == 10);
  CHECK(perf::KspStats::Bucket(1 << 20) == 11);
  CHECK(perf::KspStats::BucketName(0) == "0-0");
  CHECK(perf::KspStats::BucketName(1) == "1-1");
  CHECK(perf::KspStats::BucketName(3) == "4-7");
  CHECK(perf::KspStats::BucketName(11) == "1024+");
}

TEST_CASE("SolverStats Krylov And Eigen", "[perfstats][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  BlockTimer::Reset();
  {
    perf::KspStats stats("TestRole");
    stats.Record(10, true, 0, 1.0, 1e-10);   // Rate 0.1
    stats.Record(2, true, 0, 2.0, 2e-2);     // Rate 0.1
    stats.Record(0, true, 0, 0.0, 0.0);      // Zero-iteration solve
    stats.Record(100, false, 3, 1.0, 1e-3);  // Not converged, 3 restarts
    perf::RecordEigenSolve("Main", 7, 4, 5);
    perf::RecordEigenSolve("Main", 3, 5, 5);
  }
  BlockTimer::Finalize(comm);
  const auto &r = perf::GetResults();
  auto V = [&](const std::string &f) { return Value(r, "LinearSolver/TestRole/" + f); };
  CHECK(V("Solves") == 4.0);
  CHECK(V("Its") == 112.0);
  CHECK(V("MaxIts") == 100.0);
  CHECK(V("NotConverged") == 1.0);
  CHECK(V("NotConvergedIts") == 100.0);
  CHECK(V("ZeroItSolves") == 1.0);
  CHECK(V("Restarts") == 3.0);
  CHECK(V("ItsHistogram/0-0") == 1.0);
  CHECK(V("ItsHistogram/2-3") == 1.0);
  CHECK(V("ItsHistogram/8-15") == 1.0);
  CHECK(V("ItsHistogram/64-127") == 1.0);
  CHECK(V("ItsHistogram/1-1") == -1.0);  // Untouched buckets are not written
  CHECK_THAT(V("MaxFinalRelRes"), WithinRel(1e-2, 1e-12));
  CHECK_THAT(V("MeanReductionFactor"),
             WithinRel((0.2 + std::pow(1e-3, 0.01)) / 3.0, 1e-12));
  CHECK(Value(r, "LinearSolver/TestRole/Its", true) == 112.0);
  CHECK(Value(r, "EigenSolver/Solves") == 2.0);
  CHECK(Value(r, "EigenSolver/OuterIts") == 10.0);
  CHECK(Value(r, "EigenSolver/MaxOuterIts") == 7.0);
  CHECK(Value(r, "EigenSolver/Converged") == 9.0);
  CHECK(Value(r, "EigenSolver/Requested") == 10.0);
  CHECK(Value(r, "EigenSolver/Incomplete") == 1.0);
  BlockTimer::Reset();
}
