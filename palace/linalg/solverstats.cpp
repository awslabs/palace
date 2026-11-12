// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "solverstats.hpp"

#include <cmath>
#include <fmt/format.h>
#include "utils/perfstats.hpp"

namespace palace::perf
{

int KspStats::Bucket(int its)
{
  int b = 0;
  while (its > 0 && b < NUM_BUCKETS - 1)
  {
    its >>= 1;
    b++;
  }
  return b;
}

std::string KspStats::BucketName(int b)
{
  if (b <= 1)
  {
    return fmt::format("{0}-{0}", b);  // Not a bare number, which is a JSON array index
  }
  if (b == NUM_BUCKETS - 1)
  {
    return fmt::format("{}+", 1 << (b - 1));
  }
  return fmt::format("{}-{}", 1 << (b - 1), (1 << b) - 1);
}

KspStats::KspStats(std::string_view role)
{
  auto Path = [&](std::string_view f) { return LinearSolverPath(role, f); };
  solves = RegisterCounter(Path("Solves"));
  its = RegisterCounter(Path("Its"));
  max_its = RegisterCounter(Path("MaxIts"), Reduce::MAX, Accumulate::MAX);
  not_converged = RegisterCounter(Path("NotConverged"));
  not_converged_its = RegisterCounter(Path("NotConvergedIts"));
  zero_its = RegisterCounter(Path("ZeroItSolves"));
  restarts = RegisterCounter(Path("Restarts"));
  mean_reduction =
      RegisterCounter(Path("MeanReductionFactor"), Reduce::MAX, Accumulate::NONE);
  min_reduction = RegisterCounter(Path("MinReductionFactor"), Reduce::MIN, Accumulate::MIN);
  max_reduction = RegisterCounter(Path("MaxReductionFactor"), Reduce::MAX, Accumulate::MAX);
  max_rel_res = RegisterCounter(Path("MaxFinalRelRes"), Reduce::MAX, Accumulate::MAX);
  for (int b = 0; b < NUM_BUCKETS; b++)
  {
    histogram[b] = RegisterCounter(Path("ItsHistogram/" + BucketName(b)));
  }
}

void KspStats::Record(int n, bool converged, int nrestart, double initial_res,
                      double final_res) const
{
  Add(solves, 1);
  Add(its, n);
  SetMax(max_its, n);
  Add(histogram[Bucket(n)], 1);
  if (!converged)
  {
    Add(not_converged, 1);
    Add(not_converged_its, n);
  }
  if (n == 0)
  {
    Add(zero_its, 1);
  }
  if (nrestart > 0)
  {
    Add(restarts, nrestart);
  }
  if (initial_res > 0.0)
  {
    const double rel = final_res / initial_res;
    SetMax(max_rel_res, rel);
    if (n > 0)
    {
      // Running mean over the solves with at least one iteration.
      const double m = Get(solves) - Get(zero_its);
      const double r = std::pow(rel, 1.0 / n);
      const double mean = IsSet(mean_reduction) ? Get(mean_reduction) : 0.0;
      Set(mean_reduction, (m > 1.0) ? mean + (r - mean) / m : r);
      SetMin(min_reduction, r);
      SetMax(max_reduction, r);
    }
  }
}

void RecordEigenSolve(std::string_view role, int outer_its, int num_conv, int num_req)
{
  const std::string prefix =
      (role == "Main") ? "EigenSolver/" : fmt::format("EigenSolver/{}/", role);
  Add(prefix + "Solves", 1);
  Add(prefix + "OuterIts", outer_its);
  SetMax(RegisterCounter(prefix + "MaxOuterIts", Reduce::MAX, Accumulate::MAX), outer_its);
  Add(prefix + "Converged", num_conv);
  Add(prefix + "Requested", num_req);
  Add(prefix + "Incomplete", (num_conv < num_req) ? 1 : 0);
}

int LambdaMaxTally()
{
  return RegisterTally(fmt::format("{}/LambdaMax", CurrentRole()));
}

}  // namespace palace::perf
