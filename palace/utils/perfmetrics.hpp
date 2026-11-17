// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_PERFMETRICS_HPP
#define PALACE_UTILS_PERFMETRICS_HPP

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <deque>
#include <string>
#include <string_view>
#include <vector>
#include <mpi.h>

namespace palace::perf
{

// Statistics of a value over the processes which contribute to it.
struct Stats
{
  double min = 0.0, max = 0.0, sum = 0.0;
  int max_rank = -1;  // Lowest rank holding the maximum
  int ranks = 0;      // Number of contributing processes (0: set nowhere)
  double Avg() const { return (ranks > 0) ? sum / ranks : 0.0; }
};

// Reduce each value over the processes of comm; value i contributes on a process if
// set[i]. Collective: three MPI_Allreduce.
std::vector<Stats> ReduceStats(const std::vector<double> &value,
                               const std::vector<char> &set, MPI_Comm comm);

// Reduction across processes. COUNT: a count which all processes running the work perform
// identically (MAX; SUM over subcommunicator roots for a role on subcommunicators).
enum class Red
{
  STATS,
  COUNT,
  MAX,
  MIN,
  SUM
};

// Cumulative rule: ITERATION_ONLY metrics (sizes, settings) are not part of Cumulative.
enum class Run
{
  ACCUMULATE,
  ITERATION_ONLY
};

// Accumulation within a process, set by the handle type.
enum class Acc
{
  SUM,
  MAX,
  MIN,
  LAST
};

namespace detail
{

struct Entry
{
  std::string path;
  Acc acc;
  Red red;  // Resolved: never COUNT
  Run run;
  double iter = 0.0, run_value = 0.0;
  bool iter_set = false, run_set = false;
};

inline std::vector<Entry> entries;

inline double Combine(Acc acc, double a, double b)
{
  switch (acc)
  {
    case Acc::SUM:
      return std::clamp(a + b, -DBL_MAX, DBL_MAX);
    case Acc::MAX:
      return std::max(a, b);
    case Acc::MIN:
      return std::min(a, b);
    case Acc::LAST:
      return b;
  }
  return b;
}

// Hot path: one indexed update. Invalid ids are no-ops; NaN is ignored, inf clamped.
inline void Update(int id, double v)
{
  if (id < 0 || std::isnan(v))
  {
    return;
  }
  v = std::clamp(v, -DBL_MAX, DBL_MAX);
  auto &e = entries[id];
  e.iter = e.iter_set ? Combine(e.acc, e.iter, v) : v;
  e.iter_set = true;
}

int Register(std::string_view path, Acc acc, Red red, Run run);

// Whether a path is valid and not under a prefix owned by another writer.
bool ValidPath(std::string_view path);

// Fold the iteration values into the run values and clear them (iteration boundary).
void FoldIteration();

// Clear all values, keeping registrations (BlockTimer::Reset).
void ResetValues();

}  // namespace detail

//
// Roles: the solver a metric belongs to (LinearSolver/<Role>/...). Interned: a RoleId is
// an index into a process-wide table.
//

struct RoleId
{
  int index = 0;  // 0: Unassigned
};

class ScopedRole
{
public:
  explicit ScopedRole(std::string_view name);
  // A role on subcommunicators: its COUNT metrics are updated on the root of comm only and
  // summed over processes. Processes with comm == MPI_COMM_NULL do not contribute.
  ScopedRole(std::string_view name, MPI_Comm comm);
  explicit ScopedRole(RoleId id);
  ~ScopedRole();
  ScopedRole(const ScopedRole &) = delete;
  ScopedRole &operator=(const ScopedRole &) = delete;
};

RoleId CurrentRoleId();
const std::string &RoleName(RoleId id);
inline const std::string &CurrentRole()
{
  return RoleName(CurrentRoleId());
}

//
// Metric handles. Construction registers the metric (string lookup: at setup, or as a
// function-local static); a default-constructed handle is invalid and its updates no-ops.
// Handles stay valid for the lifetime of the process.
//

class Metric
{
protected:
  int id = -1;
  Metric() = default;
  Metric(std::string_view path, Acc acc, Red red, Run run);

public:
  bool Valid() const { return id >= 0; }
};

class SumMetric : public Metric
{
public:
  SumMetric() = default;
  SumMetric(std::string_view path, Red red, Run run = Run::ACCUMULATE)
    : Metric(path, Acc::SUM, red, run)
  {
  }
  void Add(double v) const { detail::Update(id, v); }
};

class MaxMetric : public Metric
{
public:
  MaxMetric() = default;
  MaxMetric(std::string_view path, Red red, Run run = Run::ACCUMULATE)
    : Metric(path, Acc::MAX, red, run)
  {
  }
  void Update(double v) const { detail::Update(id, v); }
};

class MinMetric : public Metric
{
public:
  MinMetric() = default;
  MinMetric(std::string_view path, Red red, Run run = Run::ACCUMULATE)
    : Metric(path, Acc::MIN, red, run)
  {
  }
  void Update(double v) const { detail::Update(id, v); }
};

class LastMetric : public Metric
{
public:
  LastMetric() = default;
  LastMetric(std::string_view path, Red red, Run run = Run::ACCUMULATE)
    : Metric(path, Acc::LAST, red, run)
  {
  }
  void Set(double v) const { detail::Update(id, v); }
};

}  // namespace palace::perf

#endif  // PALACE_UTILS_PERFMETRICS_HPP
