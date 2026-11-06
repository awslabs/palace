// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_PERFSTATS_HPP
#define PALACE_UTILS_PERFSTATS_HPP

#include <algorithm>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <nlohmann/json_fwd.hpp>
#include "utils/timer.hpp"

namespace palace::perf
{

//
// Performance statistics reported in palace.json next to the exclusive BlockTimer rows:
// inclusive timing tallies, counters and per-setup information. All values are rank-local
// and reduced once per reporting iteration in BlockTimer::Finalize. Registration is cheap
// but not free (string lookup): register once (e.g. at setup) and keep the returned id
// for use in hot paths.
//

// palace.json schema version of the performance report.
inline constexpr int SCHEMA_VERSION = 2;

// How a counter is reduced across ranks.
enum class Reduce
{
  SUM,   // Sum over ranks
  MIN,   // Minimum over ranks
  MAX,   // Maximum over ranks (use for rank-uniform values)
  STATS  // {"Min", "Avg", "Max", "MaxRank"} over the ranks which set the counter
};

// How a counter is accumulated over iterations for the Cumulative block.
enum class Accumulate
{
  SUM,
  MIN,
  MAX,
  NONE  // Not reported in Cumulative (e.g. sizes)
};

namespace detail
{

struct TallyData
{
  std::string name;
  Timer::Index parent;
  double time = 0.0, cum_time = 0.0;
  long count = 0, cum_count = 0;
  int phase_base = -1;  // For phase tallies: first of NUM_TIMINGS sub-tally ids (lazy)
};

struct CounterData
{
  std::string path;
  Reduce reduce;
  Accumulate acc;
  double value = 0.0, cum = 0.0;
  bool set = false, cum_set = false;
};

struct Registry
{
  std::vector<TallyData> tallies;
  std::vector<CounterData> counters;
  std::unordered_map<std::string, int> tally_ids, counter_ids;
  std::vector<std::string> roles;
  bool suspended = false;  // Set while BlockTimer::Finalize reduces (not tallied)
};

void AddPhaseTally(int id, double seconds, long count);

inline Registry registry;

}  // namespace detail

//
// Inclusive tallies: stackless, they do not interrupt BlockTimer rows and may overlap each
// other. Each carries the exclusive Timer row it is measured within (Parent). Reported in
// palace.json under ElapsedTime/Tallies/<name>.
//

// Register (or look up) a tally. Ids are stable for the lifetime of the process.
int RegisterTally(std::string_view name, Timer::Index parent);

// A tally measured in any phase (Parent TOTAL), additionally split by the exclusive row
// active at the time into "<name>/<ROW>" sub-tallies (Parent ROW).
int RegisterPhaseTally(std::string_view name);

inline void AddTally(int id, double seconds, long count = 1)
{
  if (id >= 0)
  {
    auto &t = detail::registry.tallies[id];
    t.time += seconds;
    t.count += count;
    if (t.phase_base != -1)
    {
      detail::AddPhaseTally(id, seconds, count);
    }
  }
}

// Time a scope into a tally (two clock reads, no system calls). A negative id is a no-op.
class ScopedTally
{
  int id;
  Timer::TimePoint t0;

public:
  explicit ScopedTally(int id) : id(id)
  {
    if (id >= 0)
    {
      t0 = Timer::Now();
    }
  }
  ~ScopedTally()
  {
    if (id >= 0)
    {
      AddTally(id, Timer::Duration(Timer::Now() - t0).count());
    }
  }
  ScopedTally(const ScopedTally &) = delete;
  ScopedTally &operator=(const ScopedTally &) = delete;
};

//
// Counters: a number written to palace.json at a slash-separated path (for example
// "LinearSolver/DivFree/Solves"). Counters are reset at the start of each iteration and a
// counter that was never set in an iteration is not written.
//

int RegisterCounter(std::string_view path, Reduce reduce = Reduce::MAX,
                    Accumulate acc = Accumulate::SUM);

inline void Add(int id, double v)
{
  auto &c = detail::registry.counters[id];
  c.value += v;
  c.set = true;
}
inline void Set(int id, double v)
{
  auto &c = detail::registry.counters[id];
  c.value = v;
  c.set = true;
}
inline void SetMax(int id, double v)
{
  auto &c = detail::registry.counters[id];
  c.value = c.set ? std::max(c.value, v) : v;
  c.set = true;
}
inline void SetMin(int id, double v)
{
  auto &c = detail::registry.counters[id];
  c.value = c.set ? std::min(c.value, v) : v;
  c.set = true;
}

// Convenience overloads registering by path (string lookup; avoid in hot paths).
inline void Add(std::string_view path, double v, Reduce reduce = Reduce::MAX,
                Accumulate acc = Accumulate::SUM)
{
  Add(RegisterCounter(path, reduce, acc), v);
}
inline void Set(std::string_view path, double v, Reduce reduce = Reduce::MAX,
                Accumulate acc = Accumulate::NONE)
{
  Set(RegisterCounter(path, reduce, acc), v);
}

//
// Per-setup information (rank-uniform values such as resolved solver choices). The global
// root's value is written at the slash-separated path as is. Cleared at each iteration;
// not part of Cumulative.
//
void SetInfo(std::string_view path, const nlohmann::json &value);

//
// Solver roles: the KSP solvers and the tallies they register are keyed by the role active
// at their construction ("Main" if none).
//
std::string_view CurrentRole();

class ScopedRole
{
public:
  explicit ScopedRole(std::string_view role);
  ~ScopedRole();
  ScopedRole(const ScopedRole &) = delete;
  ScopedRole &operator=(const ScopedRole &) = delete;
};

// Path of a linear solver statistic for the given role: the main solver's solve and
// iteration totals keep their historical names (LinearSolver/TotalSolves, TotalIts),
// everything else is under LinearSolver/<role>/<field>.
std::string LinearSolverPath(std::string_view role, std::string_view field);

//
// Reduced results, valid on all ranks after BlockTimer::Finalize.
//

struct Stats
{
  double min = 0.0, max = 0.0, sum = 0.0, avg = 0.0;
  int max_rank = -1;  // Rank (or node index) of the maximum
  int n = 0;          // Number of contributing ranks (or nodes)
};

struct TimerResults
{
  std::vector<Stats> time, count, rank_mem, node_mem;  // Per Timer::Index
  std::vector<std::string> tally_names;
  std::vector<Timer::Index> tally_parents;
  std::vector<Stats> tally_time, tally_count;
};

struct Results
{
  int iteration = 1;
  int num_nodes = 0;
  bool hwm_reset = false;   // Whether the peak RSS was reset at the iteration start
  TimerResults iter, cum;   // This iteration and iterations 1..N
  Stats peak_mem, cur_mem;  // Per-rank peak (this iteration) and current RSS
  Stats peak_node_mem, cur_node_mem;
  Stats cum_peak_mem, cum_peak_node_mem;  // Max over iterations
  // Peak memory attribution: peak = baseline (peak RSS at the iteration start) + sum of the
  // per-row growth + unattributed growth (outside any timer row).
  Stats mem_baseline, mem_unattributed, node_mem_baseline, node_mem_unattributed;
  std::vector<std::string> counter_paths;
  std::vector<Reduce> counter_reduce;
  std::vector<Accumulate> counter_acc;
  std::vector<Stats> counter, cum_counter;
};

const Results &GetResults();

// Write the performance report into the palace.json object (global root only). Removes
// keys written by previous iterations first.
void WriteJson(nlohmann::json &meta);

}  // namespace palace::perf

#endif  // PALACE_UTILS_PERFSTATS_HPP
