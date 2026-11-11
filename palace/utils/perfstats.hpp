// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_PERFSTATS_HPP
#define PALACE_UTILS_PERFSTATS_HPP

#include <algorithm>
#include <cmath>
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
  // Exclusive timer row the tally is measured within: fixed for phase tallies (TOTAL, and
  // the row of each sub-tally), otherwise the row active at each measurement (-1 until
  // measured, TOTAL once measured within several rows), this iteration and over all.
  int fixed_parent = -1;
  int row = -1, cum_row = -1;
  double time = 0.0, cum_time = 0.0;
  double wait = 0.0, cum_wait = 0.0;  // Part of time blocked in communication
  long count = 0, cum_count = 0;
  long nested = 0, cum_nested = 0;  // Scopes opened inside a scope of the same tally
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
  int paused = 0;          // ScopedPause depth: tallies and counters record nothing
};

void AddPhaseTally(int id, double seconds, long count, double wait);

// Row of a tally measured within rows a and b (-1: none): TOTAL if they differ.
inline int MergeRow(int a, int b)
{
  return (a < 0 || a == b) ? b : ((b < 0) ? a : static_cast<int>(Timer::TOTAL));
}

inline Registry registry;

// Seconds this thread has spent blocked in communication: in the MPI_Wait family (requests
// that had not completed yet, see commstats.cpp) and in Palace's blocking reductions.
inline thread_local double comm_wait = 0.0;

// Whether a ScopedTally of each tally id is open on this thread.
inline thread_local std::vector<char> open_tallies;

}  // namespace detail

//
// Inclusive tallies: stackless, they do not interrupt BlockTimer rows and may overlap each
// other. Each reports the exclusive Timer row it is measured within (Parent): the row
// active at each measurement (BlockTimer::Current(), when a ScopedTally closes), TOTAL if
// it was measured within several rows or between top-level timers. Its time includes that
// of the rows nested in its scopes. Reported in palace.json under
// ElapsedTime/Tallies/<name>, with the time imbalance over ranks and, when nonzero, the
// time blocked in communication within the tally (Wait) and the rest (Work = time - Wait),
// each reduced over ranks.
//

// Register (or look up) a tally. Ids are stable for the lifetime of the process.
int RegisterTally(std::string_view name);

// A tally measured in any phase (Parent TOTAL), additionally split by the exclusive row
// active at the time into "<name>/<ROW>" sub-tallies (Parent ROW). The sub-tallies are
// registered on first use, or all at once if eager (for tallies whose rows depend on the
// rank, to keep the registration order the same on all ranks).
int RegisterPhaseTally(std::string_view name, bool eager = false);

// Id of the sub-tally of a phase tally for an exclusive row (registered if needed).
int PhaseTallyRow(int id, Timer::Index row);

// Register the communication tallies and counters (Reduce, Barrier, Broadcast, Wait/MPI).
// Call it where all ranks reach it in the same order, before any rank-dependent
// registration; otherwise they are registered at the first communication.
void RegisterCommunicationTallies();

// Add to a tally: elapsed seconds, of which wait seconds were blocked in communication.
inline void AddTally(int id, double seconds, long count = 1, double wait = 0.0)
{
  if (id >= 0 && detail::registry.paused == 0)
  {
    auto &t = detail::registry.tallies[id];
    t.time += seconds;
    t.wait += wait;
    t.count += count;
    if (t.phase_base != -1)
    {
      detail::AddPhaseTally(id, seconds, count, wait);
    }
    else if (t.fixed_parent < 0)
    {
      t.row = detail::MergeRow(t.row, BlockTimer::Current());
    }
  }
}

// Time a scope into a tally (two clock reads, no system calls), with the communication
// wait inside it. A negative id is a no-op. As a safeguard, so is a scope nested in an open
// scope of the same tally on the same thread (for example a tallied helper which a tallied
// entry point also calls): the outermost scope counts the time once, and the nested scope
// is counted as Nested in the report.
class ScopedTally
{
  int id;
  Timer::TimePoint t0;
  double w0 = 0.0;

public:
  explicit ScopedTally(int id) : id(id)
  {
    if (id >= 0)
    {
      auto &open = detail::open_tallies;
      if (static_cast<std::size_t>(id) >= open.size())
      {
        open.resize(
            std::max(static_cast<std::size_t>(id) + 1, detail::registry.tallies.size()), 0);
      }
      if (open[id])
      {
        if (detail::registry.paused == 0)
        {
          detail::registry.tallies[id].nested++;
        }
        this->id = -1;
        return;
      }
      open[id] = 1;
      w0 = detail::comm_wait;
      t0 = Timer::Now();
    }
  }
  ~ScopedTally()
  {
    if (id >= 0)
    {
      const double seconds = Timer::Duration(Timer::Now() - t0).count();
      detail::open_tallies[id] = 0;
      AddTally(id, seconds, 1, detail::comm_wait - w0);
    }
  }
  ScopedTally(const ScopedTally &) = delete;
  ScopedTally &operator=(const ScopedTally &) = delete;
};

// Pause all reporting in a scope: nested BlockTimer rows, tallies and counters record
// nothing, and the time stays in the row active before. For verification probes, which are
// measured by their own tally around the scope instead.
class ScopedPause
{
public:
  ScopedPause()
  {
    detail::registry.paused++;
    BlockTimer::Pause(true);
  }
  ~ScopedPause()
  {
    detail::registry.paused = std::max(detail::registry.paused - 1, 0);  // See ClearPause
    BlockTimer::Pause(false);
  }
  ScopedPause(const ScopedPause &) = delete;
  ScopedPause &operator=(const ScopedPause &) = delete;
};

// End all pauses, for a failure inside a ScopedPause (see FailureReport).
inline void ClearPause()
{
  detail::registry.paused = 0;
  BlockTimer::ClearPause();
}

//
// Counters: a number written to palace.json at a slash-separated path (for example
// "LinearSolver/DivFree/Solves"). Counters are reset at the start of each iteration and a
// counter that was never set in an iteration is not written.
//

int RegisterCounter(std::string_view path, Reduce reduce = Reduce::MAX,
                    Accumulate acc = Accumulate::SUM);

inline void Add(int id, double v)
{
  if (detail::registry.paused)
  {
    return;
  }
  auto &c = detail::registry.counters[id];
  c.value += v;
  c.set = true;
}
// Local (rank, current iteration) value of a counter, and whether it was set.
inline double Get(int id)
{
  return detail::registry.counters[id].value;
}
inline bool IsSet(int id)
{
  return detail::registry.counters[id].set;
}

inline void Set(int id, double v)
{
  if (detail::registry.paused)
  {
    return;
  }
  auto &c = detail::registry.counters[id];
  c.value = v;
  c.set = true;
}
// SetMax and SetMin keep a value which is not a number (NaN) once set.
inline void SetMax(int id, double v)
{
  if (detail::registry.paused)
  {
    return;
  }
  auto &c = detail::registry.counters[id];
  if (!c.set || std::isnan(v) || (!std::isnan(c.value) && v > c.value))
  {
    c.value = v;
  }
  c.set = true;
}
inline void SetMin(int id, double v)
{
  if (detail::registry.paused)
  {
    return;
  }
  auto &c = detail::registry.counters[id];
  if (!c.set || std::isnan(v) || (!std::isnan(c.value) && v < c.value))
  {
    c.value = v;
  }
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
// Solver roles: the solvers and the tallies they register are keyed by the role active at
// their construction ("Main" if none). A sub-role "<role>/<part>" keeps apart the parts
// of one solver of a role (for example the blocks of a block preconditioner).
//
std::string_view CurrentRole();

// Whether a role is the main solver's or one of its sub-roles: verification checks only
// run for these, which run on all processes of the problem.
inline bool IsMainRole(std::string_view role)
{
  return role == "Main" || role.starts_with("Main/");
}

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
  std::vector<Stats> tally_time, tally_count, tally_wait, tally_work, tally_nested;
  std::vector<Stats> tally_active;  // Sum: number of ranks with a nonzero count
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
  // Heap in use (allocator statistics, see memory_reporting::GetHeapMemory): growth per
  // row, at the iteration start (baseline), outside any row (unattributed) and at the end
  // of the iteration. Baseline + rows + unattributed = end, per rank.
  // Unattributed also holds the growth across skipped samples; heap_coverage is the
  // fraction of each row's samples charged to it.
  std::vector<Stats> heap_growth, heap_coverage;
  Stats heap_baseline, heap_unattributed, heap_mem;
  std::vector<std::string> counter_paths;
  std::vector<Reduce> counter_reduce;
  std::vector<Accumulate> counter_acc;
  std::vector<Stats> counter, cum_counter;
};

const Results &GetResults();

// Write the performance report into the palace.json object (global root only). Removes
// keys written by previous iterations first. "Complete" is true, unless an error (see
// FailureReport) is given: then it is false and the error is written under "Error". If the
// error is not "Collective", the report holds the root process's own values and is written
// under "RootProcess" instead.
void WriteJson(nlohmann::json &meta, const nlohmann::json *error = nullptr);

}  // namespace palace::perf

#endif  // PALACE_UTILS_PERFSTATS_HPP
