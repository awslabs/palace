// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "perfstats.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <tuple>
#include <nlohmann/json.hpp>
#include "utils/communication.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/perfprobes.hpp"

namespace palace
{

using json = nlohmann::json;

namespace perf
{

namespace
{

// Reporting state outside of the registry.
struct State
{
  Results results;
  bool finalized = false;
  int iteration = 1;
  json info = json::object();
  std::set<std::string> info_paths;
  std::set<std::string> counter_paths;  // All counter paths written by the report so far

  // Peak memory bookkeeping (rank-local), see IterationPeakMemory.
  long peak_at_start = 0, cur_at_start = 0;
  double cum_peak = 0.0;
  double cum_node_peak = 0.0, last_node_peak = 0.0;
};

State state;

// Reduction entry: {min, max, sum, rank of max, number of contributions}.
constexpr int STRIDE = 5;

void PushEntry(std::vector<double> &buf, double v, bool set, int rank)
{
  if (set)
  {
    buf.insert(buf.end(), {v, v, v, static_cast<double>(rank), 1.0});
  }
  else
  {
    constexpr double inf = std::numeric_limits<double>::infinity();
    buf.insert(buf.end(), {inf, -inf, 0.0, std::numeric_limits<double>::max(), 0.0});
  }
}

void StatsReduceOp(void *in, void *inout, int *len, MPI_Datatype *)
{
  const double *a = static_cast<const double *>(in);
  double *b = static_cast<double *>(inout);
  for (int k = 0; k < *len; k++, a += STRIDE, b += STRIDE)
  {
    // A value which is not a number (NaN) wins.
    if (!std::isnan(b[0]) && (std::isnan(a[0]) || a[0] < b[0]))
    {
      b[0] = a[0];
    }
    if (!std::isnan(b[1]) &&
        (std::isnan(a[1]) || a[1] > b[1] || (a[1] == b[1] && a[3] < b[3])))
    {
      b[1] = a[1];
      b[3] = a[3];
    }
    b[2] += a[2];
    b[4] += a[4];
  }
}

std::pair<MPI_Datatype, MPI_Op> StatsTypeAndOp()
{
  static MPI_Datatype type = MPI_DATATYPE_NULL;
  static MPI_Op op = MPI_OP_NULL;
  if (type == MPI_DATATYPE_NULL)
  {
    MPI_Type_contiguous(STRIDE, MPI_DOUBLE, &type);
    MPI_Type_commit(&type);
    MPI_Op_create(&StatsReduceOp, 1, &op);
  }
  return {type, op};
}

void StatsAllreduce(std::vector<double> &buf, MPI_Comm comm)
{
  if (buf.empty())
  {
    return;
  }
  auto [type, op] = StatsTypeAndOp();
  MPI_Allreduce(MPI_IN_PLACE, buf.data(), static_cast<int>(buf.size() / STRIDE), type, op,
                comm);
}

Stats Unpack(const double *e)
{
  Stats s;
  s.n = static_cast<int>(e[4]);
  if (s.n > 0)
  {
    s.min = e[0];
    s.max = e[1];
    s.sum = e[2];
    s.avg = e[2] / s.n;
    s.max_rank = static_cast<int>(e[3]);
  }
  return s;
}

std::uint64_t Hash(std::uint64_t h, std::string_view s)
{
  for (unsigned char c : s)
  {
    h = (h ^ c) * 1099511628211ull;
  }
  return (h ^ 0xff) * 1099511628211ull;
}

std::string TallyKey(const detail::TallyData &t)
{
  return t.name;
}

std::string CounterKey(const detail::CounterData &c)
{
  return c.path + '\x1f' + std::to_string(static_cast<int>(c.reduce)) +
         std::to_string(static_cast<int>(c.acc));
}

// For each canonical (rank-consistent) entry, the local index or -1, and the keys. Fast
// path when all ranks registered the same entries in the same order (reported in
// verification mode as Verification/ConsistentRegistration).
std::tuple<std::vector<int>, std::vector<std::string>, std::vector<int>,
           std::vector<std::string>>
CanonicalOrder(MPI_Comm comm)
{
  const auto &reg = detail::registry;
  std::vector<std::string> tkeys, ckeys;
  std::uint64_t ht = 14695981039346656037ull, hc = ht;
  for (const auto &t : reg.tallies)
  {
    tkeys.push_back(TallyKey(t));
    ht = Hash(ht, tkeys.back());
  }
  for (const auto &c : reg.counters)
  {
    ckeys.push_back(CounterKey(c));
    hc = Hash(hc, ckeys.back());
  }
  std::int64_t chk[8] = {
      static_cast<std::int64_t>(tkeys.size()), static_cast<std::int64_t>(ht >> 1),
      static_cast<std::int64_t>(ckeys.size()), static_cast<std::int64_t>(hc >> 1)};
  for (int i = 0; i < 4; i++)
  {
    chk[4 + i] = -chk[i];
  }
  Mpi::GlobalMax(8, chk, comm);
  const bool consistent =
      chk[0] == -chk[4] && chk[1] == -chk[5] && chk[2] == -chk[6] && chk[3] == -chk[7];
  auto Identity = [](std::size_t n)
  {
    std::vector<int> v(n);
    for (std::size_t i = 0; i < n; i++)
    {
      v[i] = static_cast<int>(i);
    }
    return v;
  };
  if (Verification())
  {
    SetInfo("Verification/ConsistentRegistration", consistent);
  }
  if (consistent)
  {
    return {Identity(tkeys.size()), tkeys, Identity(ckeys.size()), ckeys};
  }

  // Slow path: union of all ranks' entries, sorted.
  auto Union = [&comm](const std::vector<std::string> &keys)
  {
    std::string local;
    for (const auto &k : keys)
    {
      local += k;
      local += '\0';
    }
    const int np = Mpi::Size(comm);
    int len = static_cast<int>(local.size());
    std::vector<int> lens(np), displs(np, 0);
    MPI_Allgather(&len, 1, MPI_INT, lens.data(), 1, MPI_INT, comm);
    for (int i = 1; i < np; i++)
    {
      displs[i] = displs[i - 1] + lens[i - 1];
    }
    std::string all(displs[np - 1] + lens[np - 1], '\0');
    MPI_Allgatherv(local.data(), len, MPI_CHAR, all.data(), lens.data(), displs.data(),
                   MPI_CHAR, comm);
    std::set<std::string> keyset;
    std::size_t pos = 0;
    while (pos < all.size())
    {
      auto end = all.find('\0', pos);
      keyset.insert(all.substr(pos, end - pos));
      pos = end + 1;
    }
    std::vector<std::string> canon(keyset.begin(), keyset.end());
    std::vector<int> local_idx(canon.size(), -1);
    for (std::size_t i = 0; i < canon.size(); i++)
    {
      for (std::size_t j = 0; j < keys.size(); j++)
      {
        if (keys[j] == canon[i])
        {
          local_idx[i] = static_cast<int>(j);
          break;
        }
      }
    }
    return std::make_pair(local_idx, canon);
  };
  auto [tidx, tcanon] = Union(tkeys);
  auto [cidx, ccanon] = Union(ckeys);
  return {tidx, tcanon, cidx, ccanon};
}

double AccumulateValue(Accumulate acc, double cum, bool cum_set, double v, bool set)
{
  if (!set)
  {
    return cum;
  }
  if (!cum_set)
  {
    return v;
  }
  switch (acc)
  {
    case Accumulate::SUM:
      return cum + v;
    case Accumulate::MIN:
      return (std::isnan(cum) || std::isnan(v)) ? std::numeric_limits<double>::quiet_NaN()
                                                : std::min(cum, v);
    case Accumulate::MAX:
      return (std::isnan(cum) || std::isnan(v)) ? std::numeric_limits<double>::quiet_NaN()
                                                : std::max(cum, v);
    case Accumulate::NONE:
      break;
  }
  return v;
}

// Peak RSS of this rank within the current iteration, and whether it is exact. The peak RSS
// of the process is not reset (the reset would also change the peak seen by getrusage,
// wait4 and the accounting of job schedulers): a peak of the iteration which did not raise
// it is unknown and bounded below by the RSS at the iteration boundaries and at the timer
// events where it was sampled (rss_sampled, see Timer::RssSampled).
std::pair<double, bool> IterationPeakMemory(long rss_sampled)
{
  const long peak = memory_reporting::GetPeakMemory();
  if (peak > state.peak_at_start)
  {
    return {static_cast<double>(peak), true};
  }
  return {static_cast<double>(std::max(
              {state.cur_at_start, rss_sampled, memory_reporting::GetCurrentMemory()})),
          false};
}

json Num(double v)
{
  if (std::isfinite(v) && v == std::round(v) && std::abs(v) < 9.0e15)
  {
    return static_cast<std::int64_t>(v);
  }
  return v;
}

json::json_pointer Pointer(std::string_view path)
{
  std::string p(path);
  if (p.empty() || p[0] != '/')
  {
    p.insert(p.begin(), '/');
  }
  return json::json_pointer(p);
}

// Element at a slash-separated path, creating objects as needed. Every token is an object
// key (json_pointer would create arrays for numeric tokens).
json &At(json &j, std::string_view path)
{
  json *node = &j;
  std::size_t pos = 0;
  while (pos <= path.size())
  {
    auto end = path.find('/', pos);
    if (end == std::string_view::npos)
    {
      end = path.size();
    }
    if (end > pos)
    {
      if (!node->is_object())
      {
        *node = json::object();
      }
      node = &(*node)[std::string(path.substr(pos, end - pos))];
    }
    pos = end + 1;
  }
  return *node;
}

void ErasePath(json &meta, std::string_view path)
{
  auto ptr = Pointer(path);
  if (!meta.contains(ptr))
  {
    return;
  }
  // Remove the key and any parents left empty.
  while (!ptr.empty())
  {
    auto parent = ptr.parent_pointer();
    meta[parent].erase(ptr.back());
    if (parent.empty() || !meta[parent].empty())
    {
      break;
    }
    ptr = parent;
  }
}

json StatsJson(const Stats &s, double scale = 1.0, const char *rank_key = "MaxRank")
{
  return {{"Min", Num(s.min * scale)},
          {"Avg", Num(s.avg * scale)},
          {"Max", Num(s.max * scale)},
          {rank_key, s.max_rank}};
}

json CounterJson(Reduce reduce, const Stats &s)
{
  switch (reduce)
  {
    case Reduce::SUM:
      return Num(s.sum);
    case Reduce::MIN:
      return Num(s.min);
    case Reduce::MAX:
      return Num(s.max);
    case Reduce::STATS:
      return StatsJson(s);
  }
  return nullptr;
}

// Max/avg over the ranks which ran the tally (1 is balanced), when anything was measured.
void SetImbalance(json &j, const Stats &s, const Stats &active)
{
  if (s.sum > 0.0 && active.sum > 0.0)
  {
    j["Imbalance"] = s.max * active.sum / s.sum;
  }
}

json ElapsedTimeJson(const TimerResults &r)
{
  json et;
  for (int i = 0; i < Timer::NUM_TIMINGS; i++)
  {
    std::string key(Timer::names[i]);
    et["Durations"][key] = r.time[i].avg;
    et["Min"][key] = r.time[i].min;
    et["Max"][key] = r.time[i].max;
    et["MaxRank"][key] = r.time[i].max_rank;
    et["Counts"][key] = Num(r.count[i].max);
  }
  json tallies = json::object();
  for (std::size_t i = 0; i < r.tally_names.size(); i++)
  {
    if (r.tally_count[i].max > 0)
    {
      auto t = StatsJson(r.tally_time[i]);
      t["Parent"] = Timer::names[r.tally_parents[i]];
      t["Count"] = Num(r.tally_count[i].max);
      SetImbalance(t, r.tally_time[i], r.tally_active[i]);
      if (r.tally_nested[i].max > 0)
      {
        t["Nested"] = Num(r.tally_nested[i].max);
      }
      if (r.tally_wait[i].max > 0.0)
      {
        t["Wait"] = StatsJson(r.tally_wait[i]);
        if (r.tally_work[i].max >
            1.0e-3 * r.tally_time[i].max)  // Not for (almost) pure waits
        {
          t["Work"] = StatsJson(r.tally_work[i]);
          SetImbalance(t["Work"], r.tally_work[i], r.tally_active[i]);
        }
      }
      tallies[r.tally_names[i]] = std::move(t);
    }
  }
  et["Tallies"] = std::move(tallies);
  return et;
}

json PeakJson(const Stats &s, bool node)
{
  constexpr double to_mb = 1.0 / (1024.0 * 1024.0);
  return {{"Min", s.min * to_mb},
          {"Max", s.max * to_mb},
          {"Average", s.avg * to_mb},
          {"Total", s.sum * to_mb},
          {node ? "MaxNode" : "MaxRank", s.max_rank}};
}

// Peak memory of the iteration, flagged as a lower bound if it is for any process.
json IterationPeakJson(const Stats &s, bool node, bool lower_bound)
{
  auto j = PeakJson(s, node);
  j["LowerBound"] = lower_bound;
  return j;
}

void Commit()
{
  auto &reg = detail::registry;
  for (auto &t : reg.tallies)
  {
    t.cum_time += t.time;
    t.cum_wait += t.wait;
    t.cum_count += t.count;
    t.cum_nested += t.nested;
    t.cum_row = detail::MergeRow(t.cum_row, t.row);
    t.time = t.wait = 0.0;
    t.count = t.nested = 0;
    t.row = -1;
  }
  for (auto &c : reg.counters)
  {
    c.cum = AccumulateValue(c.acc, c.cum, c.cum_set, c.value, c.set);
    c.cum_set = c.cum_set || c.set;
    c.value = 0.0;
    c.set = false;
  }
}

}  // namespace

namespace
{

int RegisterTally(std::string_view name, int fixed_parent)
{
  auto &reg = detail::registry;
  std::string key(name);
  auto it = reg.tally_ids.find(key);
  if (it != reg.tally_ids.end())
  {
    return it->second;
  }
  const int id = static_cast<int>(reg.tallies.size());
  reg.tallies.push_back({key, fixed_parent});
  reg.tally_ids.emplace(std::move(key), id);
  return id;
}

}  // namespace

int RegisterTally(std::string_view name)
{
  return RegisterTally(name, -1);
}

int RegisterPhaseTally(std::string_view name)
{
  auto &reg = detail::registry;
  const int id = RegisterTally(name, Timer::TOTAL);
  if (reg.tallies[id].phase_base < 0)
  {
    // The sub-tallies of all rows, contiguous: the sub-tally of a row is found by index.
    const int base = static_cast<int>(reg.tallies.size());
    for (int row = 0; row < Timer::NUM_TIMINGS; row++)
    {
      MFEM_VERIFY(RegisterTally(fmt::format("{}/{}", name, Timer::names[row]), row) ==
                      base + row,
                  "Phase tally rows must not be registered before the phase tally!");
    }
    reg.tallies[id].phase_base = base;
  }
  return id;
}

int RegisterCounter(std::string_view path, Reduce reduce, Accumulate acc)
{
  auto &reg = detail::registry;
  std::string key(path);
  auto it = reg.counter_ids.find(key);
  if (it != reg.counter_ids.end())
  {
    return it->second;
  }
  const int id = static_cast<int>(reg.counters.size());
  reg.counters.push_back({key, reduce, acc});
  reg.counter_ids.emplace(std::move(key), id);
  return id;
}

void SetInfo(std::string_view path, const json &value)
{
  At(state.info, path) = value;
  state.info_paths.emplace(path);
}

std::string_view CurrentRole()
{
  const auto &roles = detail::registry.roles;
  return roles.empty() ? std::string_view("Main") : std::string_view(roles.back());
}

ScopedRole::ScopedRole(std::string_view role)
{
  detail::registry.roles.emplace_back(role);
}

ScopedRole::~ScopedRole()
{
  detail::registry.roles.pop_back();
}

std::string LinearSolverPath(std::string_view role, std::string_view field)
{
  if (role == "Main" && (field == "Solves" || field == "Its"))
  {
    return (field == "Solves") ? "LinearSolver/TotalSolves" : "LinearSolver/TotalIts";
  }
  return fmt::format("LinearSolver/{}/{}", role, field);
}

const Results &GetResults()
{
  MFEM_VERIFY(state.finalized, "BlockTimer::Finalize() must be called before accessing "
                               "performance results!");
  return state.results;
}

void WriteJson(json &meta, const json *error)
{
  const auto &r = GetResults();
  constexpr double to_mb = 1.0 / (1024.0 * 1024.0);

  // Remove everything written by a previous iteration.
  for (const char *key : {"ElapsedTime", "LinearSolver", "PeakMemoryGrowthMegabytes",
                          "PeakNodeMemoryGrowthMegabytes", "PeakMemoryMegabytes",
                          "PeakNodeMemoryMegabytes", "HeapGrowthMegabytes", "HeapMegabytes",
                          "Cumulative", "Complete", "Error", "RootProcess"})
  {
    meta.erase(key);
  }
  // Counters by the union over processes of their paths, now and when last written: a
  // counter registered only on other processes is written by the root all the same.
  state.counter_paths.insert(r.counter_paths.begin(), r.counter_paths.end());
  for (const auto &c : detail::registry.counters)
  {
    state.counter_paths.insert(c.path);
  }
  for (const auto &p : state.counter_paths)
  {
    ErasePath(meta, p);
  }
  for (const auto &p : state.info_paths)
  {
    ErasePath(meta, p);
  }

  meta["SchemaVersion"] = SCHEMA_VERSION;
  meta["Complete"] = (error == nullptr);
  if (error)
  {
    meta["Error"] = *error;
    meta["Error"]["Iteration"] = r.iteration;
  }

  // The values of the root process alone (failure report without a collective) are not
  // statistics over processes: they go under RootProcess instead of the usual keys.
  json root_only;
  json &out = (error && !error->value("Collective", true)) ? root_only : meta;
  out["ElapsedTime"] = ElapsedTimeJson(r.iter);
  for (int i = 0; i < Timer::NUM_TIMINGS; i++)
  {
    std::string key(Timer::names[i]);
    const auto &rm = r.iter.rank_mem[i];
    const auto &nm = r.iter.node_mem[i];
    out["PeakMemoryGrowthMegabytes"]["Min"][key] = rm.min * to_mb;
    out["PeakMemoryGrowthMegabytes"]["Max"][key] = rm.max * to_mb;
    out["PeakMemoryGrowthMegabytes"]["Sum"][key] = rm.sum * to_mb;
    out["PeakMemoryGrowthMegabytes"]["MaxRank"][key] = rm.max_rank;
    out["PeakNodeMemoryGrowthMegabytes"]["Min"][key] = nm.min * to_mb;
    out["PeakNodeMemoryGrowthMegabytes"]["Max"][key] = nm.max * to_mb;
    out["PeakNodeMemoryGrowthMegabytes"]["Sum"][key] = nm.sum * to_mb;
    out["PeakNodeMemoryGrowthMegabytes"]["MaxNode"][key] = nm.max_rank;
  }
  for (const auto &[key, rm, nm] :
       {std::tuple{"Baseline", r.mem_baseline, r.node_mem_baseline},
        std::tuple{"Unattributed", r.mem_unattributed, r.node_mem_unattributed}})
  {
    out["PeakMemoryGrowthMegabytes"]["Min"][key] = rm.min * to_mb;
    out["PeakMemoryGrowthMegabytes"]["Max"][key] = rm.max * to_mb;
    out["PeakMemoryGrowthMegabytes"]["Sum"][key] = rm.sum * to_mb;
    out["PeakMemoryGrowthMegabytes"]["MaxRank"][key] = rm.max_rank;
    out["PeakNodeMemoryGrowthMegabytes"]["Min"][key] = nm.min * to_mb;
    out["PeakNodeMemoryGrowthMegabytes"]["Max"][key] = nm.max * to_mb;
    out["PeakNodeMemoryGrowthMegabytes"]["Sum"][key] = nm.sum * to_mb;
    out["PeakNodeMemoryGrowthMegabytes"]["MaxNode"][key] = nm.max_rank;
  }
  out["PeakMemoryMegabytes"] = IterationPeakJson(r.peak_mem, false, r.peak_lower_bound > 0);
  out["PeakNodeMemoryMegabytes"] =
      IterationPeakJson(r.peak_node_mem, true, r.peak_lower_bound > 0);
  auto &heap = out["HeapGrowthMegabytes"];
  auto HeapRow = [&heap](const std::string &key, const Stats &s)
  {
    heap["Min"][key] = s.min * to_mb;
    heap["Max"][key] = s.max * to_mb;
    heap["Sum"][key] = s.sum * to_mb;
    heap["MaxRank"][key] = s.max_rank;
  };
  for (int i = 0; i < Timer::NUM_TIMINGS; i++)
  {
    HeapRow(std::string(Timer::names[i]), r.heap_growth[i]);
  }
  HeapRow("Baseline", r.heap_baseline);
  HeapRow("Unattributed", r.heap_unattributed);
  for (int i = 0; i < Timer::NUM_TIMINGS; i++)
  {
    if (r.heap_coverage[i].n > 0)
    {
      heap["Coverage"][std::string(Timer::names[i])] = r.heap_coverage[i].min;
    }
  }
  out["HeapMegabytes"] = PeakJson(r.heap_mem, false);

  json cum;
  cum["ElapsedTime"] = ElapsedTimeJson(r.cum);
  cum["PeakMemoryMegabytes"] = PeakJson(r.cum_peak_mem, false);
  cum["PeakNodeMemoryMegabytes"] = PeakJson(r.cum_peak_node_mem, true);
  for (std::size_t i = 0; i < r.counter_paths.size(); i++)
  {
    if (r.counter[i].n > 0)
    {
      At(out, r.counter_paths[i]) = CounterJson(r.counter_reduce[i], r.counter[i]);
    }
    if (r.counter_acc[i] != Accumulate::NONE && r.cum_counter[i].n > 0)
    {
      At(cum, r.counter_paths[i]) = CounterJson(r.counter_reduce[i], r.cum_counter[i]);
    }
  }
  if (!state.info.empty())
  {
    out.merge_patch(state.info);
  }
  out["Cumulative"] = std::move(cum);
  if (&out == &root_only)
  {
    meta["RootProcess"] = std::move(root_only);
  }
}

}  // namespace perf

bool BlockTimer::IsFinalized()
{
  return perf::state.finalized;
}

void BlockTimer::Finalize(MPI_Comm comm)
{
  using namespace perf;
  auto &reg = detail::registry;

  reg.suspended = true;

  // Drain any open timers.
  while (!stack.empty())
  {
    timer.MarkTime(stack.back());
    timer.MarkMemory(stack.back());
    timer.MarkHeap(stack.back(), HeapGapRow(), true);
    stack.pop_back();
    HeapEvent(true);
  }
  timer.MarkIterationEnd();
  timer.MarkMemory(Timer::TOTAL, timer.MemoryFromStart());
  timer.MarkHeap(Timer::NUM_TIMINGS, HeapGapRow(), true);
  HeapEvent(true);
  Set(RegisterCounter("Memory/Heap/QuerySeconds", Reduce::STATS, Accumulate::SUM),
      timer.HeapQueryTime());
  Set(RegisterCounter("Memory/Heap/Samples", Reduce::MAX, Accumulate::SUM),
      static_cast<double>(timer.HeapSamples()));
  Set(RegisterCounter("Memory/Heap/SkippedSamples", Reduce::MAX, Accumulate::SUM),
      static_cast<double>(timer.HeapSkipped()));
  SetInfo("Memory/HeapSource", memory_reporting::HeapMemorySource());
  SetInfo("Communication/MPIWaitTally", detail::MpiWaitTally());
  if (detail::MpiWaitTally())
  {
    SetInfo("Communication/MPIWaitTool", detail::mpi_wait_tool);
  }

  const int rank = Mpi::Rank(comm);
  constexpr int N = Timer::NUM_TIMINGS;
  auto [tidx, tkeys, cidx, ckeys] = CanonicalOrder(comm);
  const int T = static_cast<int>(tidx.size()), C = static_cast<int>(cidx.size());
  const auto [iter_peak, peak_exact] = IterationPeakMemory(timer.RssSampled());
  const double cur_mem = static_cast<double>(memory_reporting::GetCurrentMemory());
  const double cum_peak = std::max(state.cum_peak, iter_peak);
  double phase_mem = 0.0;
  for (int i = 0; i < N; i++)
  {
    if (i != Timer::TOTAL)
    {
      phase_mem += static_cast<double>(timer.MemoryData(static_cast<Timer::Index>(i)));
    }
  }
  const double total_mem = static_cast<double>(timer.MemoryData(Timer::TOTAL));
  const double baseline_mem = iter_peak - total_mem, unattr_mem = total_mem - phase_mem;

  // Pack everything into one buffer for a single reduction.
  std::vector<double> buf;
  buf.reserve(STRIDE * (7 * N + 14 * T + 2 * C + 10));
  for (int i = 0; i < N; i++)
  {
    PushEntry(buf, timer.Data(static_cast<Timer::Index>(i)), true, rank);
  }
  for (int i = 0; i < N; i++)
  {
    PushEntry(buf, timer.Counts(static_cast<Timer::Index>(i)), true, rank);
  }
  for (int i = 0; i < N; i++)
  {
    PushEntry(buf, static_cast<double>(timer.MemoryData(static_cast<Timer::Index>(i))),
              true, rank);
  }
  for (int i = 0; i < T; i++)
  {
    PushEntry(buf, (tidx[i] >= 0) ? reg.tallies[tidx[i]].time : 0.0, true, rank);
  }
  for (int i = 0; i < T; i++)
  {
    PushEntry(buf, (tidx[i] >= 0) ? reg.tallies[tidx[i]].count : 0.0, true, rank);
  }
  for (int i = 0; i < T; i++)
  {
    const auto *t = (tidx[i] >= 0) ? &reg.tallies[tidx[i]] : nullptr;
    PushEntry(buf, t ? t->wait : 0.0, true, rank);
    PushEntry(buf, t ? t->time - t->wait : 0.0, true, rank);
    PushEntry(buf, t ? t->nested : 0.0, true, rank);
    PushEntry(buf, (t && t->count > 0) ? 1.0 : 0.0, true, rank);
    const int row = !t ? -1 : (t->fixed_parent >= 0) ? t->fixed_parent : t->row;
    PushEntry(buf, row, row >= 0, rank);
  }
  for (int i = 0; i < C; i++)
  {
    const auto *c = (cidx[i] >= 0) ? &reg.counters[cidx[i]] : nullptr;
    PushEntry(buf, c ? c->value : 0.0, c && c->set, rank);
  }
  PushEntry(buf, iter_peak, true, rank);
  PushEntry(buf, cur_mem, true, rank);
  PushEntry(buf, peak_exact ? 0.0 : 1.0, true, rank);
  for (int i = 0; i < N; i++)
  {
    PushEntry(buf, timer.CumulativeData(static_cast<Timer::Index>(i)), true, rank);
  }
  for (int i = 0; i < N; i++)
  {
    PushEntry(buf, timer.CumulativeCounts(static_cast<Timer::Index>(i)), true, rank);
  }
  for (int i = 0; i < T; i++)
  {
    const auto *t = (tidx[i] >= 0) ? &reg.tallies[tidx[i]] : nullptr;
    PushEntry(buf, t ? t->cum_time + t->time : 0.0, true, rank);
  }
  for (int i = 0; i < T; i++)
  {
    const auto *t = (tidx[i] >= 0) ? &reg.tallies[tidx[i]] : nullptr;
    PushEntry(buf, t ? t->cum_count + t->count : 0.0, true, rank);
  }
  for (int i = 0; i < T; i++)
  {
    const auto *t = (tidx[i] >= 0) ? &reg.tallies[tidx[i]] : nullptr;
    const double wait = t ? t->cum_wait + t->wait : 0.0;
    PushEntry(buf, wait, true, rank);
    PushEntry(buf, t ? t->cum_time + t->time - wait : 0.0, true, rank);
    PushEntry(buf, t ? t->cum_nested + t->nested : 0.0, true, rank);
    PushEntry(buf, (t && t->cum_count + t->count > 0) ? 1.0 : 0.0, true, rank);
    const int row = !t                       ? -1
                    : (t->fixed_parent >= 0) ? t->fixed_parent
                                             : detail::MergeRow(t->cum_row, t->row);
    PushEntry(buf, row, row >= 0, rank);
  }
  for (int i = 0; i < C; i++)
  {
    const auto *c = (cidx[i] >= 0) ? &reg.counters[cidx[i]] : nullptr;
    PushEntry(buf, c ? AccumulateValue(c->acc, c->cum, c->cum_set, c->value, c->set) : 0.0,
              c && (c->cum_set || c->set), rank);
  }
  PushEntry(buf, cum_peak, true, rank);
  PushEntry(buf, baseline_mem, true, rank);
  PushEntry(buf, unattr_mem, true, rank);
  for (int i = 0; i < N; i++)
  {
    PushEntry(buf,
              static_cast<double>((i == Timer::TOTAL)
                                      ? timer.HeapAtLastSample() - timer.HeapAtStart()
                                      : timer.HeapData(i)),
              true, rank);
  }
  PushEntry(buf, static_cast<double>(timer.HeapAtStart()), true, rank);
  PushEntry(buf, static_cast<double>(timer.HeapData(N)), true, rank);
  PushEntry(buf, static_cast<double>(timer.HeapAtLastSample()), true, rank);
  for (int i = 0; i < N; i++)
  {
    const double coverage = timer.HeapCoverage(i);
    PushEntry(buf, coverage, coverage >= 0.0, rank);
  }
  StatsAllreduce(buf, comm);

  // Node memory: sum over the ranks of each node, then statistics over nodes. A process
  // alone (failure report of the root) does not split (or replace the cached) node
  // communicators.
  static const memory_reporting::NodeComms self_comms{MPI_COMM_SELF, MPI_COMM_SELF, 0};
  const auto &nc =
      (comm == MPI_COMM_SELF) ? self_comms : memory_reporting::GetNodeComms(comm);
  std::vector<double> node_local(N + 4), node_sum(N + 4, 0.0);
  for (int i = 0; i < N; i++)
  {
    node_local[i] = static_cast<double>(timer.MemoryData(static_cast<Timer::Index>(i)));
  }
  node_local[N] = iter_peak;
  node_local[N + 1] = cur_mem;
  node_local[N + 2] = baseline_mem;
  node_local[N + 3] = unattr_mem;
  MPI_Reduce(node_local.data(), node_sum.data(), N + 4, MPI_DOUBLE, MPI_SUM, 0, nc.node);
  std::vector<double> node_buf;
  if (nc.leaders != MPI_COMM_NULL)
  {
    const int node_idx = Mpi::Rank(nc.leaders);
    for (int i = 0; i < N + 4; i++)
    {
      PushEntry(node_buf, node_sum[i], true, node_idx);
    }
    PushEntry(node_buf, std::max(state.cum_node_peak, node_sum[N]), true, node_idx);
    PushEntry(node_buf, node_sum[N], true, node_idx);  // Local node peak, not reduced
    std::vector<double> tmp(node_buf.begin(), node_buf.end() - STRIDE);
    StatsAllreduce(tmp, nc.leaders);
    std::copy(tmp.begin(), tmp.end(), node_buf.begin());
  }
  else
  {
    node_buf.resize(STRIDE * (N + 6));
  }
  MPI_Bcast(node_buf.data(), static_cast<int>(node_buf.size()), MPI_DOUBLE, 0, nc.node);

  // Unpack.
  auto &r = state.results;
  r = Results();
  r.iteration = state.iteration;

  const double *e = buf.data();
  auto Next = [&e]()
  {
    auto s = Unpack(e);
    e += STRIDE;
    return s;
  };
  auto NextN = [&Next](std::vector<Stats> &v, int n)
  {
    v.resize(n);
    for (auto &s : v)
    {
      s = Next();
    }
  };
  auto NextTallyExtras = [&Next](TimerResults &t, int n)
  {
    for (auto *v : {&t.tally_wait, &t.tally_work, &t.tally_nested, &t.tally_active})
    {
      v->resize(n);
    }
    t.tally_parents.resize(n);
    for (int i = 0; i < n; i++)
    {
      t.tally_wait[i] = Next();
      t.tally_work[i] = Next();
      t.tally_nested[i] = Next();
      t.tally_active[i] = Next();
      // The row of all processes which measured the tally, TOTAL if they differ.
      const auto row = Next();
      t.tally_parents[i] = (row.n > 0 && row.min == row.max)
                               ? static_cast<Timer::Index>(row.min)
                               : Timer::TOTAL;
    }
  };
  NextN(r.iter.time, N);
  NextN(r.iter.count, N);
  NextN(r.iter.rank_mem, N);
  NextN(r.iter.tally_time, T);
  NextN(r.iter.tally_count, T);
  NextTallyExtras(r.iter, T);
  NextN(r.counter, C);
  r.peak_mem = Next();
  r.cur_mem = Next();
  r.peak_lower_bound = static_cast<int>(Next().sum);
  NextN(r.cum.time, N);
  NextN(r.cum.count, N);
  NextN(r.cum.tally_time, T);
  NextN(r.cum.tally_count, T);
  NextTallyExtras(r.cum, T);
  NextN(r.cum_counter, C);
  r.cum_peak_mem = Next();
  r.mem_baseline = Next();
  r.mem_unattributed = Next();
  NextN(r.heap_growth, N);
  r.heap_baseline = Next();
  r.heap_unattributed = Next();
  r.heap_mem = Next();
  NextN(r.heap_coverage, N);
  e = node_buf.data();
  NextN(r.iter.node_mem, N);
  r.peak_node_mem = Next();
  r.cur_node_mem = Next();
  r.node_mem_baseline = Next();
  r.node_mem_unattributed = Next();
  r.cum_peak_node_mem = Next();
  state.last_node_peak = e[2];  // This node's own peak (sum of its ranks)
  r.num_nodes = r.peak_node_mem.n;
  r.cum.rank_mem = r.iter.rank_mem;
  r.cum.node_mem = r.iter.node_mem;

  r.iter.tally_names = tkeys;
  r.cum.tally_names = tkeys;
  reg.suspended = false;
  for (int i = 0; i < C; i++)
  {
    const auto sep = ckeys[i].rfind('\x1f');
    r.counter_paths.push_back(ckeys[i].substr(0, sep));
    r.counter_reduce.push_back(static_cast<Reduce>(ckeys[i][sep + 1] - '0'));
    r.counter_acc.push_back(static_cast<Accumulate>(ckeys[i][sep + 2] - '0'));
  }
  state.finalized = true;
}

void BlockTimer::StartIteration()
{
  using namespace perf;
  state.cum_peak = std::max(state.cum_peak, IterationPeakMemory(timer.RssSampled()).first);
  state.cum_node_peak = std::max(state.cum_node_peak, state.last_node_peak);
  Commit();
  state.info = json::object();
  state.peak_at_start = memory_reporting::GetPeakMemory();
  state.cur_at_start = memory_reporting::GetCurrentMemory();
  timer.StartIteration();
  heap_open = stack.size();
  state.iteration++;
  detail::iteration_serial++;
  state.finalized = false;
}

void BlockTimer::Reset()
{
  using namespace perf;
  timer = Timer();
  stack.clear();
  heap_open = 0;
  started = false;
  // Registrations are kept (ids may be cached by callers); values are cleared.
  for (auto &t : detail::registry.tallies)
  {
    t.time = t.cum_time = t.wait = t.cum_wait = 0.0;
    t.count = t.cum_count = t.nested = t.cum_nested = 0;
    t.row = t.cum_row = -1;
  }
  for (auto &c : detail::registry.counters)
  {
    c.value = c.cum = 0.0;
    c.set = c.cum_set = false;
  }
  detail::registry.roles.clear();
  detail::iteration_serial++;
  state = State();
}

void BlockTimer::Print(MPI_Comm comm)
{
  const auto &r = perf::GetResults();
  constexpr int w = 12;  // Data column width
  constexpr int h = 26;  // Left-hand side width
  using Row = std::tuple<std::string, std::string, std::string>;
  auto Table = [&](const std::string &title, const std::string &col1,
                   const std::string &col2, const std::string &col3, auto &&row_fn,
                   const Row *baseline = nullptr)
  {
    Mpi::Print(comm, "\n{:<{}s}{:>{}s}{:>{}s}{:>{}s}\n", title, h, col1, w, col2, w, col3,
               w);
    Mpi::Print(comm, "{}\n", std::string(h + 3 * w, '='));
    if (baseline)
    {
      const auto &[v1, v2, v3] = *baseline;
      Mpi::Print(comm, "{:<{}s}{:>{}s}{:>{}s}{:>{}s}\n", "Baseline", h, v1, w, v2, w, v3,
                 w);
    }
    for (int i = 0; i < Timer::NUM_TIMINGS; i++)
    {
      if (r.cum.count[i].max > 0 || i == Timer::TOTAL)
      {
        if (i == Timer::TOTAL)
        {
          Mpi::Print(comm, "{}\n", std::string(h + 3 * w, '-'));
        }
        auto [v1, v2, v3] = row_fn(i);
        Mpi::Print(comm, "{:<{}s}{:>{}s}{:>{}s}{:>{}s}\n", Timer::descriptions[i], h, v1, w,
                   v2, w, v3, w);
      }
    }
  };
  auto TimeRow = [](const perf::TimerResults &t)
  {
    return [&t](int i) -> Row
    {
      return {fmt::format("{:.3f}", t.time[i].min), fmt::format("{:.3f}", t.time[i].max),
              fmt::format("{:.3f}", t.time[i].avg)};
    };
  };

  // Timing table: totals over all iterations (unchanged layout).
  if (r.iteration > 1)
  {
    Mpi::Print(comm,
               "\nPerformance report after iteration {:d}: elapsed times are totals over "
               "iterations 1-{:d}, peak memory is for iteration {:d} only\n",
               r.iteration, r.iteration, r.iteration);
  }
  Table("Elapsed Time Report (s)", "Min.", "Max.", "Avg.", TimeRow(r.cum));

  // Memory table for this iteration: per-node max and total across nodes of the peak RSS
  // growth of each row, and the running total high-water mark across nodes, from the peak
  // at the start of the iteration (Baseline) to the peak of the iteration (Total).
  const auto &nb = r.node_mem_baseline;
  double hwm = nb.sum;
  const Row baseline = {memory_reporting::FormatBytes(nb.max),
                        memory_reporting::FormatBytes(nb.sum),
                        memory_reporting::FormatBytes(nb.sum)};
  Table(
      "Peak Memory", "Per-Node", "Total", "Total HWM",
      [&](int i) -> Row
      {
        const auto &nm = r.iter.node_mem[i];
        if (i != Timer::TOTAL)
        {
          hwm += nm.sum;
        }
        return {memory_reporting::FormatBytes(nm.max),
                memory_reporting::FormatBytes(nm.sum),
                memory_reporting::FormatBytes((i == Timer::TOTAL) ? nb.sum + nm.sum : hwm)};
      },
      &baseline);
  if (r.peak_lower_bound > 0)
  {
    Mpi::Print(
        comm,
        "The peak RSS did not grow in this iteration on {:d} process{}: the largest RSS "
        "sampled is a lower bound of its peak\n",
        r.peak_lower_bound, (r.peak_lower_bound > 1) ? "es" : "");
  }

  // This iteration alone, when it differs from the totals.
  if (r.iteration > 1)
  {
    Table(fmt::format("Iteration {:d} (s)", r.iteration), "Min.", "Max.", "Avg.",
          TimeRow(r.iter));
  }

  // Inclusive tallies (totals over all iterations).
  bool any = false;
  for (std::size_t i = 0; i < r.cum.tally_names.size() && !any; i++)
  {
    any = r.cum.tally_count[i].max > 0;
  }
  if (any)
  {
    constexpr int ht = 38;
    Mpi::Print(comm, "\n{:<{}s}{:>{}s}{:>{}s}{:>{}s}{:>{}s}\n", "Inclusive Tallies (s)", ht,
               "Min.", w, "Max.", w, "Avg.", w, "Count", w);
    Mpi::Print(comm, "{}\n", std::string(ht + 4 * w, '='));
    for (std::size_t i = 0; i < r.cum.tally_names.size(); i++)
    {
      if (r.cum.tally_count[i].max > 0)
      {
        const auto &t = r.cum.tally_time[i];
        Mpi::Print(comm, "{:<{}s}{:>{}.3f}{:>{}.3f}{:>{}.3f}{:>{}d}\n",
                   r.cum.tally_names[i], ht, t.min, w, t.max, w, t.avg, w,
                   static_cast<long>(r.cum.tally_count[i].max), w);
      }
    }
  }
}

}  // namespace palace
