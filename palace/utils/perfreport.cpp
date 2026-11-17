// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "perfreport.hpp"

#include <cmath>
#include <cstdint>
#include <set>
#include <string>
#include <tuple>
#include <fmt/format.h>
#include <nlohmann/json.hpp>
#include "utils/communication.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/timer.hpp"

namespace palace::perf
{

using json = nlohmann::json;

namespace
{

constexpr int N = Timer::NUM_TIMINGS;
constexpr double to_mb = 1.0 / (1024.0 * 1024.0);

struct State
{
  ReportResults results;
  std::set<std::string> written;  // Metric paths written at the top level last time
};

State &GetState()
{
  static State s;
  return s;
}

struct NodeComms
{
  MPI_Comm comm = MPI_COMM_NULL, node = MPI_COMM_NULL, leaders = MPI_COMM_NULL;
};

const NodeComms &GetNodeComms(MPI_Comm comm)
{
  static NodeComms nc = [comm]()
  {
    NodeComms c;
    c.comm = comm;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &c.node);
    MPI_Comm_split(comm, (Mpi::Rank(c.node) == 0) ? 0 : MPI_UNDEFINED, 0, &c.leaders);
    return c;
  }();
  int cmp;
  MPI_Comm_compare(comm, nc.comm, &cmp);
  MFEM_VERIFY(cmp == MPI_IDENT || cmp == MPI_CONGRUENT,
              "perf::Report must always be called on the same communicator");
  return nc;
}

// Rows which contribute on this process: used (time or count positive), and TOTAL.
std::vector<char> RowsSet(const Timer &t, bool cum)
{
  std::vector<char> s(N);
  for (int i = 0; i < N; i++)
  {
    const auto idx = static_cast<Timer::Index>(i);
    s[i] = (i == Timer::TOTAL) || (cum ? (t.RunData(idx) > 0.0 || t.RunCounts(idx) > 0)
                                       : (t.Data(idx) > 0.0 || t.Counts(idx) > 0));
  }
  return s;
}

// Fixed per-period values: N times, N counts, N memory growths, peak, baseline,
// unattributed. Memory values are also appended to node.
void PackFixed(const Timer &t, bool cum, std::vector<double> &v, std::vector<char> &s,
               std::vector<double> &node)
{
  const auto set = RowsSet(t, cum);
  for (int i = 0; i < N; i++)
  {
    const auto idx = static_cast<Timer::Index>(i);
    v.push_back(cum ? t.RunData(idx) : t.Data(idx));
    s.push_back(set[i]);
  }
  for (int i = 0; i < N; i++)
  {
    const auto idx = static_cast<Timer::Index>(i);
    v.push_back(static_cast<double>(cum ? t.RunCounts(idx) : t.Counts(idx)));
    s.push_back(set[i]);
  }
  const double peak = static_cast<double>(memory_reporting::GetPeakMemory());
  const double baseline =
      static_cast<double>(cum ? t.RunBaseline() : t.IterationBaseline());
  double rows = 0.0;
  for (int i = 0; i < N; i++)
  {
    const auto idx = static_cast<Timer::Index>(i);
    const double m = (i == Timer::TOTAL) ? 0.0
                                         : static_cast<double>(cum ? t.RunMemoryData(idx)
                                                                   : t.MemoryData(idx));
    rows += m;
    v.push_back(m);
    s.push_back(set[i] && i != Timer::TOTAL);
    node.push_back(m);
  }
  for (double m : {peak, baseline, peak - baseline - rows})
  {
    v.push_back(m);
    s.push_back(1);
    node.push_back(m);
  }
}

void UnpackFixed(const std::vector<Stats> &r, std::size_t &k, PeriodResults &p)
{
  p.time.assign(r.begin() + k, r.begin() + k + N);
  k += N;
  p.count.assign(r.begin() + k, r.begin() + k + N);
  k += N;
  p.mem.assign(r.begin() + k, r.begin() + k + N);
  k += N;
  p.peak = r[k++];
  p.baseline = r[k++];
  p.unattributed = r[k++];
}

void UnpackNode(const std::vector<Stats> &r, std::size_t &k, PeriodResults &p)
{
  p.node_mem.assign(r.begin() + k, r.begin() + k + N);
  k += N;
  p.node_peak = r[k++];
  p.node_baseline = r[k++];
  p.node_unattributed = r[k++];
}

json Num(double v)
{
  if (std::isfinite(v) && v == std::round(v) && std::abs(v) < 9.0e15)
  {
    return static_cast<std::int64_t>(v);
  }
  return v;
}

json StatsJson(const Stats &s, double scale = 1.0)
{
  return {{"Min", Num(s.min * scale)}, {"Avg", Num(s.Avg() * scale)},
          {"Max", Num(s.max * scale)}, {"Sum", Num(s.sum * scale)},
          {"MaxRank", s.max_rank},     {"Ranks", s.ranks}};
}

json NodeStatsJson(const Stats &s)
{
  return {{"Min", Num(s.min * to_mb)},
          {"Avg", Num(s.Avg() * to_mb)},
          {"Max", Num(s.max * to_mb)},
          {"Sum", Num(s.sum * to_mb)},
          {"Nodes", s.ranks}};
}

json::json_pointer Pointer(const std::string &path)
{
  return json::json_pointer("/" + path);
}

void ErasePath(json &meta, const std::string &path)
{
  auto ptr = Pointer(path);
  while (!ptr.empty() && meta.contains(ptr))
  {
    auto parent = ptr.parent_pointer();
    auto &obj = parent.empty() ? meta : meta[parent];
    obj.erase(ptr.back());
    if (!obj.empty())
    {
      break;
    }
    ptr = parent;
  }
}

void WritePeriod(json &out, const PeriodResults &p, bool cum,
                 std::set<std::string> *written)
{
  for (int i = 0; i < N; i++)
  {
    if (p.time[i].ranks > 0)
    {
      auto j = StatsJson(p.time[i]);
      j["Count"] = Num(p.count[i].max);
      out["ElapsedTime"][std::string(Timer::names[i])] = std::move(j);
    }
  }
  out["PeakMemoryMegabytes"] = StatsJson(p.peak, to_mb);
  out["PeakNodeMemoryMegabytes"] = NodeStatsJson(p.node_peak);
  auto &g = out["PeakMemoryGrowthMegabytes"];
  auto &ng = out["PeakNodeMemoryGrowthMegabytes"];
  g["Baseline"] = StatsJson(p.baseline, to_mb);
  ng["Baseline"] = NodeStatsJson(p.node_baseline);
  for (int i = 0; i < Timer::TOTAL; i++)
  {
    if (p.time[i].ranks > 0)
    {
      g[std::string(Timer::names[i])] = StatsJson(p.mem[i], to_mb);
      ng[std::string(Timer::names[i])] = NodeStatsJson(p.node_mem[i]);
    }
  }
  g["Unattributed"] = StatsJson(p.unattributed, to_mb);
  ng["Unattributed"] = NodeStatsJson(p.node_unattributed);
  const auto &metrics = GlobalMetrics();
  for (std::size_t i = 0; i < metrics.size(); i++)
  {
    const auto &m = metrics[i];
    const auto &s = p.metric[i];
    if (s.ranks == 0 || (cum && m.run == Run::ITERATION_ONLY))
    {
      continue;
    }
    json v;
    switch (m.red)
    {
      case Red::STATS:
        v = StatsJson(s);
        break;
      case Red::MAX:
        v = Num(s.max);
        break;
      case Red::MIN:
        v = Num(s.min);
        break;
      case Red::SUM:
        v = Num(s.sum);
        break;
      case Red::COUNT:
        MFEM_ABORT("Unresolved COUNT metric!");
    }
    out[Pointer(m.path)] = std::move(v);
    if (written)
    {
      written->insert(m.path);
    }
  }
}

}  // namespace

const Stats *ReportResults::Metric(std::string_view path, bool cumulative) const
{
  const auto &metrics = GlobalMetrics();
  const auto &p = cumulative ? cum : iter;
  for (std::size_t i = 0; i < metrics.size() && i < p.metric.size(); i++)
  {
    if (metrics[i].path == path)
    {
      return (p.metric[i].ranks > 0) ? &p.metric[i] : nullptr;
    }
  }
  return nullptr;
}

const ReportResults &Report(MPI_Comm comm)
{
  BlockTimer::Snapshot();
  Sync(comm);
  const auto &t = BlockTimer::GlobalTimer();
  std::vector<double> v, node;
  std::vector<char> s;
  PackFixed(t, false, v, s, node);
  PackMetrics(false, v, s);
  PackFixed(t, true, v, s, node);
  PackMetrics(true, v, s);
  const auto r = ReduceStats(v, s, comm);

  // Node memory: sums over the processes of each node, then statistics over nodes.
  const auto &nc = GetNodeComms(comm);
  std::vector<double> node_sum(node.size(), 0.0);
  MPI_Reduce(node.data(), node_sum.data(), static_cast<int>(node.size()), MPI_DOUBLE,
             MPI_SUM, 0, nc.node);
  std::vector<Stats> nr(node.size());
  if (nc.leaders != MPI_COMM_NULL)
  {
    nr = ReduceStats(node_sum, std::vector<char>(node.size(), 1), nc.leaders);
  }
  MPI_Bcast(nr.data(), static_cast<int>(nr.size() * sizeof(Stats)), MPI_BYTE, 0, nc.node);

  auto &res = GetState().results;
  res.iteration = t.Iterations();
  const std::size_t M = GlobalMetrics().size();
  std::size_t k = 0, kn = 0;
  UnpackFixed(r, k, res.iter);
  res.iter.metric.assign(r.begin() + k, r.begin() + k + M);
  k += M;
  UnpackFixed(r, k, res.cum);
  res.cum.metric.assign(r.begin() + k, r.begin() + k + M);
  UnpackNode(nr, kn, res.iter);
  UnpackNode(nr, kn, res.cum);
  return res;
}

void WriteJson(json &meta)
{
  auto &st = GetState();
  for (const char *key :
       {"MetadataVersion", "ElapsedTime", "PeakMemoryMegabytes", "PeakNodeMemoryMegabytes",
        "PeakMemoryGrowthMegabytes", "PeakNodeMemoryGrowthMegabytes", "Cumulative"})
  {
    meta.erase(key);
  }
  for (const auto &path : st.written)
  {
    ErasePath(meta, path);
  }
  st.written.clear();
  meta["MetadataVersion"] = METADATA_VERSION;
  WritePeriod(meta, st.results.iter, false, &st.written);
  json cum = json::object();
  WritePeriod(cum, st.results.cum, true, nullptr);
  meta["Cumulative"] = std::move(cum);
}

void Print(MPI_Comm comm)
{
  const auto &res = GetState().results;
  constexpr int w = 12, h = 26;
  using Row = std::tuple<std::string, std::string, std::string>;
  auto Table = [&](const std::string &title, const std::string &c1, const std::string &c2,
                   const std::string &c3, const PeriodResults &p, auto &&row_fn,
                   const Row *first = nullptr, const char *first_name = nullptr)
  {
    Mpi::Print(comm, "\n{:<{}s}{:>{}s}{:>{}s}{:>{}s}\n", title, h, c1, w, c2, w, c3, w);
    Mpi::Print(comm, "{}\n", std::string(h + 3 * w, '='));
    if (first)
    {
      const auto &[v1, v2, v3] = *first;
      Mpi::Print(comm, "{:<{}s}{:>{}s}{:>{}s}{:>{}s}\n", first_name, h, v1, w, v2, w, v3,
                 w);
    }
    for (int i = 0; i < N; i++)
    {
      if (p.time[i].ranks > 0)
      {
        if (i == Timer::TOTAL)
        {
          Mpi::Print(comm, "{}\n", std::string(h + 3 * w, '-'));
        }
        const auto [v1, v2, v3] = row_fn(i);
        Mpi::Print(comm, "{:<{}s}{:>{}s}{:>{}s}{:>{}s}\n", Timer::descriptions[i], h, v1, w,
                   v2, w, v3, w);
      }
    }
  };
  auto TimeRow = [](const PeriodResults &p)
  {
    return [&p](int i) -> Row
    {
      return {fmt::format("{:.3f}", p.time[i].min), fmt::format("{:.3f}", p.time[i].max),
              fmt::format("{:.3f}", p.time[i].Avg())};
    };
  };
  Table("Elapsed Time Report (s)", "Min.", "Max.", "Avg.", res.cum, TimeRow(res.cum));
  if (res.iteration > 1)
  {
    Table(fmt::format("Iteration {:d} (s)", res.iteration), "Min.", "Max.", "Avg.",
          res.iter, TimeRow(res.iter));
  }
  // Memory of this iteration: per-node max and total growth of each row, and the running
  // total high-water mark from the peak at the iteration start (Baseline).
  const auto &p = res.iter;
  double hwm = p.node_baseline.sum;
  const Row baseline = {memory_reporting::FormatBytes(p.node_baseline.max),
                        memory_reporting::FormatBytes(p.node_baseline.sum),
                        memory_reporting::FormatBytes(hwm)};
  Table(
      "Peak Memory", "Per-Node", "Total", "Total HWM", p,
      [&](int i) -> Row
      {
        if (i == Timer::TOTAL)
        {
          return {memory_reporting::FormatBytes(p.node_peak.max),
                  memory_reporting::FormatBytes(p.node_peak.sum),
                  memory_reporting::FormatBytes(p.node_peak.sum)};
        }
        hwm += p.node_mem[i].sum;
        return {memory_reporting::FormatBytes(p.node_mem[i].max),
                memory_reporting::FormatBytes(p.node_mem[i].sum),
                memory_reporting::FormatBytes(hwm)};
      },
      &baseline, "Baseline");
}

}  // namespace palace::perf
