// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "perfreport.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <utility>
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
  std::set<std::string> written;  // Paths written at the top level last time
};

State &GetState()
{
  static State s;
  return s;
}

// Per period: time, count and memory growth of each row, then peak, baseline and
// unattributed; memory values are also appended to node. Rows unused on this process (time
// and count zero) are unset (NaN), and so is the memory of TOTAL.
void Pack(const Timer &t, bool cum, std::vector<double> &v, std::vector<double> &node)
{
  constexpr double unset = std::numeric_limits<double>::quiet_NaN();
  const double peak = static_cast<double>(t.ReportMemory());
  const double baseline = static_cast<double>(t.Baseline(cum));
  double rows = 0.0;
  for (int i = 0; i < N; i++)
  {
    const auto r = t.Get(static_cast<Timer::Index>(i), cum);
    const bool used = (i == Timer::TOTAL) || r.time.count() > 0.0 || r.count > 0;
    const double m = static_cast<double>(r.mem);
    rows += m;
    v.push_back(used ? r.time.count() : unset);
    v.push_back(used ? static_cast<double>(r.count) : unset);
    v.push_back((used && i != Timer::TOTAL) ? m : unset);
    node.push_back(m);
  }
  for (double m : {peak, baseline, peak - baseline - rows})
  {
    v.push_back(m);
    node.push_back(m);
  }
}

json Num(double v)
{
  if (std::isfinite(v) && v == std::round(v) && std::abs(v) < 9.0e15)
  {
    return static_cast<std::int64_t>(v);
  }
  return v;
}

// Statistics over processes, or over nodes.
json StatsJson(const Stats &s, double scale = 1.0, bool node = false)
{
  json j = {{"Min", Num(s.min * scale)},
            {"Avg", Num(s.Avg() * scale)},
            {"Max", Num(s.max * scale)},
            {"Sum", Num(s.sum * scale)}};
  if (node)
  {
    j["Nodes"] = s.ranks;
  }
  else
  {
    j["MaxRank"] = s.max_rank;
    j["Ranks"] = s.ranks;
  }
  return j;
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

// Write the values of a period, recording the written paths in written if not null.
void WritePeriod(json &out, const PeriodResults &p, bool cum,
                 std::set<std::string> *written)
{
  auto put = [&](const std::string &path, json v)
  {
    out[Pointer(path)] = std::move(v);
    if (written)
    {
      written->insert(path);
    }
  };
  auto put_mem = [&](const std::string &key, const Stats &s, const Stats &node)
  {
    put("PeakMemoryGrowthMegabytes/" + key, StatsJson(s, to_mb));
    put("PeakNodeMemoryGrowthMegabytes/" + key, StatsJson(node, to_mb, true));
  };
  put("PeakMemoryMegabytes", StatsJson(p.peak, to_mb));
  put("PeakNodeMemoryMegabytes", StatsJson(p.node_peak, to_mb, true));
  put_mem("Baseline", p.baseline, p.node_baseline);
  put_mem("Unattributed", p.unattributed, p.node_unattributed);
  for (int i = 0; i < N; i++)
  {
    if (p.time[i].ranks > 0)
    {
      const std::string name(Timer::names[i]);
      auto j = StatsJson(p.time[i]);
      j["Count"] = Num(p.count[i].max);
      put("ElapsedTime/" + name, std::move(j));
      if (i != Timer::TOTAL)
      {
        put_mem(name, p.mem[i], p.node_mem[i]);
      }
    }
  }
  const auto &metrics = GlobalMetrics();
  for (std::size_t i = 0; i < metrics.size(); i++)
  {
    const auto &m = metrics[i];
    const auto &s = p.metric[i];
    if (s.ranks == 0 || (cum && m.run == Run::ITERATION_ONLY))
    {
      continue;
    }
    switch (m.red)
    {
      case Red::STATS:
        put(m.path, StatsJson(s));
        break;
      case Red::MAX:
        put(m.path, Num(s.max));
        break;
      case Red::MIN:
        put(m.path, Num(s.min));
        break;
      case Red::SUM:
        put(m.path, Num(s.sum));
        break;
      case Red::COUNT:
        MFEM_ABORT("Unresolved COUNT metric!");
    }
  }
}

}  // namespace

const Stats *ReportResults::Metric(std::string_view path, bool cumulative) const
{
  const auto &p = cumulative ? cum : iter;
  const int i = GlobalIndex(path);
  return (i >= 0 && i < static_cast<int>(p.metric.size()) && p.metric[i].ranks > 0)
             ? &p.metric[i]
             : nullptr;
}

const ReportResults &Report(MPI_Comm comm)
{
  BlockTimer::Snapshot();
  Sync(comm);
  const auto &t = BlockTimer::GlobalTimer();
  std::vector<double> v, node;
  for (bool cum : {false, true})
  {
    Pack(t, cum, v, node);
    PackMetrics(cum, v);
  }
  const auto r = ReduceStats(v, comm);

  // Node memory: sums over the processes of each node, then statistics over nodes.
  const auto &nc = memory_reporting::GetNodeComms(comm);
  std::vector<double> node_sum(node.size(), 0.0);
  MPI_Reduce(node.data(), node_sum.data(), static_cast<int>(node.size()), MPI_DOUBLE,
             MPI_SUM, 0, nc.node);
  std::vector<Stats> nr(node.size());
  if (nc.leaders != MPI_COMM_NULL)
  {
    nr = ReduceStats(node_sum, nc.leaders);
  }
  MPI_Bcast(nr.data(), static_cast<int>(nr.size() * sizeof(Stats)), MPI_BYTE, 0, nc.node);

  auto &res = GetState().results;
  res.iteration = t.Iterations();
  const std::size_t M = GlobalMetrics().size();
  std::size_t k = 0, kn = 0;
  for (auto *p : {&res.iter, &res.cum})
  {
    p->time.resize(N);
    p->count.resize(N);
    p->mem.resize(N);
    p->node_mem.resize(N);
    for (int i = 0; i < N; i++)
    {
      p->time[i] = r[k++];
      p->count[i] = r[k++];
      p->mem[i] = r[k++];
      p->node_mem[i] = nr[kn++];
    }
    for (auto [s, ns] :
         {std::pair{&p->peak, &p->node_peak}, std::pair{&p->baseline, &p->node_baseline},
          std::pair{&p->unattributed, &p->node_unattributed}})
    {
      *s = r[k++];
      *ns = nr[kn++];
    }
    p->metric.assign(r.begin() + k, r.begin() + k + M);
    k += M;
  }
  return res;
}

void WriteJson(json &meta)
{
  auto &st = GetState();
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
  auto Line =
      [&](std::string_view name, std::string_view a, std::string_view b, std::string_view c)
  { Mpi::Print(comm, "{:<{}s}{:>{}s}{:>{}s}{:>{}s}\n", name, h, a, w, b, w, c, w); };
  auto Header = [&](std::string_view title, std::string_view a, std::string_view b,
                    std::string_view c)
  {
    Mpi::Print(comm, "\n");
    Line(title, a, b, c);
    Mpi::Print(comm, "{}\n", std::string(h + 3 * w, '='));
  };
  // Row i: TOTAL after a rule.
  auto Row = [&](int i, std::string_view a, std::string_view b, std::string_view c)
  {
    if (i == Timer::TOTAL)
    {
      Mpi::Print(comm, "{}\n", std::string(h + 3 * w, '-'));
    }
    Line(Timer::descriptions[i], a, b, c);
  };
  auto Times = [&](std::string_view title, const PeriodResults &p)
  {
    Header(title, "Min.", "Max.", "Avg.");
    for (int i = 0; i < N; i++)
    {
      const auto &s = p.time[i];
      if (s.ranks > 0)
      {
        Row(i, fmt::format("{:.3f}", s.min), fmt::format("{:.3f}", s.max),
            fmt::format("{:.3f}", s.Avg()));
      }
    }
  };
  Times("Elapsed Time Report (s)", res.cum);
  if (res.iteration > 1)
  {
    Times(fmt::format("Iteration {:d} (s)", res.iteration), res.iter);
  }
  // Memory of this iteration: per-node max and total growth of each row, and the running
  // total high-water mark from the peak at the iteration start (Baseline). TOTAL: the peak.
  using memory_reporting::FormatBytes;
  const auto &p = res.iter;
  double hwm = p.node_baseline.sum;
  Header("Peak Memory", "Per-Node", "Total", "Total HWM");
  Line("Baseline", FormatBytes(p.node_baseline.max), FormatBytes(p.node_baseline.sum),
       FormatBytes(hwm));
  for (int i = 0; i < N; i++)
  {
    if (p.time[i].ranks > 0)
    {
      const auto &s = (i == Timer::TOTAL) ? p.node_peak : p.node_mem[i];
      hwm = (i == Timer::TOTAL) ? s.sum : hwm + s.sum;
      Row(i, FormatBytes(s.max), FormatBytes(s.sum), FormatBytes(hwm));
    }
  }
}

}  // namespace palace::perf
