// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "perfmetrics.hpp"

#include <array>
#include <deque>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <mfem.hpp>
#include "utils/communication.hpp"

namespace palace::perf
{

std::vector<Stats> ReduceStats(const std::vector<double> &value,
                               const std::vector<char> &set, MPI_Comm comm)
{
  MFEM_VERIFY(value.size() == set.size(), "ReduceStats: size mismatch!");
  const int n = static_cast<int>(value.size());
  const int rank = Mpi::Rank(comm);
  constexpr double inf = std::numeric_limits<double>::infinity();
  struct DoubleInt
  {
    double v;
    int r;
  };
  std::vector<double> mn(n), sm(2 * n);
  std::vector<DoubleInt> mx(n);
  for (int i = 0; i < n; i++)
  {
    mn[i] = set[i] ? value[i] : inf;
    mx[i] = {set[i] ? value[i] : -inf, rank};
    sm[i] = set[i] ? value[i] : 0.0;
    sm[n + i] = set[i] ? 1.0 : 0.0;
  }
  // Raw MPI (not Mpi::Global*), so that future tallies of Palace collectives do not
  // measure the report itself.
  MPI_Allreduce(MPI_IN_PLACE, mn.data(), n, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(MPI_IN_PLACE, mx.data(), n, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
  MPI_Allreduce(MPI_IN_PLACE, sm.data(), 2 * n, MPI_DOUBLE, MPI_SUM, comm);
  std::vector<Stats> out(n);
  for (int i = 0; i < n; i++)
  {
    const int ranks = static_cast<int>(sm[n + i]);
    if (ranks > 0)
    {
      out[i] = {mn[i], mx[i].v, sm[i], mx[i].r, ranks};
    }
  }
  return out;
}

namespace
{

// Top-level keys of palace.json owned by other writers (design section 6).
constexpr std::array<std::string_view, 11> owned_keys = {"MetadataVersion",
                                                         "ElapsedTime",
                                                         "PeakMemoryMegabytes",
                                                         "PeakNodeMemoryMegabytes",
                                                         "PeakMemoryGrowthMegabytes",
                                                         "PeakNodeMemoryGrowthMegabytes",
                                                         "Cumulative",
                                                         "Problem",
                                                         "SavedAdaptedMesh",
                                                         "Excitations",
                                                         "GitTag"};

struct LocalRegistry
{
  std::unordered_map<std::string, int> ids;
  std::vector<int> pending;  // Local ids registered since the last Sync
};

LocalRegistry &Local()
{
  static LocalRegistry r;
  return r;
}

struct RoleInfo
{
  std::string name;
  bool partitioned = false;
  bool comm_root = true;
};

struct RoleTable
{
  std::deque<RoleInfo> roles = {RoleInfo{"Unassigned"}};  // deque: stable references
  std::vector<int> stack;
};

RoleTable &Roles()
{
  static RoleTable t;
  return t;
}

int Intern(std::string_view name, bool partitioned, bool comm_root)
{
  auto &t = Roles();
  for (std::size_t i = 0; i < t.roles.size(); i++)
  {
    const auto &r = t.roles[i];
    if (r.name == name && r.partitioned == partitioned && r.comm_root == comm_root)
    {
      return static_cast<int>(i);
    }
  }
  t.roles.push_back({std::string(name), partitioned, comm_root});
  return static_cast<int>(t.roles.size() - 1);
}

}  // namespace

namespace detail
{

bool ValidPath(std::string_view path)
{
  if (path.empty() || path.front() == '/' || path.back() == '/' ||
      path.find("//") != std::string_view::npos ||
      path.find_first_of("\t\n") != std::string_view::npos)
  {
    return false;
  }
  const auto first = path.substr(0, path.find('/'));
  for (auto key : owned_keys)
  {
    if (first == key)
    {
      return false;
    }
  }
  return true;
}

int Register(std::string_view path, Acc acc, Red red, Run run)
{
  MFEM_VERIFY(red != Red::COUNT, "COUNT must be resolved before registration!");
  MFEM_VERIFY(ValidPath(path), "Invalid performance metric path \"" << path << "\"!");
  auto &local = Local();
  std::string key(path);
  if (auto it = local.ids.find(key); it != local.ids.end())
  {
    const auto &e = entries[it->second];
    MFEM_VERIFY(e.acc == acc && e.red == red && e.run == run,
                "Performance metric \"" << path << "\" registered with different rules!");
    return it->second;
  }
  const int id = static_cast<int>(entries.size());
  entries.push_back({key, acc, red, run});
  local.ids.emplace(std::move(key), id);
  local.pending.push_back(id);
  return id;
}

void FoldIteration()
{
  for (auto &e : entries)
  {
    if (e.run == Run::ACCUMULATE && e.iter_set)
    {
      e.run_value = e.run_set ? Combine(e.acc, e.run_value, e.iter) : e.iter;
      e.run_set = true;
    }
    e.iter = 0.0;
    e.iter_set = false;
  }
}

void ResetValues()
{
  for (auto &e : entries)
  {
    e.iter = e.run_value = 0.0;
    e.iter_set = e.run_set = false;
  }
}

}  // namespace detail

ScopedRole::ScopedRole(std::string_view name)
{
  Roles().stack.push_back(Intern(name, false, true));
}

ScopedRole::ScopedRole(std::string_view name, MPI_Comm comm)
{
  const bool root = (comm != MPI_COMM_NULL) && Mpi::Rank(comm) == 0;
  Roles().stack.push_back(Intern(name, true, root));
}

ScopedRole::ScopedRole(RoleId id)
{
  Roles().stack.push_back(id.index);
}

ScopedRole::~ScopedRole()
{
  Roles().stack.pop_back();
}

RoleId CurrentRoleId()
{
  const auto &s = Roles().stack;
  return {s.empty() ? 0 : s.back()};
}

const std::string &RoleName(RoleId id)
{
  return Roles().roles[id.index].name;
}

Metric::Metric(std::string_view path, Acc acc, Red red, Run run)
{
  bool drop = false;
  if (red == Red::COUNT)
  {
    const auto &role = Roles().roles[CurrentRoleId().index];
    red = role.partitioned ? Red::SUM : Red::MAX;
    drop = role.partitioned && !role.comm_root;
  }
  id = detail::Register(path, acc, red, run);
  if (drop)
  {
    id = -1;  // Registered (so the rules agree everywhere) but never set here
  }
}

namespace
{

struct GlobalRegistry
{
  std::vector<MetricInfo> list;
  std::unordered_map<std::string, int> index;
  std::vector<int> local_of_global;  // -1: not registered on this process
};

GlobalRegistry &Global()
{
  static GlobalRegistry g;
  return g;
}

std::string Rules(Acc acc, Red red, Run run)
{
  return {static_cast<char>('0' + static_cast<int>(acc)),
          static_cast<char>('0' + static_cast<int>(red)),
          static_cast<char>('0' + static_cast<int>(run))};
}

// Whether p conflicts with a path of paths: one is a component-wise prefix of the other.
bool PrefixConflict(const std::set<std::string> &paths, const std::string &p,
                    std::string &other)
{
  for (auto pos = p.find('/'); pos != std::string::npos; pos = p.find('/', pos + 1))
  {
    if (auto it = paths.find(p.substr(0, pos)); it != paths.end())
    {
      other = *it;
      return true;
    }
  }
  const auto child = p + "/";
  if (auto it = paths.lower_bound(child); it != paths.end() && it->starts_with(child))
  {
    other = *it;
    return true;
  }
  return false;
}

}  // namespace

namespace detail
{

std::vector<std::string> ValidateNew(const std::vector<MetricInfo> &global,
                                     const std::string &lines, std::string &error)
{
  std::map<std::string, std::string> rules_of;
  std::set<std::string> paths;
  for (const auto &m : global)
  {
    rules_of[m.path] = Rules(m.acc, m.red, m.run);
    paths.insert(m.path);
  }
  std::map<std::string, std::string> fresh;
  std::size_t start = 0;
  while (start < lines.size())
  {
    const auto end = lines.find('\n', start);
    const auto line = lines.substr(start, end - start);
    start = end + 1;
    const auto tab = line.find('\t');
    const auto path = line.substr(0, tab), rules = line.substr(tab + 1);
    auto known = rules_of.find(path);
    if (known == rules_of.end())
    {
      known = fresh.find(path);
      if (known == fresh.end())
      {
        fresh.emplace(path, rules);
        continue;
      }
    }
    if (known->second != rules)
    {
      error = "metric \"" + path +
              "\" is registered with different rules on different processes";
      return {};
    }
  }
  std::vector<std::string> out;
  for (const auto &[path, rules] : fresh)
  {
    std::string other;
    if (PrefixConflict(paths, path, other))
    {
      error = "metric \"" + path + "\" and \"" + other +
              "\" conflict (a metric path cannot contain another)";
      return {};
    }
    paths.insert(path);
    out.push_back(path + '\t' + rules + '\n');
  }
  return out;
}

}  // namespace detail

void Sync(MPI_Comm comm)
{
  auto &local = Local();
  auto &g = Global();
  const auto pending = std::exchange(local.pending, {});
  std::string mine;
  for (int id : pending)
  {
    const auto &e = detail::entries[id];
    mine += e.path + '\t' + Rules(e.acc, e.red, e.run) + '\n';
  }
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm);
  int len = static_cast<int>(mine.size());
  std::vector<int> lens(np, 0), displs(np, 0);
  MPI_Gather(&len, 1, MPI_INT, lens.data(), 1, MPI_INT, 0, comm);
  std::string all;
  if (rank == 0)
  {
    for (int i = 1; i < np; i++)
    {
      displs[i] = displs[i - 1] + lens[i - 1];
    }
    all.resize(displs[np - 1] + lens[np - 1]);
  }
  MPI_Gatherv(mine.data(), len, MPI_CHAR, all.data(), lens.data(), displs.data(), MPI_CHAR,
              0, comm);
  std::string error, fresh;
  if (rank == 0)
  {
    for (const auto &line : detail::ValidateNew(g.list, all, error))
    {
      fresh += line;
    }
  }
  int sizes[2] = {static_cast<int>(error.size()), static_cast<int>(fresh.size())};
  MPI_Bcast(sizes, 2, MPI_INT, 0, comm);
  error.resize(sizes[0]);
  fresh.resize(sizes[1]);
  MPI_Bcast(error.data(), sizes[0], MPI_CHAR, 0, comm);
  MPI_Bcast(fresh.data(), sizes[1], MPI_CHAR, 0, comm);
  MFEM_VERIFY(error.empty(), "Performance report: " << error << "!");
  std::size_t start = 0;
  while (start < fresh.size())
  {
    const auto end = fresh.find('\n', start);
    const auto line = fresh.substr(start, end - start);
    start = end + 1;
    const auto tab = line.find('\t');
    MetricInfo m{line.substr(0, tab), static_cast<Acc>(line[tab + 1] - '0'),
                 static_cast<Red>(line[tab + 2] - '0'),
                 static_cast<Run>(line[tab + 3] - '0')};
    g.index.emplace(m.path, static_cast<int>(g.list.size()));
    auto it = local.ids.find(m.path);
    g.local_of_global.push_back((it != local.ids.end()) ? it->second : -1);
    g.list.push_back(std::move(m));
  }
  // Metrics registered here which other processes had registered before.
  for (int id : pending)
  {
    g.local_of_global[g.index.at(detail::entries[id].path)] = id;
  }
}

const std::vector<MetricInfo> &GlobalMetrics()
{
  return Global().list;
}

void PackMetrics(bool cumulative, std::vector<double> &value, std::vector<char> &set)
{
  const auto &g = Global();
  for (std::size_t i = 0; i < g.list.size(); i++)
  {
    const int id = g.local_of_global[i];
    double v = 0.0;
    bool s = false;
    if (id >= 0)
    {
      const auto &e = detail::entries[id];
      if (!cumulative)
      {
        v = e.iter;
        s = e.iter_set;
      }
      else if (e.run == Run::ACCUMULATE)
      {
        s = e.run_set || e.iter_set;
        v = (e.run_set && e.iter_set) ? detail::Combine(e.acc, e.run_value, e.iter)
                                      : (e.run_set ? e.run_value : e.iter);
      }
    }
    value.push_back(v);
    set.push_back(s);
  }
}

}  // namespace palace::perf
