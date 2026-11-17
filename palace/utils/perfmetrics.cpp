// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "perfmetrics.hpp"

#include <array>
#include <deque>
#include <limits>
#include <string>
#include <unordered_map>
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

}  // namespace palace::perf
