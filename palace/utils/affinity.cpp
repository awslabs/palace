// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "affinity.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "utils/communication.hpp"

#if defined(__linux__)
#include <sched.h>
#endif

namespace palace::utils
{

namespace
{

// CPU bound before rebinding (or -1 if unchanged).
int rebind_original_cpu = -1;

#if defined(__linux__)

int GetEnvInt(const char *name)
{
  const char *s = std::getenv(name);
  return (s && *s) ? std::atoi(s) : -1;
}

// Parse a Linux list string like "0-95,192-287".
std::vector<int> ParseCpuList(const std::string &s)
{
  std::vector<int> cpus;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ','))
  {
    if (tok.empty() || tok == "\n")
    {
      continue;
    }
    const auto dash = tok.find('-');
    const int a = std::atoi(tok.substr(0, dash).c_str());
    const int b = (dash == std::string::npos) ? a : std::atoi(tok.substr(dash + 1).c_str());
    for (int c = a; c <= b; c++)
    {
      cpus.push_back(c);
    }
  }
  return cpus;
}

// CPU lists of the online NUMA domains with CPUs, in node order.
std::vector<std::vector<int>> GetNumaCpuLists()
{
  std::vector<std::vector<int>> domains;
  std::ifstream fo("/sys/devices/system/node/online");
  std::string online;
  if (!fo || !std::getline(fo, online))
  {
    return domains;
  }
  for (int node : ParseCpuList(online))
  {
    std::ifstream f("/sys/devices/system/node/node" + std::to_string(node) + "/cpulist");
    std::string s;
    if (f && std::getline(f, s))
    {
      auto cpus = ParseCpuList(s);
      if (!cpus.empty())
      {
        std::sort(cpus.begin(), cpus.end());
        domains.push_back(std::move(cpus));
      }
    }
  }
  return domains;
}

int GetBoundCpu()
{
  cpu_set_t set;
  CPU_ZERO(&set);
  if (sched_getaffinity(0, sizeof(set), &set) != 0 || CPU_COUNT(&set) != 1)
  {
    return -1;
  }
  for (int c = 0; c < CPU_SETSIZE; c++)
  {
    if (CPU_ISSET(c, &set))
    {
      return c;
    }
  }
  return -1;
}

bool BindToCpu(int cpu)
{
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(cpu, &set);
  return sched_setaffinity(0, sizeof(set), &set) == 0;
}

#endif

}  // namespace

void SpreadPackedProcessBinding()
{
#if defined(__linux__)
  if (const char *s = std::getenv("PALACE_DISABLE_REBIND");
      s && *s && std::string(s) != "0")
  {
    return;
  }
  int local_rank = GetEnvInt("OMPI_COMM_WORLD_LOCAL_RANK");
  int local_size = GetEnvInt("OMPI_COMM_WORLD_LOCAL_SIZE");
  if (local_rank < 0 || local_size < 0)
  {
    local_rank = GetEnvInt("MPI_LOCALRANKID");
    local_size = GetEnvInt("MPI_LOCALNRANKS");
  }
  if (local_rank < 0 || local_size < 2 || local_rank >= local_size)
  {
    return;
  }
  const int cpu = GetBoundCpu();
  if (cpu != local_rank)
  {
    return;  // Not the packed one-core-per-rank placement.
  }
  const auto domains = GetNumaCpuLists();
  const int num_domains = static_cast<int>(domains.size());
  if (num_domains < 2 || local_size < num_domains)
  {
    return;
  }

  // Ranks per domain for the current packed placement (cores 0, ..., local_size - 1) and
  // for an even block distribution over the domains.
  std::vector<int> packed(num_domains, 0);
  for (int d = 0; d < num_domains; d++)
  {
    packed[d] =
        static_cast<int>(std::count_if(domains[d].begin(), domains[d].end(),
                                       [local_size](int c) { return c < local_size; }));
  }
  const auto [min_it, max_it] = std::minmax_element(packed.begin(), packed.end());
  if (*max_it - *min_it <= 1)
  {
    return;  // Already balanced.
  }
  auto first = [&](int d)
  { return static_cast<int>((long long)d * local_size / num_domains); };
  int d = 0;
  while (local_rank >= first(d + 1))
  {
    d++;
  }
  const int i = local_rank - first(d);
  for (int k = 0; k < num_domains; k++)
  {
    if (first(k + 1) - first(k) > static_cast<int>(domains[k].size()))
    {
      return;  // Not enough cores (should not happen for packed single-core binding).
    }
  }
  const int target = domains[d][i];
  if (target != cpu && BindToCpu(target))
  {
    rebind_original_cpu = cpu;
  }
#endif
}

void CheckProcessBinding(MPI_Comm comm)
{
#if defined(__linux__)
  MPI_Comm node_comm;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
  const int node_size = Mpi::Size(node_comm);
  const int cpu = GetBoundCpu();
  std::vector<int> cpus(node_size);
  MPI_Allgather(&cpu, 1, MPI_INT, cpus.data(), 1, MPI_INT, node_comm);
  int changed = (rebind_original_cpu >= 0), any_changed = changed;
  MPI_Allreduce(MPI_IN_PLACE, &any_changed, 1, MPI_INT, MPI_MAX, node_comm);
  bool ok = true;
  if (any_changed)
  {
    // All ranks must be bound to distinct single cores.
    std::sort(cpus.begin(), cpus.end());
    ok =
        (cpus.front() >= 0) && (std::adjacent_find(cpus.begin(), cpus.end()) == cpus.end());
    if (!ok && changed)
    {
      BindToCpu(rebind_original_cpu);
      rebind_original_cpu = -1;
      changed = 0;
    }
  }
  MPI_Comm_free(&node_comm);
  int num_changed = changed;
  MPI_Allreduce(MPI_IN_PLACE, &num_changed, 1, MPI_INT, MPI_SUM, comm);
  if (num_changed > 0)
  {
    Mpi::Print(comm,
               "Rebound {:d} processes to spread the ranks of each node evenly across its "
               "NUMA domains (disable with PALACE_DISABLE_REBIND=1)\n",
               num_changed);
  }
#endif
}

}  // namespace palace::utils
