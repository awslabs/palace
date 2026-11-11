// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "memoryreporting.hpp"

#include <sys/resource.h>

#if defined(__APPLE__)
#include <mach/mach.h>
#elif defined(__linux__)
#include <algorithm>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#endif
#if defined(__GLIBC__)
#include <malloc.h>
#endif

#include "utils/communication.hpp"

namespace palace::memory_reporting
{

// To query current memory on Linux, we look at what the kernel
// reports for the process.
//
// See: man 5 proc_pid_status
// (https://man7.org/linux/man-pages/man5/proc_pid_status.5.html)
//
// We look at MaxRSS (Resident Set Size). This is not the most accurate
// way to measure memory for two reasons:
// 1. It contains memory shared with other processes,
// 2. It doesn't fully account for kernel-internal optimizations.
// For an application like Palace, MaxRSS is often accurate at percent level or
// below.
//
// For MacOS, this is very hard, but this thread suggests using phys_footprint:
// https://developer.apple.com/forums/thread/105088
//
// See: man footprint
long GetCurrentMemory()
{
#if defined(__APPLE__)
  task_vm_info_data_t vm_info;
  mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
  if (task_info(mach_task_self(), TASK_VM_INFO, reinterpret_cast<task_info_t>(&vm_info),
                &count) == KERN_SUCCESS)
  {
    return static_cast<long>(vm_info.phys_footprint);
  }
  return 0;
#elif defined(__linux__)
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line))
  {
    if (line.starts_with("VmRSS:"))
    {
      // Value is in kilobytes.
      return std::stol(line.substr(6)) * 1024;
    }
  }
  return 0;
#else
  return 0;
#endif
}

// For peak memory we rely on rusage, which is relatively portable.
// See: man 2 getrusage
// (https://man7.org/linux/man-pages/man2/getrusage.2.html for Linux)
//
// Note that Linux reports the result in kilobytes, MacOS in bytes.
long GetPeakMemory()
{
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) == 0)
  {
#if defined(__APPLE__)
    // macOS reports ru_maxrss in bytes.
    return usage.ru_maxrss;
#else
    // Linux reports ru_maxrss in kilobytes.
    return usage.ru_maxrss * 1024;
#endif
  }
  return 0;
}

void TrimHeap()
{
#if defined(__GLIBC__)
  malloc_trim(0);
#endif
}

namespace
{

long available_memory_for_testing = -1;

#if defined(__linux__)
// The first whitespace-separated word of a file, or nothing if unreadable.
std::optional<std::string> ReadWord(const std::string &file)
{
  std::ifstream in(file);
  std::string word;
  if (!(in >> word))
  {
    return std::nullopt;
  }
  return word;
}

std::optional<long long> ReadNumber(const std::string &file)
{
  const auto word = ReadWord(file);
  if (!word)
  {
    return std::nullopt;
  }
  try
  {
    return std::stoll(*word);
  }
  catch (const std::exception &)
  {
    return std::nullopt;
  }
}

// The value of a key of a memory.stat file, or nothing if unreadable.
std::optional<long long> ReadStat(const std::string &file, const std::string &key)
{
  std::ifstream in(file);
  std::string k;
  long long v;
  while (in >> k >> v)
  {
    if (k == key)
    {
      return v;
    }
  }
  return std::nullopt;
}
#endif

}  // namespace

long GetAvailableMemory()
{
  if (available_memory_for_testing >= 0)
  {
    return available_memory_for_testing;
  }
#if defined(__linux__)
  long long avail = -1;
  {
    std::ifstream meminfo("/proc/meminfo");
    std::string key, unit;
    long long value;
    while (meminfo >> key >> value)
    {
      std::getline(meminfo, unit);
      if (key == "MemAvailable:")
      {
        avail = value * 1024;
      }
    }
  }
  if (avail < 0)
  {
    return -1;
  }

  // The mounts of the cgroup hierarchies (version 2, and version 1 with the memory
  // controller), from the fields of /proc/self/mountinfo: ID, parent ID, device, root,
  // mount point, options, optional fields, "-", file system type, source, super options.
  struct Mount
  {
    std::string root, point;
    bool v2;
  };
  std::vector<Mount> mounts;
  {
    std::ifstream mountinfo("/proc/self/mountinfo");
    if (!mountinfo.is_open())
    {
      return -1;
    }
    std::string line;
    while (std::getline(mountinfo, line))
    {
      std::istringstream ss(line);
      std::vector<std::string> fields;
      for (std::string f; ss >> f;)
      {
        fields.push_back(f);
      }
      const auto sep = std::find(fields.begin(), fields.end(), "-");
      if (fields.size() < 5 || sep == fields.end() || fields.end() - sep < 4)
      {
        continue;
      }
      const auto &type = *(sep + 1), &options = *(sep + 3);
      if (type == "cgroup2")
      {
        mounts.push_back({fields[3], fields[4], true});
      }
      else if (type == "cgroup" &&
               ("," + options + ",").find(",memory,") != std::string::npos)
      {
        mounts.push_back({fields[3], fields[4], false});
      }
    }
  }

  std::ifstream cgroups("/proc/self/cgroup");
  if (!cgroups.is_open())
  {
    return -1;
  }
  std::string line;
  while (std::getline(cgroups, line))
  {
    // hierarchy-ID:controller-list:path, with an empty controller list for version 2.
    const auto c1 = line.find(':'), c2 = line.find(':', c1 + 1);
    if (c1 == std::string::npos || c2 == std::string::npos)
    {
      continue;
    }
    const std::string controllers = line.substr(c1 + 1, c2 - c1 - 1),
                      path = line.substr(c2 + 1);
    const bool v2 = controllers.empty();
    if (!v2 && ("," + controllers + ",").find(",memory,") == std::string::npos)
    {
      continue;
    }

    // The mount of the hierarchy whose root contains the group (the longest root).
    const Mount *mount = nullptr;
    for (const auto &m : mounts)
    {
      const bool contains =
          (m.root == "/") || (path.compare(0, m.root.size(), m.root) == 0 &&
                              (path.size() == m.root.size() || path[m.root.size()] == '/'));
      if (m.v2 == v2 && contains && (!mount || m.root.size() > mount->root.size()))
      {
        mount = &m;
      }
    }
    if (!mount)
    {
      // A version 2 hierarchy which is not mounted has no controllers to read (for example
      // the legacy layout of systemd); the version 1 memory controller must be readable.
      if (v2)
      {
        continue;
      }
      return -1;
    }
    if (v2)
    {
      // Without the memory controller (for example the unified hierarchy of the hybrid
      // layout of systemd), the version 2 hierarchy has no memory limits.
      std::ifstream in(mount->point + "/cgroup.controllers");
      std::string c;
      bool memory = false;
      while (in >> c)
      {
        memory = memory || (c == "memory");
      }
      if (!memory)
      {
        continue;
      }
    }

    // The group and its ancestors up to the root of the mount.
    std::string rel = (mount->root == "/") ? path : path.substr(mount->root.size());
    while (true)
    {
      const std::string dir = mount->point + ((rel == "/") ? std::string() : rel) + "/";
      const bool top = rel.empty() || rel == "/";
      std::optional<long long> limit, usage, inactive;
      if (v2)
      {
        const auto max = ReadWord(dir + "memory.max");
        if (!max)
        {
          // The root group has no limit file, and the groups below it have one where the
          // memory controller is enabled.
          std::ifstream in(dir + "cgroup.controllers");
          const bool readable = in.is_open();
          std::string c;
          bool memory = false;
          while (in >> c)
          {
            memory = memory || (c == "memory");
          }
          if (!top && (!readable || memory))
          {
            return -1;
          }
        }
        else if (*max != "max")
        {
          try
          {
            limit = std::stoll(*max);
          }
          catch (const std::exception &)
          {
            return -1;
          }
        }
        if (limit)
        {
          usage = ReadNumber(dir + "memory.current");
          inactive = ReadStat(dir + "memory.stat", "inactive_file");
        }
      }
      else
      {
        // Version 1 reports no limit as a value near the largest integer.
        limit = ReadNumber(dir + "memory.limit_in_bytes");
        usage = ReadNumber(dir + "memory.usage_in_bytes");
        inactive = ReadStat(dir + "memory.stat", "total_inactive_file");
        if (!inactive)
        {
          inactive = ReadStat(dir + "memory.stat", "inactive_file");
        }
        if (!limit)
        {
          return -1;
        }
        if (*limit >= (1LL << 60))
        {
          limit.reset();
        }
      }
      if (limit)
      {
        if (!usage || !inactive)
        {
          return -1;
        }
        avail = std::min(avail, std::max(*limit - (*usage - *inactive), 0LL));
      }
      if (top)
      {
        break;
      }
      const auto slash = rel.find_last_of('/');
      rel = (slash == 0 || slash == std::string::npos) ? "/" : rel.substr(0, slash);
    }
  }
  return static_cast<long>(avail);
#else
  return -1;
#endif
}

void SetAvailableMemoryForTesting(long bytes)
{
  available_memory_for_testing = bytes;
}

bool ResetPeakMemory()
{
#if defined(__linux__)
  std::ofstream clear_refs("/proc/self/clear_refs");
  if (!clear_refs.is_open())
  {
    return false;
  }
  clear_refs << "5";
  clear_refs.close();
  return !clear_refs.fail();
#else
  return false;
#endif
}

const NodeComms &GetNodeComms(MPI_Comm comm)
{
  static MPI_Comm cached_comm = MPI_COMM_NULL;
  static NodeComms comms;
  if (cached_comm != comm || comms.node == MPI_COMM_NULL)
  {
    // Communicators of a previous comm are leaked deliberately: freeing is collective
    // over the old communicator.
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comms.node);
    comms.node_rank = Mpi::Rank(comms.node);
    MPI_Comm_split(comm, (comms.node_rank == 0) ? 0 : MPI_UNDEFINED, 0, &comms.leaders);
    cached_comm = comm;
  }
  return comms;
}

std::string FormatBytes(double bytes)
{
  constexpr double kilobyte = 1024.0;
  constexpr double megabyte = kilobyte * 1024.0;
  constexpr double gigabyte = megabyte * 1024.0;
  constexpr double terabyte = gigabyte * 1024.0;
  if (bytes >= terabyte)
  {
    return fmt::format("{:.1f}T", bytes / terabyte);
  }
  if (bytes >= gigabyte)
  {
    return fmt::format("{:.1f}G", bytes / gigabyte);
  }
  if (bytes >= megabyte)
  {
    return fmt::format("{:.1f}M", bytes / megabyte);
  }
  return fmt::format("{:.1f}K", bytes / kilobyte);
}

namespace
{

MemoryStats ComputeStats(std::string label, long local_value, MPI_Comm comm)
{
  long val_min = local_value, val_max = local_value, val_sum = local_value;
  Mpi::GlobalMin(1, &val_min, comm);
  Mpi::GlobalMax(1, &val_max, comm);
  Mpi::GlobalSum(1, &val_sum, comm);
  return {label, val_min, val_max, val_sum, static_cast<double>(val_sum) / Mpi::Size(comm)};
}

}  // namespace

MemoryStats ComputeNodeMemoryStats(std::string label, long local_value, MPI_Comm comm)
{
  const auto &[node_comm, leaders_comm, node_rank] = GetNodeComms(comm);

  // Sum memory usage across all processes on this node. Assumes shared memory is
  // negligible.
  long node_total_memory = local_value;
  MPI_Reduce(&local_value, &node_total_memory, 1, MPI_LONG, MPI_SUM, 0, node_comm);

  // Compute statistics across node leaders.
  MemoryStats stats;
  stats.label = label;
  if (node_rank == 0 && leaders_comm != MPI_COMM_NULL)
  {
    stats.min = node_total_memory;
    stats.max = node_total_memory;
    stats.sum = node_total_memory;
    Mpi::GlobalMin(1, &stats.min, leaders_comm);
    Mpi::GlobalMax(1, &stats.max, leaders_comm);
    Mpi::GlobalSum(1, &stats.sum, leaders_comm);
    stats.avg = static_cast<double>(stats.sum) / Mpi::Size(leaders_comm);
  }

  // Broadcast results from node leaders to all ranks.
  Mpi::Broadcast(1, &stats.min, 0, node_comm);
  Mpi::Broadcast(1, &stats.max, 0, node_comm);
  Mpi::Broadcast(1, &stats.sum, 0, node_comm);
  Mpi::Broadcast(1, &stats.avg, 0, node_comm);
  return stats;
}

MemoryStats GetCurrentMemoryStats(MPI_Comm comm)
{
  return ComputeStats("current per-rank", GetCurrentMemory(), comm);
}

MemoryStats GetPeakMemoryStats(MPI_Comm comm)
{
  return ComputeStats("peak per-rank", GetPeakMemory(), comm);
}

MemoryStats GetCurrentNodeMemoryStats(MPI_Comm comm)
{
  return ComputeNodeMemoryStats("current per-node", GetCurrentMemory(), comm);
}

MemoryStats GetPeakNodeMemoryStats(MPI_Comm comm)
{
  return ComputeNodeMemoryStats("peak per-node", GetPeakMemory(), comm);
}

void PrintMemoryUsage(MPI_Comm comm, const MemoryStats &stats)
{
  // For per-node stats, we sum all processes on each node, assuming shared
  // memory is negligible.
  Mpi::Print(comm, "Estimated {} memory usage is: Min. {}, Max. {}, Avg. {}, Total {}\n",
             stats.label, FormatBytes(stats.min), FormatBytes(stats.max),
             FormatBytes(stats.avg), FormatBytes(stats.sum));
}

}  // namespace palace::memory_reporting
