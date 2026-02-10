// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "memoryreporting.hpp"

#include <sys/resource.h>

#if defined(__APPLE__)
#include <mach/mach.h>
#elif defined(__linux__)
#include <fstream>
#include <string>
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
  if (task_info(mach_task_self(), TASK_VM_INFO, (task_info_t)&vm_info, &count) ==
      KERN_SUCCESS)
  {
    return static_cast<long>(vm_info.phys_footprint);
  }
  return 0;
#elif defined(__linux__)
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line))
  {
    if (line.compare(0, 6, "VmRSS:") == 0)
    {
      // Value is in kB.
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
    // Linux reports ru_maxrss in kB.
    return usage.ru_maxrss * 1024;
#endif
  }
  return 0;
}

namespace
{

std::string FormatBytes(double bytes)
{
  constexpr double kB = 1024.0;
  constexpr double MB = kB * 1024.0;
  constexpr double GB = MB * 1024.0;
  constexpr double TB = GB * 1024.0;
  if (bytes >= TB)
  {
    return fmt::format("{:.1f}T", bytes / TB);
  }
  if (bytes >= GB)
  {
    return fmt::format("{:.1f}G", bytes / GB);
  }
  if (bytes >= MB)
  {
    return fmt::format("{:.1f}M", bytes / MB);
  }
  return fmt::format("{:.1f}K", bytes / kB);
}

MemoryStats ComputeStats(std::string label, long local_value, MPI_Comm comm)
{
  long val_min = local_value, val_max = local_value, val_sum = local_value;
  Mpi::GlobalMin(1, &val_min, comm);
  Mpi::GlobalMax(1, &val_max, comm);
  Mpi::GlobalSum(1, &val_sum, comm);
  return {label, val_min, val_max, val_sum, static_cast<double>(val_sum) / Mpi::Size(comm)};
}

}  // namespace

MemoryStats GetCurrentMemoryStats(MPI_Comm comm)
{
  return ComputeStats("Current", GetCurrentMemory(), comm);
}

MemoryStats GetPeakMemoryStats(MPI_Comm comm)
{
  return ComputeStats("Peak", GetPeakMemory(), comm);
}

void PrintMemoryUsage(MPI_Comm comm, const MemoryStats &stats)
{
  Mpi::Print(comm, "\n{} RSS: Min. {}, Max. {}, Avg. {}, Total {}\n", stats.label,
             FormatBytes(stats.min), FormatBytes(stats.max), FormatBytes(stats.avg),
             FormatBytes(stats.sum));
}

}  // namespace palace::memory_reporting
