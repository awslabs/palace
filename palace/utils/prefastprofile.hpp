// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
//
// LANE TOOL (v1prefast PART 1), NEVER MERGED: the profile image's phase marks. Every call prints ONE line to stdout from
// the calling rank: "PREFAST <rank> <epoch s> <mark> rss_mb=<VmRSS> hwm_mb=<VmHWM>" (no collective; the root's marks
// inside root-only code print from the root alone). Read by tools/prefast_marks_v1.py.
#ifndef PALACE_UTILS_PREFAST_PROFILE_HPP
#define PALACE_UTILS_PREFAST_PROFILE_HPP

#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#include <thread>
#include <mpi.h>

namespace palace::prefast
{

inline void Mark(const std::string &mark)
{
  int initialized = 0, finalized = 0, rank = -1;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized && !finalized)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
  long rss_kb = -1, hwm_kb = -1;
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line))
  {
    if (line.rfind("VmRSS:", 0) == 0)
    {
      rss_kb = std::stol(line.substr(6));
    }
    else if (line.rfind("VmHWM:", 0) == 0)
    {
      hwm_kb = std::stol(line.substr(6));
    }
  }
  const double now =
      std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch())
          .count();
  std::printf("PREFAST %d %.3f %s rss_mb=%ld hwm_mb=%ld\n", rank, now, mark.c_str(),
              rss_kb / 1024, hwm_kb / 1024);
  std::fflush(stdout);
}

// The calling process's allowed CPUs (/proc/self/status Cpus_allowed_list) and the
// hardware concurrency: whether a thread of this rank could run beside it.
inline std::string AllowedCpus()
{
  std::ifstream status("/proc/self/status");
  std::string line, list = "?";
  while (std::getline(status, line))
  {
    if (line.rfind("Cpus_allowed_list:", 0) == 0)
    {
      list = line.substr(18);
      list.erase(0, list.find_first_not_of(" \t"));
    }
  }
  return list;
}

}  // namespace palace::prefast

#endif  // PALACE_UTILS_PREFAST_PROFILE_HPP
