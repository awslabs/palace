// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <mpi.h>
#include "utils/communication.hpp"
#include "utils/perfstats.hpp"

namespace palace::perf
{

namespace detail
{

bool mpi_wait_seen = false, mpi_wait_tool = false;

namespace
{

struct CommTallies
{
  int reduce, barrier, broadcast, wait, reduce_bytes, broadcast_bytes;
  std::array<int, Timer::NUM_TIMINGS> wait_rows;  // Wait/MPI/<ROW>
};

// Communication tallies, registered together and with all their rows at once: which rows
// communicate can depend on the rank (root-only waits, collectives on sub-communicators).
const CommTallies &Tallies()
{
  static const CommTallies t = {
      RegisterPhaseTally("Reduce", true),
      RegisterPhaseTally("Barrier", true),
      RegisterPhaseTally("Broadcast", true),
      RegisterPhaseTally("Wait/MPI", true),
      RegisterCounter("Communication/ReduceBytes", Reduce::STATS, Accumulate::SUM),
      RegisterCounter("Communication/BroadcastBytes", Reduce::STATS, Accumulate::SUM),
      []()
      {
        std::array<int, Timer::NUM_TIMINGS> rows;
        const int wait = RegisterPhaseTally("Wait/MPI", true);
        for (int i = 0; i < Timer::NUM_TIMINGS; i++)
        {
          rows[i] = PhaseTallyRow(wait, static_cast<Timer::Index>(i));
        }
        return rows;
      }()};
  return t;
}

}  // namespace

void AddCollective(Collective kind, double seconds, std::size_t bytes)
{
  comm_wait += seconds;
  if (registry.suspended)
  {
    return;
  }
  const auto &t = Tallies();
  switch (kind)
  {
    case Collective::REDUCE:
      AddTally(t.reduce, seconds, 1, seconds);
      Add(t.reduce_bytes, static_cast<double>(bytes));
      break;
    case Collective::BARRIER:
      AddTally(t.barrier, seconds, 1, seconds);
      break;
    case Collective::BROADCAST:
      AddTally(t.broadcast, seconds, 1, seconds);
      Add(t.broadcast_bytes, static_cast<double>(bytes));
      break;
  }
}

void AddWait(double seconds)
{
  // Called for every blocking MPI_Wait family call: the phase tally is updated directly.
  comm_wait += seconds;
  static thread_local const bool main_thread = []()
  {
    int flag = 0;
    PMPI_Is_thread_main(&flag);
    return flag != 0;
  }();
  if (!main_thread || registry.suspended || registry.paused)
  {
    return;
  }
  const auto &t = Tallies();
  for (int id : {t.wait, t.wait_rows[BlockTimer::Current()]})
  {
    auto &tally = registry.tallies[id];
    tally.time += seconds;
    tally.wait += seconds;
    tally.count++;
  }
}

bool MpiWaitTally()
{
  static const bool intercepted = []()
  {
    mpi_wait_seen = false;
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    return mpi_wait_seen;
  }();
  return intercepted;
}

}  // namespace detail

void RegisterCommunicationTallies()
{
  detail::Tallies();
  detail::MpiWaitTally();
}

}  // namespace palace::perf
