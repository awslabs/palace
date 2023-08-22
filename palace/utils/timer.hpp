// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_TIMER_HPP
#define PALACE_UTILS_TIMER_HPP

#include <chrono>
#include <stack>
#include <string>
#include <vector>
#include "drivers/basesolver.hpp"
#include "utils/communication.hpp"

namespace palace
{

//
// Timer classes for profiling.
//

class Timer
{
public:
  using Clock = std::chrono::steady_clock;
  using Duration = std::chrono::duration<double>;
  using TimePoint = typename Clock::time_point;

  enum Index
  {
    INIT = 0,
    CONSTRUCT,
    FREQSAMPLESELECTION,
    HDMSOLVE,
    SOLVE,
    POSTPRO,
    IO,
    TOTAL,
    NUMTIMINGS
  };

  inline static const std::vector<std::string> descriptions{"Initialization",
                                                            "Operator Construction",
                                                            "Frequency Sample Selection",
                                                            "HDM Solve",
                                                            "Solve",
                                                            "Postprocessing",
                                                            "Disk IO",
                                                            "Total"};

private:
  const TimePoint start_time;
  TimePoint last_lap_time;
  std::vector<Duration> data;
  std::vector<int> counts;
  std::vector<double> data_min, data_max, data_avg;

  // Save a timing step by adding a duration, without lapping; optionally, count it.
  Duration SaveTime(Index idx, Duration time, bool count_it)
  {
    data[idx] += time;
    count_it &&counts[idx]++;
    return data[idx];
  }

public:
  Timer()
    : start_time(Now()), last_lap_time(start_time), data(NUMTIMINGS), counts(NUMTIMINGS)
  {
  }

  // Get the current time.
  static TimePoint Now() { return Clock::now(); }

  // Provide stopwatch lap split functionality.
  Duration Lap()
  {
    auto temp_time = last_lap_time;
    last_lap_time = Now();
    return last_lap_time - temp_time;
  }

  // Return the time elapsed since timer creation.
  Duration TimeFromStart() const { return Now() - start_time; }

  // Lap and record a timing step.
  Duration MarkTime(Index idx, bool count_it = true)
  {
    return SaveTime(idx, Lap(), count_it);
  }

  // Provide map-like read-only access to the timing data.
  Duration operator[](Index idx) const { return data[idx]; }

  // Provide access to the reduced timing data.
  double GetMinTime(Index idx) const { return data_min[idx]; }
  double GetMaxTime(Index idx) const { return data_max[idx]; }
  double GetAvgTime(Index idx) const { return data_avg[idx]; }

  // Return number of times timer.MarkTime(idx) or TimerBlock b(idx) was called.
  int GetCounts(Index idx) const { return counts[idx]; }

  // Reduce timing information across MPI ranks.
  void Reduce(MPI_Comm comm)
  {
    const std::size_t ntimes = data.size();
    SaveTime(TOTAL, TimeFromStart(), true);
    data_min.resize(ntimes);
    data_max.resize(ntimes);
    data_avg.resize(ntimes);

    int np = Mpi::Size(comm);
    for (std::size_t i = 0; i < ntimes; i++)
    {
      data_min[i] = data_max[i] = data_avg[i] = data[i].count();
    }

    Mpi::GlobalMin(ntimes, data_min.data(), comm);
    Mpi::GlobalMax(ntimes, data_max.data(), comm);
    Mpi::GlobalSum(ntimes, data_avg.data(), comm);

    for (std::size_t i = 0; i < ntimes; i++)
    {
      data_avg[i] /= np;
    }
  }

  // Prints timing information. We assume the data has already been reduced.
  void Print(MPI_Comm comm) const
  {
    constexpr int p = 3;   // Floating point precision
    constexpr int w = 12;  // Data column width
    constexpr int h = 26;  // Left-hand side width
    // clang-format off
    Mpi::Print(comm, "\n{:<{}s}{:>{}s}{:>{}s}{:>{}s}\n",
               "Elapsed Time Report (s)", h, "Min.", w, "Max.", w, "Avg.", w);
    // clang-format on
    Mpi::Print(comm, "{}\n", std::string(h + 3 * w, '='));
    for (int i = INIT; i < NUMTIMINGS; i++)
    {
      if (counts[i] > 0)
      {
        if (i == TOTAL)
        {
          Mpi::Print(comm, "{}\n", std::string(h + 3 * w, '-'));
        }
        // clang-format off
        Mpi::Print(comm, "{:<{}s}{:{}.{}f}{:{}.{}f}{:{}.{}f}\n",
                   descriptions[i], h,
                   data_min[i], w, p, data_max[i], w, p, data_avg[i], w, p);
        // clang-format on
      }
    }
  }
};

class BlockTimer
{
  using Index = Timer::Index;

private:
  inline static std::stack<Index> stack;
  inline static Timer timer;

public:
  BlockTimer(Index i)
  {
    // Start timing when entering the block, interrupting whatever we were timing before.
    // Take note of what we are now timing.
    (stack.empty()) ? timer.Lap() : timer.MarkTime(stack.top(), false);
    stack.push(i);
  }

  ~BlockTimer()
  {
    // When a BlockTimer is no longer in scope, record the time.
    // (check whether stack is empty in case Finalize was called already.)
    if (!stack.empty())
    {
      timer.MarkTime(stack.top());
      stack.pop();
    }
  }

  static void Finalize(MPI_Comm comm, BaseSolver &solver)
  {
    while (!stack.empty())
    {
      timer.MarkTime(stack.top());
      stack.pop();
    }
    timer.Reduce(comm);
    timer.Print(comm);
    solver.SaveMetadata(timer);
    Mpi::Print(comm, "\n");
  }
};

}  // namespace palace

#endif  // PALACE_UTILS_TIMER_HPP
