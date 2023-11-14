// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_TIMER_HPP
#define PALACE_UTILS_TIMER_HPP

#include <chrono>
#include <stack>
#include <string>
#include <vector>
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
    WAVEPORT,  // Wave port solver
    SOLVE,
    PRECONDITIONER,      // Linear solver
    COARSESOLVE,         // Linear solver
    ESTIMATION,          // Estimation
    CONSTRUCTESTIMATOR,  // Construction of estimator
    SOLVEESTIMATOR,      // Evaluation of estimator
    ADAPTATION,          // Adaptation
    REBALANCE,           // Rebalancing
    CONSTRUCTPROM,       // Adaptive frequency sweep
    SOLVEPROM,           // Adaptive frequency sweep
    POSTPRO,
    IO,
    TOTAL,
    NUMTIMINGS
  };

  // clang-format off
  inline static const std::vector<std::string> descriptions{
      "Initialization",
      "Operator Construction",
      "  Wave Ports",
      "Solve",
      "  Preconditioner",
      "  Coarse Solve",
      "Estimation",
      "  Construction",
      "  Solve",
      "Adaptation",
      "  Rebalancing",
      "PROM Construction",
      "PROM Solve",
      "Postprocessing",
      "Disk IO",
      "Total"};
  // clang-format on

private:
  const TimePoint start_time;
  TimePoint last_lap_time;
  std::vector<Duration> data;
  std::vector<int> counts;

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

  // Save a timing step by adding a duration, without lapping; optionally, count it.
  Duration SaveTime(Index idx, Duration time, bool count_it = true)
  {
    data[idx] += time;
    counts[idx] += count_it;
    return data[idx];
  }

  // Lap and record a timing step.
  Duration MarkTime(Index idx, bool count_it = true)
  {
    return SaveTime(idx, Lap(), count_it);
  }

  // Provide map-like read-only access to the timing data.
  auto Data(Index idx) const { return data[idx].count(); }

  // Return number of times timer.MarkTime(idx) or TimerBlock b(idx) was called.
  auto Counts(Index idx) const { return counts[idx]; }
};

class BlockTimer
{
  using Index = Timer::Index;

private:
  inline static Timer timer;
  inline static std::stack<Index> stack;

  // Reduce timing information across MPI ranks.
  static void Reduce(MPI_Comm comm, std::vector<double> &data_min,
                     std::vector<double> &data_max, std::vector<double> &data_avg)
  {
    data_min.resize(Timer::NUMTIMINGS);
    data_max.resize(Timer::NUMTIMINGS);
    data_avg.resize(Timer::NUMTIMINGS);
    for (int i = Timer::INIT; i < Timer::NUMTIMINGS; i++)
    {
      data_min[i] = data_max[i] = data_avg[i] = timer.Data((Timer::Index)i);
    }

    Mpi::GlobalMin(Timer::NUMTIMINGS, data_min.data(), comm);
    Mpi::GlobalMax(Timer::NUMTIMINGS, data_max.data(), comm);
    Mpi::GlobalSum(Timer::NUMTIMINGS, data_avg.data(), comm);

    const int np = Mpi::Size(comm);
    for (int i = Timer::INIT; i < Timer::NUMTIMINGS; i++)
    {
      data_avg[i] /= np;
    }
  }

public:
  BlockTimer(Index i)
  {
    // Start timing when entering the block, interrupting whatever we were timing before.
    // Take note of what we are now timing.
    stack.empty() ? timer.Lap() : timer.MarkTime(stack.top(), false);
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

  // Read-only access the static Timer object.
  static const Timer &GlobalTimer() { return timer; }

  // Print timing information after reducing the data across all processes.
  static void Print(MPI_Comm comm)
  {
    while (!stack.empty())
    {
      timer.MarkTime(stack.top());
      stack.pop();
    }
    timer.SaveTime(Timer::TOTAL, timer.TimeFromStart());

    // Reduce timing data.
    std::vector<double> data_min, data_max, data_avg;
    Reduce(comm, data_min, data_max, data_avg);

    // Print a nice table of the timing data.
    constexpr int p = 3;   // Floating point precision
    constexpr int w = 12;  // Data column width
    constexpr int h = 26;  // Left-hand side width
    // clang-format off
    Mpi::Print(comm, "\n{:<{}s}{:>{}s}{:>{}s}{:>{}s}\n",
               "Elapsed Time Report (s)", h, "Min.", w, "Max.", w, "Avg.", w);
    // clang-format on
    Mpi::Print(comm, "{}\n", std::string(h + 3 * w, '='));
    for (int i = Timer::INIT; i < Timer::NUMTIMINGS; i++)
    {
      if (timer.Counts((Timer::Index)i) > 0)
      {
        if (i == Timer::TOTAL)
        {
          Mpi::Print(comm, "{}\n", std::string(h + 3 * w, '-'));
        }
        // clang-format off
        Mpi::Print(comm, "{:<{}s}{:{}.{}f}{:{}.{}f}{:{}.{}f}\n",
                   timer.descriptions[i], h,
                   data_min[i], w, p, data_max[i], w, p, data_avg[i], w, p);
        // clang-format on
      }
    }
  }
};

}  // namespace palace

#endif  // PALACE_UTILS_TIMER_HPP
