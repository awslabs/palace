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
    MESH_PREPROCESS,       // Preprocessing mesh
    CONSTRUCT,             // Space and operator construction
    WAVE_PORT,             // Wave port solver
    KSP,                   // Linear solver
    KSP_SETUP,             // Linear solver setup
    KSP_PRECONDITIONER,    // Linear solver preconditioner
    KSP_COARSE_SOLVE,      // Linear solver coarse-level solve
    TS,                    // Time integrator
    EPS,                   // Eigenvalue problem solver
    DIV_FREE,              // Divergence-free projection
    CONSTRUCT_PROM,        // Adaptive frequency sweep offline
    SOLVE_PROM,            // Adaptive frequency sweep online
    ESTIMATION,            // Error estimation
    CONSTRUCT_ESTIMATOR,   // Construction of estimator
    SOLVE_ESTIMATOR,       // Evaluation of estimator
    ADAPTATION,            // Adaptation
    REBALANCE,             // Rebalancing
    POSTPRO,               // Solution postprocessing
    POSTPRO_FARFIELD,      // Computing far-fields
    POSTPRO_PARAVIEW,      // Paraview calculations and I/O
    POSTPRO_GRIDFUNCTION,  // MFEM gridfunction calculations and I/O
    IO,                    // Disk I/O
    TOTAL,
    NUM_TIMINGS
  };

  // clang-format off
  inline static const std::vector<std::string> descriptions{
      "Initialization",
      "  Mesh Preprocessing",
      "Operator Construction",
      "  Wave Ports",
      "Linear Solve",
      "  Setup",
      "  Preconditioner",
      "  Coarse Solve",
      "Time Stepping",
      "Eigenvalue Solve",
      "Div.-Free Projection",
      "PROM Construction",
      "PROM Solve",
      "Estimation",
      "  Construction",
      "  Solve",
      "Adaptation",
      "  Rebalancing",
      "Postprocessing",
      "  Far Fields",
      "  Paraview",
      "  Grid function",
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
    : start_time(Now()), last_lap_time(start_time), data(NUM_TIMINGS), counts(NUM_TIMINGS)
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
    return MarkTime(idx, Lap(), count_it);
  }

  // Record a timing step by adding a duration, without lapping; optionally, count it.
  Duration MarkTime(Index idx, Duration time, bool count_it = true)
  {
    if (idx == Timer::TOTAL)
    {
      data[idx] = time;
    }
    else
    {
      data[idx] += time;
    }
    counts[idx] += count_it;
    return data[idx];
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
  bool count;

  // Reduce timing information across MPI ranks.
  static void Reduce(MPI_Comm comm, std::vector<double> &data_min,
                     std::vector<double> &data_max, std::vector<double> &data_avg)
  {
    data_min.resize(Timer::NUM_TIMINGS);
    data_max.resize(Timer::NUM_TIMINGS);
    data_avg.resize(Timer::NUM_TIMINGS);
    for (int i = Timer::INIT; i < Timer::NUM_TIMINGS; i++)
    {
      data_min[i] = data_max[i] = data_avg[i] = timer.Data((Timer::Index)i);
    }

    Mpi::GlobalMin(Timer::NUM_TIMINGS, data_min.data(), comm);
    Mpi::GlobalMax(Timer::NUM_TIMINGS, data_max.data(), comm);
    Mpi::GlobalSum(Timer::NUM_TIMINGS, data_avg.data(), comm);

    const int np = Mpi::Size(comm);
    for (int i = Timer::INIT; i < Timer::NUM_TIMINGS; i++)
    {
      data_avg[i] /= np;
    }
  }

public:
  BlockTimer(Index i, bool count = true) : count(count)
  {
    // Start timing when entering the block, interrupting whatever we were timing before.
    // Take note of what we are now timing.
    if (count)
    {
      stack.empty() ? timer.Lap() : timer.MarkTime(stack.top(), false);
      stack.push(i);
    }
  }

  ~BlockTimer()
  {
    // When a BlockTimer is no longer in scope, record the time (check whether stack is
    // empty in case the timer has already been finalized).
    if (count && !stack.empty())
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
    timer.MarkTime(Timer::TOTAL, timer.TimeFromStart());

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
    for (int i = Timer::INIT; i < Timer::NUM_TIMINGS; i++)
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
