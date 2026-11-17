// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_TIMER_HPP
#define PALACE_UTILS_TIMER_HPP

#include <array>
#include <chrono>
#include <string_view>
#include <vector>
#include "utils/memoryreporting.hpp"
#include "utils/perfmetrics.hpp"

namespace palace
{

//
// Exclusive timer rows of time and peak memory growth, per iteration and over the run.
//

class Timer
{
public:
  using Clock = std::chrono::steady_clock;
  using Duration = std::chrono::duration<double>;
  using TimePoint = typename Clock::time_point;

  enum Index
  {
    STARTUP = 0,           // Process start to the first timer
    INIT,                  // Initialization
    MESH_PREPROCESS,       // Preprocessing mesh
    CONSTRUCT,             // Space and operator construction
    WAVE_PORT,             // Wave port solver
    WAVE_PORT_SOLVE,       // Reduced/exact port eigenproblem
    WAVE_PORT_FIELD,       // Port eigenvector reconstruction and distribution
    WAVE_PORT_POSTPRO,     // Port mode normalization and projection forms
    KSP,                   // Linear solver
    KSP_SETUP,             // Linear solver setup
    KSP_PRECONDITIONER,    // Linear solver preconditioner
    KSP_COARSE_SETUP,      // Linear solver coarse-level setup and factorization
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
  static constexpr auto descriptions = std::to_array<std::string_view>({
      "Startup",
      "Initialization",
      "  Mesh Preprocessing",
      "Operator Construction",
      "  Wave Ports",
      "    Eigenproblem",
      "    Field reconstruction",
      "    Mode postprocessing",
      "Linear Solve",
      "  Setup",
      "  Preconditioner",
      "  Coarse Setup",
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
      "Total"});

  // Keys in palace.json.
  static constexpr auto names = std::to_array<std::string_view>({
      "STARTUP", "INIT", "MESH_PREPROCESS", "CONSTRUCT", "WAVE_PORT", "WAVE_PORT_SOLVE",
      "WAVE_PORT_FIELD", "WAVE_PORT_POSTPRO", "KSP", "KSP_SETUP", "KSP_PRECONDITIONER",
      "KSP_COARSE_SETUP", "KSP_COARSE_SOLVE", "TS", "EPS", "DIV_FREE", "CONSTRUCT_PROM",
      "SOLVE_PROM", "ESTIMATION", "CONSTRUCT_ESTIMATOR", "SOLVE_ESTIMATOR", "ADAPTATION",
      "REBALANCE", "POSTPRO", "POSTPRO_FARFIELD", "POSTPRO_PARAVIEW",
      "POSTPRO_GRIDFUNCTION", "IO", "TOTAL"});
  // clang-format on
  static_assert(descriptions.size() == NUM_TIMINGS);
  static_assert(names.size() == NUM_TIMINGS);

private:
  TimePoint start_time, iteration_start_time, last_lap_time, report_time;
  std::vector<Duration> data, run_data;  // This iteration; previous iterations
  std::vector<long> counts, run_counts;
  std::vector<long> mem_data, run_mem_data;
  long last_memory;
  long iteration_baseline = 0;  // Process peak at the iteration start (0: first)
  int iterations = 1;

public:
  Timer()
    : start_time(Now()), iteration_start_time(start_time), last_lap_time(start_time),
      report_time(start_time), data(NUM_TIMINGS), run_data(NUM_TIMINGS),
      counts(NUM_TIMINGS, 0), run_counts(NUM_TIMINGS, 0), mem_data(NUM_TIMINGS, 0),
      run_mem_data(NUM_TIMINGS, 0), last_memory(memory_reporting::GetPeakMemory())
  {
  }

  static TimePoint Now() { return Clock::now(); }

  Duration Lap()
  {
    auto temp_time = last_lap_time;
    last_lap_time = Now();
    return last_lap_time - temp_time;
  }

  // Lap and record a timing step of a row (not TOTAL); optionally, count it.
  void MarkTime(Index idx, bool count_it = true)
  {
    data[idx] += Lap();
    counts[idx] += count_it;
  }

  long MemoryLap()
  {
    const long current = memory_reporting::GetPeakMemory();
    const long delta = current - last_memory;
    last_memory = current;
    return delta;
  }

  // Lap and record the peak memory growth of a row (not TOTAL).
  void MarkMemory(Index idx) { mem_data[idx] += MemoryLap(); }

  // Fix the time of the report: TOTAL ends there, and the next iteration starts there.
  void MarkReport() { report_time = Now(); }

  // This iteration (TOTAL: elapsed up to the report).
  double Data(Index idx) const
  {
    return (idx == TOTAL) ? Duration(report_time - iteration_start_time).count()
                          : data[idx].count();
  }
  long Counts(Index idx) const { return (idx == TOTAL) ? 1 : counts[idx]; }
  long MemoryData(Index idx) const { return mem_data[idx]; }
  long IterationBaseline() const { return iteration_baseline; }
  int Iterations() const { return iterations; }

  // Run so far, including this iteration.
  double RunData(Index idx) const
  {
    return (idx == TOTAL) ? Duration(report_time - start_time).count()
                          : (run_data[idx] + data[idx]).count();
  }
  long RunCounts(Index idx) const
  {
    return (idx == TOTAL) ? iterations : run_counts[idx] + counts[idx];
  }
  long RunMemoryData(Index idx) const { return run_mem_data[idx] + mem_data[idx]; }

  // Fold this iteration into the run and start the next one at the report time.
  void StartIteration()
  {
    for (int i = 0; i < NUM_TIMINGS; i++)
    {
      run_data[i] += data[i];
      run_counts[i] += counts[i];
      run_mem_data[i] += mem_data[i];
      data[i] = Duration::zero();
      counts[i] = 0;
      mem_data[i] = 0;
    }
    iteration_start_time = report_time;
    iteration_baseline = memory_reporting::GetPeakMemory();
    iterations++;
  }
};

class BlockTimer
{
  using Index = Timer::Index;

private:
  inline static Timer timer;
  inline static std::vector<Index> stack;
  inline static bool started = false;
  bool count;

public:
  BlockTimer(Index i, bool count_it = true) : count(count_it)
  {
    // Start timing when entering the block, interrupting whatever we were timing before.
    if (count)
    {
      if (!started)
      {
        // Everything before the first timer is startup.
        started = true;
        timer.MarkTime(Timer::STARTUP);
        timer.MarkMemory(Timer::STARTUP);
      }
      else if (stack.empty())
      {
        timer.Lap();
        timer.MemoryLap();
      }
      else
      {
        timer.MarkTime(stack.back(), false);
        timer.MarkMemory(stack.back());
      }
      stack.push_back(i);
    }
  }

  ~BlockTimer()
  {
    if (count && !stack.empty())
    {
      timer.MarkTime(stack.back());
      timer.MarkMemory(stack.back());
      stack.pop_back();
    }
  }

  static const Timer &GlobalTimer() { return timer; }

  // Before a report: account the open innermost row up to now, without closing it, and fix
  // the report time.
  static void Snapshot()
  {
    if (!stack.empty())
    {
      timer.MarkTime(stack.back(), false);
      timer.MarkMemory(stack.back());
    }
    else if (started)
    {
      timer.Lap();
      timer.MemoryLap();
    }
    timer.MarkReport();
  }

  // After a report: start the next iteration (rows and metrics).
  static void StartIteration()
  {
    timer.StartIteration();
    perf::detail::FoldIteration();
  }

  // Clear all process-wide timing state (successive runs in one process, as in the
  // regression test harness). Metric registrations are kept.
  static void Reset()
  {
    timer = Timer();
    stack.clear();
    started = false;
    perf::detail::ResetValues();
  }
};

}  // namespace palace

#endif  // PALACE_UTILS_TIMER_HPP
