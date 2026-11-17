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

  // Time, count and peak memory growth of a row.
  struct Row
  {
    Duration time{};
    long count = 0, mem = 0;
    Row &operator+=(const Row &r)
    {
      time += r.time;
      count += r.count;
      mem += r.mem;
      return *this;
    }
  };

private:
  std::array<Row, NUM_TIMINGS> iter{}, run{};  // This iteration; previous iterations
  TimePoint start_time, iteration_start_time, last_lap_time, report_time;
  long last_memory;
  long report_memory = 0;       // Process peak at the last report
  long iteration_baseline = 0;  // Process peak at the iteration start (0: first)
  long run_baseline = 0;        // Process peak at the run start (0: process start)
  int iterations = 1;

public:
  Timer()
    : start_time(Now()), iteration_start_time(start_time), last_lap_time(start_time),
      report_time(start_time), last_memory(memory_reporting::GetPeakMemory())
  {
  }

  static TimePoint Now() { return Clock::now(); }

  // Time and peak memory growth since the last lap.
  Row Lap()
  {
    const auto time = Now();
    const long memory = memory_reporting::GetPeakMemory();
    Row r{time - last_lap_time, 0, memory - last_memory};
    last_lap_time = time;
    last_memory = memory;
    return r;
  }

  // Lap and charge the step to a row (not TOTAL); optionally, count it.
  void Mark(Index idx, bool count_it = true)
  {
    auto r = Lap();
    r.count = count_it;
    iter[idx] += r;
  }

  // Fix the time of the report: TOTAL ends there, and the next iteration starts there.
  // Call right after a lap: the report time and memory are those of the last lap.
  void MarkReport()
  {
    report_time = last_lap_time;
    report_memory = last_memory;
  }

  // Peak memory at the start of the run (after a reset).
  void SetStartBaseline(long b) { iteration_baseline = run_baseline = b; }

  // A row of this iteration, or of the run so far including this iteration. TOTAL: elapsed
  // up to the report, counting iterations, no memory.
  Row Get(Index idx, bool cumulative = false) const
  {
    if (idx == TOTAL)
    {
      return {report_time - (cumulative ? start_time : iteration_start_time),
              cumulative ? iterations : 1, 0};
    }
    Row r = iter[idx];
    if (cumulative)
    {
      r += run[idx];
    }
    return r;
  }
  long Baseline(bool cumulative) const
  {
    return cumulative ? run_baseline : iteration_baseline;
  }
  long ReportMemory() const { return report_memory; }
  int Iterations() const { return iterations; }

  // Fold this iteration into the run and start the next one at the report time.
  void StartIteration()
  {
    for (int i = 0; i < NUM_TIMINGS; i++)
    {
      run[i] += iter[i];
    }
    iter = {};
    iteration_start_time = report_time;
    iteration_baseline = report_memory;
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

  // Charge the step since the last lap to the open innermost row, without closing it.
  static void Interrupt()
  {
    if (stack.empty())
    {
      timer.Lap();
    }
    else
    {
      timer.Mark(stack.back(), false);
    }
  }

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
        timer.Mark(Timer::STARTUP);
      }
      else
      {
        Interrupt();
      }
      stack.push_back(i);
    }
  }

  ~BlockTimer()
  {
    if (count && !stack.empty())
    {
      timer.Mark(stack.back());
      stack.pop_back();
    }
  }

  static const Timer &GlobalTimer() { return timer; }

  // Before a report: account the open innermost row up to now, without closing it, and fix
  // the report time.
  static void Snapshot()
  {
    Interrupt();
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
    timer.SetStartBaseline(memory_reporting::GetPeakMemory());
    stack.clear();
    started = false;
    perf::detail::ResetValues();
  }
};

}  // namespace palace

#endif  // PALACE_UTILS_TIMER_HPP
