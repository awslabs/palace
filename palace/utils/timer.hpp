// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_TIMER_HPP
#define PALACE_UTILS_TIMER_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <string>
#include <string_view>
#include <vector>
#include "utils/communication.hpp"
#include "utils/memoryreporting.hpp"

namespace palace
{

//
// Timer classes for profiling time and peak memory growth.
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
    KSP_COARSE_SETUP,      // Linear solver coarse-level setup/factorization
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
    ADAPT_MESH_OUTPUT,     // Writing the adapted mesh
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
      "  Adapted Mesh Output",
      "Postprocessing",
      "  Far Fields",
      "  Paraview",
      "  Grid function",
      "Disk IO",
      "Total"});

  // Keys used in palace.json.
  static constexpr auto names = std::to_array<std::string_view>({
      "STARTUP", "INIT", "MESH_PREPROCESS", "CONSTRUCT", "WAVE_PORT", "WAVE_PORT_SOLVE",
      "WAVE_PORT_FIELD", "WAVE_PORT_POSTPRO", "KSP", "KSP_SETUP", "KSP_PRECONDITIONER",
      "KSP_COARSE_SETUP", "KSP_COARSE_SOLVE", "TS", "EPS", "DIV_FREE", "CONSTRUCT_PROM",
      "SOLVE_PROM", "ESTIMATION", "CONSTRUCT_ESTIMATOR", "SOLVE_ESTIMATOR", "ADAPTATION",
      "REBALANCE", "ADAPT_MESH_OUTPUT", "POSTPRO", "POSTPRO_FARFIELD", "POSTPRO_PARAVIEW",
      "POSTPRO_GRIDFUNCTION", "IO", "TOTAL"});
  // clang-format on
  static_assert(descriptions.size() == NUM_TIMINGS);
  static_assert(names.size() == NUM_TIMINGS);

private:
  TimePoint start_time;
  TimePoint iteration_start_time, iteration_end_time;
  TimePoint last_lap_time;
  std::vector<Duration> data;
  std::vector<int> counts;
  long iteration_start_memory;
  long last_memory;
  std::vector<long> mem_data;

  // Totals of previous iterations (rank-local).
  std::vector<Duration> cum_data;
  std::vector<long> cum_counts;

  // Largest RSS sampled in the iteration, a lower bound of its peak: at timer events while
  // the peak RSS had not grown since the iteration start (at most every rss_interval
  // seconds, reading the RSS costs a system call and a parse) and at memory probes.
  static constexpr double rss_interval = 0.1;
  long rss_sampled = 0;
  TimePoint last_rss_sample;

public:
  Timer()
    : start_time(Now()), iteration_start_time(start_time), iteration_end_time(start_time),
      last_lap_time(start_time), data(NUM_TIMINGS), counts(NUM_TIMINGS),
      iteration_start_memory(memory_reporting::GetPeakMemory()),
      last_memory(iteration_start_memory), mem_data(NUM_TIMINGS, 0), cum_data(NUM_TIMINGS),
      cum_counts(NUM_TIMINGS, 0), last_rss_sample(start_time)
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

  // Close the current iteration's TOTAL at the current time.
  void MarkIterationEnd()
  {
    iteration_end_time = Now();
    MarkTime(TOTAL, iteration_end_time - iteration_start_time);
  }

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

  // Provide map-like read-only access to the timing data of the current iteration.
  auto Data(Index idx) const { return data[idx].count(); }

  // Return number of times timer.MarkTime(idx) or TimerBlock b(idx) was called in the
  // current iteration.
  auto Counts(Index idx) const { return counts[idx]; }

  // Totals over all iterations, including the current one.
  double CumulativeData(Index idx) const
  {
    return (idx == Timer::TOTAL) ? Duration(iteration_end_time - start_time).count()
                                 : (cum_data[idx] + data[idx]).count();
  }
  long CumulativeCounts(Index idx) const { return cum_counts[idx] + counts[idx]; }

  // Snapshot peak RSS and return delta from last snapshot.
  long MemoryLap()
  {
    long current = memory_reporting::GetPeakMemory();
    long delta = current - last_memory;
    last_memory = current;
    return delta;
  }

  // Return peak RSS growth since the start of the current iteration.
  long MemoryFromStart() const
  {
    return memory_reporting::GetPeakMemory() - iteration_start_memory;
  }

  // Lap and record a memory delta for the given phase.
  long MarkMemory(Index idx) { return MarkMemory(idx, MemoryLap()); }

  // Record a given memory delta for the given phase (without lapping).
  long MarkMemory(Index idx, long delta)
  {
    if (idx == Timer::TOTAL)
    {
      mem_data[idx] = delta;
    }
    else
    {
      mem_data[idx] += delta;
    }
    return mem_data[idx];
  }

  // Provide read-only access to the memory data (bytes) for a given phase.
  auto MemoryData(Index idx) const { return mem_data[idx]; }

  // At a timer event, sample the RSS for the lower bound of the peak of the iteration if
  // the peak RSS has not grown in the iteration (and not sampled in the last rss_interval).
  void SampleRss()
  {
    if (last_memory <= iteration_start_memory &&
        Duration(last_lap_time - last_rss_sample).count() >= rss_interval)
    {
      last_rss_sample = last_lap_time;
      NoteRss(memory_reporting::GetCurrentMemory());
    }
  }

  // Record an RSS (bytes) measured in the current iteration, for example by a memory probe.
  void NoteRss(long rss) { rss_sampled = std::max(rss_sampled, rss); }

  // Largest RSS sampled in the current iteration.
  long RssSampled() const { return rss_sampled; }

  // Fold the current iteration into the cumulative totals and restart counting from zero.
  // The new iteration starts at the end of the previous one, so that iteration totals add
  // up. The memory baseline is the peak RSS at the call.
  void StartIteration()
  {
    for (int i = 0; i < NUM_TIMINGS; i++)
    {
      if (i != TOTAL)
      {
        cum_data[i] += data[i];
      }
      cum_counts[i] += counts[i];
      data[i] = Duration::zero();
      counts[i] = 0;
      mem_data[i] = 0;
    }
    iteration_start_time = iteration_end_time;
    iteration_start_memory = last_memory = memory_reporting::GetPeakMemory();
    rss_sampled = 0;
    last_rss_sample = last_lap_time;
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
      if (stack.empty())
      {
        if (!started)
        {
          // Everything before the first timer is startup.
          started = true;
          timer.MarkTime(Timer::STARTUP);
          timer.MarkMemory(Timer::STARTUP);
        }
        else
        {
          timer.Lap();
          timer.MemoryLap();
        }
      }
      else
      {
        timer.MarkTime(stack.back(), false);
        timer.MarkMemory(stack.back(), timer.MemoryLap());
      }
      timer.SampleRss();
      stack.push_back(i);
    }
  }

  ~BlockTimer()
  {
    // When a BlockTimer is no longer in scope, record the time and memory growth.
    if (count && !stack.empty())
    {
      timer.MarkTime(stack.back());
      timer.MarkMemory(stack.back());
      timer.SampleRss();
      stack.pop_back();
    }
  }

  // Read-only access the static Timer object.
  static const Timer &GlobalTimer() { return timer; }

  // Record an RSS measured in the current iteration (for memory probes), for the lower
  // bound of its peak.
  static void NoteRss(long rss) { timer.NoteRss(rss); }

  // The exclusive row currently being timed (STARTUP before the first timer, TOTAL between
  // top-level timers).
  static Index Current()
  {
    return !stack.empty() ? stack.back() : (started ? Timer::TOTAL : Timer::STARTUP);
  }

  // Close open timers and reduce all timers, tallies, counters and memory statistics
  // across ranks (collective). Results are available on all ranks via perf::Results().
  static void Finalize(MPI_Comm comm);

  // Whether Finalize has been called since the last Reset/StartIteration.
  static bool IsFinalized();

  // Start a new reporting iteration (AMR): counts and times restart from zero. The peak RSS
  // of the process is not reset. Totals are kept for the Cumulative block.
  static void StartIteration();

  // Clear all process-wide timing and reporting state. Used by palace::Run's callers to
  // prevent statistics from bleeding between successive solves when Palace is embedded
  // (e.g. in the regression test harness).
  static void Reset();

  // Print timing and memory tables from stored reduction results.
  static void Print(MPI_Comm comm);
};

}  // namespace palace

#endif  // PALACE_UTILS_TIMER_HPP
