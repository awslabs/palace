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
    VERIFICATION,          // Opt-in verification checks
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
      "Verification",
      "Disk IO",
      "Total"});

  // Keys used in palace.json.
  static constexpr auto names = std::to_array<std::string_view>({
      "STARTUP", "INIT", "MESH_PREPROCESS", "CONSTRUCT", "WAVE_PORT", "WAVE_PORT_SOLVE",
      "WAVE_PORT_FIELD", "WAVE_PORT_POSTPRO", "KSP", "KSP_SETUP", "KSP_PRECONDITIONER",
      "KSP_COARSE_SETUP", "KSP_COARSE_SOLVE", "TS", "EPS", "DIV_FREE", "CONSTRUCT_PROM",
      "SOLVE_PROM", "ESTIMATION", "CONSTRUCT_ESTIMATOR", "SOLVE_ESTIMATOR", "ADAPTATION",
      "REBALANCE", "ADAPT_MESH_OUTPUT", "POSTPRO", "POSTPRO_FARFIELD", "POSTPRO_PARAVIEW",
      "POSTPRO_GRIDFUNCTION", "VERIFICATION", "IO", "TOTAL"});
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

  // Heap in use growth per phase (index NUM_TIMINGS: unattributed), sampled at timer
  // events. Allocator statistics can be slow, so a sample is skipped once the time spent in
  // the samples not forced exceeds a small fraction of the elapsed time. The growth since
  // the last sample taken then spans several phases and goes to a row open throughout.
  double heap_budget = 5.0e-4, heap_budget_floor = 1.0e-3;
  long iteration_start_heap = 0;
  long last_heap = 0;
  std::vector<long> heap_data;
  // Per phase: exclusive time since the last sample taken, and over the iteration, the
  // time sampled and the part of it whose growth was charged to the phase.
  std::vector<double> heap_pending, heap_time, heap_covered;
  TimePoint last_heap_event;
  bool heap_gap = false;  // A sample was skipped since the last one taken
  double heap_query_time = 0.0, budget_query_time = 0.0;
  long heap_samples = 0, heap_skipped = 0;

  // Largest RSS at the heap samples of the iteration taken while the peak RSS had not grown
  // since the iteration start: a lower bound of the peak of the iteration, if the peak RSS
  // of the process does not grow in it.
  long rss_sampled = 0;

public:
  Timer()
    : start_time(Now()), iteration_start_time(start_time), iteration_end_time(start_time),
      last_lap_time(start_time), data(NUM_TIMINGS), counts(NUM_TIMINGS),
      iteration_start_memory(memory_reporting::GetPeakMemory()),
      last_memory(iteration_start_memory), mem_data(NUM_TIMINGS, 0), cum_data(NUM_TIMINGS),
      cum_counts(NUM_TIMINGS, 0), heap_data(NUM_TIMINGS + 1, 0),
      heap_pending(NUM_TIMINGS + 1, 0.0), heap_time(NUM_TIMINGS + 1, 0.0),
      heap_covered(NUM_TIMINGS + 1, 0.0), last_heap_event(start_time)
  {
    // The first query also detects the allocator, a one-time cost (slow with a cold page
    // cache) which is not counted against the heap query budget.
    iteration_start_heap = last_heap = memory_reporting::GetHeapMemory();
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

  // Query the heap in use, counting the time spent in the report and, for samples which
  // can be skipped, against the heap query budget.
  long QueryHeap(bool budgeted = false)
  {
    const auto t0 = Now();
    const long heap = memory_reporting::GetHeapMemory();
    const double dt = Duration(Now() - t0).count();
    heap_query_time += dt;
    budget_query_time += budgeted ? dt : 0.0;
    return heap;
  }

  // At a timer event ending an interval of phase idx (NUM_TIMINGS: outside any phase),
  // sample the heap in use and charge its growth since the previous sample to idx, unless
  // over the query time budget. After skipped samples, the growth goes to gap_idx, the
  // innermost phase open since the previous sample taken. Returns whether it sampled.
  bool MarkHeap(int idx, int gap_idx, bool force = false)
  {
    const double dt = Duration(last_lap_time - last_heap_event).count();
    last_heap_event = last_lap_time;
    heap_pending[idx] += dt;
    heap_time[idx] += dt;
    if (!force &&
        budget_query_time >
            heap_budget * Duration(last_lap_time - start_time).count() + heap_budget_floor)
    {
      heap_skipped++;
      heap_gap = true;
      return false;
    }
    const long heap = QueryHeap(!force);
    heap_samples++;
    if (last_memory <= iteration_start_memory)
    {
      const auto t0 = Now();
      rss_sampled = std::max(rss_sampled, memory_reporting::GetCurrentMemory());
      const double dt = Duration(Now() - t0).count();
      heap_query_time += dt;
      budget_query_time += force ? 0.0 : dt;
    }
    const int to = heap_gap ? gap_idx : idx;
    heap_data[to] += heap - last_heap;
    heap_covered[to] += heap_pending[to];
    std::fill(heap_pending.begin(), heap_pending.end(), 0.0);
    heap_gap = false;
    last_heap = heap;
    return true;
  }

  // Fraction of the elapsed time, plus a floor in seconds, that the heap samples which can
  // be skipped may take.
  void SetHeapBudget(double fraction, double floor)
  {
    heap_budget = fraction;
    heap_budget_floor = floor;
  }

  // Heap growth (bytes) of a phase in the current iteration (NUM_TIMINGS: outside any
  // phase), the heap in use at the iteration start and at the last sample.
  long HeapData(int idx) const { return heap_data[idx]; }
  long HeapAtStart() const { return iteration_start_heap; }
  long HeapAtLastSample() const { return last_heap; }

  // Fraction of the exclusive time of a phase whose heap growth was charged to it (-1 if
  // it had none).
  double HeapCoverage(int idx) const
  {
    return (heap_time[idx] > 0.0) ? heap_covered[idx] / heap_time[idx] : -1.0;
  }

  // Largest RSS sampled in the current iteration while the peak RSS had not grown in it.
  long RssSampled() const { return rss_sampled; }

  // Heap queries of the current iteration: seconds (all queries, including forced samples,
  // the RSS samples and memory probes), samples taken and samples skipped.
  double HeapQueryTime() const { return heap_query_time; }
  long HeapSamples() const { return heap_samples; }
  long HeapSkipped() const { return heap_skipped; }

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
    std::fill(heap_data.begin(), heap_data.end(), 0);
    std::fill(heap_pending.begin(), heap_pending.end(), 0.0);
    std::fill(heap_time.begin(), heap_time.end(), 0.0);
    std::fill(heap_covered.begin(), heap_covered.end(), 0.0);
    last_heap_event = last_lap_time;
    heap_gap = false;
    heap_query_time = 0.0;
    heap_samples = heap_skipped = 0;
    rss_sampled = 0;
    iteration_start_heap = last_heap = QueryHeap();
  }
};

class BlockTimer
{
  using Index = Timer::Index;

private:
  inline static Timer timer;
  inline static std::vector<Index> stack;
  inline static bool started = false;
  // Number of stack entries which stayed open since the last heap sample taken.
  inline static std::size_t heap_open = 0;
  inline static int paused = 0;  // Pause depth (perf::ScopedPause): new timers do not count
  bool count;

  // Phase charged with the heap growth across skipped samples: the innermost one open
  // since the last sample taken (NUM_TIMINGS if none).
  static int HeapGapRow()
  {
    return (heap_open > 0) ? stack[heap_open - 1] : Timer::NUM_TIMINGS;
  }

  // After a timer event and its heap sample (taken or skipped).
  static void HeapEvent(bool taken)
  {
    heap_open = taken ? stack.size() : std::min(heap_open, stack.size());
  }

public:
  BlockTimer(Index i, bool count_it = true) : count(count_it && paused == 0)
  {
    // Start timing when entering the block, interrupting whatever we were timing before.
    if (count)
    {
      bool taken;
      if (stack.empty())
      {
        if (!started)
        {
          // Everything before the first timer is startup.
          started = true;
          timer.MarkTime(Timer::STARTUP);
          timer.MarkMemory(Timer::STARTUP);
          taken = timer.MarkHeap(Timer::STARTUP, HeapGapRow(), true);
        }
        else
        {
          timer.Lap();
          timer.MemoryLap();
          taken = timer.MarkHeap(Timer::NUM_TIMINGS, HeapGapRow());
        }
      }
      else
      {
        timer.MarkTime(stack.back(), false);
        timer.MarkMemory(stack.back(), timer.MemoryLap());
        taken = timer.MarkHeap(stack.back(), HeapGapRow());
      }
      stack.push_back(i);
      HeapEvent(taken);
    }
  }

  ~BlockTimer()
  {
    // When a BlockTimer is no longer in scope, record the time and memory growth.
    if (count && !stack.empty())
    {
      timer.MarkTime(stack.back());
      timer.MarkMemory(stack.back());
      const bool taken = timer.MarkHeap(stack.back(), HeapGapRow());
      stack.pop_back();
      HeapEvent(taken);
    }
  }

  // Read-only access the static Timer object.
  static const Timer &GlobalTimer() { return timer; }

  // Heap in use, counted in the heap query time of the report (for memory probes).
  static long QueryHeap() { return timer.QueryHeap(); }

  // Heap query budget of the samples which can be skipped (see Timer::SetHeapBudget). Reset
  // restores the default.
  static void SetHeapBudget(double fraction, double floor)
  {
    timer.SetHeapBudget(fraction, floor);
  }

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
  // of the process is not reset. Totals are kept for the Cumulative report.
  static void StartIteration();

  // Clear all process-wide timing and reporting state. Used by palace::Run's callers to
  // prevent statistics from bleeding between successive solves when Palace is embedded
  // (e.g. in the regression test harness).
  static void Reset();

  // Pause (true) or resume (false) the counting of new timers, see perf::ScopedPause.
  static void Pause(bool pause)
  {
    paused = pause ? paused + 1 : (paused > 0 ? paused - 1 : 0);
  }
  static void ClearPause() { paused = 0; }

  // Print timing and memory tables from stored reduction results.
  static void Print(MPI_Comm comm);
};

}  // namespace palace

#endif  // PALACE_UTILS_TIMER_HPP
