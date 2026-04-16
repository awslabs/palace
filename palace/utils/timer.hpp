// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_TIMER_HPP
#define PALACE_UTILS_TIMER_HPP

#include <chrono>
#include <stack>
#include <string>
#include <tuple>
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
  long start_memory;
  long last_memory;
  std::vector<long> mem_data;

public:
  Timer()
    : start_time(Now()), last_lap_time(start_time), data(NUM_TIMINGS), counts(NUM_TIMINGS),
      start_memory(memory_reporting::GetPeakMemory()), last_memory(start_memory),
      mem_data(NUM_TIMINGS, 0)
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

  // Snapshot peak RSS and return delta from last snapshot.
  long MemoryLap()
  {
    long current = memory_reporting::GetPeakMemory();
    long delta = current - last_memory;
    last_memory = current;
    return delta;
  }

  // Return peak RSS growth since timer creation.
  long MemoryFromStart() const { return memory_reporting::GetPeakMemory() - start_memory; }

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
};

class BlockTimer
{
  using Index = Timer::Index;

private:
  inline static Timer timer;
  inline static std::stack<Index> stack;
  bool count;

  // Stored reduction results (populated by Finalize).
  inline static std::vector<double> reduced_time_min;
  inline static std::vector<double> reduced_time_max;
  inline static std::vector<double> reduced_time_avg;
  inline static std::vector<double> reduced_mem_min;
  inline static std::vector<double> reduced_mem_max;
  inline static std::vector<double> reduced_mem_sum;
  inline static std::vector<double> reduced_node_mem_min;
  inline static std::vector<double> reduced_node_mem_max;
  inline static std::vector<double> reduced_node_mem_sum;
  inline static int num_nodes = 0;

  // Print a summary table with three columns. The row_fn callback produces the three
  // column values for a given timer index.
  template <typename RowFn>
  static void PrintTable(MPI_Comm comm, const std::string &title, const std::string &col1,
                         const std::string &col2, const std::string &col3, RowFn &&row_fn)
  {
    constexpr int w = 12;  // Data column width
    constexpr int h = 26;  // Left-hand side width
    // clang-format off
    Mpi::Print(comm, "\n{:<{}s}{:>{}s}{:>{}s}{:>{}s}\n",
               title, h, col1, w, col2, w, col3, w);
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
        auto [v1, v2, v3] = row_fn(i);
        // clang-format off
        Mpi::Print(comm, "{:<{}s}{:>{}s}{:>{}s}{:>{}s}\n",
                   timer.descriptions[i], h, v1, w, v2, w, v3, w);
        // clang-format on
      }
    }
  }

public:
  BlockTimer(Index i, bool count = true) : count(count)
  {
    // Start timing when entering the block, interrupting whatever we were timing before.
    if (count)
    {
      if (stack.empty())
      {
        timer.Lap();
        timer.MemoryLap();
      }
      else
      {
        timer.MarkTime(stack.top(), false);
        timer.MarkMemory(stack.top(), timer.MemoryLap());
      }
      stack.push(i);
    }
  }

  ~BlockTimer()
  {
    // When a BlockTimer is no longer in scope, record the time and memory growth.
    if (count && !stack.empty())
    {
      timer.MarkTime(stack.top());
      timer.MarkMemory(stack.top());
      stack.pop();
    }
  }

  // Read-only access the static Timer object.
  static const Timer &GlobalTimer() { return timer; }

  // Access stored per-rank memory reduction results (populated by Finalize).
  static const std::vector<double> &RankMemoryMin()
  {
    MFEM_VERIFY(IsFinalized(),
                "BlockTimer::Finalize() must be called before accessing results!");
    return reduced_mem_min;
  }
  static const std::vector<double> &RankMemoryMax()
  {
    MFEM_VERIFY(IsFinalized(),
                "BlockTimer::Finalize() must be called before accessing results!");
    return reduced_mem_max;
  }
  static const std::vector<double> &RankMemorySum()
  {
    MFEM_VERIFY(IsFinalized(),
                "BlockTimer::Finalize() must be called before accessing results!");
    return reduced_mem_sum;
  }

  // Access stored per-node memory reduction results (populated by Finalize).
  static const std::vector<double> &NodeMemoryMin()
  {
    MFEM_VERIFY(IsFinalized(),
                "BlockTimer::Finalize() must be called before accessing results!");
    return reduced_node_mem_min;
  }
  static const std::vector<double> &NodeMemoryMax()
  {
    MFEM_VERIFY(IsFinalized(),
                "BlockTimer::Finalize() must be called before accessing results!");
    return reduced_node_mem_max;
  }
  static const std::vector<double> &NodeMemorySum()
  {
    MFEM_VERIFY(IsFinalized(),
                "BlockTimer::Finalize() must be called before accessing results!");
    return reduced_node_mem_sum;
  }

  // Finalize timers and perform MPI reductions. Must be called before Print().
  static void Finalize(MPI_Comm comm)
  {
    // Drain any open timers.
    while (!stack.empty())
    {
      timer.MarkTime(stack.top());
      timer.MarkMemory(stack.top());
      stack.pop();
    }
    timer.MarkTime(Timer::TOTAL, timer.TimeFromStart());
    timer.MarkMemory(Timer::TOTAL, timer.MemoryFromStart());

    // Reduce timing data across ranks.
    reduced_time_min.resize(Timer::NUM_TIMINGS);
    reduced_time_max.resize(Timer::NUM_TIMINGS);
    reduced_time_avg.resize(Timer::NUM_TIMINGS);
    for (int i = Timer::INIT; i < Timer::NUM_TIMINGS; i++)
    {
      reduced_time_min[i] = reduced_time_max[i] = reduced_time_avg[i] =
          timer.Data((Timer::Index)i);
    }
    Mpi::GlobalMin(Timer::NUM_TIMINGS, reduced_time_min.data(), comm);
    Mpi::GlobalMax(Timer::NUM_TIMINGS, reduced_time_max.data(), comm);
    Mpi::GlobalSum(Timer::NUM_TIMINGS, reduced_time_avg.data(), comm);
    const int np = Mpi::Size(comm);
    for (int i = Timer::INIT; i < Timer::NUM_TIMINGS; i++)
    {
      reduced_time_avg[i] /= np;
    }

    // Reduce per-rank memory data across ranks.
    reduced_mem_min.resize(Timer::NUM_TIMINGS);
    reduced_mem_max.resize(Timer::NUM_TIMINGS);
    reduced_mem_sum.resize(Timer::NUM_TIMINGS);
    for (int i = Timer::INIT; i < Timer::NUM_TIMINGS; i++)
    {
      reduced_mem_min[i] = reduced_mem_max[i] = reduced_mem_sum[i] =
          static_cast<double>(timer.MemoryData((Timer::Index)i));
    }
    Mpi::GlobalMin(Timer::NUM_TIMINGS, reduced_mem_min.data(), comm);
    Mpi::GlobalMax(Timer::NUM_TIMINGS, reduced_mem_max.data(), comm);
    Mpi::GlobalSum(Timer::NUM_TIMINGS, reduced_mem_sum.data(), comm);

    // Reduce per-node memory data across nodes.
    reduced_node_mem_min.resize(Timer::NUM_TIMINGS);
    reduced_node_mem_max.resize(Timer::NUM_TIMINGS);
    reduced_node_mem_sum.resize(Timer::NUM_TIMINGS);
    for (int i = Timer::INIT; i < Timer::NUM_TIMINGS; i++)
    {
      auto stats = memory_reporting::ComputeNodeMemoryStats(
          "", timer.MemoryData((Timer::Index)i), comm);
      reduced_node_mem_min[i] = static_cast<double>(stats.min);
      reduced_node_mem_max[i] = static_cast<double>(stats.max);
      reduced_node_mem_sum[i] = static_cast<double>(stats.sum);
    }

    // Count nodes.
    MPI_Comm node_comm;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
    int node_rank = Mpi::Rank(node_comm);
    MPI_Comm_free(&node_comm);
    int is_leader = (node_rank == 0) ? 1 : 0;
    num_nodes = 0;
    MPI_Allreduce(&is_leader, &num_nodes, 1, MPI_INT, MPI_SUM, comm);
  }

  // Whether Finalize has been called.
  static bool IsFinalized() { return !reduced_time_min.empty(); }

  // Print timing and memory tables from stored reduction results.
  static void Print(MPI_Comm comm)
  {
    MFEM_VERIFY(IsFinalized(), "BlockTimer::Finalize() must be called before Print()!");
    // Timing table.
    constexpr int p = 3;
    PrintTable(comm, "Elapsed Time Report (s)", "Min.", "Max.", "Avg.",
               [&](int i) -> std::tuple<std::string, std::string, std::string>
               {
                 return {fmt::format("{:.{}f}", reduced_time_min[i], p),
                         fmt::format("{:.{}f}", reduced_time_max[i], p),
                         fmt::format("{:.{}f}", reduced_time_avg[i], p)};
               });

    // Memory table. Single node: per-rank min/max/total. Multi-node: per-node
    // min/max/total.
    const auto &m_min = (num_nodes == 1) ? reduced_mem_min : reduced_node_mem_min;
    const auto &m_max = (num_nodes == 1) ? reduced_mem_max : reduced_node_mem_max;
    const auto &m_sum = (num_nodes == 1) ? reduced_mem_sum : reduced_node_mem_sum;
    std::string title =
        (num_nodes == 1) ? "Peak Memory Growth per Rank" : "Peak Memory Growth per Node";
    PrintTable(comm, title, "Min.", "Max.", "Tot.",
               [&](int i) -> std::tuple<std::string, std::string, std::string>
               {
                 return {memory_reporting::FormatBytes(m_min[i]),
                         memory_reporting::FormatBytes(m_max[i]),
                         memory_reporting::FormatBytes(m_sum[i])};
               });
  }
};

}  // namespace palace

#endif  // PALACE_UTILS_TIMER_HPP
