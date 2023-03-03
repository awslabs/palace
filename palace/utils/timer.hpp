// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TIMER_HPP
#define PALACE_TIMER_HPP

#include <chrono>
#include <string>
#include <vector>
#include "utils/communication.hpp"

namespace palace
{

//
// A timer class for profiling.
//
class Timer
{
public:
  enum Index
  {
    INIT = 0,
    CONSTRUCT,
    SOLVE,
    ESTIMATION,
    POSTPRO,
    IO,
    TOTAL,
    NUMTIMINGS
  };

  using Clock = std::chrono::steady_clock;
  using Duration = std::chrono::duration<double>;

private:
  const typename Clock::time_point start_time;
  typename Clock::time_point last_lap_time;
  std::vector<Duration> data;
  std::vector<double> data_min, data_max, data_avg;

public:
  Timer() : start_time(Now()), last_lap_time(start_time), data((int)NUMTIMINGS) {}

  // Get current time.
  typename Clock::time_point Now() const { return Clock::now(); }

  // Stopwatch functionality.
  Duration Lap()
  {
    auto temp_time = last_lap_time;
    last_lap_time = Now();
    return last_lap_time - temp_time;
  }

  // Get time since start.
  Duration TimeFromStart() const { return Now() - start_time; }

  // Access the timing data.
  Duration &init_time = data[INIT];
  Duration &construct_time = data[CONSTRUCT];
  Duration &solve_time = data[SOLVE];
  Duration &estimation_time = data[ESTIMATION];
  Duration &postpro_time = data[POSTPRO];
  Duration &io_time = data[IO];
  Duration &total_time = data[TOTAL];

  // Access the reduced timing data.
  double GetMinTime(Index idx) const { return data_min[idx]; }
  double GetMaxTime(Index idx) const { return data_max[idx]; }
  double GetAvgTime(Index idx) const { return data_avg[idx]; }

  // Reduce timing information across MPI ranks.
  void Reduce(MPI_Comm comm)
  {
    const std::size_t ntimes = data.size();
    data[TOTAL] = TimeFromStart();
    data_min.resize(ntimes);
    data_max.resize(ntimes);
    data_avg.resize(ntimes);
    for (std::size_t i = 0; i < ntimes; i++)
    {
      data_min[i] = data_max[i] = data_avg[i] = data[i].count();
    }
    Mpi::GlobalMin(ntimes, data_min.data(), comm);
    Mpi::GlobalMax(ntimes, data_max.data(), comm);
    Mpi::GlobalSum(ntimes, data_avg.data(), comm);
    int np = Mpi::Size(comm);
    for (std::size_t i = 0; i < ntimes; i++)
    {
      data_avg[i] /= np;
    }
  }

  // Print timing information. Assumes the data has already been reduced.
  void Print(MPI_Comm comm) const
  {
    // clang-format off
    constexpr int p = 3;   // Floating point precision
    constexpr int w = 12;  // Total column width
    Mpi::Print(comm,
               "\n"
               "Elapsed Time Report (s)   {:>{}s}{:>{}s}{:>{}s}\n"
               "=========================={}\n"
               "Initialization            {:{}.{}f}{:{}.{}f}{:{}.{}f}\n"
               "Operator Construction     {:{}.{}f}{:{}.{}f}{:{}.{}f}\n"
               "Solve                     {:{}.{}f}{:{}.{}f}{:{}.{}f}\n"
               "Estimation                {:{}.{}f}{:{}.{}f}{:{}.{}f}\n"
               "Postprocessing            {:{}.{}f}{:{}.{}f}{:{}.{}f}\n"
               "Disk I/O                  {:{}.{}f}{:{}.{}f}{:{}.{}f}\n"
               "--------------------------{}\n"
               "Total Simulation          {:{}.{}f}{:{}.{}f}{:{}.{}f}\n",
               "Min.", w, "Max.", w, "Avg.", w,
               std::string(3 * w, '='),
               data_min[INIT], w, p, data_max[INIT], w, p, data_avg[INIT], w, p,
               data_min[CONSTRUCT], w, p, data_max[CONSTRUCT], w, p, data_avg[CONSTRUCT], w, p,
               data_min[SOLVE], w, p, data_max[SOLVE], w, p, data_avg[SOLVE], w, p,
               data_min[ESTIMATION], w, p, data_max[ESTIMATION], w, p, data_avg[ESTIMATION], w, p,
               data_min[POSTPRO], w, p, data_max[POSTPRO], w, p, data_avg[POSTPRO], w, p,
               data_min[IO], w, p, data_max[IO], w, p, data_avg[IO], w, p,
               std::string(3 * w, '-'),
               data_min[TOTAL], w, p, data_max[TOTAL], w, p, data_avg[TOTAL], w, p);
    // clang-format on
  }
};

}  // namespace palace

#endif  // PALACE_TIMER_HPP
