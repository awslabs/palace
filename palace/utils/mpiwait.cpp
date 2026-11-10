// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Completion of nonblocking point-to-point requests through the MPI profiling interface,
// for the performance report's Wait/MPI tally and communication wait (CMake option
// PALACE_WITH_MPI_WAIT_TALLY). Compiled into the executables only, not into libpalace.
//
// The wrappers see the C calls of the code linked into the executable (Palace and, for
// static or directly linked dependencies, the halo exchanges of MFEM and hypre and the
// sparse direct solvers), not the Fortran bindings (MUMPS) or dependencies bound to the
// MPI library directly (macOS two-level namespaces). Each wrapper calls the next definition
// of its function: a profiling tool loaded after the executable (for example with
// LD_PRELOAD) still receives every call, all of which are then timed. Without such a tool,
// requests which have already completed return from a test without reading the clock, and
// only the calls that block are timed.

#include <chrono>
#include <dlfcn.h>
#include <mpi.h>
#include "utils/communication.hpp"

namespace
{

using Clock = std::chrono::steady_clock;

// The next definition of an MPI function and whether it belongs to another tool (it
// differs from the MPI library's profiling entry point). Without dynamic symbols (static
// MPI), the profiling entry point.
template <typename F>
struct Next
{
  F f;
  bool tool;
  Next(const char *name, const char *pname, F pmpi) : f(pmpi), tool(false)
  {
    void *next = dlsym(RTLD_NEXT, name);
    if (next)
    {
      f = reinterpret_cast<F>(next);
      tool = (next != dlsym(RTLD_NEXT, pname));
    }
    palace::perf::detail::mpi_wait_tool |= tool;
  }
};

struct NextFunctions
{
  Next<decltype(&PMPI_Wait)> wait{"MPI_Wait", "PMPI_Wait", &PMPI_Wait};
  Next<decltype(&PMPI_Waitall)> waitall{"MPI_Waitall", "PMPI_Waitall", &PMPI_Waitall};
  Next<decltype(&PMPI_Waitany)> waitany{"MPI_Waitany", "PMPI_Waitany", &PMPI_Waitany};
  Next<decltype(&PMPI_Waitsome)> waitsome{"MPI_Waitsome", "PMPI_Waitsome", &PMPI_Waitsome};
};

const NextFunctions &Functions()
{
  palace::perf::detail::mpi_wait_seen = true;
  static const NextFunctions functions;
  return functions;
}

void AddWait(Clock::time_point t0)
{
  palace::perf::detail::AddWait(std::chrono::duration<double>(Clock::now() - t0).count());
}

}  // namespace

extern "C"
{
  int MPI_Wait(MPI_Request *request, MPI_Status *status)
  {
    const auto &next = Functions().wait;
    if (!next.tool)
    {
      int flag = 0;
      const int err = PMPI_Test(request, &flag, status);
      if (err != MPI_SUCCESS || flag)
      {
        return err;
      }
    }
    const auto t0 = Clock::now();
    const int err = next.f(request, status);
    AddWait(t0);
    return err;
  }

  int MPI_Waitall(int count, MPI_Request array_of_requests[],
                  MPI_Status array_of_statuses[])
  {
    const auto &next = Functions().waitall;
    if (!next.tool)
    {
      int flag = 0;
      const int err = PMPI_Testall(count, array_of_requests, &flag, array_of_statuses);
      if (err != MPI_SUCCESS || flag)
      {
        return err;
      }
    }
    const auto t0 = Clock::now();
    const int err = next.f(count, array_of_requests, array_of_statuses);
    AddWait(t0);
    return err;
  }

  int MPI_Waitany(int count, MPI_Request array_of_requests[], int *index,
                  MPI_Status *status)
  {
    const auto &next = Functions().waitany;
    if (!next.tool)
    {
      int flag = 0;
      const int err = PMPI_Testany(count, array_of_requests, index, &flag, status);
      if (err != MPI_SUCCESS || flag)
      {
        return err;
      }
    }
    const auto t0 = Clock::now();
    const int err = next.f(count, array_of_requests, index, status);
    AddWait(t0);
    return err;
  }

  int MPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount,
                   int array_of_indices[], MPI_Status array_of_statuses[])
  {
    const auto &next = Functions().waitsome;
    if (!next.tool)
    {
      const int err = PMPI_Testsome(incount, array_of_requests, outcount, array_of_indices,
                                    array_of_statuses);
      if (err != MPI_SUCCESS || *outcount != 0)
      {
        return err;
      }
    }
    const auto t0 = Clock::now();
    const int err =
        next.f(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
    AddWait(t0);
    return err;
  }
}
