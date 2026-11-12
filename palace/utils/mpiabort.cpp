// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// MPI_Abort through the MPI profiling interface, to write the failure report first (see
// perf::FailureReport, CMake option PALACE_WITH_MPI_ABORT_HOOK). Compiled into the
// executables only, not into libpalace, so that other users of the library keep their
// MPI_Abort. The symbol is weak, so that a static tool library's own MPI_Abort takes
// precedence at link time, and the wrapper calls the next definition of MPI_Abort: a
// profiling tool loaded after the executable (for example with LD_PRELOAD) still receives
// the call. Without dynamic symbols (static MPI), it calls PMPI_Abort.

#include <dlfcn.h>
#include <mpi.h>
#include <fmt/format.h>
#include "utils/failurereport.hpp"

extern "C" __attribute__((weak)) int MPI_Abort(MPI_Comm comm, int errorcode)
{
  const auto fallback = fmt::format("MPI_Abort called with error code {}", errorcode);
  palace::perf::ReportFailure({}, fallback);
  using F = decltype(&PMPI_Abort);
  void *next = dlsym(RTLD_NEXT, "MPI_Abort");
  return (next ? reinterpret_cast<F>(next) : &PMPI_Abort)(comm, errorcode);
}
