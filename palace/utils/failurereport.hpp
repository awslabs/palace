// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_FAILURE_REPORT_HPP
#define PALACE_UTILS_FAILURE_REPORT_HPP

#include <functional>
#include <string_view>
#include <mpi.h>
#include <nlohmann/json_fwd.hpp>

namespace palace::perf
{

//
// Performance report on failure. While a FailureReport is alive, a failure on the root
// process first writes the root's own performance report collected so far to palace.json
// (no communication, under "RootProcess"), flagged as incomplete ("Complete": false and an
// "Error" block), then aborts all processes. A failure on another process waits for a short
// grace period, so that the root can finish its report if it fails too (as for a check on a
// global quantity), then aborts: if the root did not fail, palace.json keeps what was
// written before. The hooks:
//   - Uncaught exceptions (std::terminate), including MFEM errors when MFEM is built with
//     exceptions (its error action is set to throw).
//   - MPI_Abort, which MFEM errors call otherwise, through a definition in the executables
//     (utils/mpiabort.cpp, CMake option PALACE_WITH_MPI_ABORT_HOOK), not in libpalace. It
//     calls the next definition of MPI_Abort, so that a profiling tool loaded at run time
//     still receives the call. It does not see the calls bound inside a shared MFEM library
//     on macOS (two-level namespace) or those of the Fortran bindings of MPI (used by MUMPS
//     and ARPACK).
// Fatal signals, for example when the process runs out of memory, are not reported. Only
// the main thread reports: on another thread the hooks abort without a report, and a thread
// which fails while another one reports waits for it to abort.
//
class FailureReport
{
public:
  // Writes palace.json from the root process's own values (reduced over MPI_COMM_SELF),
  // with the given "Error" block.
  using Writer = std::function<void(const nlohmann::json &)>;

  // grace: the time (s) a failing non-root process waits before it returns. Only one
  // FailureReport may be active at a time.
  FailureReport(MPI_Comm comm, Writer writer, double grace = 2.0);
  ~FailureReport();

  FailureReport(const FailureReport &) = delete;
  FailureReport &operator=(const FailureReport &) = delete;
};

// Handle a failure on this process, once, if a FailureReport is active (called by the
// hooks, which then abort): the root writes the failure report, another process waits for
// the grace period. An empty message stands for the error message MFEM printed just before,
// or the fallback. Returns whether this process wrote the report.
bool ReportFailure(std::string_view message, std::string_view fallback = "Unknown error");

}  // namespace palace::perf

#endif  // PALACE_UTILS_FAILURE_REPORT_HPP
