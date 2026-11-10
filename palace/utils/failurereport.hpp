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
// Performance report on failure. While a FailureReport is alive, a process which fails
// first writes the performance report collected so far to palace.json, flagged as
// incomplete ("Complete": false and an "Error" block), then aborts all processes. The
// hooks:
//   - Uncaught exceptions (std::terminate), including MFEM errors when MFEM is built with
//     exceptions (its error action is set to throw).
//   - MPI_Abort, which MFEM errors call otherwise, through a definition in libpalace which
//     calls PMPI_Abort (CMake option PALACE_WITH_MPI_ABORT_HOOK). The symbol is weak, so
//     that a static tool library defining MPI_Abort takes precedence at link time, but the
//     dynamic linker ignores weakness: with a static libpalace the definition is part of
//     the executable and comes before libraries preloaded with LD_PRELOAD (a tool's
//     MPI_Abort is then bypassed); with a shared libpalace it applies through the library
//     order (libpalace before libmpi). It does not see the calls bound inside a shared MFEM
//     library on macOS (two-level namespace) or those of the Fortran bindings of MPI (used
//     by MUMPS and ARPACK), which call PMPI_Abort.
// Fatal signals, for example when the process runs out of memory, are not reported. Only
// the main thread reports: on another thread the hooks abort without a report, and a thread
// which fails while another one reports waits for it to abort.
//
// A failing process does not know whether the others fail too, so it waits up to a timeout
// for all of them to report a failure (as for a check on a global quantity). If they do,
// the report is reduced over all processes as usual (Error/Collective true). Otherwise no
// collective is safe: the root process, if it failed, writes its own values without any
// communication under "RootProcess" (Error/Collective false), and the other processes
// write nothing. Before aborting, which terminates all processes, they wait for about twice
// the timeout more: within that time, a failure of the root is still reported (the
// processes also join a collective report if the last ones fail meanwhile), and once the
// agreement shows that the root failed, they wait until the end of that time, and at least
// the timeout, while it writes. They abort right away only once the agreement shows that
// the root did not fail. When only non-root
// processes fail, palace.json keeps what was written before the failure. The wait is part
// of the innermost timer row, and is also reported as the tally Wait/FailureReport.
//
class FailureReport
{
public:
  // Reduces the report over the given communicator (all processes, or MPI_COMM_SELF on the
  // root process alone) and writes palace.json with the given "Error" block.
  using Writer = std::function<void(MPI_Comm, const nlohmann::json &)>;

  // Collective on comm. Only one FailureReport may be active at a time.
  FailureReport(MPI_Comm comm, Writer writer, double timeout = 10.0);

  // Collective on comm, except while unwinding an exception.
  ~FailureReport();

  FailureReport(const FailureReport &) = delete;
  FailureReport &operator=(const FailureReport &) = delete;
};

// Write the failure report for an error on this process, once, if a FailureReport is
// active (called by the hooks, which then abort). An empty message stands for the error
// message MFEM printed just before. Returns whether the report was reduced over all
// processes.
bool ReportFailure(std::string_view message);

}  // namespace palace::perf

#endif  // PALACE_UTILS_FAILURE_REPORT_HPP
