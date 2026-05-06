// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVER_HPP
#define PALACE_DRIVER_HPP

#include <mpi.h>

namespace palace
{

class IoData;

// Run a Palace solve end-to-end on an already-parsed `iodata`.
//
// Preconditions:
//   - MPI has been initialised.
//   - `mfem::Device` has been constructed with the desired backend.
//   - libCEED, HYPRE, and (if built in) SLEPc/PETSc are initialised.
//   - `iodata.problem.output` is set; relative mesh paths in the config
//     are resolvable from the current working directory (or have been
//     rewritten to absolute paths by the caller).
//
// Behaviour:
//   - Creates/validates the output folder (broadcasts the resolved path
//     to all ranks).
//   - Resets the process-wide BlockTimer statics so repeated calls in
//     the same process (e.g. the regression test harness) do not
//     accumulate stale state.
//   - Runs the configured solver through SolveEstimateMarkRefine and
//     writes metadata + timings into the output folder.
//
// Mutates `iodata` in place (notably `problem.output` is normalised).
// MFEM exceptions propagate to the caller; Catch2 catches them when
// Palace is built with `MFEM_USE_EXCEPTIONS=ON`.
void Run(IoData &iodata, MPI_Comm comm, int omp_threads, const char *git_tag = nullptr);

}  // namespace palace

#endif  // PALACE_DRIVER_HPP
