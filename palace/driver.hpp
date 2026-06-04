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
void Run(IoData &iodata, MPI_Comm comm, int omp_threads, const char *git_tag = nullptr);

}  // namespace palace

#endif  // PALACE_DRIVER_HPP
