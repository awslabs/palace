// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVER_HPP
#define PALACE_DRIVER_HPP

#include <memory>
#include <string>
#include <vector>
#include <mpi.h>

namespace palace
{

class BaseSolver;
class IoData;
class Mesh;
struct SurfaceResponsePreflightOutcome;

// Run a Palace solve end-to-end on an already-parsed `iodata`.
//
// Preconditions:
//   - MPI has been initialized.
//   - `mfem::Device` has been constructed with the desired backend.
//   - libCEED, HYPRE, and (if built in) SLEPc/PETSc are initialized.
//   - The output folder (`iodata.problem.output`) has been created by the
//     caller (e.g. via MakeOutputFolder); relative mesh paths in the config
//     are resolvable from the current working directory (or have been
//     rewritten to absolute paths by the caller).
void Run(IoData &iodata, MPI_Comm comm, int omp_threads, const char *git_tag = nullptr);

// The configuration key that selects the surface-response preflight exactly as the
// command-line option --surface-response-preflight does, for front ends that pass no
// option: "Solver.SurfaceResponseCorrection.PreflightOnly", or
// "Solver.Electrostatic.ResponseCorrection.PreflightOnly" for an electrostatic problem (the
// block of the configured problem type, the one the preflight reads). Returns the key's
// path when it is true, else an empty string.
std::string SurfaceResponsePreflightKey(const IoData &iodata);

// The process exit status of a surface-response preflight that wrote its complete manifest
// and then ended because GeometryFailurePolicy FailClosed recorded geometry fallbacks.
constexpr int kSurfaceResponsePreflightFailClosedExitStatus = 3;

// Load and preprocess the configured mesh, classify automatic surface-response
// neighborhoods, write a requirements manifest, and exit without assembling or solving.
// Returns the process exit status: 0, or kSurfaceResponsePreflightFailClosedExitStatus
// (then the timing and memory metadata are not written, as before this status existed).
int RunSurfaceResponsePreflight(IoData &iodata, MPI_Comm comm, int omp_threads,
                                const char *git_tag = nullptr);

// The end of a surface-response preflight after its manifest (`outcome`, written on
// `mesh`): the resolution census of the first mesh where a solve-form run reaches it
// (CUT-PLAN step 4b), one line where that run stops before it, and the exit status of
// RunSurfaceResponsePreflight. Collective.
int FinishSurfaceResponsePreflight(IoData &iodata, const BaseSolver &solver,
                                   std::vector<std::unique_ptr<Mesh>> &mesh,
                                   const SurfaceResponsePreflightOutcome &outcome);

// Load through the normal preprocessing path (including internal-boundary cracking)
// and report the resulting global H1 size without constructing a PDE operator.
void RunMeshStatistics(IoData &iodata, MPI_Comm comm, int omp_threads,
                       const char *git_tag = nullptr);

}  // namespace palace

#endif  // PALACE_DRIVER_HPP
