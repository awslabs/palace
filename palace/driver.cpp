// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "driver.hpp"

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "drivers/basesolver.hpp"
#include "drivers/boundarymodesolver.hpp"
#include "drivers/drivensolver.hpp"
#include "drivers/eigensolver.hpp"
#include "drivers/electrostaticsolver.hpp"
#include "drivers/magnetostaticsolver.hpp"
#include "drivers/transientsolver.hpp"
#include "fem/mesh.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/outputdir.hpp"
#include "utils/timer.hpp"

namespace palace
{

void Run(IoData &iodata, MPI_Comm comm, int omp_threads, const char *git_tag)
{
  // Clear any BlockTimer state left over from a previous Run in the same
  // process. Production main() calls Run exactly once so this is a no-op
  // there; the regression test harness invokes Run many times and would
  // otherwise see cumulative timings in palace.json. See plan R2.
  BlockTimer::Reset();

  // Ensure the output folder exists, is writable, and that every rank sees
  // the same (possibly normalised) absolute/relative path string.
  MakeOutputFolder(iodata, comm);

  const bool world_root = Mpi::Root(comm);
  const int world_size = Mpi::Size(comm);

  // Construct the problem driver.
  std::unique_ptr<BaseSolver> solver = [&]() -> std::unique_ptr<BaseSolver>
  {
    switch (iodata.problem.type)
    {
      case ProblemType::DRIVEN:
        return std::make_unique<DrivenSolver>(iodata, world_root, world_size, omp_threads,
                                              git_tag);
      case ProblemType::EIGENMODE:
        return std::make_unique<EigenSolver>(iodata, world_root, world_size, omp_threads,
                                             git_tag);
      case ProblemType::ELECTROSTATIC:
        return std::make_unique<ElectrostaticSolver>(iodata, world_root, world_size,
                                                     omp_threads, git_tag);
      case ProblemType::MAGNETOSTATIC:
        return std::make_unique<MagnetostaticSolver>(iodata, world_root, world_size,
                                                     omp_threads, git_tag);
      case ProblemType::TRANSIENT:
        return std::make_unique<TransientSolver>(iodata, world_root, world_size,
                                                 omp_threads, git_tag);
      case ProblemType::BOUNDARYMODE:
        return std::make_unique<BoundaryModeSolver>(iodata, world_root, world_size,
                                                    omp_threads, git_tag);
    }
    return nullptr;
  }();
  MFEM_VERIFY(solver, "Unknown problem type in palace::Run!");

  // Load the serial mesh, apply problem-type-specific serial-stage
  // preprocessing, partition, distribute, and refine.
  std::vector<std::unique_ptr<Mesh>> mesh;
  {
    auto smesh = mesh::Load(iodata, comm);
    solver->Preprocess(iodata, smesh, comm);
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), comm));
    mesh::RefineMesh(iodata, mfem_mesh);
    Mpi::Print(comm, "\n");
    memory_reporting::PrintMemoryUsage(comm, memory_reporting::GetCurrentMemoryStats(comm));
    memory_reporting::PrintMemoryUsage(comm,
                                       memory_reporting::GetCurrentNodeMemoryStats(comm));
    for (auto &m : mfem_mesh)
    {
      mesh.push_back(std::make_unique<Mesh>(std::move(m)));
    }
  }

  // Run the problem driver.
  solver->SolveEstimateMarkRefine(mesh);

  // Print timing + peak-memory summary and record metadata.
  auto peak_mem = memory_reporting::GetPeakMemoryStats(comm);
  auto peak_node_mem = memory_reporting::GetPeakNodeMemoryStats(comm);
  Mpi::Print(comm, "\n");
  memory_reporting::PrintMemoryUsage(comm, peak_mem);
  memory_reporting::PrintMemoryUsage(comm, peak_node_mem);
  BlockTimer::Finalize(comm);
  BlockTimer::Print(comm);
  solver->SaveMetadata(BlockTimer::GlobalTimer());
  solver->SaveMetadata(peak_mem);
  solver->SaveMetadata(peak_node_mem);
  Mpi::Print(comm, "\n");
}

}  // namespace palace
