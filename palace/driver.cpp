// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "driver.hpp"

#include <fstream>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "drivers/basesolver.hpp"
#include "drivers/boundarymodesolver.hpp"
#include "drivers/drivensolver.hpp"
#include "drivers/eigensolver.hpp"
#include "drivers/electrostaticsolver.hpp"
#include "drivers/magnetostaticsolver.hpp"
#include "drivers/transientsolver.hpp"
#include "fem/mesh.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

std::unique_ptr<BaseSolver> MakeSolver(const IoData &iodata, bool root, int size,
                                       int omp_threads, const char *git_tag)
{
  switch (iodata.problem.type)
  {
    case ProblemType::DRIVEN:
      return std::make_unique<DrivenSolver>(iodata, root, size, omp_threads, git_tag);
    case ProblemType::EIGENMODE:
      return std::make_unique<EigenSolver>(iodata, root, size, omp_threads, git_tag);
    case ProblemType::ELECTROSTATIC:
      return std::make_unique<ElectrostaticSolver>(iodata, root, size, omp_threads,
                                                   git_tag);
    case ProblemType::MAGNETOSTATIC:
      return std::make_unique<MagnetostaticSolver>(iodata, root, size, omp_threads,
                                                   git_tag);
    case ProblemType::TRANSIENT:
      return std::make_unique<TransientSolver>(iodata, root, size, omp_threads, git_tag);
    case ProblemType::BOUNDARYMODE:
      return std::make_unique<BoundaryModeSolver>(iodata, root, size, omp_threads, git_tag);
  }
  return nullptr;
}

// Axisymmetric (r, z) interpretation (config::ModelData::axisymmetric): the mesh must be
// two-dimensional and lie in the half-plane x >= 0 (a tolerance of the mesh's own length
// scale absorbs roundoff of nodes generated on the axis).
void VerifyAxisymmetricMesh(mfem::ParMesh &mesh)
{
  MFEM_VERIFY(mesh.Dimension() == 2 && mesh.SpaceDimension() == 2,
              "Model.Axisymmetric requires a two-dimensional (r, z) mesh, not (dim, "
              "space_dim) = ("
                  << mesh.Dimension() << ", " << mesh.SpaceDimension() << ")!");
  mfem::Vector bbmin, bbmax;
  mesh.GetBoundingBox(bbmin, bbmax, 2);
  double extent = std::max(bbmax[0] - bbmin[0], bbmax[1] - bbmin[1]);
  double xmin = mfem::infinity();
  const auto *nodes = mesh.GetNodes();
  if (nodes)
  {
    const auto *fespace = nodes->FESpace();
    for (int i = 0; i < fespace->GetNDofs(); i++)
    {
      xmin = std::min(xmin, (*nodes)(fespace->DofToVDof(i, 0)));
    }
  }
  else
  {
    for (int i = 0; i < mesh.GetNV(); i++)
    {
      xmin = std::min(xmin, mesh.GetVertex(i)[0]);
    }
  }
  Mpi::GlobalMin(1, &xmin, mesh.GetComm());
  Mpi::GlobalMax(1, &extent, mesh.GetComm());
  MFEM_VERIFY(
      xmin >= -1.0e-12 * extent,
      "Model.Axisymmetric requires every mesh node to satisfy x = r >= 0 (found x = "
          << xmin << ")!");
}

std::vector<std::unique_ptr<Mesh>> LoadMesh(IoData &iodata, MPI_Comm comm,
                                            const BaseSolver &solver)
{
  std::vector<std::unique_ptr<Mesh>> mesh;
  BlockTimer bt(Timer::INIT);
  auto smesh = mesh::Load(iodata, comm);
  solver.Preprocess(iodata, smesh, comm);
  std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
  mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), comm));
  mesh::RefineMesh(iodata, mfem_mesh);
  Mpi::Print(comm, "\n");
  memory_reporting::PrintMemoryUsage(comm, memory_reporting::GetCurrentMemoryStats(comm));
  memory_reporting::PrintMemoryUsage(comm,
                                     memory_reporting::GetCurrentNodeMemoryStats(comm));
  for (auto &m : mfem_mesh)
  {
    if (iodata.model.axisymmetric)
    {
      VerifyAxisymmetricMesh(*m);
    }
    mesh.push_back(std::make_unique<Mesh>(std::move(m)));
    mesh.back()->SetAxisymmetric(iodata.model.axisymmetric);
  }
  if (iodata.model.axisymmetric)
  {
    Mpi::Print(comm, "\nAxisymmetric (r, z) model: x = r >= 0, every integral carries the "
                     "revolution measure 2 pi r\n");
  }
  return mesh;
}

}  // namespace

void Run(IoData &iodata, MPI_Comm comm, int omp_threads, const char *git_tag)
{
  const bool world_root = Mpi::Root(comm);
  const int world_size = Mpi::Size(comm);

  auto solver = MakeSolver(iodata, world_root, world_size, omp_threads, git_tag);
  MFEM_VERIFY(solver, "Unknown problem type in palace::Run!");

  auto mesh = LoadMesh(iodata, comm, *solver);
  solver->SolveEstimateMarkRefine(mesh);

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

void RunSurfaceResponsePreflight(IoData &iodata, MPI_Comm comm, int omp_threads,
                                 const char *git_tag)
{
  auto solver = MakeSolver(iodata, Mpi::Root(comm), Mpi::Size(comm), omp_threads, git_tag);
  MFEM_VERIFY(solver, "Unknown problem type in surface-response preflight!");
  auto mesh = LoadMesh(iodata, comm, *solver);
  MFEM_VERIFY(!mesh.empty(), "Surface-response preflight produced no mesh!");
  WriteSurfaceResponseRequirements(
      iodata, *mesh.back(),
      (fs::path(iodata.problem.output) / "surface-response-requirements.json").string());

  // Preserve the same structured timing and memory metadata as an ordinary solve so
  // geometry-only preflight benchmarks do not need to scrape formatted terminal output.
  const auto peak_mem = memory_reporting::GetPeakMemoryStats(comm);
  const auto peak_node_mem = memory_reporting::GetPeakNodeMemoryStats(comm);
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

void RunMeshStatistics(IoData &iodata, MPI_Comm comm, int omp_threads, const char *git_tag)
{
  auto solver = MakeSolver(iodata, Mpi::Root(comm), Mpi::Size(comm), omp_threads, git_tag);
  MFEM_VERIFY(solver, "Unknown problem type in mesh statistics!");
  auto mesh = LoadMesh(iodata, comm, *solver);
  MFEM_VERIFY(!mesh.empty(), "Mesh statistics preprocessing produced no mesh!");
  mfem::H1_FECollection collection(iodata.solver.order, mesh.back()->Dimension());
  mfem::ParFiniteElementSpace space(&mesh.back()->Get(), &collection);
  const auto dofs = space.GlobalTrueVSize();
  const auto elements = mesh.back()->Get().GetGlobalNE();
  if (Mpi::Root(comm))
  {
    const nlohmann::json data = {
        {"Version", 1},
        {"Scope", "Configured mesh after normal preprocessing and initial refinement"},
        {"Mesh", iodata.model.mesh},
        {"Order", iodata.solver.order},
        {"Ranks", Mpi::Size(comm)},
        {"GlobalElements", elements},
        {"H1TrueDOFs", dofs},
        {"CrackInternalBoundaryElements", iodata.model.crack_bdr_elements},
        {"RefineCrackElements", iodata.model.refine_crack_elements}};
    const auto path = fs::path(iodata.problem.output) / "mesh-statistics.json";
    std::ofstream output(path);
    MFEM_VERIFY(output, "Cannot open mesh-statistics output!");
    output << data.dump(2) << '\n';
    MFEM_VERIFY(output.good(), "Failed writing mesh-statistics output!");
    Mpi::Print(comm, "Post-preprocessing H1 DOFs: {}, elements: {}\n", dofs, elements);
  }
}

}  // namespace palace
