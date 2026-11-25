// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "driver.hpp"

#include <fstream>
#include <memory>
#include <string>
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
#include "utils/outputfile.hpp"
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

std::string SurfaceResponsePreflightKey(const IoData &iodata)
{
  const bool electrostatic = iodata.problem.type == ProblemType::ELECTROSTATIC;
  const auto &request = electrostatic ? iodata.solver.electrostatic.response_correction
                                      : iodata.solver.surface_response_correction;
  if (!request || !request->preflight_only)
  {
    return {};
  }
  return electrostatic ? "Solver.Electrostatic.ResponseCorrection.PreflightOnly"
                       : "Solver.SurfaceResponseCorrection.PreflightOnly";
}

int RunSurfaceResponsePreflight(IoData &iodata, MPI_Comm comm, int omp_threads,
                                const char *git_tag)
{
  auto solver = MakeSolver(iodata, Mpi::Root(comm), Mpi::Size(comm), omp_threads, git_tag);
  MFEM_VERIFY(solver, "Unknown problem type in surface-response preflight!");
  auto mesh = LoadMesh(iodata, comm, *solver);
  MFEM_VERIFY(!mesh.empty(), "Surface-response preflight produced no mesh!");
  const auto outcome = WriteSurfaceResponseRequirements(
      iodata, *mesh.back(),
      (fs::path(iodata.problem.output) / "surface-response-requirements.json").string());
  return FinishSurfaceResponsePreflight(iodata, *solver, mesh, outcome);
}

int FinishSurfaceResponsePreflight(IoData &iodata, const BaseSolver &solver,
                                   std::vector<std::unique_ptr<Mesh>> &mesh,
                                   const SurfaceResponsePreflightOutcome &outcome)
{
  MPI_Comm comm = mesh.back()->GetComm();
  auto &request = iodata.problem.type == ProblemType::ELECTROSTATIC
                      ? iodata.solver.electrostatic.response_correction
                      : iodata.solver.surface_response_correction;
  MFEM_VERIFY(request, "Surface-response preflight requires a response correction!");
  // CUT-PLAN step 4b (lane v1step4b; the supervisor's ruling R1): a solve-form run stops
  // before its census here; the census is not built and the exit status is the manifest's.
  const std::string refusal = SurfaceResponsePreflightCensusRefusal(outcome, *request);
  const auto policy = request->geometry_failure_policy;
  if (!outcome.fail_closed_summary.empty())
  {
    // The supervisor's rulings (b) and R2 (2026-11-25): today's FailClosed line first;
    // with GeometryFailurePolicy absent (the preflight fails closed by default while a
    // solve-form run records raw, its own default, and reaches its census) the census is
    // then built as that run builds it; an explicit FailClosed (a solve form aborts too)
    // builds none.
    std::string written;
    for (const auto &file : outcome.written)
    {
      written += (written.empty() ? "" : " and ") + file;
    }
    Mpi::Print(
        comm,
        "\nSurface-response preflight ended: GeometryFailurePolicy FailClosed with {}; "
        "written: {} (written in full); exit status {:d}\n",
        outcome.fail_closed_summary, written,
        kSurfaceResponsePreflightFailClosedExitStatus);
    if (policy)
    {
      Mpi::Print(comm, "\nResolution census: not evaluated (GeometryFailurePolicy "
                       "\"FailClosed\" is set); \"GeometryFailurePolicy\": \"RecordRaw\" "
                       "gives it\n");
    }
    else if (!refusal.empty())
    {
      Mpi::Print(comm, "\n{}\n", refusal);
    }
    else
    {
      const auto census = solver.PreflightResolutionCensus(mesh, false);
      if (census.evaluated)
      {
        Mpi::Print(
            comm,
            "The resolution census above was built as a solve-form run builds it "
            "(GeometryFailurePolicy RecordRaw, the default of a solve: raw recorded for "
            "the failed geometry, {:d} feature(s))\n",
            census.fallback_features);
      }
    }
    return kSurfaceResponsePreflightFailClosedExitStatus;
  }
  if (!refusal.empty())
  {
    Mpi::Print(comm, "\n{}\n", refusal);
  }
  else
  {
    // CUT-PLAN step 4b (lane v1step4b): the resolution census of the first mesh, as a
    // corrected run computes it before its first solve (no system matrix is assembled and
    // nothing is solved; a wave port's boundary-mode matrices are built with the space
    // operator). An explicit FailClosed (ruling R1, G3): the census operator is built under
    // RecordRaw, which takes the same branches, and a fallback it records ends the census.
    const bool fail_closed =
        policy && *policy == config::ElectrostaticSolverData::ResponseCorrectionData::
                                 GeometryFailurePolicy::FAIL_CLOSED;
    if (fail_closed)
    {
      request->geometry_failure_policy = config::ElectrostaticSolverData::
          ResponseCorrectionData::GeometryFailurePolicy::RECORD_RAW;
    }
    solver.PreflightResolutionCensus(mesh, fail_closed);
    request->geometry_failure_policy = policy;
  }

  // Preserve the same structured timing and memory metadata as an ordinary solve so
  // geometry-only preflight benchmarks do not need to scrape formatted terminal output.
  const auto peak_mem = memory_reporting::GetPeakMemoryStats(comm);
  const auto peak_node_mem = memory_reporting::GetPeakNodeMemoryStats(comm);
  Mpi::Print(comm, "\n");
  memory_reporting::PrintMemoryUsage(comm, peak_mem);
  memory_reporting::PrintMemoryUsage(comm, peak_node_mem);
  BlockTimer::Finalize(comm);
  BlockTimer::Print(comm);
  solver.SaveMetadata(BlockTimer::GlobalTimer());
  solver.SaveMetadata(peak_mem);
  solver.SaveMetadata(peak_node_mem);
  Mpi::Print(comm, "\n");
  return 0;
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
        // null when unset (the default crack set of GetMeshCrackAttributes applies).
        {"CrackInternalBoundaryElements",
         iodata.model.crack_bdr_elements ? nlohmann::json(*iodata.model.crack_bdr_elements)
                                         : nlohmann::json(nullptr)},
        {"RefineCrackElements", iodata.model.refine_crack_elements}};
    const auto path = fs::path(iodata.problem.output) / "mesh-statistics.json";
    WriteOutputFile(path, [&](std::ostream &output) { output << data.dump(2) << '\n'; });
    Mpi::Print(comm, "Post-preprocessing H1 DOFs: {}, elements: {}\n", dofs, elements);
  }
}

}  // namespace palace
