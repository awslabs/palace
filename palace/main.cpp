// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "drivers/drivensolver.hpp"
#include "drivers/eigensolver.hpp"
#include "drivers/electrostaticsolver.hpp"
#include "drivers/magnetostaticsolver.hpp"
#include "drivers/transientsolver.hpp"
#include "linalg/petsc.hpp"
#include "linalg/slepc.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

#if defined(MFEM_USE_OPENMP)
#include <omp.h>
#endif
#if defined(PALACE_GIT_COMMIT)
#include "gitversion.hpp"
#endif

using namespace palace;

void PrintBanner(MPI_Comm comm, int np, int nt, const char *git_tag)
{
  Mpi::Print(comm, "_____________     _______\n"
                   "_____   __   \\____ __   /____ ____________\n"
                   "____   /_/  /  __ ` /  /  __ ` /  ___/  _ \\\n"
                   "___   _____/  /_/  /  /  /_/  /  /__/  ___/\n"
                   "  /__/     \\___,__/__/\\___,__/\\_____\\_____/\n\n");
  if (git_tag)
  {
    Mpi::Print(comm, "Git changeset ID: {}\n", git_tag);
  }
  if (nt > 0)
  {
    Mpi::Print(comm, "Running with {:d} MPI process{} and {:d} OpenMP thread{}\n\n", np,
               (np > 1) ? "es" : "", nt, (nt > 1) ? "s" : "");
  }
  else
  {
    Mpi::Print(comm, "Running with {:d} MPI process{}\n\n", np, (np > 1) ? "es" : "");
  }
  Mpi::Barrier(comm);
}

int main(int argc, char *argv[])
{
  // Initialize MPI.
  Mpi::Init(argc, argv);
  MPI_Comm world_comm = Mpi::World();
  bool world_root = Mpi::Root(world_comm);
  int world_size = Mpi::Size(world_comm);
  Mpi::Print(world_comm, "\n");

  // Initialize timer.
  Timer timer;

  // Parse command-line options.
  std::vector<std::string_view> argv_sv(argv, argv + argc);
  bool dryrun = false;
  auto Help = [executable_path = argv_sv[0], &world_comm]()
  {
    Mpi::Print(world_comm,
               "Usage: {} [OPTIONS] CONFIG_FILE\n\n"
               "Options:\n"
               "  -h, --help           Show this help message and exit\n"
               "  -dry-run, --dry-run  Parse configuration file for errors and exit\n\n",
               executable_path.substr(executable_path.find_last_of('/') + 1));
  };
  for (int i = 1; i < argc; i++)
  {
    std::string_view argv_i = argv_sv.at(i);
    if ((argv_i == "-h") || (argv_i == "--help"))
    {
      Help();
      return 0;
    }
    if ((argv_i == "-dry-run") || (argv_i == "--dry-run"))
    {
      dryrun = true;
      continue;
    }
  }
  if (argc < 2)
  {
    Mpi::Print(world_comm, "Error: Invalid usage!\n\n");
    Help();
    return 1;
  }

  // Perform dry run: Parse configuration file for errors and exit.
  if (dryrun)
  {
    if (Mpi::Root(world_comm))
    {
      IoData iodata(argv[argc - 1], false);
    }
    Mpi::Print(world_comm, "Dry-run: No errors detected in configuration file \"{}\"\n\n",
               argv[argc - 1]);
    return 0;
  }

  // Parse configuration file.
  int num_thread = 0;
#if defined(MFEM_USE_OPENMP)
  const char *env = std::getenv("OMP_NUM_THREADS");
  if (env)
  {
    std::sscanf(env, "%d", &num_thread);
  }
  else
  {
    num_thread = 1;
    omp_set_num_threads(num_thread);
  }
#endif
#if defined(PALACE_GIT_COMMIT)
  const char *git_tag = GetGitCommit();
#else
  const char *git_tag = nullptr;
#endif
  PrintBanner(world_comm, world_size, num_thread, git_tag);
  IoData iodata(argv[1], false);

  // Initialize Hypre and PETSc, and optionally SLEPc.
  mfem::Hypre::Init();
  petsc::Initialize(argc, argv, nullptr, nullptr);
#if defined(PALACE_WITH_SLEPC)
  slepc::Initialize();
#endif
  if (PETSC_COMM_WORLD != world_comm)
  {
    Mpi::Print(world_comm, "Error: Problem during MPI initialization!\n\n");
    return 1;
  }

  // Initialize the problem driver.
  const std::unique_ptr<BaseSolver> solver = [&]() -> std::unique_ptr<BaseSolver>
  {
    switch (iodata.problem.type)
    {
      case config::ProblemData::Type::DRIVEN:
        return std::make_unique<DrivenSolver>(iodata, world_root, world_size, num_thread,
                                              git_tag);
        break;
      case config::ProblemData::Type::EIGENMODE:
        return std::make_unique<EigenSolver>(iodata, world_root, world_size, num_thread,
                                             git_tag);
        break;
      case config::ProblemData::Type::ELECTROSTATIC:
        return std::make_unique<ElectrostaticSolver>(iodata, world_root, world_size,
                                                     num_thread, git_tag);
        break;
      case config::ProblemData::Type::MAGNETOSTATIC:
        return std::make_unique<MagnetostaticSolver>(iodata, world_root, world_size,
                                                     num_thread, git_tag);
        break;
      case config::ProblemData::Type::TRANSIENT:
        return std::make_unique<TransientSolver>(iodata, world_root, world_size, num_thread,
                                                 git_tag);
        break;
      default:
        Mpi::Print(world_comm, "Error: Unsupported problem type!\n\n");
    }
    return nullptr;
  }();

  // Read the mesh from file, refine, partition, and distribute it. Then nondimensionalize
  // it and the input parameters.
  std::vector<std::unique_ptr<mfem::ParMesh>> mesh;
  mesh.push_back(std::make_unique<mfem::ParMesh>(
      mesh::ReadMesh(world_comm, iodata, false, true, true, false, timer)));
  iodata.NondimensionalizeInputs(*mesh[0]);
  mesh::RefineMesh(iodata, mesh);
  timer.init_time += timer.Lap() - timer.io_time;

  // Run the problem driver.
  int iter = 0;
  auto solver_output = solver->Solve(mesh, timer, iter++);




  timer.Reduce(world_comm);
  timer.Print(world_comm);
  solver->SaveMetadata(timer);
  Mpi::Print(world_comm, "\n");

  // Finalize PETSc.
#if defined(PALACE_WITH_SLEPC)
  slepc::Finalize();
#endif
  petsc::Finalize();

  return 0;
}
