// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <mpi.h>
#include <mfem.hpp>
#include "drivers/drivensolver.hpp"
#include "drivers/eigensolver.hpp"
#include "drivers/electrostaticsolver.hpp"
#include "drivers/magnetostaticsolver.hpp"
#include "drivers/transientsolver.hpp"
#include "fem/errorindicator.hpp"
#include "fem/libceed/utils.hpp"
#include "linalg/slepc.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

#if defined(MFEM_USE_OPENMP)
#include <omp.h>
#endif

using namespace palace;

static const char *GetPalaceGitTag()
{
#if defined(PALACE_GIT_COMMIT)
  static const char *commit = PALACE_GIT_COMMIT_ID;
#else
  static const char *commit = "UNKNOWN";
#endif
  return commit;
}

static const char *GetPalaceCeedJitSourceDir()
{
#if defined(PALACE_LIBCEED_JIT_SOURCE)
  static const char *path = PALACE_LIBCEED_JIT_SOURCE_DIR;
#else
  static const char *path = "";
#endif
  return path;
}

static int ConfigureOmp()
{
#if defined(MFEM_USE_OPENMP)
  int nt;
  const char *env = std::getenv("OMP_NUM_THREADS");
  if (env)
  {
    std::sscanf(env, "%d", &nt);
  }
  else
  {
    nt = 1;
    omp_set_num_threads(nt);
  }
  omp_set_dynamic(0);
  return nt;
#else
  return 0;
#endif
}

static int GetDeviceId(MPI_Comm comm)
{
  // Assign devices round-robin over MPI ranks if GPU support is enabled.
#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)
  MPI_Comm node_comm;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                      &node_comm);
  int node_size = Mpi::Rank(node_comm);
  MPI_Comm_free(&node_comm);
  return node_size % mfem::Device::GetNumGPU();
#else
  return 0;
#endif
}

static std::string ConfigureDeviceAndBackend(config::SolverData::Device device,
                                             const std::string &ceed_backend)
{
  // Configure
  std::string device_str, default_ceed_backend;
  switch (device)
  {
    case config::SolverData::Device::CPU:
      device_str = "cpu";
      default_ceed_backend = "/cpu/self";
      break;
    case config::SolverData::Device::GPU:
#if defined(MFEM_USE_CUDA)
      device_str = "cuda";
      default_ceed_backend = "/gpu/cuda/magma";
#elif defined(MFEM_USE_HIP)
      device_str = "hip";
      default_ceed_backend = "/gpu/hip/magma";
#else
      MFEM_ABORT(
          "Palace must be built with either CUDA or HIP support for GPU device usage!");
#endif
      break;
    case config::SolverData::Device::DEBUG:
      device_str = "cpu,debug";
      default_ceed_backend = "/cpu/self/ref";
      break;
  }
#if defined(MFEM_USE_OPENMP)
  device_str += ",omp";
#endif

  // Initialize libCEED.
  const std::string &backend =
      !ceed_backend.empty() ? ceed_backend.c_str() : default_ceed_backend.c_str();
  ceed::Initialize(backend.c_str(), GetPalaceCeedJitSourceDir());

  // Check that the provided resource matches the requested one.
  std::string ceed_resource = ceed::Print();
  if (backend.compare(0, backend.length(), ceed_resource, 0, backend.length()))
  {
    Mpi::Warning(
        "libCEED is not using the requested backend (requested \"{}\", got \"{}\")!\n",
        backend, ceed_resource);
  }

  return device_str;
}

static void PrintPalaceBanner(MPI_Comm comm)
{
  Mpi::Print(comm, "_____________     _______\n"
                   "_____   __   \\____ __   /____ ____________\n"
                   "____   /_/  /  __ ` /  /  __ ` /  ___/  _ \\\n"
                   "___   _____/  /_/  /  /  /_/  /  /__/  ___/\n"
                   "  /__/     \\___,__/__/\\___,__/\\_____\\_____/\n\n");
}

static void PrintPalaceInfo(MPI_Comm comm, int np, int nt, mfem::Device &device)
{
  if (std::strcmp(GetPalaceGitTag(), "UNKNOWN"))
  {
    Mpi::Print(comm, "Git changeset ID: {}\n", GetPalaceGitTag());
  }
  Mpi::Print(comm, "Running with {:d} MPI process{}", np, (np > 1) ? "es" : "");
  if (nt > 0)
  {
    Mpi::Print(comm, ", {:d} OpenMP thread{}", nt, (nt > 1) ? "s" : "");
  }
#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)
  int ngpu = mfem::Device::GetNumGPU();
#if defined(MFEM_USE_CUDA)
  const char *device_name = "CUDA";
#else
  const char *device_name = "HIP";
#endif
      Mpi::Print(comm, "\n{:d} detected {} device{}{}", ngpu, device_name,
                 (ngpu > 1) ? "s" : "",
                 mfem::Device::GetGPUAwareMPI() ? " (MPI is GPU aware)" : "");
#endif
  std::ostringstream resource(std::stringstream::out);
  resource << "\n";
  device.Print(resource);
  resource << "libCEED backend: " << ceed::Print();
  Mpi::Print(comm, "{}\n\n", resource.str());
  Mpi::Barrier(comm);
}

int main(int argc, char *argv[])
{
  // Initialize the timer.
  BlockTimer bt(Timer::INIT);

  // Initialize MPI.
  Mpi::Init(argc, argv);
  MPI_Comm world_comm = Mpi::World();
  bool world_root = Mpi::Root(world_comm);
  int world_size = Mpi::Size(world_comm);
  Mpi::Print(world_comm, "\n");

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
  PrintPalaceBanner(world_comm);
  IoData iodata(argv[1], false);

  // Initialize the MFEM device and configure libCEED backend.
  int omp_threads = ConfigureOmp(), device_id = GetDeviceId(world_comm);
  mfem::Device device(
      ConfigureDeviceAndBackend(iodata.solver.device, iodata.solver.ceed_backend),
      device_id);
#if defined(HYPRE_WITH_GPU_AWARE_MPI)
  device.SetGPUAwareMPI(true);
#endif

  // Initialize Hypre and, optionally, SLEPc/PETSc.
  mfem::Hypre::Init();
#if defined(PALACE_WITH_SLEPC)
  slepc::Initialize(argc, argv, nullptr, nullptr);
  if (PETSC_COMM_WORLD != world_comm)
  {
    Mpi::Print(world_comm, "Error: Problem during MPI initialization!\n\n");
    return 1;
  }
#endif

  // Initialize the problem driver.
  PrintPalaceInfo(world_comm, world_size, omp_threads, device);
  const auto solver = [&]() -> std::unique_ptr<BaseSolver>
  {
    switch (iodata.problem.type)
    {
      case config::ProblemData::Type::DRIVEN:
        return std::make_unique<DrivenSolver>(iodata, world_root, world_size, omp_threads,
                                              GetPalaceGitTag());
      case config::ProblemData::Type::EIGENMODE:
        return std::make_unique<EigenSolver>(iodata, world_root, world_size, omp_threads,
                                             GetPalaceGitTag());
      case config::ProblemData::Type::ELECTROSTATIC:
        return std::make_unique<ElectrostaticSolver>(iodata, world_root, world_size,
                                                     omp_threads, GetPalaceGitTag());
      case config::ProblemData::Type::MAGNETOSTATIC:
        return std::make_unique<MagnetostaticSolver>(iodata, world_root, world_size,
                                                     omp_threads, GetPalaceGitTag());
      case config::ProblemData::Type::TRANSIENT:
        return std::make_unique<TransientSolver>(iodata, world_root, world_size,
                                                 omp_threads, GetPalaceGitTag());
    }
    return nullptr;
  }();

  // Read the mesh from file, refine, partition, and distribute it. Then nondimensionalize
  // it and the input parameters.
  std::vector<std::unique_ptr<mfem::ParMesh>> mesh;
  mesh.push_back(mesh::ReadMesh(world_comm, iodata, false, true, true, false));
  iodata.NondimensionalizeInputs(*mesh[0]);
  mesh::RefineMesh(iodata, mesh);

  // Run the problem driver.
  solver->SolveEstimateMarkRefine(mesh);

  // Print timing summary.
  BlockTimer::Print(world_comm);
  solver->SaveMetadata(BlockTimer::GlobalTimer());
  Mpi::Print(world_comm, "\n");

  // Finalize libCEED.
  ceed::Finalize();

  // Finalize SLEPc/PETSc.
#if defined(PALACE_WITH_SLEPC)
  slepc::Finalize();
#endif

  return 0;
}
