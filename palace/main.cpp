// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>
#include <mpi.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "driver.hpp"
#include "fem/libceed/ceed.hpp"
#include "linalg/hypre.hpp"
#include "linalg/slepc.hpp"
#include "utils/communication.hpp"
#include "utils/device.hpp"
#include "utils/iodata.hpp"
#include "utils/jsonschema.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"
#include "utils/timer.hpp"

#if defined(MFEM_USE_STRUMPACK)
#include <StrumpackConfig.hpp>
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

static std::string ConfigureDevice(Device device)
{
  std::string device_str;
  switch (device)
  {
    case Device::CPU:
      device_str = "cpu";
      break;
    case Device::GPU:
#if defined(MFEM_USE_CUDA)
      device_str = "cuda";
#elif defined(MFEM_USE_HIP)
      device_str = "hip";
#else
      Mpi::Warning("Palace must be built with either CUDA or HIP support for GPU device "
                   "usage, reverting to CPU!\n");
      device_str = "cpu";
#endif
      break;
    case Device::DEBUG:
      device_str = "cpu,debug";
      break;
  }
#if defined(MFEM_USE_OPENMP)
  device_str += ",omp";
#endif
  return device_str;
}

static void ConfigureCeedBackend(const std::string &ceed_backend)
{
  // Initialize libCEED (only after MFEM device is configured).
  std::string default_ceed_backend;
  if (mfem::Device::Allows(mfem::Backend::CUDA_MASK))
  {
    default_ceed_backend = "/gpu/cuda/magma";
  }
  else if (mfem::Device::Allows(mfem::Backend::HIP_MASK))
  {
    default_ceed_backend = "/gpu/hip/magma";
  }
  else if (mfem::Device::Allows(mfem::Backend::DEBUG_DEVICE))
  {
    default_ceed_backend = "/cpu/self/ref/serial";
  }
  else
  {
    default_ceed_backend = "/cpu/self";
  }
  const std::string &backend =
      !ceed_backend.empty() ? ceed_backend.c_str() : default_ceed_backend.c_str();
  ceed::Initialize(backend.c_str(), GetPalaceCeedJitSourceDir());

  // Check that the provided resource matches the requested one.
  std::string ceed_resource = ceed::Print();
  if (backend.compare(0, backend.length(), ceed_resource, 0, backend.length()))
  {
    Mpi::Warning(
        "libCEED is not using the requested backend!\nRequested \"{}\", got \"{}\"!\n",
        backend, ceed_resource);
  }
}

static void PrintPalaceBanner(MPI_Comm comm)
{
  Mpi::Print(comm, "_____________     _______\n"
                   "_____   __   \\____ __   /____ ____________\n"
                   "____   /_/  /  __ ` /  /  __ ` /  ___/  _ \\\n"
                   "___   _____/  /_/  /  /  /_/  /  /__/  ___/\n"
                   "  /__/     \\___,__/__/\\___,__/\\_____\\_____/\n\n");
}

static void PrintPalaceInfo(MPI_Comm comm, int np, int nt, int ngpu, mfem::Device &device)
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
#if defined(MFEM_USE_CUDA)
  const char *device_name = "CUDA";
#else
  const char *device_name = "HIP";
#endif
  Mpi::Print(comm, "\nDetected {:d} {} device{}{}", ngpu, device_name,
             (ngpu != 1) ? "s" : "",
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
  // Initialize MPI.
#if defined(MFEM_USE_STRUMPACK) && \
    (defined(STRUMPACK_USE_PTSCOTCH) || defined(STRUMPACK_USE_SLATE_SCALAPACK))
  Mpi::default_thread_required = MPI_THREAD_MULTIPLE;
#endif
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
               "  --version            Show version information and exit\n"
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
    if (argv_i == "--version")
    {
      Mpi::Print(world_comm, "Palace version: {}\nSchema version: {}\n", GetPalaceGitTag(),
                 GetSchemaVersion());
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
      auto config = IoData::ParseAndValidate(argv[argc - 1]);
      IoData iodata(config, false);
      if (!std::filesystem::exists(iodata.model.mesh))
      {
        MFEM_ABORT("Unable to open mesh file \"" << iodata.model.mesh << "\"!");
      }
      // Echo the fully-resolved configuration so users can inspect every filled-in
      // default without running the simulation.
      Mpi::Print(world_comm, "\nResolved configuration:\n{}\n",
                 IoData::ConcretizeDefaults(iodata, config).dump(2));
    }
    Mpi::Print(world_comm, "Dry-run: No errors detected in configuration file \"{}\"\n\n",
               argv[argc - 1]);
    return 0;
  }

  // Parse configuration file.
  PrintPalaceBanner(world_comm);
  auto config = IoData::ParseAndValidate(argv[1]);
  IoData iodata(config, false);
  MakeOutputFolder(iodata, world_comm);

  // Write the resolved configuration to the output directory so users have a complete
  // record of every Palace decision (all defaults filled in).
  if (world_root)
  {
    iodata.WriteResolvedConfig(config, argv[1]);
  }

  // Initialise device + numerics. The BlockTimer is scoped to this block so it
  // credits the INIT phase before palace::Run starts its own timing.
  int omp_threads, ngpu;
  std::optional<mfem::Device> device;
  {
    BlockTimer bt(Timer::INIT);
    omp_threads = utils::ConfigureOmp();
    ngpu = utils::GetDeviceCount();
    device.emplace(ConfigureDevice(iodata.solver.device),
                   utils::GetDeviceId(world_comm, ngpu));
    ConfigureCeedBackend(iodata.solver.ceed_backend);
#if defined(PALACE_WITH_GPU_AWARE_MPI)
    device->SetGPUAwareMPI(true);
#endif

    // Initialize Hypre and, optionally, SLEPc/PETSc.
    hypre::Initialize();
#if defined(PALACE_WITH_SLEPC)
    slepc::Initialize(argc, argv, nullptr, nullptr);
    if (PETSC_COMM_WORLD != world_comm)
    {
      Mpi::Print(world_comm, "Error: Problem during MPI initialization!\n\n");
      return 1;
    }
#endif
  }

  // Hand off to the driver library entry point. Keeping the lifted block out of
  // main keeps the executable body thin and lets the Catch2 regression harness
  // reuse the exact same code path on an in-process IoData. See palace/driver.hpp
  // for the preconditions it expects.
  PrintPalaceInfo(world_comm, world_size, omp_threads, ngpu, *device);
  palace::Run(iodata, world_comm, omp_threads, GetPalaceGitTag());

  // Finalize libCEED.
  ceed::Finalize();

  // Finalize SLEPc/PETSc.
#if defined(PALACE_WITH_SLEPC)
  slepc::Finalize();
#endif

  return 0;
}
