// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
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
#include "utils/versions.hpp"

#if defined(MFEM_USE_STRUMPACK)
#include <StrumpackConfig.hpp>
#endif

using namespace palace;

static consteval const char *GetPalaceGitTag()
{
#if defined(PALACE_GIT_COMMIT)
  return PALACE_GIT_COMMIT_ID;
#else
  return "UNKNOWN";
#endif
}

static void PrintPalaceVersion(MPI_Comm comm)
{
  Mpi::Print(comm, "Palace version: {}\nSchema version: {}\n", GetPalaceGitTag(),
             GetSchemaVersion());
#if defined(PALACE_BUILD_ID)
  Mpi::Print(comm, "Build: {}\n", PALACE_BUILD_ID);
#endif
  Mpi::Print(comm, "\nDependencies:\n");
  for (const auto &[name, version] : GetDependencyVersions())
  {
    Mpi::Print(comm, "  {}: {}\n", name, version);
  }
}

static consteval const char *GetPalaceCeedJitSourceDir()
{
#if defined(PALACE_LIBCEED_JIT_SOURCE)
  return PALACE_LIBCEED_JIT_SOURCE_DIR;
#else
  return "";
#endif
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
  if (!ceed_resource.starts_with(backend))
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
  // The total elapsed time is measured from process start-up (the timer is a static
  // object), so credit everything before palace::Run, including MPI_Init, to INIT.
  std::optional<BlockTimer> init_timer(std::in_place, Timer::INIT);

  // Parse command-line options, before MPI_Init so that a run can read its configuration
  // file and select its device first (see utils::InitDeviceBeforeMpi). The options are
  // acted on after MPI_Init.
  std::string_view exit_option;  // The first of -h, --help, and --version, if any
  bool dryrun = false;
  for (int i = 1; i < argc; i++)
  {
    const std::string_view arg(argv[i]);
    if (exit_option.empty() && (arg == "-h" || arg == "--help" || arg == "--version"))
    {
      exit_option = arg;
    }
    dryrun = dryrun || arg == "-dry-run" || arg == "--dry-run";
  }
  // A configuration file which cannot be opened is reported after MPI_Init.
  nlohmann::json config;
  int early_device_id = -1;
  if (exit_option.empty() && !dryrun && argc >= 2 && std::ifstream(argv[1]).is_open())
  {
    config = IoData::Parse(argv[1]);
    const nlohmann::json::json_pointer device("/Solver/Device");
    if (config.contains(device) && config[device] == "GPU")
    {
      early_device_id = utils::InitDeviceBeforeMpi();
    }
  }

  // Initialize MPI.
#if defined(MFEM_USE_STRUMPACK) && \
    (defined(STRUMPACK_USE_PTSCOTCH) || defined(STRUMPACK_USE_SLATE_SCALAPACK))
  Mpi::default_thread_required = MPI_THREAD_MULTIPLE;
#endif
  {
    // MPI_Init waits for all processes to start, which on many nodes can take much longer
    // than the setup itself, so it gets its own row.
    BlockTimer bt(Timer::MPI_STARTUP);
    Mpi::Init(argc, argv);
  }
  MPI_Comm world_comm = Mpi::World();
  bool world_root = Mpi::Root(world_comm);
  int world_size = Mpi::Size(world_comm);
  Mpi::Print(world_comm, "\n");

  auto Help = [executable_path = std::string_view(argv[0]), &world_comm]()
  {
    Mpi::Print(world_comm,
               "Usage: {} [OPTIONS] CONFIG_FILE\n\n"
               "Options:\n"
               "  -h, --help           Show this help message and exit\n"
               "  --version            Show version information and exit\n"
               "  -dry-run, --dry-run  Parse configuration file for errors and exit\n\n",
               executable_path.substr(executable_path.find_last_of('/') + 1));
  };
  if (exit_option == "--version")
  {
    PrintPalaceVersion(world_comm);
    return 0;
  }
  if (!exit_option.empty())
  {
    Help();
    return 0;
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
    config = IoData::ParseAndValidate(argv[argc - 1]);
    IoData iodata(config, false);
    if (!std::filesystem::exists(iodata.model.mesh))
    {
      MFEM_ABORT("Unable to open mesh file \"" << iodata.model.mesh << "\"!");
    }
    // Echo the fully-resolved configuration so users can inspect every filled-in default
    // without running the simulation.
    Mpi::Print(world_comm, "\nResolved configuration:\n{}\n",
               IoData::ConcretizeDefaults(iodata, config).dump(2));
    Mpi::Print(world_comm, "Dry-run: No errors detected in configuration file \"{}\"\n\n",
               argv[argc - 1]);
    return 0;
  }

  // Validate the configuration file, parsed before MPI_Init if it could be opened.
  PrintPalaceBanner(world_comm);
  if (config.is_null())
  {
    config = IoData::Parse(argv[1]);
  }
  IoData::Validate(config);
  IoData iodata(config, false);
  MakeOutputFolder(iodata, world_comm);

  // Write the resolved configuration to the output directory so users have a complete
  // record of every Palace decision (all defaults filled in).
  if (world_root)
  {
    iodata.WriteResolvedConfig(config, argv[1]);
  }

  // Initialize device + numerics.
  const int omp_threads = utils::ConfigureOmp();
  const int ngpu = utils::GetDeviceCount();
  const std::string device_str = ConfigureDevice(iodata.solver.device);
  int device_id = utils::GetDeviceId(world_comm, ngpu);
  if (iodata.solver.device == Device::GPU)
  {
    // Keep a device selected before MPI_Init, but warn if it differs from the assignment by
    // MPI, for example with a launcher variable inherited from an outer launch.
    bool differs = (early_device_id >= 0 && early_device_id != device_id);
    Mpi::GlobalOr(1, &differs, world_comm);
    if (differs)
    {
      Mpi::Warning(world_comm,
                   "The device selected from the node-local rank set by the MPI launcher "
                   "differs from the node-local MPI rank on some processes!\n");
    }
    if (early_device_id >= 0)
    {
      device_id = early_device_id;
    }
  }
  utils::ConfigureMemory(device_str, device_id);
  mfem::Device device(device_str, device_id);
  ConfigureCeedBackend(iodata.solver.ceed_backend);
  device.SetGPUAwareMPI(utils::QueryGPUAwareMPI());

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

  // Hand off to the driver library entry point. Keeping the lifted block out of
  // main keeps the executable body thin and lets the Catch2 regression harness
  // reuse the exact same code path on an in-process IoData. See palace/driver.hpp
  // for the preconditions it expects. INIT must end before palace::Run: its top-level
  // timers need an empty stack, otherwise its untimed gaps are credited to INIT and
  // BlockTimer::Finalize pops this frame.
  PrintPalaceInfo(world_comm, world_size, omp_threads, ngpu, device);
  init_timer.reset();
  palace::Run(iodata, world_comm, omp_threads, GetPalaceGitTag());

  // Finalize libCEED.
  ceed::Finalize();

  // Finalize SLEPc/PETSc.
#if defined(PALACE_WITH_SLEPC)
  slepc::Finalize();
#endif

  return 0;
}
