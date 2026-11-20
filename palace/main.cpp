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
#include "utils/versions.hpp"

#if defined(MFEM_USE_STRUMPACK)
#include <StrumpackConfig.hpp>
#endif

// LAB (ordslow), not for a pull request.
#include <climits>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <malloc.h>
#include <sys/prctl.h>

using namespace palace;

namespace
{

// LAB (ordslow): the variant of a run is given in the configuration file as
// "STRUMPACKButterflyLevels": 1000 + flags (the value has no effect without HODLR
// compression), and read here before MPI_Init so that it can act at process start.
//   1: STRUMPACK prints its per-solve instrumentation (STRUMPACK_ORDSLOW=1)
//   2: no node-shared permutation
//   4: glibc malloc with fixed thresholds: no trimming, 64 MB top pad, mmap from 32 MB
//   8: no transparent huge pages for this process
//  16: STRUMPACK proportional mapping by factor size
//  64: instrumentation reduced over batches of 1024 solves instead of 64
// 128: the instrumentation waits for the receives of the redistribution one by one
int ordslow_flags = 0, ordslow_mallopt = 0, ordslow_prctl = 0;

void OrdslowLabSetup(int argc, char *argv[])
{
  if (argc < 2)
  {
    return;
  }
  std::ifstream f(argv[argc - 1]);
  if (!f)
  {
    return;
  }
  std::stringstream ss;
  ss << f.rdbuf();
  const std::string text = ss.str();
  const std::string key = "\"STRUMPACKButterflyLevels\"";
  auto pos = text.find(key);
  if (pos == std::string::npos)
  {
    return;
  }
  pos = text.find(':', pos + key.size());
  if (pos == std::string::npos)
  {
    return;
  }
  const long v = std::strtol(text.c_str() + pos + 1, nullptr, 10);
  if (v < 1000)
  {
    return;
  }
  ordslow_flags = static_cast<int>(v - 1000);
  if (ordslow_flags & 1)
  {
    // 128: also find the last source of every receiver of the right-hand side rows.
    setenv("STRUMPACK_ORDSLOW", (ordslow_flags & 128) ? "2" : "1", 1);
  }
  if (ordslow_flags & 64)
  {
    setenv("STRUMPACK_ORDSLOW_BATCH", "1024", 1);
  }
  if (ordslow_flags & 2)
  {
    setenv("PALACE_ORDSLOW_NOSHARE", "1", 1);
  }
  if (ordslow_flags & 4)
  {
    // Setting any of these also switches off the dynamic adjustment of the thresholds.
    ordslow_mallopt = mallopt(M_MMAP_THRESHOLD, 32 * 1024 * 1024) +
                      2 * mallopt(M_TRIM_THRESHOLD, -1) +
                      4 * mallopt(M_TOP_PAD, 64 * 1024 * 1024);
  }
  if (ordslow_flags & 8)
  {
    ordslow_prctl = prctl(PR_SET_THP_DISABLE, 1, 0, 0, 0);
  }
  if (ordslow_flags & 16)
  {
    setenv("PALACE_ORDSLOW_PROPMAP", "1", 1);
  }
}

}  // namespace

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

  // Initialize MPI.
#if defined(MFEM_USE_STRUMPACK) && \
    (defined(STRUMPACK_USE_PTSCOTCH) || defined(STRUMPACK_USE_SLATE_SCALAPACK))
  Mpi::default_thread_required = MPI_THREAD_MULTIPLE;
#endif
  OrdslowLabSetup(argc, argv);  // LAB (ordslow)
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
  if (ordslow_flags)
  {
    Mpi::Print(world_comm,
               "LAB ordslow: flags {:d} (1 instrumentation, 2 no shared permutation, 4 "
               "mallopt, 8 no THP, 16 mapping by factor size, 64 large batches, 128 "
               "receives one by one); mallopt "
               "results {:d} (7 = all three set), prctl result {:d}\n",
               ordslow_flags, ordslow_mallopt, ordslow_prctl);
  }

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
      PrintPalaceVersion(world_comm);
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

  // Initialize device + numerics.
  const int omp_threads = utils::ConfigureOmp();
  const int ngpu = utils::GetDeviceCount();
  const std::string device_str = ConfigureDevice(iodata.solver.device);
  const int device_id = utils::GetDeviceId(world_comm, ngpu);
  utils::ConfigureMemory(device_str, device_id);
  mfem::Device device(device_str, device_id);
  ConfigureCeedBackend(iodata.solver.ceed_backend);
#if defined(PALACE_WITH_GPU_AWARE_MPI)
  device.SetGPUAwareMPI(true);
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
