// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "device.hpp"
#include "communication.hpp"
#include "omp.hpp"

#include <mfem.hpp>
#if defined(MFEM_USE_UMPIRE)
#include <umpire/ResourceManager.hpp>
#include <umpire/strategy/QuickPool.hpp>
#endif
#if __has_include(<mpi-ext.h>) && defined(OPEN_MPI)
#include <mpi-ext.h>
#endif

namespace palace::utils
{

int GetDeviceCount()
{
#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)
  return mfem::Device::GetDeviceCount();
#else
  return 0;
#endif
}

int GetDeviceId(MPI_Comm comm, int ngpu)
{
  // Assign devices round-robin over MPI ranks if GPU support is enabled.
#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)
  MPI_Comm node_comm;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                      &node_comm);
  int node_size = Mpi::Rank(node_comm);
  MPI_Comm_free(&node_comm);
  return node_size % ngpu;
#else
  return 0;
#endif
}

bool IsSerialHost()
{
  return !mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && GetMaxThreads() == 1;
}

void ConfigureMemory(const std::string &device_str, int device_id)
{
#if defined(MFEM_USE_UMPIRE)
  // MFEM built with Umpire defaults to the Umpire host allocator (HOST_UMPIRE), which is
  // slow (mfem/mfem#5122), and to its own device pool, while hypre creates a second device
  // pool. Select plain host memory and, on GPUs, one device pool shared with hypre (see
  // hypre::Initialize). MFEM always uses its debug memory types for the debug device, so
  // there is nothing to configure for it.
  if (device_str.find("debug") != std::string::npos)
  {
    return;
  }
  // MFEM also accepts "gpu" (and "ceed-gpu", "raja-gpu", ...) as aliases for CUDA or HIP.
  if (device_str.find("cuda") != std::string::npos ||
      device_str.find("hip") != std::string::npos ||
      device_str.find("gpu") != std::string::npos)
  {
    // Build the pool explicitly on this rank's GPU ("DEVICE::<id>"; Umpire aliases plain
    // "DEVICE" to "DEVICE::0"). MFEM and hypre both reuse an existing allocator with the
    // given name and do not take ownership of it.
    auto &rm = umpire::ResourceManager::getInstance();
    if (!rm.isAllocator(UMPIRE_DEVICE_POOL_NAME))
    {
      rm.makeAllocator<umpire::strategy::QuickPool>(
          UMPIRE_DEVICE_POOL_NAME, rm.getAllocator("DEVICE::" + std::to_string(device_id)));
    }
    mfem::MemoryManager::SetUmpireDeviceAllocatorName(UMPIRE_DEVICE_POOL_NAME);
    mfem::Device::SetMemoryTypes(mfem::MemoryType::HOST, mfem::MemoryType::DEVICE_UMPIRE);
  }
  else
  {
    mfem::Device::SetMemoryTypes(mfem::MemoryType::HOST, mfem::MemoryType::HOST);
  }
#endif
}

bool QueryGPUAwareMPI()
{
  // Open MPI reports support through its CUDA and ROCm extensions, MPICH and the MPIs based
  // on it (Cray MPICH, MVAPICH) through MPIX_GPU_query_support.
  int supported = 0;
#if defined(MFEM_USE_CUDA)
#if defined(OMPI_HAVE_MPI_EXT_CUDA) && OMPI_HAVE_MPI_EXT_CUDA
  supported = MPIX_Query_cuda_support();
#elif defined(MPIX_GPU_SUPPORT_CUDA)
  MPIX_GPU_query_support(MPIX_GPU_SUPPORT_CUDA, &supported);
#endif
#elif defined(MFEM_USE_HIP)
#if defined(OMPI_HAVE_MPI_EXT_ROCM) && OMPI_HAVE_MPI_EXT_ROCM
  supported = MPIX_Query_rocm_support();
#elif defined(MPIX_GPU_SUPPORT_HIP)
  MPIX_GPU_query_support(MPIX_GPU_SUPPORT_HIP, &supported);
#endif
#endif
  return supported;
}

}  // namespace palace::utils
