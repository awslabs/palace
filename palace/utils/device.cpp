// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "device.hpp"
#include "communication.hpp"

#include <mfem.hpp>
#if defined(MFEM_USE_UMPIRE)
#include <umpire/ResourceManager.hpp>
#include <umpire/strategy/QuickPool.hpp>
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

void ConfigureMemory(const std::string &device_str, int device_id)
{
#if defined(MFEM_USE_UMPIRE)
  // MFEM built with Umpire defaults to the Umpire host allocator (HOST_UMPIRE), which is
  // slow (mfem/mfem#5122), and to its own device pool, while hypre creates a second device
  // pool. Select plain host memory and, on GPUs, one device pool shared with hypre (see
  // hypre::Initialize). The debug device needs MFEM's debug memory types, so leave it
  // alone.
  if (device_str.find("debug") != std::string::npos)
  {
    return;
  }
  if (device_str.find("cuda") != std::string::npos ||
      device_str.find("hip") != std::string::npos)
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

}  // namespace palace::utils
