// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "device.hpp"
#include "communication.hpp"
#include "omp.hpp"

#include <mfem.hpp>

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

}  // namespace palace::utils
