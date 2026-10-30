// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "device.hpp"
#include "communication.hpp"

#include <mfem.hpp>
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
