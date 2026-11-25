// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef DEVICE_HPP
#define DEVICE_HPP

#include <string>
#include "communication.hpp"

namespace palace::utils
{

// Name of the Umpire device memory pool shared by MFEM and hypre in GPU builds. MFEM stores
// the allocator name as a raw pointer, so this must have static storage duration.
inline constexpr const char *UMPIRE_DEVICE_POOL_NAME = "PALACE_DEVICE_POOL";

// Return the number of devices available.
int GetDeviceCount();

// Return the node-local rank of this process set by the MPI launcher, or -1 if none is set.
int GetLauncherLocalRank();

// Select the device from the node-local rank set by the MPI launcher and create its
// context, before MPI_Init. Return the device ID, or -1 if unavailable.
int InitDeviceBeforeMpi();

// Assign devices round-robin over MPI ranks.
int GetDeviceId(MPI_Comm comm, int ngpu);

// Return whether vector operations run on the host in a single thread (no device backend
// and one OpenMP thread).
bool IsSerialHost();

// Configure the default MFEM memory types for the given MFEM device string (for example
// "cpu" or "cuda"). For GPU devices, this creates the shared Umpire device pool on the GPU
// with the given ID. Must be called before the mfem::Device is constructed. Does nothing
// unless MFEM is built with Umpire.
void ConfigureMemory(const std::string &device_str, int device_id);

// Return whether MPI can communicate CUDA or HIP device buffers, as reported by the MPI
// library (false for CPU builds).
bool QueryGPUAwareMPI();

}  // namespace palace::utils

#endif  // DEVICE_HPP
