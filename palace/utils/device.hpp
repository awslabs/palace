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

// Assign devices round-robin over MPI ranks.
int GetDeviceId(MPI_Comm comm, int ngpu);

// Configure the default MFEM memory types for the given MFEM device string (for example
// "cpu" or "cuda"). For GPU devices, this creates the shared Umpire device pool on the GPU
// with the given ID. Must be called before the mfem::Device is constructed. Does nothing
// unless MFEM is built with Umpire.
void ConfigureMemory(const std::string &device_str, int device_id);

}  // namespace palace::utils

#endif  // DEVICE_HPP
