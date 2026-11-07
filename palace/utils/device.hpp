// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef DEVICE_HPP
#define DEVICE_HPP

#include "communication.hpp"

namespace palace::utils
{

// Return the number of devices available.
int GetDeviceCount();

// Assign devices round-robin over MPI ranks.
int GetDeviceId(MPI_Comm comm, int ngpu);

// Return whether vector operations run on the host in a single thread (no device backend
// and one OpenMP thread).
bool IsSerialHost();

}  // namespace palace::utils

#endif  // DEVICE_HPP
