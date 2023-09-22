// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fespace.hpp"

#include <uuid.h>

namespace palace
{

std::string FiniteElementSpace::GenerateId()
{
  return uuids::to_string(uuids::uuid_system_generator{}());
}

}  // namespace palace
