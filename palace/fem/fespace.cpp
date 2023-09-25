// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fespace.hpp"

#include "utils/omp.hpp"

namespace palace
{

std::size_t FiniteElementSpace::global_id = 0;

std::size_t FiniteElementSpace::GetId() const
{
  PalacePragmaOmp(critical(GetId))
  {
    if (!init || GetSequence() != prev_sequence)
    {
      id = global_id++;
      prev_sequence = GetSequence();
      init = true;
    }
  }
  return id;
}

}  // namespace palace
