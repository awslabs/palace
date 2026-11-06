// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "utils/communication.hpp"
#include "utils/perfstats.hpp"

namespace palace::perf::detail
{

void AddReduction(double seconds, std::size_t bytes)
{
  if (registry.suspended)
  {
    return;
  }
  static const int tally = RegisterPhaseTally("Reduce");
  static const int counter =
      RegisterCounter("Communication/ReduceBytes", Reduce::STATS, Accumulate::SUM);
  AddTally(tally, seconds);
  Add(counter, static_cast<double>(bytes));
}

}  // namespace palace::perf::detail
