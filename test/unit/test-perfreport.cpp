// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <thread>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/timer.hpp"

using namespace palace;
using Catch::Matchers::WithinAbs;

namespace
{
void Sleep(int ms)
{
  std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}
}  // namespace

TEST_CASE("Timer Iterations", "[perfreport][Serial]")
{
  BlockTimer::Reset();
  Sleep(5);
  {
    BlockTimer bt(Timer::IO);
    Sleep(5);
  }
  const auto &t = BlockTimer::GlobalTimer();
  CHECK(t.Counts(Timer::STARTUP) == 1);
  CHECK(t.Data(Timer::STARTUP) >= 0.004);
  // A row open at the report is accounted up to the report, and not twice.
  {
    BlockTimer bt(Timer::POSTPRO);
    Sleep(5);
    BlockTimer::Snapshot();
    const double t1 = t.Data(Timer::POSTPRO), total1 = t.Data(Timer::TOTAL);
    CHECK(t1 >= 0.004);
    CHECK(t.Counts(Timer::POSTPRO) == 0);
    BlockTimer::StartIteration();
    Sleep(5);
  }
  BlockTimer::Snapshot();
  CHECK(t.Iterations() == 2);
  CHECK(t.Counts(Timer::POSTPRO) == 1);
  CHECK(t.Data(Timer::POSTPRO) >= 0.004);
  CHECK(t.Counts(Timer::STARTUP) == 0);
  CHECK(t.RunCounts(Timer::STARTUP) == 1);
  CHECK_THAT(t.RunData(Timer::POSTPRO) - t.Data(Timer::POSTPRO), WithinAbs(0.005, 0.004));
  CHECK(t.RunCounts(Timer::TOTAL) == 2);
  BlockTimer::Reset();
}
