// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstdlib>
#include <deque>
#include <iterator>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include "test-helpers.hpp"
#include "utils/device.hpp"

using namespace palace;

TEST_CASE("GetLauncherLocalRank", "[device][Serial]")
{
  // Unset the launcher variables, restoring them at the end of the test.
  std::deque<ScopedEnvironmentVariable> environment;
  for (const char *name : utils::LAUNCHER_LOCAL_RANK_VARIABLES)
  {
    environment.emplace_back(name, nullptr);
  }
  CHECK(utils::GetLauncherLocalRank() == -1);

  // The first variable that is set wins.
  const int num_variables =
      static_cast<int>(std::size(utils::LAUNCHER_LOCAL_RANK_VARIABLES));
  for (int i = num_variables - 1; i >= 0; i--)
  {
    setenv(utils::LAUNCHER_LOCAL_RANK_VARIABLES[i], std::to_string(i).c_str(), 1);
    CHECK(utils::GetLauncherLocalRank() == i);
  }

  // A value that is not a non-negative integer gives no rank, without falling back to the
  // variables after it.
  for (const char *value : {"", "x", "-1", "2x", " 2", "99999999999"})
  {
    setenv(utils::LAUNCHER_LOCAL_RANK_VARIABLES[0], value, 1);
    CHECK(utils::GetLauncherLocalRank() == -1);
  }
}
