// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstdlib>
#include <optional>
#include <string>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "utils/device.hpp"

using namespace palace;

namespace
{

// Unsets the launcher variables read by utils::GetLauncherLocalRank for the scope of a
// test, and restores them afterwards.
class LauncherEnvironment
{
private:
  static constexpr const char *names[] = {"OMPI_COMM_WORLD_LOCAL_RANK",
                                          "MV2_COMM_WORLD_LOCAL_RANK",
                                          "MPI_LOCALRANKID",
                                          "PMI_LOCAL_RANK",
                                          "SLURM_LOCALID",
                                          "PALS_LOCAL_RANKID",
                                          "FLUX_TASK_LOCAL_ID"};
  std::vector<std::optional<std::string>> saved;

public:
  LauncherEnvironment()
  {
    for (const char *name : names)
    {
      const char *value = std::getenv(name);
      saved.push_back(value ? std::optional<std::string>(value) : std::nullopt);
      unsetenv(name);
    }
  }
  ~LauncherEnvironment()
  {
    for (std::size_t i = 0; i < saved.size(); i++)
    {
      if (saved[i])
      {
        setenv(names[i], saved[i]->c_str(), 1);
      }
      else
      {
        unsetenv(names[i]);
      }
    }
  }
};

}  // namespace

TEST_CASE("GetLauncherLocalRank", "[device][Serial]")
{
  LauncherEnvironment env;
  CHECK(utils::GetLauncherLocalRank() == -1);

  // Each launcher's variable is read, and the first one set wins.
  setenv("FLUX_TASK_LOCAL_ID", "4", 1);
  CHECK(utils::GetLauncherLocalRank() == 4);
  setenv("SLURM_LOCALID", "3", 1);
  CHECK(utils::GetLauncherLocalRank() == 3);
  setenv("OMPI_COMM_WORLD_LOCAL_RANK", "1", 1);
  CHECK(utils::GetLauncherLocalRank() == 1);

  // A value that is not a non-negative integer gives no rank, without falling back to the
  // variables after it.
  for (const char *value : {"", "x", "-1", "2x", " 2", "99999999999"})
  {
    setenv("OMPI_COMM_WORLD_LOCAL_RANK", value, 1);
    CHECK(utils::GetLauncherLocalRank() == -1);
  }
}
