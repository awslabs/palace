// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FIXTURES_HPP
#define PALACE_FIXTURES_HPP

#include <filesystem>
#include <random>
#include "utils/communication.hpp"

namespace palace::test
{

namespace fs = std::filesystem;

enum class TempDirMode
{
  PerRank,  // Each MPI rank gets its own directory
  Shared    // All ranks share one directory
};

// Fixture that creates a temporary directory for each test and cleans it up automatically.
template <TempDirMode Mode = TempDirMode::PerRank>
struct TempDirFixture
{
  fs::path temp_dir;

  TempDirFixture()
  {
    int rank = Mpi::Rank(Mpi::World());
    int random_num;

    if constexpr (Mode == TempDirMode::PerRank)
    {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<> dis(10000, 99999);
      random_num = dis(gen);
      temp_dir = fs::temp_directory_path() / ("palace_test_" + std::to_string(random_num) +
                                              "_rank" + std::to_string(rank));
      fs::create_directories(temp_dir);
    }
    else
    {
      // Rank 0 generates random number and broadcasts to all ranks
      if (rank == 0)
      {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(10000, 99999);
        random_num = dis(gen);
      }
      Mpi::Broadcast(1, &random_num, 0, Mpi::World());

      temp_dir = fs::temp_directory_path() / ("palace_test_" + std::to_string(random_num));
      if (rank == 0)
      {
        fs::create_directories(temp_dir);
      }
      Mpi::Barrier(Mpi::World());
    }
  }

  ~TempDirFixture()
  {
    if constexpr (Mode == TempDirMode::PerRank)
    {
      fs::remove_all(temp_dir);
    }
    else
    {
      Mpi::Barrier(Mpi::World());
      if (Mpi::Rank(Mpi::World()) == 0)
      {
        fs::remove_all(temp_dir);
      }
    }
  }
};

}  // namespace palace::test

#endif  // PALACE_FIXTURES_HPP
