// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"

#include <random>
#include "utils/communication.hpp"

namespace palace::test
{

PerRankTempDir::PerRankTempDir()
{
  int rank = Mpi::Rank(Mpi::World());
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(10000, 99999);
  int random_num = dis(gen);
  temp_dir = fs::temp_directory_path() /
             ("palace_test_" + std::to_string(random_num) + "_rank" + std::to_string(rank));
  fs::create_directories(temp_dir);
}

PerRankTempDir::~PerRankTempDir()
{
  fs::remove_all(temp_dir);
}

SharedTempDir::SharedTempDir()
{
  int rank = Mpi::Rank(Mpi::World());
  int random_num;

  // Rank 0 generates random number and broadcasts to all ranks.
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

SharedTempDir::~SharedTempDir()
{
  Mpi::Barrier(Mpi::World());
  if (Mpi::Rank(Mpi::World()) == 0)
  {
    fs::remove_all(temp_dir);
  }
}

}  // namespace palace::test
