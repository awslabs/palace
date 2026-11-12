// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>
#include "utils/communication.hpp"

namespace palace::test
{

namespace
{

// New uniquely named temporary directory (empty path on failure).
fs::path MakeUniqueTempDir(const std::string &prefix)
{
  std::string tmpl = (fs::temp_directory_path() / (prefix + "XXXXXX")).string();
  std::vector<char> buf(tmpl.begin(), tmpl.end());
  buf.push_back('\0');
  return mkdtemp(buf.data()) ? fs::path(buf.data()) : fs::path();
}

}  // namespace

PerRankTempDir::PerRankTempDir()
{
  int rank = Mpi::Rank(Mpi::World());
  temp_dir = MakeUniqueTempDir("palace_test_rank" + std::to_string(rank) + "_");
  if (temp_dir.empty())
  {
    throw std::runtime_error("Failed to create a temporary directory!");
  }
}

PerRankTempDir::~PerRankTempDir()
{
  fs::remove_all(temp_dir);
}

SharedTempDir::SharedTempDir()
{
  // Rank 0 creates the directory and broadcasts its path.
  std::string path;
  if (Mpi::Root(Mpi::World()))
  {
    path = MakeUniqueTempDir("palace_test_").string();
  }
  int len = static_cast<int>(path.size());
  Mpi::Broadcast(1, &len, 0, Mpi::World());
  path.resize(len);
  if (len > 0)
  {
    MPI_Bcast(path.data(), len, MPI_CHAR, 0, Mpi::World());
  }
  if (len == 0)
  {
    throw std::runtime_error("Failed to create a temporary directory!");
  }
  temp_dir = path;
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
