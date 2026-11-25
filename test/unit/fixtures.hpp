// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FIXTURES_HPP
#define PALACE_FIXTURES_HPP

#include <filesystem>
#include <string>

namespace palace::test
{

namespace fs = std::filesystem;

// New uniquely named temporary directory <prefix>XXXXXX, created by rank 0 and shared by
// all ranks of MPI_COMM_WORLD. Throws on all ranks if it cannot be created.
fs::path MakeSharedTempDir(const std::string &prefix);

// Fixture that creates a per-rank temporary directory for each test and cleans it up
// automatically.
struct PerRankTempDir
{
  fs::path temp_dir;
  PerRankTempDir();
  ~PerRankTempDir();
};

// Fixture that creates a single temporary directory shared across all MPI ranks for each
// test and cleans it up automatically.
struct SharedTempDir
{
  fs::path temp_dir;
  SharedTempDir();
  ~SharedTempDir();
};

}  // namespace palace::test

#endif  // PALACE_FIXTURES_HPP
