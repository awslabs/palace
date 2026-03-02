// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FIXTURES_HPP
#define PALACE_FIXTURES_HPP

#include <filesystem>

namespace palace::test
{

namespace fs = std::filesystem;

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
