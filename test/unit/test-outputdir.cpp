// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include "fixtures.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/outputdir.hpp"

using namespace palace;

// The core-logic cases are tagged [Serial][Parallel] so they run in both the
// single-rank and the multi-rank sweeps. In the parallel sweep they additionally
// verify that the collective (barrier / split / free) is balanced across ranks
// and does not deadlock, and that a directory created on the (shared) filesystem
// is visible to every rank afterwards.

TEST_CASE_METHOD(palace::test::SharedTempDir, "EnsureDirectory creates a new directory",
                 "[outputdir][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto dir = temp_dir / "postpro";
  REQUIRE(!fs::exists(dir));

  EnsureDirectory(dir, comm);

  CHECK(fs::is_directory(dir));
  // The directory must be writable from this rank. Use a rank-unique filename so
  // ranks sharing the directory (parallel sweep, shared filesystem) do not race
  // on the same path.
  auto written = dir / ("written_" + std::to_string(Mpi::Rank(comm)) + ".txt");
  {
    std::ofstream f(written);
    f << "ok";
  }
  CHECK(fs::is_regular_file(written));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "EnsureDirectory succeeds on a pre-existing directory",
                 "[outputdir][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto dir = temp_dir / "postpro";
  // Pre-create the directory on the (shared) filesystem before the call.
  if (Mpi::Root(comm))
  {
    fs::create_directories(dir);
  }
  Mpi::Barrier(comm);
  REQUIRE(fs::is_directory(dir));

  // Regression guard: EnsureDirectory must NOT treat an already-existing
  // directory as an error. It relies on the non-throwing create_directories
  // overload and ignores the bool return, so a pre-existing directory returns
  // cleanly rather than aborting.
  EnsureDirectory(dir, comm);

  CHECK(fs::is_directory(dir));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "EnsureDirectory preserves existing directory contents",
                 "[outputdir][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto dir = temp_dir / "postpro";
  if (Mpi::Root(comm))
  {
    fs::create_directories(dir);
    std::ofstream f(dir / "existing.txt");
    f << "keep me";
  }
  Mpi::Barrier(comm);

  EnsureDirectory(dir, comm);

  // EnsureDirectory only ensures the directory exists; it must not wipe content.
  CHECK(fs::is_regular_file(dir / "existing.txt"));
  {
    std::ifstream f(dir / "existing.txt");
    std::string content{std::istreambuf_iterator<char>(f),
                        std::istreambuf_iterator<char>()};
    CHECK(content == "keep me");
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "EnsureDirectory creates nested directories",
                 "[outputdir][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto dir = temp_dir / "a" / "b" / "c";
  REQUIRE(!fs::exists(dir));

  EnsureDirectory(dir, comm);

  CHECK(fs::is_directory(dir));
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "EnsureDirectory is idempotent",
                 "[outputdir][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto dir = temp_dir / "postpro";

  EnsureDirectory(dir, comm);
  // A second call on the now-existing directory must also succeed.
  EnsureDirectory(dir, comm);

  CHECK(fs::is_directory(dir));
}

// Node-local-filesystem simulation: with MPI_COMM_SELF every rank runs the full
// creation path independently on its own (distinct) directory, as it would when
// each node owns its filesystem. This exercises the per-node creation path on
// every rank -- which the MPI_COMM_WORLD cases above cannot on a single physical
// machine, where MPI_COMM_TYPE_SHARED groups all ranks into one node.
TEST_CASE_METHOD(palace::test::PerRankTempDir,
                 "EnsureDirectory creates a per-rank directory on MPI_COMM_SELF",
                 "[outputdir][Serial][Parallel]")
{
  auto dir = temp_dir / "postpro";
  REQUIRE(!fs::exists(dir));

  EnsureDirectory(dir, MPI_COMM_SELF);

  CHECK(fs::is_directory(dir));
  {
    std::ofstream f(dir / "written.txt");
    f << "ok";
  }
  CHECK(fs::is_regular_file(dir / "written.txt"));
}
