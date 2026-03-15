// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>
#include <catch2/catch_test_macros.hpp>
#include "drivers/basesolver.hpp"
#include "fixtures.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"

using namespace palace;

namespace
{

// Helper to create a file with given content.
void CreateFile(const fs::path &path, const std::string &content = "test")
{
  std::ofstream f(path);
  f << content;
}

// Helper to read file content.
std::string ReadFile(const fs::path &path)
{
  std::ifstream f(path);
  return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

}  // namespace

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveIteration moves files and leaves symlinks", "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // Create test files in the output directory.
  CreateFile(temp_dir / "palace.json", R"({"key": "value"})");
  CreateFile(temp_dir / "domain-E.csv", "f,E\n1.0,2.0\n");
  CreateFile(temp_dir / "port-S.csv", "f,S11\n1.0,0.5\n");

  SaveIteration(comm, temp_dir, 1, 2);

  auto iter_dir = temp_dir / "iteration01";

  SECTION("Iteration subfolder is created")
  {
    CHECK(fs::is_directory(iter_dir));
  }

  SECTION("Regular files are moved to iteration subfolder")
  {
    CHECK(fs::is_regular_file(iter_dir / "domain-E.csv"));
    CHECK(ReadFile(iter_dir / "domain-E.csv") == "f,E\n1.0,2.0\n");
    CHECK(fs::is_regular_file(iter_dir / "port-S.csv"));
    CHECK(ReadFile(iter_dir / "port-S.csv") == "f,S11\n1.0,0.5\n");
  }

  SECTION("Symlinks are left behind for moved files")
  {
    CHECK(fs::is_symlink(temp_dir / "domain-E.csv"));
    CHECK(fs::is_symlink(temp_dir / "port-S.csv"));
    // Symlinks should resolve to the moved files.
    CHECK(ReadFile(temp_dir / "domain-E.csv") == "f,E\n1.0,2.0\n");
    CHECK(ReadFile(temp_dir / "port-S.csv") == "f,S11\n1.0,0.5\n");
  }

  SECTION("Symlinks are relative")
  {
    auto target = fs::read_symlink(temp_dir / "domain-E.csv");
    CHECK(target.is_relative());
    CHECK(target == fs::path("iteration01") / "domain-E.csv");
  }

  SECTION("palace.json is copied, not moved or symlinked")
  {
    CHECK(fs::is_regular_file(temp_dir / "palace.json"));
    CHECK(!fs::is_symlink(temp_dir / "palace.json"));
    CHECK(fs::is_regular_file(iter_dir / "palace.json"));
    CHECK(ReadFile(temp_dir / "palace.json") == R"({"key": "value"})");
    CHECK(ReadFile(iter_dir / "palace.json") == R"({"key": "value"})");
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "SaveIteration handles directories",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  CreateFile(temp_dir / "palace.json", "{}");
  fs::create_directories(temp_dir / "paraview" / "driven");
  CreateFile(temp_dir / "paraview" / "driven" / "fields.vtu", "<VTK/>");

  SaveIteration(comm, temp_dir, 1, 1);

  auto iter_dir = temp_dir / "iteration1";

  SECTION("Directory is moved to iteration subfolder")
  {
    CHECK(fs::is_directory(iter_dir / "paraview" / "driven"));
    CHECK(ReadFile(iter_dir / "paraview" / "driven" / "fields.vtu") == "<VTK/>");
  }

  SECTION("Symlink is left behind for directory")
  {
    CHECK(fs::is_symlink(temp_dir / "paraview"));
    CHECK(fs::is_directory(temp_dir / "paraview" / "driven"));
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveIteration handles two consecutive iterations", "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // Simulate first iteration output.
  CreateFile(temp_dir / "palace.json", "{}");
  CreateFile(temp_dir / "domain-E.csv", "iter1");

  SaveIteration(comm, temp_dir, 1, 1);

  // After first save: domain-E.csv is a symlink to iteration1/domain-E.csv.
  CHECK(fs::is_symlink(temp_dir / "domain-E.csv"));
  CHECK(ReadFile(temp_dir / "domain-E.csv") == "iter1");

  // Simulate second iteration: solver overwrites symlink with a real file.
  fs::remove(temp_dir / "domain-E.csv");
  CreateFile(temp_dir / "domain-E.csv", "iter2");
  CHECK(!fs::is_symlink(temp_dir / "domain-E.csv"));

  SaveIteration(comm, temp_dir, 2, 1);

  SECTION("First iteration data is preserved")
  {
    CHECK(ReadFile(temp_dir / "iteration1" / "domain-E.csv") == "iter1");
  }

  SECTION("Second iteration has new data")
  {
    CHECK(ReadFile(temp_dir / "iteration2" / "domain-E.csv") == "iter2");
  }

  SECTION("Symlink now points to second iteration")
  {
    CHECK(fs::is_symlink(temp_dir / "domain-E.csv"));
    CHECK(fs::read_symlink(temp_dir / "domain-E.csv") ==
          fs::path("iteration2") / "domain-E.csv");
    CHECK(ReadFile(temp_dir / "domain-E.csv") == "iter2");
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveIteration preserves old symlinks for files not reproduced by solver",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // First iteration produces two files.
  CreateFile(temp_dir / "palace.json", "{}");
  CreateFile(temp_dir / "domain-E.csv", "iter1");
  CreateFile(temp_dir / "extra.csv", "extra_data");

  SaveIteration(comm, temp_dir, 1, 1);

  // Simulate second iteration that only produces domain-E.csv (not extra.csv).
  fs::remove(temp_dir / "domain-E.csv");
  CreateFile(temp_dir / "domain-E.csv", "iter2");

  // extra.csv is still a symlink from iteration 1.
  CHECK(fs::is_symlink(temp_dir / "extra.csv"));

  SaveIteration(comm, temp_dir, 2, 1);

  SECTION("Old symlink is preserved, still accessible")
  {
    CHECK(fs::is_symlink(temp_dir / "extra.csv"));
    CHECK(ReadFile(temp_dir / "extra.csv") == "extra_data");
  }

  SECTION("New file is moved to iteration2")
  {
    CHECK(ReadFile(temp_dir / "iteration2" / "domain-E.csv") == "iter2");
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveIteration handles dirty output directory from previous run",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // Simulate a previous run that left files in the iteration folder.
  CreateFile(temp_dir / "palace.json", "{}");
  fs::create_directories(temp_dir / "iteration1" / "paraview" / "driven");
  CreateFile(temp_dir / "iteration1" / "paraview" / "driven" / "old.vtu", "old");
  CreateFile(temp_dir / "iteration1" / "domain-E.csv", "old_data");

  // New run produces fresh output.
  CreateFile(temp_dir / "domain-E.csv", "new_data");
  fs::create_directories(temp_dir / "paraview" / "driven");
  CreateFile(temp_dir / "paraview" / "driven" / "new.vtu", "new");

  // Should not throw despite destination files existing.
  SaveIteration(comm, temp_dir, 1, 1);

  SECTION("New data replaces old in iteration folder")
  {
    CHECK(ReadFile(temp_dir / "iteration1" / "domain-E.csv") == "new_data");
    CHECK(fs::exists(temp_dir / "iteration1" / "paraview" / "driven" / "new.vtu"));
    CHECK(!fs::exists(temp_dir / "iteration1" / "paraview" / "driven" / "old.vtu"));
  }

  SECTION("Symlinks point to new data")
  {
    CHECK(fs::is_symlink(temp_dir / "domain-E.csv"));
    CHECK(ReadFile(temp_dir / "domain-E.csv") == "new_data");
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "SaveIteration skips iteration subfolders",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  CreateFile(temp_dir / "palace.json", "{}");
  CreateFile(temp_dir / "domain-E.csv", "data");

  // Create a pre-existing iteration folder.
  fs::create_directories(temp_dir / "iteration1");
  CreateFile(temp_dir / "iteration1" / "old.csv", "old");

  SaveIteration(comm, temp_dir, 2, 1);

  SECTION("Pre-existing iteration folder is untouched")
  {
    CHECK(fs::is_directory(temp_dir / "iteration1"));
    CHECK(ReadFile(temp_dir / "iteration1" / "old.csv") == "old");
  }

  SECTION("New iteration folder is created")
  {
    CHECK(fs::is_directory(temp_dir / "iteration2"));
    CHECK(ReadFile(temp_dir / "iteration2" / "domain-E.csv") == "data");
  }
}
