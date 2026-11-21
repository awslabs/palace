// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fixtures.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/outputdir.hpp"

using namespace palace;

// The core-logic cases are tagged [Serial][Parallel] so they run in both the
// single-rank and the multi-rank sweeps. In the parallel sweep they additionally
// verify that the collective (barrier / split / free) is balanced across ranks
// and does not deadlock, and that a directory created on the (shared) filesystem
// is visible to every rank afterwards. Reduce preconditions before fatal assertions
// so every rank takes the same path into the following collectives.

TEST_CASE_METHOD(palace::test::SharedTempDir, "EnsureDirectory creates a new directory",
                 "[outputdir][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto dir = temp_dir / "postpro";
  bool directory_absent = !fs::exists(dir);
  Mpi::GlobalAnd(1, &directory_absent, comm);
  REQUIRE(directory_absent);

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
  bool directory_exists = fs::is_directory(dir);
  Mpi::GlobalAnd(1, &directory_exists, comm);
  REQUIRE(directory_exists);

  // Regression guard: EnsureDirectory must NOT treat an already-existing
  // directory as an error.
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
  bool directory_absent = !fs::exists(dir);
  Mpi::GlobalAnd(1, &directory_absent, comm);
  REQUIRE(directory_absent);

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

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RemovePreviousOutput removes an existing directory and its contents",
                 "[outputdir][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto dir = temp_dir / "postpro";
  if (Mpi::Root(comm))
  {
    fs::create_directories(dir);
    std::ofstream f(dir / "stale.txt");
    f << "stale artifact";
  }
  Mpi::Barrier(comm);
  bool directory_exists = fs::is_directory(dir);
  Mpi::GlobalAnd(1, &directory_exists, comm);
  REQUIRE(directory_exists);

  RemovePreviousOutput(dir, comm);

  // The directory and everything under it must be gone on every rank.
  CHECK(!fs::exists(dir));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RemovePreviousOutput removes a symlink without following it",
                 "[outputdir][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto target = temp_dir / "target";
  auto link = temp_dir / "postpro";
  if (Mpi::Root(comm))
  {
    fs::create_directories(target);
    std::ofstream f(target / "keep.txt");
    f << "keep me";
    fs::create_directory_symlink(target, link);
  }
  Mpi::Barrier(comm);
  bool link_exists = fs::is_symlink(link);
  Mpi::GlobalAnd(1, &link_exists, comm);
  REQUIRE(link_exists);

  RemovePreviousOutput(link, comm);

  // The symlink itself must be gone, but its target's contents must survive: a symlink
  // standing in for the output directory must be unlinked, not followed and deleted.
  CHECK(!fs::exists(link));
  CHECK(fs::is_regular_file(target / "keep.txt"));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RemovePreviousOutput is a no-op on a non-existent path",
                 "[outputdir][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto dir = temp_dir / "does_not_exist";
  bool directory_absent = !fs::exists(dir);
  Mpi::GlobalAnd(1, &directory_absent, comm);
  REQUIRE(directory_absent);

  // Removing a path that was never created must not error or abort.
  RemovePreviousOutput(dir, comm);

  CHECK(!fs::exists(dir));
}

namespace
{

std::string ReadFile(const fs::path &path)
{
  std::ifstream f(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

// A device on which every write fails, as on a full disk, while the open succeeds. Only
// write it in place (internal::WriteFile, or appending): a writer which replaces the file
// would rename a temporary file over the device when run as root.
const fs::path full_device = "/dev/full";

}  // namespace

TEST_CASE("Output files report failed writes", "[outputdir][Serial]")
{
  using Catch::Matchers::ContainsSubstring;
  if (!fs::is_character_file(full_device))
  {
    SKIP("No " << full_device << " to test a failed write");
  }
  // The data are buffered, so the write fails when the file is closed, with the reason.
  CHECK_THROWS_WITH(internal::WriteFile(full_device, [](std::ostream &os) { os << "x"; }),
                    ContainsSubstring("Could not write output file") &&
                        ContainsSubstring("No space left on device"));
  // A write larger than the buffer fails before the close.
  const std::string data(1 << 20, 'x');
  CHECK_THROWS_WITH(internal::WriteFile(full_device, [&](std::ostream &os) { os << data; }),
                    ContainsSubstring("Could not write output file"));
  CHECK_THROWS_WITH(WriteOutputFile(
                        full_device, [](std::ostream &os) { os << "x"; }, std::ios::app),
                    ContainsSubstring("No space left on device"));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "WriteOutputFile replaces or appends to a file, and aborts on a failure",
                 "[outputdir][Serial]")
{
  using Catch::Matchers::ContainsSubstring;
  auto path = temp_dir / "file.txt", target = temp_dir / "target.txt";
  WriteOutputFile(path, [](std::ostream &os) { os << "a\n"; });
  WriteOutputFile(path, [](std::ostream &os) { os << "b\n"; }, std::ios::app);
  CHECK(ReadFile(path) == "a\nb\n");

  // A failed write keeps the previous file, and leaves no temporary file.
  CHECK_THROWS_WITH(WriteOutputFile(path,
                                    [](std::ostream &os)
                                    {
                                      os << "c\n";
                                      throw std::runtime_error("write_func failed");
                                    }),
                    ContainsSubstring("Could not write output file \"" + path.string() +
                                      "\": write_func failed"));
  CHECK(ReadFile(path) == "a\nb\n");
  std::vector<std::string> files;
  for (const auto &entry : fs::directory_iterator(temp_dir))
  {
    files.push_back(entry.path().filename().string());
  }
  CHECK(files == std::vector<std::string>{"file.txt"});

  // A symlink, like those SaveIteration leaves, is replaced rather than followed.
  WriteOutputFile(target, [](std::ostream &os) { os << "keep\n"; });
  fs::remove(path);
  fs::create_symlink(target, path);
  WriteOutputFile(path, [](std::ostream &os) { os << "d\n"; });
  CHECK(!fs::is_symlink(path));
  CHECK(ReadFile(path) == "d\n");
  CHECK(ReadFile(target) == "keep\n");

  CHECK_THROWS_WITH(WriteOutputFile(temp_dir / "missing" / "file.txt",
                                    [](std::ostream &os) { os << "x"; }),
                    ContainsSubstring("Could not open output file \"" +
                                      (temp_dir / "missing" / "file.txt").string() +
                                      "\"") &&
                        ContainsSubstring("No such file or directory"));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "WriteRootOutputFile writes compressed files, and aborts on every rank",
                 "[outputdir][Serial][Parallel]")
{
  using Catch::Matchers::ContainsSubstring;
  MPI_Comm comm = MPI_COMM_WORLD;
  auto path = temp_dir / "file.gz";
  // More than the zlib stream buffers, and a flush within, which ends a gzip member.
  const std::string data(3 << 20, 'z');
  WriteRootOutputFile(
      path, comm, [&](std::ostream &os) { os << data << std::flush << "end"; }, true);
  {
    mfem::ifgzstream f(path.string());
    std::string read{std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
    CHECK(read == data + "end");
  }
#if defined(MFEM_USE_ZLIB)
  CHECK(fs::file_size(path) < data.size());
#endif

  CHECK_THROWS_WITH(WriteRootOutputFile(temp_dir / "missing" / "file.txt", comm,
                                        [](std::ostream &os) { os << "x"; }),
                    ContainsSubstring("Writing output files failed on MPI rank 0") &&
                        ContainsSubstring("Could not open output file"));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "WriteRankOutputFile failures abort on every rank when one rank fails",
                 "[outputdir][Serial][Parallel]")
{
  using Catch::Matchers::ContainsSubstring;
  MPI_Comm comm = MPI_COMM_WORLD;
  const int rank = Mpi::Rank(comm), last = Mpi::Size(comm) - 1;
  auto path = temp_dir / ("file_" + std::to_string(rank) + ".txt");
  std::string error;
  WriteRankOutputFile(path, [&](std::ostream &os) { os << rank; }, error);
  CHECK(error.empty());
  AbortOnOutputFileErrors(comm, error);
  CHECK(ReadFile(path) == std::to_string(rank));

  // Errors of write_func itself are not caught.
  CHECK_THROWS_AS(
      WriteRankOutputFile(
          path, [](std::ostream &) { throw std::runtime_error("write_func"); }, error),
      std::runtime_error);
  CHECK(error.empty());

  // The last rank cannot open its file. write_func still runs on it, so a write_func which
  // communicates does not deadlock the other ranks.
  auto bad_path = rank == last ? temp_dir / "missing" / "file.txt" : path;
  WriteRankOutputFile(
      bad_path,
      [&](std::ostream &os)
      {
        int one = 1;
        Mpi::GlobalSum(1, &one, comm);
        os << one;
      },
      error);
  CHECK(error.empty() == (rank != last));
  CHECK_THROWS_WITH(AbortOnOutputFileErrors(comm, error),
                    ContainsSubstring("Writing output files failed on MPI rank " +
                                      std::to_string(last)) &&
                        ContainsSubstring("Could not open output file"));
}
