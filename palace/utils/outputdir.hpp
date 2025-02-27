// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_OUTPUTDIR_HPP
#define PALACE_UTILS_OUTPUTDIR_HPP

#include <fmt/format.h>
#include <fmt/os.h>
#include "communication.hpp"
#include "filesystem.hpp"
#include "iodata.hpp"
#include "timer.hpp"

namespace palace
{

inline void MakeOutputFolder(IoData &iodata, MPI_Comm &comm)
{
  BlockTimer bt(Timer::IO);
  // Validate and make folder on root
  auto root = Mpi::Root(comm);
  auto &output_str = iodata.problem.output;
  if (root)
  {
    MFEM_VERIFY(!output_str.empty(),
                fmt::format("Invalid output directory, got empty string \"\"."))
    // Remove any trailing "/" to get folder name
    if (output_str.back() == '/')
    {
      output_str.erase(output_str.end() - 1);
    }
    // Resolve canonical path (no ".", etc)
    auto output_path = fs::path(output_str);
    // Make folder if it does not exist
    if (!fs::exists(output_path))
    {
      MFEM_VERIFY(fs::create_directories(output_path),
                  fmt::format("Error std::filesystem could not create a directory at {}",
                              output_path.string()));
    }
    else
    {
      MFEM_VERIFY(fs::is_directory(output_path),
                  fmt::format("Output path already exists but is not a directory: {}",
                              output_path.string()));
      if (!fs::is_empty(output_path))
      {
        Mpi::Warning("Output folder is not empty; program will overwrite content! ({})",
                     output_path.string());
      }
    }
    // Ensure we can write to folder by making test file
    {
      fs::path tmp_ = output_path / "tmp_test_file.txt";
      auto file_buf = fmt::output_file(
          tmp_.string(), fmt::file::WRONLY | fmt::file::CREATE | fmt::file::TRUNC);
      file_buf.print("Test Print");
      file_buf.close();
      MFEM_VERIFY(
          fs::exists(tmp_) && fs::is_regular_file(tmp_),
          fmt::format("Error creating test file in output folder: {}", tmp_.string()));
      fs::remove(tmp_);
    }
    output_str = output_path.string();
  }

  // Broadcast new output_str to all ranks
  if (Mpi::Size(comm) > 1)
  {
    int str_len = static_cast<int>(output_str.size());
    if (root)
    {
      MFEM_VERIFY(output_str.size() == std::size_t(str_len),
                  "Overflow in stringbuffer size!");
    }
    Mpi::Broadcast(1, &str_len, 0, comm);
    output_str.resize(str_len);
    Mpi::Broadcast(str_len, output_str.data(), 0, comm);
  }
}

}  // namespace palace
#endif  // PALACE_UTILS_OUTPUTDIR_HPP