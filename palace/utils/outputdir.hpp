// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_OUTPUTDIR_HPP
#define PALACE_UTILS_OUTPUTDIR_HPP

#include <system_error>
#include <fmt/format.h>
#include "communication.hpp"
#include "filesystem.hpp"
#include "iodata.hpp"
#include "timer.hpp"

namespace palace
{

// Ensure an output directory exists on every node's filesystem. This is
// node-local-filesystem aware: on a shared filesystem only the global root creates the
// directory (all other creations hit an already-existing directory harmlessly), while on
// a node-local filesystem each node's root process creates the node's own copy. The code
// never classifies the filesystem topology; it only guarantees the invariant "the
// directory exists here."
inline void EnsureDirectory(const fs::path &dir, MPI_Comm comm)
{
  BlockTimer bt(Timer::IO);
  // Global root creates the directory. Use the non-throwing overload and ignore both the
  // error code and the bool return: neither is a reliable success signal (the bool is
  // false for an already-existing directory, and a benign EEXIST race may populate the
  // error code). A genuine write failure surfaces at the first real write.
  if (Mpi::Root(comm))
  {
    std::error_code ec;
    fs::create_directories(dir, ec);
  }
  Mpi::Barrier(comm);

  // Split off a node-local communicator so one process per node can ensure the directory
  // exists on that node's filesystem.
  MPI_Comm node_comm = MPI_COMM_NULL;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                      &node_comm);

  // On a node-local filesystem this is the real create. On a shared filesystem it
  // is a harmless EEXIST no-op that, on typical network filesystems, also refreshes
  // this node's metadata cache so ranks here are less likely to hit a stale
  // "not found" before their first write (not guaranteed by POSIX).
  if (Mpi::Root(node_comm) && !Mpi::Root(comm))
  {
    std::error_code ec;
    fs::create_directories(dir, ec);
  }

  MPI_Comm_free(&node_comm);
  Mpi::Barrier(comm);
}

inline void MakeOutputFolder(IoData &iodata, MPI_Comm &comm)
{
  BlockTimer bt(Timer::IO);
  // Validate output directory name and warn about overwrites on root.
  auto root = Mpi::Root(comm);
  auto &output_str = iodata.problem.output;
  if (root)
  {
    MFEM_VERIFY(!output_str.empty(),
                fmt::format("Invalid output directory, got empty string \"\"."))
    // Remove any trailing "/" to get folder name.
    if (output_str.back() == '/')
    {
      output_str.erase(output_str.end() - 1);
    }
    auto output_path = fs::path(output_str);
    // Warn if the folder already exists and is not empty; content will be overwritten.
    if (fs::exists(output_path))
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
    output_str = output_path.string();
  }

  // Broadcast new output_str to all ranks.
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

  // Create the output directory (node-local aware).
  EnsureDirectory(fs::path(output_str), comm);
}

}  // namespace palace
#endif  // PALACE_UTILS_OUTPUTDIR_HPP
