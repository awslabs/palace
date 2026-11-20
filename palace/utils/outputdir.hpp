// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_OUTPUTDIR_HPP
#define PALACE_UTILS_OUTPUTDIR_HPP

#include <deque>
#include <exception>
#include <fstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>
#include <fmt/format.h>
#include "communication.hpp"
#include "filesystem.hpp"
#include "iodata.hpp"
#include "timer.hpp"

namespace palace
{

namespace internal
{

inline void AbortOnNodeFilesystemError(MPI_Comm comm, const std::string &operation,
                                       const std::string &local_error)
{
  int error_rank = local_error.empty() ? Mpi::Size(comm) : Mpi::Rank(comm);
  Mpi::GlobalMin(1, &error_rank, comm);
  if (error_rank == Mpi::Size(comm))
  {
    return;
  }

  std::string error = Mpi::Rank(comm) == error_rank ? local_error : std::string();
  int error_size = static_cast<int>(error.size());
  Mpi::Broadcast(1, &error_size, error_rank, comm);
  error.resize(error_size);
  if (error_size > 0)
  {
    Mpi::Broadcast(error_size, error.data(), error_rank, comm);
  }
  MFEM_ABORT(operation << " failed on MPI rank " << error_rank << ": " << error);
}

}  // namespace internal

// Apply a filesystem operation on the global root and propagate any error collectively so
// other ranks do not continue into a mismatched MPI control path.
template <typename Func>
inline void ApplyOnRootFilesystem(MPI_Comm comm, const std::string &operation, Func &&func)
{
  std::string local_error;
  if (Mpi::Root(comm))
  {
    try
    {
      func();
    }
    catch (const std::exception &e)
    {
      local_error = e.what();
    }
    catch (...)
    {
      local_error = "Unknown filesystem error";
    }
  }
  internal::AbortOnNodeFilesystemError(comm, operation, local_error);
}

// Apply an operation once to each node's view of the filesystem. The global root runs
// first so that shared filesystems are updated before the remaining node roots inspect
// them. Filesystem exceptions are reported only after every rank has completed the same
// MPI collectives, avoiding partial-rank failures and deadlocks.
template <typename Func>
inline void ApplyOnEachNodeFilesystem(MPI_Comm comm, const std::string &operation,
                                      Func &&func)
{
  std::string local_error;
  auto apply_here = [&]()
  {
    try
    {
      func();
    }
    catch (const std::exception &e)
    {
      local_error = e.what();
    }
    catch (...)
    {
      local_error = "Unknown filesystem error";
    }
  };

  if (Mpi::Root(comm))
  {
    apply_here();
  }
  Mpi::Barrier(comm);

  MPI_Comm node_comm = MPI_COMM_NULL;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                      &node_comm);
  if (Mpi::Root(node_comm) && !Mpi::Root(comm))
  {
    apply_here();
  }
  MPI_Comm_free(&node_comm);

  internal::AbortOnNodeFilesystemError(comm, operation, local_error);
}

// Apply a non-idempotent operation exactly once to each distinct filesystem view. The
// global root attempts the claim first; later node roots sharing its NFS/SMB mount observe
// the claim and skip the operation, while roots with node-local storage each acquire their
// own claim. Directory creation is a portable atomic create-if-absent operation, avoiding
// platform-specific file-lock APIs. The caller must first create the claim's parent
// directory on every filesystem view.
template <typename Func>
inline void ApplyOnceOnEachNodeFilesystem(MPI_Comm comm, const fs::path &claim_path,
                                          const std::string &operation,
                                          Func &&operation_func)
{
  ApplyOnEachNodeFilesystem(comm, "Removing stale claims for " + operation,
                            [&]() { fs::remove_all(claim_path); });

  ApplyOnEachNodeFilesystem(comm, operation,
                            [&]()
                            {
                              const bool acquired_claim = fs::create_directory(claim_path);
                              if (acquired_claim)
                              {
                                operation_func();
                              }
                            });

  ApplyOnEachNodeFilesystem(comm, "Cleaning up " + operation,
                            [&]() { fs::remove_all(claim_path); });
}

// Ensure an output directory exists on every node's filesystem. This is
// node-local-filesystem aware: on a shared filesystem only the global root creates the
// directory (all other creations hit an already-existing directory harmlessly), while on
// a node-local filesystem each node's root process creates the node's own copy. The code
// never classifies the filesystem topology; it only guarantees the invariant "the
// directory exists here."
inline void EnsureDirectory(const fs::path &dir, MPI_Comm comm)
{
  BlockTimer bt(Timer::IO);
  ApplyOnEachNodeFilesystem(
      comm, fmt::format("Creating output directory \"{}\"", dir.string()),
      [&]()
      {
        fs::create_directories(dir);
        if (!fs::is_directory(dir))
        {
          throw fs::filesystem_error("Path is not a directory", dir,
                                     std::make_error_code(std::errc::not_a_directory));
        }
      });
}

// Remove a previous run's output root (a directory, or a symlink standing in for one that
// was left by an earlier save) on every node's filesystem, so fresh output does not mix
// with stale artifacts. This is the cleanup counterpart of EnsureDirectory and mirrors
// its node-local awareness: on a node-local (non-shared) filesystem each node holds its
// own copy of the previous run's output, so removal runs once per node rather than only
// on the global root.
inline void RemovePreviousOutput(const fs::path &dir, MPI_Comm comm)
{
  BlockTimer bt(Timer::IO);
  // std::filesystem::remove_all removes a directory symlink itself without following it.
  ApplyOnEachNodeFilesystem(comm,
                            fmt::format("Removing previous output \"{}\"", dir.string()),
                            [&]() { fs::remove_all(dir); });
}

// Replace a root-written output file without following a symlink left by SaveIteration.
// The postprocessing directory is Palace-owned, so any existing file, directory, or
// symlink at the path can be removed. Open and write failures are propagated collectively
// to keep all ranks on the same MPI control path.
//
// With compress = true the file is written through mfem::ofgzstream, which gzip-compresses
// the bytes when MFEM is built with zlib (MFEM_USE_ZLIB) and writes them verbatim
// otherwise.
template <typename Func>
inline void WriteRootOutputFile(const fs::path &path, MPI_Comm comm, Func &&write_func,
                                bool compress = false)
{
  BlockTimer bt(Timer::IO);
  ApplyOnRootFilesystem(
      comm, fmt::format("Writing output file \"{}\"", path.string()),
      [&]()
      {
        fs::remove_all(path);
        if (compress)
        {
          mfem::ofgzstream stream(path.string(), true);
          write_func(stream);
          return;
        }
        std::ofstream stream(path);
        if (!stream.is_open())
        {
          throw fs::filesystem_error("Could not open output file", path,
                                     std::make_error_code(std::errc::io_error));
        }
        write_func(stream);
        stream.close();
        if (stream.fail())
        {
          throw fs::filesystem_error("Could not write output file", path,
                                     std::make_error_code(std::errc::io_error));
        }
      });
}

// Write a file made of the parts of all processes, part index by part index and process by
// process within each, without gathering it: the root writes, as for WriteRootOutputFile,
// and asks the other processes for their parts, a few ahead. Collective; each part must be
// smaller than 2 GiB.
inline void WriteRootOutputFileParts(const fs::path &path, MPI_Comm comm,
                                     const std::vector<std::string> &parts)
{
  BlockTimer bt(Timer::IO);
  const int P = Mpi::Size(comm), S = static_cast<int>(parts.size());
  std::vector<long long> size(S), sizes(Mpi::Root(comm) ? S * P : 0);
  for (int s = 0; s < S; s++)
  {
    size[s] = static_cast<long long>(parts[s].size());
  }
  MPI_Gather(size.data(), S, MPI_LONG_LONG, sizes.data(), S, MPI_LONG_LONG, 0, comm);
  MPI_Comm wcomm;
  MPI_Comm_dup(comm, &wcomm);
  std::string error;
  if (Mpi::Root(comm))
  {
    std::error_code ec;
    fs::remove_all(path, ec);
    std::ofstream fo(path, std::ios::binary);
    std::deque<std::pair<std::string, MPI_Request>> pending;  // at most 16 parts, 256 MB
    std::size_t bytes = 0;
    int next = 0;  // the next part to ask for, as s * P + p
    auto Ask = [&]()
    {
      for (; next < S * P && (pending.empty() || (pending.size() < 16 &&
                                                  bytes + sizes[(next % P) * S + next / P] <
                                                      (std::size_t(1) << 28)));
           next++)
      {
        const int s = next / P, p = next % P, n = static_cast<int>(sizes[p * S + s]);
        if (p > 0 && n > 0)
        {
          auto &[buf, req] = pending.emplace_back(std::string(n, '\0'), MPI_REQUEST_NULL);
          const char go = 0;
          MPI_Send(&go, 1, MPI_BYTE, p, s, wcomm);
          MPI_Irecv(buf.data(), n, MPI_BYTE, p, s, wcomm, &req);
          bytes += n;
        }
      }
    };
    for (int s = 0; s < S; s++)
    {
      for (int p = 0; p < P; p++)
      {
        if (sizes[p * S + s] > 0 && p == 0)
        {
          fo.write(parts[s].data(), parts[s].size());
        }
        else if (sizes[p * S + s] > 0)
        {
          Ask();
          auto &[buf, req] = pending.front();
          MPI_Wait(&req, MPI_STATUS_IGNORE);
          fo.write(buf.data(), buf.size());
          bytes -= buf.size();
          pending.pop_front();
        }
      }
    }
    fo.close();
    error = fo.fail() ? "Could not write output file" : "";
  }
  else
  {
    for (int s = 0; s < S; s++)
    {
      if (size[s] > 0)
      {
        char go;
        MPI_Recv(&go, 1, MPI_BYTE, 0, s, wcomm, MPI_STATUS_IGNORE);
        MPI_Send(parts[s].data(), static_cast<int>(size[s]), MPI_BYTE, 0, s, wcomm);
      }
    }
  }
  MPI_Comm_free(&wcomm);
  internal::AbortOnNodeFilesystemError(
      comm, fmt::format("Writing output file \"{}\"", path.string()), error);
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
  auto output_path = fs::path(output_str);
  EnsureDirectory(output_path, comm);

  // Fail fast if the output directory is not writable, with a clear message at setup rather
  // than an obscure failure deep in the solve. Checked on every node's filesystem, since a
  // node-local (non-shared) filesystem may be writable on one node but not another. Each
  // node root writes and removes a rank-unique probe file so concurrent checks on a shared
  // filesystem do not race on the same path.
  ApplyOnEachNodeFilesystem(
      comm, fmt::format("Verifying output directory \"{}\" is writable", output_str),
      [&]()
      {
        auto tmp = output_path / fmt::format("tmp_test_file_{}.txt", Mpi::Rank(comm));
        {
          std::ofstream f(tmp);
          f << "Test Print";
        }
        if (!fs::is_regular_file(tmp))
        {
          throw fs::filesystem_error("Could not create test file in output folder", tmp,
                                     std::make_error_code(std::errc::permission_denied));
        }
        fs::remove(tmp);
      });
}

}  // namespace palace
#endif  // PALACE_UTILS_OUTPUTDIR_HPP
