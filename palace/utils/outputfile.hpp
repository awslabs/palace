// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_OUTPUTFILE_HPP
#define PALACE_UTILS_OUTPUTFILE_HPP

#include <cerrno>
#include <exception>
#include <fstream>
#include <ios>
#include <random>
#include <stdexcept>
#include <string>
#include <system_error>
#include <mfem/general/error.hpp>
#include "utils/filesystem.hpp"

namespace palace
{

// Every output file is written by write_func(std::ostream &) through internal::WriteFile,
// which checks the open, the write, and the close: WriteOutputFile here for files written
// by a single process, and WriteRankOutputFile and WriteRootOutputFile (outputdir.hpp) for
// files written by the processes of a collective call.

namespace internal
{

// Run func, returning the message of any exception it throws, or an empty string if it
// throws none.
template <typename Func>
inline std::string CatchError(Func &&func)
{
  try
  {
    func();
  }
  catch (const std::exception &e)
  {
    return e.what();
  }
  catch (...)
  {
    return "Unknown filesystem error";
  }
  return {};
}

// A failure to open, write, close, or replace an output file, as opposed to an error of
// the write_func writing it.
struct OutputFileError : std::runtime_error
{
  OutputFileError(const std::string &what, const fs::path &path, int err)
    : std::runtime_error(what + " \"" + path.string() + "\"" +
                         (err ? ": " + std::generic_category().message(err) : ""))
  {
  }
};

// Write a file with write_func, throwing OutputFileError if it cannot be opened, written,
// or closed. The check comes after the close, since the stream buffers data until then,
// and some failures, like a full disk or an exceeded quota on a network filesystem, are
// only reported by the close. The reason (errno) is reported when the open or the close
// fails with one. write_func runs even if the file does not open, writing to a stream that
// discards its output, so it may communicate with other processes running it. Errors name
// the file name, by default the path.
template <typename Func>
inline void WriteFile(const fs::path &path, Func &&write_func,
                      std::ios::openmode mode = std::ios::trunc, const fs::path &name = {})
{
  errno = 0;
  std::ofstream stream(path, std::ios::out | std::ios::binary | mode);
  const bool opened = stream.is_open();
  const int open_errno = errno;
  write_func(stream);
  errno = 0;
  stream.close();
  if (!opened)
  {
    throw OutputFileError("Could not open output file", name.empty() ? path : name,
                          open_errno);
  }
  if (stream.fail())
  {
    throw OutputFileError("Could not write output file", name.empty() ? path : name, errno);
  }
}

// Replace a file with one written as by WriteFile, through a temporary file renamed over
// it: readers never see a partial file, a failed write keeps the previous file, and a
// symlink at the path is replaced rather than followed. The temporary file has a random
// name, so processes replacing the same file, or a file left by a killed process, do not
// collide.
template <typename Func>
inline void ReplaceFile(const fs::path &path, Func &&write_func)
{
  auto tmp = path;
  tmp += "." + std::to_string(std::random_device{}()) + ".tmp";
  std::error_code ec;
  try
  {
    WriteFile(tmp, write_func, std::ios::trunc, path);
  }
  catch (...)
  {
    fs::remove(tmp, ec);
    throw;
  }
  fs::rename(tmp, path, ec);
  if (ec)
  {
    const int err = ec.value();
    fs::remove(tmp, ec);
    throw OutputFileError("Could not replace output file", path, err);
  }
}

// Run func, which writes the file at path, returning the message of any exception it
// throws, naming the file, or an empty string if it throws none.
template <typename Func>
inline std::string CatchFileError(const fs::path &path, Func &&func)
{
  try
  {
    func();
  }
  catch (const OutputFileError &e)
  {
    return e.what();
  }
  catch (const std::exception &e)
  {
    return "Could not write output file \"" + path.string() + "\": " + e.what();
  }
  catch (...)
  {
    return "Could not write output file \"" + path.string() + "\"";
  }
  return {};
}

}  // namespace internal

// Write a file on this process, aborting if it cannot be opened, written, or closed. The
// file is replaced as by internal::ReplaceFile, or appended to with mode = std::ios::app.
// Not collective: for files written by a single process, like the tables and metadata the
// root writes.
template <typename Func>
inline void WriteOutputFile(const fs::path &path, Func &&write_func,
                            std::ios::openmode mode = std::ios::trunc)
{
  const auto error =
      internal::CatchFileError(path,
                               [&]()
                               {
                                 if (mode & std::ios::app)
                                 {
                                   internal::WriteFile(path, write_func, mode);
                                 }
                                 else
                                 {
                                   internal::ReplaceFile(path, write_func);
                                 }
                               });
  if (!error.empty())
  {
    MFEM_ABORT(error);
  }
}

}  // namespace palace

#endif  // PALACE_UTILS_OUTPUTFILE_HPP
