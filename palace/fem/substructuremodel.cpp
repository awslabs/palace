// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuremodel.hpp"

#include <algorithm>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <system_error>
#include <fmt/format.h>
#include "utils/communication.hpp"
#include "utils/outputdir.hpp"

namespace palace::substructuring
{

namespace
{

// The tag and the size of a record.
constexpr long long kRecordHead = sizeof(std::uint32_t) + sizeof(long long);

}  // namespace

ModelFile ModelFile::Read(const std::string &path, MPI_Comm comm)
{
  int head[2] = {-1, 0};        // version, incomplete
  std::vector<long long> flat;  // tag, position and size of each record
  if (Mpi::Root(comm))
  {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (f)
    {
      const long long end = f.tellg();
      f.seekg(0);
      std::uint32_t mv[2] = {0, 0};
      head[0] = 0;
      if (f.read(reinterpret_cast<char *>(mv), sizeof(mv)) && mv[0] == kMagic)
      {
        head[0] = static_cast<int>(mv[1]);
        long long pos = sizeof(mv);
        while (pos < end)
        {
          std::uint32_t tag = 0;
          long long size = -1;
          f.seekg(pos);
          f.read(reinterpret_cast<char *>(&tag), sizeof(tag));
          f.read(reinterpret_cast<char *>(&size), sizeof(size));
          if (!f || size < 0 || size > end - pos - kRecordHead)
          {
            head[1] = 1;
            break;
          }
          flat.insert(flat.end(), {static_cast<long long>(tag), pos + kRecordHead, size});
          pos += kRecordHead + size;
        }
      }
    }
  }
  Mpi::Broadcast(2, head, 0, comm);
  int n = static_cast<int>(flat.size());
  Mpi::Broadcast(1, &n, 0, comm);
  flat.resize(n);
  Mpi::Broadcast(n, flat.data(), 0, comm);
  ModelFile file;
  file.version = head[0];
  file.incomplete = head[1];
  for (int i = 0; i < n; i += 3)
  {
    file.records.push_back({static_cast<std::uint32_t>(flat[i]), flat[i + 1], flat[i + 2]});
  }
  return file;
}

const ModelFile::Record *ModelFile::Find(std::uint32_t tag) const
{
  for (auto r = records.rbegin(); r != records.rend(); ++r)
  {
    if (r->tag == tag)
    {
      return &*r;
    }
  }
  return nullptr;
}

std::vector<char> ModelFile::ReadPayload(const std::string &path, const Record &r,
                                         MPI_Comm comm)
{
  std::vector<char> payload;
  long long size = 0;
  if (Mpi::Root(comm))
  {
    std::ifstream f(path, std::ios::binary);
    payload.resize(r.size);
    f.seekg(r.pos);
    f.read(payload.data(), r.size);
    size = f ? r.size : 0;
  }
  Mpi::Broadcast(1, &size, 0, comm);
  payload.resize(size);
  for (long long i = 0; i < size; i += std::numeric_limits<int>::max())
  {
    const int m =
        static_cast<int>(std::min<long long>(std::numeric_limits<int>::max(), size - i));
    MPI_Bcast(payload.data() + i, m, MPI_BYTE, 0, comm);
  }
  return payload;
}

void ModelFile::Create(const std::string &path, MPI_Comm comm)
{
  ApplyOnRootFilesystem(comm, fmt::format("Writing the substructuring model \"{}\"", path),
                        [&]()
                        {
                          std::ofstream f(path, std::ios::binary | std::ios::trunc);
                          const std::uint32_t mv[2] = {kMagic, kVersion};
                          f.write(reinterpret_cast<const char *>(mv), sizeof(mv));
                          f.close();
                          if (f.fail())
                          {
                            throw std::filesystem::filesystem_error(
                                "Could not write the model file", path,
                                std::make_error_code(std::errc::io_error));
                          }
                        });
}

void ModelFile::CheckWritable(const std::string &path, MPI_Comm comm)
{
  ApplyOnRootFilesystem(comm, fmt::format("Writing the substructuring model \"{}\"", path),
                        [&]()
                        {
                          std::ofstream f(path, std::ios::binary | std::ios::app);
                          if (!f)
                          {
                            throw std::filesystem::filesystem_error(
                                "Could not open the model file", path,
                                std::make_error_code(std::errc::io_error));
                          }
                        });
}

void ModelFile::Append(const std::string &path, MPI_Comm comm, std::uint32_t tag,
                       long long size, const std::function<void(std::ofstream &)> &write,
                       const std::function<void()> &others)
{
  if (!Mpi::Root(comm) && others)
  {
    others();
  }
  ApplyOnRootFilesystem(
      comm, fmt::format("Writing the substructuring model \"{}\"", path),
      [&]()
      {
        std::ofstream f(path, std::ios::binary | std::ios::app);
        f.write(reinterpret_cast<const char *>(&tag), sizeof(tag));
        f.write(reinterpret_cast<const char *>(&size), sizeof(size));
        const long long start = f.tellp();
        write(f);
        const long long written = f.tellp() - start;
        f.close();
        if (f.fail())
        {
          throw std::filesystem::filesystem_error(
              "Could not write the model file", path,
              std::make_error_code(std::errc::io_error));
        }
        if (written != size)
        {
          throw std::runtime_error(fmt::format(
              "A record of the model holds {:d} bytes instead of {:d}", written, size));
        }
      });
}

}  // namespace palace::substructuring
