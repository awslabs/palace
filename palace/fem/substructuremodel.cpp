// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuremodel.hpp"

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <fmt/format.h>
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"
#include "utils/outputfile.hpp"

namespace palace::substructuring
{

std::string SaveModelPath(const IoData &iodata)
{
  const auto &sub = iodata.solver.substructuring;
  if (!sub || sub->Online() || !sub->save_model)
  {
    return {};
  }
  return (std::filesystem::path(iodata.problem.output) / "environment.model").string();
}

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
      const long long end = static_cast<long long>(f.tellg());
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
                          internal::ReplaceFile(
                              path,
                              [](std::ostream &f)
                              {
                                const std::uint32_t mv[2] = {kMagic, kVersion};
                                f.write(reinterpret_cast<const char *>(mv), sizeof(mv));
                              });
                        });
}

void ModelFile::CheckWritable(const std::string &path, MPI_Comm comm)
{
  ApplyOnRootFilesystem(
      comm, fmt::format("Writing the substructuring model \"{}\"", path),
      [&]() { internal::WriteFile(path, [](std::ostream &) {}, std::ios::app); });
}

void ModelFile::Append(const std::string &path, MPI_Comm comm, std::uint32_t tag,
                       long long size, const std::function<void(std::ostream &)> &write,
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
        long long written = 0;
        internal::WriteFile(
            path,
            [&](std::ostream &f)
            {
              f.write(reinterpret_cast<const char *>(&tag), sizeof(tag));
              f.write(reinterpret_cast<const char *>(&size), sizeof(size));
              const long long start = static_cast<long long>(f.tellp());
              write(f);
              written = static_cast<long long>(f.tellp()) - start;
            },
            std::ios::app);
        if (written != size)
        {
          throw std::runtime_error(fmt::format(
              "A record of the model holds {:d} bytes instead of {:d}", written, size));
        }
      });
}

void InterfaceRecords::Append(const std::string &path, MPI_Comm comm) const
{
  MFEM_VERIFY(signatures.size() == static_cast<std::size_t>(n) * width,
              "Invalid interface records of a substructuring model!");
  const int head[2] = {n, width};
  ModelFile::Append(
      path, comm, kSignaturesTag,
      static_cast<long long>(sizeof(head) + sizeof(double) * signatures.size()),
      [&](std::ostream &f)
      {
        f.write(reinterpret_cast<const char *>(head), sizeof(head));
        f.write(reinterpret_cast<const char *>(signatures.data()),
                static_cast<std::streamsize>(sizeof(double) * signatures.size()));
      });
  ModelFile::Append(path, comm, kLengthTag, sizeof(double), [&](std::ostream &f)
                    { f.write(reinterpret_cast<const char *>(&length), sizeof(double)); });
}

bool InterfaceRecords::Read(const std::string &path, const ModelFile &file, MPI_Comm comm)
{
  const auto *s = file.Find(kSignaturesTag), *l = file.Find(kLengthTag);
  int head[2] = {0, 0};
  if (!s || !l || s->size < static_cast<long long>(sizeof(head)) ||
      l->size != static_cast<long long>(sizeof(double)))
  {
    return false;
  }
  const std::vector<char> p = ModelFile::ReadPayload(path, *s, comm),
                          q = ModelFile::ReadPayload(path, *l, comm);
  if (p.size() < sizeof(head) || q.size() != sizeof(double))
  {
    return false;
  }
  std::memcpy(head, p.data(), sizeof(head));
  // The sizes, in floating point, against the record before anything is allocated.
  if (head[0] < 0 || head[1] <= 0 ||
      static_cast<double>(head[0]) * head[1] * sizeof(double) !=
          static_cast<double>(p.size() - sizeof(head)))
  {
    return false;
  }
  n = head[0];
  width = head[1];
  signatures.resize(static_cast<std::size_t>(n) * width);
  std::memcpy(signatures.data(), p.data() + sizeof(head),
              sizeof(double) * signatures.size());
  std::memcpy(&length, q.data(), sizeof(double));
  return true;
}

double InterfaceRecords::ReadLength(const std::string &path, MPI_Comm comm)
{
  const auto file = ModelFile::Read(path, comm);
  const auto *l = file.Find(kLengthTag);
  double length = 0.0;
  if (l && l->size == static_cast<long long>(sizeof(double)))
  {
    const std::vector<char> q = ModelFile::ReadPayload(path, *l, comm);
    if (q.size() == sizeof(double))
    {
      std::memcpy(&length, q.data(), sizeof(double));
    }
  }
  return length;
}

}  // namespace palace::substructuring
