// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuremodel.hpp"

#include <algorithm>
#include <complex>
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

std::array<long long, 2> TriangleEntry(long long n, long long e)
{
  long long lo = 0, hi = std::max(n, 1LL);
  while (hi - lo > 1)
  {
    const long long mid = (lo + hi) / 2;
    (TriangleColumnStart(n, mid) <= e ? lo : hi) = mid;
  }
  return {lo + (e - TriangleColumnStart(n, lo)), lo};
}

std::array<long long, 2> TriangleShare(long long n, MPI_Comm comm)
{
  const long long total = TriangleSize(n), p = Mpi::Size(comm), r = Mpi::Rank(comm);
  return {total * r / p, total * (r + 1) / p};
}

std::vector<long long> TriangleShareSizes(long long n, MPI_Comm comm)
{
  const long long total = TriangleSize(n), p = Mpi::Size(comm);
  std::vector<long long> size(p);
  for (long long r = 0; r < p; r++)
  {
    size[r] = total * (r + 1) / p - total * r / p;
  }
  return size;
}

template <typename T>
void ToRoot(MPI_Comm comm, const std::vector<T> &x, const std::vector<long long> &size,
            std::size_t chunk, const std::function<void(const T *, std::size_t)> &write)
{
  const auto type = mpi::DataType<T>();
  const std::size_t per = std::max<std::size_t>(1, chunk);
  if (!Mpi::Root(comm))
  {
    for (std::size_t i = 0; i < x.size(); i += per)
    {
      MPI_Send(x.data() + i, static_cast<int>(std::min(per, x.size() - i)), type, 0, 0,
               comm);
    }
    return;
  }
  // The chunks of another process as it sends them, until it has sent its part.
  std::vector<T> buf;
  for (std::size_t q = 0; q < size.size(); q++)
  {
    for (std::size_t i = 0; i < static_cast<std::size_t>(size[q]);)
    {
      const std::size_t n = std::min(per, static_cast<std::size_t>(size[q]) - i);
      const T *y = x.data() + i;
      if (q > 0)
      {
        buf.resize(n);
        MPI_Recv(buf.data(), static_cast<int>(n), type, static_cast<int>(q), 0, comm,
                 MPI_STATUS_IGNORE);
        y = buf.data();
      }
      write(y, n);
      i += n;
    }
  }
}

template void ToRoot<double>(MPI_Comm, const std::vector<double> &,
                             const std::vector<long long> &, std::size_t,
                             const std::function<void(const double *, std::size_t)> &);

template <typename T>
std::vector<T> FromRoot(MPI_Comm comm, const std::vector<long long> &size,
                        std::size_t chunk,
                        const std::function<void(long long, long long, T *)> &read)
{
  const auto type = mpi::DataType<T>();
  const std::size_t per = std::max<std::size_t>(1, chunk);
  std::vector<T> x(static_cast<std::size_t>(size[Mpi::Rank(comm)]));
  if (!Mpi::Root(comm))
  {
    for (std::size_t i = 0; i < x.size(); i += per)
    {
      MPI_Recv(x.data() + i, static_cast<int>(std::min(per, x.size() - i)), type, 0, 0,
               comm, MPI_STATUS_IGNORE);
    }
    return x;
  }
  std::vector<T> buf;
  long long e = 0;
  for (std::size_t q = 0; q < size.size(); q++)
  {
    for (std::size_t i = 0; i < static_cast<std::size_t>(size[q]);)
    {
      const std::size_t n = std::min(per, static_cast<std::size_t>(size[q]) - i);
      if (q == 0)
      {
        read(e, static_cast<long long>(n), x.data() + i);
      }
      else
      {
        buf.resize(n);
        read(e, static_cast<long long>(n), buf.data());
        MPI_Send(buf.data(), static_cast<int>(n), type, static_cast<int>(q), 0, comm);
      }
      e += static_cast<long long>(n);
      i += n;
    }
  }
  return x;
}

template std::vector<std::complex<double>> FromRoot<std::complex<double>>(
    MPI_Comm, const std::vector<long long> &, std::size_t,
    const std::function<void(long long, long long, std::complex<double> *)> &);

template <typename T>
void RedistributeTriangle(MPI_Comm comm, long long n,
                          const std::function<void(long long, long long, T *)> &read,
                          const std::vector<std::vector<std::pair<int, double>>> &to,
                          const std::function<int(int, int)> &owner,
                          const std::function<std::size_t(int, int)> &index,
                          std::vector<T> &y)
{
  const int nranks = Mpi::Size(comm);
  const auto [e0, e1] = TriangleShare(n, comm);
  long long rounds = (e1 - e0 + kTriangleChunk - 1) / kTriangleChunk;
  Mpi::GlobalMax(1, &rounds, comm);

  auto [a, b] = TriangleEntry(n, e0);

  const auto type = mpi::DataType<T>();
  std::vector<T> buf, svals, rvals;
  std::vector<std::vector<int>> idx(nranks);
  std::vector<std::vector<T>> vals(nranks);
  std::vector<int> sidx, ridx, scnt(nranks), rcnt(nranks), sdsp(nranks), rdsp(nranks);
  auto add = [&](int g, int h, T v)
  {
    const int q = owner(g, h);
    idx[q].push_back(g);
    idx[q].push_back(h);
    vals[q].push_back(v);
  };
  for (long long round = 0, e = e0; round < rounds; round++)
  {
    const long long e_end = std::min(e1, e + kTriangleChunk);
    buf.resize(static_cast<std::size_t>(e_end - e));
    if (!buf.empty())
    {
      read(e, e_end - e, buf.data());
    }
    for (int q = 0; q < nranks; q++)
    {
      idx[q].clear();
      vals[q].clear();
    }
    for (const T v : buf)
    {
      if (to.empty())
      {
        add(static_cast<int>(a), static_cast<int>(b), v);
        if (a != b)
        {
          add(static_cast<int>(b), static_cast<int>(a), v);
        }
      }
      else
      {
        for (const auto &[g, cg] : to[a])
        {
          for (const auto &[h, ch] : to[b])
          {
            add(g, h, (cg * ch) * v);
            if (a != b)
            {
              add(h, g, (cg * ch) * v);
            }
          }
        }
      }
      if (++a == n)
      {
        a = ++b;
      }
    }
    e = e_end;

    // The entries to the processes holding them, which add them up.
    for (int q = 0; q < nranks; q++)
    {
      scnt[q] = static_cast<int>(vals[q].size());
    }
    MPI_Alltoall(scnt.data(), 1, MPI_INT, rcnt.data(), 1, MPI_INT, comm);
    sdsp[0] = rdsp[0] = 0;
    for (int q = 1; q < nranks; q++)
    {
      sdsp[q] = sdsp[q - 1] + scnt[q - 1];
      rdsp[q] = rdsp[q - 1] + rcnt[q - 1];
    }
    svals.clear();
    sidx.clear();
    for (int q = 0; q < nranks; q++)
    {
      svals.insert(svals.end(), vals[q].begin(), vals[q].end());
      sidx.insert(sidx.end(), idx[q].begin(), idx[q].end());
    }
    const int nrecv = rdsp[nranks - 1] + rcnt[nranks - 1];
    rvals.resize(nrecv);
    ridx.resize(2 * static_cast<std::size_t>(nrecv));
    MPI_Alltoallv(svals.data(), scnt.data(), sdsp.data(), type, rvals.data(), rcnt.data(),
                  rdsp.data(), type, comm);
    for (int q = 0; q < nranks; q++)
    {
      scnt[q] *= 2;
      sdsp[q] *= 2;
      rcnt[q] *= 2;
      rdsp[q] *= 2;
    }
    MPI_Alltoallv(sidx.data(), scnt.data(), sdsp.data(), MPI_INT, ridx.data(), rcnt.data(),
                  rdsp.data(), MPI_INT, comm);
    for (int k = 0; k < nrecv; k++)
    {
      y[index(ridx[2 * k], ridx[2 * k + 1])] += rvals[k];
    }
  }
}

template void RedistributeTriangle<double>(
    MPI_Comm, long long, const std::function<void(long long, long long, double *)> &,
    const std::vector<std::vector<std::pair<int, double>>> &,
    const std::function<int(int, int)> &, const std::function<std::size_t(int, int)> &,
    std::vector<double> &);
template void RedistributeTriangle<std::complex<double>>(
    MPI_Comm, long long,
    const std::function<void(long long, long long, std::complex<double> *)> &,
    const std::vector<std::vector<std::pair<int, double>>> &,
    const std::function<int(int, int)> &, const std::function<std::size_t(int, int)> &,
    std::vector<std::complex<double>> &);

}  // namespace palace::substructuring
