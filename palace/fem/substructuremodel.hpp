// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SUBSTRUCTURE_MODEL_HPP
#define PALACE_FEM_SUBSTRUCTURE_MODEL_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include <mpi.h>

namespace palace
{

class IoData;

}  // namespace palace

namespace palace::substructuring
{

// The model which an offline substructuring run saves, environment.model in its output
// directory, or empty: an online run, or an offline one with "SaveModel": false.
std::string SaveModelPath(const IoData &iodata);

// The file of a saved substructuring model: a magic number and the version of the
// container, then records [tag][size][payload], with a tag of four characters and the size
// of the payload in bytes. The records are appended as the model is computed. A record that
// runs past the end of the file, as one being written when a run stopped, is incomplete,
// and is dropped with the rest of the file. A later record with a tag supersedes an earlier
// one, and a reader skips the tags it does not know. Rank 0 reads and writes the file.
struct ModelFile
{
  static constexpr std::uint32_t kMagic = 0x42555350;  // "PSUB"
  static constexpr std::uint32_t kVersion = 1;

  // The tag of four characters.
  static constexpr std::uint32_t Tag(const char (&s)[5])
  {
    return static_cast<std::uint32_t>(static_cast<unsigned char>(s[0])) |
           static_cast<std::uint32_t>(static_cast<unsigned char>(s[1])) << 8 |
           static_cast<std::uint32_t>(static_cast<unsigned char>(s[2])) << 16 |
           static_cast<std::uint32_t>(static_cast<unsigned char>(s[3])) << 24;
  }

  struct Record
  {
    std::uint32_t tag = 0;
    long long pos = 0, size = 0;  // of the payload, in bytes
  };

  // The version of the container of the file (-1: the file cannot be read, 0: it is not a
  // model), whether it ends in an incomplete record, and its complete records.
  int version = -1;
  bool incomplete = false;
  std::vector<Record> records;

  // The records of a file, read on rank 0 and broadcast. Collective.
  static ModelFile Read(const std::string &path, MPI_Comm comm);

  // The last record with the tag, or null.
  const Record *Find(std::uint32_t tag) const;

  // The payload of a record, read on rank 0 and broadcast (empty if the read fails).
  // Collective.
  static std::vector<char> ReadPayload(const std::string &path, const Record &r,
                                       MPI_Comm comm);

  // A new file with only the magic number and the version. Collective.
  static void Create(const std::string &path, MPI_Comm comm);

  // Stop every rank unless rank 0 can write the file, which it opens for appending (an
  // existing file is kept). Collective.
  static void CheckWritable(const std::string &path, MPI_Comm comm);

  // Append a record of the given size: rank 0 writes its tag and size, then its payload
  // with write, while the other ranks run others (sending rank 0 their parts of the
  // payload, for example). A stream that fails writes nothing, so write still receives
  // them. Stops every rank if rank 0 cannot write the record or write does not write size
  // bytes. Collective.
  static void Append(const std::string &path, MPI_Comm comm, std::uint32_t tag,
                     long long size, const std::function<void(std::ostream &)> &write,
                     const std::function<void()> &others = {});
};

// The lower triangle of a symmetric n x n matrix by columns, as the records of a model
// hold it (column b holds the rows a >= b): its number of entries, and the first entry of
// column b.
constexpr long long TriangleSize(long long n)
{
  return n * (n + 1) / 2;
}
constexpr long long TriangleColumnStart(long long n, long long b)
{
  return b * n - b * (b - 1) / 2;
}

// The row a and column b of entry e of the triangle of an n x n matrix.
std::array<long long, 2> TriangleEntry(long long n, long long e);

// The share [e0, e1) of this process of the entries of the triangle of an n x n matrix, an
// equal one of the processes of comm, and the sizes of the shares of all of them.
std::array<long long, 2> TriangleShare(long long n, MPI_Comm comm);
std::vector<long long> TriangleShareSizes(long long n, MPI_Comm comm);

// The triangle of an n x n matrix, each process of comm holding its share of it
// (TriangleShare), added onto a target layout through a map of its DOFs: read(e, m, x)
// fills x with the m entries from e of the share of this process, at most kTriangleChunk at
// a time, and each entry S_ab adds c_ga c_hb S_ab to the entries (g, h) and, for a > b,
// (h, g) of the target, for the (g, c_ga) of to[a] and the (h, c_hb) of to[b] (with to
// empty, g = a and h = b with coefficients 1). The process owner(g, h) holds the target
// entry (g, h), at index(g, h) of its local part y. Collective.
inline constexpr long long kTriangleChunk = 1LL << 18;
template <typename T>
void RedistributeTriangle(MPI_Comm comm, long long n,
                          const std::function<void(long long, long long, T *)> &read,
                          const std::vector<std::vector<std::pair<int, double>>> &to,
                          const std::function<int(int, int)> &owner,
                          const std::function<std::size_t(int, int)> &index,
                          std::vector<T> &y);

// To rank 0: the values x of every process (size[q] of process q) in rank order, which
// write(x, n) takes n at a time at x on rank 0, in chunks of at most chunk values, so that
// it holds one chunk of another process's at a time. It receives them all even if writing
// fails, and the other processes only send theirs. Collective.
template <typename T>
void ToRoot(MPI_Comm comm, const std::vector<T> &x, const std::vector<long long> &size,
            std::size_t chunk, const std::function<void(const T *, std::size_t)> &write);

// From rank 0: the values of every process (size[q] of process q) in rank order, which
// read(e, m, x) gives on rank 0, filling x with the m values from e of all of them, in
// chunks of at most chunk values, so that it holds one chunk at a time. Returns those of
// this process. Collective.
template <typename T>
std::vector<T> FromRoot(MPI_Comm comm, const std::vector<long long> &size,
                        std::size_t chunk,
                        const std::function<void(long long, long long, T *)> &read);

// The records of the interface every substructuring model holds: the signatures of its DOFs
// (n x width, row-major, see TrueDofSignatures) and the characteristic length of the run
// that saved it (model.Lc, which an online run takes, see BaseSolver::Preprocess).
struct InterfaceRecords
{
  static constexpr std::uint32_t kSignaturesTag = ModelFile::Tag("SIGS"),
                                 kLengthTag = ModelFile::Tag("LENC");
  int n = 0, width = 0;
  std::vector<double> signatures;
  double length = 0.0;

  // Append them to a model file (see ModelFile::Append). Collective.
  void Append(const std::string &path, MPI_Comm comm) const;

  // Read them from the records of a file: false if one is missing or does not hold its
  // sizes. Collective.
  bool Read(const std::string &path, const ModelFile &file, MPI_Comm comm);

  // The characteristic length of a model file, from its length record alone (0 if it holds
  // none). Collective.
  static double ReadLength(const std::string &path, MPI_Comm comm);
};

}  // namespace palace::substructuring

#endif  // PALACE_FEM_SUBSTRUCTURE_MODEL_HPP
