// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SUBSTRUCTURE_MODEL_HPP
#define PALACE_FEM_SUBSTRUCTURE_MODEL_HPP

#include <cstdint>
#include <fstream>
#include <functional>
#include <string>
#include <vector>
#include <mpi.h>

namespace palace::substructuring
{

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
                     long long size, const std::function<void(std::ofstream &)> &write,
                     const std::function<void()> &others = {});
};

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
