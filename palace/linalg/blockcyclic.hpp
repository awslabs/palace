// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_BLOCK_CYCLIC_HPP
#define PALACE_LINALG_BLOCK_CYCLIC_HPP

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <vector>
#include <mpi.h>

namespace palace
{

//
// The 2D block cyclic distribution of a dense n x n matrix (ScaLAPACK's): blocks of
// mb x mb dealt cyclically over a grid of nprow x npcol processes, the first ranks of a
// communicator, rank p at (p / npcol, p % npcol), which is the layout of MUMPS's
// distributed Schur complement. The local part of a process is column-major, with leading
// dimension LeadingDim().
//
class BlockCyclicGrid
{
public:
  // A grid for an n x n matrix of nproc processes (at most those of comm), as square as
  // their number allows. By default the grid is of all processes, unless the matrix is so
  // small that some would hold less than kMinLocal x kMinLocal entries, which only adds
  // latency. Without ScaLAPACK, the grid is the first process alone.
  BlockCyclicGrid(MPI_Comm comm, int n, int nproc = 0, int mb = 64);

  static constexpr int kMinLocal = 512;

  MPI_Comm GetComm() const { return comm; }
  int Size() const { return n; }
  int BlockSize() const { return mb; }
  int NumProcessRows() const { return nprow; }
  int NumProcessCols() const { return npcol; }

  // This process: whether it is in the grid, its place there, and its local part.
  bool InGrid() const { return myrow >= 0; }
  int MyRow() const { return myrow; }
  int MyCol() const { return mycol; }
  int LocalRows() const { return mloc; }
  int LocalCols() const { return nloc; }
  int LeadingDim() const { return std::max(1, mloc); }
  std::size_t LocalSize() const
  {
    return InGrid() ? static_cast<std::size_t>(LeadingDim()) * std::max(1, nloc) : 0;
  }

  // The global row and column of a local row and column of this process.
  int GlobalRow(int l) const { return GlobalIndex(l, mb, myrow, nprow); }
  int GlobalCol(int l) const { return GlobalIndex(l, mb, mycol, npcol); }

  // In one dimension of n indices in blocks of b over np processes: the number held by
  // process p (ScaLAPACK's NUMROC), and the global index of its local index l.
  static int LocalCount(int n, int b, int p, int np);
  static int GlobalIndex(int l, int b, int p, int np)
  {
    return ((l / b) * np + p) * b + l % b;
  }

  // The rank of the process at (pr, pc).
  int Rank(int pr, int pc) const { return pr * npcol + pc; }

  // The local parts of an n x n matrix on rank 0 (column-major, or row-major), which it
  // sends to each process of the grid in turn. Collective.
  std::vector<std::complex<double>>
  ScatterFromRoot(const std::vector<std::complex<double>> &full,
                  bool row_major = false) const;

private:
  MPI_Comm comm;
  int n, mb, nprow = 1, npcol = 1, myrow = -1, mycol = -1, mloc = 0, nloc = 0;
};

//
// The LU factorization with partial pivoting of a dense complex matrix distributed on a
// BlockCyclicGrid, in place, and solves with it: ScaLAPACK's, or LAPACK's on a grid of one
// process. Collective on the communicator of the grid.
//
class BlockCyclicLU
{
public:
  explicit BlockCyclicLU(const BlockCyclicGrid &grid);
  ~BlockCyclicLU();
  BlockCyclicLU(const BlockCyclicLU &) = delete;
  BlockCyclicLU &operator=(const BlockCyclicLU &) = delete;

  // Factor A (the local part of this process), which then holds the factors: it must
  // outlive the solves.
  void Factor(std::vector<std::complex<double>> &A);

  // Solve for nb right-hand sides on rank 0 (n x nb, column-major), in place.
  void Solve(std::vector<std::complex<double>> &b, int nb) const;

private:
  const BlockCyclicGrid &grid;
  int system = -1, ctxt = -1;  // BLACS system handle and grid context
  std::array<int, 9> desc{};   // ScaLAPACK descriptor of the factors
  std::vector<int> ipiv;
  const std::vector<std::complex<double>> *factors = nullptr;
};

}  // namespace palace

#endif  // PALACE_LINALG_BLOCK_CYCLIC_HPP
