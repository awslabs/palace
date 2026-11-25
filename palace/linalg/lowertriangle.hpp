// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_LOWER_TRIANGLE_HPP
#define PALACE_LINALG_LOWER_TRIANGLE_HPP

#include <array>
#include <complex>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class BlockCyclicGrid;

// Lower-triangle COO pattern (1-based global indices, one entry per position, columns
// ascending in each row) of the local rows of the union of the given symmetric matrices
// (same row distribution; null ones skipped). The entries of each matrix in pattern order
// are returned in vals (empty for a null one), and the first pattern entry of each local
// row in row_ptr (n_loc + 1).
void LowerTrianglePattern(const std::vector<const mfem::HypreParMatrix *> &A,
                          std::vector<int> &irn, std::vector<int> &jcn,
                          std::vector<std::vector<double>> &vals,
                          std::vector<int> &row_ptr);

// The values on their owners of the local true DOFs of 1-based global columns: v[q] =
// value[cols[q] - 1 - tstart] of the rank owning column cols[q]. Collective.
std::vector<HYPRE_BigInt> OwnerValues(MPI_Comm comm, HYPRE_BigInt tstart,
                                      const std::vector<HYPRE_BigInt> &value,
                                      const std::vector<int> &cols);

// val += a X on a pattern from LowerTrianglePattern, for a matrix X with entries inside the
// pattern (or zero).
void AddToPattern(const mfem::HypreParMatrix &X, double a, const std::vector<int> &row_ptr,
                  const std::vector<int> &jcn, std::vector<double> &val);
void AddToPattern(const mfem::HypreParMatrix &X, std::complex<double> a,
                  const std::vector<int> &row_ptr, const std::vector<int> &jcn,
                  std::vector<std::complex<double>> &val);

// The rows and columns of a symmetric operator of local true DOFs it keeps (the others
// pinned), numbered by rank (the system): the lower-triangle pattern of the union of its
// parts (same row distribution; null ones skipped) on the kept DOFs, as COO entries in
// system indices (1-based), with the entries of each part on it, so that the system's
// values are a combination of its parts. Collective.
class SideSystem
{
public:
  SideSystem(const std::vector<const mfem::HypreParMatrix *> &parts,
             const std::vector<char> &keep);

  // The system indices (0-based) of global true DOFs (replicated, in their order; -1 for
  // one not kept), each rank mapping those it owns. Collective.
  std::vector<HYPRE_BigInt> Indices(const std::vector<HYPRE_BigInt> &tdofs) const;

  // Σ_p coef[p] parts[p] on the pattern (null parts skipped).
  template <typename T>
  std::vector<T> Values(const std::vector<T> &coef) const;

  // val += a X on the entries of the system, for a matrix X of the rows of the parts with
  // entries inside the pattern of their union (or zero): a part with new values.
  void AddPart(const mfem::HypreParMatrix &X, double a, std::vector<double> &val) const;
  void AddPart(const mfem::HypreParMatrix &X, std::complex<double> a,
               std::vector<std::complex<double>> &val) const;

  // The local true DOF of each local row of the system, and the global size of the system.
  std::vector<int> rows;
  HYPRE_BigInt n_glob = 0;
  // The pattern, its entries in system indices, and the entries of each part on it.
  std::vector<int> irn, jcn;
  std::vector<std::vector<double>> parts;

private:
  MPI_Comm comm;
  HYPRE_BigInt tstart = 0;
  std::vector<HYPRE_BigInt> index;  // the system index of each local true DOF (-1: pinned)
  // The pattern of the union of the parts on the local true DOFs (see LowerTrianglePattern)
  // and the place in it of each entry of the system.
  std::vector<int> row_ptr, tdof_jcn, kept;

  template <typename T>
  void AddPartT(const mfem::HypreParMatrix &X, T a, std::vector<T> &val) const;
};

// The lower triangle of a symmetric n x n matrix by columns, as the records of a
// substructuring model hold it (column b holds the rows a >= b): its number of entries, and
// the first entry of column b.
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

// The most entries of a triangle that a process handles at a time when it moves them.
inline constexpr long long kTriangleChunk = 1LL << 18;

// The triangle of an n x n matrix, each process of comm holding its share of it
// (TriangleShare), added onto a target layout through a map of its DOFs: read(e, m, x)
// fills x with the m entries from e of the share of this process, at most kTriangleChunk at
// a time, and each entry S_ab adds c_ga c_hb S_ab to the entries (g, h) and, for a > b,
// (h, g) of the target, for the (g, c_ga) of to[a] and the (h, c_hb) of to[b] (with to
// empty, g = a and h = b with coefficients 1). The process owner(g, h) holds the target
// entry (g, h), at index(g, h) of its local part y. Collective.
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

// The share of this process (TriangleShare) of the lower triangle of a symmetric matrix on
// a grid, from the local parts of the processes in it (the entries above the diagonal
// ignored). Collective.
template <typename T>
std::vector<T> ShareFromGrid(const BlockCyclicGrid &grid, const std::vector<T> &local);

}  // namespace palace

#endif  // PALACE_LINALG_LOWER_TRIANGLE_HPP
