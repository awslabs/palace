// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lowertriangle.hpp"

#include <algorithm>
#include <complex>
#include <limits>
#include <numeric>
#include <utility>
#include "linalg/blockcyclic.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// f(i, j, a) for the entries of local row i of A with global column j <= global row.
template <typename F>
void ForEachLowerEntry(const mfem::HypreParMatrix &A, int i, F &&f)
{
  auto *parcsr = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(A);
  const hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(parcsr);
  const hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(parcsr);
  const HYPRE_BigInt *cmap = hypre_ParCSRMatrixColMapOffd(parcsr);
  const HYPRE_BigInt row = hypre_ParCSRMatrixFirstRowIndex(parcsr) + i;
  const HYPRE_BigInt col0 = hypre_ParCSRMatrixFirstColDiag(parcsr);
  for (HYPRE_Int k = diag->i[i]; k < diag->i[i + 1]; k++)
  {
    const HYPRE_BigInt j = col0 + diag->j[k];
    if (j <= row)
    {
      f(row, j, diag->data[k]);
    }
  }
  for (HYPRE_Int k = offd->i ? offd->i[i] : 0; offd->i && k < offd->i[i + 1]; k++)
  {
    const HYPRE_BigInt j = cmap[offd->j[k]];
    if (j <= row)
    {
      f(row, j, offd->data[k]);
    }
  }
}

// The values on their owners of the local true DOFs of 1-based global columns: v[q] =
// value[cols[q] - 1 - tstart] of the rank owning column cols[q]. Collective.
std::vector<HYPRE_BigInt> OwnerValues(MPI_Comm comm, HYPRE_BigInt tstart,
                                      const std::vector<HYPRE_BigInt> &value,
                                      const std::vector<int> &cols)
{
  const int nranks = Mpi::Size(comm);
  const HYPRE_BigInt tend = tstart + static_cast<HYPRE_BigInt>(value.size());
  std::vector<HYPRE_BigInt> starts(nranks);
  MPI_Allgather(&tstart, 1, HYPRE_MPI_BIG_INT, starts.data(), 1, HYPRE_MPI_BIG_INT, comm);
  auto owner = [&](HYPRE_BigInt j)
  {
    return static_cast<int>(std::upper_bound(starts.begin(), starts.end(), j) -
                            starts.begin()) -
           1;
  };
  std::vector<std::vector<HYPRE_BigInt>> req(nranks);
  for (int c : cols)
  {
    const HYPRE_BigInt j = c - 1;
    if (j < tstart || j >= tend)
    {
      req[owner(j)].push_back(j);
    }
  }
  std::vector<int> scnt(nranks), sdisp(nranks, 0), rcnt(nranks), rdisp(nranks, 0);
  std::vector<HYPRE_BigInt> sbuf;
  for (int r = 0; r < nranks; r++)
  {
    std::sort(req[r].begin(), req[r].end());
    req[r].erase(std::unique(req[r].begin(), req[r].end()), req[r].end());
    scnt[r] = static_cast<int>(req[r].size());
    sdisp[r] = static_cast<int>(sbuf.size());
    sbuf.insert(sbuf.end(), req[r].begin(), req[r].end());
  }
  MPI_Alltoall(scnt.data(), 1, MPI_INT, rcnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    rdisp[r] = rdisp[r - 1] + rcnt[r - 1];
  }
  std::vector<HYPRE_BigInt> rbuf(rdisp[nranks - 1] + rcnt[nranks - 1]);
  MPI_Alltoallv(sbuf.data(), scnt.data(), sdisp.data(), HYPRE_MPI_BIG_INT, rbuf.data(),
                rcnt.data(), rdisp.data(), HYPRE_MPI_BIG_INT, comm);
  for (auto &j : rbuf)
  {
    j = value[j - tstart];
  }
  std::vector<HYPRE_BigInt> ans(sbuf.size());
  MPI_Alltoallv(rbuf.data(), rcnt.data(), rdisp.data(), HYPRE_MPI_BIG_INT, ans.data(),
                scnt.data(), sdisp.data(), HYPRE_MPI_BIG_INT, comm);
  std::vector<HYPRE_BigInt> v(cols.size());
  for (std::size_t q = 0; q < cols.size(); q++)
  {
    const HYPRE_BigInt j = cols[q] - 1;
    if (j >= tstart && j < tend)
    {
      v[q] = value[j - tstart];
    }
    else
    {
      const int r = owner(j);
      v[q] = ans[sdisp[r] +
                 (std::lower_bound(req[r].begin(), req[r].end(), j) - req[r].begin())];
    }
  }
  return v;
}

}  // namespace

void LowerTrianglePattern(const std::vector<const mfem::HypreParMatrix *> &A,
                          std::vector<int> &irn, std::vector<int> &jcn,
                          std::vector<std::vector<double>> &vals, std::vector<int> &row_ptr)
{
  const mfem::HypreParMatrix *A0 = nullptr;
  for (const auto *X : A)
  {
    if (X)
    {
      MFEM_VERIFY(!A0 || (X->Height() == A0->Height() &&
                          X->GetRowStarts()[0] == A0->GetRowStarts()[0]),
                  "Matrices of a pattern must have the same row distribution!");
      A0 = A0 ? A0 : X;
      X->HostRead();
    }
  }
  MFEM_VERIFY(A0, "LowerTrianglePattern needs a matrix!");
  MFEM_VERIFY(A0->GetGlobalNumRows() <= std::numeric_limits<int>::max(),
              "LowerTrianglePattern needs fewer than 2^31 unknowns (32-bit indices)!");
  const int n = A0->Height();
  const HYPRE_BigInt row0 = A0->GetRowStarts()[0];
  // Two passes over the rows (count, then fill) to allocate the pattern exactly.
  std::vector<std::pair<HYPRE_BigInt, int>> row;  // (column, matrix)
  std::vector<double> rv;
  auto gather = [&](int i)
  {
    row.clear();
    rv.clear();
    for (std::size_t p = 0; p < A.size(); p++)
    {
      if (A[p])
      {
        ForEachLowerEntry(*A[p], i,
                          [&](HYPRE_BigInt, HYPRE_BigInt j, double a)
                          {
                            row.emplace_back(j, static_cast<int>(p));
                            rv.push_back(a);
                          });
      }
    }
  };
  std::vector<int> perm;
  auto sorted = [&]()
  {
    perm.resize(row.size());
    for (std::size_t q = 0; q < perm.size(); q++)
    {
      perm[q] = static_cast<int>(q);
    }
    std::sort(perm.begin(), perm.end(),
              [&](int a, int b) { return row[a].first < row[b].first; });
  };
  row_ptr.assign(n + 1, 0);
  for (int i = 0; i < n; i++)
  {
    gather(i);
    sorted();
    int m = 0;
    for (std::size_t q = 0; q < perm.size(); q++)
    {
      m += (q == 0 || row[perm[q]].first != row[perm[q - 1]].first);
    }
    row_ptr[i + 1] = row_ptr[i] + m;
  }
  const std::size_t nnz = row_ptr[n];
  irn.assign(nnz, 0);
  jcn.assign(nnz, 0);
  vals.assign(A.size(), {});
  for (std::size_t p = 0; p < A.size(); p++)
  {
    if (A[p])
    {
      vals[p].assign(nnz, 0.0);
    }
  }
  for (int i = 0; i < n; i++)
  {
    gather(i);
    sorted();
    int pos = row_ptr[i] - 1;
    for (std::size_t q = 0; q < perm.size(); q++)
    {
      const auto [j, p] = row[perm[q]];
      if (q == 0 || j != row[perm[q - 1]].first)
      {
        pos++;
        irn[pos] = static_cast<int>(row0 + i + 1);
        jcn[pos] = static_cast<int>(j + 1);
      }
      vals[p][pos] += rv[perm[q]];
    }
  }
  for (const auto *X : A)
  {
    if (X)
    {
      X->HypreRead();
    }
  }
}

SideSystem::SideSystem(const std::vector<const mfem::HypreParMatrix *> &parts_in,
                       const std::vector<char> &keep)
{
  const mfem::HypreParMatrix *A0 = nullptr;
  for (const auto *X : parts_in)
  {
    A0 = A0 ? A0 : X;
  }
  MFEM_VERIFY(A0, "SideSystem needs a matrix!");
  comm = A0->GetComm();
  const int nt = A0->Height();
  MFEM_VERIFY(keep.size() == static_cast<std::size_t>(nt),
              "SideSystem needs one flag per local true DOF!");
  tstart = A0->GetRowStarts()[0];
  index.assign(nt, -1);
  for (int i = 0; i < nt; i++)
  {
    if (keep[i])
    {
      rows.push_back(i);
    }
  }
  HYPRE_BigInt nloc = static_cast<HYPRE_BigInt>(rows.size()), off = 0;
  MPI_Exscan(&nloc, &off, 1, HYPRE_MPI_BIG_INT, MPI_SUM, comm);
  off = Mpi::Root(comm) ? 0 : off;
  MPI_Allreduce(&nloc, &n_glob, 1, HYPRE_MPI_BIG_INT, MPI_SUM, comm);
  for (std::size_t q = 0; q < rows.size(); q++)
  {
    index[rows[q]] = off + static_cast<HYPRE_BigInt>(q);
  }

  // The pattern without the entries of pinned rows and columns.
  std::vector<int> t_irn, t_jcn, t_ptr;
  std::vector<std::vector<double>> vals;
  LowerTrianglePattern(parts_in, t_irn, t_jcn, vals, t_ptr);
  const auto col_index = OwnerValues(comm, tstart, index, t_jcn);
  parts.assign(parts_in.size(), {});
  for (int i = 0; i < nt; i++)
  {
    for (int q = t_ptr[i]; index[i] >= 0 && q < t_ptr[i + 1]; q++)
    {
      if (col_index[q] >= 0)
      {
        irn.push_back(static_cast<int>(index[i] + 1));
        jcn.push_back(static_cast<int>(col_index[q] + 1));
        for (std::size_t p = 0; p < parts.size(); p++)
        {
          if (!vals[p].empty())
          {
            parts[p].push_back(vals[p][q]);
          }
        }
      }
    }
  }
}

std::vector<HYPRE_BigInt> SideSystem::Indices(const std::vector<HYPRE_BigInt> &tdofs) const
{
  const HYPRE_BigInt tend = tstart + static_cast<HYPRE_BigInt>(index.size());
  std::vector<HYPRE_BigInt> out(tdofs.size(), -1);
  for (std::size_t q = 0; q < tdofs.size(); q++)
  {
    if (tdofs[q] >= tstart && tdofs[q] < tend)
    {
      out[q] = index[tdofs[q] - tstart];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, out.data(), static_cast<int>(out.size()), HYPRE_MPI_BIG_INT,
                MPI_MAX, comm);
  return out;
}

template <typename T>
std::vector<T> SideSystem::Values(const std::vector<T> &coef) const
{
  MFEM_VERIFY(coef.size() == parts.size(), "SideSystem: one coefficient per part!");
  std::vector<T> v(irn.size(), T(0.0));
  for (std::size_t p = 0; p < parts.size(); p++)
  {
    for (std::size_t q = 0; q < parts[p].size(); q++)
    {
      v[q] += coef[p] * parts[p][q];
    }
  }
  return v;
}

template std::vector<double> SideSystem::Values(const std::vector<double> &) const;
template std::vector<std::complex<double>>
SideSystem::Values(const std::vector<std::complex<double>> &) const;

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
template void ToRoot<std::complex<double>>(
    MPI_Comm, const std::vector<std::complex<double>> &, const std::vector<long long> &,
    std::size_t, const std::function<void(const std::complex<double> *, std::size_t)> &);

namespace
{

// Entries (g, h, v) to the processes holding them, owner(g, h), which add v at index(g,
// h) of y: this process gives its entries in rounds, produce(round, add) calling add(g,
// h, v) for those of a round, as many rounds as the process with the most. Collective.
template <typename T, typename Produce>
void SendEntries(MPI_Comm comm, long long rounds, Produce &&produce,
                 const std::function<int(int, int)> &owner,
                 const std::function<std::size_t(int, int)> &index, std::vector<T> &y)
{
  const int nranks = Mpi::Size(comm);
  Mpi::GlobalMax(1, &rounds, comm);
  const auto type = mpi::DataType<T>();
  std::vector<T> svals, rvals;
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
  for (long long round = 0; round < rounds; round++)
  {
    for (int q = 0; q < nranks; q++)
    {
      idx[q].clear();
      vals[q].clear();
    }
    produce(round, add);
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

}  // namespace

template <typename T>
void RedistributeTriangle(MPI_Comm comm, long long n,
                          const std::function<void(long long, long long, T *)> &read,
                          const std::vector<std::vector<std::pair<int, double>>> &to,
                          const std::function<int(int, int)> &owner,
                          const std::function<std::size_t(int, int)> &index,
                          std::vector<T> &y)
{
  const auto range = TriangleShare(n, comm);
  const long long e0 = range[0], e1 = range[1];
  const auto first = TriangleEntry(n, e0);
  long long a = first[0], b = first[1];
  std::vector<T> buf;
  SendEntries<T>(
      comm, (e1 - e0 + kTriangleChunk - 1) / kTriangleChunk,
      [&](long long round, auto &add)
      {
        const long long e = e0 + round * kTriangleChunk;
        buf.resize(
            static_cast<std::size_t>(std::max(0LL, std::min(e1, e + kTriangleChunk) - e)));
        if (!buf.empty())
        {
          read(e, static_cast<long long>(buf.size()), buf.data());
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
      },
      owner, index, y);
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

template <typename T>
std::vector<T> ShareFromGrid(const BlockCyclicGrid &grid, const std::vector<T> &local)
{
  MPI_Comm comm = grid.GetComm();
  const long long n = grid.Size();
  const int rank = Mpi::Rank(comm);
  // The first entry of the share of each process, to find the owner of an entry.
  const auto size = TriangleShareSizes(n, comm);
  std::vector<long long> start(size.size() + 1, 0);
  std::partial_sum(size.begin(), size.end(), start.begin() + 1);
  std::vector<T> y(static_cast<std::size_t>(start[rank + 1] - start[rank]), T(0.0));
  auto entry = [n](int g, int h) { return TriangleColumnStart(n, h) + (g - h); };
  // The local entries on or below the diagonal, in rounds of local columns.
  const int mloc = grid.LocalRows(), nloc = grid.LocalCols();
  const std::size_t lld = grid.LeadingDim();
  const long long cols = std::max(1LL, kTriangleChunk / std::max(1, mloc));
  SendEntries<T>(
      comm, (nloc + cols - 1) / cols,
      [&](long long round, auto &add)
      {
        for (long long jl = round * cols;
             jl < std::min<long long>(nloc, (round + 1) * cols); jl++)
        {
          const int h = grid.GlobalCol(static_cast<int>(jl));
          for (int il = 0; il < mloc; il++)
          {
            const int g = grid.GlobalRow(il);
            if (g >= h)
            {
              add(g, h, local[il + lld * static_cast<std::size_t>(jl)]);
            }
          }
        }
      },
      [&](int g, int h)
      {
        return static_cast<int>(std::upper_bound(start.begin(), start.end(), entry(g, h)) -
                                start.begin() - 1);
      },
      [&](int g, int h) { return static_cast<std::size_t>(entry(g, h) - start[rank]); }, y);
  return y;
}

template std::vector<double> ShareFromGrid(const BlockCyclicGrid &,
                                           const std::vector<double> &);
template std::vector<std::complex<double>>
ShareFromGrid(const BlockCyclicGrid &, const std::vector<std::complex<double>> &);

}  // namespace palace
