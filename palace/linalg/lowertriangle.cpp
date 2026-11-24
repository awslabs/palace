// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lowertriangle.hpp"

#include <algorithm>
#include <limits>
#include <utility>
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

namespace
{

template <typename T>
void AddToPatternT(const mfem::HypreParMatrix &X, T a, const std::vector<int> &row_ptr,
                   const std::vector<int> &jcn, std::vector<T> &val)
{
  X.HostRead();
  for (int i = 0; i < X.Height(); i++)
  {
    const auto first = jcn.begin() + row_ptr[i], last = jcn.begin() + row_ptr[i + 1];
    ForEachLowerEntry(X, i,
                      [&](HYPRE_BigInt, HYPRE_BigInt j, double x)
                      {
                        const auto it =
                            std::lower_bound(first, last, static_cast<int>(j + 1));
                        if (it != last && *it == j + 1)
                        {
                          val[it - jcn.begin()] += a * x;
                        }
                        else
                        {
                          MFEM_VERIFY(x == 0.0, "Entry outside of the sparsity pattern!");
                        }
                      });
  }
  X.HypreRead();
}

}  // namespace

void AddToPattern(const mfem::HypreParMatrix &X, double a, const std::vector<int> &row_ptr,
                  const std::vector<int> &jcn, std::vector<double> &val)
{
  AddToPatternT(X, a, row_ptr, jcn, val);
}

void AddToPattern(const mfem::HypreParMatrix &X, std::complex<double> a,
                  const std::vector<int> &row_ptr, const std::vector<int> &jcn,
                  std::vector<std::complex<double>> &val)
{
  AddToPatternT(X, a, row_ptr, jcn, val);
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
  row_ptr.assign(nt + 1, 0);
  for (int i = 0; i < nt; i++)
  {
    for (int q = t_ptr[i]; index[i] >= 0 && q < t_ptr[i + 1]; q++)
    {
      if (col_index[q] >= 0)
      {
        tdof_jcn.push_back(t_jcn[q]);
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
    row_ptr[i + 1] = static_cast<int>(tdof_jcn.size());
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

}  // namespace palace
