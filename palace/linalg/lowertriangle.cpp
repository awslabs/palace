// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lowertriangle.hpp"

#include <algorithm>
#include <utility>

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

}  // namespace palace
