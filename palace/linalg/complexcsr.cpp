// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "complexcsr.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <string_view>
#include <utility>
#include "linalg/rap.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

const mfem::HypreParMatrix *GetAssembled(const Operator *op)
{
  if (!op)
  {
    return nullptr;
  }
  if (const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(op))
  {
    return hA;
  }
  const auto *PtAP = dynamic_cast<const ParOperator *>(op);
  MFEM_VERIFY(PtAP, "Complex-valued sparse direct solvers require HypreParMatrix or "
                    "ParOperator real and imaginary parts!");
  return &PtAP->ParallelAssemble();
}

}  // namespace

ComplexLocalCSR::ComplexLocalCSR(const ComplexOperator &op)
{
  const auto *hAr = GetAssembled(op.Real());
  const auto *hAi = GetAssembled(op.Imag());
  MFEM_VERIFY(hAr || hAi,
              "Empty ComplexOperator for a complex-valued sparse direct solver!");
  const auto &hA = hAr ? *hAr : *hAi;
  MFEM_VERIFY(hA.GetGlobalNumRows() == hA.GetGlobalNumCols(),
              "Complex-valued sparse direct solvers require a square matrix!");
  n_glob = hA.GetGlobalNumRows();
  first_row = hA.GetRowStarts()[0];
  n_loc = hA.Height();
  MFEM_VERIFY(!hAr || !hAi ||
                  (hAr->Height() == hAi->Height() &&
                   hAr->GetRowStarts()[0] == hAi->GetRowStarts()[0]),
              "Real and imaginary parts with different row distributions!");

  // The diagonal and off-diagonal blocks of the parts are read in place, and the rows are
  // counted before they are stored, so that the matrix is the only allocation of its size.
  const std::pair<const mfem::HypreParMatrix *, std::complex<double>> parts[] = {
      {hAr, {1.0, 0.0}}, {hAi, {0.0, 1.0}}};
  for (const auto &[part, s] : parts)
  {
    if (part)
    {
      part->HostRead();
    }
  }

  // Entries of a local row sorted by global column index, with the entries of the two parts
  // in the same column added.
  std::vector<std::pair<HYPRE_BigInt, std::complex<double>>> row;
  auto MergeRow = [&](int i)
  {
    row.clear();
    for (const auto &[part, s] : parts)
    {
      if (!part)
      {
        continue;
      }
      auto *A = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(*part);
      const hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(A);
      const hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
      const HYPRE_BigInt first_col = hypre_ParCSRMatrixFirstColDiag(A);
      const HYPRE_BigInt *col_map = hypre_ParCSRMatrixColMapOffd(A);
      for (HYPRE_Int k = hypre_CSRMatrixI(diag)[i]; k < hypre_CSRMatrixI(diag)[i + 1]; k++)
      {
        row.emplace_back(first_col + hypre_CSRMatrixJ(diag)[k],
                         s * hypre_CSRMatrixData(diag)[k]);
      }
      for (HYPRE_Int k = hypre_CSRMatrixI(offd)[i]; k < hypre_CSRMatrixI(offd)[i + 1]; k++)
      {
        row.emplace_back(col_map[hypre_CSRMatrixJ(offd)[k]],
                         s * hypre_CSRMatrixData(offd)[k]);
      }
    }
    std::sort(row.begin(), row.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
    std::size_t n = 0;
    for (const auto &entry : row)
    {
      if (n > 0 && row[n - 1].first == entry.first)
      {
        row[n - 1].second += entry.second;
      }
      else
      {
        row[n++] = entry;
      }
    }
    row.resize(n);
  };
  I.assign(n_loc + 1, 0);
  for (int i = 0; i < n_loc; i++)
  {
    MergeRow(i);
    MFEM_VERIFY(row.size() <= static_cast<std::size_t>(
                                  std::numeric_limits<HYPRE_BigInt>::max() - I[i]),
                "Too many local nonzeros for the row offsets of a complex-valued sparse "
                "direct solver!");
    I[i + 1] = I[i] + static_cast<HYPRE_BigInt>(row.size());
  }
  J.resize(I[n_loc]);
  V.resize(I[n_loc]);
  for (int i = 0; i < n_loc; i++)
  {
    MergeRow(i);
    for (std::size_t k = 0; k < row.size(); k++)
    {
      J[I[i] + k] = row[k].first;
      V[I[i] + k] = row[k].second;
    }
  }
  for (const auto &[part, s] : parts)
  {
    if (part)
    {
      part->HypreRead();
    }
  }
}

bool ComplexLocalCSR::SamePattern(MPI_Comm comm, std::size_t &hash) const
{
  auto Hash = [](const std::vector<HYPRE_BigInt> &v)
  {
    return std::hash<std::string_view>{}(
        {reinterpret_cast<const char *>(v.data()), v.size() * sizeof(HYPRE_BigInt)});
  };
  const std::size_t hash_I = Hash(I), hash_J = Hash(J);
  const std::size_t previous = std::exchange(
      hash, hash_I ^ (hash_J + 0x9e3779b97f4a7c15ull + (hash_I << 6) + (hash_I >> 2)));
  int same = (hash == previous);
  Mpi::GlobalMin(1, &same, comm);
  return same;
}

bool HasZeroImag(const ComplexOperator &op)
{
  const auto *hAi = GetAssembled(op.Imag());
  if (!hAi)
  {
    return true;
  }

  // The test is exact (not a number counts as nonzero), so that the real part alone is the
  // same matrix as the operator.
  hAi->HostRead();
  auto *A = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(*hAi);
  bool nonzero = false;
  for (const hypre_CSRMatrix *block :
       {hypre_ParCSRMatrixDiag(A), hypre_ParCSRMatrixOffd(A)})
  {
    const auto *data = hypre_CSRMatrixData(block);
    nonzero = nonzero || std::any_of(data, data + hypre_CSRMatrixNumNonzeros(block),
                                     [](double v) { return v != 0.0; });
  }
  hAi->HypreRead();
  Mpi::GlobalOr(1, &nonzero, hAi->GetComm());
  return !nonzero;
}

void ReleaseParallelAssemble(const ComplexOperator &op)
{
  for (const auto *part : {op.Real(), op.Imag()})
  {
    if (const auto *PtAP = dynamic_cast<const ParOperator *>(part))
    {
      PtAP->StealParallelAssemble();
    }
  }
}

}  // namespace palace
