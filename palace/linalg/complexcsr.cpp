// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "complexcsr.hpp"

#include <algorithm>
#include <utility>
#include "linalg/rap.hpp"

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

// Local rows of a HypreParMatrix with global column indices, on the host.
struct LocalRows
{
  hypre_CSRMatrix *csr = nullptr;
  explicit LocalRows(const mfem::HypreParMatrix *A)
  {
    if (A)
    {
      csr = hypre_MergeDiagAndOffd(
          (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(*A));
      hypre_CSRMatrixMigrate(csr, HYPRE_MEMORY_HOST);
    }
  }
  ~LocalRows()
  {
    if (csr)
    {
      hypre_CSRMatrixDestroy(csr);
    }
  }
};

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

  {
    LocalRows re(hAr), im(hAi);
    I.assign(n_loc + 1, 0);
    std::vector<std::pair<HYPRE_BigInt, std::complex<double>>> row;
    for (int i = 0; i < n_loc; i++)
    {
      row.clear();
      for (auto [part, s] : {std::pair{&re, std::complex<double>(1.0, 0.0)},
                             std::pair{&im, std::complex<double>(0.0, 1.0)}})
      {
        if (!part->csr)
        {
          continue;
        }
        const HYPRE_Int *Ip = hypre_CSRMatrixI(part->csr);
        const HYPRE_BigInt *Jp = hypre_CSRMatrixBigJ(part->csr);
        const double *Dp = hypre_CSRMatrixData(part->csr);
        for (HYPRE_Int k = Ip[i]; k < Ip[i + 1]; k++)
        {
          row.emplace_back(Jp[k], s * Dp[k]);
        }
      }
      std::sort(row.begin(), row.end(),
                [](const auto &a, const auto &b) { return a.first < b.first; });
      for (const auto &[j, v] : row)
      {
        if (static_cast<int>(J.size()) > I[i] && J.back() == j)
        {
          V.back() += v;
        }
        else
        {
          J.push_back(j);
          V.push_back(v);
        }
      }
      I[i + 1] = static_cast<int>(J.size());
    }
  }

  // The solvers copy the matrix: release the parallel assembled parts (as MfemWrapperSolver
  // does for the sparse direct solvers).
  for (const auto *part : {op.Real(), op.Imag()})
  {
    if (const auto *PtAP = dynamic_cast<const ParOperator *>(part))
    {
      PtAP->StealParallelAssemble();
    }
  }
}

}  // namespace palace
