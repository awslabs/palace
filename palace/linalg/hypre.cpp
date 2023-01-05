// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hypre.hpp"

#include <map>

namespace mfem
{

mfem::HypreParMatrix *
HypreParMatrixFromBlocks(mfem::Array2D<const mfem::HypreParMatrix *> &blocks,
                         mfem::Array2D<double> *blockCoeff)
{
  mfem::Array2D<mfem::HypreParMatrix *> blocks_without_const(blocks.NumRows(),
                                                             blocks.NumCols());
  for (int i = 0; i < blocks.NumRows(); i++)
  {
    for (int j = 0; j < blocks.NumCols(); j++)
    {
      blocks_without_const(i, j) = const_cast<mfem::HypreParMatrix *>(blocks(i, j));
    }
  }
  return HypreParMatrixFromBlocks(blocks_without_const, blockCoeff);
}

}  // namespace mfem

namespace palace::hypre
{

void hypreParCSREliminateRowsCols(hypre_ParCSRMatrix *A, const mfem::Array<int> &rows_cols,
                                  hypre::DiagonalPolicy diag_policy, HYPRE_Complex diag,
                                  bool ignore_rows)
{
  hypre_error_flag = 0;

  hypre_CSRMatrix *A_diag = hypre_ParCSRMatrixDiag(A);
  HYPRE_Real *A_diag_a = hypre_CSRMatrixData(A_diag);
  HYPRE_Int *A_diag_i = hypre_CSRMatrixI(A_diag);
  HYPRE_Int *A_diag_j = hypre_CSRMatrixJ(A_diag);
  HYPRE_Int ncols_A_diag = hypre_CSRMatrixNumCols(A_diag);

  hypre_CSRMatrix *A_offd = hypre_ParCSRMatrixOffd(A);
  HYPRE_Real *A_offd_a = hypre_CSRMatrixData(A_offd);
  HYPRE_Int *A_offd_i = hypre_CSRMatrixI(A_offd);
  HYPRE_Int *A_offd_j = hypre_CSRMatrixJ(A_offd);
  HYPRE_Int ncols_A_offd = hypre_CSRMatrixNumCols(A_offd);

  HYPRE_BigInt *col_map_offd_A = hypre_ParCSRMatrixColMapOffd(A);
  HYPRE_Int *marker_offd = nullptr;

  HYPRE_BigInt first_row = hypre_ParCSRMatrixFirstRowIndex(A);
  HYPRE_Int nrows_local = hypre_CSRMatrixNumRows(A_diag);

  HYPRE_Int i, j, k, nnz_diag, nnz_offd, A_diag_i_i, A_offd_i_i;

  // Get markers for columns of the diagonal and off-diagonal matrix to eliminate
  // (from mfem::internal::hypre_ParCSRMatrixEliminateAAe).
  HYPRE_Int *eliminate_diag_col, *eliminate_offd_col;
  {
    hypre_ParCSRCommHandle *comm_handle;
    hypre_ParCSRCommPkg *comm_pkg;
    HYPRE_Int num_sends, *int_buf_data;
    HYPRE_Int index, start;

    eliminate_diag_col = mfem_hypre_CTAlloc_host(HYPRE_Int, ncols_A_diag);
    eliminate_offd_col = mfem_hypre_CTAlloc_host(HYPRE_Int, ncols_A_offd);

    // Make sure A has a communication package.
    comm_pkg = hypre_ParCSRMatrixCommPkg(A);
    if (!comm_pkg)
    {
      hypre_MatvecCommPkgCreate(A);
      comm_pkg = hypre_ParCSRMatrixCommPkg(A);
    }

    // Which of the local rows are to be eliminated.
    for (i = 0; i < ncols_A_diag; i++)
    {
      eliminate_diag_col[i] = 0;
    }
    for (i = 0; i < rows_cols.Size(); i++)
    {
      eliminate_diag_col[rows_cols[i]] = 1;
    }

    // Use a Matvec communication pattern to find (in eliminate_col) which of the local offd
    // columns are to be eliminated.
    num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
    int_buf_data = mfem_hypre_CTAlloc_host(
        HYPRE_Int, hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
    index = 0;
    for (i = 0; i < num_sends; i++)
    {
      start = hypre_ParCSRCommPkgSendMapStart(comm_pkg, i);
      for (j = start; j < hypre_ParCSRCommPkgSendMapStart(comm_pkg, i + 1); j++)
      {
        k = hypre_ParCSRCommPkgSendMapElmt(comm_pkg, j);
        int_buf_data[index++] = eliminate_diag_col[k];
      }
    }
    comm_handle =
        hypre_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, eliminate_offd_col);

    // Finish the communication.
    hypre_ParCSRCommHandleDestroy(comm_handle);

    mfem_hypre_TFree_host(int_buf_data);
  }

  marker_offd = hypre_CTAlloc(HYPRE_Int, ncols_A_offd, HYPRE_MEMORY_HOST);

  nnz_diag = nnz_offd = A_diag_i_i = A_offd_i_i = 0;
  for (i = 0; i < nrows_local; i++)
  {
    // Drop eliminated entries in the diagonal block.
    for (j = A_diag_i_i; j < A_diag_i[i + 1]; j++)
    {
      HYPRE_Int col = A_diag_j[j];
      HYPRE_Complex val = A_diag_a[j];
      if ((!ignore_rows && eliminate_diag_col[i]) || eliminate_diag_col[col])
      {
        // Always keep the diagonal entry (even if it is 0).
        if (!ignore_rows && i == col)
        {
          if (diag_policy == DiagonalPolicy::USER)
          {
            val = diag;
          }
          else if (diag_policy == DiagonalPolicy::ONE)
          {
            val = 1.0;
          }
          else if (diag_policy == DiagonalPolicy::ZERO)
          {
            val = 0.0;
          }
          // else (diag_policy == DiagonalPolicy::KEEP)
        }
        else
        {
          continue;
        }
      }
      A_diag_j[nnz_diag] = col;
      A_diag_a[nnz_diag] = val;
      nnz_diag++;
    }

    // Drop eliminated entries in the off-diagonal block.
    for (j = A_offd_i_i; j < A_offd_i[i + 1]; j++)
    {
      HYPRE_Int col = A_offd_j[j];
      HYPRE_Complex val = A_offd_a[j];
      if ((!ignore_rows && eliminate_diag_col[i]) || eliminate_offd_col[col])
      {
        // In normal cases: diagonal entry should not appear in A_offd (but this can still
        // be possible).
        if (!ignore_rows && i + first_row == col_map_offd_A[col])
        {
          if (diag_policy == DiagonalPolicy::USER)
          {
            val = diag;
          }
          else if (diag_policy == DiagonalPolicy::ONE)
          {
            val = 1.0;
          }
          else if (diag_policy == DiagonalPolicy::ZERO)
          {
            val = 0.0;
          }
          // else (diag_policy == DiagonalPolicy::KEEP)
        }
        else
        {
          continue;
        }
      }
      if (marker_offd[col] == 0)
      {
        marker_offd[col] = 1;
      }
      A_offd_j[nnz_offd] = col;
      A_offd_a[nnz_offd] = val;
      nnz_offd++;
    }
    A_diag_i_i = A_diag_i[i + 1];
    A_offd_i_i = A_offd_i[i + 1];
    A_diag_i[i + 1] = nnz_diag;
    A_offd_i[i + 1] = nnz_offd;
  }

  mfem_hypre_TFree_host(eliminate_offd_col);
  mfem_hypre_TFree_host(eliminate_diag_col);

  hypre_CSRMatrixNumNonzeros(A_diag) = nnz_diag;
  hypre_CSRMatrixNumNonzeros(A_offd) = nnz_offd;
  hypre_ParCSRMatrixSetNumNonzeros(A);
  hypre_ParCSRMatrixDNumNonzeros(A) = (HYPRE_Real)hypre_ParCSRMatrixNumNonzeros(A);

  for (i = 0, k = 0; i < ncols_A_offd; i++)
  {
    if (marker_offd[i])
    {
      col_map_offd_A[k] = col_map_offd_A[i];
      marker_offd[i] = k++;
    }
  }
  hypre_CSRMatrixNumCols(A_offd) = k;  // ncols_A_offd = k
  for (i = 0; i < nnz_offd; i++)
  {
    A_offd_j[i] = marker_offd[A_offd_j[i]];
  }

  hypre_TFree(marker_offd, HYPRE_MEMORY_HOST);

  if (hypre_ParCSRMatrixCommPkg(A))
  {
    hypre_MatvecCommPkgDestroy(hypre_ParCSRMatrixCommPkg(A));
  }
  hypre_MatvecCommPkgCreate(A);

  MFEM_VERIFY(!hypre_error_flag,
              "HYPRE error encountered: Error code = " << hypre_error_flag << "!");
}

void hypreParCSREliminateRowsColsv2(hypre_ParCSRMatrix *A,
                                    const mfem::Array<int> &rows_cols,
                                    hypre::DiagonalPolicy diag_policy, HYPRE_Complex diag,
                                    bool ignore_rows)
{
  hypre_error_flag = 0;

  hypre_CSRMatrix *A_diag = hypre_ParCSRMatrixDiag(A);
  HYPRE_Int ncols_A_diag = hypre_CSRMatrixNumRows(A_diag);

  hypre_CSRMatrix *A_offd = hypre_ParCSRMatrixOffd(A);
  HYPRE_Int ncols_A_offd = hypre_CSRMatrixNumCols(A_offd);

  const auto n_ess_dofs = rows_cols.Size();
  const auto ess_dofs_d =
      rows_cols.GetMemory().Read(mfem::GetHypreMemoryClass(), n_ess_dofs);

  // Start communication to figure out which columns need to be eliminated in the
  // off-diagonal block.
  hypre_ParCSRCommHandle *comm_handle;
  HYPRE_Int *int_buf_data, *eliminate_row, *eliminate_col;
  {
    eliminate_row = mfem_hypre_CTAlloc(HYPRE_Int, ncols_A_diag);
    eliminate_col = mfem_hypre_CTAlloc(HYPRE_Int, ncols_A_offd);

    // Make sure A has a communication package.
    hypre_ParCSRCommPkg *comm_pkg = hypre_ParCSRMatrixCommPkg(A);
    if (!comm_pkg)
    {
      hypre_MatvecCommPkgCreate(A);
      comm_pkg = hypre_ParCSRMatrixCommPkg(A);
    }

    // Which of the local rows are to be eliminated?
    MFEM_HYPRE_FORALL(i, ncols_A_diag, { eliminate_row[i] = 0; });
    MFEM_HYPRE_FORALL(i, n_ess_dofs, { eliminate_row[ess_dofs_d[i]] = 1; });

    // Use a matvec communication pattern to find (in eliminate_col) which of the local offd
    // columns are to be eliminated.
    HYPRE_Int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
    HYPRE_Int int_buf_sz = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
    int_buf_data = mfem_hypre_CTAlloc(HYPRE_Int, int_buf_sz);

    HYPRE_Int *send_map_elmts;
#if defined(HYPRE_USING_GPU)
    hypre_ParCSRCommPkgCopySendMapElmtsToDevice(comm_pkg);
    send_map_elmts = hypre_ParCSRCommPkgDeviceSendMapElmts(comm_pkg);
#else
    send_map_elmts = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
#endif
    MFEM_HYPRE_FORALL(i, int_buf_sz, {
      auto k = send_map_elmts[i];
      int_buf_data[i] = eliminate_row[k];
    });

#if defined(HYPRE_USING_GPU)
    // Try to use device-aware MPI for the communication if available.
    comm_handle =
        hypre_ParCSRCommHandleCreate_v2(11, comm_pkg, HYPRE_MEMORY_DEVICE, int_buf_data,
                                        HYPRE_MEMORY_DEVICE, eliminate_col);
#else
    comm_handle = hypre_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, eliminate_col);
#endif
  }

  // Eliminate rows and columns in the diagonal block.
  if (ignore_rows)
  {
    const auto nrows_A_diag = hypre_CSRMatrixNumRows(A_diag);
    const auto I = A_diag->i;
    const auto J = A_diag->j;
    auto data = A_diag->data;
    MFEM_HYPRE_FORALL(i, nrows_A_diag, {
      for (auto j = I[i]; j < I[i + 1]; j++)
      {
        data[j] *= 1 - eliminate_row[J[j]];
      }
    });
  }
  else
  {
    const auto I = A_diag->i;
    const auto J = A_diag->j;
    auto data = A_diag->data;
    MFEM_HYPRE_FORALL(i, n_ess_dofs, {
      const auto idof = ess_dofs_d[i];
      for (auto j = I[idof]; j < I[idof + 1]; j++)
      {
        const auto jdof = J[j];
        if (jdof == idof)
        {
          if (diag_policy == DiagonalPolicy::USER)
          {
            data[j] = diag;
          }
          else if (diag_policy == DiagonalPolicy::ONE)
          {
            data[j] = 1.0;
          }
          else if (diag_policy == DiagonalPolicy::ZERO)
          {
            data[j] = 0.0;
          }
          // else (diag_policy == DiagonalPolicy::KEEP)
        }
        else
        {
          data[j] = 0.0;
          for (auto k = I[jdof]; k < I[jdof + 1]; k++)
          {
            if (J[k] == idof)
            {
              data[k] = 0.0;
              break;
            }
          }
        }
      }
    });
  }

  // Eliminate rows in the off-diagonal block.
  if (!ignore_rows)
  {
    const auto I = A_offd->i;
    auto data = A_offd->data;
    MFEM_HYPRE_FORALL(i, n_ess_dofs, {
      const auto idof = ess_dofs_d[i];
      for (auto j = I[idof]; j < I[idof + 1]; j++)
      {
        data[j] = 0.0;
      }
    });
  }

  // Wait for MPI communication to finish.
  hypre_ParCSRCommHandleDestroy(comm_handle);
  mfem_hypre_TFree(int_buf_data);
  mfem_hypre_TFree(eliminate_row);

  // Eliminate columns in the off-diagonal block.
  {
    const auto nrows_A_offd = hypre_CSRMatrixNumRows(A_offd);
    const auto I = A_offd->i;
    const auto J = A_offd->j;
    auto data = A_offd->data;
    MFEM_HYPRE_FORALL(i, nrows_A_offd, {
      for (auto j = I[i]; j < I[i + 1]; j++)
      {
        data[j] *= 1 - eliminate_col[J[j]];
      }
    });
  }

  mfem_hypre_TFree(eliminate_col);

  MFEM_VERIFY(!hypre_error_flag,
              "HYPRE error encountered: Error code = " << hypre_error_flag << "!");
}

hypre_ParCSRMatrix *hypreParCSREliminateRowsWithCols(hypre_ParCSRMatrix *A,
                                                     const mfem::Array<int> &rows)
{
  hypre_error_flag = 0;

  HYPRE_Int nrows_local = hypre_ParCSRMatrixNumRows(A);
  HYPRE_Int ncols_local = hypre_ParCSRMatrixNumCols(A);

  HYPRE_Int *diag_rows_bc, *offd_rows_bc;

  hypre_ParCSRMatrix *At, *B;

  HYPRE_Int i, j, k;

  diag_rows_bc = mfem_hypre_CTAlloc_host(HYPRE_Int, nrows_local);

  // Which of the local rows are to be eliminated.
  for (i = 0; i < rows.Size(); i++)
  {
    diag_rows_bc[rows[i]] = 1;
  }

  hypre_ParCSRMatrixTranspose(A, &At, 1);
  hypre_MatvecCommPkgCreate(At);

  // Use a Matvec communication pattern to find which of the rows connected to local columns
  // are to be eliminated.
  {
    hypre_ParCSRCommHandle *comm_handle;
    hypre_ParCSRCommPkg *comm_pkg;
    HYPRE_Int num_sends, *int_buf_data;
    HYPRE_Int index, start;

    offd_rows_bc = mfem_hypre_TAlloc_host(
        HYPRE_Int, hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(At)));

    comm_pkg = hypre_ParCSRMatrixCommPkg(At);
    num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
    int_buf_data = mfem_hypre_TAlloc_host(
        HYPRE_Int, hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
    index = 0;
    for (i = 0; i < num_sends; i++)
    {
      start = hypre_ParCSRCommPkgSendMapStart(comm_pkg, i);
      for (j = start; j < hypre_ParCSRCommPkgSendMapStart(comm_pkg, i + 1); j++)
      {
        k = hypre_ParCSRCommPkgSendMapElmt(comm_pkg, j);
        int_buf_data[index++] = diag_rows_bc[k];
      }
    }
    comm_handle = hypre_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, offd_rows_bc);

    // Finish the communication.
    hypre_ParCSRCommHandleDestroy(comm_handle);

    hypre_TFree(int_buf_data, HYPRE_MEMORY_HOST);
  }

  // Eliminate the columns of the original matrix (rows in the transposed matrix).
  {
    hypre_CSRMatrix *Atd = hypre_ParCSRMatrixDiag(At);
    HYPRE_Real *AtdA = hypre_CSRMatrixData(Atd);
    HYPRE_Int *AtdI = hypre_CSRMatrixI(Atd);
    HYPRE_Int *AtdJ = hypre_CSRMatrixJ(Atd);

    hypre_CSRMatrix *Ato = hypre_ParCSRMatrixOffd(At);
    HYPRE_Real *AtoA = hypre_CSRMatrixData(Ato);
    HYPRE_Int *AtoI = hypre_CSRMatrixI(Ato);
    HYPRE_Int *AtoJ = hypre_CSRMatrixJ(Ato);

    HYPRE_Int elim;

    for (i = 0; i < ncols_local; i++)
    {
      // A column is eliminated if it has a nonzero in an eliminated row.
      elim = 0;
      for (j = AtdI[i]; j < AtdI[i + 1]; j++)
      {
        if (diag_rows_bc[AtdJ[j]])
        {
          elim = 1;
          break;
        }
      }
      if (!elim && AtoI)
      {
        for (j = AtoI[i]; j < AtoI[i + 1]; j++)
        {
          if (offd_rows_bc[AtoJ[j]])
          {
            elim = 1;
            break;
          }
        }
      }
      if (elim)
      {
        for (j = AtdI[i]; j < AtdI[i + 1]; j++)
        {
          // if (!diag_rows_bc[AtdJ[j]])
          AtdA[j] = 0.0;
        }
        if (AtoI)
        {
          for (j = AtoI[i]; j < AtoI[i + 1]; j++)
          {
            // if (!offd_rows_bc[AtoJ[j]])
            AtoA[j] = 0.0;
          }
        }
      }
    }
  }

  hypre_TFree(diag_rows_bc, HYPRE_MEMORY_HOST);
  hypre_TFree(offd_rows_bc, HYPRE_MEMORY_HOST);

  // Create as a new matrix.
  hypre_ParCSRMatrixTranspose(At, &B, 1);
  hypre_MatvecCommPkgCreate(B);
  hypre_ParCSRMatrixDestroy(At);

  MFEM_VERIFY(!hypre_error_flag,
              "HYPRE error encountered: Error code = " << hypre_error_flag << "!");
  return B;
}

hypre_ParCSRMatrix *hypreParCSREliminateColsWithRows(hypre_ParCSRMatrix *A,
                                                     const mfem::Array<int> &cols)
{
  hypre_error_flag = 0;

  HYPRE_Int nrows_local = hypre_ParCSRMatrixNumRows(A);

  HYPRE_Int *diag_cols_bc, *offd_cols_bc;

  hypre_ParCSRMatrix *B;

  HYPRE_Int i, j, k;

  diag_cols_bc = mfem_hypre_CTAlloc_host(HYPRE_Int, nrows_local);

  // Which of the local columns are to be eliminated.
  for (i = 0; i < cols.Size(); i++)
  {
    diag_cols_bc[cols[i]] = 1;
  }

  // Clone the original matrix.
  B = hypre_ParCSRMatrixClone(A, 1);
  hypre_MatvecCommPkgCreate(B);

  // Use a Matvec communication pattern to find which of the off-diagonal columns are to be
  // eliminated.
  {
    hypre_ParCSRCommHandle *comm_handle;
    hypre_ParCSRCommPkg *comm_pkg;
    HYPRE_Int num_sends, *int_buf_data;
    HYPRE_Int index, start;

    offd_cols_bc = mfem_hypre_TAlloc_host(
        HYPRE_Int, hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(B)));

    comm_pkg = hypre_ParCSRMatrixCommPkg(B);
    num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
    int_buf_data = mfem_hypre_TAlloc_host(
        HYPRE_Int, hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
    index = 0;
    for (i = 0; i < num_sends; i++)
    {
      start = hypre_ParCSRCommPkgSendMapStart(comm_pkg, i);
      for (j = start; j < hypre_ParCSRCommPkgSendMapStart(comm_pkg, i + 1); j++)
      {
        k = hypre_ParCSRCommPkgSendMapElmt(comm_pkg, j);
        int_buf_data[index++] = diag_cols_bc[k];
      }
    }
    comm_handle = hypre_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, offd_cols_bc);

    // Finish the communication.
    hypre_ParCSRCommHandleDestroy(comm_handle);

    hypre_TFree(int_buf_data, HYPRE_MEMORY_HOST);
  }

  // Eliminate the rows of the original matrix.
  {
    hypre_CSRMatrix *Bd = hypre_ParCSRMatrixDiag(B);
    HYPRE_Real *BdA = hypre_CSRMatrixData(Bd);
    HYPRE_Int *BdI = hypre_CSRMatrixI(Bd);
    HYPRE_Int *BdJ = hypre_CSRMatrixJ(Bd);

    hypre_CSRMatrix *Bo = hypre_ParCSRMatrixOffd(B);
    HYPRE_Real *BoA = hypre_CSRMatrixData(Bo);
    HYPRE_Int *BoI = hypre_CSRMatrixI(Bo);
    HYPRE_Int *BoJ = hypre_CSRMatrixJ(Bo);

    HYPRE_Int elim;

    for (i = 0; i < nrows_local; i++)
    {
      // A column is eliminated if it has a nonzero in an eliminated row.
      elim = 0;
      for (j = BdI[i]; j < BdI[i + 1]; j++)
      {
        if (diag_cols_bc[BdJ[j]])
        {
          elim = 1;
          break;
        }
      }
      if (!elim && BoI)
      {
        for (j = BoI[i]; j < BoI[i + 1]; j++)
        {
          if (offd_cols_bc[BoJ[j]])
          {
            elim = 1;
            break;
          }
        }
      }
      if (elim)
      {
        for (j = BdI[i]; j < BdI[i + 1]; j++)
        {
          // if (!diag_cols_bc[BdJ[j]])
          BdA[j] = 0.0;
        }
        if (BoI)
        {
          for (j = BoI[i]; j < BoI[i + 1]; j++)
          {
            // if (!offd_cols_bc[BoJ[j]])
            BoA[j] = 0.0;
          }
        }
      }
    }
  }

  hypre_TFree(diag_cols_bc, HYPRE_MEMORY_HOST);
  hypre_TFree(offd_cols_bc, HYPRE_MEMORY_HOST);

  MFEM_VERIFY(!hypre_error_flag,
              "HYPRE error encountered: Error code = " << hypre_error_flag << "!");
  return B;
}

void hypreParCSRCopy(hypre_ParCSRMatrix *A, hypre_ParCSRMatrix *B)
{
  hypre_error_flag = 0;

  hypre_CSRMatrix *A_diag = hypre_ParCSRMatrixDiag(A);
  HYPRE_Real *A_diag_a = hypre_CSRMatrixData(A_diag);
  HYPRE_Int *A_diag_i = hypre_CSRMatrixI(A_diag);
  HYPRE_Int *A_diag_j = hypre_CSRMatrixJ(A_diag);
  HYPRE_Int ncols_A_diag = hypre_CSRMatrixNumCols(A_diag);

  hypre_CSRMatrix *A_offd = hypre_ParCSRMatrixOffd(A);
  HYPRE_Real *A_offd_a = hypre_CSRMatrixData(A_offd);
  HYPRE_Int *A_offd_i = hypre_CSRMatrixI(A_offd);
  HYPRE_Int *A_offd_j = hypre_CSRMatrixJ(A_offd);

  HYPRE_BigInt *col_map_offd_A = hypre_ParCSRMatrixColMapOffd(A);

  hypre_CSRMatrix *B_diag = hypre_ParCSRMatrixDiag(B);
  HYPRE_Real *B_diag_a = hypre_CSRMatrixData(B_diag);
  HYPRE_Int *B_diag_i = hypre_CSRMatrixI(B_diag);
  HYPRE_Int *B_diag_j = hypre_CSRMatrixJ(B_diag);
  HYPRE_Int ncols_B_diag = hypre_CSRMatrixNumCols(B_diag);

  hypre_CSRMatrix *B_offd = hypre_ParCSRMatrixOffd(B);
  HYPRE_Real *B_offd_a = hypre_CSRMatrixData(B_offd);
  HYPRE_Int *B_offd_i = hypre_CSRMatrixI(B_offd);
  HYPRE_Int *B_offd_j = hypre_CSRMatrixJ(B_offd);

  HYPRE_BigInt *col_map_offd_B = hypre_ParCSRMatrixColMapOffd(B);

  HYPRE_Int i, j, pos;

  HYPRE_BigInt first_row = hypre_ParCSRMatrixFirstRowIndex(A);
  HYPRE_Int nrows_local = hypre_CSRMatrixNumRows(A_diag);
  MFEM_VERIFY(first_row == hypre_ParCSRMatrixFirstRowIndex(B) &&
                  nrows_local == hypre_CSRMatrixNumRows(B_diag) &&
                  ncols_A_diag == ncols_B_diag,
              "Invalid mismatch in matrix sizes/distribution!");

  // Copy the diagonal block A => B.
  {
    HYPRE_Int *marker = mfem_hypre_CTAlloc_host(HYPRE_Int, ncols_A_diag);
    for (j = 0; j < ncols_A_diag; j++)
    {
      marker[j] = -1;
    }

    for (i = 0; i < nrows_local; i++)
    {
      for (j = A_diag_i[i]; j < A_diag_i[i + 1]; j++)
      {
        marker[A_diag_j[j]] = j;
      }

      for (j = B_diag_i[i]; j < B_diag_i[i + 1]; j++)
      {
        // Skip entries not in sparsity pattern of B to copy. All columns of B are marked in
        // the array because sparsity(B) ⊆ sparsity(A).
        pos = marker[B_diag_j[j]];
        MFEM_VERIFY(pos >= A_diag_i[i],
                    "Found nonzero entry of B in copy which is not in A!");
        B_diag_a[j] = A_diag_a[pos];
      }
    }
    mfem_hypre_TFree_host(marker);
  }

  // Copy the off-diagonal block A => B.
  {
    for (i = 0; i < nrows_local; i++)
    {
      std::map<HYPRE_BigInt, HYPRE_Int> marker;
      // std::unordered_map<HYPRE_BigInt, HYPRE_Int> marker;
      for (j = A_offd_i[i]; j < A_offd_i[i + 1]; j++)
      {
        marker.insert(std::make_pair(col_map_offd_A[A_offd_j[j]], j));
      }

      for (j = B_offd_i[i]; j < B_offd_i[i + 1]; j++)
      {
        auto it = marker.find(col_map_offd_B[B_offd_j[j]]);
        MFEM_VERIFY(it != marker.end(),
                    "Found nonzero entry of B in copy which is not in A!");
        pos = it->second;
        B_offd_a[j] = A_offd_a[pos];
      }
    }
  }

  MFEM_VERIFY(!hypre_error_flag,
              "HYPRE error encountered: Error code = " << hypre_error_flag << "!");
}

void hypreParCSRRowSums(hypre_ParCSRMatrix *A, mfem::Vector &rowsums)
{
  hypre_error_flag = 0;

  hypre_CSRMatrix *A_diag = hypre_ParCSRMatrixDiag(A);
  HYPRE_Real *A_diag_a = hypre_CSRMatrixData(A_diag);
  HYPRE_Int *A_diag_i = hypre_CSRMatrixI(A_diag);

  hypre_CSRMatrix *A_offd = hypre_ParCSRMatrixOffd(A);
  HYPRE_Real *A_offd_a = hypre_CSRMatrixData(A_offd);
  HYPRE_Int *A_offd_i = hypre_CSRMatrixI(A_offd);

  HYPRE_Int nrows_local = hypre_CSRMatrixNumRows(A_diag);

  HYPRE_Int i, j;
  HYPRE_Real rowsum;

  for (i = 0; i < nrows_local; i++)
  {
    rowsum = 0.0;
    for (j = A_diag_i[i]; j < A_diag_i[i + 1]; j++)
    {
      rowsum += std::abs(A_diag_a[j]);
    }
    for (j = A_offd_i[i]; j < A_offd_i[i + 1]; j++)
    {
      rowsum += std::abs(A_offd_a[j]);
    }
    rowsums(i) = rowsum;
  }

  MFEM_VERIFY(!hypre_error_flag,
              "HYPRE error encountered: Error code = " << hypre_error_flag << "!");
}

void hypreParCSRInfNorm(hypre_ParCSRMatrix *Ar, hypre_ParCSRMatrix *Ai, HYPRE_Real *norm)
{
  hypre_error_flag = 0;

  MPI_Comm comm = hypre_ParCSRMatrixComm(Ar);

  hypre_CSRMatrix *Ar_diag = hypre_ParCSRMatrixDiag(Ar);
  HYPRE_Real *Ar_diag_a = hypre_CSRMatrixData(Ar_diag);
  HYPRE_Int *Ar_diag_i = hypre_CSRMatrixI(Ar_diag);
  HYPRE_Int *Ar_diag_j = hypre_CSRMatrixJ(Ar_diag);
  HYPRE_Int ncols_Ar_diag = hypre_CSRMatrixNumCols(Ar_diag);

  hypre_CSRMatrix *Ar_offd = hypre_ParCSRMatrixOffd(Ar);
  HYPRE_Real *Ar_offd_a = hypre_CSRMatrixData(Ar_offd);
  HYPRE_Int *Ar_offd_i = hypre_CSRMatrixI(Ar_offd);
  HYPRE_Int *Ar_offd_j = hypre_CSRMatrixJ(Ar_offd);

  HYPRE_BigInt *col_map_offd_Ar = hypre_ParCSRMatrixColMapOffd(Ar);

  hypre_CSRMatrix *Ai_diag = hypre_ParCSRMatrixDiag(Ai);
  HYPRE_Real *Ai_diag_a = hypre_CSRMatrixData(Ai_diag);
  HYPRE_Int *Ai_diag_i = hypre_CSRMatrixI(Ai_diag);
  HYPRE_Int *Ai_diag_j = hypre_CSRMatrixJ(Ai_diag);
  HYPRE_Int ncols_Ai_diag = hypre_CSRMatrixNumCols(Ai_diag);

  hypre_CSRMatrix *Ai_offd = hypre_ParCSRMatrixOffd(Ai);
  HYPRE_Real *Ai_offd_a = hypre_CSRMatrixData(Ai_offd);
  HYPRE_Int *Ai_offd_i = hypre_CSRMatrixI(Ai_offd);
  HYPRE_Int *Ai_offd_j = hypre_CSRMatrixJ(Ai_offd);

  HYPRE_BigInt *col_map_offd_Ai = hypre_ParCSRMatrixColMapOffd(Ai);

  HYPRE_Int *marker_diag;

  HYPRE_BigInt first_row = hypre_ParCSRMatrixFirstRowIndex(Ar);
  HYPRE_Int nrows_local = hypre_CSRMatrixNumRows(Ar_diag);
  MFEM_VERIFY(first_row == hypre_ParCSRMatrixFirstRowIndex(Ai) &&
                  nrows_local == hypre_CSRMatrixNumRows(Ai_diag) &&
                  ncols_Ar_diag == ncols_Ai_diag,
              "Invalid mismatch in matrix sizes/distribution!");

  HYPRE_Int i, j, pos;
  HYPRE_Real rowsum, maxsum = 0.0;

  // We assume the sparsity of the imaginary part is a subset of the real part. Entries
  // outside the sparsity of the real part will be ignored for the calculation of matrix
  // norm.
  marker_diag = mfem_hypre_CTAlloc_host(HYPRE_Int, ncols_Ai_diag);
  for (j = 0; j < ncols_Ai_diag; j++)
  {
    marker_diag[j] = -1;
  }

  for (i = 0; i < nrows_local; i++)
  {
    rowsum = 0.0;

    // Diagonal part
    for (j = Ai_diag_i[i]; j < Ai_diag_i[i + 1]; j++)
    {
      marker_diag[Ai_diag_j[j]] = j;
    }

    for (j = Ar_diag_i[i]; j < Ar_diag_i[i + 1]; j++)
    {
      pos = marker_diag[Ar_diag_j[j]];
      if (pos >= Ai_diag_i[i])
      {
        // Column entry is nonzero in both Ar and Ai.
        rowsum += std::hypot(Ar_diag_a[j], Ai_diag_a[pos]);
      }
      else
      {
        rowsum += std::abs(Ar_diag_a[j]);
      }
    }

    // Off-diagonal part
    std::map<HYPRE_BigInt, HYPRE_Int> marker_offd;
    // std::unordered_map<HYPRE_BigInt, HYPRE_Int> marker_offd;
    for (j = Ai_offd_i[i]; j < Ai_offd_i[i + 1]; j++)
    {
      marker_offd.insert(std::make_pair(col_map_offd_Ai[Ai_offd_j[j]], j));
    }

    for (j = Ar_offd_i[i]; j < Ar_offd_i[i + 1]; j++)
    {
      auto it = marker_offd.find(col_map_offd_Ar[Ar_offd_j[j]]);
      if (it != marker_offd.end())
      {
        // Column entry is nonzero in both Ar and Ai.
        pos = it->second;
        rowsum += std::hypot(Ar_offd_a[j], Ai_offd_a[pos]);
      }
      else
      {
        rowsum += std::abs(Ar_offd_a[j]);
      }
    }

    maxsum = std::max(maxsum, rowsum);
  }

  mfem_hypre_TFree_host(marker_diag);

  MPI_Allreduce(&maxsum, norm, 1, HYPRE_MPI_REAL, MPI_MAX, comm);

  MFEM_VERIFY(!hypre_error_flag,
              "HYPRE error encountered: Error code = " << hypre_error_flag << "!");
}

}  // namespace palace::hypre
