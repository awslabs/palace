// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "solver.hpp"

#include "linalg/mumps.hpp"
#include "linalg/rap.hpp"

#include <fstream> // for matrix export?

namespace palace
{

void ExportHypreParMatrixToMM(const mfem::HypreParMatrix &A, const std::string &filename)
{
  MPI_Comm comm = A.GetComm();
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // Get local CSR data
  hypre_ParCSRMatrix *parcsr = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(A);
  hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(parcsr);
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(parcsr);

  HYPRE_BigInt *col_map = hypre_ParCSRMatrixColMapOffd(parcsr);
  HYPRE_BigInt first_row = hypre_ParCSRMatrixFirstRowIndex(parcsr);
  HYPRE_BigInt glob_rows = hypre_ParCSRMatrixGlobalNumRows(parcsr);
  HYPRE_BigInt glob_cols = hypre_ParCSRMatrixGlobalNumCols(parcsr);

  // Merge diag and offd for this rank
  hypre_CSRMatrix *merged = hypre_MergeDiagAndOffd(parcsr);
  hypre_CSRMatrixMigrate(merged, HYPRE_MEMORY_HOST);

  HYPRE_Int n_loc = hypre_CSRMatrixNumRows(merged);
  HYPRE_Int *I = hypre_CSRMatrixI(merged);
  HYPRE_BigInt *J = hypre_CSRMatrixBigJ(merged);
  double *data = hypre_CSRMatrixData(merged);
  HYPRE_Int nnz_loc = I[n_loc];

  // Gather global NNZ
  HYPRE_BigInt nnz_glob;
  HYPRE_BigInt nnz_loc_big = nnz_loc;
  MPI_Allreduce(&nnz_loc_big, &nnz_glob, 1, HYPRE_MPI_BIG_INT, MPI_SUM, comm);

  // Write sequentially (rank 0 writes header, then each rank appends)
  for (int r = 0; r < size; r++)
  {
    MPI_Barrier(comm);
    if (r == rank)
    {
      std::ofstream file(filename, r == 0 ? std::ios::out : std::ios::app);
      if (r == 0)
      {
        file << "%%MatrixMarket matrix coordinate real general\n";
        file << glob_rows << " " << glob_cols << " " << nnz_glob << "\n";
      }
      for (HYPRE_Int i = 0; i < n_loc; i++)
      {
        for (HYPRE_Int j = I[i]; j < I[i + 1]; j++)
        {
          // Matrix Market is 1-indexed
          file << (first_row + i + 1) << " " << (J[j] + 1) << " " << data[j] << "\n";
        }
      }
    }
  }

  hypre_CSRMatrixDestroy(merged);
}

void ExportVectorToBinary(const mfem::Vector &v, const std::string &filename, MPI_Comm comm)
{
  int rank, size;
  rank = Mpi::Rank(comm);
  size = Mpi::Size(comm);

  // Gather all data to rank 0
  int local_size = v.Size();
  std::vector<int> sizes(size), displs(size);
  MPI_Gather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, comm);

  int global_size = 0;
  if (rank == 0)
  {
    for (int i = 0; i < size; i++)
    {
      displs[i] = global_size;
      global_size += sizes[i];
    }
  }

  std::vector<double> global_data(rank == 0 ? global_size : 0);
  const double *local_data = v.HostRead();
  MPI_Gatherv(local_data, local_size, MPI_DOUBLE,
              global_data.data(), sizes.data(), displs.data(), MPI_DOUBLE, 0, comm);

  if (rank == 0)
  {
    std::ofstream file(filename, std::ios::binary);
    file.write(reinterpret_cast<const char*>(&global_size), sizeof(int));
    file.write(reinterpret_cast<const char*>(global_data.data()), global_size * sizeof(double));
  }
}

template <>
void MfemWrapperSolver<Operator>::SetOperator(const Operator &op)
{
  // Operator is always assembled as a HypreParMatrix.
  if (const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(&op))
  {
    pc->SetOperator(*hA);
  }
  else
  {
    const auto *PtAP = dynamic_cast<const ParOperator *>(&op);
    MFEM_VERIFY(PtAP,
                "MfemWrapperSolver must be able to construct a HypreParMatrix operator!");
    pc->SetOperator(!save_assembled ? *PtAP->StealParallelAssemble()
                                    : PtAP->ParallelAssemble());
  }
  this->height = op.Height();
  this->width = op.Width();
}

template <>
void MfemWrapperSolver<ComplexOperator>::SetOperator(const ComplexOperator &op)
{
  // Assemble the real and imaginary parts, then add.
  // XX TODO: Test complex matrix assembly if coarse solve supports it.
  const mfem::HypreParMatrix *hAr = dynamic_cast<const mfem::HypreParMatrix *>(op.Real());
  const mfem::HypreParMatrix *hAi = dynamic_cast<const mfem::HypreParMatrix *>(op.Imag());
  const ParOperator *PtAPr = nullptr, *PtAPi = nullptr;
  if (op.Real() && !hAr)
  {
    PtAPr = dynamic_cast<const ParOperator *>(op.Real());
    MFEM_VERIFY(PtAPr,
                "MfemWrapperSolver must be able to construct a HypreParMatrix operator!");
    hAr = &PtAPr->ParallelAssemble();
  }
  if (op.Imag() && !hAi)
  {
    PtAPi = dynamic_cast<const ParOperator *>(op.Imag());
    MFEM_VERIFY(PtAPi,
                "MfemWrapperSolver must be able to construct a HypreParMatrix operator!");
    hAi = &PtAPi->ParallelAssemble();
  }
  if (hAr && hAi)
  {
    if (complex_matrix)
    {
      // A = [Ar, Ai]
      //     [Ai, -Ar]
      // We solve A [xr; -xi] = [br; bi]
      mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
      mfem::Array2D<double> block_coeffs(2, 2);
      blocks(0, 0) = hAr;
      blocks(0, 1) = hAi;
      blocks(1, 0) = hAi;
      blocks(1, 1) = hAr;
      block_coeffs(0, 0) = 1.0;
      block_coeffs(0, 1) = 1.0;
      block_coeffs(1, 0) = 1.0;
      block_coeffs(1, 1) = -1.0;
      A.reset(mfem::HypreParMatrixFromBlocks(blocks, &block_coeffs));
    }
    else
    {
      // A = Ar + Ai.
      A.reset(mfem::Add(1.0, *hAr, 1.0, *hAi));
      ExportHypreParMatrixToMM(*A, "matrix.mtx");
      Mpi::Print("Exported matrix to matrix.mtx\n");
    }
    if (PtAPr)
    {
      PtAPr->StealParallelAssemble();
    }
    if (PtAPi)
    {
      PtAPi->StealParallelAssemble();
    }
    if (drop_small_entries)
    {
      DropSmallEntries();
    }
    pc->SetOperator(*A);
    if (!save_assembled)
    {
      A.reset();
    }
  }
  else if (hAr)
  {
    if (drop_small_entries)
    {
      A = std::make_unique<mfem::HypreParMatrix>(*hAr);
      DropSmallEntries();
      pc->SetOperator(*A);
    }
    else
    {
      pc->SetOperator(*hAr);
    }
    if (PtAPr && !save_assembled)
    {
      PtAPr->StealParallelAssemble();
    }
  }
  else if (hAi)
  {
    if (drop_small_entries)
    {
      A = std::make_unique<mfem::HypreParMatrix>(*hAi);
      DropSmallEntries();
      pc->SetOperator(*A);
    }
    else
    {
      pc->SetOperator(*hAi);
    }
    if (PtAPi && !save_assembled)
    {
      PtAPi->StealParallelAssemble();
    }
  }
  else
  {
    MFEM_ABORT("Empty ComplexOperator for MfemWrapperSolver!");
  }
  this->height = op.Height();
  this->width = op.Width();
}

template <>
void MfemWrapperSolver<Operator>::Mult(const Vector &x, Vector &y) const
{
  pc->Mult(x, y);
}

template <>
void MfemWrapperSolver<ComplexOperator>::Mult(const ComplexVector &x,
                                              ComplexVector &y) const
{
  if (pc->Height() == x.Size())
  {
    static bool exported_static = false;
    if (!exported_static)
    {
      // Check if RHS is non-zero
      double norm = 0.0;
      const double *data = x.Real().HostRead();
      for (int i = 0; i < x.Real().Size(); i++) norm += data[i] * data[i];
      if (norm > 1e-30)
      {
        ExportVectorToBinary(x.Real(), "rhs.bin", MPI_COMM_WORLD);
        Mpi::Print("Exported RHS to rhs.bin\n");
        exported_static = true;
      }
    }
    mfem::Array<const Vector *> X(2);
    mfem::Array<Vector *> Y(2);
    X[0] = &x.Real();
    X[1] = &x.Imag();
    Y[0] = &y.Real();
    Y[1] = &y.Imag();
    pc->ArrayMult(X, Y);
  }
  else
  {
    const int Nx = x.Size(), Ny = y.Size();
    Vector X(2 * Nx), Y(2 * Ny), yr, yi;
    X.UseDevice(true);
    Y.UseDevice(true);
    yr.UseDevice(true);
    yi.UseDevice(true);
    linalg::SetSubVector(X, 0, x.Real());
    linalg::SetSubVector(X, Nx, x.Imag());
    pc->Mult(X, Y);
    Y.ReadWrite();
    yr.MakeRef(Y, 0, Ny);
    yi.MakeRef(Y, Ny, Ny);
    // [yr; yi] is the complex conjugate of the solution
    yi *= -1.0;
    y.Real() = yr;
    y.Imag() = yi;
  }
}

template <typename OperType>
void MfemWrapperSolver<OperType>::DropSmallEntries()
{
  const auto nnz_before = A->NNZ();
  A->DropSmallEntries(std::pow(std::numeric_limits<double>::epsilon(), 2));
  const auto nnz_after = A->NNZ();
#if defined(MFEM_USE_MUMPS)
  if (auto *mumps = dynamic_cast<MumpsSolver *>(pc.get()))
  {
    if (reorder_reuse && (num_dropped_entries != 0) &&
        (num_dropped_entries != (nnz_before - nnz_after)))
    {
      // MUMPS errors out if there are any changes to the symmetry pattern after the first
      // factorization so we don't reuse the reordering if the number of dropped entries has
      // changed.
      mumps->SetReorderReuse(false);
    }
    else if (reorder_reuse && (num_dropped_entries == (nnz_before - nnz_after)))
    {
      // Reuse the column ordering if the number of dropped entries has not changed.
      mumps->SetReorderReuse(true);
    }
  }
#endif
  num_dropped_entries = nnz_before - nnz_after;
  Mpi::Print(" Dropping {} small entries in sparse matrix out of {} ({:.1f}%)\n",
             num_dropped_entries, nnz_before,
             (double)(num_dropped_entries) / nnz_before * 100.0);
}

}  // namespace palace
