// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "solver.hpp"

#include "linalg/rap.hpp"

namespace palace
{

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

bool IsSymmetric(const mfem::HypreParMatrix& A, double tol = 1e-12) {
    MPI_Comm comm = A.GetComm();

    // Check if square
    if (A.GetGlobalNumRows() != A.GetGlobalNumCols()) return false;

    // Create A^T
    mfem::HypreParMatrix* AT = A.Transpose();

    // Compute A - A^T
    mfem::HypreParMatrix* diff = mfem::Add(1.0, A, -1.0, *AT);

    // Get Frobenius norm of difference
    double local_norm_sq = 0.0;
    hypre_ParCSRMatrix* parcsr_diff = *diff;
    hypre_CSRMatrix* diag = hypre_ParCSRMatrixDiag(parcsr_diff);
    hypre_CSRMatrix* offd = hypre_ParCSRMatrixOffd(parcsr_diff);

    // Sum squares of diagonal block entries
    HYPRE_Real* diag_data = hypre_CSRMatrixData(diag);
    HYPRE_Int diag_nnz = hypre_CSRMatrixNumNonzeros(diag);
    for (HYPRE_Int i = 0; i < diag_nnz; i++) {
        local_norm_sq += diag_data[i] * diag_data[i];
    }

    // Sum squares of off-diagonal block entries
    HYPRE_Real* offd_data = hypre_CSRMatrixData(offd);
    HYPRE_Int offd_nnz = hypre_CSRMatrixNumNonzeros(offd);
    for (HYPRE_Int i = 0; i < offd_nnz; i++) {
        local_norm_sq += offd_data[i] * offd_data[i];
    }

    // Global reduction
    double global_norm_sq;
    MPI_Allreduce(&local_norm_sq, &global_norm_sq, 1, MPI_DOUBLE, MPI_SUM, comm);

    // Cleanup
    delete AT;
    delete diff;
    std::cout << "sqrt(global_norm_sq): " << sqrt(global_norm_sq) << "\n";
    return sqrt(global_norm_sq) <= tol;
}

bool IsSymmetricSerial(const mfem::HypreParMatrix& A, double tol = 1e-12) {
    // Get CSR data (assuming single processor)
    hypre_ParCSRMatrix* parcsr_A = A;
    hypre_CSRMatrix* csr = hypre_ParCSRMatrixDiag(parcsr_A);

    HYPRE_Int* I = hypre_CSRMatrixI(csr);
    HYPRE_Int* J = hypre_CSRMatrixJ(csr);
    HYPRE_Real* data = hypre_CSRMatrixData(csr);
    HYPRE_Int nrows = hypre_CSRMatrixNumRows(csr);
    HYPRE_Real max_diff = 0.0;
    // Check each entry A(i,j) against A(j,i)
    bool ret = true;
    for (HYPRE_Int i = 0; i < nrows; i++) {
        for (HYPRE_Int k = I[i]; k < I[i+1]; k++) {
            HYPRE_Int j = J[k];
            HYPRE_Real a_ij = data[k];

            // Find A(j,i)
            HYPRE_Real a_ji = 0.0;
            bool found = false;
            for (HYPRE_Int l = I[j]; l < I[j+1]; l++) {
                if (J[l] == i) {
                    a_ji = data[l];
                    found = true;
                    break;
                }
            }
            max_diff = std::max(max_diff, std::abs(a_ij - a_ji) / std::abs(a_ij));
            // Check symmetry
            if (std::abs(a_ij - a_ji) > tol * std::abs(a_ij))
            {
              ret = false;
              //return false;
            }
        }
    }
    std::cout << "max_diff: " << max_diff << "\n";
    return ret;
    //return true;
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
      // A = [Ar, -Ai]
      //     [Ai,  Ar]
      mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
      mfem::Array2D<double> block_coeffs(2, 2);
      blocks(0, 0) = hAr;
      blocks(0, 1) = hAi;
      blocks(1, 0) = hAi;
      blocks(1, 1) = hAr;
      block_coeffs(0, 0) = 1.0;
      block_coeffs(0, 1) = -1.0;
      block_coeffs(1, 0) = 1.0;
      block_coeffs(1, 1) = 1.0;
      A.reset(mfem::HypreParMatrixFromBlocks(blocks, &block_coeffs));
    }
    else
    {
      // A = Ar + Ai.
      A.reset(mfem::Add(1.0, *hAr, 1.0, *hAi));
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
      std::cout << "before NNZ: " << A->NNZ() << " IsSymmetric(A): " << IsSymmetric(*A.get()) << " SerialVersion: " << IsSymmetricSerial(*A.get()) << "\n";
      A->DropSmallEntries(std::numeric_limits<double>::epsilon());
      std::cout << "after NNZ: " << A->NNZ() << " IsSymmetric(A): " << IsSymmetric(*A.get()) << " SerialVersion: " << IsSymmetricSerial(*A.get()) << "\n";
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
      A->DropSmallEntries(std::numeric_limits<double>::epsilon());
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
      A->DropSmallEntries(std::numeric_limits<double>::epsilon());
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
    y.Real() = yr;
    y.Imag() = yi;
  }
}

}  // namespace palace
