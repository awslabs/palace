// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "solver.hpp"
#include <mfem.hpp> // for PetscPar?
#include <petsc.h>
#include "linalg/mumps.hpp"
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

template <>
void MfemWrapperSolver<ComplexOperator>::SetOperator(const ComplexOperator &op)
{
  // Assemble the real and imaginary parts, then add.
  // XX TODO: Test complex matrix assembly if coarse solve supports it.
  const mfem::HypreParMatrix *hAr = dynamic_cast<const mfem::HypreParMatrix *>(op.Real());
  const mfem::HypreParMatrix *hAi = dynamic_cast<const mfem::HypreParMatrix *>(op.Imag());
  const ParOperator *PtAPr = nullptr, *PtAPi = nullptr;
  //Operator::Type tid = Operator::PETSC_MATIS; // test
  //mfem::OperatorHandle Arh(tid), Aih(tid); //test
  if (op.Real() && !hAr)
  {
    Mpi::Print("calling ParallelAssemble()\n");
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
    /*
    std::cout << "solver.cpp L78\n";
PetscOptionsSetValue(NULL, "-matis_localmat_type", "aij");
    //auto pAr = std::make_unique<mfem::PetscParMatrix>(hAr->GetComm(), op.Real(), Operator::PETSC_MATIS);
    //auto pAi = std::make_unique<mfem::PetscParMatrix>(hAi->GetComm(), op.Imag(), Operator::PETSC_MATIS);
    auto pAr = std::make_unique<mfem::PetscParMatrix>(hAr->GetComm(), hAr, Operator::PETSC_MATIS);
    auto pAi = std::make_unique<mfem::PetscParMatrix>(hAi->GetComm(), hAi, Operator::PETSC_MATIS);
    std::cout << "solver.cpp L83\n";
    mfem::Array<int> block_offsets(3); // number of variables + 1
    block_offsets[0] = 0;
    block_offsets[1] = op.Real()->Height();
    block_offsets[2] = op.Imag()->Height();
    block_offsets.PartialSum();
    std::cout << "solver.cpp L89\n";
    block_op = std::make_unique<mfem::BlockOperator>(block_offsets);
std::cout << "solver.cpp L91\n";
    block_op->SetBlock(0, 0, pAr.get(), 1.0);
    block_op->SetBlock(0, 1, pAi.get(), -1.0);
    block_op->SetBlock(1, 0, pAi.get(), 1.0);
    block_op->SetBlock(1, 1, pAr.get(), 1.0);
    //block_op->SetBlock(0, 0, const_cast<mfem::Operator*>(op.Real()), 1.0);
    //block_op->SetBlock(0, 1, const_cast<mfem::Operator*>(op.Imag()), -1.0);
    //block_op->SetBlock(1, 0, const_cast<mfem::Operator*>(op.Imag()), 1.0);
    //block_op->SetBlock(1, 1, const_cast<mfem::Operator*>(op.Real()), 1.0);
std::cout << "solver.cpp L100\n";
    test_op = std::make_unique<mfem::PetscParMatrix>(hAr->GetComm(), block_op.get(), Operator::PETSC_MATIS);
    std::cout << "solver.cpp L102\n";
    //auto *petsc_mat = dynamic_cast<mfem::PetscParMatrix*>(test_op.get());
    //std::cout << "solver.cpp L108\n";
    //if (petsc_mat) {
    //  std::cout << "solver.cpp L110\n";
      //Mat mat = petsc_mat->A;
    //  std::cout << "solver.cpp L112\n";
    //MatConvert(mat, MATIS, MAT_INPLACE_MATRIX, &mat);
    //std::cout << "solver.cpp L114\n";
    //}
std::cout << "solver.cpp L116\n";
    //pc->SetOperator(*test_op);
    pc->SetOperator(*test_op);
    std::cout << "solver.cpp L118\n";
*/
int rank = Mpi::Rank(hAr->GetComm());
const HYPRE_BigInt *row_starts = hAr->RowPart();
const HYPRE_BigInt *col_starts = hAr->ColPart();

std::cout << "Rank " << rank << " row range: [" << row_starts[rank]
          << ", " << row_starts[rank+1] << ")" << std::endl;

// Check if there are off-diagonal entries (communication between ranks)
hypre_ParCSRMatrix *hypre_A = (hypre_ParCSRMatrix*)(*hAr);
HYPRE_Int num_cols_offd = hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(hypre_A));
std::cout << "Rank " << rank << " off-diagonal cols: " << num_cols_offd << std::endl;


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
    //pc->SetOperator(*block_op);
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
