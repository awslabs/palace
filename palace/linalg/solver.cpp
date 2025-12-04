// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "solver.hpp"
// #include <mfem.hpp> // for PetscPar? in solver.hpp?
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
  mfem::PetscParMatrix *pAr = NULL, *pAi = NULL;
  const ParOperator *PtAPr = nullptr, *PtAPi = nullptr;
  //Operator::Type tid = Operator::PETSC_MATIS; // test
  //mfem::OperatorHandle Arh(tid), Aih(tid); //test
  //if (op.Real() && !hAr)
  if (op.Real() && !pAr)
  {
    Mpi::Print("op.Real calling ParallelAssemble()\n");
    PtAPr = dynamic_cast<const ParOperator *>(op.Real());
    MFEM_VERIFY(PtAPr,
                "MfemWrapperSolver must be able to construct a HypreParMatrix operator!");
    //hAr = &PtAPr->ParallelAssemble();
    std::cout << "calling PtApr->PetscParallelAssemble()\n";
    pAr = &PtAPr->PetscParallelAssemble();
  }
  //if (op.Imag() && !hAi)
  if (op.Imag() && !pAi)
  {
    Mpi::Print("op.Imag calling ParallelAssemble()\n");
    PtAPi = dynamic_cast<const ParOperator *>(op.Imag());
    MFEM_VERIFY(PtAPi,
                "MfemWrapperSolver must be able to construct a HypreParMatrix operator!");
    //hAi = &PtAPi->ParallelAssemble();
    std::cout << "calling PtApi->PetscParallelAssemble()\n";
    pAi = &PtAPi->PetscParallelAssemble();
  }
  //if (hAr && hAi)
  if (pAr && pAi)
  {
    std::cout << "if pAr && pAi\n";
    std::cout << "pAr height: " << pAr->Height() << ", width: " << pAr->Width() << std::endl;
    std::cout << "pAi height: " << pAi->Height() << ", width: " << pAi->Width() << std::endl;
    if (complex_matrix)
    {
      // A = [Ar, -Ai]
      //     [Ai,  Ar]
      /*
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
      */
      std::cout << "block offsets\n";
      mfem::Array<int> block_offsets(3); // number of variables + 1
      block_offsets[0] = 0;
      block_offsets[1] = op.Real()->Height(); // check this?
      std::cout << "block_offsets[1]: " << block_offsets[1] << "\n";
      block_offsets[2] = op.Imag()->Height();
      std::cout << "block_offsets[2]: " << block_offsets[2] << "\n";
      block_offsets.PartialSum();
      block_op = std::make_unique<mfem::BlockOperator>(block_offsets);
      std::cout << "setblock 0, 0\n";
      block_op->SetBlock(0, 0, pAr, 1.0);
      std::cout << "setblock 0, 1\n";
      block_op->SetBlock(0, 1, pAi, -1.0);
      std::cout << "setblock 1, 0\n";
      block_op->SetBlock(1, 0, pAi, 1.0);
      std::cout << "setblock 1, 1\n";
      block_op->SetBlock(1, 1, pAr, 1.0);
      std::cout << "create pA from block_op\n";
      pA = std::make_unique<mfem::PetscParMatrix>(pAr->GetComm(), block_op.get(), Operator::PETSC_MATIS);
    }
    else
    {
      // A = Ar + Ai.
      //A.reset(mfem::Add(1.0, *hAr, 1.0, *hAi));
      std::cout << "create pA from pAr\n";
      pA = std::make_unique<mfem::PetscParMatrix>(*pAr);
      // Create pA as a copy of pAr
Mat matAr = pAr->GetMat();
Mat matAi = pAi->GetMat();
Mat matSum;

MatDuplicate(matAr, MAT_COPY_VALUES, &matSum);
MatAXPY(matSum, 1.0, matAi, DIFFERENT_NONZERO_PATTERN); // compare with SAME_NONZERO_PATTERN?

// Create PetscParMatrix from the sum
pA = std::make_unique<mfem::PetscParMatrix>(matSum, mfem::Operator::PETSC_MATIS);

      //pA.reset(mfem::Add(1.0, *pAr, 1.0, *pAi));
    }

    if (PtAPr)
    {
      //PtAPr->StealParallelAssemble();
      PtAPr->StealPetscParallelAssemble();
    }
    if (PtAPi)
    {
      //PtAPi->StealParallelAssemble();
      PtAPi->StealPetscParallelAssemble();
    }
    if (drop_small_entries)
    {
      DropSmallEntries();
    }
    //pc->SetOperator(*block_op);
    //pc->SetOperator(*A);
    std::cout << "pc->SetOperator(*pA)\n";
    pc->SetOperator(*pA);
    if (!save_assembled)
    {
      A.reset();
      pA.reset();
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
