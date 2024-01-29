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

template <>
void MfemWrapperSolver<ComplexOperator>::SetOperator(const ComplexOperator &op)
{
  // Assemble the real and imaginary parts, then add.
  // XX TODO: Test complex matrix assembly if coarse solve supports it
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
    A.reset(mfem::Add(1.0, *hAr, 1.0, *hAi));
    if (PtAPr)
    {
      PtAPr->StealParallelAssemble();
    }
    if (PtAPi)
    {
      PtAPi->StealParallelAssemble();
    }
    pc->SetOperator(*A);
    if (!save_assembled)
    {
      A.reset();
    }
  }
  else if (hAr)
  {
    pc->SetOperator(*hAr);
    if (PtAPr && !save_assembled)
    {
      PtAPr->StealParallelAssemble();
    }
  }
  else if (hAi)
  {
    pc->SetOperator(*hAi);
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
  mfem::Array<const Vector *> X(2);
  mfem::Array<Vector *> Y(2);
  X[0] = &x.Real();
  X[1] = &x.Imag();
  Y[0] = &y.Real();
  Y[1] = &y.Imag();
  pc->ArrayMult(X, Y);
}

}  // namespace palace
