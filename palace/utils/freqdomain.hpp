// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_FREQ_DOMAIN_HPP
#define PALACE_UTILS_FREQ_DOMAIN_HPP

#include <mfem.hpp>
#include "linalg/petsc.hpp"
#include "utils/mfemoperators.hpp"

namespace palace::utils
{

//
// Some utility methods for frequency domain problems.
//

// Convinience method for constructing a the frequency domain matrix-vector product with the
// operator K + iω C - ω² M + A2(ω).
inline std::unique_ptr<petsc::PetscParMatrix> GetSystemMatrixShell(
    double omega, const petsc::PetscParMatrix &K, const petsc::PetscParMatrix &M,
    const petsc::PetscParMatrix *C = nullptr, const petsc::PetscParMatrix *A2 = nullptr)
{
  constexpr auto ExtractReal = petsc::PetscParMatrix::ExtractStructure::REAL;
  constexpr auto ExtractImag = petsc::PetscParMatrix::ExtractStructure::IMAGINARY;
  auto Ar = std::make_unique<SumOperator>(K.GetNumRows(), K.GetNumCols());
  auto Ai = std::make_unique<SumOperator>(K.GetNumRows(), K.GetNumCols());
  if (K.HasReal())
  {
    Ar->AddOperator(*K.GetOperator(ExtractReal));
  }
  if (K.HasImag())
  {
    Ai->AddOperator(*K.GetOperator(ExtractImag));
  }
  if (M.HasReal())
  {
    Ar->AddOperator(*M.GetOperator(ExtractReal), -omega * omega);
  }
  if (M.HasImag())
  {
    Ai->AddOperator(*M.GetOperator(ExtractImag), -omega * omega);
  }
  if (C)
  {
    if (C->HasReal())
    {
      Ai->AddOperator(*C->GetOperator(ExtractReal), omega);
    }
    if (C->HasImag())
    {
      Ar->AddOperator(*C->GetOperator(ExtractImag), -omega);
    }
  }
  if (A2)
  {
    if (A2->HasReal())
    {
      Ar->AddOperator(*A2->GetOperator(ExtractReal));
    }
    if (A2->HasImag())
    {
      Ai->AddOperator(*A2->GetOperator(ExtractImag));
    }
  }
  auto A =
      std::make_unique<petsc::PetscShellMatrix>(K.GetComm(), std::move(Ar), std::move(Ai));
  A->SetSymmetric();
  return A;
}

}  // namespace palace::utils

#endif  // PALACE_UTILS_FREQ_DOMAIN_HPP
