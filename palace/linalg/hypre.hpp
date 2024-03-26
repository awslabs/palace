// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_HYPRE_HPP
#define PALACE_LINALG_HYPRE_HPP

#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace::hypre
{

// Helper function to initialize HYPRE and control use of GPU at runtime. This will call
// HYPRE_SetMemoryLocation and HYPRE_SetExecutionPolicy to match the mfem::Device
// configuration.
inline void Initialize()
{
  mfem::Hypre::Init();
  // HYPRE_SetSpGemmUseCusparse(1);  // MFEM sets to zero, so leave as is for now
}

//
// Wrapper class for HYPRE's hypre_Vector, which can alias an mfem::Vector object for use
// with HYPRE.
//
class HypreVector
{
private:
  hypre_Vector *vec;

public:
  HypreVector(hypre_Vector *vec = nullptr);
  HypreVector(const Vector &x);
  ~HypreVector();

  auto Size() const { return hypre_VectorSize(vec); }

  void Update(const Vector &x);

  operator hypre_Vector *() const { return vec; }
};

//
// Wrapper class for HYPRE's hypre_CSRMatrix, an alternative to mfem::SparseMatrix with
// increased functionality from HYPRE.
//
class HypreCSRMatrix : public palace::Operator
{
private:
  hypre_CSRMatrix *mat;
  mfem::Array<HYPRE_Int> data_I, data_J;
  bool hypre_own_I;

public:
  HypreCSRMatrix(int h, int w, int nnz);
  HypreCSRMatrix(hypre_CSRMatrix *mat);
  HypreCSRMatrix(const mfem::SparseMatrix &m);
  ~HypreCSRMatrix();

  auto NNZ() const { return hypre_CSRMatrixNumNonzeros(mat); }

  const auto *GetI() const { return hypre_CSRMatrixI(mat); }
  auto *GetI() { return hypre_CSRMatrixI(mat); }
  const auto *GetJ() const { return hypre_CSRMatrixJ(mat); }
  auto *GetJ() { return hypre_CSRMatrixJ(mat); }
  const auto *GetData() const { return hypre_CSRMatrixData(mat); }
  auto *GetData() { return hypre_CSRMatrixData(mat); }

  void AssembleDiagonal(Vector &diag) const override;

  void Mult(const Vector &x, Vector &y) const override;

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override;

  void MultTranspose(const Vector &x, Vector &y) const override;

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override;

  operator hypre_CSRMatrix *() const { return mat; }
};

}  // namespace palace::hypre

#endif  // PALACE_LINALG_HYPRE_HPP
