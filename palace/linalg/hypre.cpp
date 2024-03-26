// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hypre.hpp"

namespace palace::hypre
{

HypreVector::HypreVector(hypre_Vector *vec) : vec(vec) {}

HypreVector::HypreVector(const Vector &x) : vec(nullptr)
{
  Update(x);
}

HypreVector::~HypreVector()
{
  hypre_SeqVectorDestroy(vec);
}

void HypreVector::Update(const Vector &x)
{
  const HYPRE_Int N = x.Size();
  if (!vec)
  {
    vec = hypre_SeqVectorCreate(N);
    hypre_SeqVectorSetDataOwner(vec, 0);
    hypre_VectorData(vec) = const_cast<double *>(x.Read());
    hypre_SeqVectorInitialize(vec);
  }
  else
  {
    hypre_SeqVectorSetSize(vec, N);
    hypre_VectorData(vec) = const_cast<double *>(x.Read());
  }
}

HypreCSRMatrix::HypreCSRMatrix(int h, int w, int nnz)
  : palace::Operator(h, w), hypre_own_I(true)
{
  mat = hypre_CSRMatrixCreate(h, w, nnz);
  hypre_CSRMatrixInitialize(mat);
}

HypreCSRMatrix::HypreCSRMatrix(hypre_CSRMatrix *mat) : mat(mat), hypre_own_I(true)
{
  height = hypre_CSRMatrixNumRows(mat);
  width = hypre_CSRMatrixNumCols(mat);
}

HypreCSRMatrix::HypreCSRMatrix(const mfem::SparseMatrix &m)
  : palace::Operator(m.Height(), m.Width()), hypre_own_I(false)
{
  const int nnz = m.NumNonZeroElems();
  mat = hypre_CSRMatrixCreate(height, width, nnz);
  hypre_CSRMatrixSetDataOwner(mat, 0);
  hypre_CSRMatrixData(mat) = const_cast<double *>(m.ReadData());
#if !defined(HYPRE_BIGINT)
  hypre_CSRMatrixI(mat) = const_cast<int *>(m.ReadI());
  hypre_CSRMatrixJ(mat) = const_cast<int *>(m.ReadJ());
#else
  data_I.SetSize(height);
  data_J.SetSize(nnz);
  {
    const auto *I = m.ReadI();
    const auto *J = m.ReadJ();
    auto *DI = data_I.Write();
    auto *DJ = data_J.Write();
    mfem::forall(height, [=] MFEM_HOST_DEVICE(int i) { DI[i] = I[i]; });
    mfem::forall(nnz, [=] MFEM_HOST_DEVICE(int i) { DJ[i] = J[i]; });
  }
#endif
  hypre_CSRMatrixInitialize(mat);
}

HypreCSRMatrix::~HypreCSRMatrix()
{
  if (!hypre_own_I)
  {
    hypre_CSRMatrixI(mat) = nullptr;
  }
  hypre_CSRMatrixDestroy(mat);
}

void HypreCSRMatrix::AssembleDiagonal(Vector &diag) const
{
  diag.SetSize(height);
  hypre_CSRMatrixExtractDiagonal(mat, diag.Write(), 0);
}

namespace
{

static HypreVector X, Y;

}  // namespace

void HypreCSRMatrix::Mult(const Vector &x, Vector &y) const
{
  X.Update(x);
  Y.Update(y);
  hypre_CSRMatrixMatvec(1.0, mat, X, 0.0, Y);
}

void HypreCSRMatrix::AddMult(const Vector &x, Vector &y, const double a) const
{
  X.Update(x);
  Y.Update(y);
  hypre_CSRMatrixMatvec(a, mat, X, 1.0, Y);
}

void HypreCSRMatrix::MultTranspose(const Vector &x, Vector &y) const
{
  X.Update(x);
  Y.Update(y);
  hypre_CSRMatrixMatvecT(1.0, mat, X, 0.0, Y);
}

void HypreCSRMatrix::AddMultTranspose(const Vector &x, Vector &y, const double a) const
{
  X.Update(x);
  Y.Update(y);
  hypre_CSRMatrixMatvecT(a, mat, X, 1.0, Y);
}

}  // namespace palace::hypre
