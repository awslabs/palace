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

HypreCSRMatrix::HypreCSRMatrix(int h, int w, int nnz) : palace::Operator(h, w)
{
  mat = hypre_CSRMatrixCreate(h, w, nnz);
  hypre_CSRMatrixInitialize(mat);
}

HypreCSRMatrix::HypreCSRMatrix(hypre_CSRMatrix *mat) : mat(mat)
{
  height = hypre_CSRMatrixNumRows(mat);
  width = hypre_CSRMatrixNumCols(mat);
}

HypreCSRMatrix::HypreCSRMatrix(mfem::SparseMatrix &m)
  : palace::Operator(m.Height(), m.Width())
{
  // Initialize the internal hypre matrix.
  mat = hypre_CSRMatrixCreate(height, width, m.NumNonZeroElems());

  // Copy the I, J and Data vectors.
  hypre_CSRMatrixI(mat) = m.GetI();
  hypre_CSRMatrixJ(mat) = m.GetJ();
  hypre_CSRMatrixData(mat) = m.GetData();

  // Match the memory class of the sparse matrix.
  hypre_CSRMatrixMemoryLocation(mat) = (m.GetMemoryClass() == mfem::MemoryClass::DEVICE)
                                           ? HYPRE_MEMORY_DEVICE
                                           : HYPRE_MEMORY_HOST;
}

HypreCSRMatrix::HypreCSRMatrix(mfem::SparseMatrix &&m) : HypreCSRMatrix(m)
{
  // Given m is an rvalue, additionally take ownership of the pointers.
  m.SetGraphOwner(false);
  m.SetDataOwner(false);
  hypre_CSRMatrixOwnsData(mat) = true;
}

HypreCSRMatrix::~HypreCSRMatrix()
{
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
