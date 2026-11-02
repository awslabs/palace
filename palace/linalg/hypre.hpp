// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_HYPRE_HPP
#define PALACE_LINALG_HYPRE_HPP

#include <vector>
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

//
// Wrapper for a HypreParMatrix prolongation from true DOFs to local DOFs (the conforming
// prolongation of a nonconforming space, or of a Nedelec space on tetrahedra), which
// applies it and its transpose on the host with the identity rows (most of them) as plain
// copies. Results are identical to those of hypre's products with the wrapped matrix.
//
class HypreProlongation : public palace::Operator
{
private:
  const mfem::HypreParMatrix &P;

  // Forward product: runs of identity rows y[row + k] = x[col + k], the remaining rows of
  // the diagonal block, and the rows with entries in the off-diagonal block.
  std::vector<int> run_row, run_col, run_len, diag_rows, offd_rows;

  // Transpose product: the transposes of the diagonal and off-diagonal blocks (with rows
  // sorted like hypre_CSRMatrixTranspose), runs of identity rows of the transposed diagonal
  // block, and its remaining rows.
  std::vector<int> diagT_i, diagT_j, offdT_i, offdT_j;
  std::vector<double> diagT_a, offdT_a;
  std::vector<int> runT_row, runT_col, runT_len, diagT_rows, offdT_rows;

  // Communication buffers.
  mutable std::vector<double> send_buf, recv_buf;
  mutable std::vector<MPI_Request> requests;
  mutable std::vector<double> packed, unpacked;

  // Exchange n vectors at once using the communication pattern of the matrix, forward
  // (job = 1: data indexed by the send map are sent, off-diagonal columns are received) or
  // transposed (job = 2). The n parts of each array are stored as consecutive blocks.
  void ExchangeStart(int n, int job, const double *send_data, int send_n, int recv_n) const;
  void ExchangeFinish(int n, int job, double *recv_data, int recv_n) const;
  void Exchange2(int job, const double *send_data, int send_n, double *recv_data,
                 int recv_n) const
  {
    ExchangeStart(2, job, send_data, send_n, recv_n);
    ExchangeFinish(2, job, recv_data, recv_n);
  }

public:
  HypreProlongation(const mfem::HypreParMatrix &P);

  const mfem::HypreParMatrix &GetMatrix() const { return P; }

  void Mult(const Vector &x, Vector &y) const override;

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override
  {
    P.AddMult(x, y, a);
  }

  void MultTranspose(const Vector &x, Vector &y) const override;

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override
  {
    P.AddMultTranspose(x, y, a);
  }

  // Products for the real and imaginary parts of a complex vector at once (each identical
  // to the product with a single vector).
  void Mult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi) const;
  void MultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi) const;

  // Products for n vectors stored one after the other in x and y (each identical to the
  // product with a single vector), with one exchange of off-process entries for all.
  void MultBatch(int n, const Vector &x, Vector &y) const;
  void MultTransposeBatch(int n, const Vector &x, Vector &y) const;
};

}  // namespace palace::hypre

#endif  // PALACE_LINALG_HYPRE_HPP
