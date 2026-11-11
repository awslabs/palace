// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_HYPRE_HPP
#define PALACE_LINALG_HYPRE_HPP

#include <cstdint>
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

  // Compact copy of the remaining rows of the diagonal block (row r is diag_rows[r]). The
  // values of these rows and of the transposed diagonal block are stored as indices into a
  // table of the distinct values of the diagonal block when there are at most 65536 of
  // them (diag_v, diagT_v), and as values otherwise (diag_a, diagT_a).
  std::vector<int> diag_ia, diag_ja;
  std::vector<double> diag_a, table;
  std::vector<std::uint16_t> diag_v, diagT_v;

  // Communication buffers.
  mutable std::vector<double> send_buf, recv_buf;
  mutable std::vector<MPI_Request> requests;
  mutable std::vector<double> packed, unpacked;

  // Exchange two vectors at once using the communication pattern of the matrix, forward
  // (job = 1: data indexed by the send map are sent, off-diagonal columns are received) or
  // transposed (job = 2). The two parts of each array are stored as consecutive blocks.
  void Exchange2(int job, const double *send_data, int send_n, double *recv_data,
                 int recv_n) const;

  // The remaining rows of the diagonal block (diag_rows), the listed rows of the transposed
  // diagonal block, and a single row of it, for one or two vectors.
  void DiagRowSums(const double *X, double *Y) const;
  void DiagRowSums(const double *XR, const double *XI, double *YR, double *YI) const;
  void DiagTRowSums(const std::vector<int> &rows, const double *X, double *Y) const;
  void DiagTRowSums(const std::vector<int> &rows, const double *XR, const double *XI,
                    double *YR, double *YI) const;
  double DiagTRowSum(int t, const double *X) const;
  void DiagTRowSum(int t, const double *XR, const double *XI, double &tr, double &ti) const;

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

  friend class MaskedProlongation;
};

//
// Products with a HypreProlongation P around an operator with eliminated essential true
// dofs: P x̃, where x̃ is x with zeros at the essential true dofs, without forming x̃, and
// transposed products added to a true-dof vector, y += a ỹ, where ỹ is Pᵀ x except at the
// essential true dofs, where it is zero or given (the diagonal policies of ParOperator),
// without forming ỹ. The results are identical to those of the products with the explicit
// vectors, up to the signs of zeros.
//
class MaskedProlongation
{
private:
  const HypreProlongation &P;

  // The essential true dofs (sorted and unique), and as given.
  std::vector<int> list, given;

  // Forward product: identity runs of the unlisted columns, identity rows of listed columns
  // (zero), and the other rows of the diagonal block without the entries of the listed
  // columns (with values or value indices as HypreProlongation), and the positions of the
  // send map at listed true dofs (sent as zero).
  std::vector<int> run_row, run_col, run_len, zero_rows, rows, ia, ja, send_zero;
  std::vector<double> a;
  std::vector<std::uint16_t> v;

  // Transposed product: identity runs and other rows of the transposed diagonal block for
  // the unlisted true dofs without contributions from other processes (added to the result
  // as they are computed), and the unlisted true dofs with such contributions (shared),
  // with their local L-dof if they are identity rows (or -1), and the slot in the shared
  // buffer of each send map entry (-1 for listed true dofs).
  std::vector<int> runT_row, runT_col, runT_len, rowsT, shared, shared_ldof, send_slot;
  mutable std::vector<double> sbuf;

  // Products without L-vectors for local operators applied directly to the true dofs and
  // the extra dofs (see LocalDofMap): the map of the local dofs, the extra index of each
  // compact other row, zero row and row with off-diagonal entries; the identity rows of
  // true dofs with an earlier identity row (rows, columns and extra indices); the sources
  // of the entries of the transposed diagonal block (-1 for the identity row of the true
  // dof itself, the extra index of the row otherwise) and of the transposed off-diagonal
  // block; the unlisted rows of the transposed diagonal block with other entries.
  std::vector<int> lmap, rows_x, zero_x, offd_x, dup_rows, dup_cols, dup_x, diagT_x,
      offdT_x, depT_rows;
  int num_extra = 0;
  std::uint64_t key;

  void ForwardLocal(const double *X, double *Y) const;
  void ForwardLocal(const double *XR, const double *XI, double *YR, double *YI) const;

public:
  MaskedProlongation(const HypreProlongation &P, const mfem::Array<int> &list);

  // Whether the essential true dofs are the given ones.
  bool Matches(const mfem::Array<int> &dofs) const;

  const HypreProlongation &GetProlongation() const { return P; }

  void Mult(const Vector &x, Vector &y) const;
  void Mult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi) const;

  // y += a ỹ with ỹ = Pᵀ x, except at the essential true dofs, where ỹ = z (zero for a null
  // z), for one or two (real and imaginary parts, with a real coefficient) vectors.
  void AddMultTranspose(const Vector &x, Vector &y, double a, const Vector *z) const;
  void AddMultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                        double a, const Vector *zr, const Vector *zi) const;

  // Products for local operators which read and write the true dofs directly: the local
  // dofs which are identity rows of unlisted columns are the true dofs map[l] ≥ 0, and the
  // others ("extra" dofs, the hanging, ghost and essential identity rows) are the entries
  // -1 - map[l] of separate vectors of NumExtra() values. ExtraMult gives the extra values
  // of P x̃ (x̃ as for Mult). Given the contributions y of a local operator accumulated at
  // the identity rows into a true-dof vector and the contributions yx at the extra dofs,
  // TransposeComplete turns y into Pᵀ ỹ (with ỹ the local vector of these contributions)
  // except at the essential true dofs, whose entries are left unspecified. The results are
  // those of Mult and of Pᵀ, up to the signs of zeros.
  const std::vector<int> &LocalDofMap() const { return lmap; }
  int NumExtra() const { return num_extra; }
  std::uint64_t Key() const { return key; }
  const std::vector<int> &EssentialTrueDofs() const { return list; }
  void ExtraMult(const Vector &x, Vector &xx) const;
  void ExtraMult(const Vector &xr, const Vector &xi, Vector &xxr, Vector &xxi) const;
  void TransposeComplete(Vector &y, const Vector &yx) const;
  void TransposeComplete(Vector &yr, Vector &yi, const Vector &yxr,
                         const Vector &yxi) const;
};

}  // namespace palace::hypre

#endif  // PALACE_LINALG_HYPRE_HPP
