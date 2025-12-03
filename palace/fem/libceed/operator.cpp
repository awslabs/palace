// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "operator.hpp"

#include <numeric>
#include <ceed/backend.h>
#include <mfem.hpp>
#include <mfem/general/forall.hpp>
#include "fem/fespace.hpp"
#include "linalg/hypre.hpp"
#include "utils/omp.hpp"

namespace palace::ceed
{

Operator::Operator(int h, int w) : palace::Operator(h, w)
{
  const std::size_t nt = internal::GetCeedObjects().size();
  op.resize(nt, nullptr);
  op_t.resize(nt, nullptr);
  u.resize(nt, nullptr);
  v.resize(nt, nullptr);
  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed = ceed::internal::GetCeedObjects()[utils::GetThreadNum()];
    CeedOperator loc_op, loc_op_t;
    CeedVector loc_u, loc_v;
    PalaceCeedCall(ceed, CeedOperatorCreateComposite(ceed, &loc_op));
    PalaceCeedCall(ceed, CeedOperatorCreateComposite(ceed, &loc_op_t));
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, width, &loc_u));
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, height, &loc_v));
    op[id] = loc_op;
    op_t[id] = loc_op_t;
    u[id] = loc_u;
    v[id] = loc_v;
  }
  temp.UseDevice(true);
}

Operator::~Operator()
{
  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedOperatorDestroy(&op[id]));
    PalaceCeedCall(ceed, CeedOperatorDestroy(&op_t[id]));
    PalaceCeedCall(ceed, CeedVectorDestroy(&u[id]));
    PalaceCeedCall(ceed, CeedVectorDestroy(&v[id]));
  }
}

void Operator::AddSubOperator(CeedOperator sub_op, CeedOperator sub_op_t)
{
  // This should be called from within a OpenMP parallel region.
  const int id = utils::GetThreadNum();
  MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
              "Out of bounds access for thread number " << id << "!");
  Ceed ceed;
  PalaceCeedCallBackend(CeedOperatorGetCeed(sub_op, &ceed));
  CeedSize l_in, l_out;
  PalaceCeedCall(ceed, CeedOperatorGetActiveVectorLengths(sub_op, &l_in, &l_out));
  MFEM_VERIFY((l_in < 0 || mfem::internal::to_int(l_in) == width) &&
                  (l_out < 0 || mfem::internal::to_int(l_out) == height),
              "Dimensions mismatch for CeedOperator!");
  PalaceCeedCall(ceed, CeedOperatorCompositeAddSub(op[id], sub_op));
  PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op));
  if (sub_op_t)
  {
    Ceed ceed_t;
    PalaceCeedCallBackend(CeedOperatorGetCeed(sub_op_t, &ceed_t));
    MFEM_VERIFY(ceed_t == ceed, "Ceed context mismatch for transpose CeedOperator!");
    CeedSize l_in_t, l_out_t;
    PalaceCeedCall(ceed, CeedOperatorGetActiveVectorLengths(sub_op_t, &l_in_t, &l_out_t));
    MFEM_VERIFY(l_in_t == l_out && l_out_t == l_in,
                "Dimensions mismatch for transpose CeedOperator!");
    PalaceCeedCall(ceed, CeedOperatorCompositeAddSub(op_t[id], sub_op_t));
    PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op_t));
  }
}

void Operator::Finalize()
{
  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedOperatorCheckReady(op[id]));
    PalaceCeedCall(ceed, CeedOperatorCheckReady(op_t[id]));
  }
}

void Operator::DestroyAssemblyData() const
{
  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(op[id]));
  }
}

void Operator::AssembleDiagonal(Vector &diag) const
{
  Ceed ceed;
  CeedMemType mem;
  MFEM_VERIFY(diag.Size() == height, "Invalid size for diagonal vector!");
  diag = 0.0;
  PalaceCeedCallBackend(CeedOperatorGetCeed(op[0], &ceed));
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem == CEED_MEM_DEVICE)
  {
    mem = CEED_MEM_HOST;
  }
  auto *diag_data = diag.ReadWrite(mem == CEED_MEM_DEVICE);

  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedVectorSetArray(v[id], mem, CEED_USE_POINTER, diag_data));
    PalaceCeedCall(
        ceed, CeedOperatorLinearAssembleAddDiagonal(op[id], v[id], CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(ceed, CeedVectorTakeArray(v[id], mem, nullptr));
    PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(op[id]));
  }
}

namespace
{

inline void CeedAddMult(const std::vector<CeedOperator> &op,
                        const std::vector<CeedVector> &u, const std::vector<CeedVector> &v,
                        const Vector &x, Vector &y)
{
  Ceed ceed;
  CeedMemType mem;
  PalaceCeedCallBackend(CeedOperatorGetCeed(op[0], &ceed));
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem == CEED_MEM_DEVICE)
  {
    mem = CEED_MEM_HOST;
  }
  const auto *x_data = x.Read(mem == CEED_MEM_DEVICE);
  auto *y_data = y.ReadWrite(mem == CEED_MEM_DEVICE);

  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedVectorSetArray(u[id], mem, CEED_USE_POINTER,
                                            const_cast<CeedScalar *>(x_data)));
    PalaceCeedCall(ceed, CeedVectorSetArray(v[id], mem, CEED_USE_POINTER, y_data));
    PalaceCeedCall(ceed,
                   CeedOperatorApplyAdd(op[id], u[id], v[id], CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(ceed, CeedVectorTakeArray(u[id], mem, nullptr));
    PalaceCeedCall(ceed, CeedVectorTakeArray(v[id], mem, nullptr));
  }
}

}  // namespace

void Operator::Mult(const Vector &x, Vector &y) const
{
  y = 0.0;
  CeedAddMult(op, u, v, x, y);
  if (dof_multiplicity.Size() > 0)
  {
    y *= dof_multiplicity;
  }
}

void Operator::AddMult(const Vector &x, Vector &y, const double a) const
{
  MFEM_VERIFY(a == 1.0, "ceed::Operator::AddMult only supports coefficient = 1.0!");
  if (dof_multiplicity.Size() > 0)
  {
    temp.SetSize(height);
    temp = 0.0;
    CeedAddMult(op, u, v, x, temp);
    {
      const auto *d_dof_multiplicity = dof_multiplicity.Read();
      const auto *d_temp = temp.Read();
      auto *d_y = y.ReadWrite();
      mfem::forall(height, [=] MFEM_HOST_DEVICE(int i)
                   { d_y[i] += d_dof_multiplicity[i] * d_temp[i]; });
    }
  }
  else
  {
    CeedAddMult(op, u, v, x, y);
  }
}

void Operator::MultTranspose(const Vector &x, Vector &y) const
{
  y = 0.0;
  AddMultTranspose(x, y);
}

void Operator::AddMultTranspose(const Vector &x, Vector &y, const double a) const
{
  MFEM_VERIFY(a == 1.0,
              "ceed::Operator::AddMultTranspose only supports coefficient = 1.0!");
  if (dof_multiplicity.Size() > 0)
  {
    temp.SetSize(height);
    {
      const auto *d_dof_multiplicity = dof_multiplicity.Read();
      const auto *d_x = x.Read();
      auto *d_temp = temp.Write();
      mfem::forall(height, [=] MFEM_HOST_DEVICE(int i)
                   { d_temp[i] = d_dof_multiplicity[i] * d_x[i]; });
    }
    CeedAddMult(op_t, v, u, temp, y);
  }
  else
  {
    CeedAddMult(op_t, v, u, x, y);
  }
}

namespace
{

int CeedInternalCallocArray(size_t n, size_t unit, void *p)
{
  *(void **)p = calloc(n, unit);
  MFEM_ASSERT(!n || !unit || *(void **)p,
              "calloc failed to allocate " << n << " members of size " << unit << "!");
  return 0;
}

int CeedInternalFree(void *p)
{
  free(*(void **)p);
  *(void **)p = nullptr;
  return 0;
}

#define CeedInternalCalloc(n, p) CeedInternalCallocArray((n), sizeof(**(p)), p)

void CeedOperatorAssembleCOO(Ceed ceed, CeedOperator op, bool skip_zeros, CeedSize *nnz,
                             CeedInt **rows, CeedInt **cols, CeedVector *vals,
                             CeedMemType *mem)
{
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, mem));

  // Assemble sparsity pattern (rows, cols are always host pointers).
  PalaceCeedCall(ceed, CeedOperatorLinearAssembleSymbolic(op, nnz, rows, cols));

  // Assemble values.
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, *nnz, vals));
  PalaceCeedCall(ceed, CeedOperatorLinearAssemble(op, *vals));

  // Filter out zero entries. For now, eliminating zeros happens all on the host.
  // std::cout << "  Operator full assembly (COO) has " << *nnz << " NNZ";
  if (skip_zeros && *nnz > 0)
  {
    // XX TODO: Use Thrust for this (thrust::copy_if and thrust::zip_iterator)
    CeedInt *new_rows, *new_cols;
    PalaceCeedCall(ceed, CeedInternalCalloc(*nnz, &new_rows));
    PalaceCeedCall(ceed, CeedInternalCalloc(*nnz, &new_cols));

    CeedVector new_vals;
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, *nnz, &new_vals));

    CeedSize q = 0;
    const CeedScalar *vals_array;
    CeedScalar *new_vals_array;
    PalaceCeedCall(ceed, CeedVectorGetArrayRead(*vals, CEED_MEM_HOST, &vals_array));
    PalaceCeedCall(ceed, CeedVectorGetArrayWrite(new_vals, CEED_MEM_HOST, &new_vals_array));
    for (CeedSize k = 0; k < *nnz; k++)
    {
      if (vals_array[k] != 0.0)
      {
        new_rows[q] = (*rows)[k];
        new_cols[q] = (*cols)[k];
        new_vals_array[q] = vals_array[k];
        q++;
      }
    }
    PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(*vals, &vals_array));
    PalaceCeedCall(ceed, CeedVectorRestoreArray(new_vals, &new_vals_array));

    PalaceCeedCall(ceed, CeedInternalFree(rows));
    PalaceCeedCall(ceed, CeedInternalFree(cols));
    PalaceCeedCall(ceed, CeedVectorDestroy(vals));

    *nnz = q;
    *rows = new_rows;
    *cols = new_cols;
    *vals = new_vals;

    // std::cout << " (new NNZ after removal: " << *nnz << ")";
  }
  // std::cout << "\n";
}

std::unique_ptr<hypre::HypreCSRMatrix> OperatorCOOtoCSR(Ceed ceed, CeedInt m, CeedInt n,
                                                        CeedSize nnz, CeedInt *rows,
                                                        CeedInt *cols, CeedVector vals,
                                                        CeedMemType mem, bool set)
{
  // Preallocate CSR memory on host (like PETSc's MatSetValuesCOO). Check for overflow for
  // large nonzero counts.
  const int nnz_int = mfem::internal::to_int(nnz);
  mfem::Array<int> I(m + 1), J(nnz_int), perm(nnz_int), Jmap(nnz_int + 1);
  I = 0;
  for (int k = 0; k < nnz_int; k++)
  {
    perm[k] = k;
  }
  std::sort(perm.begin(), perm.end(),
            [&](const int &i, const int &j) { return (rows[i] < rows[j]); });

  int q = -1;  // True nnz index
  for (int k = 0; k < nnz_int;)
  {
    // Sort column entries in the row.
    const int row = rows[perm[k]];
    const int start = k;
    while (k < nnz_int && rows[perm[k]] == row)
    {
      k++;
    }
    std::sort(perm.begin() + start, perm.begin() + k,
              [&](const int &i, const int &j) { return (cols[i] < cols[j]); });

    q++;
    I[row + 1] = 1;
    J[q] = cols[perm[start]];
    Jmap[q + 1] = 1;
    for (int p = start + 1; p < k; p++)
    {
      if (cols[perm[p]] != cols[perm[p - 1]])
      {
        // New nonzero.
        q++;
        I[row + 1]++;
        J[q] = cols[perm[p]];
        Jmap[q + 1] = 1;
      }
      else
      {
        Jmap[q + 1]++;
      }
    }
  }
  PalaceCeedCall(ceed, CeedInternalFree(&rows));
  PalaceCeedCall(ceed, CeedInternalFree(&cols));

  // Finalize I, Jmap.
  const int nnz_new = q + 1;
  I[0] = 0;
  for (int i = 0; i < m; i++)
  {
    I[i + 1] += I[i];
  }
  Jmap[0] = 0;
  for (int k = 0; k < nnz_new; k++)
  {
    Jmap[k + 1] += Jmap[k];
  }

  // Construct and fill the final CSR matrix. On GPU, MFEM and Hypre share the same memory
  // space. On CPU, the inner nested OpenMP loop (if enabled in MFEM) should be ignored.
  auto mat = std::make_unique<hypre::HypreCSRMatrix>(m, n, nnz_new);
  {
    const auto *d_I_old = I.Read();
    auto *d_I = mat->GetI();
    mfem::forall(m + 1, [=] MFEM_HOST_DEVICE(int i) { d_I[i] = d_I_old[i]; });
  }
  {
    const auto *d_J_old = J.Read();
    auto *d_J = mat->GetJ();
    mfem::forall(nnz_new, [=] MFEM_HOST_DEVICE(int k) { d_J[k] = d_J_old[k]; });
  }
  {
    auto FillValues = [&](const double *vals_array)
    {
      const auto *d_perm = perm.Read();
      const auto *d_Jmap = Jmap.Read();
      auto *d_A = mat->GetData();
      if (set)
      {
        mfem::forall(nnz_new, [=] MFEM_HOST_DEVICE(int k)
                     { d_A[k] = vals_array[d_perm[d_Jmap[k]]]; });
      }
      else
      {
        mfem::forall(nnz_new,
                     [=] MFEM_HOST_DEVICE(int k)
                     {
                       double sum = 0.0;
                       for (int p = d_Jmap[k]; p < d_Jmap[k + 1]; p++)
                       {
                         sum += vals_array[d_perm[p]];
                       }
                       d_A[k] = sum;
                     });
      }
    };
    Ceed ceed;
    const CeedScalar *vals_array;
    PalaceCeedCallBackend(CeedVectorGetCeed(vals, &ceed));
    PalaceCeedCall(ceed, CeedVectorGetArrayRead(vals, mem, &vals_array));
    if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem != CEED_MEM_DEVICE)
    {
      // Copy values to device before filling.
      Vector d_vals(nnz_int);
      {
        auto *d_vals_array = d_vals.HostWrite();
        PalacePragmaOmp(parallel for schedule(static))
        for (int k = 0; k < nnz_int; k++)
        {
          d_vals_array[k] = vals_array[k];
        }
      }
      FillValues(d_vals.Read());
    }
    else
    {
      // No copy required.
      FillValues(vals_array);
    }
    PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(vals, &vals_array));
    PalaceCeedCall(ceed, CeedVectorDestroy(&vals));
  }

  return mat;
}

}  // namespace

std::unique_ptr<hypre::HypreCSRMatrix> CeedOperatorFullAssemble(const Operator &op,
                                                                bool skip_zeros, bool set)
{
  // Assemble operators on each thread.
  std::vector<std::unique_ptr<hypre::HypreCSRMatrix>> loc_mat(op.Size());
  PalacePragmaOmp(parallel if (op.Size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.Size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));

    // Check if the operator is empty, otherwise assemble.
    CeedInt nsub_ops;
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op[id], &nsub_ops));
    if (nsub_ops == 0)
    {
      loc_mat[id] = std::make_unique<hypre::HypreCSRMatrix>(op.Height(), op.Width(), 0);
    }
    else
    {
      // First, get matrix on master thread in COO format, with rows/cols always on host
      // and vals potentially on the device. Process skipping zeros if desired.
      CeedSize nnz;
      CeedInt *rows, *cols;
      CeedVector vals;
      CeedMemType mem;
      CeedOperatorAssembleCOO(ceed, op[id], skip_zeros, &nnz, &rows, &cols, &vals, &mem);
      PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(op[id]));

      // Convert COO to CSR (on each thread). The COO memory is free'd internally.
      loc_mat[id] =
          OperatorCOOtoCSR(ceed, op.Height(), op.Width(), nnz, rows, cols, vals, mem, set);
    }
  }

  // Add CSR matrix objects from each thread (HYPRE's hypre_CSRMatrixAdd uses threads
  // internally as available). We have to scale the duplicated nonzeros when set = true.
  auto mat = std::move(loc_mat[0]);
  std::unique_ptr<hypre::HypreCSRMatrix> b_mat;
  if (set && op.Size() > 1)
  {
    b_mat = std::make_unique<hypre::HypreCSRMatrix>(hypre_CSRMatrixClone(*mat, 0));
    hypre_CSRMatrixSetConstantValues(*b_mat, 1.0);
    for (std::size_t id = 1; id < op.Size(); id++)
    {
      hypre_CSRMatrix *b_loc_mat = hypre_CSRMatrixClone(*loc_mat[id], 0);
      hypre_CSRMatrixSetConstantValues(b_loc_mat, 1.0);
      b_mat = std::make_unique<hypre::HypreCSRMatrix>(
          hypre_CSRMatrixAdd(1.0, *b_mat, 1.0, b_loc_mat));
      hypre_CSRMatrixDestroy(b_loc_mat);
    }
  }
  for (std::size_t id = 1; id < op.Size(); id++)
  {
    mat = std::make_unique<hypre::HypreCSRMatrix>(
        hypre_CSRMatrixAdd(1.0, *mat, 1.0, *loc_mat[id]));
  }
  if (set && op.Size() > 1)
  {
    const auto *d_b_data = b_mat->GetData();
    auto *d_data = mat->GetData();
    mfem::forall(mat->NNZ(),
                 [=] MFEM_HOST_DEVICE(int i) { d_data[i] *= 1.0 / d_b_data[i]; });
  }

  return mat;
}

std::unique_ptr<Operator> CeedOperatorCoarsen(const Operator &op_fine,
                                              const FiniteElementSpace &fespace_coarse)
{
  auto SingleOperatorCoarsen =
      [&fespace_coarse](Ceed ceed, CeedOperator op_fine, CeedOperator *op_coarse)
  {
    CeedBasis basis_fine;
    CeedElemTopology geom;
    PalaceCeedCall(ceed, CeedOperatorGetActiveBasis(op_fine, &basis_fine));
    PalaceCeedCall(ceed, CeedBasisGetTopology(basis_fine, &geom));

    const auto &geom_data =
        fespace_coarse.GetMesh().GetCeedGeomFactorData(ceed).at(GetMfemTopology(geom));
    CeedElemRestriction restr_coarse = fespace_coarse.GetCeedElemRestriction(
        ceed, GetMfemTopology(geom), geom_data.indices);
    CeedBasis basis_coarse = fespace_coarse.GetCeedBasis(ceed, GetMfemTopology(geom));

    PalaceCeedCall(ceed, CeedOperatorMultigridLevelCreate(op_fine, nullptr, restr_coarse,
                                                          basis_coarse, op_coarse, nullptr,
                                                          nullptr));
    PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(*op_coarse));
  };

  // Initialize the coarse operator.
  auto op_coarse = std::make_unique<SymmetricOperator>(fespace_coarse.GetVSize(),
                                                       fespace_coarse.GetVSize());

  // Assemble the coarse operator by coarsening each sub-operator (over threads, geometry
  // types, integrators) of the original fine operator.
  PalacePragmaOmp(parallel if (op_fine.Size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op_fine.Size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op_fine[id], &ceed));
    {
      Ceed ceed_parent;
      PalaceCeedCall(ceed, CeedGetParent(ceed, &ceed_parent));
      if (ceed_parent)
      {
        ceed = ceed_parent;
      }
    }
    CeedInt nsub_ops_fine;
    CeedOperator *sub_ops_fine;
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op_fine[id], &nsub_ops_fine));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op_fine[id], &sub_ops_fine));
    for (CeedInt k = 0; k < nsub_ops_fine; k++)
    {
      CeedOperator sub_op_coarse;
      SingleOperatorCoarsen(ceed, sub_ops_fine[k], &sub_op_coarse);
      op_coarse->AddSubOperator(sub_op_coarse);  // Sub-operator owned by ceed::Operator
    }
  }

  // Finalize the operator (call CeedOperatorCheckReady).
  op_coarse->Finalize();

  return op_coarse;
}

}  // namespace palace::ceed
