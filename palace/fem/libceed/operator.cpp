// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "operator.hpp"

#include <numeric>
#include <ceed/backend.h>
#include <mfem/general/forall.hpp>
#include "fem/fespace.hpp"
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
  temp.UseDevice(true);
}

Operator::~Operator()
{
  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    Ceed ceed;
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(id < op.size() && op[id],
                "Out of bounds access for thread number " << id << "!");
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedOperatorDestroy(&op[id]));
    PalaceCeedCall(ceed, CeedOperatorDestroy(&op_t[id]));
    PalaceCeedCall(ceed, CeedVectorDestroy(&u[id]));
    PalaceCeedCall(ceed, CeedVectorDestroy(&v[id]));
  }
}

void Operator::AddOper(CeedOperator op_, CeedOperator op_t_)
{
  Ceed ceed;
  CeedSize l_in, l_out;
  CeedVector loc_u, loc_v;
  PalaceCeedCallBackend(CeedOperatorGetCeed(op_, &ceed));
  PalaceCeedCall(ceed, CeedOperatorGetActiveVectorLengths(op_, &l_in, &l_out));
  MFEM_VERIFY((l_in == 0 && l_out == 0) || (mfem::internal::to_int(l_in) == width &&
                                            mfem::internal::to_int(l_out) == height),
              "Dimensions mismatch for CeedOperator!");
  if (op_t_)
  {
    CeedSize l_in_t, l_out_t;
    PalaceCeedCall(ceed, CeedOperatorGetActiveVectorLengths(op_t_, &l_in_t, &l_out_t));
    MFEM_VERIFY((l_in_t == 0 && l_out_t == 0) || (l_in_t == l_out && l_out_t == l_in),
                "Dimensions mismatch for transpose CeedOperator!");
  }
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, l_in, &loc_u));
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, l_out, &loc_v));

  const int id = utils::GetThreadNum();
  MFEM_ASSERT(id < op.size(), "Out of bounds access for thread number " << id << "!");
  op[id] = op_;
  op_t[id] = op_t_;
  u[id] = loc_u;
  v[id] = loc_v;
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
    Ceed ceed;
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(id < op.size() && op[id],
                "Out of bounds access for thread number " << id << "!");
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedVectorSetArray(v[id], mem, CEED_USE_POINTER, diag_data));
    PalaceCeedCall(
        ceed, CeedOperatorLinearAssembleAddDiagonal(op[id], v[id], CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(ceed, CeedVectorTakeArray(v[id], mem, nullptr));
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
    Ceed ceed;
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(id < op.size() && op[id],
                "Out of bounds access for thread number " << id << "!");
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

void CeedOperatorAssembleCOORemoveZeros(Ceed ceed, CeedSize *nnz, CeedInt **rows,
                                        CeedInt **cols, CeedVector *vals, CeedMemType *mem)
{
  // Filter out zero entries. For now, eliminating zeros happens all on the host.
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

  *rows = new_rows;
  *cols = new_cols;
  *vals = new_vals;
  *nnz = q;
}

void CeedOperatorAssembleCOO(const Operator &op, bool skip_zeros, CeedSize *nnz,
                             CeedInt **rows, CeedInt **cols, CeedVector *vals,
                             CeedMemType *mem)
{
  const std::size_t nt = op.Size();
  if (nt == 0)
  {
    *nnz = 0;
    *rows = nullptr;
    *cols = nullptr;
    *vals = nullptr;
    *mem = CEED_MEM_HOST;
    return;
  }

  Ceed ceed;
  CeedScalar *vals_array;
  std::vector<CeedSize> loc_nnz(nt), loc_offsets(nt + 1);
  std::vector<CeedInt *> loc_rows(nt), loc_cols(nt);
  std::vector<CeedVector> loc_vals(nt);

  PalaceCeedCallBackend(CeedOperatorGetCeed(op[0], &ceed));
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, mem));
  if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    *mem = CEED_MEM_HOST;
  }

  PalacePragmaOmp(parallel if (nt > 1))
  {
    Ceed ceed;
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(id < op.Size() && op[id],
                "Out of bounds access for thread number " << id << "!");
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));

    // Assemble sparsity pattern (rows, cols are always host pointers).
    PalaceCeedCall(ceed, CeedOperatorLinearAssembleSymbolic(op[id], &loc_nnz[id],
                                                            &loc_rows[id], &loc_cols[id]));

    // Assemble values.
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, loc_nnz[id], &loc_vals[id]));
    PalaceCeedCall(ceed, CeedOperatorLinearAssemble(op[id], loc_vals[id]));
  }

  loc_offsets[0] = 0;
  std::inclusive_scan(loc_nnz.begin(), loc_nnz.end(), loc_offsets.begin() + 1);
  *nnz = loc_offsets.back();
  if (nt == 1)
  {
    // Assemble values.
    *rows = loc_rows[0];
    *cols = loc_cols[0];
    *vals = loc_vals[0];
  }
  else
  {
    // Global assembly.
    PalaceCeedCall(ceed, CeedInternalCalloc(*nnz, rows));
    PalaceCeedCall(ceed, CeedInternalCalloc(*nnz, cols));
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, *nnz, vals));
    PalaceCeedCall(ceed, CeedVectorGetArrayWrite(*vals, *mem, &vals_array));

    PalacePragmaOmp(parallel if (nt > 1))
    {
      const int id = utils::GetThreadNum();
      MFEM_ASSERT(id < op.Size() && op[id],
                  "Out of bounds access for thread number " << id << "!");
      const auto start = loc_offsets[id];
      const auto end = loc_offsets[id + 1];
      for (auto k = start; k < end; k++)
      {
        (*rows)[k] = loc_rows[id][k - start];
        (*cols)[k] = loc_cols[id][k - start];
      }

      // The CeedVector is on only on device when MFEM is also using the device. This is
      // also correctly a non-OpenMP loop when executed on the CPU (OpenMP is handled in
      // outer scope above).
      Ceed ceed;
      const CeedScalar *loc_vals_array;
      PalaceCeedCallBackend(CeedVectorGetCeed(loc_vals[id], &ceed));
      PalaceCeedCall(ceed, CeedVectorGetArrayRead(loc_vals[id], *mem, &loc_vals_array));
      mfem::forall_switch(*mem == CEED_MEM_DEVICE, end - start,
                          [=] MFEM_HOST_DEVICE(int k)
                          { vals_array[k + start] = loc_vals_array[k]; });
      PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(loc_vals[id], &loc_vals_array));
      PalaceCeedCall(ceed, CeedInternalFree(&loc_rows[id]));
      PalaceCeedCall(ceed, CeedInternalFree(&loc_cols[id]));
      PalaceCeedCall(ceed, CeedVectorDestroy(&loc_vals[id]));
    }

    PalaceCeedCall(ceed, CeedVectorRestoreArray(*vals, &vals_array));
  }

  // std::cout << "  Operator full assembly (COO) has " << *nnz << " NNZ";
  if (skip_zeros && *nnz > 0)
  {
    CeedOperatorAssembleCOORemoveZeros(ceed, nnz, rows, cols, vals, mem);
    // std::cout << " (new NNZ after removal: " << *nnz << ")";
  }
  // std::cout << "\n";
}

}  // namespace

std::unique_ptr<mfem::SparseMatrix> CeedOperatorFullAssemble(const Operator &op,
                                                             bool skip_zeros, bool set)
{
  // First, get matrix on master thread in COO format, withs rows/cols always on host and
  // vals potentially on the device. Process skipping zeros if desired.
  CeedSize nnz;
  CeedInt *rows, *cols;
  CeedVector vals;
  CeedMemType mem;
  CeedOperatorAssembleCOO(op, skip_zeros, &nnz, &rows, &cols, &vals, &mem);

  // Preallocate CSR memory on host (like PETSc's MatSetValuesCOO).
  auto mat = std::make_unique<mfem::SparseMatrix>();
  mat->OverrideSize(op.Height(), op.Width());
  mat->GetMemoryI().New(op.Height() + 1);
  auto *I = mat->GetI();
  mfem::Array<int> J(nnz), perm(nnz), Jmap(nnz + 1);

  for (int i = 0; i < op.Height() + 1; i++)
  {
    I[i] = 0;
  }
  for (int k = 0; k < nnz; k++)
  {
    perm[k] = k;
  }
  std::sort(perm.begin(), perm.end(),
            [&](const int &i, const int &j) { return (rows[i] < rows[j]); });

  int q = -1;  // True nnz index
  for (int k = 0; k < nnz;)
  {
    // Sort column entries in the row.
    const int row = rows[perm[k]];
    const int start = k;
    while (k < nnz && rows[perm[k]] == row)
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
  const int nnz_new = q + 1;

  // Finalize I, Jmap.
  I[0] = 0;
  for (int i = 0; i < op.Height(); i++)
  {
    I[i + 1] += I[i];
  }
  Jmap[0] = 0;
  for (int k = 0; k < nnz_new; k++)
  {
    Jmap[k + 1] += Jmap[k];
  }

  mat->GetMemoryJ().New(nnz_new, mat->GetMemoryJ().GetMemoryType());
  mat->GetMemoryData().New(nnz_new, mat->GetMemoryJ().GetMemoryType());
  {
    // This always executes on the device.
    const auto *d_J_old = J.Read();
    auto *d_J = mfem::Write(mat->GetMemoryJ(), nnz_new);
    mfem::forall(nnz_new, [=] MFEM_HOST_DEVICE(int k) { d_J[k] = d_J_old[k]; });
  }

  // Fill the values (also always on device).
  if (vals)
  {
    auto FillValues = [&](const double *vals_array)
    {
      const auto *d_perm = perm.Read();
      const auto *d_Jmap = Jmap.Read();
      auto *d_A = mfem::Write(mat->GetMemoryData(), nnz_new);
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
      Vector d_vals(nnz);
      d_vals.UseDevice(true);
      {
        auto *d_vals_array = d_vals.HostWrite();
        PalacePragmaOmp(parallel for schedule(static))
        for (int k = 0; k < nnz; k++)
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
    PalaceCeedCall(ceed, CeedInternalFree(&rows));
    PalaceCeedCall(ceed, CeedInternalFree(&cols));
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
  };

  // Initialize the coarse operator.
  auto op_coarse = std::make_unique<SymmetricOperator>(fespace_coarse.GetVSize(),
                                                       fespace_coarse.GetVSize());

  // Assemble the coarse operator by coarsening each sub-operator (over threads, geometry
  // types, integrators) of the original fine operator.
  PalacePragmaOmp(parallel if (op_fine.Size() > 1))
  {
    Ceed ceed;
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(id < op_fine.Size() && op_fine[id],
                "Out of bounds access for thread number " << id << "!");
    PalaceCeedCallBackend(CeedOperatorGetCeed(op_fine[id], &ceed));
    {
      Ceed ceed_parent;
      PalaceCeedCall(ceed, CeedGetParent(ceed, &ceed_parent));
      if (ceed_parent)
      {
        ceed = ceed_parent;
      }
    }

    // Initialize the composite operator on each thread.
    CeedOperator loc_op;
    PalaceCeedCall(ceed, CeedCompositeOperatorCreate(ceed, &loc_op));

    bool composite;
    PalaceCeedCall(ceed, CeedOperatorIsComposite(op_fine[id], &composite));
    if (composite)
    {
      CeedInt nloc_ops_fine;
      CeedOperator *loc_ops_fine;
      PalaceCeedCall(ceed, CeedCompositeOperatorGetNumSub(op_fine[id], &nloc_ops_fine));
      PalaceCeedCall(ceed, CeedCompositeOperatorGetSubList(op_fine[id], &loc_ops_fine));
      for (CeedInt k = 0; k < nloc_ops_fine; k++)
      {
        CeedOperator sub_op;
        SingleOperatorCoarsen(ceed, loc_ops_fine[k], &sub_op);
        PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op, sub_op));
        PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op));
      }
    }
    else
    {
      CeedOperator sub_op;
      SingleOperatorCoarsen(ceed, op_fine[id], &sub_op);
      PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op, sub_op));
      PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op));
    }
    PalaceCeedCall(ceed, CeedOperatorCheckReady(loc_op));
    op_coarse->AddOper(loc_op);  // Thread-safe
  }

  return op_coarse;
}

}  // namespace palace::ceed
