// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "superlu.hpp"

#if defined(MFEM_USE_SUPERLU)

#include <algorithm>
#include "utils/communication.hpp"

namespace palace
{

namespace
{

int GetNpDep(int np, bool use_3d)
{
  // Return heuristic choice of 3D processor grid depth based on communicator
  // size. Performance doesn't matter here.
  if (!use_3d)
  {
    return 1;
  }
  else
  {
    int npdep = static_cast<int>(std::pow(2, std::floor(std::log2(std::cbrt(np)))));
    while (npdep > 1 && np % npdep > 0)
    {
      npdep /= 2;
    }
    return npdep;
  }
}

}  // namespace

SuperLUSolver::SuperLUSolver(MPI_Comm comm, SymbolicFactorization reorder, bool use_3d,
                             bool, int print)
  : mfem::Solver(), comm(comm), solver(comm, GetNpDep(Mpi::Size(comm), use_3d))
{
  // Configure the solver.
  if (print > 1)
  {
    if (solver.npdep_ > 1)
    {
      Mpi::Print(comm, " SuperLUSolver: Using 3D processor grid {:d} x {:d} x {:d}\n",
                 solver.nprow_, solver.npcol_, solver.npdep_);
    }
    else
    {
      Mpi::Print(comm, " SuperLUSolver: Using 2D processor grid {:d} x {:d}\n",
                 solver.nprow_, solver.npcol_);
    }
  }
  solver.SetPrintStatistics(print > 1);
  solver.SetEquilibriate(false);
  solver.SetReplaceTinyPivot(false);
  switch (reorder)
  {
    case SymbolicFactorization::METIS:
      solver.SetColumnPermutation(mfem::superlu::METIS_AT_PLUS_A);
      break;
    case SymbolicFactorization::PARMETIS:
      solver.SetColumnPermutation(mfem::superlu::PARMETIS);
      break;
    case SymbolicFactorization::AMD:
    case SymbolicFactorization::RCM:
      solver.SetColumnPermutation(mfem::superlu::MMD_AT_PLUS_A);
      break;
    case SymbolicFactorization::SCOTCH:
    case SymbolicFactorization::PTSCOTCH:
    case SymbolicFactorization::PORD:
    case SymbolicFactorization::DEFAULT:
      // Should have good default.
      break;
  }
  // solver.SetRowPermutation(mfem::superlu::NOROWPERM);
  solver.SetIterativeRefine(mfem::superlu::NOREFINE);
  solver.SetSymmetricPattern(true);  // Always symmetric sparsity pattern

  // Relaxed supernodes of up to 200 columns instead of the SuperLU_DIST default of 30 give
  // fewer, larger blocks, and many fewer messages in the distributed triangular solves,
  // which are latency bound when SuperLU_DIST is the coarse solver.
  solver.SetRelaxedSupernodeSize(200);
}

void SuperLUSolver::SetOperator(const Operator &op)
{
  MFEM_VERIFY(!failed, "SuperLUSolver cannot be used after a failed factorization!");

  // This is very similar to the MFEM SuperLURowLocMatrix from a HypreParMatrix but avoids
  // using the communicator from the Hypre matrix in the case that the solver is
  // constructed on a different communicator.
  const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(&op);
  MFEM_VERIFY(hA && hA->GetGlobalNumRows() == hA->GetGlobalNumCols(),
              "SuperLUSolver requires a square HypreParMatrix operator!");
  auto *parcsr = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(*hA);

  // The native solver is bound to one layout (it rejects other local sizes on its own).
  const HYPRE_BigInt new_layout[3] = {hA->Height(), hypre_ParCSRMatrixFirstRowIndex(parcsr),
                                      hypre_ParCSRMatrixGlobalNumRows(parcsr)};
  bool changed = (layout[0] >= 0 && !std::equal(layout, layout + 3, new_layout));
  Mpi::GlobalOr(1, &changed, comm);
  MFEM_VERIFY(!changed, "SuperLUSolver does not support a change of the matrix layout!");
  std::copy(new_layout, new_layout + 3, layout);

  hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(parcsr);
  hypre_CSRMatrixMigrate(csr, HYPRE_MEMORY_HOST);

  // Create the SuperLURowLocMatrix by taking the internal data from a hypre_CSRMatrix.
  HYPRE_BigInt glob_n = hypre_ParCSRMatrixGlobalNumRows(parcsr);
  HYPRE_BigInt first_row = hypre_ParCSRMatrixFirstRowIndex(parcsr);
  HYPRE_Int n_loc = hypre_CSRMatrixNumRows(csr);
  HYPRE_Int *I = hypre_CSRMatrixI(csr);
  HYPRE_BigInt *J = hypre_CSRMatrixBigJ(csr);
  double *data = hypre_CSRMatrixData(csr);

  // Remove the stored exact zeros from this host copy (the coarse operators keep
  // zero-valued entries in their sparsity pattern). With them the factorization was
  // inaccurate for nonsingular matrices, without them accurate: the operator is the same,
  // the pattern of the factorization is not. Every other entry keeps its value and order.
  HYPRE_Int nnz = 0;
  for (HYPRE_Int i = 0, p = 0; i < n_loc; i++)
  {
    for (const HYPRE_Int end = I[i + 1]; p < end; p++)
    {
      if (data[p] != 0.0)
      {
        J[nnz] = J[p];
        data[nnz++] = data[p];
      }
    }
    I[i + 1] = nnz;
  }
  hypre_CSRMatrixNumNonzeros(csr) = nnz;

  // We need to save the matrix because SuperLU does not copy the input matrix. Also clean
  // up the Hypre data structure once we are done with it.
#if !defined(HYPRE_BIGINT)
  pending = std::make_unique<mfem::SuperLURowLocMatrix>(comm, n_loc, first_row, glob_n,
                                                        glob_n, I, J, data);
#else
  int n_loc_int = static_cast<int>(n_loc);
  MFEM_ASSERT(n_loc == (HYPRE_Int)n_loc_int,
              "Overflow error for local sparse matrix size!");
  mfem::Array<int> II(n_loc_int + 1);
  for (int i = 0; i <= n_loc_int; i++)
  {
    II[i] = static_cast<int>(I[i]);
    MFEM_ASSERT(I[i] == (HYPRE_Int)II[i], "Overflow error for local sparse matrix index!");
  }
  pending = std::make_unique<mfem::SuperLURowLocMatrix>(comm, n_loc_int, first_row, glob_n,
                                                        glob_n, II.HostRead(), J, data);
#endif
  height = pending->Height();
  width = pending->Width();
  hypre_CSRMatrixDestroy(csr);
}

}  // namespace palace

#endif
