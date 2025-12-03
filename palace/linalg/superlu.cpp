// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "superlu.hpp"

#if defined(MFEM_USE_SUPERLU)

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
    int npdep = (int)std::pow(2, std::floor(std::log2(std::cbrt(np))));
    while (npdep > 1 && np % npdep > 0)
    {
      npdep /= 2;
    }
    return npdep;
  }
}

}  // namespace

SuperLUSolver::SuperLUSolver(MPI_Comm comm, SymbolicFactorization reorder, bool use_3d,
                             bool reorder_reuse, int print)
  : mfem::Solver(), comm(comm), A(nullptr), solver(comm, GetNpDep(Mpi::Size(comm), use_3d)),
    reorder_reuse(reorder_reuse)
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
}

void SuperLUSolver::SetOperator(const Operator &op)
{
  // For repeated factorizations, always reuse the sparsity pattern.
  if (A && reorder_reuse)
  {
    solver.SetFact(mfem::superlu::SamePattern_SameRowPerm);
  }

  // This is very similar to the MFEM SuperLURowLocMatrix from a HypreParMatrix but avoids
  // using the communicator from the Hypre matrix in the case that the solver is
  // constructed on a different communicator.
  const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(&op);
  MFEM_VERIFY(hA && hA->GetGlobalNumRows() == hA->GetGlobalNumCols(),
              "SuperLUSolver requires a square HypreParMatrix operator!");
  auto *parcsr = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(*hA);
  hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(parcsr);
  hypre_CSRMatrixMigrate(csr, HYPRE_MEMORY_HOST);

  // Create the SuperLURowLocMatrix by taking the internal data from a hypre_CSRMatrix.
  HYPRE_BigInt glob_n = hypre_ParCSRMatrixGlobalNumRows(parcsr);
  HYPRE_BigInt first_row = hypre_ParCSRMatrixFirstRowIndex(parcsr);
  HYPRE_Int n_loc = hypre_CSRMatrixNumRows(csr);
  HYPRE_Int *I = hypre_CSRMatrixI(csr);
  HYPRE_BigInt *J = hypre_CSRMatrixBigJ(csr);
  double *data = hypre_CSRMatrixData(csr);

  // We need to save A because SuperLU does not copy the input matrix. Also clean up the
  // Hypre data structure once we are done with it.
#if !defined(HYPRE_BIGINT)
  A = std::make_unique<mfem::SuperLURowLocMatrix>(comm, n_loc, first_row, glob_n, glob_n, I,
                                                  J, data);
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
  A = std::make_unique<mfem::SuperLURowLocMatrix>(comm, n_loc_int, first_row, glob_n,
                                                  glob_n, II.HostRead(), J, data);
#endif
  solver.SetOperator(*A);
  height = solver.Height();
  width = solver.Width();
  hypre_CSRMatrixDestroy(csr);
}

}  // namespace palace

#endif
