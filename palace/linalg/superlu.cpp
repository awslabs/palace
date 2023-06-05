// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "superlu.hpp"

#if defined(MFEM_USE_SUPERLU)

#include "linalg/rap.hpp"
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

SuperLUSolver::SuperLUSolver(MPI_Comm comm, int sym_fact_type, bool use_3d, int print)
  : mfem::Solver(), comm(comm), A(nullptr), solver(comm, GetNpDep(Mpi::Size(comm), use_3d))
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
  // solver.SetEquilibriate(false);
  // solver.SetReplaceTinyPivot(false);
  if (sym_fact_type == 2)
  {
    solver.SetColumnPermutation(mfem::superlu::PARMETIS);
  }
  else if (sym_fact_type == 1)
  {
    solver.SetColumnPermutation(mfem::superlu::METIS_AT_PLUS_A);
  }
  else
  {
    // Use default
  }
  // solver.SetRowPermutation(mfem::superlu::NOROWPERM);
  solver.SetIterativeRefine(mfem::superlu::NOREFINE);
  solver.SetSymmetricPattern(true);  // Always symmetric sparsity pattern
}

void SuperLUSolver::SetOperator(const Operator &op)
{
  // For repeated factorizations, always reuse the sparsity pattern. This is very similar to
  // the MFEM SuperLURowLocMatrix from a HypreParMatrix but avoids using the communicator
  // from the Hypre matrix in the case that the solver is constructed on a different
  // communicator.
  if (A)
  {
    solver.SetFact(mfem::superlu::SamePattern_SameRowPerm);
  }
  const mfem::HypreParMatrix *hypA;
  const auto *PtAP = dynamic_cast<const ParOperator *>(&op);
  if (PtAP)
  {
    hypA = &const_cast<ParOperator *>(PtAP)->ParallelAssemble();
  }
  else
  {
    hypA = dynamic_cast<const mfem::HypreParMatrix *>(&op);
    MFEM_VERIFY(hypA, "SuperLUSolver requires a HypreParMatrix operator!");
  }
  hypre_ParCSRMatrix *parcsr =
      (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(*hypA);
  hypA->HostRead();
  hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(parcsr);
  hypA->HypreRead();

  // Create the SuperLURowLocMatrix by taking the internal data from a hypre_CSRMatrix.
  HYPRE_Int n_loc = csr->num_rows;
  HYPRE_BigInt first_row = parcsr->first_row_index;
  HYPRE_Int *I = csr->i;
  HYPRE_BigInt *J = csr->big_j;
  double *data = csr->data;

  // We need to save A because SuperLU does not copy the input matrix. Also clean up the
  // Hypre data structure once we are done with it.
#if !defined(HYPRE_BIGINT)
  A = std::make_unique<mfem::SuperLURowLocMatrix>(comm, n_loc, first_row,
                                                  hypA->GetGlobalNumRows(),
                                                  hypA->GetGlobalNumCols(), I, J, data);
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
  A = std::make_unique<mfem::SuperLURowLocMatrix>(comm, n_loc_int, first_row,
                                                  hypA->GetGlobalNumRows(),
                                                  hypA->GetGlobalNumCols(), II, J, data);
#endif
  solver.SetOperator(*A);
  height = solver.Height();
  width = solver.Width();
  hypre_CSRMatrixDestroy(csr);
}

}  // namespace palace

#endif
