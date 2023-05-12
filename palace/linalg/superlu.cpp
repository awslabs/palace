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

SuperLUSolver::SuperLUSolver(MPI_Comm comm, int sym_fact_type, bool use_3d, int print_lvl)
  : mfem::Solver(), solver(comm, GetNpDep(Mpi::Size(comm), use_3d))
{
  // Configure the solver.
  if (print_lvl > 1)
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
  solver.SetPrintStatistics(print_lvl > 1);
  solver.SetEquilibriate(false);
  solver.SetReplaceTinyPivot(false);
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
  solver.SetRowPermutation(mfem::superlu::NOROWPERM);
  solver.SetIterativeRefine(mfem::superlu::NOREFINE);
  solver.SetSymmetricPattern(true);  // Always symmetric sparsity pattern
}

void SuperLUSolver::SetOperator(const mfem::Operator &op)
{
  // We need to save A because SuperLU does not copy the input matrix. For repeated
  // factorizations, always reuse the sparsity pattern.
  if (A)
  {
    solver.SetFact(mfem::superlu::SamePattern_SameRowPerm);
  }
  auto *PtAP = const_cast<ParOperator *>(dynamic_cast<const ParOperator *>(&op));
  MFEM_VERIFY(PtAP, "SuperLUSolver requires a ParOperator operator!");
  A = std::make_unique<mfem::SuperLURowLocMatrix>(PtAP->ParallelAssemble());

  // Set up base class.
  solver.SetOperator(*A);
}

}  // namespace palace

#endif
