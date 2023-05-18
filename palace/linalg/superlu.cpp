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

SuperLUSolver::SuperLUSolver(MPI_Comm comm, config::LinearSolverData::SymFactType reorder,
                             bool use_3d, int print)
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
  solver.SetEquilibriate(false);
  solver.SetReplaceTinyPivot(false);
  if (reorder == config::LinearSolverData::SymFactType::METIS)
  {
    solver.SetColumnPermutation(mfem::superlu::METIS_AT_PLUS_A);
  }
  else if (reorder == config::LinearSolverData::SymFactType::PARMETIS)
  {
    solver.SetColumnPermutation(mfem::superlu::PARMETIS);
  }
  else
  {
    // Use default
  }
  solver.SetRowPermutation(mfem::superlu::NOROWPERM);
  solver.SetIterativeRefine(mfem::superlu::NOREFINE);
  solver.SetSymmetricPattern(true);  // Always symmetric sparsity pattern
}

void SuperLUSolver::SetOperator(const Operator &op)
{
  // We need to save A because SuperLU does not copy the input matrix. For repeated
  // factorizations, always reuse the sparsity pattern.
  if (A)
  {
    solver.SetFact(mfem::superlu::SamePattern_SameRowPerm);
  }
  const auto *PtAP = dynamic_cast<const ParOperator *>(&op);
  if (PtAP)
  {
    A = std::make_unique<mfem::SuperLURowLocMatrix>(
        const_cast<ParOperator *>(PtAP)->ParallelAssemble());
  }
  else
  {
    A = std::make_unique<mfem::SuperLURowLocMatrix>(op);
  }
  solver.SetOperator(*A);
  height = solver.Height();
  width = solver.Width();
}

}  // namespace palace

#endif
