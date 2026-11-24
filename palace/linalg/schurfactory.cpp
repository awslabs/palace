// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "schurfactory.hpp"

#include <array>
#include <complex>
#include "linalg/directschur.hpp"
#include "linalg/mumpsschur.hpp"

namespace palace
{

namespace
{

bool Built(LinearSolver type)
{
  switch (type)
  {
#if defined(MFEM_USE_STRUMPACK)
    case LinearSolver::STRUMPACK:
      return true;
#endif
#if defined(MFEM_USE_MUMPS)
    case LinearSolver::MUMPS:
      return true;
#endif
#if defined(MFEM_USE_SUPERLU)
    case LinearSolver::SUPERLU:
      return true;
#endif
    default:
      return false;
  }
}

}  // namespace

LinearSolver DirectSolverFor(LinearSolver requested, Factorization kind)
{
  if (Built(requested))
  {
    return requested;
  }
  const auto order =
      (kind == Factorization::SYSTEM)
          ? std::array{LinearSolver::STRUMPACK, LinearSolver::MUMPS, LinearSolver::SUPERLU}
          : std::array{LinearSolver::MUMPS, LinearSolver::STRUMPACK, LinearSolver::SUPERLU};
  for (const auto type : order)
  {
    if (Built(type))
    {
      return type;
    }
  }
  MFEM_ABORT("Substructuring needs a sparse direct solver: STRUMPACK, MUMPS or "
             "SuperLU_DIST!");
  return LinearSolver::DEFAULT;
}

const char *DirectSolverName(LinearSolver type)
{
  switch (type)
  {
    case LinearSolver::STRUMPACK:
      return "STRUMPACK";
    case LinearSolver::MUMPS:
      return "MUMPS";
    default:
      return "SuperLU_DIST";
  }
}

template <typename T>
std::unique_ptr<SchurFactorization<T>>
MakeSchurFactorization(LinearSolver type, SymbolicFactorization ordering, MPI_Comm comm,
                       HYPRE_BigInt n_glob, int n_loc,
                       typename SchurFactorization<T>::Coo &&coo,
                       const std::vector<HYPRE_BigInt> &schur_vars, std::vector<int> rows,
                       bool spd, bool refactor)
{
#if defined(MFEM_USE_MUMPS)
  if (type == LinearSolver::MUMPS)
  {
    return std::make_unique<MumpsSchurSolverT<T>>(comm, n_glob, n_loc, std::move(coo),
                                                  schur_vars, kBlrTol, refactor, spd,
                                                  std::move(rows), ordering);
  }
#else
  MFEM_CONTRACT_VAR(spd);
  MFEM_CONTRACT_VAR(refactor);
#endif
  return std::make_unique<DirectSchurSolverT<T>>(
      comm, n_glob, n_loc, std::move(coo), schur_vars, std::move(rows), type, ordering);
}

template std::unique_ptr<SchurFactorization<double>>
MakeSchurFactorization<double>(LinearSolver, SymbolicFactorization, MPI_Comm, HYPRE_BigInt,
                               int, SchurFactorization<double>::Coo &&,
                               const std::vector<HYPRE_BigInt> &, std::vector<int>, bool,
                               bool);
template std::unique_ptr<SchurFactorization<std::complex<double>>>
MakeSchurFactorization<std::complex<double>>(
    LinearSolver, SymbolicFactorization, MPI_Comm, HYPRE_BigInt, int,
    SchurFactorization<std::complex<double>>::Coo &&, const std::vector<HYPRE_BigInt> &,
    std::vector<int>, bool, bool);

}  // namespace palace
