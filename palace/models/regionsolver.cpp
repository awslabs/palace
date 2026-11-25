// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "regionsolver.hpp"

#include <complex>
#include <type_traits>
#include <utility>
#include "linalg/mumpsschur.hpp"
#include "linalg/schurfactory.hpp"

namespace palace
{

template <typename T>
RegionSolver<T>::RegionSolver(const BlockCyclicGrid &grid, LinearSolver type,
                              SymbolicFactorization ordering, HYPRE_BigInt n_glob,
                              int n_loc, typename SchurFactorization<T>::Coo &&A,
                              const std::vector<HYPRE_BigInt> &gamma, std::vector<int> rows,
                              const std::function<void(std::vector<T> &)> &add_interface)
  : grid(grid)
{
  MFEM_VERIFY(static_cast<int>(gamma.size()) == grid.Size(),
              "The interface of a region solver is not that of its grid!");
  const bool has_gamma = (grid.Size() > 0);
  factorization = MakeSchurFactorization<T>(
      type, ordering, grid.GetComm(), n_glob, n_loc, std::move(A), gamma, std::move(rows),
      std::is_same_v<T, double>, false, has_gamma ? &grid : nullptr);
  if (!has_gamma)
  {
    return;
  }
  std::vector<T> *C = &dense;
#if defined(MFEM_USE_MUMPS)
  if (auto *mumps = dynamic_cast<MumpsSchurSolverT<T> *>(factorization.get()))
  {
    C = &mumps->SchurBlocks();
  }
  else
#endif
  {
    dense = grid.ScatterFromRoot(factorization->ReleaseSchur());
  }
  add_interface(*C);
  factors = std::make_unique<BlockCyclicFactorization<T>>(grid);
  factors->Factor(*C);
}

template <typename T>
void RegionSolver<T>::Mult(const std::vector<const VecType *> &X,
                           const std::vector<VecType *> &Y) const
{
  if (X.empty())
  {
    return;
  }
  if (!factors)
  {
    factorization->SolveInternal(X, Y);
    return;
  }
  std::vector<T> u;
  factorization->Reduce(X, u);
  factors->Solve(u, static_cast<int>(X.size()));
  factorization->Expand(u, Y);
}

template class RegionSolver<double>;
template class RegionSolver<std::complex<double>>;

}  // namespace palace
