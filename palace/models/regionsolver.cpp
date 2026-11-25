// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "regionsolver.hpp"

#include <algorithm>
#include <complex>
#include <type_traits>
#include <utility>
#include "linalg/lowertriangle.hpp"
#include "linalg/mumpsschur.hpp"
#include "linalg/rap.hpp"
#include "linalg/schurfactory.hpp"
#include "models/environmentoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

template <typename T>
RegionSolver<T>::RegionSolver(const BlockCyclicGrid &grid, LinearSolver type,
                              SymbolicFactorization ordering, HYPRE_BigInt n_glob,
                              int n_loc, typename SchurFactorization<T>::Coo &&A,
                              const std::vector<HYPRE_BigInt> &gamma, std::vector<int> rows,
                              const std::function<void(std::vector<T> &)> &add_interface,
                              bool refactor)
  : grid(grid)
{
  MFEM_VERIFY(static_cast<int>(gamma.size()) == grid.Size(),
              "The interface of a region solver is not that of its grid!");
  const bool has_gamma = (grid.Size() > 0);
  factorization = MakeSchurFactorization<T>(
      type, ordering, grid.GetComm(), n_glob, n_loc, std::move(A), gamma, std::move(rows),
      std::is_same_v<T, double>, refactor, has_gamma ? &grid : nullptr);
  if (has_gamma)
  {
    FactorInterface(add_interface);
  }
}

template <typename T>
void RegionSolver<T>::Refactor(const std::function<void(std::vector<T> &)> &add_interface)
{
  factorization->Refactor();
  if (grid.Size() > 0)
  {
    FactorInterface(add_interface);
  }
}

template <typename T>
void RegionSolver<T>::FactorInterface(
    const std::function<void(std::vector<T> &)> &add_interface)
{
  std::vector<T> *C = &dense;
#if defined(MFEM_USE_MUMPS)
  if (auto *mumps = dynamic_cast<MumpsSchurSolverT<T> *>(factorization.get()))
  {
    C = &mumps->SchurBlocks();
  }
  else
#endif
  {
    dense = factorization->ReleaseSchurBlocks(grid);
  }
  add_interface(*C);
  if (!factors)
  {
    factors = std::make_unique<BlockCyclicFactorization<T>>(grid);
  }
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

DrivenRegionSolver::DrivenRegionSolver(const IoData &iodata, SpaceOperator &space_op,
                                       EnvironmentOperator &env)
  : space_op(space_op), env(env),
    type(SchurSolverFor(iodata.solver.substructuring->direct_solver)),
    ordering(iodata.solver.linear.sym_factorization)
{
}

DrivenRegionSolver::~DrivenRegionSolver() = default;

void DrivenRegionSolver::SetOperator(const ComplexOperator &op)
{
  height = op.Height();
  width = op.Width();
  auto assembled = [](const Operator *part) -> const mfem::HypreParMatrix *
  { return part ? &dynamic_cast<const ParOperator &>(*part).ParallelAssemble() : nullptr; };
  const auto *Ar = assembled(op.Real()), *Ai = assembled(op.Imag());
  // S_E(ω) on the grid, at the frequency of the matrix.
  auto add = [this](std::vector<std::complex<double>> &C)
  { env.AddSchurBlocks(env.PreconditionerFrequency(), C); };
  if (region)
  {
    // The new values in the same pattern.
    auto &val = region->Values();
    std::fill(val.begin(), val.end(), 0.0);
    if (Ar)
    {
      side->AddPart(*Ar, std::complex<double>(1.0, 0.0), val);
    }
    if (Ai)
    {
      side->AddPart(*Ai, std::complex<double>(0.0, 1.0), val);
    }
    region->Refactor(add);
    return;
  }

  // The system: the true DOFs which are not eliminated (Dirichlet or pinned), with Γ in
  // interface order as the Schur variables.
  const auto &fes = space_op.GetNDSpace().Get();
  MPI_Comm comm = fes.GetComm();
  const int nt = fes.GetTrueVSize();
  const auto &dbc = space_op.GetEliminatedTrueDofs();
  eliminated.assign(dbc.HostRead(), dbc.HostRead() + dbc.Size());
  std::vector<char> keep(nt, 1);
  for (int i : eliminated)
  {
    keep[i] = 0;
  }
  side =
      std::make_unique<SideSystem>(std::vector<const mfem::HypreParMatrix *>{Ar, Ai}, keep);
  const HYPRE_BigInt tstart = fes.GetMyTDofOffset();
  const auto &gamma = env.InterfaceTrueDofs();
  std::vector<HYPRE_BigInt> mine(gamma.size()), tdofs(env.Grid().Size());
  for (std::size_t a = 0; a < gamma.size(); a++)
  {
    mine[a] = tstart + gamma[a];
  }
  MPI_Allgatherv(mine.data(), static_cast<int>(mine.size()), HYPRE_MPI_BIG_INT,
                 tdofs.data(), env.InterfaceCounts().data(), env.InterfaceOffsets().data(),
                 HYPRE_MPI_BIG_INT, comm);
  auto val = side->Values(std::vector<std::complex<double>>{
      std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 1.0)});
  const int n_loc = static_cast<int>(side->rows.size());
  region = std::make_unique<RegionSolver<std::complex<double>>>(
      env.Grid(), type, ordering, side->n_glob, n_loc,
      SchurFactorization<std::complex<double>>::Coo{std::move(side->irn),
                                                    std::move(side->jcn), std::move(val)},
      side->Indices(tdofs), side->rows, add, true);
  Mpi::Print(comm, " Region solver: {}, dense LU on a {:d} x {:d} process grid\n",
             SchurFactorizationName(type), env.Grid().NumProcessRows(),
             env.Grid().NumProcessCols());
}

void DrivenRegionSolver::Mult(const ComplexVector &x, ComplexVector &y) const
{
  region->Mult({&x}, {&y});
  // The eliminated DOFs have identity rows.
  const double *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
  double *yr = y.Real().HostReadWrite(), *yi = y.Imag().HostReadWrite();
  for (int i : eliminated)
  {
    yr[i] = xr[i];
    yi[i] = xi[i];
  }
}

std::unique_ptr<ComplexKspSolver> MakeDrivenKspSolver(const IoData &iodata,
                                                      SpaceOperator &space_op)
{
  if (auto *env = space_op.GetEnvironment())
  {
    // FGMRES, whose solution combines the preconditioned vectors of its residual: near a
    // resonance of the region, where the region solver is nearly singular, GMRES applies it
    // to a combination of its basis, which loses digits.
    config::LinearSolverData linear = iodata.solver.linear;
    linear.krylov_solver = KrylovSolver::FGMRES;
    linear.pc_side = PreconditionerSide::DEFAULT;
    return std::make_unique<ComplexKspSolver>(
        linear, iodata.problem.verbose, space_op.GetComm(),
        std::make_unique<DrivenRegionSolver>(iodata, space_op, *env));
  }
  return std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                            &space_op.GetH1Spaces());
}

}  // namespace palace
