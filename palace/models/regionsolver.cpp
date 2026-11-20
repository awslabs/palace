// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "regionsolver.hpp"

#include "linalg/lowertriangle.hpp"
#include "linalg/mumpsschur.hpp"
#include "linalg/rap.hpp"
#include "models/drivensubstructure.hpp"
#include "models/environmentoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

RegionSolver::RegionSolver(const IoData &iodata, SpaceOperator &space_op,
                           EnvironmentOperator &env)
  : space_op(space_op), env(env),
    type(SchurSolverType(iodata.solver.substructuring->direct_solver)),
    ordering(iodata.solver.linear.sym_factorization)
{
}

RegionSolver::~RegionSolver() = default;

void RegionSolver::Setup(const mfem::HypreParMatrix *Ar, const mfem::HypreParMatrix *Ai)
{
  // The system: the true DOFs which are not eliminated (Dirichlet or pinned), numbered by
  // rank.
  const auto &fes = space_op.GetNDSpace().Get();
  MPI_Comm comm = fes.GetComm();
  const int nt = fes.GetTrueVSize();
  const HYPRE_BigInt tstart = fes.GetMyTDofOffset();
  const auto &dbc = space_op.GetNDDbcTDofLists().back();
  eliminated.assign(dbc.HostRead(), dbc.HostRead() + dbc.Size());
  std::vector<HYPRE_BigInt> sys(nt, 0);
  for (int i : eliminated)
  {
    sys[i] = -1;
  }
  for (int i = 0; i < nt; i++)
  {
    if (sys[i] == 0)
    {
      rows.push_back(i);
    }
  }
  HYPRE_BigInt nloc = static_cast<HYPRE_BigInt>(rows.size()), off = 0, n_sys = 0;
  MPI_Exscan(&nloc, &off, 1, HYPRE_MPI_BIG_INT, MPI_SUM, comm);
  off = Mpi::Root(comm) ? 0 : off;
  MPI_Allreduce(&nloc, &n_sys, 1, HYPRE_MPI_BIG_INT, MPI_SUM, comm);
  for (std::size_t q = 0; q < rows.size(); q++)
  {
    sys[rows[q]] = off + static_cast<HYPRE_BigInt>(q);
  }

  // The Schur variables: Γ, in interface order.
  const auto &gamma = env.InterfaceTrueDofs();
  std::vector<HYPRE_BigInt> mine(gamma.size());
  for (std::size_t a = 0; a < gamma.size(); a++)
  {
    mine[a] = sys[gamma[a]];
    MFEM_ASSERT(mine[a] >= 0, "An interface DOF of the region is eliminated!");
  }
  std::vector<HYPRE_BigInt> schur_vars(env.Grid().Size());
  MPI_Allgatherv(mine.data(), static_cast<int>(mine.size()), HYPRE_MPI_BIG_INT,
                 schur_vars.data(), env.InterfaceCounts().data(),
                 env.InterfaceOffsets().data(), HYPRE_MPI_BIG_INT, comm);

  // The lower triangle of the matrix on the system, from the pattern of its local rows.
  std::vector<int> irn;
  std::vector<std::vector<double>> vals;
  LowerTrianglePattern({Ar, Ai}, irn, jcn, vals, row_ptr);
  const auto col_sys = OwnerValues(comm, tstart, sys, jcn);
  SchurFactorization<std::complex<double>>::Coo coo;
  for (int i = 0; i < nt; i++)
  {
    for (int q = row_ptr[i]; sys[i] >= 0 && q < row_ptr[i + 1]; q++)
    {
      if (col_sys[q] >= 0)
      {
        kept.push_back(q);
        coo.irn.push_back(static_cast<int>(sys[i] + 1));
        coo.jcn.push_back(static_cast<int>(col_sys[q] + 1));
        coo.val.emplace_back(vals[0].empty() ? 0.0 : vals[0][q],
                             vals[1].empty() ? 0.0 : vals[1][q]);
      }
    }
  }
  factorization = MakeSchurFactorization(type, ordering, comm, n_sys, std::move(coo),
                                         schur_vars, rows, &env.Grid());
  lu = std::make_unique<BlockCyclicLU>(env.Grid());
}

void RegionSolver::SetOperator(const ComplexOperator &op)
{
  height = op.Height();
  width = op.Width();
  auto assembled = [](const Operator *part) -> const mfem::HypreParMatrix *
  { return part ? &dynamic_cast<const ParOperator &>(*part).ParallelAssemble() : nullptr; };
  const auto *Ar = assembled(op.Real()), *Ai = assembled(op.Imag());
  if (!factorization)
  {
    Setup(Ar, Ai);
  }
  else
  {
    // The new values in the same pattern.
    std::vector<std::complex<double>> val(jcn.size(), 0.0);
    if (Ar)
    {
      AddToPattern(*Ar, std::complex<double>(1.0, 0.0), row_ptr, jcn, val);
    }
    if (Ai)
    {
      AddToPattern(*Ai, std::complex<double>(0.0, 1.0), row_ptr, jcn, val);
    }
    auto &v = factorization->Values();
    for (std::size_t e = 0; e < kept.size(); e++)
    {
      v[e] = val[kept[e]];
    }
    factorization->Refactor();
  }

  // S_R + S_E(ω) on the grid, at the frequency of the matrix, factored.
#if defined(MFEM_USE_MUMPS)
  if (auto *mumps =
          dynamic_cast<MumpsSchurSolverT<std::complex<double>> *>(factorization.get()))
  {
    C = &mumps->SchurBlocks();
  }
  else
#endif
  {
    dense = env.Grid().ScatterFromRoot(factorization->ReleaseSchur());
    C = &dense;
  }
  env.AddSchurBlocks(env.PreconditionerFrequency(), *C);
  lu->Factor(*C);
}

void RegionSolver::Mult(const ComplexVector &x, ComplexVector &y) const
{
  std::vector<std::complex<double>> u;
  factorization->Reduce({&x}, u);
  lu->Solve(u, 1);
  factorization->Expand(u, {&y});
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
    return std::make_unique<ComplexKspSolver>(
        iodata.solver.linear, iodata.problem.verbose, space_op.GetComm(),
        std::make_unique<RegionSolver>(iodata, space_op, *env));
  }
  return std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                            &space_op.GetH1Spaces());
}

}  // namespace palace
