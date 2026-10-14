// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "strumpackcomplex.hpp"

#if defined(MFEM_USE_STRUMPACK)

#include <algorithm>
#include <utility>
#include <StrumpackSparseSolverMPIDist.hpp>
#include <StrumpackSparseSolverMixedPrecisionMPIDist.hpp>
#include "linalg/rap.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Options as for StrumpackSolver: the factorization is used as a preconditioner or direct
// solver, without compression.
template <typename Options>
void ConfigureOptions(Options &opts, SymbolicFactorization reorder, int print)
{
  opts.set_verbose(print > 1);
  opts.set_Krylov_solver(strumpack::KrylovSolver::DIRECT);
  opts.set_matching(strumpack::MatchingJob::NONE);
  switch (reorder)
  {
    case SymbolicFactorization::PARMETIS:
      opts.set_reordering_method(strumpack::ReorderingStrategy::PARMETIS);
      break;
    case SymbolicFactorization::SCOTCH:
      opts.set_reordering_method(strumpack::ReorderingStrategy::SCOTCH);
      break;
    case SymbolicFactorization::PTSCOTCH:
      opts.set_reordering_method(strumpack::ReorderingStrategy::PTSCOTCH);
      break;
    case SymbolicFactorization::AMD:
      opts.set_reordering_method(strumpack::ReorderingStrategy::AMD);
      break;
    case SymbolicFactorization::RCM:
      opts.set_reordering_method(strumpack::ReorderingStrategy::RCM);
      break;
    case SymbolicFactorization::METIS:
    case SymbolicFactorization::PORD:
    case SymbolicFactorization::DEFAULT:
      opts.set_reordering_method(strumpack::ReorderingStrategy::METIS);
      break;
  }
  opts.set_compression(strumpack::CompressionType::NONE);
  opts.disable_gpu();
}

const mfem::HypreParMatrix *GetAssembled(const Operator *op)
{
  if (!op)
  {
    return nullptr;
  }
  if (const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(op))
  {
    return hA;
  }
  const auto *PtAP = dynamic_cast<const ParOperator *>(op);
  MFEM_VERIFY(PtAP, "StrumpackComplexSolver requires HypreParMatrix or ParOperator parts!");
  return &PtAP->ParallelAssemble();
}

}  // namespace

StrumpackComplexSolver::StrumpackComplexSolver(MPI_Comm comm, SymbolicFactorization reorder,
                                               bool single, int refine_its, int print)
  : Solver<ComplexOperator>(), comm(comm), reorder(reorder), single(single),
    refine_its(refine_its), print(print)
{
}

StrumpackComplexSolver::~StrumpackComplexSolver() = default;

void StrumpackComplexSolver::SetOperator(const ComplexOperator &op)
{
  const auto *hAr = GetAssembled(op.Real());
  const auto *hAi = GetAssembled(op.Imag());
  MFEM_VERIFY(hAr || hAi, "Empty ComplexOperator for StrumpackComplexSolver!");
  const auto &hA = hAr ? *hAr : *hAi;

  // Local rows of the real and imaginary parts, with global column indices, merged into a
  // complex-valued matrix (the sparsity pattern is assumed to be symmetric).
  struct LocalCSR
  {
    hypre_CSRMatrix *csr = nullptr;
    ~LocalCSR()
    {
      if (csr)
      {
        hypre_CSRMatrixDestroy(csr);
      }
    }
  } re, im;
  auto Merge = [](const mfem::HypreParMatrix *A)
  {
    hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(
        (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(*A));
    hypre_CSRMatrixMigrate(csr, HYPRE_MEMORY_HOST);
    return csr;
  };
  if (hAr)
  {
    re.csr = Merge(hAr);
  }
  if (hAi)
  {
    im.csr = Merge(hAi);
  }
  const int n_loc = hA.Height();
  std::vector<int> I(n_loc + 1, 0), J;
  std::vector<std::complex<double>> V;
  std::vector<std::pair<HYPRE_BigInt, std::complex<double>>> row;
  for (int i = 0; i < n_loc; i++)
  {
    row.clear();
    for (auto [part, s] : {std::pair{&re, std::complex<double>(1.0, 0.0)},
                           std::pair{&im, std::complex<double>(0.0, 1.0)}})
    {
      if (!part->csr)
      {
        continue;
      }
      const HYPRE_Int *Ip = hypre_CSRMatrixI(part->csr);
      const HYPRE_BigInt *Jp = hypre_CSRMatrixBigJ(part->csr);
      const double *Dp = hypre_CSRMatrixData(part->csr);
      for (HYPRE_Int k = Ip[i]; k < Ip[i + 1]; k++)
      {
        row.emplace_back(Jp[k], s * Dp[k]);
      }
    }
    std::sort(row.begin(), row.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
    for (std::size_t k = 0; k < row.size(); k++)
    {
      if (!J.empty() && static_cast<int>(J.size()) > I[i] &&
          J.back() == static_cast<int>(row[k].first))
      {
        V.back() += row[k].second;
      }
      else
      {
        J.push_back(static_cast<int>(row[k].first));
        V.push_back(row[k].second);
      }
    }
    I[i + 1] = static_cast<int>(J.size());
  }

  // Row distribution (assumed partition of the rows over the processes).
  const int n_proc = Mpi::Size(comm);
  std::vector<int> dist(n_proc + 1, 0);
  int first_row = static_cast<int>(hA.GetRowStarts()[0]);
  MPI_Allgather(&first_row, 1, MPI_INT, dist.data(), 1, MPI_INT, comm);
  dist[n_proc] = static_cast<int>(hA.GetGlobalNumRows());

  // The solver copies the matrix: release the parallel assembled parts (as
  // MfemWrapperSolver does for the sparse direct solvers).
  for (const auto *part : {op.Real(), op.Imag()})
  {
    if (const auto *PtAP = dynamic_cast<const ParOperator *>(part))
    {
      PtAP->StealParallelAssemble();
    }
  }

  // Reuse the reordering for a matrix with the same sparsity pattern on all processes.
  int same = (solver_d || solver_s) && I == I_prev && J == J_prev;
  Mpi::GlobalMin(1, &same, comm);
  if (!single)
  {
    if (!same)
    {
      solver_d =
          std::make_unique<strumpack::SparseSolverMPIDist<std::complex<double>, int>>(
              comm, print > 1);
      ConfigureOptions(solver_d->options(), reorder, print);
      solver_d->set_distributed_csr_matrix(n_loc, I.data(), J.data(), V.data(), dist.data(),
                                           true);
      MFEM_VERIFY(solver_d->reorder() == strumpack::ReturnCode::SUCCESS,
                  "STRUMPACK reordering failed!");
    }
    else
    {
      solver_d->update_matrix_values(n_loc, I.data(), J.data(), V.data(), dist.data(),
                                     true);
    }
    MFEM_VERIFY(solver_d->factor() == strumpack::ReturnCode::SUCCESS,
                "STRUMPACK factorization failed!");
  }
  else
  {
    strumpack::CSRMatrixMPI<std::complex<double>, int> A(
        n_loc, I.data(), J.data(), V.data(), dist.data(), strumpack::MPIComm(comm), true);
    if (!same)
    {
      solver_s = std::make_unique<strumpack::SparseSolverMixedPrecisionMPIDist<
          std::complex<float>, std::complex<double>, int>>(comm, print > 1);
      ConfigureOptions(solver_s->solver().options(), reorder, print);
      ConfigureOptions(solver_s->options(), reorder, print);
      if (refine_its > 0)
      {
        solver_s->options().set_Krylov_solver(strumpack::KrylovSolver::REFINE);
        solver_s->options().set_maxit(refine_its + 1);
        solver_s->options().set_rel_tol(1.0e-12);
        solver_s->options().set_abs_tol(0.0);
      }
      solver_s->set_matrix(A);
      MFEM_VERIFY(solver_s->reorder() == strumpack::ReturnCode::SUCCESS,
                  "STRUMPACK reordering failed!");
    }
    else
    {
      solver_s->update_matrix_values(A);
    }
    MFEM_VERIFY(solver_s->factor() == strumpack::ReturnCode::SUCCESS,
                "STRUMPACK factorization failed!");
  }
  I_prev = std::move(I);
  J_prev = std::move(J);
  b.resize(n_loc);
  x.resize(n_loc);
  this->height = op.Height();
  this->width = op.Width();
}

void StrumpackComplexSolver::Mult(const ComplexVector &xv, ComplexVector &yv) const
{
  const int n = static_cast<int>(b.size());
  const auto *xr = xv.Real().HostRead(), *xi = xv.Imag().HostRead();
  for (int i = 0; i < n; i++)
  {
    b[i] = {xr[i], xi[i]};
  }
  if (solver_d)
  {
    solver_d->solve(b.data(), x.data(), false);
  }
  else
  {
    solver_s->solve(b.data(), x.data(), false);
  }
  auto *yr = yv.Real().HostWrite(), *yi = yv.Imag().HostWrite();
  for (int i = 0; i < n; i++)
  {
    yr[i] = x[i].real();
    yi[i] = x[i].imag();
  }
}

}  // namespace palace

#endif
