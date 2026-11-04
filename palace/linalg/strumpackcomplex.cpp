// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "strumpackcomplex.hpp"

#if defined(MFEM_USE_STRUMPACK)

#include <limits>
#include <utility>
#include <StrumpackSparseSolverMPIDist.hpp>
#include "linalg/complexcsr.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Options as for StrumpackSolver: the factorization is used as a preconditioner or direct
// solver, without compression.
void ConfigureOptions(strumpack::SPOptions<std::complex<double>> &opts,
                      SymbolicFactorization reorder, int print)
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

}  // namespace

StrumpackComplexSolver::StrumpackComplexSolver(MPI_Comm comm, SymbolicFactorization reorder,
                                               int print)
  : Solver<ComplexOperator>(), comm(comm), reorder(reorder), print(print)
{
}

StrumpackComplexSolver::~StrumpackComplexSolver() = default;

void StrumpackComplexSolver::SetOperator(const ComplexOperator &op)
{
  // Local rows of the complex-valued matrix (the sparsity pattern is assumed to be
  // symmetric), and row distribution (assumed partition of the rows over the processes).
  ComplexLocalCSR csr(op);
  const int n_loc = csr.n_loc;
  MFEM_VERIFY(csr.n_glob <= std::numeric_limits<int>::max(),
              "StrumpackComplexSolver does not support 64-bit global indices!");
  std::vector<int> I = std::move(csr.I), J(csr.J.begin(), csr.J.end());
  std::vector<std::complex<double>> V = std::move(csr.V);
  const int n_proc = Mpi::Size(comm);
  std::vector<int> dist(n_proc + 1, 0);
  int first_row = static_cast<int>(csr.first_row);
  MPI_Allgather(&first_row, 1, MPI_INT, dist.data(), 1, MPI_INT, comm);
  dist[n_proc] = static_cast<int>(csr.n_glob);

  // Reuse the reordering for a matrix with the same sparsity pattern on all processes.
  int same = solver && I == I_prev && J == J_prev;
  Mpi::GlobalMin(1, &same, comm);
  if (!same)
  {
    solver = std::make_unique<strumpack::SparseSolverMPIDist<std::complex<double>, int>>(
        comm, print > 1);
    ConfigureOptions(solver->options(), reorder, print);
    solver->set_distributed_csr_matrix(n_loc, I.data(), J.data(), V.data(), dist.data(),
                                       true);
    MFEM_VERIFY(solver->reorder() == strumpack::ReturnCode::SUCCESS,
                "STRUMPACK reordering failed!");
  }
  else
  {
    solver->update_matrix_values(n_loc, I.data(), J.data(), V.data(), dist.data(), true);
  }
  MFEM_VERIFY(solver->factor() == strumpack::ReturnCode::SUCCESS,
              "STRUMPACK factorization failed!");
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
  solver->solve(b.data(), x.data(), false);
  auto *yr = yv.Real().HostWrite(), *yi = yv.Imag().HostWrite();
  for (int i = 0; i < n; i++)
  {
    yr[i] = x[i].real();
    yi[i] = x[i].imag();
  }
}

}  // namespace palace

#endif
