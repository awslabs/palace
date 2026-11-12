// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "strumpackcomplex.hpp"

#if defined(MFEM_USE_STRUMPACK)

#include <utility>
#include <StrumpackSparseSolverMPIDist.hpp>
#include "linalg/complexcsr.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Options as for StrumpackSolver without compression: the factorization is used as a
// preconditioner or direct solver.
void ConfigureOptions(strumpack::SPOptions<std::complex<double>> &opts,
                      SymbolicFactorization reorder, int print)
{
  opts.set_verbose(print > 1);
  opts.set_Krylov_solver(strumpack::KrylovSolver::DIRECT);
  opts.set_matching(strumpack::MatchingJob::NONE);
  switch (reorder)
  {
    case SymbolicFactorization::METIS:
      opts.set_reordering_method(strumpack::ReorderingStrategy::METIS);
      break;
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
    case SymbolicFactorization::PORD:
    case SymbolicFactorization::DEFAULT:
      // Should have good default.
      break;
  }
  opts.set_compression(strumpack::CompressionType::NONE);
  opts.disable_gpu();
}

}  // namespace

StrumpackComplexSolver::StrumpackComplexSolver(MPI_Comm comm, SymbolicFactorization reorder,
                                               bool reorder_reuse, int print)
  : Solver<ComplexOperator>(), comm(comm), reorder(reorder), reorder_reuse(reorder_reuse),
    print(print)
{
}

StrumpackComplexSolver::~StrumpackComplexSolver() = default;

void StrumpackComplexSolver::SetOperator(const ComplexOperator &op)
{
  // Local rows of the complex-valued matrix, and row distribution (assumed partition of the
  // rows over the processes). The sparsity pattern is not declared symmetric, as for
  // StrumpackSolver: hypre's parallel matrix products can store explicit zeros without
  // their transposes on small partitions, and STRUMPACK symmetrizes the pattern itself.
  ComplexLocalCSR csr(op);
  const HYPRE_BigInt n_loc = csr.n_loc;
  std::vector<HYPRE_BigInt> I(csr.I.begin(), csr.I.end()), J = std::move(csr.J);
  std::vector<std::complex<double>> V = std::move(csr.V);
  const int n_proc = Mpi::Size(comm);
  std::vector<HYPRE_BigInt> dist(n_proc + 1);
  Mpi::Allgather(1, &csr.first_row, dist.data(), comm);
  dist[n_proc] = csr.n_glob;

  // Reuse the reordering for a matrix with the same sparsity pattern on all processes.
  int same = solver && reorder_reuse && I == I_prev && J == J_prev;
  Mpi::GlobalMin(1, &same, comm);
  if (!same)
  {
    solver = std::make_unique<
        strumpack::SparseSolverMPIDist<std::complex<double>, HYPRE_BigInt>>(comm,
                                                                            print > 1);
    ConfigureOptions(solver->options(), reorder, print);
    solver->set_distributed_csr_matrix(n_loc, I.data(), J.data(), V.data(), dist.data(),
                                       false);
    const auto ret = solver->reorder();
    MFEM_VERIFY(ret == strumpack::ReturnCode::SUCCESS,
                "STRUMPACK reordering failed with return code " << static_cast<int>(ret)
                                                                << "!");
  }
  else
  {
    solver->update_matrix_values(n_loc, I.data(), J.data(), V.data(), dist.data(), false);
  }
  const auto ret = solver->factor();
  MFEM_VERIFY(ret == strumpack::ReturnCode::SUCCESS,
              "STRUMPACK factorization failed with return code " << static_cast<int>(ret)
                                                                 << "!");
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
  const auto ret = solver->solve(b.data(), x.data(), false);
  MFEM_VERIFY(ret == strumpack::ReturnCode::SUCCESS,
              "STRUMPACK solve failed with return code " << static_cast<int>(ret) << "!");
  auto *yr = yv.Real().HostWrite(), *yi = yv.Imag().HostWrite();
  for (int i = 0; i < n; i++)
  {
    yr[i] = x[i].real();
    yi[i] = x[i].imag();
  }
}

}  // namespace palace

#endif
