// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "strumpackcomplex.hpp"

#if defined(MFEM_USE_STRUMPACK)

#include <StrumpackSparseSolverMPIDist.hpp>
#include "linalg/complexcsr.hpp"
#include "linalg/strumpack.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Options as for StrumpackSolver: the factorization is used as a preconditioner or direct
// solver.
void ConfigureOptions(strumpack::SPOptions<std::complex<double>> &opts,
                      SymbolicFactorization reorder, int print)
{
  opts.set_verbose(print > 1);
  opts.set_Krylov_solver(strumpack::KrylovSolver::DIRECT);
  opts.set_matching(strumpack::MatchingJob::NONE);
  ConfigureStrumpackReordering(opts, reorder);
  opts.disable_gpu();
}

}  // namespace

StrumpackComplexSolver::StrumpackComplexSolver(MPI_Comm comm, SymbolicFactorization reorder,
                                               SparseCompression compression, double lr_tol,
                                               int butterfly_l, int lossy_prec,
                                               bool drop_small_entries, bool reorder_reuse,
                                               int print)
  : Solver<ComplexOperator>(), comm(comm), reorder(reorder), compression(compression),
    lr_tol(lr_tol), butterfly_l(butterfly_l), lossy_prec(lossy_prec),
    drop_small_entries(drop_small_entries), reorder_reuse(reorder_reuse), print(print)
{
}

StrumpackComplexSolver::~StrumpackComplexSolver() = default;

void StrumpackComplexSolver::SetOperator(const ComplexOperator &op)
{
  bool same;
  {
    // Local rows of the complex-valued matrix (freed before the factorization), and row
    // distribution (assumed partition of the rows over the processes). The sparsity pattern
    // is not declared symmetric, as for StrumpackSolver: hypre's parallel matrix products
    // can store explicit zeros without their transposes on small partitions, and STRUMPACK
    // symmetrizes the pattern itself.
    ComplexLocalCSR csr(op, drop_small_entries);

    // STRUMPACK copies the matrix: release the parallel assembled parts (as
    // MfemWrapperSolver does for the sparse direct solvers).
    ReleaseParallelAssemble(op);
    const int n_proc = Mpi::Size(comm);
    std::vector<HYPRE_BigInt> dist(n_proc + 1);
    Mpi::Allgather(1, &csr.first_row, dist.data(), comm);
    dist[n_proc] = csr.n_glob;

    // Reuse the reordering for a matrix with the same sparsity pattern on all processes
    // (when small entries are dropped, the sparsity pattern can change with the values).
    same = csr.SamePattern(comm, pattern) && solver && reorder_reuse;
    if (!same)
    {
      solver = std::make_unique<
          strumpack::SparseSolverMPIDist<std::complex<double>, HYPRE_BigInt>>(comm,
                                                                              print > 1);
      ConfigureOptions(solver->options(), reorder, print);
      ConfigureStrumpackCompression(solver->options(), compression, lr_tol, butterfly_l,
                                    lossy_prec);
      EnableStrumpackSharedPermutation(solver->options());
      solver->set_distributed_csr_matrix(csr.n_loc, csr.I.data(), csr.J.data(),
                                         csr.V.data(), dist.data(), false);
    }
    else
    {
      solver->update_matrix_values(csr.n_loc, csr.I.data(), csr.J.data(), csr.V.data(),
                                   dist.data(), false);
    }
    b.resize(csr.n_loc);
    x.resize(csr.n_loc);
  }
  if (!same)
  {
    const auto ret = solver->reorder();
    MFEM_VERIFY(ret == strumpack::ReturnCode::SUCCESS,
                "STRUMPACK reordering failed with return code " << static_cast<int>(ret)
                                                                << "!");
  }
  const auto ret = solver->factor();
  MFEM_VERIFY(ret == strumpack::ReturnCode::SUCCESS,
              "STRUMPACK factorization failed with return code " << static_cast<int>(ret)
                                                                 << "!");
  this->height = op.Height();
  this->width = op.Width();
}

void StrumpackComplexSolver::Mult(const ComplexVector &xv, ComplexVector &yv) const
{
  xv.Get(b.data(), static_cast<int>(b.size()), false);
  const auto ret = solver->solve(b.data(), x.data(), false);
  MFEM_VERIFY(ret == strumpack::ReturnCode::SUCCESS,
              "STRUMPACK solve failed with return code " << static_cast<int>(ret) << "!");
  yv.Set(x.data(), static_cast<int>(x.size()), false);
}

}  // namespace palace

#endif
