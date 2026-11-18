// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "strumpackschur.hpp"

#if defined(PALACE_WITH_STRUMPACK_SCHUR)

#include <algorithm>
#include <climits>
#include "utils/communication.hpp"

namespace palace
{

StrumpackSchurSolver::StrumpackSchurSolver(const mfem::HypreParMatrix &A,
                                           const std::vector<HYPRE_BigInt> &schur_vars,
                                           SymbolicFactorization ordering)
  : comm(A.GetComm()), n_loc(A.Height()), n_schur(static_cast<int>(schur_vars.size())),
    solver(comm, false), rows(comm, A.GetRowStarts()[0], A.Height(), schur_vars)
{
  MFEM_VERIFY(A.GetGlobalNumRows() <= INT_MAX,
              "STRUMPACK Schur complement needs fewer than 2^31 unknowns!");
  // The Schur complement needs neither matching nor the GPU (and a sequential graph
  // reordering without compression).
  solver.options().set_matching(strumpack::MatchingJob::NONE);
  solver.options().disable_gpu();
  switch (ordering)
  {
    case SymbolicFactorization::METIS:
    case SymbolicFactorization::PARMETIS:
      solver.options().set_reordering_method(strumpack::ReorderingStrategy::METIS);
      break;
    case SymbolicFactorization::SCOTCH:
    case SymbolicFactorization::PTSCOTCH:
      solver.options().set_reordering_method(strumpack::ReorderingStrategy::SCOTCH);
      break;
    case SymbolicFactorization::AMD:
      solver.options().set_reordering_method(strumpack::ReorderingStrategy::AMD);
      break;
    case SymbolicFactorization::RCM:
      solver.options().set_reordering_method(strumpack::ReorderingStrategy::RCM);
      break;
    case SymbolicFactorization::PORD:
    case SymbolicFactorization::DEFAULT:
      break;
  }

  // The local rows with global columns: the block-row distributed CSR STRUMPACK takes (it
  // copies the entries).
  auto *parcsr = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(A);
  A.HostRead();
  hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(parcsr);
  A.HypreRead();
  const int nnz = csr->i[n_loc];
#if MFEM_HYPRE_VERSION >= 21600
  const std::vector<int> col(csr->big_j, csr->big_j + nnz);
#else
  const std::vector<int> col(csr->j, csr->j + nnz);
#endif
  const std::vector<int> ptr(csr->i, csr->i + n_loc + 1);
  const int nranks = Mpi::Size(comm);
  std::vector<int> dist(nranks + 1, 0);
  MPI_Allgather(&n_loc, 1, MPI_INT, dist.data() + 1, 1, MPI_INT, comm);
  for (int r = 0; r < nranks; r++)
  {
    dist[r + 1] += dist[r];
  }
  solver.set_distributed_csr_matrix(n_loc, ptr.data(), col.data(), csr->data, dist.data(),
                                    true);
  hypre_CSRMatrixDestroy(csr);

  const std::vector<int> vars(schur_vars.begin(), schur_vars.end());
  MFEM_VERIFY(solver.set_schur_variables(n_schur, vars.data()) ==
                  strumpack::ReturnCode::SUCCESS,
              "STRUMPACK rejected the Schur variables!");
  MFEM_VERIFY(solver.factor() == strumpack::ReturnCode::SUCCESS,
              "STRUMPACK partial factorization failed!");
  strumpack::DenseMatrix<double> S;
  MFEM_VERIFY(solver.get_schur_complement(S) == strumpack::ReturnCode::SUCCESS,
              "STRUMPACK Schur complement failed!");
  if (Mpi::Root(comm))
  {
    schur.resize(static_cast<std::size_t>(n_schur) * n_schur);
    for (int j = 0; j < n_schur; j++)
    {
      for (int i = 0; i < n_schur; i++)
      {
        schur[i + static_cast<std::size_t>(n_schur) * j] = S(i, j);
      }
    }
  }
}

std::vector<double> StrumpackSchurSolver::ReleaseSchurRows()
{
  std::vector<double> S = rows.ScatterSchur(schur);
  std::vector<double>().swap(schur);
  return S;
}

void StrumpackSchurSolver::SolveInternal(const std::vector<const mfem::Vector *> &X,
                                         const std::vector<mfem::Vector *> &Y)
{
  // The expansion of zero values of the Schur variables: x_I = A_II^-1 b_I, x_S = 0.
  const int nb = static_cast<int>(X.size()), ld = std::max(n_loc, 1);
  std::vector<double> b(static_cast<std::size_t>(ld) * nb), x(b.size(), 0.0);
  for (int k = 0; k < nb; k++)
  {
    const double *xk = X[k]->HostRead();
    std::copy(xk, xk + n_loc, b.begin() + static_cast<std::size_t>(ld) * k);
  }
  MFEM_VERIFY(solver.schur_expand(nb, b.data(), ld, x.data(), ld) ==
                  strumpack::ReturnCode::SUCCESS,
              "STRUMPACK internal solve failed!");
  for (int k = 0; k < nb; k++)
  {
    Y[k]->SetSize(n_loc);
    std::copy(x.begin() + static_cast<std::size_t>(ld) * k,
              x.begin() + static_cast<std::size_t>(ld) * k + n_loc, Y[k]->HostWrite());
  }
}

void StrumpackSchurSolver::Condense(const std::vector<const mfem::Vector *> &X,
                                    std::vector<double> &g)
{
  // The condensations in the rows of the Schur variables (by the forward elimination),
  // gathered on rank 0.
  const int nb = static_cast<int>(X.size()), ld = std::max(n_loc, 1);
  g.assign(Mpi::Root(comm) ? static_cast<std::size_t>(n_schur) * nb : 0, 0.0);
  if (nb == 0)
  {
    return;
  }
  std::vector<double> b(static_cast<std::size_t>(ld) * nb), c(b.size(), 0.0);
  for (int k = 0; k < nb; k++)
  {
    const double *xk = X[k]->HostRead();
    std::copy(xk, xk + n_loc, b.begin() + static_cast<std::size_t>(ld) * k);
  }
  MFEM_VERIFY(solver.schur_condense(nb, b.data(), ld, c.data(), ld) ==
                  strumpack::ReturnCode::SUCCESS,
              "STRUMPACK condensation failed!");
  const std::vector<int> &row = rows.Rows();
  const int m = static_cast<int>(row.size());
  std::vector<double> loc(static_cast<std::size_t>(m) * nb);
  for (int k = 0; k < nb; k++)
  {
    for (int r = 0; r < m; r++)
    {
      loc[static_cast<std::size_t>(k) * m + r] =
          c[static_cast<std::size_t>(ld) * k + row[r]];
    }
  }
  rows.Gather(loc, nb, g.data());
}

}  // namespace palace

#endif
