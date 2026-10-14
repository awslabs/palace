// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mumpscomplex.hpp"

#if defined(MFEM_USE_MUMPS)

#include <algorithm>
#include <complex>
#include <limits>
#include <vector>
#include "linalg/complexcsr.hpp"
#include "utils/communication.hpp"

// MUMPS complex double precision header.
#include "zmumps_c.h"

#if MFEM_MUMPS_VERSION < 530
#error "MumpsComplexSolver requires MUMPS 5.3.0 or later (distributed right-hand sides)!"
#endif

namespace palace
{

struct MumpsComplexSolver::Data
{
  ZMUMPS_STRUC_C id{};
  bool initialized = false;

  // Local rows, global index of the first local row, local solution entries, and first
  // global row of each process.
  int n_loc = 0, row_start = 0, n_sol = 0;
  std::vector<int> row_starts;

  // Sparsity pattern of the factored matrix (1-based global row and column indices of the
  // local entries), for the reuse of the analysis.
  std::vector<MUMPS_INT> I, J;

  // Distributed right-hand side (in the row distribution of the matrix) and distributed
  // solution (in the distribution of MUMPS).
  std::vector<MUMPS_INT> irhs_loc, isol_loc;
  std::vector<mumps_double_complex> rhs_loc, sol_loc;

  // Redistribution of the solution to the row distribution of the matrix, set up on the
  // first solve after a factorization: for each local solution entry, its position in the
  // send buffer or -1 - (local row) for a local row, and the local rows of the received
  // entries (counts and displacements in doubles).
  bool sol_map = false;
  std::vector<int> send_pos, recv_rows, send_count, send_displ, recv_count, recv_displ;
  std::vector<double> send_buf, recv_buf;

  // Accessors with the indexing of the MUMPS documentation.
  MUMPS_INT &icntl(int i) { return id.icntl[i - 1]; }
  double &cntl(int i) { return id.cntl[i - 1]; }
  MUMPS_INT info(int i) const { return id.info[i - 1]; }
  MUMPS_INT infog(int i) const { return id.infog[i - 1]; }

  void Call(int job)
  {
    id.job = job;
    zmumps_c(&id);
  }

  void Destroy()
  {
    if (initialized)
    {
      Call(-2);
      initialized = false;
    }
  }

  void SetUpSolutionMap(MPI_Comm comm);
};

void MumpsComplexSolver::Data::SetUpSolutionMap(MPI_Comm comm)
{
  const int n_proc = static_cast<int>(row_starts.size()), rank = Mpi::Rank(comm);
  std::vector<int> owner(n_sol), count(n_proc, 0), rcount(n_proc), sdispl(n_proc + 1, 0),
      rdispl(n_proc + 1, 0);
  for (int k = 0; k < n_sol; k++)
  {
    const int g = isol_loc[k] - 1;
    owner[k] = static_cast<int>(std::upper_bound(row_starts.begin(), row_starts.end(), g) -
                                row_starts.begin() - 1);
    if (owner[k] != rank)
    {
      count[owner[k]]++;
    }
  }
  MPI_Alltoall(count.data(), 1, MPI_INT, rcount.data(), 1, MPI_INT, comm);
  for (int p = 0; p < n_proc; p++)
  {
    sdispl[p + 1] = sdispl[p] + count[p];
    rdispl[p + 1] = rdispl[p] + rcount[p];
  }
  send_pos.resize(n_sol);
  std::vector<int> send_rows(sdispl[n_proc]), offset(sdispl.begin(), sdispl.end() - 1);
  for (int k = 0; k < n_sol; k++)
  {
    const int g = isol_loc[k] - 1;
    if (owner[k] == rank)
    {
      send_pos[k] = -1 - (g - row_start);
    }
    else
    {
      send_pos[k] = offset[owner[k]]++;
      send_rows[send_pos[k]] = g;
    }
  }
  recv_rows.resize(rdispl[n_proc]);
  MPI_Alltoallv(send_rows.data(), count.data(), sdispl.data(), MPI_INT, recv_rows.data(),
                rcount.data(), rdispl.data(), MPI_INT, comm);
  for (auto &i : recv_rows)
  {
    i -= row_start;
  }

  // Values are exchanged as pairs of doubles.
  send_count.resize(n_proc);
  send_displ.resize(n_proc);
  recv_count.resize(n_proc);
  recv_displ.resize(n_proc);
  for (int p = 0; p < n_proc; p++)
  {
    send_count[p] = 2 * count[p];
    send_displ[p] = 2 * sdispl[p];
    recv_count[p] = 2 * rcount[p];
    recv_displ[p] = 2 * rdispl[p];
  }
  send_buf.resize(2 * sdispl[n_proc]);
  recv_buf.resize(2 * rdispl[n_proc]);
  sol_map = true;
}

MumpsComplexSolver::MumpsComplexSolver(MPI_Comm comm, MatrixSymmetry sym,
                                       SymbolicFactorization reorder, double blr_tol,
                                       bool reorder_reuse, int print)
  : Solver<ComplexOperator>(), comm(comm), sym(sym), reorder(reorder), blr_tol(blr_tol),
    reorder_reuse(reorder_reuse), print(print), data(std::make_unique<Data>())
{
}

MumpsComplexSolver::~MumpsComplexSolver()
{
  data->Destroy();
}

void MumpsComplexSolver::SetOperator(const ComplexOperator &op)
{
  // Local entries of the complex-valued matrix in coordinate format with 1-based global
  // indices, only of the lower triangle for symmetric matrices (complex symmetric, not
  // Hermitian).
  const bool symmetric = (sym != MatrixSymmetry::UNSYMMETRIC);
  std::vector<MUMPS_INT> I, J;
  std::vector<mumps_double_complex> V;
  HYPRE_BigInt n_glob;
  int n_loc, row_start;
  {
    ComplexLocalCSR csr(op);
    MFEM_VERIFY(csr.n_glob <= std::numeric_limits<MUMPS_INT>::max(),
                "MumpsComplexSolver does not support 64-bit global indices!");
    n_glob = csr.n_glob;
    n_loc = csr.n_loc;
    row_start = static_cast<int>(csr.first_row);
    const std::size_t nnz = csr.J.size();
    I.reserve(symmetric ? nnz / 2 + n_loc : nnz);
    J.reserve(I.capacity());
    V.reserve(I.capacity());
    for (int i = 0; i < n_loc; i++)
    {
      const MUMPS_INT gi = row_start + i + 1;
      for (int k = csr.I[i]; k < csr.I[i + 1]; k++)
      {
        const auto gj = static_cast<MUMPS_INT>(csr.J[k] + 1);
        if (symmetric && gj > gi)
        {
          continue;
        }
        I.push_back(gi);
        J.push_back(gj);
        V.push_back({csr.V[k].real(), csr.V[k].imag()});
      }
    }
  }

  // Analysis, unless the analysis of a previous matrix with the same sparsity pattern is
  // reused.
  auto &d = *data;
  int same = d.initialized && reorder_reuse && d.id.n == static_cast<MUMPS_INT>(n_glob) &&
             I == d.I && J == d.J;
  Mpi::GlobalMin(1, &same, comm);
  if (!same)
  {
    d.Destroy();
    d.id = ZMUMPS_STRUC_C{};
    d.id.sym = symmetric ? 2 : 0;
    d.id.par = 1;
    d.id.comm_fortran = static_cast<MUMPS_INT>(MPI_Comm_c2f(comm));
    d.Call(-1);
    MFEM_VERIFY(d.infog(1) >= 0,
                "MUMPS initialization failed (INFOG(1) = " << d.infog(1) << ")!");
    d.initialized = true;

    // Parameters as for MumpsSolver (mfem::MUMPSSolver), with a distributed matrix,
    // right-hand side, and solution.
    d.icntl(1) = 6;
    d.icntl(2) = 0;
    d.icntl(3) = 6;
    d.icntl(4) = print;
    d.icntl(5) = 0;
    d.icntl(9) = 1;
    d.icntl(10) = 0;
    d.icntl(11) = 0;
    d.icntl(13) = 0;
    d.icntl(14) = 20;
    d.icntl(16) = 0;
    d.icntl(18) = 3;
    d.icntl(19) = 0;
    d.icntl(20) = 10;
    d.icntl(21) = 1;
    d.icntl(22) = 0;
    d.icntl(23) = 0;
    switch (reorder)
    {
      case SymbolicFactorization::METIS:
        d.icntl(28) = 1;
        d.icntl(7) = 5;
        break;
      case SymbolicFactorization::PARMETIS:
        d.icntl(28) = 2;
        d.icntl(29) = 2;
        break;
      case SymbolicFactorization::SCOTCH:
        d.icntl(28) = 1;
        d.icntl(7) = 3;
        break;
      case SymbolicFactorization::PTSCOTCH:
        d.icntl(28) = 2;
        d.icntl(29) = 1;
        break;
      case SymbolicFactorization::PORD:
        d.icntl(28) = 1;
        d.icntl(7) = 4;
        break;
      case SymbolicFactorization::AMD:
      case SymbolicFactorization::RCM:
        d.icntl(28) = 1;
        d.icntl(7) = 0;
        break;
      case SymbolicFactorization::DEFAULT:
        d.icntl(28) = 0;
        d.icntl(7) = 7;
        d.icntl(29) = 0;
        break;
    }
    if (blr_tol > 0.0)
    {
      d.icntl(35) = 1;
      d.cntl(7) = blr_tol;
    }

    d.id.n = static_cast<MUMPS_INT>(n_glob);
    d.id.nnz_loc = static_cast<MUMPS_INT8>(I.size());
    d.id.irn_loc = I.data();
    d.id.jcn_loc = J.data();
    d.id.a_loc = V.data();
    d.Call(1);
    MFEM_VERIFY(d.infog(1) >= 0,
                "MUMPS analysis failed (INFOG(1) = " << d.infog(1) << ")!");
  }
  else
  {
    d.id.irn_loc = I.data();
    d.id.jcn_loc = J.data();
    d.id.a_loc = V.data();
  }

  // Numerical factorization, with more workspace if the estimate is insufficient.
  constexpr int mem_relax_lim = 200;
  while (true)
  {
    d.Call(2);
    if (d.infog(1) >= 0)
    {
      break;
    }
    MFEM_VERIFY(d.infog(1) == -8 || d.infog(1) == -9,
                "MUMPS numerical factorization failed (INFOG(1) = " << d.infog(1) << ")!");
    d.icntl(14) += 20;
    MFEM_VERIFY(d.icntl(14) <= mem_relax_lim,
                "Memory relaxation limit reached for MUMPS factorization!");
    if (print > 0)
    {
      Mpi::Print(comm, " Re-running MUMPS factorization with memory relaxation {:d}\n",
                 d.icntl(14));
    }
  }

  // The values are not needed after the factorization (no iterative refinement), the
  // sparsity pattern is kept for the reuse of the analysis.
  d.I = std::move(I);
  d.J = std::move(J);
  d.id.irn_loc = d.I.data();
  d.id.jcn_loc = d.J.data();
  d.id.a_loc = nullptr;

  // Storage for the distributed right-hand side and solution.
  d.n_loc = n_loc;
  d.row_start = row_start;
  d.row_starts.resize(Mpi::Size(comm));
  MPI_Allgather(&row_start, 1, MPI_INT, d.row_starts.data(), 1, MPI_INT, comm);
  d.n_sol = d.info(23);
  d.irhs_loc.resize(std::max(n_loc, 1));
  for (int i = 0; i < n_loc; i++)
  {
    d.irhs_loc[i] = row_start + i + 1;
  }
  d.rhs_loc.resize(std::max(n_loc, 1));
  d.isol_loc.resize(std::max(d.n_sol, 1));
  d.sol_loc.resize(std::max(d.n_sol, 1));
  d.id.nrhs = 1;
  d.id.nloc_rhs = n_loc;
  d.id.lrhs_loc = std::max(n_loc, 1);
  d.id.irhs_loc = d.irhs_loc.data();
  d.id.rhs_loc = d.rhs_loc.data();
  d.id.lsol_loc = d.n_sol;
  d.id.isol_loc = d.isol_loc.data();
  d.id.sol_loc = d.sol_loc.data();
  d.sol_map = false;
  this->height = op.Height();
  this->width = op.Width();
}

void MumpsComplexSolver::Mult(const ComplexVector &x, ComplexVector &y) const
{
  auto &d = *data;
  const auto *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
  for (int i = 0; i < d.n_loc; i++)
  {
    d.rhs_loc[i] = {xr[i], xi[i]};
  }
  d.Call(3);
  MFEM_VERIFY(d.infog(1) >= 0, "MUMPS solve failed (INFOG(1) = " << d.infog(1) << ")!");

  // Redistribute the solution (the indices of the local solution entries are set by the
  // first solve after a factorization).
  if (!d.sol_map)
  {
    d.SetUpSolutionMap(comm);
  }
  auto *yr = y.Real().HostWrite(), *yi = y.Imag().HostWrite();
  for (int k = 0; k < d.n_sol; k++)
  {
    const int pos = d.send_pos[k];
    if (pos < 0)
    {
      yr[-1 - pos] = d.sol_loc[k].r;
      yi[-1 - pos] = d.sol_loc[k].i;
    }
    else
    {
      d.send_buf[2 * pos] = d.sol_loc[k].r;
      d.send_buf[2 * pos + 1] = d.sol_loc[k].i;
    }
  }
  MPI_Alltoallv(d.send_buf.data(), d.send_count.data(), d.send_displ.data(), MPI_DOUBLE,
                d.recv_buf.data(), d.recv_count.data(), d.recv_displ.data(), MPI_DOUBLE,
                comm);
  for (std::size_t k = 0; k < d.recv_rows.size(); k++)
  {
    yr[d.recv_rows[k]] = d.recv_buf[2 * k];
    yi[d.recv_rows[k]] = d.recv_buf[2 * k + 1];
  }
}

}  // namespace palace

#endif
