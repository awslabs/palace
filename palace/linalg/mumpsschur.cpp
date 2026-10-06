// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mumpsschur.hpp"

#if defined(MFEM_USE_MUMPS)

#include <algorithm>
#include <cmath>
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Local rows of A as 1-based COO, lower triangle only (symmetric storage: MUMPS sums
// duplicates).
void LowerTriangleCOO(const mfem::HypreParMatrix &A, std::vector<MUMPS_INT> &irn,
                      std::vector<MUMPS_INT> &jcn, std::vector<double> &val)
{
  irn.clear();
  jcn.clear();
  val.clear();
  auto *parcsr = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(A);
  A.HostRead();
  hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(parcsr);
  const HYPRE_Int *Ip = csr->i;
#if MFEM_HYPRE_VERSION >= 21600
  const HYPRE_BigInt *Jp = csr->big_j;
#else
  const HYPRE_Int *Jp = csr->j;
#endif
  const HYPRE_BigInt row0 = parcsr->first_row_index;
  for (int i = 0; i < A.Height(); i++)
  {
    for (HYPRE_Int k = Ip[i]; k < Ip[i + 1]; k++)
    {
      const HYPRE_BigInt ii = row0 + i + 1, jj = static_cast<HYPRE_BigInt>(Jp[k]) + 1;
      if (ii >= jj)
      {
        irn.push_back(static_cast<MUMPS_INT>(ii));
        jcn.push_back(static_cast<MUMPS_INT>(jj));
        val.push_back(csr->data[k]);
      }
    }
  }
  hypre_CSRMatrixDestroy(csr);
}

}  // namespace

MumpsSchurSolver::MumpsSchurSolver(const mfem::HypreParMatrix &A,
                                   const std::vector<HYPRE_BigInt> &schur_vars,
                                   double blr_tol, bool serial, bool refactor)
  : comm(A.GetComm()), serial(serial), refactor(refactor), n_glob(A.GetGlobalNumRows()),
    n_loc(A.Height()), n_schur(static_cast<int>(schur_vars.size())), blr_tol(blr_tol)
{
  MPI_Comm_rank(comm, &rank);
  active = !serial || rank == 0;
  int nranks;
  MPI_Comm_size(comm, &nranks);
  row_cnt.assign(nranks, 0);
  row_disp.assign(nranks, 0);
  MPI_Allgather(&n_loc, 1, MPI_INT, row_cnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    row_disp[r] = row_disp[r - 1] + row_cnt[r - 1];
  }

  LowerTriangleCOO(A, irn, jcn, val);
  MFEM_VERIFY(!(refactor && serial), "MumpsSchurSolver: Refactor needs distributed input!");
  if (blr_tol > 0.0)
  {
    // The BLR dropping parameter CNTL(7) is absolute, and the Schur option excludes
    // MUMPS's own scaling: scale the matrix by its largest diagonal entry so the
    // tolerance acts as a relative one (the Schur and the solves are rescaled exactly
    // below).
    double dmax = 0.0;
    for (std::size_t q = 0; q < irn.size(); q++)
    {
      if (irn[q] == jcn[q])
      {
        dmax = std::max(dmax, std::abs(val[q]));
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &dmax, 1, MPI_DOUBLE, MPI_MAX, comm);
    scale = (dmax > 0.0) ? dmax : 1.0;
    for (auto &v : val)
    {
      v /= scale;
    }
  }
  if (serial)
  {
    // Gather the entries on rank 0, which factors alone (centralized assembled input).
    static_assert(sizeof(MUMPS_INT) == sizeof(int), "MUMPS_INT must be a 32-bit int!");
    int nnz_loc = static_cast<int>(irn.size());
    std::vector<int> cnt(nranks), disp(nranks, 0);
    MPI_Gather(&nnz_loc, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, comm);
    std::size_t total = 0;
    for (int r = 0; rank == 0 && r < nranks; r++)
    {
      disp[r] = static_cast<int>(total);
      total += cnt[r];
    }
    std::vector<MUMPS_INT> irn_g(total), jcn_g(total);
    std::vector<double> val_g(total);
    MPI_Gatherv(irn.data(), nnz_loc, MPI_INT, irn_g.data(), cnt.data(), disp.data(),
                MPI_INT, 0, comm);
    MPI_Gatherv(jcn.data(), nnz_loc, MPI_INT, jcn_g.data(), cnt.data(), disp.data(),
                MPI_INT, 0, comm);
    MPI_Gatherv(val.data(), nnz_loc, MPI_DOUBLE, val_g.data(), cnt.data(), disp.data(),
                MPI_DOUBLE, 0, comm);
    irn.swap(irn_g);
    jcn.swap(jcn_g);
    val.swap(val_g);
  }
  for (HYPRE_BigInt v : schur_vars)
  {
    listvar.push_back(static_cast<MUMPS_INT>(v + 1));
  }

  auto icntl = [&](int i) -> MUMPS_INT & { return id.icntl[i - 1]; };
  if (active)
  {
    id.sym = 2;  // symmetric (general)
    id.par = 1;  // the host takes part in the factorization
    id.comm_fortran = static_cast<MUMPS_INT>(MPI_Comm_c2f(serial ? MPI_COMM_SELF : comm));
    id.job = -1;
    dmumps_c(&id);
    icntl(1) = -1;  // silence errors / diagnostics / global info / printing
    icntl(2) = -1;
    icntl(3) = -1;
    icntl(4) = 0;
    icntl(5) = 0;                // assembled input
    icntl(18) = serial ? 0 : 3;  // centralized or distributed matrix entries
    // Complete Schur (if any Schur variables), 2D block cyclic on a 1 x 1 grid: centralized
    // on rank 0.
    icntl(19) = (n_schur > 0) ? 3 : 0;
    icntl(28) = 1;   // sequential analysis (the Schur option excludes parallel analysis)
    icntl(7) = 5;    // METIS ordering
    icntl(20) = 0;   // dense, centralized right-hand sides
    icntl(21) = 0;   // centralized solution
    icntl(26) = 0;   // solve the internal problem (Schur variables held at 0)
    icntl(14) = 50;  // workspace relaxation (%); raised on a workspace failure below
    if (blr_tol > 0.0)
    {
      icntl(35) = 2;     // BLR in both factorization and solve (compressed factors)
      icntl(36) = 1;     // UCFS variant: compress earlier, fewer operations
      id.cntl[0] = 0.0;  // CNTL(1): no numerical pivoting (the environment operator is SPD)
      id.cntl[6] = blr_tol;  // CNTL(7): dropping parameter (relative, after the scaling)
    }
    id.n = static_cast<MUMPS_INT>(n_glob);
    if (serial)
    {
      id.nnz = static_cast<MUMPS_INT8>(irn.size());
      id.irn = irn.data();
      id.jcn = jcn.data();
      id.a = val.data();
    }
    else
    {
      id.nnz_loc = static_cast<MUMPS_INT8>(irn.size());
      id.irn_loc = irn.data();
      id.jcn_loc = jcn.data();
      id.a_loc = val.data();
    }
    id.size_schur = n_schur;
    id.listvar_schur = listvar.data();
    id.nprow = 1;
    id.npcol = 1;
    id.mblock = 64;
    id.nblock = 64;
    id.job = 1;  // analysis
    dmumps_c(&id);
  }
  Check("analysis");
  if (rank == 0 && n_schur > 0)
  {
    id.schur_lld = std::max<MUMPS_INT>(1, id.schur_mloc);
    schur.assign(static_cast<std::size_t>(id.schur_lld) *
                     std::max<MUMPS_INT>(1, id.schur_nloc),
                 0.0);
    id.schur = schur.data();
  }
  Factor();
}

void MumpsSchurSolver::Factor()
{
  // Factorization (+ Schur). The workspace estimate from the analysis can be too small
  // (the dense Schur root sits on one rank): on a workspace failure (INFOG(1) = -8, -9,
  // -20) double the relaxation and retry, as the MUMPS user guide recommends.
  for (int attempt = 0; active; attempt++)
  {
    id.job = 2;
    dmumps_c(&id);
    const MUMPS_INT err = id.infog[0];
    if ((err == -8 || err == -9 || err == -20) && attempt < 5)
    {
      id.icntl[13] *= 2;  // ICNTL(14)
      continue;
    }
    break;
  }
  Check("factorization");
  // The factors do not need the input entries (no iterative refinement), unless they are
  // refactored with new values.
  if (!refactor)
  {
    irn = {};
    jcn = {};
    val = {};
  }
  if (blr_tol > 0.0)
  {
    for (auto &v : schur)
    {
      v *= scale;  // Schur(A / s) = S / s
    }
    // INFOG(29/35): theoretical / effective factor entries; RINFOG(3/14): operations.
    auto big = [](MUMPS_INT v) { return v >= 0 ? static_cast<double>(v) : -1.0e6 * v; };
    Mpi::Print(
        comm,
        " MUMPS BLR (tol = {:.1e}): factor entries {:.1f}%, operations {:.1f}% of full "
        "rank\n",
        blr_tol, 100.0 * big(id.infog[34]) / big(id.infog[28]),
        100.0 * id.rinfog[13] / id.rinfog[2]);
  }
}

void MumpsSchurSolver::Refactor(const mfem::HypreParMatrix &A)
{
  MFEM_VERIFY(refactor, "MumpsSchurSolver was not set up for refactorization!");
  std::vector<MUMPS_INT> irn_new, jcn_new;
  std::vector<double> val_new;
  LowerTriangleCOO(A, irn_new, jcn_new, val_new);
  int same = (irn_new == irn && jcn_new == jcn);
  MPI_Allreduce(MPI_IN_PLACE, &same, 1, MPI_INT, MPI_MIN, comm);
  MFEM_VERIFY(same,
              "MumpsSchurSolver::Refactor needs the sparsity pattern of the analysis!");
  val.swap(val_new);
  for (auto &v : val)
  {
    v /= scale;
  }
  id.a_loc = val.data();
  Factor();
}

MumpsSchurSolver::~MumpsSchurSolver()
{
  if (active)
  {
    id.job = -2;
    dmumps_c(&id);
  }
}

void MumpsSchurSolver::SolveInternal(const std::vector<const mfem::Vector *> &X,
                                     const std::vector<mfem::Vector *> &Y)
{
  const int n = static_cast<int>(X.size());
  constexpr int B = 32;
  for (int c0 = 0; c0 < n; c0 += B)
  {
    const int nb = std::min(B, n - c0);
    if (rank == 0)
    {
      rhs.assign(static_cast<std::size_t>(n_glob) * nb, 0.0);
    }
    for (int k = 0; k < nb; k++)
    {
      const mfem::Vector &x = *X[c0 + k];
      x.HostRead();
      MPI_Gatherv(x.GetData(), n_loc, MPI_DOUBLE,
                  rank == 0 ? rhs.data() + static_cast<std::size_t>(k) * n_glob : nullptr,
                  row_cnt.data(), row_disp.data(), MPI_DOUBLE, 0, comm);
    }
    if (rank == 0)
    {
      id.rhs = rhs.data();
      id.nrhs = nb;
      id.lrhs = static_cast<MUMPS_INT>(n_glob);
    }
    if (active)
    {
      id.job = 3;
      dmumps_c(&id);
    }
    Check("solve");
    if (rank == 0 && scale != 1.0)
    {
      for (auto &v : rhs)
      {
        v /= scale;  // (A / s)^-1 b = s A^-1 b
      }
    }
    for (int k = 0; k < nb; k++)
    {
      mfem::Vector &y = *Y[c0 + k];
      y.SetSize(n_loc);
      MPI_Scatterv(rank == 0 ? rhs.data() + static_cast<std::size_t>(k) * n_glob : nullptr,
                   row_cnt.data(), row_disp.data(), MPI_DOUBLE, y.HostWrite(), n_loc,
                   MPI_DOUBLE, 0, comm);
    }
  }
}

void MumpsSchurSolver::Check(const char *phase) const
{
  // INFOG is global (identical on the ranks taking part); serially, rank 0 has it.
  int info[2] = {id.infog[0], id.infog[1]};
  if (serial)
  {
    MPI_Bcast(info, 2, MPI_INT, 0, comm);
  }
  MFEM_VERIFY(info[0] >= 0, "MUMPS " << phase << " failed: INFOG(1) = " << info[0]
                                     << ", INFOG(2) = " << info[1]);
}

}  // namespace palace

#endif
