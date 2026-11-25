// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mumpsschur.hpp"

#if defined(MFEM_USE_MUMPS)

#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Free a vector's storage (v = {} empties it but keeps its capacity).
template <typename V>
void Release(V &v)
{
  V().swap(v);
}

// The number of the n indices of a block-cyclic dimension (blocks of b over np processes,
// from process 0) held by process me (ScaLAPACK's NUMROC).
int LocalCount(int n, int b, int me, int np)
{
  const int blocks = n / b, extra = blocks % np;
  return (blocks / np) * b + ((me < extra) ? b : (me == extra) ? n % b : 0);
}

// The global index of the local index l of process me in that dimension.
int GlobalIndex(int l, int b, int me, int np)
{
  return ((l / b) * np + me) * b + l % b;
}

}  // namespace

namespace
{

template <typename T>
typename MumpsSchurSolverT<T>::Coo LowerTriangle(const mfem::HypreParMatrix &A)
{
  typename MumpsSchurSolverT<T>::Coo coo;
  std::vector<std::vector<double>> vals;
  std::vector<int> row_ptr;
  LowerTrianglePattern({&A}, coo.irn, coo.jcn, vals, row_ptr);
  if constexpr (std::is_same_v<T, double>)
  {
    coo.val = std::move(vals[0]);
  }
  else
  {
    coo.val.assign(vals[0].begin(), vals[0].end());
  }
  return coo;
}

}  // namespace

template <typename T>
MumpsSchurSolverT<T>::MumpsSchurSolverT(const mfem::HypreParMatrix &A,
                                        const std::vector<HYPRE_BigInt> &schur_vars,
                                        double blr_tol, bool refactor, bool spd,
                                        SymbolicFactorization ordering, int block)
  : MumpsSchurSolverT(A.GetComm(), A.GetGlobalNumRows(), A.Height(), LowerTriangle<T>(A),
                      schur_vars, blr_tol, refactor, spd, {}, ordering, block)
{
}

template <typename T>
MumpsSchurSolverT<T>::MumpsSchurSolverT(MPI_Comm comm, HYPRE_BigInt n_glob, int n_loc,
                                        Coo &&A,
                                        const std::vector<HYPRE_BigInt> &schur_vars,
                                        double blr_tol, bool refactor, bool spd,
                                        std::vector<int> rows,
                                        SymbolicFactorization ordering, int block,
                                        int grid_rows, int grid_cols)
  : comm(comm), refactor(refactor), spd(spd), ordering(ordering), n_glob(n_glob),
    n_loc(n_loc), n_schur(static_cast<int>(schur_vars.size())), nprow(grid_rows),
    npcol(grid_cols), mb(block), gather(comm, n_loc, schur_vars), rows(std::move(rows)),
    irn(std::move(A.irn)), jcn(std::move(A.jcn)), val(std::move(A.val)), blr_tol(blr_tol)
{
  // Parent vectors if any rank gives rows (a rank with no local rows gives none).
  int mine = !this->rows.empty();
  MPI_Allreduce(&mine, &parent, 1, MPI_INT, MPI_MAX, comm);
  MFEM_VERIFY(!parent || static_cast<int>(this->rows.size()) == n_loc,
              "MumpsSchurSolver: wrong number of parent rows!");
  Init(schur_vars);
}

template <typename T>
void MumpsSchurSolverT<T>::Init(const std::vector<HYPRE_BigInt> &schur_vars)
{
  MPI_Comm_rank(comm, &rank);
  int nranks;
  MPI_Comm_size(comm, &nranks);
  if (nprow <= 0 || npcol <= 0)
  {
    nprow = static_cast<int>(std::sqrt(static_cast<double>(nranks)));
    npcol = nranks / nprow;
  }
  MFEM_VERIFY(nprow * npcol <= nranks, "MumpsSchurSolver: grid larger than the ranks!");
  row_cnt.assign(nranks, 0);
  row_disp.assign(nranks, 0);
  MPI_Allgather(&n_loc, 1, MPI_INT, row_cnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    row_disp[r] = row_disp[r - 1] + row_cnt[r - 1];
  }

  if (blr_tol > 0.0)
  {
    // The BLR dropping parameter CNTL(7) is absolute, and the Schur option excludes
    // MUMPS's own scaling: scale the matrix by its largest diagonal entry so the
    // tolerance acts as a relative one (the Schur and the solves are rescaled exactly
    // below), at each factorization.
    for (std::size_t q = 0; q < irn.size(); q++)
    {
      if (irn[q] == jcn[q])
      {
        diag_pos.push_back(q);
      }
    }
    scale = DiagonalScale();
    for (auto &v : val)
    {
      v /= scale;
    }
  }
  // The rows of the right-hand sides this rank passes to MUMPS (1-based): its own.
  irhs_loc.resize(std::max(n_loc, 1), 1);
  for (int i = 0; i < n_loc; i++)
  {
    irhs_loc[i] = static_cast<MUMPS_INT>(row_disp[rank] + i + 1);
  }
  for (HYPRE_BigInt v : schur_vars)
  {
    listvar.push_back(static_cast<MUMPS_INT>(v + 1));
  }

  auto icntl = [&](int i) -> MUMPS_INT & { return id.icntl[i - 1]; };
  id.sym = (spd && !kComplex) ? 1 : 2;  // positive definite (no pivoting), or symmetric
  id.par = 1;                           // the host takes part in the factorization
  id.comm_fortran = static_cast<MUMPS_INT>(MPI_Comm_c2f(comm));
  id.job = -1;
  Call();
  icntl(1) = -1;  // silence errors / diagnostics / global info / printing
  icntl(2) = -1;
  icntl(3) = -1;
  icntl(4) = 0;
  icntl(5) = 0;   // assembled input
  icntl(18) = 3;  // distributed matrix entries
  // Complete Schur (if any Schur variables), 2D block cyclic on a near-square grid of the
  // ranks (MUMPS symmetrizes it with blocks of the same rows and columns).
  icntl(19) = (n_schur > 0) ? 3 : 0;
  // Fill-reducing ordering: sequential (ICNTL(7)), or parallel (ICNTL(29)) for ParMETIS
  // and PT-SCOTCH if there are no Schur variables.
  const bool parallel = n_schur == 0 && (ordering == SymbolicFactorization::PARMETIS ||
                                         ordering == SymbolicFactorization::PTSCOTCH);
  icntl(28) = parallel ? 2 : 1;
  switch (ordering)
  {
    case SymbolicFactorization::PARMETIS:
      icntl(7) = 5;  // METIS
      icntl(29) = 2;
      break;
    case SymbolicFactorization::PTSCOTCH:
      icntl(7) = 3;  // SCOTCH
      icntl(29) = 1;
      break;
    case SymbolicFactorization::SCOTCH:
      icntl(7) = 3;
      break;
    case SymbolicFactorization::PORD:
      icntl(7) = 4;
      break;
    case SymbolicFactorization::AMD:
    case SymbolicFactorization::RCM:
      icntl(7) = 0;
      break;
    case SymbolicFactorization::METIS:
    case SymbolicFactorization::DEFAULT:
      icntl(7) = 5;
      break;
  }
  icntl(14) = 50;  // workspace relaxation (%); raised on a workspace failure below
  if (blr_tol > 0.0)
  {
    icntl(35) = 2;  // BLR in both factorization and solve (compressed factors)
    icntl(36) = 1;  // UCFS variant: compress earlier, fewer operations
    if (spd)
    {
      id.cntl[0] = 0.0;  // CNTL(1): no numerical pivoting
    }
    id.cntl[6] = blr_tol;  // CNTL(7): dropping parameter (relative, after the scaling)
  }
  id.n = static_cast<MUMPS_INT>(n_glob);
  id.nnz_loc = static_cast<MUMPS_INT8>(irn.size());
  id.irn_loc = irn.data();
  id.jcn_loc = jcn.data();
  id.a_loc = Entries(val.data());
  id.size_schur = n_schur;
  id.listvar_schur = listvar.data();
  id.nprow = nprow;
  id.npcol = npcol;
  id.mblock = mb;
  id.nblock = mb;
  id.job = 1;  // analysis
  Call();
  Check("analysis");
  if (InGrid() && n_schur > 0)
  {
    id.schur_lld = std::max<MUMPS_INT>(1, id.schur_mloc);
  }
  Factor();
}

template <typename T>
std::vector<T> MumpsSchurSolverT<T>::ReleaseSchurBlocks(const BlockCyclicGrid &grid)
{
  MFEM_VERIFY(grid.Size() == n_schur && grid.NumProcessRows() == nprow &&
                  grid.NumProcessCols() == npcol && grid.BlockSize() == mb,
              "MumpsSchurSolver: the grid is not that of its Schur complement!");
  std::vector<T> S;
  S.swap(schur);
  return S;
}

template <typename T>
std::vector<T> MumpsSchurSolverT<T>::ReleaseSchurRows()
{
  // S is 2D block cyclic (local parts column-major) on the grid of the first nprow npcol
  // ranks, rank p at (p / npcol, p % npcol). Each round sends
  // one local block column of every process to the owners of its rows, column by column and
  // in a column the rows of each owner in turn, which place them by the same traversal.
  const int nranks = static_cast<int>(row_cnt.size()), ns = n_schur,
            m = static_cast<int>(gather.Rows().size());
  std::vector<T> S(static_cast<std::size_t>(m) * ns);
  if (ns == 0)
  {
    return S;
  }
  // The owner of each Schur row and its place among the rows of its owner, and the places
  // of the rows of each process row owned here.
  std::vector<int> owner(ns), slot(ns), count(nranks, 0);
  for (int p = 0; p < ns; p++)
  {
    owner[p] = static_cast<int>(
                   std::upper_bound(row_disp.begin(), row_disp.end(), listvar[p] - 1) -
                   row_disp.begin()) -
               1;
    slot[p] = count[owner[p]]++;
  }
  std::vector<std::vector<int>> mine(nprow);
  for (int pr = 0; pr < nprow; pr++)
  {
    for (int l = 0; l < LocalCount(ns, mb, pr, nprow); l++)
    {
      const int g = GlobalIndex(l, mb, pr, nprow);
      if (owner[g] == rank)
      {
        mine[pr].push_back(slot[g]);
      }
    }
  }
  const int ngrid = nprow * npcol;
  const bool grid = InGrid();
  const int myrow = grid ? rank / npcol : 0, mycol = grid ? rank % npcol : 0,
            mloc = grid ? LocalCount(ns, mb, myrow, nprow) : 0;
  const std::size_t lld = std::max(1, mloc);
  std::vector<int> dest(mloc);
  for (int l = 0; l < mloc; l++)
  {
    dest[l] = owner[GlobalIndex(l, mb, myrow, nprow)];
  }
  std::vector<int> scnt(nranks), sdsp(nranks, 0), rcnt(nranks, 0), rdsp(nranks, 0), at;
  std::vector<T> sbuf, rbuf;
  const int rounds = (LocalCount(ns, mb, 0, npcol) + mb - 1) / mb;
  for (int k = 0; k < rounds; k++)
  {
    // The local columns [k mb, k mb + width) of the processes of a grid column.
    auto width = [&](int pc)
    { return std::clamp(LocalCount(ns, mb, pc, npcol) - k * mb, 0, mb); };
    const int w = grid ? width(mycol) : 0;
    std::fill(scnt.begin(), scnt.end(), 0);
    for (int l = 0; l < mloc; l++)
    {
      scnt[dest[l]] += w;
    }
    for (int p = 0; p < ngrid; p++)
    {
      rcnt[p] = static_cast<int>(mine[p / npcol].size()) * width(p % npcol);
    }
    for (int r = 1; r < nranks; r++)
    {
      sdsp[r] = sdsp[r - 1] + scnt[r - 1];
      rdsp[r] = rdsp[r - 1] + rcnt[r - 1];
    }
    sbuf.resize(sdsp[nranks - 1] + scnt[nranks - 1]);
    rbuf.resize(rdsp[nranks - 1] + rcnt[nranks - 1]);
    at = sdsp;
    for (int j = 0; j < w; j++)
    {
      const T *col = schur.data() + lld * (k * mb + j);
      for (int l = 0; l < mloc; l++)
      {
        sbuf[at[dest[l]]++] = col[l];
      }
    }
    MPI_Alltoallv(sbuf.data(), scnt.data(), sdsp.data(), mpi::DataType<T>(), rbuf.data(),
                  rcnt.data(), rdsp.data(), mpi::DataType<T>(), comm);
    for (int p = 0; p < ngrid; p++)
    {
      const int pr = p / npcol, pc = p % npcol;
      std::size_t q = rdsp[p];
      for (int j = 0; j < width(pc); j++)
      {
        const int c = GlobalIndex(k * mb + j, mb, pc, npcol);
        for (int r : mine[pr])
        {
          S[static_cast<std::size_t>(r) * ns + c] = rbuf[q++];
        }
      }
    }
  }
  std::vector<T>().swap(schur);
  return S;
}

template <typename T>
void MumpsSchurSolverT<T>::Call()
{
  if constexpr (kComplex)
  {
    zmumps_c(&id);
  }
  else
  {
    dmumps_c(&id);
  }
}

template <typename T>
void MumpsSchurSolverT<T>::Factor()
{
  // Factorization (+ Schur). The workspace estimates of the analysis can be too small in
  // parallel (the dense Schur root sits on one rank). If the integer workspace or the
  // buffers are (INFOG(1) = -8, -17, -20), double their relaxation ICNTL(14) and retry. If
  // the real workspace is (-9), its shortfall can be in the main workspace, which MUMPS
  // sizes from the estimate relaxed by ICNTL(14), or in the allocations it makes within the
  // memory allowed per process (ICNTL(23), in MB, which MUMPS reads on each process when it
  // is set on more than the host), such as that of the dense Schur root, which a
  // refactorization with ICNTL(23) set can exceed (-19): double the relaxation and the
  // memory allowed, from what each process allocated (INFO(16)) with twice their shortfall
  // (INFO(2) entries) more on the processes that ran short, and retry.
  if (InGrid() && n_schur > 0 && schur.empty())
  {
    // (Again after ReleaseSchurRows.)
    schur.assign(static_cast<std::size_t>(id.schur_lld) *
                     std::max<MUMPS_INT>(1, id.schur_nloc),
                 T(0.0));
    id.schur = Entries(schur.data());
  }
  for (int attempt = 0;; attempt++)
  {
    id.job = 2;
    Call();
    const MUMPS_INT err = id.infog[0];
    if (attempt >= 5 || (err != -8 && err != -9 && err != -17 && err != -19 && err != -20))
    {
      break;
    }
    Mpi::Print(comm,
               " MUMPS factorization short of workspace (INFOG(1) = {:d}, INFOG(2) = "
               "{:d}): retrying\n",
               err, id.infog[1]);
    id.icntl[13] *= 2;  // ICNTL(14)
    if (err != -9 && err != -19)
    {
      continue;
    }
    double short_mb = 0.0;
    if (id.info[0] == -9 || id.info[0] == -19)
    {
      const double entries = (id.info[1] > 0) ? id.info[1] : -1.0e6 * id.info[1];
      short_mb = entries * sizeof(T) / 1.0e6;
    }
    id.icntl[22] = 2 * std::max<MUMPS_INT>(id.icntl[22], id.info[15]) +  // ICNTL(23)
                   static_cast<MUMPS_INT>(std::ceil(2.0 * short_mb)) + 1;
  }
  Check("factorization");
  // The factors do not need the input entries (no iterative refinement), unless they are
  // refactored with new values.
  if (!refactor)
  {
    Release(irn);
    Release(jcn);
    Release(val);
  }
  if (blr_tol > 0.0)
  {
    for (auto &v : schur)
    {
      v *= scale;  // Schur(A / s) = S / s
    }
    if (!factored)
    {
      // INFOG(21/22), (18/19), (16/17): memory used, allocated and estimated (full rank),
      // max per rank / total (MB); INFOG(29/35): theoretical / effective factor entries.
      auto big = [](MUMPS_INT v) { return v >= 0 ? static_cast<double>(v) : -1.0e6 * v; };
      const int nranks = static_cast<int>(row_cnt.size());
      Mpi::Print(
          comm,
          " MUMPS BLR (tol = {:.1e}): factor entries {:.1f}% of full rank, memory "
          "{:.2f} GB used of {:.2f} GB allocated ({:.2f} GB estimated) on {:d} rank{} "
          "({:.2f} / {:.2f} GB max per rank)\n",
          blr_tol, 100.0 * big(id.infog[34]) / big(id.infog[28]), id.infog[21] / 1024.0,
          id.infog[18] / 1024.0, id.infog[16] / 1024.0, nranks, (nranks > 1) ? "s" : "",
          id.infog[20] / 1024.0, id.infog[17] / 1024.0);
    }
  }
  factored = true;
}

template <typename T>
double MumpsSchurSolverT<T>::DiagonalScale() const
{
  double dmax = 0.0;
  for (std::size_t q : diag_pos)
  {
    dmax = std::max(dmax, static_cast<double>(std::abs(val[q])));
  }
  MPI_Allreduce(MPI_IN_PLACE, &dmax, 1, MPI_DOUBLE, MPI_MAX, comm);
  return (dmax > 0.0) ? dmax : 1.0;
}

template <typename T>
void MumpsSchurSolverT<T>::Refactor()
{
  MFEM_VERIFY(refactor, "MumpsSchurSolver was not set up for refactorization!");
  if (blr_tol > 0.0)
  {
    scale = DiagonalScale();
  }
  for (auto &v : val)
  {
    v /= scale;
  }
  id.a_loc = Entries(val.data());
  Factor();
}

template <typename T>
MumpsSchurSolverT<T>::~MumpsSchurSolverT()
{
  id.job = -2;
  Call();
}

template <typename T>
void MumpsSchurSolverT<T>::SetRhs(const std::vector<const VecType *> &X)
{
  // The right-hand sides distributed by rows (ICNTL(20) = 10): each rank passes its rows,
  // column-major.
  const int nb = static_cast<int>(X.size());
  std::vector<T> loc(static_cast<std::size_t>(n_loc) * nb);
  for (int k = 0; k < nb; k++)
  {
    T *l = loc.data() + static_cast<std::size_t>(k) * n_loc;
    const VecType &x = *X[k];
    if constexpr (kComplex)
    {
      const double *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
      for (int i = 0; i < n_loc; i++)
      {
        const int r = parent ? rows[i] : i;
        l[i] = (r >= 0) ? T(xr[r], xi[r]) : T(0.0);
      }
    }
    else
    {
      const double *xr = x.HostRead();
      for (int i = 0; i < n_loc; i++)
      {
        const int r = parent ? rows[i] : i;
        l[i] = (r >= 0) ? xr[r] : 0.0;
      }
    }
  }
  loc.resize(static_cast<std::size_t>(std::max(n_loc, 1)) * nb);
  rhs_loc.swap(loc);
  id.icntl[19] = 10;  // ICNTL(20): distributed right-hand sides
  id.nrhs = nb;
  id.nloc_rhs = n_loc;
  id.lrhs_loc = std::max(n_loc, 1);
  id.rhs_loc = Entries(rhs_loc.data());
  id.irhs_loc = irhs_loc.data();
}

template <typename T>
void MumpsSchurSolverT<T>::SetSolution(int nb)
{
  // Room for the solution distributed as MUMPS leaves it (ICNTL(21) = 1): INFO(23) rows on
  // each rank after the factorization.
  const int ns = std::max<MUMPS_INT>(id.info[22], 1);
  sol_loc.assign(static_cast<std::size_t>(ns) * nb, T(0.0));
  isol_loc.assign(ns, 0);
  id.icntl[20] = 1;  // ICNTL(21): distributed solution
  id.nrhs = nb;
  id.lsol_loc = ns;
  id.sol_loc = Entries(sol_loc.data());
  id.isol_loc = isol_loc.data();
}

template <typename T>
void MumpsSchurSolverT<T>::ScatterSolution(const std::vector<VecType *> &Y)
{
  // The entries of MUMPS's distributed solution, with the Schur variables (0 for the
  // internal problem, u_S for an expansion), to the ranks owning their rows, unscaled:
  // (A / s)^-1 b = s A^-1 b.
  const int nb = static_cast<int>(Y.size()), nranks = static_cast<int>(row_cnt.size());
  const int ns = id.info[22];
  std::vector<int> owner(ns), scnt(nranks, 0), rcnt(nranks), sdsp(nranks, 0),
      rdsp(nranks, 0);
  for (int q = 0; q < ns; q++)
  {
    const int g = isol_loc[q] - 1;
    owner[q] = static_cast<int>(std::upper_bound(row_disp.begin(), row_disp.end(), g) -
                                row_disp.begin()) -
               1;
    scnt[owner[q]]++;
  }
  MPI_Alltoall(scnt.data(), 1, MPI_INT, rcnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    sdsp[r] = sdsp[r - 1] + scnt[r - 1];
    rdsp[r] = rdsp[r - 1] + rcnt[r - 1];
  }
  const int nrecv = rdsp.back() + rcnt.back();
  std::vector<int> sidx(ns), ridx(nrecv), pos(sdsp);
  std::vector<T> sval(static_cast<std::size_t>(ns) * nb);
  for (int q = 0; q < ns; q++)
  {
    const int r = owner[q], e = pos[r]++;
    sidx[e] = isol_loc[q] - 1 - row_disp[r];
    for (int k = 0; k < nb; k++)
    {
      sval[static_cast<std::size_t>(e) * nb + k] =
          sol_loc[static_cast<std::size_t>(k) * id.lsol_loc + q] / scale;
    }
  }
  Release(sol_loc);
  Release(isol_loc);
  MPI_Alltoallv(sidx.data(), scnt.data(), sdsp.data(), MPI_INT, ridx.data(), rcnt.data(),
                rdsp.data(), MPI_INT, comm);
  for (int r = 0; r < nranks; r++)
  {
    scnt[r] *= nb;
    sdsp[r] *= nb;
    rcnt[r] *= nb;
    rdsp[r] *= nb;
  }
  std::vector<T> rval(static_cast<std::size_t>(nrecv) * nb);
  MPI_Alltoallv(sval.data(), scnt.data(), sdsp.data(), mpi::DataType<T>(), rval.data(),
                rcnt.data(), rdsp.data(), mpi::DataType<T>(), comm);
  Release(sval);
  for (int k = 0; k < nb; k++)
  {
    VecType &y = *Y[k];
    if (!parent)
    {
      y.SetSize(n_loc);
    }
    y = 0.0;
    if constexpr (kComplex)
    {
      double *yr = y.Real().HostReadWrite(), *yi = y.Imag().HostReadWrite();
      for (int e = 0; e < nrecv; e++)
      {
        const int r = parent ? rows[ridx[e]] : ridx[e];
        if (r >= 0)
        {
          const T v = rval[static_cast<std::size_t>(e) * nb + k];
          yr[r] = v.real();
          yi[r] = v.imag();
        }
      }
    }
    else
    {
      double *yr = y.HostReadWrite();
      for (int e = 0; e < nrecv; e++)
      {
        const int r = parent ? rows[ridx[e]] : ridx[e];
        if (r >= 0)
        {
          yr[r] = rval[static_cast<std::size_t>(e) * nb + k];
        }
      }
    }
  }
}

template <typename T>
void MumpsSchurSolverT<T>::SolveInternal(const std::vector<const VecType *> &X,
                                         const std::vector<VecType *> &Y)
{
  MFEM_VERIFY(!reduced, "MumpsSchurSolver: a Reduce is pending its Expand!");
  const int n = static_cast<int>(X.size());
  constexpr int B = 32;
  for (int c0 = 0; c0 < n; c0 += B)
  {
    const int nb = std::min(B, n - c0);
    SetRhs({X.begin() + c0, X.begin() + c0 + nb});
    SetSolution(nb);
    id.icntl[25] = 0;  // ICNTL(26): the internal problem (Schur variables held at 0)
    id.job = 3;
    Call();
    Check("solve");
    Release(rhs_loc);
    ScatterSolution({Y.begin() + c0, Y.begin() + c0 + nb});
  }
}

template <typename T>
void MumpsSchurSolverT<T>::Reduce(const std::vector<const VecType *> &B,
                                  std::vector<T> &red)
{
  // The forward elimination is kept for Expand (MUMPS's RHSINTR), so the whole batch is
  // condensed in one solve. The condensation is independent of the scaling of A.
  MFEM_VERIFY(n_schur > 0 && !reduced, "MumpsSchurSolver: invalid Reduce!");
  const int nb = static_cast<int>(B.size());
  SetRhs(B);
  SetSolution(nb);
  if (rank == 0)
  {
    redrhs.assign(static_cast<std::size_t>(n_schur) * std::max(nb, 1), T(0.0));
    id.redrhs = Entries(redrhs.data());
    id.lredrhs = n_schur;
  }
  id.icntl[25] = 1;  // ICNTL(26): condensation onto the Schur variables
  id.job = 3;
  Call();
  Check("condensation");
  Release(rhs_loc);
  Release(sol_loc);
  Release(isol_loc);
  if (rank == 0)
  {
    red.swap(redrhs);
    red.resize(static_cast<std::size_t>(n_schur) * nb);
  }
  reduced = nb;
}

template <typename T>
void MumpsSchurSolverT<T>::Expand(const std::vector<T> &u, const std::vector<VecType *> &X)
{
  MFEM_VERIFY(reduced && static_cast<int>(X.size()) == reduced,
              "MumpsSchurSolver: Expand must follow Reduce of the same batch!");
  const int nb = reduced;
  if (rank == 0)
  {
    // The Schur part of the solution of the scaled system (A / s) x = b is s u_S.
    redrhs.resize(static_cast<std::size_t>(n_schur) * nb);
    for (std::size_t q = 0; q < redrhs.size(); q++)
    {
      redrhs[q] = u[q] * scale;
    }
    id.redrhs = Entries(redrhs.data());
    id.lredrhs = n_schur;
  }
  SetSolution(nb);
  // The right-hand sides are those of the condensation (kept by MUMPS).
  id.icntl[19] = 0;  // ICNTL(20)
  id.rhs = nullptr;
  id.icntl[25] = 2;  // ICNTL(26): expansion from the Schur part of the solution
  id.job = 3;
  Call();
  Check("expansion");
  Release(redrhs);
  reduced = 0;
  ScatterSolution(X);
}

template <typename T>
void MumpsSchurSolverT<T>::Condense(const std::vector<const VecType *> &X,
                                    std::vector<T> &g)
{
  // A Reduce that no Expand follows (the next solve replaces MUMPS's RHSINTR).
  Reduce(X, g);
  reduced = 0;
}

template <typename T>
void MumpsSchurSolverT<T>::Check(const char *phase) const
{
  // INFOG is global (identical on all ranks).
  const int info[2] = {id.infog[0], id.infog[1]};
  MFEM_VERIFY(info[0] >= 0, "MUMPS " << phase << " failed: INFOG(1) = " << info[0]
                                     << ", INFOG(2) = " << info[1] << "!");
}

template class MumpsSchurSolverT<double>;
template class MumpsSchurSolverT<std::complex<double>>;

}  // namespace palace

#endif
