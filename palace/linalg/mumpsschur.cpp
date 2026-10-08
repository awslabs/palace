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

// f(i, j, a) for the entries of local row i of A with global column j <= global row.
template <typename F>
void ForEachLowerEntry(const mfem::HypreParMatrix &A, int i, F &&f)
{
  auto *parcsr = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(A);
  const hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(parcsr);
  const hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(parcsr);
  const HYPRE_BigInt *cmap = hypre_ParCSRMatrixColMapOffd(parcsr);
  const HYPRE_BigInt row = hypre_ParCSRMatrixFirstRowIndex(parcsr) + i;
  const HYPRE_BigInt col0 = hypre_ParCSRMatrixFirstColDiag(parcsr);
  for (HYPRE_Int k = diag->i[i]; k < diag->i[i + 1]; k++)
  {
    const HYPRE_BigInt j = col0 + diag->j[k];
    if (j <= row)
    {
      f(row, j, diag->data[k]);
    }
  }
  for (HYPRE_Int k = offd->i ? offd->i[i] : 0; offd->i && k < offd->i[i + 1]; k++)
  {
    const HYPRE_BigInt j = cmap[offd->j[k]];
    if (j <= row)
    {
      f(row, j, offd->data[k]);
    }
  }
}

template <typename T>
MPI_Datatype MpiType()
{
  return std::is_same_v<T, double> ? MPI_DOUBLE : MPI_C_DOUBLE_COMPLEX;
}

}  // namespace

void LowerTrianglePattern(const std::vector<const mfem::HypreParMatrix *> &A,
                          std::vector<MUMPS_INT> &irn, std::vector<MUMPS_INT> &jcn,
                          std::vector<std::vector<double>> &vals, std::vector<int> &row_ptr)
{
  const mfem::HypreParMatrix *A0 = nullptr;
  for (const auto *X : A)
  {
    if (X)
    {
      MFEM_VERIFY(!A0 || (X->Height() == A0->Height() &&
                          X->GetRowStarts()[0] == A0->GetRowStarts()[0]),
                  "Matrices of a pattern must have the same row distribution!");
      A0 = A0 ? A0 : X;
      X->HostRead();
    }
  }
  MFEM_VERIFY(A0, "LowerTrianglePattern needs a matrix!");
  const int n = A0->Height();
  const HYPRE_BigInt row0 = A0->GetRowStarts()[0];
  // Two passes over the rows (count, then fill) to allocate the pattern exactly.
  std::vector<std::pair<HYPRE_BigInt, int>> row;  // (column, matrix)
  std::vector<double> rv;
  auto gather = [&](int i)
  {
    row.clear();
    rv.clear();
    for (std::size_t p = 0; p < A.size(); p++)
    {
      if (A[p])
      {
        ForEachLowerEntry(*A[p], i,
                          [&](HYPRE_BigInt, HYPRE_BigInt j, double a)
                          {
                            row.emplace_back(j, static_cast<int>(p));
                            rv.push_back(a);
                          });
      }
    }
  };
  std::vector<int> perm;
  auto sorted = [&]()
  {
    perm.resize(row.size());
    for (std::size_t q = 0; q < perm.size(); q++)
    {
      perm[q] = static_cast<int>(q);
    }
    std::sort(perm.begin(), perm.end(),
              [&](int a, int b) { return row[a].first < row[b].first; });
  };
  row_ptr.assign(n + 1, 0);
  for (int i = 0; i < n; i++)
  {
    gather(i);
    sorted();
    int m = 0;
    for (std::size_t q = 0; q < perm.size(); q++)
    {
      m += (q == 0 || row[perm[q]].first != row[perm[q - 1]].first);
    }
    row_ptr[i + 1] = row_ptr[i] + m;
  }
  const std::size_t nnz = row_ptr[n];
  irn.assign(nnz, 0);
  jcn.assign(nnz, 0);
  vals.assign(A.size(), {});
  for (std::size_t p = 0; p < A.size(); p++)
  {
    if (A[p])
    {
      vals[p].assign(nnz, 0.0);
    }
  }
  for (int i = 0; i < n; i++)
  {
    gather(i);
    sorted();
    int pos = row_ptr[i] - 1;
    for (std::size_t q = 0; q < perm.size(); q++)
    {
      const auto [j, p] = row[perm[q]];
      if (q == 0 || j != row[perm[q - 1]].first)
      {
        pos++;
        irn[pos] = static_cast<MUMPS_INT>(row0 + i + 1);
        jcn[pos] = static_cast<MUMPS_INT>(j + 1);
      }
      vals[p][pos] += rv[perm[q]];
    }
  }
}

namespace
{

template <typename T>
void AddToPatternT(const mfem::HypreParMatrix &X, T a, const std::vector<int> &row_ptr,
                   const std::vector<MUMPS_INT> &jcn, std::vector<T> &val)
{
  X.HostRead();
  for (int i = 0; i < X.Height(); i++)
  {
    const auto first = jcn.begin() + row_ptr[i], last = jcn.begin() + row_ptr[i + 1];
    ForEachLowerEntry(X, i,
                      [&](HYPRE_BigInt, HYPRE_BigInt j, double x)
                      {
                        const auto it =
                            std::lower_bound(first, last, static_cast<MUMPS_INT>(j + 1));
                        if (it != last && *it == j + 1)
                        {
                          val[it - jcn.begin()] += a * x;
                        }
                        else
                        {
                          MFEM_VERIFY(x == 0.0, "Entry outside of the sparsity pattern!");
                        }
                      });
  }
}

}  // namespace

void AddToPattern(const mfem::HypreParMatrix &X, double a, const std::vector<int> &row_ptr,
                  const std::vector<MUMPS_INT> &jcn, std::vector<double> &val)
{
  AddToPatternT(X, a, row_ptr, jcn, val);
}

void AddToPattern(const mfem::HypreParMatrix &X, std::complex<double> a,
                  const std::vector<int> &row_ptr, const std::vector<MUMPS_INT> &jcn,
                  std::vector<std::complex<double>> &val)
{
  AddToPatternT(X, a, row_ptr, jcn, val);
}

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
                                        double blr_tol, int procs, bool refactor, bool spd)
  : MumpsSchurSolverT(A.GetComm(), A.GetGlobalNumRows(), A.Height(), LowerTriangle<T>(A),
                      schur_vars, blr_tol, procs, refactor, spd)
{
}

template <typename T>
MumpsSchurSolverT<T>::MumpsSchurSolverT(MPI_Comm comm, HYPRE_BigInt n_glob, int n_loc,
                                        Coo &&A,
                                        const std::vector<HYPRE_BigInt> &schur_vars,
                                        double blr_tol, int procs, bool refactor, bool spd,
                                        std::vector<int> rows)
  : comm(comm), procs(procs), refactor(refactor), spd(spd), n_glob(n_glob), n_loc(n_loc),
    n_schur(static_cast<int>(schur_vars.size())), rows(std::move(rows)),
    irn(std::move(A.irn)), jcn(std::move(A.jcn)), val(std::move(A.val)), blr_tol(blr_tol)
{
  MFEM_VERIFY(this->rows.empty() || static_cast<int>(this->rows.size()) == n_loc,
              "MumpsSchurSolver: wrong number of parent rows!");
  Init(schur_vars);
}

template <typename T>
void MumpsSchurSolverT<T>::Init(const std::vector<HYPRE_BigInt> &schur_vars)
{
  MPI_Comm_rank(comm, &rank);
  int nranks;
  MPI_Comm_size(comm, &nranks);
  const int stride = (procs > 0 && procs < nranks) ? (nranks + procs - 1) / procs : 1;
  procs = (nranks + stride - 1) / stride;
  active = (rank % stride == 0);
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
    // below).
    double dmax = 0.0;
    for (std::size_t q = 0; q < irn.size(); q++)
    {
      if (irn[q] == jcn[q])
      {
        dmax = std::max(dmax, static_cast<double>(std::abs(val[q])));
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &dmax, 1, MPI_DOUBLE, MPI_MAX, comm);
    scale = (dmax > 0.0) ? dmax : 1.0;
    for (auto &v : val)
    {
      v /= scale;
    }
  }
  if (stride > 1)
  {
    // Groups of stride consecutive ranks, whose first rank factors: the pattern of the
    // group's entries on it.
    MPI_Comm_split(comm, active ? 0 : MPI_UNDEFINED, rank, &sub);
    MPI_Comm_split(comm, rank / stride, rank, &grp);
    int ng, nnz_loc = static_cast<int>(irn.size());
    MPI_Comm_size(grp, &ng);
    grp_cnt.assign(active ? ng : 0, 0);
    grp_disp.assign(active ? ng : 0, 0);
    MPI_Gather(&nnz_loc, 1, MPI_INT, grp_cnt.data(), 1, MPI_INT, 0, grp);
    for (int r = 1; r < static_cast<int>(grp_cnt.size()); r++)
    {
      grp_disp[r] = grp_disp[r - 1] + grp_cnt[r - 1];
    }
    const std::size_t total = active ? grp_disp.back() + grp_cnt.back() : 0;
    std::vector<MUMPS_INT> irn_g(total), jcn_g(total);
    MPI_Gatherv(irn.data(), nnz_loc, MPI_INT, irn_g.data(), grp_cnt.data(), grp_disp.data(),
                MPI_INT, 0, grp);
    MPI_Gatherv(jcn.data(), nnz_loc, MPI_INT, jcn_g.data(), grp_cnt.data(), grp_disp.data(),
                MPI_INT, 0, grp);
    irn.swap(irn_g);
    jcn.swap(jcn_g);
  }
  GatherEntries();
  for (HYPRE_BigInt v : schur_vars)
  {
    listvar.push_back(static_cast<MUMPS_INT>(v + 1));
  }

  auto icntl = [&](int i) -> MUMPS_INT & { return id.icntl[i - 1]; };
  if (active)
  {
    id.sym = 2;  // symmetric (general)
    id.par = 1;  // the host takes part in the factorization
    id.comm_fortran =
        static_cast<MUMPS_INT>(MPI_Comm_c2f(sub != MPI_COMM_NULL ? sub : comm));
    id.job = -1;
    Call();
    icntl(1) = -1;  // silence errors / diagnostics / global info / printing
    icntl(2) = -1;
    icntl(3) = -1;
    icntl(4) = 0;
    icntl(5) = 0;   // assembled input
    icntl(18) = 3;  // distributed matrix entries
    // Complete Schur (if any Schur variables), 2D block cyclic on a 1 x 1 grid: centralized
    // on rank 0.
    icntl(19) = (n_schur > 0) ? 3 : 0;
    icntl(28) = 1;   // sequential analysis (the Schur option excludes parallel analysis)
    icntl(7) = 5;    // METIS ordering
    icntl(20) = 0;   // dense, centralized right-hand sides
    icntl(21) = 0;   // centralized solution
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
    id.a_loc = Entries(grp != MPI_COMM_NULL ? grp_val.data() : val.data());
    id.size_schur = n_schur;
    id.listvar_schur = listvar.data();
    id.nprow = 1;
    id.npcol = 1;
    id.mblock = 64;
    id.nblock = 64;
    id.job = 1;  // analysis
    Call();
  }
  Check("analysis");
  if (rank == 0 && n_schur > 0)
  {
    id.schur_lld = std::max<MUMPS_INT>(1, id.schur_mloc);
    schur.assign(static_cast<std::size_t>(id.schur_lld) *
                     std::max<MUMPS_INT>(1, id.schur_nloc),
                 T(0.0));
    id.schur = Entries(schur.data());
  }
  Factor();
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
  // Factorization (+ Schur). The workspace estimate from the analysis can be too small
  // (the dense Schur root sits on one rank): on a workspace failure (INFOG(1) = -8, -9,
  // -20) raise the relaxation and retry, as the MUMPS user guide recommends. The missing
  // workspace (INFOG(2)) shrinks linearly with the relaxation ICNTL(14), whose base can be
  // small next to the Schur root, so the relaxation is extrapolated from two attempts.
  int r_prev = -1, m_prev = 0;
  for (int attempt = 0; active; attempt++)
  {
    id.job = 2;
    Call();
    const MUMPS_INT err = id.infog[0];
    if ((err == -8 || err == -9 || err == -20) && attempt < 8)
    {
      const int r = id.icntl[13], m = id.infog[1];  // ICNTL(14), INFOG(2)
      int r_next = 2 * r;
      if (r_prev >= 0 && m > 0 && m_prev > m)
      {
        const double slope = static_cast<double>(m_prev - m) / (r - r_prev);
        r_next = std::max(r_next, r + static_cast<int>(std::ceil(1.25 * m / slope)) + 10);
      }
      r_prev = r;
      m_prev = m;
      id.icntl[13] = r_next;
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
    grp_val = {};
  }
  if (blr_tol > 0.0)
  {
    for (auto &v : schur)
    {
      v *= scale;  // Schur(A / s) = S / s
    }
    if (!factored)
    {
      // INFOG(21/22): effective memory, max per rank / total (MB); INFOG(29/35):
      // theoretical / effective factor entries.
      auto big = [](MUMPS_INT v) { return v >= 0 ? static_cast<double>(v) : -1.0e6 * v; };
      Mpi::Print(comm,
                 " MUMPS BLR (tol = {:.1e}): factor entries {:.1f}% of full rank, memory "
                 "{:.2f} GB on {:d} rank{} ({:.2f} GB max per rank)\n",
                 blr_tol, 100.0 * big(id.infog[34]) / big(id.infog[28]),
                 id.infog[21] / 1024.0, procs, (procs > 1) ? "s" : "",
                 id.infog[20] / 1024.0);
    }
  }
  factored = true;
}

template <typename T>
void MumpsSchurSolverT<T>::Refactor()
{
  MFEM_VERIFY(refactor, "MumpsSchurSolver was not set up for refactorization!");
  for (auto &v : val)
  {
    v /= scale;
  }
  GatherEntries();
  if (active)
  {
    id.a_loc = Entries(grp != MPI_COMM_NULL ? grp_val.data() : val.data());
  }
  Factor();
}

template <typename T>
MumpsSchurSolverT<T>::~MumpsSchurSolverT()
{
  if (active)
  {
    id.job = -2;
    Call();
  }
  for (MPI_Comm *c : {&sub, &grp})
  {
    if (*c != MPI_COMM_NULL)
    {
      MPI_Comm_free(c);
    }
  }
}

template <typename T>
void MumpsSchurSolverT<T>::GatherEntries()
{
  // The values of the group's entries on its factoring rank.
  if (grp != MPI_COMM_NULL)
  {
    grp_val.resize(active ? grp_disp.back() + grp_cnt.back() : 0);
    MPI_Gatherv(val.data(), static_cast<int>(val.size()), MpiType<T>(), grp_val.data(),
                grp_cnt.data(), grp_disp.data(), MpiType<T>(), 0, grp);
  }
}

template <typename T>
void MumpsSchurSolverT<T>::Gather(const std::vector<const VecType *> &X)
{
  // The right-hand sides on rank 0 (n_glob x |X|, column-major).
  const int nb = static_cast<int>(X.size());
  std::vector<T> loc(n_loc);
  if (rank == 0)
  {
    rhs.assign(static_cast<std::size_t>(n_glob) * nb, T(0.0));
  }
  for (int k = 0; k < nb; k++)
  {
    const VecType &x = *X[k];
    if constexpr (kComplex)
    {
      const double *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
      for (int i = 0; i < n_loc; i++)
      {
        const int r = rows.empty() ? i : rows[i];
        loc[i] = {xr[r], xi[r]};
      }
    }
    else
    {
      const double *xr = x.HostRead();
      for (int i = 0; i < n_loc; i++)
      {
        loc[i] = xr[rows.empty() ? i : rows[i]];
      }
    }
    MPI_Gatherv(loc.data(), n_loc, MpiType<T>(),
                rank == 0 ? rhs.data() + static_cast<std::size_t>(k) * n_glob : nullptr,
                row_cnt.data(), row_disp.data(), MpiType<T>(), 0, comm);
  }
  if (rank == 0)
  {
    id.rhs = Entries(rhs.data());
    id.nrhs = nb;
    id.lrhs = static_cast<MUMPS_INT>(n_glob);
  }
}

template <typename T>
void MumpsSchurSolverT<T>::Scatter(const std::vector<VecType *> &Y)
{
  // The solutions from rank 0, unscaled: (A / s)^-1 b = s A^-1 b.
  std::vector<T> loc(n_loc);
  for (std::size_t k = 0; k < Y.size(); k++)
  {
    MPI_Scatterv(rank == 0 ? rhs.data() + k * static_cast<std::size_t>(n_glob) : nullptr,
                 row_cnt.data(), row_disp.data(), MpiType<T>(), loc.data(), n_loc,
                 MpiType<T>(), 0, comm);
    VecType &y = *Y[k];
    if (rows.empty())
    {
      y.SetSize(n_loc);
    }
    else
    {
      y = 0.0;
    }
    if constexpr (kComplex)
    {
      double *yr = y.Real().HostReadWrite(), *yi = y.Imag().HostReadWrite();
      for (int i = 0; i < n_loc; i++)
      {
        const int r = rows.empty() ? i : rows[i];
        yr[r] = loc[i].real() / scale;
        yi[r] = loc[i].imag() / scale;
      }
    }
    else
    {
      double *yr = y.HostReadWrite();
      for (int i = 0; i < n_loc; i++)
      {
        yr[rows.empty() ? i : rows[i]] = loc[i] / scale;
      }
    }
  }
  rhs = {};
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
    Gather({X.begin() + c0, X.begin() + c0 + nb});
    if (active)
    {
      id.icntl[25] = 0;  // ICNTL(26): the internal problem (Schur variables held at 0)
      id.job = 3;
      Call();
    }
    Check("solve");
    Scatter({Y.begin() + c0, Y.begin() + c0 + nb});
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
  Gather(B);
  if (rank == 0)
  {
    redrhs.assign(static_cast<std::size_t>(n_schur) * std::max(nb, 1), T(0.0));
    id.redrhs = Entries(redrhs.data());
    id.lredrhs = n_schur;
  }
  if (active)
  {
    id.icntl[25] = 1;  // ICNTL(26): condensation onto the Schur variables
    id.job = 3;
    Call();
  }
  Check("condensation");
  rhs = {};
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
    rhs.assign(static_cast<std::size_t>(n_glob) * nb, T(0.0));
    id.rhs = Entries(rhs.data());
    id.nrhs = nb;
    id.lrhs = static_cast<MUMPS_INT>(n_glob);
  }
  if (active)
  {
    id.icntl[25] = 2;  // ICNTL(26): expansion from the Schur part of the solution
    id.job = 3;
    Call();
  }
  Check("expansion");
  redrhs = {};
  reduced = 0;
  Scatter(X);
}

template <typename T>
void MumpsSchurSolverT<T>::Check(const char *phase) const
{
  // INFOG is global (identical on the ranks taking part, rank 0 among them).
  int info[2] = {id.infog[0], id.infog[1]};
  if (grp != MPI_COMM_NULL)
  {
    MPI_Bcast(info, 2, MPI_INT, 0, comm);
  }
  MFEM_VERIFY(info[0] >= 0, "MUMPS " << phase << " failed: INFOG(1) = " << info[0]
                                     << ", INFOG(2) = " << info[1]);
}

template class MumpsSchurSolverT<double>;
template class MumpsSchurSolverT<std::complex<double>>;

}  // namespace palace

#endif
