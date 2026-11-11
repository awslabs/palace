// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "directschur.hpp"

#include <algorithm>
#include <utility>
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// The kinds of the entries: of A_II, A_IS and A_SS.
constexpr int kII = 0, kIS = 1, kSS = 2;

// Right-hand sides per solve with A_II (a fixed number: MFEM's SuperLU_DIST wrapper keeps
// that of its first solve).
constexpr int kBatch = 32;

}  // namespace

DirectSchurSolver::DirectSchurSolver(MPI_Comm comm, HYPRE_BigInt n_glob, int n_loc, Coo &&A,
                                     const std::vector<HYPRE_BigInt> &schur_vars,
                                     std::vector<int> local_rows,
                                     std::unique_ptr<mfem::Solver> &&direct_solver)
  : comm(comm), rank(Mpi::Rank(comm)), n_loc(n_loc),
    n_schur(static_cast<int>(schur_vars.size())), rows(std::move(local_rows)),
    val(std::move(A.val)), solver(std::move(direct_solver))
{
  MFEM_VERIFY(solver, "DirectSchurSolver needs a sparse direct solver!");
  MFEM_VERIFY(rows.size() == static_cast<std::size_t>(n_loc),
              "DirectSchurSolver needs the vector index of each local row!");
  const int nranks = Mpi::Size(comm);
  std::vector<HYPRE_BigInt> row_disp(nranks + 1, 0);
  {
    std::vector<long long> cnt(nranks);
    const long long m = n_loc;
    MPI_Allgather(&m, 1, MPI_LONG_LONG, cnt.data(), 1, MPI_LONG_LONG, comm);
    for (int r = 0; r < nranks; r++)
    {
      row_disp[r + 1] = row_disp[r] + static_cast<HYPRE_BigInt>(cnt[r]);
    }
    MFEM_VERIFY(row_disp.back() == n_glob, "DirectSchurSolver: wrong number of rows!");
  }

  // The Schur variables by index, with their positions; A_II's rows and columns are the
  // other indices, in order.
  std::vector<std::pair<HYPRE_BigInt, int>> sorted(n_schur);
  for (int p = 0; p < n_schur; p++)
  {
    sorted[p] = {schur_vars[p], p};
  }
  std::sort(sorted.begin(), sorted.end());
  std::vector<HYPRE_BigInt> s_index(n_schur);
  std::vector<int> s_pos(n_schur);
  for (int p = 0; p < n_schur; p++)
  {
    s_index[p] = sorted[p].first;
    s_pos[p] = sorted[p].second;
  }
  auto pos = [&](HYPRE_BigInt g)
  {
    const auto it = std::lower_bound(s_index.begin(), s_index.end(), g);
    return (it != s_index.end() && *it == g) ? s_pos[it - s_index.begin()] : -1;
  };
  auto interior = [&](HYPRE_BigInt g)
  { return g - (std::lower_bound(s_index.begin(), s_index.end(), g) - s_index.begin()); };
  auto owner = [&](HYPRE_BigInt g)
  {
    return static_cast<int>(std::upper_bound(row_disp.begin(), row_disp.end(), g) -
                            row_disp.begin()) -
           1;
  };
  schur_pos.assign(n_loc, -1);
  int_row.assign(n_loc, -1);
  for (int q = 0; q < n_loc; q++)
  {
    const int p = pos(row_disp[rank] + q);
    if (p >= 0)
    {
      schur_pos[q] = p;
    }
    else
    {
      int_row[q] = n_int++;
    }
  }
  int_off = interior(row_disp[rank]);
  n_int_glob = n_glob - n_schur;

  // Each entry and its transpose go to the rank of their row, those of A_SS to rank 0 (an
  // entry of A_SI is the transpose of one of A_IS).
  std::vector<std::vector<long long>> out(nranks);
  std::vector<std::vector<int>> src(nranks);
  auto add = [&](HYPRE_BigInt i, HYPRE_BigInt j, int e)
  {
    const int pi = pos(i), pj = pos(j);
    if (pi < 0)
    {
      const int r = owner(i);
      out[r].insert(out[r].end(), {(pj < 0) ? kII : kIS, i, (pj < 0) ? j : pj});
      src[r].push_back(e);
    }
    else if (pj >= 0)
    {
      out[0].insert(out[0].end(), {kSS, pi, pj});
      src[0].push_back(e);
    }
  };
  for (std::size_t e = 0; e < A.irn.size(); e++)
  {
    const HYPRE_BigInt i = A.irn[e] - 1, j = A.jcn[e] - 1;
    add(i, j, static_cast<int>(e));
    if (i != j)
    {
      add(j, i, static_cast<int>(e));
    }
  }
  send_cnt.assign(nranks, 0);
  send_dsp.assign(nranks, 0);
  recv_cnt.assign(nranks, 0);
  recv_dsp.assign(nranks, 0);
  std::vector<long long> sidx;
  for (int r = 0; r < nranks; r++)
  {
    send_cnt[r] = static_cast<int>(src[r].size());
    send_dsp[r] = (r == 0) ? 0 : send_dsp[r - 1] + send_cnt[r - 1];
    send_src.insert(send_src.end(), src[r].begin(), src[r].end());
    sidx.insert(sidx.end(), out[r].begin(), out[r].end());
  }
  MPI_Alltoall(send_cnt.data(), 1, MPI_INT, recv_cnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    recv_dsp[r] = recv_dsp[r - 1] + recv_cnt[r - 1];
  }
  const int nrecv = recv_dsp.back() + recv_cnt.back();
  std::vector<long long> ridx(3 * static_cast<std::size_t>(nrecv));
  {
    std::vector<int> c3(nranks), d3(nranks), rc3(nranks), rd3(nranks);
    for (int r = 0; r < nranks; r++)
    {
      c3[r] = 3 * send_cnt[r];
      d3[r] = 3 * send_dsp[r];
      rc3[r] = 3 * recv_cnt[r];
      rd3[r] = 3 * recv_dsp[r];
    }
    MPI_Alltoallv(sidx.data(), c3.data(), d3.data(), MPI_LONG_LONG, ridx.data(), rc3.data(),
                  rd3.data(), MPI_LONG_LONG, comm);
  }

  // The local patterns of A_II and A_IS (CSR, columns sorted), and the slots of the
  // received entries in their values or in A_SS.
  recv_kind.resize(nrecv);
  recv_slot.resize(nrecv);
  std::vector<std::pair<std::pair<int, long long>, int>> ii, is;
  for (int p = 0; p < nrecv; p++)
  {
    const int kind = static_cast<int>(ridx[3 * p]);
    const long long a = ridx[3 * p + 1], b = ridx[3 * p + 2];
    recv_kind[p] = kind;
    if (kind == kSS)
    {
      recv_slot[p] = static_cast<int>(a + static_cast<long long>(n_schur) * b);
    }
    else
    {
      const int li = static_cast<int>(interior(a) - int_off);
      (kind == kII ? ii : is).push_back({{li, (kind == kII) ? interior(b) : b}, p});
    }
  }
  auto pattern = [&](auto &entries, std::vector<int> &ptr, auto &col)
  {
    std::sort(entries.begin(), entries.end());
    ptr.assign(n_int + 1, 0);
    col.clear();
    for (std::size_t q = 0; q < entries.size(); q++)
    {
      const auto &[ij, p] = entries[q];
      if (q == 0 || ij != entries[q - 1].first)
      {
        col.push_back(
            static_cast<typename std::decay_t<decltype(col)>::value_type>(ij.second));
        ptr[ij.first + 1]++;
      }
      recv_slot[p] = static_cast<int>(col.size()) - 1;
    }
    for (int i = 0; i < n_int; i++)
    {
      ptr[i + 1] += ptr[i];
    }
  };
  pattern(ii, ii_ptr, ii_col);
  pattern(is, is_ptr, is_col);
  Factor();
}

DirectSchurSolver::~DirectSchurSolver() = default;

void DirectSchurSolver::Distribute()
{
  std::vector<T> sv(send_src.size()), rv(recv_kind.size());
  for (std::size_t q = 0; q < send_src.size(); q++)
  {
    sv[q] = val[send_src[q]];
  }
  MPI_Alltoallv(sv.data(), send_cnt.data(), send_dsp.data(), MPI_C_DOUBLE_COMPLEX,
                rv.data(), recv_cnt.data(), recv_dsp.data(), MPI_C_DOUBLE_COMPLEX, comm);
  ii_val.assign(ii_col.size(), 0.0);
  is_val.assign(is_col.size(), 0.0);
  ss_val.assign((rank == 0) ? static_cast<std::size_t>(n_schur) * n_schur : 0, 0.0);
  for (std::size_t p = 0; p < rv.size(); p++)
  {
    auto &v = (recv_kind[p] == kII) ? ii_val : (recv_kind[p] == kIS) ? is_val : ss_val;
    v[recv_slot[p]] += rv[p];
  }
}

void DirectSchurSolver::Factor()
{
  Distribute();
  {
    // A_II in its real equivalent form [Ar, Ai; Ai, -Ar] (local rows: those of Ar, then of
    // Ai, so that a local vector is [Re x; -Im x] for the solution x).
    std::vector<double> re(ii_val.size()), im(ii_val.size());
    for (std::size_t q = 0; q < ii_val.size(); q++)
    {
      re[q] = ii_val[q].real();
      im[q] = ii_val[q].imag();
    }
    HYPRE_BigInt part[2] = {int_off, int_off + n_int};
    mfem::HypreParMatrix Ar(comm, n_int, n_int_glob, n_int_glob, ii_ptr.data(),
                            ii_col.data(), re.data(), part, part);
    mfem::HypreParMatrix Ai(comm, n_int, n_int_glob, n_int_glob, ii_ptr.data(),
                            ii_col.data(), im.data(), part, part);
    mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
    mfem::Array2D<double> coef(2, 2);
    blocks(0, 0) = blocks(1, 1) = &Ar;
    blocks(0, 1) = blocks(1, 0) = &Ai;
    coef(0, 0) = coef(0, 1) = coef(1, 0) = 1.0;
    coef(1, 1) = -1.0;
    A2.reset(mfem::HypreParMatrixFromBlocks(blocks, &coef));
  }
  solver->SetOperator(*A2);

  // S = A_SS - A_SI A_II^-1 A_IS, for a batch of the Schur variables at a time.
  schur = std::move(ss_val);
  for (int c0 = 0; c0 < n_schur; c0 += kBatch)
  {
    const int nb = std::min(kBatch, n_schur - c0);
    std::vector<T> b(static_cast<std::size_t>(n_int) * nb, 0.0), x;
    for (int i = 0; i < n_int; i++)
    {
      for (int q = is_ptr[i]; q < is_ptr[i + 1]; q++)
      {
        if (is_col[q] >= c0 && is_col[q] < c0 + nb)
        {
          b[i + static_cast<std::size_t>(n_int) * (is_col[q] - c0)] = is_val[q];
        }
      }
    }
    SolveInterior(b, x, nb);
    AddSchurProduct(x, nb, -1.0,
                    (rank == 0) ? schur.data() + static_cast<std::size_t>(n_schur) * c0
                                : nullptr);
  }
}

void DirectSchurSolver::Refactor()
{
  MFEM_VERIFY(!reduced, "DirectSchurSolver: Refactor between Reduce and Expand!");
  Factor();
}

std::vector<DirectSchurSolver::T> DirectSchurSolver::ReleaseSchur()
{
  std::vector<T> s;
  s.swap(schur);
  return s;
}

void DirectSchurSolver::SolveInterior(const std::vector<T> &b, std::vector<T> &x,
                                      int nb) const
{
  x.assign(static_cast<std::size_t>(n_int) * nb, 0.0);
  std::vector<Vector> xv(kBatch), yv(kBatch);
  mfem::Array<const Vector *> X(kBatch);
  mfem::Array<Vector *> Y(kBatch);
  for (int k = 0; k < kBatch; k++)
  {
    xv[k].SetSize(2 * n_int);
    yv[k].SetSize(2 * n_int);
    X[k] = &xv[k];
    Y[k] = &yv[k];
  }
  for (int c0 = 0; c0 < nb; c0 += kBatch)
  {
    const int m = std::min(kBatch, nb - c0);
    for (int k = 0; k < kBatch; k++)
    {
      double *p = xv[k].HostWrite();
      for (int i = 0; i < n_int; i++)
      {
        const T v = (k < m) ? b[i + static_cast<std::size_t>(n_int) * (c0 + k)] : T(0.0);
        p[i] = v.real();
        p[n_int + i] = v.imag();
      }
    }
    solver->ArrayMult(X, Y);
    for (int k = 0; k < m; k++)
    {
      const double *p = yv[k].HostRead();
      for (int i = 0; i < n_int; i++)
      {
        x[i + static_cast<std::size_t>(n_int) * (c0 + k)] = T(p[i], -p[n_int + i]);
      }
    }
  }
}

void DirectSchurSolver::AddSchurProduct(const std::vector<T> &x, int nb, double sign,
                                        T *y) const
{
  // A_SI x = A_IS^T x (A is symmetric), summed over the ranks onto rank 0.
  std::vector<T> loc(static_cast<std::size_t>(n_schur) * nb, 0.0), sum;
  for (int i = 0; i < n_int; i++)
  {
    for (int q = is_ptr[i]; q < is_ptr[i + 1]; q++)
    {
      for (int k = 0; k < nb; k++)
      {
        loc[is_col[q] + static_cast<std::size_t>(n_schur) * k] +=
            is_val[q] * x[i + static_cast<std::size_t>(n_int) * k];
      }
    }
  }
  if (rank == 0)
  {
    sum.resize(loc.size());
  }
  MPI_Reduce(loc.data(), sum.data(), static_cast<int>(loc.size()), MPI_C_DOUBLE_COMPLEX,
             MPI_SUM, 0, comm);
  if (rank == 0)
  {
    for (std::size_t q = 0; q < sum.size(); q++)
    {
      y[q] += sign * sum[q];
    }
  }
}

void DirectSchurSolver::Reduce(const std::vector<const ComplexVector *> &B,
                               std::vector<T> &red)
{
  MFEM_VERIFY(n_schur > 0 && !reduced, "DirectSchurSolver: invalid Reduce!");
  const int nb = static_cast<int>(B.size());
  std::vector<T> b(static_cast<std::size_t>(n_int) * nb, 0.0),
      bs(static_cast<std::size_t>(n_schur) * nb, 0.0);
  for (int k = 0; k < nb; k++)
  {
    const double *br = B[k]->Real().HostRead(), *bi = B[k]->Imag().HostRead();
    for (int q = 0; q < n_loc; q++)
    {
      const int r = rows[q];
      if (r < 0)
      {
        continue;
      }
      const T v(br[r], bi[r]);
      if (int_row[q] >= 0)
      {
        b[int_row[q] + static_cast<std::size_t>(n_int) * k] = v;
      }
      else
      {
        bs[schur_pos[q] + static_cast<std::size_t>(n_schur) * k] = v;
      }
    }
  }
  SolveInterior(b, y_last, nb);
  // red = b_S - A_SI A_II^-1 b_I, with each b_S from the rank of its row.
  red.assign((rank == 0) ? bs.size() : 0, 0.0);
  MPI_Reduce(bs.data(), red.data(), static_cast<int>(bs.size()), MPI_C_DOUBLE_COMPLEX,
             MPI_SUM, 0, comm);
  AddSchurProduct(y_last, nb, -1.0, red.data());
  reduced = nb;
}

void DirectSchurSolver::Expand(const std::vector<T> &u,
                               const std::vector<ComplexVector *> &X)
{
  MFEM_VERIFY(reduced && static_cast<int>(X.size()) == reduced,
              "DirectSchurSolver: Expand must follow Reduce of the same batch!");
  const int nb = reduced;
  std::vector<T> us(static_cast<std::size_t>(n_schur) * nb);
  if (rank == 0)
  {
    MFEM_VERIFY(u.size() == us.size(), "DirectSchurSolver: wrong size of u_S!");
    us = u;
  }
  MPI_Bcast(us.data(), static_cast<int>(us.size()), MPI_C_DOUBLE_COMPLEX, 0, comm);
  // x_I = A_II^-1 b_I - A_II^-1 A_IS u_S.
  std::vector<T> r(static_cast<std::size_t>(n_int) * nb, 0.0), z;
  for (int i = 0; i < n_int; i++)
  {
    for (int q = is_ptr[i]; q < is_ptr[i + 1]; q++)
    {
      for (int k = 0; k < nb; k++)
      {
        r[i + static_cast<std::size_t>(n_int) * k] +=
            is_val[q] * us[is_col[q] + static_cast<std::size_t>(n_schur) * k];
      }
    }
  }
  SolveInterior(r, z, nb);
  for (int k = 0; k < nb; k++)
  {
    ComplexVector &x = *X[k];
    x = 0.0;
    double *xr = x.Real().HostReadWrite(), *xi = x.Imag().HostReadWrite();
    for (int q = 0; q < n_loc; q++)
    {
      const int row = rows[q];
      if (row < 0)
      {
        continue;
      }
      const T v = (int_row[q] >= 0)
                      ? y_last[int_row[q] + static_cast<std::size_t>(n_int) * k] -
                            z[int_row[q] + static_cast<std::size_t>(n_int) * k]
                      : us[schur_pos[q] + static_cast<std::size_t>(n_schur) * k];
      xr[row] = v.real();
      xi[row] = v.imag();
    }
  }
  std::vector<T>().swap(y_last);
  reduced = 0;
}

}  // namespace palace
