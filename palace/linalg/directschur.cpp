// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "directschur.hpp"

#include <algorithm>
#include <numeric>
#include <utility>
#include "linalg/superlu.hpp"
#include "utils/communication.hpp"

#if defined(MFEM_USE_STRUMPACK)
#include <StrumpackSparseSolverMPIDist.hpp>
#endif

namespace palace
{

// The blocks of the matrix on the local rows of A_II (rows [off, off + n) of its n_glob,
// numbered by rank): A_II (CSR with global columns), A_IS (CSR, by the positions of the
// Schur variables) and, on rank 0, A_SS (dense, column-major, with ss_slots its structural
// entries).
struct DirectSchurBlocks
{
  using T = std::complex<double>;
  int n, n_schur;
  HYPRE_BigInt off, n_glob;
  const std::vector<int> &ii_ptr;
  const std::vector<HYPRE_BigInt> &ii_col;
  const std::vector<T> &ii_val;
  const std::vector<int> &is_ptr, &is_col;
  const std::vector<T> &is_val;
  const std::vector<int> &ss_slots;
  const std::vector<T> &ss_val;
};

// The factorization of A_II or, partially, of the whole matrix, which then also gives S on
// rank 0 (Factor returns true).
class DirectSchurFactorization
{
public:
  using T = std::complex<double>;
  virtual ~DirectSchurFactorization() = default;
  virtual bool Factor(const DirectSchurBlocks &A, std::vector<T> &schur) = 0;
  // x = A_II^-1 b for nb local column-major blocks (collective).
  virtual void Solve(const std::vector<T> &b, std::vector<T> &x, int nb) const = 0;
};

namespace
{

// The kinds of the entries: of A_II, A_IS and A_SS.
constexpr int kII = 0, kIS = 1, kSS = 2;

// Right-hand sides per solve with A_II (a fixed number: MFEM's SuperLU_DIST wrapper keeps
// that of its first solve).
constexpr int kBatch = 32;

#if defined(MFEM_USE_STRUMPACK)
// STRUMPACK, in complex arithmetic, with the options of Palace's STRUMPACK solver and its
// ordering reused when refactored.
class StrumpackFactorization : public DirectSchurFactorization
{
public:
  explicit StrumpackFactorization(MPI_Comm comm) : comm(comm), sp(comm, false)
  {
    auto &opt = sp.options();
    opt.set_Krylov_solver(strumpack::KrylovSolver::DIRECT);
    opt.set_matching(strumpack::MatchingJob::NONE);
    opt.set_reordering_method(strumpack::ReorderingStrategy::METIS);
    opt.set_compression(strumpack::CompressionType::NONE);
    opt.set_verbose(false);
  }

  bool Factor(const DirectSchurBlocks &A, std::vector<T> &) override
  {
    n_loc = A.n;
    cols.assign(A.ii_col.begin(), A.ii_col.end());
    if (dist.empty())
    {
      const int nranks = Mpi::Size(comm);
      std::vector<int> cnt(nranks);
      MPI_Allgather(&n_loc, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
      dist.assign(nranks + 1, 0);
      for (int r = 0; r < nranks; r++)
      {
        dist[r + 1] = dist[r] + cnt[r];
      }
      sp.set_distributed_csr_matrix(n_loc, A.ii_ptr.data(), cols.data(), A.ii_val.data(),
                                    dist.data(), true);
      MFEM_VERIFY(sp.reorder() == strumpack::ReturnCode::SUCCESS,
                  "STRUMPACK reordering failed!");
    }
    else
    {
      sp.update_matrix_values(n_loc, A.ii_ptr.data(), cols.data(), A.ii_val.data(),
                              dist.data(), true);
    }
    MFEM_VERIFY(sp.factor() == strumpack::ReturnCode::SUCCESS,
                "STRUMPACK factorization failed!");
    return false;
  }

  void Solve(const std::vector<T> &b, std::vector<T> &x, int nb) const override
  {
    x.resize(b.size());
    const int ld = std::max(n_loc, 1);
    MFEM_VERIFY(sp.solve(nb, b.data(), ld, x.data(), ld) == strumpack::ReturnCode::SUCCESS,
                "STRUMPACK solve failed!");
  }

private:
  MPI_Comm comm;
  int n_loc = 0;
  mutable strumpack::SparseSolverMPIDist<T, int> sp;
  std::vector<int> dist, cols;
};
#endif

#if defined(PALACE_WITH_STRUMPACK_SCHUR)
// STRUMPACK's partial factorization of the whole matrix, with its Schur complement
// interface: numbered with the Schur variables first (rows of rank 0), then A_II, it gives
// S, and the expansions of zero values of the Schur variables are the solves with A_II.
class StrumpackSchurFactorization : public DirectSchurFactorization
{
public:
  explicit StrumpackSchurFactorization(MPI_Comm comm) : comm(comm), sp(comm, false)
  {
    // The Schur complement needs no matching (and a sequential graph reordering without
    // compression, STRUMPACK's defaults).
    sp.options().set_matching(strumpack::MatchingJob::NONE);
  }

  bool Factor(const DirectSchurBlocks &A, std::vector<T> &schur) override
  {
    if (dist.empty())
    {
      Pattern(A);
    }
    std::vector<T> si(si_pos.size());
    MPI_Gatherv(A.is_val.data(), static_cast<int>(A.is_val.size()), MPI_C_DOUBLE_COMPLEX,
                si.data(), si_cnt.data(), si_dsp.data(), MPI_C_DOUBLE_COMPLEX, 0, comm);
    val.assign(col.size(), 0.0);
    for (std::size_t e = 0; e < si.size(); e++)
    {
      val[si_pos[e]] = si[e];
    }
    for (std::size_t k = 0; k < ss_pos.size(); k++)
    {
      val[ss_pos[k]] = A.ss_val[A.ss_slots[k]];
    }
    for (int i = 0; i < n; i++)
    {
      auto v = std::copy(A.is_val.begin() + A.is_ptr[i], A.is_val.begin() + A.is_ptr[i + 1],
                         val.begin() + ptr[r0 + i]);
      std::copy(A.ii_val.begin() + A.ii_ptr[i], A.ii_val.begin() + A.ii_ptr[i + 1], v);
    }
    if (first)
    {
      sp.set_distributed_csr_matrix(n_rows, ptr.data(), col.data(), val.data(), dist.data(),
                                    true);
      std::vector<int> vars(n_schur);
      std::iota(vars.begin(), vars.end(), 0);
      MFEM_VERIFY(sp.set_schur_variables(n_schur, vars.data()) ==
                      strumpack::ReturnCode::SUCCESS,
                  "STRUMPACK rejected the Schur variables!");
      first = false;
    }
    else
    {
      sp.update_matrix_values(n_rows, ptr.data(), col.data(), val.data(), dist.data(),
                              true);
    }
    MFEM_VERIFY(sp.factor() == strumpack::ReturnCode::SUCCESS,
                "STRUMPACK partial factorization failed!");
    strumpack::DenseMatrix<T> S;
    MFEM_VERIFY(sp.get_schur_complement(S) == strumpack::ReturnCode::SUCCESS,
                "STRUMPACK Schur complement failed!");
    schur.resize(static_cast<std::size_t>(r0) * r0);
    for (int j = 0; j < r0; j++)
    {
      for (int i = 0; i < r0; i++)
      {
        schur[i + static_cast<std::size_t>(r0) * j] = S(i, j);
      }
    }
    return true;
  }

  void Solve(const std::vector<T> &b, std::vector<T> &x, int nb) const override
  {
    const int ld = std::max(n_rows, 1);
    std::vector<T> bf(static_cast<std::size_t>(ld) * nb, 0.0), xf(bf.size(), 0.0);
    for (int k = 0; k < nb; k++)
    {
      std::copy(b.begin() + static_cast<std::size_t>(n) * k,
                b.begin() + static_cast<std::size_t>(n) * (k + 1),
                bf.begin() + static_cast<std::size_t>(ld) * k + r0);
    }
    MFEM_VERIFY(sp.schur_expand(nb, bf.data(), ld, xf.data(), ld) ==
                    strumpack::ReturnCode::SUCCESS,
                "STRUMPACK solve failed!");
    x.resize(static_cast<std::size_t>(n) * nb);
    for (int k = 0; k < nb; k++)
    {
      std::copy(xf.begin() + static_cast<std::size_t>(ld) * k + r0,
                xf.begin() + static_cast<std::size_t>(ld) * k + r0 + n,
                x.begin() + static_cast<std::size_t>(n) * k);
    }
  }

private:
  // The rows of the Schur variables on rank 0 (the structural entries of A_SS, and those of
  // A_SI, the transposes of the A_IS entries of every rank), and those of A_II (A_IS, then
  // A_II), with the slots of the values.
  void Pattern(const DirectSchurBlocks &A)
  {
    const int rank = Mpi::Rank(comm), nranks = Mpi::Size(comm);
    n = A.n;
    n_schur = A.n_schur;
    r0 = (rank == 0) ? n_schur : 0;
    n_rows = r0 + n;
    std::vector<int> cnt(nranks);
    MPI_Allgather(&n_rows, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    dist.assign(nranks + 1, 0);
    for (int r = 0; r < nranks; r++)
    {
      dist[r + 1] = dist[r] + cnt[r];
    }
    std::vector<int> loc;  // (position, row) of the local A_IS entries
    for (int i = 0; i < n; i++)
    {
      for (int q = A.is_ptr[i]; q < A.is_ptr[i + 1]; q++)
      {
        loc.insert(loc.end(), {A.is_col[q], n_schur + static_cast<int>(A.off) + i});
      }
    }
    const int m = static_cast<int>(loc.size()) / 2;
    si_cnt.assign((rank == 0) ? nranks : 0, 0);
    si_dsp.assign(si_cnt.size(), 0);
    MPI_Gather(&m, 1, MPI_INT, si_cnt.data(), 1, MPI_INT, 0, comm);
    std::vector<int> cnt2(si_cnt.size()), dsp2(si_cnt.size());
    for (std::size_t r = 0; r < si_cnt.size(); r++)
    {
      si_dsp[r] = (r == 0) ? 0 : si_dsp[r - 1] + si_cnt[r - 1];
      cnt2[r] = 2 * si_cnt[r];
      dsp2[r] = 2 * si_dsp[r];
    }
    std::vector<int> all((rank == 0) ? 2 * (si_dsp.back() + si_cnt.back()) : 0);
    MPI_Gatherv(loc.data(), 2 * m, MPI_INT, all.data(), cnt2.data(), dsp2.data(), MPI_INT,
                0, comm);
    ptr.assign(n_rows + 1, 0);
    col.clear();
    if (rank == 0)
    {
      // The entries of each row of a Schur variable: (column, k) for the k-th structural
      // entry of A_SS (k < 0, as -1 - k) and the k-th gathered one of A_SI.
      std::vector<std::vector<std::pair<int, int>>> srow(n_schur);
      for (std::size_t k = 0; k < A.ss_slots.size(); k++)
      {
        srow[A.ss_slots[k] % n_schur].push_back(
            {A.ss_slots[k] / n_schur, -1 - static_cast<int>(k)});
      }
      for (std::size_t k = 0; k < all.size() / 2; k++)
      {
        srow[all[2 * k]].push_back({all[2 * k + 1], static_cast<int>(k)});
      }
      ss_pos.resize(A.ss_slots.size());
      si_pos.resize(all.size() / 2);
      for (int p = 0; p < n_schur; p++)
      {
        std::sort(srow[p].begin(), srow[p].end());
        for (const auto &[c, k] : srow[p])
        {
          ((k < 0) ? ss_pos[-1 - k] : si_pos[k]) = static_cast<int>(col.size());
          col.push_back(c);
        }
        ptr[p + 1] = static_cast<int>(col.size());
      }
    }
    // A_IS's columns (the Schur positions) come before A_II's (after the Schur variables).
    for (int i = 0; i < n; i++)
    {
      col.insert(col.end(), A.is_col.begin() + A.is_ptr[i],
                 A.is_col.begin() + A.is_ptr[i + 1]);
      for (int q = A.ii_ptr[i]; q < A.ii_ptr[i + 1]; q++)
      {
        col.push_back(n_schur + static_cast<int>(A.ii_col[q]));
      }
      ptr[r0 + i + 1] = static_cast<int>(col.size());
    }
  }

  MPI_Comm comm;
  int n = 0, n_schur = 0, r0 = 0, n_rows = 0;
  bool first = true;
  std::vector<int> dist, ptr, col, si_cnt, si_dsp, si_pos, ss_pos;
  std::vector<T> val;
  mutable strumpack::SparseSolverMPIDist<T, int> sp;
};
#endif

#if defined(MFEM_USE_SUPERLU)
// SuperLU_DIST, in the real equivalent form [Ar, Ai; Ai, -Ar] of A_II (local rows: those of
// Ar, then of Ai, so that a local vector is [Re x; -Im x] for the solution x).
class RealFormFactorization : public DirectSchurFactorization
{
public:
  explicit RealFormFactorization(MPI_Comm comm)
    : comm(comm), solver(std::make_unique<SuperLUSolver>(
                      comm, SymbolicFactorization::DEFAULT, false, true, 0))
  {
  }

  bool Factor(const DirectSchurBlocks &A, std::vector<T> &) override
  {
    n_loc = A.n;
    std::vector<double> re(A.ii_val.size()), im(A.ii_val.size());
    for (std::size_t q = 0; q < A.ii_val.size(); q++)
    {
      re[q] = A.ii_val[q].real();
      im[q] = A.ii_val[q].imag();
    }
    HYPRE_BigInt part[2] = {A.off, A.off + A.n};
    mfem::HypreParMatrix Ar(comm, A.n, A.n_glob, A.n_glob, A.ii_ptr.data(), A.ii_col.data(),
                            re.data(), part, part);
    mfem::HypreParMatrix Ai(comm, A.n, A.n_glob, A.n_glob, A.ii_ptr.data(), A.ii_col.data(),
                            im.data(), part, part);
    mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
    mfem::Array2D<double> coef(2, 2);
    blocks(0, 0) = blocks(1, 1) = &Ar;
    blocks(0, 1) = blocks(1, 0) = &Ai;
    coef(0, 0) = coef(0, 1) = coef(1, 0) = 1.0;
    coef(1, 1) = -1.0;
    A2.reset(mfem::HypreParMatrixFromBlocks(blocks, &coef));
    solver->SetOperator(*A2);
    return false;
  }

  void Solve(const std::vector<T> &b, std::vector<T> &x, int nb) const override
  {
    x.assign(static_cast<std::size_t>(n_loc) * nb, 0.0);
    std::vector<Vector> xv(kBatch), yv(kBatch);
    mfem::Array<const Vector *> X(kBatch);
    mfem::Array<Vector *> Y(kBatch);
    for (int k = 0; k < kBatch; k++)
    {
      xv[k].SetSize(2 * n_loc);
      yv[k].SetSize(2 * n_loc);
      X[k] = &xv[k];
      Y[k] = &yv[k];
    }
    for (int c0 = 0; c0 < nb; c0 += kBatch)
    {
      const int m = std::min(kBatch, nb - c0);
      for (int k = 0; k < kBatch; k++)
      {
        double *p = xv[k].HostWrite();
        for (int i = 0; i < n_loc; i++)
        {
          const T v = (k < m) ? b[i + static_cast<std::size_t>(n_loc) * (c0 + k)] : T(0.0);
          p[i] = v.real();
          p[n_loc + i] = v.imag();
        }
      }
      solver->ArrayMult(X, Y);
      for (int k = 0; k < m; k++)
      {
        const double *p = yv[k].HostRead();
        for (int i = 0; i < n_loc; i++)
        {
          x[i + static_cast<std::size_t>(n_loc) * (c0 + k)] = T(p[i], -p[n_loc + i]);
        }
      }
    }
  }

private:
  MPI_Comm comm;
  int n_loc = 0;
  std::unique_ptr<mfem::Solver> solver;
  std::unique_ptr<mfem::HypreParMatrix> A2;
};
#endif

}  // namespace

DirectSchurSolver::DirectSchurSolver(MPI_Comm comm, HYPRE_BigInt n_glob, int n_loc, Coo &&A,
                                     const std::vector<HYPRE_BigInt> &schur_vars,
                                     std::vector<int> local_rows, LinearSolver type)
  : comm(comm), rank(Mpi::Rank(comm)), n_loc(n_loc),
    n_schur(static_cast<int>(schur_vars.size())), rows(std::move(local_rows)),
    val(std::move(A.val))
{
  switch (type)
  {
#if defined(MFEM_USE_STRUMPACK)
    case LinearSolver::STRUMPACK:
#if defined(PALACE_WITH_STRUMPACK_SCHUR)
      if (n_schur > 0)
      {
        fact = std::make_unique<StrumpackSchurFactorization>(comm);
        break;
      }
#endif
      fact = std::make_unique<StrumpackFactorization>(comm);
      break;
#endif
#if defined(MFEM_USE_SUPERLU)
    case LinearSolver::SUPERLU:
      fact = std::make_unique<RealFormFactorization>(comm);
      break;
#endif
    default:
      MFEM_ABORT("DirectSchurSolver factors with STRUMPACK or SuperLU_DIST, and Palace is "
                 "not built with the one requested!");
  }
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
      ss_slots.push_back(recv_slot[p]);
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
  std::sort(ss_slots.begin(), ss_slots.end());
  ss_slots.erase(std::unique(ss_slots.begin(), ss_slots.end()), ss_slots.end());
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
  if (fact->Factor({n_int, n_schur, int_off, n_int_glob, ii_ptr, ii_col, ii_val, is_ptr,
                    is_col, is_val, ss_slots, ss_val},
                   schur))
  {
    return;
  }

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
  fact->Solve(b, x, nb);
}

void DirectSchurSolver::SolveInternal(const std::vector<const ComplexVector *> &X,
                                      const std::vector<ComplexVector *> &Y)
{
  const int nb = static_cast<int>(X.size());
  std::vector<T> b(static_cast<std::size_t>(n_int) * nb, 0.0), x;
  for (int k = 0; k < nb; k++)
  {
    const double *br = X[k]->Real().HostRead(), *bi = X[k]->Imag().HostRead();
    for (int q = 0; q < n_loc; q++)
    {
      if (rows[q] >= 0 && int_row[q] >= 0)
      {
        b[int_row[q] + static_cast<std::size_t>(n_int) * k] = T(br[rows[q]], bi[rows[q]]);
      }
    }
  }
  SolveInterior(b, x, nb);
  for (int k = 0; k < nb; k++)
  {
    ComplexVector &y = *Y[k];
    y = 0.0;
    double *yr = y.Real().HostReadWrite(), *yi = y.Imag().HostReadWrite();
    for (int q = 0; q < n_loc; q++)
    {
      if (rows[q] >= 0 && int_row[q] >= 0)
      {
        const T v = x[int_row[q] + static_cast<std::size_t>(n_int) * k];
        yr[rows[q]] = v.real();
        yi[rows[q]] = v.imag();
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
