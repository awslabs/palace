// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hypre.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include "linalg/fusedvector.hpp"

namespace palace::hypre
{

HypreVector::HypreVector(hypre_Vector *vec) : vec(vec) {}

HypreVector::HypreVector(const Vector &x) : vec(nullptr)
{
  Update(x);
}

HypreVector::~HypreVector()
{
  hypre_SeqVectorDestroy(vec);
}

void HypreVector::Update(const Vector &x)
{
  const HYPRE_Int N = x.Size();
  if (!vec)
  {
    vec = hypre_SeqVectorCreate(N);
    hypre_SeqVectorSetDataOwner(vec, 0);
    hypre_VectorData(vec) = const_cast<double *>(x.Read());
    hypre_SeqVectorInitialize(vec);
  }
  else
  {
    hypre_SeqVectorSetSize(vec, N);
    hypre_VectorData(vec) = const_cast<double *>(x.Read());
  }
}

HypreCSRMatrix::HypreCSRMatrix(int h, int w, int nnz)
  : palace::Operator(h, w), hypre_own_I(true)
{
  mat = hypre_CSRMatrixCreate(h, w, nnz);
  hypre_CSRMatrixInitialize(mat);
}

HypreCSRMatrix::HypreCSRMatrix(hypre_CSRMatrix *mat) : mat(mat), hypre_own_I(true)
{
  height = hypre_CSRMatrixNumRows(mat);
  width = hypre_CSRMatrixNumCols(mat);
}

HypreCSRMatrix::HypreCSRMatrix(const mfem::SparseMatrix &m)
  : palace::Operator(m.Height(), m.Width()), hypre_own_I(false)
{
  const int nnz = m.NumNonZeroElems();
  mat = hypre_CSRMatrixCreate(height, width, nnz);
  hypre_CSRMatrixSetDataOwner(mat, 0);
  hypre_CSRMatrixData(mat) = const_cast<double *>(m.ReadData());
#if !defined(HYPRE_BIGINT)
  hypre_CSRMatrixI(mat) = const_cast<int *>(m.ReadI());
  hypre_CSRMatrixJ(mat) = const_cast<int *>(m.ReadJ());
#else
  data_I.SetSize(height);
  data_J.SetSize(nnz);
  {
    const auto *I = m.ReadI();
    const auto *J = m.ReadJ();
    auto *DI = data_I.Write();
    auto *DJ = data_J.Write();
    mfem::forall(height, [=] MFEM_HOST_DEVICE(int i) { DI[i] = I[i]; });
    mfem::forall(nnz, [=] MFEM_HOST_DEVICE(int i) { DJ[i] = J[i]; });
  }
#endif
  hypre_CSRMatrixInitialize(mat);
}

HypreCSRMatrix::~HypreCSRMatrix()
{
  if (!hypre_own_I)
  {
    hypre_CSRMatrixI(mat) = nullptr;
  }
  hypre_CSRMatrixDestroy(mat);
}

void HypreCSRMatrix::AssembleDiagonal(Vector &diag) const
{
  diag.SetSize(height);
  hypre_CSRMatrixExtractDiagonal(mat, diag.Write(), 0);
}

namespace
{

static HypreVector X, Y;

}  // namespace

void HypreCSRMatrix::Mult(const Vector &x, Vector &y) const
{
  X.Update(x);
  Y.Update(y);
  hypre_CSRMatrixMatvec(1.0, mat, X, 0.0, Y);
}

void HypreCSRMatrix::AddMult(const Vector &x, Vector &y, const double a) const
{
  X.Update(x);
  Y.Update(y);
  hypre_CSRMatrixMatvec(a, mat, X, 1.0, Y);
}

void HypreCSRMatrix::MultTranspose(const Vector &x, Vector &y) const
{
  X.Update(x);
  Y.Update(y);
  hypre_CSRMatrixMatvecT(1.0, mat, X, 0.0, Y);
}

void HypreCSRMatrix::AddMultTranspose(const Vector &x, Vector &y, const double a) const
{
  X.Update(x);
  Y.Update(y);
  hypre_CSRMatrixMatvecT(a, mat, X, 1.0, Y);
}

namespace
{

// Transpose of a CSR block, with the entries of each row in increasing column order (as
// hypre_CSRMatrixTranspose).
void TransposeCSR(int nrows, int ncols, const HYPRE_Int *ia, const HYPRE_Int *ja,
                  const HYPRE_Complex *a, std::vector<int> &iat, std::vector<int> &jat,
                  std::vector<double> &at)
{
  iat.assign(ncols + 1, 0);
  for (int k = 0; k < ia[nrows]; k++)
  {
    iat[ja[k] + 1]++;
  }
  for (int j = 0; j < ncols; j++)
  {
    iat[j + 1] += iat[j];
  }
  jat.resize(ia[nrows]);
  at.resize(ia[nrows]);
  std::vector<int> pos(iat.begin(), iat.end() - 1);
  for (int i = 0; i < nrows; i++)
  {
    for (int k = ia[i]; k < ia[i + 1]; k++)
    {
      const int p = pos[ja[k]]++;
      jat[p] = i;
      at[p] = a[k];
    }
  }
}

// Identity rows (a single unit entry) of a CSR matrix as runs y[row + k] = x[col + k], and
// the other rows. Rows listed in skip (sorted) are always other rows.
template <typename IntType, typename ScalarType>
void FindIdentityRuns(int nrows, const IntType *ia, const IntType *ja, const ScalarType *a,
                      const std::vector<int> &skip, std::vector<int> &run_row,
                      std::vector<int> &run_col, std::vector<int> &run_len,
                      std::vector<int> &rows)
{
  std::size_t s = 0;
  for (int i = 0; i < nrows; i++)
  {
    const bool skipped = (s < skip.size() && skip[s] == i);
    s += skipped;
    if (!skipped && ia[i + 1] - ia[i] == 1 && a[ia[i]] == 1.0)
    {
      const int j = ja[ia[i]];
      if (!run_len.empty() && run_row.back() + run_len.back() == i &&
          run_col.back() + run_len.back() == j)
      {
        run_len.back()++;
      }
      else
      {
        run_row.push_back(i);
        run_col.push_back(j);
        run_len.push_back(1);
      }
    }
    else
    {
      rows.push_back(i);
    }
  }
}

// y[row + k] = 0 + x[col + k] over the runs (the result of hypre for a unit entry).
void CopyRuns(const std::vector<int> &run_row, const std::vector<int> &run_col,
              const std::vector<int> &run_len, const double *X, double *Y)
{
  for (std::size_t r = 0; r < run_len.size(); r++)
  {
    const double *__restrict__ Xr = X + run_col[r];
    double *__restrict__ Yr = Y + run_row[r];
    const int len = run_len[r];
    for (int k = 0; k < len; k++)
    {
      Yr[k] = 0.0 + Xr[k];
    }
  }
}

// For each listed row of a transposed block, y[i] = Σ a x accumulated with fused multiply-
// adds in increasing order of the original rows, as hypre's transposed (scatter) CSR
// product on the host. The values at[k] are given by a functor.
template <typename Values>
void FmaRowSumsV(const std::vector<int> &rows, const std::vector<int> &iat,
                 const std::vector<int> &jat, const Values &at, const double *X, double *Y)
{
  for (const int i : rows)
  {
    double t = 0.0;
    for (int k = iat[i]; k < iat[i + 1]; k++)
    {
      t = linalg::fused::Fma(at(k), X[jat[k]], t);
    }
    Y[i] = t;
  }
}

template <typename Values>
void FmaRowSumsV(const std::vector<int> &rows, const std::vector<int> &iat,
                 const std::vector<int> &jat, const Values &at, const double *XR,
                 const double *XI, double *YR, double *YI)
{
  for (const int i : rows)
  {
    double tr = 0.0, ti = 0.0;
    for (int k = iat[i]; k < iat[i + 1]; k++)
    {
      const double ak = at(k);
      tr = linalg::fused::Fma(ak, XR[jat[k]], tr);
      ti = linalg::fused::Fma(ak, XI[jat[k]], ti);
    }
    YR[i] = tr;
    YI[i] = ti;
  }
}

void FmaRowSums(const std::vector<int> &rows, const std::vector<int> &iat,
                const std::vector<int> &jat, const std::vector<double> &at, const double *X,
                double *Y)
{
  FmaRowSumsV(rows, iat, jat, [&at](int k) { return at[k]; }, X, Y);
}

void CopyRuns(const std::vector<int> &run_row, const std::vector<int> &run_col,
              const std::vector<int> &run_len, const double *XR, const double *XI,
              double *YR, double *YI)
{
  CopyRuns(run_row, run_col, run_len, XR, YR);
  CopyRuns(run_row, run_col, run_len, XI, YI);
}

void FmaRowSums(const std::vector<int> &rows, const std::vector<int> &iat,
                const std::vector<int> &jat, const std::vector<double> &at,
                const double *XR, const double *XI, double *YR, double *YI)
{
  FmaRowSumsV(rows, iat, jat, [&at](int k) { return at[k]; }, XR, XI, YR, YI);
}

hypre_ParCSRCommPkg *GetCommPkg(hypre_ParCSRMatrix *A, int vecstride)
{
  if (!hypre_ParCSRMatrixCommPkg(A))
  {
    hypre_MatvecCommPkgCreate(A);
  }
  hypre_ParCSRCommPkg *comm_pkg = hypre_ParCSRMatrixCommPkg(A);
  hypre_ParCSRCommPkgUpdateVecStarts(comm_pkg, 1, vecstride, 1);
  return comm_pkg;
}

}  // namespace

HypreProlongation::HypreProlongation(const mfem::HypreParMatrix &P)
  : palace::Operator(P.Height(), P.Width()), P(P)
{
  hypre_ParCSRMatrix *A = P;
  hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(A), *offd = hypre_ParCSRMatrixOffd(A);
  const HYPRE_Int *Di = hypre_CSRMatrixI(diag), *Dj = hypre_CSRMatrixJ(diag);
  const HYPRE_Complex *Da = hypre_CSRMatrixData(diag);
  const HYPRE_Int *Oi = hypre_CSRMatrixI(offd), *Oj = hypre_CSRMatrixJ(offd);
  const HYPRE_Complex *Oa = hypre_CSRMatrixData(offd);
  const int nrows = hypre_CSRMatrixNumRows(diag), ncols = hypre_CSRMatrixNumCols(diag);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  if (ncols_offd > 0)
  {
    for (int i = 0; i < nrows; i++)
    {
      if (Oi[i + 1] > Oi[i])
      {
        offd_rows.push_back(i);
      }
    }
    TransposeCSR(nrows, ncols_offd, Oi, Oj, Oa, offdT_i, offdT_j, offdT_a);
    offdT_rows.resize(ncols_offd);
    for (int c = 0; c < ncols_offd; c++)
    {
      offdT_rows[c] = c;
    }
  }
  FindIdentityRuns(nrows, Di, Dj, Da, offd_rows, run_row, run_col, run_len, diag_rows);
  TransposeCSR(nrows, ncols, Di, Dj, Da, diagT_i, diagT_j, diagT_a);
  FindIdentityRuns(ncols, diagT_i.data(), diagT_j.data(), diagT_a.data(), {}, runT_row,
                   runT_col, runT_len, diagT_rows);

  // Compact copy of the remaining rows of the diagonal block, and the values as indices
  // into the table of the distinct values (by bit pattern) when it is small enough.
  diag_ia.reserve(diag_rows.size() + 1);
  diag_ia.push_back(0);
  for (const int i : diag_rows)
  {
    for (int k = Di[i]; k < Di[i + 1]; k++)
    {
      diag_ja.push_back(Dj[k]);
      diag_a.push_back(Da[k]);
    }
    diag_ia.push_back(static_cast<int>(diag_ja.size()));
  }
  auto Bits = [](double a)
  {
    std::uint64_t b;
    std::memcpy(&b, &a, sizeof(b));
    return b;
  };
  std::vector<std::uint64_t> bits(diagT_a.size());
  for (std::size_t k = 0; k < diagT_a.size(); k++)
  {
    bits[k] = Bits(diagT_a[k]);
  }
  std::sort(bits.begin(), bits.end());
  bits.erase(std::unique(bits.begin(), bits.end()), bits.end());
  if (bits.size() <= 65536)
  {
    table.resize(bits.size());
    for (std::size_t k = 0; k < bits.size(); k++)
    {
      std::memcpy(&table[k], &bits[k], sizeof(double));
    }
    auto Index = [&](double a)
    {
      return static_cast<std::uint16_t>(
          std::lower_bound(bits.begin(), bits.end(), Bits(a)) - bits.begin());
    };
    diag_v.resize(diag_a.size());
    for (std::size_t k = 0; k < diag_a.size(); k++)
    {
      diag_v[k] = Index(diag_a[k]);
    }
    diagT_v.resize(diagT_a.size());
    for (std::size_t k = 0; k < diagT_a.size(); k++)
    {
      diagT_v[k] = Index(diagT_a[k]);
    }
    std::vector<double>().swap(diag_a);
    std::vector<double>().swap(diagT_a);
  }
}

void HypreProlongation::DiagRowSums(const double *X, double *Y) const
{
  if (!table.empty())
  {
    linalg::fused::CompactCsrRowSums(diag_rows, diag_ia, diag_ja, diag_v.data(),
                                     table.data(), X, Y);
  }
  else
  {
    linalg::fused::CompactCsrRowSums(diag_rows, diag_ia, diag_ja, diag_a.data(), X, Y);
  }
}

void HypreProlongation::DiagRowSums(const double *XR, const double *XI, double *YR,
                                    double *YI) const
{
  if (!table.empty())
  {
    linalg::fused::CompactCsrRowSums(diag_rows, diag_ia, diag_ja, diag_v.data(),
                                     table.data(), XR, XI, YR, YI);
  }
  else
  {
    linalg::fused::CompactCsrRowSums(diag_rows, diag_ia, diag_ja, diag_a.data(), XR, XI, YR,
                                     YI);
  }
}

void HypreProlongation::DiagTRowSums(const std::vector<int> &rows, const double *X,
                                     double *Y) const
{
  if (!table.empty())
  {
    const std::uint16_t *v = diagT_v.data();
    const double *tab = table.data();
    FmaRowSumsV(rows, diagT_i, diagT_j, [v, tab](int k) { return tab[v[k]]; }, X, Y);
  }
  else
  {
    FmaRowSums(rows, diagT_i, diagT_j, diagT_a, X, Y);
  }
}

void HypreProlongation::DiagTRowSums(const std::vector<int> &rows, const double *XR,
                                     const double *XI, double *YR, double *YI) const
{
  if (!table.empty())
  {
    const std::uint16_t *v = diagT_v.data();
    const double *tab = table.data();
    FmaRowSumsV(
        rows, diagT_i, diagT_j, [v, tab](int k) { return tab[v[k]]; }, XR, XI, YR, YI);
  }
  else
  {
    FmaRowSums(rows, diagT_i, diagT_j, diagT_a, XR, XI, YR, YI);
  }
}

double HypreProlongation::DiagTRowSum(int t, const double *X) const
{
  double s = 0.0;
  if (!table.empty())
  {
    for (int k = diagT_i[t]; k < diagT_i[t + 1]; k++)
    {
      s = linalg::fused::Fma(table[diagT_v[k]], X[diagT_j[k]], s);
    }
  }
  else
  {
    for (int k = diagT_i[t]; k < diagT_i[t + 1]; k++)
    {
      s = linalg::fused::Fma(diagT_a[k], X[diagT_j[k]], s);
    }
  }
  return s;
}

void HypreProlongation::DiagTRowSum(int t, const double *XR, const double *XI, double &tr,
                                    double &ti) const
{
  tr = 0.0;
  ti = 0.0;
  if (!table.empty())
  {
    for (int k = diagT_i[t]; k < diagT_i[t + 1]; k++)
    {
      const double ak = table[diagT_v[k]];
      tr = linalg::fused::Fma(ak, XR[diagT_j[k]], tr);
      ti = linalg::fused::Fma(ak, XI[diagT_j[k]], ti);
    }
  }
  else
  {
    for (int k = diagT_i[t]; k < diagT_i[t + 1]; k++)
    {
      tr = linalg::fused::Fma(diagT_a[k], XR[diagT_j[k]], tr);
      ti = linalg::fused::Fma(diagT_a[k], XI[diagT_j[k]], ti);
    }
  }
}

void HypreProlongation::Mult(const Vector &x, Vector &y) const
{
  // Same operations as hypre_ParCSRMatrixMatvec, y = P x.
  hypre_ParCSRMatrix *A = P;
  hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(A), *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *X = x.HostRead();
  double *Y = y.HostWrite();

  // Start the exchange of the off-process entries of x.
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, x.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  send_buf.resize(send_end);
  recv_buf.resize(ncols_offd);
  for (int i = hypre_ParCSRCommPkgSendMapStart(comm_pkg, 0); i < send_end; i++)
  {
    send_buf[i] = X[send_map[i]];
  }
  hypre_ParCSRCommHandle *comm_handle = hypre_ParCSRCommHandleCreate_v2(
      1, comm_pkg, HYPRE_MEMORY_HOST, send_buf.data(), HYPRE_MEMORY_HOST, recv_buf.data());

  // Local part.
  CopyRuns(run_row, run_col, run_len, X, Y);
  DiagRowSums(X, Y);

  // Off-process part.
  hypre_ParCSRCommHandleDestroy(comm_handle);
  if (ncols_offd > 0)
  {
    linalg::fused::CsrRowSums(offd_rows, hypre_CSRMatrixI(offd), hypre_CSRMatrixJ(offd),
                              hypre_CSRMatrixData(offd), recv_buf.data(), Y, true);
  }
}

void HypreProlongation::MultTranspose(const Vector &x, Vector &y) const
{
  // Same operations as hypre_ParCSRMatrixMatvecT, y = Pᵀ x.
  hypre_ParCSRMatrix *A = P;
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *X = x.HostRead();
  double *Y = y.HostWrite();

  // Off-process contributions, sent to their owners.
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, y.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  recv_buf.resize(ncols_offd);
  send_buf.resize(send_end);
  if (ncols_offd > 0)
  {
    FmaRowSums(offdT_rows, offdT_i, offdT_j, offdT_a, X, recv_buf.data());
  }
  hypre_ParCSRCommHandle *comm_handle = hypre_ParCSRCommHandleCreate_v2(
      2, comm_pkg, HYPRE_MEMORY_HOST, recv_buf.data(), HYPRE_MEMORY_HOST, send_buf.data());

  // Local part.
  CopyRuns(runT_row, runT_col, runT_len, X, Y);
  DiagTRowSums(diagT_rows, X, Y);

  // Add received contributions.
  hypre_ParCSRCommHandleDestroy(comm_handle);
  for (int i = hypre_ParCSRCommPkgSendMapStart(comm_pkg, 0); i < send_end; i++)
  {
    Y[send_map[i]] += send_buf[i];
  }
}

void HypreProlongation::Exchange2(int job, const double *send_data, int send_n,
                                  double *recv_data, int recv_n) const
{
  // Forward (job 1): send x[send_map] to the send processes, receive the off-diagonal
  // columns from the receive processes. Transposed (job 2): the reverse. The messages for
  // both parts are concatenated per neighbor.
  hypre_ParCSRMatrix *A = P;
  hypre_ParCSRCommPkg *comm_pkg = hypre_ParCSRMatrixCommPkg(A);
  MPI_Comm comm = hypre_ParCSRCommPkgComm(comm_pkg);
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int num_recvs = hypre_ParCSRCommPkgNumRecvs(comm_pkg);
  const HYPRE_Int *send_procs = hypre_ParCSRCommPkgSendProcs(comm_pkg);
  const HYPRE_Int *recv_procs = hypre_ParCSRCommPkgRecvProcs(comm_pkg);
  const HYPRE_Int *send_starts = hypre_ParCSRCommPkgSendMapStarts(comm_pkg);
  const HYPRE_Int *recv_starts = hypre_ParCSRCommPkgRecvVecStarts(comm_pkg);
  const int ns = (job == 1) ? num_sends : num_recvs,
            nr = (job == 1) ? num_recvs : num_sends;
  const HYPRE_Int *sprocs = (job == 1) ? send_procs : recv_procs;
  const HYPRE_Int *rprocs = (job == 1) ? recv_procs : send_procs;
  const HYPRE_Int *sstarts = (job == 1) ? send_starts : recv_starts;
  const HYPRE_Int *rstarts = (job == 1) ? recv_starts : send_starts;
  constexpr int tag = 25431;
  requests.resize(ns + nr);
  send_buf.resize(2 * std::size_t(send_n));
  recv_buf.resize(2 * std::size_t(recv_n));
  for (int p = 0; p < nr; p++)
  {
    const int b = rstarts[p], n = rstarts[p + 1] - b;
    MPI_Irecv(recv_buf.data() + 2 * b, 2 * n, MPI_DOUBLE, rprocs[p], tag, comm,
              &requests[p]);
  }
  for (int p = 0; p < ns; p++)
  {
    const int b = sstarts[p], n = sstarts[p + 1] - b;
    double *buf = send_buf.data() + 2 * b;
    std::copy(send_data + b, send_data + b + n, buf);
    std::copy(send_data + send_n + b, send_data + send_n + b + n, buf + n);
    MPI_Isend(buf, 2 * n, MPI_DOUBLE, sprocs[p], tag, comm, &requests[nr + p]);
  }
  MPI_Waitall(ns + nr, requests.data(), MPI_STATUSES_IGNORE);
  for (int p = 0; p < nr; p++)
  {
    const int b = rstarts[p], n = rstarts[p + 1] - b;
    const double *buf = recv_buf.data() + 2 * b;
    std::copy(buf, buf + n, recv_data + b);
    std::copy(buf + n, buf + 2 * n, recv_data + recv_n + b);
  }
}

void HypreProlongation::Mult(const Vector &xr, const Vector &xi, Vector &yr,
                             Vector &yi) const
{
  hypre_ParCSRMatrix *A = P;
  hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(A), *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *XR = xr.HostRead(), *XI = xi.HostRead();
  double *YR = yr.HostWrite(), *YI = yi.HostWrite();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, xr.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  packed.resize(2 * std::size_t(send_end));
  unpacked.resize(2 * std::size_t(ncols_offd));
  for (int i = 0; i < send_end; i++)
  {
    packed[i] = XR[send_map[i]];
    packed[send_end + i] = XI[send_map[i]];
  }

  // Local part, then off-process part (rows with off-diagonal entries are not identity
  // rows, and their local part is complete before the update).
  CopyRuns(run_row, run_col, run_len, XR, XI, YR, YI);
  DiagRowSums(XR, XI, YR, YI);
  Exchange2(1, packed.data(), send_end, unpacked.data(), ncols_offd);
  if (ncols_offd > 0)
  {
    linalg::fused::CsrRowSums(offd_rows, hypre_CSRMatrixI(offd), hypre_CSRMatrixJ(offd),
                              hypre_CSRMatrixData(offd), unpacked.data(),
                              unpacked.data() + ncols_offd, YR, YI, true);
  }
}

void HypreProlongation::MultTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                      Vector &yi) const
{
  hypre_ParCSRMatrix *A = P;
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *XR = xr.HostRead(), *XI = xi.HostRead();
  double *YR = yr.HostWrite(), *YI = yi.HostWrite();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, yr.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  packed.resize(2 * std::size_t(ncols_offd));
  unpacked.resize(2 * std::size_t(send_end));
  if (ncols_offd > 0)
  {
    FmaRowSums(offdT_rows, offdT_i, offdT_j, offdT_a, XR, XI, packed.data(),
               packed.data() + ncols_offd);
  }
  CopyRuns(runT_row, runT_col, runT_len, XR, XI, YR, YI);
  DiagTRowSums(diagT_rows, XR, XI, YR, YI);
  Exchange2(2, packed.data(), ncols_offd, unpacked.data(), send_end);
  for (int i = 0; i < send_end; i++)
  {
    YR[send_map[i]] += unpacked[i];
    YI[send_map[i]] += unpacked[send_end + i];
  }
}

namespace
{

// y += a v as the contracted updates of Vector::Add and of ComplexVector::AXPY with a real
// coefficient.
inline void AddScaled(double a, double v, double &y)
{
  y = linalg::fused::Fma(a, v, y);
}

inline void AddScaled(double a, double vr, double vi, double &yr, double &yi)
{
  yr = linalg::fused::Fma(a, vr, yr);
  yi = linalg::fused::Fma(a, vi, yi);
}

}  // namespace

MaskedProlongation::MaskedProlongation(const HypreProlongation &P,
                                       const mfem::Array<int> &dofs)
  : P(P)
{
  const int *d = dofs.HostRead();
  given.assign(d, d + dofs.Size());
  list = given;
  std::sort(list.begin(), list.end());
  list.erase(std::unique(list.begin(), list.end()), list.end());
  hypre_ParCSRMatrix *A = P.P;
  hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(A);
  const int ncols = hypre_CSRMatrixNumCols(diag);
  std::vector<char> listed(ncols, 0);
  for (const int t : list)
  {
    MFEM_VERIFY(t >= 0 && t < ncols, "Invalid essential true dof for MaskedProlongation!");
    listed[t] = 1;
  }

  // Forward product: identity runs split at the listed columns, which give zero rows.
  for (std::size_t r = 0; r < P.run_len.size(); r++)
  {
    const int row = P.run_row[r], col = P.run_col[r], len = P.run_len[r];
    for (int k = 0; k < len;)
    {
      if (listed[col + k])
      {
        zero_rows.push_back(row + k++);
        continue;
      }
      int e = k;
      while (e < len && !listed[col + e])
      {
        e++;
      }
      run_row.push_back(row + k);
      run_col.push_back(col + k);
      run_len.push_back(e - k);
      k = e;
    }
  }
  ia.push_back(0);
  const bool use_table = !P.table.empty();
  for (std::size_t r = 0; r < P.diag_rows.size(); r++)
  {
    rows.push_back(P.diag_rows[r]);
    for (int k = P.diag_ia[r]; k < P.diag_ia[r + 1]; k++)
    {
      if (!listed[P.diag_ja[k]])
      {
        ja.push_back(P.diag_ja[k]);
        if (use_table)
        {
          v.push_back(P.diag_v[k]);
        }
        else
        {
          a.push_back(P.diag_a[k]);
        }
      }
    }
    ia.push_back(static_cast<int>(ja.size()));
  }
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, P.Width());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  for (int i = 0; i < send_end; i++)
  {
    if (listed[send_map[i]])
    {
      send_zero.push_back(i);
    }
  }

  // Transposed product: the unlisted true dofs with contributions from other processes are
  // completed in a buffer, the others are added to the result directly.
  std::vector<int> slot(ncols, -1);
  send_slot.resize(send_end);
  for (int i = 0; i < send_end; i++)
  {
    const int t = send_map[i];
    if (listed[t])
    {
      send_slot[i] = -1;
      continue;
    }
    if (slot[t] < 0)
    {
      slot[t] = static_cast<int>(shared.size());
      shared.push_back(t);
    }
    send_slot[i] = slot[t];
  }
  std::vector<int> ident(ncols, -1);
  for (std::size_t r = 0; r < P.runT_len.size(); r++)
  {
    for (int k = 0; k < P.runT_len[r]; k++)
    {
      ident[P.runT_row[r] + k] = P.runT_col[r] + k;
    }
  }
  shared_ldof.resize(shared.size());
  for (std::size_t s = 0; s < shared.size(); s++)
  {
    shared_ldof[s] = ident[shared[s]];
  }
  for (std::size_t r = 0; r < P.runT_len.size(); r++)
  {
    const int row = P.runT_row[r], col = P.runT_col[r], len = P.runT_len[r];
    for (int k = 0; k < len;)
    {
      if (listed[row + k] || slot[row + k] >= 0)
      {
        k++;
        continue;
      }
      int e = k;
      while (e < len && !listed[row + e] && slot[row + e] < 0)
      {
        e++;
      }
      runT_row.push_back(row + k);
      runT_col.push_back(col + k);
      runT_len.push_back(e - k);
      k = e;
    }
  }
  for (const int t : P.diagT_rows)
  {
    if (!listed[t] && slot[t] < 0)
    {
      rowsT.push_back(t);
    }
  }

  // Local dof map for products without L-vectors: the identity rows of unlisted columns are
  // their true dofs, the other rows are numbered as extra dofs in increasing order.
  // A true dof with several identity rows keeps the first one, the others are extra dofs
  // (with the copied value).
  const int nrows = hypre_CSRMatrixNumRows(diag);
  lmap.assign(nrows, -1);
  std::vector<char> owned(ncols, 0);
  for (std::size_t r = 0; r < run_len.size(); r++)
  {
    for (int k = 0; k < run_len[r]; k++)
    {
      const int l = run_row[r] + k, t = run_col[r] + k;
      if (!owned[t])
      {
        owned[t] = 1;
        lmap[l] = t;
      }
      else
      {
        dup_rows.push_back(l);
        dup_cols.push_back(t);
      }
    }
  }
  for (int i = 0; i < nrows; i++)
  {
    if (lmap[i] < 0)
    {
      lmap[i] = -1 - num_extra++;
    }
  }
  auto Extra = [this](int i)
  {
    MFEM_ASSERT(lmap[i] < 0, "Unexpected identity row for an extra dof!");
    return -1 - lmap[i];
  };
  rows_x.resize(rows.size());
  for (std::size_t r = 0; r < rows.size(); r++)
  {
    rows_x[r] = Extra(rows[r]);
  }
  zero_x.resize(zero_rows.size());
  for (std::size_t r = 0; r < zero_rows.size(); r++)
  {
    zero_x[r] = Extra(zero_rows[r]);
  }
  offd_x.resize(P.offd_rows.size());
  for (std::size_t r = 0; r < P.offd_rows.size(); r++)
  {
    offd_x[r] = Extra(P.offd_rows[r]);
  }
  dup_x.resize(dup_rows.size());
  for (std::size_t r = 0; r < dup_rows.size(); r++)
  {
    dup_x[r] = Extra(dup_rows[r]);
  }
  diagT_x.resize(P.diagT_j.size());
  for (std::size_t k = 0; k < P.diagT_j.size(); k++)
  {
    const int i = P.diagT_j[k];
    diagT_x[k] = (lmap[i] >= 0) ? -1 : -1 - lmap[i];
  }
  offdT_x.resize(P.offdT_j.size());
  for (std::size_t k = 0; k < P.offdT_j.size(); k++)
  {
    offdT_x[k] = Extra(P.offdT_j[k]);
  }
  for (const int t : P.diagT_rows)
  {
    if (!listed[t])
    {
      depT_rows.push_back(t);
    }
  }
  static std::uint64_t next_key = 0;
  key = ++next_key;
}

bool MaskedProlongation::Matches(const mfem::Array<int> &dofs) const
{
  const int *d = dofs.HostRead();
  if (static_cast<std::size_t>(dofs.Size()) == given.size() &&
      std::equal(given.begin(), given.end(), d))
  {
    return true;
  }
  std::vector<int> sorted(d, d + dofs.Size());
  std::sort(sorted.begin(), sorted.end());
  sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
  return sorted == list;
}

void MaskedProlongation::ForwardLocal(const double *X, double *Y) const
{
  CopyRuns(run_row, run_col, run_len, X, Y);
  for (const int i : zero_rows)
  {
    Y[i] = 0.0;
  }
  if (!P.table.empty())
  {
    linalg::fused::CompactCsrRowSums(rows, ia, ja, v.data(), P.table.data(), X, Y);
  }
  else
  {
    linalg::fused::CompactCsrRowSums(rows, ia, ja, a.data(), X, Y);
  }
}

void MaskedProlongation::ForwardLocal(const double *XR, const double *XI, double *YR,
                                      double *YI) const
{
  CopyRuns(run_row, run_col, run_len, XR, XI, YR, YI);
  for (const int i : zero_rows)
  {
    YR[i] = 0.0;
    YI[i] = 0.0;
  }
  if (!P.table.empty())
  {
    linalg::fused::CompactCsrRowSums(rows, ia, ja, v.data(), P.table.data(), XR, XI, YR,
                                     YI);
  }
  else
  {
    linalg::fused::CompactCsrRowSums(rows, ia, ja, a.data(), XR, XI, YR, YI);
  }
}

void MaskedProlongation::Mult(const Vector &x, Vector &y) const
{
  hypre_ParCSRMatrix *A = P.P;
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *X = x.HostRead();
  double *Y = y.HostWrite();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, x.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  auto &send_buf = P.send_buf;
  auto &recv_buf = P.recv_buf;
  send_buf.resize(send_end);
  recv_buf.resize(ncols_offd);
  for (int i = 0; i < send_end; i++)
  {
    send_buf[i] = X[send_map[i]];
  }
  for (const int i : send_zero)
  {
    send_buf[i] = 0.0;
  }
  hypre_ParCSRCommHandle *comm_handle = hypre_ParCSRCommHandleCreate_v2(
      1, comm_pkg, HYPRE_MEMORY_HOST, send_buf.data(), HYPRE_MEMORY_HOST, recv_buf.data());
  ForwardLocal(X, Y);
  hypre_ParCSRCommHandleDestroy(comm_handle);
  if (ncols_offd > 0)
  {
    linalg::fused::CsrRowSums(P.offd_rows, hypre_CSRMatrixI(offd), hypre_CSRMatrixJ(offd),
                              hypre_CSRMatrixData(offd), recv_buf.data(), Y, true);
  }
}

void MaskedProlongation::Mult(const Vector &xr, const Vector &xi, Vector &yr,
                              Vector &yi) const
{
  hypre_ParCSRMatrix *A = P.P;
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *XR = xr.HostRead(), *XI = xi.HostRead();
  double *YR = yr.HostWrite(), *YI = yi.HostWrite();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, xr.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  auto &packed = P.packed;
  auto &unpacked = P.unpacked;
  packed.resize(2 * std::size_t(send_end));
  unpacked.resize(2 * std::size_t(ncols_offd));
  for (int i = 0; i < send_end; i++)
  {
    packed[i] = XR[send_map[i]];
    packed[send_end + i] = XI[send_map[i]];
  }
  for (const int i : send_zero)
  {
    packed[i] = 0.0;
    packed[send_end + i] = 0.0;
  }
  ForwardLocal(XR, XI, YR, YI);
  P.Exchange2(1, packed.data(), send_end, unpacked.data(), ncols_offd);
  if (ncols_offd > 0)
  {
    linalg::fused::CsrRowSums(P.offd_rows, hypre_CSRMatrixI(offd), hypre_CSRMatrixJ(offd),
                              hypre_CSRMatrixData(offd), unpacked.data(),
                              unpacked.data() + ncols_offd, YR, YI, true);
  }
}

void MaskedProlongation::AddMultTranspose(const Vector &x, Vector &y, const double a,
                                          const Vector *z) const
{
  if (a == 0.0)
  {
    return;
  }
  hypre_ParCSRMatrix *A = P.P;
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *X = x.HostRead();
  double *Y = y.HostReadWrite();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, y.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  auto &send_buf = P.send_buf;
  auto &recv_buf = P.recv_buf;
  recv_buf.resize(ncols_offd);
  send_buf.resize(send_end);
  if (ncols_offd > 0)
  {
    FmaRowSums(P.offdT_rows, P.offdT_i, P.offdT_j, P.offdT_a, X, recv_buf.data());
  }
  hypre_ParCSRCommHandle *comm_handle = hypre_ParCSRCommHandleCreate_v2(
      2, comm_pkg, HYPRE_MEMORY_HOST, recv_buf.data(), HYPRE_MEMORY_HOST, send_buf.data());

  // Local part, added to the result for the true dofs without contributions from other
  // processes.
  for (std::size_t r = 0; r < runT_len.size(); r++)
  {
    const double *__restrict__ Xr = X + runT_col[r];
    double *__restrict__ Yr = Y + runT_row[r];
    const int len = runT_len[r];
    for (int k = 0; k < len; k++)
    {
      AddScaled(a, 0.0 + Xr[k], Yr[k]);
    }
  }
  for (const int t : rowsT)
  {
    AddScaled(a, P.DiagTRowSum(t, X), Y[t]);
  }
  sbuf.resize(shared.size());
  for (std::size_t s = 0; s < shared.size(); s++)
  {
    const int t = shared[s], l = shared_ldof[s];
    sbuf[s] = (l >= 0) ? 0.0 + X[l] : P.DiagTRowSum(t, X);
  }

  // Received contributions, then the shared true dofs and the essential ones.
  hypre_ParCSRCommHandleDestroy(comm_handle);
  for (int i = 0; i < send_end; i++)
  {
    if (send_slot[i] >= 0)
    {
      sbuf[send_slot[i]] += send_buf[i];
    }
  }
  for (std::size_t s = 0; s < shared.size(); s++)
  {
    AddScaled(a, sbuf[s], Y[shared[s]]);
  }
  if (z)
  {
    const double *Z = z->HostRead();
    for (const int t : list)
    {
      AddScaled(a, Z[t], Y[t]);
    }
  }
  else
  {
    // ỹ = 0 at the essential true dofs: y += a 0 (as Vector::Add).
    for (const int t : list)
    {
      AddScaled(a, 0.0, Y[t]);
    }
  }
}

void MaskedProlongation::AddMultTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                          Vector &yi, const double a, const Vector *zr,
                                          const Vector *zi) const
{
  hypre_ParCSRMatrix *A = P.P;
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *XR = xr.HostRead(), *XI = xi.HostRead();
  double *YR = yr.HostReadWrite(), *YI = yi.HostReadWrite();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, yr.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  auto &packed = P.packed;
  auto &unpacked = P.unpacked;
  packed.resize(2 * std::size_t(ncols_offd));
  unpacked.resize(2 * std::size_t(send_end));
  if (ncols_offd > 0)
  {
    FmaRowSums(P.offdT_rows, P.offdT_i, P.offdT_j, P.offdT_a, XR, XI, packed.data(),
               packed.data() + ncols_offd);
  }
  for (std::size_t r = 0; r < runT_len.size(); r++)
  {
    const double *__restrict__ XRr = XR + runT_col[r];
    const double *__restrict__ XIr = XI + runT_col[r];
    double *__restrict__ YRr = YR + runT_row[r];
    double *__restrict__ YIr = YI + runT_row[r];
    const int len = runT_len[r];
    for (int k = 0; k < len; k++)
    {
      AddScaled(a, 0.0 + XRr[k], 0.0 + XIr[k], YRr[k], YIr[k]);
    }
  }
  for (const int t : rowsT)
  {
    double tr, ti;
    P.DiagTRowSum(t, XR, XI, tr, ti);
    AddScaled(a, tr, ti, YR[t], YI[t]);
  }
  sbuf.resize(2 * shared.size());
  const std::size_t ns = shared.size();
  for (std::size_t s = 0; s < ns; s++)
  {
    const int t = shared[s], l = shared_ldof[s];
    if (l >= 0)
    {
      sbuf[s] = 0.0 + XR[l];
      sbuf[ns + s] = 0.0 + XI[l];
    }
    else
    {
      P.DiagTRowSum(t, XR, XI, sbuf[s], sbuf[ns + s]);
    }
  }
  P.Exchange2(2, packed.data(), ncols_offd, unpacked.data(), send_end);
  for (int i = 0; i < send_end; i++)
  {
    if (send_slot[i] >= 0)
    {
      sbuf[send_slot[i]] += unpacked[i];
      sbuf[ns + send_slot[i]] += unpacked[send_end + i];
    }
  }
  for (std::size_t s = 0; s < ns; s++)
  {
    AddScaled(a, sbuf[s], sbuf[ns + s], YR[shared[s]], YI[shared[s]]);
  }
  if (zr)
  {
    const double *ZR = zr->HostRead(), *ZI = zi->HostRead();
    for (const int t : list)
    {
      AddScaled(a, ZR[t], ZI[t], YR[t], YI[t]);
    }
  }
  else
  {
    // ỹ = 0 at the essential true dofs: y += a 0 (as ComplexVector::AXPY).
    for (const int t : list)
    {
      AddScaled(a, 0.0, 0.0, YR[t], YI[t]);
    }
  }
}

void MaskedProlongation::ExtraMult(const Vector &x, Vector &xx) const
{
  hypre_ParCSRMatrix *A = P.P;
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *X = x.HostRead();
  double *XX = xx.HostWrite();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, x.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  auto &send_buf = P.send_buf;
  auto &recv_buf = P.recv_buf;
  send_buf.resize(send_end);
  recv_buf.resize(ncols_offd);
  for (int i = 0; i < send_end; i++)
  {
    send_buf[i] = X[send_map[i]];
  }
  for (const int i : send_zero)
  {
    send_buf[i] = 0.0;
  }
  hypre_ParCSRCommHandle *comm_handle = hypre_ParCSRCommHandleCreate_v2(
      1, comm_pkg, HYPRE_MEMORY_HOST, send_buf.data(), HYPRE_MEMORY_HOST, recv_buf.data());
  for (const int k : zero_x)
  {
    XX[k] = 0.0;
  }
  for (std::size_t r = 0; r < dup_x.size(); r++)
  {
    XX[dup_x[r]] = 0.0 + X[dup_cols[r]];
  }
  if (!P.table.empty())
  {
    linalg::fused::CompactCsrRowSums(rows_x, ia, ja, v.data(), P.table.data(), X, XX);
  }
  else
  {
    linalg::fused::CompactCsrRowSums(rows_x, ia, ja, a.data(), X, XX);
  }
  hypre_ParCSRCommHandleDestroy(comm_handle);
  if (ncols_offd > 0)
  {
    linalg::fused::CsrRowSums(P.offd_rows, offd_x, hypre_CSRMatrixI(offd),
                              hypre_CSRMatrixJ(offd), hypre_CSRMatrixData(offd),
                              recv_buf.data(), XX, true);
  }
}

void MaskedProlongation::ExtraMult(const Vector &xr, const Vector &xi, Vector &xxr,
                                   Vector &xxi) const
{
  hypre_ParCSRMatrix *A = P.P;
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *XR = xr.HostRead(), *XI = xi.HostRead();
  double *XXR = xxr.HostWrite(), *XXI = xxi.HostWrite();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, xr.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  auto &packed = P.packed;
  auto &unpacked = P.unpacked;
  packed.resize(2 * std::size_t(send_end));
  unpacked.resize(2 * std::size_t(ncols_offd));
  for (int i = 0; i < send_end; i++)
  {
    packed[i] = XR[send_map[i]];
    packed[send_end + i] = XI[send_map[i]];
  }
  for (const int i : send_zero)
  {
    packed[i] = 0.0;
    packed[send_end + i] = 0.0;
  }
  for (const int k : zero_x)
  {
    XXR[k] = 0.0;
    XXI[k] = 0.0;
  }
  for (std::size_t r = 0; r < dup_x.size(); r++)
  {
    XXR[dup_x[r]] = 0.0 + XR[dup_cols[r]];
    XXI[dup_x[r]] = 0.0 + XI[dup_cols[r]];
  }
  if (!P.table.empty())
  {
    linalg::fused::CompactCsrRowSums(rows_x, ia, ja, v.data(), P.table.data(), XR, XI, XXR,
                                     XXI);
  }
  else
  {
    linalg::fused::CompactCsrRowSums(rows_x, ia, ja, a.data(), XR, XI, XXR, XXI);
  }
  P.Exchange2(1, packed.data(), send_end, unpacked.data(), ncols_offd);
  if (ncols_offd > 0)
  {
    linalg::fused::CsrRowSums(P.offd_rows, offd_x, hypre_CSRMatrixI(offd),
                              hypre_CSRMatrixJ(offd), hypre_CSRMatrixData(offd),
                              unpacked.data(), unpacked.data() + ncols_offd, XXR, XXI,
                              true);
  }
}

namespace
{

// Fused multiply-add row sums of a transposed block whose entry k reads its value from
// x[src[k]], or from self when src[k] < 0, with the values at[k] given by a functor.
template <typename Values>
inline double SourceRowSum(int b, int e, const std::vector<int> &src, const Values &at,
                           double self, const double *X)
{
  double t = 0.0;
  for (int k = b; k < e; k++)
  {
    t = linalg::fused::Fma(at(k), (src[k] < 0) ? self : X[src[k]], t);
  }
  return t;
}

template <typename Values>
inline void SourceRowSum(int b, int e, const std::vector<int> &src, const Values &at,
                         double selfr, double selfi, const double *XR, const double *XI,
                         double &tr, double &ti)
{
  tr = 0.0;
  ti = 0.0;
  for (int k = b; k < e; k++)
  {
    const double ak = at(k);
    const bool s = src[k] < 0;
    tr = linalg::fused::Fma(ak, s ? selfr : XR[src[k]], tr);
    ti = linalg::fused::Fma(ak, s ? selfi : XI[src[k]], ti);
  }
}

}  // namespace

void MaskedProlongation::TransposeComplete(Vector &y, const Vector &yx) const
{
  hypre_ParCSRMatrix *A = P.P;
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  double *Y = y.HostReadWrite();
  const double *YX = yx.HostRead();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, y.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  auto &send_buf = P.send_buf;
  auto &recv_buf = P.recv_buf;
  recv_buf.resize(ncols_offd);
  send_buf.resize(send_end);
  const auto offd_a = [this](int k) { return P.offdT_a[k]; };
  for (int c = 0; c < ncols_offd; c++)
  {
    recv_buf[c] = SourceRowSum(P.offdT_i[c], P.offdT_i[c + 1], offdT_x, offd_a, 0.0, YX);
  }
  hypre_ParCSRCommHandle *comm_handle = hypre_ParCSRCommHandleCreate_v2(
      2, comm_pkg, HYPRE_MEMORY_HOST, recv_buf.data(), HYPRE_MEMORY_HOST, send_buf.data());
  if (!P.table.empty())
  {
    const std::uint16_t *vv = P.diagT_v.data();
    const double *tab = P.table.data();
    const auto at = [vv, tab](int k) { return tab[vv[k]]; };
    for (const int t : depT_rows)
    {
      Y[t] = SourceRowSum(P.diagT_i[t], P.diagT_i[t + 1], diagT_x, at, Y[t], YX);
    }
  }
  else
  {
    const auto at = [this](int k) { return P.diagT_a[k]; };
    for (const int t : depT_rows)
    {
      Y[t] = SourceRowSum(P.diagT_i[t], P.diagT_i[t + 1], diagT_x, at, Y[t], YX);
    }
  }
  hypre_ParCSRCommHandleDestroy(comm_handle);
  for (int i = 0; i < send_end; i++)
  {
    Y[send_map[i]] += send_buf[i];
  }
}

void MaskedProlongation::TransposeComplete(Vector &yr, Vector &yi, const Vector &yxr,
                                           const Vector &yxi) const
{
  hypre_ParCSRMatrix *A = P.P;
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  double *YR = yr.HostReadWrite(), *YI = yi.HostReadWrite();
  const double *YXR = yxr.HostRead(), *YXI = yxi.HostRead();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, yr.Size());
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  auto &packed = P.packed;
  auto &unpacked = P.unpacked;
  packed.resize(2 * std::size_t(ncols_offd));
  unpacked.resize(2 * std::size_t(send_end));
  const auto offd_a = [this](int k) { return P.offdT_a[k]; };
  for (int c = 0; c < ncols_offd; c++)
  {
    SourceRowSum(P.offdT_i[c], P.offdT_i[c + 1], offdT_x, offd_a, 0.0, 0.0, YXR, YXI,
                 packed[c], packed[ncols_offd + c]);
  }
  if (!P.table.empty())
  {
    const std::uint16_t *vv = P.diagT_v.data();
    const double *tab = P.table.data();
    const auto at = [vv, tab](int k) { return tab[vv[k]]; };
    for (const int t : depT_rows)
    {
      SourceRowSum(P.diagT_i[t], P.diagT_i[t + 1], diagT_x, at, YR[t], YI[t], YXR, YXI,
                   YR[t], YI[t]);
    }
  }
  else
  {
    const auto at = [this](int k) { return P.diagT_a[k]; };
    for (const int t : depT_rows)
    {
      SourceRowSum(P.diagT_i[t], P.diagT_i[t + 1], diagT_x, at, YR[t], YI[t], YXR, YXI,
                   YR[t], YI[t]);
    }
  }
  P.Exchange2(2, packed.data(), ncols_offd, unpacked.data(), send_end);
  for (int i = 0; i < send_end; i++)
  {
    YR[send_map[i]] += unpacked[i];
    YI[send_map[i]] += unpacked[send_end + i];
  }
}

}  // namespace palace::hypre
