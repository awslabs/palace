// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hypre.hpp"

#include <algorithm>
#include <cmath>
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
// product on the host.
void FmaRowSums(const std::vector<int> &rows, const std::vector<int> &iat,
                const std::vector<int> &jat, const std::vector<double> &at, const double *X,
                double *Y)
{
  for (const int i : rows)
  {
    double t = 0.0;
    for (int k = iat[i]; k < iat[i + 1]; k++)
    {
      t = linalg::fused::Fma(at[k], X[jat[k]], t);
    }
    Y[i] = t;
  }
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
  for (const int i : rows)
  {
    double tr = 0.0, ti = 0.0;
    for (int k = iat[i]; k < iat[i + 1]; k++)
    {
      tr = linalg::fused::Fma(at[k], XR[jat[k]], tr);
      ti = linalg::fused::Fma(at[k], XI[jat[k]], ti);
    }
    YR[i] = tr;
    YI[i] = ti;
  }
}

// Same as FmaRowSums, for n vectors at once (each with the same operations), stored with
// strides ldx and ldy.
void FmaRowSums(const std::vector<int> &rows, const std::vector<int> &iat,
                const std::vector<int> &jat, const std::vector<double> &at, int n,
                const double *X, int ldx, double *Y, int ldy)
{
  constexpr int max_n = 32;
  if (n > max_n)
  {
    for (int c = 0; c < n; c++)
    {
      FmaRowSums(rows, iat, jat, at, X + std::size_t(c) * ldx, Y + std::size_t(c) * ldy);
    }
    return;
  }
  double t[max_n];
  for (const int i : rows)
  {
    for (int c = 0; c < n; c++)
    {
      t[c] = 0.0;
    }
    for (int k = iat[i]; k < iat[i + 1]; k++)
    {
      const double ak = at[k];
      const double *Xk = X + jat[k];
      for (int c = 0; c < n; c++)
      {
        t[c] = linalg::fused::Fma(ak, Xk[std::size_t(c) * ldx], t[c]);
      }
    }
    for (int c = 0; c < n; c++)
    {
      Y[std::size_t(c) * ldy + i] = t[c];
    }
  }
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
  linalg::fused::CsrRowSums(diag_rows, hypre_CSRMatrixI(diag), hypre_CSRMatrixJ(diag),
                            hypre_CSRMatrixData(diag), X, Y, false);

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
  FmaRowSums(diagT_rows, diagT_i, diagT_j, diagT_a, X, Y);

  // Add received contributions.
  hypre_ParCSRCommHandleDestroy(comm_handle);
  for (int i = hypre_ParCSRCommPkgSendMapStart(comm_pkg, 0); i < send_end; i++)
  {
    Y[send_map[i]] += send_buf[i];
  }
}

void HypreProlongation::ExchangeStart(int n, int job, const double *send_data, int send_n,
                                      int recv_n) const
{
  // Forward (job 1): send x[send_map] to the send processes, receive the off-diagonal
  // columns from the receive processes. Transposed (job 2): the reverse. The messages for
  // the n parts are concatenated per neighbor.
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
  send_buf.resize(n * std::size_t(send_n));
  recv_buf.resize(n * std::size_t(recv_n));
  for (int p = 0; p < nr; p++)
  {
    const int b = rstarts[p], m = rstarts[p + 1] - b;
    MPI_Irecv(recv_buf.data() + std::size_t(n) * b, n * m, MPI_DOUBLE, rprocs[p], tag, comm,
              &requests[p]);
  }
  for (int p = 0; p < ns; p++)
  {
    const int b = sstarts[p], m = sstarts[p + 1] - b;
    double *buf = send_buf.data() + std::size_t(n) * b;
    for (int k = 0; k < n; k++)
    {
      const double *data = send_data + std::size_t(k) * send_n + b;
      std::copy(data, data + m, buf + std::size_t(k) * m);
    }
    MPI_Isend(buf, n * m, MPI_DOUBLE, sprocs[p], tag, comm, &requests[nr + p]);
  }
}

void HypreProlongation::ExchangeFinish(int n, int job, double *recv_data, int recv_n) const
{
  hypre_ParCSRMatrix *A = P;
  hypre_ParCSRCommPkg *comm_pkg = hypre_ParCSRMatrixCommPkg(A);
  const int nr = (job == 1) ? hypre_ParCSRCommPkgNumRecvs(comm_pkg)
                            : hypre_ParCSRCommPkgNumSends(comm_pkg);
  const HYPRE_Int *rstarts = (job == 1) ? hypre_ParCSRCommPkgRecvVecStarts(comm_pkg)
                                        : hypre_ParCSRCommPkgSendMapStarts(comm_pkg);
  MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
  for (int p = 0; p < nr; p++)
  {
    const int b = rstarts[p], m = rstarts[p + 1] - b;
    const double *buf = recv_buf.data() + std::size_t(n) * b;
    for (int k = 0; k < n; k++)
    {
      std::copy(buf + std::size_t(k) * m, buf + std::size_t(k + 1) * m,
                recv_data + std::size_t(k) * recv_n + b);
    }
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
  linalg::fused::CsrRowSums(diag_rows, hypre_CSRMatrixI(diag), hypre_CSRMatrixJ(diag),
                            hypre_CSRMatrixData(diag), XR, XI, YR, YI, false);
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
  FmaRowSums(diagT_rows, diagT_i, diagT_j, diagT_a, XR, XI, YR, YI);
  Exchange2(2, packed.data(), ncols_offd, unpacked.data(), send_end);
  for (int i = 0; i < send_end; i++)
  {
    YR[send_map[i]] += unpacked[i];
    YI[send_map[i]] += unpacked[send_end + i];
  }
}

void HypreProlongation::MultBatch(int n, const Vector &x, Vector &y) const
{
  const int w = Width(), h = Height();
  MFEM_ASSERT(x.Size() == n * w && y.Size() == n * h,
              "Invalid vector sizes for HypreProlongation::MultBatch!");
  hypre_ParCSRMatrix *A = P;
  hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(A), *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *X = x.HostRead();
  double *Y = y.HostWrite();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, w);
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  packed.resize(n * std::size_t(send_end));
  unpacked.resize(n * std::size_t(ncols_offd));
  for (int k = 0; k < n; k++)
  {
    const double *Xk = X + std::size_t(k) * w;
    double *buf = packed.data() + std::size_t(k) * send_end;
    for (int i = 0; i < send_end; i++)
    {
      buf[i] = Xk[send_map[i]];
    }
  }
  ExchangeStart(n, 1, packed.data(), send_end, ncols_offd);

  // Local part, then off-process part (rows with off-diagonal entries are not identity
  // rows, and their local part is complete before the update).
  for (int k = 0; k < n; k++)
  {
    CopyRuns(run_row, run_col, run_len, X + std::size_t(k) * w, Y + std::size_t(k) * h);
  }
  linalg::fused::CsrRowSums(diag_rows, hypre_CSRMatrixI(diag), hypre_CSRMatrixJ(diag),
                            hypre_CSRMatrixData(diag), n, X, w, Y, h, false);
  ExchangeFinish(n, 1, unpacked.data(), ncols_offd);
  if (ncols_offd > 0)
  {
    linalg::fused::CsrRowSums(offd_rows, hypre_CSRMatrixI(offd), hypre_CSRMatrixJ(offd),
                              hypre_CSRMatrixData(offd), n, unpacked.data(), ncols_offd, Y,
                              h, true);
  }
}

void HypreProlongation::MultTransposeBatch(int n, const Vector &x, Vector &y) const
{
  const int w = Width(), h = Height();
  MFEM_ASSERT(x.Size() == n * h && y.Size() == n * w,
              "Invalid vector sizes for HypreProlongation::MultTransposeBatch!");
  hypre_ParCSRMatrix *A = P;
  hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(A);
  const int ncols_offd = hypre_CSRMatrixNumCols(offd);
  const double *X = x.HostRead();
  double *Y = y.HostWrite();
  hypre_ParCSRCommPkg *comm_pkg = GetCommPkg(A, w);
  const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
  const int send_end = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
  const HYPRE_Int *send_map = hypre_ParCSRCommPkgSendMapElmts(comm_pkg);
  packed.resize(n * std::size_t(ncols_offd));
  unpacked.resize(n * std::size_t(send_end));
  if (ncols_offd > 0)
  {
    FmaRowSums(offdT_rows, offdT_i, offdT_j, offdT_a, n, X, h, packed.data(), ncols_offd);
  }
  ExchangeStart(n, 2, packed.data(), ncols_offd, send_end);
  for (int k = 0; k < n; k++)
  {
    CopyRuns(runT_row, runT_col, runT_len, X + std::size_t(k) * h, Y + std::size_t(k) * w);
  }
  FmaRowSums(diagT_rows, diagT_i, diagT_j, diagT_a, n, X, h, Y, w);
  ExchangeFinish(n, 2, unpacked.data(), send_end);
  for (int k = 0; k < n; k++)
  {
    double *Yk = Y + std::size_t(k) * w;
    const double *buf = unpacked.data() + std::size_t(k) * send_end;
    for (int i = 0; i < send_end; i++)
    {
      Yk[send_map[i]] += buf[i];
    }
  }
}

}  // namespace palace::hypre
