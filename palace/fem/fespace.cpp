// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fespace.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <Eigen/Dense>
#include "fem/bilinearform.hpp"
#include "fem/brokenspace.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/restriction.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"

namespace palace
{

CeedBasis FiniteElementSpace::GetCeedBasis(Ceed ceed, mfem::Geometry::Type geom) const
{
  auto it = basis.find(ceed);
  MFEM_ASSERT(it != basis.end(), "Unknown Ceed context in GetCeedBasis!");
  auto &basis_map = it->second;
  auto basis_it = basis_map.find(geom);
  if (basis_it != basis_map.end())
  {
    return basis_it->second;
  }
  return basis_map.emplace(geom, BuildCeedBasis(*this, ceed, geom)).first->second;
}

CeedElemRestriction
FiniteElementSpace::GetCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                           const ceed::CeedGeomFactorData &data,
                                           ceed::CeedElementSubset subset) const
{
  const bool use_active = subset == ceed::CeedElementSubset::Active &&
                          data.active_indices.size() < data.indices.size();
  MFEM_ASSERT(!use_active || !data.active_indices.empty(),
              "Empty active subset in GetCeedElemRestriction!");
  auto &restr_map = use_active ? active_restr : restr;
  auto it = restr_map.find(ceed);
  MFEM_ASSERT(it != restr_map.end(), "Unknown Ceed context in GetCeedElemRestriction!");
  auto &geom_restr = it->second;
  auto restr_it = geom_restr.find(geom);
  if (restr_it != geom_restr.end())
  {
    return restr_it->second;
  }

  std::vector<int> active_indices;
  if (use_active)
  {
    active_indices.reserve(data.active_indices.size());
    for (auto i : data.active_indices)
    {
      MFEM_ASSERT(i >= 0 && static_cast<std::size_t>(i) < data.indices.size(),
                  "Invalid element subset index!");
      active_indices.push_back(data.indices[i]);
    }
  }
  const auto &indices = use_active ? active_indices : data.indices;
  if (broken)
  {
    return geom_restr
        .emplace(geom, BuildBrokenCeedElemRestriction(ceed, geom, indices, false, false))
        .first->second;
  }
  return geom_restr.emplace(geom, BuildCeedElemRestriction(*this, ceed, geom, indices))
      .first->second;
}

CeedElemRestriction
FiniteElementSpace::GetInterpCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                                 const ceed::CeedGeomFactorData &data) const
{
  const mfem::FiniteElement &fe = *GetFEColl().FiniteElementForGeometry(geom);
  if (!HasUniqueInterpRestriction(fe))
  {
    return GetCeedElemRestriction(ceed, geom, data);
  }
  auto it = interp_restr.find(ceed);
  MFEM_ASSERT(it != interp_restr.end(),
              "Unknown Ceed context in GetInterpCeedElemRestriction!");
  auto &restr_map = it->second;
  auto restr_it = restr_map.find(geom);
  if (restr_it != restr_map.end())
  {
    return restr_it->second;
  }
  // A broken space reads the copies of the side of each element (the element DOF
  // overrides), also for the discrete interpolators between broken spaces.
  return restr_map
      .emplace(geom,
               broken
                   ? BuildBrokenCeedElemRestriction(ceed, geom, data.indices, true, false)
                   : BuildCeedElemRestriction(*this, ceed, geom, data.indices, true, false))
      .first->second;
}

CeedElemRestriction FiniteElementSpace::GetInterpRangeCeedElemRestriction(
    Ceed ceed, mfem::Geometry::Type geom, const ceed::CeedGeomFactorData &data) const
{
  const mfem::FiniteElement &fe = *GetFEColl().FiniteElementForGeometry(geom);
  if (!HasUniqueInterpRangeRestriction(fe))
  {
    return GetInterpCeedElemRestriction(ceed, geom, data);
  }
  auto it = interp_range_restr.find(ceed);
  MFEM_ASSERT(it != interp_range_restr.end(),
              "Unknown Ceed context in GetInterpRangeCeedElemRestriction!");
  auto &restr_map = it->second;
  auto restr_it = restr_map.find(geom);
  if (restr_it != restr_map.end())
  {
    return restr_it->second;
  }
  return restr_map
      .emplace(geom,
               broken
                   ? BuildBrokenCeedElemRestriction(ceed, geom, data.indices, true, true)
                   : BuildCeedElemRestriction(*this, ceed, geom, data.indices, true, true))
      .first->second;
}

CeedElemRestriction FiniteElementSpace::BuildBrokenCeedElemRestriction(
    Ceed ceed, mfem::Geometry::Type geom, const std::vector<int> &indices, bool is_interp,
    bool is_interp_range) const
{
  MFEM_VERIFY(broken && mfem::Geometry::Dimension[geom] == GetParMesh().Dimension(),
              "Only domain element restrictions are available for a broken space!");
  CeedElemRestriction val;
  ceed::ElementDofRemap remap;
  remap.offsets = broken->override_offsets.data();
  remap.local = broken->override_local.data();
  remap.ldof = broken->override_ldof.data();
  remap.l_size = broken->vsize;
  ceed::InitRestriction(Get(), indices, false, is_interp, is_interp_range, ceed, &val,
                        &remap);
  return val;
}

bool FiniteElementSpace::GetRestrictionWeights(Vector &w) const
{
  MFEM_VERIFY(broken, "Restriction weights are only used for broken spaces!");
  // An interpolated field (an interpolant of a field of a conforming space) has at every
  // L-DOF i read by a local element the value (P t)_i of its true DOFs t. A true DOF is
  // thus determined by any row of P with a single nonzero entry (an unconstrained DOF or a
  // copy, or a DOF of a hanging entity constrained by a single master DOF), and two true
  // DOFs by two rows of P whose nonzeros are on them and form a nonsingular block T (the
  // DOFs of a triangular face shared with another process, whose orientation P
  // transforms, or the constrained DOFs of a hanging face). Only the L-DOFs read by a local
  // element hold values: the copy block of a process can hold copies which only the
  // elements of other processes read. With the blocks B of disjoint rows chosen on each
  // process, q = T_B⁻ᵀ 1 on their rows, c = Pᵀ q counts the determinations of each true DOF
  // over all processes, and R = Pᵀ W with W = T_B⁻ᵀ diag(1 / c) T_B⁻¹ on the blocks (the
  // diagonal weights w for the single rows, the pairs for the others, with w = 1 there).
  const auto &P = *broken->P.P;
  hypre_ParCSRMatrix *hP = const_cast<mfem::HypreParMatrix &>(P);
  const hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(hP),
                        *offd = hypre_ParCSRMatrixOffd(hP);
  const HYPRE_Int *di = hypre_CSRMatrixI(diag), *dj = hypre_CSRMatrixJ(diag),
                  *oi = hypre_CSRMatrixI(offd), *oj = hypre_CSRMatrixJ(offd);
  const HYPRE_Real *dv = hypre_CSRMatrixData(diag), *ov = hypre_CSRMatrixData(offd);
  const HYPRE_BigInt *cmap = hypre_ParCSRMatrixColMapOffd(hP);
  const HYPRE_BigInt col0 = hypre_ParCSRMatrixFirstColDiag(hP);
  const int n = broken->vsize;
  MFEM_VERIFY(hypre_CSRMatrixNumRows(diag) == n, "Unexpected broken prolongation size!");
  std::vector<char> read(n, 0);
  {
    mfem::Array<int> dofs;
    mfem::DofTransformation dof_trans;
    for (int e = 0; e < GetMesh().GetNE(); e++)
    {
      GetElementDofs(e, dofs, dof_trans);
      for (int d : dofs)
      {
        read[(d >= 0) ? d : -1 - d] = 1;
      }
    }
  }
  // The nonzero entries of a row (global columns).
  auto Nonzeros = [&](int i, std::array<HYPRE_BigInt, 3> &col, std::array<double, 3> &val)
  {
    int m = 0;
    for (HYPRE_Int k = di[i]; k < di[i + 1]; k++)
    {
      if (dv[k] != 0.0)
      {
        if (m < 3)
        {
          col[m] = col0 + dj[k];
          val[m] = dv[k];
        }
        m++;
      }
    }
    for (HYPRE_Int k = (oi ? oi[i] : 0); k < (oi ? oi[i + 1] : 0); k++)
    {
      if (ov[k] != 0.0)
      {
        if (m < 3)
        {
          col[m] = cmap[oj[k]];
          val[m] = ov[k];
        }
        m++;
      }
    }
    return m;
  };
  // Rows with one nonzero (flagged) and with two (grouped by their sorted columns; few:
  // shared faces and hanging entities).
  std::vector<char> single(n, 0);
  std::map<std::array<HYPRE_BigInt, 2>, std::vector<int>> two;
  for (int i = 0; i < n; i++)
  {
    if (!read[i])
    {
      continue;
    }
    std::array<HYPRE_BigInt, 3> col;
    std::array<double, 3> val;
    const int m = Nonzeros(i, col, val);
    if (m == 1)
    {
      single[i] = 1;
    }
    else if (m == 2)
    {
      two[{std::min(col[0], col[1]), std::max(col[0], col[1])}].push_back(i);
    }
  }
  auto SingleValue = [&](int i)
  {
    std::array<HYPRE_BigInt, 3> col;
    std::array<double, 3> val;
    Nonzeros(i, col, val);
    return val[0];
  };
  // The single rows on the columns of the groups with an odd number of rows (candidates to
  // complete their last row to a pair).
  std::map<HYPRE_BigInt, std::vector<int>> one;
  for (const auto &[cols, rows] : two)
  {
    if (rows.size() % 2)
    {
      one[cols[0]];
      one[cols[1]];
    }
  }
  if (!one.empty())
  {
    for (int i = 0; i < n; i++)
    {
      if (single[i])
      {
        std::array<HYPRE_BigInt, 3> col;
        std::array<double, 3> val;
        Nonzeros(i, col, val);
        auto it = one.find(col[0]);
        if (it != one.end())
        {
          it->second.push_back(i);
        }
      }
    }
  }
  // The 2 x 2 blocks: two rows of a group of two-entry rows, or a leftover row of a group
  // with a single-entry row on one of its columns. T has the columns in increasing global
  // order.
  auto &pairs = broken->restriction_pairs;
  auto &blocks = broken->restriction_blocks;
  pairs.clear();
  blocks.clear();
  std::vector<char> used(n, 0);
  std::vector<std::array<double, 4>> inverse;  // T⁻¹ (row major), per pair.
  Vector q(n), c(broken->tsize);
  q = 0.0;
  auto Row = [&](int i, const std::array<HYPRE_BigInt, 2> &cols, double *t)
  {
    std::array<HYPRE_BigInt, 3> col;
    std::array<double, 3> val;
    const int m = Nonzeros(i, col, val);
    t[0] = t[1] = 0.0;
    for (int k = 0; k < m; k++)
    {
      t[col[k] == cols[0] ? 0 : 1] = val[k];
    }
  };
  auto AddPair = [&](int i0, int i1, const std::array<HYPRE_BigInt, 2> &cols)
  {
    double T[4];
    Row(i0, cols, T);
    Row(i1, cols, T + 2);
    const double det = T[0] * T[3] - T[1] * T[2],
                 scale = std::max(
                     {std::abs(T[0]), std::abs(T[1]), std::abs(T[2]), std::abs(T[3])});
    if (!(std::abs(det) > 1.0e-8 * scale * scale))
    {
      return false;
    }
    const std::array<double, 4> Ti = {T[3] / det, -T[1] / det, -T[2] / det, T[0] / det};
    pairs.push_back(i0);
    pairs.push_back(i1);
    inverse.push_back(Ti);
    q(i0) = Ti[0] + Ti[2];  // T⁻ᵀ (1, 1).
    q(i1) = Ti[1] + Ti[3];
    used[i0] = used[i1] = 1;
    return true;
  };
  for (const auto &[cols, rows] : two)
  {
    int left = -1;
    for (int i : rows)
    {
      if (left < 0)
      {
        left = i;
      }
      else if (AddPair(left, i, cols))
      {
        left = -1;
      }
    }
    if (left >= 0)
    {
      for (const auto col : cols)
      {
        auto it = one.find(col);
        if (it == one.end())
        {
          continue;
        }
        bool done = false;
        for (int i : it->second)
        {
          if (!used[i] && AddPair(i, left, cols))
          {
            done = true;
            break;
          }
        }
        if (done)
        {
          break;
        }
      }
    }
  }
  for (int i = 0; i < n; i++)
  {
    if (single[i] && !used[i])
    {
      q(i) = 1.0 / SingleValue(i);
    }
  }
  P.MultTranspose(q, c);
  // Every true DOF must be determined somewhere.
  int missing = 0;
  for (int t = 0; t < c.Size(); t++)
  {
    c(t) = std::round(c(t));
    missing += (c(t) < 0.5);
  }
  Mpi::GlobalSum(1, &missing, GetComm());
  // The true DOFs undetermined by single rows and pairs (copies of master DOFs read on a
  // side only through constraint rows with several masters, on nonconforming meshes) are
  // recovered from the constraint rows of their hanging entities (see
  // BuildHangingRestriction).
  int unresolved = missing;
  if (missing > 0)
  {
    if (!broken->hanging)
    {
      BuildHangingRestriction(read, c);
    }
    unresolved = broken->hanging->unresolved;
  }
  else
  {
    broken->hanging.reset();
  }
  broken->restriction_undetermined = unresolved;
  broken->restriction_hanging = missing - unresolved;
  if (unresolved > 0 && !broken->allow_undetermined)
  {
    pairs.clear();
    blocks.clear();
    broken->hanging.reset();
    w.SetSize(0);
    return false;
  }
  // The counts at the rows: P c = val c_t at the single rows and T (c_t0, c_t1) at the
  // pairs.
  w.SetSize(n);
  P.Mult(c, w);
  for (std::size_t k = 0; k < inverse.size(); k++)
  {
    const int i0 = pairs[2 * k], i1 = pairs[2 * k + 1];
    const auto &Ti = inverse[k];
    const double c0 = Ti[0] * w(i0) + Ti[1] * w(i1), c1 = Ti[2] * w(i0) + Ti[3] * w(i1);
    // W = T⁻ᵀ diag(1 / c0, 1 / c1) T⁻¹.
    const double d0 = 1.0 / std::round(c0), d1 = 1.0 / std::round(c1);
    blocks.push_back(Ti[0] * d0 * Ti[0] + Ti[2] * d1 * Ti[2]);
    blocks.push_back(Ti[0] * d0 * Ti[1] + Ti[2] * d1 * Ti[3]);
    blocks.push_back(Ti[1] * d0 * Ti[0] + Ti[3] * d1 * Ti[2]);
    blocks.push_back(Ti[1] * d0 * Ti[1] + Ti[3] * d1 * Ti[3]);
  }
  for (int i = 0; i < n; i++)
  {
    // W = 1 / (val² c_t) at a single row, (P c)_i = val c_t.
    if (single[i] && !used[i])
    {
      const double val = SingleValue(i);
      w(i) = 1.0 / (val * val * std::round(w(i) / val));
    }
    else
    {
      w(i) = 0.0;
    }
  }
  for (int i : pairs)
  {
    w(i) = 1.0;  // W is applied by ApplyRestrictionPairs.
  }
  if (broken->hanging)
  {
    for (int i : broken->hanging->rows)
    {
      w(i) = 1.0;  // The constraint rows are read by BrokenRestrictionMult.
    }
  }
  return true;
}

void FiniteElementSpace::ApplyRestrictionPairs(Vector &lx) const
{
  MFEM_ASSERT(broken, "Restriction pairs are only used for broken spaces!");
  const auto &pairs = broken->restriction_pairs;
  const auto &blocks = broken->restriction_blocks;
  if (pairs.empty())
  {
    return;
  }
  auto *x = lx.HostReadWrite();
  for (std::size_t k = 0; k < pairs.size() / 2; k++)
  {
    const int i0 = pairs[2 * k], i1 = pairs[2 * k + 1];
    const double a = x[i0], b = x[i1];
    x[i0] = blocks[4 * k + 0] * a + blocks[4 * k + 1] * b;
    x[i1] = blocks[4 * k + 2] * a + blocks[4 * k + 3] * b;
  }
}

namespace
{

// The partition array of a distributed set with n local entries, as for MFEM's
// HypreParMatrix (2 entries with HYPRE's assumed partition, the offsets of all processes
// otherwise), and the global size.
std::vector<HYPRE_BigInt> CompactOffsets(MPI_Comm comm, int n, HYPRE_BigInt &global)
{
  const int rank = Mpi::Rank(comm), nprocs = Mpi::Size(comm);
  std::vector<HYPRE_BigInt> all(nprocs + 1, 0);
  HYPRE_BigInt local = n;
  std::vector<HYPRE_BigInt> counts(nprocs);
  MPI_Allgather(&local, 1, HYPRE_MPI_BIG_INT, counts.data(), 1, HYPRE_MPI_BIG_INT, comm);
  for (int p = 0; p < nprocs; p++)
  {
    all[p + 1] = all[p] + counts[p];
  }
  global = all[nprocs];
  if (HYPRE_AssumedPartitionCheck())
  {
    return {all[rank], all[rank + 1]};
  }
  return all;
}

// A HypreParMatrix with the given local rows (CSR with global columns) and partitions.
std::unique_ptr<mfem::HypreParMatrix>
MakeParMatrix(MPI_Comm comm, const std::vector<HYPRE_BigInt> &row_part,
              HYPRE_BigInt global_rows, const HYPRE_BigInt *col_part,
              HYPRE_BigInt global_cols, std::vector<int> &I, std::vector<HYPRE_BigInt> &J,
              std::vector<double> &V)
{
  const int nrows = static_cast<int>(I.size()) - 1;
  if (J.empty())
  {
    // Valid (non-null) arrays for the constructor.
    J.push_back(0);
    V.push_back(0.0);
  }
  return std::make_unique<mfem::HypreParMatrix>(
      comm, nrows, global_rows, global_cols, I.data(), J.data(), V.data(),
      const_cast<HYPRE_BigInt *>(row_part.data()), const_cast<HYPRE_BigInt *>(col_part));
}

}  // namespace

void FiniteElementSpace::BuildHangingRestriction(const std::vector<char> &read,
                                                 const Vector &c) const
{
  // An interpolated field y = P t has at every constraint row i read by a local element
  // y_i = sum_j P_ij t_j. The true DOFs which no single row or pair of rows determines
  // (u: c < 1/2) are the copies of master DOFs which the elements of their side read only
  // through such rows: with the others (d) known, y_c - P_cd t_d = P_cu t_u on the rows
  // with an entry on u, and P_cu has full column rank on each cluster of u (connected by
  // the rows), since the constraints embed the parent space into the children's space.
  // The rows of all processes are gathered on the first process, which solves each
  // cluster's least-squares problem (SVD; columns not in the row space of their cluster's
  // rows remain unresolved) and returns the coefficients of the pseudo-inverse on each
  // process's rows (setup only: the clusters are the hanging entities of the interior
  // boundaries, a few thousand DOFs).
  MFEM_VERIFY(broken, "Hanging restriction requires a broken space!");
  auto &H = *(broken->hanging = std::make_unique<BrokenData::HangingRestriction>());
  const MPI_Comm comm = GetComm();
  const int rank = Mpi::Rank(comm), nprocs = Mpi::Size(comm);
  const auto &P = *broken->P.P;
  hypre_ParCSRMatrix *hP = const_cast<mfem::HypreParMatrix &>(P);
  const hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(hP),
                        *offd = hypre_ParCSRMatrixOffd(hP);
  const HYPRE_Int *di = hypre_CSRMatrixI(diag), *dj = hypre_CSRMatrixJ(diag),
                  *oi = hypre_CSRMatrixI(offd), *oj = hypre_CSRMatrixJ(offd);
  const HYPRE_Real *dv = hypre_CSRMatrixData(diag), *ov = hypre_CSRMatrixData(offd);
  const HYPRE_BigInt *cmap = hypre_ParCSRMatrixColMapOffd(hP);
  const int num_cols_offd = hypre_CSRMatrixNumCols(offd);
  const HYPRE_BigInt col0 = hypre_ParCSRMatrixFirstColDiag(hP);
  const int n = broken->vsize, tsize = broken->tsize;

  // Undetermined flags of the local true DOFs and of the off-diagonal columns of P.
  std::vector<char> undet(tsize), ghost_undet(num_cols_offd, 0);
  for (int t = 0; t < tsize; t++)
  {
    undet[t] = (c(t) < 0.5);
  }
  {
    if (!hypre_ParCSRMatrixCommPkg(hP))
    {
      hypre_MatvecCommPkgCreate(hP);
    }
    hypre_ParCSRCommPkg *comm_pkg = hypre_ParCSRMatrixCommPkg(hP);
    const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
    const int send_size = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
    std::vector<HYPRE_Complex> send(send_size), recv(num_cols_offd);
    for (int k = 0; k < send_size; k++)
    {
      send[k] = undet[hypre_ParCSRCommPkgSendMapElmt(comm_pkg, k)];
    }
    auto *handle = hypre_ParCSRCommHandleCreate(1, comm_pkg, send.data(), recv.data());
    hypre_ParCSRCommHandleDestroy(handle);
    for (int k = 0; k < num_cols_offd; k++)
    {
      ghost_undet[k] = (recv[k] > 0.5);
    }
  }

  // The rows read by local elements with an entry on an undetermined column, with these
  // entries (global columns).
  std::vector<int> rows, row_nnz;
  std::vector<std::int64_t> cols;
  std::vector<double> vals;
  for (int i = 0; i < n; i++)
  {
    if (!read[i])
    {
      continue;
    }
    int m = 0;
    for (HYPRE_Int k = di[i]; k < di[i + 1]; k++)
    {
      if (dv[k] != 0.0 && undet[dj[k]])
      {
        cols.push_back(col0 + dj[k]);
        vals.push_back(dv[k]);
        m++;
      }
    }
    for (HYPRE_Int k = (oi ? oi[i] : 0); k < (oi ? oi[i + 1] : 0); k++)
    {
      if (ov[k] != 0.0 && ghost_undet[oj[k]])
      {
        cols.push_back(cmap[oj[k]]);
        vals.push_back(ov[k]);
        m++;
      }
    }
    if (m > 0)
    {
      rows.push_back(i);
      row_nnz.push_back(m);
    }
  }

  // Gather the rows on the first process.
  auto Gather =
      [&](const auto &local, auto &all, MPI_Datatype type, std::vector<int> &counts)
  {
    int count = static_cast<int>(local.size());
    counts.assign(nprocs, 0);
    MPI_Gather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
    std::vector<int> displs(nprocs, 0);
    for (int p = 1; p < nprocs; p++)
    {
      displs[p] = displs[p - 1] + counts[p - 1];
    }
    if (rank == 0)
    {
      all.resize(displs[nprocs - 1] + counts[nprocs - 1]);
    }
    MPI_Gatherv(local.data(), count, type, all.data(), counts.data(), displs.data(), type,
                0, comm);
  };
  std::vector<int> all_nnz, row_counts, entry_counts;
  std::vector<std::int64_t> all_cols;
  std::vector<double> all_vals;
  Gather(row_nnz, all_nnz, MPI_INT, row_counts);
  Gather(cols, all_cols, MPI_INT64_T, entry_counts);
  Gather(vals, all_vals, MPI_DOUBLE, entry_counts);

  // On the first process: the clusters (rows connected by their undetermined columns),
  // their pseudo-inverses, and the coefficients per process (compact local row, global
  // column, coefficient).
  std::vector<std::vector<std::int64_t>> send_idx(nprocs);
  std::vector<std::vector<double>> send_coef(nprocs);
  int counts[2] = {0, 0};  // Recovered, unresolved.
  if (rank == 0)
  {
    const int nr = static_cast<int>(all_nnz.size());
    std::vector<int> row_rank(nr), row_local(nr), row_start(nr + 1, 0);
    for (int p = 0, r = 0; p < nprocs; p++)
    {
      for (int k = 0; k < row_counts[p]; k++, r++)
      {
        row_rank[r] = p;
        row_local[r] = k;
      }
    }
    for (int r = 0; r < nr; r++)
    {
      row_start[r + 1] = row_start[r] + all_nnz[r];
    }
    std::vector<std::int64_t> ucols(all_cols);
    std::sort(ucols.begin(), ucols.end());
    ucols.erase(std::unique(ucols.begin(), ucols.end()), ucols.end());
    auto Index = [&](std::int64_t g)
    {
      return static_cast<int>(std::lower_bound(ucols.begin(), ucols.end(), g) -
                              ucols.begin());
    };
    std::vector<int> parent(ucols.size());
    for (std::size_t k = 0; k < parent.size(); k++)
    {
      parent[k] = static_cast<int>(k);
    }
    std::function<int(int)> Find = [&](int k)
    { return (parent[k] == k) ? k : (parent[k] = Find(parent[k])); };
    for (int r = 0; r < nr; r++)
    {
      const int a = Find(Index(all_cols[row_start[r]]));
      for (int k = row_start[r] + 1; k < row_start[r + 1]; k++)
      {
        const int b = Find(Index(all_cols[k]));
        if (a != b)
        {
          parent[b] = a;
        }
      }
    }
    std::map<int, std::vector<int>> cluster_rows, cluster_cols;
    for (int r = 0; r < nr; r++)
    {
      cluster_rows[Find(Index(all_cols[row_start[r]]))].push_back(r);
    }
    for (std::size_t k = 0; k < ucols.size(); k++)
    {
      cluster_cols[Find(static_cast<int>(k))].push_back(static_cast<int>(k));
    }
    for (const auto &[root, crows] : cluster_rows)
    {
      const auto &ccols = cluster_cols[root];
      const int m = static_cast<int>(crows.size()), k = static_cast<int>(ccols.size());
      std::unordered_map<int, int> local_col;
      for (int j = 0; j < k; j++)
      {
        local_col[ccols[j]] = j;
      }
      Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, k);
      for (int r = 0; r < m; r++)
      {
        for (int e = row_start[crows[r]]; e < row_start[crows[r] + 1]; e++)
        {
          A(r, local_col[Index(all_cols[e])]) += all_vals[e];
        }
      }
      Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
      const auto &s = svd.singularValues();
      int rk = 0;
      while (rk < s.size() && s(rk) > 1.0e-10 * s(0))
      {
        rk++;
      }
      // X = V_r Σ_r⁻¹ U_rᵀ (k x m); column j is resolved when e_j is in the row space of A
      // ((V_r V_rᵀ)_jj = 1), then x_j = (X y)_j exactly for y = A x.
      const Eigen::MatrixXd Vr = svd.matrixV().leftCols(rk);
      const Eigen::MatrixXd X = Vr * s.head(rk).cwiseInverse().asDiagonal() *
                                svd.matrixU().leftCols(rk).transpose();
      for (int j = 0; j < k; j++)
      {
        if (Vr.row(j).squaredNorm() < 1.0 - 1.0e-8)
        {
          counts[1]++;
          continue;
        }
        counts[0]++;
        for (int r = 0; r < m; r++)
        {
          const double x = X(j, r);
          if (x == 0.0)
          {
            continue;
          }
          const int p = row_rank[crows[r]];
          send_idx[p].push_back(row_local[crows[r]]);
          send_idx[p].push_back(ucols[ccols[j]]);
          send_coef[p].push_back(x);
        }
      }
    }
  }
  Mpi::Broadcast(2, counts, 0, comm);
  H.recovered = counts[0];
  H.unresolved = counts[1];

  // Return the coefficients of each process's rows.
  auto Scatter = [&](const auto &send, auto &recv, MPI_Datatype type)
  {
    std::vector<int> scounts(nprocs, 0), displs(nprocs, 0);
    std::vector<typename std::decay_t<decltype(recv)>::value_type> buf;
    if (rank == 0)
    {
      for (int p = 0; p < nprocs; p++)
      {
        scounts[p] = static_cast<int>(send[p].size());
        displs[p] = (p > 0) ? displs[p - 1] + scounts[p - 1] : 0;
        buf.insert(buf.end(), send[p].begin(), send[p].end());
      }
    }
    int count = 0;
    MPI_Scatter(scounts.data(), 1, MPI_INT, &count, 1, MPI_INT, 0, comm);
    recv.resize(count);
    MPI_Scatterv(buf.data(), scounts.data(), displs.data(), type, recv.data(), count, type,
                 0, comm);
  };
  std::vector<std::int64_t> idx;
  std::vector<double> coef;
  Scatter(send_idx, idx, MPI_INT64_T);
  Scatter(send_coef, coef, MPI_DOUBLE);

  // The compact rows: the rows with coefficients, in increasing order.
  std::vector<std::vector<std::pair<HYPRE_BigInt, double>>> row_coef(rows.size());
  for (std::size_t e = 0; e < coef.size(); e++)
  {
    row_coef[idx[2 * e]].emplace_back(idx[2 * e + 1], coef[e]);
  }
  std::vector<int> CtI(1, 0), PcI(1, 0);
  std::vector<HYPRE_BigInt> CtJ, PcJ;
  std::vector<double> CtV, PcV;
  for (std::size_t r = 0; r < rows.size(); r++)
  {
    if (row_coef[r].empty())
    {
      continue;
    }
    const int i = rows[r];
    H.rows.push_back(i);
    for (const auto &[j, x] : row_coef[r])
    {
      CtJ.push_back(j);
      CtV.push_back(x);
    }
    CtI.push_back(static_cast<int>(CtJ.size()));
    for (HYPRE_Int k = di[i]; k < di[i + 1]; k++)
    {
      if (dv[k] != 0.0 && !undet[dj[k]])
      {
        PcJ.push_back(col0 + dj[k]);
        PcV.push_back(dv[k]);
      }
    }
    for (HYPRE_Int k = (oi ? oi[i] : 0); k < (oi ? oi[i + 1] : 0); k++)
    {
      if (ov[k] != 0.0 && !ghost_undet[oj[k]])
      {
        PcJ.push_back(cmap[oj[k]]);
        PcV.push_back(ov[k]);
      }
    }
    PcI.push_back(static_cast<int>(PcJ.size()));
  }
  HYPRE_BigInt global_rows;
  const auto row_part = CompactOffsets(comm, static_cast<int>(H.rows.size()), global_rows);
  const HYPRE_BigInt *col_part = broken->P.col_offsets.GetData();
  H.Ct = MakeParMatrix(comm, row_part, global_rows, col_part, broken->global_tsize, CtI,
                       CtJ, CtV);
  H.Pc = MakeParMatrix(comm, row_part, global_rows, col_part, broken->global_tsize, PcI,
                       PcJ, PcV);

  // M = P Pcᵀ, on its nonzero rows (the L-DOFs affected by the determined true DOFs of the
  // constraint rows).
  {
    std::unique_ptr<mfem::HypreParMatrix> PcT(H.Pc->Transpose());
    std::unique_ptr<mfem::HypreParMatrix> PPcT(mfem::ParMult(&P, PcT.get()));
    PcT.reset();
    hypre_ParCSRMatrix *hM = *PPcT;
    const hypre_CSRMatrix *Md = hypre_ParCSRMatrixDiag(hM),
                          *Mo = hypre_ParCSRMatrixOffd(hM);
    const HYPRE_Int *Mdi = hypre_CSRMatrixI(Md), *Mdj = hypre_CSRMatrixJ(Md),
                    *Moi = hypre_CSRMatrixI(Mo), *Moj = hypre_CSRMatrixJ(Mo);
    const HYPRE_Real *Mdv = hypre_CSRMatrixData(Md), *Mov = hypre_CSRMatrixData(Mo);
    const HYPRE_BigInt *Mcmap = hypre_ParCSRMatrixColMapOffd(hM);
    const HYPRE_BigInt Mcol0 = hypre_ParCSRMatrixFirstColDiag(hM);
    std::vector<int> MI(1, 0);
    std::vector<HYPRE_BigInt> MJ;
    std::vector<double> MV;
    for (int i = 0; i < hypre_CSRMatrixNumRows(Md); i++)
    {
      const std::size_t nnz = MJ.size();
      for (HYPRE_Int k = Mdi[i]; k < Mdi[i + 1]; k++)
      {
        if (Mdv[k] != 0.0)
        {
          MJ.push_back(Mcol0 + Mdj[k]);
          MV.push_back(Mdv[k]);
        }
      }
      for (HYPRE_Int k = (Moi ? Moi[i] : 0); k < (Moi ? Moi[i + 1] : 0); k++)
      {
        if (Mov[k] != 0.0)
        {
          MJ.push_back(Mcmap[Moj[k]]);
          MV.push_back(Mov[k]);
        }
      }
      if (MJ.size() > nnz)
      {
        H.affected.push_back(i);
        MI.push_back(static_cast<int>(MJ.size()));
      }
    }
    PPcT.reset();
    HYPRE_BigInt global_affected;
    const auto affected_part =
        CompactOffsets(comm, static_cast<int>(H.affected.size()), global_affected);
    H.M = MakeParMatrix(comm, affected_part, global_affected, row_part.data(), global_rows,
                        MI, MJ, MV);
  }
  H.yc.SetSize(static_cast<int>(H.rows.size()));
  H.mc.SetSize(static_cast<int>(H.affected.size()));
  H.yc.UseDevice(false);
  H.mc.UseDevice(false);
}

void FiniteElementSpace::BrokenRestrictionMult(Vector &ly, Vector &ty) const
{
  MFEM_ASSERT(broken, "Broken restriction requires a broken space!");
  // ty = (I - Ctᵀ Pc) Pᵀ Z ly + Ctᵀ S ly, with Z the 2 x 2 blocks of W and the constraint
  // rows S zeroed (the diagonal of W is applied to ly by the interpolator).
  ApplyRestrictionPairs(ly);
  const auto *H = broken->hanging.get();
  if (H && !H->rows.empty())
  {
    auto *x = ly.HostReadWrite();
    auto *yc = H->yc.HostWrite();
    for (std::size_t k = 0; k < H->rows.size(); k++)
    {
      yc[k] = x[H->rows[k]];
      x[H->rows[k]] = 0.0;
    }
  }
  broken->P.P->MultTranspose(ly, ty);
  if (H)
  {
    H->Pc->Mult(-1.0, ty, 1.0, H->yc);
    H->Ct->MultTranspose(1.0, H->yc, 1.0, ty);
  }
}

void FiniteElementSpace::BrokenRestrictionMultTranspose(const Vector &ty, Vector &ly) const
{
  MFEM_ASSERT(broken, "Broken restriction requires a broken space!");
  // ly = Zᵀ (P ty - (P Pcᵀ) Ct ty) + Sᵀ Ct ty.
  broken->P.P->Mult(ty, ly);
  const auto *H = broken->hanging.get();
  if (H)
  {
    H->Ct->Mult(ty, H->yc);
    H->M->Mult(H->yc, H->mc);
    auto *x = ly.HostReadWrite();
    const auto *mc = H->mc.HostRead();
    for (std::size_t k = 0; k < H->affected.size(); k++)
    {
      x[H->affected[k]] -= mc[k];
    }
  }
  ApplyRestrictionPairs(ly);  // W is symmetric.
  if (H)
  {
    auto *x = ly.HostReadWrite();
    const auto *yc = H->yc.HostRead();
    for (std::size_t k = 0; k < H->rows.size(); k++)
    {
      x[H->rows[k]] = yc[k];
    }
  }
}

void FiniteElementSpace::ResetCeedObjects()
{
  for (auto &[ceed, basis_map] : basis)
  {
    for (auto &[key, val] : basis_map)
    {
      PalaceCeedCall(ceed, CeedBasisDestroy(&val));
    }
  }
  auto DestroyRestrictions = [](auto &restrictions)
  {
    for (auto &[ceed, restr_map] : restrictions)
    {
      for (auto &[key, val] : restr_map)
      {
        PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
      }
    }
  };
  DestroyRestrictions(restr);
  DestroyRestrictions(active_restr);
  DestroyRestrictions(interp_restr);
  DestroyRestrictions(interp_range_restr);
  basis.clear();
  restr.clear();
  active_restr.clear();
  interp_restr.clear();
  interp_range_restr.clear();
  for (std::size_t i = 0; i < ceed::internal::GetCeedObjects().size(); i++)
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[i];
    basis.emplace(ceed, ceed::GeometryObjectMap<CeedBasis>());
    restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
    active_restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
    interp_restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
    interp_range_restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
  }
}

FiniteElementSpace::FiniteElementSpace(FiniteElementSpace &fespace, const CrackSides &sides)
  : fespace(fespace.Get()), mesh(fespace.GetMesh()), aux_fespace(nullptr)
{
  MFEM_VERIFY(!fespace.IsBroken(), "Finite element space is already broken!");
  ResetCeedObjects();
  tx.UseDevice(true);
  lx.UseDevice(true);
  ly.UseDevice(true);
  InitBroken(sides);
  broken->unbroken = &fespace;
  broken->sides = &sides;
}

namespace
{

// Exchange arrays of 64-bit integers between all processes (send[p] is sent to process p),
// returning the concatenation of the arrays received, in the order of the sending
// processes, with their sizes. Collective.
std::vector<std::int64_t> ExchangeAll(MPI_Comm comm,
                                      const std::vector<std::vector<std::int64_t>> &send,
                                      std::vector<int> &recv_counts)
{
  const int nprocs = Mpi::Size(comm);
  std::vector<int> send_counts(nprocs), send_displs(nprocs, 0), recv_displs(nprocs, 0);
  recv_counts.assign(nprocs, 0);
  for (int p = 0; p < nprocs; p++)
  {
    send_counts[p] = static_cast<int>(send[p].size());
  }
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);
  for (int p = 1; p < nprocs; p++)
  {
    send_displs[p] = send_displs[p - 1] + send_counts[p - 1];
    recv_displs[p] = recv_displs[p - 1] + recv_counts[p - 1];
  }
  std::vector<std::int64_t> send_buf,
      recv_buf(recv_displs[nprocs - 1] + recv_counts[nprocs - 1]);
  send_buf.reserve(send_displs[nprocs - 1] + send_counts[nprocs - 1]);
  for (int p = 0; p < nprocs; p++)
  {
    send_buf.insert(send_buf.end(), send[p].begin(), send[p].end());
  }
  MPI_Alltoallv(send_buf.data(), send_counts.data(), send_displs.data(), MPI_INT64_T,
                recv_buf.data(), recv_counts.data(), recv_displs.data(), MPI_INT64_T, comm);
  return recv_buf;
}

}  // namespace

void FiniteElementSpace::InitBroken(const CrackSides &sides)
{
  const mfem::ParFiniteElementSpace &fespace = Get();
  MFEM_VERIFY(fespace.GetVDim() == 1,
              "Broken finite element spaces are only supported for vdim = 1!");
  MFEM_VERIFY(!fespace.IsVariableOrder(),
              "Broken finite element spaces are not supported for variable order spaces!");
  const auto &pmesh = GetParMesh();
  MPI_Comm comm = pmesh.GetComm();
  const int ne = pmesh.GetNE(), rank = Mpi::Rank(comm), nprocs = Mpi::Size(comm);
  MFEM_VERIFY(sides.split.size() == static_cast<std::size_t>(ne) &&
                  sides.version.size() == static_cast<std::size_t>(ne) &&
                  sides.side.size() == static_cast<std::size_t>(ne),
              "Invalid interior boundary sides for broken finite element space!");
  auto data = std::make_unique<BrokenData>();
  const int vsize = fespace.GetVSize(), tsize = fespace.GetTrueVSize();
  mfem::HypreParMatrix &P = *fespace.Dof_TrueDof_Matrix();
  hypre_ParCSRMatrix *hP = P;

  // Partition of the true DOFs, to find the process owning each one.
  const HYPRE_BigInt col_begin = hypre_ParCSRMatrixFirstColDiag(hP);
  std::vector<HYPRE_BigInt> col_starts(nprocs + 1);
  MPI_Allgather(&col_begin, 1, HYPRE_MPI_BIG_INT, col_starts.data(), 1, HYPRE_MPI_BIG_INT,
                comm);
  col_starts[nprocs] = hypre_ParCSRMatrixGlobalNumCols(hP);
  auto Owner = [&](HYPRE_BigInt t)
  {
    return static_cast<int>(std::upper_bound(col_starts.begin(), col_starts.end(), t) -
                            col_starts.begin()) -
           1;
  };

  // Rows of the prolongation: the global true DOF of each entry (diagonal entries, then
  // off-diagonal ones), and the true DOF of the rows which are unit vectors (L-DOFs of
  // unconstrained entities), or -1.
  std::vector<int> row_ptr(vsize + 1, 0);
  std::vector<HYPRE_BigInt> row_col, unit_tdof(vsize, -1), cmap;
  {
    P.HostRead();
    mfem::SparseMatrix diag, offd;
    HYPRE_BigInt *offd_cmap;
    P.GetDiag(diag);
    P.GetOffd(offd, offd_cmap);
    cmap.assign(offd_cmap, offd_cmap + offd.Width());
    for (int i = 0; i < vsize; i++)
    {
      int nnz = 0;
      double val = 0.0;
      HYPRE_BigInt col = -1;
      auto Add = [&](HYPRE_BigInt g, double a)
      {
        row_col.push_back(g);
        if (a != 0.0)
        {
          nnz++;
          val = a;
          col = g;
        }
      };
      if (diag.Height() > 0)
      {
        for (int k = diag.GetI()[i]; k < diag.GetI()[i + 1]; k++)
        {
          Add(col_begin + diag.GetJ()[k], diag.GetData()[k]);
        }
      }
      if (offd.Height() > 0)
      {
        for (int k = offd.GetI()[i]; k < offd.GetI()[i + 1]; k++)
        {
          Add(cmap[offd.GetJ()[k]], offd.GetData()[k]);
        }
      }
      row_ptr[i + 1] = static_cast<int>(row_col.size());
      if (nnz == 1 && std::abs(val) == 1.0)
      {
        unit_tdof[i] = col;
      }
    }
    P.HypreRead();
  }

  auto IsSplit = [&](int e, int b) { return b >= 0 && ((sides.split[e] >> b) & 1u); };

  // Versions of the true DOFs read by the elements with each side label, from the elements
  // with split entities: for the L-DOFs of unconstrained entities, the version of the
  // entity (first priority), and for the L-DOFs of hanging entities, the version of the
  // hanging entity for all true DOFs of its constraint, by increasing dimension of the
  // entity carrying the hanging entity. The master entity of a hanging entity has its
  // version, but not necessarily the entities of lower dimension of its closure (at a
  // junction of interior boundaries), which are read by unconstrained or hanging entities
  // carried by them, with the same side label. The records (true DOF, side label, priority,
  // version) are resolved on the processes owning the true DOFs.
  mfem::Array<int> dofs, entities;
  mfem::DofTransformation dof_trans;
  std::vector<std::array<std::int64_t, 4>> records;
  {
    std::vector<std::array<std::int64_t, 4>> local;
    for (int e = 0; e < ne; e++)
    {
      if (!sides.split[e])
      {
        continue;
      }
      fespace.GetElementDofs(e, dofs, dof_trans);
      fem::GetElementDofEntities(fespace, e, entities);
      for (int j = 0; j < dofs.Size(); j++)
      {
        const int b = entities[j];
        if (!IsSplit(e, b))
        {
          continue;
        }
        const int l = (dofs[j] >= 0) ? dofs[j] : -1 - dofs[j], v = sides.GetVersion(e, b);
        if (unit_tdof[l] >= 0)
        {
          local.push_back({unit_tdof[l], sides.side[e], -1, v});
        }
        else
        {
          for (int k = row_ptr[l]; k < row_ptr[l + 1]; k++)
          {
            local.push_back({row_col[k], sides.side[e], sides.GetCarrier(e, b), v});
          }
        }
      }
    }
    std::sort(local.begin(), local.end());
    local.erase(std::unique(local.begin(), local.end()), local.end());
    std::vector<std::vector<std::int64_t>> send(nprocs);
    for (const auto &r : local)
    {
      const int p = Owner(r[0]);
      if (p == rank)
      {
        records.push_back(r);
      }
      else
      {
        send[p].insert(send[p].end(), r.begin(), r.end());
      }
    }
    std::vector<int> recv_counts;
    const auto recv = ExchangeAll(comm, send, recv_counts);
    for (std::size_t i = 0; i + 3 < recv.size(); i += 4)
    {
      records.push_back({recv[i], recv[i + 1], recv[i + 2], recv[i + 3]});
    }
  }

  // Resolve the version of each owned true DOF for each side label (the record of first
  // priority), and the number of versions of each owned true DOF. The true DOFs of split
  // entities are those with a record from an unconstrained entity (an entity with master
  // DOFs is unconstrained for the elements on its unrefined side): the constraints of
  // hanging entities also involve the true DOFs of entities which are not split (at the
  // free edge of an interior boundary, for example), whose records are discarded.
  std::vector<std::array<std::int64_t, 3>> resolved;  // (true DOF, side label, version)
  std::vector<int> num_versions(tsize, 1);
  int num_conflicts = 0;
  {
    std::sort(records.begin(), records.end());
    std::vector<char> split_tdof(tsize, 0);
    for (const auto &r : records)
    {
      split_tdof[r[0] - col_begin] |= (r[2] < 0);
    }
    for (std::size_t i = 0; i < records.size();)
    {
      std::size_t j = i + 1;
      const bool split = split_tdof[records[i][0] - col_begin];
      while (j < records.size() && records[j][0] == records[i][0] &&
             records[j][1] == records[i][1])
      {
        num_conflicts +=
            (split && records[j][2] == records[i][2] && records[j][3] != records[i][3]);
        j++;
      }
      if (split)
      {
        resolved.push_back({records[i][0], records[i][1], records[i][3]});
        auto &n = num_versions[records[i][0] - col_begin];
        n = std::max(n, static_cast<int>(records[i][3]) + 1);
      }
      i = j;
    }
  }

  // Number of versions of the true DOFs of the off-diagonal columns of P.
  std::vector<int> ghost_versions(cmap.size(), 1);
  {
    if (!hypre_ParCSRMatrixCommPkg(hP))
    {
      hypre_MatvecCommPkgCreate(hP);
    }
    hypre_ParCSRCommPkg *comm_pkg = hypre_ParCSRMatrixCommPkg(hP);
    const int num_sends = hypre_ParCSRCommPkgNumSends(comm_pkg);
    const int send_size = hypre_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
    std::vector<HYPRE_Complex> send(send_size), recv(cmap.size());
    for (int k = 0; k < send_size; k++)
    {
      send[k] = num_versions[hypre_ParCSRCommPkgSendMapElmt(comm_pkg, k)];
    }
    auto *handle = hypre_ParCSRCommHandleCreate(1, comm_pkg, send.data(), recv.data());
    hypre_ParCSRCommHandleDestroy(handle);
    for (std::size_t c = 0; c < cmap.size(); c++)
    {
      ghost_versions[c] = static_cast<int>(std::lround(recv[c]));
    }
  }
  auto NumVersions = [&](HYPRE_BigInt g)
  {
    if (g >= col_begin && g < col_begin + tsize)
    {
      return num_versions[g - col_begin];
    }
    return ghost_versions[std::lower_bound(cmap.begin(), cmap.end(), g) - cmap.begin()];
  };

  // The versions needed by each element for the true DOFs with several versions in the
  // constraints of its L-DOFs which are not read directly (L-DOFs of hanging entities, or
  // with a version from other entities), for its side label, queried from the processes
  // owning the true DOFs.
  auto IsDirect = [&](int e, int b, int l) { return IsSplit(e, b) && unit_tdof[l] >= 0; };
  std::vector<std::array<std::int64_t, 3>> needed;  // (true DOF, side label, version)
  {
    std::vector<std::array<std::int64_t, 2>> keys;
    for (int e = 0; e < ne; e++)
    {
      fespace.GetElementDofs(e, dofs, dof_trans);
      if (sides.split[e])
      {
        fem::GetElementDofEntities(fespace, e, entities);
      }
      for (int j = 0; j < dofs.Size(); j++)
      {
        const int l = (dofs[j] >= 0) ? dofs[j] : -1 - dofs[j];
        if (sides.split[e] && IsDirect(e, entities[j], l))
        {
          continue;
        }
        for (int k = row_ptr[l]; k < row_ptr[l + 1]; k++)
        {
          if (NumVersions(row_col[k]) > 1)
          {
            keys.push_back({row_col[k], sides.side[e]});
          }
        }
      }
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    auto Lookup = [&](std::int64_t t, std::int64_t side) -> std::int64_t
    {
      const std::array<std::int64_t, 3> key = {t, side, -1};
      auto it = std::lower_bound(resolved.begin(), resolved.end(), key);
      return (it != resolved.end() && (*it)[0] == t && (*it)[1] == side) ? (*it)[2] : -1;
    };
    std::vector<std::vector<std::int64_t>> send(nprocs);
    std::vector<std::vector<std::array<std::int64_t, 2>>> queried(nprocs);
    for (const auto &key : keys)
    {
      const int p = Owner(key[0]);
      if (p == rank)
      {
        needed.push_back({key[0], key[1], Lookup(key[0], key[1])});
      }
      else
      {
        send[p].insert(send[p].end(), key.begin(), key.end());
        queried[p].push_back(key);
      }
    }
    std::vector<int> recv_counts;
    const auto queries = ExchangeAll(comm, send, recv_counts);
    std::vector<std::vector<std::int64_t>> replies(nprocs);
    for (int p = 0, pos = 0; p < nprocs; p++)
    {
      for (int i = 0; i < recv_counts[p]; i += 2, pos += 2)
      {
        replies[p].push_back(Lookup(queries[pos], queries[pos + 1]));
      }
    }
    const auto answers = ExchangeAll(comm, replies, recv_counts);
    for (int p = 0, pos = 0; p < nprocs; p++)
    {
      for (const auto &key : queried[p])
      {
        needed.push_back({key[0], key[1], answers[pos++]});
      }
    }
    std::sort(needed.begin(), needed.end());
  }

  // Each element DOF reads the original L-DOF when all true DOFs of its row are read in
  // their original version, and otherwise a copy of the L-DOF with the versions of its true
  // DOFs, shared by all element DOFs reading the same versions.
  std::vector<int> copy_ldofs, copy_offsets(1, 0);
  std::vector<std::uint8_t> copy_versions, ver;
  std::unordered_map<std::string, int> copy_index;
  int num_unresolved = 0;
  data->override_offsets.resize(ne + 1, 0);
  for (int e = 0; e < ne; e++)
  {
    data->override_offsets[e] = static_cast<int>(data->override_local.size());
    fespace.GetElementDofs(e, dofs, dof_trans);
    if (sides.split[e])
    {
      fem::GetElementDofEntities(fespace, e, entities);
    }
    for (int j = 0; j < dofs.Size(); j++)
    {
      const int l = (dofs[j] >= 0) ? dofs[j] : -1 - dofs[j];
      const int n = row_ptr[l + 1] - row_ptr[l];
      ver.assign(n, 0);
      bool any = false;
      if (sides.split[e] && IsDirect(e, entities[j], l))
      {
        const int v = sides.GetVersion(e, entities[j]);
        for (int k = 0; k < n && v > 0; k++)
        {
          if (row_col[row_ptr[l] + k] == unit_tdof[l])
          {
            ver[k] = static_cast<std::uint8_t>(v);
            any = true;
          }
        }
      }
      else
      {
        for (int k = 0; k < n; k++)
        {
          const HYPRE_BigInt g = row_col[row_ptr[l] + k];
          if (NumVersions(g) <= 1)
          {
            continue;
          }
          const std::array<std::int64_t, 3> key = {g, sides.side[e], -1};
          const auto it = std::lower_bound(needed.begin(), needed.end(), key);
          if (it == needed.end() || (*it)[0] != g || (*it)[1] != sides.side[e] ||
              (*it)[2] < 0 || (*it)[2] >= NumVersions(g))
          {
            num_unresolved++;
          }
          else if ((*it)[2] > 0)
          {
            ver[k] = static_cast<std::uint8_t>((*it)[2]);
            any = true;
          }
        }
      }
      if (!any)
      {
        continue;
      }
      std::string key(reinterpret_cast<const char *>(&l), sizeof(l));
      key.append(reinterpret_cast<const char *>(ver.data()), ver.size());
      auto [it, inserted] =
          copy_index.try_emplace(std::move(key), static_cast<int>(copy_ldofs.size()));
      if (inserted)
      {
        copy_ldofs.push_back(l);
        copy_versions.insert(copy_versions.end(), ver.begin(), ver.end());
        copy_offsets.push_back(static_cast<int>(copy_versions.size()));
      }
      data->override_local.push_back(j);
      data->override_ldof.push_back(vsize + it->second);
    }
  }
  data->override_offsets[ne] = static_cast<int>(data->override_local.size());
  {
    // Inconsistent versions from records of the same priority, or versions which could not
    // be resolved, would give wrong broken spaces.
    int counts[2] = {num_conflicts, num_unresolved};
    Mpi::GlobalSum(2, counts, comm);
    MFEM_VERIFY(counts[0] == 0 && counts[1] == 0,
                "Inconsistent (" << counts[0] << ") or missing (" << counts[1]
                                 << ") versions of the DOFs of interior boundaries for a "
                                    "broken finite element space!");
  }

  // The broken prolongation.
  data->vsize = vsize + static_cast<int>(copy_ldofs.size());
  fem::BuildBrokenProlongation(pmesh, P, copy_ldofs, copy_offsets, copy_versions,
                               num_versions, data->P);
  data->tsize = data->P.P->Width();
  data->global_vsize = data->P.P->GetGlobalNumRows();
  data->global_tsize = data->P.P->GetGlobalNumCols();

  // For the interpolators: a side label of the elements reading each version of each owned
  // true DOF, from the resolved records for the split true DOFs, and from the local
  // elements reading the true DOF (all of which read its only version) for the others.
  data->version_label.assign(data->tsize, -1);
  for (const auto &r : resolved)
  {
    auto &label = data->version_label[data->P.GetVersionTrueDof(
        static_cast<int>(r[0] - col_begin), static_cast<int>(r[2]))];
    label = (label < 0) ? r[1] : label;
  }
  for (int e = 0; e < ne; e++)
  {
    fespace.GetElementDofs(e, dofs, dof_trans);
    for (auto d : dofs)
    {
      const int l = (d >= 0) ? d : -1 - d;
      const HYPRE_BigInt t = unit_tdof[l];
      if (t >= col_begin && t < col_begin + tsize && num_versions[t - col_begin] == 1)
      {
        auto &label = data->version_label[t - col_begin];
        label = (label < 0) ? sides.side[e] : label;
      }
    }
  }
  {
    int num_unlabeled = static_cast<int>(
        std::count(data->version_label.begin(), data->version_label.end(), -1));
    Mpi::GlobalSum(1, &num_unlabeled, comm);
    MFEM_VERIFY(num_unlabeled == 0, num_unlabeled
                                        << " true DOF(s) of a broken finite element space "
                                           "are read by no element!");
  }
  data->resolved = std::move(resolved);
  data->col_starts = std::move(col_starts);
  broken = std::move(data);
}

void FiniteElementSpace::GetElementDofs(int e, mfem::Array<int> &dofs,
                                        mfem::DofTransformation &dof_trans) const
{
  Get().GetElementDofs(e, dofs, dof_trans);
  if (!broken)
  {
    return;
  }
  for (int k = broken->override_offsets[e]; k < broken->override_offsets[e + 1]; k++)
  {
    const int j = broken->override_local[k];
    dofs[j] = (dofs[j] >= 0) ? broken->override_ldof[k] : -1 - broken->override_ldof[k];
  }
}

void FiniteElementSpace::GetBdrElementDofs(int be, mfem::Array<int> &dofs,
                                           mfem::DofTransformation &dof_trans) const
{
  Get().GetBdrElementDofs(be, dofs, dof_trans);
  if (!broken)
  {
    return;
  }
  // The DOFs of the boundary element are among those of the adjacent element, which the
  // broken space reads in the versions of its side.
  int e, info;
  GetParMesh().GetBdrElementAdjacentElement(be, e, info);
  mfem::Array<int> edofs, bdofs;
  mfem::DofTransformation edof_trans;
  Get().GetElementDofs(e, edofs, edof_trans);
  GetElementDofs(e, bdofs, edof_trans);
  for (auto &d : dofs)
  {
    const int l = (d >= 0) ? d : -1 - d;
    int j = 0;
    while (j < edofs.Size() && ((edofs[j] >= 0) ? edofs[j] : -1 - edofs[j]) != l)
    {
      j++;
    }
    MFEM_ASSERT(j < edofs.Size(),
                "Boundary element DOF not found in its adjacent element!");
    const int lb = (bdofs[j] >= 0) ? bdofs[j] : -1 - bdofs[j];
    d = (d >= 0) ? lb : -1 - lb;
  }
}

void FiniteElementSpace::GetEssentialTrueDofs(const mfem::Array<int> &bdr_marker,
                                              mfem::Array<int> &ess_tdof_list) const
{
  if (!broken)
  {
    Get().GetEssentialTrueDofs(bdr_marker, ess_tdof_list);
    return;
  }
  // Mark the L-DOFs of the boundary elements in the versions of their adjacent elements. A
  // true DOF is essential when a marked L-DOF reads it, on any process (|P|ᵀ sums the marks
  // of every process, where the unbroken space synchronizes them across processes).
  const auto &mesh = GetParMesh();
  Vector marker(broken->vsize), tmarker(broken->tsize);
  marker.UseDevice(false);
  tmarker.UseDevice(false);
  marker = 0.0;
  mfem::Array<int> dofs;
  mfem::DofTransformation dof_trans;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (!bdr_marker[mesh.GetBdrAttribute(be) - 1])
    {
      continue;
    }
    GetBdrElementDofs(be, dofs, dof_trans);
    for (auto d : dofs)
    {
      marker((d >= 0) ? d : -1 - d) = 1.0;
    }
  }
  broken->P.P->AbsMultTranspose(1.0, marker, 0.0, tmarker);
  ess_tdof_list.SetSize(0);
  for (int t = 0; t < broken->tsize; t++)
  {
    if (tmarker(t) > 0.0)
    {
      ess_tdof_list.Append(t);
    }
  }
}

std::unique_ptr<mfem::HypreParMatrix>
FiniteElementSpace::BuildVersionMap(VersionMap type) const
{
  MFEM_VERIFY(broken, "Version maps are only available for broken spaces!");
  const int tsize = Get().GetTrueVSize();
  auto S = std::make_unique<mfem::SparseMatrix>(tsize, broken->tsize);
  for (int t = 0; t < tsize; t++)
  {
    const int n = GetNumVersions(t);
    for (int v = 0; v < ((type == VersionMap::FIRST) ? 1 : n); v++)
    {
      S->Set(t, GetVersionTrueDof(t, v), (type == VersionMap::SIDES && n == 1) ? 2.0 : 1.0);
    }
  }
  return linalg::BlockDiagonalParMatrix(GetComm(), Get().GlobalTrueVSize(),
                                        broken->global_tsize, Get().GetTrueDofOffsets(),
                                        broken->P.col_offsets.GetData(), std::move(S));
}

std::vector<HYPRE_BigInt> FiniteElementSpace::GetVersionedTrueDofs(
    const std::vector<std::array<std::int64_t, 2>> &keys) const
{
  MFEM_VERIFY(broken, "Versioned true DOFs are only available for broken spaces!");
  MPI_Comm comm = GetComm();
  const int rank = Mpi::Rank(comm), nprocs = Mpi::Size(comm);
  const auto &col_starts = broken->col_starts;
  const HYPRE_BigInt col_begin = col_starts[rank],
                     col_begin_b =
                         broken->P.col_offsets[HYPRE_AssumedPartitionCheck() ? 0 : rank];
  auto Lookup = [&](std::int64_t t, std::int64_t label) -> HYPRE_BigInt
  {
    const int tl = static_cast<int>(t - col_begin);
    if (GetNumVersions(tl) == 1)
    {
      return col_begin_b + tl;
    }
    const std::array<std::int64_t, 3> key = {t, label, -1};
    auto it = std::lower_bound(broken->resolved.begin(), broken->resolved.end(), key);
    if (it == broken->resolved.end() || (*it)[0] != t || (*it)[1] != label)
    {
      return -1;
    }
    return col_begin_b + GetVersionTrueDof(tl, static_cast<int>((*it)[2]));
  };

  // Queries to the processes owning the true DOFs, and their replies in the order of the
  // keys.
  std::vector<HYPRE_BigInt> out(keys.size(), -1);
  std::vector<std::vector<std::int64_t>> send(nprocs);
  std::vector<std::vector<std::size_t>> queried(nprocs);
  for (std::size_t i = 0; i < keys.size(); i++)
  {
    const int p = static_cast<int>(
                      std::upper_bound(col_starts.begin(), col_starts.end(), keys[i][0]) -
                      col_starts.begin()) -
                  1;
    if (p == rank)
    {
      out[i] = Lookup(keys[i][0], keys[i][1]);
    }
    else
    {
      send[p].insert(send[p].end(), keys[i].begin(), keys[i].end());
      queried[p].push_back(i);
    }
  }
  std::vector<int> recv_counts;
  const auto queries = ExchangeAll(comm, send, recv_counts);
  std::vector<std::vector<std::int64_t>> replies(nprocs);
  for (int p = 0, pos = 0; p < nprocs; p++)
  {
    for (int i = 0; i < recv_counts[p]; i += 2, pos += 2)
    {
      replies[p].push_back(Lookup(queries[pos], queries[pos + 1]));
    }
  }
  const auto answers = ExchangeAll(comm, replies, recv_counts);
  for (int p = 0, pos = 0; p < nprocs; p++)
  {
    for (auto i : queried[p])
    {
      out[i] = answers[pos++];
    }
  }
  return out;
}

std::unique_ptr<mfem::HypreParMatrix>
FiniteElementSpace::BuildBrokenInterpolator(const mfem::HypreParMatrix &G,
                                            const FiniteElementSpace &trial_fespace,
                                            const FiniteElementSpace &test_fespace)
{
  MFEM_VERIFY(trial_fespace.broken && test_fespace.broken &&
                  (&trial_fespace.GetMesh() != &test_fespace.GetMesh() ||
                   trial_fespace.broken->sides == test_fespace.broken->sides),
              "Interpolators between broken spaces require both spaces to be broken with "
              "the same interior boundary sides!");
  const auto &test = *test_fespace.broken;
  MFEM_VERIFY(G.Height() == test_fespace.Get().GetTrueVSize() &&
                  G.GetGlobalNumCols() == trial_fespace.Get().GlobalTrueVSize(),
              "Invalid interpolator between the underlying spaces of broken spaces!");
  MPI_Comm comm = test_fespace.GetComm();
  G.HostRead();
  mfem::SparseMatrix diag, offd;
  HYPRE_BigInt *cmap;
  G.GetDiag(diag);
  G.GetOffd(offd, cmap);
  const HYPRE_BigInt col_begin = hypre_ParCSRMatrixFirstColDiag((hypre_ParCSRMatrix *)G);
  const int nrows = G.Height();
  // The nonzero entries of a row (an assembled interpolator stores zeros for the DOFs of
  // the element which are not in the closure of the entity of the interpolated DOF, which
  // elements on one side of the entity need not share).
  auto ForEachEntry = [&](int t, auto &&f)
  {
    if (diag.Height() > 0)
    {
      for (int k = diag.GetI()[t]; k < diag.GetI()[t + 1]; k++)
      {
        if (diag.GetData()[k] != 0.0)
        {
          f(col_begin + diag.GetJ()[k], diag.GetData()[k]);
        }
      }
    }
    if (offd.Height() > 0)
    {
      for (int k = offd.GetI()[t]; k < offd.GetI()[t + 1]; k++)
      {
        if (offd.GetData()[k] != 0.0)
        {
          f(cmap[offd.GetJ()[k]], offd.GetData()[k]);
        }
      }
    }
  };

  // The versions of the trial DOFs read by the elements reading each version of each test
  // DOF, from the trial space.
  std::vector<std::array<std::int64_t, 2>> keys;
  for (int t = 0; t < nrows; t++)
  {
    for (int v = 0; v < test_fespace.GetNumVersions(t); v++)
    {
      const std::int64_t label = test.version_label[test_fespace.GetVersionTrueDof(t, v)];
      ForEachEntry(t, [&](HYPRE_BigInt s, double) { keys.push_back({s, label}); });
    }
  }
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  const auto cols = trial_fespace.GetVersionedTrueDofs(keys);
  int num_missing = static_cast<int>(std::count(cols.begin(), cols.end(), -1));
  Mpi::GlobalSum(1, &num_missing, comm);
  MFEM_VERIFY(num_missing == 0, num_missing << " versioned true DOF(s) not found for the "
                                               "interpolator between broken spaces!");

  // Rows of the broken interpolator, the original true DOFs then the copies.
  std::vector<std::vector<std::pair<HYPRE_BigInt, double>>> rows(test.tsize);
  for (int t = 0; t < nrows; t++)
  {
    for (int v = 0; v < test_fespace.GetNumVersions(t); v++)
    {
      const int r = test_fespace.GetVersionTrueDof(t, v);
      const std::int64_t label = test.version_label[r];
      ForEachEntry(t,
                   [&](HYPRE_BigInt s, double a)
                   {
                     const std::array<std::int64_t, 2> key = {s, label};
                     const auto it = std::lower_bound(keys.begin(), keys.end(), key);
                     rows[r].emplace_back(cols[it - keys.begin()], a);
                   });
    }
  }
  G.HypreRead();  // Return G to hypre's memory space after reading it on the host
  std::vector<int> I(test.tsize + 1, 0);
  std::vector<HYPRE_BigInt> J;
  std::vector<double> data;
  for (int r = 0; r < test.tsize; r++)
  {
    for (const auto &[j, a] : rows[r])
    {
      J.push_back(j);
      data.push_back(a);
    }
    I[r + 1] = static_cast<int>(J.size());
  }
  auto Gb = std::make_unique<mfem::HypreParMatrix>(
      comm, test.tsize, test.global_tsize, trial_fespace.broken->global_tsize, I.data(),
      J.data(), data.data(), test_fespace.GetTrueDofOffsets(),
      trial_fespace.GetTrueDofOffsets());
  hypre_ParCSRMatrixSetNumNonzeros(*Gb);
  return Gb;
}

CeedBasis FiniteElementSpace::BuildCeedBasis(const mfem::FiniteElementSpace &fespace,
                                             Ceed ceed, mfem::Geometry::Type geom)
{
  // Find the appropriate integration rule for the element.
  mfem::IsoparametricTransformation T;
  const mfem::FiniteElement *fe_nodal =
      fespace.GetMesh()->GetNodalFESpace()->FEColl()->FiniteElementForGeometry(geom);
  if (!fe_nodal)
  {
    fe_nodal =
        fespace.GetMesh()->GetNodalFESpace()->FEColl()->TraceFiniteElementForGeometry(geom);
  }
  T.SetFE(fe_nodal);
  const int q_order = fem::DefaultIntegrationOrder::Get(T);
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(geom, q_order);

  // Build the libCEED basis.
  CeedBasis val;
  const mfem::FiniteElement *fe = fespace.FEColl()->FiniteElementForGeometry(geom);
  if (!fe)
  {
    fe = fespace.FEColl()->TraceFiniteElementForGeometry(geom);
  }
  const int vdim = fespace.GetVDim();
  ceed::InitBasis(*fe, ir, vdim, ceed, &val);
  return val;
}

CeedElemRestriction FiniteElementSpace::BuildCeedElemRestriction(
    const mfem::FiniteElementSpace &fespace, Ceed ceed, mfem::Geometry::Type geom,
    const std::vector<int> &indices, bool is_interp, bool is_interp_range)
{
  // Construct the libCEED element restriction for this element type.
  CeedElemRestriction val;
  const bool use_bdr = (mfem::Geometry::Dimension[geom] != fespace.GetMesh()->Dimension());
  ceed::InitRestriction(fespace, indices, use_bdr, is_interp, is_interp_range, ceed, &val);
  return val;
}

const Operator &FiniteElementSpace::BuildDiscreteInterpolator() const
{
  // Allow finite element spaces to be swapped in their order (intended as deriv(aux) ->
  // primal). G is always partially assembled.
  const int dim = Dimension();
  const bool forward =
      (GetFEColl().GetMapType(dim) == aux_fespace->GetFEColl().GetDerivMapType(dim));
  const bool swap = !forward && (aux_fespace->GetFEColl().GetMapType(dim) ==
                                 GetFEColl().GetDerivMapType(dim));
  MFEM_VERIFY(!swap, "Incorrect order for primal/auxiliary (test/trial) spaces in discrete "
                     "interpolator construction!");
  MFEM_VERIFY(forward, "Unsupported trial/test FE spaces for FiniteElementSpace discrete "
                       "interpolator!");
  const FiniteElementSpace &trial_fespace = !swap ? *aux_fespace : *this;
  const FiniteElementSpace &test_fespace = !swap ? *this : *aux_fespace;
  if (broken || aux_fespace->IsBroken())
  {
    MFEM_VERIFY(
        broken && aux_fespace->IsBroken(),
        "Discrete interpolators require both or neither of the spaces to be broken!");
    // Interpolators between broken spaces are partially assembled as for unbroken spaces,
    // with the element restrictions of the broken spaces (each element reads the copies of
    // its side) and the restriction of the broken test space (see GetRestrictionWeights),
    // unless some copy of the test space is determined by none of the rows of its
    // prolongation: then they are assembled from the one between the underlying spaces.
    // The 2D scalar curl uses MFEM's native assembly, which knows no broken spaces.
    const bool native =
        (trial_fespace.GetFEColl().GetMapType(dim) == mfem::FiniteElement::H_CURL &&
         test_fespace.GetFEColl().GetMapType(dim) == mfem::FiniteElement::INTEGRAL);
    Vector w;
    if (native || !test_fespace.GetRestrictionWeights(w))
    {
      const auto *G_u = dynamic_cast<const ParOperator *>(
          &GetUnbroken().GetDiscreteInterpolator(aux_fespace->GetUnbroken()));
      MFEM_VERIFY(G_u, "Unexpected discrete interpolator type for the underlying spaces of "
                       "broken spaces!");
      G = BuildBrokenInterpolator(G_u->ParallelAssemble(), trial_fespace, test_fespace);
      return *G;
    }
  }
  const auto aux_map_type = trial_fespace.GetFEColl().GetMapType(dim);
  const auto primal_map_type = test_fespace.GetFEColl().GetMapType(dim);
  if (aux_map_type == mfem::FiniteElement::VALUE &&
      primal_map_type == mfem::FiniteElement::H_CURL)
  {
    // Discrete gradient interpolator.
    DiscreteLinearOperator interp(trial_fespace, test_fespace);
    interp.AddDomainInterpolator<GradientInterpolator>();
    G = std::make_unique<ParOperator>(interp.PartialAssemble(), trial_fespace, test_fespace,
                                      true);
  }
  else if (aux_map_type == mfem::FiniteElement::H_CURL &&
           primal_map_type == mfem::FiniteElement::H_DIV)
  {
    // Discrete curl interpolator (3D: H(curl) → H(div)).
    DiscreteLinearOperator interp(trial_fespace, test_fespace);
    interp.AddDomainInterpolator<CurlInterpolator>();
    G = std::make_unique<ParOperator>(interp.PartialAssemble(), trial_fespace, test_fespace,
                                      true);
  }
  else if (aux_map_type == mfem::FiniteElement::H_CURL &&
           primal_map_type == mfem::FiniteElement::INTEGRAL)
  {
    // Discrete curl interpolator (2D: H(curl) → L2, scalar curl). Uses MFEM's native
    // assembly because libCEED does not support partial assembly for this operator type.
    auto *trial_pfes = const_cast<mfem::ParFiniteElementSpace *>(&trial_fespace.Get());
    auto *test_pfes = const_cast<mfem::ParFiniteElementSpace *>(&test_fespace.Get());
    mfem::DiscreteLinearOperator interp(trial_pfes, test_pfes);
    interp.AddDomainInterpolator(new mfem::CurlInterpolator);
    interp.Assemble();
    interp.Finalize();
    G = std::make_unique<ParOperator>(std::unique_ptr<mfem::SparseMatrix>(interp.LoseMat()),
                                      trial_fespace, test_fespace, true);
  }
  else if (aux_map_type == mfem::FiniteElement::H_DIV &&
           primal_map_type == mfem::FiniteElement::INTEGRAL)
  {
    // Discrete divergence interpolator.
    DiscreteLinearOperator interp(trial_fespace, test_fespace);
    interp.AddDomainInterpolator<DivergenceInterpolator>();
    G = std::make_unique<ParOperator>(interp.PartialAssemble(), trial_fespace, test_fespace,
                                      true);
  }
  else
  {
    MFEM_ABORT(
        "Unsupported trial/test FE spaces for FiniteElementSpace discrete interpolator!");
  }

  return *G;
}

namespace
{

// Transfer between the true DOFs of a space on a coarse mesh and on the mesh refined from
// it, as R_f T P_c with T assembled element by element from the coarse-to-fine
// interpolation, including the DOF transformations. MFEM's assembled refinement matrix is
// not used: for ND order 2 on tetrahedra it does not commute with the discrete gradient.
std::unique_ptr<mfem::HypreParMatrix>
BuildTrueTransferMatrix(const FiniteElementSpace &coarse, const FiniteElementSpace &fine)
{
  const mfem::ParFiniteElementSpace &cfes = coarse.Get(), &ffes = fine.Get();
  MFEM_VERIFY(!coarse.IsBroken() && !fine.IsBroken() && ffes.GetVDim() == 1 &&
                  !ffes.IsVariableOrder() && ffes.FEColl() == cfes.FEColl(),
              "Unsupported spaces for the assembled transfer between meshes!");
  auto &fmesh = const_cast<mfem::ParMesh &>(fine.GetParMesh());
  const auto &cmesh = coarse.GetParMesh();
  MFEM_VERIFY(
      fmesh.GetLastOperation() == mfem::Mesh::REFINE,
      "The fine mesh of an assembled transfer must have been refined from the coarse "
      "mesh!");
  const mfem::CoarseFineTransformations &cf = fmesh.GetRefinementTransforms();
  MFEM_VERIFY(cf.embeddings.Size() == fmesh.GetNE(),
              "Invalid coarse-to-fine transformations for the assembled transfer!");
  auto T = std::make_unique<mfem::SparseMatrix>(ffes.GetVSize(), cfes.GetVSize());
  mfem::Array<int> fdofs, cdofs;
  mfem::DofTransformation fdt, cdt;
  mfem::IsoparametricTransformation isotr;
  mfem::DenseMatrix lP, M;
  mfem::Vector v, w;
  for (int k = 0; k < fmesh.GetNE(); k++)
  {
    const auto &emb = cf.embeddings[k];
    const auto geom = fmesh.GetElementGeometry(k);
    MFEM_VERIFY(emb.parent >= 0 && emb.parent < cmesh.GetNE() &&
                    cmesh.GetElementGeometry(emb.parent) == geom,
                "Unsupported refinement for the assembled transfer between meshes!");
    const mfem::FiniteElement &fe = *ffes.GetFE(k);
    isotr.SetIdentityTransformation(geom);
    isotr.SetPointMat(cf.point_matrices[geom](emb.matrix));
    fe.GetLocalInterpolation(isotr, lP);
    ffes.GetElementDofs(k, fdofs, fdt);
    cfes.GetElementDofs(emb.parent, cdofs, cdt);
    const int n = fe.GetDof();
    M.SetSize(n);
    v.SetSize(n);
    w.SetSize(n);
    for (int j = 0; j < n; j++)
    {
      v = 0.0;
      v(j) = 1.0;
      if (!cdt.IsIdentity())
      {
        cdt.InvTransformPrimal(v);
      }
      lP.Mult(v, w);
      if (!fdt.IsIdentity())
      {
        fdt.TransformPrimal(w);
      }
      M.SetCol(j, w);
    }
    T->SetSubMatrix(fdofs, cdofs, M);
  }
  auto T_par = linalg::BlockDiagonalParMatrix(fine.GetComm(), ffes.GlobalVSize(),
                                              cfes.GlobalVSize(), ffes.GetDofOffsets(),
                                              cfes.GetDofOffsets(), std::move(T));
  mfem::HypreParMatrix R(fine.GetComm(), ffes.GlobalTrueVSize(), ffes.GlobalVSize(),
                         ffes.GetTrueDofOffsets(), ffes.GetDofOffsets(),
                         const_cast<mfem::SparseMatrix *>(ffes.GetRestrictionMatrix()));
  hypre_ParCSRMatrix *TP = hypre_ParCSRMatMat(*T_par, *cfes.Dof_TrueDof_Matrix());
  auto RTP = std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatMat(R, TP), true);
  hypre_ParCSRMatrixDestroy(TP);
  hypre_ParCSRMatrixSetNumNonzeros(*RTP);
  return RTP;
}

}  // namespace

const Operator &FiniteElementSpaceHierarchy::BuildProlongationAtLevel(std::size_t l) const
{
  Vector w;
  // P is always partially assembled.
  MFEM_VERIFY(l + 1 < GetNumLevels(),
              "Can only construct a finite element space prolongation with more than one "
              "space in the hierarchy!");
  if (&fespaces[l]->GetMesh() != &fespaces[l + 1]->GetMesh())
  {
    if (fespaces[l]->IsBroken() || fespaces[l + 1]->IsBroken())
    {
      // h-prolongation between broken spaces, assembled from the transfer between the
      // underlying spaces on their true DOFs. The sides of the fine mesh must be inherited
      // from the coarse mesh (see Mesh::InheritCrackSides) for the side labels of the two
      // spaces to agree.
      P[l] = FiniteElementSpace::BuildBrokenInterpolator(
          *BuildTrueTransferMatrix(fespaces[l]->GetUnbroken(),
                                   fespaces[l + 1]->GetUnbroken()),
          *fespaces[l], *fespaces[l + 1]);
    }
    else
    {
      P[l] = std::make_unique<ParOperator>(
          std::make_unique<mfem::TransferOperator>(*fespaces[l], *fespaces[l + 1]),
          *fespaces[l], *fespaces[l + 1], true);
    }
  }
  else if ((fespaces[l]->IsBroken() || fespaces[l + 1]->IsBroken()) &&
           !fespaces[l + 1]->GetRestrictionWeights(w))
  {
    // Prolongation between broken spaces when some copy of the fine space is determined by
    // none of the rows of its prolongation: assembled from the one between the underlying
    // spaces.
    MFEM_VERIFY(fespaces[l]->IsBroken() && fespaces[l + 1]->IsBroken(),
                "Prolongation requires both or neither of the spaces to be broken!");
    const auto &coarse = fespaces[l]->GetUnbroken(), &fine = fespaces[l + 1]->GetUnbroken();
    DiscreteLinearOperator p(coarse, fine);
    p.AddDomainInterpolator<IdentityInterpolator>();
    ParOperator P_u(p.PartialAssemble(), coarse, fine, true);
    P[l] = FiniteElementSpace::BuildBrokenInterpolator(P_u.ParallelAssemble(), *fespaces[l],
                                                       *fespaces[l + 1]);
  }
  else
  {
    // Partially assembled, also between broken spaces (with the element restrictions of the
    // broken spaces, each element reading the copies of its side, and the restriction of
    // the broken fine space, see FiniteElementSpace::GetRestrictionWeights): no matrix.
    MFEM_VERIFY(fespaces[l]->IsBroken() == fespaces[l + 1]->IsBroken(),
                "Prolongation requires both or neither of the spaces to be broken!");
    DiscreteLinearOperator p(*fespaces[l], *fespaces[l + 1]);
    p.AddDomainInterpolator<IdentityInterpolator>();
    P[l] = std::make_unique<ParOperator>(p.PartialAssemble(), *fespaces[l],
                                         *fespaces[l + 1], true);
  }

  return *P[l];
}

}  // namespace palace
