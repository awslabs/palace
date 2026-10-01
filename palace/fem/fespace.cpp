// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fespace.hpp"

#include <algorithm>
#include <array>
#include <string>
#include <unordered_map>
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
                                           const std::vector<int> &indices) const
{
  auto it = restr.find(ceed);
  MFEM_ASSERT(it != restr.end(), "Unknown Ceed context in GetCeedElemRestriction!");
  auto &restr_map = it->second;
  auto restr_it = restr_map.find(geom);
  if (restr_it != restr_map.end())
  {
    return restr_it->second;
  }
  if (broken)
  {
    CeedElemRestriction val;
    MFEM_VERIFY(mfem::Geometry::Dimension[geom] == GetParMesh().Dimension(),
                "Only domain element restrictions are available for a broken space!");
    ceed::ElementDofRemap remap;
    remap.offsets = broken->override_offsets.data();
    remap.local = broken->override_local.data();
    remap.ldof = broken->override_ldof.data();
    remap.l_size = broken->vsize;
    ceed::InitRestriction(Get(), indices, false, false, false, ceed, &val, &remap);
    return restr_map.emplace(geom, val).first->second;
  }
  return restr_map.emplace(geom, BuildCeedElemRestriction(*this, ceed, geom, indices))
      .first->second;
}

CeedElemRestriction
FiniteElementSpace::GetInterpCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                                 const std::vector<int> &indices) const
{
  const mfem::FiniteElement &fe = *GetFEColl().FiniteElementForGeometry(geom);
  if (!HasUniqueInterpRestriction(fe))
  {
    return GetCeedElemRestriction(ceed, geom, indices);
  }
  MFEM_VERIFY(!broken, "Interpolation restrictions are not available for a broken space!");
  auto it = interp_restr.find(ceed);
  MFEM_ASSERT(it != interp_restr.end(),
              "Unknown Ceed context in GetInterpCeedElemRestriction!");
  auto &restr_map = it->second;
  auto restr_it = restr_map.find(geom);
  if (restr_it != restr_map.end())
  {
    return restr_it->second;
  }
  return restr_map
      .emplace(geom, BuildCeedElemRestriction(*this, ceed, geom, indices, true, false))
      .first->second;
}

CeedElemRestriction
FiniteElementSpace::GetInterpRangeCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                                      const std::vector<int> &indices) const
{
  const mfem::FiniteElement &fe = *GetFEColl().FiniteElementForGeometry(geom);
  if (!HasUniqueInterpRangeRestriction(fe))
  {
    return GetInterpCeedElemRestriction(ceed, geom, indices);
  }
  MFEM_VERIFY(!broken, "Interpolation restrictions are not available for a broken space!");
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
      .emplace(geom, BuildCeedElemRestriction(*this, ceed, geom, indices, true, true))
      .first->second;
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
  for (auto &[ceed, restr_map] : restr)
  {
    for (auto &[key, val] : restr_map)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
    }
  }
  for (auto &[ceed, restr_map] : interp_restr)
  {
    for (auto &[key, val] : restr_map)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
    }
  }
  for (auto &[ceed, restr_map] : interp_range_restr)
  {
    for (auto &[key, val] : restr_map)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
    }
  }
  basis.clear();
  restr.clear();
  interp_restr.clear();
  interp_range_restr.clear();
  for (std::size_t i = 0; i < ceed::internal::GetCeedObjects().size(); i++)
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[i];
    basis.emplace(ceed, ceed::GeometryObjectMap<CeedBasis>());
    restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
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
          MFEM_ASSERT(it != needed.end() && (*it)[0] == g && (*it)[1] == sides.side[e],
                      "Missing version query for a broken space true DOF!");
          if ((*it)[2] < 0)
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
    int counts[2] = {num_conflicts, num_unresolved};
    Mpi::GlobalSum(2, counts, comm);
    if (counts[0] > 0 || counts[1] > 0)
    {
      Mpi::Warning(comm,
                   "Inconsistent ({:d}) or missing ({:d}) versions of degrees of freedom "
                   "of interior boundaries for a broken finite element space!\n",
                   counts[0], counts[1]);
    }
  }

  // The broken prolongation.
  data->vsize = vsize + static_cast<int>(copy_ldofs.size());
  fem::BuildBrokenProlongation(pmesh, P, copy_ldofs, copy_offsets, copy_versions,
                               num_versions, data->P);
  data->tsize = data->P.P->Width();
  data->global_vsize = data->P.P->GetGlobalNumRows();
  data->global_tsize = data->P.P->GetGlobalNumCols();
  broken = std::move(data);
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

const Operator &FiniteElementSpaceHierarchy::BuildProlongationAtLevel(std::size_t l) const
{
  // P is always partially assembled.
  MFEM_VERIFY(l + 1 < GetNumLevels(),
              "Can only construct a finite element space prolongation with more than one "
              "space in the hierarchy!");
  if (&fespaces[l]->GetMesh() != &fespaces[l + 1]->GetMesh())
  {
    P[l] = std::make_unique<ParOperator>(
        std::make_unique<mfem::TransferOperator>(*fespaces[l], *fespaces[l + 1]),
        *fespaces[l], *fespaces[l + 1], true);
  }
  else
  {
    DiscreteLinearOperator p(*fespaces[l], *fespaces[l + 1]);
    p.AddDomainInterpolator<IdentityInterpolator>();
    P[l] = std::make_unique<ParOperator>(p.PartialAssemble(), *fespaces[l],
                                         *fespaces[l + 1], true);
  }

  return *P[l];
}

}  // namespace palace
