// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "sliverpatch.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/libceed/operator.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

inline void Cross(const double *a, const double *b, double *c)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

inline double Dot(const double *a, const double *b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Marker of the local rows of residual updates which are skipped (essential true dofs).
constexpr int kSkip = std::numeric_limits<int>::min();

// Patches of at least this many true dofs are assigned to processes so as to balance their
// dense work (see the constructor of SliverPatchSolver).
constexpr std::size_t kBalanceSize = 256;

// Process owning a global true dof, from the first true dof of each process (and the
// global size at the end).
inline int OwnerRank(const std::vector<long long> &starts, long long g)
{
  return static_cast<int>(std::upper_bound(starts.begin(), starts.end(), g) -
                          starts.begin()) -
         1;
}

// Greedy merge of candidate patches (sorted true dofs) in order: a candidate merges with
// the patches it overlaps while the result has at most max_size true dofs, otherwise only
// its true dofs not yet in a patch form a new patch, so that the patches are disjoint. A
// merged patch belongs to the owner of its first candidate, and counts their slivers.
struct MergedPatches
{
  std::vector<std::vector<long long>> tdofs;
  std::vector<int> owner, slivers;
};

MergedPatches MergePatches(const std::vector<std::vector<long long>> &cand,
                           const std::vector<int> &cand_owner,
                           const std::vector<int> &cand_slivers, int max_size)
{
  std::vector<std::vector<long long>> groups;
  std::vector<int> owner, slivers;
  std::vector<char> live;
  std::unordered_map<long long, int> tdof_group;
  for (std::size_t c = 0; c < cand.size(); c++)
  {
    std::vector<int> overlap;
    for (long long g : cand[c])
    {
      if (auto it = tdof_group.find(g); it != tdof_group.end())
      {
        overlap.push_back(it->second);
      }
    }
    std::sort(overlap.begin(), overlap.end());
    overlap.erase(std::unique(overlap.begin(), overlap.end()), overlap.end());
    std::size_t merged_size = cand[c].size();
    for (int o : overlap)
    {
      merged_size += groups[o].size();
    }
    std::vector<long long> merged;
    int own = cand_owner[c], num = cand_slivers[c];
    if (overlap.empty() || merged_size <= static_cast<std::size_t>(max_size))
    {
      merged = cand[c];
      for (int o : overlap)
      {
        merged.insert(merged.end(), groups[o].begin(), groups[o].end());
        own = std::min(own, owner[o]);
        num += slivers[o];
        groups[o].clear();
        live[o] = 0;
      }
      std::sort(merged.begin(), merged.end());
      merged.erase(std::unique(merged.begin(), merged.end()), merged.end());
    }
    else
    {
      for (long long g : cand[c])
      {
        if (tdof_group.find(g) == tdof_group.end())
        {
          merged.push_back(g);
        }
      }
    }
    if (merged.empty() || merged.size() > static_cast<std::size_t>(max_size))
    {
      continue;
    }
    const int id = static_cast<int>(groups.size());
    for (long long g : merged)
    {
      tdof_group[g] = id;
    }
    groups.push_back(std::move(merged));
    owner.push_back(own);
    slivers.push_back(num);
    live.push_back(1);
  }
  MergedPatches result;
  for (std::size_t o = 0; o < groups.size(); o++)
  {
    if (live[o])
    {
      result.tdofs.push_back(std::move(groups[o]));
      result.owner.push_back(owner[o]);
      result.slivers.push_back(slivers[o]);
    }
  }
  return result;
}

}  // namespace

double MinDihedralAngle(const double *v0, const double *v1, const double *v2,
                        const double *v3)
{
  // Outward normals (not normalized) of the faces opposite to each vertex. The dihedral
  // angle at the edge shared by the faces opposite to vertices k and l is π minus the angle
  // between their outward normals.
  const double *v[4] = {v0, v1, v2, v3};
  double n[4][3];
  for (int k = 0; k < 4; k++)
  {
    const double *a = v[(k + 1) % 4], *b = v[(k + 2) % 4], *c = v[(k + 3) % 4];
    const double e1[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]},
                 e2[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]},
                 h[3] = {v[k][0] - a[0], v[k][1] - a[1], v[k][2] - a[2]};
    Cross(e1, e2, n[k]);
    if (Dot(n[k], h) > 0.0)
    {
      n[k][0] = -n[k][0];
      n[k][1] = -n[k][1];
      n[k][2] = -n[k][2];
    }
  }
  double min_angle = std::numbers::pi;
  for (int k = 0; k < 4; k++)
  {
    for (int l = k + 1; l < 4; l++)
    {
      double c[3];
      Cross(n[k], n[l], c);
      min_angle = std::min(min_angle, std::atan2(std::sqrt(Dot(c, c)), -Dot(n[k], n[l])));
    }
  }
  return min_angle * 180.0 / std::numbers::pi;
}

bool SliverPatchesSupported()
{
  return !mfem::Device::Allows(mfem::Backend::DEVICE_MASK | mfem::Backend::OMP_MASK);
}

std::vector<int> FindSliverElements(const mfem::ParMesh &mesh, double min_angle,
                                    SliverStatistics *stats)
{
  std::vector<int> slivers;
  long long counts[5] = {};
  double local_min = 180.0;
  if (mesh.Dimension() == 3 && mesh.SpaceDimension() == 3)
  {
    mfem::Array<int> verts;
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      if (mesh.GetElementGeometry(e) != mfem::Geometry::TETRAHEDRON)
      {
        continue;
      }
      mesh.GetElementVertices(e, verts);
      const double angle =
          MinDihedralAngle(mesh.GetVertex(verts[0]), mesh.GetVertex(verts[1]),
                           mesh.GetVertex(verts[2]), mesh.GetVertex(verts[3]));
      counts[0]++;
      counts[1] += (angle < 1.0);
      counts[2] += (angle < 2.0);
      counts[3] += (angle < 5.0);
      if (angle < min_angle)
      {
        slivers.push_back(e);
        counts[4]++;
      }
      local_min = std::min(local_min, angle);
    }
  }
  if (stats)
  {
    Mpi::GlobalSum(5, counts, mesh.GetComm());
    Mpi::GlobalMin(1, &local_min, mesh.GetComm());
    stats->num_tets = counts[0];
    stats->num_below_1 = counts[1];
    stats->num_below_2 = counts[2];
    stats->num_below_5 = counts[3];
    stats->num_slivers = counts[4];
    stats->min_angle = local_min;
  }
  return slivers;
}

template <OperatorType OperType>
struct SliverPatchSolver<OperType>::Factors
{
  // Symmetric diagonal scaling of the patch blocks (for each entry of loc), and the LU
  // factorizations of the scaled blocks, real or complex. Patches whose block is singular
  // are skipped.
  std::vector<double> scale;
  std::vector<Eigen::PartialPivLU<Eigen::MatrixXd>> lu_real;
  std::vector<Eigen::PartialPivLU<Eigen::MatrixXcd>> lu_complex;
  std::vector<char> valid;
  bool complex = false;
};

template <OperatorType OperType>
SliverPatchSolver<OperType>::SliverPatchSolver(const FiniteElementSpace &fespace,
                                               const std::vector<int> &slivers, int layers,
                                               const mfem::Array<int> *dbc_tdof_list,
                                               int max_size, bool update_residual)
  : update_residual(update_residual)
{
  ptr.assign(1, 0);
  // The conditions are the same on all processes, which take part in the collective
  // operations below together.
  if (fespace.IsBroken() || fespace.Dimension() != 3 || fespace.GetVDim() != 1 ||
      !SliverPatchesSupported())
  {
    return;
  }
  MPI_Comm space_comm = fespace.GetComm();
  const int rank = Mpi::Rank(space_comm), nranks = Mpi::Size(space_comm);
  const int nl = fespace.GetVSize(), nt = fespace.GetTrueVSize();
  const mfem::HypreParMatrix *P = fespace.Dof_TrueDof_Matrix();
  P->HostRead();
  Pd = std::make_unique<mfem::SparseMatrix>();
  Po = std::make_unique<mfem::SparseMatrix>();
  HYPRE_BigInt *cmap;
  P->GetDiag(*Pd);
  P->GetOffd(*Po, cmap);
  offd_cols.assign(cmap, cmap + Po->Width());
  my_start = static_cast<long long>(fespace.Get().GetMyTDofOffset());
  starts.resize(nranks + 1);
  Mpi::Allgather(1, &my_start, starts.data(), space_comm);
  starts[nranks] = static_cast<long long>(fespace.Get().GlobalTrueVSize());

  // Number of L-dofs of all processes with a nonzero prolongation entry for each owned true
  // dof (1 for the true dofs which only this process sees).
  mfem::Vector count(nt);
  {
    mfem::HypreParMatrix Q(*P);
    Q.HostReadWrite();
    hypre_ParCSRMatrix *q = Q;
    for (hypre_CSRMatrix *m : {hypre_ParCSRMatrixDiag(q), hypre_ParCSRMatrixOffd(q)})
    {
      const HYPRE_Int nnz =
          hypre_CSRMatrixI(m) ? hypre_CSRMatrixI(m)[hypre_CSRMatrixNumRows(m)] : 0;
      double *a = hypre_CSRMatrixData(m);
      for (HYPRE_Int k = 0; k < nnz; k++)
      {
        a[k] = (a[k] != 0.0) ? 1.0 : 0.0;
      }
    }
    mfem::Vector ones(nl);
    ones = 1.0;
    Q.MultTranspose(ones, count);
  }
  std::vector<char> ess(nt, 0);
  if (dbc_tdof_list)
  {
    for (int t : *dbc_tdof_list)
    {
      ess[t] = 1;
    }
  }
  if (update_residual)
  {
    essential = ess;
  }

  // The global true dof g and the weight w of an unconstrained L-dof (a single ±1
  // prolongation entry).
  auto Unconstrained = [&](int l, long long &g, double &w)
  {
    const int nd = Pd->RowSize(l), no = Po->RowSize(l);
    if (nd + no != 1)
    {
      return false;
    }
    if (nd)
    {
      g = my_start + Pd->GetRowColumns(l)[0];
      w = Pd->GetRowEntries(l)[0];
    }
    else
    {
      g = offd_cols[Po->GetRowColumns(l)[0]];
      w = Po->GetRowEntries(l)[0];
    }
    return std::abs(std::abs(w) - 1.0) < 1.0e-12;
  };

  // Candidate patch of each local sliver: the global true dofs of the unconstrained L-dofs
  // of its own entities (layers = 0), or of the local entities incident to its vertices.
  const auto &mesh = fespace.GetParMesh();
  const auto &fes = fespace.Get();
  std::vector<int> v2e_ptr, v2e;
  mfem::Array<int> verts, edges, faces, ori, dofs, ev;
  if (layers > 0 && !slivers.empty())
  {
    v2e_ptr.assign(mesh.GetNV() + 1, 0);
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      mesh.GetElementVertices(e, verts);
      for (int v : verts)
      {
        v2e_ptr[v + 1]++;
      }
    }
    std::partial_sum(v2e_ptr.begin(), v2e_ptr.end(), v2e_ptr.begin());
    v2e.resize(v2e_ptr.back());
    std::vector<int> pos(v2e_ptr.begin(), v2e_ptr.end() - 1);
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      mesh.GetElementVertices(e, verts);
      for (int v : verts)
      {
        v2e[pos[v]++] = e;
      }
    }
  }
  std::vector<long long> mine;  // For each sliver: number of true dofs, then the true dofs
  std::vector<int> star;
  std::vector<long long> cand;
  std::vector<int> excl[2];
  for (int sliver : slivers)
  {
    mesh.GetElementVertices(sliver, verts);
    auto HasSliverVertex = [&](const mfem::Array<int> &v)
    {
      for (int w : v)
      {
        for (int u : verts)
        {
          if (w == u)
          {
            return true;
          }
        }
      }
      return false;
    };
    star.assign(1, sliver);
    if (layers > 0)
    {
      for (int v : verts)
      {
        star.insert(star.end(), v2e.begin() + v2e_ptr[v], v2e.begin() + v2e_ptr[v + 1]);
      }
      std::sort(star.begin(), star.end());
      star.erase(std::unique(star.begin(), star.end()), star.end());
    }
    cand.clear();
    auto AddDofs = [&]()
    {
      for (int d : dofs)
      {
        const int l = (d >= 0) ? d : -1 - d;
        long long g;
        double w;
        if (!Unconstrained(l, g, w))
        {
          excl[0].push_back(l);
        }
        else if (g >= starts[rank] && g < starts[rank + 1] && ess[g - my_start])
        {
          excl[1].push_back(l);
        }
        else
        {
          cand.push_back(g);
        }
      }
    };
    for (int v : verts)
    {
      fes.GetVertexDofs(v, dofs);
      AddDofs();
    }
    for (int e : star)
    {
      fes.GetElementInteriorDofs(e, dofs);
      AddDofs();
      mesh.GetElementEdges(e, edges, ori);
      for (int edge : edges)
      {
        mesh.GetEdgeVertices(edge, ev);
        if (e == sliver || HasSliverVertex(ev))
        {
          fes.GetEdgeInteriorDofs(edge, dofs);
          AddDofs();
        }
      }
      mesh.GetElementFaces(e, faces, ori);
      for (int face : faces)
      {
        mesh.GetFaceVertices(face, ev);
        if (e == sliver || HasSliverVertex(ev))
        {
          fes.GetFaceInteriorDofs(face, dofs);
          AddDofs();
        }
      }
    }
    std::sort(cand.begin(), cand.end());
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
    if (!cand.empty())
    {
      mine.push_back(static_cast<long long>(cand.size()));
      mine.insert(mine.end(), cand.begin(), cand.end());
    }
  }
  for (int k = 0; k < 2; k++)
  {
    std::sort(excl[k].begin(), excl[k].end());
    num_excluded[k] = std::unique(excl[k].begin(), excl[k].end()) - excl[k].begin();
  }

  // Merge the candidate patches of this process: the small interior ones (all their true
  // dofs owned by this process and without L-dofs on other processes) stay local, and those
  // at the interfaces between processes or with at least kBalanceSize true dofs are
  // gathered from all processes and merged again (in the order of the processes), with the
  // same result on all processes. The gathered patches come first, numbered globally, then
  // the interior patches of this process.
  std::vector<std::vector<long long>> local_cand;
  for (std::size_t k = 0; k < mine.size();)
  {
    const auto n = static_cast<std::size_t>(mine[k]);
    local_cand.emplace_back(mine.begin() + k + 1, mine.begin() + k + 1 + n);
    k += n + 1;
  }
  std::vector<int> local_owner(local_cand.size(), rank),
      local_slivers(local_cand.size(), 1);
  auto local = MergePatches(local_cand, local_owner, local_slivers, max_size);
  const double *cnt = count.HostRead();
  std::vector<long long> iface;  // For each interface patch: size, slivers, true dofs
  std::vector<std::size_t> interior;
  for (std::size_t o = 0; o < local.tdofs.size(); o++)
  {
    bool inner = true;
    for (long long g : local.tdofs[o])
    {
      inner =
          inner && (g >= starts[rank] && g < starts[rank + 1] && cnt[g - my_start] == 1.0);
    }
    if (inner && local.tdofs[o].size() < kBalanceSize)
    {
      interior.push_back(o);
    }
    else
    {
      iface.push_back(static_cast<long long>(local.tdofs[o].size()));
      iface.push_back(local.slivers[o]);
      iface.insert(iface.end(), local.tdofs[o].begin(), local.tdofs[o].end());
    }
  }
  long long num_any = static_cast<long long>(local.tdofs.size());
  Mpi::GlobalSum(1, &num_any, space_comm);
  if (num_any == 0)
  {
    return;
  }
  int my_count = static_cast<int>(iface.size());
  std::vector<int> counts(nranks), displs(nranks + 1, 0);
  Mpi::Allgather(1, &my_count, counts.data(), space_comm);
  std::partial_sum(counts.begin(), counts.end(), displs.begin() + 1);
  std::vector<long long> all(displs.back());
  Mpi::Allgatherv(my_count, iface.data(), all.data(), counts.data(), displs.data(),
                  space_comm);

  // Leave out the essential true dofs of other processes (their owners flag them).
  std::vector<int> flag(all.size(), 0);
  for (int q = 0; q < nranks; q++)
  {
    for (int k = displs[q]; k < displs[q + 1];)
    {
      const int n = static_cast<int>(all[k]);
      for (int i = k + 2; i < k + 2 + n; i++)
      {
        if (all[i] >= starts[rank] && all[i] < starts[rank + 1] && ess[all[i] - my_start])
        {
          flag[i] = 1;
        }
      }
      k += n + 2;
    }
  }
  if (!flag.empty())
  {
    Mpi::GlobalMax(static_cast<int>(flag.size()), flag.data(), space_comm);
  }
  std::vector<std::vector<long long>> iface_cand;
  std::vector<int> iface_owner, iface_slivers;
  for (int q = 0; q < nranks; q++)
  {
    for (int k = displs[q]; k < displs[q + 1];)
    {
      const int n = static_cast<int>(all[k]);
      std::vector<long long> c;
      for (int i = k + 2; i < k + 2 + n; i++)
      {
        if (!flag[i])
        {
          c.push_back(all[i]);
        }
      }
      if (!c.empty())
      {
        iface_cand.push_back(std::move(c));
        iface_owner.push_back(q);
        iface_slivers.push_back(static_cast<int>(all[k + 1]));
      }
      k += n + 2;
    }
  }
  auto merged = MergePatches(iface_cand, iface_owner, iface_slivers, max_size);

  // The dense work of an application, n² for a patch of n true dofs, is that of its owner:
  // each large gathered patch, from the largest, stays with the process owning most of its
  // true dofs unless that raises the work of the process above the largest n² or the mean
  // work (the best possible maximum), and goes to the process with the least work so far
  // otherwise, so that clusters of slivers in a few processes do not set the cost of the
  // patch solves, without moving the patches away from their true dofs where not needed.
  // The assignment is the same on all processes.
  {
    double my_work = 0.0;
    for (std::size_t o : interior)
    {
      const auto n = static_cast<double>(local.tdofs[o].size());
      my_work += n * n;
    }
    std::vector<double> work(nranks);
    Mpi::Allgather(1, &my_work, work.data(), space_comm);
    std::vector<std::size_t> large;
    for (std::size_t o = 0; o < merged.tdofs.size(); o++)
    {
      const auto n = static_cast<double>(merged.tdofs[o].size());
      if (merged.tdofs[o].size() >= kBalanceSize)
      {
        large.push_back(o);
      }
      else
      {
        work[merged.owner[o]] += n * n;
      }
    }
    std::stable_sort(large.begin(), large.end(), [&](std::size_t a, std::size_t b)
                     { return merged.tdofs[a].size() > merged.tdofs[b].size(); });
    double target = std::accumulate(work.begin(), work.end(), 0.0);
    for (std::size_t o : large)
    {
      const auto n = static_cast<double>(merged.tdofs[o].size());
      target += n * n;
    }
    target /= nranks;
    if (!large.empty())
    {
      const auto n = static_cast<double>(merged.tdofs[large.front()].size());
      target = std::max(target, n * n);
    }
    using Load = std::pair<double, int>;
    std::priority_queue<Load, std::vector<Load>, std::greater<Load>> least;
    for (int q = 0; q < nranks; q++)
    {
      least.emplace(work[q], q);
    }
    for (std::size_t o : large)
    {
      const auto n = static_cast<double>(merged.tdofs[o].size());
      std::map<int, int> count;
      for (long long g : merged.tdofs[o])
      {
        count[OwnerRank(starts, g)]++;
      }
      int home = count.begin()->first;
      for (const auto &[q, c] : count)
      {
        home = (c > count[home]) ? q : home;
      }
      int q = home;
      if (work[home] + n * n > target)
      {
        // The least loaded process (entries of the queue can be stale: loads only grow).
        while (least.top().first != work[least.top().second])
        {
          const int r = least.top().second;
          least.pop();
          least.emplace(work[r], r);
        }
        q = least.top().second;
        least.pop();
      }
      merged.owner[o] = q;
      work[q] += n * n;
      least.emplace(work[q], q);
    }
  }
  std::vector<std::vector<long long>> groups = std::move(merged.tdofs);
  std::vector<int> group_owner = std::move(merged.owner),
                   group_slivers = std::move(merged.slivers);
  for (std::size_t o : interior)
  {
    groups.push_back(std::move(local.tdofs[o]));
    group_owner.push_back(rank);
    group_slivers.push_back(local.slivers[o]);
  }
  std::vector<char> group_live(groups.size(), 1);

  // The final patches, the patch and position of each patch true dof, and the patch
  // entries of the L-dofs of this process.
  std::unordered_map<long long, std::pair<int, int>> tdof_patch;
  for (std::size_t o = 0; o < groups.size(); o++)
  {
    if (!group_live[o])
    {
      continue;
    }
    const int id = static_cast<int>(patch_owner.size());
    for (std::size_t i = 0; i < groups[o].size(); i++)
    {
      // Only the true dofs of the L-dofs of this process (owned, or of other processes).
      const long long g = groups[o][i];
      if ((g >= starts[rank] && g < starts[rank + 1]) ||
          std::binary_search(offd_cols.begin(), offd_cols.end(), g))
      {
        tdof_patch[g] = {id, static_cast<int>(i)};
      }
    }
    patch_owner.push_back(group_owner[o]);
    patch_size.push_back(static_cast<int>(groups[o].size()));
    if (group_owner[o] == rank)
    {
      num_slivers += group_slivers[o];
    }
  }
  const int num_patches = static_cast<int>(patch_owner.size());
  ldof_ent_ptr.assign(nl + 1, 0);
  for (int pass = 0; pass < 2; pass++)
  {
    for (int l = 0; l < nl; l++)
    {
      auto AddEntry = [&](long long g, double w)
      {
        if (auto it = tdof_patch.find(g); it != tdof_patch.end())
        {
          if (pass == 0)
          {
            ldof_ent_ptr[l + 1]++;
          }
          else
          {
            const int k = ldof_ent_ptr[l]++;
            ldof_ent_patch[k] = it->second.first;
            ldof_ent_pos[k] = it->second.second;
            ldof_ent_w[k] = w;
          }
        }
      };
      const int *cd = Pd->GetRowColumns(l);
      const double *wd = Pd->GetRowEntries(l);
      for (int c = 0; c < Pd->RowSize(l); c++)
      {
        AddEntry(my_start + cd[c], wd[c]);
      }
      const int *co = Po->GetRowColumns(l);
      const double *wo = Po->GetRowEntries(l);
      for (int c = 0; c < Po->RowSize(l); c++)
      {
        AddEntry(offd_cols[co[c]], wo[c]);
      }
    }
    if (pass == 0)
    {
      std::partial_sum(ldof_ent_ptr.begin(), ldof_ent_ptr.end(), ldof_ent_ptr.begin());
      ldof_ent_patch.resize(ldof_ent_ptr.back());
      ldof_ent_pos.resize(ldof_ent_ptr.back());
      ldof_ent_w.resize(ldof_ent_ptr.back());
    }
    else
    {
      // Restore the offsets shifted by the second pass.
      for (int l = nl; l > 0; l--)
      {
        ldof_ent_ptr[l] = ldof_ent_ptr[l - 1];
      }
      ldof_ent_ptr[0] = 0;
    }
  }
  if (ldof_ent_ptr.back() > 0)
  {
    marker.assign(nl, 0);
    for (int l = 0; l < nl; l++)
    {
      marker[l] = (ldof_ent_ptr[l + 1] > ldof_ent_ptr[l]) ? 2 : 0;
    }
    if (update_residual)
    {
      // The other L-dofs of the elements with patch L-dofs are coupled to the patches.
      for (int e = 0; e < mesh.GetNE(); e++)
      {
        fes.GetElementDofs(e, dofs);
        bool touches = false;
        for (int d : dofs)
        {
          touches = touches || (marker[(d >= 0) ? d : -1 - d] == 2);
        }
        if (touches)
        {
          for (int d : dofs)
          {
            char &m = marker[(d >= 0) ? d : -1 - d];
            m = std::max(m, char(1));
          }
        }
      }
    }
  }
  else
  {
    ldof_ent_ptr.clear();
  }

  // Patches owned by this process: local true dofs, and slots for the residual values of
  // the true dofs of other processes, ordered by process, patch and position. The other
  // processes send the values of their true dofs of the patches of this process in the
  // same order.
  local_patch.assign(num_patches, -1);
  std::vector<std::tuple<int, int, int>> remote;       // (process, patch, position in loc)
  std::vector<std::tuple<int, int, int, int>> helper;  // (owner, patch, position, tdof)
  for (std::size_t o = 0, id = 0; o < groups.size(); o++)
  {
    if (!group_live[o])
    {
      continue;
    }
    const int k = static_cast<int>(id++);
    if (patch_owner[k] == rank)
    {
      local_patch[k] = NumPatches();
      for (std::size_t i = 0; i < groups[o].size(); i++)
      {
        const long long g = groups[o][i];
        gtdofs.push_back(g);
        if (g >= starts[rank] && g < starts[rank + 1])
        {
          loc.push_back(static_cast<int>(g - my_start));
        }
        else
        {
          remote.emplace_back(OwnerRank(starts, g), k, static_cast<int>(loc.size()));
          loc.push_back(-1);
        }
      }
      ptr.push_back(static_cast<int>(loc.size()));
    }
    else
    {
      for (std::size_t i = 0; i < groups[o].size(); i++)
      {
        const long long g = groups[o][i];
        if (g >= starts[rank] && g < starts[rank + 1])
        {
          helper.emplace_back(patch_owner[k], k, static_cast<int>(i),
                              static_cast<int>(g - my_start));
        }
      }
    }
  }
  num_remote = static_cast<long long>(remote.size());
  std::sort(remote.begin(), remote.end());
  recv_ptr.assign(1, 0);
  for (std::size_t s = 0; s < remote.size(); s++)
  {
    const auto [q, k, i] = remote[s];
    if (s == 0 || q != std::get<0>(remote[s - 1]))
    {
      recv_ranks.push_back(q);
      recv_ptr.push_back(recv_ptr.back());
    }
    loc[i] = -1 - static_cast<int>(s);
    recv_ptr.back()++;
  }
  std::sort(helper.begin(), helper.end());
  send_ptr.assign(1, 0);
  for (std::size_t s = 0; s < helper.size(); s++)
  {
    const auto [o, k, i, t] = helper[s];
    if (s == 0 || o != std::get<0>(helper[s - 1]))
    {
      send_ranks.push_back(o);
      send_ptr.push_back(send_ptr.back());
    }
    send_idx.push_back(t);
    send_ptr.back()++;
  }
  MPI_Comm_dup(space_comm, &comm);
}

template <OperatorType OperType>
SliverPatchSolver<OperType>::~SliverPatchSolver()
{
  if (comm != MPI_COMM_NULL)
  {
    MPI_Comm_free(&comm);
  }
}

template <OperatorType OperType>
int SliverPatchSolver<OperType>::MaxPatchSize() const
{
  int n = 0;
  for (int p = 0; p < NumPatches(); p++)
  {
    n = std::max(n, ptr[p + 1] - ptr[p]);
  }
  return n;
}

template <OperatorType OperType>
double SliverPatchSolver<OperType>::PatchWork() const
{
  double work = 0.0;
  for (int p = 0; p < NumPatches(); p++)
  {
    const double n = ptr[p + 1] - ptr[p];
    work += n * n;
  }
  return work;
}

template <OperatorType OperType>
bool SliverPatchSolver<OperType>::SetOperator(const OperType &op)
{
  factors.reset();
  active = false;
  if (comm == MPI_COMM_NULL)
  {
    return true;
  }

  // Local real and imaginary parts of the operator.
  const Operator *parts[2] = {nullptr, nullptr};
  int supported = 1;
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    const auto *par_op = dynamic_cast<const ComplexParOperator *>(&op);
    supported = (par_op != nullptr);
    if (par_op)
    {
      parts[0] = par_op->LocalOperator().Real();
      parts[1] = par_op->LocalOperator().Imag();
    }
  }
  else
  {
    const auto *par_op = dynamic_cast<const ParOperator *>(&op);
    supported = (par_op != nullptr);
    if (par_op)
    {
      parts[0] = &par_op->LocalOperator();
    }
  }
  for (const auto *part : parts)
  {
    if (part && !dynamic_cast<const ceed::Operator *>(part) &&
        !dynamic_cast<const hypre::HypreCSRMatrix *>(part))
    {
      supported = 0;
    }
  }

  // The contributions of the element matrices of this process to the patch blocks, and
  // with residual updates to the couplings of the patch true dofs (columns) to other true
  // dofs (rows): (patch, global row, column) -> value, sent to the owners of the patches.
  // Without residual updates, only the entries between the L-dofs of the same patch.
  std::map<std::tuple<int, long long, int>, std::array<double, 2>> entries;
  for (int k = 0; k < 2; k++)
  {
    if (!parts[k] || ldof_ent_ptr.empty() || !supported)
    {
      continue;
    }
    auto AddEntry = [&](int i, int j, double a)
    {
      for (int ej = ldof_ent_ptr[j]; ej < ldof_ent_ptr[j + 1]; ej++)
      {
        const int q = ldof_ent_patch[ej], pj = ldof_ent_pos[ej];
        const double wa = ldof_ent_w[ej] * a;
        if (!update_residual)
        {
          for (int ei = ldof_ent_ptr[i]; ei < ldof_ent_ptr[i + 1]; ei++)
          {
            if (ldof_ent_patch[ei] == q)
            {
              entries[{q, -1 - ldof_ent_pos[ei], pj}][k] += ldof_ent_w[ei] * wa;
            }
          }
          continue;
        }
        // All true dofs of row i, through the prolongation.
        const int *cd = Pd->GetRowColumns(i);
        const double *wd = Pd->GetRowEntries(i);
        for (int c = 0; c < Pd->RowSize(i); c++)
        {
          entries[{q, my_start + cd[c], pj}][k] += wd[c] * wa;
        }
        const int *co = Po->GetRowColumns(i);
        const double *wo = Po->GetRowEntries(i);
        for (int c = 0; c < Po->RowSize(i); c++)
        {
          entries[{q, offd_cols[co[c]], pj}][k] += wo[c] * wa;
        }
      }
    };
    if (const auto *ceed_op = dynamic_cast<const ceed::Operator *>(parts[k]))
    {
      if (!ceed::AssembleMarkedEntries(*ceed_op, marker, AddEntry))
      {
        supported = 0;
      }
    }
    else if (const auto *csr = dynamic_cast<const hypre::HypreCSRMatrix *>(parts[k]))
    {
      hypre_CSRMatrix *mat = *csr, *host = nullptr;
      if (hypre_CSRMatrixMemoryLocation(mat) != HYPRE_MEMORY_HOST)
      {
        host = hypre_CSRMatrixClone_v2(mat, 1, HYPRE_MEMORY_HOST);
        mat = host;
      }
      const HYPRE_Int *I = hypre_CSRMatrixI(mat), *J = hypre_CSRMatrixJ(mat);
      const double *a = hypre_CSRMatrixData(mat);
      for (int l = 0; l < static_cast<int>(marker.size()); l++)
      {
        if (!marker[l])
        {
          continue;
        }
        for (HYPRE_Int c = I[l]; c < I[l + 1]; c++)
        {
          if (marker[J[c]] == 2)
          {
            AddEntry(l, J[c], a[c]);
          }
        }
      }
      if (host)
      {
        hypre_CSRMatrixDestroy(host);
      }
    }
  }
  Mpi::GlobalMin(1, &supported, comm);
  if (!supported)
  {
    return false;
  }

  // Send the entries to the owners of the patches as (patch, row, column, real, imaginary)
  // records: the row is a global true dof, or -1 - (position in the patch) without
  // residual updates.
  const int nranks = Mpi::Size(comm), rank = Mpi::Rank(comm);
  std::vector<int> scount(nranks, 0), rcount(nranks), sdispl(nranks + 1, 0),
      rdispl(nranks + 1, 0);
  for (const auto &[key, val] : entries)
  {
    scount[patch_owner[std::get<0>(key)]] += 5;
  }
  MPI_Alltoall(scount.data(), 1, MPI_INT, rcount.data(), 1, MPI_INT, comm);
  std::partial_sum(scount.begin(), scount.end(), sdispl.begin() + 1);
  std::partial_sum(rcount.begin(), rcount.end(), rdispl.begin() + 1);
  std::vector<double> sendv(sdispl.back()), recvv(rdispl.back());
  {
    std::vector<int> pos(sdispl.begin(), sdispl.end() - 1);
    for (const auto &[key, val] : entries)
    {
      const int o = patch_owner[std::get<0>(key)];
      double *rec = sendv.data() + pos[o];
      rec[0] = std::get<0>(key);
      rec[1] = static_cast<double>(std::get<1>(key));
      rec[2] = std::get<2>(key);
      rec[3] = val[0];
      rec[4] = val[1];
      pos[o] += 5;
    }
  }
  MPI_Alltoallv(sendv.data(), scount.data(), sdispl.data(), MPI_DOUBLE, recvv.data(),
                rcount.data(), rdispl.data(), MPI_DOUBLE, comm);
  entries.clear();

  // Blocks of the owned patches (column-major), and with residual updates the rows of the
  // owned patches outside the patch: (patch, global row) -> entries (column, value).
  block_off.assign(NumPatches() + 1, 0);
  for (int p = 0; p < NumPatches(); p++)
  {
    const std::size_t n = ptr[p + 1] - ptr[p];
    block_off[p + 1] = block_off[p] + n * n;
  }
  for (auto &b : block_val)
  {
    b.assign(block_off.back(), 0.0);
  }
  std::unordered_map<long long, std::pair<int, int>> owned_pos;  // Patch true dofs
  if (update_residual)
  {
    for (int p = 0; p < NumPatches(); p++)
    {
      for (int k = ptr[p]; k < ptr[p + 1]; k++)
      {
        owned_pos[gtdofs[k]] = {p, k - ptr[p]};
      }
    }
  }
  std::map<std::pair<int, long long>, std::vector<std::array<double, 3>>> coupling;
  for (std::size_t r = 0; r + 4 < recvv.size(); r += 5)
  {
    const int p = local_patch[static_cast<int>(recvv[r])];
    MFEM_ASSERT(p >= 0, "Received a sliver patch entry for another process!");
    const long long row = static_cast<long long>(recvv[r + 1]);
    const int n = ptr[p + 1] - ptr[p], pj = static_cast<int>(recvv[r + 2]);
    int pi = -1;
    if (row < 0)
    {
      pi = static_cast<int>(-1 - row);
    }
    else if (auto it = owned_pos.find(row); it != owned_pos.end() && it->second.first == p)
    {
      pi = it->second.second;
    }
    if (pi >= 0)
    {
      block_val[0][block_off[p] + std::size_t(pj) * n + pi] += recvv[r + 3];
      block_val[1][block_off[p] + std::size_t(pj) * n + pi] += recvv[r + 4];
    }
    else
    {
      coupling[{p, row}].push_back({double(pj), recvv[r + 3], recvv[r + 4]});
    }
  }

  // Rows of the residual updates of the owned patches (the patch rows, then the coupled
  // ones by global true dof), with the entries of the coupled rows summed by column, and
  // the update messages: each process receives the global true dofs of its rows from the
  // owners of the patches.
  up_ptr.assign(1, 0);
  up_rows.clear();
  cp_ptr.assign(1, 0);
  cp_col.clear();
  cp_val[0].clear();
  cp_val[1].clear();
  up_send_ranks.clear();
  up_send_ptr.assign(1, 0);
  up_recv_ranks.clear();
  up_recv_ptr.assign(1, 0);
  up_recv_idx.clear();
  if (update_residual)
  {
    std::vector<std::tuple<int, int, long long>>
        remote_rows;  // (rank, index in up_rows, row)
    auto AddRow = [&](long long g)
    {
      if (g >= starts[rank] && g < starts[rank + 1])
      {
        // Essential rows are decoupled from the patches (eliminated operator).
        const int t = static_cast<int>(g - my_start);
        up_rows.push_back(essential[t] ? kSkip : t);
      }
      else
      {
        remote_rows.emplace_back(OwnerRank(starts, g), static_cast<int>(up_rows.size()), g);
        up_rows.push_back(-1);
      }
    };
    auto it = coupling.begin();
    for (int p = 0; p < NumPatches(); p++)
    {
      for (int k = ptr[p]; k < ptr[p + 1]; k++)
      {
        AddRow(gtdofs[k]);
      }
      for (; it != coupling.end() && it->first.first == p; ++it)
      {
        AddRow(it->first.second);
        std::map<int, std::array<double, 2>> sum;
        for (const auto &e : it->second)
        {
          auto &v = sum[static_cast<int>(e[0])];
          v[0] += e[1];
          v[1] += e[2];
        }
        for (const auto &[c, v] : sum)
        {
          cp_col.push_back(c);
          cp_val[0].push_back(v[0]);
          cp_val[1].push_back(v[1]);
        }
        cp_ptr.push_back(static_cast<int>(cp_col.size()));
      }
      up_ptr.push_back(static_cast<int>(up_rows.size()));
    }
    std::sort(remote_rows.begin(), remote_rows.end());
    std::vector<int> ucount(nranks, 0), urcount(nranks), udispl(nranks + 1, 0),
        urdispl(nranks + 1, 0);
    std::vector<long long> ug;
    for (std::size_t k = 0; k < remote_rows.size(); k++)
    {
      const auto [q, i, g] = remote_rows[k];
      if (k == 0 || q != std::get<0>(remote_rows[k - 1]))
      {
        up_send_ranks.push_back(q);
        up_send_ptr.push_back(up_send_ptr.back());
      }
      up_rows[i] = -1 - static_cast<int>(k);
      up_send_ptr.back()++;
      ucount[q]++;
      ug.push_back(g);
    }
    MPI_Alltoall(ucount.data(), 1, MPI_INT, urcount.data(), 1, MPI_INT, comm);
    std::partial_sum(ucount.begin(), ucount.end(), udispl.begin() + 1);
    std::partial_sum(urcount.begin(), urcount.end(), urdispl.begin() + 1);
    std::vector<long long> urg(urdispl.back());
    MPI_Alltoallv(ug.data(), ucount.data(), udispl.data(), MPI_LONG_LONG, urg.data(),
                  urcount.data(), urdispl.data(), MPI_LONG_LONG, comm);
    for (int q = 0; q < nranks; q++)
    {
      if (urcount[q] > 0)
      {
        up_recv_ranks.push_back(q);
        up_recv_ptr.push_back(urdispl[q + 1]);
      }
    }
    for (long long g : urg)
    {
      MFEM_ASSERT(g >= starts[rank] && g < starts[rank + 1],
                  "Received a residual update row of another process!");
      const int t = static_cast<int>(g - my_start);
      up_recv_idx.push_back(essential[t] ? -1 : t);
    }
  }
  std::vector<double>(&blocks)[2] = block_val;
  const auto &offsets = block_off;

  // Factor the symmetrically scaled blocks, in real arithmetic when the imaginary part
  // vanishes on all patches of this process.
  auto f = std::make_unique<Factors>();
  f->complex =
      std::any_of(blocks[1].begin(), blocks[1].end(), [](double a) { return a != 0.0; });
  f->scale.resize(loc.size());
  f->valid.assign(NumPatches(), 0);
  if (f->complex)
  {
    f->lu_complex.resize(NumPatches());
  }
  else
  {
    f->lu_real.resize(NumPatches());
  }
  int all_valid = 1;
  for (int p = 0; p < NumPatches(); p++)
  {
    const int n = ptr[p + 1] - ptr[p];
    Eigen::Map<const Eigen::MatrixXd> Ar(blocks[0].data() + offsets[p], n, n),
        Ai(blocks[1].data() + offsets[p], n, n);
    double *s = f->scale.data() + ptr[p];
    for (int i = 0; i < n; i++)
    {
      const double d = std::hypot(Ar(i, i), Ai(i, i));
      s[i] = (d > 0.0 && std::isfinite(d)) ? 1.0 / std::sqrt(d) : 1.0;
    }
    Eigen::Map<const Eigen::VectorXd> S(s, n);
    double rcond;
    if (f->complex)
    {
      Eigen::MatrixXcd A(n, n);
      A.real() = S.asDiagonal() * Ar * S.asDiagonal();
      A.imag() = S.asDiagonal() * Ai * S.asDiagonal();
      f->lu_complex[p].compute(A);
      rcond = f->lu_complex[p].rcond();
    }
    else
    {
      f->lu_real[p].compute(S.asDiagonal() * Ar * S.asDiagonal());
      rcond = f->lu_real[p].rcond();
    }
    f->valid[p] = (std::isfinite(rcond) && rcond > 1.0e-14);
    all_valid = all_valid && f->valid[p];
  }
  factors = std::move(f);
  active = (NumPatches() > 0 || !send_ranks.empty() || !up_recv_ranks.empty());
  Mpi::GlobalMin(1, &all_valid, comm);
  return all_valid;
}

template <OperatorType OperType>
void SliverPatchSolver<OperType>::AddMultUpdate(VecType &r, VecType &y) const
{
  MFEM_VERIFY(
      update_residual,
      "SliverPatchSolver residual updates require construction with update_residual!");
  Apply(r, y, &r);
}

template <OperatorType OperType>
void SliverPatchSolver<OperType>::Apply(const VecType &r, VecType &y,
                                        VecType *r_update) const
{
  if (!active)
  {
    return;
  }
  const auto &f = *factors;
  constexpr bool cplx = std::is_same_v<OperType, ComplexOperator>;
  constexpr int nv = cplx ? 2 : 1;
  const double *R[2] = {nullptr, nullptr};
  double *Y[2] = {nullptr, nullptr};
  if constexpr (cplx)
  {
    R[0] = r.Real().HostRead();
    R[1] = r.Imag().HostRead();
    Y[0] = y.Real().HostReadWrite();
    Y[1] = y.Imag().HostReadWrite();
  }
  else
  {
    R[0] = r.HostRead();
    Y[0] = y.HostReadWrite();
  }
  const bool update = (r_update != nullptr);

  // Buffers: residuals received for the owned patches, then the corrections received for
  // the local true dofs of other patches; the residual values sent to the patch owners,
  // then the corrections returned for the received residuals. With residual updates, the
  // updates sent and received in their own buffers.
  const int nrecv = recv_ptr.back(), nsend = send_ptr.back();
  const int nup_s = update ? up_send_ptr.back() : 0,
            nup_r = update ? up_recv_ptr.back() : 0;
  rbuf.resize(std::size_t(nrecv + nsend) * nv);
  sbuf.assign(std::size_t(nsend + nrecv) * nv, 0.0);
  ubuf_s.assign(std::size_t(nup_s) * nv, 0.0);
  ubuf_r.resize(std::size_t(nup_r) * nv);
  reqs.resize(2 * (recv_ranks.size() + send_ranks.size()) +
              (update ? (up_send_ranks.size() + up_recv_ranks.size()) : 0));
  std::size_t nreq = 0;
  constexpr int tag_res = 1, tag_cor = 2, tag_up = 3;
  for (std::size_t i = 0; i < recv_ranks.size(); i++)
  {
    MPI_Irecv(rbuf.data() + std::size_t(recv_ptr[i]) * nv,
              (recv_ptr[i + 1] - recv_ptr[i]) * nv, MPI_DOUBLE, recv_ranks[i], tag_res,
              comm, &reqs[nreq++]);
  }
  const std::size_t nreq_res = nreq;
  for (std::size_t i = 0; i < send_ranks.size(); i++)
  {
    MPI_Irecv(rbuf.data() + std::size_t(nrecv + send_ptr[i]) * nv,
              (send_ptr[i + 1] - send_ptr[i]) * nv, MPI_DOUBLE, send_ranks[i], tag_cor,
              comm, &reqs[nreq++]);
  }
  if (update)
  {
    for (std::size_t i = 0; i < up_recv_ranks.size(); i++)
    {
      MPI_Irecv(ubuf_r.data() + std::size_t(up_recv_ptr[i]) * nv,
                (up_recv_ptr[i + 1] - up_recv_ptr[i]) * nv, MPI_DOUBLE, up_recv_ranks[i],
                tag_up, comm, &reqs[nreq++]);
    }
  }
  const std::size_t nreq_cor = nreq;
  for (std::size_t i = 0; i < send_ranks.size(); i++)
  {
    for (int k = send_ptr[i]; k < send_ptr[i + 1]; k++)
    {
      for (int v = 0; v < nv; v++)
      {
        sbuf[std::size_t(k) * nv + v] = R[v][send_idx[k]];
      }
    }
    MPI_Isend(sbuf.data() + std::size_t(send_ptr[i]) * nv,
              (send_ptr[i + 1] - send_ptr[i]) * nv, MPI_DOUBLE, send_ranks[i], tag_res,
              comm, &reqs[nreq++]);
  }
  MPI_Waitall(static_cast<int>(nreq_res), reqs.data(), MPI_STATUSES_IGNORE);

  // Solve on the owned patches.
  auto Residual = [&](int k, int v)
  { return (loc[k] >= 0) ? R[v][loc[k]] : rbuf[std::size_t(-1 - loc[k]) * nv + v]; };
  delta.assign(loc.size() * nv, 0.0);
  for (int p = 0; p < NumPatches(); p++)
  {
    if (!f.valid[p])
    {
      continue;
    }
    const int n = ptr[p + 1] - ptr[p];
    const double *s = f.scale.data() + ptr[p];
    double *d = delta.data() + std::size_t(ptr[p]) * nv;
    if (f.complex)
    {
      // Only complex operators have complex blocks.
      Eigen::VectorXcd b(n);
      for (int i = 0; i < n; i++)
      {
        b(i) = std::complex<double>(s[i] * Residual(ptr[p] + i, 0),
                                    s[i] * Residual(ptr[p] + i, nv - 1));
      }
      const Eigen::VectorXcd x = f.lu_complex[p].solve(b);
      for (int i = 0; i < n; i++)
      {
        d[i * nv] = s[i] * x(i).real();
        if constexpr (cplx)
        {
          d[i * nv + 1] = s[i] * x(i).imag();
        }
      }
    }
    else
    {
      Eigen::MatrixXd b(n, nv);
      for (int i = 0; i < n; i++)
      {
        for (int v = 0; v < nv; v++)
        {
          b(i, v) = s[i] * Residual(ptr[p] + i, v);
        }
      }
      const Eigen::MatrixXd x = f.lu_real[p].solve(b);
      for (int i = 0; i < n; i++)
      {
        for (int v = 0; v < nv; v++)
        {
          d[i * nv + v] = s[i] * x(i, v);
        }
      }
    }
  }
  for (std::size_t k = 0; k < loc.size(); k++)
  {
    for (int v = 0; v < nv; v++)
    {
      if (loc[k] >= 0)
      {
        Y[v][loc[k]] += delta[k * nv + v];
      }
      else
      {
        sbuf[std::size_t(nsend - 1 - loc[k]) * nv + v] = delta[k * nv + v];
      }
    }
  }

  // Residual updates u = A[rows, patch] δ: the patch rows with the block, the others with
  // their couplings, for the local rows and the update messages.
  double *RU[2] = {nullptr, nullptr};
  if (update)
  {
    if constexpr (cplx)
    {
      RU[0] = r_update->Real().HostReadWrite();
      RU[1] = r_update->Imag().HostReadWrite();
    }
    else
    {
      RU[0] = r_update->HostReadWrite();
    }
    int c = 0;  // Coupled row index (CSR)
    for (int p = 0; p < NumPatches(); p++)
    {
      const int n = ptr[p + 1] - ptr[p];
      const double *d = delta.data() + std::size_t(ptr[p]) * nv;
      const double *br = block_val[0].data() + block_off[p],
                   *bi = block_val[1].data() + block_off[p];
      for (int k = up_ptr[p]; k < up_ptr[p + 1]; k++)
      {
        double u[2] = {0.0, 0.0};
        auto Add = [&](int j, double ar, double ai)
        {
          if constexpr (cplx)
          {
            u[0] += ar * d[j * nv] - ai * d[j * nv + 1];
            u[1] += ar * d[j * nv + 1] + ai * d[j * nv];
          }
          else
          {
            u[0] += ar * d[j * nv];
          }
        };
        const int i = k - up_ptr[p];
        if (i < n)
        {
          for (int j = 0; j < n; j++)
          {
            Add(j, br[std::size_t(j) * n + i], bi[std::size_t(j) * n + i]);
          }
        }
        else
        {
          for (int e = cp_ptr[c]; e < cp_ptr[c + 1]; e++)
          {
            Add(cp_col[e], cp_val[0][e], cp_val[1][e]);
          }
          c++;
        }
        for (int v = 0; v < nv; v++)
        {
          if (up_rows[k] >= 0)
          {
            RU[v][up_rows[k]] -= u[v];
          }
          else if (up_rows[k] != kSkip)
          {
            ubuf_s[std::size_t(-1 - up_rows[k]) * nv + v] = u[v];
          }
        }
      }
    }
  }

  // Return the corrections (and the updates) of the true dofs of other processes, and add
  // those received.
  for (std::size_t i = 0; i < recv_ranks.size(); i++)
  {
    MPI_Isend(sbuf.data() + std::size_t(nsend + recv_ptr[i]) * nv,
              (recv_ptr[i + 1] - recv_ptr[i]) * nv, MPI_DOUBLE, recv_ranks[i], tag_cor,
              comm, &reqs[nreq++]);
  }
  if (update)
  {
    for (std::size_t i = 0; i < up_send_ranks.size(); i++)
    {
      MPI_Isend(ubuf_s.data() + std::size_t(up_send_ptr[i]) * nv,
                (up_send_ptr[i + 1] - up_send_ptr[i]) * nv, MPI_DOUBLE, up_send_ranks[i],
                tag_up, comm, &reqs[nreq++]);
    }
  }
  MPI_Waitall(static_cast<int>(nreq_cor - nreq_res), reqs.data() + nreq_res,
              MPI_STATUSES_IGNORE);
  for (int k = 0; k < nsend; k++)
  {
    for (int v = 0; v < nv; v++)
    {
      Y[v][send_idx[k]] += rbuf[std::size_t(nrecv + k) * nv + v];
    }
  }
  if (update)
  {
    for (int k = 0; k < nup_r; k++)
    {
      if (up_recv_idx[k] < 0)
      {
        continue;
      }
      for (int v = 0; v < nv; v++)
      {
        RU[v][up_recv_idx[k]] -= ubuf_r[std::size_t(k) * nv + v];
      }
    }
  }
  MPI_Waitall(static_cast<int>(nreq - nreq_cor), reqs.data() + nreq_cor,
              MPI_STATUSES_IGNORE);
}

template class SliverPatchSolver<Operator>;
template class SliverPatchSolver<ComplexOperator>;

}  // namespace palace
