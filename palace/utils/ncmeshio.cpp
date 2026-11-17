// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ncmeshio.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#include <mfem/mesh/ncmesh_tables.hpp>
#if defined(MFEM_USE_ZLIB)
#include <zlib.h>
#endif
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/outputdir.hpp"

namespace palace::mesh
{

namespace
{

// Read access to the refinement tree of an NCMesh, which MFEM keeps protected: pointers to
// members are formed in a derived class, which C++ allows, and applied to the mesh.
//
// The writer below also relies on how MFEM's serial reader numbers a nonconforming mesh,
// which the unit test checks against a gathered save: the leaves come depth first from
// the roots in file order (NCMesh::CollectLeafElements); the vertices are numbered top
// level first, by node id, then in the order in which the leaves first reach them
// (NCMesh::UpdateVertices); the edges and faces in the order in which the elements first
// reach them (Mesh::GetElementToEdgeTable, Mesh::GetElementToFaceTable), and each face
// is oriented as the first element sees it (Mesh::GenerateFaces).
struct NCMeshData : public mfem::NCMesh
{
  static const auto &Elements(const mfem::NCMesh &m) { return m.*&NCMeshData::elements; }
  static const auto &Nodes(const mfem::NCMesh &m) { return m.*&NCMeshData::nodes; }
  static const auto &Faces(const mfem::NCMesh &m) { return m.*&NCMeshData::faces; }
  static const auto &Leaves(const mfem::NCMesh &m) { return m.*&NCMeshData::leaf_elements; }
  static const auto &RootState(const mfem::NCMesh &m) { return m.*&NCMeshData::root_state; }
  static auto &Coordinates(mfem::NCMesh &m) { return m.*&NCMeshData::coordinates; }
  static bool Scaled(const mfem::NCMesh &m) { return m.*&NCMeshData::using_scaling; }
  static bool Isotropic(const mfem::NCMesh &m) { return m.*&NCMeshData::Iso; }
  static const auto &Geom(int geom)
  {
    if (!GI[geom].initialized)
    {
      GI[geom].InitGeom(mfem::Geometry::Type(geom));
    }
    return GI[geom];
  }
  static int EdgeId(const mfem::NCMesh &m, const Element &el, int j)
  {
    const auto &ev = GI[el.Geom()].edges[j];
    return Nodes(m).FindId(el.node[ev[0]], el.node[ev[1]]);
  }
  static int FaceId(const mfem::NCMesh &m, const Element &el, int j)
  {
    const auto &fv = GI[el.Geom()].faces[j];
    return Faces(m).FindId(el.node[fv[0]], el.node[fv[1]], el.node[fv[2]], el.node[fv[3]]);
  }
};

// The same for the vertices of a mesh with nodes.
struct MeshData : public mfem::Mesh
{
  static void UpdateVertices(mfem::Mesh &m)
  {
    (m.*&MeshData::SetVerticesFromNodes)(m.GetNodes());
  }
};

// Send lists of trivially copyable records to every process. Returns the records received,
// ordered by source process, and their offsets per source (P + 1).
template <typename T>
std::vector<T> ExchangeRecords(MPI_Comm comm, const std::vector<std::vector<T>> &send,
                               std::vector<int> &recv_displs)
{
  static_assert(std::is_trivially_copyable_v<T>);
  const int P = static_cast<int>(send.size());
  std::vector<int> send_counts(P), send_displs(P + 1, 0), recv_counts(P);
  for (int p = 0; p < P; p++)
  {
    send_counts[p] = static_cast<int>(send[p].size());
    send_displs[p + 1] = send_displs[p] + send_counts[p];
  }
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);
  recv_displs.assign(P + 1, 0);
  std::partial_sum(recv_counts.begin(), recv_counts.end(), recv_displs.begin() + 1);
  std::vector<T> sbuf, rbuf(recv_displs[P]);
  sbuf.reserve(send_displs[P]);
  for (const auto &s : send)
  {
    sbuf.insert(sbuf.end(), s.begin(), s.end());
  }
  MPI_Datatype type;
  MPI_Type_contiguous(sizeof(T), MPI_BYTE, &type);
  MPI_Type_commit(&type);
  MPI_Alltoallv(sbuf.data(), send_counts.data(), send_displs.data(), type, rbuf.data(),
                recv_counts.data(), recv_displs.data(), type, comm);
  MPI_Type_free(&type);
  return rbuf;
}

// Exclusive prefix sums over the processes, and totals, of a few counts.
std::pair<std::vector<long long>, std::vector<long long>>
ExscanCounts(MPI_Comm comm, std::vector<long long> counts)
{
  std::vector<long long> offsets(counts.size(), 0), totals = counts;
  MPI_Exscan(counts.data(), offsets.data(), static_cast<int>(counts.size()), MPI_LONG_LONG,
             MPI_SUM, comm);
  if (Mpi::Root(comm))
  {
    std::fill(offsets.begin(), offsets.end(), 0);
  }
  Mpi::GlobalSum(static_cast<int>(totals.size()), totals.data(), comm);
  return {offsets, totals};
}

// For each root of the refinement tree: the own leaves under it, which come consecutively
// in the mesh, in the order of the roots; whether all its leaves are own; the process with
// most of its leaves; the index of its first leaf in the order of the processes and their
// leaves; and its number of leaves. The home process of a root on several processes finds
// the last three.
struct RootLeaves
{
  int leaves = 0;
  bool whole = true;
  int target = -1;
  long long first = 0, total = 0;
};

std::vector<RootLeaves> GetRootLeaves(const mfem::ParMesh &mesh, bool &split)
{
  MPI_Comm comm = mesh.GetComm();
  const int rank = Mpi::Rank(comm), P = Mpi::Size(comm);
  const auto &E = NCMeshData::Elements(*mesh.pncmesh);
  std::function<void(int, RootLeaves &)> Count = [&](int e, RootLeaves &r)
  {
    const auto &el = E[e];
    for (int j = 0; el.ref_type && j < mfem::NCMesh::MaxElemChildren && el.child[j] >= 0;
         j++)
    {
      Count(el.child[j], r);
    }
    r.leaves += (!el.ref_type && el.rank == rank);
    r.whole = r.whole && (el.ref_type || el.rank == rank);
  };
  std::vector<RootLeaves> roots(NCMeshData::RootState(*mesh.pncmesh).Size());
  long long first = ExscanCounts(comm, {mesh.GetNE()}).first[0];
  split = false;
  for (std::size_t i = 0; i < roots.size(); first += roots[i++].leaves)
  {
    Count(static_cast<int>(i), roots[i]);
    roots[i].target = rank;
    roots[i].first = first;
    roots[i].total = roots[i].leaves;
    split = split || (roots[i].leaves && !roots[i].whole);
  }
  Mpi::GlobalOr(1, &split, comm);
  if (!split)
  {
    return roots;
  }
  std::vector<std::vector<std::array<long long, 4>>> share(P);  // root, process, leaves,
  for (std::size_t i = 0; i < roots.size(); i++)                // first leaf
  {
    if (roots[i].leaves && !roots[i].whole)
    {
      share[i % P].push_back({(long long)i, rank, roots[i].leaves, roots[i].first});
    }
  }
  std::vector<int> displs;
  const auto in = ExchangeRecords(comm, share, displs);
  struct Best
  {
    long long rank = -1, leaves = 0, first = std::numeric_limits<long long>::max(),
              total = 0;
  };
  std::unordered_map<long long, Best> best;  // the most leaves, then the lowest process
  for (const auto &s : in)
  {
    auto &b = best[s[0]];
    if (std::tie(s[2], b.rank) > std::tie(b.leaves, s[1]))
    {
      std::tie(b.rank, b.leaves) = std::tie(s[1], s[2]);
    }
    b.first = std::min(b.first, s[3]);
    b.total += s[2];
  }
  std::vector<std::vector<std::array<long long, 3>>> reply(P);
  for (int p = 0; p < P; p++)
  {
    for (int i = displs[p]; i < displs[p + 1]; i++)
    {
      const auto &b = best[in[i][0]];
      reply[p].push_back({b.rank, b.first, b.total});
    }
  }
  const auto replies = ExchangeRecords(comm, reply, displs);
  for (std::size_t i = 0; i < roots.size(); i++)
  {
    if (roots[i].leaves && !roots[i].whole)
    {
      const auto &r = replies[displs[i % P]++];
      roots[i].target = static_cast<int>(r[0]);
      std::tie(roots[i].first, roots[i].total) = std::tie(r[1], r[2]);
    }
  }
  return roots;
}

// The numbers in the file of the vertices of a mesh with whole roots (see
// PrintParallelNCMesh): the flags of the vertices, edges and faces on the leaves of lower
// processes, which belong to them (1, 2 and 4); the own top-level and other vertices, in
// the order in which the own leaves reach them; and the numbers of the vertices, those of
// other processes from their owner, the leader of their group. Every vertex numbered is a
// corner of an own leaf, or a parent of an own vertex, which is a corner of an own element
// since the roots are whole. Collective.
struct FileVertices
{
  std::vector<char> skip;
  std::vector<long long> vid, vtot;
  std::vector<int> own[2];
  bool ok = true;
};

FileVertices NumberFileVertices(const mfem::ParMesh &mesh)
{
  using NC = NCMeshData;
  MPI_Comm comm = mesh.GetComm();
  const int rank = Mpi::Rank(comm), dim = mesh.Dimension(), ne = mesh.GetNE();
  const mfem::NCMesh &nc = *mesh.pncmesh;
  const auto &E = NC::Elements(nc);
  const auto &nodes = NC::Nodes(nc);
  const auto &leaf = NC::Leaves(nc);
  FileVertices fv;
  auto &skip = fv.skip;
  auto &vid = fv.vid;
  auto &own = fv.own;
  skip.assign(std::max(nodes.NumIds(), NC::Faces(nc).NumIds()), 0);
  for (int i = ne; i < ne + nc.GetNGhostElements(); i++)
  {
    const auto &el = E[leaf[i]];
    const auto &gi = NC::Geom(el.Geom());
    if (el.rank > rank)
    {
      continue;
    }
    for (int j = 0; j < gi.nv; j++)
    {
      skip[el.node[j]] |= 1;
    }
    for (int j = 0; j < gi.ne; j++)
    {
      skip[NC::EdgeId(nc, el, j)] |= 2;
    }
    for (int j = 0; dim == 3 && j < gi.nf; j++)
    {
      skip[NC::FaceId(nc, el, j)] |= 4;
    }
  }
  vid.assign(nodes.NumIds(), -1);
  for (int e = 0; e < ne; e++)
  {
    const auto &el = E[leaf[e]];
    for (int j = 0; j < NC::Geom(el.Geom()).nv; j++)
    {
      const int n = el.node[j];
      if (!(skip[n] & 1) && vid[n] == -1)
      {
        const bool top = (nodes[n].p1 == nodes[n].p2);
        vid[n] = own[!top].size();
        own[!top].push_back(n);
      }
    }
  }
  const auto [voff, vtot] =
      ExscanCounts(comm, {(long long)own[0].size(), (long long)own[1].size()});
  fv.vtot = vtot;
  for (int k = 0; k < 2; k++)
  {
    for (int n : own[k])
    {
      vid[n] += (k ? vtot[0] : 0) + voff[k];
    }
  }
  mfem::Array<int> group(mesh.GetNV()), number(mesh.GetNV());
  group = 0;
  number = -1;
  for (int g = 1; g < mesh.GetNGroups(); g++)
  {
    for (int i = 0; i < mesh.GroupNVertices(g); i++)
    {
      group[mesh.GroupVertex(g, i)] = g;
    }
  }
  for (int k = 0; k < 2; k++)
  {
    for (int n : own[k])
    {
      fv.ok = fv.ok && mesh.gtopo.IAmMaster(group[nodes[n].vert_index]);
      number[nodes[n].vert_index] = static_cast<int>(vid[n]);
    }
  }
  mfem::GroupCommunicator gc(mesh.gtopo);
  gc.Create(group);
  gc.Bcast(number);
  auto Number = [&](int n)
  {
    const int v = nodes[n].vert_index;
    fv.ok = fv.ok && v >= 0 && v < mesh.GetNV() && number[v] >= 0;
    vid[n] = fv.ok ? number[v] : -1;
  };
  for (int e = 0; e < ne; e++)
  {
    const auto &el = E[leaf[e]];
    for (int j = 0; j < NC::Geom(el.Geom()).nv; j++)
    {
      Number(el.node[j]);
    }
  }
  for (int n : own[1])
  {
    Number(nodes[n].p1);
    Number(nodes[n].p2);
  }
  return fv;
}

// The DOFs of the nodes of the own vertices (top-level, then the others), edges, faces and
// elements of a mesh, in the order of the file: each from the leaf which reaches it first,
// in its own (reference) order, the DOFs of an edge from its lower vertex number.
std::array<std::vector<int>, 5> FileDofs(const mfem::ParMesh &mesh, const FileVertices &fv)
{
  using NC = NCMeshData;
  using Geometry = mfem::Geometry;
  const mfem::NCMesh &nc = *mesh.pncmesh;
  const auto &E = NC::Elements(nc);
  const auto &nodes = NC::Nodes(nc);
  const auto &leaf = NC::Leaves(nc);
  const auto *fes = mesh.GetNodes()->FESpace();
  const auto *fec = fes->FEColl();
  const int nvd = fec->DofForGeometry(Geometry::POINT);
  const int ned = fec->DofForGeometry(Geometry::SEGMENT);
  auto skip = fv.skip;
  std::array<std::vector<int>, 5> dofs;
  mfem::Array<int> edofs;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const auto &el = E[leaf[e]];
    const auto &gi = NC::Geom(el.Geom());
    fes->GetElementDofs(e, edofs);
    int k = 0;
    auto Add = [&](std::vector<int> &list, int n, bool add, bool reverse = false)
    {
      for (int i = 0; add && i < n; i++)
      {
        list.push_back(edofs[k + (reverse ? n - 1 - i : i)]);
      }
      k += n;
    };
    for (int j = 0; j < gi.nv; j++)
    {
      const int n = el.node[j];
      Add(dofs[nodes[n].p1 != nodes[n].p2], nvd, !(skip[n] & 1));
      skip[n] |= 1;
    }
    for (int j = 0; j < gi.ne; j++)
    {
      const int m = NC::EdgeId(nc, el, j);
      Add(dofs[2], ned, !(skip[m] & 2),
          fv.vid[el.node[gi.edges[j][0]]] > fv.vid[el.node[gi.edges[j][1]]]);
      skip[m] |= 2;
    }
    for (int j = 0; mesh.Dimension() == 3 && j < gi.nf; j++)
    {
      const int f = NC::FaceId(nc, el, j);
      const auto fgeom = (gi.nfv[j] == 3) ? Geometry::TRIANGLE : Geometry::SQUARE;
      Add(dofs[3], fec->DofForGeometry(fgeom), !(skip[f] & 4));
      skip[f] |= 4;
    }
    Add(dofs[4], edofs.Size() - k, true);
  }
  return dofs;
}

// The process of a root of the refinement tree in MFEM's default partition of the leaves,
// when its leaves are [first, first + leaves) of all `total`: that of its middle leaf.
int RootProcess(long long first, long long leaves, long long total, int P)
{
  return static_cast<int>(std::min<long long>(P - 1, (first + leaves / 2) * P / total));
}

}  // namespace

void RebalanceRootAligned(mfem::ParMesh &mesh, double tol)
{
  MPI_Comm comm = mesh.GetComm();
  const int P = Mpi::Size(comm);
  bool split;
  const auto roots = GetRootLeaves(mesh, split);
  const long long total = mesh.GetGlobalNE();
  mfem::Array<int> part(mesh.GetNE());
  std::vector<long long> load(P, 0);
  for (std::size_t i = 0, e = 0; i < roots.size(); e += roots[i++].leaves)
  {
    if (roots[i].leaves)
    {
      const int r = RootProcess(roots[i].first, roots[i].total, total, P);
      std::fill_n(part.begin() + e, roots[i].leaves, r);
      load[r] += roots[i].leaves;
    }
  }
  // Accept the snapped partition when its largest load is within half of the tolerance of
  // the mean, and its largest over its smallest (as for MaximumImbalance) within the
  // tolerance, so that no process is empty and the next rebalance is not due at once.
  long long max_load = 0, min_load;
  MPI_Reduce_scatter_block(load.data(), &max_load, 1, MPI_LONG_LONG, MPI_SUM, comm);
  min_load = max_load;
  Mpi::GlobalMax(1, &max_load, comm);
  Mpi::GlobalMin(1, &min_load, comm);
  if (max_load <= (1.0 + 0.5 * (tol - 1.0)) * total / P && max_load <= tol * min_load)
  {
    mesh.Rebalance(part);
  }
  else
  {
    mesh.Rebalance();
  }
}

void ScaleNCMeshCoordinates(mfem::Mesh &mesh, double length_scale)
{
  if (mesh.ncmesh && !mesh.GetNodes())
  {
    for (auto &x : NCMeshData::Coordinates(*mesh.ncmesh))
    {
      x *= length_scale;
    }
  }
}

// The roots of the refinement tree are listed process by process (a root whose leaves are
// on several processes first moves whole to one of them, in a copy of the mesh). The serial
// reader then orders its elements process by process, and numbers its vertices, edges,
// faces and nodes in the order in which the elements first reach them (see NCMeshData). So
// each entity belongs to the lowest process with a leaf on it, which the ghost layer tells,
// and each process numbers its own entities; the numbers of the vertices of other processes
// come from their owner, the leader of the vertex's group in the ParMesh.
bool PrintParallelNCMesh(const mfem::ParMesh &orig, double length_scale,
                         const std::string &path, bool compress, MeshEntityCounts *counts,
                         std::string *reason)
{
  using Geometry = mfem::Geometry;
  using NC = NCMeshData;
  MPI_Comm comm = orig.GetComm();
  const int rank = Mpi::Rank(comm), P = Mpi::Size(comm), dim = orig.Dimension();
  const auto *nodes0 = orig.GetNodes();
  const auto *fec0 = nodes0 ? nodes0->FESpace()->FEColl() : nullptr;
  std::string why;
  if (!orig.pncmesh || orig.NURBSext || (dim != 2 && dim != 3) || NC::Scaled(*orig.pncmesh))
  {
    why = "the mesh type";
  }
  else if (dim == 3 && !NC::Isotropic(*orig.pncmesh))
  {
    why = "anisotropic refinement";
  }
  else if (nodes0 && (nodes0->FESpace()->IsVariableOrder() ||
                      !(dynamic_cast<const mfem::H1_FECollection *>(fec0) ||
                        dynamic_cast<const mfem::L2_FECollection *>(fec0))))
  {
    why = "the nodes";
  }
  auto Fail = [&](const std::string &w)
  {
    if (reason)
    {
      *reason = w;
    }
    return false;
  };
  bool ok = why.empty();
  Mpi::GlobalAnd(1, &ok, comm);
  if (!ok)
  {
    return Fail(why.empty() ? "the mesh type on another process" : why);
  }

  // 1. Every root must belong to one process: a root on several processes moves, in a copy
  //    of the mesh, to the process with most of its leaves.
  bool split;
  auto roots = GetRootLeaves(orig, split);
  std::unique_ptr<mfem::ParMesh> copy;
  if (split)
  {
    mfem::Array<int> part(orig.GetNE());
    for (std::size_t i = 0, e = 0; i < roots.size(); e += roots[i++].leaves)
    {
      std::fill_n(part.begin() + e, roots[i].leaves, roots[i].target);
    }
    copy = std::make_unique<mfem::ParMesh>(orig);
    copy->Rebalance(part);
    roots = GetRootLeaves(*copy, split);
    if (split)
    {
      return Fail("roots still split after the move");
    }
  }
  const mfem::ParMesh &mesh = copy ? *copy : orig;
  const int ne = mesh.GetNE();
  const mfem::NCMesh &nc = *mesh.pncmesh;
  const auto &E = NC::Elements(nc);
  const auto &nodes = NC::Nodes(nc);
  const auto &faces = NC::Faces(nc);
  const auto &leaf = NC::Leaves(nc);
  const auto *x = mesh.GetNodes();
  const auto *fes = x ? x->FESpace() : nullptr;

  // 2. and 3. The numbers of the vertices in the file (see NumberFileVertices).
  auto fv = NumberFileVertices(mesh);
  ok = fv.ok;
  auto &skip = fv.skip;
  const auto &vid = fv.vid;
  const auto &own = fv.own;
  const auto &vtot = fv.vtot;

  // 4. The file numbers of the elements: the roots, then the others, depth first.
  std::vector<int> sub;
  std::function<void(int)> Sub = [&](int e)
  {
    for (int j = 0;
         E[e].ref_type && j < mfem::NCMesh::MaxElemChildren && E[e].child[j] >= 0; j++)
    {
      sub.push_back(E[e].child[j]);
      Sub(E[e].child[j]);
    }
  };
  std::vector<int> root;
  for (std::size_t i = 0; i < roots.size(); i++)
  {
    if (roots[i].leaves)
    {
      root.push_back(static_cast<int>(i));
      Sub(root.back());
    }
  }
  const auto [eoff, etot] =
      ExscanCounts(comm, {(long long)root.size(), (long long)sub.size()});
  std::vector<long long> fid(E.Size(), -1);  // the number of an element in the file
  for (std::size_t k = 0; k < root.size(); k++)
  {
    fid[root[k]] = eoff[0] + k;
  }
  for (std::size_t k = 0; k < sub.size(); k++)
  {
    fid[sub[k]] = etot[0] + eoff[1] + k;
  }
  Mpi::GlobalAnd(1, &ok, comm);
  if (!ok || etot[0] + etot[1] >= std::numeric_limits<int>::max() ||
      vtot[0] + vtot[1] >= std::numeric_limits<int>::max())
  {
    return Fail(ok ? "more than 2^31 elements or vertices" : "the vertex groups");
  }

  // 5. The text of the file. Each section is printed in turn by all processes (the root
  //    also prints the headers), into segments: the root elements, the other elements, the
  //    boundary header, the boundary, the vertex parents, the root states, then the
  //    coordinates or the nodes, and the end. The nodes are the DOFs of the own vertices
  //    (top-level first), edges and faces, in the order in which the own leaves reach them,
  //    then of the elements; each comes from the element which reaches it first, in its own
  //    (reference) order, edge DOFs running from the lower vertex number.
  enum Segment
  {
    ROOTS,
    ELEMENTS,
    BOUNDARY_HEADER,
    BOUNDARY,
    PARENTS,
    STATES,
    GEOMETRY
  };
  const bool byvdim = fes && fes->GetOrdering() == mfem::Ordering::byVDIM;
  const int vdim = fes ? fes->GetVDim() : mesh.SpaceDimension();
  const int nseg = GEOMETRY + (!fes ? 1 : (byvdim ? 5 : 5 * vdim)) + 1;  // with the end
  std::vector<std::ostringstream> seg(nseg);
  for (auto &os : seg)
  {
    os << std::scientific;
    os.precision(std::numeric_limits<double>::max_digits10);
  }
  const mfem::NCMesh::NCList *face_list = nullptr;
  if (counts && dim == 3)
  {
    mesh.pncmesh->GetEdgeList();  // ParNCMesh::BuildFaceList may need the edge owners
    face_list = &mesh.pncmesh->GetFaceList();
  }
  const auto dofs = fes ? FileDofs(mesh, fv) : std::array<std::vector<int>, 5>{};
  std::array<long long, 2 * Geometry::NumGeom + 1> counted{};  // cells, true faces, bdr
  for (int e = 0; e < ne; e++)
  {
    const auto &el = E[leaf[e]];
    const auto &gi = NC::Geom(el.Geom());
    counted[el.Geom()]++;
    for (int j = 0; j < gi.nf; j++)  // the boundary as NCMesh::PrintBoundary, faces in 3D
    {
      const int f = NC::FaceId(nc, el, j);
      if (faces[f].Boundary())
      {
        constexpr int nfv2geom[5] = {Geometry::INVALID, Geometry::POINT, Geometry::SEGMENT,
                                     Geometry::TRIANGLE, Geometry::SQUARE};
        seg[BOUNDARY] << faces[f].attribute << ' ' << nfv2geom[gi.nfv[j]];
        for (int i = 0; i < gi.nfv[j]; i++)
        {
          seg[BOUNDARY] << ' ' << vid[el.node[gi.faces[j][i * (dim == 2 ? 2 : 1)]]];
        }
        seg[BOUNDARY] << '\n';
        counted[2 * Geometry::NumGeom]++;
      }
      if (dim == 3)
      {
        const auto fgeom = (gi.nfv[j] == 3) ? Geometry::TRIANGLE : Geometry::SQUARE;
        counted[Geometry::NumGeom + fgeom] += !(skip[f] & 4) && face_list &&
                                              face_list->GetMeshIdType(faces[f].index) !=
                                                  mfem::NCMesh::NCList::MeshIdType::SLAVE;
        skip[f] |= 4;
      }
    }
  }
  Mpi::GlobalSum(static_cast<int>(counted.size()), counted.data(), comm);
  auto Line = [&](std::ostream &os, int e)
  {
    const auto &el = E[e];  // as NCMesh::Print, all the children of a refined pyramid too
    os << "0 " << el.attribute << ' ' << int(el.geom) << ' ' << int(el.ref_type);
    for (int j = 0; el.ref_type && j < mfem::NCMesh::MaxElemChildren && el.child[j] >= 0;
         j++)
    {
      os << ' ' << fid[el.child[j]];
    }
    for (int j = 0; !el.ref_type && j < mfem::NCMesh::MaxElemNodes && el.node[j] >= 0; j++)
    {
      os << ' ' << vid[el.node[j]];
    }
    os << '\n';
  };
  if (Mpi::Root(comm))
  {
    seg[ROOTS] << "MFEM NC mesh v1.0\n\n# NCMesh supported geometry types:\n"
                  "# SEGMENT     = 1\n# TRIANGLE    = 2\n# SQUARE      = 3\n"
                  "# TETRAHEDRON = 4\n# CUBE        = 5\n# PRISM       = 6\n"
                  "# PYRAMID     = 7\n\ndimension\n"
               << dim
               << "\n\nrank\n0\n\n# rank attr geom ref_type nodes/children\nelements\n"
               << etot[0] + etot[1] << "\n";
    seg[BOUNDARY_HEADER] << "\n# attr geom nodes\nboundary\n"
                         << counted[2 * Geometry::NumGeom] << "\n";
    seg[PARENTS] << "\n# vert_id p1 p2\nvertex_parents\n" << vtot[1] << "\n";
    seg[STATES] << "\n# root element orientation\nroot_state\n" << etot[0] << "\n";
    if (fes)
    {
      seg[GEOMETRY] << "\n# mesh curvature GridFunction\nnodes\n";
      fes->Save(seg[GEOMETRY]);
      seg[GEOMETRY] << '\n';
    }
    else
    {
      seg[GEOMETRY] << "\n# top-level node coordinates\ncoordinates\n"
                    << vtot[0] << '\n'
                    << vdim << '\n';
    }
    seg[nseg - 1] << "\nmfem_mesh_end\n";
  }
  for (int r : root)
  {
    Line(seg[ROOTS], r);
    seg[STATES] << NC::RootState(nc)[r] << '\n';
  }
  for (int e : sub)
  {
    Line(seg[ELEMENTS], e);
  }
  for (int n : own[1])
  {
    seg[PARENTS] << vid[n] << ' ' << vid[nodes[n].p1] << ' ' << vid[nodes[n].p2] << '\n';
  }
  for (std::size_t i = 0; !fes && i < own[0].size(); i++)
  {
    const double *v = mesh.GetVertex(nodes[own[0][i]].vert_index);
    for (int c = 0; c < vdim; c++)
    {
      seg[GEOMETRY] << v[c] * length_scale << (c + 1 < vdim ? ' ' : '\n');
    }
  }
  for (int c = 0; fes && c < (byvdim ? 1 : vdim); c++)
  {
    for (int cat = 0; cat < 5; cat++)
    {
      for (int d : dofs[cat])
      {
        for (int i = (byvdim ? 0 : c); i < (byvdim ? vdim : c + 1); i++)
        {
          seg[GEOMETRY + 5 * c + cat]
              << mfem::ZeroSubnormal((*x)(fes->DofToVDof(d, i)) * length_scale)
              << (byvdim && i + 1 < vdim ? ' ' : '\n');
        }
      }
    }
  }

  // 6. The root writes the parts of all processes in file order.
  std::vector<std::string> part(nseg);
  for (int s = 0; s < nseg; s++)
  {
    part[s] = std::move(seg[s]).str();
#if defined(MFEM_USE_ZLIB)
    if (compress && !part[s].empty())  // one gzip member per part
    {
      std::ostringstream z;
      zstr::ostream(z).write(part[s].data(), part[s].size());
      part[s] = std::move(z).str();
    }
#endif
    ok = ok && (Mpi::Root(comm) || part[s].size() < std::numeric_limits<int>::max());
  }
  Mpi::GlobalAnd(1, &ok, comm);
  if (!ok)
  {
    return Fail("a part of the file of 2 GiB or more");
  }
  WriteRootOutputFileParts(path, comm, part);
  if (counts && Mpi::Root(comm))
  {
    *counts = MeshEntityCounts{};
    counts->dim = dim;
    for (int g = 0; g < Geometry::NumGeom; g++)
    {
      if (counted[g])
      {
        counts->cells[Geometry::Type(g)] = counted[g];
      }
      if (counted[Geometry::NumGeom + g])
      {
        counts->true_faces[Geometry::Type(g)] = counted[Geometry::NumGeom + g];
      }
    }
    counts->domain_attributes = {mesh.attributes.begin(), mesh.attributes.end()};
    counts->boundary_attributes = {mesh.bdr_attributes.begin(), mesh.bdr_attributes.end()};
  }
  return true;
}

// What ReadNCMeshCoarse keeps for PartitionNCMeshFile: the process of each root (in file
// order, the order of the coarse elements), the refinement types of the elements of the
// trees of the own roots (level by level over these roots), and for the nodes, their order
// and ordering and the lines of the file which start in this process's slice of its bytes.
struct NCMeshFile
{
  std::vector<int> part;
  std::string tree, text;
  std::vector<std::size_t> start;  // of each own line in text, then the end
  std::vector<long long> first;    // the first line of each process (P + 1)
  long long g0 = 0, values = 0, ndofs = 0;
  int order = 0, vdim = 0;
  bool byvdim = true;
  bool Owns(long long g) const { return g >= g0 && g < g0 + long(start.size()) - 1; }
  std::string_view Line(long long g) const
  {
    const auto i = g - g0;
    return std::string_view(text).substr(start[i], start[i + 1] - 1 - start[i]);
  }
};

namespace
{

#if defined(MFEM_USE_ZLIB)
// The root decompresses a gzip-compressed nonconforming mesh and sends each process the
// lines which start in its slice of the compressed bytes (from the compressed offset after
// each read), in messages of at most 64 MiB. Collective.
bool DecompressFileLines(MPI_Comm comm, const std::string &path, std::string &text)
{
  constexpr int more = 1, last = 2;  // the tags of the messages
  constexpr int chunk = 1 << 16, message = 1 << 26;
  const int P = Mpi::Size(comm);
  bool ok = true;
  if (!Mpi::Root(comm))
  {
    for (MPI_Status s = {}; s.MPI_TAG != last;)
    {
      int n;
      MPI_Probe(0, MPI_ANY_TAG, comm, &s);
      MPI_Get_count(&s, MPI_CHAR, &n);
      text.resize(text.size() + n);
      MPI_Recv(text.data() + text.size() - n, n, MPI_CHAR, 0, s.MPI_TAG, comm,
               MPI_STATUS_IGNORE);
    }
  }
  else
  {
    std::ifstream fi(path, std::ios::binary | std::ios::ate);
    const long long size = fi ? static_cast<long long>(fi.tellg()) : 0;
    gzFile f = gzopen(path.c_str(), "rb");
    std::string buf(17, '\0');  // the bytes of process r not sent yet
    ok = f && gzbuffer(f, 1 << 20) == 0 && gzread(f, buf.data(), 17) == 17 &&
         buf == "MFEM NC mesh v1.0";
    int r = 0;
    auto Send = [&](std::size_t n, int tag)
    {
      if (r == 0)
      {
        text.append(buf, 0, n);
      }
      else
      {
        MPI_Send(buf.data(), static_cast<int>(n), MPI_CHAR, r, tag, comm);
      }
      buf.erase(0, n);
      r += (tag == last);
    };
    for (int n = 1; ok && n > 0;)
    {
      const auto n0 = buf.size();
      buf.resize(n0 + chunk);
      n = gzread(f, buf.data() + n0, chunk);
      ok = (n >= 0);
      buf.resize(n0 + std::max(n, 0));
      // When a read reaches the end of the slice of r, its lines end at the next end of
      // line.
      const long long in = gzoffset(f);
      for (std::size_t from = n0, eol; r + 1 < P && in >= size * (r + 1) / P &&
                                       (eol = buf.find('\n', from)) != std::string::npos;
           from = 0)
      {
        Send(eol + 1, last);
      }
      if (buf.size() >= message && (r + 1 == P || in < size * (r + 1) / P))
      {
        Send(buf.size(), more);
      }
    }
    int err = Z_ERRNO;
    if (f)
    {
      gzerror(f, &err);
      gzclose_r(f);
    }
    ok = ok && (err == Z_OK);
    while (r < P)  // the rest, then nothing
    {
      Send(buf.size(), last);
    }
  }
  Mpi::GlobalAnd(1, &ok, comm);
  return ok;
}
#endif

// Read the lines of a file which start in this process's slice of its bytes (with enough of
// the next bytes to end them), and number them. A compressed file is decompressed by the
// root. Collective.
bool ReadFileLines(MPI_Comm comm, const std::string &path, bool compressed,
                   NCMeshFile &file)
{
  const int rank = Mpi::Rank(comm), P = Mpi::Size(comm);
  bool ok = true;
  std::size_t p = 0, end = 0;  // the own lines start in [p, end)
  if (compressed)
  {
#if defined(MFEM_USE_ZLIB)
    ok = DecompressFileLines(comm, path, file.text);
#endif
    end = file.text.size();
  }
  else
  {
    std::ifstream fi(path, std::ios::binary | std::ios::ate);
    const long long size = fi ? static_cast<long long>(fi.tellg()) : 0;
    const long long a = size * rank / P, b = size * (rank + 1) / P;
    const long long lo = std::max(a - 1, 0LL), hi = std::min(b + (1LL << 16), size);
    file.text.assign(std::max(hi - lo, 0LL), '\0');
    fi.seekg(lo);
    fi.read(file.text.data(), file.text.size());
    ok = static_cast<bool>(fi);
    p = (a > 0) ? std::min(file.text.find('\n'), file.text.size() - 1) + 1 : 0;
    end = std::max(b - lo, 0LL);
  }
  while (ok && p < end)
  {
    const auto eol = file.text.find('\n', p);
    ok = (eol != std::string::npos);
    file.start.push_back(p);
    p = eol + 1;
  }
  file.start.push_back(p);
  long long n = static_cast<long long>(file.start.size()) - 1;
  file.first.assign(P + 1, 0);
  MPI_Allgather(&n, 1, MPI_LONG_LONG, file.first.data() + 1, 1, MPI_LONG_LONG, comm);
  std::partial_sum(file.first.begin(), file.first.end(), file.first.begin());
  file.g0 = file.first[rank];
  Mpi::GlobalAnd(1, &ok, comm);
  return ok;
}

// The numbers of a line (all of it, or false).
template <typename T>
bool LineNumbers(std::string_view line, std::vector<T> &v)
{
  v.clear();
  const char *p = line.data(), *end = p + line.size();
  for (T x; p < end; p += (p < end && *p == ' '))
  {
    const char *q = p;
    if constexpr (std::is_floating_point_v<T>)
    {
      // Not std::from_chars, which libc++ only has for floating point from LLVM 20. The
      // line ends at an end of line, so strtod ends there at the latest, but for spaces.
      char *r;
      x = std::strtod(p, &r);
      q = r;
    }
    else
    {
      const auto r = std::from_chars(p, end, x);
      q = (r.ec == std::errc()) ? r.ptr : p;
    }
    if (q == p || q > end)
    {
      return false;
    }
    v.push_back(x);
    p = q;
  }
  return !v.empty();
}

// The records of all processes, in their order, on the root.
template <typename T>
std::vector<T> GatherRecords(MPI_Comm comm, std::vector<T> records)
{
  std::vector<std::vector<T>> send(Mpi::Size(comm));
  send[0] = std::move(records);
  std::vector<int> displs;
  return ExchangeRecords(comm, send, displs);
}

// The child of a refined element which holds a corner of it, or a face, as MFEM derefines
// it (NCMesh::DerefineElement).
int DerefinementChild(int geom, int ref_type, int i, bool face)
{
  using Geometry = mfem::Geometry;
  const int r = ref_type - 1;
  switch (geom)
  {
    case Geometry::TRIANGLE:
      return i;
    case Geometry::TETRAHEDRON:
      return face ? (i + 1) & 3 : i;
    case Geometry::SQUARE:
      return mfem::quad_deref_table[r][i + 4 * face];
    case Geometry::CUBE:
      return mfem::hex_deref_table[r][i + 8 * face];
    case Geometry::PRISM:
      return mfem::prism_deref_table[r][i + 6 * face];
    case Geometry::PYRAMID:
      return mfem::pyramid_deref_table[r][i + 5 * face];
  }
  return -1;
}

// Set the nodes of a mesh refined as a file to the nodes of the file, divided by the length
// scale: each process asks the processes which read them for the lines of its own DOFs, in
// the order of the file (see FileDofs), and the others follow from the true DOFs.
// Collective.
bool ReadFileNodes(mfem::ParMesh &mesh, const NCMeshFile &file, double length_scale)
{
  MPI_Comm comm = mesh.GetComm();
  const int P = Mpi::Size(comm);
  auto &x = static_cast<mfem::ParGridFunction &>(*mesh.GetNodes());
  const auto fv = NumberFileVertices(mesh);
  const auto dofs = FileDofs(mesh, fv);
  std::vector<long long> n(5);
  for (int c = 0; c < 5; c++)
  {
    n[c] = static_cast<long long>(dofs[c].size());
  }
  const auto [off, tot] = ExscanCounts(comm, n);
  bool ok = fv.ok && std::accumulate(tot.begin(), tot.end(), 0LL) == file.ndofs;
  for (int v : fv.own[0])  // the top-level vertices are numbered as in the file
  {
    ok = ok && fv.vid[v] == v;
  }
  Mpi::GlobalAnd(1, &ok, comm);
  if (!ok)
  {
    return false;
  }
  struct Piece
  {
    int p, cat, c;
    long long line, m, i;
  };
  std::vector<Piece> pieces;  // the own lines of each component and category, by process
  std::vector<std::vector<std::array<long long, 2>>> ask(P);
  for (int c = 0; c < (file.byvdim ? 1 : file.vdim); c++)
  {
    for (long long cat = 0, base = 0; cat < 5; base += tot[cat++])
    {
      long long line = file.values + c * file.ndofs + base + off[cat];
      for (long long i = 0; i < n[cat]; i += pieces.back().m, line += pieces.back().m)
      {
        const int p =
            static_cast<int>(std::upper_bound(file.first.begin(), file.first.end(), line) -
                             file.first.begin() - 1);
        pieces.push_back(
            {p, int(cat), c, line, std::min(n[cat] - i, file.first[p + 1] - line), i});
        ask[p].push_back({line, pieces.back().m});
      }
    }
  }
  std::vector<int> displs;
  const auto asked = ExchangeRecords(comm, ask, displs);
  const int nval = file.byvdim ? file.vdim : 1;  // per line
  std::vector<std::vector<double>> answer(P);
  std::vector<double> v;
  for (int p = 0; p < P; p++)
  {
    for (int i = displs[p]; i < displs[p + 1]; i++)
    {
      for (long long g = asked[i][0]; g < asked[i][0] + asked[i][1]; g++)
      {
        ok = ok && file.Owns(g) && LineNumbers(file.Line(g), v) && int(v.size()) == nval;
        answer[p].insert(answer[p].end(), ok ? v.begin() : v.end(), v.end());
        answer[p].resize(answer[p].size() + (ok ? 0 : nval));
      }
    }
  }
  const auto values = ExchangeRecords(comm, answer, displs);
  std::vector<int> cursor(displs.begin(), displs.end() - 1);
  auto *fes = x.ParFESpace();
  x = std::numeric_limits<double>::quiet_NaN();
  for (const auto &q : pieces)
  {
    for (long long i = q.i; i < q.i + q.m; i++)
    {
      for (int k = 0; k < nval; k++)
      {
        x(fes->DofToVDof(dofs[q.cat][i], file.byvdim ? k : q.c)) =
            values[cursor[q.p]++] / length_scale;
      }
    }
  }
  mfem::Vector tv;
  x.GetTrueDofs(tv);
  for (int i = 0; i < tv.Size(); i++)
  {
    ok = ok && std::isfinite(tv(i));
  }
  x.SetFromTrueDofs(tv);
  MeshData::UpdateVertices(mesh);
  Mpi::GlobalAnd(1, &ok, comm);
  return ok;
}

}  // namespace

// Every process reads the lines which start in its slice of the file (the root
// decompresses a compressed file and sends them). The roots, whose trees are consecutive
// in the file (PrintParallelNCMesh lists the elements of each root depth first), are
// partitioned into consecutive ranges with about the same number of leaves, and each
// process gets the elements of its roots: their corners come from the leaves through the
// children which hold them, as MFEM derefines, and so do the attributes of their faces,
// looked up in the boundary of the file. The root prints the coarse mesh, with the
// top-level vertex coordinates (the first vertex DOFs of the nodes).
std::shared_ptr<NCMeshFile> ReadNCMeshCoarse(MPI_Comm comm, const std::string &path,
                                             bool enabled, std::string &coarse)
{
  using NC = NCMeshData;
  using Geometry = mfem::Geometry;
  const int rank = Mpi::Rank(comm), P = Mpi::Size(comm);
  coarse.clear();
  auto Fail = []() { return std::shared_ptr<NCMeshFile>(); };
  if (!enabled)
  {
    return Fail();
  }
  auto kept = std::make_shared<NCMeshFile>();
  auto &file = *kept;
  int format = 0;  // 1: a nonconforming mesh, 2: compressed (its header is checked later)
  if (Mpi::Root(comm))
  {
    std::ifstream fi(path, std::ios::binary);
    std::string head(17, '\0');
    fi.read(head.data(), head.size());
    format = (fi && head == "MFEM NC mesh v1.0") ? 1 : 0;
#if defined(MFEM_USE_ZLIB)
    format = (head.compare(0, 2, "\x1f\x8b") == 0) ? 2 : format;
#endif
  }
  Mpi::Broadcast(1, &format, 0, comm);
  if (format == 0 || !ReadFileLines(comm, path, format == 2, file))
  {
    return Fail();
  }
  bool ok = true;

  // 1. The sections: the line of each keyword, then the numbers after them.
  constexpr std::string_view keys[] = {"dimension",    "elements",    "boundary",
                                       "root_state",   "coordinates", "nodes",
                                       "mfem_mesh_end"};
  std::array<long long, 7> at;
  at.fill(-1);
  for (long long g = file.g0; g < file.first[rank + 1]; g++)
  {
    for (int k = 0; k < 7; k++)
    {
      at[k] = (file.Line(g) == keys[k]) ? g : at[k];
    }
  }
  Mpi::GlobalMax(7, at.data(), comm);
  const auto [dim_at, el_at, bd_at, rs_at, co_at, nd_at, end_at] = at;
  auto Int = [&](long long g, std::string_view prefix = "")
  {
    std::vector<long long> v;
    const auto line = file.Owns(g) ? file.Line(g) : std::string_view();
    return (line.substr(0, prefix.size()) == prefix &&
            LineNumbers(line.substr(prefix.size()), v))
               ? v[0]
               : -1LL;
  };
  const auto fec = file.Owns(nd_at + 2) ? file.Line(nd_at + 2) : std::string_view();
  std::array<long long, 8> num = {Int(dim_at + 1),
                                  Int(el_at + 1),
                                  Int(bd_at + 1),
                                  Int(rs_at + 1),
                                  (co_at >= 0) ? Int(co_at + 2) : Int(nd_at + 3, "VDim: "),
                                  Int(nd_at + 4, "Ordering: "),
                                  (fec.substr(0, 28) == "FiniteElementCollection: H1_")
                                      ? Int(nd_at + 2, fec.substr(0, fec.rfind('P') + 1))
                                      : -1,
                                  (co_at >= 0) ? Int(co_at + 1) : 0};
  Mpi::GlobalMax(8, num.data(), comm);
  const long long dim = num[0], nel = num[1], nbd = num[2], R = num[3], vdim = num[4],
                  ordering = num[5], order = num[6], ncoord = num[7];
  file.order = (nd_at >= 0) ? static_cast<int>(order) : 0;
  file.vdim = static_cast<int>(vdim);
  file.byvdim = (ordering == mfem::Ordering::byVDIM);
  file.values = nd_at + 6;
  file.ndofs = (end_at - 1 - file.values) / (file.byvdim ? 1 : std::max(vdim, 1LL));
  if (!((dim == 2 || dim == 3) && nel > 0 && R > 0 && R <= nel &&
        nel < std::numeric_limits<int>::max() && nbd >= 0 && end_at > rs_at && vdim > 0 &&
        (co_at >= 0) != (nd_at >= 0) && (nd_at < 0 || (file.order > 0 && file.ndofs > 0))))
  {
    return Fail();
  }

  // 2. The elements of the own lines (rank, attribute, geometry, type, then the nodes or
  //    the children), and the first child of each refined root: the other elements of the
  //    trees of the roots come in the order of the roots. The number of leaves of each
  //    root.
  struct Elem
  {
    int id, attr, geom, ref;
    std::array<int, mfem::NCMesh::MaxElemChildren> ids;
  };
  std::vector<Elem> lines;
  std::vector<int> start(R + 1, -1), leaves(R, 0), v;
  start[R] = static_cast<int>(nel);
  for (long long g = std::max(file.g0, el_at + 2);
       g < std::min(file.first[rank + 1], el_at + 2 + nel); g++)
  {
    ok = ok && LineNumbers(file.Line(g), v) && v.size() >= 5 &&
         v.size() <= 4 + mfem::NCMesh::MaxElemChildren && v[2] >= Geometry::TRIANGLE &&
         v[2] <= Geometry::PYRAMID && (v[2] >= Geometry::TETRAHEDRON) == (dim == 3) &&
         (dim == 3 ? (v[3] == 0 || v[3] == 7) : (v[3] >= 0 && v[3] <= 3)) &&
         (v[3] || int(v.size()) - 4 == Geometry::NumVerts[v[2]]);
    if (ok)
    {
      Elem &el = lines.emplace_back(Elem{int(g - el_at - 2), v[1], v[2], v[3], {}});
      el.ids.fill(-1);
      std::copy(v.begin() + 4, v.end(), el.ids.begin());
      if (el.id < R && el.ref)
      {
        start[el.id] = el.ids[0];
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, start.data(), static_cast<int>(R + 1), MPI_INT, MPI_MAX,
                comm);
  for (long long k = R - 1; k >= 0; k--)
  {
    start[k] = (start[k] >= 0) ? start[k] : start[k + 1];
  }
  ok = ok && start[0] == R && std::is_sorted(start.begin(), start.end());
  auto RootOf = [&](int e)
  {
    return (e < R) ? e
                   : static_cast<int>(std::upper_bound(start.begin(), start.end(), e) -
                                      start.begin() - 1);
  };
  for (const auto &el : lines)
  {
    const int r = RootOf(el.id);
    ok = ok && r >= 0 && r < R;
    leaves[ok ? r : 0] += (el.ref == 0);
  }
  long long nlines = static_cast<long long>(lines.size());
  MPI_Allreduce(MPI_IN_PLACE, leaves.data(), static_cast<int>(R), MPI_INT, MPI_SUM, comm);
  Mpi::GlobalSum(1, &nlines, comm);
  ok = ok && nlines == nel;
  Mpi::GlobalAnd(1, &ok, comm);
  if (!ok)
  {
    return Fail();
  }

  // 3. The partition of the roots, and the elements of the own roots.
  const long long total = std::accumulate(leaves.begin(), leaves.end(), 0LL);
  file.part.resize(R);
  for (long long k = 0, f = 0; k < R; f += leaves[k++])
  {
    file.part[k] = RootProcess(f, leaves[k], total, P);
  }
  std::vector<std::vector<Elem>> send(P);
  for (const auto &el : lines)
  {
    send[file.part[RootOf(el.id)]].push_back(el);
  }
  lines.clear();
  std::vector<int> displs;
  const auto recv = ExchangeRecords(comm, send, displs);
  send.clear();
  const auto own = std::equal_range(file.part.begin(), file.part.end(), rank);
  const int k0 = int(own.first - file.part.begin()),
            k1 = int(own.second - file.part.begin());
  const int s0 = start[k0], s1 = start[k1];
  auto Index = [&](int e)  // in the elements of the own roots
  {
    return (e < R) ? ((e >= k0 && e < k1) ? e - k0 : -1)
                   : ((e >= s0 && e < s1) ? (k1 - k0) + (e - s0) : -1);
  };
  std::vector<Elem> loc((k1 - k0) + (s1 - s0));
  for (const auto &el : recv)
  {
    const int j = Index(el.id);
    ok = ok && j >= 0;
    if (j >= 0)
    {
      loc[j] = el;
    }
  }
  ok = ok && recv.size() == loc.size();

  // 4. The trees of the own roots, level by level, and their corners and faces, from their
  //    leaves; each face is looked up, by its leaf vertices, in the boundary of the file.
  std::vector<char> seen(loc.size(), 0);
  std::vector<int> level(k1 - k0), next;
  std::iota(level.begin(), level.end(), 0);
  std::fill_n(seen.begin(), k1 - k0, 1);
  for (; ok && !level.empty(); level.swap(next))
  {
    next.clear();
    for (int e : level)
    {
      file.tree += static_cast<char>('0' + loc[e].ref);
      for (int j = 0;
           ok && loc[e].ref && j < mfem::NCMesh::MaxElemChildren && loc[e].ids[j] >= 0; j++)
      {
        const int c = Index(loc[e].ids[j]);
        ok = c >= 0 && !seen[c] && loc[c].attr == loc[e].attr;  // as MFEM refines
        seen[std::max(c, 0)] = 1;
        next.push_back(c);
      }
    }
  }
  ok = ok && std::all_of(seen.begin(), seen.end(), [](char s) { return s; });
  auto Leaf = [&](int e, int i, bool face)  // the leaf which holds a corner or a face
  {
    while (ok && loc[e].ref)
    {
      const int c = DerefinementChild(loc[e].geom, loc[e].ref, i, face);
      const int j = (c >= 0) ? Index(loc[e].ids[c]) : -1;
      ok = j >= 0 && loc[j].geom == loc[e].geom;
      e = ok ? j : e;
    }
    return e;
  };
  // A face: its vertices, sorted, and for a boundary face its attribute; for a root face
  // looked up, -1, and its corners and its geometry.
  struct Face
  {
    std::array<int, 4> v, c;
    int attr, geom;
  };
  constexpr int nfv2geom[5] = {Geometry::INVALID, Geometry::POINT, Geometry::SEGMENT,
                               Geometry::TRIANGLE, Geometry::SQUARE};
  struct Root
  {
    int attr, geom;
    std::array<int, 8> v;
  };
  std::vector<Root> roots;
  std::vector<std::vector<Face>> faces(P);
  long long top = -1, prev = -1, local_top = -1;  // the largest top-level vertex numbers
  for (int k = 0; ok && k < k1 - k0; k++)
  {
    const auto &gi = NC::Geom(loc[k].geom);
    Root &root = roots.emplace_back(Root{loc[k].attr, loc[k].geom, {}});
    for (int i = 0; ok && i < gi.nv; i++)
    {
      root.v[i] = loc[Leaf(k, i, false)].ids[i];
      local_top = std::max<long long>(local_top, root.v[i]);
    }
    for (int f = 0; ok && f < gi.nf; f++)
    {
      const auto &leaf = loc[Leaf(k, f, true)];
      Face q{{-1, -1, -1, -1}, {-1, -1, -1, -1}, -1, nfv2geom[gi.nfv[f]]};
      for (int i = 0; i < gi.nfv[f]; i++)
      {
        q.v[i] = leaf.ids[gi.faces[f][i * (dim == 2 ? 2 : 1)]];
        q.c[i] = root.v[gi.faces[f][i * (dim == 2 ? 2 : 1)]];
      }
      std::sort(q.v.begin(), q.v.begin() + gi.nfv[f]);
      faces[std::max(q.v[0], 0) % P].push_back(q);
    }
  }
  for (long long g = std::max(file.g0, bd_at + 2);
       g < std::min(file.first[rank + 1], bd_at + 2 + nbd); g++)
  {
    // The attribute, the geometry, and the vertices.
    ok = ok && LineNumbers(file.Line(g), v) && v.size() >= 4 && v.size() <= 6;
    if (ok)
    {
      Face b{{-1, -1, -1, -1}, {}, v[0], 0};
      std::copy(v.begin() + 2, v.end(), b.v.begin());
      std::sort(b.v.begin(), b.v.begin() + (v.size() - 2));
      faces[std::max(b.v[0], 0) % P].push_back(b);
    }
  }
  const auto in = ExchangeRecords(comm, faces, displs);
  faces.clear();
  std::map<std::array<int, 4>, int> attribute;
  std::vector<Face> boundary;
  for (const auto &f : in)
  {
    if (f.attr >= 0)
    {
      attribute[f.v] = f.attr;
    }
  }
  for (const auto &f : in)
  {
    if (const auto it = attribute.find(f.v); f.attr < 0 && it != attribute.end())
    {
      boundary.push_back({f.c, {}, it->second, f.geom});
    }
  }

  // 5. The top-level vertices are numbered in the order in which the roots reach them, as
  //    PrintParallelNCMesh numbers them, checked only where the nodes depend on it.
  MPI_Exscan(&local_top, &prev, 1, MPI_LONG_LONG, MPI_MAX, comm);
  prev = (rank == 0) ? -1 : prev;
  for (const auto &root : roots)
  {
    long long fresh = 0, high = prev;
    for (int i = 0; i < Geometry::NumVerts[root.geom]; i++)
    {
      fresh += (root.v[i] > prev);
      high = std::max<long long>(high, root.v[i]);
    }
    ok = ok && (file.order == 0 || fresh == high - prev);
    prev = high;
  }
  top = local_top;
  Mpi::GlobalMax(1, &top, comm);
  top++;
  ok = ok && top <= (nd_at >= 0 ? file.ndofs : ncoord);

  // 6. The root states and the coordinates of the top-level vertices of the own lines, and
  //    the coarse mesh, printed on the root and sent to all.
  std::vector<int> states;
  std::vector<double> x, w;
  for (long long g = std::max(file.g0, rs_at + 2);
       g < std::min(file.first[rank + 1], rs_at + 2 + R); g++)
  {
    ok = ok && LineNumbers(file.Line(g), v) && v.size() == 1;
    states.push_back(ok ? v[0] : 0);
  }
  const long long c0 = (nd_at >= 0) ? file.values : co_at + 3;
  for (long long c = 0; c < ((nd_at >= 0 && !file.byvdim) ? vdim : 1); c++)
  {
    const long long a = c0 + c * file.ndofs;
    for (long long g = std::max(file.g0, a); g < std::min(file.first[rank + 1], a + top);
         g++)
    {
      ok = ok && LineNumbers(file.Line(g), w) &&
           long(w.size()) == ((nd_at >= 0 && !file.byvdim) ? 1 : vdim);
      x.insert(x.end(), w.begin(), w.end());
    }
  }
  roots = GatherRecords(comm, std::move(roots));
  boundary = GatherRecords(comm, std::move(boundary));
  states = GatherRecords(comm, std::move(states));
  x = GatherRecords(comm, std::move(x));
  Mpi::GlobalAnd(1, &ok, comm);
  if (!ok)
  {
    return Fail();
  }
  if (Mpi::Root(comm))
  {
    std::ostringstream os;
    os.precision(std::numeric_limits<double>::max_digits10);
    os << "MFEM NC mesh v1.0\n\ndimension\n" << dim << "\n\nelements\n" << R << '\n';
    for (const auto &root : roots)
    {
      os << "0 " << root.attr << ' ' << root.geom << " 0";
      for (int i = 0; i < Geometry::NumVerts[root.geom]; i++)
      {
        os << ' ' << root.v[i];
      }
      os << '\n';
    }
    os << "\nboundary\n" << boundary.size() << '\n';
    for (const auto &b : boundary)
    {
      os << b.attr << ' ' << b.geom;
      for (int i = 0; i < Geometry::NumVerts[b.geom]; i++)
      {
        os << ' ' << b.v[i];
      }
      os << '\n';
    }
    os << "\nroot_state\n" << R << '\n';
    for (int s : states)
    {
      os << s << '\n';
    }
    os << "\ncoordinates\n" << top << '\n' << vdim << '\n';
    for (long long i = 0; i < top * vdim; i++)  // byNODES: one component after the other
    {
      os << x[(nd_at < 0 || file.byvdim) ? i : (i % vdim) * top + i / vdim]
         << ((i + 1) % vdim ? ' ' : '\n');
    }
    os << "\nmfem_mesh_end\n";
    coarse = std::move(os).str();
  }
  std::int64_t size = static_cast<std::int64_t>(coarse.size());
  Mpi::Broadcast(1, &size, 0, comm);
  coarse.resize(size);
  Mpi::BroadcastLarge(size, coarse.data(), 0, comm);
  if (file.order == 0)
  {
    std::string().swap(file.text);
    std::vector<std::size_t>().swap(file.start);
  }
  return kept;
}

std::unique_ptr<mfem::ParMesh> PartitionNCMeshFile(MPI_Comm comm, NCMeshFile &file,
                                                   mfem::Mesh &coarse, double length_scale,
                                                   double tol)
{
  using NC = NCMeshData;
  bool ok = coarse.ncmesh && file.part.size() == static_cast<std::size_t>(coarse.GetNE());
  Mpi::GlobalAnd(1, &ok, comm);
  if (!ok)  // for example, if the preprocessing removed elements
  {
    return nullptr;
  }
  auto mesh = std::make_unique<mfem::ParMesh>(comm, coarse, file.part.data());

  // The refinements of the own roots, level by level.
  const int k0 = static_cast<int>(
      std::lower_bound(file.part.begin(), file.part.end(), Mpi::Rank(comm)) -
      file.part.begin());
  std::vector<int> level(mesh->GetNE()), next;
  for (int i = 0; i < mesh->GetNE(); i++)
  {
    level[i] = NC::Leaves(*mesh->pncmesh)[i];  // the roots, before the refinements
    ok = ok && level[i] == k0 + i;
  }
  const char iso = (1 << mesh->Dimension()) - 1;  // MFEM's isotropic refinement type
  std::size_t t = 0;
  for (bool more = true; more; level.swap(next))
  {
    const auto &E = NC::Elements(*mesh->pncmesh);
    mfem::Array<mfem::Refinement> refinements;
    for (int e : level)
    {
      ok = ok && t < file.tree.size() && !E[e].ref_type;
      const char type = ok ? file.tree[t++] - '0' : 0;
      if (type)  // isotropic as Refinement::XYZ, so that MFEM keeps the mesh isotropic
      {
        refinements.Append(
            mfem::Refinement(E[e].index, type == iso ? char(mfem::Refinement::XYZ) : type));
      }
    }
    more = refinements.Size() > 0;
    Mpi::GlobalOr(1, &more, comm);
    if (more)
    {
      mesh->GeneralRefinement(refinements, 1, 0);
    }
    next.clear();
    for (int e : level)
    {
      const auto &el = NC::Elements(*mesh->pncmesh)[e];
      for (int j = 0; el.ref_type && j < mfem::NCMesh::MaxElemChildren && el.child[j] >= 0;
           j++)
      {
        next.push_back(el.child[j]);
      }
    }
  }
  ok = ok && t == file.tree.size() &&
       mesh->GetNE() == std::count(file.tree.begin(), file.tree.end(), '0');
  Mpi::GlobalAnd(1, &ok, comm);
  if (ok && file.order > 0)  // the nodes of the file
  {
    mesh->SetCurvature(file.order, false, -1,
                       file.byvdim ? mfem::Ordering::byVDIM : mfem::Ordering::byNODES);
    ok = ReadFileNodes(*mesh, file, length_scale);
  }
  file = {};
  if (!ok)
  {
    return nullptr;
  }
  int ne_min = mesh->GetNE(), ne_max = ne_min;
  Mpi::GlobalMin(1, &ne_min, comm);
  Mpi::GlobalMax(1, &ne_max, comm);
  if (ne_max > tol * ne_min)
  {
    mesh->Rebalance();
  }
  return mesh;
}

}  // namespace palace::mesh
