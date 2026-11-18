// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ncmeshio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
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
  static const auto &Geom(int geom) { return GI[geom]; }
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

}  // namespace

void RebalanceRootAligned(mfem::ParMesh &mesh, double tol, const mfem::Array<int> *base)
{
  MPI_Comm comm = mesh.GetComm();
  const int P = Mpi::Size(comm);
  bool split;
  const auto roots = GetRootLeaves(mesh, split);
  const long long total = mesh.GetGlobalNE();
  mfem::Array<int> part(mesh.GetNE());
  for (std::size_t i = 0, e = 0; i < roots.size(); e += roots[i++].leaves)
  {
    const int n = roots[i].leaves;
    if (n && !base)
    {
      const int r = static_cast<int>((roots[i].first + roots[i].total / 2) * P / total);
      std::fill_n(part.begin() + e, n, r);
    }
    else if (n && roots[i].whole)
    {
      // The process which gets the most leaves of the root, then the lowest.
      std::map<int, int> votes;
      std::for_each(base->begin() + e, base->begin() + e + n, [&](int r) { votes[r]++; });
      const int r = std::max_element(votes.begin(), votes.end(),
                                     [](const auto &a, const auto &b)
                                     { return a.second < b.second; })
                        ->first;
      std::fill_n(part.begin() + e, n, r);
    }
    else if (n)
    {
      std::copy_n(base->begin() + e, n, part.begin() + e);
    }
  }
  std::vector<long long> load(P, 0);
  for (int r : part)
  {
    load[r]++;
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
  else if (base)
  {
    mesh.Rebalance(*base);
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

  // 2. The vertices, edges and faces on the leaves of lower processes, which belong to them
  //    (flags 1, 2 and 4), and the file numbers of the own vertices, top-level vertices
  //    first, in the order in which the own leaves reach them.
  std::vector<char> skip(std::max(nodes.NumIds(), faces.NumIds()), 0);
  auto EdgeId = [&](const auto &el, int j)
  {
    const auto &ev = NC::Geom(el.Geom()).edges[j];
    return nodes.FindId(el.node[ev[0]], el.node[ev[1]]);
  };
  auto FaceId = [&](const auto &el, int j)
  {
    const auto &fv = NC::Geom(el.Geom()).faces[j];
    return faces.FindId(el.node[fv[0]], el.node[fv[1]], el.node[fv[2]], el.node[fv[3]]);
  };
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
      skip[EdgeId(el, j)] |= 2;
    }
    for (int j = 0; dim == 3 && j < gi.nf; j++)
    {
      skip[FaceId(el, j)] |= 4;
    }
  }
  std::vector<long long> vid(nodes.NumIds(), -1);  // the number of a vertex in the file
  std::vector<int> own[2];                         // own top-level and other vertices
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
  for (int k = 0; k < 2; k++)
  {
    for (int n : own[k])
    {
      vid[n] += (k ? vtot[0] : 0) + voff[k];
    }
  }

  // 3. The numbers of the vertices of other processes, from their owner, the leader of
  //    their group. Every vertex printed here is a vertex of the mesh: the corners of the
  //    own leaves, and the parents of the own vertices, which are corners of an own
  //    element since the roots are whole.
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
      ok = ok && mesh.gtopo.IAmMaster(group[nodes[n].vert_index]);
      number[nodes[n].vert_index] = static_cast<int>(vid[n]);
    }
  }
  mfem::GroupCommunicator gc(mesh.gtopo);
  gc.Create(group);
  gc.Bcast(number);
  auto Number = [&](int n)
  {
    const int v = nodes[n].vert_index;
    ok = ok && v >= 0 && v < mesh.GetNV() && number[v] >= 0;
    vid[n] = ok ? number[v] : -1;
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
  const auto *fec = fes ? fes->FEColl() : nullptr;
  const int nvd = fec ? fec->DofForGeometry(Geometry::POINT) : 0;
  const int ned = fec ? fec->DofForGeometry(Geometry::SEGMENT) : 0;
  const mfem::NCMesh::NCList *face_list = nullptr;
  if (counts && dim == 3)
  {
    mesh.pncmesh->GetEdgeList();  // ParNCMesh::BuildFaceList may need the edge owners
    face_list = &mesh.pncmesh->GetFaceList();
  }
  std::array<std::vector<int>, 5> dofs;
  std::array<long long, 2 * Geometry::NumGeom + 1> counted{};  // cells, true faces, bdr
  mfem::Array<int> edofs;
  for (int e = 0; e < ne; e++)
  {
    const auto &el = E[leaf[e]];
    const auto &gi = NC::Geom(el.Geom());
    counted[el.Geom()]++;
    if (fes)
    {
      fes->GetElementDofs(e, edofs);
    }
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
      const int m = EdgeId(el, j);
      Add(dofs[2], ned, !(skip[m] & 2),
          vid[el.node[gi.edges[j][0]]] > vid[el.node[gi.edges[j][1]]]);
      skip[m] |= 2;
    }
    for (int j = 0; j < gi.nf; j++)  // the boundary as NCMesh::PrintBoundary, faces in 3D
    {
      const int f = FaceId(el, j);
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
        Add(dofs[3], fec ? fec->DofForGeometry(fgeom) : 0, !(skip[f] & 4));
        counted[Geometry::NumGeom + fgeom] += !(skip[f] & 4) && face_list &&
                                              face_list->GetMeshIdType(faces[f].index) !=
                                                  mfem::NCMesh::NCList::MeshIdType::SLAVE;
        skip[f] |= 4;
      }
    }
    Add(dofs[4], edofs.Size() - k, true);
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

}  // namespace palace::mesh
