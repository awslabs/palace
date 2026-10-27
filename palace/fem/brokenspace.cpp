// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "brokenspace.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include "fem/fespace.hpp"
#include "utils/communication.hpp"

namespace palace::fem
{
namespace
{

using Index = HYPRE_BigInt;
int Unsigned(int d)
{
  return d >= 0 ? d : -1 - d;
}

// Rendezvous records, not replicated mesh data. All counts are checked before MPI's
// int-sized Alltoallv interface is used. Empty ranks take part in every exchange.
template <typename T>
std::vector<std::vector<T>> Exchange(MPI_Comm comm, const std::vector<std::vector<T>> &send)
{
  const int np = Mpi::Size(comm);
  std::vector<int> sc(np), rc(np), sd(np), rd(np);
  std::size_t total = 0;
  for (int p = 0; p < np; p++)
  {
    total += sizeof(T) * send[p].size();
  }
  int valid = total <= static_cast<std::size_t>(std::numeric_limits<int>::max());
  Mpi::GlobalMin(1, &valid, comm);
  MFEM_VERIFY(valid, "Broken-space exchange exceeds MPI's count range!");
  for (int p = 0, pos = 0; p < np; p++)
  {
    sd[p] = pos;
    sc[p] = sizeof(T) * send[p].size();
    pos += sc[p];
  }
  MPI_Alltoall(sc.data(), 1, MPI_INT, rc.data(), 1, MPI_INT, comm);
  total = std::accumulate(rc.begin(), rc.end(), std::size_t(0));
  valid = total <= static_cast<std::size_t>(std::numeric_limits<int>::max());
  Mpi::GlobalMin(1, &valid, comm);
  MFEM_VERIFY(valid, "Broken-space receive exceeds MPI's count range!");
  for (int p = 1; p < np; p++)
  {
    rd[p] = rd[p - 1] + rc[p - 1];
  }
  std::vector<unsigned char> sb(std::accumulate(sc.begin(), sc.end(), 0)), rb(total);
  for (int p = 0; p < np; p++)
  {
    if (sc[p])
    {
      std::memcpy(sb.data() + sd[p], send[p].data(), sc[p]);
    }
  }
  MPI_Alltoallv(sb.data(), sc.data(), sd.data(), MPI_BYTE, rb.data(), rc.data(), rd.data(),
                MPI_BYTE, comm);
  std::vector<std::vector<T>> recv(np);
  for (int p = 0; p < np; p++)
  {
    recv[p].resize(rc[p] / sizeof(T));
    if (rc[p])
    {
      std::memcpy(recv[p].data(), rb.data() + rd[p], rc[p]);
    }
  }
  return recv;
}

struct Entry
{
  Index col;
  double value;
};
using Row = std::vector<Entry>;

std::vector<Row> ReadRows(const mfem::HypreParMatrix &P)
{
  P.HostRead();
  mfem::SparseMatrix diag, offd;
  Index *cmap;
  P.GetDiag(diag);
  P.GetOffd(offd, cmap);
  const auto begin = hypre_ParCSRMatrixFirstColDiag((hypre_ParCSRMatrix *)P);
  std::vector<Row> rows(P.Height());
  for (int l = 0; l < P.Height(); l++)
  {
    for (int k = diag.GetI()[l]; k < diag.GetI()[l + 1]; k++)
    {
      if (diag.GetData()[k] != 0.0)
      {
        rows[l].push_back({begin + diag.GetJ()[k], diag.GetData()[k]});
      }
    }
    for (int k = offd.GetI()[l]; k < offd.GetI()[l + 1]; k++)
    {
      if (offd.GetData()[k] != 0.0)
      {
        rows[l].push_back({cmap[offd.GetJ()[k]], offd.GetData()[k]});
      }
    }
  }
  P.HypreRead();
  return rows;
}

// Exact names for current NC entities. Node parent keys, not coordinates or
// ancestor-element sides, identify vertices across rank-local NC numbering.
// Serialize the ancestry DAG once per requested node: no exponential tree expansion.
class NodeNames
{
  const mfem::NCMesh &nc;
  std::map<std::pair<int, int>, int> comparison;
  std::map<int, std::vector<std::int64_t>> names;
  int Compare(int a, int b)
  {
    if (a == b)
    {
      return 0;
    }
    const auto key = std::make_pair(a, b);
    auto found = comparison.find(key);
    if (found != comparison.end())
    {
      return found->second;
    }
    const auto &x = nc.GetNode(a), &y = nc.GetNode(b);
    const bool xr = x.p1 == x.p2, yr = y.p1 == y.p2;
    int result;
    if (xr != yr)
    {
      result = xr ? -1 : 1;
    }
    else if (xr)
    {
      result = (x.p1 > y.p1) - (x.p1 < y.p1);
    }
    else
    {
      int x1 = x.p1, x2 = x.p2, y1 = y.p1, y2 = y.p2;
      if (Compare(x1, x2) > 0)
      {
        std::swap(x1, x2);
      }
      if (Compare(y1, y2) > 0)
      {
        std::swap(y1, y2);
      }
      result = Compare(x1, y1);
      if (!result)
      {
        result = Compare(x2, y2);
      }
    }
    comparison[key] = result;
    comparison[{b, a}] = -result;
    return result;
  }

public:
  explicit NodeNames(const mfem::NCMesh &nc) : nc(nc) {}
  const std::vector<std::int64_t> &Name(int node)
  {
    auto found = names.find(node);
    if (found != names.end())
    {
      return found->second;
    }
    std::set<int> seen;
    std::vector<int> pending{node};
    for (std::size_t i = 0; i < pending.size(); i++)
    {
      const int id = pending[i];
      if (!seen.insert(id).second)
      {
        continue;
      }
      const auto &n = nc.GetNode(id);
      if (n.p1 != n.p2)
      {
        pending.push_back(n.p1);
        pending.push_back(n.p2);
      }
    }
    std::vector<int> ordered(seen.begin(), seen.end());
    std::sort(ordered.begin(), ordered.end(),
              [&](int a, int b) { return Compare(a, b) < 0; });
    std::map<int, int> index;
    for (std::size_t i = 0; i < ordered.size(); i++)
    {
      index[ordered[i]] = i;
    }
    std::vector<std::int64_t> name{static_cast<std::int64_t>(ordered.size()), index[node]};
    for (auto id : ordered)
    {
      const auto &n = nc.GetNode(id);
      if (n.p1 == n.p2)
      {
        name.insert(name.end(), {0, n.p1});
      }
      else
      {
        int a = n.p1, b = n.p2;
        if (Compare(a, b) > 0)
        {
          std::swap(a, b);
        }
        name.insert(name.end(), {1, index[a], index[b]});
      }
    }
    return names.emplace(node, std::move(name)).first->second;
  }
  std::vector<std::int64_t> Entity(int dim, const std::vector<int> &vertices)
  {
    std::vector<std::vector<std::int64_t>> sorted;
    for (auto v : vertices)
    {
      sorted.push_back(Name(v));
    }
    std::sort(sorted.begin(), sorted.end());
    std::vector<std::int64_t> key{dim, static_cast<std::int64_t>(sorted.size())};
    for (const auto &v : sorted)
    {
      key.push_back(v.size());
      key.insert(key.end(), v.begin(), v.end());
    }
    return key;
  }
};

void CoordinateEntityDofs(const mfem::ParFiniteElementSpace &space, int e, int dim,
                          int local, mfem::Array<int> &dofs)
{
  mfem::Array<int> entities, orientation;
  if (dim == 1)
  {
    space.GetParMesh()->GetElementEdges(e, entities, orientation);
    space.GetEdgeDofs(entities[local], dofs);
  }
  else
  {
    space.GetParMesh()->GetElementFaces(e, entities, orientation);
    space.GetFaceDofs(entities[local], dofs);
  }
}

// Include children touching a master only at an interior edge/vertex, even
// when RT has no read there or a high-order nodal basis vanishes there.
void CompleteCarrierStars(mfem::ParFiniteElementSpace &space, const std::vector<Row> &base,
                          std::vector<std::set<Index>> &support)
{
  auto &mesh = *space.GetParMesh();
  if (!mesh.Nonconforming())
  {
    return;
  }
  auto &nc = *mesh.ncmesh;
  const auto comm = mesh.GetComm();
  const int np = Mpi::Size(comm);
  NodeNames names(nc);
  std::vector<std::vector<int>> vertex_star(mesh.GetNV());
  std::map<std::pair<int, int>, std::vector<int>> edge_star;
  std::vector<int> vertex_node(mesh.GetNV(), -1);
  mfem::Array<int> vertices, edges, orientation, dofs;
  for (const auto &v : nc.GetVertexList().conforming)
  {
    if (v.index >= 0 && v.index < mesh.GetNV())
    {
      vertex_node[v.index] = nc.GetElement(v.element).node[v.local];
    }
  }
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementVertices(e, vertices);
    for (auto v : vertices)
    {
      vertex_star[v].push_back(e);
    }
    mesh.GetElementEdges(e, edges, orientation);
    for (auto edge : edges)
    {
      mfem::Array<int> v;
      mesh.GetEdgeVertices(edge, v);
      edge_star[{std::min(v[0], v[1]), std::max(v[0], v[1])}].push_back(e);
    }
  }
  auto Home = [&](const auto &key)
  {
    std::uint64_t h = 14695981039346656037ULL;
    for (auto v : key)
    {
      h = (h ^ static_cast<std::uint64_t>(v)) * 1099511628211ULL;
    }
    return h % np;  // Hashing routes requests only; full keys decide identity.
  };
  for (int dim = 1; dim < mesh.Dimension(); dim++)
  {
    const auto &list = dim == 1 ? nc.GetEdgeList() : nc.GetFaceList();
    std::vector<std::set<int>> incidence(list.masters.Size());
    for (int i = 0; i < list.masters.Size(); i++)
    {
      const auto &master = list.masters[i];
      auto Interior = [&](double x, double y)
      {
        if (master.Geom() == mfem::Geometry::SEGMENT)
        {
          return x > 0.0 && x < 1.0;
        }
        if (master.Geom() == mfem::Geometry::TRIANGLE)
        {
          return x > 0.0 && y > 0.0 && x + y < 1.0;
        }
        return x > 0.0 && y > 0.0 && x < 1.0 && y < 1.0;
      };
      auto Add = [&](const auto &elements)
      { incidence[i].insert(elements.begin(), elements.end()); };
      for (int k = master.slaves_begin; k < master.slaves_end; k++)
      {
        const auto &slave = list.slaves[k];
        mfem::DenseMatrix points;
        list.OrientedPointMatrix(slave, points);
        std::vector<int> v;
        if (dim == 1 || slave.index < 0)
        {
          const int edge = Unsigned(slave.index);
          if (edge >= mesh.GetNEdges())
          {
            continue;
          }
          mesh.GetEdgeVertices(edge, vertices);
          v = {std::min(vertices[0], vertices[1]), std::max(vertices[0], vertices[1])};
        }
        else
        {
          const auto &leaf = nc.GetElement(slave.element);
          if (leaf.rank != Mpi::Rank(comm) || leaf.index >= mesh.GetNE())
          {
            continue;
          }
          mesh.GetElementVertices(leaf.index, vertices);
          const int *fv = mesh.GetElement(leaf.index)->GetFaceVertices(slave.local);
          for (int j = 0; j < mfem::Geometry::NumVerts[slave.Geom()]; j++)
          {
            v.push_back(vertices[fv[j]]);
          }
        }
        for (int j = 0; j < static_cast<int>(v.size()); j++)
        {
          if (Interior(points(0, j), dim == 1 ? 0.0 : points(1, j)))
          {
            Add(vertex_star[v[j]]);
          }
        }
        const int nedge = v.size() == 2 ? 1 : v.size();
        for (int j = 0; j < nedge; j++)
        {
          const int next = (j + 1) % v.size();
          const double x = 0.5 * (points(0, j) + points(0, next));
          const double y = dim == 1 ? 0.0 : 0.5 * (points(1, j) + points(1, next));
          if (Interior(x, y))
          {
            Add(edge_star.at({std::min(v[j], v[next]), std::max(v[j], v[next])}));
          }
        }
      }
    }
    std::vector<std::vector<std::int64_t>> send(np);
    std::set<std::vector<std::int64_t>> offered;
    auto Offer = [&](int e, int local, const std::vector<int> &nodes)
    {
      auto key = names.Entity(dim, nodes);
      if (!offered.insert(key).second)
      {
        return;
      }
      CoordinateEntityDofs(space, e, dim, local, dofs);
      std::set<Index> cols;
      for (auto d : dofs)
      {
        for (const auto &entry : base[Unsigned(d)])
        {
          cols.insert(entry.col);
        }
      }
      auto &packet = send[Home(key)];
      packet.insert(packet.end(), {0, -1, static_cast<std::int64_t>(key.size())});
      packet.insert(packet.end(), key.begin(), key.end());
      packet.push_back(cols.size());
      packet.insert(packet.end(), cols.begin(), cols.end());
    };
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      mesh.GetElementVertices(e, vertices);
      const auto *el = mesh.GetElement(e);
      const int count = dim == 1 ? el->GetNEdges() : el->GetNFaces();
      for (int k = 0; k < count; k++)
      {
        const int *local = dim == 1 ? el->GetEdgeVertices(k) : el->GetFaceVertices(k);
        const int n = dim == 1 ? 2 : el->GetNFaceVertices(k);
        std::vector<int> nodes;
        for (int j = 0; j < n; j++)
        {
          nodes.push_back(vertex_node[vertices[local[j]]]);
        }
        Offer(e, k, nodes);
      }
    }
    for (int i = 0; i < list.masters.Size(); i++)
    {
      if (incidence[i].empty())
      {
        continue;
      }
      const auto &m = list.masters[i];
      const auto &leaf = nc.GetElement(m.element);
      std::unique_ptr<mfem::Element> el;
      switch (leaf.Geom())
      {
        case mfem::Geometry::TRIANGLE:
          el = std::make_unique<mfem::Triangle>();
          break;
        case mfem::Geometry::SQUARE:
          el = std::make_unique<mfem::Quadrilateral>();
          break;
        case mfem::Geometry::TETRAHEDRON:
          el = std::make_unique<mfem::Tetrahedron>();
          break;
        case mfem::Geometry::CUBE:
          el = std::make_unique<mfem::Hexahedron>();
          break;
        case mfem::Geometry::PRISM:
          el = std::make_unique<mfem::Wedge>();
          break;
        case mfem::Geometry::PYRAMID:
          el = std::make_unique<mfem::Pyramid>();
          break;
        default:
          MFEM_ABORT("Unsupported NC entity geometry!");
      }
      const int *local =
          dim == 1 ? el->GetEdgeVertices(m.local) : el->GetFaceVertices(m.local);
      const int n = dim == 1 ? 2 : el->GetNFaceVertices(m.local);
      std::vector<int> nodes;
      for (int j = 0; j < n; j++)
      {
        nodes.push_back(leaf.node[local[j]]);
      }
      auto key = names.Entity(dim, nodes);
      auto &packet = send[Home(key)];
      packet.insert(packet.end(), {1, i, static_cast<std::int64_t>(key.size())});
      packet.insert(packet.end(), key.begin(), key.end());
      packet.push_back(0);
    }
    struct Query
    {
      int rank, id;
      std::vector<std::int64_t> key;
    };
    std::vector<Query> queries;
    std::map<std::vector<std::int64_t>, std::set<Index>> registry;
    const auto received = Exchange(comm, send);
    for (int p = 0; p < np; p++)
    {
      const auto &packet = received[p];
      for (std::size_t pos = 0; pos < packet.size();)
      {
        const int query = packet.at(pos++), id = packet.at(pos++),
                  length = packet.at(pos++);
        std::vector<std::int64_t> key(packet.begin() + pos, packet.begin() + pos + length);
        pos += length;
        const int n = packet.at(pos++);
        if (query)
        {
          queries.push_back({p, id, std::move(key)});
        }
        else
        {
          registry[key].insert(packet.begin() + pos, packet.begin() + pos + n);
        }
        pos += n;
      }
    }
    send.assign(np, {});
    int missing = 0;
    for (const auto &q : queries)
    {
      const auto found = registry.find(q.key);
      missing += found == registry.end();
      auto &packet = send[q.rank];
      packet.insert(packet.end(),
                    {q.id, found == registry.end()
                               ? 0
                               : static_cast<std::int64_t>(found->second.size())});
      if (found != registry.end())
      {
        packet.insert(packet.end(), found->second.begin(), found->second.end());
      }
    }
    for (const auto &packet : Exchange(comm, send))
    {
      for (std::size_t pos = 0; pos < packet.size();)
      {
        const int id = packet.at(pos++), n = packet.at(pos++);
        for (auto e : incidence[id])
        {
          support[e].insert(packet.begin() + pos, packet.begin() + pos + n);
        }
        pos += n;
      }
    }
    Mpi::GlobalSum(1, &missing, comm);
    MFEM_VERIFY(missing == 0,
                missing << " current NC master entities have no owned trace!");
  }
}

struct StarRecord
{
  Index t, a, b;
  int rank, cut;
};
struct Component
{
  Index t, label, column;
  int owner;
};
struct Read
{
  Index t, element, column;
};

}  // namespace

std::unique_ptr<BrokenDofLayout>
BuildBrokenDofLayout(mfem::ParFiniteElementSpace &space,
                     const std::vector<int> &break_attributes)
{
  auto &mesh = *space.GetParMesh();
  const auto comm = mesh.GetComm();
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm), ne = mesh.GetNE();
  MFEM_VERIFY(space.GetVDim() == 1 && !space.IsVariableOrder(),
              "Broken spaces require vdim = 1 and uniform polynomial order!");
  space.ExchangeFaceNbrData();
  auto base = ReadRows(*space.Dof_TrueDof_Matrix());
  const auto begin =
      hypre_ParCSRMatrixFirstColDiag((hypre_ParCSRMatrix *)*space.Dof_TrueDof_Matrix());
  std::vector<Index> starts(np + 1);
  MPI_Allgather(&begin, 1, HYPRE_MPI_BIG_INT, starts.data(), 1, HYPRE_MPI_BIG_INT, comm);
  starts[np] = space.GlobalTrueVSize();
  auto Owner = [&](Index t)
  { return std::upper_bound(starts.begin(), starts.end(), t) - starts.begin() - 1; };

  // Current element identity and marked local facets, exchanged in precisely MFEM's
  // face-neighbour ordering. A coarse marked facet marks each slave patch; a marked
  // slave marks only that patch. No coordinate matching or root ancestry is involved.
  mfem::L2_FECollection p0(0, mesh.Dimension());
  mfem::ParFiniteElementSpace header_space(&mesh, &p0, 2);
  mfem::ParGridFunction header(&header_space);
  Index offset = 0, local_ne = ne;
  MPI_Exscan(&local_ne, &offset, 1, HYPRE_MPI_BIG_INT, MPI_SUM, comm);
  if (rank == 0)
  {
    offset = 0;
  }
  // HYPRE_BigInt can be 32-bit; the double exactness bound must not shift that type.
  MFEM_VERIFY(std::int64_t(offset) + ne < (std::int64_t(1) << 53),
              "Element labels exceed exact double range!");
  std::vector<int> marked(ne, 0);
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (std::binary_search(break_attributes.begin(), break_attributes.end(),
                           mesh.GetBdrAttribute(be)))
    {
      int e, info;
      mesh.GetBdrElementAdjacentElement(be, e, info);
      marked[e] |= 1 << mesh.DecodeFaceInfoLocalIndex(info);
    }
  }
  mfem::Array<int> dofs, positions;
  mfem::DofTransformation dt;
  mfem::Vector pair(2);
  for (int e = 0; e < ne; e++)
  {
    header_space.GetElementVDofs(e, dofs);
    pair(0) = offset + e;
    pair(1) = marked[e];
    header.SetSubVector(dofs, pair);
  }
  header.ExchangeFaceNbrData();
  std::vector<Index> nbr_id(mesh.GetNFaceNeighborElements());
  std::vector<int> nbr_mark(nbr_id.size());
  for (int e = 0; e < static_cast<int>(nbr_id.size()); e++)
  {
    header_space.GetFaceNbrElementVDofs(e, dofs, dt);
    header.FaceNbrData().GetSubVector(dofs, pair);
    nbr_id[e] = static_cast<Index>(pair(0));
    nbr_mark[e] = static_cast<int>(pair(1));
  }

  std::vector<std::vector<StarRecord>> send(np);
  std::vector<std::set<Index>> element_read(ne);
  std::vector<std::vector<std::set<Index>>> element_trace(ne);
  for (int e = 0; e < ne; e++)
  {
    space.GetElementDofs(e, dofs, dt);
    auto &read = element_read[e];
    for (auto d : dofs)
    {
      for (const auto &entry : base[Unsigned(d)])
      {
        read.insert(entry.col);
      }
    }
    const auto geom = mesh.GetElementGeometry(e);
    const int nf = mesh.Dimension() == 3 ? mesh.GetElement(e)->GetNFaces()
                                         : mfem::Geometry::NumEdges[geom];
    element_trace[e].resize(nf);
    for (int f = 0; f < nf; f++)
    {
      CoordinateEntityDofs(space, e, mesh.Dimension() - 1, f, positions);
      for (auto d : positions)
      {
        for (const auto &entry : base[Unsigned(d)])
        {
          element_trace[e][f].insert(entry.col);
        }
      }
    }
  }
  const auto actual_read = element_read;
  CompleteCarrierStars(space, base, element_read);
  for (int e = 0; e < ne; e++)
  {
    for (auto t : element_read[e])
    {
      send[Owner(t)].push_back({t, offset + e, -1, actual_read[e].count(t) ? rank : np, 0});
    }
  }
  // A master face coefficient is read by every child touching that face, even
  // when its normal/tangential trace vanishes on the face between two children.
  // Join common entity reads, NOT just nonzeros of that intermediate face's P.
  // The latter falsely gives RT one master copy per child and breaks curl grad.
  std::vector<std::vector<Index>> support(np);
  for (int fn = 0; fn < mesh.GetNFaceNeighbors(); fn++)
  {
    auto &packet = support[mesh.GetFaceNbrRank(fn)];
    for (int k = 0; k < mesh.send_face_nbr_elements.RowSize(fn); k++)
    {
      const int e = mesh.send_face_nbr_elements.GetRow(fn)[k];
      const auto &read = element_read[e];
      packet.push_back(read.size());
      packet.insert(packet.end(), read.begin(), read.end());
      packet.push_back(element_trace[e].size());
      for (const auto &trace : element_trace[e])
      {
        packet.push_back(trace.size());
        packet.insert(packet.end(), trace.begin(), trace.end());
      }
    }
  }
  auto ghost_support = Exchange(comm, support);
  std::vector<std::set<Index>> nbr_read(mesh.GetNFaceNeighborElements());
  std::vector<std::vector<std::set<Index>>> nbr_trace(nbr_read.size());
  for (int fn = 0; fn < mesh.GetNFaceNeighbors(); fn++)
  {
    const auto &packet = ghost_support[mesh.GetFaceNbrRank(fn)];
    std::size_t pos = 0;
    for (int e = mesh.face_nbr_elements_offset[fn];
         e < mesh.face_nbr_elements_offset[fn + 1]; e++)
    {
      const int n = packet.at(pos++);
      nbr_read[e].insert(packet.begin() + pos, packet.begin() + pos + n);
      pos += n;
      const int nf = packet.at(pos++);
      nbr_trace[e].resize(nf);
      for (auto &trace : nbr_trace[e])
      {
        const int size = packet.at(pos++);
        trace.insert(packet.begin() + pos, packet.begin() + pos + size);
        pos += size;
      }
    }
    MFEM_VERIFY(pos == packet.size(), "Invalid neighbour entity support packet!");
  }
  for (int f = 0; f < mesh.GetNumFacesWithGhost(); f++)
  {
    const auto face = mesh.GetFaceInformation(f);
    if (!face.IsLocal() && !face.IsShared())
    {
      continue;  // Exterior and master faces: connectivity is on slave patches.
    }
    const int e = face.element[0].index, other = face.element[1].index;
    const bool remote = face.IsShared();
    const Index a = offset + e, b = remote ? nbr_id[other] : offset + other;
    const int cut =
        ((marked[e] >> face.element[0].local_face_id) & 1) |
        (((remote ? nbr_mark[other] : marked[other]) >> face.element[1].local_face_id) & 1);
    const auto &left = element_read[e],
               &right = remote ? nbr_read[other] : element_read[other];
    std::vector<Index> common;
    std::set_intersection(left.begin(), left.end(), right.begin(), right.end(),
                          std::back_inserter(common));
    const auto &trace1 = element_trace[e][face.element[0].local_face_id];
    const auto &trace2 =
        (remote ? nbr_trace[other] : element_trace[other])[face.element[1].local_face_id];
    for (auto t : common)
    {
      // Only closure coefficients of the marked patch acquire copies. Other
      // shared reads may arise from hanging constraints on a different facet.
      const int touched = cut && trace1.count(t) && trace2.count(t);
      // A marked face blocks every common carrier incidence, even if this
      // particular basis coefficient vanishes in that patch's actual trace.
      send[Owner(t)].push_back({t, a, b, rank, cut ? 1 + touched : 0});
    }
  }

  std::map<Index, std::vector<StarRecord>> stars;
  for (const auto &part : Exchange(comm, send))
  {
    for (const auto &record : part)
    {
      stars[record.t].push_back(record);
    }
  }
  // Incident reads include slave edge/vertex rows in the master entity's chart.
  // Components are computed at the old owner, not on an incomplete
  // one-layer ghost star. Unmarked entities retain one copy even on non-manifolds.
  std::map<std::pair<Index, Index>, Component> components;
  std::vector<std::vector<Read>> reads(np);
  std::vector<std::vector<Component>> proposed(np);
  int missing = 0;
  for (const auto &[t, records] : stars)
  {
    std::map<Index, int> readers;
    bool touched = false;
    for (const auto &r : records)
    {
      if (r.b < 0)
      {
        readers.emplace(r.a, r.rank);
      }
      touched |= r.cut == 2;
    }
    std::map<Index, Index> parent;
    for (const auto &[e, p] : readers)
    {
      parent[e] = e;
    }
    auto Find = [&](Index e)
    {
      while (parent.at(e) != e)
      {
        e = parent[e] = parent.at(parent.at(e));
      }
      return e;
    };
    for (const auto &r : records)
    {
      if (r.b >= 0 && !r.cut)
      {
        if (!readers.count(r.a) || !readers.count(r.b))
        {
          missing++;
          continue;
        }
        const auto a = Find(r.a), b = Find(r.b);
        parent[std::max(a, b)] = std::min(a, b);
      }
    }
    std::map<Index, int> owners;
    for (const auto &[e, p] : readers)
    {
      const Index label = touched ? Find(e) : readers.begin()->first;
      auto [it, inserted] = owners.emplace(label, p);
      it->second = std::min(it->second, p);
      if (p < np)
      {
        reads[p].push_back({t, e, label});
      }
    }
    for (const auto &[label, p] : owners)
    {
      if (p == np)
      {
        continue;
      }  // A virtual-only component is not an unknown.
      // Preserve conforming-first ordering where ownership permits it. This flag
      // is temporary; numbering below replaces it with the actual global column.
      Component c{t, label, label == owners.begin()->first ? 0 : 1, p};
      components.emplace(std::make_pair(t, label), c);
      proposed[p].push_back(c);
    }
  }
  Mpi::GlobalSum(1, &missing, comm);
  MFEM_VERIFY(missing == 0,
              "Face trace coefficient missing from " << missing << " incident stars!");

  std::vector<Component> owned;
  for (const auto &part : Exchange(comm, proposed))
  {
    owned.insert(owned.end(), part.begin(), part.end());
  }
  std::sort(owned.begin(), owned.end(),
            [&](const auto &a, const auto &b)
            {
              return std::tuple(Owner(a.t), a.column, a.t, a.label) <
                     std::tuple(Owner(b.t), b.column, b.t, b.label);
            });
  auto out = std::make_unique<BrokenDofLayout>();
  Index nt = owned.size();
  mfem::Array<Index> *col_partition[] = {&out->col_offsets};
  mesh.GenerateOffsets(1, &nt, col_partition);
  const auto col_begin = out->col_offsets[HYPRE_AssumedPartitionCheck() ? 0 : rank];
  proposed.assign(np, {});
  for (std::size_t i = 0; i < owned.size(); i++)
  {
    owned[i].column = col_begin + i;
    proposed[Owner(owned[i].t)].push_back(owned[i]);
  }
  for (const auto &part : Exchange(comm, proposed))
  {
    for (const auto &c : part)
    {
      components.at({c.t, c.label}).column = c.column;
    }
  }
  for (auto &part : reads)
  {
    for (auto &r : part)
    {
      r.column = components.at({r.t, r.column}).column;
    }
  }
  std::map<std::pair<Index, Index>, Index> columns;
  for (const auto &part : Exchange(comm, reads))
  {
    for (const auto &r : part)
    {
      columns.emplace(std::make_pair(r.t, r.element), r.column);
    }
  }

  // Reuse a coordinate row for its first actual read; append only distinct side reads.
  // Dormant MFEM master rows have no element read and remain zero, rather than inventing
  // an ownership vote or an interpolation anchor which element-local I cannot write.
  std::vector<Row> rows(base.size());
  out->row_origin.resize(base.size());
  std::iota(out->row_origin.begin(), out->row_origin.end(), 0);
  std::vector<char> used(base.size(), 0);
  std::map<std::pair<int, std::vector<Index>>, int> row_map;
  out->offsets.push_back(0);
  for (int e = 0; e < ne; e++)
  {
    space.GetElementDofs(e, dofs, dt);
    for (auto d : dofs)
    {
      const int l = Unsigned(d);
      Row row = base[l];
      std::vector<Index> key;
      for (auto &entry : row)
      {
        entry.col = columns.at({entry.col, offset + e});
        key.push_back(entry.col);
      }
      auto [it, inserted] = row_map.emplace(std::make_pair(l, key), -1);
      if (inserted)
      {
        if (!used[l])
        {
          it->second = l;
          rows[l] = std::move(row);
          used[l] = 1;
        }
        else
        {
          it->second = rows.size();
          rows.push_back(std::move(row));
          out->row_origin.push_back(l);
        }
      }
      out->dofs.push_back(d >= 0 ? it->second : -1 - it->second);
    }
    out->offsets.push_back(out->dofs.size());
  }

  // Exchange the neighbour's actual broken rows, not a global side flip. MFEM owns
  // pairing, subface incidence, element order, signs and DOF transformations.
  std::vector<std::vector<std::int64_t>> packets(np);
  for (int fn = 0; fn < mesh.GetNFaceNeighbors(); fn++)
  {
    auto &packet = packets[mesh.GetFaceNbrRank(fn)];
    for (int k = 0; k < mesh.send_face_nbr_elements.RowSize(fn); k++)
    {
      const int e = mesh.send_face_nbr_elements.GetRow(fn)[k];
      packet.push_back(out->offsets[e + 1] - out->offsets[e]);
      for (int j = out->offsets[e]; j < out->offsets[e + 1]; j++)
      {
        const int d = out->dofs[j], l = Unsigned(d);
        packet.insert(packet.end(), {d, static_cast<std::int64_t>(rows[l].size())});
        for (const auto &entry : rows[l])
        {
          packet.push_back(entry.col);
          packet.push_back(std::bit_cast<std::int64_t>(entry.value));
        }
      }
    }
  }
  const auto received = Exchange(comm, packets);
  std::map<std::pair<int, int>, int> ghost_rows;
  out->nbr_offsets.push_back(0);
  for (int fn = 0; fn < mesh.GetNFaceNeighbors(); fn++)
  {
    const int p = mesh.GetFaceNbrRank(fn);
    const auto &packet = received[p];
    std::size_t pos = 0;
    for (int e = mesh.face_nbr_elements_offset[fn];
         e < mesh.face_nbr_elements_offset[fn + 1]; e++)
    {
      const int n = packet.at(pos++);
      for (int j = 0; j < n; j++)
      {
        const int d = packet.at(pos++), nnz = packet.at(pos++);
        auto [it, inserted] =
            ghost_rows.emplace(std::make_pair(p, Unsigned(d)), rows.size());
        Row row;
        for (int k = 0; k < nnz; k++)
        {
          const Index col = packet.at(pos++);
          const double value = std::bit_cast<double>(packet.at(pos++));
          row.push_back({col, value});
        }
        if (inserted)
        {
          out->neighbour_rows.push_back({it->second, p, Unsigned(d)});
          rows.push_back(std::move(row));
          out->row_origin.push_back(-1);
        }
        out->nbr_dofs.push_back(d >= 0 ? it->second : -1 - it->second);
      }
      out->nbr_offsets.push_back(out->nbr_dofs.size());
    }
    MFEM_VERIFY(pos == packet.size(), "Invalid broken face-neighbour packet!");
  }

  Index nl = rows.size();
  mfem::Array<Index> *row_partition[] = {&out->row_offsets};
  mesh.GenerateOffsets(1, &nl, row_partition);
  std::vector<int> I(1, 0);
  std::vector<Index> J;
  std::vector<double> values;
  for (const auto &row : rows)
  {
    for (const auto &entry : row)
    {
      J.push_back(entry.col);
      values.push_back(entry.value);
    }
    I.push_back(J.size());
  }
  const int last = HYPRE_AssumedPartitionCheck() ? 2 : np;
  out->P = std::make_unique<mfem::HypreParMatrix>(
      comm, static_cast<int>(nl), out->row_offsets[last], out->col_offsets[last], I.data(),
      J.data(), values.data(), out->row_offsets.GetData(), out->col_offsets.GetData());
  return out;
}

namespace
{
std::vector<int> TracePositions(const mfem::FiniteElementCollection &fec,
                                const mfem::Element &el, int order, int dim, int local)
{
  const int nv = el.GetNVertices(), ne = el.GetNEdges();
  const int vd = fec.GetNumDof(mfem::Geometry::POINT, order);
  const int ed = fec.GetNumDof(mfem::Geometry::SEGMENT, order);
  const int *vertices = dim == 1 ? el.GetEdgeVertices(local) : el.GetFaceVertices(local);
  const int n = dim == 1 ? 2 : el.GetNFaceVertices(local);
  std::set<int> closure(vertices, vertices + n);
  std::vector<int> positions;
  for (auto v : closure)
  {
    for (int j = 0; j < vd; j++)
    {
      positions.push_back(v * vd + j);
    }
  }
  for (int e = 0; e < ne; e++)
  {
    const auto *ev = el.GetEdgeVertices(e);
    if ((dim == 1 && e == local) ||
        (dim == 2 && closure.count(ev[0]) && closure.count(ev[1])))
    {
      for (int j = 0; j < ed; j++)
      {
        positions.push_back(nv * vd + e * ed + j);
      }
    }
  }
  if (dim == 2)
  {
    int offset = nv * vd + ne * ed;
    for (int f = 0; f <= local; f++)
    {
      const auto geom =
          el.GetNFaceVertices(f) == 3 ? mfem::Geometry::TRIANGLE : mfem::Geometry::SQUARE;
      const int fd = fec.GetNumDof(geom, order);
      if (f == local)
      {
        for (int j = 0; j < fd; j++)
        {
          positions.push_back(offset + j);
        }
      }
      offset += fd;
    }
  }
  return positions;
}
}  // namespace

std::unique_ptr<mfem::SparseMatrix> BuildBrokenRestriction(const FiniteElementSpace &space)
{
  const auto &P = *space.Dof_TrueDof_Matrix();
  const auto rows = ReadRows(P);
  const auto begin = hypre_ParCSRMatrixFirstColDiag((hypre_ParCSRMatrix *)P);
  const int nt = P.Width();
  auto R = std::make_unique<mfem::SparseMatrix>(nt, P.Height());
  std::vector<char> selected(nt, 0), active(P.Height(), 0);
  mfem::Array<int> dofs, positions;
  mfem::DofTransformation dt;
  for (int e = 0; e < space.GetMesh().GetNE(); e++)
  {
    space.GetElementDofs(e, dofs, dt);
    for (auto d : dofs)
    {
      active[Unsigned(d)] = 1;
    }
  }
  auto Select = [&](int l)
  {
    if (rows[l].size() == 1)
    {
      const auto &entry = rows[l][0];
      const auto t = entry.col - begin;
      if (t >= 0 && t < nt && !selected[t])
      {
        R->Set(t, l, 1.0 / entry.value);
        selected[t] = 1;
      }
    }
  };
  for (int l = 0; l < P.Height(); l++)
  {
    if (active[l])
    {
      Select(l);
    }
  }

  // Hanging constraints and parallel ND orientations can leave a master coefficient
  // without a singleton P row, so even signed/scaled row selection cannot invert P.
  // A child's complete trace is injective. Its inverse includes the closure,
  // since face rows alone also read edge/vertex
  // coefficients. Parallel ND orientation blocks are recovered by the same rule.
  // The SVD is a setup-only local inverse, not an iterative/global solve or average.
  auto Recover = [&](std::vector<int> trace)
  {
    std::sort(trace.begin(), trace.end());
    trace.erase(std::unique(trace.begin(), trace.end()), trace.end());
    std::vector<Index> cols;
    for (auto l : trace)
    {
      for (const auto &entry : rows[l])
      {
        cols.push_back(entry.col);
      }
    }
    std::sort(cols.begin(), cols.end());
    cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
    const int m = trace.size(), n = cols.size();
    if (!n || m < n ||
        std::none_of(cols.begin(), cols.end(), [&](auto c)
                     { return c >= begin && c < begin + nt && !selected[c - begin]; }))
    {
      return;
    }
    mfem::DenseMatrix block(m, n);
    block = 0.0;
    for (int r = 0; r < m; r++)
    {
      for (const auto &entry : rows[trace[r]])
      {
        block(r, std::lower_bound(cols.begin(), cols.end(), entry.col) - cols.begin()) =
            entry.value;
      }
    }
    mfem::DenseMatrix original(block);
    mfem::DenseMatrixSVD svd(block, 'S', 'S');
    svd.Eval(block);
    const auto &s = svd.Singularvalues();
    if (s(n - 1) <= 1.0e-12 * s(0))
    {
      return;
    }
    const auto &U = svd.LeftSingularvectors(), &Vt = svd.RightSingularvectors();
    mfem::DenseMatrix inverse(n, m);
    inverse = 0.0;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < m; j++)
      {
        for (int k = 0; k < n; k++)
        {
          inverse(i, j) += Vt(k, i) * U(j, k) / s(k);
        }
      }
    }
    mfem::DenseMatrix identity(n);
    mfem::Mult(inverse, original, identity);
    double error = 0.0;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        error = std::max(error, std::abs(identity(i, j) - (i == j)));
      }
    }
    if (error > 1.0e-11)
    {
      return;
    }
    for (int i = 0; i < n; i++)
    {
      const auto t = cols[i] - begin;
      if (t >= 0 && t < nt && !selected[t])
      {
        for (int j = 0; j < m; j++)
        {
          R->Set(t, trace[j], inverse(i, j));
        }
        selected[t] = 1;
      }
    }
  };
  const auto &mesh = space.GetParMesh();
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    space.GetElementDofs(e, dofs, dt);
    const auto geom = mesh.GetElementGeometry(e);
    for (int dim = 1; dim < mesh.Dimension(); dim++)
    {
      const int count =
          dim == 1 ? mfem::Geometry::NumEdges[geom] : mesh.GetElement(e)->GetNFaces();
      for (int entity = 0; entity < count; entity++)
      {
        space.GetElementEntityDofs(e, dim, entity, positions);
        std::vector<int> trace;
        for (auto d : positions)
        {
          trace.push_back(Unsigned(d));
        }
        Recover(std::move(trace));
      }
    }
  }
  int missing = std::count(selected.begin(), selected.end(), 0);
  Mpi::GlobalSum(1, &missing, space.GetComm());
  if (missing)
  {
    // Lowest reading rank need not own a full slave face: a tet can touch it
    // only along a hanging edge. Its stored face-neighbour L-rows can supply
    // direct reads or a complete injective trace.
    for (int e = 0; e < mesh.GetNFaceNeighborElements(); e++)
    {
      space.GetFaceNbrElementDofs(e, dofs, dt);
      // A halo singleton can recover an owned coefficient even when its containing
      // trace also has unrelated columns which prevent inversion of the whole block.
      for (auto d : dofs)
      {
        Select(Unsigned(d));
      }
      const auto &el = *mesh.face_nbr_elements[e];
      for (int dim = 1; dim < mesh.Dimension(); dim++)
      {
        const int count = dim == 1 ? el.GetNEdges() : el.GetNFaces();
        for (int entity = 0; entity < count; entity++)
        {
          std::vector<int> trace;
          for (auto j : TracePositions(space.GetFEColl(), el, space.GetMaxElementOrder(),
                                       dim, entity))
          {
            trace.push_back(Unsigned(dofs[j]));
          }
          Recover(std::move(trace));
        }
      }
    }
  }
  missing = std::count(selected.begin(), selected.end(), 0);
  Mpi::GlobalSum(1, &missing, space.GetComm());
  MFEM_VERIFY(missing == 0,
              missing << " owned broken DOF(s) have no injective local trace!");
  R->Finalize();
  return R;
}

std::unique_ptr<mfem::HypreParMatrix> BuildTrueInterpolator(const FiniteElementSpace &trial,
                                                            const FiniteElementSpace &test,
                                                            const BrokenDofLayout &layout,
                                                            const mfem::SparseMatrix &local)
{
  const auto comm = test.GetComm();
  mfem::HypreParMatrix I(comm, test.GlobalVSize(), trial.GlobalVSize(),
                         test.GetDofOffsets(), trial.GetDofOffsets(),
                         const_cast<mfem::SparseMatrix *>(&local));
  std::unique_ptr<mfem::HypreParMatrix> IP(
      new mfem::HypreParMatrix(hypre_ParCSRMatMat(I, *trial.Dof_TrueDof_Matrix()), true));
  auto rows = ReadRows(*IP);
  const auto &mesh = test.GetParMesh();
  std::vector<std::vector<std::int64_t>> packets(Mpi::Size(comm));
  mfem::Array<int> dofs;
  mfem::DofTransformation dt;
  for (int fn = 0; fn < mesh.GetNFaceNeighbors(); fn++)
  {
    std::set<int> read;
    for (int k = 0; k < mesh.send_face_nbr_elements.RowSize(fn); k++)
    {
      test.GetElementDofs(mesh.send_face_nbr_elements.GetRow(fn)[k], dofs, dt);
      for (auto d : dofs)
      {
        read.insert(Unsigned(d));
      }
    }
    auto &packet = packets[mesh.GetFaceNbrRank(fn)];
    for (auto l : read)
    {
      packet.insert(packet.end(), {l, static_cast<std::int64_t>(rows[l].size())});
      for (const auto &entry : rows[l])
      {
        packet.push_back(entry.col);
        packet.push_back(std::bit_cast<std::int64_t>(entry.value));
      }
    }
  }
  const auto received = Exchange(comm, packets);
  std::map<std::pair<int, int>, Row> neighbours;
  for (int p = 0; p < Mpi::Size(comm); p++)
  {
    const auto &packet = received[p];
    for (std::size_t pos = 0; pos < packet.size();)
    {
      const int l = packet.at(pos++), n = packet.at(pos++);
      auto &row = neighbours[{p, l}];
      for (int k = 0; k < n; k++)
      {
        const Index col = packet.at(pos++);
        row.push_back({col, std::bit_cast<double>(packet.at(pos++))});
      }
    }
  }
  for (const auto &origin : layout.neighbour_rows)
  {
    rows[origin[0]] = neighbours.at({origin[1], origin[2]});
  }
  std::vector<int> indices(1, 0);
  std::vector<Index> columns;
  std::vector<double> values;
  for (const auto &row : rows)
  {
    for (const auto &entry : row)
    {
      columns.push_back(entry.col);
      values.push_back(entry.value);
    }
    indices.push_back(columns.size());
  }
  mfem::HypreParMatrix complete(
      comm, test.GetVSize(), test.GlobalVSize(), trial.GlobalTrueVSize(), indices.data(),
      columns.data(), values.data(), test.GetDofOffsets(), trial.GetTrueDofOffsets());
  mfem::HypreParMatrix R(comm, test.GlobalTrueVSize(), test.GlobalVSize(),
                         test.GetTrueDofOffsets(), test.GetDofOffsets(),
                         const_cast<mfem::SparseMatrix *>(test.GetRestrictionMatrix()));
  auto result =
      std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatMat(R, complete), true);
  hypre_ParCSRMatrixSetNumNonzeros(*result);
  return result;
}

}  // namespace palace::fem
