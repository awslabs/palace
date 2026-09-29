// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "brokenspace.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include "utils/communication.hpp"

namespace palace
{

namespace
{

enum EntityType : int
{
  VERTEX = 0,
  EDGE = 1,
  FACE = 2
};

// Build a communicator over the shared entities of the given type, with nslots consecutive
// L-dofs per shared entity, and return the local index of each shared entity in the order
// of the communicator L-dofs. In 2D, faces are the mesh edges. Shared entities are listed
// in the same order by all processes of a group, which makes the slot layout consistent.
std::unique_ptr<mfem::GroupCommunicator> BuildSharedEntityComm(mfem::ParMesh &mesh,
                                                               EntityType type, int nslots,
                                                               std::vector<int> &entities)
{
  const int dim = mesh.Dimension();
  const bool edges = (type == EDGE || (type == FACE && dim == 2));
  auto Count = [&](int g)
  {
    if (type == VERTEX)
    {
      return mesh.GroupNVertices(g);
    }
    if (edges)
    {
      return mesh.GroupNEdges(g);
    }
    return mesh.GroupNTriangles(g) + mesh.GroupNQuadrilaterals(g);
  };
  auto Entity = [&](int g, int i)
  {
    if (type == VERTEX)
    {
      return mesh.GroupVertex(g, i);
    }
    if (edges)
    {
      return mesh.GroupEdge(g, i);
    }
    const int ntri = mesh.GroupNTriangles(g);
    return (i < ntri) ? mesh.GroupTriangle(g, i) : mesh.GroupQuadrilateral(g, i - ntri);
  };

  auto gc = std::make_unique<mfem::GroupCommunicator>(mesh.gtopo);
  mfem::Table &group_ldof = gc->GroupLDofTable();
  const int ngroups = mesh.GetNGroups();
  group_ldof.MakeI(ngroups);
  for (int g = 1; g < ngroups; g++)
  {
    group_ldof.AddColumnsInRow(g, Count(g) * nslots);
  }
  group_ldof.MakeJ();
  entities.clear();
  for (int g = 1; g < ngroups; g++)
  {
    for (int i = 0; i < Count(g); i++)
    {
      const int p = static_cast<int>(entities.size());
      entities.push_back(Entity(g, i));
      for (int s = 0; s < nslots; s++)
      {
        group_ldof.AddConnection(g, p * nslots + s);
      }
    }
  }
  group_ldof.ShiftUpI();
  gc->Finalize();
  return gc;
}

// Reduce values over all processes sharing each entity, and broadcast the result.
template <typename T>
void ReduceShared(const mfem::GroupCommunicator &gc, const std::vector<int> &entities,
                  int nslots, std::vector<T> &buf,
                  void (*Op)(mfem::GroupCommunicator::OpData<T>))
{
  MFEM_ASSERT(buf.size() == entities.size() * nslots,
              "Unexpected buffer size for shared entity reduction!");
  gc.Reduce<T>(buf.data(), Op);
  gc.Bcast<T>(buf.data());
}

// Reduce a per-entity value over the processes sharing each entity (one slot per entity).
template <typename T>
void ReduceSharedEntities(const mfem::GroupCommunicator &gc,
                          const std::vector<int> &entities, std::vector<T> &values,
                          void (*Op)(mfem::GroupCommunicator::OpData<T>))
{
  std::vector<T> buf(entities.size());
  for (std::size_t p = 0; p < entities.size(); p++)
  {
    buf[p] = values[entities[p]];
  }
  ReduceShared(gc, entities, 1, buf, Op);
  for (std::size_t p = 0; p < entities.size(); p++)
  {
    values[entities[p]] = buf[p];
  }
}

}  // namespace

bool CrackSides::Any() const
{
  return std::any_of(copy.begin(), copy.end(), [](auto bits) { return bits != 0; });
}

namespace mesh
{

CrackSides ComputeCrackSides(const mfem::ParMesh &cmesh,
                             const mfem::Array<int> &bdr_attr_marker)
{
  // Many mfem::ParMesh accessors for shared entities are not const.
  auto &mesh = const_cast<mfem::ParMesh &>(cmesh);
  MPI_Comm comm = mesh.GetComm();
  const int dim = mesh.Dimension();
  MFEM_VERIFY(dim == 2 || dim == 3,
              "Interior boundary sides are only supported for 2D and 3D meshes!");
  const int ne = mesh.GetNE(), nf = mesh.GetNumFaces(), nv = mesh.GetNV();
  const int nedges = (dim == 3) ? mesh.GetNEdges() : nf;
  CrackSides sides;
  sides.copy.assign(ne, 0);
  sides.split.assign(ne, 0);
  MFEM_VERIFY(!HasHangingEntities(mesh),
              "Interior boundary sides can only be discovered on a mesh without hanging "
              "entities!");

  // Faces shared with other processes (interior faces, the neighbor element is remote).
  std::vector<int> shared_faces;
  auto gc_face = BuildSharedEntityComm(mesh, FACE, 1, shared_faces);
  std::vector<char> is_shared_face(nf, 0);
  for (auto f : shared_faces)
  {
    is_shared_face[f] = 1;
  }

  // Interior faces carrying one of the marked boundary attributes. The boundary element of
  // a face on a process boundary may only exist on one side, so agree on shared faces.
  std::vector<int> sheet_face(nf, 0);
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int attr = mesh.GetBdrAttribute(be);
    if (attr > 0 && attr <= bdr_attr_marker.Size() && bdr_attr_marker[attr - 1])
    {
      sheet_face[mesh.GetBdrElementFaceIndex(be)] = 1;
    }
  }
  ReduceSharedEntities(*gc_face, shared_faces, sheet_face,
                       mfem::GroupCommunicator::Max<int>);
  int num_sheet_faces = 0;
  for (int f = 0; f < nf; f++)
  {
    if (sheet_face[f])
    {
      int e1, e2;
      mesh.GetFaceElements(f, &e1, &e2);
      if (!is_shared_face[f] && e2 < 0)
      {
        sheet_face[f] = 0;  // Exterior boundary, nothing to break
      }
      else
      {
        num_sheet_faces++;
      }
    }
  }
  Mpi::GlobalSum(1, &num_sheet_faces, comm);
  if (num_sheet_faces == 0)
  {
    return sides;
  }

  // Vertices and edges of interior boundary faces (in 2D, the edges are the faces). An
  // entity on a process boundary may belong to interior boundary faces of another process
  // only, so agree on shared entities.
  std::vector<int> sheet_vert(nv, 0), sheet_edge(nedges, 0);
  mfem::Array<int> fv, fe, fo;
  for (int f = 0; f < nf; f++)
  {
    if (!sheet_face[f])
    {
      continue;
    }
    mesh.GetFaceVertices(f, fv);
    for (auto v : fv)
    {
      sheet_vert[v] = 1;
    }
    if (dim == 3)
    {
      mesh.GetFaceEdges(f, fe, fo);
      for (auto edge : fe)
      {
        sheet_edge[edge] = 1;
      }
    }
    else
    {
      sheet_edge[f] = 1;
    }
  }
  std::vector<int> shared_verts, shared_edges;
  auto gc_vert = BuildSharedEntityComm(mesh, VERTEX, 1, shared_verts);
  ReduceSharedEntities(*gc_vert, shared_verts, sheet_vert,
                       mfem::GroupCommunicator::Max<int>);
  std::unique_ptr<mfem::GroupCommunicator> gc_edge;
  if (dim == 3)
  {
    gc_edge = BuildSharedEntityComm(mesh, EDGE, 1, shared_edges);
    ReduceSharedEntities(*gc_edge, shared_edges, sheet_edge,
                         mfem::GroupCommunicator::Max<int>);
  }
  const auto &gc_edge_ref = (dim == 3) ? *gc_edge : *gc_face;
  const auto &shared_edges_ref = (dim == 3) ? shared_edges : shared_faces;

  // Collect (element, local entity) pairs for all interior boundary entities. Each such
  // pair is a node of a union-find structure: nodes of the same entity are joined when
  // their elements share a face containing the entity which is not an interior boundary
  // face.
  std::vector<char> touched(ne, 0);
  {
    std::unique_ptr<mfem::Table> vert_to_elem(mesh.GetVertexToElementTable());
    for (int v = 0; v < nv; v++)
    {
      if (sheet_vert[v])
      {
        const int *elems = vert_to_elem->GetRow(v);
        for (int k = 0; k < vert_to_elem->RowSize(v); k++)
        {
          touched[elems[k]] = 1;
        }
      }
    }
  }
  std::vector<int> node_offsets(ne + 1, 0), node_elem, node_bit, node_type, node_ent;
  {
    mfem::Array<int> ev, ee, eo, ef, efo;
    auto AddNode = [&](int e, int bit, int type, int ent)
    {
      node_elem.push_back(e);
      node_bit.push_back(bit);
      node_type.push_back(type);
      node_ent.push_back(ent);
    };
    for (int e = 0; e < ne; e++)
    {
      node_offsets[e] = static_cast<int>(node_elem.size());
      if (!touched[e])
      {
        continue;
      }
      const auto geom = mesh.GetElementGeometry(e);
      const int nv_e = mfem::Geometry::NumVerts[geom];
      const int ne_e = mfem::Geometry::NumEdges[geom];
      mesh.GetElementVertices(e, ev);
      for (int i = 0; i < ev.Size(); i++)
      {
        if (sheet_vert[ev[i]])
        {
          AddNode(e, i, VERTEX, ev[i]);
        }
      }
      mesh.GetElementEdges(e, ee, eo);
      for (int i = 0; i < ee.Size(); i++)
      {
        if (sheet_edge[ee[i]])
        {
          AddNode(e, nv_e + i, EDGE, ee[i]);
        }
      }
      if (dim == 3)
      {
        mesh.GetElementFaces(e, ef, efo);
        for (int i = 0; i < ef.Size(); i++)
        {
          if (sheet_face[ef[i]])
          {
            AddNode(e, nv_e + ne_e + i, FACE, ef[i]);
          }
        }
      }
      MFEM_VERIFY(nv_e + ne_e + ((dim == 3) ? ef.Size() : 0) <= 32,
                  "Too many local entities for an interior boundary side bitmask!");
    }
    node_offsets[ne] = static_cast<int>(node_elem.size());
  }
  auto FindNode = [&](int e, int type, int ent)
  {
    for (int k = node_offsets[e]; k < node_offsets[e + 1]; k++)
    {
      if (node_type[k] == type && node_ent[k] == ent)
      {
        return k;
      }
    }
    MFEM_ABORT("Unable to locate interior boundary entity in element " << e << "!");
    return -1;
  };

  // Global element numbers (as mfem::ParMesh::GetGlobalElementNum, which is collective on
  // first use while not all processes have interior boundary entities).
  long long elem_offset = 0;
  {
    long long ne_local = ne;
    MPI_Exscan(&ne_local, &elem_offset, 1, MPI_LONG_LONG, MPI_SUM, comm);
    if (Mpi::Rank(comm) == 0)
    {
      elem_offset = 0;  // MPI_Exscan leaves the result on the first process undefined
    }
  }
  auto GlobalElement = [elem_offset](int e)
  { return static_cast<double>(elem_offset + e); };

  // Union-find where each component carries the smallest global element number of its
  // members as its label.
  const int nn = static_cast<int>(node_elem.size());
  std::vector<int> parent(nn);
  std::iota(parent.begin(), parent.end(), 0);
  std::vector<double> label(nn);
  for (int k = 0; k < nn; k++)
  {
    label[k] = GlobalElement(node_elem[k]);
  }
  auto Find = [&parent](int k)
  {
    while (parent[k] != k)
    {
      parent[k] = parent[parent[k]];
      k = parent[k];
    }
    return k;
  };
  auto Union = [&](int a, int b)
  {
    a = Find(a);
    b = Find(b);
    if (a != b)
    {
      parent[a] = b;
      label[b] = std::min(label[a], label[b]);
    }
  };
  for (int f = 0; f < nf; f++)
  {
    if (sheet_face[f] || is_shared_face[f])
    {
      continue;
    }
    int e1, e2;
    mesh.GetFaceElements(f, &e1, &e2);
    if (e1 < 0 || e2 < 0 || !touched[e1] || !touched[e2])
    {
      continue;
    }
    mesh.GetFaceVertices(f, fv);
    for (auto v : fv)
    {
      if (sheet_vert[v])
      {
        Union(FindNode(e1, VERTEX, v), FindNode(e2, VERTEX, v));
      }
    }
    if (dim == 3)
    {
      mesh.GetFaceEdges(f, fe, fo);
      for (auto edge : fe)
      {
        if (sheet_edge[edge])
        {
          Union(FindNode(e1, EDGE, edge), FindNode(e2, EDGE, edge));
        }
      }
    }
  }

  // Join components across process boundaries: exchange the labels of the entities of each
  // shared face which is not an interior boundary face, until no label changes. The entity
  // slots of a face are ordered by global vertex numbers, identically on both processes.
  if (Mpi::Size(comm) > 1)
  {
    // Globally unique vertex identifiers: (rank, local index) of the owning process, which
    // is the master of the group of a shared vertex.
    std::vector<double> vert_id(nv);
    for (int v = 0; v < nv; v++)
    {
      vert_id[v] = std::ldexp(static_cast<double>(Mpi::Rank(comm)), 32) + v;
    }
    {
      std::vector<double> buf(shared_verts.size());
      for (std::size_t p = 0; p < shared_verts.size(); p++)
      {
        buf[p] = vert_id[shared_verts[p]];
      }
      gc_vert->Bcast<double>(buf.data());
      for (std::size_t p = 0; p < shared_verts.size(); p++)
      {
        vert_id[shared_verts[p]] = buf[p];
      }
    }
    auto GlobalVertex = [&](int v) { return vert_id[v]; };
    constexpr int max_face_verts = 4;
    const int nslots = (dim == 3) ? 2 * max_face_verts : 2;
    std::vector<int> faces;
    auto gc_slots = BuildSharedEntityComm(mesh, FACE, nslots, faces);
    std::vector<int> slot_node(faces.size() * nslots, -1);
    mfem::Array<int> edge_verts;
    for (std::size_t p = 0; p < faces.size(); p++)
    {
      const int f = faces[p];
      if (sheet_face[f])
      {
        continue;
      }
      int e1, e2;
      mesh.GetFaceElements(f, &e1, &e2);
      if (e1 < 0 || !touched[e1])
      {
        continue;
      }
      mesh.GetFaceVertices(f, fv);
      std::vector<std::pair<double, int>> verts;
      for (auto v : fv)
      {
        verts.emplace_back(GlobalVertex(v), v);
      }
      std::sort(verts.begin(), verts.end());
      for (std::size_t s = 0; s < verts.size(); s++)
      {
        if (sheet_vert[verts[s].second])
        {
          slot_node[p * nslots + s] = FindNode(e1, VERTEX, verts[s].second);
        }
      }
      if (dim == 3)
      {
        mesh.GetFaceEdges(f, fe, fo);
        std::vector<std::pair<std::array<double, 2>, int>> edges;
        for (auto edge : fe)
        {
          mesh.GetEdgeVertices(edge, edge_verts);
          const auto g0 = GlobalVertex(edge_verts[0]), g1 = GlobalVertex(edge_verts[1]);
          edges.push_back({{std::min(g0, g1), std::max(g0, g1)}, edge});
        }
        std::sort(edges.begin(), edges.end());
        for (std::size_t s = 0; s < edges.size(); s++)
        {
          if (sheet_edge[edges[s].second])
          {
            slot_node[p * nslots + max_face_verts + s] =
                FindNode(e1, EDGE, edges[s].second);
          }
        }
      }
    }
    std::vector<double> buf(slot_node.size());
    while (true)
    {
      for (std::size_t k = 0; k < slot_node.size(); k++)
      {
        buf[k] = (slot_node[k] >= 0) ? label[Find(slot_node[k])] : mfem::infinity();
      }
      ReduceShared(*gc_slots, faces, nslots, buf, mfem::GroupCommunicator::Min<double>);
      int changed = 0;
      for (std::size_t k = 0; k < slot_node.size(); k++)
      {
        if (slot_node[k] >= 0)
        {
          const int r = Find(slot_node[k]);
          if (buf[k] < label[r])
          {
            label[r] = buf[k];
            changed = 1;
          }
        }
      }
      Mpi::GlobalMax(1, &changed, comm);
      if (!changed)
      {
        break;
      }
    }
  }

  // Per-entity reductions over the nodes of an entity on all processes. The value function
  // returns NaN for nodes to skip.
  struct EntityInfo
  {
    int nent;
    const mfem::GroupCommunicator *gc;
    const std::vector<int> *shared;
  };
  const std::array<EntityInfo, 3> entity_info = {
      EntityInfo{nv, gc_vert.get(), &shared_verts},
      EntityInfo{nedges, &gc_edge_ref, &shared_edges_ref},
      EntityInfo{(dim == 3) ? nf : 0, gc_face.get(), &shared_faces}};
  auto EntityReduce = [&](auto &&value, bool use_max)
  {
    std::array<std::vector<double>, 3> out;
    for (int type = VERTEX; type <= FACE; type++)
    {
      const auto &info = entity_info[type];
      out[type].assign(info.nent, use_max ? -mfem::infinity() : mfem::infinity());
    }
    for (int k = 0; k < nn; k++)
    {
      const double v = value(k);
      if (!std::isnan(v))
      {
        auto &x = out[node_type[k]][node_ent[k]];
        x = use_max ? std::max(x, v) : std::min(x, v);
      }
    }
    for (int type = VERTEX; type <= FACE; type++)
    {
      const auto &info = entity_info[type];
      if (type != FACE || dim == 3)
      {
        ReduceSharedEntities(*info.gc, *info.shared, out[type],
                             use_max ? mfem::GroupCommunicator::Max<double>
                                     : mfem::GroupCommunicator::Min<double>);
      }
    }
    return out;
  };

  // An entity is split when the elements around it form more than one group.
  std::vector<double> comp(nn);
  for (int k = 0; k < nn; k++)
  {
    comp[k] = label[Find(k)];
  }
  const auto comp_min = EntityReduce([&](int k) { return comp[k]; }, false);
  const auto comp_max = EntityReduce([&](int k) { return comp[k]; }, true);
  auto Split = [&](int type, int ent)
  { return comp_min[type][ent] != comp_max[type][ent]; };

  // The choice of the group keeping the original DOFs has to be consistent between an
  // entity and the entities of its closure, since a nonconforming mesh constrains the DOFs
  // of hanging entities on an interior boundary face or edge to those of the face or edge
  // and its closure. Therefore, groups are identified with regions: sets of elements on the
  // same side of a connected interior boundary, joined through faces which contain a split
  // entity. The region with the smallest label keeps the original DOFs. Only at junctions
  // of three or more regions this can still be inconsistent.
  std::vector<int> rparent(ne);
  std::iota(rparent.begin(), rparent.end(), 0);
  std::vector<double> rlabel(ne);
  for (int e = 0; e < ne; e++)
  {
    rlabel[e] = GlobalElement(e);
  }
  auto RFind = [&rparent](int e)
  {
    while (rparent[e] != e)
    {
      rparent[e] = rparent[rparent[e]];
      e = rparent[e];
    }
    return e;
  };
  auto HasSplitEntity = [&](int f)
  {
    if (sheet_face[f])
    {
      return false;
    }
    mesh.GetFaceVertices(f, fv);
    for (auto v : fv)
    {
      if (sheet_vert[v] && Split(VERTEX, v))
      {
        return true;
      }
    }
    if (dim == 3)
    {
      mesh.GetFaceEdges(f, fe, fo);
      for (auto edge : fe)
      {
        if (sheet_edge[edge] && Split(EDGE, edge))
        {
          return true;
        }
      }
    }
    return false;
  };
  for (int f = 0; f < nf; f++)
  {
    if (is_shared_face[f] || !HasSplitEntity(f))
    {
      continue;
    }
    int e1, e2;
    mesh.GetFaceElements(f, &e1, &e2);
    if (e1 < 0 || e2 < 0)
    {
      continue;
    }
    const int r1 = RFind(e1), r2 = RFind(e2);
    if (r1 != r2)
    {
      rparent[r1] = r2;
      rlabel[r2] = std::min(rlabel[r1], rlabel[r2]);
    }
  }
  if (Mpi::Size(comm) > 1)
  {
    std::vector<int> slot_elem(shared_faces.size(), -1);
    for (std::size_t p = 0; p < shared_faces.size(); p++)
    {
      const int f = shared_faces[p];
      if (HasSplitEntity(f))
      {
        int e1, e2;
        mesh.GetFaceElements(f, &e1, &e2);
        slot_elem[p] = e1;
      }
    }
    std::vector<double> buf(shared_faces.size());
    while (true)
    {
      for (std::size_t p = 0; p < slot_elem.size(); p++)
      {
        buf[p] = (slot_elem[p] >= 0) ? rlabel[RFind(slot_elem[p])] : mfem::infinity();
      }
      ReduceShared(*gc_face, shared_faces, 1, buf, mfem::GroupCommunicator::Min<double>);
      int changed = 0;
      for (std::size_t p = 0; p < slot_elem.size(); p++)
      {
        if (slot_elem[p] >= 0)
        {
          const int r = RFind(slot_elem[p]);
          if (buf[p] < rlabel[r])
          {
            rlabel[r] = buf[p];
            changed = 1;
          }
        }
      }
      Mpi::GlobalMax(1, &changed, comm);
      if (!changed)
      {
        break;
      }
    }
  }
  std::vector<double> region(nn);
  for (int k = 0; k < nn; k++)
  {
    region[k] = rlabel[RFind(node_elem[k])];
  }

  // For each split entity, the base group is the one in the region with the smallest label
  // (falling back to the group with the smallest label if that region contains several of
  // its groups). Elements of all other groups read the copy.
  const auto region_min = EntityReduce([&](int k) { return region[k]; }, false);
  auto BaseCandidate = [&](int k)
  { return (region[k] == region_min[node_type[k]][node_ent[k]]) ? comp[k] : std::nan(""); };
  const auto base_min = EntityReduce(BaseCandidate, false);
  const auto base_max = EntityReduce(BaseCandidate, true);
  for (int k = 0; k < nn; k++)
  {
    const int type = node_type[k], ent = node_ent[k];
    if (!Split(type, ent))
    {
      continue;
    }
    sides.split[node_elem[k]] |= (std::uint32_t(1) << node_bit[k]);
    const double base = (base_min[type][ent] == base_max[type][ent]) ? base_min[type][ent]
                                                                     : comp_min[type][ent];
    if (comp[k] != base)
    {
      sides.copy[node_elem[k]] |= (std::uint32_t(1) << node_bit[k]);
    }
  }
  return sides;
}

bool HasHangingEntities(const mfem::ParMesh &mesh)
{
  int hanging = 0;
  if (mesh.Nonconforming())
  {
    // The lists include the ghost layer, so it is enough to check them locally. In 2D, the
    // faces are in the edge list (and the face list is empty).
    auto &ncmesh = *mesh.ncmesh;
    hanging = (ncmesh.GetEdgeList().slaves.Size() > 0 ||
               (mesh.Dimension() == 3 && ncmesh.GetFaceList().slaves.Size() > 0));
  }
  Mpi::GlobalMax(1, &hanging, mesh.GetComm());
  return hanging;
}

namespace
{

// Reference element topology for the parent of a refined element: vertex coordinates and
// the vertex lists of edges and faces, in MFEM's local ordering.
struct ReferenceTopology
{
  std::vector<std::array<double, 3>> verts;
  std::vector<std::vector<int>> edges, faces;

  // Bitmask of the closure of each local entity (vertices, edges, faces), in the local
  // entity numbering of interior boundary sides.
  std::vector<std::uint32_t> closure;
};

const ReferenceTopology &GetReferenceTopology(mfem::Geometry::Type geom)
{
  static std::array<std::unique_ptr<ReferenceTopology>, mfem::Geometry::NumGeom> topo;
  auto &t = topo[geom];
  if (!t)
  {
    t = std::make_unique<ReferenceTopology>();
    const mfem::IntegrationRule &ir = *mfem::Geometries.GetVertices(geom);
    for (int i = 0; i < ir.GetNPoints(); i++)
    {
      const auto &ip = ir.IntPoint(i);
      t->verts.push_back({ip.x, ip.y, ip.z});
    }
    // Elements of the mesh may come from a memory pool (mfem::Mesh::NewElement), so use
    // separately allocated ones.
    std::unique_ptr<mfem::Element> el;
    switch (geom)
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
        MFEM_ABORT("Unsupported element geometry for interior boundary sides!");
    }
    for (int j = 0; j < el->GetNEdges(); j++)
    {
      const int *ev = el->GetEdgeVertices(j);
      t->edges.push_back({ev[0], ev[1]});
    }
    if (mfem::Geometry::Dimension[geom] == 3)
    {
      for (int k = 0; k < el->GetNFaces(); k++)
      {
        const int *fv = el->GetFaceVertices(k);
        t->faces.emplace_back(fv, fv + el->GetNFaceVertices(k));
      }
    }
    const int nv = static_cast<int>(t->verts.size()),
              ne = static_cast<int>(t->edges.size());
    for (int i = 0; i < nv; i++)
    {
      t->closure.push_back(std::uint32_t(1) << i);
    }
    for (int j = 0; j < ne; j++)
    {
      t->closure.push_back((std::uint32_t(1) << (nv + j)) | t->closure[t->edges[j][0]] |
                           t->closure[t->edges[j][1]]);
    }
    for (std::size_t k = 0; k < t->faces.size(); k++)
    {
      const auto &fv = t->faces[k];
      std::uint32_t mask = std::uint32_t(1) << (nv + ne + k);
      for (int j = 0; j < ne; j++)
      {
        // The edges of a face are those with both vertices on the face.
        if (std::find(fv.begin(), fv.end(), t->edges[j][0]) != fv.end() &&
            std::find(fv.begin(), fv.end(), t->edges[j][1]) != fv.end())
        {
          mask |= t->closure[nv + j];
        }
      }
      t->closure.push_back(mask);
    }
  }
  return *t;
}

using Point = std::array<double, 3>;

bool SamePoint(const Point &p, const Point &q)
{
  constexpr double tol = 1.0e-10;
  return std::abs(p[0] - q[0]) < tol && std::abs(p[1] - q[1]) < tol &&
         std::abs(p[2] - q[2]) < tol;
}

Point Sub(const Point &p, const Point &q)
{
  return {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
}

Point Cross(const Point &a, const Point &b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

double Dot(const Point &a, const Point &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Point on the (closed) segment [a, b].
bool OnSegment(const Point &p, const Point &a, const Point &b)
{
  constexpr double tol = 1.0e-10;
  const auto ab = Sub(b, a), ap = Sub(p, a);
  const auto c = Cross(ab, ap);
  if (Dot(c, c) > tol * tol * Dot(ab, ab))
  {
    return false;
  }
  const double t = Dot(ap, ab) / Dot(ab, ab);
  return t > -tol && t < 1.0 + tol;
}

// Point on the plane of a (planar) reference element face. Since the point lies in the
// reference element, it then lies on the face.
bool OnFacePlane(const Point &p, const std::vector<Point> &face)
{
  constexpr double tol = 1.0e-10;
  const auto n = Cross(Sub(face[1], face[0]), Sub(face[2], face[0]));
  return std::abs(Dot(n, Sub(p, face[0]))) < tol * std::sqrt(Dot(n, n));
}

}  // namespace

namespace
{

// See InheritCrackSides, with the geometries of the fine elements.
CrackSides InheritCrackSides(int dim, const std::vector<mfem::Geometry::Type> &fine_geoms,
                             const mfem::CoarseFineTransformations &cf,
                             const CrackSides &coarse_sides)
{
  const int ne = static_cast<int>(fine_geoms.size());
  MFEM_VERIFY(cf.embeddings.Size() >= ne,
              "Invalid coarse-to-fine transformations for interior boundary sides!");

  // The embeddings give the geometry of the fine element, but the point matrices are those
  // of the parent geometry, which differ for the tetrahedra of a refined pyramid: a parent
  // is a pyramid if any of its fine elements is, and otherwise has their geometry.
  std::vector<mfem::Geometry::Type> coarse_geoms(coarse_sides.copy.size(),
                                                 mfem::Geometry::INVALID);
  for (int e = 0; e < ne; e++)
  {
    const int parent = cf.embeddings[e].parent;
    MFEM_VERIFY(parent >= 0 && static_cast<std::size_t>(parent) < coarse_geoms.size(),
                "Invalid parent element for interior boundary sides!");
    if (coarse_geoms[parent] != mfem::Geometry::PYRAMID)
    {
      coarse_geoms[parent] = fine_geoms[e];
    }
  }

  CrackSides sides;
  sides.copy.assign(ne, 0);
  sides.split.assign(ne, 0);
  std::vector<Point> pts, sub, face_pts;
  for (int e = 0; e < ne; e++)
  {
    const auto &emb = cf.embeddings[e];
    const std::uint32_t coarse_copy = coarse_sides.copy[emb.parent];
    const std::uint32_t coarse_split = coarse_sides.split[emb.parent];
    if (!coarse_copy && !coarse_split)
    {
      continue;
    }
    const auto coarse_geom = coarse_geoms[emb.parent];
    const auto &ct = GetReferenceTopology(coarse_geom);
    const int nv_c = static_cast<int>(ct.verts.size());
    const int ne_c = static_cast<int>(ct.edges.size());
    const auto &ft = GetReferenceTopology(fine_geoms[e]);
    const int nv_f = static_cast<int>(ft.verts.size());
    const int ne_f = static_cast<int>(ft.edges.size());

    // The vertices of the fine element are the first columns of its point matrix (which
    // has as many columns as the parent geometry has vertices).
    MFEM_VERIFY(static_cast<int>(emb.matrix) < cf.point_matrices[coarse_geom].SizeK() &&
                    cf.point_matrices[coarse_geom].SizeJ() >= nv_f,
                "Unexpected refinement point matrix for interior boundary sides!");
    const mfem::DenseMatrix &pm = cf.point_matrices[coarse_geom](emb.matrix);
    pts.resize(nv_f);
    for (int i = 0; i < nv_f; i++)
    {
      pts[i] = {pm(0, i), (pm.Height() > 1) ? pm(1, i) : 0.0,
                (pm.Height() > 2) ? pm(2, i) : 0.0};
    }

    // Local index of the smallest entity of the parent containing all given points (the
    // vertices of a fine element entity), or -1 when the entity is interior to the parent.
    auto ParentEntity = [&](const std::vector<Point> &x) -> int
    {
      if (x.size() == 1)
      {
        for (int i = 0; i < nv_c; i++)
        {
          if (SamePoint(x[0], ct.verts[i]))
          {
            return i;
          }
        }
      }
      for (int j = 0; j < ne_c; j++)
      {
        const auto &a = ct.verts[ct.edges[j][0]], &b = ct.verts[ct.edges[j][1]];
        if (std::all_of(x.begin(), x.end(),
                        [&](const Point &p) { return OnSegment(p, a, b); }))
        {
          return nv_c + j;
        }
      }
      for (std::size_t k = 0; k < ct.faces.size(); k++)
      {
        face_pts.clear();
        for (auto v : ct.faces[k])
        {
          face_pts.push_back(ct.verts[v]);
        }
        if (std::all_of(x.begin(), x.end(),
                        [&](const Point &p) { return OnFacePlane(p, face_pts); }))
        {
          return nv_c + ne_c + static_cast<int>(k);
        }
      }
      return -1;
    };

    // Set the flags of the fine element entity with local index b.
    auto Inherit = [&](const std::vector<Point> &x, int b)
    {
      const int pb = ParentEntity(x);
      if (pb < 0)
      {
        return;  // Interior entities are never constrained
      }
      if ((coarse_split >> pb) & 1u)
      {
        sides.split[e] |= (std::uint32_t(1) << b);
        sides.copy[e] |= (((coarse_copy >> pb) & 1u) << b);
      }
      else if (coarse_copy & ct.closure[pb])
      {
        // Not split, but possibly constrained to DOFs of split entities which are read as
        // copies by the parent. Reading copies has no effect for unconstrained entities.
        sides.copy[e] |= (std::uint32_t(1) << b);
      }
    };

    for (int i = 0; i < nv_f; i++)
    {
      sub.assign(1, pts[i]);
      Inherit(sub, i);
    }
    if (dim > 1)
    {
      for (int j = 0; j < ne_f; j++)
      {
        sub.assign({pts[ft.edges[j][0]], pts[ft.edges[j][1]]});
        Inherit(sub, nv_f + j);
      }
    }
    if (dim == 3)
    {
      for (std::size_t k = 0; k < ft.faces.size(); k++)
      {
        sub.clear();
        for (auto v : ft.faces[k])
        {
          sub.push_back(pts[v]);
        }
        Inherit(sub, nv_f + ne_f + static_cast<int>(k));
      }
    }
  }
  return sides;
}

std::vector<mfem::Geometry::Type> GetElementGeometries(const mfem::Mesh &mesh)
{
  std::vector<mfem::Geometry::Type> geoms(mesh.GetNE());
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    geoms[e] = mesh.GetElementGeometry(e);
  }
  return geoms;
}

}  // namespace

CrackSides InheritCrackSides(const mfem::Mesh &fine_mesh,
                             const mfem::CoarseFineTransformations &cf,
                             const CrackSides &coarse_sides)
{
  return InheritCrackSides(fine_mesh.Dimension(), GetElementGeometries(fine_mesh), cf,
                           coarse_sides);
}

bool ReconstructCrackSides(const mfem::ParMesh &mesh,
                           const mfem::Array<int> &bdr_attr_marker, CrackSides &sides)
{
  MPI_Comm comm = mesh.GetComm();
  const int dim = mesh.Dimension(), ne = mesh.GetNE();

  // MFEM does not support the derefinement of anisotropic refinements in 3D, and its
  // derefinement transformations do not support pyramids (refined into pyramids and
  // tetrahedra: the point matrices are computed for the geometry of the fine element). An
  // isotropic refinement reduces the element size by 8 in 3D.
  int supported = mesh.Nonconforming();
  for (int e = 0; supported && e < ne; e++)
  {
    const int depth = mesh.pncmesh->GetElementDepth(e);
    supported = (mesh.GetElementGeometry(e) != mfem::Geometry::PYRAMID) &&
                (dim < 3 || (depth <= 10 && mesh.pncmesh->GetElementSizeReduction(e) ==
                                                (1 << (3 * depth))));
  }
  Mpi::GlobalMin(1, &supported, comm);
  if (!supported)
  {
    return false;
  }

  // Gather a copy of the mesh on the root process, keeping track of the process and local
  // index of each element (the element ordering after rebalancing is not specified). The
  // copy needs its own nodes, which would otherwise be shared with (and modified for) the
  // original mesh.
  mfem::ParMesh gmesh(mesh, true);
  std::vector<int> orig_rank, orig_index;
  {
    mfem::L2_FECollection fec(0, dim);
    mfem::ParFiniteElementSpace fespace(&gmesh, &fec, 2);
    mfem::ParGridFunction id(&fespace);
    auto *h_id = id.HostWrite();
    for (int e = 0; e < ne; e++)
    {
      h_id[e] = Mpi::Rank(comm);
      h_id[ne + e] = e;
    }
    mfem::Array<int> partition(ne);
    partition = 0;
    gmesh.Rebalance(partition);
    fespace.Update();
    id.Update();
    const int gne = gmesh.GetNE();
    const auto *h_gid = id.HostRead();
    orig_rank.resize(gne);
    orig_index.resize(gne);
    for (int e = 0; e < gne; e++)
    {
      orig_rank[e] = static_cast<int>(std::lround(h_gid[e]));
      orig_index[e] = static_cast<int>(std::lround(h_gid[gne + e]));
    }
  }

  // Derefine until there are no hanging entities, one level of the refinement hierarchy at
  // a time, keeping the fine element geometries and derefinement transformations.
  struct Level
  {
    std::vector<mfem::Geometry::Type> fine_geoms;
    mfem::CoarseFineTransformations cf;
  };
  std::vector<Level> levels;
  while (HasHangingEntities(gmesh))
  {
    Level level;
    level.fine_geoms = GetElementGeometries(gmesh);
    mfem::Vector zero(gmesh.GetNE());
    zero = 0.0;
    if (!gmesh.DerefineByError(zero, 1.0))
    {
      return false;  // The coarsest mesh of the hierarchy has hanging entities
    }
    level.cf = gmesh.pncmesh->GetDerefinementTransforms();
    levels.push_back(std::move(level));
  }

  // Discover the sides on the coarsest mesh and inherit them through the refinements.
  CrackSides gsides = ComputeCrackSides(gmesh, bdr_attr_marker);
  for (auto it = levels.rbegin(); it != levels.rend(); ++it)
  {
    gsides = InheritCrackSides(dim, it->fine_geoms, it->cf, gsides);
  }

  // Send the sides of each element back to its process.
  const int num_procs = Mpi::Size(comm), gne = static_cast<int>(orig_rank.size());
  std::vector<int> counts(num_procs, 0), displs(num_procs, 0);
  for (int e = 0; e < gne; e++)
  {
    counts[orig_rank[e]] += 3;
  }
  for (int p = 1; p < num_procs; p++)
  {
    displs[p] = displs[p - 1] + counts[p - 1];
  }
  std::vector<std::uint32_t> send(3 * static_cast<std::size_t>(gne)), recv(3 * ne);
  {
    std::vector<int> pos(displs);
    for (int e = 0; e < gne; e++)
    {
      auto *entry = send.data() + pos[orig_rank[e]];
      entry[0] = static_cast<std::uint32_t>(orig_index[e]);
      entry[1] = gsides.copy[e];
      entry[2] = gsides.split[e];
      pos[orig_rank[e]] += 3;
    }
  }
  int recv_count = 0;
  MPI_Scatter(counts.data(), 1, MPI_INT, &recv_count, 1, MPI_INT, 0, comm);
  MFEM_VERIFY(recv_count == 3 * ne,
              "Unexpected element count when reconstructing interior boundary sides!");
  MPI_Scatterv(send.data(), counts.data(), displs.data(), MPI_UINT32_T, recv.data(), 3 * ne,
               MPI_UINT32_T, 0, comm);
  sides.copy.assign(ne, 0);
  sides.split.assign(ne, 0);
  for (int i = 0; i < ne; i++)
  {
    const auto e = recv[3 * i];
    sides.copy[e] = recv[3 * i + 1];
    sides.split[e] = recv[3 * i + 2];
  }
  return true;
}

}  // namespace mesh

namespace fem
{

void GetElementDofEntities(const mfem::FiniteElementSpace &fespace, int e,
                           mfem::Array<int> &entities)
{
  // Mirrors the DOF layout of mfem::FiniteElementSpace::GetElementDofs: vertex DOFs, edge
  // DOFs, face DOFs (3D), and interior DOFs.
  const mfem::Mesh &mesh = *fespace.GetMesh();
  const mfem::FiniteElementCollection &fec = *fespace.FEColl();
  const int dim = mesh.Dimension();
  const auto geom = mesh.GetElementGeometry(e);
  const int order = fespace.GetElementOrder(e);
  const int nv_e = mfem::Geometry::NumVerts[geom];
  const int ne_e = (dim > 1) ? mfem::Geometry::NumEdges[geom] : 0;
  const int nvd = fec.GetNumDof(mfem::Geometry::POINT, order);
  const int ned = (dim > 1) ? fec.GetNumDof(mfem::Geometry::SEGMENT, order) : 0;
  const int nbd = (dim > 0) ? fec.GetNumDof(geom, order) : 0;
  mfem::Array<int> faces, orients;
  if (dim > 2 && fec.HasFaceDofs(geom, order))
  {
    mesh.GetElementFaces(e, faces, orients);
  }
  entities.SetSize(0);
  entities.Reserve(nvd * nv_e + ned * ne_e + nbd);
  for (int i = 0; i < nv_e; i++)
  {
    for (int j = 0; j < nvd; j++)
    {
      entities.Append(i);
    }
  }
  for (int i = 0; i < ne_e; i++)
  {
    for (int j = 0; j < ned; j++)
    {
      entities.Append(nv_e + i);
    }
  }
  for (int i = 0; i < faces.Size(); i++)
  {
    const int nfd = fec.GetNumDof(mesh.GetFaceGeometry(faces[i]), order);
    for (int j = 0; j < nfd; j++)
    {
      entities.Append(nv_e + ne_e + i);
    }
  }
  for (int j = 0; j < nbd; j++)
  {
    entities.Append(-1);
  }
}

}  // namespace fem

namespace
{

// Device kernels are in free functions: CUDA does not allow extended device lambdas in
// private member functions.

// y[idx[i]] = x[i], i = 0, ..., n - 1.
void Scatter(bool use_dev, const mfem::Array<int> &idx, const double *x, double *y)
{
  const auto *I = idx.Read(use_dev);
  mfem::forall_switch(use_dev, idx.Size(), [=] MFEM_HOST_DEVICE(int i) { y[I[i]] = x[i]; });
}

// y[i] = a x[idx[i]] + b y[i], i = 0, ..., n - 1 (y is not read for b = 0).
void Gather(bool use_dev, const mfem::Array<int> &idx, const double *x, double a, double b,
            double *y)
{
  const auto *I = idx.Read(use_dev);
  if (b == 0.0)
  {
    mfem::forall_switch(use_dev, idx.Size(),
                        [=] MFEM_HOST_DEVICE(int i) { y[i] = a * x[I[i]]; });
  }
  else
  {
    mfem::forall_switch(use_dev, idx.Size(),
                        [=] MFEM_HOST_DEVICE(int i) { y[i] = a * x[I[i]] + b * y[i]; });
  }
}

// x[idx[i]] = 0, i = 0, ..., n - 1.
void ZeroIndexed(bool use_dev, const mfem::Array<int> &idx, double *x)
{
  const auto *I = idx.Read(use_dev);
  mfem::forall_switch(use_dev, idx.Size(), [=] MFEM_HOST_DEVICE(int i) { x[I[i]] = 0.0; });
}

}  // namespace

BrokenProlongation::BrokenProlongation(const Operator &P,
                                       const mfem::Array<int> &copy_ldofs_,
                                       const mfem::Array<int> &split_tdofs_)
  : Operator(P.Height() + copy_ldofs_.Size(), P.Width() + split_tdofs_.Size()), P(P),
    hP(dynamic_cast<const mfem::HypreParMatrix *>(&P)), vsize(P.Height()), tsize(P.Width())
{
  copy_ldofs = copy_ldofs_;
  split_tdofs = split_tdofs_;
  tx.SetSize(tsize);
  lx.SetSize(vsize);
  tx.UseDevice(true);
  lx.UseDevice(true);
}

void BrokenProlongation::Mult(const Vector &x, Vector &y) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Invalid vector sizes for BrokenProlongation::Mult!");
  const Vector xt(const_cast<Vector &>(x), 0, tsize);
  Vector yl(y, 0, vsize);
  P.Mult(xt, yl);

  // Copies: y_c = (P z)[copy], with z = x with the split true DOFs replaced by their
  // copies. The second application of P is collective, so it is done on every process even
  // when there are no local copies.
  const bool use_dev = x.UseDevice() || y.UseDevice();
  tx = xt;
  Scatter(use_dev, split_tdofs, x.Read(use_dev) + tsize, tx.ReadWrite(use_dev));
  P.Mult(tx, lx);
  Gather(use_dev, copy_ldofs, lx.Read(use_dev), 1.0, 0.0, y.ReadWrite(use_dev) + vsize);
}

void BrokenProlongation::AddCopyTranspose(const Vector &x, double a, double b, Vector &y,
                                          bool abs) const
{
  const bool use_dev = x.UseDevice() || y.UseDevice();
  lx = 0.0;
  Scatter(use_dev, copy_ldofs, x.Read(use_dev) + vsize, lx.ReadWrite(use_dev));
  if (abs && hP)
  {
    hP->AbsMultTranspose(1.0, lx, 0.0, tx);
  }
  else
  {
    // The prolongation of a conforming space (not a HypreParMatrix) has only nonnegative
    // entries.
    P.MultTranspose(lx, tx);
  }

  // Split true DOFs: to the copies; unsplit ones (from constrained copy rows): to y_t.
  Gather(use_dev, split_tdofs, tx.Read(use_dev), a, b, y.ReadWrite(use_dev) + tsize);
  ZeroIndexed(use_dev, split_tdofs, tx.ReadWrite(use_dev));
  Vector yt(y, 0, tsize);
  yt.Add(a, tx);
}

void BrokenProlongation::MultTranspose(const Vector &x, Vector &y) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Invalid vector sizes for BrokenProlongation::MultTranspose!");
  const Vector xl(const_cast<Vector &>(x), 0, vsize);
  Vector yt(y, 0, tsize);
  P.MultTranspose(xl, yt);
  AddCopyTranspose(x, 1.0, 0.0, y, false);
}

void BrokenProlongation::AbsMultTranspose(double a, const Vector &x, double b,
                                          Vector &y) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Invalid vector sizes for BrokenProlongation::AbsMultTranspose!");
  const Vector xl(const_cast<Vector &>(x), 0, vsize);
  Vector yt(y, 0, tsize);
  if (hP)
  {
    hP->AbsMultTranspose(a, xl, b, yt);
  }
  else
  {
    // The prolongation of a conforming space (not a HypreParMatrix) has only nonnegative
    // entries.
    Vector t(tsize);
    t.UseDevice(true);
    P.MultTranspose(xl, t);
    if (b == 0.0)
    {
      yt = 0.0;
    }
    else
    {
      yt *= b;
    }
    yt.Add(a, t);
  }
  AddCopyTranspose(x, a, b, y, true);
}

}  // namespace palace
