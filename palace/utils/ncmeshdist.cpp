// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ncmeshdist.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "utils/communication.hpp"
#include "utils/geodata.hpp"

namespace palace::mesh
{

namespace
{

// Append-only byte buffer and its reader, for the serialized mesh data.
class ByteWriter
{
public:
  std::vector<char> buf;

  template <typename T>
  void Put(const T &v)
  {
    const char *p = reinterpret_cast<const char *>(&v);
    buf.insert(buf.end(), p, p + sizeof(T));
  }

  template <typename T>
  void Put(const T *v, std::size_t n)
  {
    const char *p = reinterpret_cast<const char *>(v);
    buf.insert(buf.end(), p, p + n * sizeof(T));
  }
};

class ByteReader
{
  const char *p, *end;

public:
  explicit ByteReader(const std::vector<char> &buf)
    : p(buf.data()), end(buf.data() + buf.size())
  {
  }

  template <typename T>
  void Get(T *v, std::size_t n = 1)
  {
    MFEM_VERIFY(p + n * sizeof(T) <= end, "Corrupt serialized mesh data!");
    std::memcpy(v, p, n * sizeof(T));
    p += n * sizeof(T);
  }

  template <typename T>
  T Get()
  {
    T v;
    Get(&v);
    return v;
  }
};

// A serial mesh with only elements, boundary elements and vertices (no edges, faces or
// other topology tables), and the NCMesh of its elements as the coarse (root) elements.
// It provides exactly what mfem::ParMesh needs to construct a nonconforming parallel mesh
// from a serial one: the NCMesh, the element count, meshgen and the attribute lists.
class RootElementMesh : public mfem::Mesh
{
public:
  RootElementMesh(int dim, int sdim, int nv, int ne, int nbe)
    : mfem::Mesh(dim, nv, ne, nbe, sdim)
  {
  }

  void AddRootElement(int geom, int attr, const int *v, bool bdr)
  {
    mfem::Element *el = NewElement(geom);
    el->SetVertices(v);
    el->SetAttribute(attr);
    bdr ? AddBdrElement(el) : AddElement(el);
  }

  void BuildNCMesh()
  {
    SetMeshGen();
    SetAttributes();
    ncmesh = new mfem::NCMesh(this);
  }
};

// A mesh which only holds an NCMesh (no elements, vertices or topology tables), with the
// element count, mesh generator and attribute lists the mfem::ParMesh constructor takes
// from a nonconforming serial mesh.
class NCMeshShell : public mfem::Mesh
{
public:
  NCMeshShell(mfem::NCMesh *nc, int ne, int mgen, const std::vector<int> &attr,
              const std::vector<int> &bdr_attr)
  {
    Dim = nc->Dimension();
    spaceDim = nc->SpaceDimension();
    NumOfElements = ne;
    meshgen = mgen;
    attributes.SetSize(static_cast<int>(attr.size()));
    std::copy(attr.begin(), attr.end(), attributes.begin());
    bdr_attributes.SetSize(static_cast<int>(bdr_attr.size()));
    std::copy(bdr_attr.begin(), bdr_attr.end(), bdr_attributes.begin());
    ncmesh = nc;
  }

  ~NCMeshShell() override
  {
    // There are no element objects to free.
    NumOfElements = 0;
  }
};

// Serialized root mesh: elements and boundary elements (geometry, attribute, vertices),
// vertex coordinates and, for a curved mesh, the description of its nodal space.
struct RootMeshData
{
  int dim = 0, sdim = 0, nv = 0;
  std::vector<int> geom, attr, verts, bgeom, battr, bverts;
  std::vector<double> coords;
  int curved = 0, vdim = 0, ordering = 0;
  std::string fec;

  void SetCurved(const mfem::Mesh &m)
  {
    const mfem::GridFunction *nodes = m.GetNodes();
    curved = (nodes != nullptr);
    if (curved)
    {
      fec = nodes->FESpace()->FEColl()->Name();
      vdim = nodes->FESpace()->GetVDim();
      ordering = nodes->FESpace()->GetOrdering();
    }
  }

  void Pack(ByteWriter &w) const
  {
    w.Put(dim);
    w.Put(sdim);
    w.Put(nv);
    for (const auto *v : {&geom, &attr, &verts, &bgeom, &battr, &bverts})
    {
      w.Put<std::int64_t>(v->size());
      w.Put(v->data(), v->size());
    }
    w.Put(coords.data(), coords.size());
    w.Put(curved);
    w.Put(vdim);
    w.Put(ordering);
    w.Put<int>(fec.size());
    w.Put(fec.data(), fec.size());
  }

  void Unpack(ByteReader &r)
  {
    r.Get(&dim);
    r.Get(&sdim);
    r.Get(&nv);
    for (auto *v : {&geom, &attr, &verts, &bgeom, &battr, &bverts})
    {
      v->resize(r.Get<std::int64_t>());
      r.Get(v->data(), v->size());
    }
    coords.resize(static_cast<std::size_t>(nv) * sdim);
    r.Get(coords.data(), coords.size());
    r.Get(&curved);
    r.Get(&vdim);
    r.Get(&ordering);
    fec.resize(r.Get<int>());
    r.Get(fec.data(), fec.size());
  }

  // Build the root-element mesh and its NCMesh, releasing the serialized arrays.
  std::unique_ptr<RootElementMesh> Build()
  {
    auto rmesh = std::make_unique<RootElementMesh>(dim, sdim, nv, int(geom.size()),
                                                   int(bgeom.size()));
    for (int i = 0, off = 0; i < int(geom.size()); i++)
    {
      rmesh->AddRootElement(geom[i], attr[i], verts.data() + off, false);
      off += mfem::Geometry::NumVerts[geom[i]];
    }
    for (int i = 0, off = 0; i < int(bgeom.size()); i++)
    {
      rmesh->AddRootElement(bgeom[i], battr[i], bverts.data() + off, true);
      off += mfem::Geometry::NumVerts[bgeom[i]];
    }
    for (int i = 0; i < nv; i++)
    {
      rmesh->AddVertex(coords.data() + static_cast<std::size_t>(i) * sdim);
    }
    for (auto *v : {&geom, &attr, &verts, &bgeom, &battr, &bverts})
    {
      std::vector<int>().swap(*v);
    }
    std::vector<double>().swap(coords);
    rmesh->BuildNCMesh();
    return rmesh;
  }
};

// Send/receive a buffer whose length may exceed INT_MAX, in chunks.
constexpr std::int64_t MAX_CHUNK = std::numeric_limits<int>::max() / 8;

template <typename T>
void SendLarge(const T *buf, std::int64_t n, int dest, MPI_Comm comm)
{
  for (std::int64_t off = 0; off < n; off += MAX_CHUNK)
  {
    MPI_Send(buf + off, static_cast<int>(std::min(MAX_CHUNK, n - off)), mpi::DataType<T>(),
             dest, 0, comm);
  }
}

template <typename T>
void RecvLarge(T *buf, std::int64_t n, int src, MPI_Comm comm)
{
  for (std::int64_t off = 0; off < n; off += MAX_CHUNK)
  {
    MPI_Recv(buf + off, static_cast<int>(std::min(MAX_CHUNK, n - off)), mpi::DataType<T>(),
             src, 0, comm, MPI_STATUS_IGNORE);
  }
}

// Broadcast a byte buffer from the node root.
void BroadcastBuffer(std::vector<char> &buf, MPI_Comm node_comm)
{
  std::int64_t len = static_cast<std::int64_t>(buf.size());
  Mpi::Broadcast(1, &len, 0, node_comm);
  buf.resize(len);
  Mpi::BroadcastLarge(len, buf.data(), 0, node_comm);
}

// Every process of a node sends its node root a list of keys and receives the values the
// node root produces for them with fill(key, values) (which appends to values), in order.
template <typename T, typename F>
std::vector<T> PullFromNodeRoot(MPI_Comm node_comm, const std::vector<std::int64_t> &keys,
                                F &&fill)
{
  std::vector<T> values;
  if (Mpi::Root(node_comm))
  {
    std::vector<std::int64_t> rkeys;
    for (int k = 1; k < Mpi::Size(node_comm); k++)
    {
      std::int64_t n;
      MPI_Recv(&n, 1, MPI_INT64_T, k, 0, node_comm, MPI_STATUS_IGNORE);
      rkeys.resize(n);
      RecvLarge(rkeys.data(), n, k, node_comm);
      values.clear();
      for (auto key : rkeys)
      {
        fill(key, values);
      }
      n = static_cast<std::int64_t>(values.size());
      MPI_Send(&n, 1, MPI_INT64_T, k, 0, node_comm);
      SendLarge(values.data(), n, k, node_comm);
    }
    values.clear();
    values.shrink_to_fit();
    for (auto key : keys)
    {
      fill(key, values);
    }
  }
  else
  {
    std::int64_t n = static_cast<std::int64_t>(keys.size());
    MPI_Send(&n, 1, MPI_INT64_T, 0, 0, node_comm);
    SendLarge(keys.data(), n, 0, node_comm);
    MPI_Recv(&n, 1, MPI_INT64_T, 0, 0, node_comm, MPI_STATUS_IGNORE);
    values.resize(n);
    RecvLarge(values.data(), n, 0, node_comm);
  }
  return values;
}

// Append the element-wise node values of element i of the serial mesh, as the
// mfem::ParMesh constructor copies them.
void AppendElementNodes(const mfem::Mesh &smesh, int i, std::vector<double> &values)
{
  const mfem::GridFunction &nodes = *smesh.GetNodes();
  mfem::Array<int> vdofs;
  nodes.FESpace()->GetElementVDofs(i, vdofs);
  mfem::Vector el_nodes;
  nodes.GetSubVector(vdofs, el_nodes);
  values.insert(values.end(), el_nodes.begin(), el_nodes.end());
}

// Set the curved mesh nodes of the parallel mesh from element-wise values in its local
// element order, and the vertex coordinates from the nodes (as the mfem::ParMesh
// constructor does with mfem::Mesh::SetVerticesFromNodes).
void SetElementNodes(mfem::ParMesh &pmesh, const RootMeshData &data,
                     const std::vector<double> &values)
{
  auto *nfec = mfem::FiniteElementCollection::New(data.fec.c_str());
  auto *pfes = new mfem::ParFiniteElementSpace(&pmesh, nfec, data.vdim, data.ordering);
  auto *pnodes = new mfem::ParGridFunction(pfes);
  pnodes->MakeOwner(nfec);
  mfem::Array<int> vdofs;
  std::size_t off = 0;
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    pfes->GetElementVDofs(e, vdofs);
    MFEM_VERIFY(off + vdofs.Size() <= values.size(),
                "Mismatch in the received curved mesh node values!");
    pnodes->SetSubVector(vdofs, const_cast<double *>(values.data()) + off);
    off += vdofs.Size();
  }
  MFEM_VERIFY(off == values.size(), "Mismatch in the received curved mesh node values!");
  pmesh.NewNodes(*pnodes, true);
  const int nv = pmesh.GetNV(), sdim = pmesh.SpaceDimension();
  mfem::Vector vert(sdim * nv), vals;
  for (int d = 0; d < sdim; d++)
  {
    pnodes->GetNodalValues(vals, d + 1);
    for (int i = 0; i < nv; i++)
    {
      vert(d * nv + i) = vals(i);
    }
  }
  pmesh.SetVertices(vert);
}

// An element of the given geometry with the vertex numbering of the reference element (not
// from the element memory pool of a mfem::Mesh, which mfem::Mesh::NewElement may use).
std::unique_ptr<mfem::Element> ReferenceElement(mfem::Geometry::Type geom)
{
  switch (geom)
  {
    case mfem::Geometry::TRIANGLE:
      return std::make_unique<mfem::Triangle>();
    case mfem::Geometry::SQUARE:
      return std::make_unique<mfem::Quadrilateral>();
    case mfem::Geometry::TETRAHEDRON:
      return std::make_unique<mfem::Tetrahedron>();
    case mfem::Geometry::CUBE:
      return std::make_unique<mfem::Hexahedron>();
    case mfem::Geometry::PRISM:
      return std::make_unique<mfem::Wedge>();
    case mfem::Geometry::PYRAMID:
      return std::make_unique<mfem::Pyramid>();
    default:
      MFEM_ABORT("Unsupported element geometry for nonconforming mesh distribution!");
      return nullptr;
  }
}

// Number of children of a refined element of the refinement tree.
template <typename Element>
int NumChildren(const Element &el)
{
  constexpr int max_children = sizeof(el.child) / sizeof(el.child[0]);
  int n = 0;
  while (n < max_children && el.child[n] >= 0)
  {
    n++;
  }
  return n;
}

// Permute the lines of a section of an MFEM nonconforming mesh text (the n lines after the
// line with its name and the line with n): line k becomes the old line order[k]. A missing
// optional section is left out.
void PermuteSectionLines(std::string &text, const std::string &name,
                         const std::vector<int> &order, bool required)
{
  const std::string key = "\n" + name + "\n";
  const std::size_t pos = text.find(key);
  MFEM_VERIFY(pos != std::string::npos || !required,
              "Missing section " << name << " in the nonconforming mesh!");
  if (pos == std::string::npos)
  {
    return;
  }
  std::size_t p = pos + key.size();
  std::size_t eol = text.find('\n', p);
  const long long n = std::stoll(text.substr(p, eol - p));
  MFEM_VERIFY(n == static_cast<long long>(order.size()),
              "Unexpected size of section " << name << " in the nonconforming mesh!");
  std::vector<std::size_t> start(order.size() + 1);
  p = eol + 1;
  for (std::size_t i = 0; i < order.size(); i++)
  {
    start[i] = p;
    eol = text.find('\n', p);
    MFEM_VERIFY(eol != std::string::npos,
                "Truncated section " << name << " in the nonconforming mesh!");
    p = eol + 1;
  }
  start[order.size()] = p;
  std::string lines;
  lines.reserve(p - start[0]);
  for (const int i : order)
  {
    lines.append(text, start[i], start[i + 1] - start[i]);
  }
  text.replace(start[0], p - start[0], lines);
}

}  // namespace

std::unique_ptr<mfem::ParMesh> DistributeNCMesh(MPI_Comm comm,
                                                std::unique_ptr<mfem::Mesh> &smesh,
                                                const int *partitioning, bool local_reorder)
{
  MPI_Comm node_comm;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                      &node_comm);
  const bool node_root = Mpi::Root(node_comm);
  MFEM_VERIFY(node_root == bool(smesh) && (!smesh || (partitioning && smesh->ncmesh)),
              "The serial nonconforming mesh and partitioning must be held by the first "
              "process of each node!");

  // The node root sends the NCMesh of the serial mesh in the MFEM nonconforming mesh
  // format, without the curved mesh nodes, which is what every process would otherwise
  // parse from the serialized serial mesh (so that the NCMesh, and the meshes refined from
  // it, are the same), plus the element count, mesh generator and attribute lists the
  // mfem::ParMesh constructor takes from the serial mesh, and the partitioning.
  //
  // With local_reorder, the root elements are numbered along a Hilbert curve through their
  // centroids (the same order on every node), so that the elements of each process, which
  // are the leaves of its root elements in the order of the roots (and the vertices, edges
  // and faces numbered from them), are local in space for any partition: the partition is
  // unchanged, and every process keeps the serial index of each root element (origin).
  RootMeshData data;
  std::string nc_text;
  std::vector<int> part, attr, bdr_attr, origin;
  int meshgen = 0;
  {
    std::vector<char> buf;
    if (node_root)
    {
      const int ne = smesh->GetNE(), sdim = smesh->SpaceDimension();
      if (local_reorder)
      {
        std::vector<double> xc(static_cast<std::size_t>(ne) * sdim, 0.0);
        for (int e = 0; e < ne; e++)
        {
          const mfem::Element *el = smesh->GetElement(e);
          const int nv = el->GetNVertices();
          const int *v = el->GetVertices();
          for (int j = 0; j < nv; j++)
          {
            const double *x = smesh->GetVertex(v[j]);
            for (int d = 0; d < sdim; d++)
            {
              xc[static_cast<std::size_t>(e) * sdim + d] += x[d] / nv;
            }
          }
        }
        origin = HilbertOrder(xc.data(), ne, sdim);
      }
      std::ostringstream os;
      os << std::scientific;
      os.precision(std::numeric_limits<double>::max_digits10);
      smesh->ncmesh->Print(os);
      if (smesh->GetNodes())
      {
        // The NCMesh of a curved mesh has no vertex coordinates. Its space dimension is
        // set when the serial mesh reads its nodes, so give the vertex coordinates here
        // instead (they are removed by mfem::Mesh::NewNodes).
        os << "\ncoordinates\n"
           << smesh->GetNV() << "\n"
           << smesh->SpaceDimension() << "\n";
        for (int i = 0; i < smesh->GetNV(); i++)
        {
          const double *x = smesh->GetVertex(i);
          for (int d = 0; d < smesh->SpaceDimension(); d++)
          {
            os << (d ? " " : "") << x[d];
          }
          os << "\n";
        }
      }
      nc_text = os.str();
      if (!origin.empty())
      {
        PermuteSectionLines(nc_text, "elements", origin, true);
        PermuteSectionLines(nc_text, "root_state", origin, false);
      }
      data.SetCurved(*smesh);
      ByteWriter w;
      data.Pack(w);
      w.Put<int>(smesh->MeshGenerator());
      w.Put<std::int64_t>(smesh->attributes.Size());
      w.Put(smesh->attributes.GetData(), smesh->attributes.Size());
      w.Put<std::int64_t>(smesh->bdr_attributes.Size());
      w.Put(smesh->bdr_attributes.GetData(), smesh->bdr_attributes.Size());
      w.Put<std::int64_t>(ne);
      if (origin.empty())
      {
        w.Put(partitioning, ne);
      }
      else
      {
        std::vector<int> p(ne);
        for (int k = 0; k < ne; k++)
        {
          p[k] = partitioning[origin[k]];
        }
        w.Put(p.data(), ne);
      }
      w.Put<std::int64_t>(origin.size());
      w.Put(origin.data(), origin.size());
      w.Put<std::int64_t>(nc_text.size());
      w.Put(nc_text.data(), nc_text.size());
      std::string().swap(nc_text);
      buf = std::move(w.buf);
    }
    BroadcastBuffer(buf, node_comm);
    ByteReader r(buf);
    data.Unpack(r);
    meshgen = r.Get<int>();
    attr.resize(r.Get<std::int64_t>());
    r.Get(attr.data(), attr.size());
    bdr_attr.resize(r.Get<std::int64_t>());
    r.Get(bdr_attr.data(), bdr_attr.size());
    part.resize(r.Get<std::int64_t>());
    r.Get(part.data(), part.size());
    origin.resize(r.Get<std::int64_t>());
    r.Get(origin.data(), origin.size());
    nc_text.resize(r.Get<std::int64_t>());
    r.Get(nc_text.data(), nc_text.size());
  }

  // Every process parses the NCMesh, then builds the ParMesh.
  std::unique_ptr<mfem::ParMesh> pmesh;
  {
    std::istringstream is(std::move(nc_text));
    std::string header;
    std::getline(is, header);
    MFEM_VERIFY(header == "MFEM NC mesh v1.0", "Unexpected nonconforming mesh format!");
    int curved = 0, is_nc = 0;
    auto *nc = new mfem::NCMesh(is, 10, curved, is_nc);
    NCMeshShell shell(nc, static_cast<int>(part.size()), meshgen, attr, bdr_attr);
    pmesh = std::make_unique<mfem::ParMesh>(comm, shell, part.data());
  }

  // Curved mesh nodes: each process receives from its node root the node values of its
  // elements, in increasing serial element order (the local element order of the
  // ParMesh), as the mfem::ParMesh constructor would copy them.
  if (data.curved)
  {
    const int rank = Mpi::Rank(comm);
    std::vector<std::int64_t> keys;
    keys.reserve(pmesh->GetNE());
    for (std::size_t i = 0; i < part.size(); i++)
    {
      if (part[i] == rank)
      {
        keys.push_back(origin.empty() ? i : origin[i]);
      }
    }
    std::vector<int>().swap(part);
    auto values = PullFromNodeRoot<double>(
        node_comm, keys, [&](std::int64_t i, std::vector<double> &v)
        { AppendElementNodes(*smesh, static_cast<int>(i), v); });
    SetElementNodes(*pmesh, data, values);
  }
  MPI_Comm_free(&node_comm);
  SetRootElementOrigin(*pmesh, std::move(origin));
  return pmesh;
}

std::unique_ptr<mfem::ParMesh> DistributeNCMeshTree(MPI_Comm comm,
                                                    std::unique_ptr<mfem::Mesh> &smesh,
                                                    const int *partitioning)
{
  MPI_Comm node_comm;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                      &node_comm);
  const bool node_root = Mpi::Root(node_comm);
  MFEM_VERIFY(node_root == bool(smesh) && (!smesh || (partitioning && smesh->ncmesh)),
              "The serial nonconforming mesh and partitioning must be held by the first "
              "process of each node!");
  const int rank = Mpi::Rank(comm);
  const double t0 = MPI_Wtime();

  // On the node root: the root elements (the coarse mesh of the refinement tree), their
  // boundary faces, the refinement tree in depth-first order (per element, the refinement
  // type and, if refined, the number of children), and the serial element index of each
  // leaf in that order. The replay reproduces the mesh if the root elements form a
  // conforming mesh (their corners are top-level nodes), the refinement is isotropic, and
  // every boundary element lies on a face of a root element, whose attribute the
  // refinement passes on. Otherwise (for example, for a mesh rebuilt with hanging vertices
  // in its root elements), return null for the caller to distribute the serial mesh.
  RootMeshData data;
  std::vector<std::uint8_t> tree;
  std::vector<int> root_part, serial_leaf;
  std::int64_t num_leaves = 0;
  {
    std::vector<char> buf;
    const char *unsupported = nullptr;
    if (node_root)
    {
      mfem::NCMesh &nc = *smesh->ncmesh;
      const int nroots = nc.GetNumRootElements();
      data.dim = smesh->Dimension();
      data.sdim = smesh->SpaceDimension();
      data.SetCurved(*smesh);
      const int iso_ref_type = (data.dim == 3) ? 7 : 3;

      // Tree, leaves and the node of each serial mesh vertex.
      std::vector<int> node_of_vertex(smesh->GetNV(), -1), leaves_per_root(nroots, 0);
      serial_leaf.reserve(smesh->GetNE());
      std::vector<int> stack;
      mfem::Array<int> v;
      for (int r = 0; r < nroots; r++)
      {
        stack.push_back(r);
        while (!stack.empty())
        {
          const int e = stack.back();
          stack.pop_back();
          const auto &el = nc.GetElement(e);
          if (el.ref_type)
          {
            if (el.ref_type != iso_ref_type)
            {
              unsupported = "anisotropic refinement";
            }
            const int nch = NumChildren(el);
            tree.push_back(static_cast<std::uint8_t>(el.ref_type));
            tree.push_back(static_cast<std::uint8_t>(nch));
            for (int k = nch - 1; k >= 0; k--)
            {
              stack.push_back(el.child[k]);
            }
          }
          else
          {
            tree.push_back(0);
            MFEM_VERIFY(el.index >= 0 && el.index < smesh->GetNE(),
                        "Unexpected leaf element of the nonconforming mesh!");
            serial_leaf.push_back(el.index);
            leaves_per_root[r]++;
            smesh->GetElementVertices(el.index, v);
            for (int k = 0; k < v.Size(); k++)
            {
              node_of_vertex[v[k]] = el.node[k];
            }
          }
        }
      }
      num_leaves = static_cast<std::int64_t>(serial_leaf.size());
      MFEM_VERIFY(num_leaves == smesh->GetNE(),
                  "Leaf count of the refinement tree does not match the mesh!");

      // Root elements, with the top-level nodes as vertices.
      mfem::Array<int> corners;
      int max_node = -1;
      for (int r = 0; r < nroots; r++)
      {
        const auto &el = nc.GetElement(r);
        nc.GetElementCornerNodes(r, corners);
        data.geom.push_back(el.Geom());
        data.attr.push_back(el.GetAttribute());
        for (int k = 0; k < corners.Size(); k++)
        {
          const auto &nd = nc.GetNode(corners[k]);
          if (nd.p1 != nd.p2)
          {
            unsupported = "root elements with hanging vertices";
          }
          data.verts.push_back(corners[k]);
          max_node = std::max(max_node, corners[k]);
        }
      }
      data.nv = max_node + 1;
      data.coords.assign(static_cast<std::size_t>(data.nv) * data.sdim, 0.0);
      for (int t : data.verts)
      {
        const int vi = nc.GetNode(t).vert_index;
        MFEM_VERIFY(vi >= 0, "Root element corner is not a mesh vertex!");
        std::copy(smesh->GetVertex(vi), smesh->GetVertex(vi) + data.sdim,
                  data.coords.begin() + static_cast<std::size_t>(t) * data.sdim);
      }

      // Boundary elements of the root elements: the root face of each boundary element.
      // A boundary segment in 2D is an edge of an edge of a root element, found from the
      // parents of its midpoint nodes. In 3D, the parent faces are followed up to a face
      // made of top-level nodes.
      struct RootFace
      {
        int attr;
        std::array<int, 4> nodes;
        bool found;
      };
      std::map<std::array<int, 4>, RootFace> root_faces;
      constexpr int max_levels = 64;
      for (int b = 0; b < smesh->GetNBE() && !unsupported; b++)
      {
        smesh->GetBdrElementVertices(b, v);
        std::array<int, 4> fn = {-1, -1, -1, -1};
        if (v.Size() == 2)
        {
          int n0 = node_of_vertex[v[0]], n1 = node_of_vertex[v[1]];
          for (int level = 0;; level++)
          {
            const auto &a = nc.GetNode(n0), &c = nc.GetNode(n1);
            if (c.p1 != c.p2 && (c.p1 == n0 || c.p2 == n0))
            {
              n1 = (c.p1 == n0) ? c.p2 : c.p1;
            }
            else if (a.p1 != a.p2 && (a.p1 == n1 || a.p2 == n1))
            {
              n0 = (a.p1 == n1) ? a.p2 : a.p1;
            }
            else
            {
              break;
            }
            if (level == max_levels)
            {
              unsupported = "boundary elements outside the edges of the root elements";
              break;
            }
          }
          fn = {n0, n0, n1, n1};
        }
        else
        {
          MFEM_VERIFY(v.Size() == 3 || v.Size() == 4,
                      "Unsupported boundary element for nonconforming mesh distribution!");
          for (int k = 0; k < v.Size(); k++)
          {
            fn[k] = node_of_vertex[v[k]];
          }
          for (int level = 0;; level++)
          {
            const auto prev = fn;
            if (nc.ParentFaceNodes(fn) == -1)
            {
              break;
            }
            if (fn == prev || level == max_levels)
            {
              unsupported = "boundary elements outside the faces of the root elements";
              break;
            }
          }
        }
        auto key = fn;
        std::sort(key.begin(), key.end());
        const int battr = smesh->GetBdrAttribute(b);
        auto [it, inserted] = root_faces.try_emplace(key, RootFace{battr, fn, false});
        if (!inserted && it->second.attr != battr)
        {
          unsupported = "boundary elements with different attributes on a root face";
        }
      }

      // Every root boundary face must be a face (an edge in 2D) of a root element.
      std::unique_ptr<mfem::Element> ref;
      for (int r = 0; r < nroots && !unsupported; r++)
      {
        const auto geom = nc.GetElement(r).Geom();
        if (!ref || ref->GetGeometryType() != geom)
        {
          ref = ReferenceElement(geom);
        }
        nc.GetElementCornerNodes(r, corners);
        const int nf = (data.dim == 3) ? ref->GetNFaces() : ref->GetNEdges();
        for (int f = 0; f < nf; f++)
        {
          std::array<int, 4> key = {-1, -1, -1, -1};
          if (data.dim == 3)
          {
            const int *fv = ref->GetFaceVertices(f);
            for (int k = 0; k < ref->GetNFaceVertices(f); k++)
            {
              key[k] = corners[fv[k]];
            }
          }
          else
          {
            const int *ev = ref->GetEdgeVertices(f);
            key = {corners[ev[0]], corners[ev[0]], corners[ev[1]], corners[ev[1]]};
          }
          std::sort(key.begin(), key.end());
          if (auto it = root_faces.find(key); it != root_faces.end())
          {
            it->second.found = true;
          }
        }
      }
      for (const auto &[key, rf] : root_faces)
      {
        if (!rf.found && !unsupported)
        {
          unsupported = "boundary elements outside the faces of the root elements";
        }
        const auto &fn = rf.nodes;
        if (fn[0] == fn[1] && fn[2] == fn[3])
        {
          data.bgeom.push_back(mfem::Geometry::SEGMENT);
          data.bverts.insert(data.bverts.end(), {fn[0], fn[2]});
        }
        else if (fn[3] < 0)
        {
          data.bgeom.push_back(mfem::Geometry::TRIANGLE);
          data.bverts.insert(data.bverts.end(), {fn[0], fn[1], fn[2]});
        }
        else
        {
          data.bgeom.push_back(mfem::Geometry::SQUARE);
          data.bverts.insert(data.bverts.end(), fn.begin(), fn.end());
        }
        data.battr.push_back(rf.attr);
      }

      // Partition for the replay of the refinement: each root element goes to the process
      // which owns most of its leaves in the given partitioning, so that the replay works
      // on compact subdomains and the final rebalancing moves few elements.
      root_part.resize(nroots);
      {
        std::size_t off = 0;
        std::map<int, int> count;
        for (int r = 0; r < nroots; r++)
        {
          count.clear();
          int best = -1, best_count = 0;
          for (int k = 0; k < leaves_per_root[r]; k++)
          {
            const int p = partitioning[serial_leaf[off + k]];
            const int c = ++count[p];
            if (c > best_count || (c == best_count && p < best))
            {
              best = p;
              best_count = c;
            }
          }
          root_part[r] = best;
          off += leaves_per_root[r];
        }
      }

      ByteWriter w;
      data.Pack(w);
      w.Put(num_leaves);
      w.Put<std::int64_t>(root_part.size());
      w.Put(root_part.data(), root_part.size());
      w.Put<std::int64_t>(tree.size());
      w.Put(tree.data(), tree.size());
      buf = std::move(w.buf);
    }
    int replayable = (unsupported == nullptr);
    MPI_Allreduce(MPI_IN_PLACE, &replayable, 1, MPI_INT, MPI_MIN, comm);
    if (!replayable)
    {
      Mpi::Print(comm, " Nonconforming mesh with {}: distributing the serial mesh\n",
                 unsupported ? unsupported : "an unsupported refinement tree");
      MPI_Comm_free(&node_comm);
      return nullptr;
    }
    BroadcastBuffer(buf, node_comm);
    ByteReader r(buf);
    data.Unpack(r);
    r.Get(&num_leaves);
    root_part.resize(r.Get<std::int64_t>());
    r.Get(root_part.data(), root_part.size());
    tree.resize(r.Get<std::int64_t>());
    r.Get(tree.data(), tree.size());
  }
  const int nroots = static_cast<int>(root_part.size());

  // Every process builds the root-element NCMesh, then the ParMesh.
  std::unique_ptr<mfem::ParMesh> pmesh;
  {
    auto rmesh = data.Build();
    pmesh = std::make_unique<mfem::ParMesh>(comm, *rmesh, root_part.data());
  }
  Mpi::Print(comm, " Distributed {:d} root elements ({:.1f} s)\n", nroots,
             MPI_Wtime() - t0);

  // For the root elements of this process, the position of each tree element in the
  // depth-first order of its subtree, the end of its subtree, and the leaf count before it.
  std::map<int, std::size_t> root_index;
  std::vector<std::int64_t> root_tree_off, root_leaf_off;
  std::vector<std::vector<std::int32_t>> subtree_end, leaf_before;
  {
    std::size_t pos = 0;
    std::int64_t leaves = 0;
    std::vector<std::int32_t> stack;
    for (int r = 0; r < nroots; r++)
    {
      // Length of the subtree of root r.
      const std::size_t start = pos;
      const std::int64_t leaf_start = leaves;
      std::int64_t pending = 1;
      while (pending > 0)
      {
        MFEM_VERIFY(pos < tree.size(), "Corrupt refinement tree!");
        pending--;
        if (tree[pos])
        {
          pending += tree[pos + 1];
          pos += 2;
        }
        else
        {
          pos++;
          leaves++;
        }
      }
      if (root_part[r] != rank)
      {
        continue;
      }
      root_index[r] = root_tree_off.size();
      root_tree_off.push_back(start);
      root_leaf_off.push_back(leaf_start);
      const std::size_t len = pos - start;
      auto &end = subtree_end.emplace_back(len, 0);
      auto &lb = leaf_before.emplace_back(len, 0);
      std::int32_t nleaf = 0;
      // Subtree ends: an element's subtree ends where the subtree of its last child ends.
      // Walk backwards keeping a stack of completed subtree ends.
      stack.clear();
      std::vector<std::int32_t> starts;
      for (std::size_t p = 0; p < len;)
      {
        starts.push_back(static_cast<std::int32_t>(p));
        p += tree[start + p] ? 2 : 1;
      }
      for (auto it = starts.rbegin(); it != starts.rend(); ++it)
      {
        const std::int32_t p = *it;
        if (tree[start + p])
        {
          std::int32_t e = p + 2;
          for (int k = 0; k < tree[start + p + 1]; k++)
          {
            e = stack.back();
            stack.pop_back();
          }
          end[p] = e;
        }
        else
        {
          end[p] = p + 1;
        }
        stack.push_back(end[p]);
      }
      for (std::int32_t p : starts)
      {
        lb[p] = nleaf;
        nleaf += (tree[start + p] == 0);
      }
    }
  }
  // Position in the tree of each local leaf: walk up to the root element recording the
  // child indices, then down the subtree.
  auto Locate = [&](std::vector<std::size_t> &tree_pos, std::vector<std::int64_t> *leaf_id)
  {
    const mfem::NCMesh &nc = *pmesh->ncmesh;
    const int ne = pmesh->GetNE();
    std::vector<int> leaf_elem(ne, -1);
    for (int e = 0; e < nc.GetNumElements(); e++)
    {
      const auto &el = nc.GetElement(e);
      if (el.IsLeaf() && el.rank == rank && el.index >= 0 && el.index < ne)
      {
        leaf_elem[el.index] = e;
      }
    }
    tree_pos.resize(ne);
    if (leaf_id)
    {
      leaf_id->resize(ne);
    }
    std::vector<int> path;
    for (int j = 0; j < ne; j++)
    {
      int e = leaf_elem[j];
      MFEM_VERIFY(e >= 0, "Local element not found in the refinement tree!");
      path.clear();
      while (nc.GetElement(e).parent >= 0)
      {
        const int parent = nc.GetElement(e).parent;
        const auto &pel = nc.GetElement(parent);
        int k = 0;
        while (pel.child[k] != e)
        {
          k++;
        }
        path.push_back(k);
        e = parent;
      }
      auto it = root_index.find(e);
      MFEM_VERIFY(it != root_index.end(), "Local element of an unexpected root element!");
      const std::size_t ri = it->second, start = root_tree_off[ri];
      const auto &end = subtree_end[ri];
      std::int32_t p = 0;
      for (auto k = path.rbegin(); k != path.rend(); ++k)
      {
        MFEM_VERIFY(tree[start + p] && *k < tree[start + p + 1],
                    "Mesh element is not in the refinement tree!");
        p += 2;
        for (int c = 0; c < *k; c++)
        {
          p = end[p];
        }
      }
      tree_pos[j] = start + p;
      if (leaf_id)
      {
        MFEM_VERIFY(tree[start + p] == 0,
                    "Mesh element is not a leaf of the refinement tree!");
        (*leaf_id)[j] = root_leaf_off[ri] + leaf_before[ri][p];
      }
    }
  };

  // Replay the refinement, level by level.
  std::vector<std::size_t> tree_pos;
  for (int level = 0;; level++)
  {
    Locate(tree_pos, nullptr);
    mfem::Array<mfem::Refinement> refs;
    for (int j = 0; j < pmesh->GetNE(); j++)
    {
      if (tree[tree_pos[j]])
      {
        refs.Append(mfem::Refinement(j, tree[tree_pos[j]]));
      }
    }
    long long nref = refs.Size();
    Mpi::GlobalSum(1, &nref, comm);
    if (nref == 0)
    {
      break;
    }
    pmesh->GeneralRefinement(refs, 1, 0);
    Mpi::Print(comm,
               " Replayed refinement level {:d}: {:d} elements refined, {:d} elements "
               "({:.1f} s)\n",
               level + 1, nref, pmesh->GetGlobalNE(), MPI_Wtime() - t0);
  }
  MFEM_VERIFY(pmesh->GetGlobalNE() == num_leaves,
              "Replayed refinement does not reproduce the nonconforming mesh!");

  // Target partition (that of the serial mesh) and curved mesh nodes of the local leaves,
  // from the node root.
  std::vector<std::int64_t> leaf_id;
  Locate(tree_pos, &leaf_id);
  std::vector<std::size_t>().swap(tree_pos);
  auto target =
      PullFromNodeRoot<int>(node_comm, leaf_id, [&](std::int64_t g, std::vector<int> &v)
                            { v.push_back(partitioning[serial_leaf[g]]); });
  if (data.curved)
  {
    auto values = PullFromNodeRoot<double>(
        node_comm, leaf_id, [&](std::int64_t g, std::vector<double> &v)
        { AppendElementNodes(*smesh, serial_leaf[g], v); });
    SetElementNodes(*pmesh, data, values);
    Mpi::Print(comm, " Copied the curved mesh nodes ({:.1f} s)\n", MPI_Wtime() - t0);
  }
  MPI_Comm_free(&node_comm);
  {
    mfem::Array<int> partition(target.data(), static_cast<int>(target.size()));
    pmesh->Rebalance(partition);
  }
  Mpi::Print(comm, " Rebalanced the nonconforming mesh ({:.1f} s)\n", MPI_Wtime() - t0);
  SetRootElementOrigin(*pmesh, {});  // Root elements in the input order
  return pmesh;
}

}  // namespace palace::mesh
