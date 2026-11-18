// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

// Tests of the rebalancing of conformal meshes: the element graph, the PT-Scotch partition,
// and RebalanceMesh with the gather and the parallel methods.

namespace palace
{
using json = nlohmann::json;

namespace
{

// Unit cube of tetrahedra (or unit square of triangles) with two element attributes, the
// exterior boundary attributes of MakeCartesian, and interior boundary elements (attribute
// 7) on the faces between the two attributes, every other one oriented towards the second.
mfem::Mesh MakeMesh(int dim)
{
  mfem::Mesh mesh = (dim == 3)
                        ? mfem::Mesh::MakeCartesian3D(4, 3, 2, mfem::Element::TETRAHEDRON)
                        : mfem::Mesh::MakeCartesian2D(6, 5, mfem::Element::TRIANGLE, true);
  mfem::Vector c;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementCenter(e, c);
    mesh.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
  }
  mesh.SetAttributes();
  int count = 0;
  for (int f = 0; f < mesh.GetNumFaces(); f++)
  {
    int e1, e2;
    mesh.GetFaceElements(f, &e1, &e2);
    if (e2 < 0 || mesh.GetAttribute(e1) == mesh.GetAttribute(e2))
    {
      continue;
    }
    mfem::Array<int> v;
    mesh.GetFaceVertices(f, v);
    if (count++ % 2)
    {
      std::swap(v[0], v[1]);
    }
    mfem::Element *be =
        mesh.NewElement((dim == 3) ? mfem::Geometry::TRIANGLE : mfem::Geometry::SEGMENT);
    be->SetVertices(v);
    be->SetAttribute(7);
    mesh.AddBdrElement(be);
  }
  mesh.FinalizeTopology(false);
  mesh.Finalize(true);
  mesh.SetAttributes();
  return mesh;
}

// Refine the elements with a center within `r` of `p`, conformally.
template <typename MeshType>
void RefineBall(MeshType &mesh, double x, double r)
{
  mfem::Array<int> marked;
  mfem::Vector c;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementCenter(e, c);
    double d2 = 0.0;
    for (int i = 0; i < c.Size(); i++)
    {
      d2 += (c(i) - x) * (c(i) - x);
    }
    if (d2 < r * r)
    {
      marked.Append(e);
    }
  }
  mesh.GeneralRefinement(marked, 0);
}

// Global counts of a parallel mesh: elements, vertices, edges, faces (3D), boundary
// elements of each attribute, volume.
struct Counts
{
  long long ne, nv, nedges, nfaces;
  std::vector<long long> nbe;
  double volume;
};

bool operator==(const Counts &a, const Counts &b)
{
  return a.ne == b.ne && a.nv == b.nv && a.nedges == b.nedges && a.nfaces == b.nfaces &&
         a.nbe == b.nbe && std::abs(a.volume - b.volume) <= 1.0e-12 * std::abs(b.volume);
}

Counts GetCounts(mfem::ParMesh &mesh)
{
  const int dim = mesh.Dimension();
  mfem::H1_FECollection h1(1, dim);
  mfem::ND_FECollection nd(1, dim);
  mfem::RT_FECollection rt(0, dim);
  Counts c;
  c.ne = mesh.GetGlobalNE();
  c.nv = mfem::ParFiniteElementSpace(&mesh, &h1).GlobalTrueVSize();
  c.nedges = mfem::ParFiniteElementSpace(&mesh, &nd).GlobalTrueVSize();
  c.nfaces = (dim == 3) ? mfem::ParFiniteElementSpace(&mesh, &rt).GlobalTrueVSize() : 0;
  c.nbe.assign(mesh.bdr_attributes.Max() + 1, 0);
  for (int b = 0; b < mesh.GetNBE(); b++)
  {
    c.nbe[mesh.GetBdrAttribute(b)]++;
  }
  Mpi::GlobalSum(static_cast<int>(c.nbe.size()), c.nbe.data(), mesh.GetComm());
  c.volume = 0.0;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    c.volume += mesh.GetElementVolume(e);
  }
  Mpi::GlobalSum(1, &c.volume, mesh.GetComm());
  return c;
}

Counts GetCounts(mfem::Mesh &mesh)
{
  const int dim = mesh.Dimension();
  Counts c;
  c.ne = mesh.GetNE();
  c.nv = mesh.GetNV();
  c.nedges = mesh.GetNEdges();
  c.nfaces = (dim == 3) ? mesh.GetNFaces() : 0;
  c.nbe.assign(mesh.bdr_attributes.Max() + 1, 0);
  for (int b = 0; b < mesh.GetNBE(); b++)
  {
    c.nbe[mesh.GetBdrAttribute(b)]++;
  }
  c.volume = 0.0;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    c.volume += mesh.GetElementVolume(e);
  }
  return c;
}

// All the elements of a parallel mesh, on the root process, sorted: attribute, refinement
// flag (tetrahedra) and vertex coordinates in the vertex order of the element, which gives
// its refinement edge.
std::vector<std::vector<double>> GatherElements(mfem::ParMesh &mesh)
{
  std::vector<double> local;
  const int nv = mesh.Dimension() + 1, sdim = mesh.SpaceDimension();
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const mfem::Element *el = mesh.GetElement(e);
    local.push_back(el->GetAttribute());
    local.push_back((el->GetType() == mfem::Element::TETRAHEDRON)
                        ? static_cast<const mfem::Tetrahedron *>(el)->GetRefinementFlag()
                        : 0);
    for (int j = 0; j < nv; j++)
    {
      for (int d = 0; d < sdim; d++)
      {
        local.push_back(std::round(mesh.GetVertex(el->GetVertices()[j])[d] * 1.0e9));
      }
    }
  }
  MPI_Comm comm = mesh.GetComm();
  const int size = Mpi::Size(comm), n = static_cast<int>(local.size());
  std::vector<int> counts(size), displs(size + 1, 0);
  MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
  for (int i = 0; i < size; i++)
  {
    displs[i + 1] = displs[i] + counts[i];
  }
  std::vector<double> all(Mpi::Root(comm) ? displs[size] : 0);
  MPI_Gatherv(local.data(), n, MPI_DOUBLE, all.data(), counts.data(), displs.data(),
              MPI_DOUBLE, 0, comm);
  const std::size_t m = 2 + nv * sdim;
  std::vector<std::vector<double>> elems;
  for (std::size_t i = 0; i < all.size(); i += m)
  {
    elems.emplace_back(all.begin() + i, all.begin() + i + m);
  }
  std::sort(elems.begin(), elems.end());
  return elems;
}

// The same elements, with the vertex order and the refinement flag forgotten.
std::vector<std::vector<double>> WithoutOrder(const std::vector<std::vector<double>> &elems,
                                              int dim)
{
  std::vector<std::vector<double>> out;
  const int nv = dim + 1,
            sdim = static_cast<int>((elems.empty() ? 2 : elems[0].size() - 2) / nv);
  for (const auto &x : elems)
  {
    std::vector<std::vector<double>> verts;
    for (int j = 0; j < nv; j++)
    {
      verts.emplace_back(x.begin() + 2 + j * sdim, x.begin() + 2 + (j + 1) * sdim);
    }
    std::sort(verts.begin(), verts.end());
    std::vector<double> y = {x[0]};
    for (const auto &v : verts)
    {
      y.insert(y.end(), v.begin(), v.end());
    }
    out.push_back(y);
  }
  std::sort(out.begin(), out.end());
  return out;
}

// The partition of a serial mesh which puts most elements on the first process.
std::vector<int> UnbalancedPartition(mfem::Mesh &mesh, int size)
{
  std::vector<int> part(mesh.GetNE(), 0);
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    part[e] = (e % 7 == 0) ? (e / 7) % size : 0;
  }
  return part;
}

double Imbalance(const mfem::ParMesh &mesh)
{
  int min = mesh.GetNE(), max = mesh.GetNE();
  Mpi::GlobalMin(1, &min, mesh.GetComm());
  Mpi::GlobalMax(1, &max, mesh.GetComm());
  return static_cast<double>(max) / min;
}

IoData MakeIoData(RebalanceMethod method)
{
  IoData iodata(Units(1.0, 1.0));
  iodata.model.refinement.maximum_imbalance = 1.0;
  iodata.model.refinement.rebalance_method = method;
  return iodata;
}

}  // namespace

TEST_CASE("ElementGraph", "[rebalance][Serial][Parallel]")
{
  // The element graph of a distributed mesh is the dual graph of the serial mesh: one
  // vertex per element in the order of the processes, one edge per interior face.
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);
  for (int dim : {2, 3})
  {
    mfem::Mesh serial = MakeMesh(dim);
    std::vector<int> part = UnbalancedPartition(serial, size);
    mfem::ParMesh mesh(comm, serial, part.data());
    std::vector<long long> xadj, adjncy;
    const long long offset = mesh::GetElementGraph(mesh, xadj, adjncy);
    REQUIRE(xadj.size() == static_cast<std::size_t>(mesh.GetNE() + 1));
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      CHECK(std::is_sorted(adjncy.begin() + xadj[e], adjncy.begin() + xadj[e + 1]));
    }

    // Global number of each serial element: the ParMesh constructor keeps the serial order
    // on each process.
    std::vector<long long> global(serial.GetNE());
    {
      long long g = 0;
      for (int p = 0; p < size; p++)
      {
        for (int e = 0; e < serial.GetNE(); e++)
        {
          if (part[e] == p)
          {
            global[e] = g++;
          }
        }
      }
      long long first = -1;
      for (int e = 0; e < serial.GetNE(); e++)
      {
        if (part[e] == Mpi::Rank(comm) && first < 0)
        {
          first = global[e];
        }
      }
      CHECK((mesh.GetNE() == 0 || first == offset));
    }
    std::vector<int> serial_of_global(serial.GetNE());
    for (int e = 0; e < serial.GetNE(); e++)
    {
      serial_of_global[global[e]] = e;
    }
    bool ok = true;
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      const int se = serial_of_global[offset + e];
      std::vector<long long> expected;
      for (int f = 0; f < serial.GetNumFaces(); f++)
      {
        int e1, e2;
        serial.GetFaceElements(f, &e1, &e2);
        if (e2 >= 0 && (e1 == se || e2 == se))
        {
          expected.push_back(global[(e1 == se) ? e2 : e1]);
        }
      }
      std::sort(expected.begin(), expected.end());
      ok = ok &&
           std::equal(expected.begin(), expected.end(), adjncy.begin() + xadj[e],
                      adjncy.begin() + xadj[e + 1]) &&
           static_cast<long long>(expected.size()) == xadj[e + 1] - xadj[e];
    }
    CHECK(ok);
  }
}

TEST_CASE("RebalancePartitioning", "[rebalance][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);
  mfem::Mesh serial = MakeMesh(3);
  RefineBall(serial, 0.3, 0.3);
  std::vector<int> part = UnbalancedPartition(serial, size);
  mfem::ParMesh mesh(comm, serial, part.data());
  mfem::Array<int> partitioning;
#if defined(PALACE_WITH_SCOTCH)
  // Default groups of partitioning processes, one partitioning process, all processes.
  for (long long elements_per_process : {100000LL, 1000000000LL, 0LL})
  {
    REQUIRE(mesh::GetRebalancePartitioning(mesh, partitioning, elements_per_process));
    REQUIRE(partitioning.Size() == mesh.GetNE());
    std::vector<int> counts(size, 0);
    for (int p : partitioning)
    {
      REQUIRE((p >= 0 && p < size));
      counts[p]++;
    }
    Mpi::GlobalSum(size, counts.data(), comm);
    const long long ne = mesh.GetGlobalNE();
    CHECK(*std::min_element(counts.begin(), counts.end()) > 0);
    CHECK(*std::max_element(counts.begin(), counts.end()) <=
          std::ceil(1.05 * static_cast<double>(ne) / size));

    // The partition does not depend on the run.
    mfem::Array<int> again;
    REQUIRE(mesh::GetRebalancePartitioning(mesh, again, elements_per_process));
    bool same = (again.Size() == partitioning.Size());
    for (int e = 0; same && e < again.Size(); e++)
    {
      same = (again[e] == partitioning[e]);
    }
    CHECK(same);
  }

  // A partition leaving a process without elements is refused: a single tetrahedron.
  if (size > 1)
  {
    mfem::Mesh single(3, 4, 1, 0, 3);
    const double v[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (const auto &x : v)
    {
      single.AddVertex(x);
    }
    const int tet[4] = {0, 1, 2, 3};
    single.AddTet(tet);
    single.FinalizeTetMesh(1, 1, true);
    std::vector<int> part1(1, 0);
    mfem::ParMesh one(comm, single, part1.data());
    CHECK_FALSE(mesh::GetRebalancePartitioning(one, partitioning));
  }
#else
  CHECK_FALSE(mesh::GetRebalancePartitioning(mesh, partitioning));
#endif
}

TEST_CASE("GetGlobalEntityCounts", "[rebalance][Serial][Parallel]")
{
  // The counts from the shared entity groups are those of the serial mesh, and the true
  // sizes of the lowest-order H1, ND and RT spaces, whatever the distribution.
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);
  for (int dim : {2, 3})
  {
    mfem::Mesh serial = MakeMesh(dim);
    RefineBall(serial, 0.3, 0.3);
    const Counts counts = GetCounts(serial);
    std::vector<long long> expected = {counts.ne, counts.nv, counts.nedges, counts.nfaces};
    expected.insert(expected.end(), counts.nbe.begin() + 1, counts.nbe.end());
    std::vector<int> cyclic(serial.GetNE());
    for (int e = 0; e < serial.GetNE(); e++)
    {
      cyclic[e] = e % size;
    }
    std::vector<std::vector<int>> parts = {UnbalancedPartition(serial, size), cyclic};
    for (auto &part : parts)
    {
      mfem::ParMesh pmesh(comm, serial, part.data());
      CHECK(mesh::GetGlobalEntityCounts(pmesh) == expected);
      const Counts pcounts = GetCounts(pmesh);
      CHECK(pcounts.nv == expected[1]);
      CHECK(pcounts.nedges == expected[2]);
      CHECK(pcounts.nfaces == expected[3]);
    }
  }
}

TEST_CASE("RebalanceMesh conformal", "[rebalance][Serial][Parallel]")
{
  // The same unbalanced mesh, rebalanced with the gather method and with the parallel
  // method: same mesh (counts, boundary elements of each attribute, volume), balanced. Both
  // methods mark the rebalanced mesh again: its vertex orders and refinement flags are
  // those of the serial mesh marked again, and a later refinement gives the mesh of the
  // same refinement of that serial mesh. On one process there is no rebalance.
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);
  for (int dim : {2, 3})
  {
    mfem::Mesh serial = MakeMesh(dim);
    RefineBall(serial, 0.3, 0.3);
    RefineBall(serial, 0.35, 0.2);
    std::vector<int> part = UnbalancedPartition(serial, size);
    const Counts counts = GetCounts(serial);
    mfem::Mesh marked(serial);
    marked.Finalize(true);
    std::vector<std::vector<double>> marked_elements;
    {
      mfem::ParMesh pmarked(comm, marked, part.data());
      marked_elements = GatherElements(pmarked);
    }
    std::vector<RebalanceMethod> methods = {RebalanceMethod::GATHER};
#if defined(PALACE_WITH_SCOTCH)
    methods.push_back(RebalanceMethod::PARALLEL);
#endif
    for (auto method : methods)
    {
      auto pmesh = std::make_unique<mfem::ParMesh>(comm, serial, part.data());
      const auto before = GatherElements(*pmesh);
      const double ratio = mesh::RebalanceMesh(MakeIoData(method), pmesh);
      CHECK(GetCounts(*pmesh) == counts);
      const auto after = GatherElements(*pmesh);
      CHECK(WithoutOrder(after, dim) == WithoutOrder(before, dim));
      if (size > 1)
      {
        CHECK(ratio > 1.0);
        CHECK(Imbalance(*pmesh) <= 1.2);
      }
      if (size > 1)
      {
        // The vertex orders and refinement flags of a new marking, which differ in 3D from
        // those the mesh had.
        CHECK(after == marked_elements);
        if (dim == 3 && Mpi::Root(comm))
        {
          CHECK(after != before);
        }
      }
      else
      {
        CHECK(after == before);
      }

      // Refine the rebalanced mesh, and the serial mesh (marked again, if the mesh was
      // rebalanced) the same way.
      mfem::Mesh serial_refined((size > 1) ? marked : serial);
      RefineBall(*pmesh, 0.6, 0.25);
      RefineBall(serial_refined, 0.6, 0.25);
      const Counts refined = GetCounts(*pmesh), serial_counts = GetCounts(serial_refined);
      CHECK(refined == serial_counts);
      CHECK(std::abs(refined.volume - serial_counts.volume) <= 1.0e-12);
    }
  }
}

TEST_CASE("RebalanceMesh conformal fallback", "[rebalance][Serial][Parallel]")
{
  // A conformal hexahedral mesh is rebalanced on the root process whatever the method.
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);
  mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(8, 6, 4, mfem::Element::HEXAHEDRON);
  std::vector<int> part(serial.GetNE(), 0);
  part.back() = size - 1;
  const Counts counts = GetCounts(serial);
  for (auto method : {RebalanceMethod::GATHER, RebalanceMethod::PARALLEL})
  {
    auto pmesh = std::make_unique<mfem::ParMesh>(comm, serial, part.data());
    mesh::RebalanceMesh(MakeIoData(method), pmesh);
    CHECK(GetCounts(*pmesh) == counts);
    if (size > 1)
    {
      CHECK(Imbalance(*pmesh) <= 1.2);
    }
  }
}

}  // namespace palace
