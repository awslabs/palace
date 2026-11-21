// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/geodata_impl.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

// Tests of the rebalancing of conformal meshes: the element graph gathered on the root, and
// RebalanceMesh, which moves the elements in parallel, and what it reports.

namespace palace
{
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

// The numbers of times the three timers of the rebalance were entered.
std::array<int, 3> RebalanceTimerCounts()
{
  const auto &timer = BlockTimer::GlobalTimer();
  return {static_cast<int>(timer.Counts(Timer::REBALANCE_GATHER)),
          static_cast<int>(timer.Counts(Timer::REBALANCE_PARTITION)),
          static_cast<int>(timer.Counts(Timer::REBALANCE_MOVE))};
}

// True if the result is the same on all processes.
bool SameOnAllProcesses(const mesh::RebalanceResult &result, MPI_Comm comm)
{
  double v[3] = {result.ratio_before, result.ratio_after, double(result.moved)};
  double vmin[3] = {v[0], v[1], v[2]}, vmax[3] = {v[0], v[1], v[2]};
  Mpi::GlobalMin(3, vmin, comm);
  Mpi::GlobalMax(3, vmax, comm);
  return std::equal(vmin, vmin + 3, vmax);
}

IoData MakeIoData()
{
  IoData iodata(Units(1.0, 1.0));
  iodata.model.refinement.maximum_imbalance = 1.0;
  return iodata;
}

// A mesh cracked along the plane x = 0.5: the elements and boundary elements beyond it get
// their own copies of the vertices on it, so that no face joins the two sides.
mfem::Mesh Crack(const mfem::Mesh &orig)
{
  const int nv = orig.GetNV();
  auto Beyond = [&](const mfem::Element &el)
  {
    double x = 0.0;
    for (int j = 0; j < el.GetNVertices(); j++)
    {
      x += orig.GetVertex(el.GetVertices()[j])[0];
    }
    return x / el.GetNVertices() > 0.5;
  };
  std::vector<int> copy(nv, -1), source(nv);
  std::iota(source.begin(), source.end(), 0);
  for (int e = 0; e < orig.GetNE(); e++)
  {
    const mfem::Element &el = *orig.GetElement(e);
    for (int j = 0; j < el.GetNVertices() && Beyond(el); j++)
    {
      const int v = el.GetVertices()[j];
      if (copy[v] < 0 && std::abs(orig.GetVertex(v)[0] - 0.5) < 1.0e-12)
      {
        copy[v] = static_cast<int>(source.size());
        source.push_back(v);
      }
    }
  }
  mfem::Mesh mesh(orig.Dimension(), static_cast<int>(source.size()), orig.GetNE(),
                  orig.GetNBE(), orig.SpaceDimension());
  for (int v : source)
  {
    mesh.AddVertex(orig.GetVertex(v));
  }
  auto Copy = [&](const mfem::Element &el)
  {
    mfem::Element *copied = el.Duplicate(&mesh);
    int *v = copied->GetVertices();
    for (int j = 0; j < copied->GetNVertices() && Beyond(el); j++)
    {
      v[j] = (copy[v[j]] >= 0) ? copy[v[j]] : v[j];
    }
    return copied;
  };
  for (int e = 0; e < orig.GetNE(); e++)
  {
    mesh.AddElement(Copy(*orig.GetElement(e)));
  }
  for (int b = 0; b < orig.GetNBE(); b++)
  {
    mesh.AddBdrElement(Copy(*orig.GetBdrElement(b)));
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  return mesh;
}

}  // namespace

TEST_CASE("GatherElementGraph conforming", "[rebalance][Serial][Parallel]")
{
  // The graph of a conforming mesh gathered on the root is the dual graph of the serial
  // mesh, in the global numbering: one vertex per element, one edge in each direction per
  // face between two elements, also across processes, the periodic boundary or none
  // across a crack.
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);
  for (std::string_view name : {"triangles", "tetrahedra", "periodic", "cracked"})
  {
    CAPTURE(name);
    mfem::Mesh serial;
    if (name == "triangles")
    {
      serial = MakeMesh(2);
    }
    else if (name == "tetrahedra")
    {
      serial = MakeMesh(3);
    }
    else if (name == "periodic")
    {
      auto orig = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::TETRAHEDRON);
      std::vector<mfem::Vector> translations = {mfem::Vector({1.0, 0.0, 0.0}),
                                                mfem::Vector({0.0, 1.0, 0.0})};
      serial =
          mfem::Mesh::MakePeriodic(orig, orig.CreatePeriodicVertexMapping(translations));
    }
    else
    {
      serial = Crack(mfem::Mesh::MakeCartesian3D(4, 2, 2, mfem::Element::TETRAHEDRON));
      CHECK(serial.GetNV() == 5 * 3 * 3 + 3 * 3);
    }
    std::vector<int> part = UnbalancedPartition(serial, size);
    mfem::ParMesh pmesh(comm, serial, part.data());
    std::vector<idx_t> xadj, adjncy;
    std::vector<int> rank_ne, rank_node;
    REQUIRE(mesh::GatherElementGraph(pmesh, xadj, adjncy, rank_ne, rank_node));
    if (!Mpi::Root(comm))
    {
      continue;
    }

    // The ParMesh constructor keeps the serial order of the elements on each process.
    std::vector<idx_t> global(serial.GetNE());
    std::vector<int> count(size, 0);
    idx_t g = 0;
    for (int p = 0; p < size; p++)
    {
      for (int e = 0; e < serial.GetNE(); e++)
      {
        if (part[e] == p)
        {
          global[e] = g++;
          count[p]++;
        }
      }
    }
    CHECK(rank_ne == count);
    std::vector<std::set<idx_t>> expected(serial.GetNE());
    int across = 0;
    for (int f = 0; f < serial.GetNumFaces(); f++)
    {
      int e1, e2;
      serial.GetFaceElements(f, &e1, &e2);
      if (e2 >= 0 && e1 != e2)
      {
        expected[global[e1]].insert(global[e2]);
        expected[global[e2]].insert(global[e1]);
        across += (part[e1] != part[e2]);
      }
    }
    CHECK((size == 1 || across > 0));
    REQUIRE(xadj.size() == expected.size() + 1);
    bool same = true;
    for (std::size_t v = 0; v < expected.size(); v++)
    {
      const std::vector<idx_t> row(adjncy.begin() + xadj[v], adjncy.begin() + xadj[v + 1]);
      same = same && (row == std::vector<idx_t>(expected[v].begin(), expected[v].end()));
    }
    CHECK(same);
  }
}

TEST_CASE("RebalanceMesh conformal", "[rebalance][Serial][Parallel]")
{
  // An unbalanced mesh, with the marks its refinements left as in the AMR loop, is
  // rebalanced in place: the same mesh object and the same mesh (counts, boundary elements
  // of each attribute, volume), balanced. The elements move with their vertex orders and
  // refinement flags, so a later refinement gives the mesh of the same refinement of the
  // serial mesh. Each of the three timers of the rebalance runs once. On one process
  // nothing moves.
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);
  for (int dim : {2, 3})
  {
    mfem::Mesh serial = MakeMesh(dim);
    RefineBall(serial, 0.3, 0.3);
    RefineBall(serial, 0.35, 0.2);
    std::vector<int> part = UnbalancedPartition(serial, size);
    const Counts counts = GetCounts(serial);
    auto pmesh = std::make_unique<mfem::ParMesh>(comm, serial, part.data());
    const mfem::ParMesh *object = pmesh.get();
    const auto before = GatherElements(*pmesh);
    const auto timers = RebalanceTimerCounts();
    const auto result = mesh::RebalanceMesh(MakeIoData(), pmesh);
    CHECK(pmesh.get() == object);
    CHECK(GetCounts(*pmesh) == counts);
    const auto after = GatherElements(*pmesh);
    CHECK(after == before);
    CHECK(SameOnAllProcesses(result, comm));
    CHECK(result.moved == (size > 1));
    if (size > 1)
    {
      CHECK(result.ratio_before > 1.0);
      CHECK(Imbalance(*pmesh) <= 1.2);
      CHECK(result.ratio_after == Imbalance(*pmesh));
      CHECK(RebalanceTimerCounts() ==
            std::array{timers[0] + 1, timers[1] + 1, timers[2] + 1});
      CHECK(mesh::RebalanceMessage(result, 1.0).starts_with("Rebalanced mesh:"));
    }
    else
    {
      CHECK(RebalanceTimerCounts() == timers);
    }

    // Refine the rebalanced mesh and the serial mesh the same way.
    RefineBall(*pmesh, 0.6, 0.25);
    RefineBall(serial, 0.6, 0.25);
    const Counts refined = GetCounts(*pmesh), serial_counts = GetCounts(serial);
    CHECK(refined == serial_counts);
    CHECK(std::abs(refined.volume - serial_counts.volume) <= 1.0e-12);
  }
}

TEST_CASE("RebalanceMesh conformal without a partition", "[rebalance][Serial][Parallel]")
{
  // One tetrahedron, on the first process: no partition gives every process an element, so
  // none comes back. The mesh keeps its distribution, and the result says so, the same on
  // all processes. Only the gather and the partition run, not the move.
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);
  mfem::Mesh serial(3, 4, 1, 4);
  const double x[4][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  for (const auto &v : x)
  {
    serial.AddVertex(v);
  }
  const int tet[4] = {0, 1, 2, 3},
            faces[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
  serial.AddTet(tet);
  for (const auto &f : faces)
  {
    serial.AddBdrTriangle(f);
  }
  serial.FinalizeTopology();
  serial.Finalize(true);
  std::vector<int> part = {0};
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, serial, part.data());
  const mfem::ParMesh *object = pmesh.get();
  const int ne = pmesh->GetNE();
  const auto before = GatherElements(*pmesh);
  const auto timers = RebalanceTimerCounts();
  const auto result = mesh::RebalanceMesh(MakeIoData(), pmesh);
  CHECK(pmesh.get() == object);
  CHECK(pmesh->GetNE() == ne);
  CHECK(GatherElements(*pmesh) == before);
  CHECK_FALSE(result.moved);
  CHECK(result.ratio_after == result.ratio_before);
  CHECK(SameOnAllProcesses(result, comm));
  if (size > 1)
  {
    CHECK(result.ratio_before > 1.0);
    CHECK(RebalanceTimerCounts() == std::array{timers[0] + 1, timers[1] + 1, timers[2]});
    const auto message = mesh::RebalanceMessage(result, 1.0);
    CHECK(message.starts_with("Mesh not rebalanced:"));
    CHECK(message.find("the distribution is kept") != std::string::npos);
  }
  else
  {
    CHECK(RebalanceTimerCounts() == timers);
  }
}

}  // namespace palace
