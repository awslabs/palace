// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <numbers>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/interpolator.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

using json = nlohmann::json;
namespace fs = std::filesystem;

namespace
{

// A tetrahedral mesh of the unit cube with its elements (and so its vertices) in a random
// order, written to a file in the given directory (one of its own for each test case and
// process: ctest runs the serial and parallel versions of the test cases at the same time).
fs::path WriteShuffledCube(const fs::path &dir, int n)
{
  auto mesh = mfem::Mesh::MakeCartesian3D(n, n, n, mfem::Element::TETRAHEDRON);
  mfem::Array<int> ordering(mesh.GetNE());
  std::iota(ordering.begin(), ordering.end(), 0);
  std::mt19937 gen(12345);
  std::shuffle(ordering.begin(), ordering.end(), gen);
  mesh.ReorderElements(ordering);
  const auto path = dir / ("palace-test-localreorder-" + std::to_string(n) + ".mesh");
  if (Mpi::Root(Mpi::World()))
  {
    std::ofstream fo(path);
    fo.precision(17);
    mesh.Print(fo);
  }
  Mpi::Barrier(Mpi::World());
  return path;
}

IoData MakeIoData(const fs::path &mesh_path, bool local_reorder, int max_its)
{
  json config = {{"Problem", {{"Type", "Eigenmode"}, {"Verbose", 0}}},
                 {"Model",
                  {{"Mesh", mesh_path.string()},
                   {"LocalReorder", local_reorder},
                   {"Refinement", {{"MaxIts", max_its}, {"MaximumImbalance", 1.0}}}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                 {"Boundaries", {{"PEC", {{"Attributes", {1, 2, 3, 4, 5, 6}}}}}},
                 {"Solver", {{"Eigenmode", {{"Target", 1.0}}}}}};
  return IoData(config, false);
}

// The centroids of the local elements, in their order.
std::vector<std::array<double, 3>> Centroids(mfem::ParMesh &mesh)
{
  std::vector<std::array<double, 3>> xc(mesh.GetNE());
  mfem::Vector c;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementCenter(e, c);
    xc[e] = {c(0), c(1), c(2)};
  }
  return xc;
}

// Number of processes whose sets of element centroids differ between two meshes.
int DifferentParts(mfem::ParMesh &a, mfem::ParMesh &b)
{
  auto Sorted = [](mfem::ParMesh &mesh)
  {
    auto xc = Centroids(mesh);
    for (auto &x : xc)
    {
      for (auto &v : x)
      {
        v = std::round(v * 1.0e9) * 1.0e-9;
      }
    }
    std::ranges::sort(xc);
    return xc;
  };
  int differ = (Sorted(a) != Sorted(b));
  Mpi::GlobalSum(1, &differ, a.GetComm());
  return differ;
}

// Sum of the distances between the centroids of consecutive local elements, for a
// nonconforming mesh only between consecutive elements of different root elements (the
// leaves of a root element are consecutive, in an order which does not depend on the
// order of the root elements).
double PathLength(mfem::ParMesh &mesh)
{
  const auto xc = Centroids(mesh);
  mfem::CoarseFineTransformations transforms;
  if (mesh.ncmesh)
  {
    mesh.ncmesh->GetRootTransforms(transforms);
  }
  double length = 0.0;
  for (std::size_t e = 1; e < xc.size(); e++)
  {
    if (!mesh.ncmesh ||
        transforms.embeddings[e].parent != transforms.embeddings[e - 1].parent)
    {
      length += std::hypot(xc[e][0] - xc[e - 1][0], xc[e][1] - xc[e - 1][1],
                           xc[e][2] - xc[e - 1][2]);
    }
  }
  Mpi::GlobalSum(1, &length, mesh.GetComm());
  return length;
}

// The squared norm in the mass matrix of a Nedelec space of order 2 of the interpolant of a
// smooth field, and the number of true dofs.
std::pair<double, HYPRE_BigInt> MassNorm(mfem::ParMesh &mesh)
{
  mfem::ND_FECollection fec(2, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &fec);
  mfem::VectorFunctionCoefficient f(3,
                                    [](const mfem::Vector &x, mfem::Vector &v)
                                    {
                                      v(0) = std::sin(3.0 * x(1)) + x(2);
                                      v(1) = std::cos(2.0 * x(0) * x(2));
                                      v(2) = x(0) * x(1) - x(2);
                                    });
  mfem::ParGridFunction gf(&fespace);
  gf.ProjectCoefficient(f);
  mfem::ParBilinearForm m(&fespace);
  m.AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
  m.Assemble();
  m.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> M(m.ParallelAssemble());
  std::unique_ptr<mfem::HypreParVector> x(gf.ParallelProject());
  mfem::HypreParVector y(*x);
  M->Mult(*x, y);
  return {mfem::InnerProduct(*x, y), fespace.GlobalTrueVSize()};
}

}  // namespace

TEST_CASE("HilbertOrder", "[localreorder][geodata][Serial]")
{
  // The cell centers of a 2^b grid, in a random order, are ordered along the Hilbert curve:
  // consecutive points are face neighbors.
  constexpr int b = 3, m = 1 << b;
  std::vector<std::array<int, 3>> cells;
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < m; j++)
    {
      for (int k = 0; k < m; k++)
      {
        cells.push_back({i, j, k});
      }
    }
  }
  std::mt19937 gen(7);
  std::shuffle(cells.begin(), cells.end(), gen);
  std::vector<double> x;
  for (const auto &c : cells)
  {
    for (int d = 0; d < 3; d++)
    {
      x.push_back(2.0 + 0.5 * (c[d] + 0.5));  // Arbitrary offset and scaling
    }
  }
  const int n = static_cast<int>(cells.size());
  const auto order = mesh::HilbertOrder(x.data(), n, 3);
  REQUIRE(static_cast<int>(order.size()) == n);
  std::vector<int> sorted(order);
  std::ranges::sort(sorted);
  std::vector<int> iota(n);
  std::iota(iota.begin(), iota.end(), 0);
  REQUIRE(sorted == iota);
  for (int k = 1; k < n; k++)
  {
    int dist = 0;
    for (int d = 0; d < 3; d++)
    {
      dist += std::abs(cells[order[k]][d] - cells[order[k - 1]][d]);
    }
    CHECK(dist == 1);
  }

  // Ties keep their order.
  const std::vector<double> same(12, 1.0);
  CHECK(mesh::HilbertOrder(same.data(), 4, 3) == std::vector<int>{0, 1, 2, 3});
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Local element reordering keeps the partition",
                 "[localreorder][geodata][Serial][Parallel]")
{
  const auto mesh_path = WriteShuffledCube(temp_dir, 6);
  for (const int max_its : {0, 1})
  {
    // Conforming distribution without adaptive mesh refinement, nonconforming with it.
    auto iodata_a = MakeIoData(mesh_path, false, max_its);
    auto iodata_b = MakeIoData(mesh_path, true, max_its);
    auto a = mesh::ReadMesh(iodata_a, Mpi::World());
    auto b = mesh::ReadMesh(iodata_b, Mpi::World());
    REQUIRE(a->Nonconforming() == (max_its > 0));
    REQUIRE(b->Nonconforming() == (max_its > 0));
    CHECK(a->GetGlobalNE() == b->GetGlobalNE());
    CHECK(DifferentParts(*a, *b) == 0);
    CHECK(PathLength(*b) < 0.5 * PathLength(*a));
    if (max_its == 0)
    {
      // The vertices keep their numbering (the orientations of the edges and faces).
      int differ = (a->GetNV() != b->GetNV());
      for (int v = 0; !differ && v < a->GetNV(); v++)
      {
        for (int d = 0; d < 3; d++)
        {
          differ |= (a->GetVertex(v)[d] != b->GetVertex(v)[d]);
        }
      }
      Mpi::GlobalSum(1, &differ, Mpi::World());
      CHECK(differ == 0);
    }
    {
      const auto [norm_a, size_a] = MassNorm(*a);
      const auto [norm_b, size_b] = MassNorm(*b);
      CHECK(size_a == size_b);
      CHECK(std::abs(norm_a - norm_b) <= 1.0e-12 * std::abs(norm_a));
    }
    if (max_its == 0)
    {
      continue;
    }

    // Nonconforming refinement and rebalancing give the same partition: the rebalancing
    // numbers the elements in the order of the input mesh.
    for (int it = 0; it < 2; it++)
    {
      for (auto *pmesh : {a.get(), b.get()})
      {
        mfem::Array<int> marked;
        mfem::Vector c;
        for (int e = 0; e < pmesh->GetNE(); e++)
        {
          pmesh->GetElementCenter(e, c);
          if (c(0) + 0.5 * c(1) < 0.4 + 0.2 * it)
          {
            marked.Append(e);
          }
        }
        pmesh->GeneralRefinement(marked, -1, 1);
      }
      CHECK(a->GetGlobalNE() == b->GetGlobalNE());
      CHECK(mesh::GetRootElementOrigin(*a) == nullptr);
      CHECK(mesh::GetRootElementOrigin(*b) != nullptr);
      mesh::RebalanceMesh(iodata_a, a);
      mesh::RebalanceMesh(iodata_b, b);
      CHECK(DifferentParts(*a, *b) == 0);
      CHECK(PathLength(*b) < 0.5 * PathLength(*a));
      const auto [norm_a, size_a] = MassNorm(*a);
      const auto [norm_b, size_b] = MassNorm(*b);
      CHECK(size_a == size_b);
      CHECK(std::abs(norm_a - norm_b) <= 1.0e-12 * std::abs(norm_a));
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Probes on element boundaries do not depend on the element order",
                 "[localreorder][Serial][Parallel]")
{
  // A piecewise constant field which identifies each element by its centroid, probed at a
  // vertex, an edge midpoint, a face centroid and an interior point of the shuffled cube:
  // with the reordering, a point on the boundary between elements takes the value of the
  // containing element first in the input (serial) mesh, whatever the partition.
  const int n = 6;
  const auto mesh_path = WriteShuffledCube(temp_dir, n);
  auto Field = [](const mfem::Vector &c)
  { return 1.0 + c(0) + std::numbers::sqrt2 * c(1) + std::numbers::pi * c(2); };
  const std::vector<std::array<double, 3>> points = {{3.0 / n, 2.0 / n, 4.0 / n},
                                                     {2.5 / n, 2.5 / n, 1.0 / n},
                                                     {1.0 / n, 3.4 / n, 4.2 / n},
                                                     {0.21, 0.37, 0.53}};

  // The expected values, from the serial mesh: the first element containing each point.
  std::vector<double> expected;
  {
    mfem::Mesh smesh(mesh_path.string());
    mfem::Vector c, x(3);
    for (const auto &pt : points)
    {
      std::copy(pt.begin(), pt.end(), x.begin());
      double value = 0.0;
      for (int e = 0; e < smesh.GetNE(); e++)
      {
        mfem::InverseElementTransformation inv(smesh.GetElementTransformation(e));
        inv.SetElementTol(1.0e-10);
        mfem::IntegrationPoint ip;
        if (inv.Transform(x, ip) == mfem::InverseElementTransformation::Inside)
        {
          smesh.GetElementCenter(e, c);
          value = Field(c);
          break;
        }
      }
      expected.push_back(value);
    }
  }
  std::map<int, config::ProbeData> probes;
  for (std::size_t k = 0; k < points.size(); k++)
  {
    probes[static_cast<int>(k) + 1].center = points[k];
  }
  const Units units(1.0, 1.0);
  for (const int max_its : {0, 1})
  {
    auto iodata = MakeIoData(mesh_path, true, max_its);
    Mesh mesh(mesh::ReadMesh(iodata, Mpi::World()));
    REQUIRE(!mesh::GetElementInputKeys(mesh.Get()).empty());
    mfem::L2_FECollection fec(0, 3);
    FiniteElementSpace fespace(mesh, &fec);
    GridFunction U(fespace);
    mfem::Vector c;
    mfem::Array<int> dofs;
    for (int e = 0; e < mesh.Get().GetNE(); e++)
    {
      mesh.Get().GetElementCenter(e, c);
      fespace.Get().GetElementDofs(e, dofs);
      U.Real()(dofs[0]) = Field(c);
    }
    InterpolationOperator op(probes, units, fespace);
    const auto values = op.ProbeField(U);
    REQUIRE(values.size() == points.size());
    for (std::size_t k = 0; k < points.size(); k++)
    {
      CHECK(std::abs(values[k].real() - expected[k]) <= 1.0e-12 * std::abs(expected[k]));
    }
  }
}

}  // namespace palace
