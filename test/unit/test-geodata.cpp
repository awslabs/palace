// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"
#include "test-helpers.hpp"

#include "fem/fespace.hpp"
#include "utils/geodata.hpp"
#include "utils/geodata_impl.hpp"
#include "utils/iodata.hpp"

#include "fem/interpolator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/edgedistance.hpp"
#include "utils/filesystem.hpp"
#include "utils/metaledge.hpp"

namespace palace
{
using json = nlohmann::json;
namespace fs = std::filesystem;
using namespace Catch::Matchers;

namespace
{

std::unique_ptr<mfem::Mesh> TwoEdgeFanTetMesh()
{
  auto mesh = std::make_unique<mfem::Mesh>(3, 10, 4, 0, 3);

  mesh->AddVertex(0.0, 0.0, 0.0);
  mesh->AddVertex(1.0, 0.0, 0.0);
  mesh->AddVertex(0.0, 1.0, 0.0);
  mesh->AddVertex(0.0, 0.0, 1.0);
  mesh->AddVertex(0.0, 0.0, -1.0);
  mesh->AddVertex(10.0, 0.0, 0.0);
  mesh->AddVertex(11.0, 0.0, 0.0);
  mesh->AddVertex(10.0, 1.0, 0.0);
  mesh->AddVertex(10.0, 0.0, 1.0);
  mesh->AddVertex(10.0, 0.0, -1.0);

  mesh->AddTet(0, 1, 2, 3, 1);
  mesh->AddTet(0, 2, 1, 4, 1);
  mesh->AddTet(5, 6, 7, 8, 1);
  mesh->AddTet(5, 7, 6, 9, 1);
  mesh->FinalizeTopology();
  return mesh;
}

std::unique_ptr<mfem::Mesh> CrackedSquareSheetTetMesh(bool add_truncation = false)
{
  auto mesh = std::make_unique<mfem::Mesh>(3, 8, 8, 8, 3);

  mesh->AddVertex(0.0, 0.0, 0.0);
  mesh->AddVertex(1.0, 0.0, 0.0);
  mesh->AddVertex(1.0, 1.0, 0.0);
  mesh->AddVertex(0.0, 1.0, 0.0);
  mesh->AddVertex(0.5, 0.5, 0.0);
  mesh->AddVertex(0.5, 0.5, 0.0);
  mesh->AddVertex(0.5, 0.5, 1.0);
  mesh->AddVertex(0.5, 0.5, -1.0);

  for (int i = 0; i < 4; i++)
  {
    const int j = (i + 1) % 4;
    mesh->AddTet(i, j, 4, 6, 2);
    mesh->AddTet(i, 5, j, 7, 1);
    mesh->AddBdrTriangle(i, j, 4, 5);
    mesh->AddBdrTriangle(i, 5, j, 5);
    if (add_truncation && i == 0)
    {
      mesh->AddBdrTriangle(i, j, 6, 6);
      mesh->AddBdrTriangle(i, 7, j, 6);
    }
  }
  mesh->FinalizeTopology(false);
  return mesh;
}

std::unique_ptr<mfem::ParMesh> RoundedIslandSheetHexMesh(int in_plane_elements,
                                                         bool high_order, bool rounded)
{
  constexpr double center = 0.5;
  constexpr double half_width = 0.25;
  constexpr double radius = 0.125;
  constexpr double tolerance = 1.0e-12;

  mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(in_plane_elements, 2, in_plane_elements,
                                                  mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  for (int face = 0; face < serial.GetNumFaces(); face++)
  {
    int element1, element2;
    serial.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    serial.GetFaceVertices(face, vertices);
    bool on_plane = true;
    double xmin = 1.0, xmax = 0.0, zmin = 1.0, zmax = 0.0;
    for (const int vertex : vertices)
    {
      const double *point = serial.GetVertex(vertex);
      on_plane = on_plane && std::abs(point[1] - center) < tolerance;
      xmin = std::min(xmin, point[0]);
      xmax = std::max(xmax, point[0]);
      zmin = std::min(zmin, point[2]);
      zmax = std::max(zmax, point[2]);
    }
    if (on_plane && xmin >= center - half_width - tolerance &&
        xmax <= center + half_width + tolerance &&
        zmin >= center - half_width - tolerance && zmax <= center + half_width + tolerance)
    {
      serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
      serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
    }
  }
  serial.FinalizeTopology();
  serial.Finalize();

  if (high_order)
  {
    serial.SetCurvature(2);
  }
  if (rounded)
  {
    serial.Transform(
        [&](const mfem::Vector &input, mfem::Vector &output)
        {
          output = input;
          if (std::abs(input[1] - center) > tolerance)
          {
            return;
          }
          for (const double sign_x : {-1.0, 1.0})
          {
            for (const double sign_z : {-1.0, 1.0})
            {
              const double corner_x = center + sign_x * half_width;
              const double corner_z = center + sign_z * half_width;
              const double circle_x = corner_x - sign_x * radius;
              const double circle_z = corner_z - sign_z * radius;
              const double local_x = sign_x * (input[0] - circle_x);
              const double local_z = sign_z * (input[2] - circle_z);
              if (local_x < -tolerance || local_x > radius + tolerance ||
                  local_z < -tolerance || local_z > radius + tolerance)
              {
                continue;
              }

              double angle;
              if (std::abs(input[2] - corner_z) <= tolerance)
              {
                angle = 0.5 * std::acos(-1.0) - 0.25 * std::acos(-1.0) * local_x / radius;
              }
              else if (std::abs(input[0] - corner_x) <= tolerance)
              {
                angle = 0.25 * std::acos(-1.0) * local_z / radius;
              }
              else
              {
                continue;
              }
              output[0] = circle_x + sign_x * radius * std::cos(angle);
              output[2] = circle_z + sign_z * radius * std::sin(angle);
              return;
            }
          }
        });
    if (high_order)
    {
      for (int d = 0; d < serial.SpaceDimension(); d++)
      {
        mfem::Vector values;
        serial.GetNodes()->GetNodalValues(values, d + 1);
        for (int vertex = 0; vertex < serial.GetNV(); vertex++)
        {
          serial.GetVertex(vertex)[d] = values[vertex];
        }
      }
    }
  }
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

bool ElementContainsEdge(const mfem::Element &el, int v0, int v1)
{
  bool has_v0 = false, has_v1 = false;
  const int *verts = el.GetVertices();
  for (int i = 0; i < el.GetNVertices(); i++)
  {
    has_v0 = has_v0 || verts[i] == v0;
    has_v1 = has_v1 || verts[i] == v1;
  }
  return has_v0 && has_v1;
}

void CheckVertex(const mfem::Mesh &mesh, int v, const std::array<double, 3> &coord)
{
  const double *actual = mesh.GetVertex(v);
  for (int d = 0; d < 3; d++)
  {
    CAPTURE(v, d);
    CHECK_THAT(actual[d], WithinAbs(coord[d], 1.0e-12));
  }
}

}  // namespace

// TODO: Add this test when we can access MFEM_DATA_PATH from Spack
// This requires MFEM to move to a CMake-based build system for spack.
// Waiting on https://github.com/spack/spack-packages/pull/143
//
// TEST_CASE("ReferenceTetrahedron", "[geodata][Serial]")
// {
//   auto tet1 = SingleTetMesh();
//   auto ref_tet_path = fs::path(MFEM_DATA_PATH) / "ref-tetrahedron.mesh";
//   mfem::Mesh tet2(ref_tet_path.string());

//   REQUIRE(tet1.GetNV() == tet2.GetNV());
//   REQUIRE(tet1.GetNE() == tet2.GetNE());
//   for (int i = 0; i < tet1.GetNV(); i++)
//   {
//     REQUIRE_THAT(tet1.GetVertex(i)[0], WithinAbs(tet2.GetVertex(i)[0], 1.0e-12));
//     REQUIRE_THAT(tet1.GetVertex(i)[1], WithinAbs(tet2.GetVertex(i)[1], 1.0e-12));
//     REQUIRE_THAT(tet1.GetVertex(i)[2], WithinAbs(tet2.GetVertex(i)[2], 1.0e-12));
//   }
//   for (int i = 0; i < tet1.GetNE(); i++)
//   {
//     for (int j = 0; j < 4; j++)
//     {
//       CHECK(tet1.GetElement(i)->GetVertices()[j] ==
//       tet2.GetElement(i)->GetVertices()[j]);
//     }
//   }
//   for (int i = 0; i < tet1.GetNBE(); i++)
//   {
//     for (int j = 0; j < 3; j++)
//     {
//       CHECK(tet1.GetBdrElement(i)->GetVertices()[j] ==
//               tet2.GetBdrElement(i)->GetVertices()[j]);
//     }
//   }
// }

TEST_CASE("TetToHex", "[geodata][Serial]")
{

  mfem::Mesh single_tet(SingleTetMesh());

  SECTION("Linear")
  {
    int order = 1;
    single_tet.EnsureNodes();
    single_tet.SetCurvature(order);
    auto four_hex = mesh::MeshTetToHex(single_tet);
    CHECK(four_hex.GetNE() == 4);

    // DOFs are added in vert -> edge -> face order, based on local ordering in vertex,
    // which has dofs (3,2,1,0).
    const std::vector<std::array<double, 3>> global_dof_vals{
        {0.0, 0.0, 0.0},              // in elem 3
        {1.0, 0.0, 0.0},              // in elem 2
        {0.0, 1.0, 0.0},              // in elem 1
        {0.0, 0.0, 1.0},              // in elem 0
        {0.0, 0.5, 0.5},              // 3 -> 2
        {0.5, 0.0, 0.5},              // 3 -> 1
        {0.0, 0.0, 0.5},              // 3 -> 0
        {0.5, 0.5, 0.0},              // 2 -> 1
        {0.0, 0.5, 0.0},              // 2 -> 0
        {0.5, 0.0, 0.0},              // 1 -> 0
        {1.0 / 3, 1.0 / 3, 0.0},      // opp 3
        {1.0 / 3, 0.0, 1.0 / 3},      // opp 2
        {0.0, 1.0 / 3, 1.0 / 3},      // opp 1
        {1.0 / 3, 1.0 / 3, 1.0 / 3},  // opp 0
        {0.25, 0.25, 0.25}};
    // From drawing out dof diagram.
    const std::vector<std::array<int, 8>> elem_dofs{{3, 4, 13, 5, 6, 12, 14, 11},
                                                    {2, 7, 13, 4, 8, 10, 14, 12},
                                                    {1, 5, 13, 7, 9, 11, 14, 10},
                                                    {0, 9, 10, 8, 6, 11, 14, 12}};

    for (int i = 0; i < 14; i++)
      for (int j = 0; j < 3; j++)
      {
        // margin(1e-12) for comparing zeros.
        CHECK_THAT((*four_hex.GetNodes())(j + 3 * i),
                   WithinAbs(global_dof_vals[i][j], 1e-12));
      }

    mfem::Vector vdof_vals, col;
    four_hex.GetNodes()->GetElementDofValues(0, vdof_vals);
    REQUIRE(four_hex.GetNodes()->FESpace()->GetOrdering() == mfem::Ordering::byVDIM);
    REQUIRE(vdof_vals.Size() == 3 * 8);
    mfem::DenseMatrix vdof_vals_mat(vdof_vals.GetData(), 8, 3);

    auto check_mat =
        [&col, &global_dof_vals, &vdof_vals_mat](const std::array<int, 8> &verts)
    {
      for (std::size_t i = 0; i < 3; i++)
      {
        vdof_vals_mat.GetColumn(i, col);
        for (int j = 0; j < col.Size(); j++)
        {
          CAPTURE(i, j, global_dof_vals[verts[j]][i], col(j));
          // margin(1e-12) for comparing zeros.
          CHECK_THAT(col(j), WithinAbs(global_dof_vals[verts[j]][i], 1e-12));
        }
      }
    };

    for (int i = 0; i < 4; i++)
    {
      four_hex.GetNodes()->GetElementDofValues(i, vdof_vals);
      check_mat(elem_dofs[i]);
    }
  }
#if defined(PALACE_WITH_GSLIB)
  SECTION("UniformSampler")
  {
    // Use GSLIB to find a mapping for the single tet mesh, with some randomly perturbed
    // position data.
    int order = GENERATE(2, 3);
    const int sdim = 3;
    // Create linear meshes, to copy the nodes to.
    mfem::Mesh linear_single_tet(single_tet);
    auto linear_four_hex = mesh::MeshTetToHex(linear_single_tet);
    linear_single_tet.EnsureNodes();
    linear_four_hex.EnsureNodes();
    REQUIRE(linear_single_tet.GetNodes());
    REQUIRE(linear_single_tet.GetNodes()->FESpace()->GetMaxElementOrder() == 1);
    REQUIRE(linear_four_hex.GetNodes());
    REQUIRE(linear_four_hex.GetNodes()->FESpace()->GetMaxElementOrder() == 1);

    single_tet.EnsureNodes();
    single_tet.SetCurvature(order);

    // Randomly perturb the non-vertex data with positive values (ensures non-zeros in later
    // comparison).
    for (int i = 0; i < single_tet.GetNodes()->Size(); i++)
    {
      (*single_tet.GetNodes())(i) += 0.05 * (1.0 + static_cast<double>(rand()) / RAND_MAX);
    }

    auto four_hex = mesh::MeshTetToHex(single_tet);
    REQUIRE(four_hex.GetNE() == 4);

    // Helper to generate a set of locations to sample tet at.
    // n_sample is number of points along an edge, thus >= 2.
    // In byVDIM ordering [x1,y1,z1,x2,y2,z2,...]
    auto gen_samples = [](int n_sample)
    {
      REQUIRE(n_sample >= 2);
      mfem::Vector xyz_samples(3 * (n_sample * (n_sample + 1) * (n_sample + 2) / 6));
      int o = 0;
      for (double k = 0; k < n_sample; k++)
        for (double j = 0; j < n_sample - k; j++)
          for (double i = 0; i < n_sample - j - k; i++)
          {
            xyz_samples(o) = i / (n_sample - 1);
            xyz_samples(o + 1) = j / (n_sample - 1);
            xyz_samples(o + 2) = k / (n_sample - 1);
            o++;
          }
      return xyz_samples;
    };

    // Uniform sampling over the tet as "xyz" coords of the original reference tet.
    auto xyz_samples = gen_samples(order + 2);

    // Create FiniteElementSpace on the linear meshes, with the same dofs from the higher
    // mesh nodes, then sample the node functions using coordinates from the linear meshes.
    // These should be equal to each other, as the sample points correspond to the original
    // reference space on the tet.
    const auto &tet_FESpace = single_tet.GetNodes()->FESpace();
    const auto &hex_FESpace = four_hex.GetNodes()->FESpace();
    mfem::FiniteElementSpace linear_tet_FESpace(&linear_single_tet, tet_FESpace->FEColl(),
                                                sdim, tet_FESpace->GetOrdering());
    mfem::FiniteElementSpace linear_hex_FESpace(&linear_four_hex, hex_FESpace->FEColl(),
                                                sdim, hex_FESpace->GetOrdering());
    mfem::GridFunction tet_nodes_on_linear_tet(&linear_tet_FESpace);
    mfem::GridFunction hex_nodes_on_linear_hex(&linear_hex_FESpace);
    REQUIRE(tet_nodes_on_linear_tet.Size() == single_tet.GetNodes()->Size());
    REQUIRE(hex_nodes_on_linear_hex.Size() == four_hex.GetNodes()->Size());
    tet_nodes_on_linear_tet = *single_tet.GetNodes();
    hex_nodes_on_linear_hex = *four_hex.GetNodes();

    mfem::Vector tet_vals(xyz_samples.Size()), hex_vals(xyz_samples.Size());
    fem::InterpolateFunction(xyz_samples, tet_nodes_on_linear_tet, tet_vals,
                             tet_FESpace->GetOrdering());
    fem::InterpolateFunction(xyz_samples, hex_nodes_on_linear_hex, hex_vals,
                             hex_FESpace->GetOrdering());
    const auto *T = tet_vals.HostRead();
    const auto *H = hex_vals.HostRead();
    for (int i = 0; i < tet_vals.Size(); i++)
    {
      CHECK_THAT(T[i], WithinAbs(H[i], 1e-9));
    }
  }
#endif
}

TEST_CASE("LocalEdgeSplit", "[geodata][Serial]")
{
  auto mesh = TwoEdgeFanTetMesh();

  REQUIRE(mesh->GetNV() == 10);
  REQUIRE(mesh->GetNE() == 4);
  REQUIRE(mesh->GetNBE() == 12);

  const std::vector<std::pair<int, int>> split_edges{{5, 6}, {0, 1}};
  REQUIRE(mesh::LocalEdgeSplit(mesh, split_edges) == 2);

  CHECK(mesh->GetNV() == 12);
  CHECK(mesh->GetNE() == 8);
  CHECK(mesh->GetNBE() == 16);

  CheckVertex(*mesh, 10, {0.5, 0.0, 0.0});
  CheckVertex(*mesh, 11, {10.5, 0.0, 0.0});

  for (int e = 0; e < mesh->GetNE(); e++)
  {
    CAPTURE(e);
    CHECK(!ElementContainsEdge(*mesh->GetElement(e), 0, 1));
    CHECK(!ElementContainsEdge(*mesh->GetElement(e), 5, 6));
  }
  for (int be = 0; be < mesh->GetNBE(); be++)
  {
    CAPTURE(be);
    CHECK(!ElementContainsEdge(*mesh->GetBdrElement(be), 0, 1));
    CHECK(!ElementContainsEdge(*mesh->GetBdrElement(be), 5, 6));

    int f, o, e1, e2;
    mesh->GetBdrElementFace(be, &f, &o);
    mesh->GetFaceElements(f, &e1, &e2);
    CAPTURE(f, e1, e2);
    CHECK((e1 >= 0) != (e2 >= 0));
  }

  CHECK(mesh->CheckElementOrientation(false) == 0);
  CHECK(mesh->CheckBdrElementOrientation(false) == 0);
}

TEST_CASE("Boundary edge extraction", "[geodata][Parallel]")
{
  SECTION("2D endpoints")
  {
    auto serial_mesh =
        mfem::Mesh::MakeCartesian2D(4, 2, mfem::Element::TRIANGLE, false, 2.0, 1.0);
    mfem::ParMesh mesh(Mpi::World(), serial_mesh);
    auto marker = mesh::AttrToMarker(mesh.bdr_attributes.Max(), std::vector<int>{1});
    auto edges = mesh::GetBoundaryEdgeSegments(mesh, marker);

    REQUIRE(edges.size() == 2);
    for (const auto &edge : edges)
    {
      CHECK(edge.p0 == edge.p1);
      CHECK_THAT(edge.p0[1], WithinAbs(0.0, 1.0e-12));
      CHECK_THAT(edge.p0[2], WithinAbs(0.0, 1.0e-12));
    }
    CHECK_THAT(edges[0].p0[0] + edges[1].p0[0], WithinAbs(2.0, 1.0e-12));
  }

  SECTION("3D perimeter")
  {
    auto serial_mesh =
        mfem::Mesh::MakeCartesian3D(2, 3, 1, mfem::Element::TETRAHEDRON, 2.0, 3.0, 1.0);
    mfem::ParMesh mesh(Mpi::World(), serial_mesh);
    auto marker = mesh::AttrToMarker(mesh.bdr_attributes.Max(), std::vector<int>{1});
    auto edges = mesh::GetBoundaryEdgeSegments(mesh, marker);

    double length = 0.0;
    for (const auto &edge : edges)
    {
      double length_squared = 0.0;
      for (int d = 0; d < 3; d++)
      {
        const double delta = edge.p1[d] - edge.p0[d];
        length_squared += delta * delta;
      }
      length += std::sqrt(length_squared);
      CHECK_THAT(edge.p0[2], WithinAbs(0.0, 1.0e-12));
      CHECK_THAT(edge.p1[2], WithinAbs(0.0, 1.0e-12));
    }
    CHECK_THAT(length, WithinAbs(10.0, 1.0e-12));
  }
}

TEST_CASE("Boundary element edge extraction", "[geodata][Parallel]")
{
  SECTION("2D boundary segments")
  {
    auto serial_mesh =
        mfem::Mesh::MakeCartesian2D(4, 2, mfem::Element::TRIANGLE, false, 2.0, 1.0);
    mfem::ParMesh mesh(Mpi::World(), serial_mesh);
    auto marker = mesh::AttrToMarker(mesh.bdr_attributes.Max(), std::vector<int>{1});
    auto segments = mesh::GetBoundaryElementEdgeSegments(mesh, marker);

    REQUIRE(segments.size() == 4);
    double length = 0.0;
    for (const auto &segment : segments)
    {
      length += std::abs(segment.p1[0] - segment.p0[0]);
      CHECK_THAT(segment.p0[1], WithinAbs(0.0, 1.0e-12));
      CHECK_THAT(segment.p1[1], WithinAbs(0.0, 1.0e-12));
    }
    CHECK_THAT(length, WithinAbs(2.0, 1.0e-12));
  }

  SECTION("3D boundary face edges")
  {
    auto serial_mesh =
        mfem::Mesh::MakeCartesian3D(2, 3, 1, mfem::Element::TETRAHEDRON, 2.0, 3.0, 1.0);
    mfem::ParMesh mesh(Mpi::World(), serial_mesh);
    auto marker = mesh::AttrToMarker(mesh.bdr_attributes.Max(), std::vector<int>{1});
    auto segments = mesh::GetBoundaryElementEdgeSegments(mesh, marker);

    REQUIRE(!segments.empty());
    for (const auto &segment : segments)
    {
      CHECK_THAT(segment.p0[2], WithinAbs(0.0, 1.0e-12));
      CHECK_THAT(segment.p1[2], WithinAbs(0.0, 1.0e-12));
    }
  }
}

TEST_CASE("Edge distance tree and correction-tube marking", "[geodata][Serial]")
{
  mesh::BoundaryEdgeSegment segment{{2.0, 1.0, 0.0}, {2.0, 3.0, 0.0}};
  auto tree =
      std::make_shared<EdgeDistanceTree>(std::vector<mesh::BoundaryEdgeSegment>{segment});

  mfem::Vector point(2);
  point[0] = 3.0;
  point[1] = 2.0;
  CHECK_THAT(tree->DistanceSquared(point), WithinAbs(1.0, 1.0e-12));
  point[0] = 2.0;
  point[1] = 4.0;
  CHECK_THAT(tree->DistanceSquared(point), WithinAbs(1.0, 1.0e-12));

  auto serial_mesh =
      mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::QUADRILATERAL, false, 4.0, 4.0);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);

  const EdgeRefinementContext marking_context{tree, 0.2, 0.5, 0.25, 0.0};
  auto marked = MarkEdgeRefinementElements(mesh, {marking_context});
  CHECK(marked.Size() == 8);

  mfem::Vector indicators(mesh.GetNE());
  indicators = 1.0;
  const EdgeRefinementContext weighting_context{tree, 1.3, 0.5, 2.6, 0.25};
  CHECK(WeightEdgeCoreIndicators(mesh, {weighting_context}, indicators) == 4);
  int downweighted = 0;
  for (int i = 0; i < indicators.Size(); i++)
  {
    const double value = indicators[i];
    CHECK((value == 1.0 || value == 0.25));
    downweighted += (value == 0.25);
  }
  CHECK(downweighted == 4);
}

TEST_CASE("Polarized edge energy frame", "[geodata][Serial]")
{
  const mesh::BoundaryEdgeSegment segment{{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}};
  const auto top_frame = BuildEdgeFrame(segment, {0.0, 0.0, 1.0}, 3);
  const auto bottom_frame = BuildEdgeFrame(segment, {0.0, 0.0, -1.0}, 3);

  double point_data[] = {1.0, 0.5, 0.25};
  double field_data[] = {2.0, 3.0, 5.0};
  double displacement_data[] = {7.0, 11.0, 13.0};
  mfem::Vector point(point_data, 3);
  mfem::Vector field(field_data, 3);
  mfem::Vector displacement(displacement_data, 3);

  const auto top =
      GetPolarizedEdgeEnergyDensity(point, segment, top_frame, field, displacement);
  const auto reversed =
      GetPolarizedEdgeEnergyDensity(point, segment, bottom_frame, field, displacement);
  CHECK_THAT(top[0], WithinAbs(0.5 * 5.0 * 13.0, 1.0e-12));
  CHECK_THAT(top[1], WithinAbs(0.5 * 3.0 * 11.0, 1.0e-12));
  CHECK_THAT(top[2], WithinAbs(0.5 * 2.0 * 7.0, 1.0e-12));
  CHECK_THAT(top[3] + top[4] + top[5], WithinAbs(0.0, 1.0e-12));
  CHECK_THAT(top[0] + top[1] + top[2],
             WithinAbs(0.5 * (2.0 * 7.0 + 3.0 * 11.0 + 5.0 * 13.0), 1.0e-12));
  for (std::size_t component = 0; component < 3; component++)
  {
    CHECK_THAT(reversed[component], WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(reversed[3 + component], WithinAbs(top[component], 1.0e-12));
  }
}

TEST_CASE("Manual polarized edge extraction ignores process-normal seams",
          "[geodata][Serial]")
{
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  const std::array<double, 3> process_normal{0.0, 0.0, 1.0};

  int attribute = 0;
  std::vector<mesh::BoundaryEdgeSegment> unfiltered;
  for (const int candidate : serial_mesh.bdr_attributes)
  {
    auto marker = mesh::BdrAttrToMarker(mesh, std::vector<int>{candidate}, true);
    auto segments = mesh::GetBoundaryEdgeSegments(mesh, marker);
    const bool has_parallel =
        std::any_of(segments.begin(), segments.end(), [&](const auto &segment)
                    { return std::abs(segment.p1[2] - segment.p0[2]) > 0.5; });
    const bool has_transverse =
        std::any_of(segments.begin(), segments.end(),
                    [&](const auto &segment)
                    {
                      return std::hypot(segment.p1[0] - segment.p0[0],
                                        segment.p1[1] - segment.p0[1]) > 0.5;
                    });
    if (has_parallel && has_transverse)
    {
      attribute = candidate;
      unfiltered = std::move(segments);
      break;
    }
  }
  REQUIRE(attribute > 0);

  const auto tree = BuildEdgeDistanceTree(mesh, {attribute}, {}, process_normal);
  REQUIRE(tree->Size() < unfiltered.size());
  for (std::size_t i = 0; i < tree->Size(); i++)
  {
    const auto &segment = tree->GetSegment(i);
    CHECK(std::hypot(segment.p1[0] - segment.p0[0], segment.p1[1] - segment.p0[1]) > 0.0);
    CHECK_NOTHROW(BuildEdgeFrame(segment, process_normal, 3));
  }
}

TEST_CASE("Explicit edge segment exclusion (EdgeExcludeSegments)", "[geodata][Serial]")
{
  // Decision 282 rule B5: a spatial coupon's FOREIGN edges are removed from the
  // edge-distance perimeter by explicit segments. On the unit cube's bottom face (z = 0)
  // the perimeter is the four unit edges; excluding the edge x = 0 by its exact segment
  // drops exactly the mesh segments on it (any subdivision), an empty list reproduces the
  // tree, a tolerance miss keeps it, and excluding every edge fails closed.
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(2, 2, 1, mfem::Element::TETRAHEDRON);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  int bottom = 0;
  for (const int candidate : serial_mesh.bdr_attributes)
  {
    auto marker = mesh::BdrAttrToMarker(mesh, std::vector<int>{candidate}, true);
    auto segments = mesh::GetBoundaryEdgeSegments(mesh, marker);
    if (!segments.empty() &&
        std::all_of(segments.begin(), segments.end(), [](const auto &segment)
                    { return segment.p0[2] == 0.0 && segment.p1[2] == 0.0; }))
    {
      bottom = candidate;
      break;
    }
  }
  REQUIRE(bottom > 0);
  const auto full = BuildEdgeDistanceTree(mesh, {bottom}, {});
  const auto identity = BuildEdgeDistanceTree(mesh, {bottom}, {}, std::nullopt, {}, 0.0);
  REQUIRE(identity->Size() == full->Size());
  auto OnLeftEdge = [](const mesh::BoundaryEdgeSegment &segment)
  { return segment.p0[0] == 0.0 && segment.p1[0] == 0.0; };
  auto CountLeft = [&](const EdgeDistanceTree &tree)
  {
    std::size_t count = 0;
    for (std::size_t i = 0; i < tree.Size(); i++)
    {
      count += OnLeftEdge(tree.GetSegment(i)) ? 1 : 0;
    }
    return count;
  };
  const std::size_t left = CountLeft(*full);
  REQUIRE(left >= 2);  // the 2 x 2 subdivision splits the unit edge
  const std::vector<std::array<double, 6>> exclude = {{0.0, 0.0, 0.0, 0.0, 1.0, 0.0}};
  const auto excluded =
      BuildEdgeDistanceTree(mesh, {bottom}, {}, std::nullopt, exclude, 1.0e-3);
  CHECK(excluded->Size() == full->Size() - left);
  CHECK(CountLeft(*excluded) == 0);
  mfem::Vector probe(3);
  probe[0] = 0.05;
  probe[1] = 0.5;
  probe[2] = 0.0;
  CHECK_THAT(std::sqrt(full->DistanceSquared(probe)), WithinAbs(0.05, 1.0e-12));
  CHECK_THAT(std::sqrt(excluded->DistanceSquared(probe)), WithinAbs(0.5, 1.0e-12));
  // A listed segment 1e-2 away from the edge (beyond the tolerance) excludes nothing.
  const std::vector<std::array<double, 6>> miss = {{0.01, 0.0, 0.0, 0.01, 1.0, 0.0}};
  CHECK(BuildEdgeDistanceTree(mesh, {bottom}, {}, std::nullopt, miss, 1.0e-3)->Size() ==
        full->Size());
  // Half of the edge excludes only the mesh segments inside that half.
  const std::vector<std::array<double, 6>> half = {{0.0, 0.0, 0.0, 0.0, 0.5, 0.0}};
  CHECK(BuildEdgeDistanceTree(mesh, {bottom}, {}, std::nullopt, half, 1.0e-3)->Size() ==
        full->Size() - left / 2);
  const std::vector<std::array<double, 6>> all = {{0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                                                  {1.0, 0.0, 0.0, 1.0, 1.0, 0.0},
                                                  {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                                                  {0.0, 1.0, 0.0, 1.0, 1.0, 0.0}};
  CHECK_THROWS(BuildEdgeDistanceTree(mesh, {bottom}, {}, std::nullopt, all, 1.0e-3));
  CHECK_THROWS(BuildEdgeDistanceTree(mesh, {bottom}, {}, std::nullopt, exclude, 0.0));
}

TEST_CASE("Sagitta-aware edge segment exclusion (EdgeExcludeSegmentSagittas)",
          "[geodata][Serial]")
{
  // Decision 657 (a) (the 32b0 device-validity diagnosis, decision 656 (3)): a coupon mesh
  // builds a Chain: false context ARC as the true arc while EdgeExcludeSegments list the
  // device's CHORDS, so the plain rule B5 (ends + midpoint within 1e-3 R of the nearest
  // listed segment) leaves most of the arc's perimeter in the within-R tree. With the
  // chords' sagitta listed per segment (EdgeExcludeSegmentSagittas) the tolerance against
  // a chord is tolerance + sagitta and the whole arc is removed; a listed segment of
  // sagitta 0 keeps the plain rule bitwise (the straight-foreign models unchanged).
  // Geometry: a quarter circle of radius 1 followed by a straight piece along its tangent.
  constexpr double pi = 3.14159265358979323846;
  const double half_angle = pi / 8.0;  // two 45-degree chords
  const double sagitta = 1.0 - std::cos(half_angle);
  const double tolerance = 1.0e-3;
  auto ArcPoint = [](double angle)
  { return std::array<double, 3>{std::cos(angle), std::sin(angle), 0.0}; };
  // The perimeter: 40 segments on the true arc (0 -> 90 degrees) + 10 on the straight
  // piece x = 1, y in [-1, 0] (the tangent continuation at the arc's start) + 4 of another
  // conductor's edge far away (never excluded).
  std::vector<mesh::BoundaryEdgeSegment> perimeter;
  constexpr int arc_segments = 40;
  for (int i = 0; i < arc_segments; i++)
  {
    perimeter.push_back({ArcPoint(0.5 * pi * i / arc_segments),
                         ArcPoint(0.5 * pi * (i + 1) / arc_segments)});
  }
  for (int i = 0; i < 10; i++)
  {
    perimeter.push_back({{1.0, -0.1 * i, 0.0}, {1.0, -0.1 * (i + 1), 0.0}});
  }
  for (int i = 0; i < 4; i++)
  {
    perimeter.push_back({{3.0, 0.25 * i, 0.0}, {3.0, 0.25 * (i + 1), 0.0}});
  }
  const std::size_t total = perimeter.size();
  auto ArcLength = [&]()
  {
    double length = 0.0;
    for (int i = 0; i < arc_segments; i++)
    {
      length += 2.0 * std::sin(0.25 * pi / arc_segments);
    }
    return length;
  };
  // The exclusion list: the two chords of the arc (sagitta s each) + the straight piece.
  const auto a = ArcPoint(0.0), m = ArcPoint(0.25 * pi), b = ArcPoint(0.5 * pi);
  const std::vector<std::array<double, 6>> exclude = {{a[0], a[1], 0.0, m[0], m[1], 0.0},
                                                      {m[0], m[1], 0.0, b[0], b[1], 0.0},
                                                      {1.0, 0.0, 0.0, 1.0, -1.0, 0.0}};
  const double chords_length = 2.0 * 2.0 * std::sin(half_angle) + 1.0;

  SECTION("the plain rule keeps the arc (the 32b0 defect reproduced)")
  {
    auto segments = perimeter;
    EdgeExclusionRecord record;
    ExcludeCoincidentEdgeSegments(exclude, {}, tolerance, 3, segments, &record);
    // The straight piece is removed exactly; no arc segment has all three points within
    // 1e-3 of a chord (the first off-vertex node already deviates ~0.015).
    CHECK(segments.size() == total - 10);
    CHECK(record.exclude_segments_count == 3);
    CHECK_THAT(record.exclude_segments_length, WithinRel(chords_length, 1.0e-12));
    CHECK(record.perimeter_count == total);
    CHECK(record.excluded_perimeter_count == 10);
    CHECK_THAT(record.excluded_perimeter_length, WithinRel(1.0, 1.0e-12));
    CHECK(record.max_sagitta == 0.0);
    CHECK(record.tolerance == tolerance);
    CHECK(record.excluded_perimeter_length / record.exclude_segments_length < 0.5);
    // The three listed pieces form one chain (1, -1) -> (1, 0) -> arc -> (0, 1); (1, 0) is
    // a joint, not an end, so the first arc segment (adjacent to the removed piece) is NOT
    // credited; the end (0, 1) is contained by the kept last arc segment (the same length),
    // the end (1, -1) by nothing kept. The straight piece is one run.
    CHECK_THAT(record.straddling_perimeter_length,
               WithinRel(2.0 * std::sin(0.25 * pi / arc_segments), 1.0e-12));
    CHECK(record.removed_runs == 1);
    CHECK(record.listed_chains == 1);
    // Explicit zero sagittas are the plain rule.
    auto zero = perimeter;
    ExcludeCoincidentEdgeSegments(exclude, {0.0, 0.0, 0.0}, tolerance, 3, zero);
    REQUIRE(zero.size() == segments.size());
    for (std::size_t i = 0; i < zero.size(); i++)
    {
      CHECK(zero[i].p0 == segments[i].p0);
      CHECK(zero[i].p1 == segments[i].p1);
    }
  }

  SECTION("the sagitta-aware rule removes the whole arc and nothing else")
  {
    auto segments = perimeter;
    EdgeExclusionRecord record;
    ExcludeCoincidentEdgeSegments(exclude, {sagitta, sagitta, 0.0}, tolerance, 3, segments,
                                  &record);
    // Every arc segment (its ends on the circle within the sagitta of its chord, its
    // midpoint inside) and the straight piece are removed; the far edge stays. The first
    // arc segment is nearest to the STRAIGHT piece at its midpoint and far end (beyond
    // 1e-3 of it) and is found on the chord by the per-segment test.
    CHECK(segments.size() == 4);
    for (const auto &segment : segments)
    {
      CHECK(segment.p0[0] == 3.0);
      CHECK(segment.p1[0] == 3.0);
    }
    CHECK(record.excluded_perimeter_count == arc_segments + 10);
    CHECK_THAT(record.excluded_perimeter_length, WithinRel(ArcLength() + 1.0, 1.0e-12));
    CHECK_THAT(record.max_sagitta, WithinRel(sagitta, 1.0e-12));
    // Coverage = (arc length + 1) / (chords + 1) >= 1: the arc is longer than its chords.
    const double coverage =
        record.excluded_perimeter_length / record.exclude_segments_length;
    CHECK(coverage >= 1.0);
    CHECK(coverage < 1.03);
    // Nothing kept contains the chain's ends ((1, -1) and (0, 1)); one run of arc + piece.
    CHECK(record.straddling_perimeter_length == 0.0);
    CHECK(record.removed_runs == 1);
    CHECK(record.listed_chains == 1);
  }

  SECTION("the end allowance: a listed end inside a perimeter segment keeps it whole")
  {
    // Decision 665: a straight piece listed from y = 0.02 to y = 0.98 along a perimeter of
    // 20 segments of length 0.05 on x = 1 (y in [0, 1]). The two end segments have one
    // point 0.02 beyond the listed end (> 1e-3) and are kept whole; the 18 inside are
    // removed: the raw coverage reads 0.9 / 0.96 = 0.9375 (< 0.99 with a CORRECT
    // exclusion, the census's absolute deficit per listed end), the kept segments sharing
    // a vertex with a removed one are exactly the two straddlers (0.1), and the end
    // allowance coverage reads 1.0 / 0.96 >= 0.99.
    std::vector<mesh::BoundaryEdgeSegment> segments;
    for (int i = 0; i < 20; i++)
    {
      segments.push_back({{1.0, 0.05 * i, 0.0}, {1.0, 0.05 * (i + 1), 0.0}});
    }
    segments.push_back({{3.0, 0.0, 0.0}, {3.0, 1.0, 0.0}});
    const std::vector<std::array<double, 6>> piece = {{1.0, 0.02, 0.0, 1.0, 0.98, 0.0}};
    EdgeExclusionRecord record;
    ExcludeCoincidentEdgeSegments(piece, {0.0}, tolerance, 3, segments, &record);
    CHECK(segments.size() == 3);
    CHECK(record.excluded_perimeter_count == 18);
    CHECK_THAT(record.excluded_perimeter_length, WithinRel(0.9, 1.0e-12));
    CHECK_THAT(record.straddling_perimeter_length, WithinRel(0.1, 1.0e-12));
    CHECK(record.removed_runs == 1);
    CHECK(record.listed_chains == 1);
    CHECK(record.excluded_perimeter_length / record.exclude_segments_length < 0.99);
    CHECK((record.excluded_perimeter_length + record.straddling_perimeter_length) /
              record.exclude_segments_length >=
          0.99);
    // The same piece listed through the automatic-edge path records the same lengths.
    MetalEdgeGeometry geometry;
    for (int i = 0; i <= 20; i++)
    {
      MetalEdgeVertex vertex;
      vertex.coordinate = {1.0, 0.05 * i, 0.0};
      geometry.vertices.push_back(vertex);
    }
    std::vector<std::size_t> indices;
    for (int i = 0; i < 20; i++)
    {
      MetalEdgeSegment segment;
      segment.vertices = {static_cast<std::size_t>(i), static_cast<std::size_t>(i + 1)};
      geometry.segments.push_back(segment);
      indices.push_back(i);
    }
    EdgeExclusionRecord automatic;
    ExcludeCoincidentMetalEdgeSegmentIndices(geometry, piece, {0.0}, tolerance, 3, indices,
                                             &automatic);
    CHECK(indices == std::vector<std::size_t>{0, 19});
    CHECK_THAT(automatic.excluded_perimeter_length, WithinRel(0.9, 1.0e-12));
    CHECK_THAT(automatic.straddling_perimeter_length, WithinRel(0.1, 1.0e-12));
    CHECK(automatic.removed_runs == 1);
    CHECK(automatic.listed_chains == 1);
  }

  SECTION("the end allowance credits one kept segment per listed end only")
  {
    // Decision 667 (1): a kept segment is credited only when it CONTAINS a listed chain's
    // end point (within the tolerance), at most one per listed end, so a fragmented
    // exclusion (holes of kept perimeter inside the listed piece) gains nothing from its
    // runs' neighbours. The confirmation's 7 patterns on a perimeter of 20 segments of 0.05
    // on x = 1 (y in [0, 1]) with the piece listed over all of them: the same masks and
    // values as the census-v3 tool's test (test_exclusion_coverage_census_v3.py).
    std::vector<mesh::BoundaryEdgeSegment> line;
    for (int i = 0; i < 20; i++)
    {
      line.push_back({{1.0, 0.05 * i, 0.0}, {1.0, 0.05 * (i + 1), 0.0}});
    }
    const std::vector<std::array<double, 6>> whole = {{1.0, 0.0, 0.0, 1.0, 1.0, 0.0}};
    const std::vector<std::array<double, 6>> halves = {{1.0, 0.0, 0.0, 1.0, 0.5, 0.0},
                                                       {1.0, 0.5, 0.0, 1.0, 1.0, 0.0}};
    auto Read = [&](const std::vector<int> &kept_indices,
                    const std::vector<std::array<double, 6>> &listed)
    {
      std::vector<bool> excluded(line.size(), true);
      for (const int i : kept_indices)
      {
        excluded[i] = false;
      }
      EdgeExclusionRecord record;
      for (const auto &s : listed)
      {
        record.exclude_segments_length += std::abs(s[4] - s[1]);
      }
      for (std::size_t i = 0; i < line.size(); i++)
      {
        record.excluded_perimeter_length += excluded[i] ? 0.05 : 0.0;
      }
      RecordEdgeExclusionEnds(line, excluded, listed, tolerance, record);
      return record;
    };
    auto Allowance = [](const EdgeExclusionRecord &r)
    {
      return (r.excluded_perimeter_length + r.straddling_perimeter_length) /
             r.exclude_segments_length;
    };
    auto Passes = [&](const EdgeExclusionRecord &r)
    {
      const double raw = r.excluded_perimeter_length / r.exclude_segments_length;
      return raw <= 1.01 && std::max(raw, Allowance(r)) >= 0.99;
    };
    // correct: all 20 removed.
    auto correct = Read({}, whole);
    CHECK(correct.straddling_perimeter_length == 0.0);
    CHECK(correct.removed_runs == 1);
    CHECK(correct.listed_chains == 1);
    CHECK_THAT(Allowance(correct), WithinRel(1.0, 1.0e-12));
    CHECK(Passes(correct));
    // alternating: 0, 2, ..., 18 removed; the kept segment 19 contains the end y = 1.
    auto alternating = Read({1, 3, 5, 7, 9, 11, 13, 15, 17, 19}, whole);
    CHECK_THAT(alternating.straddling_perimeter_length, WithinRel(0.05, 1.0e-12));
    CHECK(alternating.removed_runs == 10);
    CHECK(alternating.listed_chains == 1);
    CHECK_THAT(Allowance(alternating), WithinRel(0.55, 1.0e-12));
    CHECK_FALSE(Passes(alternating));
    // 1-segment holes at 3, 7, 11, 15.
    auto holes1 = Read({3, 7, 11, 15}, whole);
    CHECK(holes1.straddling_perimeter_length == 0.0);
    CHECK(holes1.removed_runs == 5);
    CHECK_THAT(Allowance(holes1), WithinRel(0.8, 1.0e-12));
    CHECK_FALSE(Passes(holes1));
    // three 2-segment holes.
    auto holes2 = Read({3, 4, 9, 10, 15, 16}, whole);
    CHECK(holes2.straddling_perimeter_length == 0.0);
    CHECK(holes2.removed_runs == 4);
    CHECK_THAT(Allowance(holes2), WithinRel(0.7, 1.0e-12));
    CHECK_FALSE(Passes(holes2));
    // two 3-segment holes.
    auto holes3 = Read({4, 5, 6, 13, 14, 15}, whole);
    CHECK(holes3.straddling_perimeter_length == 0.0);
    CHECK(holes3.removed_runs == 3);
    CHECK_THAT(Allowance(holes3), WithinRel(0.7, 1.0e-12));
    CHECK_FALSE(Passes(holes3));
    // one 4-segment hole.
    auto hole4 = Read({8, 9, 10, 11}, whole);
    CHECK(hole4.straddling_perimeter_length == 0.0);
    CHECK(hole4.removed_runs == 2);
    CHECK_THAT(Allowance(hole4), WithinRel(0.8, 1.0e-12));
    CHECK_FALSE(Passes(hole4));
    // 32b0-like: 2 segments at each chord vertex of two listed halves (one chain).
    std::vector<int> kept_32b0;
    for (int i = 2; i < 18; i++)
    {
      if (i != 9 && i != 10)
      {
        kept_32b0.push_back(i);
      }
    }
    auto like_32b0 = Read(kept_32b0, halves);
    CHECK(like_32b0.straddling_perimeter_length == 0.0);
    CHECK(like_32b0.removed_runs == 3);
    CHECK(like_32b0.listed_chains == 1);
    CHECK_THAT(Allowance(like_32b0), WithinRel(0.3, 1.0e-12));
    CHECK_FALSE(Passes(like_32b0));
    // One kept segment per listed end: the piece (1, 0.5) -> (1, 1) with segments 9 AND 10
    // kept (both contain the end y = 0.5) credits one segment.
    const std::vector<std::array<double, 6>> upper = {{1.0, 0.5, 0.0, 1.0, 1.0, 0.0}};
    std::vector<int> kept_two_at_end = {9, 10};
    for (int i = 0; i < 9; i++)
    {
      kept_two_at_end.push_back(i);
    }
    auto two_at_end = Read(kept_two_at_end, upper);
    CHECK_THAT(two_at_end.straddling_perimeter_length, WithinRel(0.05, 1.0e-12));
    CHECK(two_at_end.removed_runs == 1);
    // A kept segment credited once: the piece (1, 0.51) -> (1, 0.54) lies inside segment
    // 10, nothing is removed, both ends are in that one segment.
    const std::vector<std::array<double, 6>> inside = {{1.0, 0.51, 0.0, 1.0, 0.54, 0.0}};
    std::vector<int> all_kept(20);
    std::iota(all_kept.begin(), all_kept.end(), 0);
    auto once = Read(all_kept, inside);
    CHECK_THAT(once.straddling_perimeter_length, WithinRel(0.05, 1.0e-12));
    CHECK(once.removed_runs == 0);
    CHECK(once.listed_chains == 1);
    // The chain count joins listed ends within the tolerance: the two halves offset by
    // 5e-4 at their joint are still one chain with 2 ends; offset by 2e-3, two chains.
    const std::vector<std::array<double, 6>> nearly = {{1.0, 0.0, 0.0, 1.0, 0.5, 0.0},
                                                       {1.0, 0.5005, 0.0, 1.0, 1.0, 0.0}};
    CHECK(Read({}, nearly).listed_chains == 1);
    const std::vector<std::array<double, 6>> apart = {{1.0, 0.0, 0.0, 1.0, 0.5, 0.0},
                                                      {1.0, 0.502, 0.0, 1.0, 1.0, 0.0}};
    CHECK(Read({}, apart).listed_chains == 2);
  }

  SECTION("a sagitta widens the band of ITS segment only")
  {
    // A straight perimeter 2e-3 off the straight listed piece is kept whatever the arc
    // chords' sagitta: the straight segment's own band stays 1e-3.
    std::vector<mesh::BoundaryEdgeSegment> segments = {
        {{1.002, -0.2, 0.0}, {1.002, -0.6, 0.0}},
        {{1.002, -0.6, 0.0}, {1.002, -1.0, 0.0}},
        {{3.0, 0.0, 0.0}, {3.0, 1.0, 0.0}}};
    auto kept = segments;
    ExcludeCoincidentEdgeSegments(exclude, {sagitta, sagitta, 0.0}, tolerance, 3, kept);
    CHECK(kept.size() == 3);
    // Listed with a sagitta of its own, the straight piece's band widens and the offset
    // perimeter is removed.
    auto widened = segments;
    ExcludeCoincidentEdgeSegments(exclude, {sagitta, sagitta, 0.005}, tolerance, 3,
                                  widened);
    CHECK(widened.size() == 1);
  }

  SECTION("a sagitta list of the wrong length fails closed")
  {
    auto segments = perimeter;
    CHECK_THROWS(ExcludeCoincidentEdgeSegments(exclude, {sagitta}, tolerance, 3, segments));
  }
}

TEST_CASE("Boundary edge exclusion", "[geodata][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 2, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  const int dim = par_mesh->Dimension();
  Mesh palace_mesh(std::move(par_mesh));

  mfem::H1_FECollection h1_fec(1, dim);
  mfem::ND_FECollection nd_fec(1, dim);
  FiniteElementSpace h1_fespace(palace_mesh, &h1_fec);
  FiniteElementSpace nd_fespace(palace_mesh, &nd_fec);

  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC, palace_mesh);

  config::InterfaceDielectricData dielectric;
  dielectric.attributes = {1};
  dielectric.edge_attributes = {1};
  dielectric.edge_exclude_attributes = {1};
  dielectric.edge_distances = {0.1};
  config::BoundaryPostData postpro;
  postpro.dielectric.try_emplace(1, std::move(dielectric));

  CHECK_THROWS(SurfacePostOperator(postpro, ProblemType::ELECTROSTATIC, mat_op, h1_fespace,
                                   nd_fespace));
}

TEST_CASE("Flux recovery rejects cracked interfaces", "[geodata][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(1, 1, mfem::Element::TRIANGLE));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  const int dim = par_mesh->Dimension();
  Mesh palace_mesh(std::move(par_mesh));

  mfem::H1_FECollection h1_fec(1, dim);
  mfem::ND_FECollection nd_fec(1, dim);
  FiniteElementSpace h1_fespace(palace_mesh, &h1_fec);
  FiniteElementSpace nd_fespace(palace_mesh, &nd_fec);

  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC, palace_mesh);

  config::InterfaceDielectricData dielectric;
  dielectric.attributes = {1};
  dielectric.type = InterfaceDielectric::MA;
  dielectric.flux_recovery = true;
  config::BoundaryPostData postpro;
  postpro.dielectric.try_emplace(1, std::move(dielectric));

  // The legacy key: the interface attribute is in the mesh crack set.
  config::BoundaryData boundaries;
  boundaries.cracked_attributes = {1};
  CHECK_THROWS_WITH(
      SurfacePostOperator(postpro, ProblemType::ELECTROSTATIC, mat_op, h1_fespace,
                          nd_fespace, {}, &boundaries),
      Catch::Matchers::ContainsSubstring("\"FluxRecovery\" is not supported"));
}

TEST_CASE("Boundary edge extraction on cracked 3D sheet", "[geodata][Serial][Parallel]")
{
  auto serial_mesh = CrackedSquareSheetTetMesh();
  mfem::ParMesh mesh(Mpi::World(), *serial_mesh);
  auto marker = mesh::AttrToMarker(mesh.bdr_attributes.Max(), std::vector<int>{5});
  auto edges = mesh::GetBoundaryEdgeSegments(mesh, marker);

  double length = 0.0;
  for (const auto &edge : edges)
  {
    CAPTURE(edge.p0[0], edge.p0[1], edge.p1[0], edge.p1[1]);
    double length_squared = 0.0;
    for (int d = 0; d < 3; d++)
    {
      const double delta = edge.p1[d] - edge.p0[d];
      length_squared += delta * delta;
    }
    length += std::sqrt(length_squared);
    CHECK_THAT(edge.p0[2], WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(edge.p1[2], WithinAbs(0.0, 1.0e-12));
    const bool on_boundary =
        (std::abs(edge.p0[0]) < 1.0e-12 && std::abs(edge.p1[0]) < 1.0e-12) ||
        (std::abs(edge.p0[0] - 1.0) < 1.0e-12 && std::abs(edge.p1[0] - 1.0) < 1.0e-12) ||
        (std::abs(edge.p0[1]) < 1.0e-12 && std::abs(edge.p1[1]) < 1.0e-12) ||
        (std::abs(edge.p0[1] - 1.0) < 1.0e-12 && std::abs(edge.p1[1] - 1.0) < 1.0e-12);
    CHECK(on_boundary);
  }
  REQUIRE(edges.size() == 4);
  CHECK_THAT(length, WithinAbs(4.0, 1.0e-12));
}

TEST_CASE("Automatic metal edge extraction and classification",
          "[geodata][metaledge][Serial][Parallel]")
{
  auto serial_mesh = CrackedSquareSheetTetMesh(true);
  mfem::ParMesh mesh(Mpi::World(), *serial_mesh);

  auto MakeDielectric = [](InterfaceDielectric type)
  {
    config::InterfaceDielectricData dielectric;
    dielectric.type = type;
    dielectric.attributes = {5};
    return dielectric;
  };

  config::BoundaryData boundaries;
  boundaries.pec.attributes = {5};
  boundaries.postpro.dielectric.emplace(11, MakeDielectric(InterfaceDielectric::SA));
  boundaries.postpro.dielectric.emplace(12, MakeDielectric(InterfaceDielectric::MS));
  boundaries.postpro.dielectric.emplace(13, MakeDielectric(InterfaceDielectric::MA));

  // Unit square sheet: the joint rule at R = 0.1 (the geometry has no joints).
  MetalSurfaceExtraction surface = JointNoiseExtraction(0.1);
  surface.classify_components = true;
  surface.retain_faces = true;
  auto geometry = ExtractMetalEdgeGeometry(mesh, boundaries, surface);
  REQUIRE(geometry.components == 1);
  REQUIRE(geometry.physical_components == 1);
  REQUIRE(geometry.metal_components == 1);
  REQUIRE(!geometry.surface_faces.empty());
  CHECK(std::all_of(geometry.surface_faces.begin(), geometry.surface_faces.end(),
                    [](const auto &face)
                    { return face.component == 0 && face.vertices.size() >= 3; }));
  REQUIRE(geometry.vertices.size() == 4);
  REQUIRE(geometry.segments.size() == 4);
  int truncated_vertices = 0;
  int physical_endpoints = 0;
  for (const auto &vertex : geometry.vertices)
  {
    CHECK(vertex.type == MetalEdgeVertexType::CORNER);
    CHECK(vertex.segments.size() == 2);
    truncated_vertices += vertex.on_truncation_boundary;
    physical_endpoints +=
        vertex.physical_type == std::optional{MetalEdgeVertexType::ENDPOINT};
  }
  CHECK(truncated_vertices == 2);
  CHECK(physical_endpoints == 2);
  int truncation_segments = 0;
  for (const auto &segment : geometry.segments)
  {
    CHECK(segment.component == 0);
    CHECK(segment.metal_component == 0);
    CHECK(segment.metal_attributes == std::vector<int>{5});
    REQUIRE(segment.conditions.size() == 1);
    CHECK(segment.conditions[0].type == MetalBoundaryConditionType::PEC);
    CHECK(segment.conditions[0].index == 0);
    CHECK(segment.sa_interfaces == std::vector<int>{11});
    CHECK(segment.ms_interfaces == std::vector<int>{12});
    CHECK(segment.ma_interfaces == std::vector<int>{13});
    if (segment.type == MetalEdgeSegmentType::TRUNCATION)
    {
      truncation_segments++;
      CHECK(segment.physical_component == -1);
      CHECK(segment.truncation_attributes == std::vector<int>{6});
    }
    else
    {
      CHECK(segment.physical_component == 0);
      CHECK(segment.truncation_attributes.empty());
    }
  }
  CHECK(truncation_segments == 1);
  CHECK(geometry.physical_chains == 3);
  std::set<int> physical_chains;
  for (const auto &segment : geometry.segments)
  {
    if (segment.type == MetalEdgeSegmentType::PHYSICAL)
    {
      CHECK(segment.physical_chain >= 0);
      physical_chains.insert(segment.physical_chain);
    }
    else
    {
      CHECK(segment.physical_chain == -1);
    }
  }
  CHECK(physical_chains.size() == 3);

  const auto physical_segments =
      GetInterfaceMetalEdgeSegmentIndices(geometry, 11, InterfaceDielectric::SA);
  const auto process_normals = BuildMetalEdgeProcessNormals(
      mesh, geometry, physical_segments,
      [](int material_attribute) { return material_attribute == 2 ? 1.0 : 0.0; });
  REQUIRE(process_normals.size() == 3);
  for (const auto &normal : process_normals)
  {
    CHECK_THAT(normal[0], WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(normal[1], WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(normal[2], WithinAbs(1.0, 1.0e-12));
  }
  auto process_tree = BuildEdgeDistanceTree(geometry, physical_segments, process_normals);
  REQUIRE(process_tree->HasProcessNormals());
  REQUIRE(process_tree->HasMetadata());
  CHECK(process_tree->GetProcessNormal(0) == process_normals[0]);
  int metadata_endpoints = 0;
  for (std::size_t i = 0; i < process_tree->Size(); i++)
  {
    const auto &metadata = process_tree->GetMetadata(i);
    CHECK(metadata.component == 0);
    CHECK(metadata.chain >= 0);
    for (const int type : metadata.vertex_types)
    {
      CHECK((type == static_cast<int>(MetalEdgeVertexType::CORNER) ||
             type == static_cast<int>(MetalEdgeVertexType::ENDPOINT)));
      metadata_endpoints += type == static_cast<int>(MetalEdgeVertexType::ENDPOINT);
    }
  }
  CHECK(metadata_endpoints == 2);

  const auto fallback_normals = BuildMetalEdgeProcessNormals(
      mesh, geometry, physical_segments, [](int) { return 1.0; },
      std::array<double, 3>{0.0, 0.0, -1.0});
  for (const auto &normal : fallback_normals)
  {
    CHECK_THAT(normal[2], WithinAbs(-1.0, 1.0e-12));
  }

  surface.retain_faces = false;
  geometry = ExtractMetalEdgeGeometry(mesh, boundaries, surface);
  CHECK(geometry.metal_components == 1);
  CHECK(geometry.surface_faces.empty());
  CHECK(std::all_of(geometry.segments.begin(), geometry.segments.end(),
                    [](const auto &segment) { return segment.metal_component == 0; }));

  boundaries.pec.attributes.clear();
  config::ImpedanceData impedance;
  impedance.attributes = {5};
  boundaries.impedance.push_back(std::move(impedance));
  geometry = ExtractMetalEdgeGeometry(mesh, boundaries, JointNoiseExtraction(0.1));
  REQUIRE(geometry.segments.size() == 4);
  for (const auto &segment : geometry.segments)
  {
    REQUIRE(segment.conditions.size() == 1);
    CHECK(segment.conditions[0].type == MetalBoundaryConditionType::IMPEDANCE);
    CHECK(segment.conditions[0].index == 0);
  }
  const auto impedance_segments =
      GetInterfaceMetalEdgeSegmentIndices(geometry, 11, InterfaceDielectric::SA);
  const auto impedance_normals = BuildMetalEdgeProcessNormals(
      mesh, geometry, impedance_segments,
      [](int material_attribute) { return material_attribute == 2 ? 1.0 : 0.0; });
  CHECK(impedance_normals == process_normals);
}

TEST_CASE("Automatic metal edge extraction samples high-order rounded edges",
          "[geodata][metaledge][Serial][Parallel]")
{
  config::BoundaryData boundaries;
  boundaries.pec.attributes = {9};
  config::InterfaceDielectricData dielectric;
  dielectric.type = InterfaceDielectric::SA;
  dielectric.attributes = {9};
  boundaries.postpro.dielectric.emplace(1, std::move(dielectric));

  auto Extract = [&](int elements, bool high_order, bool rounded, double radius = 0.125)
  {
    auto mesh = RoundedIslandSheetHexMesh(elements, high_order, rounded);
    // The geometric joint rule (kJointNoiseSagittaOverRadius, USER decision 121 (B)) at
    // R = the fillet radius (0.125) unless given: a fillet joint of tens of degrees on
    // chords of ~0.2 rho implies a sagitta (c / 2) tan(t / 4) far below 0.05 R and is a
    // REGULAR joint of the extraction's chain (the identification's arc rule reads the
    // non-collinear joints regardless of their type).
    MetalSurfaceExtraction surface = JointNoiseExtraction(radius);
    surface.retain_faces = true;
    auto geometry = ExtractMetalEdgeGeometry(*mesh, boundaries, surface);
    return std::pair{std::move(mesh), std::move(geometry)};
  };
  auto [coarse_mesh, coarse] = Extract(8, true, true);
  mesh::MeshEdgeSegmentCache edge_cache(*coarse_mesh);
  const auto &cached_edge = edge_cache.Get(0);
  CHECK(&cached_edge == &edge_cache.Get(0));
  const auto uncached_edge = mesh::GetMeshEdgeSegments(*coarse_mesh, 0);
  REQUIRE(cached_edge.size() == uncached_edge.size());
  for (std::size_t i = 0; i < cached_edge.size(); i++)
  {
    CHECK(cached_edge[i].p0 == uncached_edge[i].p0);
    CHECK(cached_edge[i].p1 == uncached_edge[i].p1);
  }
  auto [fine_mesh, fine] = Extract(16, false, true);
  auto [sharp_mesh, sharp] = Extract(8, false, false);

  auto CountPhysicalVertexTypes =
      [](const MetalEdgeGeometry &geometry, MetalEdgeVertexType type)
  {
    return std::count_if(geometry.vertices.begin(), geometry.vertices.end(),
                         [&](const auto &vertex)
                         { return vertex.physical_type == std::optional{type}; });
  };
  auto PhysicalLength = [](const MetalEdgeGeometry &geometry)
  {
    double length = 0.0;
    for (const auto &segment : geometry.segments)
    {
      if (segment.type != MetalEdgeSegmentType::PHYSICAL)
      {
        continue;
      }
      const auto &p0 = geometry.vertices[segment.vertices[0]].coordinate;
      const auto &p1 = geometry.vertices[segment.vertices[1]].coordinate;
      double length_squared = 0.0;
      for (int d = 0; d < 3; d++)
      {
        length_squared += (p1[d] - p0[d]) * (p1[d] - p0[d]);
      }
      length += std::sqrt(length_squared);
    }
    return length;
  };

  CHECK(coarse.physical_components == 1);
  CHECK(fine.physical_components == 1);
  // The perimeter's vertex classification uses the geometric joint noise rule
  // (kJointNoiseSagittaOverRadius = 0.05, USER decision 121 (B)): at R = rho every joint of
  // the sampled rounded corners (2-3 deg per sampled sub-segment of 0.006: implied sagitta
  // ~3e-5 = 2e-4 rho) is a REGULAR joint, the rounded loop is ONE chain of no corners; at
  // R = rho / 250 (threshold 2.5e-5) the same joints imply sagittas above 0.05 R and are
  // CORNER vertices and chain breaks (one chain per corner on a closed loop) — the former
  // 1 deg angular reading of decision 117(4), which made every sampled fillet joint a
  // corner whatever its size.
  const auto coarse_corners = CountPhysicalVertexTypes(coarse, MetalEdgeVertexType::CORNER);
  const auto fine_corners = CountPhysicalVertexTypes(fine, MetalEdgeVertexType::CORNER);
  CHECK(coarse_corners == 0);
  CHECK(fine_corners == 0);
  CHECK(coarse.physical_chains == 1);
  CHECK(fine.physical_chains == 1);
  {
    auto [small_mesh, small_R] = Extract(8, true, true, 0.0005);
    const auto small_corners =
        CountPhysicalVertexTypes(small_R, MetalEdgeVertexType::CORNER);
    CHECK(small_corners >= 8);
    CHECK(small_R.physical_chains == small_corners);
    CHECK(small_R.segments.size() == coarse.segments.size());
  }
  CHECK(coarse.segments.size() > fine.segments.size());

  const double exact_perimeter =
      4.0 * (2.0 * 0.25 - 2.0 * 0.125) + 2.0 * std::acos(-1.0) * 0.125;
  CHECK_THAT(PhysicalLength(coarse), WithinAbs(exact_perimeter, 2.0e-3));
  CHECK_THAT(PhysicalLength(fine), WithinAbs(exact_perimeter, 6.0e-3));
  CHECK_THAT(PhysicalLength(coarse), WithinAbs(PhysicalLength(fine), 6.0e-3));

  // The retained faces are rank-local (metaledge.cpp, "(7) Retained faces"): a rank of a
  // parallel run may hold no metal face, so the preconditions and checks on them are
  // collective (a rank-local REQUIRE ends the case on some ranks only and leaves the others
  // in the next collective).
  const auto comm = coarse_mesh->GetComm();
  std::array<long long, 2> surface_face_counts = {
      static_cast<long long>(coarse.surface_faces.size()),
      static_cast<long long>(fine.surface_faces.size())};
  Mpi::GlobalSum(2, surface_face_counts.data(), comm);
  REQUIRE(surface_face_counts[0] > 0);
  REQUIRE(surface_face_counts[1] > 0);
  bool coarse_faces_in_component =
      std::all_of(coarse.surface_faces.begin(), coarse.surface_faces.end(),
                  [](const auto &face) { return face.component == 0; });
  bool fine_faces_are_quadrilaterals =
      std::all_of(fine.surface_faces.begin(), fine.surface_faces.end(), [](const auto &face)
                  { return face.component == 0 && face.vertices.size() == 4; });
  bool coarse_faces_high_order =
      std::any_of(coarse.surface_faces.begin(), coarse.surface_faces.end(),
                  [](const auto &face) { return face.vertices.size() > 8; });
  Mpi::GlobalAnd(1, &coarse_faces_in_component, comm);
  Mpi::GlobalAnd(1, &fine_faces_are_quadrilaterals, comm);
  Mpi::GlobalOr(1, &coarse_faces_high_order, comm);
  CHECK(coarse_faces_in_component);
  CHECK(fine_faces_are_quadrilaterals);
  CHECK(coarse_faces_high_order);

  CHECK(sharp.physical_components == 1);
  CHECK(sharp.physical_chains == 4);
  CHECK(CountPhysicalVertexTypes(sharp, MetalEdgeVertexType::CORNER) == 4);

  const auto coarse_segments =
      GetInterfaceMetalEdgeSegmentIndices(coarse, 1, InterfaceDielectric::SA);
  REQUIRE(coarse_segments.size() == coarse.segments.size());
  const auto process_normals = BuildMetalEdgeProcessNormals(
      *coarse_mesh, coarse, coarse_segments, [](int) { return 1.0; },
      std::array<double, 3>{0.0, 1.0, 0.0});
  REQUIRE(process_normals.size() == coarse_segments.size());
  for (const auto &normal : process_normals)
  {
    CHECK_THAT(normal[0], WithinAbs(0.0, 1.0e-10));
    CHECK_THAT(normal[1], WithinAbs(1.0, 1.0e-10));
    CHECK_THAT(normal[2], WithinAbs(0.0, 1.0e-10));
  }
  const auto gap_directions =
      BuildMetalEdgeGapDirections(*coarse_mesh, coarse, coarse_segments, process_normals);
  REQUIRE(gap_directions.size() == coarse_segments.size());
  for (const auto &gap : gap_directions)
  {
    CHECK_THAT(gap[1], WithinAbs(0.0, 1.0e-10));
    CHECK_THAT(std::hypot(gap[0], gap[2]), WithinAbs(1.0, 1.0e-10));
  }
}

TEST_CASE("Automatic edge distance to nonregular vertex crosses mesh segments",
          "[geodata][metaledge][Serial][Parallel]")
{
  MetalEdgeGeometry geometry;
  geometry.vertices.resize(4);
  geometry.segments.resize(3);
  for (std::size_t vertex = 0; vertex < geometry.vertices.size(); vertex++)
  {
    geometry.vertices[vertex].coordinate = {static_cast<double>(vertex), 0.0, 0.0};
    geometry.vertices[vertex].physical_type =
        vertex == 0 || vertex + 1 == geometry.vertices.size()
            ? MetalEdgeVertexType::ENDPOINT
            : MetalEdgeVertexType::REGULAR;
  }
  for (std::size_t segment = 0; segment < geometry.segments.size(); segment++)
  {
    geometry.segments[segment].vertices = {segment, segment + 1};
    geometry.segments[segment].physical_component = 0;
    geometry.segments[segment].physical_chain = 0;
    geometry.vertices[segment].segments.push_back(segment);
    geometry.vertices[segment + 1].segments.push_back(segment);
  }

  const auto tree = BuildEdgeDistanceTree(geometry, {0, 1, 2});
  CHECK(tree->GetMetadata(0).vertex_distances == std::array<double, 2>{0.0, 1.0});
  CHECK(tree->GetMetadata(1).vertex_distances == std::array<double, 2>{1.0, 1.0});
  CHECK(tree->GetMetadata(2).vertex_distances == std::array<double, 2>{1.0, 0.0});

  mfem::Vector point(3);
  point = 0.0;
  point[0] = 1.5;
  const auto middle = tree->Nearest(point);
  REQUIRE(middle.segment == 1);
  CHECK_THAT(tree->DistanceAlongEdgeToNonregularVertex(point, middle.segment),
             WithinAbs(1.5, 1.0e-12));

  point[0] = 0.25;
  const auto endpoint = tree->Nearest(point);
  REQUIRE(endpoint.segment == 0);
  CHECK_THAT(tree->DistanceAlongEdgeToNonregularVertex(point, endpoint.segment),
             WithinAbs(0.25, 1.0e-12));

  const EdgeDistanceTree manual_tree(
      std::vector<mesh::BoundaryEdgeSegment>{{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}});
  CHECK(std::isinf(manual_tree.DistanceAlongEdgeToNonregularVertex(point, 0)));

  for (auto &vertex : geometry.vertices)
  {
    vertex.physical_type = MetalEdgeVertexType::REGULAR;
  }
  const auto regular_loop_tree = BuildEdgeDistanceTree(geometry, {0, 1, 2});
  CHECK(std::isinf(
      regular_loop_tree->DistanceAlongEdgeToNonregularVertex(point, endpoint.segment)));
}

TEST_CASE("Automatic metal edge extraction on 3D CPW",
          "[geodata][metaledge][Serial][Parallel]")
{
  // P3-23: both crack forms (the legacy key of SD-801; upstream's default = the key unset).
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/cpw3d_surface/cpw3d_surface_validation_thin.json";
  IoData iodata(config_path.c_str(), false);
  iodata.model.mesh = fs::path(PALACE_TEST_DATA_DIR) / "mesh/cpw3d-surface-nc.msh";
  test::ApplyCrackPolicy(iodata, crack_policy);  // explicit crack policy (P3-23)
  iodata.model.refinement.max_it = 0;
  auto mesh = mesh::ReadMesh(iodata, Mpi::World());

  for (const auto &[index, dielectric] : iodata.boundaries.postpro.dielectric)
  {
    (void)index;
    CHECK(dielectric.automatic_edges);
    CHECK(dielectric.edge_attributes.empty());
    CHECK(dielectric.edge_exclude_attributes.empty());
    CHECK_FALSE(dielectric.edge_frame_normal);
  }

  const auto geometry = ExtractMetalEdgeGeometry(
      *mesh, iodata.boundaries, JointNoiseExtractionFor(iodata.boundaries));
  REQUIRE(geometry.components == 3);
  REQUIRE(geometry.segments.size() == 92);
  int sa_segments = 0;
  int physical_segments = 0;
  int truncation_segments = 0;
  std::set<int> truncation_attributes;
  for (const auto &segment : geometry.segments)
  {
    REQUIRE(segment.conditions.size() == 1);
    CHECK(segment.conditions[0].type == MetalBoundaryConditionType::PEC);
    CHECK((segment.metal_attributes == std::vector<int>{1} ||
           segment.metal_attributes == std::vector<int>{2}));
    CHECK(segment.ms_interfaces == std::vector<int>{2});
    CHECK(segment.ma_interfaces == std::vector<int>{3});
    sa_segments += !segment.sa_interfaces.empty();
    physical_segments += segment.type == MetalEdgeSegmentType::PHYSICAL;
    truncation_segments += segment.type == MetalEdgeSegmentType::TRUNCATION;
    truncation_attributes.insert(segment.truncation_attributes.begin(),
                                 segment.truncation_attributes.end());
  }
  CHECK(sa_segments > 0);
  CHECK(geometry.physical_components == 4);
  CHECK(geometry.physical_chains == 4);
  CHECK(physical_segments == 8);
  CHECK(truncation_segments == 84);
  CHECK(truncation_attributes == std::set<int>{4, 5, 6});

  const auto sa_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry, 1, InterfaceDielectric::SA);
  const auto ms_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry, 2, InterfaceDielectric::MS);
  const auto ma_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry, 3, InterfaceDielectric::MA);
  CHECK(sa_indices.size() == 8);
  CHECK(ms_indices.size() == 8);
  CHECK(ma_indices.size() == 8);
  const auto sa_tree = BuildEdgeDistanceTree(geometry, sa_indices);
  CHECK(sa_tree->Size() == 8);
  CHECK(sa_tree->HasMetadata());
  for (std::size_t i = 0; i < sa_tree->Size(); i++)
  {
    CHECK(sa_tree->GetMetadata(i).component >= 0);
    CHECK(sa_tree->GetMetadata(i).chain >= 0);
  }
  const auto process_normals =
      BuildMetalEdgeProcessNormals(*mesh, geometry, sa_indices, [](int material_attribute)
                                   { return material_attribute == 2 ? 1.0 : 0.0; });
  for (const auto &normal : process_normals)
  {
    CHECK_THAT(normal[0], WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(normal[1], WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(normal[2], WithinAbs(0.0, 1.0e-12));
  }
  const auto gap_directions =
      BuildMetalEdgeGapDirections(*mesh, geometry, sa_indices, process_normals);
  REQUIRE(gap_directions.size() == sa_indices.size());
  for (std::size_t i = 0; i < gap_directions.size(); i++)
  {
    const auto &segment = geometry.segments[sa_indices[i]];
    const double x = 0.5 * (geometry.vertices[segment.vertices[0]].coordinate[0] +
                            geometry.vertices[segment.vertices[1]].coordinate[0]);
    const double expected_x = (x < 46.0 || (x > 62.0 && x < 78.0)) ? 1.0 : -1.0;
    CHECK_THAT(gap_directions[i][0], WithinAbs(expected_x, 1.0e-12));
    CHECK_THAT(gap_directions[i][1], WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(gap_directions[i][2], WithinAbs(0.0, 1.0e-12));
  }

  const auto contexts = BuildEdgeRefinementContexts(*mesh, iodata.boundaries);
  REQUIRE(contexts.size() == 1);
  CHECK(contexts[0].distance_tree->Size() == 8);

  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(mesh)));
  LaplaceOperator laplace_op(iodata, meshes);
  CHECK_NOTHROW((PostOperator<ProblemType::ELECTROSTATIC>(iodata, laplace_op)));
}

TEST_CASE("Automatic metal edge chains survive local refinement",
          "[geodata][metaledge][Serial][Parallel]")
{
  // P3-23: both crack forms (the legacy key of SD-801; upstream's default = the key unset).
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/cpw3d_surface/cpw3d_surface_validation_thin.json";
  IoData iodata(config_path.c_str(), false);
  iodata.model.mesh = fs::path(PALACE_TEST_DATA_DIR) / "mesh/cpw3d-surface-nc.msh";
  test::ApplyCrackPolicy(iodata, crack_policy);  // explicit crack policy (P3-23)
  iodata.model.refinement.max_it = 1;
  auto mesh = mesh::ReadMesh(iodata, Mpi::World());

  const auto coarse = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries,
                                               JointNoiseExtractionFor(iodata.boundaries));
  REQUIRE(coarse.physical_chains == 4);

  std::set<int> adjacent_elements;
  for (int be = 0; be < mesh->GetNBE(); be++)
  {
    const int attribute = mesh->GetBdrAttribute(be);
    if (attribute != 1 && attribute != 2)
    {
      continue;
    }
    int element, info;
    mesh->GetBdrElementAdjacentElement(be, element, info);
    adjacent_elements.insert(element);
  }
  // The boundary elements are rank-local: the precondition before the collective
  // GeneralRefinement is a global count (a rank with no metal boundary element marks none).
  long long adjacent_element_count = static_cast<long long>(adjacent_elements.size());
  Mpi::GlobalSum(1, &adjacent_element_count, mesh->GetComm());
  REQUIRE(adjacent_element_count > 0);
  mfem::Array<int> marked_elements(static_cast<int>(adjacent_elements.size()));
  std::copy(adjacent_elements.begin(), adjacent_elements.end(), marked_elements.begin());
  mesh->GeneralRefinement(marked_elements);

  const auto refined = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries,
                                                JointNoiseExtractionFor(iodata.boundaries));
  CHECK(refined.physical_components == coarse.physical_components);
  CHECK(refined.physical_chains == coarse.physical_chains);
  CHECK(refined.segments.size() > coarse.segments.size());

  auto GetPhysicalChainLengths = [](const MetalEdgeGeometry &geometry)
  {
    std::vector<double> lengths(geometry.physical_chains);
    for (const auto &segment : geometry.segments)
    {
      if (segment.type != MetalEdgeSegmentType::PHYSICAL)
      {
        continue;
      }
      const auto &p0 = geometry.vertices[segment.vertices[0]].coordinate;
      const auto &p1 = geometry.vertices[segment.vertices[1]].coordinate;
      double length_squared = 0.0;
      for (int d = 0; d < 3; d++)
      {
        length_squared += (p1[d] - p0[d]) * (p1[d] - p0[d]);
      }
      lengths[segment.physical_chain] += std::sqrt(length_squared);
    }
    std::sort(lengths.begin(), lengths.end());
    return lengths;
  };
  const auto coarse_lengths = GetPhysicalChainLengths(coarse);
  const auto refined_lengths = GetPhysicalChainLengths(refined);
  REQUIRE(refined_lengths.size() == coarse_lengths.size());
  for (std::size_t i = 0; i < coarse_lengths.size(); i++)
  {
    CHECK_THAT(refined_lengths[i], WithinAbs(coarse_lengths[i], 1.0e-10));
  }

  auto CountPhysicalVertexTypes = [](const MetalEdgeGeometry &geometry)
  {
    std::array<int, 4> counts{};
    for (const auto &vertex : geometry.vertices)
    {
      if (vertex.physical_type)
      {
        counts[static_cast<std::size_t>(*vertex.physical_type)]++;
      }
    }
    return counts;
  };
  const auto coarse_vertex_types = CountPhysicalVertexTypes(coarse);
  const auto refined_vertex_types = CountPhysicalVertexTypes(refined);
  for (std::size_t i = static_cast<std::size_t>(MetalEdgeVertexType::CORNER);
       i < coarse_vertex_types.size(); i++)
  {
    CHECK(refined_vertex_types[i] == coarse_vertex_types[i]);
  }
}

// Root-level perimeter (decision 899, OI-CTX-ARC-ADAPT): on a nonconforming mesh the
// identification's extraction merges the refinement-generated perimeter vertices out, so
// the segments, vertices and numbering are the root level's (the initial mesh of the run
// that refined it), and the frame builders resolve a leaf mesh edge to the root segment
// containing it: both crack copies of a sheet score the one root segment whatever their
// depth. On the unrefined mesh the option changes nothing (no refinement-generated vertex).
TEST_CASE("Root-level metal perimeter on a locally refined NC mesh",
          "[geodata][metaledge][Serial][Parallel]")
{
  // P3-23: both crack forms (the legacy key of SD-801; upstream's default = the key unset).
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/cpw3d_surface/cpw3d_surface_validation_thin.json";
  IoData iodata(config_path.c_str(), false);
  iodata.model.mesh = fs::path(PALACE_TEST_DATA_DIR) / "mesh/cpw3d-surface-nc.msh";
  test::ApplyCrackPolicy(iodata, crack_policy);  // explicit crack policy (P3-23)
  iodata.model.refinement.max_it = 1;
  auto mesh = mesh::ReadMesh(iodata, Mpi::World());
  REQUIRE(mesh->Nonconforming());

  MetalSurfaceExtraction leaf = JointNoiseExtractionFor(iodata.boundaries);
  MetalSurfaceExtraction root = leaf;
  root.root_level_perimeter = true;

  // Depth 0: the census is empty and the root-level extraction is the leaf extraction.
  const auto coarse = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries, leaf);
  const auto coarse_root = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries, root);
  REQUIRE(coarse.physical_chains == 4);
  CHECK(coarse_root.refinement_points.empty());
  CHECK(coarse_root.refinement_vertices_kept == 0);
  REQUIRE(coarse_root.segments.size() == coarse.segments.size());
  REQUIRE(coarse_root.vertices.size() == coarse.vertices.size());
  for (std::size_t i = 0; i < coarse.segments.size(); i++)
  {
    CHECK(coarse_root.segments[i].vertices == coarse.segments[i].vertices);
    CHECK(coarse_root.segments[i].type == coarse.segments[i].type);
  }

  // Refine the elements on ONE side of the two PEC sheets (the vacuum-side neighbours of
  // the sheet faces): the crack copies of the sheets are then at different depths and the
  // leaf perimeter is split at the refinement midpoints. P3-23: the neighbours are read
  // through the sheet face (GetFaceElements), not the boundary element's adjacent element:
  // on the cracked mesh each copy's face has the one element of its side (the same set as
  // GetBdrElementAdjacentElement); on the uncracked mesh (upstream's default) a sheet face
  // is interior and its boundary element names one neighbour only, the other side's
  // element is the face's second (local) neighbour.
  std::set<int> adjacent_elements;
  for (int be = 0; be < mesh->GetNBE(); be++)
  {
    const int attribute = mesh->GetBdrAttribute(be);
    if (attribute != 1 && attribute != 2)
    {
      continue;
    }
    int element1, element2;
    mesh->GetFaceElements(mesh->GetBdrElementFaceIndex(be), &element1, &element2);
    for (const int element : {element1, element2})
    {
      if (element >= 0 && mesh->GetAttribute(element) == 2)  // the vacuum side only
      {
        adjacent_elements.insert(element);
      }
    }
  }
  int marked_count = static_cast<int>(adjacent_elements.size());
  Mpi::GlobalSum(1, &marked_count, Mpi::World());
  REQUIRE(marked_count > 0);  // a rank may hold no vacuum-side element at the sheets
  mfem::Array<int> marked_elements(static_cast<int>(adjacent_elements.size()));
  std::copy(adjacent_elements.begin(), adjacent_elements.end(), marked_elements.begin());
  mesh->GeneralRefinement(marked_elements);

  const auto refined = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries, leaf);
  REQUIRE(refined.segments.size() > coarse.segments.size());
  const auto refined_root = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries, root);
  CHECK_FALSE(refined_root.refinement_points.empty());
  CHECK(refined_root.refinement_vertices_kept == 0);
  CHECK(refined_root.physical_chains == coarse.physical_chains);
  CHECK(refined_root.metal_components == coarse.metal_components);

  // The root-level segments ARE the unrefined mesh's: count, numbering, endpoints, types,
  // attributes and interfaces.
  REQUIRE(refined_root.segments.size() == coarse.segments.size());
  REQUIRE(refined_root.vertices.size() == coarse.vertices.size());
  for (std::size_t i = 0; i < coarse.segments.size(); i++)
  {
    const auto &a = coarse.segments[i];
    const auto &b = refined_root.segments[i];
    CHECK(a.vertices == b.vertices);
    CHECK(a.type == b.type);
    CHECK(a.metal_attributes == b.metal_attributes);
    CHECK(a.sa_interfaces == b.sa_interfaces);
    CHECK(a.ms_interfaces == b.ms_interfaces);
    CHECK(a.ma_interfaces == b.ma_interfaces);
    CHECK(a.metal_component == b.metal_component);
    CHECK(a.physical_chain == b.physical_chain);
    // The merged segment's face normals (the union of its children's leaf faces') are the
    // root edge's within the non-planar-metal classification's parallel tolerance: every
    // normal on either side has a parallel one on the other (face_count is a leaf-level
    // diagnostic without a consumer and is not compared).
    auto Parallel = [](const std::vector<std::array<double, 3>> &normals,
                       const std::vector<std::array<double, 3>> &others)
    {
      return std::all_of(normals.begin(), normals.end(),
                         [&](const std::array<double, 3> &n)
                         {
                           return std::any_of(others.begin(), others.end(),
                                              [&](const std::array<double, 3> &m)
                                              {
                                                return std::abs(n[0] * m[0] + n[1] * m[1] +
                                                                n[2] * m[2]) > 1.0 - 1.0e-8;
                                              });
                         });
    };
    REQUIRE(!a.face_normals.empty());
    REQUIRE(!b.face_normals.empty());
    CHECK(Parallel(a.face_normals, b.face_normals));
    CHECK(Parallel(b.face_normals, a.face_normals));
  }
  for (std::size_t v = 0; v < coarse.vertices.size(); v++)
  {
    for (int d = 0; d < 3; d++)
    {
      CHECK_THAT(refined_root.vertices[v].coordinate[d],
                 WithinAbs(coarse.vertices[v].coordinate[d], 1.0e-12));
    }
    CHECK(refined_root.vertices[v].segments == coarse.vertices[v].segments);
    CHECK(refined_root.vertices[v].type == coarse.vertices[v].type);
    CHECK(refined_root.vertices[v].physical_type == coarse.vertices[v].physical_type);
  }
  // Every refinement point lies strictly inside its segment.
  for (const auto &point : refined_root.refinement_points)
  {
    REQUIRE(point.segment < refined_root.segments.size());
    const auto &segment = refined_root.segments[point.segment];
    const auto &p0 = refined_root.vertices[segment.vertices[0]].coordinate;
    const auto &p1 = refined_root.vertices[segment.vertices[1]].coordinate;
    double length_squared = 0.0, projection = 0.0, distance_squared = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length_squared += (p1[d] - p0[d]) * (p1[d] - p0[d]);
      projection += (point.coordinate[d] - p0[d]) * (p1[d] - p0[d]);
    }
    const double t = projection / length_squared;
    CHECK(t > 1.0e-9);
    CHECK(t < 1.0 - 1.0e-9);
    for (int d = 0; d < 3; d++)
    {
      const double delta = point.coordinate[d] - p0[d] - t * (p1[d] - p0[d]);
      distance_squared += delta * delta;
    }
    CHECK(distance_squared < 1.0e-20);
  }

  // The frame builders on the refined mesh with the root geometry: every physical segment
  // sees both crack copies (no ambiguous process side) and reads the unrefined mesh's
  // process normal and gap direction; the leaf geometry on the same mesh leaves the
  // coarse copies' segments one-sided (the OI-CTX-ARC-ADAPT mechanism (A)).
  auto PhysicalSegments = [](const MetalEdgeGeometry &geometry)
  {
    std::vector<std::size_t> indices;
    for (std::size_t i = 0; i < geometry.segments.size(); i++)
    {
      if (geometry.segments[i].type == MetalEdgeSegmentType::PHYSICAL)
      {
        indices.push_back(i);
      }
    }
    return indices;
  };
  auto Score = [](int attribute) { return attribute == 1 ? 1.0 : 2.0; };
  auto coarse_mesh = mesh::ReadMesh(iodata, Mpi::World());
  const auto coarse_physical = PhysicalSegments(coarse);
  std::vector<bool> coarse_ambiguous, root_ambiguous, leaf_ambiguous;
  const auto coarse_normals = BuildMetalEdgeProcessNormals(
      *coarse_mesh, coarse, coarse_physical, Score, std::nullopt, &coarse_ambiguous);
  const auto coarse_gaps =
      BuildMetalEdgeGapDirections(*coarse_mesh, coarse, coarse_physical, coarse_normals);
  const auto root_normals = BuildMetalEdgeProcessNormals(
      *mesh, refined_root, coarse_physical, Score, std::nullopt, &root_ambiguous);
  const auto root_gaps =
      BuildMetalEdgeGapDirections(*mesh, refined_root, coarse_physical, root_normals);
  REQUIRE(root_normals.size() == coarse_normals.size());
  CHECK(std::count(coarse_ambiguous.begin(), coarse_ambiguous.end(), true) == 0);
  CHECK(std::count(root_ambiguous.begin(), root_ambiguous.end(), true) == 0);
  for (std::size_t i = 0; i < coarse_normals.size(); i++)
  {
    for (int d = 0; d < 3; d++)
    {
      CHECK_THAT(root_normals[i][d], WithinAbs(coarse_normals[i][d], 1.0e-12));
      CHECK_THAT(root_gaps[i][d], WithinAbs(coarse_gaps[i][d], 1.0e-12));
    }
  }
  const auto leaf_physical = PhysicalSegments(refined);
  BuildMetalEdgeProcessNormals(*mesh, refined, leaf_physical, Score, std::nullopt,
                               &leaf_ambiguous);
  // P3-23 (MEASURED): mechanism (A) needs two crack copies of a sheet refined to different
  // depths. Uncracked (upstream's default), the sheet is one set of interior faces that the
  // one-sided refinement splits for both sides at once: no leaf segment is one-sided, the
  // leaf frame is unambiguous as the root one.
  const auto leaf_ambiguous_count =
      std::count(leaf_ambiguous.begin(), leaf_ambiguous.end(), true);
  if (crack_policy == test::CrackPolicy::Legacy)
  {
    CHECK(leaf_ambiguous_count > 0);
  }
  else
  {
    CHECK(leaf_ambiguous_count == 0);
  }
}

TEST_CASE("Automatic metal edge extraction on 3D transmon",
          "[geodata][metaledge][Serial][Parallel]")
{
  // P3-23: both crack forms (the legacy key of SD-801; upstream's default = the key unset).
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/transmon/transmon_surface_coarse.json";
  IoData iodata(config_path.c_str(), false);
  iodata.model.mesh = config_path.parent_path() / iodata.model.mesh;
  test::ApplyCrackPolicy(iodata, crack_policy);  // explicit crack policy (P3-23)
  iodata.model.refinement.max_it = 0;
  auto mesh = mesh::ReadMesh(iodata, Mpi::World());

  for (const auto &[index, dielectric] : iodata.boundaries.postpro.dielectric)
  {
    (void)index;
    CHECK(dielectric.automatic_edges);
    CHECK(dielectric.edge_attributes.empty());
    CHECK(dielectric.edge_exclude_attributes.empty());
    CHECK_FALSE(dielectric.edge_frame_normal);
  }

  const auto geometry = ExtractMetalEdgeGeometry(
      *mesh, iodata.boundaries, JointNoiseExtractionFor(iodata.boundaries));
  REQUIRE_FALSE(geometry.Empty());
  int physical_segments = 0;
  int truncation_segments = 0;
  int fold_segments = 0;
  int nonmanifold_segments = 0, port_segments = 0;
  int embedded_segments = 0;
  int corners = 0;
  double fold_length = 0.0, embedded_length = 0.0;  // P3-23: the form-invariant measure
  auto SegmentLength = [&geometry](const MetalEdgeSegment &segment)
  {
    const auto &p0 = geometry.vertices[segment.vertices[0]].coordinate;
    const auto &p1 = geometry.vertices[segment.vertices[1]].coordinate;
    return std::sqrt((p1[0] - p0[0]) * (p1[0] - p0[0]) + (p1[1] - p0[1]) * (p1[1] - p0[1]) +
                     (p1[2] - p0[2]) * (p1[2] - p0[2]));
  };
  int sa_segments = 0;
  std::set<int> truncation_attributes;
  for (const auto &segment : geometry.segments)
  {
    CHECK(segment.metal_attributes == std::vector<int>{5});
    REQUIRE(segment.conditions.size() == 1);
    CHECK(segment.conditions[0].type == MetalBoundaryConditionType::PEC);
    CHECK(segment.ms_interfaces == std::vector<int>{2});
    CHECK(segment.ma_interfaces == std::vector<int>{3});
    physical_segments += segment.type == MetalEdgeSegmentType::PHYSICAL;
    truncation_segments += segment.type == MetalEdgeSegmentType::TRUNCATION;
    fold_segments += segment.type == MetalEdgeSegmentType::FOLD;
    nonmanifold_segments += segment.type == MetalEdgeSegmentType::NONMANIFOLD;
    port_segments += segment.type == MetalEdgeSegmentType::PORT;
    embedded_segments += segment.side_attributes.size() == 1;
    if (segment.type == MetalEdgeSegmentType::FOLD)
    {
      fold_length += SegmentLength(segment);
    }
    if (segment.side_attributes.size() == 1)
    {
      embedded_length += SegmentLength(segment);
    }
    sa_segments += !segment.sa_interfaces.empty();
    truncation_attributes.insert(segment.truncation_attributes.begin(),
                                 segment.truncation_attributes.end());
    if (!segment.sa_interfaces.empty())
    {
      CHECK(segment.type == MetalEdgeSegmentType::PHYSICAL);
    }
  }
  for (const auto &vertex : geometry.vertices)
  {
    corners += vertex.physical_type == std::optional{MetalEdgeVertexType::CORNER};
  }
  CAPTURE(geometry.components, geometry.physical_components, geometry.segments.size(),
          physical_segments, truncation_segments, fold_segments, nonmanifold_segments,
          embedded_segments, corners, sa_segments, truncation_attributes);
  // Phase 3 of the identification fix: the airbridge (attribute 5, vacuum on both sides)
  // is part of the perimeter (the crack-independent extraction no longer cancels a sheet
  // with one material on both sides): 40 one-sided edges (24 deck + 16 leg edges, single
  // side attribute), 16 deck-leg folds and 8 nonmanifold foot segments join the 3,088
  // planar physical and 58 truncation segments; 4 more perimeter components, 24 more
  // chains and 16 more corners (the deck's and legs' rectangles).
  // Decision 82(5): ports are not metal — the 6 feedline edges bordering the two lumped
  // port faces are PORT segments (cuts): the feedline centre conductor's perimeter splits
  // into 2 more physical components (the cut pieces), 6 fewer chains, and the 12 port-end
  // vertices are no longer corners.
  // USER decision 117(4) / 121 (B): the vertex classification uses the geometric joint
  // noise rule (kJointNoiseSagittaOverRadius at the config's edge distance R = 1.9 um)
  // instead of the 30 deg corner class: the polyline joints of the CPW bends (11-16 deg on
  // 7 um chords: implied sagitta 0.1-0.13 R) and rounded features (60 sharp corners before)
  // are CORNER vertices and chain breaks of the extraction (74 chains -> 460; the 1 deg
  // angular threshold of 117(4) read 461 chains / 447 corners: one joint whose implied
  // sagitta on its shorter piece is below 0.05 R is a straight continuation now); the
  // identification's arc rule absorbs the joints of fitted arcs and merges their chains.
  CHECK(geometry.components == 10);
  CHECK(geometry.physical_components == 15);
  CHECK(physical_segments + truncation_segments + fold_segments + nonmanifold_segments +
            port_segments ==
        static_cast<int>(geometry.segments.size()));
  CHECK(physical_segments == 3122);
  CHECK(port_segments == 6);
  // P3-23 (MEASURED, both forms): the 16 deck-leg folds of the legacy form are the 8 fold
  // mesh edges of the uncracked mesh, each split in two by the legacy crack's bisection
  // (Load, Model.RefineCrackElements: an element with an edge interior to the crack between
  // two seam vertices is refined); the fold length, 80, and the one-sided length, 496, are
  // the same in both forms, and so is everything else checked here.
  const bool legacy = crack_policy == test::CrackPolicy::Legacy;
  CHECK(fold_segments == (legacy ? 16 : 8));
  CHECK(nonmanifold_segments == 8);
  // The feet also touch the ground sheet faces.
  CHECK(embedded_segments == 40 + (legacy ? 16 : 8));
  CHECK_THAT(fold_length, WithinRel(80.0, 1.0e-12));
  CHECK_THAT(embedded_length, WithinRel(496.0, 1.0e-12));
  CHECK(geometry.physical_chains == 460);
  CHECK(truncation_segments == 58);
  CHECK(truncation_attributes == std::set<int>{3});
  CHECK(corners == 446);
  CHECK(sa_segments == 3082);
  // Conductors by metal connectivity: ground plane with the airbridge, feedline centre
  // conductor (between the lumped ports), island.
  CHECK(geometry.metal_components == 3);
  CHECK_THAT(geometry.layer_normal[2], WithinAbs(1.0, 1.0e-9));

  const auto sa_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry, 1, InterfaceDielectric::SA);
  const auto ms_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry, 2, InterfaceDielectric::MS);
  const auto ma_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry, 3, InterfaceDielectric::MA);
  CHECK(sa_indices.size() == 3082);
  CHECK(ms_indices.size() == 3122);
  CHECK(ma_indices.size() == 3122);

  // The 40 airbridge edges have vacuum on both sides: their process side is ambiguous and
  // (without a fallback) taken from the dominant normal component; every planar edge on the
  // substrate keeps its material-derived +z normal.
  std::vector<bool> ambiguous;
  const auto process_normals = BuildMetalEdgeProcessNormals(
      *mesh, geometry, ms_indices, [](int material_attribute)
      { return material_attribute == 2 ? 1.0 : 0.0; }, std::nullopt, &ambiguous);
  REQUIRE(process_normals.size() == ms_indices.size());
  REQUIRE(ambiguous.size() == ms_indices.size());
  int ambiguous_count = 0;
  for (std::size_t i = 0; i < ms_indices.size(); i++)
  {
    ambiguous_count += ambiguous[i];
    if (ambiguous[i])
    {
      continue;
    }
    const auto &normal = process_normals[i];
    CHECK_THAT(normal[0], WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(normal[1], WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(normal[2], WithinAbs(1.0, 1.0e-12));
  }
  CHECK(ambiguous_count == 40);
}

TEST_CASE("Automatic metal edge extraction treats a junction face as the open gap",
          "[geodata][metaledge][Serial][Parallel]")
{
  // Decision 876 (5) (DESIGN v3h/maxwell-ports A1.4 as fixed by decision 875): a chip's
  // junction is a small, purely reactive LumpedPort (R 0, L_J / C_J) whose face sits
  // between the junction leads of the island and the ground. The ES reference models the
  // same face as the open gap, so the identification treats a JUNCTION FACE like the gap it
  // replaced: the bordering lead-end segments stay PHYSICAL SA edges (no Port cut), and the
  // Maxwell run with the port, the ES comparator with the face listed in its SA attributes,
  // and the record's open gap see the same perimeter. An R 50 port and an oversize reactive
  // port (the face extent at or above the smallest EdgeDistances) stay Port cuts.
  //
  // Fixture: a hex box with interior faces on the plane z = 0.5 (cells 0.1 x 0.1): the
  // ground sheet (attribute 11, y in [0.1, 0.4]) and the pad (attribute 12, x in [0.3,
  // 0.6], y in [0.5, 0.8]) are PEC; the junction square (attribute 14) is the single cell
  // x in [0.4, 0.5], y in [0.4, 0.5] with one edge on the ground, one on the pad and two on
  // the gap (the chip form); the gap sheet (attribute 13) is the rest of the plane (the box
  // faces keep MakeCartesian3D's attributes 1-6).
  constexpr int n = 10;
  constexpr double h = 1.0 / n;
  constexpr double tolerance = 1.0e-12;
  mfem::Mesh serial =
      mfem::Mesh::MakeCartesian3D(n, n, 2, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  for (int e = 0; e < serial.GetNE(); e++)
  {
    mfem::Vector center(3);
    serial.GetElementCenter(e, center);
    serial.SetAttribute(e, center[2] < 0.5 ? 1 : 2);
  }
  for (int face = 0; face < serial.GetNumFaces(); face++)
  {
    int element1, element2;
    serial.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    serial.GetFaceVertices(face, vertices);
    bool on_plane = true;
    double xmin = 1.0, xmax = 0.0, ymin = 1.0, ymax = 0.0;
    for (const int vertex : vertices)
    {
      const double *point = serial.GetVertex(vertex);
      on_plane = on_plane && std::abs(point[2] - 0.5) < tolerance;
      xmin = std::min(xmin, point[0]);
      xmax = std::max(xmax, point[0]);
      ymin = std::min(ymin, point[1]);
      ymax = std::max(ymax, point[1]);
    }
    if (!on_plane)
    {
      continue;
    }
    auto Inside = [&](double x0, double x1, double y0, double y1)
    {
      return xmin >= x0 - tolerance && xmax <= x1 + tolerance && ymin >= y0 - tolerance &&
             ymax <= y1 + tolerance;
    };
    int attribute = 13;
    if (Inside(0.1, 0.9, 0.1, 0.4))
    {
      attribute = 11;
    }
    else if (Inside(0.3, 0.6, 0.5, 0.8))
    {
      attribute = 12;
    }
    else if (Inside(0.4, 0.5, 0.4, 0.5))
    {
      attribute = 14;
    }
    serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
    serial.SetBdrAttribute(serial.GetNBE() - 1, attribute);
  }
  serial.FinalizeTopology();
  serial.Finalize();
  mfem::ParMesh mesh(Mpi::World(), serial);

  // The junction square's extent (its bounding-box diagonal) is h sqrt(2) = 0.1414.
  const double junction_extent = h * std::sqrt(2.0);
  auto MakeBoundaries = [&](const std::vector<int> &sa_attributes,
                            const std::vector<double> &edge_distances,
                            std::optional<config::LumpedPortData> port)
  {
    config::BoundaryData boundaries;
    boundaries.pec.attributes = {11, 12};
    config::InterfaceDielectricData dielectric;
    dielectric.type = InterfaceDielectric::SA;
    dielectric.attributes = sa_attributes;
    dielectric.automatic_edges = true;
    dielectric.edge_distances = edge_distances;
    boundaries.postpro.dielectric.emplace(1, std::move(dielectric));
    if (port)
    {
      boundaries.lumpedport.emplace(1, std::move(*port));
    }
    return boundaries;
  };
  auto MakePort = [](double R, double L, double C)
  {
    config::LumpedPortData port;
    port.R = R;
    port.L = L;
    port.C = C;
    config::internal::ElementData element;
    element.direction = {0.0, 1.0, 0.0};
    element.attributes = {14};
    port.elements.push_back(std::move(element));
    return port;
  };
  // A lead-end segment borders the junction square: on the ground's edge y = 0.4 or the
  // pad's edge y = 0.5, with x in [0.4, 0.5].
  auto BordersJunction =
      [&](const MetalEdgeGeometry &geometry, const MetalEdgeSegment &segment)
  {
    const auto &p0 = geometry.vertices[segment.vertices[0]].coordinate;
    const auto &p1 = geometry.vertices[segment.vertices[1]].coordinate;
    const bool y_on_square =
        (std::abs(p0[1] - 0.4) < tolerance || std::abs(p0[1] - 0.5) < tolerance) &&
        std::abs(p1[1] - p0[1]) < tolerance;
    const bool x_on_square = std::min(p0[0], p1[0]) >= 0.4 - tolerance &&
                             std::max(p0[0], p1[0]) <= 0.5 + tolerance;
    return y_on_square && x_on_square;
  };
  // The perimeter as a canonical, partition-independent record: every segment's end
  // coordinates with its type and interface memberships, and every vertex's cut flags.
  struct SegmentRecord
  {
    std::array<double, 3> p0, p1;
    MetalEdgeSegmentType type;
    std::vector<int> sa_interfaces;
    std::vector<int> port_attributes;
    std::vector<int> truncation_attributes;
    bool on_port_boundary;
    bool operator==(const SegmentRecord &other) const
    {
      return p0 == other.p0 && p1 == other.p1 && type == other.type &&
             sa_interfaces == other.sa_interfaces &&
             port_attributes == other.port_attributes &&
             truncation_attributes == other.truncation_attributes &&
             on_port_boundary == other.on_port_boundary;
    }
  };
  auto Perimeter = [&](const MetalEdgeGeometry &geometry)
  {
    std::vector<SegmentRecord> records;
    for (const auto &segment : geometry.segments)
    {
      SegmentRecord record;
      record.p0 = geometry.vertices[segment.vertices[0]].coordinate;
      record.p1 = geometry.vertices[segment.vertices[1]].coordinate;
      if (record.p1 < record.p0)
      {
        std::swap(record.p0, record.p1);
      }
      record.type = segment.type;
      record.sa_interfaces = segment.sa_interfaces;
      record.port_attributes = segment.port_attributes;
      record.truncation_attributes = segment.truncation_attributes;
      record.on_port_boundary = geometry.vertices[segment.vertices[0]].on_port_boundary ||
                                geometry.vertices[segment.vertices[1]].on_port_boundary;
      records.push_back(std::move(record));
    }
    std::sort(records.begin(), records.end(), [](const auto &a, const auto &b)
              { return std::tie(a.p0, a.p1) < std::tie(b.p0, b.p1); });
    return records;
  };
  struct Counts
  {
    int physical = 0, port = 0, truncation = 0, port_vertices = 0, bordering = 0,
        bordering_physical_sa = 0;
  };
  auto Count = [&](const MetalEdgeGeometry &geometry)
  {
    Counts counts;
    for (const auto &segment : geometry.segments)
    {
      counts.physical += segment.type == MetalEdgeSegmentType::PHYSICAL;
      counts.port += segment.type == MetalEdgeSegmentType::PORT;
      counts.truncation += segment.type == MetalEdgeSegmentType::TRUNCATION;
      if (BordersJunction(geometry, segment))
      {
        counts.bordering++;
        counts.bordering_physical_sa += segment.type == MetalEdgeSegmentType::PHYSICAL &&
                                        segment.sa_interfaces == std::vector<int>{1};
      }
    }
    for (const auto &vertex : geometry.vertices)
    {
      counts.port_vertices += vertex.on_port_boundary;
    }
    return counts;
  };
  const MetalSurfaceExtraction surface = JointNoiseExtraction(0.3);

  // (0) The open gap of the ES reference and of the ES same-mesh comparator: the junction
  // face is listed in the SA interface's attributes, no LumpedPort.
  const auto control =
      ExtractMetalEdgeGeometry(mesh, MakeBoundaries({13, 14}, {0.3}, {}), surface);
  const auto control_counts = Count(control);
  REQUIRE(control.segments.size() == 2 * (8 + 3) + 2 * (3 + 3));
  CHECK(control_counts.physical == static_cast<int>(control.segments.size()));
  CHECK(control_counts.port == 0);
  CHECK(control_counts.truncation == 0);
  CHECK(control_counts.port_vertices == 0);
  CHECK(control_counts.bordering == 2);
  CHECK(control_counts.bordering_physical_sa == 2);
  CHECK(control.metal_components == 2);
  CHECK(control.physical_chains == 8);
  for (const auto &segment : control.segments)
  {
    CHECK(segment.sa_interfaces == std::vector<int>{1});
  }

  // (a) The Maxwell run: the junction face is a reactive LumpedPort below the smallest
  // EdgeDistances -> the same perimeter as the open gap.
  {
    REQUIRE(junction_extent < 0.3);
    const auto geometry = ExtractMetalEdgeGeometry(
        mesh, MakeBoundaries({13}, {0.3}, MakePort(0.0, 1.0e-9, 1.0e-15)), surface);
    const auto counts = Count(geometry);
    CHECK(counts.port == 0);
    CHECK(counts.port_vertices == 0);
    CHECK(counts.truncation == 0);
    CHECK(counts.bordering == 2);
    CHECK(counts.bordering_physical_sa == 2);
    CHECK(geometry.physical_chains == control.physical_chains);
    CHECK(Perimeter(geometry) == Perimeter(control));
  }
  // The excitation flag and an L-only / C-only junction play no part.
  {
    auto port = MakePort(0.0, 1.0e-9, 0.0);
    port.excitation = 1;
    const auto geometry = ExtractMetalEdgeGeometry(
        mesh, MakeBoundaries({13}, {0.3}, std::move(port)), surface);
    CHECK(Perimeter(geometry) == Perimeter(control));
    const auto capacitive = ExtractMetalEdgeGeometry(
        mesh, MakeBoundaries({13}, {0.3}, MakePort(0.0, 0.0, 1.0e-15)), surface);
    CHECK(Perimeter(capacitive) == Perimeter(control));
  }

  // (b) An R 50 port on the same face is a Port cut (decision 82(5), unchanged): the two
  // bordering lead-end segments are PORT with no SA membership, their 4 vertices port cuts.
  auto CheckPortCut = [&](const MetalEdgeGeometry &geometry)
  {
    const auto counts = Count(geometry);
    CHECK(counts.port == 2);
    CHECK(counts.port_vertices == 4);
    CHECK(counts.truncation == 0);
    CHECK(counts.bordering == 2);
    CHECK(counts.bordering_physical_sa == 0);
    CHECK(counts.physical == static_cast<int>(geometry.segments.size()) - 2);
    CHECK(geometry.physical_chains == control.physical_chains + 2);
    for (const auto &segment : geometry.segments)
    {
      if (segment.type == MetalEdgeSegmentType::PORT)
      {
        CHECK(BordersJunction(geometry, segment));
        CHECK(segment.port_attributes == std::vector<int>{14});
        CHECK(segment.sa_interfaces.empty());
      }
      else
      {
        CHECK(segment.sa_interfaces == std::vector<int>{1});
      }
    }
    CHECK_FALSE(Perimeter(geometry) == Perimeter(control));
  };
  CheckPortCut(ExtractMetalEdgeGeometry(
      mesh, MakeBoundaries({13}, {0.3}, MakePort(50.0, 0.0, 0.0)), surface));
  CheckPortCut(ExtractMetalEdgeGeometry(
      mesh, MakeBoundaries({13}, {0.3}, MakePort(50.0, 1.0e-9, 1.0e-15)), surface));

  // (c) A reactive port whose extent is not below the smallest EdgeDistances is a Port cut
  // (the size bound), as is one with no EdgeDistances configured (no radius to compare).
  REQUIRE(junction_extent > 0.1);
  CheckPortCut(ExtractMetalEdgeGeometry(
      mesh, MakeBoundaries({13}, {0.1}, MakePort(0.0, 1.0e-9, 1.0e-15)), surface));
  CheckPortCut(ExtractMetalEdgeGeometry(
      mesh, MakeBoundaries({13}, {}, MakePort(0.0, 1.0e-9, 1.0e-15)), surface));
}

TEST_CASE("Automatic metal edge classification is partition independent",
          "[geodata][metaledge][Serial][Parallel]")
{
  // P3-23: both crack forms (the legacy key of SD-801; upstream's default = the key unset).
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  // The metal perimeter and its edge frames must be a function of the global mesh only.
  // Palace's interior-boundary cracking displaces the duplicated crack vertices by
  // CrackDisplacementFactor x h_min, so coincident copies of a perimeter vertex carry
  // ulp-different coordinates; the canonical geometry gathered from all ranks must not
  // depend on which rank contributes which copy. Every rank loads the preprocessed serial
  // mesh, classifies it on MPI_COMM_SELF as the reference, and compares bit for bit
  // against two different distributions of the same mesh over MPI_COMM_WORLD.
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/transmon/transmon_surface_coarse.json";
  IoData iodata(config_path.c_str(), false);
  iodata.model.mesh = config_path.parent_path() / iodata.model.mesh;
  test::ApplyCrackPolicy(iodata, crack_policy);  // explicit crack policy (P3-23)
  iodata.model.refinement.max_it = 0;
  auto serial_mesh = mesh::Load(iodata, MPI_COMM_SELF);
  REQUIRE(serial_mesh);
  const int ne = serial_mesh->GetNE();
  const int size = Mpi::Size(Mpi::World());

  struct Classification
  {
    std::vector<mesh::BoundaryEdgeSegment> perimeter;
    MetalEdgeGeometry geometry;
    std::vector<std::array<double, 3>> process_normals;
    std::vector<std::array<double, 3>> gap_directions;
  };
  auto Classify = [&](mfem::ParMesh &mesh)
  {
    Classification result;
    auto marker = mesh::BdrAttrToMarker(mesh, std::vector<int>{5}, true);
    result.perimeter = mesh::GetBoundaryEdgeSegments(mesh, marker);
    MetalSurfaceExtraction surface = JointNoiseExtractionFor(iodata.boundaries);
    surface.classify_components = true;
    surface.retain_faces = true;
    result.geometry = ExtractMetalEdgeGeometry(mesh, iodata.boundaries, surface);
    const auto ms_indices =
        GetInterfaceMetalEdgeSegmentIndices(result.geometry, 2, InterfaceDielectric::MS);
    result.process_normals = BuildMetalEdgeProcessNormals(
        mesh, result.geometry, ms_indices,
        [](int material_attribute) { return material_attribute == 2 ? 1.0 : 0.0; });
    result.gap_directions = BuildMetalEdgeGapDirections(mesh, result.geometry, ms_indices,
                                                        result.process_normals);
    return result;
  };
  auto CompareToReference = [&](const Classification &reference,
                                const Classification &candidate, const char *label)
  {
    CAPTURE(label);
    REQUIRE(candidate.perimeter.size() == reference.perimeter.size());
    REQUIRE(candidate.geometry.vertices.size() == reference.geometry.vertices.size());
    REQUIRE(candidate.geometry.segments.size() == reference.geometry.segments.size());
    REQUIRE(candidate.process_normals.size() == reference.process_normals.size());
    REQUIRE(candidate.gap_directions.size() == reference.gap_directions.size());
    CHECK(candidate.geometry.components == reference.geometry.components);
    CHECK(candidate.geometry.physical_components == reference.geometry.physical_components);
    CHECK(candidate.geometry.physical_chains == reference.geometry.physical_chains);
    CHECK(candidate.geometry.metal_components == reference.geometry.metal_components);
    std::size_t perimeter_mismatches = 0, vertex_mismatches = 0, segment_mismatches = 0,
                normal_mismatches = 0, gap_mismatches = 0;
    for (std::size_t i = 0; i < reference.perimeter.size(); i++)
    {
      perimeter_mismatches += candidate.perimeter[i].p0 != reference.perimeter[i].p0 ||
                              candidate.perimeter[i].p1 != reference.perimeter[i].p1;
    }
    for (std::size_t i = 0; i < reference.geometry.vertices.size(); i++)
    {
      const auto &a = candidate.geometry.vertices[i];
      const auto &b = reference.geometry.vertices[i];
      vertex_mismatches += a.coordinate != b.coordinate || a.segments != b.segments ||
                           a.type != b.type || a.physical_type != b.physical_type ||
                           a.on_truncation_boundary != b.on_truncation_boundary;
    }
    for (std::size_t i = 0; i < reference.geometry.segments.size(); i++)
    {
      const auto &a = candidate.geometry.segments[i];
      const auto &b = reference.geometry.segments[i];
      segment_mismatches +=
          a.vertices != b.vertices || a.component != b.component ||
          a.physical_component != b.physical_component ||
          a.physical_chain != b.physical_chain || a.metal_component != b.metal_component ||
          a.type != b.type || a.metal_attributes != b.metal_attributes ||
          a.sa_interfaces != b.sa_interfaces || a.ms_interfaces != b.ms_interfaces ||
          a.ma_interfaces != b.ma_interfaces ||
          a.truncation_attributes != b.truncation_attributes;
    }
    for (std::size_t i = 0; i < reference.process_normals.size(); i++)
    {
      normal_mismatches += candidate.process_normals[i] != reference.process_normals[i];
      gap_mismatches += candidate.gap_directions[i] != reference.gap_directions[i];
    }
    CAPTURE(perimeter_mismatches, vertex_mismatches, segment_mismatches, normal_mismatches,
            gap_mismatches);
    CHECK(perimeter_mismatches == 0);
    CHECK(vertex_mismatches == 0);
    CHECK(segment_mismatches == 0);
    CHECK(normal_mismatches == 0);
    CHECK(gap_mismatches == 0);
  };

  Classification reference;
  {
    mfem::ParMesh mesh(MPI_COMM_SELF, *serial_mesh);
    reference = Classify(mesh);
  }
  REQUIRE(reference.geometry.segments.size() > 3000);
  REQUIRE(reference.process_normals.size() ==
          3122);  // 3,082 planar + 40 airbridge edges (6 port cuts)

  // Round-robin element distribution: every rank owns crack copies from everywhere.
  std::vector<int> round_robin(ne);
  for (int e = 0; e < ne; e++)
  {
    round_robin[e] = e % size;
  }
  {
    mfem::ParMesh mesh(Mpi::World(), *serial_mesh, round_robin.data());
    CompareToReference(reference, Classify(mesh), "round-robin partition");
  }

  // Contiguous blocks in reversed rank order: the gathered contribution order differs
  // from both the serial order and the round-robin order.
  const int block = std::max(1, (ne + size - 1) / size);
  std::vector<int> reversed_blocks(ne);
  for (int e = 0; e < ne; e++)
  {
    reversed_blocks[e] = size - 1 - std::min(size - 1, e / block);
  }
  {
    mfem::ParMesh mesh(Mpi::World(), *serial_mesh, reversed_blocks.data());
    CompareToReference(reference, Classify(mesh), "reversed block partition");
  }
}

TEST_CASE("Boundary edge extraction ignores coincident crack copies", "[geodata][Serial]")
{
  auto serial_mesh =
      mfem::Mesh::MakeCartesian2D(4, 2, mfem::Element::TRIANGLE, false, 2.0, 1.0);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  const int original_nbe = mesh.GetNBE();
  mfem::Array<int> first_vertices;
  for (int be = 0; be < original_nbe; be++)
  {
    if (mesh.GetBdrAttribute(be) != 1)
    {
      continue;
    }
    mesh.AddBdrElement(mesh.GetBdrElement(be)->Duplicate(&mesh));
    if (first_vertices.IsEmpty())
    {
      mesh.GetBdrElementVertices(be, first_vertices);
    }
  }
  mfem::Vector midpoint(2);
  for (int d = 0; d < 2; d++)
  {
    midpoint[d] =
        0.5 * (mesh.GetVertex(first_vertices[0])[d] + mesh.GetVertex(first_vertices[1])[d]);
  }
  const int mid = mesh.AddVertex(midpoint);
  mesh.AddBdrElement(new mfem::Segment(first_vertices[0], mid, 1));
  mesh.AddBdrElement(new mfem::Segment(mid, first_vertices[1], 1));

  auto marker = mesh::AttrToMarker(mesh.bdr_attributes.Max(), std::vector<int>{1});
  auto edges = mesh::GetBoundaryEdgeSegments(mesh, marker);

  REQUIRE(edges.size() == 2);
  CHECK_THAT(edges[0].p0[0] + edges[1].p0[0], WithinAbs(2.0, 1.0e-12));
}

TEST_CASE("Boundary edge extraction on cracked CPW mesh", "[geodata][Serial]")
{
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/cpw2d/cpw2d_thin_electrostatic.json";
  IoData iodata(config_path.c_str(), false);
  iodata.model.mesh = config_path.parent_path() / iodata.model.mesh;
  iodata.model.crack_bdr_elements = true;  // explicit crack policy: legacy (cracked)
  auto mesh = mesh::ReadMesh(iodata, Mpi::World());

  auto marker = mesh::AttrToMarker(mesh->bdr_attributes.Max(), std::vector<int>{1, 2});
  auto edges = mesh::GetBoundaryEdgeSegments(*mesh, marker);
  std::vector<double> x;
  x.reserve(edges.size());
  for (const auto &edge : edges)
  {
    CHECK(edge.p0 == edge.p1);
    CHECK_THAT(edge.p0[1], WithinAbs(0.0, 1.0e-12));
    x.push_back(edge.p0[0]);
  }
  std::sort(x.begin(), x.end());

  REQUIRE(x.size() == 6);
  const std::array<double, 6> expected = {0.0, 500.0, 503.5, 518.5, 522.0, 1022.0};
  for (std::size_t i = 0; i < x.size(); i++)
  {
    CHECK_THAT(x[i], WithinAbs(expected[i], 1.0e-12));
  }
}

// P3-23 (the sibling of the case above for upstream's default crack form, SL8-RESOLUTION s4
// (a)): the same 2D CPW metal line read with Model.CrackInternalBoundaryElements unset (the
// Ground / Terminal sheets 1 / 2 are then NOT cracked) and, in the same case, with the
// legacy key: the premise (the cracked mesh carries the duplicated sheet vertices, the
// uncracked one does not) and the extraction's result, order-free: the same six metal end
// points on y = 0.
TEST_CASE("Boundary edge extraction on uncracked CPW mesh", "[geodata][Serial]")
{
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/cpw2d/cpw2d_thin_electrostatic.json";
  auto EdgePoints = [&](test::CrackPolicy crack_policy, int &vertices)
  {
    IoData iodata(config_path.c_str(), false);
    iodata.model.mesh = config_path.parent_path() / iodata.model.mesh;
    test::ApplyCrackPolicy(iodata, crack_policy);
    auto mesh = mesh::ReadMesh(iodata, Mpi::World());
    // MEASURED: the legacy key cracks the Ground (PEC) sheet 2 only: Terminal attributes
    // are not in BoundaryData::attributes, the legacy crack set, so the Terminal sheet 1 is
    // not cracked in either form.
    CHECK(std::set<int>(iodata.boundaries.cracked_attributes.begin(),
                        iodata.boundaries.cracked_attributes.end()) ==
          (crack_policy == test::CrackPolicy::Legacy ? std::set<int>{2} : std::set<int>{}));
    vertices = mesh->GetNV();
    auto marker = mesh::AttrToMarker(mesh->bdr_attributes.Max(), std::vector<int>{1, 2});
    std::vector<std::array<double, 2>> points;
    for (const auto &edge : mesh::GetBoundaryEdgeSegments(*mesh, marker))
    {
      CHECK(edge.p0 == edge.p1);
      points.push_back({edge.p0[0], edge.p0[1]});
    }
    std::sort(points.begin(), points.end());
    return points;
  };
  int cracked_vertices = 0, uncracked_vertices = 0;
  const auto cracked = EdgePoints(test::CrackPolicy::Legacy, cracked_vertices);
  const auto uncracked = EdgePoints(test::CrackPolicy::Default, uncracked_vertices);
  CHECK(uncracked_vertices < cracked_vertices);
  REQUIRE(uncracked.size() == cracked.size());
  REQUIRE(uncracked.size() == 6);
  for (std::size_t i = 0; i < uncracked.size(); i++)
  {
    CHECK_THAT(uncracked[i][0], WithinAbs(cracked[i][0], 1.0e-12));
    CHECK_THAT(uncracked[i][1], WithinAbs(0.0, 1.0e-12));
  }
}

TEST_CASE("Boundary edge extraction on NC 3D CPW mesh", "[geodata][Serial][Parallel]")
{
  // P3-23: both crack forms (the legacy key of SD-801; upstream's default = the key unset).
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/cpw3d_surface/cpw3d_surface_validation_thin.json";
  auto ExtractEdges = [&](bool use_amr)
  {
    IoData iodata(config_path.c_str(), false);
    iodata.model.mesh = fs::path(PALACE_TEST_DATA_DIR) / "mesh/cpw3d-surface-nc.msh";
    test::ApplyCrackPolicy(iodata, crack_policy);  // explicit crack policy (P3-23)
    iodata.model.refinement.max_it = use_amr ? 1 : 0;
    auto mesh = mesh::ReadMesh(iodata, Mpi::World());
    auto marker = mesh::AttrToMarker(mesh->bdr_attributes.Max(), std::vector<int>{1, 2});
    return mesh::GetBoundaryEdgeSegments(*mesh, marker);
  };

  auto conforming_edges = ExtractEdges(false);
  auto nc_edges = ExtractEdges(true);
  auto TotalLength = [](const auto &edges)
  {
    double length = 0.0;
    for (const auto &edge : edges)
    {
      double length_squared = 0.0;
      for (int d = 0; d < 3; d++)
      {
        const double delta = edge.p1[d] - edge.p0[d];
        length_squared += delta * delta;
      }
      length += std::sqrt(length_squared);
    }
    return length;
  };

  CAPTURE(conforming_edges.size(), nc_edges.size(), TotalLength(conforming_edges),
          TotalLength(nc_edges));
  REQUIRE(conforming_edges.size() == 92);
  CHECK_THAT(TotalLength(conforming_edges), WithinAbs(320.0, 1.0e-10));
  for (const auto &nc_edge : nc_edges)
  {
    bool found = false;
    for (const auto &edge : conforming_edges)
    {
      double direct = 0.0, reverse = 0.0;
      for (int d = 0; d < 3; d++)
      {
        direct = std::max(direct, std::abs(nc_edge.p0[d] - edge.p0[d]));
        direct = std::max(direct, std::abs(nc_edge.p1[d] - edge.p1[d]));
        reverse = std::max(reverse, std::abs(nc_edge.p0[d] - edge.p1[d]));
        reverse = std::max(reverse, std::abs(nc_edge.p1[d] - edge.p0[d]));
      }
      found = found || std::min(direct, reverse) < 1.0e-10;
    }
    CAPTURE(nc_edge.p0[0], nc_edge.p0[1], nc_edge.p0[2], nc_edge.p1[0], nc_edge.p1[1],
            nc_edge.p1[2]);
    CHECK(found);
  }
  CHECK(nc_edges.size() == conforming_edges.size());
  CHECK_THAT(TotalLength(nc_edges), WithinAbs(TotalLength(conforming_edges), 1.0e-10));
}

TEST_CASE("Boundary edge extraction on NC mesh", "[geodata][Serial]")
{
  SECTION("2D endpoints")
  {
    auto serial_mesh =
        mfem::Mesh::MakeCartesian2D(4, 2, mfem::Element::TRIANGLE, false, 2.0, 1.0);
    serial_mesh.EnsureNCMesh(true);
    mfem::ParMesh mesh(Mpi::World(), serial_mesh);
    mfem::Array<int> marked_elements(1);
    marked_elements[0] = 0;
    mesh.GeneralRefinement(marked_elements);

    REQUIRE(mesh.Nonconforming());
    auto marker = mesh::AttrToMarker(mesh.bdr_attributes.Max(), std::vector<int>{1});
    auto edges = mesh::GetBoundaryEdgeSegments(mesh, marker);

    REQUIRE(edges.size() == 2);
    for (const auto &edge : edges)
    {
      CHECK(edge.p0 == edge.p1);
      CHECK_THAT(edge.p0[1], WithinAbs(0.0, 1.0e-12));
      CHECK_THAT(edge.p0[2], WithinAbs(0.0, 1.0e-12));
    }
    CHECK_THAT(edges[0].p0[0] + edges[1].p0[0], WithinAbs(2.0, 1.0e-12));
  }

  SECTION("3D perimeter")
  {
    auto serial_mesh =
        mfem::Mesh::MakeCartesian3D(2, 2, 1, mfem::Element::HEXAHEDRON, 2.0, 2.0, 1.0);
    serial_mesh.EnsureNCMesh(true);
    mfem::ParMesh mesh(Mpi::World(), serial_mesh);
    mfem::Array<int> marked_elements(1);
    marked_elements[0] = 0;
    mesh.GeneralRefinement(marked_elements);

    REQUIRE(mesh.Nonconforming());
    auto marker = mesh::AttrToMarker(mesh.bdr_attributes.Max(), std::vector<int>{1});
    auto edges = mesh::GetBoundaryEdgeSegments(mesh, marker);
    double length = 0.0;
    for (const auto &edge : edges)
    {
      double length_squared = 0.0;
      for (int d = 0; d < 3; d++)
      {
        const double delta = edge.p1[d] - edge.p0[d];
        length_squared += delta * delta;
      }
      length += std::sqrt(length_squared);
    }
    CHECK_THAT(length, WithinAbs(8.0, 1.0e-12));
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Mesh loader removes boundary elements which are not faces of the mesh",
                 "[geodata][Serial]")
{
  // Gmsh mesh of the unit square, with two triangles split along the diagonal from (0, 0)
  // to (1, 1), and four boundary elements on its boundary. Three more boundary elements are
  // not an edge of any triangle: the other diagonal of the square (attribute 2), the edge
  // from (2, 0) to (2, 1) whose vertices are not vertices of any triangle (attribute 3),
  // and the edge from (1, 0) to (2, 0) with the attribute 1 of the square boundary. A
  // second order mesh has a curved edge from (0, 0) to (1, 0) through (0.5, -0.1), which
  // adds 2/3 * 0.1 to the area of the square.
  const int order = GENERATE(1, 2);
  const bool clean_unused_elements = GENERATE(false, true);
  CAPTURE(order, clean_unused_elements);
  std::vector<std::array<double, 2>> nodes = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0},
                                              {0.0, 1.0}, {2.0, 0.0}, {2.0, 1.0}};
  std::vector<std::array<int, 3>> edge_nodes;  // Vertices and second order node of edges
  auto EdgeNode = [&](int v0, int v1)
  {
    for (const auto &[w0, w1, node] : edge_nodes)
    {
      if ((w0 == v0 && w1 == v1) || (w0 == v1 && w1 == v0))
      {
        return node;
      }
    }
    edge_nodes.push_back({v0, v1, static_cast<int>(nodes.size())});
    nodes.push_back({0.5 * (nodes[v0][0] + nodes[v1][0]),
                     0.5 * (nodes[v0][1] + nodes[v1][1]) - (v0 + v1 == 1 ? 0.1 : 0.0)});
    return edge_nodes.back()[2];
  };
  struct GmshElement
  {
    int type, attr;
    std::vector<int> nodes;
  };
  std::vector<GmshElement> elements;
  for (auto [v0, v1, attr] : std::vector<std::array<int, 3>>{
           {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1}, {1, 3, 2}, {4, 5, 3}, {1, 4, 1}})
  {
    elements.push_back({order == 1 ? 1 : 8, attr,
                        order == 1 ? std::vector<int>{v0, v1}
                                   : std::vector<int>{v0, v1, EdgeNode(v0, v1)}});
  }
  for (auto [v0, v1, v2] : std::vector<std::array<int, 3>>{{0, 1, 2}, {0, 2, 3}})
  {
    elements.push_back({order == 1 ? 2 : 9, 1,
                        order == 1 ? std::vector<int>{v0, v1, v2}
                                   : std::vector<int>{v0, v1, v2, EdgeNode(v0, v1),
                                                      EdgeNode(v1, v2), EdgeNode(v2, v0)}});
  }
  const auto mesh_file = temp_dir / "mesh.msh";
  {
    std::ofstream fo(mesh_file);
    fo << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n" << nodes.size() << "\n";
    for (std::size_t i = 0; i < nodes.size(); i++)
    {
      fo << i + 1 << " " << nodes[i][0] << " " << nodes[i][1] << " 0\n";
    }
    fo << "$EndNodes\n$Elements\n" << elements.size() << "\n";
    for (std::size_t i = 0; i < elements.size(); i++)
    {
      fo << i + 1 << " " << elements[i].type << " 2 " << elements[i].attr << " "
         << elements[i].attr;
      for (auto v : elements[i].nodes)
      {
        fo << " " << v + 1;
      }
      fo << "\n";
    }
    fo << "$EndElements\n";
  }

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.mesh = mesh_file.string();
  iodata.model.clean_unused_elements = clean_unused_elements;
  iodata.domains.attributes = {1};
  auto mesh = mesh::Load(iodata, Mpi::World());
  REQUIRE(mesh);

  // The boundary elements which are not edges of the mesh are removed, together with the
  // vertices of only these boundary elements.
  CHECK(mesh->GetNE() == 2);
  CHECK(mesh->GetNV() == 4);
  REQUIRE(mesh->GetNBE() == 4);
  REQUIRE(mesh->bdr_attributes.Size() == 1);
  CHECK(mesh->bdr_attributes[0] == 1);
  double length = 0.0;
  for (int be = 0; be < mesh->GetNBE(); be++)
  {
    CAPTURE(be);
    int e1, e2;
    REQUIRE(mesh->GetBdrElementFaceIndex(be) >= 0);
    mesh->GetFaceElements(mesh->GetBdrElementFaceIndex(be), &e1, &e2);
    CHECK(e1 >= 0);
    CHECK(e2 < 0);
    auto *T = mesh->GetBdrElementTransformation(be);
    const auto &ir = mfem::IntRules.Get(T->GetGeometryType(), 10);
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      T->SetIntPoint(&ir.IntPoint(q));
      length += ir.IntPoint(q).weight * T->Weight();
    }
  }
  CHECK(mesh->CheckBdrElementOrientation(false) == 0);

  // The geometry of the mesh is unchanged.
  REQUIRE(mesh->GetNodes());
  CHECK(mesh->GetNodes()->FESpace()->GetMaxElementOrder() == order);
  double area = 0.0;
  for (int e = 0; e < mesh->GetNE(); e++)
  {
    area += mesh->GetElementVolume(e);
  }
  CHECK_THAT(area, WithinRel(order == 1 ? 1.0 : 1.0 + 0.2 / 3.0, 1.0e-12));
  const double curved_length = 0.5 * std::sqrt(1.16) + std::asinh(0.4) / 0.8;
  CHECK_THAT(length, WithinRel(order == 1 ? 4.0 : 3.0 + curved_length, 1.0e-8));
}

TEST_CASE("PeriodicGmsh", "[geodata][Serial]")
{
  auto torus_path = fs::path(PALACE_TEST_DATA_DIR) / "mesh" / "periodic-torus-sector.msh";
  std::ifstream fi(torus_path.string());
  std::unique_ptr<mfem::Mesh> mesh = std::make_unique<mfem::Mesh>(fi, false, false, true);

  json boundaries = {
      {"Periodic",
       {{"BoundaryPairs", {{{"DonorAttributes", {1}}, {"ReceiverAttributes", {2}}}}}}}};
  config::BoundaryData boundary_torus(boundaries);

  for (const auto &data : boundary_torus.periodic.boundary_pairs)
  {
    auto periodic_mapping = mesh::DeterminePeriodicVertexMapping(mesh, data);
    REQUIRE(periodic_mapping.empty());
  }

  // Round-trip the periodic mesh through the nonconformal (AMR) save format. The NCMesh
  // conversion stores a single boundary attribute per face, so the donor/receiver pair
  // sharing each seam face collapses to one attribute; the reloaded mesh must still be
  // detected as already periodic.
  mesh->EnsureNCMesh(true);
  std::stringstream buffer;
  mesh->Print(buffer);
  buffer.seekg(0);
  mesh = std::make_unique<mfem::Mesh>(buffer, false, true, false);
  REQUIRE((mesh->bdr_attributes.Find(1) >= 0) != (mesh->bdr_attributes.Find(2) >= 0));
  REQUIRE(mesh::DeterminePeriodicVertexMapping(
              mesh, boundary_torus.periodic.boundary_pairs.front())
              .empty());
}

TEST_CASE("Interior boundary mesh cracking", "[geodata][Serial]")
{
  // Unit cube hexahedral mesh with interior boundary elements (attribute 7) on z = 0.5,
  // spanning the domain (exterior boundary attributes 1-6).
  auto base = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON);
  std::vector<int> sheet_faces;
  mfem::Array<int> fv;
  for (int f = 0; f < base.GetNumFaces(); f++)
  {
    int e1, e2;
    base.GetFaceElements(f, &e1, &e2);
    base.GetFaceVertices(f, fv);
    bool on_sheet = (e2 >= 0);
    for (auto v : fv)
    {
      on_sheet = on_sheet && std::abs(base.GetVertex(v)[2] - 0.5) < 1.0e-12;
    }
    if (on_sheet)
    {
      sheet_faces.push_back(f);
    }
  }
  REQUIRE(sheet_faces.size() == 4);
  mfem::Mesh mesh(3, base.GetNV(), base.GetNE(),
                  base.GetNBE() + static_cast<int>(sheet_faces.size()));
  for (int v = 0; v < base.GetNV(); v++)
  {
    mesh.AddVertex(base.GetVertex(v));
  }
  for (int e = 0; e < base.GetNE(); e++)
  {
    mesh.AddElement(base.GetElement(e)->Duplicate(&mesh));
  }
  for (int be = 0; be < base.GetNBE(); be++)
  {
    mesh.AddBdrElement(base.GetBdrElement(be)->Duplicate(&mesh));
  }
  for (auto f : sheet_faces)
  {
    auto *el = base.GetFace(f)->Duplicate(&mesh);
    el->SetAttribute(7);
    mesh.AddBdrElement(el);
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  mesh.SetAttributes();
  const auto mesh_path = fs::temp_directory_path() / "palace-test-interior-sheet.mesh";
  {
    std::ofstream fo(mesh_path);
    mesh.Print(fo);
  }
  const int nv = mesh.GetNV(), nv_cut = nv + 9;  // 3 x 3 vertices on the sheet

  auto LoadMesh = [&](const json &boundaries, const json &model, const fs::path &path)
  {
    json config = {{"Problem", {{"Type", "Eigenmode"}, {"Verbose", 0}}},
                   {"Model", {{"Mesh", path.string()}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", boundaries},
                   {"Solver", {{"Eigenmode", {{"Target", 1.0}}}}}};
    config["Model"].update(model);
    IoData iodata(config, false);
    return mesh::Load(iodata, Mpi::World());
  };
  auto NumVertices = [&](const json &boundaries, const json &model)
  { return LoadMesh(boundaries, model, mesh_path)->GetNV(); };
  const json none = json::object();
  const json legacy = {{"CrackInternalBoundaryElements", true}};
  const json off = {{"CrackInternalBoundaryElements", false}};
  auto Impedance = [](json crack)
  {
    json data = {{"Attributes", {7}}, {"Ls", 1.0e-12}};
    if (!crack.is_null())
    {
      data["Crack"] = crack;
    }
    return json{{"Impedance", {data}}};
  };
  auto Conductivity = [](json crack)
  {
    json data = {{"Attributes", {7}}, {"Conductivity", 1.0e7}};
    if (!crack.is_null())
    {
      data["Crack"] = crack;
    }
    return json{{"Conductivity", {data}}};
  };
  auto RationalImpedance = [](json crack)
  {
    json data = {{"Attributes", {7}}, {"Numerator", {1.0}}, {"Denominator", {1.0}}};
    if (!crack.is_null())
    {
      data["Crack"] = crack;
    }
    return json{{"RationalImpedance", {data}}};
  };

  // Sheets with a surface impedance: thin film (not cracked) or thick conductor (cracked).
  CHECK(NumVertices(Impedance(nullptr), none) == nv);
  CHECK(NumVertices(Impedance(false), none) == nv);
  CHECK(NumVertices(Impedance(true), none) == nv_cut);
  CHECK(NumVertices(Conductivity(nullptr), none) == nv);
  CHECK(NumVertices(Conductivity(true), none) == nv_cut);
  CHECK(NumVertices(RationalImpedance(nullptr), none) == nv);
  CHECK(NumVertices(RationalImpedance(true), none) == nv_cut);

  // Perfect conductors are never cracked, one-sided boundary conditions always are.
  CHECK(NumVertices({{"PEC", {{"Attributes", {7}}}}}, none) == nv);
  CHECK(NumVertices({{"PMC", {{"Attributes", {7}}}}}, none) == nv_cut);

  // Deprecated global option. Disabling mesh cracking contradicts boundaries which request
  // it.
  CHECK(NumVertices({{"PEC", {{"Attributes", {7}}}}}, legacy) == nv_cut);
  CHECK(NumVertices(Impedance(false), legacy) == nv_cut);
  CHECK(NumVertices(Impedance(true), legacy) == nv_cut);
  CHECK(NumVertices(Conductivity(nullptr), off) == nv);
  CHECK(NumVertices(Impedance(false), off) == nv);
  CHECK(NumVertices({{"PMC", {{"Attributes", {7}}}}}, off) == nv);
  CHECK_THROWS(NumVertices(Impedance(true), off));
  CHECK_THROWS(NumVertices(Conductivity(true), off));
  CHECK_THROWS(NumVertices(RationalImpedance(true), off));

  // Each face of a cracked impedance boundary has the given surface impedance, except with
  // the deprecated global option, which splits it between the two faces.
  auto SplitAttributes = [&](const json &boundaries, const json &model)
  {
    json config = {{"Problem", {{"Type", "Eigenmode"}, {"Verbose", 0}}},
                   {"Model", {{"Mesh", mesh_path.string()}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                   {"Boundaries", boundaries},
                   {"Solver", {{"Eigenmode", {{"Target", 1.0}}}}}};
    config["Model"].update(model);
    IoData iodata(config, false);
    mesh::ReadMesh(iodata, Mpi::World());
    REQUIRE(iodata.boundaries.cracked_attributes == std::unordered_set<int>{7});
    return iodata.boundaries.split_impedance_attributes;
  };
  CHECK(SplitAttributes(Impedance(true), none).empty());
  CHECK(SplitAttributes(RationalImpedance(true), none).empty());
  CHECK(SplitAttributes(Impedance(false), legacy) == std::unordered_set<int>{7});
  CHECK(SplitAttributes(Impedance(true), legacy) == std::unordered_set<int>{7});

  // Refinement boxes make the hexahedral mesh nonconformal, which cannot be cracked.
  const json box = {{"Refinement",
                     {{"Boxes",
                       {{{"Levels", 1},
                         {"BoundingBoxMin", {0.0, 0.0, 0.0}},
                         {"BoundingBoxMax", {0.5, 0.5, 1.0}}}}}}}};
  CHECK_THROWS(NumVertices(Impedance(true), box));
  CHECK_THROWS(NumVertices({{"PMC", {{"Attributes", {7}}}}}, box));
  CHECK_NOTHROW(NumVertices(Impedance(false), box));
  CHECK_NOTHROW(NumVertices(
      {{"PMC", {{"Attributes", {7}}}}},
      {{"CrackInternalBoundaryElements", false}, {"Refinement", box["Refinement"]}}));

  // A nonconformal input mesh, for example from a previous adaptive simulation, cannot be
  // cracked: it must already be cracked along the boundaries which are to be cracked.
  // Refine the element below the sheet next to the origin, such that it has hanging faces.
  auto WriteNonconformal = [](mfem::Mesh &&m, const fs::path &path)
  {
    m.EnsureNCMesh(true);
    mfem::Array<int> marked;
    mfem::Vector c;
    for (int e = 0; e < m.GetNE(); e++)
    {
      m.GetElementCenter(e, c);
      if (c(0) < 0.5 && c(1) < 0.5 && c(2) < 0.5)
      {
        marked.Append(e);
      }
    }
    REQUIRE(marked.Size() == 1);
    m.GeneralRefinement(marked, 1, 0);
    std::ofstream fo(path);
    fo.precision(17);
    m.Print(fo);
    return m.GetNV();
  };
  const auto nc_path = fs::temp_directory_path() / "palace-test-interior-sheet-nc.mesh";
  const auto nc_cut_path =
      fs::temp_directory_path() / "palace-test-interior-sheet-nc-cracked.mesh";
  WriteNonconformal(mfem::Mesh(mesh), nc_path);
  const int nv_nc_cut = WriteNonconformal(
      std::move(*LoadMesh(Impedance(true), none, mesh_path)), nc_cut_path);
  CHECK_THROWS(LoadMesh(Impedance(true), none, nc_path));
  CHECK_THROWS(LoadMesh(Conductivity(true), none, nc_path));
  CHECK_THROWS(LoadMesh({{"PMC", {{"Attributes", {7}}}}}, none, nc_path));
  CHECK_THROWS(LoadMesh({{"PEC", {{"Attributes", {7}}}}}, legacy, nc_path));
  CHECK_NOTHROW(LoadMesh(Impedance(false), none, nc_path));
  CHECK_NOTHROW(LoadMesh({{"PEC", {{"Attributes", {7}}}}}, none, nc_path));
  CHECK_NOTHROW(LoadMesh({{"PMC", {{"Attributes", {7}}}}}, off, nc_path));
  CHECK(LoadMesh(Impedance(true), none, nc_cut_path)->GetNV() == nv_nc_cut);
  CHECK(LoadMesh({{"PMC", {{"Attributes", {7}}}}}, none, nc_cut_path)->GetNV() ==
        nv_nc_cut);
  fs::remove(nc_path);
  fs::remove(nc_cut_path);
  fs::remove(mesh_path);
}

TEST_CASE("Default conductivity sheet warning", "[geodata][Serial]")
{
  // Unit cube hexahedral mesh with interior boundary elements (attribute 7) on z = 0.5,
  // spanning the domain (exterior boundary attributes 1-6). One exterior element and one
  // element on z = 0.5 have attribute 8, which is both exterior and interior.
  auto base = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON);
  std::vector<int> sheet_faces;
  mfem::Array<int> fv;
  for (int f = 0; f < base.GetNumFaces(); f++)
  {
    int e1, e2;
    base.GetFaceElements(f, &e1, &e2);
    base.GetFaceVertices(f, fv);
    bool on_sheet = (e2 >= 0);
    for (auto v : fv)
    {
      on_sheet = on_sheet && std::abs(base.GetVertex(v)[2] - 0.5) < 1.0e-12;
    }
    if (on_sheet)
    {
      sheet_faces.push_back(f);
    }
  }
  REQUIRE(sheet_faces.size() == 4);
  mfem::Mesh mesh(3, base.GetNV(), base.GetNE(),
                  base.GetNBE() + static_cast<int>(sheet_faces.size()));
  for (int v = 0; v < base.GetNV(); v++)
  {
    mesh.AddVertex(base.GetVertex(v));
  }
  for (int e = 0; e < base.GetNE(); e++)
  {
    mesh.AddElement(base.GetElement(e)->Duplicate(&mesh));
  }
  for (int be = 0; be < base.GetNBE(); be++)
  {
    auto *el = base.GetBdrElement(be)->Duplicate(&mesh);
    if (be == 0)
    {
      el->SetAttribute(8);
    }
    mesh.AddBdrElement(el);
  }
  for (auto f : sheet_faces)
  {
    auto *el = base.GetFace(f)->Duplicate(&mesh);
    el->SetAttribute(f == sheet_faces.back() ? 8 : 7);
    mesh.AddBdrElement(el);
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  mesh.SetAttributes();

  config::ConductivityData cond;
  cond.sigma = 5.8e7;
  cond.attributes = {7};
  auto Sheets = [&](const config::ConductivityData &data)
  {
    return mesh::GetDefaultConductivitySheetAttributes(
        mesh, std::vector<config::ConductivityData>{data});
  };

  // Interior boundary without "Crack", also with a "Thickness".
  CHECK(Sheets(cond) == std::vector<int>{7});
  config::ConductivityData film = cond;
  film.h = 1.0e-3;
  CHECK(Sheets(film) == std::vector<int>{7});

  // Exterior boundaries, and attributes with an exterior element, are not reported.
  cond.attributes = {1};
  CHECK(Sheets(cond).empty());
  cond.attributes = {8};
  CHECK(Sheets(cond).empty());
  cond.attributes = {1, 7, 8};
  CHECK(Sheets(cond) == std::vector<int>{7});
  cond.attributes = {7, 12};  // No boundary elements with attribute 12.
  CHECK(Sheets(cond) == std::vector<int>{7});
  cond.attributes = {7};

  // An explicit "Crack", either way, is an informed choice.
  config::ConductivityData given = cond;
  given.crack_given = true;
  CHECK(Sheets(given).empty());
  given.crack = true;
  CHECK(Sheets(given).empty());

  // Nothing is reported without conductivity boundaries.
  CHECK(mesh::GetDefaultConductivitySheetAttributes(mesh,
                                                    std::vector<config::ConductivityData>{})
            .empty());
}

// MFEM identifies a face by its three smallest vertex numbers. A cracked quadrilateral with
// three vertices on the crack front (at a corner of the crack) merges with its twin and the
// crack closes there, unless the front vertices are numbered last.
TEST_CASE("CrackQuadrilateralCorners", "[geodata][Serial]")
{
  // Hexahedral mesh with an interior boundary (attribute 7) made of the four quadrilaterals
  // of z = 0.5 around the center of the domain: each has three vertices on the crack front
  // and one (the center) in its interior.
  auto base = mfem::Mesh::MakeCartesian3D(4, 4, 2, mfem::Element::HEXAHEDRON);
  std::vector<int> sheet_faces;
  mfem::Array<int> fv;
  for (int f = 0; f < base.GetNumFaces(); f++)
  {
    int e1, e2;
    base.GetFaceElements(f, &e1, &e2);
    base.GetFaceVertices(f, fv);
    bool on_sheet = (e2 >= 0);
    for (auto v : fv)
    {
      const double *x = base.GetVertex(v);
      on_sheet = on_sheet && std::abs(x[2] - 0.5) < 1.0e-12 && x[0] > 0.25 - 1.0e-12 &&
                 x[0] < 0.75 + 1.0e-12 && x[1] > 0.25 - 1.0e-12 && x[1] < 0.75 + 1.0e-12;
    }
    if (on_sheet)
    {
      sheet_faces.push_back(f);
    }
  }
  REQUIRE(sheet_faces.size() == 4);
  mfem::Mesh mesh(3, base.GetNV(), base.GetNE(),
                  base.GetNBE() + static_cast<int>(sheet_faces.size()));
  for (int v = 0; v < base.GetNV(); v++)
  {
    mesh.AddVertex(base.GetVertex(v));
  }
  for (int e = 0; e < base.GetNE(); e++)
  {
    mesh.AddElement(base.GetElement(e)->Duplicate(&mesh));
  }
  for (int be = 0; be < base.GetNBE(); be++)
  {
    mesh.AddBdrElement(base.GetBdrElement(be)->Duplicate(&mesh));
  }
  for (auto f : sheet_faces)
  {
    auto *el = base.GetFace(f)->Duplicate(&mesh);
    el->SetAttribute(7);
    mesh.AddBdrElement(el);
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  mesh.SetAttributes();
  const auto mesh_path = fs::temp_directory_path() / "palace-test-crack-quad-corners.mesh";
  {
    std::ofstream fo(mesh_path);
    mesh.Print(fo);
  }

  // Crack the mesh along the interior boundary. Only the center vertex is duplicated.
  json config = {
      {"Problem", {{"Type", "Eigenmode"}, {"Verbose", 0}}},
      {"Model", {{"Mesh", mesh_path.string()}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
      {"Boundaries",
       {{"Impedance", {{{"Attributes", {7}}, {"Ls", 1.0e-12}, {"Crack", true}}}}}},
      {"Solver", {{"Eigenmode", {{"Target", 1.0}}}}}};
  IoData iodata(config, false);
  auto cracked = mesh::Load(iodata, Mpi::World());
  REQUIRE(cracked);
  CHECK(cracked->GetNV() == mesh.GetNV() + 1);

  // Every boundary element of the crack must border a single element, on its own side.
  // Without the renumbering, the quadrilateral whose center vertex has the largest number
  // of its four merges with its twin, and its 2 boundary elements lie on the same face.
  int n_crack = 0, n_merged = 0;
  for (int be = 0; be < cracked->GetNBE(); be++)
  {
    if (cracked->GetBdrAttribute(be) != 7)
    {
      continue;
    }
    int e1, e2;
    cracked->GetFaceElements(cracked->GetBdrElementFaceIndex(be), &e1, &e2);
    n_crack++;
    n_merged += (e2 >= 0);
  }
  CHECK(n_crack == 8);
  CHECK(n_merged == 0);
  fs::remove(mesh_path);
}

TEST_CASE("Default IOData", "[iodata][Serial]")
{
  config::MaterialData material;
  material.attributes = {1};

  config::PeriodicBoundaryData periodic;

  // Pull from the mfem externals data folder.
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON);
  mfem::ParMesh pmesh(Mpi::World(), mesh);

  MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC, pmesh);

  REQUIRE(mat_op.HasLossTangent() == false);
}

}  // namespace palace
