// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "test-helpers.hpp"

#include "fem/fespace.hpp"
#include "utils/geodata.hpp"
#include "utils/geodata_impl.hpp"

#include "fem/interpolator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/edgedistance.hpp"
#include "utils/filesystem.hpp"
#include "utils/iodata.hpp"
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

TEST_CASE("TwoDimensionalDiagonalSquarePort", "[geodata][Serial]")
{
  std::vector<Eigen::Vector3d> vertices{
      {-0.19942181818181828, -0.5838274545454543, 0},
      {-0.19926108900667502, -0.5836667480978528, 0},
      {-0.19926061925486355, -0.5838279242308415, -8.56131891164518e-17},
      {-0.19926014950305207, -0.5839891003638303, 0},
      {-0.19910035983153176, -0.5835060416502512, 0},
      {-0.1990998900797203, -0.5836672177832399, -8.56131891164518e-17},
      {-0.1990994203279088, -0.5838283939162288, -1.3698110258632289e-16},
      {-0.19909895057609733, -0.5839895700492176, -8.56131891164518e-17},
      {-0.19909848082428586, -0.5841507461822064, 0},
      {-0.1989396306563885, -0.5833453352026496, 0},
      {-0.19893916090457703, -0.5835065113356385, -8.56131891164518e-17},
      {-0.19893869115276555, -0.5836676874686273, -1.3698110258632289e-16},
      {-0.19893822140095407, -0.583828863601616, -1.5410374040961326e-16},
      {-0.1989377516491426, -0.583990039734605, -1.3698110258632289e-16},
      {-0.19893728189733112, -0.5841512158675937, -8.56131891164518e-17},
      {-0.19893681214551964, -0.5843123920005826, 0},
      {-0.19877890148124525, -0.5831846287550482, 0},
      {-0.19877843172943377, -0.5833458048880369, -8.56131891164518e-17},
      {-0.1987779619776223, -0.5835069810210257, -1.3698110258632289e-16},
      {-0.19877749222581081, -0.5836681571540145, -1.5410374040961326e-16},
      {-0.19877702247399934, -0.5838293332870034, -1.3698110258632289e-16},
      {-0.19877655272218786, -0.5839905094199922, -1.5410374040961326e-16},
      {-0.19877608297037638, -0.584151685552981, -1.3698110258632289e-16},
      {-0.1987756132185649, -0.5843128616859699, -8.56131891164518e-17},
      {-0.19877514346675343, -0.5844740378189587, 0},
      {-0.19861817611093396, -0.5830239261117406, 0},
      {-0.1986153576000651, -0.5846338010914915, 0},
      {-0.19845746797462832, -0.583506036157253, -8.56131891164518e-17},
      {-0.19845745074062265, -0.582863223468433, 0},
      {-0.19845699822281687, -0.5836672122902418, -1.3698110258632289e-16},
      {-0.19845652969381722, -0.5831855725288573, -8.56131891164518e-17},
      {-0.19845652847100537, -0.5838283884232307, -1.5410374040961326e-16},
      {-0.1984560587191939, -0.5839895645562194, -1.3698110258632289e-16},
      {-0.1984555889673824, -0.5841507406892084, -1.5410374040961326e-16},
      {-0.19845557173337677, -0.5847935643640245, 0},
      {-0.19845511921557096, -0.5843119168221971, -1.3698110258632289e-16},
      {-0.19845464946375949, -0.5844730929551859, -8.56131891164518e-17},
      {-0.19829672537031137, -0.5827025208251255, 0},
      {-0.19829580432350594, -0.5830248698855497, -8.56131891164518e-17},
      {-0.19829578586668845, -0.5849533276365575, 0},
      {-0.19829486359707116, -0.5846328562277189, -8.56131891164518e-17},
      {-0.1981360344680114, -0.583827443559458, -1.3698110258632289e-16},
      {-0.1981360000000001, -0.5851130909090906, 0},
      {-0.1981360000000001, -0.5825418181818179, 0},
      {-0.19813556471619992, -0.5839886196924468, -1.5410374040961324e-16},
      {-0.19813509618720032, -0.5835069799310624, -1.3698110258632289e-16},
      {-0.19813509496438844, -0.5841497958254356, -1.3698110258632289e-16},
      {-0.19813507895319465, -0.5828641672422421, -8.56131891164518e-17},
      {-0.19813507773038283, -0.5847926195002519, -8.56131891164518e-17},
      {-0.198134625212577, -0.5843109719584244, -1.5410374040961324e-16},
      {-0.19813415790638922, -0.5831865163026665, -1.3698110258632289e-16},
      {-0.19813415546076552, -0.5844721480914132, -1.3698110258632289e-16},
      {-0.19797529186369447, -0.5849523827727849, -8.56131891164518e-17},
      {-0.1979743695940772, -0.5846319113639462, -1.3698110258632289e-16},
      {-0.19797435358288337, -0.5827034645989346, -8.56131891164518e-17},
      {-0.19797343253607796, -0.5830258136593589, -1.3698110258632289e-16},
      {-0.19781460096139447, -0.584148850961663, -1.5410374040961326e-16},
      {-0.19781458372738883, -0.5847916746364793, -1.3698110258632289e-16},
      {-0.197814131209583, -0.5843100270946517, -1.3698110258632289e-16},
      {-0.19781366268058337, -0.5838283873332673, -1.5410374040961326e-16},
      {-0.19781366145777152, -0.5844712032276406, -1.5410374040961326e-16},
      {-0.19781272439977232, -0.5835079237048715, -1.5410374040961324e-16},
      {-0.19781270716576665, -0.5828651110160513, -1.3698110258632289e-16},
      {-0.19781178611896122, -0.5831874600764757, -1.5410374040961326e-16},
      {-0.1976538755910832, -0.5846309665001737, -1.5410374040961326e-16},
      {-0.1976510607486499, -0.583026757433168, -1.5410374040961326e-16},
      {-0.19749316745477752, -0.584470258363868, -1.3698110258632289e-16},
      {-0.19749222917396644, -0.5841497947354722, -1.3698110258632289e-16},
      {-0.19749129089315537, -0.5838293311070765, -1.3698110258632289e-16},
      {-0.19749035261234427, -0.5835088674786807, -1.3698110258632289e-16},
      {-0.1974894143315332, -0.5831884038502848, -1.3698110258632289e-16},
      {-0.1973324665001741, -0.5843095574092645, -1.3698110258632289e-16},
      {-0.19733152821936303, -0.5839890937808687, -1.3698110258632289e-16},
      {-0.19733058993855196, -0.583668630152473, 1.3698110258632289e-16},
      {-0.19732965165774086, -0.5833481665240772, -1.3698110258632289e-16},
      {-0.19717176554557067, -0.5841488564546611, -1.3698110258632289e-16},
      {-0.1971708272647596, -0.5838283928262653, -1.3698110258632289e-16},
      {-0.1971698889839485, -0.5835079291978695, -1.3698110258632289e-16},
      {-0.1970110645909672, -0.5839881555000577, -1.3698110258632289e-16},
      {-0.19701012631015613, -0.5836676918716619, -1.3698110258632289e-16},
      {-0.19685036363636374, -0.5838274545454543, -1.3698110258632289e-16}};

  auto comm = Mpi::World();
  auto box = mesh::BoundingBoxFromPointCloud(comm, vertices, 0);

  // True box at 45 degrees
  auto invsqrt2 = 1.0 / std::sqrt(2);
  std::array<double, 3> ax0{invsqrt2, -invsqrt2, 0.0}, ax1{invsqrt2, invsqrt2, 0.0};

  // Find the bounding points from knowing its at 45.
  auto inf = std::numeric_limits<double>::infinity();
  double min_x = inf, min_y = inf, max_x = -inf, max_y = -inf;
  for (const auto &v : vertices)
  {
    min_x = min_x > v(0) ? v(0) : min_x;
    min_y = min_y > v(1) ? v(1) : min_y;
    max_x = max_x < v(0) ? v(0) : max_x;
    max_y = max_y < v(1) ? v(1) : max_y;
  }
  auto length_x = (max_x - min_x) * invsqrt2;
  auto length_y = (max_y - min_y) * invsqrt2;

  CHECK_THAT(length_x, WithinAbs(length_y, 1e-6));

  auto length = (length_x + length_y) / 2;
  auto lengths = box.Lengths();
  CHECK_THAT(lengths[0], WithinAbs(length, 1e-6));
  CHECK_THAT(lengths[1], WithinAbs(length, 1e-6));
  CHECK_THAT(lengths[2], WithinRel(0.0));

  auto normals = box.Normals();
  CHECK_THAT(normals(0, 0), WithinAbs(ax0[0], 1e-4));
  CHECK_THAT(normals(1, 0), WithinAbs(ax0[1], 1e-4));
  CHECK_THAT(normals(2, 0), WithinAbs(ax0[2], 1e-4));
  CHECK_THAT(normals(0, 1), WithinAbs(ax1[0], 1e-4));
  CHECK_THAT(normals(1, 1), WithinAbs(ax1[1], 1e-4));
  CHECK_THAT(normals(2, 1), WithinAbs(ax1[2], 1e-4));
  CHECK(box.planar);
}

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
      (*single_tet.GetNodes())(i) += 0.05 * (1.0 + (double)rand() / RAND_MAX);
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
    for (int i = 0; i < tet_vals.Size(); i++)
    {
      CHECK_THAT(tet_vals(i), WithinAbs(hex_vals(i), 1e-9));
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
  const std::unordered_set<int> cracked_attributes = {1};

  CHECK_THROWS(SurfacePostOperator(postpro, ProblemType::ELECTROSTATIC, mat_op, h1_fespace,
                                   nd_fespace, &cracked_attributes));
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

  MetalSurfaceExtraction surface;
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
  geometry = ExtractMetalEdgeGeometry(mesh, boundaries);
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

  auto Extract = [&](int elements, bool high_order, bool rounded)
  {
    auto mesh = RoundedIslandSheetHexMesh(elements, high_order, rounded);
    MetalSurfaceExtraction surface;
    surface.retain_faces = true;
    auto geometry = ExtractMetalEdgeGeometry(*mesh, boundaries, surface);
    return std::pair{std::move(mesh), std::move(geometry)};
  };
  auto [coarse_mesh, coarse] = Extract(8, true, true);
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
  CHECK(coarse.physical_chains == 1);
  CHECK(fine.physical_chains == 1);
  CHECK(CountPhysicalVertexTypes(coarse, MetalEdgeVertexType::CORNER) == 0);
  CHECK(CountPhysicalVertexTypes(fine, MetalEdgeVertexType::CORNER) == 0);
  CHECK(coarse.segments.size() > fine.segments.size());

  const double exact_perimeter =
      4.0 * (2.0 * 0.25 - 2.0 * 0.125) + 2.0 * std::acos(-1.0) * 0.125;
  CHECK_THAT(PhysicalLength(coarse), WithinAbs(exact_perimeter, 2.0e-3));
  CHECK_THAT(PhysicalLength(fine), WithinAbs(exact_perimeter, 6.0e-3));
  CHECK_THAT(PhysicalLength(coarse), WithinAbs(PhysicalLength(fine), 6.0e-3));

  REQUIRE(!coarse.surface_faces.empty());
  REQUIRE(!fine.surface_faces.empty());
  CHECK(std::all_of(coarse.surface_faces.begin(), coarse.surface_faces.end(),
                    [](const auto &face) { return face.component == 0; }));
  CHECK(std::all_of(fine.surface_faces.begin(), fine.surface_faces.end(),
                    [](const auto &face)
                    { return face.component == 0 && face.vertices.size() == 4; }));
  CHECK(std::any_of(coarse.surface_faces.begin(), coarse.surface_faces.end(),
                    [](const auto &face) { return face.vertices.size() > 8; }));

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
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/cpw3d_surface/cpw3d_surface_validation_thin.json";
  IoData iodata(config_path.c_str(), false);
  iodata.model.mesh = fs::path(PALACE_TEST_DATA_DIR) / "mesh/cpw3d-surface-nc.msh";
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

  const auto geometry = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries);
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
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/cpw3d_surface/cpw3d_surface_validation_thin.json";
  IoData iodata(config_path.c_str(), false);
  iodata.model.mesh = fs::path(PALACE_TEST_DATA_DIR) / "mesh/cpw3d-surface-nc.msh";
  iodata.model.refinement.max_it = 1;
  auto mesh = mesh::ReadMesh(iodata, Mpi::World());

  const auto coarse = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries);
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
  REQUIRE_FALSE(adjacent_elements.empty());
  mfem::Array<int> marked_elements(static_cast<int>(adjacent_elements.size()));
  std::copy(adjacent_elements.begin(), adjacent_elements.end(), marked_elements.begin());
  mesh->GeneralRefinement(marked_elements);

  const auto refined = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries);
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

TEST_CASE("Automatic metal edge extraction on 3D transmon",
          "[geodata][metaledge][Serial][Parallel]")
{
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/transmon/transmon_surface_coarse.json";
  IoData iodata(config_path.c_str(), false);
  iodata.model.mesh = config_path.parent_path() / iodata.model.mesh;
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

  const auto geometry = ExtractMetalEdgeGeometry(*mesh, iodata.boundaries);
  REQUIRE_FALSE(geometry.Empty());
  int physical_segments = 0;
  int truncation_segments = 0;
  int corners = 0;
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
          physical_segments, truncation_segments, corners, sa_segments,
          truncation_attributes);
  CHECK(geometry.components == 6);
  CHECK(geometry.physical_components == 5);
  CHECK(physical_segments + truncation_segments ==
        static_cast<int>(geometry.segments.size()));
  CHECK(physical_segments > 0);
  CHECK(geometry.physical_chains == 56);
  CHECK(truncation_segments > 0);
  CHECK(truncation_attributes == std::set<int>{3});
  CHECK(corners == 56);
  CHECK(sa_segments > 0);

  const auto sa_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry, 1, InterfaceDielectric::SA);
  const auto ms_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry, 2, InterfaceDielectric::MS);
  const auto ma_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry, 3, InterfaceDielectric::MA);
  CHECK(sa_indices.size() == 3082);
  CHECK(ms_indices.size() == 3088);
  CHECK(ma_indices.size() == 3088);

  const auto process_normals =
      BuildMetalEdgeProcessNormals(*mesh, geometry, ms_indices, [](int material_attribute)
                                   { return material_attribute == 2 ? 1.0 : 0.0; });
  REQUIRE(process_normals.size() == ms_indices.size());
  for (const auto &normal : process_normals)
  {
    CHECK_THAT(normal[0], WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(normal[1], WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(normal[2], WithinAbs(1.0, 1.0e-12));
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

TEST_CASE("Boundary edge extraction on NC 3D CPW mesh", "[geodata][Serial][Parallel]")
{
  const auto config_path = fs::path(__FILE__).parent_path().parent_path().parent_path() /
                           "examples/cpw3d_surface/cpw3d_surface_validation_thin.json";
  auto ExtractEdges = [&](bool use_amr)
  {
    IoData iodata(config_path.c_str(), false);
    iodata.model.mesh = fs::path(PALACE_TEST_DATA_DIR) / "mesh/cpw3d-surface-nc.msh";
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
