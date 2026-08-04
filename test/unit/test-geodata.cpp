// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "test-helpers.hpp"

#include "utils/geodata.hpp"
#include "utils/geodata_impl.hpp"
#include "utils/iodata.hpp"

#include "fem/interpolator.hpp"
#include "fem/singulardofs.hpp"
#include "fem/singularfeatures.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/filesystem.hpp"

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

mfem::Mesh RectangleInternalSheetMesh()
{
  mfem::Mesh mesh(3, 6, 4, 2, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(1.0, 1.0, 0.0);
  mesh.AddVertex(0.5, 0.5, 1.0);
  mesh.AddVertex(0.5, 0.5, -1.0);
  mesh.AddTet(0, 1, 2, 4, 1);
  mesh.AddTet(0, 2, 1, 5, 1);
  mesh.AddTet(1, 3, 2, 4, 1);
  mesh.AddTet(1, 2, 3, 5, 1);
  mesh.AddBdrTriangle(0, 1, 2, 7);
  mesh.AddBdrTriangle(1, 3, 2, 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh ThreeInternalSheetMesh()
{
  mfem::Mesh mesh(3, 15, 6, 3, 3);
  for (int component = 0; component < 3; component++)
  {
    const double x = 3.0 * component;
    mesh.AddVertex(x, 0.0, 0.0);
    mesh.AddVertex(x + 1.0, 0.0, 0.0);
    mesh.AddVertex(x, 1.0, 0.0);
    mesh.AddVertex(x + 0.3, 0.3, 1.0);
    mesh.AddVertex(x + 0.3, 0.3, -1.0);
    const int v = 5 * component;
    mesh.AddTet(v, v + 1, v + 2, v + 3, 1);
    mesh.AddTet(v, v + 2, v + 1, v + 4, 1);
    mesh.AddBdrTriangle(v, v + 1, v + 2, 7 + component);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
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

TEST_CASE("Singular internal sheets remain conforming during mesh loading",
          "[geodata][singularfeatures][Serial]")
{
  const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
  const fs::path mesh_path =
      fs::temp_directory_path() / fmt::format("palace-singular-sheet-{}.mesh", nonce);
  {
    auto mesh = RectangleInternalSheetMesh();
    std::ofstream output(mesh_path);
    REQUIRE(output.good());
    output.precision(std::numeric_limits<double>::max_digits10);
    mesh.Print(output);
  }

  auto MakeConfig = [&mesh_path](bool singular)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", mesh_path.string()}, {"AddInterfaceBoundaryElements", false}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries", {{"PEC", {{"Attributes", {7}}}}}},
        {"Solver", json::object()}};
    if (singular)
    {
      config["Solver"]["SingularElements"] = {{"Attributes", {7}}};
    }
    return config;
  };

  {
    IoData iodata(MakeConfig(true), false);
    auto mesh = mesh::Load(iodata, Mpi::World());
    REQUIRE(mesh);
    CHECK(mesh->GetNV() == 6);
    CHECK(mesh->GetNBE() == 2);
    CHECK(iodata.boundaries.cracked_attributes.count(7) == 0);
    const auto features = fem::singular::ExtractSerialSheetFeatures(*mesh, {7});
    CHECK_FALSE(features.Empty());
  }
  {
    IoData iodata(MakeConfig(false), false);
    auto mesh = mesh::Load(iodata, Mpi::World());
    REQUIRE(mesh);
    CHECK(mesh->GetNV() > 6);
    CHECK(mesh->GetNBE() == 4);
    CHECK(iodata.boundaries.cracked_attributes.count(7) == 1);
  }
  CHECK(fs::remove(mesh_path));
}

TEST_CASE("Mesh loading cracks only ordinary internal sheets",
          "[geodata][singularfeatures][Serial]")
{
  const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
  const fs::path mesh_path =
      fs::temp_directory_path() / fmt::format("palace-selective-crack-{}.mesh", nonce);
  {
    auto mesh = ThreeInternalSheetMesh();
    std::ofstream output(mesh_path);
    REQUIRE(output.good());
    output.precision(std::numeric_limits<double>::max_digits10);
    mesh.Print(output);
  }

  json config = {
      {"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", mesh_path.string()}, {"AddInterfaceBoundaryElements", false}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
      {"Boundaries",
       {{"PEC", {{"Attributes", {7, 8}}}},
        {"LumpedPort",
         {{{"Index", 1},
           {"Attributes", {9}},
           {"R", 50.0},
           {"Direction", "+Z"},
           {"Excitation", true}}}}}},
      {"Solver", {{"SingularElements", {{"Attributes", {7}}}}}}};
  IoData iodata(config, false);
  auto mesh = mesh::Load(iodata, Mpi::World());
  REQUIRE(mesh);
  CAPTURE(mesh->GetNV(), mesh->GetNBE(), iodata.boundaries.cracked_attributes);
  CHECK(iodata.boundaries.cracked_attributes.size() == 1);
  CHECK(iodata.boundaries.cracked_attributes.count(8) == 1);
  CHECK(iodata.boundaries.cracked_attributes.count(7) == 0);
  CHECK(iodata.boundaries.cracked_attributes.count(9) == 0);
  CHECK(mesh->GetNV() > 15);
  CHECK(mesh->GetNBE() == 4);

  // The singular sheet remains two-sided and conforming after the unrelated sheet cracks.
  const auto features = fem::singular::ExtractSerialSheetFeatures(*mesh, {7});
  CHECK_FALSE(features.Empty());
  REQUIRE(features.sheet_faces.size() == 1);
  CHECK(features.sheet_faces.front().boundary_attribute == 7);
  CHECK(fs::remove(mesh_path));
}

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

TEST_CASE("Mesh partition transports exact source serial entity IDs", "[geodata][Parallel]")
{
  if (Mpi::Size(Mpi::World()) == 1)
  {
    SUCCEED("Source-entity transport is exercised by the [Parallel] test run.");
    return;
  }
  REQUIRE(Mpi::Size(Mpi::World()) == 2);

  std::unique_ptr<mfem::Mesh> serial_mesh;
  if (Mpi::Root(Mpi::World()))
  {
    serial_mesh = std::make_unique<mfem::Mesh>(3, 6, 4, 0, 3);
    serial_mesh->AddVertex(0.0, 0.0, 0.0);
    serial_mesh->AddVertex(1.0, 0.0, 0.0);
    serial_mesh->AddVertex(0.0, 1.0, 0.0);
    serial_mesh->AddVertex(1.0, 1.0, 0.0);
    serial_mesh->AddVertex(0.5, 0.5, 1.0);
    serial_mesh->AddVertex(0.5, 0.5, -1.0);
    serial_mesh->AddTet(0, 1, 2, 4, 1);
    serial_mesh->AddTet(0, 2, 1, 5, 1);
    serial_mesh->AddTet(1, 3, 2, 4, 1);
    serial_mesh->AddTet(1, 2, 3, 5, 1);
    serial_mesh->FinalizeTopology();
    serial_mesh->Finalize(true, false);
  }

  IoData iodata(Units(1.0, 1.0));
  mesh::PartitionMetadata metadata;
  auto parallel_mesh =
      mesh::Partition(iodata, std::move(serial_mesh), Mpi::World(), &metadata);

  REQUIRE(metadata.source_element_ids.size() ==
          static_cast<std::size_t>(parallel_mesh->GetNE()));
  REQUIRE(metadata.source_vertex_ids.size() ==
          static_cast<std::size_t>(parallel_mesh->GetNV()));
  CHECK(std::is_sorted(metadata.source_element_ids.begin(),
                       metadata.source_element_ids.end()));
  CHECK(
      std::is_sorted(metadata.source_vertex_ids.begin(), metadata.source_vertex_ids.end()));

  const std::array<std::array<std::int64_t, 4>, 4> source_vertices{
      std::array<std::int64_t, 4>{0, 1, 2, 4}, std::array<std::int64_t, 4>{0, 1, 2, 5},
      std::array<std::int64_t, 4>{1, 2, 3, 4}, std::array<std::int64_t, 4>{1, 2, 3, 5}};
  std::array<int, 4> source_element_owners{};
  for (int local_element = 0; local_element < parallel_mesh->GetNE(); local_element++)
  {
    const auto source_element = metadata.source_element_ids[local_element];
    REQUIRE(source_element >= 0);
    REQUIRE(source_element < static_cast<std::int64_t>(source_vertices.size()));
    source_element_owners[source_element]++;

    const auto &element = *parallel_mesh->GetElement(local_element);
    std::array<std::int64_t, 4> mapped_vertices;
    for (int i = 0; i < 4; i++)
    {
      mapped_vertices[i] = metadata.source_vertex_ids[element.GetVertices()[i]];
    }
    std::sort(mapped_vertices.begin(), mapped_vertices.end());
    CHECK(mapped_vertices == source_vertices[source_element]);
  }
  Mpi::GlobalSum(static_cast<int>(source_element_owners.size()),
                 source_element_owners.data(), Mpi::World());
  for (int owners : source_element_owners)
  {
    CHECK(owners == 1);
  }
}

TEST_CASE("Nonconforming mesh partition retains source vertex identities",
          "[geodata][singularelements][Parallel]")
{
  if (Mpi::Size(Mpi::World()) == 1)
  {
    SUCCEED("Nonconforming source identities are exercised by the parallel test run.");
    return;
  }
  REQUIRE(Mpi::Size(Mpi::World()) == 2);

  const std::array<std::array<double, 3>, 6> source_coordinates{
      std::array<double, 3>{0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0},
      {0.5, 0.5, 1.0},
      {0.5, 0.5, -1.0}};
  std::unique_ptr<mfem::Mesh> serial_mesh;
  if (Mpi::Root(Mpi::World()))
  {
    serial_mesh = std::make_unique<mfem::Mesh>(3, 6, 4, 0, 3);
    for (const auto &coordinate : source_coordinates)
    {
      serial_mesh->AddVertex(coordinate.data());
    }
    serial_mesh->AddTet(0, 1, 2, 4, 1);
    serial_mesh->AddTet(0, 2, 1, 5, 1);
    serial_mesh->AddTet(1, 3, 2, 4, 1);
    serial_mesh->AddTet(1, 2, 3, 5, 1);
    serial_mesh->FinalizeTopology();
    serial_mesh->Finalize(true, false);
  }

  IoData iodata(Units(1.0, 1.0));
  iodata.model.refinement.max_it = 1;
  iodata.model.refinement.nonconformal = true;
  mesh::PartitionMetadata metadata;
  auto parallel_mesh =
      mesh::Partition(iodata, std::move(serial_mesh), Mpi::World(), &metadata);
  REQUIRE(parallel_mesh->Nonconforming());
  REQUIRE(metadata.source_vertex_ids.size() ==
          static_cast<std::size_t>(parallel_mesh->GetNV()));

  std::array<int, source_coordinates.size()> found{};
  for (int vertex = 0; vertex < parallel_mesh->GetNV(); vertex++)
  {
    const auto source = metadata.source_vertex_ids[vertex];
    REQUIRE(source >= 0);
    REQUIRE(source < static_cast<std::int64_t>(source_coordinates.size()));
    found[source]++;
    const auto *coordinate = parallel_mesh->GetVertex(vertex);
    for (int d = 0; d < 3; d++)
    {
      CHECK(coordinate[d] == source_coordinates[source][d]);
    }
  }
  Mpi::GlobalSum(static_cast<int>(found.size()), found.data(), Mpi::World());
  for (int count : found)
  {
    CHECK(count > 0);
  }
}

TEST_CASE("Mesh rebalancing transports persistent source identities",
          "[geodata][singularelements][Parallel]")
{
  if (Mpi::Size(Mpi::World()) == 1)
  {
    SUCCEED("Rebalancing identity transport is exercised by the parallel test run.");
    return;
  }
  REQUIRE(Mpi::Size(Mpi::World()) == 2);

  const bool nonconforming = GENERATE(false, true);
  CAPTURE(nonconforming);
  auto serial_mesh =
      mfem::Mesh::MakeCartesian2D(4, 2, mfem::Element::TRIANGLE, true, 2.0, 1.0);
  if (nonconforming)
  {
    serial_mesh.EnsureNCMesh(true);
  }
  const int serial_vertices = serial_mesh.GetNV();
  const int serial_elements = serial_mesh.GetNE();
  std::vector<std::array<double, 2>> serial_coordinates(serial_vertices);
  for (int vertex = 0; vertex < serial_vertices; vertex++)
  {
    const auto *coordinate = serial_mesh.GetVertex(vertex);
    serial_coordinates[vertex] = {coordinate[0], coordinate[1]};
  }

  std::vector<int> initial_partition(serial_elements, 0);
  auto parallel_mesh =
      std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh, initial_partition.data());
  REQUIRE(parallel_mesh->GetNE() == (Mpi::Root(Mpi::World()) ? serial_elements : 0));

  mesh::PartitionMetadata metadata;
  if (nonconforming)
  {
    metadata.source_vertex_ids.assign(parallel_mesh->GetNV(), -1);
    const auto &ncmesh = *parallel_mesh->ncmesh;
    for (int node = 0; node < ncmesh.GetNumNodes(); node++)
    {
      const auto &record = ncmesh.GetNode(node);
      if (record.HasVertex() && record.vert_index >= 0 &&
          record.vert_index < parallel_mesh->GetNV())
      {
        metadata.source_vertex_ids[record.vert_index] = node;
      }
    }
  }
  else
  {
    metadata.source_vertex_ids = fem::singular::MapPartitionedSerialVertexIds(
        serial_mesh, *parallel_mesh, initial_partition.data());
  }
  metadata.source_element_ids.resize(parallel_mesh->GetNE());
  std::iota(metadata.source_element_ids.begin(), metadata.source_element_ids.end(), 0);

  IoData iodata(Units(1.0, 1.0));
  iodata.model.refinement.maximum_imbalance = 1.01;
  const double initial_imbalance = mesh::RebalanceMesh(iodata, parallel_mesh, &metadata);
  CHECK(std::isinf(initial_imbalance));
  REQUIRE(parallel_mesh->GetNE() > 0);
  REQUIRE(metadata.source_vertex_ids.size() ==
          static_cast<std::size_t>(parallel_mesh->GetNV()));
  REQUIRE(metadata.source_element_ids.size() ==
          static_cast<std::size_t>(parallel_mesh->GetNE()));

  int min_elements = parallel_mesh->GetNE();
  int max_elements = parallel_mesh->GetNE();
  Mpi::GlobalMin(1, &min_elements, Mpi::World());
  Mpi::GlobalMax(1, &max_elements, Mpi::World());
  CHECK(min_elements > 0);
  CHECK(double(max_elements) / min_elements <= 1.01);

  std::map<std::array<double, 2>, std::int64_t> local_id_by_coordinate;
  for (int vertex = 0; vertex < parallel_mesh->GetNV(); vertex++)
  {
    const auto *coordinate = parallel_mesh->GetVertex(vertex);
    const std::array<double, 2> point{coordinate[0], coordinate[1]};
    const auto source = metadata.source_vertex_ids[vertex];
    REQUIRE(source >= 0);
    CHECK(local_id_by_coordinate.emplace(point, source).second);
    if (!nonconforming)
    {
      REQUIRE(source < serial_vertices);
      CHECK(point == serial_coordinates[source]);
    }
  }

  const int local_vertices = parallel_mesh->GetNV();
  std::vector<int> counts(Mpi::Size(Mpi::World()));
  Mpi::Allgather(1, &local_vertices, counts.data(), Mpi::World());
  std::vector<int> offsets(counts.size());
  std::partial_sum(counts.begin(), counts.end() - 1, offsets.begin() + 1);
  const int total_vertices = offsets.back() + counts.back();
  std::vector<double> local_coordinates(2 * local_vertices);
  for (int vertex = 0; vertex < local_vertices; vertex++)
  {
    const auto *coordinate = parallel_mesh->GetVertex(vertex);
    local_coordinates[2 * vertex] = coordinate[0];
    local_coordinates[2 * vertex + 1] = coordinate[1];
  }
  std::vector<int> coordinate_counts(counts), coordinate_offsets(offsets);
  for (int &count : coordinate_counts)
  {
    count *= 2;
  }
  for (int &offset : coordinate_offsets)
  {
    offset *= 2;
  }
  std::vector<double> global_coordinates(2 * total_vertices);
  std::vector<std::int64_t> global_ids(total_vertices);
  Mpi::Allgatherv(2 * local_vertices, local_coordinates.data(), global_coordinates.data(),
                  coordinate_counts.data(), coordinate_offsets.data(), Mpi::World());
  Mpi::Allgatherv(local_vertices, metadata.source_vertex_ids.data(), global_ids.data(),
                  counts.data(), offsets.data(), Mpi::World());

  std::map<std::array<double, 2>, std::int64_t> global_id_by_coordinate;
  int shared_vertices = 0;
  for (int vertex = 0; vertex < total_vertices; vertex++)
  {
    const std::array<double, 2> point{global_coordinates[2 * vertex],
                                      global_coordinates[2 * vertex + 1]};
    const auto [record, inserted] =
        global_id_by_coordinate.emplace(point, global_ids[vertex]);
    if (!inserted)
    {
      shared_vertices++;
      CHECK(record->second == global_ids[vertex]);
    }
  }
  CHECK(shared_vertices > 0);
  CHECK(global_id_by_coordinate.size() == static_cast<std::size_t>(serial_vertices));
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
