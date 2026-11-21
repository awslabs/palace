// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <fstream>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fixtures.hpp"
#include "test-helpers.hpp"

#include "utils/geodata.hpp"
#include "utils/geodata_impl.hpp"
#include "utils/iodata.hpp"

#include "fem/interpolator.hpp"
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
